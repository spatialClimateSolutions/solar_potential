# ===================================================================
# MULTIFAMILY HOUSING ANALYSIS SCRIPT
# ===================================================================
#
# PURPOSE:
# This script processes geospatial data to identify potential development
# sites (Ground Layer) and multifamily building footprints (Rooftop Layer).
# It performs the following major steps:
#   1.  Loads and filters housing data to isolate multifamily properties.
#   2.  Generates a "Ground Layer" of potentially developable land by
#       analyzing impervious surfaces and excluding existing buildings,
#       schools, and public lands.
#   3.  Generates a "Rooftop Layer" by identifying building footprints
#       associated with multifamily housing units.
#   4.  Refines and cleans both layers through clustering, simplification,
#       and area-based filtering.
#   5.  Merges the ground layer with an additional impervious surface dataset
#       for further refinement.
#   6.  Generates a summary CSV file aggregating the results by census tract.
#   7.  Creates interactive HTML maps to visualize the results.
#
# INSTRUCTIONS:
# - Ensure all input files are available at the paths specified in the
#   CONFIGURATION section.
# - The script is designed to be run sequentially from top to bottom.
#
# ===================================================================

# CELL 1: IMPORTS
# ===================================================================
# Import all necessary libraries for the analysis.
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from shapely.geometry import shape, mapping
import os
import numpy as np
from scipy.ndimage import binary_opening
from sklearn.cluster import DBSCAN
import plotly.express as px
import gc

# ===================================================================
# CELL 2: CONFIGURATION
# ===================================================================
# Update all file paths in this section before running the script.
print("--- Step 1: Loading Configuration ---")

# --- General ---
CRS = "EPSG:3310"  # California Albers (projected CRS for accurate measurements)

# --- Input Files ---
HOUSING_CSV_PATH = "housing_classified_zillow1.csv"
IMPERVIOUS_RASTER_PATH = "Annual_NLCD_FctImp_2024_CU_C1V1.tif"
ALL_BUILDINGS_PATH = "california_building_footprints_complete.gpkg"
SCHOOLS_GDB_PATH = "CSCD_2021.gdb"
PUBLIC_LANDS_CPAD_PATH = "CPAD_2024b_Units.shp"
PUBLIC_LANDS_CCED_PATH = "CCED_2024b_Release.shp"
CENSUS_TRACTS_PATH = "tl_2024_06_tract.shp"
IMPERVIOUS_SURFACE_PATH_2 = "LC24_EVT_250_filtered_flat_points_OPTIMIZED.tif"

# A list of all polygon-based school layers to load and merge
SCHOOLS_LAYER_NAMES = [
    "University_lands", "Schools_Current_Stacked", "School_Property",
    "NonPublicK12School_EducationOwnedLands", "CA_Community_Colleges"
]

# --- Output Directories and Files ---
OUTPUT_DIR_GROUND = "output_ground_layer"
OUTPUT_DIR_ROOFTOP = "output_rooftop_layer"
OUTPUT_DIR_MAPS = "interactive_plotly_maps"

# Intermediate and Final File Paths
FINAL_GROUND_LAYER_PATH = os.path.join(OUTPUT_DIR_GROUND, "final_multifamily_ground_layer.gpkg")
FINAL_ROOFTOP_LAYER_PATH = os.path.join(OUTPUT_DIR_ROOFTOP, "final_multifamily_rooftop_layer.gpkg")
CLUSTERED_ROOFTOP_PATH_SIMPLIFIED = os.path.join(OUTPUT_DIR_ROOFTOP, "final_multifamily_rooftop_layer_clustered_simplified.gpkg")
FINAL_CLEANED_GROUND_LAYER = os.path.join(OUTPUT_DIR_GROUND, "final_multifamily_ground_layer_cleaned.gpkg")
OUTPUT_SUMMARY_CSV = "final_summary_by_tract.csv"
OUTPUT_MERGED_IMPERVIOUS_PATH = os.path.join(OUTPUT_DIR_GROUND, 'final_ground_layer_with_impervious_surfaces.gpkg')
FINAL_GROUND_LAYER_OPTIMIZED = os.path.join(OUTPUT_DIR_GROUND, "final_multifamily_ground_layer_optimized.gpkg")

# Create output directories if they don't exist
for path in [OUTPUT_DIR_GROUND, OUTPUT_DIR_ROOFTOP, OUTPUT_DIR_MAPS]:
    if not os.path.exists(path):
        os.makedirs(path)

print("Configuration loaded and output directories are ready.")

# ===================================================================
# CELL 3: LOAD AND PREPARE MULTIFAMILY HOUSING DATA
# ===================================================================
print("\n--- Step 2: Preparing Multifamily Data ---")
print(f"Loading housing data from {HOUSING_CSV_PATH}...")

housing_df = pd.read_csv(HOUSING_CSV_PATH)
multifamily_df = housing_df[housing_df['single_or_multi'].str.strip() == 'Multi'].copy()
print(f"Filtered to {len(multifamily_df)} 'Multi' family properties.")

# --- Data Cleaning: Ensure lon and lat are numeric ---
multifamily_df['lon'] = pd.to_numeric(multifamily_df['lon'], errors='coerce')
multifamily_df['lat'] = pd.to_numeric(multifamily_df['lat'], errors='coerce')
original_count = len(multifamily_df)
multifamily_df.dropna(subset=['lon', 'lat'], inplace=True)
cleaned_count = len(multifamily_df)
if original_count > cleaned_count:
    print(f"Removed {original_count - cleaned_count} rows with invalid coordinates.")

# Check if 'ID' column exists, if not, create a placeholder
if 'ID' not in multifamily_df.columns:
    multifamily_df['ID'] = range(len(multifamily_df))
    print("Warning: 'ID' column not found. Created a placeholder ID.")

multifamily_points = gpd.GeoDataFrame(
    multifamily_df,
    geometry=gpd.points_from_xy(multifamily_df.lon, multifamily_df.lat),
    crs="EPSG:4326"
).to_crs(CRS)

# ===================================================================
# PART 1: GENERATE MULTIFAMILY GROUND LAYER
# ===================================================================
print("\n--- Step 3: Ground Layer Generation ---")

# --- CONFIGURATION FOR CHUNKING ---
CHUNK_SIZE_POINTS = 50000  # Number of points to process in each chunk

# Step 3.1: Process multifamily points in chunks to manage memory
print(f"Processing {len(multifamily_points)} points in chunks of {CHUNK_SIZE_POINTS}...")
total_chunks = (len(multifamily_points) + CHUNK_SIZE_POINTS - 1) // CHUNK_SIZE_POINTS
all_non_developed_polygons = []

for i in range(total_chunks):
    start_idx = i * CHUNK_SIZE_POINTS
    end_idx = min((i + 1) * CHUNK_SIZE_POINTS, len(multifamily_points))
    chunk_points = multifamily_points.iloc[start_idx:end_idx]
    
    print(f"\n--- Processing Chunk {i+1}/{total_chunks} ---")

    # Step 3.1a: Buffer Multifamily Points to Create Area of Interest (AOI) for the chunk
    print("Step 3.1a: Buffering multifamily points to create AOI...")
    aoi = chunk_points.geometry.buffer(250).union_all()
    print("AOI for chunk created successfully.")

    # Step 3.1b: Clip, Clean, and Vectorize Impervious Land from NLCD Raster
    print("\nStep 3.1b: Vectorizing impervious land from raster for chunk...")
    geoms = []
    try:
        with rasterio.open(IMPERVIOUS_RASTER_PATH) as src:
            aoi_gdf = gpd.GeoDataFrame(geometry=[aoi], crs=CRS).to_crs(src.crs)
            if aoi_gdf.empty or aoi_gdf.geometry.iloc[0].is_empty:
                print("Skipping empty or invalid AOI for this chunk.")
                continue
            
            geom_mapping = [mapping(aoi_gdf.geometry.iloc[0])]
            out_image, out_transform = mask(src, geom_mapping, crop=True)

            binary_image = (out_image[0] >= 1) & (out_image[0] <= 20)
            structure = np.ones((3, 3), dtype=bool)
            cleaned_image = binary_opening(binary_image, structure=structure)
            cleaned_image = cleaned_image.astype(np.uint8)

            for geom, val in shapes(cleaned_image, mask=cleaned_image, transform=out_transform):
                if val == 1:
                    geoms.append(shape(geom))

        if geoms:
            chunk_polygons = gpd.GeoDataFrame(geometry=geoms, crs=src.crs).to_crs(CRS)
            all_non_developed_polygons.append(chunk_polygons)
            print(f"Found {len(chunk_polygons)} non-developed land patches in chunk.")
        else:
            print("No non-developed land patches found in this chunk.")

    except Exception as e:
        print(f"  ⚠️ Error processing chunk {i+1}: {e}")
    
    # Clean up memory
    del chunk_points, aoi, geoms
    gc.collect()

# Concatenate results from all chunks
if not all_non_developed_polygons:
    raise ValueError("No non-developed polygons were generated. Check input data and AOI creation.")

non_developed_polygons = pd.concat(all_non_developed_polygons, ignore_index=True)
print(f"\n--- Total of {len(non_developed_polygons)} non-developed land patches found across all chunks. ---")


# Step 3.3: Exclude Schools and Public Lands
print("\nStep 3.3: Excluding institutional and public lands...")
all_school_layers = []
for layer in SCHOOLS_LAYER_NAMES:
    try:
        school_layer = gpd.read_file(SCHOOLS_GDB_PATH, layer=layer).to_crs(CRS)
        all_school_layers.append(school_layer)
    except Exception as e:
        print(f"- WARNING: Could not load school layer '{layer}'. Error: {e}")

schools = pd.concat(all_school_layers, ignore_index=True)
cpad = gpd.read_file(PUBLIC_LANDS_CPAD_PATH).to_crs(CRS)
cced = gpd.read_file(PUBLIC_LANDS_CCED_PATH).to_crs(CRS)
all_exclusions_gdf = pd.concat([schools, cpad, cced], ignore_index=True)

ground_layer_intermediate = non_developed_polygons.overlay(
    all_exclusions_gdf, how='difference', keep_geom_type=True
)
ground_layer_intermediate = ground_layer_intermediate[ground_layer_intermediate.is_valid]
print("Exclusion of institutional lands complete.")

# Step 3.4: Exclude Building Buffers
print("\nStep 3.4: Excluding building footprints...")
all_buildings = gpd.read_file(ALL_BUILDINGS_PATH).to_crs(CRS)
all_buildings['geometry'] = all_buildings.geometry.buffer(5)

final_ground_layer = ground_layer_intermediate.overlay(
    all_buildings, how='difference', keep_geom_type=True
)
final_ground_layer = final_ground_layer[final_ground_layer.is_valid]
final_ground_layer.to_file(FINAL_GROUND_LAYER_PATH, driver="GPKG")
print(f"--- Initial Ground Layer saved to: {FINAL_GROUND_LAYER_PATH} ---")

# ===================================================================
# PART 2: GENERATE MULTIFAMILY ROOFTOP LAYER
# ===================================================================
print("\n--- Step 4: Rooftop Layer Generation ---")
print("Step 4.1: Identifying intersecting building footprints...")

nearby_buildings = gpd.sjoin_nearest(
    all_buildings,
    multifamily_points[['ID', 'geometry']],
    max_distance=15,
    how="inner"
)
unique_multifamily_rooftops = nearby_buildings.drop_duplicates(subset=['geometry'])
unique_multifamily_rooftops.to_file(FINAL_ROOFTOP_LAYER_PATH, driver="GPKG")
print(f"Identified {len(unique_multifamily_rooftops)} unique multifamily building footprints.")
print(f"--- Initial Rooftop Layer saved to: {FINAL_ROOFTOP_LAYER_PATH} ---")

# ===================================================================
# PART 3: REFINE AND CLEAN LAYERS
# ===================================================================
print("\n--- Step 5: Refining and Cleaning Layers ---")

# Step 5.1: Cluster and Simplify Rooftop Layer
print("Step 5.1: Clustering and simplifying rooftop layer...")
rooftops_gdf = gpd.read_file(FINAL_ROOFTOP_LAYER_PATH)
rooftops_gdf['centroid'] = rooftops_gdf.geometry.centroid
centroids = np.array(list(rooftops_gdf['centroid'].apply(lambda p: (p.x, p.y))))

dbscan = DBSCAN(eps=100, min_samples=5)
clusters = dbscan.fit_predict(centroids)
rooftops_gdf['cluster'] = clusters

clustered_rooftops = rooftops_gdf[rooftops_gdf['cluster'] != -1].copy()

# Ensure 'ID' is in the dataframe before selection
if 'ID' not in clustered_rooftops.columns:
    # If ID was lost, we can try to recover it or just proceed without it
    # For now, let's create a placeholder if it's missing to avoid an error
    print("Warning: 'ID' column was lost during processing. Creating a placeholder.")
    clustered_rooftops['ID'] = range(len(clustered_rooftops))


final_rooftops_minimal = clustered_rooftops[['ID', 'geometry', 'cluster']].copy()
final_rooftops_minimal['geometry'] = final_rooftops_minimal.geometry.simplify(tolerance=2, preserve_topology=True)
final_rooftops_minimal['cluster'] = final_rooftops_minimal['cluster'].astype('int32')
final_rooftops_minimal.to_file(CLUSTERED_ROOFTOP_PATH_SIMPLIFIED, driver="GPKG")
print(f"Simplified and clustered rooftop layer saved to: {CLUSTERED_ROOFTOP_PATH_SIMPLIFIED}")

# Step 5.2: Filter and Clean Ground Layer
print("\nStep 5.2: Filtering and cleaning ground layer...")
ground_gdf = gpd.read_file(FINAL_GROUND_LAYER_PATH)
ground_gdf['area_sqm'] = ground_gdf.geometry.area
min_area_sqm = 100
filtered_ground_layer = ground_gdf[ground_gdf['area_sqm'] >= min_area_sqm].copy()

exploded_gdf = filtered_ground_layer.explode(index_parts=False)
exploded_gdf['area_sqm'] = exploded_gdf.geometry.area
final_gdf = exploded_gdf[exploded_gdf['area_sqm'] >= min_area_sqm].copy()
final_gdf = final_gdf.drop(columns=['area_sqm'])
final_gdf.to_file(FINAL_CLEANED_GROUND_LAYER, driver="GPKG")
print(f"Final cleaned ground layer saved to: {FINAL_CLEANED_GROUND_LAYER}")

# ===================================================================
# PART 4: MERGE GROUND LAYER WITH ADDITIONAL IMPERVIOUS SURFACE DATA
# ===================================================================
print("\n--- Step 6: Merging Ground Layer with Secondary Impervious Surface Data ---")
CHUNK_SIZE = 50000
ground_gdf = gpd.read_file(FINAL_CLEANED_GROUND_LAYER).to_crs(CRS)
filtered_polygons = []
total_chunks = (len(ground_gdf) + CHUNK_SIZE - 1) // CHUNK_SIZE

for chunk_idx in range(total_chunks):
    start_idx = chunk_idx * CHUNK_SIZE
    end_idx = min((chunk_idx + 1) * CHUNK_SIZE, len(ground_gdf))
    print(f"Processing chunk {chunk_idx + 1}/{total_chunks}...")
    chunk_gdf = ground_gdf.iloc[start_idx:end_idx].copy()
    valid_polygons = []
    try:
        with rasterio.open(IMPERVIOUS_SURFACE_PATH_2) as src:
            chunk_raster_crs = chunk_gdf.to_crs(src.crs)
            for idx, row in chunk_raster_crs.iterrows():
                if not row.geometry.is_valid: continue
                out_image, _ = mask(src, [row.geometry], crop=True, nodata=src.nodata)
                valid_pixels = out_image[0][out_image[0] != src.nodata]
                if len(valid_pixels) > 0:
                    original_row = ground_gdf.iloc[idx].copy()
                    original_row['has_impervious_overlap'] = True
                    valid_polygons.append(original_row)
        filtered_polygons.extend(valid_polygons)
    except Exception as e:
        print(f"  ⚠️ Error processing chunk {chunk_idx + 1}: {e}")
    del chunk_gdf
    gc.collect()

if filtered_polygons:
    merged_gdf_impervious = gpd.GeoDataFrame(filtered_polygons, crs=CRS)
    merged_gdf_impervious.to_file(OUTPUT_MERGED_IMPERVIOUS_PATH, driver='GPKG')
    print(f"Merged ground layer with impervious surfaces saved to {OUTPUT_MERGED_IMPERVIOUS_PATH}")
else:
    print("No polygons found with impervious surface overlap.")

# ===================================================================
# PART 5: OPTIMIZE FINAL GEOMETRIES
# ===================================================================
print("\n--- Step 7: Optimizing Final Layer Geometries for Size ---")

# Step 7.1: Optimize Ground Layer
print("Optimizing final ground layer...")
try:
    ground_to_optimize = gpd.read_file(OUTPUT_MERGED_IMPERVIOUS_PATH)
    ground_minimal_gdf = ground_to_optimize[['geometry']].copy()
    ground_minimal_gdf['geometry'] = ground_minimal_gdf.geometry.simplify(tolerance=2.0, preserve_topology=True)
    ground_minimal_gdf.to_file(FINAL_GROUND_LAYER_OPTIMIZED, driver="GPKG")
    print(f"Optimized ground layer saved to: {FINAL_GROUND_LAYER_OPTIMIZED}")
except FileNotFoundError:
    print(f"Could not find {OUTPUT_MERGED_IMPERVIOUS_PATH} to optimize. Skipping.")
except Exception as e:
    print(f"An error occurred during ground layer optimization: {e}")

# Note: Rooftop layer was already simplified in Step 5.1.

# ===================================================================
# PART 6: GENERATE SUMMARY CSV BY CENSUS TRACT
# ===================================================================
print("\n--- Step 8: Generating Summary CSV by Census Tract ---")
try:
    ground_gdf = gpd.read_file(FINAL_GROUND_LAYER_OPTIMIZED).to_crs(CRS)
    rooftop_gdf = gpd.read_file(CLUSTERED_ROOFTOP_PATH_SIMPLIFIED).to_crs(CRS)
    tracts_gdf = gpd.read_file(CENSUS_TRACTS_PATH).to_crs(CRS)

    ground_gdf['area_sqm'] = ground_gdf.geometry.area
    rooftop_gdf['area_sqm'] = rooftop_gdf.geometry.area

    ground_tract_join = gpd.sjoin(ground_gdf, tracts_gdf, how="inner", predicate="intersects")
    rooftop_tract_join = gpd.sjoin(rooftop_gdf, tracts_gdf, how="inner", predicate="intersects")

    ground_summary = ground_tract_join.groupby('GEOID').agg(
        ground_total_area_sqm=('area_sqm', 'sum'),
        ground_polygon_count=('area_sqm', 'count')
    ).reset_index()
    rooftop_summary = rooftop_tract_join.groupby('GEOID').agg(
        rooftop_total_area_sqm=('area_sqm', 'sum'),
        rooftop_polygon_count=('area_sqm', 'count')
    ).reset_index()

    final_summary_df = pd.merge(
        ground_summary, rooftop_summary, on='GEOID', how='outer'
    ).fillna(0).rename(columns={'GEOID': 'census_tract_id'})

    final_summary_df.to_csv(OUTPUT_SUMMARY_CSV, index=False)
    print(f"Summary CSV saved to: {OUTPUT_SUMMARY_CSV}")
except FileNotFoundError as e:
    print(f"Error generating summary CSV. A file was not found: {e}")
except Exception as e:
    print(f"An unexpected error occurred during summary generation: {e}")

# ===================================================================
# PART 7: GENERATE INTERACTIVE MAPS
# ===================================================================
print("\n--- Step 9: Generating Interactive Plotly Maps ---")
try:
    summary_df = pd.read_csv(OUTPUT_SUMMARY_CSV, dtype={'census_tract_id': str})
    tracts_gdf = gpd.read_file(CENSUS_TRACTS_PATH, dtype={'GEOID': str}).to_crs("EPSG:4326")

    merged_gdf = tracts_gdf.merge(
        summary_df, left_on='GEOID', right_on='census_tract_id', how='left'
    ).fillna(0)

    maps_to_create = {
        'rooftop_polygon_count': 'Number of Rooftop Polygons',
        'ground_polygon_count': 'Number of Ground Polygons',
        'rooftop_total_area_sqm': 'Total Rooftop Area (sqm)',
        'ground_total_area_sqm': 'Total Ground Area (sqm)'
    }

    filtered_gdf = merged_gdf[(merged_gdf['rooftop_polygon_count'] > 0) | (merged_gdf['ground_polygon_count'] > 0)].copy()

    if not filtered_gdf.empty:
        data_bounds = filtered_gdf.total_bounds
        center_lat = (data_bounds[1] + data_bounds[3]) / 2
        center_lon = (data_bounds[0] + data_bounds[2]) / 2
        zoom_level = 8

        for column, title in maps_to_create.items():
            fig = px.choropleth_mapbox(
                filtered_gdf,
                geojson=filtered_gdf.geometry,
                locations=filtered_gdf.index,
                color=column,
                center={"lat": center_lat, "lon": center_lon},
                mapbox_style="carto-positron",
                zoom=zoom_level,
                opacity=0.8,
                hover_name="GEOID",
                labels={column: title},
                title=f'<b>{title} by Census Tract</b>',
                color_continuous_scale="Viridis"
            )
            fig.update_layout(margin={"r":0,"t":40,"l":0,"b":0})
            output_path = os.path.join(OUTPUT_DIR_MAPS, f'{column}_interactive_map.html')
            fig.write_html(output_path)
            print(f"  - Saved interactive map to {output_path}")
    else:
        print("No data found to map after merging with census tracts.")

except FileNotFoundError as e:
    print(f"Error generating maps. A file was not found: {e}")
except Exception as e:
    print(f"An unexpected error occurred during map generation: {e}")

# ===================================================================
# SCRIPT COMPLETE
# ===================================================================
print("\n🎉 Processing complete!")
