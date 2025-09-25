import geopandas as gpd
import pandas as pd
import os

# --- 1. Setup File Paths and Load Data ---
# A good practice in Python is to define paths as variables.
GROUND_PATH = "../data/final_multifamily_ground_layer_optimized.gpkg"
ROOF_PATH = "../data/final_multifamily_rooftop_layer_clustered_simplified.gpkg"
OUTPUT_DIR = "../derived/"
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "ground_joined.shp")

print("Loading data...")
ground_gdf = gpd.read_file(GROUND_PATH)
roof_gdf = gpd.read_file(ROOF_PATH)

# --- 2. Calculate the Center of Each Cluster and Create Buffers ---
# This section is equivalent to your dplyr pipe chain for creating roof_buffered.

print("Calculating cluster centers and creating buffers...")
# Get the centroid of each individual rooftop
roof_centroids = roof_gdf.copy()
roof_centroids['geometry'] = roof_centroids.geometry.centroid

# Group by cluster and calculate the mean center point for each cluster
# This is equivalent to `group_by(cluster) %>% summarise(x=mean(x), y=mean(y))`
cluster_centers = roof_centroids.dissolve(by='cluster', aggfunc='mean')
cluster_centers = cluster_centers.centroid

# Convert the cluster center points back to a GeoDataFrame
# We need to preserve the cluster ID, which is now in the index
cluster_centers_gdf = gpd.GeoDataFrame(geometry=cluster_centers, crs=roof_gdf.crs)
cluster_centers_gdf = cluster_centers_gdf.reset_index() # Move 'cluster' from index to column

# Buffer these cluster centers by 150 meters
cluster_buffers_gdf = cluster_centers_gdf.copy()
cluster_buffers_gdf['geometry'] = cluster_buffers_gdf.geometry.buffer(150)

# --- 3. Filter Ground Polygons Intersecting the Buffers ---
# This is equivalent to `st_intersects` and the `lengths(idx) > 0` filter.
# A spatial join is the most efficient way to do this in geopandas.

print("Filtering ground polygons within the buffer zones...")
# Perform an 'inner' spatial join to keep only ground polygons that intersect a buffer
sjoin_result_gdf  = gpd.sjoin(ground_gdf, cluster_buffers_gdf, how="inner", predicate="intersects")
is_first_occurrence = ~sjoin_result_gdf.index.duplicated(keep='first')

ground_filtered_gdf = sjoin_result_gdf[is_first_occurrence].drop(columns=['index_right'])

print(f"Saving final joined layer to {OUTPUT_PATH}...")
# Create the output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
ground_filtered_gdf.to_file(OUTPUT_PATH)

print("Operation complete.")