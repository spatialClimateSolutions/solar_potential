"""
process_california_final.py
---------------------------

Purpose:
    This script extracts and compiles building footprint data for California
    from large-scale, tile-based datasets of the United States. It identifies 
    relevant data tiles via quadkeys, downloads the associated files, 
    extracts building footprints within the California boundary, 
    and combines them into a single GeoPackage output.

Inputs:
    1. links.csv
       - A CSV file containing URLs to building footprint data.
       - Must include columns for location and file URL.

    2. tl_2024_us_state.shp
       - A U.S. states shapefile, required to identify the California boundary.

Outputs:
    - california_building_footprints_complete.gpkg
      A GeoPackage file containing all unique California building footprints.

Dependencies:
    - Python 3.8+
    - pandas
    - geopandas
    - shapely
    - requests
    - pyquadkey2
    - gzip
    - fiona (for shapefile support)

Usage:
    1. Update the configuration section (file paths, shapefile location, 
       output file name) before running.
    2. Run as a standalone script:
           python process_california_final.py
    3. The script will:
       - Load the California boundary
       - Parse the input CSV for U.S. building footprint URLs
       - Filter URLs to tiles overlapping California
       - Download the relevant files
       - Extract, clean, and combine building footprints
       - Save the final dataset to GeoPackage

"""

import os
import requests
import pandas as pd
import geopandas as gpd
import warnings
import gzip
from urllib.parse import urlparse
from pyquadkey2 import quadkey
from shapely.geometry import box

# Suppress non-critical warnings
warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

# --- Configuration (verify before running) ---

# 1. The name of the CSV file containing all the URLs.
URL_CSV_FILE = 'links.csv'

# 2. The folder containing the U.S. states shapefile.
ALL_STATES_BOUNDARY_FILE = 'tl_2024_us_state'
SHAPEFILE_FOLDER = 'tl_2024_us_state'
SHAPEFILE_NAME = 'tl_2024_us_state.shp'

# 3. The name for the final output file.
OUTPUT_FILE = 'california_building_footprints_complete.gpkg'

# --- End of Configuration ---

# --- Script Functions ---
DOWNLOAD_DIR = 'building_data_ca_only'


def get_california_boundary():
    """Loads the California boundary from the shapefile and returns its geometry."""
    print("Loading California boundary...")
    shapefile_full_path = os.path.join(SHAPEFILE_FOLDER, SHAPEFILE_NAME)
    try:
        ca_boundary = gpd.read_file(shapefile_full_path).to_crs("EPSG:4326")
        print("Boundary loaded successfully.")
        return ca_boundary
    except Exception as e:
        print(f"Error reading shapefile: {e}")
        return None


def download_california_files(urls_to_download):
    """Downloads files from URLs into the local directory."""
    if not os.path.exists(DOWNLOAD_DIR):
        os.makedirs(DOWNLOAD_DIR)

    print(f"\nDownloading {len(urls_to_download)} relevant data files...")

    for i, url in enumerate(urls_to_download):
        filename = os.path.basename(urlparse(url).path)
        file_save_path = os.path.join(DOWNLOAD_DIR, filename)

        if not os.path.exists(file_save_path):
            print(f"({i+1}/{len(urls_to_download)}) Downloading: {filename}")
            try:
                response = requests.get(url, stream=True)
                response.raise_for_status()
                with open(file_save_path, 'wb') as f_out:
                    f_out.write(response.content)
            except Exception as e:
                print(f" -> Could not download {filename}. Error: {e}")
        else:
            print(f"({i+1}/{len(urls_to_download)}) Already exists: {filename}")


def process_local_files(ca_boundary):
    """Processes downloaded files and extracts California building footprints."""
    all_ca_buildings = []

    try:
        files_to_process = [f for f in os.listdir(DOWNLOAD_DIR) if f.endswith('.csv.gz')]
        if not files_to_process:
            print(f"Error: No .csv.gz files found in '{DOWNLOAD_DIR}'.")
            return
    except FileNotFoundError:
        print(f"Error: The directory '{DOWNLOAD_DIR}' does not exist.")
        return

    print(f"\nProcessing {len(files_to_process)} downloaded files...")

    for i, filename in enumerate(files_to_process):
        filepath = os.path.join(DOWNLOAD_DIR, filename)

        try:
            with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                gdf = gpd.read_file(f, driver='GeoJSONSeq')

            if gdf.empty:
                continue

            # Assign CRS to building footprint data
            gdf.crs = "EPSG:4326"

            # Select only buildings intersecting the California boundary
            buildings_in_ca = gpd.sjoin(gdf, ca_boundary, how="inner", predicate="intersects")

            if not buildings_in_ca.empty:
                print(f"({i+1}/{len(files_to_process)}) Found {len(buildings_in_ca):,} buildings in {filename}")
                all_ca_buildings.append(buildings_in_ca)

        except Exception as e:
            print(f"Could not process {filename}. Error: {e}")

    if all_ca_buildings:
        print("\nCombining all building footprints into a single dataset...")
        final_gdf = pd.concat(all_ca_buildings, ignore_index=True)

        # Keep only geometry
        final_gdf = final_gdf[['geometry']]

        # Remove duplicates
        final_gdf = final_gdf.drop_duplicates(subset='geometry').reset_index(drop=True)

        print(f"Saving {len(final_gdf):,} total building footprints to {OUTPUT_FILE}...")
        final_gdf.to_file(OUTPUT_FILE, driver='GPKG')
        print(f"Completed. Saved {len(final_gdf):,} California building footprints.")
    else:
        print("No buildings were found. Script finished without producing data.")


def filter_urls_for_boundary(url_list, boundary_gdf):
    """Filters a list of URLs to only those whose tiles intersect the California boundary."""
    print("\nDetermining relevant data tiles for California...")
    relevant_urls = []
    total_urls = len(url_list)
    boundary_unary = boundary_gdf.union_all()

    parsed_count = 0
    error_count = 0

    for i, url in enumerate(url_list):
        if i % 100 == 0:  # Progress update every 100 URLs
            print(f"\rChecking URL {i+1}/{total_urls}", end="", flush=True)
        try:
            if 'quadkey=' not in url:
                raise ValueError("No quadkey parameter found in URL")

            qk_str = url.split('quadkey=')[1].split('/')[0]
            qk = quadkey.QuadKey(qk_str)

            tile_result = qk.to_tile()
            if len(tile_result) == 2:
                if isinstance(tile_result[0], tuple):
                    (tile_x, tile_y), zoom = tile_result
                else:
                    tile_x, tile_y = tile_result
                    zoom = len(qk_str)
            else:
                tile_x, tile_y, zoom = tile_result

            import math

            def tile_to_lat_lon(x, y, z):
                n = 2.0 ** z
                lon_deg = x / n * 360.0 - 180.0
                lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * y / n)))
                lat_deg = math.degrees(lat_rad)
                return lat_deg, lon_deg

            # Get bounds of the tile
            min_lat, min_lon = tile_to_lat_lon(tile_x, tile_y + 1, zoom)
            max_lat, max_lon = tile_to_lat_lon(tile_x + 1, tile_y, zoom)

            tile_bbox = box(min_lon, min_lat, max_lon, max_lat)
            parsed_count += 1

            if tile_bbox.intersects(boundary_unary):
                relevant_urls.append(url)
        except Exception:
            error_count += 1
            continue

    print(f"\nSuccessfully parsed {parsed_count} quadkeys ({error_count} errors).")
    print(f"Identified {len(relevant_urls)} relevant data tiles for California.")
    return relevant_urls


# --- Main Execution ---
if __name__ == "__main__":
    california_boundary = get_california_boundary()

    if california_boundary is not None:
        try:
            df = pd.read_csv(URL_CSV_FILE, on_bad_lines='warn')
            df.columns = df.columns.str.strip()

            # Identify location and URL columns
            location_col = None
            url_col = None

            for col in df.columns:
                if 'location' in col.lower():
                    location_col = col
                if 'url' in col.lower() or 'link' in col.lower():
                    url_col = col

            if location_col is None or url_col is None:
                print(f"Error: Could not find location and URL columns. Available columns: {list(df.columns)}")
                exit()

            print(f"Using '{location_col}' as location column and '{url_col}' as URL column")

            # Extract U.S. URLs
            us_variations = ['UnitedStates', 'United States', 'US', 'USA', 'united_states', 'united states']
            us_urls = []

            for variation in us_variations:
                matches = df[df[location_col].str.strip().str.contains(variation, case=False, na=False)]
                if not matches.empty:
                    us_urls = matches[url_col].tolist()
                    print(f"Found {len(us_urls)} URLs for United States")
                    break

            if not us_urls:
                print("Error: No URLs for United States found.")
                exit()
        except FileNotFoundError:
            print(f"Error: The file '{URL_CSV_FILE}' was not found.")
            exit()
        except Exception as e:
            print(f"Error reading CSV: {e}")
            exit()

        if us_urls:
            relevant_urls = filter_urls_for_boundary(us_urls, california_boundary)
            if relevant_urls:
                download_california_files(relevant_urls)
                process_local_files(california_boundary)
            else:
                print("No data tiles intersect with California.")
