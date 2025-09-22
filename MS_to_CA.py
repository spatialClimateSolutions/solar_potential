# -*- coding: utf-8 -*-
"""
MS_to_CA_footprints.py

Description:
    This script automates the acquisition and processing of building footprint
    data for the state of California from a global dataset. It performs the
    following key operations:
    1.  Reads a list of data URLs from a standard CSV file.
    2.  Filters the list to identify data tiles relevant to the United States.
    3.  Downloads the relevant compressed GeoJSON-L (.csv.gz) files in parallel
        to maximize efficiency.
    4.  Processes each downloaded file in parallel, extracting building
        footprints that intersect with the California state boundary.
    5.  Aggregates the results, removes duplicate geometries, and saves the
        final, comprehensive dataset to a single GeoPackage file.


Inputs:
    - links.csv: A CSV file containing download URLs for the building
      footprint data. It contains 'Location', 'QuadKey', 'Url', 'Size',
      and 'UploadDate' columns in standard CSV format.
    - tl_2024_us_state/tl_2024_us_state.shp: A TIGER/Line shapefile of U.S.
      state boundaries, used to define the area of interest for California.

Outputs:
    - california_building_footprints_complete.gpkg: A GeoPackage file
      containing the aggregated and deduplicated building footprint geometries
      for California.

"""

# --- Standard Library Imports ---
import os
import gzip
import re
import warnings
import concurrent.futures
from urllib.parse import urlparse

# --- Third-Party Library Imports ---
import pandas as pd
import geopandas as gpd
import requests
from tqdm import tqdm

# --- Global Configuration ---

# Suppress a specific, non-critical UserWarning from GeoPandas to declutter output.
warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

# --- File and Directory Paths ---
# NOTE: Verify these paths before execution.
URL_SOURCE_FILE = 'links.csv'
SHAPEFILE_DIR = 'tl_2024_us_state'
SHAPEFILE_NAME = 'tl_2024_us_state.shp'
DOWNLOAD_DIR = 'building_data_us'
OUTPUT_GPKG = 'california_building_footprints_complete.gpkg'

# --- Data Filtering Parameters ---
CALIFORNIA_STATE_NAME = 'California'


# --- Helper Functions for Parallel Processing ---

def download_single_file(url_and_quadkey):
    """
    Downloads a single data file if it does not already exist locally.

    A unique filename is generated using the file's quadkey to prevent
    naming collisions from different data sources with identical filenames.

    Args:
        url_and_quadkey (tuple): A tuple containing the download URL (str)
                                 and the associated QuadKey (str).

    Returns:
        str: An error message if the download fails.
        None: If the download is successful or the file already exists.
    """
    url, quadkey = url_and_quadkey
    try:
        original_filename = os.path.basename(urlparse(url).path)
        unique_filename = f"quadkey_{quadkey}_{original_filename}"
        file_save_path = os.path.join(DOWNLOAD_DIR, unique_filename)

        if not os.path.exists(file_save_path):
            response = requests.get(url, stream=True, timeout=60)
            response.raise_for_status()  # Raise an exception for bad status codes
            with open(file_save_path, 'wb') as f_out:
                f_out.write(response.content)
        return None
    except requests.exceptions.RequestException as e:
        return f"Failed to download {url}: {e}"


def process_single_file(task):
    """
    Reads a gzipped GeoJSON-L file, extracts building footprints within
    the California boundary, and reprojects them.

    Args:
        task (tuple): A tuple containing the file path (str) and the
                      California state boundary GeoDataFrame.

    Returns:
        GeoDataFrame: A GeoDataFrame containing geometries of buildings within
                      California, or None if no buildings are found or an
                      error occurs.
    """
    filepath, ca_boundary_utm = task
    try:
        # The source files are gzipped text files where each line is a GeoJSON object.
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            # Explicitly use the 'fiona' engine, as it correctly handles the
            # GeoJSON-L (newline-delimited GeoJSON) format.
            gdf = gpd.read_file(f, engine='fiona')

        if gdf.empty:
            return None

        # Project the building data to match the California boundary's CRS.
        gdf = gdf.set_crs("EPSG:4326", allow_override=True)
        gdf = gdf.to_crs(ca_boundary_utm.crs)

        # Perform a spatial join to find buildings that intersect with the boundary.
        buildings_in_ca = gpd.sjoin(gdf, ca_boundary_utm, how="inner", predicate="intersects")

        if not buildings_in_ca.empty:
            # Return only the essential geometry column to conserve memory.
            return buildings_in_ca[['geometry']]
        return None
    except Exception as e:
        # Log a warning if a file is malformed or cannot be processed.
        print(f"\n[Warning] Could not process {os.path.basename(filepath)}. Error: {e}\n")
        return None


# --- Main Workflow Functions ---

def get_california_boundary():
    """
    Loads the US states shapefile and extracts the boundary for California.

    The boundary is re-projected to a suitable UTM zone (NAD83 / UTM zone 10N)
    to ensure accurate spatial analysis and metric units.

    Returns:
        GeoDataFrame: A GeoDataFrame containing the California boundary, or None
                      if an error occurs.
    """
    print("Loading California state boundary...")
    shapefile_path = os.path.join(SHAPEFILE_DIR, SHAPEFILE_NAME)
    try:
        states_gdf = gpd.read_file(shapefile_path)
        ca_boundary = states_gdf[states_gdf['NAME'] == CALIFORNIA_STATE_NAME].copy()

        if ca_boundary.empty:
            print(f"Error: State '{CALIFORNIA_STATE_NAME}' not found in shapefile.")
            return None

        # Reproject to a CRS suitable for California (e.g., NAD83 / UTM zone 10N).
        # This improves the accuracy of geometric operations.
        ca_boundary_utm = ca_boundary.to_crs("EPSG:26910")

        print("California boundary loaded and projected successfully.")
        return ca_boundary_utm
    except Exception as e:
        print(f"Error reading or processing shapefile: {e}")
        return None


def download_files_parallel(df):
    """
    Manages the parallel download of data files listed in the DataFrame.

    Args:
        df (DataFrame): A pandas DataFrame containing 'Url' and 'QuadKey' columns.
    """
    if df.empty:
        print("No data files to download.")
        return

    if not os.path.exists(DOWNLOAD_DIR):
        os.makedirs(DOWNLOAD_DIR)

    print(f"\nInitiating download of {len(df)} files using parallel processing...")
    url_quadkey_pairs = list(zip(df['Url'], df['QuadKey']))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(download_single_file, url_quadkey_pairs), total=len(url_quadkey_pairs)))

    errors = [res for res in results if res]
    if errors:
        print(f"\nEncountered {len(errors)} download errors.")
        print("See the first 5 examples below:")
        for error in errors[:5]:
            print(f" - {error}")
    print("\nDownload process completed.")


def process_files_parallel(ca_boundary_utm):
    """
    Manages the parallel processing of downloaded files to extract footprints.

    Args:
        ca_boundary_utm (GeoDataFrame): The projected California state boundary.
    """
    try:
        all_files = os.listdir(DOWNLOAD_DIR)
        files_to_process = [os.path.join(DOWNLOAD_DIR, f) for f in all_files if f.endswith('.csv.gz')]
        if not files_to_process:
            print(f"Error: No .csv.gz files found in the download directory '{DOWNLOAD_DIR}'.")
            return
    except FileNotFoundError:
        print(f"Error: The download directory '{DOWNLOAD_DIR}' does not exist.")
        return

    print(f"\nProcessing {len(files_to_process)} downloaded files...")
    tasks = [(filepath, ca_boundary_utm) for filepath in files_to_process]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(process_single_file, tasks), total=len(tasks)))

    all_ca_gdfs = [gdf for gdf in results if gdf is not None]

    print("\n--- Processing Summary ---")
    print(f"  - Files processed: {len(files_to_process)}")
    print(f"  - Files containing California buildings: {len(all_ca_gdfs)}")

    if not all_ca_gdfs:
        print("\nWarning: No building footprints were found within the California boundary.")
        return

    # --- Data Aggregation and Export ---
    print("\nCombining all building footprints into a single dataset...")
    final_gdf = pd.concat(all_ca_gdfs, ignore_index=True)
    initial_count = len(final_gdf)
    print(f"Found {initial_count:,} total building footprints before deduplication.")

    final_gdf = final_gdf.drop_duplicates(subset='geometry').reset_index(drop=True)
    final_count = len(final_gdf)
    print(f"Removed {initial_count - final_count:,} duplicate entries.")

    print(f"Saving {final_count:,} unique building footprints to {OUTPUT_GPKG}...")
    final_gdf.to_file(OUTPUT_GPKG, driver='GPKG')
    print(f"Operation complete. Final dataset saved to '{OUTPUT_GPKG}'.")

    bounds = final_gdf.total_bounds
    print("\n--- Final Dataset Extent (UTM Zone 10N) ---")
    print(f"  - Easting (X):  {bounds[0]:.2f} to {bounds[2]:.2f}")
    print(f"  - Northing (Y): {bounds[1]:.2f} to {bounds[3]:.2f}")


def main():
    """
    Main function to execute the full data acquisition and processing workflow.
    """
    ca_boundary = get_california_boundary()
    if ca_boundary is None:
        print("Halting execution because the California boundary could not be loaded.")
        return

    try:
        df = pd.read_csv(URL_SOURCE_FILE)
        
        # Ensure column names are lowercase and stripped of whitespace
        df.columns = [col.strip().lower() for col in df.columns]

        required_cols = {'location', 'url', 'quadkey'}
        if not required_cols.issubset(df.columns):
            print(f"Error: CSV must contain {required_cols}.")
            print(f"Available columns: {list(df.columns)}")
            return

        # Filter for United States entries
        us_mask = df['location'].str.strip().str.lower() == 'unitedstates'
        target_df = df[us_mask].copy()

        if target_df.empty:
            print("No data for 'UnitedStates' found in the URL file.")
            return

        target_df.dropna(subset=['url'], inplace=True)
        target_df = target_df[target_df['url'].str.startswith('http')]
        target_df.rename(columns={'url': 'Url', 'quadkey': 'QuadKey'}, inplace=True)
        print(f"Identified {len(target_df)} files for 'UnitedStates' to be processed.")

    except Exception as e:
        print(f"An error occurred while reading or filtering the URL file: {e}")
        return

    if not target_df.empty:
        download_files_parallel(target_df)
        process_files_parallel(ca_boundary)
    else:
        print("No valid URLs were found to process.")


# --- Script Execution ---
if __name__ == "__main__":
    main()
