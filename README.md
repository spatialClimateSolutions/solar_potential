# Multifamily Housing Geospatial Analysis

## Overview

This project contains a Python script (`multifamily_analysis.py`) designed to perform a comprehensive geospatial analysis for identifying multifamily housing development potential. The script processes various geospatial datasets to generate two primary outputs:

1.  **Ground Layer**: Identifies potentially developable land parcels located near existing multifamily housing units.
2.  **Rooftop Layer**: Extracts the building footprints of existing multifamily housing structures.

The script automates the entire workflow, including data loading, cleaning, spatial analysis, layer refinement, and the generation of summary statistics and interactive visualizations.

## Features

-   **Data Ingestion**: Loads data from various sources including CSV files, raster images (GeoTIFF), and vector files (Shapefile, GeoPackage).
-   **Area of Interest (AOI) Generation**: Creates a processing AOI by buffering existing multifamily housing locations.
-   **Impervious Surface Analysis**: Vectorizes impervious surfaces from NLCD raster data to identify non-developed areas.
-   **Exclusion Analysis**: Excludes non-developable areas such as existing buildings, schools, and public lands.
-   **Rooftop Identification**: Identifies multifamily building footprints using spatial joins with a comprehensive building footprint dataset.
-   **Data Cleaning & Refinement**:
    -   Clusters rooftop polygons to identify housing complexes.
    -   Simplifies geometries to optimize file size and performance.
    -   Filters small and irrelevant polygons.
-   **Data Aggregation**: Merges final layers with a secondary impervious surface dataset for further refinement.
-   **Summarization**: Aggregates the final ground and rooftop data by census tract and generates a summary CSV file.
-   **Visualization**: Creates interactive choropleth maps in HTML format to visualize the results.

## Requirements

To run this script, you will need Python 3 and the following libraries. You can install them using pip:

```bash
pip install geopandas pandas rasterio shapely numpy scipy scikit-learn plotly
```

## How to Use

### 1. Configuration

Before running the script, you must configure the file paths in the `CONFIGURATION` section at the top of `multifamily_analysis.py`. Update the paths to point to the correct locations of your input files and desired output directories.

**Key Input Files:**
- `HOUSING_CSV_PATH`: CSV file with housing data, including latitude/longitude and a 'single_or_multi' column.
- `IMPERVIOUS_RASTER_PATH`: National Land Cover Database (NLCD) impervious surface raster.
- `ALL_BUILDINGS_PATH`: Comprehensive building footprint data for your region (e.g., Microsoft Building Footprints).
- `SCHOOLS_GDB_PATH`: Geodatabase containing school property polygons.
- `PUBLIC_LANDS_CPAD_PATH` & `PUBLIC_LANDS_CCED_PATH`: Shapefiles for public and conserved lands.
- `CENSUS_TRACTS_PATH`: Shapefile containing census tract boundaries for your area of analysis.
- `IMPERVIOUS_SURFACE_PATH_2`: A secondary impervious surface raster for final merging.

### 2. Running the Script

Once the paths are configured, you can execute the script from your terminal:

```bash
python multifamily_analysis.py
```

The script will print progress updates to the console as it moves through the various processing steps.

## Outputs

The script will generate the following files in the specified output directories:

-   **Optimized Ground Layer**:
    -   `output_ground_layer/final_multifamily_ground_layer_optimized.gpkg`: A GeoPackage file containing the final, cleaned, and size-optimized polygons of potentially developable land.

-   **Clustered Rooftop Layer**:
    -   `output_rooftop_layer/final_multifamily_rooftop_layer_clustered_simplified.gpkg`: A GeoPackage file with the clustered and simplified building footprints of multifamily properties.

-   **Summary Data**:
    -   `final_summary_by_tract.csv`: A CSV file containing aggregated statistics (e.g., total area, polygon count) for both the ground and rooftop layers, summarized by census tract.

-   **Interactive Maps**:
    -   `interactive_plotly_maps/`: This directory will contain several `.html` files. Each file is an interactive map visualizing a different metric (e.g., 'Number of Rooftop Polygons', 'Total Ground Area') by census tract.
