
-----

# Multifamily Housing Geospatial Analysis

## Overview

This project contains a two-part Python workflow designed to perform a comprehensive geospatial analysis for identifying multifamily housing development potential in California.

1.  **Part 1 (`MS_to_CA.py`)**: Extracts a complete California building footprint dataset from Microsoft's global tile-based data.
2.  **Part 2 (`multifamily_analysis.py`)**: Uses the extracted footprints and other geospatial data to generate two primary outputs:
      * **Ground Layer**: Identifies potentially developable land parcels near existing multifamily housing.
      * **Rooftop Layer**: Identifies and clusters the building footprints of existing multifamily housing.

The scripts automate the entire workflow, from data acquisition and cleaning to spatial analysis, layer refinement, and the generation of summary statistics and interactive visualizations.

-----

## Data Download & Sources

All required input data and the final outputs from the scripts can be downloaded from the following link:

https://drive.google.com/drive/folders/1_ncTtJZp6hUCiBmsNqRjZYxZn7POeGBt?usp=sharing

-----

### Data Sources

  * **Impervious Raster**: National Land Cover Database (NLCD) Impervious Surface.
  * **Multifamily Housing Data (CA)**: Zillow housing data
  * **Microsoft Building Footprints (Global)**: The script uses the most recent global dataset, as the dedicated U.S. dataset is outdated.
  * **California Schools/Colleges (cscd\_2021)**: California School Campus Database.
  * **California Public Lands**: California Protected Areas Database (CPAD) and California Conservation Easement Database (CCED).
  * **National Resources Inventory (NRI) Data**: Census tract shapefiles and tables.
  * **Protected Areas Database of the United States (PADUS)**: USGS protected areas data for California.

-----

## How to Use: A Two-Step Process

This analysis requires running two scripts in sequence. First, you will generate the California building footprint data. Second, you will use that data to perform the multifamily housing analysis.

### **Step 1: Extract California Building Footprints**

This initial step prepares the necessary building footprint data for the main analysis.

**Script**: `MS_to_CA.py`

This script automates the extraction and compilation of building footprint data for California by sourcing it from Microsoft's large-scale global dataset. It identifies the data tiles that overlap with California, downloads them, extracts the building footprints, and merges them into a single file.

#### **Inputs & Outputs**

  * **Required Inputs**:
      * `links.csv`: A CSV file containing the download URLs for the global building footprint data tiles. You can get the latest version (called `dataset_links.csv`) from the [Microsoft Building Footprints GitHub](https://github.com/microsoft/GlobalMLBuildingFootprints?tab=readme-ov-file).
      * `tl_2024_us_state.shp`: A shapefile containing the boundaries of all U.S. states.
  * **Primary Output**:
      * `california_building_footprints_complete.gpkg`: A GeoPackage file containing all unique building footprints located within California.

#### **Running the Script**

1.  **Configure Paths**: Open `MS_to_CA.py` and verify the file paths in the `Configuration` section.
2.  **Execute**: Run the script from your terminal:
    ```bash
    python MS_to_CA.py
    ```

-----

### **Step 2: Generate Multifamily Ground and Rooftop Layers**

After generating the California building footprints, you can proceed with the main analysis.

**Script**: `multifamily_analysis.py`

This script takes the building footprints and other datasets to produce the final ground and rooftop layers.

#### **Configuration**

Before running, you must configure the file paths in the `CONFIGURATION` section at the top of `multifamily_analysis.py`.

  * **Crucially, set `ALL_BUILDINGS_PATH` to point to the output from Step 1**:
    `'california_building_footprints_complete.gpkg'`

#### **Running the Script**

Once the paths are configured, execute the script from your terminal:

```bash
python multifamily_analysis.py
```

#### **Analysis Workflow**

  * **Multifamily Ground Layer Generation**:

      * Identifies potentially developable land near multifamily housing through a multi-stage process of buffering, impervious surface analysis (1-20% imperviousness), and exclusion of institutional lands and existing buildings.
      * The layer is further refined against a secondary high-resolution impervious raster and cleaned by removing small polygons (\< 100 m²).

  * **Multifamily Rooftop Layer Generation**:

      * Identifies building footprints associated with multifamily housing via a spatial join (`sjoin_nearest`) within a 15-meter radius.
      * Applies the DBSCAN clustering algorithm (`eps=100m`, `min_samples=5`) to group buildings into dense communities.

  * **Data Aggregation and Visualization**:

      * Summarizes the final ground and rooftop polygon counts and areas by census tract, saving the results to a CSV.
      * Generates interactive Plotly choropleth maps to visualize the results, saved as HTML files.

-----

## Final Outputs

After completing both steps, the following files will be generated:

  * **California Building Footprints**:

      * `california_building_footprints_complete.gpkg`: The comprehensive building footprint dataset for California (from Step 1).

  * **Optimized Ground Layer**:

      * `output_ground_layer/final_multifamily_ground_layer_optimized.gpkg`: A GeoPackage file containing the final, cleaned polygons of potentially developable land.

  * **Clustered Rooftop Layer**:

      * `output_rooftop_layer/final_multifamily_rooftop_layer_clustered_simplified.gpkg`: A GeoPackage file with the clustered and simplified building footprints of multifamily properties.

  * **Summary Data**:

      * `final_summary_by_tract.csv`: A CSV file containing aggregated statistics (total area, polygon count) for both layers, summarized by census tract.

  * **Interactive Maps**:

      * `interactive_plotly_maps/`: A directory containing several `.html` files, each an interactive map visualizing a different metric by census tract.

-----

## Requirements

To run both scripts, you will need Python 3 and the following libraries. You can install them all using pip:

```bash
pip install geopandas pandas rasterio shapely numpy scipy scikit-learn plotly requests pyquadkey2 fiona
```
