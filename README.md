# MS-263 Final Project
## Project Summary

This data project will investigate the interaction between regional and upwelling winds, and their relative impact on the rate of turbulence dissipation for a small upwelling bay (San Luis Obispo Bay) in the California Current System.

The repository includes:
- A Jupyter Notebook that performs structured analysis of oceanographic microstructure data using Python functions.
- A complete archive of MATLAB functions and ODAS-based analysis used during the project.
- All code and supporting files are located in the `python_analysis/` and `MATLAB_analysis_MS-263/` folders.
- Due to size limitations, some raw data files are hosted externally on Google Drive.

This analysis was developed for MS-263 (Spring 2025) and is designed to be run on systems with either Jupyter Notebooks or MATLAB installed.

### Google drive data access: insert link

## Steps needed to run the analysis code on another computer
To run Jupyter Notebook containing the final project, copy and paste the python_analysis folder into Juptyter Notebooks/Jupyter Labs. All data needed for running the notebook are included in the folder and any necessary data_paths are correct as long as the structure of the folder remains unchanged.
## Location of data
### Processed/cleaned data:
Located in `MATLAB_analysis_MS-263/data/`
### Raw ODAS .p files and NetCDFs:
Available from Google Drive: need link
### Python-parsed or reformatted datasets:
Located in `python_analysis/` or downloaded from the Drive folder above

## Source
Data were collected during the 2024 San Luis Obispo Bay field campaign as part of the SUFEX project (insert link to website):
Instrumentation:
* MicroCTD
* NOAA MET buoys
* Shipboard MET data

## Dependencies
### Python:
* numpy
* pandas
* matplotlib
* cartopy
* cmocean
* Final_Project_Functions.py (provided in repo)
### MATLAB:
* ODAS Library (from Rockland Scientific) — included in odas_functions/
* TEOS-10 GSW Toolbox — included in Imported_Functions/
* Custom MATLAB scripts in:
     - My_functions/
     - RKW_Functions/

