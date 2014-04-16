egm4
====

R script for processing output files from PP-Systems EGM-4 IRGA

You can point this script to a directory of output files downloaded from the EGM-4, an Infrared Gas Analyzer made by PP-Systems, Inc. It will process them: identifying outliers, computing gas fluxes, summarizing files, and writing outputs and logs.

The script is written for R 3.0.3 and doesn't use anything particularly exotic. Required packages are listed in the script.

*NOTE* This script doesn't currently read in plot-specific mass (or area) values, so the raw flux values it produces won't be accurate.

Ben Bond-Lamberty
bondlamberty@pnnl.gov


