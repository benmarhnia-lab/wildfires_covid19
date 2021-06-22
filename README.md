#### wildfires_covid19 ####

##Data sources

Raw datasets used in code provided are available at the following link: https://www.dropbox.com/sh/ywa94ldyrlcltg2/AADEA9aNUj2b7bVw_xYSyg2ca?dl=0
This includes the following data:

Daily COVID-19 cases and death data by County
CDC COVID tracker
Link: https://covid.cdc.gov/covid-data-tracker/
Dataset name in dropbox: covid_confirmed_usafacts.csv & covid_deaths_usafacts.csv

Wildfire smoke data
National Oceanic and Atmospheric Administration
Office of Satellite and Product Operations
Hazard Mapping System Fire and Smoke Product
Link: https://www.ospo.noaa.gov/Products/land/hms.html
Dataset name in dropbox: smk_intsct_county_area_2020.csv

Mobility data
U.S. Department of Transportation 
Bureau of Transportation Statistics
Link: https://data.bts.gov/Research-and-Statistics/Trips-by-Distance/w96p-f2qv
Dataset name in dropbox: Trips_by_Distance.csv

## Code 

File 1: 1_Compile the data_61021.do
STATA do file used to compile all data sources for analysis
Updated: June 10th, 2021

File 2: 2_R_Data diagnostics and revision_61021.R
R code used for data diagnostics and further preparation for analysis
Updated: June 10th, 2021

File 3: 3_Data_analysis_61021.do
STATA do file used for analysis and visualization of results 
Updated: June 10th, 2021

