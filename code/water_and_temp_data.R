library(tidyverse)
library(here)


# ---------------------------------- WARNING ----------------------------------------
# Some of the data that you have obtained from this U.S. Geological Survey database
# may not have received Director's approval. Any such data values are qualified
# as provisional and are subject to revision. Provisional data are released on the
# condition that neither the USGS nor the United States Government may be held liable
# for any damages resulting from its use.
#
# Additional info: https://help.waterdata.usgs.gov/policies/provisional-data-statement
#
# File-format description:  https://help.waterdata.usgs.gov/faq/about-tab-delimited-output
# Automated-retrieval info: https://help.waterdata.usgs.gov/faq/automated-retrievals
#
# Contact:   gs-w_waterdata_support@usgs.gov
# retrieved: 2023-04-07 08:22:26 EDT       (sdww02)
#
# Data for the following 1 site(s) are contained in this file
#    USGS 11446500 AMERICAN R A FAIR OAKS CA
# -----------------------------------------------------------------------------------
#
# Data provided for site 11446500
#            TS   parameter     statistic     Description
#         10973       00010     00001     Temperature, water, degrees Celsius (Maximum)
#         10974       00010     00002     Temperature, water, degrees Celsius (Minimum)
#         10975       00010     00008     Temperature, water, degrees Celsius (Median)
#         10976       00010     00011     Temperature, water, degrees Celsius (Instantaneous)
#         10977       00060     00003     Discharge, cubic feet per second (Mean)
#        234322       00010     00003     Temperature, water, degrees Celsius (Mean)
#
# Data-value qualification codes included in this output:
#        
#     A  Approved for publication -- Processing and review completed.
#     P  Provisional data subject to revision.
#     R  Records for these data have been revised. [https://waterdata.usgs.gov/usa/nwis/revision/?site_no=11446500&ts_ids=10974,10973,234322]
#     e  Value has been estimated.
#     4  Statistic computed from less than expected number of instantaneous values for the period
# 


dt <- fread(input = glue::glue("grep -v '5s' {here('data', 'dv.tsv')}"), header=TRUE, drop=c("agency_cd", "site_no"))

labels <- c(
  "datetime",
  "max_temp", 
  "max_temp_cd", 
  "min_temp", 
  "min_temp_cd", 
  "median_temp", 
  "median_temp_cd", 
  "inst_temp", 
  "inst_temp_cd", 
  "discharge_cfs",  
  "discharge_cfs_cd", 
  "mean_temp", 
  "mean_temp_cd"
  )

dt %>% 
  setnames(old = names(.), new = labels)

daily_data <- dt[year(datetime) == 2019, c("datetime", "mean_temp", "discharge_cfs")] %>% na.omit()

daily_data %>% 
  ggplot(aes(datetime, discharge_cfs)) +
  geom_line()
