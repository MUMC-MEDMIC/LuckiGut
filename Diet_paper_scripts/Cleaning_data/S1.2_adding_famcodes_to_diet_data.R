# libraries ====
library(tidyverse)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
baku_f_I <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
load("Input/4_timepoints_general_4m_to_14m.RData") #list, with six elements (timepoints)
load("Input/5_other_variables.RData") #non-dietary variables

# adding family codes ====
timepoints_general_4m_to_14m[[1]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[2]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[3]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[4]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[5]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[6]]$Family_ID <- other_variables$ï..kindcode

# saving data ====
#save(timepoints_general_4m_to_14m, file = "S1.2_timepoints_general_4m_to_14m_famcodes.RData")

