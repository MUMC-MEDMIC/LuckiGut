# load libraries ====
library(tidyverse)
library(haven)
library(dplyr)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
# read in .sav files from SPSS, read as data frame to change into normal format
envir_raw = as.data.frame(read_sav("Input/10_Mergedfile_tm1223_210506_v4_2_SEL2CR_Garden_PlayGround_Sand.sav"))
#read in environmental data for the 6 missing families
missing_env = read.delim("Input/11_missing_seven_families_environmental_data.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
month9 = as.data.frame(read_sav("Input/12_step5_rawdata_clean_9m_210506.sav"))
month11 = as.data.frame(read_sav("Input/13_step4_rawdata_BASICFILEvarnameadjusted_11M_220403.sav"))
month14 = as.data.frame(read_sav("Input/14_step4_BASICFILE1001_1223_220324_14M.sav"))

# cleaning data ====
# renaming columns 
names(envir_raw)[names(envir_raw) == "kindcode"] <- "Family_ID"
#switching/changing to labels
envir_label <- as_factor(envir_raw)
#change all non-numeric values to characters
envir_label <- envir_label %>%
  mutate_if(~ !is.numeric(.), as.character)
#change all "NA" values to missing values NA in character columns
envir_label <- envir_label %>%
  mutate_if(is.character, ~ifelse(. == "NA", NA, .))
#change all "" values to missing values NA in numeric columns
envir_label <- envir_label %>%
  mutate_if(is.character, ~ifelse(. == "", NA, .))
#change all "missing" values to missing values NA in numeric columns
envir_label <- envir_label %>%
  mutate_if(is.character, ~ifelse(. == "missing", NA, .))

rm(envir_raw)

# renaming columns ====
names(envir_label)[names(envir_label) == "kindcode"] <- "Family_ID"
names(envir_label)[names(envir_label) == "ch4mq21_sand_playground"] <- "SandPlayground_4m"
names(envir_label)[names(envir_label) == "ch5mq21_sand_playground"] <- "SandPlayground_5m"
names(envir_label)[names(envir_label) == "ch4mq21_grass_garden"] <- "GrassGarden_4m"
names(envir_label)[names(envir_label) == "ch5mq21_grass_garden"] <- "GrassGarden_5m"
names(envir_label)[names(envir_label) == "ch4mq21_grass_playground"] <- "GrassPlayground_4m"
names(envir_label)[names(envir_label) == "ch5mq21_grass_playground"] <- "GrassPlayground_5m"
names(envir_label)[names(envir_label) == "ch4mq21_sand_garden"] <- "SandGarden_4m"
names(envir_label)[names(envir_label) == "ch5mq21_sand_garden"] <- "SandGarden_5m"
names(envir_label)[names(envir_label) == "bqq26_stables"] <- "stables"
names(envir_label)[names(envir_label) == "ch4mq21_farm"] <- "farm_4m"
names(envir_label)[names(envir_label) == "ch5mq21_farm"] <- "farm_5m"
names(envir_label)[names(envir_label) == "bqq26_stables_contact"] <- "stables_contact"

# removing families with no metagenomic data ====
envir_filtered <- filter(envir_label,!Family_ID %in% c(1051,1056,1089)) 
#removing families to avoid doubles. These two had only one question answered in MM file, the answer was transferred to missing file. These families were checked manually and added to missing file.
envir_filtered_2 <- filter(envir_filtered,!Family_ID %in% c(1222,1223)) 

# adding info for missing TP 9, 11 and 14 months ====
names(month9)[names(month9) == "kindcode"] <- "Family_ID"
names(month11)[names(month11) == "kindcode"] <- "Family_ID"
names(month14)[names(month14) == "kindcode"] <- "Family_ID"

#switching/changing to labels
month9 <- as_factor(month9)
month11 <- as_factor(month11)
month14 <- as_factor(month14)

#picking environmental columns
cols_with_envir_9 <- names(month9)[grepl("22", names(month9))] 
cols_with_envir_11 <- names(month11)[grepl("22", names(month11))] 
cols_with_envir_14 <- names(month14)[grepl("39", names(month14))]

#selecting them into one dataframe
missing_TP_enviroment_9 <- month9[, c(cols_with_envir_9, "Family_ID")]
missing_TP_enviroment_11 <- month11[, c(cols_with_envir_11, "Family_ID")]
missing_TP_enviroment_14 <- month14[, c(cols_with_envir_14, "Family_ID")]
rm(cols_with_envir_9,cols_with_envir_11,cols_with_envir_14, month9,month11,month14)

#renaming the columns
names(missing_TP_enviroment_9)[names(missing_TP_enviroment_9) == "ch9mq22_farm"] <- "farm_9m"
names(missing_TP_enviroment_9)[names(missing_TP_enviroment_9) == "ch9mq22_sand_garden"] <- "SandGarden_9m"
names(missing_TP_enviroment_9)[names(missing_TP_enviroment_9) == "ch9mq22_sand_playground"] <- "SandPlayground_9m"
names(missing_TP_enviroment_9)[names(missing_TP_enviroment_9) == "ch9mq22_grass_garden"] <- "GrassGarden_9m"
names(missing_TP_enviroment_9)[names(missing_TP_enviroment_9) == "ch9mq22_grass_playground"] <- "GrassPlayground_9m"

names(missing_TP_enviroment_11)[names(missing_TP_enviroment_11) == "ch11mq22_farm"] <- "farm_11m"
names(missing_TP_enviroment_11)[names(missing_TP_enviroment_11) == "ch11mq22_sand_garden"] <- "SandGarden_11m"
names(missing_TP_enviroment_11)[names(missing_TP_enviroment_11) == "ch11mq22_sand_playground"] <- "SandPlayground_11m"
names(missing_TP_enviroment_11)[names(missing_TP_enviroment_11) == "ch11mq22_grass_garden"] <- "GrassGarden_11m"
names(missing_TP_enviroment_11)[names(missing_TP_enviroment_11) == "ch11mq22_grass_playground"] <- "GrassPlayground_11m"

names(missing_TP_enviroment_14)[names(missing_TP_enviroment_14) == "CH14Mq39_farm"] <- "farm_14m"
names(missing_TP_enviroment_14)[names(missing_TP_enviroment_14) == "CH14Mq39_sand_garden"] <- "SandGarden_14m"
names(missing_TP_enviroment_14)[names(missing_TP_enviroment_14) == "CH14Mq39_sand_playground"] <- "SandPlayground_14m"
names(missing_TP_enviroment_14)[names(missing_TP_enviroment_14) == "CH14Mq39_grass_garden"] <- "GrassGarden_14m"
names(missing_TP_enviroment_14)[names(missing_TP_enviroment_14) == "CH14Mq39_grass_playground"] <- "GrassPlayground_14m"

# merging together TP ====
missing_TP_enviroment_9_11 <- full_join(missing_TP_enviroment_9,missing_TP_enviroment_11,by = "Family_ID")
missing_TP_enviroment_9_11_14 <- full_join(missing_TP_enviroment_9_11,missing_TP_enviroment_14,by = "Family_ID")
rm(missing_TP_enviroment_9,missing_TP_enviroment_11,missing_TP_enviroment_14,missing_TP_enviroment_9_11)

# merging main data provided by MM with missing families of 6, that MM did not include ====
#removing duplicate families from main file, as they are added manually to missing_env file
missing_TP_enviroment_9_11_14 <- filter(missing_TP_enviroment_9_11_14,!Family_ID %in% c(1222,1225))
#merging files
envir_all_TP <- full_join(envir_filtered_2, missing_TP_enviroment_9_11_14,by = "Family_ID")
envir_all_TP_all_family = rbind(envir_all_TP,missing_env)
rm(envir_label,envir_filtered,envir_filtered_2, missing_env,envir_all_TP,missing_TP_enviroment_9_11_14)

# cleaning data ====
#making all Yes -> yes, No -> no
envir_all_TP_all_family[envir_all_TP_all_family == "Yes"] <- "yes"
envir_all_TP_all_family[envir_all_TP_all_family == "No"] <- "no"

envir_all_TP_all_family[envir_all_TP_all_family == "2"] <- "yes"
envir_all_TP_all_family[envir_all_TP_all_family == "1"] <- "no"
envir_all_TP_all_family[envir_all_TP_all_family == "0"] <- NA

#uniting all yes contacts with stables to one column
envir_all_TP_all_family <- envir_all_TP_all_family %>%
  mutate(stables_yn = if_else(
    condition =  .$stables == "no", 
    true = "no", 
    false = if_else(condition =  is.na(stables), 
                    true = NA, 
                    false = "yes")))

#renaming numbers into animals, according to spss file
envir_all_TP_all_family$stables_contact <- gsub('1', 'stables', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('2', 'horse(s)', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('3', 'Cow(s)', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('4', 'goats(s)', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('5', 'sheep', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('6', 'pig(s)', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('7', 'poultry', envir_all_TP_all_family$stables_contact)
envir_all_TP_all_family$stables_contact <- gsub('8', 'pigeons', envir_all_TP_all_family$stables_contact)

#if stables are no, then contact with animal is no, not NA
envir_all_TP_all_family <- envir_all_TP_all_family %>%
  mutate(stables_contact = if_else(
    condition =  .$stables_yn == "no", 
    true = "no", 
    false = if_else(condition =  is.na(stables_yn), 
                    true = NA, 
                    false = .$stables_contact)))

#creating a new column for each of the animals
envir_all_TP_all_family$stables_contact_horse <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("horse", envir_all_TP_all_family$stables_contact), "yes", "no"))
envir_all_TP_all_family$stables_contact_cow <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("Cow", envir_all_TP_all_family$stables_contact), "yes", "no"))
envir_all_TP_all_family$stables_contact_goat <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("goat", envir_all_TP_all_family$stables_contact), "yes", "no"))
envir_all_TP_all_family$stables_contact_sheep <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("sheep", envir_all_TP_all_family$stables_contact), "yes", "no"))
envir_all_TP_all_family$stables_contact_pig <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("pig", envir_all_TP_all_family$stables_contact), "yes", "no"))
envir_all_TP_all_family$stables_contact_poultry <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("poultry", envir_all_TP_all_family$stables_contact), "yes", "no"))
envir_all_TP_all_family$stables_contact_pigeons <- ifelse(is.na(envir_all_TP_all_family$stables_contact), NA, ifelse(grepl("pigeons", envir_all_TP_all_family$stables_contact), "yes", "no"))

# writing data ====
#write.table (envir_all_TP_all_family, "S1.3_environmental_cleaned_data.txt", row.names=FALSE,sep = "\t")

# check ====
test <- envir_all_TP_all_family %>% filter(stables_contact_sheep == "no")
colnames(envir_all_TP_all_family)


