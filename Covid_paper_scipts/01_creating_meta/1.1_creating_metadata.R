# load libraries ====
# the lab metadata was taken from CD excel file L:\KLON\MMB-LAB\RESEARCH\3-Ongoing-research-projects\Theme_3-Microbiome_and_Resistome\LucKi Gut Study\5 Vriezers\Feces Aliquots\Diepvrieslocaties_LucKi Gut_13-06_2022
# the seq meta plate_Latvia and box_Latvia variables were constructed from files L:\KLON\MMB-LAB\RESEARCH\3-Ongoing-research-projects\Theme_3-Microbiome_and_Resistome\LucKi Gut Study\7 DNA-isolatie\DNA Isolation MMHP\DNA isolation_library preparation _MMHP_John Penders
# data for the last 6 families (1222,1223,1225,1226,1227,1228) were collected by hand from pdfs from uni L drive
# cleaning fecal dates - SPSS file. It was created by JP from SPSS file sampling_questionnaire_date_combinatiefile_final_191106_tm1152_cleaning_1153_onwards_JP14032023

library(tidyverse)
library(haven)

# read in data ====
#loading raw meta files
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")
raw_labmeta_full = read.delim("Input/1_raw_labmeta.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
raw_seqmeta_full = read.delim("Input/2_raw_seqmeta.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
raw_fecal_dates = as.data.frame(read_sav("Input/3_sampling_questionnaire_date_combinatiefile_final_191106_tm1152_cleaning_1153_onwards_JP14032023_ED200323.sav"))
table_to_check_merging = as.data.frame(read_sav("Input/4_230104_gebdat_cleaned_148records_tm1239_SEL_ED.sav"))
basic = as.data.frame(read_sav("Input/5_Mergedfile_0-14m_230113_part1_tm1221_SEL1_ED.sav"))
allergy = as.data.frame(read_sav("Input/6_Mergedfile_0-14m_230123_part1_tm1221_SEL2_ED.sav"))
pet_birth = as.data.frame(read_sav("Input/7_step5_clean_geboortelijst_210308_birth_pets.sav"))
pet_6m = as.data.frame(read_sav("Input/8_step5_clean_6m_210322_pets.sav"))

# cleaning data ====
# replace "-" symbols with "_" in sample names
raw_seqmeta_full$MMHP_SampleID <- str_replace_all(raw_seqmeta_full$MMHP_SampleID,'-','_')
# removing spaces form the column
raw_labmeta_full$Sample_aliquote <- str_replace_all(raw_labmeta_full$Sample_aliquote,' ','')

# removing duplicate sample MMHP_UM_LCK_780 -> N = 930
raw_labmeta <- subset(raw_labmeta_full, MMHP_SampleID !='MMHP_UM_LCK_780')
#make sure sample names (MMHP_UM_LCK_) and column names are unique

# counting number of infant samples per family and total number per family and adding them to lab meta data file
n_infant_samples <- raw_labmeta %>% subset(Individual == "infant") %>% count(Family_ID)
names(n_infant_samples)[names(n_infant_samples) == "n"] <- "n_of_infant_samples"
n_total_samples <- raw_labmeta %>% count(Family_ID)
names(n_total_samples)[names(n_total_samples) == "n"] <- "n_of_samples"
samples = full_join(n_infant_samples, n_total_samples, by = "Family_ID")
lab_premeta = full_join(raw_labmeta, samples, by = "Family_ID")
rm(n_infant_samples,n_total_samples,raw_labmeta_full,samples,raw_labmeta)

#creating sets
fathers <- lab_premeta %>% subset(Individual == "father") %>% count(Family_ID)
names(fathers)[names(fathers) == "n"] <- "n_of_fathers"
mothers <- lab_premeta %>% subset(Individual == "mother") %>% count(Family_ID)
names(mothers)[names(mothers) == "n"] <- "n_of_mothers"
infants <- lab_premeta %>% subset(Individual == "infant") %>% count(Family_ID)
names(infants)[names(infants) == "n"] <- "n_of_infant"
sets = full_join(mothers, fathers, by = "Family_ID")
sets = full_join(sets, infants, by = "Family_ID")
sets <- sets %>% 
  mutate(set = if_else(
    condition = is.na(n_of_fathers)&is.na(n_of_mothers),
    true = "solo",
    false = if_else(is.na(n_of_fathers), "duo", "trio"),
  ))
sets <- select(sets,Family_ID,set)
lab_premeta = full_join(lab_premeta, sets, by = "Family_ID")
rm(fathers,mothers,infants,sets)

# merge seq and lab meta ====
premeta = inner_join(lab_premeta, raw_seqmeta_full, by = "MMHP_SampleID")
rm(lab_premeta,raw_seqmeta_full)

# creating fecal and dob dates file ====
# rename variables
names(raw_fecal_dates)[names(raw_fecal_dates) == "bq_birthdate_ch_new"] <- "DOB_infant"
names(raw_fecal_dates)[names(raw_fecal_dates) == "kindcode"] <- "Family_ID"

# selecting only those columns to be used
raw_fecal_dates = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch1w_date_f_sample_cleaned_JP14032023,
                                                     ch4w_date_f_sample_cleaned_JP14032023,ch8w_date_f_sample_cleaned_JP14032023,
                                                     ch4m_date_f_sample_cleaned_JP14032023,ch5m_date_f_sample_cleaned_JP14032023,
                                                     ch6m_date_f_sample_cleaned_JP14032023,ch9m_date_f_sample_cleaned_JP14032023,
                                                     ch11m_date_f_sample_cleaned_JP14032023,ch14m_date_f_sample_cleaned_JP14032023))

# selecting specific time points in lab metadata
week_1 = subset(premeta, Age_individual == "1-2 weeks")
week_4 = subset(premeta, Age_individual == "4 weeks")
week_8 = subset(premeta, Age_individual == "8 weeks")
month_4 = subset(premeta, Age_individual == "4 months")
month_5 = subset(premeta, Age_individual == "5 months")
month_6 = subset(premeta, Age_individual == "6 months")
month_9 = subset(premeta, Age_individual == "9 months")
month_11 = subset(premeta, Age_individual == "11 months")
month_14 = subset(premeta, Age_individual == "14 months")
mother = subset(premeta, Age_individual == "mother")
father = subset(premeta, Age_individual == "father")

# selecting columns in raw fecal data by each time point 
d_week_1 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch1w_date_f_sample_cleaned_JP14032023))
names(d_week_1)[names(d_week_1) == "ch1w_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_week_4 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch4w_date_f_sample_cleaned_JP14032023))
names(d_week_4)[names(d_week_4) == "ch4w_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_week_8 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch8w_date_f_sample_cleaned_JP14032023))
names(d_week_8)[names(d_week_8) == "ch8w_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_month_4 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch4m_date_f_sample_cleaned_JP14032023))
names(d_month_4)[names(d_month_4) == "ch4m_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_month_5 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch5m_date_f_sample_cleaned_JP14032023))
names(d_month_5)[names(d_month_5) == "ch5m_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_month_6 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch6m_date_f_sample_cleaned_JP14032023))
names(d_month_6)[names(d_month_6) == "ch6m_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_month_9 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch9m_date_f_sample_cleaned_JP14032023))
names(d_month_9)[names(d_month_9) == "ch9m_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_month_11 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch11m_date_f_sample_cleaned_JP14032023))
names(d_month_11)[names(d_month_11) == "ch11m_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_month_14 = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings,ch14m_date_f_sample_cleaned_JP14032023))
names(d_month_14)[names(d_month_14) == "ch14m_date_f_sample_cleaned_JP14032023"] <- "Sample_date"
d_mother = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings))
d_mother <- cbind(d_mother, Sample_date=NA)
d_father = subset(raw_fecal_dates, select = c(Family_ID,DOB_infant,twins,siblings))
d_father <- cbind(d_father, Sample_date=NA)

#merging fecal dates and the rest of meta by time point
meta_1w = left_join(week_1,d_week_1, by="Family_ID") 
meta_4w = left_join(week_4,d_week_4, by="Family_ID") 
meta_8w = left_join(week_8,d_week_8, by="Family_ID") 
meta_4m = left_join(month_4,d_month_4, by="Family_ID") 
meta_5m = left_join(month_5,d_month_5, by="Family_ID") 
meta_6m = left_join(month_6,d_month_6, by="Family_ID") 
meta_9m = left_join(month_9,d_month_9, by="Family_ID") 
meta_11m = left_join(month_11,d_month_11, by="Family_ID") 
meta_14m = left_join(month_14,d_month_14, by="Family_ID")
meta_mother = left_join(mother,d_mother, by="Family_ID")
meta_father = left_join(father,d_father, by="Family_ID")
rm(d_father,d_month_11,d_month_14,d_month_4,d_month_5,d_month_6,d_month_9,d_mother,d_week_1,d_week_4,d_week_8)
rm(father,month_11,month_14,month_4,month_5,month_6,month_9,mother,week_1,week_4,week_8)

# merging into one table and cleaning ====
premeta_fecal = rbind(meta_1w,meta_4w,meta_8w,meta_4m,meta_5m,meta_6m,meta_9m,meta_11m,meta_14m,meta_mother,meta_father)
rm(meta_11m,meta_14m,meta_1w,meta_4m,meta_4w,meta_5m,meta_6m,meta_8w,meta_9m,meta_mother,meta_father)
rm(premeta,raw_fecal_dates)

#test <- premeta_fecal[is.na(premeta_fecal$Sample_date),]
#test2 <- subset(test, Age_individual == "14 months")
#plyr::count(test2$Family_ID)

# creating meta from questioner files ====
# renaming columns 
names(allergy)[names(allergy) == "kindcode"] <- "Family_ID"

names(table_to_check_merging)[names(table_to_check_merging) == "kindcode"] <- "Family_ID"
names(table_to_check_merging)[names(table_to_check_merging) == "bq_birthdate_ch_new"] <- "DOB_check"

names(basic)[names(basic) == "ch1wq1_bristol_type_f"] <- "bristol_1w"
names(basic)[names(basic) == "ch4wq1_bristol_type_f"] <- "bristol_4w"
names(basic)[names(basic) == "ch8wq1_bristol_type_f"] <- "bristol_8w"
names(basic)[names(basic) == "ch4mq1_bristol_type"] <- "bristol_4m"
names(basic)[names(basic) == "ch5mq1_bristol_type"] <- "bristol_5m"
names(basic)[names(basic) == "ch6mq1_bristol_type"] <- "bristol_6m"
names(basic)[names(basic) == "ch9mq1_bristol_type"] <- "bristol_9m"
names(basic)[names(basic) == "ch11mq1_bristol_type"] <- "bristol_11m"
names(basic)[names(basic) == "CH14Mq1_Bristol_type"] <- "bristol_14m"

names(basic)[names(basic) == "AB_child_1w"] <- "inf_ab_1w"
names(basic)[names(basic) == "ch4wq4_ab_freq"] <- "inf_ab_4w"
names(basic)[names(basic) == "ch8wq4_ab_freq"] <- "inf_ab_8w"
names(basic)[names(basic) == "ch4mq5_ab_freq"] <- "inf_ab_4m"
names(basic)[names(basic) == "ch5mq5_ab_freq"] <- "inf_ab_5m"
names(basic)[names(basic) == "ch6mq16_ab_freq"] <- "inf_ab_6m"
names(basic)[names(basic) == "ch9mq5_ab_freq"] <- "inf_ab_9m"
names(basic)[names(basic) == "ch11mq5_ab_freq"] <- "inf_ab_11m"
names(basic)[names(basic) == "CH14Mq19_AB_freq"] <- "inf_ab_14m"

names(basic)[names(basic) == "kindcode"] <- "Family_ID"
names(basic)[names(basic) == "CH14Mq36_smoke"] <- "smoke_14m"
names(basic)[names(basic) == "bqq3_gender"] <- "sex"
names(basic)[names(basic) == "bqq45_siblings"] <- "older_siblings"
names(basic)[names(basic) == "bqq13_delivery_place"] <- "delivery_place"
names(basic)[names(basic) == "bqq15_hospital"] <- "hospital"
names(basic)[names(basic) == "bqq14_delivery_type"] <- "delivery_type"
names(basic)[names(basic) == "bqq1_birthweight"] <- "birthweight"
names(basic)[names(basic) == "ch14mq37_daycare"] <- "daycare_14m"
names(basic)[names(basic) == "bqq2_pregn_weeks"] <- "pregn_weeks"
names(basic)[names(basic) == "bqq23_smoke"] <- "smoke_birth"
names(basic)[names(basic) == "ch6mq31_smoke"] <- "smoke_6m"

# removing columns we won't need
basic_cleaned = subset(basic, select = -c(CH1w_version,ch8w_version,CH5m_version,CH9m_version,CH14M_version,
                                          CH4w_version,CH4m_version,CH6m_version,CH11m_version,BQ_version))
# removing rows we won't need, as we need families only up to 1228 and only the one that were metagenomicly sequenced
check <- filter(table_to_check_merging,!Family_ID %in% c(1239,1238,1237,1236,1232,1230,1229,1051,1056,1089))
allergy_cleaned <- filter(allergy,!Family_ID %in% c(1051,1056,1089))
basic_cleaned <- filter(basic_cleaned,!Family_ID %in% c(1051,1056,1089))
rm(basic,table_to_check_merging,allergy)


# joining tables
basic_cleaned_allergy = full_join(basic_cleaned, allergy_cleaned, by = "Family_ID")
basic_cleaned_allergy_check = full_join(basic_cleaned_allergy, check, by = "Family_ID")
# removing 6 missing families as no data for them
quest_with_missings <- filter(basic_cleaned_allergy_check,!Family_ID %in% c(1222,1223,1225,1226,1227,1228))
rm(basic_cleaned_allergy,basic_cleaned_allergy_check,check,allergy_cleaned,basic_cleaned)
#premeta2 will should contain 138 families, as family 1224 lacks completely any meta data.

#turning labelled meta into factors with levels, as SPSS provides labels
labelled_quest_with_missings <- as_factor(quest_with_missings)

# adding pet data ====
names(pet_birth)[names(pet_birth) == "kindcode"] <- "Family_ID"
names(pet_birth)[names(pet_birth) == "bqq25_duringpregn_text"] <- "animal_pregnancy"
names(pet_birth)[names(pet_birth) == "bqq25_thismoment_text"] <- "animal_after_pregnancy"
names(pet_birth)[names(pet_birth) == "Fur_pets"] <- "Fur_pets_check"
names(pet_6m)[names(pet_6m) == "kindcode"] <- "Family_ID"
names(pet_6m)[names(pet_6m) == "ch6mq33_pets_home_text"] <- "pets_type_home_6m"
names(pet_6m)[names(pet_6m) == "ch6mq33_pets_daycare_text"] <- "pets_type_daycare_6m"
pet_6m <- select(pet_6m,Family_ID,pets_type_home_6m,pets_type_daycare_6m)
pets <- full_join(pet_birth,pet_6m, by = "Family_ID")
pets <- as_factor(pets) #make labels into levels
pets <- filter(pets,!Family_ID %in% c(1051,1056,1089)) #these families have no metagenomic data
pets <- filter(pets,!Family_ID %in% c(1222,1223)) #these families are in missing database
labelled_quest_with_missings_pets = full_join(labelled_quest_with_missings, pets, by = "Family_ID")
rm(pet_6m,pet_birth,pets)

# writing down created files ====
# saving meta file as .txt file to change format as tab delimited for R
#write.table (labelled_quest_with_missings_pets, "S1.1_labelled_quest_with_missings.txt", row.names=FALSE,sep = "\t")
#write.table (premeta_fecal, "S1.1_premeta_fecal.txt", row.names=FALSE,sep = "\t")















