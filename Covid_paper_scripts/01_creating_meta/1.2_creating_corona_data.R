# load in libraries ====
library(tidyverse)
library(haven)
library(dplyr)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
corona_raw = as.data.frame(read_sav("Input/9_corona_n44_240109.sav"))

# renaming columns ==== 
names(corona_raw)[names(corona_raw) == "kindcode"] <- "Family_ID"
#switching/changing to labels
corona_label <- as_factor(corona_raw)
#change all non-numeric values to characters
corona_label <- corona_label %>%
  mutate_if(~ !is.numeric(.), as.character)
#change all "NA" values to missing values NA in character columns
corona_label <- corona_label %>%
  mutate_if(is.character, ~ifelse(. == "NA", NA, .))
#change all "" values to missing values NA in numeric columns
corona_label <- corona_label %>%
  mutate_if(is.character, ~ifelse(. == "", NA, .))

# saving table ====
#write.table (corona_label, "S1.2_corona_uncleaned_data.txt", row.names=FALSE,sep = "\t")
rm(corona_raw)

# cleaning data MOTHER+GENERAL ====
#cleaning definitely typos
corona_label$corm_othhome_in_chnr[grepl("Q",corona_label$corm_othhome_in_chnr)]<-"0"
corona_label$natr_incorona_freq[grepl("3,5",corona_label$natr_incorona_freq)]<-"3.5"
corona_label$natr_precorona_freq[grepl("1,5",corona_label$natr_precorona_freq)]<-"1.5"
#works in uni -> no children there
corona_label[corona_label$Family_ID =="1222", "corm_pspace_chinr"] <- "0"
#renaming level in one variable: natr_postcor_apprec
corona_label["natr_postcor_apprec"][corona_label["natr_postcor_apprec"] == "yes, more appreciation"] <- "yes"
#columns to replace NAs for no change "same" of nature visits during corona, as same is already used in the data
columns_to_replace_NAs <- c("natr_change_garden", "natr_change_park","natr_change_woods", "natr_change_rural", "natr_change_beach")
for (col in columns_to_replace_NAs) {
  corona_label[[col]][is.na(corona_label[[col]])] <- "same"
}
rm(columns_to_replace_NAs)

#replacing "<maak uw keuze>" in whole df with NAs
corona_label <- corona_label %>%
  mutate_all(~ replace(., . == "<maak uw keuze>", NA))

# removing empty columns with only NAs for the answers
corona_label <- corona_label %>% select(-c(natr_change_othname1,natr_change_othname2,natr_change_othfreq1,natr_change_othfreq2,`otherchildren#`,household_othchild))
#renaming variables in selected columns

corona_label$household_adults_MM <- corona_label$household_adults
corona_label$household_adults[grepl("1,2,3",corona_label$household_adults)]<-"mother, father, other"
corona_label$household_adults[grepl("1,2",corona_label$household_adults)]<-"mother, father"
corona_label$household_adults[grepl("1,3",corona_label$household_adults)]<-"mother, other"
# this "niet van toepassing (want mijn kind is dan niet buiten (of op dezelfde plek) als het bezoek)" was under number 6, will be replaced with "always"
corona_label$cor_homeout_dist_cha[grepl("6",corona_label$cor_homeout_dist_cha)]<-"always"
corona_label$cor_homeout_dist_chc[grepl("6",corona_label$cor_homeout_dist_chc)]<-"always"
# this "NA kind is dan niet in zelfde ruimte" will be changed to "always" in all df
corona_label <- corona_label %>%
  mutate_all(~ replace(., . == "NA kind is dan niet in zelfde ruimte", "always"))
#if days in nature 0 -> min in nature also 0 from NA
corona_label <- corona_label %>%
  mutate(natr_precorona_min = if_else(
    condition =  .$natr_precorona_freq == "0", 
    true = "0", 
    false = natr_precorona_min))
corona_label <- corona_label %>%
  mutate(natr_incorona_min = if_else(
    condition =  .$natr_incorona_freq == "0", 
    true = "0", 
    false = natr_incorona_min))
#fixing time according to the pattern of answering, no change in nature pre/post corona 7 minutes is typo
corona_label[corona_label$Family_ID=="1221", "natr_precorona_min"] <- "60"
#fixing 1 minute with average number, as impossible to be 1 minute for families 1157 and 1198
num_vector_precorona_min <- as.numeric(corona_label$natr_precorona_min)
mean(num_vector_precorona_min,na.rm = TRUE) # 37.74359 -> we will use 40 minutes
corona_label <- corona_label %>%
  mutate(natr_precorona_min = if_else(
    condition =  .$natr_precorona_min == "1", 
    true = "40", 
    false = natr_precorona_min))
rm(num_vector_precorona_min)
#the same families had 1 also during corona, we will take the same number 40, as there was no change in nature
corona_label[corona_label$Family_ID=="1157", "natr_incorona_min"] <- "40"
corona_label[corona_label$Family_ID=="1198", "natr_incorona_min"] <- "40"
#siblings
sibling_columns <- c("household_osibs", "household_ysibs")
for (col in sibling_columns) {
  corona_label[[col]][is.na(corona_label[[col]])] <- "no"
}
sibling_columns2 <- c("oldersibs#","youngersibs#")
for (col in sibling_columns2) {
  corona_label[[col]][is.na(corona_label[[col]])] <- "0"
}
#renaming the columns as names are problematic
names(corona_label)[names(corona_label) == "oldersibs#"] <- "oldersibs_num"
names(corona_label)[names(corona_label) == "youngersibs#"] <- "youngersibs_num"
rm(sibling_columns, sibling_columns2)

#if there was no contact then NA are filled with "no contact" or NA will be replaced with 0 adults/infants MOTHER/INFANT
columns_to_replace <- c("cor_home_in_adultnr", "cor_home_in_childnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$cor_home_in == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("cor_hom_in_dist_self", "cor_hom_in_dist_cha","cor_hom_in_dist_chic")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$cor_home_in == "no contact", "no contact", x)
})
rm(columns_to_replace)
columns_to_replace <- c("cor_home_out_adultnr", "cor_home_out_childnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$cor_home_out == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("cor_homout_dist_self", "cor_homeout_dist_cha","cor_homeout_dist_chc")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$cor_home_out == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corm_work_adultnr", "corm_work_childnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_work == "I'm not working at the moment", "0", ifelse(corona_label$corm_work == "I work from home (100%)", "0",x))
})
rm(columns_to_replace)
columns_to_replace <- c("corm_work_dist_self", "corm_work_cap","corm_work_screen")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_work == "I'm not working at the moment", "no contact", ifelse(corona_label$corm_work == "I work from home (100%)", "no contact",x))
})
rm(columns_to_replace)
#family 1228 mother is not working at the moment
columns_to_replace <- c("corm_work_hospital", "corm_work_nothospit","corm_work_daycare","corm_work_primschool","corm_work_secschool",	"corm_work_bus",
                        "corm_work_otherocc")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$Family_ID == "1228", "no", x)
})
rm(columns_to_replace) 

columns_to_replace <- c("corm_pspace_in_adunr", "corm_pspace_chinr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_publspace_in == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corm_pspace_in_dist", "corm_pspace_in_cap")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_publspace_in == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corm_pspace_out_adnr", "corm_pspace_out_chnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_publicspace_out == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corm_pspace_out_dist", "corm_pspace_out_cap")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_publicspace_out == "no contact", "no contact", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corm_othhome_in_adnr", "corm_othhome_in_chnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_othhome_in == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corm_othhome_in_dist", "corm_othhome_in_cap")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_othhome_in == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corm_othhom_out_adnr", "corm_othhom_out_chnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_othhome_out == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corm_othhom_out_dist", "corm_othhome_out_cap")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_othhome_out == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corm_transp_adnr", "corm_trans_chnr")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_transp == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corm_transp_dist", "corm_transp_cap")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_transp == "no contact", "no contact", x)
})
rm(columns_to_replace) 

#if no covid test was done -> no test MOTHER
columns_to_replace <- c("corm_covid19_testpos")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corm_covid19_test == "no", "no test", x)
})
rm(columns_to_replace) 

# cleaning data FATHER ====
columns_to_replace <- c("corf_work_adultnr1", "corf_work_childnr1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_work1 == "I'm not working at the moment", "0", ifelse(corona_label$corf_work1 == "I work from home (100%)", "0",x))
})
rm(columns_to_replace)
columns_to_replace <- c("corf_work_dist_self1", "corf_work_cap1","corf_work_screen1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_work1 == "I'm not working at the moment", "no contact", ifelse(corona_label$corf_work1 == "I work from home (100%)", "no contact",x))
})
rm(columns_to_replace)
columns_to_replace <- c("corf_pspace_in_adun1", "corf_pspace_chinr1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_publspace_in1 == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corf_pspace_in_dist1", "corf_pspace_in_cap1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_publspace_in1 == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corf_pspace_out_adn1", "corf_pspace_out_chn1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_publicspace_ou1 == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corf_pspace_out_dis1", "corf_pspace_out_cap1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_publicspace_ou1 == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corf_othhome_in_adn1", "corf_othhome_in_chn1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_othhome_in1 == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corf_othhome_in_dis1", "corf_othhome_in_cap1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_othhome_in1 == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corf_othhom_out_adn1", "corf_othhom_out_chn1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_othhome_out1 == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corf_othhom_out_dis1", "corf_othhome_out_ca1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_othhome_out1 == "no contact", "no contact", x)
})
rm(columns_to_replace) 
columns_to_replace <- c("corf_transp_adnr1", "corf_trans_chnr1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_transp1 == "no contact", "0", x)
})
rm(columns_to_replace)
columns_to_replace <- c("corf_transp_dist1", "corf_transp_cap1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_transp1 == "no contact", "no contact", x)
})
rm(columns_to_replace) 
#if no covid test was done -> no test MOTHER
columns_to_replace <- c("corf_covid19_testpo1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(corona_label$corf_covid19_test1 == "no", "no test", x)
})
rm(columns_to_replace) 

# PARENTS ====
# changing NAs of infants into 0, if adults is given, and infants was left blank. Probably was 0 but they did not put it. 
corona_label <- corona_label %>%
  mutate(cor_home_in_childnr = if_else(
    condition =  is.na(cor_home_in_childnr)&!is.na(cor_home_in_adultnr), 
    true = "0", 
    false = cor_home_in_childnr))
corona_label <- corona_label %>%
  mutate(cor_home_out_childnr = if_else(
    condition =  is.na(cor_home_out_childnr)&!is.na(cor_home_out_adultnr), 
    true = "0", 
    false = cor_home_out_childnr))

corona_label <- corona_label %>%
  mutate(corm_pspace_chinr = if_else(
    condition =  is.na(corm_pspace_chinr)&!is.na(corm_pspace_in_adunr), 
    true = "0", 
    false = corm_pspace_chinr))

corona_label <- corona_label %>%
  mutate(corm_pspace_out_chnr = if_else(
    condition =  is.na(corm_pspace_out_chnr)&!is.na(corm_pspace_out_adnr), 
    true = "0", 
    false = corm_pspace_out_chnr))

corona_label <- corona_label %>%
  mutate(corm_othhome_in_chnr = if_else(
    condition =  is.na(corm_othhome_in_chnr)&!is.na(corm_othhome_in_adnr), 
    true = "0", 
    false = corm_othhome_in_chnr))

corona_label <- corona_label %>%
  mutate(corm_othhom_out_chnr = if_else(
    condition =  is.na(corm_othhom_out_chnr)&!is.na(corm_othhom_out_adnr), 
    true = "0", 
    false = corm_othhom_out_chnr))

corona_label <- corona_label %>%
  mutate(corf_work_childnr1 = if_else(
    condition =  is.na(corf_work_childnr1)&!is.na(corf_work_adultnr1), 
    true = "0", 
    false = corf_work_childnr1))

corona_label <- corona_label %>%
  mutate(corf_pspace_chinr1 = if_else(
    condition =  is.na(corf_pspace_chinr1)&!is.na(corf_pspace_in_adun1), 
    true = "0", 
    false = corf_pspace_chinr1))

corona_label <- corona_label %>%
  mutate(corf_pspace_out_chn1 = if_else(
    condition =  is.na(corf_pspace_out_chn1)&!is.na(corf_pspace_out_adn1), 
    true = "0", 
    false = corf_pspace_out_chn1))

corona_label <- corona_label %>%
  mutate(corf_pspace_out_chn1 = if_else(
    condition =  is.na(corf_pspace_out_chn1)&!is.na(corf_pspace_out_adn1), 
    true = "0", 
    false = corf_pspace_out_chn1))

corona_label <- corona_label %>%
  mutate(corf_othhome_in_chn1 = if_else(
    condition =  is.na(corf_othhome_in_chn1)&!is.na(corf_othhome_in_adn1), 
    true = "0", 
    false = corf_othhome_in_chn1))

corona_label <- corona_label %>%
  mutate(corf_othhom_out_chn1 = if_else(
    condition =  is.na(corf_othhom_out_chn1)&!is.na(corf_othhom_out_adn1), 
    true = "0", 
    false = corf_othhom_out_chn1))

corona_label <- corona_label %>%
  mutate(corf_trans_chnr1 = if_else(
    condition =  is.na(corf_trans_chnr1)&!is.na(corf_transp_adnr1), 
    true = "0", 
    false = corf_trans_chnr1))

#replacing 0 values with NA
columns_to_replace <- c("corf_work1","corf_work_dist_self1", "corf_work_cap1","corf_work_screen1", "corf_work_hospital1",	
                        "corf_work_nothospit1",	"corf_work_daycare1",	"corf_work_primschoo1",	"corf_work_secschool1",
                        "corf_work_bus1",	"corf_work_otherocc1","corf_publspace_in1","corf_pspace_in_dist1", "corf_pspace_in_cap1",
                        "corf_publicspace_ou1", "corf_pspace_out_dis1", "corf_pspace_out_cap1","corf_othhome_in1",
                        "corf_othhome_in_dis1", "corf_othhome_in_cap1","corf_othhome_out1","corf_othhom_out_dis1",
                        "corf_othhome_out_ca1","corf_transp1","corf_transp_dist1", "corf_transp_cap1","corf_hands_water1",
                        "corf_hands_watersoa1",	"corf_hands_desinf1",	"corf_fever1",	"corf_shortofbreath1",	"corf_nosmelltaste1",
                        "corf_tired1",	"corf_jointpain1",	"corf_headache1",	"corf_sorethroat1",	"corf_diarrhea1",	"corf_vomit1",
                        "corf_cough1",	"corf_sneeze1",	"corf_nose1",	"corf_eyes1",	"corf_covid19_infect1",	"corf_covid19_test1", "corf_covid19_testpo1")
corona_label[columns_to_replace] <- lapply(corona_label[columns_to_replace], function(x) {
  ifelse(x == "0", NA, x)
})
rm(columns_to_replace)

# creating new variables ====
#covid symptoms FATHER
columns_to_run <- c("corf_fever1", "corf_shortofbreath1", "corf_nosmelltaste1", "corf_tired1","corf_jointpain1","corf_headache1",
                    "corf_sorethroat1","corf_diarrhea1","corf_vomit1","corf_cough1","corf_sneeze1" ,"corf_nose1", "corf_eyes1")

corona_label <- corona_label %>%
  mutate(f_covid_symptoms = rowSums(select(.,all_of(columns_to_run)) == "yes"))
rm(columns_to_run)

#covid symptoms MOTHER
columns_to_run <- c("corm_fever", "corm_shortofbreath", "corm_nosmelltaste", "corm_tired","corm_jointpain","corm_headache",
                    "corm_sorethroat","corm_diarrhea","corm_vomit","corm_cough","corm_sneeze" ,"corm_nose", "corm_eyes")
corona_label <- corona_label %>%
  mutate(m_covid_symptoms = rowSums(select(.,all_of(columns_to_run)) == "yes"))
rm(columns_to_run)

#renaming columns to original names, removing ending "_duplicate"
corona_label <- corona_label %>% rename_with(~gsub("_duplicate", "", .), contains("_duplicate"))

plyr::count(corona_label$corf_nosmelltaste1)
grep("_duplicate", names(corona_label), value = TRUE)

#time in total in nature before and during corona per week
corona_label$during_corona_time_nature_min_per_week <-as.numeric(corona_label$natr_incorona_freq)*as.numeric(corona_label$natr_incorona_min)
corona_label$before_corona_time_nature_min_per_week <-as.numeric(corona_label$natr_precorona_freq)*as.numeric(corona_label$natr_precorona_min)

# Removing families without metagenomic data ====
corona_label <- filter(corona_label,!Family_ID %in% c(1229,1232,1235,1236,1237,1238,1239,1240)) 

# histogram ====
hist_data <- hist(corona_label$before_corona_time_nature_min_per_week, plot = FALSE)
proportions <- hist_data$counts / length(corona_label$before_corona_time_nature_min_per_week)
barplot(proportions, names.arg = hist_data$breaks[-length(hist_data$breaks)], 
        xlab = "Time in nature before corona", ylab = "Proportions", main = "Histogram")

# writing data ====
#write.table (corona_label, "S1.2_corona_cleaned_data.txt", row.names=FALSE,sep = "\t")

