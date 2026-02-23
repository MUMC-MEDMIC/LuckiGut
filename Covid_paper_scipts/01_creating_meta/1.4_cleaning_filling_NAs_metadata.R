# load libraries ====
#meta file was created in 1_creating_metadata.R
#additional data was cleaned in 2a_corona_data.R and 2b_environmental_data.R
#Filling NAs in metadata and regrouping variables

library(tidyverse)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# load in data ====
#merging main quest file and missing families together into one quest file
#missing file was created from pdf from uni L drive, for family 1224 from Ldot as family has no pdfs
missing = read.delim("Input/15_missing_seven_families.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
labelled_quest_with_missings = read.delim("Output/Files/S1.1_labelled_quest_with_missings.txt",sep = "\t",header = TRUE, na.strings=c("","NA","missing","unknown"))
premeta_fecal = read.delim("Output/Files/S1.1_premeta_fecal.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
fam112 = read.delim("L:/KLON/MMB-LAB/RESEARCH/3-Ongoing-research-projects/Theme_3-Microbiome_and_Resistome/LucKi Gut Study/Michal_Skawinski_project/Generated_files/11_clusters_solids.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
samples389 = read.delim("L:/KLON/MMB-LAB/RESEARCH/3-Ongoing-research-projects/Theme_3-Microbiome_and_Resistome/LucKi Gut Study/Michal_Skawinski_project/Generated_files/10_microbial_age_df.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# merging data together ====
quest = rbind(labelled_quest_with_missings,missing)
rm(missing,labelled_quest_with_missings)

# regrouping variables ====

# ****************************************************
#                      delivery
# ****************************************************
#type of delivery, finding all vaginal deliveries and group them like one
quest$delivery_type_MM <- quest$delivery_type
quest$delivery_type[grepl("vaginal",quest$delivery_type)]<-"vaginal"

# ****************************************************
#                      infant's antibiotics
# ****************************************************
#merging never and No together for 1 week antibiotics
quest$inf_ab_1w_MM <- quest$inf_ab_1w
quest$inf_ab_1w[grepl("No",quest$inf_ab_1w)]<-"never"
quest$inf_ab_1w[grepl("Yes",quest$inf_ab_1w)]<-"1 time"

#renaming "2" antibiotics to "more than 1" antibiotics, as the rest TP have "more than 1".
quest$inf_ab_6m_MM <- quest$inf_ab_6m
quest$inf_ab_6m[grepl("2",quest$inf_ab_6m)]<-"more than 1 time"

#1 time antibiotics and more than 1 time are united as subgroups are too small, other TP has no more than 1 time answers
quest$inf_ab_9m_MM <- quest$inf_ab_9m
quest$inf_ab_11m_MM <- quest$inf_ab_11m
quest$inf_ab_14m_MM <- quest$inf_ab_14m
# regrouping as yes for all antibiotic intake for all TP
quest$inf_ab_1w[grepl("time",quest$inf_ab_1w)]<-"yes"
quest$inf_ab_4w[grepl("time",quest$inf_ab_4w)]<-"yes"
quest$inf_ab_8w[grepl("time",quest$inf_ab_8w)]<-"yes"
quest$inf_ab_4m[grepl("time",quest$inf_ab_4m)]<-"yes"
quest$inf_ab_5m[grepl("time",quest$inf_ab_5m)]<-"yes"
quest$inf_ab_6m[grepl("time",quest$inf_ab_6m)]<-"yes"
quest$inf_ab_9m[grepl("time",quest$inf_ab_9m)]<-"yes"
quest$inf_ab_11m[grepl("time",quest$inf_ab_11m)]<-"yes"
quest$inf_ab_14m[grepl("time",quest$inf_ab_14m)]<-"yes"

plyr::count(quest$inf_ab_11m) 
test <- quest %>% select(c(inf_ab_11m,Family_ID))

# ****************************************************
#                      mother's antibiotics
# ****************************************************
# regrouping mother's antibiotics
quest$mat_AB_use_final_MM <- quest$mat_AB_use_final
quest$mat_AB_use_final[grepl("no|after delivery only",quest$mat_AB_use_final)]<-"No"
quest$mat_AB_use_final[grepl("(peri)",quest$mat_AB_use_final)]<-"Yes, peri partum"
quest$mat_AB_use_final[grepl("during pregnancy only|during pregnancy and after",quest$mat_AB_use_final)]<-"Yes, during pregnancy"
quest$mat_AB_use_final[grepl("partus",quest$mat_AB_use_final)]<-"Yes, pregnancy and delivery"

#creating variables for linear regression yes/no for AB at pregnancy and delivery
quest$mat_AB_use_final_preg[grepl("pregnancy",quest$mat_AB_use_final)]<-"yes"
quest$mat_AB_use_final_preg[grepl("partum",quest$mat_AB_use_final)]<-"no"
quest$mat_AB_use_final_preg[grepl("No",quest$mat_AB_use_final)]<-"no"
quest$mat_AB_use_final_delivery[grepl("partum",quest$mat_AB_use_final)]<-"yes"
quest$mat_AB_use_final_delivery[grepl("delivery",quest$mat_AB_use_final)]<-"yes"
quest$mat_AB_use_final_delivery[grepl("during pregnancy",quest$mat_AB_use_final)]<-"no"
quest$mat_AB_use_final_delivery[grepl("No",quest$mat_AB_use_final)]<-"no"

# ****************************************************
#                      family history of allergy
# ****************************************************
#famhistory regrouping
quest$FamHistory_MM <- quest$FamHistory
quest <- quest %>% select(-c(FamHistory))
quest$atopy_mother[grepl("mother",quest$FamHistory_MM)]<-"yes"
quest$atopy_mother[grepl("only father",quest$FamHistory_MM)]<-"no"
quest$atopy_mother[grepl("No",quest$FamHistory_MM)]<-"no"
quest$atopy_father[grepl("father",quest$FamHistory_MM)]<-"yes"
quest$atopy_father[grepl("only mother",quest$FamHistory_MM)]<-"no"
quest$atopy_father[grepl("No",quest$FamHistory_MM)]<-"no"

# ****************************************************
#                      daycare
# ****************************************************
quest$daycare_6m_MM <- quest$daycare_6m
quest$daycare_6m[grepl("daycare",quest$daycare_6m)]<-"yes"
quest$daycare_6m[grepl("guest",quest$daycare_6m)]<-"yes"
quest$daycare_6m[grepl("No",quest$daycare_6m)]<-"no"
quest$daycare_14m_MM <- quest$daycare_14m
quest$daycare_14m[grepl("Yes",quest$daycare_14m)]<-"yes"
quest$daycare_14m[grepl("No",quest$daycare_14m)]<-"no"

# ****************************************************
#                      bristol
# ****************************************************
#cleaning bristol as has small and big letters in the variables
test <- quest %>% select(Family_ID, matches("bristol"))

quest <- quest %>% mutate(across(matches("bristol"), ~ gsub("type", "Type", .)))
plyr::count(quest$bristol_1w)
plyr::count(quest$bristol_4w)
plyr::count(quest$bristol_8w)
plyr::count(quest$bristol_4m)
plyr::count(quest$bristol_5m)
plyr::count(quest$bristol_6m)
plyr::count(quest$bristol_9m)
plyr::count(quest$bristol_11m)
plyr::count(quest$bristol_14m)

# ****************************************************
#                      smoking before pregnancy
# ****************************************************
quest$smoke_before_preg[grepl("before",quest$smoke_birth)]<-"yes"
quest$smoke_before_preg[grepl("Yes",quest$smoke_birth)]<-"yes"
quest$smoke_before_preg[grepl("never",quest$smoke_birth)]<-"no"

# ****************************************************
#                      allergy
# ****************************************************
#replacing Yes with yes and No with no
# Define the columns to modify
columns_to_modify <- c("AtopicEczema_6m", "AtopicEczema_14m",
                       "TotalEczema_6m", "TotalEczema_14m")

# Replace "Yes" with "yes" and "No" with "no" in specified columns
quest <- quest %>%
  mutate(across(all_of(columns_to_modify), 
                ~ recode(., "Yes" = "yes", "No" = "no")))
rm(columns_to_modify)

# ****************************************************
#                      probiotics
# ****************************************************
#replacing Yes with yes and No with no
# Define the columns to modify
columns_to_modify <- c("probiotics_1w", "probiotics_4w", "probiotics_8w",
                       "probiotics_4m", "probiotics_5m", "probiotics_6m",
                       "probiotics_9m", "probiotics_11m", "probiotics_14m")

# Replace "Yes" with "yes" and "No" with "no" in specified columns
quest <- quest %>%
  mutate(across(all_of(columns_to_modify), 
                ~ recode(., "Yes" = "yes", "No" = "no")))


# filling NAs in variables and cleaning ====
# ****************************************************
#                      pets
# ****************************************************
quest$furrypet_indoors_6m_MM <- quest$furrypet_indoors_6m
quest$furrypet_indoors_14m_MM <- quest$furrypet_indoors_14m

plyr::count(quest$Fur_pets_pregn)
plyr::count(quest$furrypet_indoors_6m) 
plyr::count(quest$furrypet_indoors_14m) 

# if preg no and 6m no -> 14m is no
quest <- quest %>%
  mutate(furrypet_indoors_14m = if_else(
    condition =  is.na(furrypet_indoors_14m)&.$Fur_pets_pregn == "No"&.$furrypet_indoors_6m == "No (furry) pets", 
    true = "No (furry) pets", 
    false = furrypet_indoors_14m))

# if preg no and 14m no -> 6m is no
quest <- quest %>%
  mutate(furrypet_indoors_6m = if_else(
    condition =  is.na(furrypet_indoors_6m)&.$Fur_pets_pregn == "No"&.$furrypet_indoors_14m == "No (furry) pets", 
    true = "No (furry) pets", 
    false = furrypet_indoors_6m))

# if preg yes and 14m yes -> 6m is yes home, preg has 2 options for pets at home
quest <- quest %>%
  mutate(furrypet_indoors_6m = if_else(
    condition = is.na(furrypet_indoors_6m) & str_detect(Fur_pets_pregn, "furry pets in pregnancy") & furrypet_indoors_14m == "Yes, in the home (no daycare or no furry pet at daycare)",  
    true = "Yes, in the home (no daycare or no furry pet at daycare)",  
    false = furrypet_indoors_6m))

# if preg yes and 6m yes -> 14m is yes home, preg has 2 options for pets at home
quest <- quest %>%
  mutate(furrypet_indoors_14m = if_else(
    condition = is.na(furrypet_indoors_14m) & str_detect(Fur_pets_pregn, "furry pets in pregnancy") & furrypet_indoors_6m == "Yes, in the home (no daycare or no furry pet at daycare)",  
    true = "Yes, in the home (no daycare or no furry pet at daycare)",  
    false = furrypet_indoors_14m))

# if preg yes and 6m yes -> 14m is yes home and daycare, preg has 2 options for pets at home
quest <- quest %>%
  mutate(furrypet_indoors_14m = if_else(
    condition = is.na(furrypet_indoors_14m) & str_detect(Fur_pets_pregn, "furry pets in pregnancy") & furrypet_indoors_6m == "Yes, in the home and at daycare",  
    true =  "Yes, in the home and at daycare", 
    false = furrypet_indoors_14m))

# if preg yes and 14m yes -> 6m is yes home and daycare, preg has 2 options for pets at home
quest <- quest %>%
  mutate(furrypet_indoors_6m = if_else(
    condition = is.na(furrypet_indoors_6m) & str_detect(Fur_pets_pregn, "furry pets in pregnancy") & furrypet_indoors_14m == "Yes, in the home and at daycare",  
    true =  "Yes, in the home and at daycare", 
    false = furrypet_indoors_6m))

# if preg no and 6m yes -> 14m is yes daycare
quest <- quest %>%
  mutate(furrypet_indoors_14m = if_else(
    condition =  is.na(furrypet_indoors_14m)&.$Fur_pets_pregn == "No"
    &.$furrypet_indoors_6m == "Yes at daycare (not at home)", 
    true = "Yes at daycare (not at home)", 
    false = furrypet_indoors_14m))

plyr::count(quest$Fur_pets_pregn)
plyr::count(quest$furrypet_indoors_6m) 
plyr::count(quest$furrypet_indoors_14m) 

#cleaning pets at 14 months
quest$furrypet_indoors_14m_MM_NA_filled <- quest$furrypet_indoors_14m
quest$furrypet_indoors_14m[grepl("in the home",quest$furrypet_indoors_14m)]<-"at home"
quest$furrypet_indoors_14m[grepl("No",quest$furrypet_indoors_14m)]<-"no"
quest$furrypet_indoors_14m[grepl("daycare",quest$furrypet_indoors_14m)]<-"at daycare"
#cleaning pets at 6 months
quest$furrypet_indoors_6m_MM_NA_filled <- quest$furrypet_indoors_6m
quest$furrypet_indoors_6m[grepl("in the home",quest$furrypet_indoors_6m)]<-"at home"
quest$furrypet_indoors_6m[grepl("No",quest$furrypet_indoors_6m)]<-"no"
quest$furrypet_indoors_6m[grepl("daycare",quest$furrypet_indoors_6m)]<-"at daycare"
#cleaning pets at pregnancy
quest$Fur_pets_pregn_MM <- quest$Fur_pets_pregn
quest$Fur_pets_pregn_yn[grepl("Yes",quest$Fur_pets_pregn)]<-"yes"
quest$Fur_pets_pregn_yn[grepl("No",quest$Fur_pets_pregn)]<-"no"
quest <- quest %>% select(-c(Fur_pets_pregn))

plyr::count(quest$Fur_pets_pregn_yn)
plyr::count(quest$furrypet_indoors_6m) 
plyr::count(quest$furrypet_indoors_14m) 


#only pets at home 6m
quest <- quest %>%
  mutate(furrypet_home_yn_6m = if_else(
    condition =  .$furrypet_indoors_6m_MM_NA_filled == "No (furry) pets", 
    true = "no", 
    false = if_else(condition =  .$furrypet_indoors_6m_MM_NA_filled == "Yes at daycare (not at home)", 
                    true = "no", 
                    false = if_else(condition =  .$furrypet_indoors_6m_MM_NA_filled == "Yes, in the home and at daycare", 
                                    true = "yes", 
                                    false = if_else(condition =  .$furrypet_indoors_6m_MM_NA_filled == "Yes, in the home (no daycare or no furry pet at daycare)", 
                                                    true = "yes", 
                                                    false = furrypet_indoors_6m_MM_NA_filled)))))
#only at home at 14 months
quest <- quest %>%
  mutate(furrypet_home_yn_14m = if_else(
    condition =  .$furrypet_indoors_14m_MM_NA_filled == "No (furry) pets", 
    true = "no", 
    false = if_else(condition =  .$furrypet_indoors_14m_MM_NA_filled == "Yes at daycare (not at home)", 
                    true = "no", 
                    false = if_else(condition =  .$furrypet_indoors_14m_MM_NA_filled == "Yes, in the home and at daycare", 
                                    true = "yes", 
                                    false = if_else(condition =  .$furrypet_indoors_14m_MM_NA_filled == "Yes, in the home (no daycare or no furry pet at daycare)", 
                                                    true = "yes", 
                                                    false = furrypet_indoors_14m_MM_NA_filled)))))
#if no at birth, no pets at 6m
quest <- quest %>%
  mutate(furrypet_home_yn_6m = if_else(
    condition = is.na(furrypet_home_yn_6m)&.$Fur_pets_pregn_yn == "no",
    true = "no",
    false = if_else(is.na(furrypet_home_yn_6m)&.$Fur_pets_pregn_yn == "yes"&is.na(furrypet_home_yn_14m),
                    true = "yes",
                    false = furrypet_home_yn_6m)))

# pet type dog, cat
quest[quest$Family_ID=="1083", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1102", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1126", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1128", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1147", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1150", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1172", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1218", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1223", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1226", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1227", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1041", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1042", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1123", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1216", "animal_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1013", "animal_pregnancy"] <- "Hond(en)\nKat(ten)"
quest[quest$Family_ID=="1026", "animal_pregnancy"] <- "Hond(en)\nKippen in de tuin"
quest[quest$Family_ID=="1081", "animal_pregnancy"] <- "Hond(en)\nKat(ten)\nKnaagdier(en)"
quest[quest$Family_ID=="1198", "animal_pregnancy"] <- "Hond(en)\nKnaagdier(en)"
quest[quest$Family_ID=="1209", "animal_pregnancy"] <- "Hond(en)\nKnaagdier(en)"
quest[quest$Family_ID=="1017", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1022", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1030", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1034", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1036", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1048", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1054", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1062", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1072", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1090", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1099", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1115", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1124", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1134", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1146", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1155", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1154", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1173", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1183", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1213", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1221", "animal_pregnancy"] <- "Kat(ten)"
quest[quest$Family_ID=="1214", "animal_pregnancy"] <- "Hond(en)\nKat(ten)"
quest[quest$Family_ID=="1137", "animal_pregnancy"] <- "Hond(en)\nKnaagdier(en)\nkonijn"
quest[quest$Family_ID=="1086", "animal_pregnancy"] <- "Hond(en)\nkonijnen"
quest[quest$Family_ID=="1067", "animal_pregnancy"] <- "Hond(en)\nVogel(s)"
quest[quest$Family_ID=="1082", "animal_pregnancy"] <- "Kat(ten)\nKnaagdier(en)"
quest[quest$Family_ID=="1053", "animal_pregnancy"] <- "Kat(ten)\nKnaagdier(en)"
quest[quest$Family_ID=="1096", "animal_pregnancy"] <- "Kat(ten)\nKnaagdier(en)"
quest[quest$Family_ID=="1001", "animal_pregnancy"] <- "Kat(ten)\nVis"
quest[quest$Family_ID=="1003", "animal_pregnancy"] <- "Knaagdier(en)"
quest[quest$Family_ID=="1138", "animal_pregnancy"] <- "Kat(ten)\nKnaagdier(en)"
quest[quest$Family_ID=="1200", "animal_pregnancy"] <- "Knaagdier(en)"
quest[quest$Family_ID=="1215", "animal_pregnancy"] <- "Knaagdier(en)"
quest[quest$Family_ID=="1211", "animal_pregnancy"] <- "Knaagdier(en)"


quest[quest$Family_ID=="1013", "animal_after_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1226", "animal_after_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1227", "animal_after_pregnancy"] <- "Hond(en)"
quest[quest$Family_ID=="1138", "animal_after_pregnancy"] <- "Kat(ten)"

quest[quest$Family_ID=="1223", "pets_type_home_6m"] <- "Hond(en)"
quest[quest$Family_ID=="1227", "pets_type_home_6m"] <- "Hond(en)"

#dogs and cats at pregnancy
quest$dog_home_pregnancy <- quest$animal_pregnancy
quest$dog_home_pregnancy[grepl("Hond",quest$dog_home_pregnancy)]<-"yes"
quest <- quest %>%
  mutate(dog_home_pregnancy = if_else(
    condition = is.na(dog_home_pregnancy)&.$Fur_pets_pregn_yn == "no",
    true = "no",
    false = if_else(condition = .$dog_home_pregnancy == "yes",
                    true = "yes",
                    false = if_else(condition = is.na(dog_home_pregnancy),
                                    true = NA,
                                    false = "no"
                    ))))

quest$cat_home_pregnancy <- quest$animal_pregnancy
quest$cat_home_pregnancy[grepl("Kat",quest$cat_home_pregnancy)]<-"yes"
quest <- quest %>%
  mutate(cat_home_pregnancy = if_else(
    condition = is.na(cat_home_pregnancy)&.$Fur_pets_pregn_yn == "no",
    true = "no",
    false = if_else(condition = .$cat_home_pregnancy == "yes",
                    true = "yes",
                    false = if_else(condition = is.na(cat_home_pregnancy),
                                    true = NA,
                                    false = "no"
                    ))))
#dogs and cats at 6m
quest$dog_home_6m <- quest$pets_type_home_6m
quest$dog_home_6m[grepl("Hond",quest$dog_home_6m)]<-"yes"
quest$dog_home_6m[grepl("hond",quest$dog_home_6m)]<-"yes"
quest <- quest %>%
  mutate(dog_home_6m = if_else(
    condition = is.na(dog_home_6m)&.$furrypet_home_yn_6m == "no",
    true = "no",
    false = if_else(condition = .$dog_home_6m == "yes",
                    true = "yes",
                    false = if_else(condition = is.na(dog_home_6m),
                                    true = NA,
                                    false = "no"
                    ))))
quest$cat_home_6m <- quest$pets_type_home_6m
quest$cat_home_6m[grepl("kat",quest$cat_home_6m)]<-"yes"
quest <- quest %>%
  mutate(cat_home_6m = if_else(
    condition = is.na(cat_home_6m)&.$furrypet_home_yn_6m == "no",
    true = "no",
    false = if_else(condition = .$cat_home_6m == "yes",
                    true = "yes",
                    false = if_else(condition = is.na(cat_home_6m),
                                    true = NA,
                                    false = "no"
                    ))))

# ****************************************************
#                      pregnancy weeks
# ****************************************************
#gestational age
quest$pregn_weeks_MM <- quest$pregn_weeks
quest[quest$Family_ID=="1224", "pregn_weeks"] <- round(mean(quest$pregn_weeks,na.rm = TRUE),digits = 0)

# ****************************************************
#                      birth weight
# ****************************************************
quest$birthweight_MM <- quest$birthweight
quest[quest$Family_ID=="1224", "birthweight"] <- round(mean(quest$birthweight,na.rm = TRUE),digits = 0)

# ****************************************************
#                      maternal weight gain
# ****************************************************
quest$mat_weightgain_MM <- quest$mat_weightgain
quest[quest$Family_ID=="1224", "mat_weightgain"] <- round(mean(quest$mat_weightgain,na.rm = TRUE),digits = 0)

# ****************************************************
#                      solids
# ****************************************************
#age when solids were introduced
quest$ageintrosolids_MM <- quest$ageintrosolids

quest[quest$Family_ID=="1088", "ageintrosolids"] <- 16 #as mother put 4 weeks as the earliest solid food introduction, must be confusion with 4 months, which is 16 weeks

mean <- round(mean(quest$ageintrosolids,na.rm = TRUE),digits = 0)
quest[quest$Family_ID=="1224", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1066", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1119", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1137", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1156", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1165", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1214", "ageintrosolids"] <- mean
quest[quest$Family_ID=="1215", "ageintrosolids"] <- mean

#test <- subset(quest, is.na(ageintrosolids))
#plyr::count(test$Family_ID)

# ****************************************************
#                      siblings
# ****************************************************
#siblings in lucki, those infants that have a sibling in lucki but it was not metagenomicly sequenced -> renamed as no siblings in lucki
quest$siblings_MM <- quest$siblings
quest$siblings[grepl("11",quest$siblings)]<-"no"
quest$siblings[grepl("10",quest$siblings)]<-"no"
#renaming the variables
names(quest)[names(quest) == "twins"] <- "twins_in_lucki"
names(quest)[names(quest) == "siblings"] <- "siblings_in_lucki"
#regrouping older siblings
quest$older_siblings_MM <- quest$older_siblings
quest$older_siblings[grepl("0",quest$older_siblings_MM)]<-"no"
quest$older_siblings[grepl("1",quest$older_siblings_MM)]<-"yes"
quest$older_siblings[grepl("2",quest$older_siblings_MM)]<-"yes"
quest$older_siblings[grepl("3",quest$older_siblings_MM)]<-"yes"
#filling NAs fo 1224 as no metadata is available for them
quest[quest$Family_ID=="1224", "twins_in_lucki"] <- "no"
quest[quest$Family_ID=="1224", "siblings_in_lucki"] <- "no"
#creating a new variable siblings in lucki+twins in lucki
quest <- quest %>%
  mutate(twins.siblings_in_lucki = if_else(
    condition = .$twins_in_lucki == "no"&.$siblings_in_lucki == "no", 
    true = "no", 
    false = "yes"))#creating new variable
#check_output <- quest %>% group_by(twins.siblings_in_lucki,twins_in_lucki,siblings_in_lucki) %>% tally() 

# ****************************************************
#                      breastfeeding and formula
# ****************************************************
test <- quest %>% select(Family_ID, matches(c("breastfeeding","formulafeeding")))

#if BF or FF were NA, but before and after info is available -> filling the values that was before AND after
# week 1
quest$breastfeeding_1_2w_MM <- quest$breastfeeding_1_2w
quest$formulafeeding_1_2w_MM <- quest$formulafeeding_1_2w
test <- quest %>% filter(is.na(breastfeeding_1_2w))
quest <- quest %>%
  mutate(breastfeeding_1_2w = if_else(
    is.na(breastfeeding_1_2w) & breastfeeding_birth == breastfeeding_4w,  # Condition: breastfeeding_1_2w is NA and breastfeeding_birth equals breastfeeding_4w
    breastfeeding_birth,  # If condition is TRUE, compute breastfeeding_1_2w as breastfeeding_birth
    breastfeeding_1_2w))   # Otherwise, leave NA

quest <- quest %>%
  mutate(formulafeeding_1_2w = if_else(
    is.na(formulafeeding_1_2w) & formulafeeding_birth == formulafeeding_4w, 
    formulafeeding_birth, 
    formulafeeding_1_2w))

plyr::count(quest$breastfeeding_1_2w)
plyr::count(quest$breastfeeding_1_2w)

# week 4
plyr::count(quest$formulafeeding_4w)
plyr::count(quest$breastfeeding_4w)
quest$formulafeeding_4w_MM <- quest$formulafeeding_4w
quest$breastfeeding_4w_MM <- quest$breastfeeding_4w
quest <- quest %>%
  mutate(formulafeeding_4w = if_else(
    is.na(formulafeeding_4w) & formulafeeding_1_2w == formulafeeding_8w, 
    formulafeeding_1_2w,
    formulafeeding_4w))
quest <- quest %>%
  mutate(breastfeeding_4w = if_else(
    is.na(breastfeeding_4w) & breastfeeding_1_2w == breastfeeding_8w, 
    breastfeeding_1_2w,
    breastfeeding_4w))

# week 8
plyr::count(quest$breastfeeding_8w)
plyr::count(quest$formulafeeding_8w)
quest$formulafeeding_8w_MM <- quest$formulafeeding_8w
quest$breastfeeding_8w_MM <- quest$breastfeeding_8w
quest <- quest %>%
  mutate(formulafeeding_8w = if_else(
    is.na(formulafeeding_8w) & formulafeeding_4w == formulafeeding_4m, 
    formulafeeding_4w,
    formulafeeding_8w))

quest <- quest %>%
  mutate(breastfeeding_8w = if_else(
    is.na(breastfeeding_8w) & breastfeeding_4w == breastfeeding_4m, 
    breastfeeding_4w,
    breastfeeding_8w))

# month 4
quest$breastfeeding_4m_MM <- quest$breastfeeding_4m
quest$formulafeeding_4m_MM <- quest$formulafeeding_4m
plyr::count(quest$breastfeeding_4m)
plyr::count(quest$formulafeeding_4m)
quest <- quest %>%
  mutate(formulafeeding_4m = if_else(
    is.na(formulafeeding_4m) & formulafeeding_8w == formulafeeding_5m, 
    formulafeeding_8w,
    formulafeeding_4m))

quest <- quest %>%
  mutate(breastfeeding_4m = if_else(
    is.na(breastfeeding_4m) & breastfeeding_8w == breastfeeding_5m, 
    breastfeeding_8w,
    breastfeeding_4m))

# month 5
quest$breastfeeding_5m_MM <- quest$breastfeeding_5m
quest$formulafeeding_5m_MM <- quest$formulafeeding_5m
plyr::count(quest$breastfeeding_5m)
plyr::count(quest$formulafeeding_5m)

quest <- quest %>%
  mutate(formulafeeding_5m = if_else(
    is.na(formulafeeding_5m) & formulafeeding_4m == formulafeeding_6m, 
    formulafeeding_4m,
    formulafeeding_5m))
quest <- quest %>%
  mutate(breastfeeding_5m = if_else(
    is.na(breastfeeding_5m) & breastfeeding_4m == breastfeeding_6m, 
    breastfeeding_4m,
    breastfeeding_5m))

# month 6
quest$breastfeeding_6m_MM <- quest$breastfeeding_6m
quest$formulafeeding_6m_MM <- quest$formulafeeding_6m
plyr::count(quest$breastfeeding_6m)
plyr::count(quest$formulafeeding_6m)

quest <- quest %>%
  mutate(formulafeeding_6m = if_else(
    is.na(formulafeeding_6m) & formulafeeding_5m == formulafeeding_9m, 
    formulafeeding_5m,
    formulafeeding_6m))
quest <- quest %>%
  mutate(breastfeeding_6m = if_else(
    is.na(breastfeeding_6m) & breastfeeding_5m == breastfeeding_9m, 
    breastfeeding_5m,
    breastfeeding_6m))

# month 9
quest$breastfeeding_9m_MM <- quest$breastfeeding_9m
quest$formulafeeding_9m_MM <- quest$formulafeeding_9m
plyr::count(quest$breastfeeding_9m)
plyr::count(quest$formulafeeding_9m)
quest <- quest %>%
  mutate(formulafeeding_9m = if_else(
    is.na(formulafeeding_9m) & formulafeeding_6m == formulafeeding_11m, 
    formulafeeding_6m,
    formulafeeding_9m))
quest <- quest %>%
  mutate(breastfeeding_9m = if_else(
    is.na(breastfeeding_9m) & breastfeeding_6m == breastfeeding_11m, 
    breastfeeding_6m,
    breastfeeding_9m))

# month 11
quest$breastfeeding_11m_MM <- quest$breastfeeding_11m
quest$formulafeeding_11m_MM <- quest$formulafeeding_11m
plyr::count(quest$breastfeeding_11m)
plyr::count(quest$formulafeeding_11m)
quest <- quest %>%
  mutate(formulafeeding_11m = if_else(
    is.na(formulafeeding_11m) & formulafeeding_9m == formulafeeding_14m, 
    formulafeeding_9m,
    formulafeeding_11m))
quest <- quest %>%
  mutate(breastfeeding_11m = if_else(
    is.na(breastfeeding_11m) & breastfeeding_9m == breastfeeding_14m, 
    breastfeeding_9m,
    breastfeeding_11m))

# 14 months
quest$breastfeeding_14m_MM <- quest$breastfeeding_14m
quest$formulafeeding_14m_MM <- quest$formulafeeding_14m

#if BF was NA, but before  is "no" -> filling the values with no
test <- quest %>% select(Family_ID, matches(c("breastfeeding","formulafeeding"))) %>% filter(if_any(everything(), is.na))
quest <- quest %>%
  mutate(breastfeeding_5m = if_else(
    is.na(breastfeeding_5m) & breastfeeding_4m == "No", 
    breastfeeding_4m,
    breastfeeding_5m))
quest <- quest %>%
  mutate(breastfeeding_6m = if_else(
    is.na(breastfeeding_6m) & breastfeeding_5m == "No", 
    breastfeeding_5m,
    breastfeeding_6m))
quest <- quest %>%
  mutate(breastfeeding_9m = if_else(
    is.na(breastfeeding_9m) & breastfeeding_6m == "No", 
    breastfeeding_6m,
    breastfeeding_9m))
quest <- quest %>%
  mutate(breastfeeding_11m = if_else(
    is.na(breastfeeding_11m) & breastfeeding_9m == "No", 
    breastfeeding_9m,
    breastfeeding_11m))
quest <- quest %>%
  mutate(breastfeeding_14m = if_else(
    is.na(breastfeeding_14m) & breastfeeding_11m == "No", 
    breastfeeding_11m,
    breastfeeding_14m))

#if FF was NA, but before  is "no" -> filling the values with no only for 14 months
test <- quest %>% select(Family_ID, matches(c("breastfeeding","formulafeeding"))) %>% filter(if_any(everything(), is.na))
quest <- quest %>%
  mutate(formulafeeding_14m = if_else(
    is.na(formulafeeding_14m) & formulafeeding_11m == "No", 
    formulafeeding_11m,
    formulafeeding_14m))

# the rest is filled manually by logically deducting, only for thos that have samples available
test <- quest %>% select(Family_ID, matches(c("breastfeeding","formulafeeding"))) %>% filter(if_any(everything(), is.na))
quest[quest$Family_ID=="1172", "formulafeeding_1_2w"] <- "No"
quest[quest$Family_ID=="1027", "formulafeeding_4w"] <- "Yes" #net
quest[quest$Family_ID=="1088", "formulafeeding_4w"] <- "Yes"
quest[quest$Family_ID=="1156", "breastfeeding_8w"] <- "Yes"
quest[quest$Family_ID=="1156", "formulafeeding_8w"] <- "Yes"
quest[quest$Family_ID=="1027", "formulafeeding_8w"] <- "Yes"#net
quest[quest$Family_ID=="1214", "breastfeeding_8w"] <- "Yes"
quest[quest$Family_ID=="1214", "formulafeeding_8w"] <- "No"
quest[quest$Family_ID=="1198", "breastfeeding_4m"] <- "Yes"
quest[quest$Family_ID=="1198", "formulafeeding_4m"] <- "Yes"
quest[quest$Family_ID=="1034", "breastfeeding_5m"] <- "Yes"
quest[quest$Family_ID=="1034", "formulafeeding_5m"] <- "No"
quest[quest$Family_ID=="1077", "breastfeeding_5m"] <- "Yes"
quest[quest$Family_ID=="1077", "formulafeeding_5m"] <- "No"
quest[quest$Family_ID=="1225", "breastfeeding_5m"] <- "No"
quest[quest$Family_ID=="1225", "formulafeeding_5m"] <- "Yes"
quest[quest$Family_ID=="1034", "breastfeeding_6m"] <- "Yes"
quest[quest$Family_ID=="1034", "formulafeeding_6m"] <- "No"
quest[quest$Family_ID=="1173", "formulafeeding_6m"] <- "Yes"
quest[quest$Family_ID=="1117", "formulafeeding_9m"] <- "Yes"
quest[quest$Family_ID=="1194", "formulafeeding_11m"] <- "No"
quest[quest$Family_ID=="1115", "breastfeeding_14m"] <- "No"
quest[quest$Family_ID=="1115", "formulafeeding_14m"] <- "No"
quest[quest$Family_ID=="1155", "breastfeeding_14m"] <- "No"
quest[quest$Family_ID=="1194", "formulafeeding_14m"] <- "No"

#replacing all No with no and Yes with yes
quest <- quest %>% mutate(across(everything(), ~ ifelse(. == "Yes", "yes", ifelse(. == "No", "no", .))))

# merging fecal_premeta and questioners ====
meta_not_fixed = full_join(premeta_fecal,quest, by = "Family_ID",multiple = "all")

#removing extra columns
meta <- meta_not_fixed %>% select(-c(twins,siblings))
rm(premeta_fecal,quest)

# creating new variables in meta ====
# calculating age in weeks and days
meta <- meta %>%
  mutate(Age_weeks = if_else(.$Age_individual == "1-2 weeks", "2", 
                             if_else(.$Age_individual == "4 weeks", "4",
                                     if_else(.$Age_individual == "8 weeks", "8",
                                             if_else(.$Age_individual == "4 months", "16",
                                                     if_else(.$Age_individual == "5 months", "20",
                                                             if_else(.$Age_individual == "6 months", "24",
                                                                     if_else(.$Age_individual == "9 months", "36",
                                                                             if_else(.$Age_individual == "11 months", "44",
                                                                                     if_else(.$Age_individual == "14 months", "56",
                                                                                             NA))))))))))
meta$Age_weeks <- as.numeric(meta$Age_weeks)
meta$Age_days = ifelse(meta$Age_individual == "1-2 weeks",meta$Age_weeks*7,
                       ifelse(meta$Age_individual == "4 weeks", meta$Age_weeks*7,
                              ifelse(meta$Age_individual == "8 weeks", meta$Age_weeks*7,
                                     ifelse(meta$Age_individual == "4 months", meta$Age_weeks*7,
                                            ifelse(meta$Age_individual == "5 months", meta$Age_weeks*7,
                                                   ifelse(meta$Age_individual == "6 months", meta$Age_weeks*7,
                                                          ifelse(meta$Age_individual == "9 months", meta$Age_weeks*7,
                                                                 ifelse(meta$Age_individual == "11 months", meta$Age_weeks*7,
                                                                        ifelse(meta$Age_individual == "14 months", meta$Age_weeks*7,
                                                                               NA)))))))))
meta$Age_days <- as.numeric(meta$Age_days)

# turning dates into date format
meta$DOB_infant <- as.Date(meta$DOB_infant, "%Y-%m-%d")
meta$Sample_date <- as.Date(meta$Sample_date, "%Y-%m-%d")

# dealing with dates ====
#adding infant sample dates for 6 missing families from pdf questionnaires from uni L drive
meta[meta$Sample_aliquote %in% "1222.1.1", "Sample_date"] <- dmy("15/01/2021") 
meta[meta$Sample_aliquote %in% "1225.1.1", "Sample_date"] <- dmy("24/03/2021") 
meta[meta$Sample_aliquote %in% "1226.1.1", "Sample_date"] <- dmy("06/05/2021")
meta[meta$Sample_aliquote %in% "1227.1.1", "Sample_date"] <- dmy("28/05/2021")
meta[meta$Sample_aliquote %in% "1228.1.1", "Sample_date"] <- dmy("04/12/2021")

meta[meta$Sample_aliquote %in% "1222.1.2", "Sample_date"] <- dmy("26/01/2021") 
meta[meta$Sample_aliquote %in% "1223.1.2", "Sample_date"] <- dmy("07/12/2020") 
meta[meta$Sample_aliquote %in% "1225.1.2", "Sample_date"] <- dmy("08/04/2021") 
meta[meta$Sample_aliquote %in% "1226.1.2", "Sample_date"] <- dmy("06/05/2021")
meta[meta$Sample_aliquote %in% "1227.1.2", "Sample_date"] <- dmy("06/06/2021")
meta[meta$Sample_aliquote %in% "1228.1.2", "Sample_date"] <- dmy("15/12/2021")

meta[meta$Sample_aliquote %in% "1222.1.3", "Sample_date"] <- dmy("22/02/2021") 
meta[meta$Sample_aliquote %in% "1223.1.3", "Sample_date"] <- dmy("04/01/2021") 
meta[meta$Sample_aliquote %in% "1225.1.3", "Sample_date"] <- dmy("03/05/2021") 
meta[meta$Sample_aliquote %in% "1226.1.3", "Sample_date"] <- dmy("01/06/2021")
meta[meta$Sample_aliquote %in% "1227.1.3", "Sample_date"] <- dmy("12/07/2021")
meta[meta$Sample_aliquote %in% "1228.1.3", "Sample_date"] <- dmy("16/01/2022")

meta[meta$Sample_aliquote %in% "1222.1.4", "Sample_date"] <- dmy("24/04/2021") 
meta[meta$Sample_aliquote %in% "1223.1.4", "Sample_date"] <- dmy("13/03/2021") 
meta[meta$Sample_aliquote %in% "1226.1.4", "Sample_date"] <- dmy("07/08/2021")
meta[meta$Sample_aliquote %in% "1227.1.4", "Sample_date"] <- dmy("08/09/2021")
meta[meta$Sample_aliquote %in% "1228.1.4", "Sample_date"] <- dmy("13/04/2022")

meta[meta$Sample_aliquote %in% "1222.1.5", "Sample_date"] <- dmy("30/05/2021") 
meta[meta$Sample_aliquote %in% "1223.1.5", "Sample_date"] <- dmy("12/04/2021")
meta[meta$Sample_aliquote %in% "1226.1.5", "Sample_date"] <- dmy("07/09/2021")
meta[meta$Sample_aliquote %in% "1228.1.5", "Sample_date"] <- dmy("14/04/2022")

meta[meta$Sample_aliquote %in% "1222.1.6", "Sample_date"] <- dmy("30/06/2021") 
meta[meta$Sample_aliquote %in% "1223.1.6", "Sample_date"] <- dmy("16/05/2021")

meta[meta$Sample_aliquote %in% "1222.1.7", "Sample_date"] <- dmy("30/09/2021") 
meta[meta$Sample_aliquote %in% "1223.1.7", "Sample_date"] <- dmy("10/08/2021")

meta[meta$Sample_aliquote %in% "1222.1.8", "Sample_date"] <- dmy("01/12/2021") 

meta[meta$Sample_aliquote %in% "1222.1.9", "Sample_date"] <- dmy("05/03/2022")

meta[meta$Family_ID  %in% "1224", "DOB_infant"] <- dmy("19/12/2020") #adding DOB for infant 1224 from ldot

# calculating theoretical dates not working!
meta <- meta %>% 
  mutate(sample_date_theoretical = DOB_infant + Age_days)

meta$Sample_date_filled <- meta$Sample_date
#Sample_date_filled if does not have real date, will be filled with theoretical values
meta <- meta %>% 
  mutate(Sample_date_filled = if_else(
    condition = is.na(Sample_date_filled),
    true = sample_date_theoretical,
    false = Sample_date_filled
  ))

# creating covid parameter =====
meta <- meta %>%
  mutate(covid_270220 = if_else(.$Sample_date_filled >= "2020-02-27", "yes", "no"))
meta <- meta %>% 
  mutate(covid_270220 = if_else(
    condition = is.na(covid_270220),
    true = if_else(.$DOB_infant >= "2020-02-27", "yes", "no"),
    false = covid_270220
  ))

# checking how milk and formula are after filling ====

#test <- subset(meta, is.na(breastfeeding_14m))
#test2 <- subset(test, Age_individual == "14 months") 
#plyr::count(test2$Family_ID)


#test <- subset(meta, is.na(breastfeeding_14m_MM))
#test2 <- subset(test, Age_individual == "14 months") 
#plyr::count(test2$Family_ID)

# checking that merged was ok from different dataframes ====
#test <- meta %>% select(Family_ID,DOB_check, DOB_infant)
#test$comapre_DOB = ifelse(test$DOB_infant == test$DOB_check, 'equal','Not equal') # DOB are equal
#rm(test)

# adding MS diet clusters to the data ====
samples389$used_in_MS_analysis <- 'yes'
names(fam112)[names(fam112) == "kindcode"] <- "Family_ID"
names(fam112)[names(fam112) == "class"] <- "diet_class"
names(samples389)[names(samples389) == "kindcode"] <- "Sample_aliquote"
#selecting only those two columns
fam112 <- fam112 %>% select(Family_ID, diet_class)
samples389 <- samples389 %>% select(Sample_aliquote, used_in_MS_analysis)
diet <- full_join(meta,fam112, by = "Family_ID")
meta_diet <- full_join(diet,samples389, by = "Sample_aliquote")
#checking merge
plyr::count(meta_diet$diet_class)
uniq = meta_diet[!duplicated(meta_diet$Family_ID),]
plyr::count(uniq$Family_ID)
rm(meta,meta_not_fixed,diet)
#changing diet class to binary factor
meta_diet <- meta_diet %>%
  mutate(diet_class_1 = case_when(
    diet_class == "1"  ~ "yes",
    diet_class == "2"  ~ "no",
    diet_class == "3"  ~ "no",
    is.na(diet_class) ~ NA
  ))
meta_diet <- meta_diet %>%
  mutate(diet_class_2 = case_when(
    diet_class == "1"  ~ "no",
    diet_class == "2"  ~ "yes",
    diet_class == "3"  ~ "no",
    is.na(diet_class) ~ NA
  ))
meta_diet <- meta_diet %>%
  mutate(diet_class_3 = case_when(
    diet_class == "1"  ~ "no",
    diet_class == "2"  ~ "no",
    diet_class == "3"  ~ "yes",
    is.na(diet_class) ~ NA
  ))

plyr::count(meta_diet$diet_class_3)

# writing the table ====
#write.table (meta_diet, "S1.4_full_meta.txt", row.names=FALSE,sep = "\t")
