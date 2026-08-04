# libraries ====
library(tidyverse)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ps_diet <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
freezer <- read.delim("Input/7_Diepvrieslocaties_LucKi Gut_14-01_2025.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
diet_weight <- read.delim("Input/12_4mon_scfa_fecal_weight.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
non_diet_weight <- read.delim("Input/13_non_diet_weight.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# selecting samples codes from phyloseq ====
meta <- as.data.frame(as.matrix(ps_diet@sam_data))
meta <- meta %>% filter(Age_individual == "4 months")
families_diet_4mon <- unique(sort(meta$Family_ID))

rm(meta,ps_diet)

# freezer data cleaning ====
#filtering for 4 months
freezer_4mon <- freezer %>% filter(Time.point == "4 months")
#selecting needed columns
freezer_4mon <- freezer_4mon %>% select(all_of(c("SampleID","Sample.nr", "Sample.weight..mg.", "Box", "Place","Comments.and.purpose","CB.child.code")))

#renaming columns
names(freezer_4mon)[names(freezer_4mon) == "Sample.nr"] <- "aliquot"
names(freezer_4mon)[names(freezer_4mon) == "Sample.weight..mg."] <- "original_weight"
names(freezer_4mon)[names(freezer_4mon) == "Comments.and.purpose"] <- "comments_upon_1st_freezing"
names(freezer_4mon)[names(freezer_4mon) == "CB.child.code"] <- "comments_upon_1st_freezing2"
freezer_4mon$aliquotID <- paste(freezer_4mon$SampleID, freezer_4mon$aliquot, sep = "_")

#taking family codes out 
freezer_4mon$backup <- freezer_4mon$SampleID
freezer_4mon <- freezer_4mon %>% separate("backup",into = c("Family_ID", "musor1", "musor2"),sep = "\\.")
freezer_4mon <- freezer_4mon %>% select(!(all_of(c("musor1", "musor2"))))

#taking all family codes
families_all <- unique(sort(freezer_4mon$Family_ID))

#filtering into diet and non diet groups
freezer_4mon_diet <- freezer_4mon %>% filter(Family_ID %in% families_diet_4mon)
freezer_4mon_nondiet <- freezer_4mon %>% filter(!(Family_ID %in% families_diet_4mon))

samples_4mon_nondiet <- unique(sort(freezer_4mon_nondiet$SampleID))

#sanity check
tp <- unique(sort(freezer_4mon_diet$Family_ID)) #should be equal to families_diet_4mon  
tp <- unique(sort(freezer_4mon_nondiet$Family_ID)) #should be equal to 50 or families_all-families_diet_4mon

rm(freezer,freezer_4mon,tp,families_diet_4mon)

# cleaning weight data ====
diet_weight$diet_paper <- "yes"
diet_weight$aliquotID <- paste(diet_weight$SampleID, diet_weight$Sample.nr, sep = "_")
names(diet_weight)[names(diet_weight) == "weight_sample.tube_2026"] <- "weight_sample_tube_2026"
colnames(diet_weight)
diet_weight <- diet_weight %>% select(all_of(c("SampleID","aliquotID","weight_sample_tube_2026","cap_markings","sample_status","diet_paper", "Place")))

non_diet_weight$diet_paper <- "no"
colnames(non_diet_weight)
non_diet_weight <- non_diet_weight %>% select(all_of(c("SampleID","aliquotID","weight_sample_tube_2026","cap_markings","sample_status","diet_paper", "Place")))

# merging data together ====
freezer <- rbind(freezer_4mon_diet, freezer_4mon_nondiet)
test <- plyr::count(freezer$aliquotID)
test <- test %>% filter(freq > 1) # 1157.1.4_1 is a problem

weight <- rbind(diet_weight,non_diet_weight)
test2 <- plyr::count(weight$aliquotID)
test2 <- test2 %>% filter(freq > 1)  # 1157.1.4_1 is a problem

# removing problematic sample
place <- freezer %>% dplyr::filter(!(is.na(comments_upon_1st_freezing2)))
place <- place$Place #should be 11

freezer <- freezer %>% dplyr::filter(is.na(comments_upon_1st_freezing2))
weight <- weight %>% dplyr::filter(!(aliquotID == "1157.1.4_1" & Place == place))

all <- full_join(weight,freezer, by = "aliquotID")

rm(place, test, test2, freezer_4mon_nondiet, freezer_4mon_diet, diet_weight,non_diet_weight, freezer, weight)

# sanity check 
colnames(all)
all(all$Place.x == all$Place.y, na.rm = FALSE) # should be TRUE
all(all$SampleID.x == all$SampleID.y, na.rm = FALSE) # should be TRUE

all <- all %>%
  select(-Place.y, -SampleID.y) %>%
  rename(
    Place = Place.x,
    SampleID = SampleID.x)

# creating new columns and cleaning new table ====
#adding info about weight
average_empty_tube_weight <- mean(1.5223,1.5346,1.5382,1.5276,1.5309,1.5334,1.5341,1.5329,1.5363,1.5307)
all$empty_tube_g <- as.numeric(average_empty_tube_weight)
all$fecal_weight_2026 <- as.numeric(all$weight_sample_tube_2026) - all$empty_tube_g

# replacing info with NA if sample is absent
colnames(all)
plyr::count(all$cap_markings)
all <- all %>%
  mutate(
    cap_markings = if_else(
      sample_status == "absent" & cap_markings == "no",
      NA_character_,
      cap_markings
    ),
    weight_sample_tube_2026 = if_else(
      sample_status == "absent" & weight_sample_tube_2026 == "no",
      NA_character_,
      weight_sample_tube_2026))

rm(average_empty_tube_weight)

# cleaning comment section 
plyr::count(all$comments_upon_1st_freezing)
test <- all %>% filter(grepl("16S", comments_upon_1st_freezing, ignore.case = TRUE))%>% filter(grepl("Ontdooid", comments_upon_1st_freezing, ignore.case = TRUE))
all <- all %>%
  mutate(
    comments_upon_1st_freezing = recode(
      comments_upon_1st_freezing,
      "ontdooid voor NMR spec" = "Ontdooid voor NMR spec",
      "Opgebruikt MMHP" = "Opgebruikt voor MMHP",
      "further aliquoted" = "Further aliquoted",
      "Ontdooid aliquoten 16S seq " = "Ontdooid aliquoten 16S seq",
      "16S sequencing" = "Opgebruikt voor 16S seq"))

all <- all %>%
  mutate(
    comments_upon_1st_freezing2 = if_else(
      comments_upon_1st_freezing %in% c(
        "1, 3, 4 aliquote do not exist",
        "4 and 2 aliquotes do not exsist"
      ),
      comments_upon_1st_freezing,
      comments_upon_1st_freezing2
    ),
    comments_upon_1st_freezing = if_else(
      comments_upon_1st_freezing %in% c(
        "1, 3, 4 aliquote do not exist",
        "4 and 2 aliquotes do not exsist"
      ),
      NA_character_,
      comments_upon_1st_freezing
    )
  )

rm(test)

# adding about buffer in the samples
plyr::count(all$comments_upon_1st_freezing)
samples_with_buffer <- c("1228.1.4_1","1236.1.4_1","1232.1.4_1","1239.1.4_1","1235.1.4_1","1240.1.4_1")
all <- all %>% mutate(comments_upon_1st_freezing = if_else(
  is.na(comments_upon_1st_freezing) & aliquotID %in% samples_with_buffer,
  "buffer added before freezing",
  comments_upon_1st_freezing))

test <- all %>% filter(aliquotID %in% samples_with_buffer)
test <- test %>% filter(comments_upon_1st_freezing == "Opgebruikt voor 16S seq")
test <- test$aliquotID

all <- all %>%
  mutate(
    comments_upon_1st_freezing = if_else(
      aliquotID %in% test,
      if_else(
        is.na(comments_upon_1st_freezing),
        "buffer added before freezing",
        paste0(comments_upon_1st_freezing, " ; buffer added before freezing")
      ),
      comments_upon_1st_freezing))

rm(test,samples_with_buffer)

# adding counts of aliquots present ====
plyr::count(all$aliquot)
all <- all %>%
  group_by(Family_ID) %>%
  mutate(aliquot_n_theory = n()) %>%
  ungroup()

plyr::count(all$sample_status)
all <- all %>%
  group_by(Family_ID) %>%
  mutate(aliquot_n_present = sum(sample_status == "present")) %>%
  ungroup()

plyr::count(all$aliquot)
all <- all %>%
  group_by(Family_ID) %>%
  mutate(aliquot_n_ali1_4 = sum(sample_status == "present" & aliquot < 5)) %>% 
  ungroup()

plyr::count(all$cap_markings)
all <- all %>%
  group_by(Family_ID) %>%
  mutate(aliquot_n_noMarkings = sum(sample_status == "present" & aliquot < 5 & cap_markings == "no")) %>% 
  ungroup()

test <- all %>% filter(sample_status == "present" & aliquot < 5 & cap_markings == "no")
plyr::count(test$comments_upon_1st_freezing)
all <- all %>%
  group_by(Family_ID) %>%
  mutate(aliquot_n_noComments = sum(sample_status == "present" & aliquot < 5 & cap_markings == "no" &
                                      is.na(comments_upon_1st_freezing))) %>% 
  ungroup()

test <- all %>% filter(sample_status == "present" & aliquot < 5 & cap_markings == "no"&
                         is.na(comments_upon_1st_freezing))

# first preselection ====
# samples must be present, have to have no markings on the cap, and no comments in the freezer map, also only 
#looking at aliquots 1-4.
rm(test)

all$preselection <- NA

#preselection 1, samples with 3 aliquotes
all <- all %>%
  group_by(Family_ID) %>%
  mutate(has_match = any(
      is.na(comments_upon_1st_freezing) &
        sample_status == "present" &
        aliquot_n_noComments == 3),
    min_aliquot = if_else(
      has_match,
      min(aliquot[is.na(comments_upon_1st_freezing) &
          sample_status == "present" &
          aliquot_n_noComments == 3
      ], na.rm = TRUE),NA_real_),
    preselection = case_when(
      has_match &
        is.na(comments_upon_1st_freezing) &
        sample_status == "present" &
        aliquot_n_noComments == 3 &
        aliquot == min_aliquot ~ "yes",
      has_match ~ "no",
      TRUE ~ preselection)) %>%
  select(-has_match, -min_aliquot) %>%
  ungroup()

#sanity check
test <- all %>% filter(aliquot_n_noComments == 3) %>% select(all_of(c("Family_ID", "aliquotID","aliquot_n_noComments",
                                                                      "preselection")))

fam1 <- unique(sort(test$Family_ID))
test2 <- all %>% filter(Family_ID %in% fam1) # should have the same obs as test

selected_families <- fam1

#preselection 2, samples with 2 aliquotes
all <- all %>%
  group_by(Family_ID) %>%
  mutate(has_match = any(
    is.na(comments_upon_1st_freezing) &
      sample_status == "present" &
      aliquot_n_noComments == 2),
    min_aliquot = if_else(
      has_match,
      min(aliquot[is.na(comments_upon_1st_freezing) &
                    sample_status == "present" &
                    aliquot_n_noComments == 2
      ], na.rm = TRUE),NA_real_),
    preselection = case_when(
      has_match &
        is.na(comments_upon_1st_freezing) &
        sample_status == "present" &
        aliquot_n_noComments == 2 &
        aliquot == min_aliquot ~ "yes",
      has_match ~ "no",
      TRUE ~ preselection)) %>%
  select(-has_match, -min_aliquot) %>%
  ungroup()

#sanity check
test <- all %>% filter(aliquot_n_noComments == 2) %>% select(all_of(c("Family_ID", "aliquotID","aliquot_n_noComments",
                                                                      "preselection")))

fam1 <- unique(sort(test$Family_ID))
test2 <- all %>% filter(Family_ID %in% fam1) # should have the same obs as test

selected_families <- c(selected_families, fam1)

#preselection 3, samples with 1 aliquot
all <- all %>%
  group_by(Family_ID) %>%
  mutate(has_match = any(
    is.na(comments_upon_1st_freezing) &
      sample_status == "present" &
      aliquot_n_noComments == 1),
    min_aliquot = if_else(
      has_match,
      min(aliquot[is.na(comments_upon_1st_freezing) &
                    sample_status == "present" &
                    aliquot_n_noComments == 1
      ], na.rm = TRUE),NA_real_),
    preselection = case_when(
      has_match &
        is.na(comments_upon_1st_freezing) &
        sample_status == "present" &
        aliquot_n_noComments == 1 &
        aliquot == min_aliquot ~ "yes, last aliquot",
      has_match ~ "no",
      TRUE ~ preselection)) %>%
  select(-has_match, -min_aliquot) %>%
  ungroup()

#sanity check
test <- all %>% filter(aliquot_n_noComments == 1) %>% select(all_of(c("Family_ID", "aliquotID","aliquot_n_noComments",
                                                                      "preselection")))

fam1 <- unique(sort(test$Family_ID))
test2 <- all %>% filter(Family_ID %in% fam1) # should have the same obs as test

selected_families <- c(selected_families, fam1)

plyr::count(all$preselection)

rm(test,test2,fam1)

# preselection 4
test <- all %>% filter(aliquot > 4) %>% filter(!(Family_ID %in% selected_families))
bigaliquot <- unique(sort(test$Family_ID))
test <- all %>% filter(Family_ID %in% bigaliquot) 
all <- all %>%
  mutate(
    preselection = if_else(
      Family_ID %in% bigaliquot,
      "only 5-7 unreliable aliquots left",
      preselection))

rm(test,bigaliquot)

#preselectino 5
plyr::count(all$preselection)
all$preselection[is.na(all$preselection)] <- "no usable aliquotes left"

# indicating weight ====
all <- all %>%
  mutate(storage_maastricht = case_when(
      is.na(fecal_weight_2026) ~ NA_character_,
      fecal_weight_2026 >= 0.2 ~ "enough",
      fecal_weight_2026 < 0.2 ~ "not enough"))

tempo <- all %>% filter(preselection == "yes, last aliquot") %>% filter(storage_maastricht == "not enough")
plyr::count(tempo$storage_maastricht)

# narrowing down the samples with only 1 sample left ====
rm(tempo,selected_families)
plyr::count(all$preselection)

tempo <- all %>% filter(preselection == "yes, last aliquot")
tp <- unique(sort(tempo$Family_ID))
tempo <- all %>% filter(Family_ID %in% tp)
plyr::count(tempo$aliquot)

# adding columns for extraction ====
all$V_MQ_added_ml <- NA
all$V_of_aliquote_1_mikrol <- NA
all$V_of_aliquote_2_mikrol <- NA

#already extracted sampels
all <- all %>%
  mutate(
    V_MQ_added_ml = case_when(
      aliquotID == "1143.1.4_2" & is.na(V_MQ_added_ml) ~ 0.72,
      aliquotID == "1149.1.4_2" & is.na(V_MQ_added_ml) ~ 0.73,
      TRUE ~ V_MQ_added_ml))

all <- all %>%
  mutate(
    V_of_aliquote_1_mikrol = case_when(
      aliquotID == "1143.1.4_2" & is.na(V_of_aliquote_1_mikrol) ~ 0.3,
      aliquotID == "1149.1.4_2" & is.na(V_of_aliquote_1_mikrol) ~ 0.3,
      TRUE ~ V_of_aliquote_1_mikrol))

all <- all %>%
  mutate(
    V_of_aliquote_2_mikrol = case_when(
      aliquotID == "1143.1.4_2" & is.na(V_of_aliquote_2_mikrol) ~ 0.25,
      aliquotID == "1149.1.4_2" & is.na(V_of_aliquote_2_mikrol) ~ 0.3,
      TRUE ~ V_of_aliquote_2_mikrol))

# saving results ====
#write.table (all, "S1.3_4mon_scfa_freezer_preselection.txt", row.names=FALSE,sep = "\t")