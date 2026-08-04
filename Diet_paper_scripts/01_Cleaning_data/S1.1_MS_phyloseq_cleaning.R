# libraries ====
library(microViz)
library(tidyverse)
library(phyloseq)
library(haven) # for spss
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in phyloseq data and break it into pieces ====
spss_file = as.data.frame(read_sav("Input/1_Mergedfile_0-14m_221130_part1_tm1221_SEL_MS.sav"))
microbial_age_df <- read.csv(file = "Input/2_MS_original_10_microbial_age_df.csv")
ms <- readRDS("Input/3_MS_original_9_phyloseq_object.rds")

otu <- as.data.frame(as(otu_table(ms), "matrix"))
sam <- as(sample_data(ms), "data.frame")
taxa <- as.data.frame(as(tax_table(ms), "matrix"))
rm(ms)

# cleaning sam table ====
names(sam)[names(sam) == "kindcode"] <- "Sample_aliquote"
names(sam)[names(sam) == "tp"] <- "timepoint"
sam$Age_individual <- sam$timepoint
# Replace underscores with spaces in 'Age_individual'
sam$Age_individual <- gsub("m", " months", sam$Age_individual)
sam$Age_individual<- factor(sam$Age_individual, 
                                    levels = c("4 months",
                                               "5 months",
                                               "6 months",
                                               "9 months",
                                               "11 months",
                                               "14 months"
                                    ))

sam$keys <- rownames(sam)

# adding microbial maturity, background data ====
names(microbial_age_df)[names(microbial_age_df) == "kindcode"] <- "Sample_aliquote"
names(microbial_age_df)[names(microbial_age_df) == "gender"] <- "sex"
names(microbial_age_df)[names(microbial_age_df) == "id"] <- "Family_ID"
col_keep <- c("Sample_aliquote","Family_ID","sex","siblings","birth_place","delivery_type","ageintrosolids",
              "birth_weigth","duration_breastfeeding","ChronoAge","MicrobialAge","spline","med","sd","RelativeMaturity",
              "maz") #the rest of the columns have matching data
microbial_age_df <- microbial_age_df %>% select(all_of(col_keep))
sam <-full_join(sam,microbial_age_df, by = "Sample_aliquote")
rm(col_keep,microbial_age_df)

# adding feeding data ====
names(spss_file)[names(spss_file) == "kindcode"] <- "Family_ID"
spss_file <- spss_file %>% select("breastfeeding_4m","breastfeeding_5m","breastfeeding_6m","breastfeeding_9m","breastfeeding_11m","breastfeeding_14m",
                                  "formulafeeding_4m", "formulafeeding_5m", "formulafeeding_6m", "formulafeeding_9m", "formulafeeding_11m","formulafeeding_14m",
                                  "Family_ID") %>% as.matrix() %>% as.data.frame()
family_ps <- unique(sort(sam$Family_ID))
family_spss <- unique(sort(spss_file$Family_ID))
missing_families <- setdiff(family_spss,family_ps)

spss_file <- dplyr::filter(spss_file, !(Family_ID %in% missing_families))

sam <-full_join(sam,spss_file, by = "Family_ID")

cols <- c("breastfeeding_4m","breastfeeding_5m","breastfeeding_6m","breastfeeding_9m","breastfeeding_11m","breastfeeding_14m",
          "formulafeeding_4m", "formulafeeding_5m", "formulafeeding_6m", "formulafeeding_9m", "formulafeeding_11m","formulafeeding_14m")
recode_vals <- c("0" = "no", "1" = "yes")

sam[cols] <- lapply(sam[cols], function(x) {
  x <- as.character(x)
  x[x %in% names(recode_vals)] <- recode_vals[x[x %in% names(recode_vals)]]
  x
})

rm(family_ps,family_spss,missing_families,spss_file)

# sanity check ====
x <- sort(sam$ageintrosolids.x)
y <- sort(sam$ageintrosolids.y)
setdiff(x,y) # other common columns are the same, except for age solid introduced
check <- sam %>% select(Family_ID,ageintrosolids.x,ageintrosolids.y)
colnames(sam)

#checking with original data
baku_f_I <- readRDS("H:/Penders_lab/Theoretical_work/Paper_1/Back-up_files/Generated_files/S4_107b_phyloseq_object_bacteria_filtered_infant.rds")
baku_f_I <- baku_f_I %>% subset_samples(used_in_MS_analysis == "yes")
real_meta <- as.data.frame(as.matrix(baku_f_I@sam_data))
m <- sort(real_meta$ageintrosolids)
setdiff(x,y)
setdiff(m,y)
setdiff(x,m)
#vector x from ps MS original object is incorrect, no idea what is it, keeping vector y from microbial age
names(sam)[names(sam) == "ageintrosolids.y"] <- "ageintrosolids"
sam$ageintrosolids.x <- NULL
rm(m,x,y,baku_f_I)

# adding word variables to data ====
#sex
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","sex")))
names(toadd)[names(toadd) == "sex"] <- "sex_word"
sam <- full_join(sam,toadd, by = "Sample_aliquote")
plyr::count(sam$bqq3_gender)
plyr::count(sam$sex)
plyr::count(sam$sex_word)
sam$bqq3_gender <- NULL
test <- sam %>% select(all_of(c("Sample_aliquote","sex","sex_word")))
sam$sex[grepl("1",sam$sex)]<-"0"
sam$sex[grepl("2",sam$sex)]<-"1"
rm(toadd,test)

#delivery mode
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","delivery_type_MM","delivery_type")))
names(toadd)[names(toadd) == "delivery_type_MM"] <- "delivery_type_word"
names(toadd)[names(toadd) == "delivery_type"] <- "delivery_type_bi_word"
sam <- full_join(sam,toadd, by = "Sample_aliquote")
plyr::count(sam$delivery_type)
plyr::count(sam$delivery_type_word)
plyr::count(sam$bqq14_delivery_type)
sam$bqq14_delivery_type <- NULL
sam$delivery_type_bi <- sam$delivery_type_word
sam$delivery_type_bi[grepl("vaginal",sam$delivery_type_bi)]<-"0"
sam$delivery_type_bi[grepl("caesarian section",sam$delivery_type_bi)]<-"1"
plyr::count(sam$delivery_type_bi)

rm(toadd)

#delivery place
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","delivery_place")))
names(toadd)[names(toadd) == "delivery_place"] <- "delivery_place_word"
sam <- full_join(sam,toadd, by = "Sample_aliquote")
plyr::count(sam$delivery_place_word)
plyr::count(sam$bqq13_delivery_place)
plyr::count(sam$birth_place)
sam$birth_place <- NULL
names(sam)[names(sam) == "bqq13_delivery_place"] <- "delivery_place"
plyr::count(sam$delivery_place)
sam$delivery_place <- as.character(sam$delivery_place)
sam$delivery_place[grepl("2",sam$delivery_place)]<-"0"
sam$delivery_place[grepl("1",sam$delivery_place)]<-"1"
rm(toadd)

#diet
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","diet_class","diet_class_1","diet_class_2","diet_class_3")))
sam <- full_join(sam,toadd, by = "Sample_aliquote")
plyr::count(sam$diet_class.x)
plyr::count(sam$diet_class.y)
test <- sam %>% select(all_of(c("Sample_aliquote","diet_class.x","diet_class.y")))
names(sam)[names(sam) == "diet_class.x"] <- "diet_class"
sam$diet_class.y <- NULL

plyr::count(sam$diet_class_3)
plyr::count(sam$CLASS_3)
sam$diet_class_3 <- NULL
names(sam)[names(sam) == "CLASS_3"] <- "diet_class_3"

plyr::count(sam$diet_class_2)
plyr::count(sam$CLASS_2)
sam$diet_class_2 <- NULL
names(sam)[names(sam) == "CLASS_2"] <- "diet_class_2"

sam$diet_class_1[grepl("no",sam$diet_class_1)]<-"0"
sam$diet_class_1[grepl("yes",sam$diet_class_1)]<-"1"
plyr::count(sam$diet_class_1)
plyr::count(sam$diet_class)

rm(toadd,test)
sort(colnames(sam))

#furrypet_indoors_6m
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","furrypet_indoors_6m")))
names(toadd)[names(toadd) == "furrypet_indoors_6m"] <- "furrypet_indoors_6m_word_ED"
sam <- full_join(sam,toadd, by = "Sample_aliquote")
names(sam)[names(sam) == "furrypet_indoors_6m"] <- "furrypet_indoors_MS_6m"
sam$furrypet_indoors_MS_6m_word <- as.character(sam$furrypet_indoors_MS_6m)
sam$furrypet_indoors_MS_6m_word[grepl("0",sam$furrypet_indoors_MS_6m_word)]<-"no"
sam$furrypet_indoors_MS_6m_word[grepl("1",sam$furrypet_indoors_MS_6m_word)]<-"home daycare"
sam$furrypet_indoors_MS_6m_word[grepl("2",sam$furrypet_indoors_MS_6m_word)]<-"home"
sam$furrypet_indoors_MS_6m_word[grepl("3",sam$furrypet_indoors_MS_6m_word)]<-"daycare"
plyr::count(sam$furrypet_indoors_6m_word_ED)
plyr::count(sam$furrypet_indoors_6m_MS)
plyr::count(sam$furrypet_indoors_MS_6m_word)
test <- sam %>% select(all_of(c("furrypet_indoors_6m_MS","furrypet_indoors_6m_word_ED","furrypet_indoors_MS_6m_word")))
rm(toadd,test)

#siblings
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","older_siblings")))
names(toadd)[names(toadd) == "older_siblings"] <- "older_siblings_ED"
sam <- full_join(sam,toadd, by = "Sample_aliquote")
names(sam)[names(sam) == "siblings"] <- "number_of_older_siblings_MS"
plyr::count(sam$number_of_older_siblings_MS)
plyr::count(sam$older_siblings_ED)
test <- sam %>% select(all_of(c("older_siblings_ED","number_of_older_siblings_MS")))
rm(toadd,test)

#antibiotics
toadd <- real_meta %>% select(all_of(c("Sample_aliquote","inf_ab_4m","inf_ab_5m","inf_ab_6m","inf_ab_9m","inf_ab_11m","inf_ab_14m")))
names(toadd)[names(toadd) == "inf_ab_4m"] <- "inf_ab_ED_4m"
names(toadd)[names(toadd) == "inf_ab_5m"] <- "inf_ab_ED_5m"
names(toadd)[names(toadd) == "inf_ab_6m"] <- "inf_ab_ED_6m"
names(toadd)[names(toadd) == "inf_ab_9m"] <- "inf_ab_ED_9m"
names(toadd)[names(toadd) == "inf_ab_11m"] <- "inf_ab_ED_11m"
names(toadd)[names(toadd) == "inf_ab_14m"] <- "inf_ab_ED_14m"
sam <- full_join(sam,toadd, by = "Sample_aliquote")
rm(toadd)

names(sam)[names(sam) == "ch4mq5_ab_freq"] <- "ch4mq5_ab_freq_MS_4m"
names(sam)[names(sam) == "ch5mq5_ab_freq"] <- "ch5mq5_ab_freq_MS_5m"
names(sam)[names(sam) == "ch6mq16_ab_freq"] <- "ch6mq16_ab_freq_MS_6m"

antibio <- c("ch4mq5_ab_freq_MS_4m", "ch5mq5_ab_freq_MS_5m", "ch6mq16_ab_freq_MS_6m")
recode_vals <- c("1" = "never", "2" = "1 time", "3" = "more than 1 time")

sam[antibio] <- lapply(sam[antibio], function(x) {
  x <- as.character(x)
  x[x %in% names(recode_vals)] <- recode_vals[x[x %in% names(recode_vals)]]
  x
})

# picking up keys ====
keys <- sam %>% select(keys,Sample_aliquote)
codes <- real_meta %>% select(MMHP_SampleID,Sample_aliquote)
codes <- full_join(codes, keys, by = "Sample_aliquote")
codes$keys <- NULL
sam <- full_join(codes, sam, by = "Sample_aliquote")
rownames(sam) <- sam$Sample_aliquote
rm(codes)

# cleaning otu ====
#change the names of the rows in otu from whole string to only species names to make it easy to analyse
taxa_raw <- otu
otu <- tibble::rownames_to_column(otu, "taxa")
otu[c('k','p','c','o','f','g','x')]<- str_split_fixed(otu$taxa, "\\|", n = 7)
otu[c('z','s')]<- str_split_fixed(otu$x, "__", n = 2)
rownames(otu) <- otu$s
otu = subset(otu, select = -c(k,p,c,o,f,g,z,x,s,taxa))

otu <- as.data.frame(t(otu))
otu$keys <- rownames(otu)
otu <- full_join(keys, otu, by = "keys")
rownames(otu) <- otu$Sample_aliquote
otu <- otu %>%  select(-c(keys,Sample_aliquote))
otu <- as.data.frame(t(otu))
rm(keys)

# cleaning taxa ====
rownames(taxa) <- taxa$Species

# combining back to phyloseq ====
#change the names of the rows in taxa and otu from whole string to only species names to make it easy to analyse
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)
#transform data to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(sam)
rm(taxa_raw)
#and put them in one object
phylo_object_bacteria<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)
#check if ok
sample_sums(phylo_object_bacteria) 
sample_names(phylo_object_bacteria) 
rank_names(phylo_object_bacteria)        
rm(phylo_OTU,phylo_TAX,otu_mat,phylo_samples,tax_mat)
rm(taxa,sam,otu)
#adding observed richness to the data, chao1 had the same value sas observed_richness
phylo_object_bacteria <- ps_calc_richness(
  phylo_object_bacteria,
  "Species",
  index = "observed",
  detection = 0
)

# saving phyloseq ====
#saveRDS(phylo_object_bacteria, 'S1.1_cleaned_phyloseq_object_MS.rds')
