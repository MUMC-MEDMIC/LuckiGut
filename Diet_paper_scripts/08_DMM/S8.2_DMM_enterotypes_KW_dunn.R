# libraries ====
library(microViz)
library(tidyverse)
library(phyloseq)
library(FSA)
library(colorRamp2)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# import data ====
baku = read.delim("Input/11_DMM_top30_driving_bacterial_species_MS.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
dmm = read.delim("Input/10_DMM_cluster_per_sample_from_MS.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
baku_f_I <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# cleaning phyloseq data ====
keys_sample_name <- as.data.frame(as.matrix(baku_f_I@sam_data))
keys_sample_name <- keys_sample_name  %>%  select (MMHP_SampleID,Family_ID, Age_individual, diet_class,Sample_aliquote)

rank_viz <- "Species"
baku_f_I_rel <- baku_f_I %>% tax_transform("compositional", rank = rank_viz)
reads_rel <- as.data.frame(as.matrix(baku_f_I_rel@otu_table))
rm(baku_f_I_rel,rank_viz)

# add DMM to data ====
dmm$Sample_aliquote <- rownames(dmm)
dmm <- dmm %>% mutate(DMM_Cluster = recode(DMM_Cluster, "5" = 5, "3" = 2, "6" = 3, "4" = 6, "1" = 4, "2" = 1))
combined_df <- full_join(dmm,keys_sample_name, by = "Sample_aliquote")
rm(keys_sample_name,dmm)
colnames(combined_df)
#select only 30 bakus from DMM
baku$bacteria <- rownames(baku)
bacteria_names <- unique(baku$bacteria)
bacteria_names30 <- bacteria_names[1:30]

#filter otu table for 30 baku
reads_rel$bacteria <- rownames(reads_rel)
reads_rel <- reads_rel %>% filter(bacteria %in% bacteria_names30)
reads_rel$bacteria <- NULL

reads_rel <- reads_rel %>% t() %>% as.matrix() %>% as.data.frame()
reads_rel$Sample_aliquote <- rownames(reads_rel)
combined_df_otu_rel <- full_join(combined_df,reads_rel, by = "Sample_aliquote")

# relative ====
combined_df_otu_rel <- full_join(combined_df,reads_rel, by = "Sample_aliquote")
abundant_species_enterotypes_rel <- combined_df_otu_rel %>%  select(all_of(c(bacteria_names30,"DMM_Cluster")))
abundant_species_enterotypes_rel_long <- abundant_species_enterotypes_rel %>% gather(key = "key", value = "value", -DMM_Cluster)
abundant_species_enterotypes_rel_wide <- abundant_species_enterotypes_rel_long %>%
  group_by(key, DMM_Cluster) %>% summarise(value = mean(value)) %>%
  pivot_wider(names_from = DMM_Cluster, values_from = value)

#write.table (abundant_species_enterotypes_rel_wide, "S8.2_DMM_top30_driving_bacterial_species_rel_abundance.txt", row.names=FALSE,sep = "\t")

rm(combined_df,reads_rel,baku,abundant_species_enterotypes_rel,abundant_species_enterotypes_rel_long,abundant_species_enterotypes_rel_wide,bacteria_names)

# Kurska wallis test ====
#spesify levels of clusters
combined_df_otu_rel$DMM_Cluster <- as.character(combined_df_otu_rel$DMM_Cluster)
combined_df_otu_rel$DMM_Cluster <- ordered(combined_df_otu_rel$DMM_Cluster,
                         levels = c("1","2","3","4","5","6"))

# Loop over each bacterium name in the vector
# Initialize an empty data frame to store results
results_KW <- data.frame(bacteria = character(),
                         h_statistics = numeric(), #H statistic is asymptotically chi-square distributed, that is why it is chi-squared in R
                         df = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each bacterium name in the vector
for (bacterium_name in bacteria_names30) {
  # Create the formula for the Kruskal-Wallis test dynamically
  formula <- as.formula(paste(bacterium_name, "~ DMM_Cluster"))
  
  # Perform the Kruskal-Wallis test
  result <- kruskal.test(formula, data = combined_df_otu_rel)
  
  # Add the results to the data frame
  results_KW <- rbind(results_KW, data.frame(bacteria = bacterium_name,
                                             h_statistics = result$statistic,
                                             df = result$parameter,
                                             p_value = result$p.value,
                                             stringsAsFactors = FALSE))
}

#marking significant values with * for adj p values
results_KW <- results_KW %>% 
  mutate(sig_p_value = if_else(
    condition = p_value < 0.05, 
    true = "*", 
    false = ""
  ))

results_KW_sig <- results_KW %>% filter(sig_p_value == "*")
bacteria_names_after_KW <- results_KW_sig$bacteria 

#saving table, groupping food items per groups
#write.table (results_KW, "S17_23_DMM_enterotype_results_KW.txt", row.names=FALSE,sep = "\t")

rm(result,bacterium_name,formula,results_KW,results_KW_sig,bacteria_names30)

# Dunn test ====
combined_df_otu_rel$DMM_Cluster <- as.integer(combined_df_otu_rel$DMM_Cluster)
combined_df_otu_rel$DMM_Cluster <- as.factor(combined_df_otu_rel$DMM_Cluster)

# Initialize an empty list to store results
dunn_results <- list()

# Loop over each bacterium name in the vector
for (bacterium_name in bacteria_names_after_KW) {
  # Create the formula for the Dunn's test dynamically
  formula <- as.formula(paste(bacterium_name, "~ DMM_Cluster"))
  
  # Perform Dunn's test with adjustment method "bh"
  result <- dunnTest(formula, data = combined_df_otu_rel, method = "bh")
  
  # Store the result in the list
  dunn_results[[bacterium_name]] <- result
}

# Create an empty list to store the results
results_dunn <- list()

# Loop through each bacteria name in bacteria_names_after_KW
for (bacteria_name in bacteria_names_after_KW) {
  # Replace "Bifidobacterium_longum" with the current bacteria name
  result_df <- dunn_results[[bacteria_name]]$res
  
  # Add the result data frame to the list with the corresponding bacteria name
  results_dunn[[bacteria_name]] <- result_df
}

# Merge the list of data frames into a single data frame
results_dunn <- bind_rows(results_dunn, .id = "bacteria")

results_dunn <- results_dunn %>% 
  mutate(sig_p_adj = if_else(
    condition = P.adj < 0.05, 
    true = "*", 
    false = ""
  ))

results_dunn <- results_dunn %>% mutate(Comparison = gsub("-", "_", Comparison))

#write.table (results_dunn, "S8.2_DMM_enterotype_dunn.txt", row.names=FALSE,sep = "\t")



