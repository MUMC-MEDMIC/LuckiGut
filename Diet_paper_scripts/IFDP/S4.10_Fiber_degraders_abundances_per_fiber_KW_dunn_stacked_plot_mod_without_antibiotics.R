# LOAD LIBRARIES ====
library(phyloseq)
library(microViz)
library(patchwork)
library(ggh4x)
library(FSA)
library(ggpubr)
library(tidyverse)
library(svglite)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
spearman = read.delim("Output/Files/S4.9_baku_vs_fiber_per_TP_spearman_adjust_without_antibiotics.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
baku_f_I <- readRDS("Output/Files/S1.6_cleaned_phyloseq_object_MS_without_antibiotics.rds")

# cleaning the data ====
#marking significant values with * for adj p values
spearman <- spearman %>% 
  mutate(sig_padj = if_else(
    condition = p_value_adj < 0.05, 
    true = "*", 
    false = ""
  ))

spearman <- spearman %>% mutate(bacteria = gsub(" ", "_", bacteria))

# filtering according to conditions ====
spearman_corr <- spearman %>% filter(Correlation_coefficient_r > 0.3) #filtering weak and negative correlations left n=411
spearman_corr_pval_05 <- spearman_corr %>% filter(p_value_adj < 0.05) #filtering weak and negative correlations left n=251
rm(spearman,spearman_corr)

#counting at how many TP fiber-baku pairs can be detected
spearman_corr_pval_05 <- spearman_corr_pval_05 %>% mutate(fib_baku = paste(bacteria, fiber, sep = "@"))
needed_pairs <- spearman_corr_pval_05 %>% group_by(fib_baku) %>% tally() %>% filter(n >= 2) #filtering out those that only encountered once
needed_pairs <- needed_pairs %>% separate(fib_baku, into = c("bacteria", "fiber"), sep = "@")
unik_baku <- sort(unique(needed_pairs$bacteria)) # 19 uniq bacteria
unik_fiber <- sort(unique(needed_pairs$fiber)) # 15 uniq fibers

needed_pairs <- spearman_corr_pval_05 %>% group_by(fib_baku) %>% tally() %>% filter(n >= 2) 
needed_pairs <- needed_pairs$fib_baku

spearman_corr_pval_05 <- spearman_corr_pval_05 %>% filter(fib_baku %in% needed_pairs)
baku_vector <- unique(spearman_corr_pval_05$bacteria)
rm(needed_pairs)

# cleaning phyloseq data ====
#subsetting per TP and transforming data 
rank_viz <- "Species"
mon4 <- baku_f_I %>% subset_samples(Age_individual == "4 months") %>% tax_transform("compositional", rank = rank_viz)

#taking otu table out
reads_4m <- as.data.frame(mon4@otu_table)
rm(mon4)

#making baku into column
reads_4m$bacteria <- rownames(reads_4m) 

#turning into long format
reads_4m <- reads_4m %>% pivot_longer(cols = -c(bacteria), names_to = "Sample_aliquote", values_to = "read")

#removing not needed bacteria
filtered_reads <- reads_4m %>% filter(bacteria %in% baku_vector)

#creating keys
keys_sample_name <- as.data.frame(as.matrix(baku_f_I@sam_data))
keys_sample_name <- keys_sample_name  %>%  select (MMHP_SampleID,Family_ID, Age_individual, diet_class,Sample_aliquote)
rm(reads_4m,baku_vector)

# connecting fiber to baku ====
count_of_fb <- spearman_corr_pval_05 %>% select (fiber,bacteria) %>% group_by(fiber,bacteria) %>% count()

#if this code produces only 1 element in the list, restart R. should be equal to the fibers n
#collecting all bacteria for each fiber
bacteria_per_fiber <- count_of_fb %>%
  group_by(fiber) %>%
  summarise(bacteria_list = list(bacteria)) %>%
  tibble::deframe()

rm(count_of_fb,spearman_corr_pval_05)

#summirising as a table
max_length <- max(sapply(bacteria_per_fiber, length))
pad_with_na <- function(x, max_length) {
  length(x) <- max_length
  return(x)
}
list <- lapply(bacteria_per_fiber, pad_with_na, max_length)
list <- as.data.frame(list)

#reshaping data into different formats
list_long <- list %>% pivot_longer(
  cols = everything(),   # Pivot all columns
  names_to = "fibers",     # Name of the new column that will store column names
  values_to = "bacteria"    # Name of the new column that will store values
)

list_long <- list_long %>% filter(!is.na(bacteria))

list_wider <- list_long %>%
  mutate(present = "x") %>%
  pivot_wider(
    names_from = fibers,
    values_from = present,
    values_fn = ~ "x",
    values_fill = NA
  )

list_wider <- list_wider %>% mutate(bacteria = gsub("_", " ", bacteria))

#write.table (list_wider, "S7_29_bacteria_in_selected_fibers.txt", row.names=FALSE,sep = "\t")
rm(list,max_length,pad_with_na,list_wider,list_long)

# counting for 15 fibers read count for all degrades ====
# Function to calculate sum of read values for each bacteria in the list
calculate_sum_reads <- function(fiber_name, bacteria_list, reads_data) {
  # Filter reads_data to include only bacteria from the filtered list
  filtered_reads <- reads_data %>%
    filter(bacteria %in% bacteria_list)
  
  # Summarize the sum of read values for each MMHP_SampleID
  summarized_reads <- filtered_reads %>%
    group_by(Sample_aliquote) %>%
    summarise(total_read = sum(read))
  
  # Return summarized reads and fiber name
  return(list(fiber_name = fiber_name, summarized_reads = summarized_reads))
}

# Iterate over each element in bacteria_per_fiber list and calculate sum of reads
result_list <- lapply(names(bacteria_per_fiber), function(fiber_name) {
  calculate_sum_reads(fiber_name, bacteria_per_fiber[[fiber_name]], filtered_reads)
})
fiber_names <- names(bacteria_per_fiber)
names(result_list) <- fiber_names
rm(fiber_names)

# adding TP to the sum of the bakus ====
# Function to join each element of result_list with df using MMHP_code
# Convert each element of result_list into a dataframe and combine them
result_df_list <- lapply(result_list, function(element) {
  fiber_name <- element$fiber_name
  summarized_reads <- element$summarized_reads
  df <- as.data.frame(summarized_reads)
  df$fiber_name <- fiber_name
  return(df)
})

# Combine all dataframes into a single dataframe
combined_df <- do.call(rbind, result_df_list)
check <- sort(unique(combined_df$Sample_aliquote)) # 84 samples had these fiber-baku pairs
check2 <- sort(unique(keys_sample_name$Sample_aliquote)) #we have 389 samples altogether

combined_df <- full_join(combined_df,keys_sample_name, by = "Sample_aliquote")
check3 <- sort(unique(combined_df$Sample_aliquote)) # should be 389
rm(check,check2,check3)

rm(keys_sample_name,result_df_list,result_list,calculate_sum_reads)

# cleaning data
combined_df$Age_individual <- gsub(" months", "", combined_df$Age_individual)
combined_df$Age_individual <- as.character(combined_df$Age_individual)
combined_df$Age_individual <- as.numeric(combined_df$Age_individual)

# data to wide format ====
wide_format <- combined_df %>% pivot_wider(names_from = fiber_name, values_from = total_read)
fiber_names15 <- unik_fiber

wide_format$diet_class <- paste0("diet_", wide_format$diet_class)
wide_format$diet_class <- as.factor(wide_format$diet_class)

mon4 <- wide_format %>%  subset(Age_individual == "4")
rm(wide_format)

# Kruskal wallis test ====
# specify levels of clusters
# Loop over each bacterium name in the vector

# 4 months
# Initialize an empty data frame to store results
results_KW <- data.frame(fiber = character(),
                         h_statistics = numeric(),#H statistic is asymptotically chi-square distributed, that is why it is chi-squared in R
                         df = numeric(),
                         p_value = numeric(),
                         method = character(),
                         stringsAsFactors = FALSE)

# Loop over each bacterium name in the vector
for (fiber_name in fiber_names15) {
  # Create the formula for the Kruskal-Wallis test dynamically
  formula <- as.formula(paste(fiber_name, "~ diet_class"))
  
  # Perform the Kruskal-Wallis test
  result <- kruskal.test(formula, data = mon4)
  
  # Add the results to the data frame
  results_KW <- rbind(results_KW, data.frame(fiber = fiber_name,
                                             h_statistics = result$statistic,
                                             df = result$parameter,
                                             p_value = result$p.value,
                                             metod = result$method,
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
fiber_after_KW_4mon <- results_KW_sig$fiber 
results_KW$age <- "4"
results_KW_4mon <- results_KW

rm(result,formula,results_KW,results_KW_sig)

#write.table (results_KW_4mon, "S4.10_fiber_vs_diet_results_KW_4months_without_antibiotics.txt", row.names=FALSE,sep = "\t")

# Dunn test ====
# 4 months
# Initialize an empty list to store results
dunn_results <- list()

# Loop over each bacterium name in the vector
for (fiber_name in fiber_after_KW_4mon) {
  # Create the formula for the Dunn's test dynamically
  formula <- as.formula(paste(fiber_name, "~ diet_class"))
  
  # Perform Dunn's test with adjustment method "bh"
  result <- dunnTest(formula, data = mon4, method = "bh")
  
  # Store the result in the list
  dunn_results[[fiber_name]] <- result
}

# Create an empty list to store the results
results_dunn <- list()

# Loop through each bacteria name in bacteria_names_after_KW
for (fiber_name in fiber_after_KW_4mon) {
  #pick results to a separate df
  result_df <- dunn_results[[fiber_name]]$res
  
  # Add the result data frame to the list with the corresponding fiber name
  results_dunn[[fiber_name]] <- result_df
}

# Merge the list of data frames into a single data frame
results_dunn <- bind_rows(results_dunn, .id = "fiber")

results_dunn <- results_dunn %>% mutate(Comparison = gsub("-", "@", Comparison))
results_dunn$age <- "4"
results_dunn_mon4 <- results_dunn

rm(results_dunn,dunn_results,mon4,fiber_after_KW_4mon)

results_dunn_mon4 <- results_dunn_mon4 %>% 
  mutate(sig_p_adj = if_else(
    condition = P.adj < 0.05, 
    true = "*", 
    false = ""
  ))

#write.table (results_dunn_mon4, "S4.10_fiber_vs_diet_results_dunn_4months_withou_antibiotics.txt", row.names=FALSE,sep = "\t")

rm(result_df,result,fiber_name)

# selecting needed baku species ====
rm(fiber_names15,formula,unik_baku,unik_fiber)

sig <- results_dunn_mon4 %>% filter(P.adj < 0.05)
fibers <- sort(unique(sig$fiber))
fibers

baku_arabinoxylan <- bacteria_per_fiber$Arabinoxylan
baku_rhamnogalacturonan <- bacteria_per_fiber$Rhamnogalacturonan
baku_xylan <- bacteria_per_fiber$Xylan
rm(fibers)

#subsetting per TP, changing counts
#selecting only needed meta data
#processing data, this will keep original reads for each sample
baku_f_I_filtered <- baku_f_I %>% ps_select(diet_class,Age_individual,MMHP_SampleID) 
rm(baku_f_I)

#processing data, this will keep original reads for each sample
mon4_arabinoxylan <- baku_f_I_filtered %>% subset_samples(Age_individual == "4 months") %>% #subsetting for needed age
  phyloseq::merge_samples(group = "MMHP_SampleID", fun = mean) %>% #merges samples based on seq name
  psmelt() %>% #melting ps into long df
  group_by(Sample) %>% #groupping based on a specified grouping variable, sample names variable is named Sample, check colnames
  mutate(Species_groupped = ifelse(Species %in% baku_arabinoxylan, Species, "Other"),  #renaming all other bakus, if they are not in arabinoxylan vector
         Abundance_groupped = ifelse(Species_groupped == "Other", sum(Abundance[Species_groupped == "Other"]), Abundance)) %>% #counting sum for other bacteria
  distinct(diet_class, Species_groupped, .keep_all = TRUE) #removing duplicates of others

mon4_arabinoxylan$Species_groupped <- gsub("_", " ", mon4_arabinoxylan$Species_groupped)

mon4_rhamnogal <- baku_f_I_filtered %>% subset_samples(Age_individual == "4 months") %>%
  phyloseq::merge_samples(group = "MMHP_SampleID", fun = mean) %>% 
  psmelt() %>% 
  group_by(Sample) %>% 
  mutate(Species_groupped = ifelse(Species %in% baku_rhamnogalacturonan, Species, "Other"),  
         Abundance_groupped = ifelse(Species_groupped == "Other", sum(Abundance[Species_groupped == "Other"]), Abundance)) %>% 
  distinct(diet_class, Species_groupped, .keep_all = TRUE) 

mon4_rhamnogal$Species_groupped <- gsub("_", " ", mon4_rhamnogal$Species_groupped)

mon4_xylan <- baku_f_I_filtered %>% subset_samples(Age_individual == "4 months") %>%
  phyloseq::merge_samples(group = "MMHP_SampleID", fun = mean) %>%
  psmelt() %>% 
  group_by(Sample) %>% 
  mutate(Species_groupped = ifelse(Species %in% baku_xylan, Species, "Other"),  
         Abundance_groupped = ifelse(Species_groupped == "Other", sum(Abundance[Species_groupped == "Other"]), Abundance)) %>% 
  distinct(diet_class, Species_groupped, .keep_all = TRUE) 

mon4_xylan$Species_groupped <- gsub("_", " ", mon4_xylan$Species_groupped)

# cleaning for plotting ====
mon4_arabinoxylan <- mon4_arabinoxylan %>% filter(Species_groupped != "Other")
mon4_rhamnogal <- mon4_rhamnogal %>% filter(Species_groupped != "Other")
mon4_xylan <- mon4_xylan %>% filter(Species_groupped != "Other")

# adding real relative abundances to the scirpt ====
relative_abundace <- baku_f_I_filtered %>% tax_transform("compositional", rank = rank_viz) %>% 
  phyloseq::merge_samples(group = "MMHP_SampleID", fun = mean) %>% 
  psmelt()
relative_abundace <- relative_abundace %>% rename(Relative_abundance = Abundance)
relative_abundace <- relative_abundace %>% mutate(merging_column = paste(OTU, Sample, sep = "@"))
relative_abundace_merge <- relative_abundace %>% select(merging_column, Relative_abundance)

#creating common merged column
mon4_arabinoxylan <- mon4_arabinoxylan %>% mutate(merging_column = paste(OTU, Sample, sep = "@"))
mon4_rhamnogal <- mon4_rhamnogal %>% mutate(merging_column = paste(OTU, Sample, sep = "@"))
mon4_xylan <- mon4_xylan %>% mutate(merging_column = paste(OTU, Sample, sep = "@"))

#merging relative and count abundaces 
mon4_arabinoxylan <- left_join(mon4_arabinoxylan,relative_abundace_merge, by = "merging_column")
mon4_rhamnogal <- left_join(mon4_rhamnogal,relative_abundace_merge, by = "merging_column")
mon4_xylan <- left_join(mon4_xylan,relative_abundace_merge, by = "merging_column")

rm(relative_abundace,relative_abundace_merge,baku_f_I_filtered,bacteria_per_fiber,combined_df,filtered_reads)

# filtering dunn test ====
sig$unpack <- sig$Comparison
sig <- sig %>% separate(unpack, into = c("group1", "group2"), sep = "@")
sig <- sig %>% mutate(group1 = str_replace(group1, "diet_", ""))
sig <- sig %>% mutate(group2 = str_replace(group2, "diet_", ""))
sig$group1 <- gsub(" ", "", sig$group1)
sig$group2 <- gsub(" ", "", sig$group2)
sig$group1 <- as.numeric(sig$group1)
sig$group2 <- as.numeric(sig$group2)
sig$y.position <- 4

sig_xyl <- sig %>% filter(fiber == "Xylan")
sig_rha <- sig %>% filter(fiber == "Rhamnogalacturonan")
sig_arab <- sig %>% filter(fiber == "Arabinoxylan")

rm(results_KW_4mon,sig,results_dunn_mon4,rank_viz,baku_xylan,baku_arabinoxylan,baku_rhamnogalacturonan)

# plotting ====
species_colors_5 <- c("#DBD3D1","#F5B8B6","#A56266","#34686C","#84869B")
species_colors_6 <- c("#DBD3D1","#F5B8B6","#A56266","#34686C","#84869B","#003856")
species_colors_8 <- c("#8B9D77","#A1B4A5","#F5B8B6", "#34686C","#B3B2A5","#C5D3E3","#D7CFA1","#D3B8C5")
species_colors_14 <- c("#8B9D77","#A1B4A5","#DBD3D1","#F5B8B6","#A56266","#22272A", "#34686C", "#B3B2A5","#C5D3E3","#003856","#84869B","#849E8A","#D7CFA1","#D3B8C5")

plot_mon4_rhamnogal <- ggplot(mon4_rhamnogal, aes(x = diet_class, y = Relative_abundance)) +
  geom_bar(aes(fill = Species_groupped),stat = "identity", position = "stack") +
  scale_fill_manual(values = species_colors_6) +
  labs(x = "Diet Class", y = "Relative abundance")  + ggtitle("Rhamnogalacturonan at 4 months")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, margin = margin(b = 1)), #grid size
        plot.title = element_text(size = 12))
plot_mon4_rhamnogal <- plot_mon4_rhamnogal + stat_pvalue_manual(sig_rha, label = "sig_p_adj",y.position = 3.8, step.increase = 0.05,tip.length = 0)

plot_mon4_arabinoxylan_with <- ggplot(mon4_arabinoxylan, aes(x = diet_class, y = Relative_abundance)) +
  geom_bar(aes(fill = Species_groupped),stat = "identity", position = "stack") +
  scale_fill_manual(values = species_colors_14) +
  labs(x = "Diet Class", y = "Relative abundance")  + ggtitle("Arabinoxylan at 4 months")+
  theme_minimal()+
  labs(fill = "Degraders")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "cm"),        # Shrink key boxes
        legend.key.height = unit(0.5, "cm"),     # Make thinner
        legend.key.width = unit(0.5, "cm"),      # Make narrower
        legend.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(size = 12))+
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 1)))
plot_mon4_arabinoxylan_with <-  plot_mon4_arabinoxylan_with + stat_pvalue_manual(sig_arab, label = "sig_p_adj",y.position = 4.8, step.increase = 0.001,tip.length = 0)
  
plot_mon4_xylan <- ggplot(mon4_xylan, aes(x = diet_class, y = Relative_abundance)) +
  geom_bar(aes(fill = Species_groupped),stat = "identity", position = "stack") +
  scale_fill_manual(values = species_colors_8) +
  labs(x = "Diet Class", y = NULL)  + ggtitle("Xylan at 4 months")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 1)), #grid size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12))
plot_mon4_xylan <- plot_mon4_xylan + stat_pvalue_manual(sig_xyl, label = "sig_p_adj",y.position = 3.7, step.increase = 0.05,tip.length = 0)

barplot_with <- (plot_mon4_rhamnogal + plot_mon4_xylan) / (plot_mon4_arabinoxylan_with)+
  plot_layout() & 
  theme(plot.margin = margin(1, 0, 2, 0)) #top, right, bottom, left
barplot_with

# saving figure for the paper ====
ggsave(
  filename   = "S4.10_ifdp_degraders_abundance_per_fiber_without_antibiotics.tiff",
  plot       = barplot_with,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

