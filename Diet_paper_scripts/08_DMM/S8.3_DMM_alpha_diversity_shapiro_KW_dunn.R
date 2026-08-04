# LOAD LIBRARIES ====
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(FSA)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ps <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
#in original MS calculations the DMM
dmm = read.delim("Input/10_DMM_cluster_per_sample_from_MS.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning and renaming clusters ====
dmm$Sample_aliquote <- rownames(dmm)
dmm <- dmm %>% mutate(DMM_Cluster = recode(DMM_Cluster, "5" = 5, "3" = 2, "6" = 3, "4" = 6, "1" = 4, "2" = 1))
dmm$Family_ID <- as.integer(str_sub(dmm$Sample_aliquote, start = 1, end = 4))

# merging data ====
ps <- ps %>% ps_join(dmm, by = "Sample_aliquote")
df_dmm_diversity <- rbind(data.frame(shannon = vegan::diversity(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1))@otu_table)), dmm = 1),
                          data.frame(shannon = vegan::diversity(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2))@otu_table)), dmm = 2),
                          data.frame(shannon = vegan::diversity(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3))@otu_table)), dmm = 3),
                          data.frame(shannon = vegan::diversity(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4))@otu_table), index = "shannon"),
                                     specnr = vegan::specnumber(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4))@otu_table)), dmm = 4),
                          data.frame(shannon = vegan::diversity(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5))@otu_table), index = "shannon"),
                                     specnr = vegan::specnumber(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5))@otu_table)), dmm = 5),
                          data.frame(shannon = vegan::diversity(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6))@otu_table), index = "shannon"),
                                     specnr = vegan::specnumber(t((ps %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6))@otu_table)), dmm = 6))
df_dmm_diversity$Sample_aliquote <- rownames(df_dmm_diversity)
df_dmm_diversity$Family_ID <- as.character(str_sub(df_dmm_diversity$Sample_aliquote, start = 1, end = 4))
df_dmm_diversity$dmm <- as.factor(df_dmm_diversity$dmm)

full<- as.data.frame(as.matrix(ps@sam_data))
sel_full <- full %>%  select(diet_class, Sample_aliquote)
df_dmm_diversity <- sel_full %>% full_join(df_dmm_diversity, by = "Sample_aliquote")

rm(ps,dmm,full)

# calculating statistics ====
# checking normality for shannon ====
# Create a list to store the results
shannon_list <- list()
# Loop through the treatment groups
for (i in 1:6) {
  # Subset the data for the current treatment group
  cluster_data <- subset(df_dmm_diversity, dmm == as.character(i))$shannon
  
  # Assign the subset data to the list with a dynamic name
  shannon_list[[paste0("cluster", i)]] <- cluster_data
}

# Create an empty data frame to store the Shapiro-Wilk test results
shapiro_shannon_results <- data.frame(
  cluster = character(),
  W_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the list and perform the Shapiro-Wilk test
for (i in 1:length(shannon_list)) {
  cluster_name <- names(shannon_list)[i]
  cluster_data <- shannon_list[[i]]
  
  # Perform the Shapiro-Wilk test
  shapiro_test <- shapiro.test(cluster_data)
  
  # Store the results in the data frame
  shapiro_shannon_results <- rbind(shapiro_shannon_results, data.frame(
    cluster = cluster_name,
    W_statistic = shapiro_test$statistic,
    p_value = shapiro_test$p.value
  ))
}

rm(shapiro_test,cluster_data,cluster_name)

# checking normality for species number ====
# Create a list to store the results
richness_list <- list()
# Loop through the treatment groups
for (i in 1:6) {
  # Subset the data for the current treatment group
  cluster_data <- subset(df_dmm_diversity, dmm == as.character(i))$specnr
  
  # Assign the subset data to the list with a dynamic name
  richness_list[[paste0("cluster", i)]] <- cluster_data
}

# Create an empty data frame to store the Shapiro-Wilk test results
shapiro_richness_results <- data.frame(
  cluster = character(),
  W_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the list and perform the Shapiro-Wilk test
for (i in 1:length(richness_list)) {
  cluster_name <- names(richness_list)[i]
  cluster_data <- richness_list[[i]]
  
  # Perform the Shapiro-Wilk test
  shapiro_test <- shapiro.test(cluster_data)
  
  # Store the results in the data frame
  shapiro_richness_results <- rbind(shapiro_richness_results, data.frame(
    cluster = cluster_name,
    W_statistic = shapiro_test$statistic,
    p_value = shapiro_test$p.value
  ))
}

rm(shapiro_test,cluster_data,cluster_name)

# merging results ====
shapiro_shannon_results$alpha_diversity <- 'shannon'
shapiro_richness_results$alpha_diversity <- 'richness'
shapiro_all <- rbind(shapiro_shannon_results,shapiro_richness_results)

shapiro_all <- shapiro_all %>% 
  mutate(sig_padj = if_else(
    condition = p_value >= 0.05, 
    true = "normally distributed", 
    false = "not normally distributed"
  ))

#write.table (shapiro_all, "S18_35_DMM_alpha_diversity_shapiro.txt", row.names=FALSE,sep = "\t")

# plotting ====
ggplot(df_dmm_diversity, aes(x=shannon)) + geom_histogram()+facet_wrap(~dmm)
ggplot(df_dmm_diversity, aes(x=specnr)) + geom_histogram()+facet_wrap(~dmm)

# KW test ====
kw_test_shannon <- kruskal.test(shannon ~ as.factor(dmm), data = df_dmm_diversity)
kw_test_specnr <- kruskal.test(specnr ~ as.factor(dmm), data = df_dmm_diversity)

# Create a data frame with the results
results <- data.frame(
  Test = c("Shannon", "Specnr"),
  h_statistic = c(kw_test_shannon$statistic, kw_test_specnr$statistic),
  P_Value = c(kw_test_shannon$p.value, kw_test_specnr$p.value)
)

#marking significant values with * for adj p values
results <- results %>% 
  mutate(sig_p_value = if_else(
    condition = P_Value < 0.05, 
    true = "*", 
    false = ""
  ))

#write.table (results, "S8.3_alpha_diversity_DMM_KW.txt", row.names=FALSE,sep = "\t")

# dunn posthoc ====
# Perform Dunn's test for 'shannon' across different 'dmm' groups
dunn_test_shannon <- dunnTest(shannon ~ as.factor(dmm), data = df_dmm_diversity, method = "bh")

# Convert results to data frame
dunn_shannon_df <- data.frame(
  Comparison = dunn_test_shannon$res$Comparison,
  Z = dunn_test_shannon$res$Z,
  p_value = dunn_test_shannon$res$P.unadj,
  p_value_adj = dunn_test_shannon$res$P.adj
)

# Perform Dunn's test for 'specnr' across different 'dmm' groups
dunn_test_specnr <- dunnTest(specnr ~ as.factor(dmm), data = df_dmm_diversity, method = "bh")

# Convert results to data frame
dunn_specnr_df <- data.frame(
  Comparison = dunn_test_specnr$res$Comparison,
  Z = dunn_test_specnr$res$Z,
  p_value = dunn_test_specnr$res$P.unadj,
  p_value_adj = dunn_test_specnr$res$P.adj
)

# Combine results into a single data frame
dunn_results <- rbind(
  cbind(Test = "Shannon", dunn_shannon_df),
  cbind(Test = "Specnr", dunn_specnr_df)
)

#marking significant values with * for adj p values
dunn_results <- dunn_results %>% 
  mutate(sig_p_adj = if_else(
    condition = p_value_adj < 0.05, 
    true = "*", 
    false = "ns"
  ))

dunn_results <- dunn_results %>% separate(Comparison, into = c("group1", "group2"), sep = " - ")

#write.table (dunn_results, "S8.3_alpha_diversity_DMM_dunn.txt", row.names=FALSE,sep = "\t")

# cleaning for the figure ====
dunn_results_filtered <- dunn_results %>%
  filter((group1 == "1" & group2 == "2") |
      (group1 == "2" & group2 == "3")  |
      (group1 == "3" & group2 == "4")  |
  (group1 == "4" & group2 == "5")      |
  (group1 == "5" & group2 == "6"))

# Supplementary figure 8 ====
shannon <- dunn_results_filtered %>% filter(Test == "Shannon") %>% filter(sig_p_adj != "ns")
shannon_plot <- ggplot(data = df_dmm_diversity, aes(x=dmm, y=shannon, col = dmm)) + 
  geom_boxplot() +
  labs(x = "DMM cluster", y = "Shannon diversity") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  #scale_x_discrete(labels = c(1:6)) +
  scale_color_manual(values = c("#796248","#a47852","#cd9c58", "#e6cd98","#f8edcf", "#d4ede8")) +
  theme(legend.position = "none")

plot1 <- shannon_plot + stat_pvalue_manual(shannon, 
                                           y.position = 2.8, step.increase = 0.05,
                                           label = "sig_p_adj",hide.ns = TRUE)

richness <- dunn_results_filtered %>% filter(Test == "Specnr")%>% filter(sig_p_adj != "ns")
richness_plot <- ggplot(data = df_dmm_diversity, aes(x=dmm, y=specnr, col = dmm)) + 
  geom_boxplot() +
  labs(x = "DMM cluster", y = "Species richenss") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  #scale_x_discrete(labels = c(1:6)) +
  scale_color_manual(values = c("#796248","#a47852","#cd9c58", "#e6cd98","#f8edcf", "#d4ede8")) +
  theme(legend.position = "none")

plot2 <- richness_plot + stat_pvalue_manual(richness, 
                                            y.position = 95, step.increase = 0.05,
                                            label = "sig_p_adj",hide.ns = TRUE)
all <- (plot1+plot2)

# saving figure for the paper ====
ggsave(
  filename   = "S8.3_DMM_alpha_diversity.tiff",
  plot       = all,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

