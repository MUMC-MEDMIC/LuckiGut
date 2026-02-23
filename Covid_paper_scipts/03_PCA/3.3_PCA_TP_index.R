# LOAD LIBRARIES ====
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index_mh <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# data cleaning ====
index_mh$Family_ID <- rownames(index_mh)
fam_mh <- sort(index_mh$Family_ID)

infant_covid_yes <- infant %>%ps_filter(covid_270220 == "yes") 
meta <- as.data.frame(as.matrix(infant_covid_yes@sam_data))
meta <- full_join(index_mh,meta, by = "Family_ID",multiple = "all")

index_mean <- mean(index_mh$index_m_h)
index_mh <- index_mh %>% mutate(index_bi = ifelse(index_m_h > index_mean, "above", "below"))
index_mh$Family_ID <- as.numeric(index_mh$Family_ID)

updated_phyloseq <- ps_join(infant_covid_yes,index_mh,by = "Family_ID")

#removing not needed TP
tp <- c("5 months","6 months","9 months","11 months","14 months")
meta <- meta %>% filter(Age_individual %in% tp)
meta$index_m_h <- as.numeric(meta$index_m_h)

col_needed <- c("PC1","PC2")

keys <- meta %>% select(MMHP_SampleID,Family_ID)

# taking PCA scores out ====
infant_covid_yes <- infant_covid_yes %>%ps_filter(Age_individual == "5 months" | Age_individual == "6 months" | Age_individual == "9 months"| Age_individual == "11 months"| Age_individual == "14 months") 

sample_data_df <- as.data.frame(sample_data(infant))# Extract the sample data
sample_data_df$Family_ID <- as.character(sample_data_df$Family_ID)# Convert column 'Z' from integer to character
sample_data(infant_covid_yes) <- sample_data(sample_data_df) # Update the sample data in the phyloseq object
rm(sample_data_df)

# 5 months 
pca <- infant_covid_yes %>% ps_filter(Age_individual == "5 months") %>%  tax_transform("clr", rank = "Species") %>% ord_calc(method = "PCA")
ord_obj <- pca@ord
scores5 <- as.data.frame(vegan::scores(ord_obj, display = "sites")) # Get the scores (site scores)
#loadings <- as.data.frame(vegan::scores(ord_obj, display = "species")) # Get the loadings (species scores)
rm(ord_obj,pca)

# 6 months 
pca <- infant_covid_yes %>% ps_filter(Age_individual == "6 months") %>%  tax_transform("clr", rank = "Species") %>% ord_calc(method = "PCA")
ord_obj <- pca@ord
scores6 <- as.data.frame(vegan::scores(ord_obj, display = "sites")) # Get the scores (site scores)
#loadings <- as.data.frame(vegan::scores(ord_obj, display = "species")) # Get the loadings (species scores)
rm(ord_obj,pca)

# 9 months 
pca <- infant_covid_yes %>% ps_filter(Age_individual == "9 months") %>% tax_transform("clr", rank = "Species") %>% ord_calc(method = "PCA")
ord_obj <- pca@ord
scores9 <- as.data.frame(vegan::scores(ord_obj, display = "sites")) # Get the scores (site scores)
rm(ord_obj,pca)

# 11 months 
pca <- infant_covid_yes %>% ps_filter(Age_individual == "11 months") %>% tax_transform("clr", rank = "Species") %>% ord_calc(method = "PCA")
ord_obj <- pca@ord
scores11 <- as.data.frame(vegan::scores(ord_obj, display = "sites")) # Get the scores (site scores)
rm(ord_obj,pca)

# 14 months 
pca <- infant_covid_yes %>% ps_filter(Age_individual == "14 months") %>% tax_transform("clr", rank = "Species") %>% ord_calc(method = "PCA")
ord_obj <- pca@ord
scores14 <- as.data.frame(vegan::scores(ord_obj, display = "sites")) # Get the scores (site scores)
rm(ord_obj,pca)

# adding family code ====
scores5$MMHP_SampleID <- rownames(scores5)
scores6$MMHP_SampleID <- rownames(scores6)
scores9$MMHP_SampleID <- rownames(scores9)
scores11$MMHP_SampleID <- rownames(scores11)
scores14$MMHP_SampleID <- rownames(scores14)

scores5 <- left_join(scores5,keys, by = "MMHP_SampleID")
scores5 <- scores5[!is.na(scores5$Family_ID), ]
scores6 <- left_join(scores6,keys, by = "MMHP_SampleID")
scores6 <- scores6[!is.na(scores6$Family_ID), ]
scores9 <- left_join(scores9,keys, by = "MMHP_SampleID")
scores9 <- scores9[!is.na(scores9$Family_ID), ]
scores11 <- left_join(scores11,keys, by = "MMHP_SampleID")
scores11 <- scores11[!is.na(scores11$Family_ID), ]
scores14 <- left_join(scores14,keys, by = "MMHP_SampleID")
scores14 <- scores14[!is.na(scores14$Family_ID), ]

scores_all <- rbind(scores5,scores6,scores9,scores11,scores14)

scores5$MMHP_SampleID <- NULL
scores14$MMHP_SampleID <- NULL
scores6$MMHP_SampleID <- NULL
scores9$MMHP_SampleID <- NULL
scores11$MMHP_SampleID <- NULL
scores_all$Family_ID <- NULL

# separating index by TP ====
meta5 <- meta %>% filter(Age_individual == "5 months")
meta6 <- meta %>% filter(Age_individual == "6 months")
meta9 <- meta %>% filter(Age_individual == "9 months")
meta11 <- meta %>% filter(Age_individual == "11 months")
meta14 <- meta %>% filter(Age_individual == "14 months")

# correlations m-h ====
#5 months
index <- meta5$index_m_h
cor_test_results <- list()

# Loop through each column in scores5 that is in col_needed
for (beta_column in col_needed) {
  # Check if the column exists in scores5
  if (beta_column %in% colnames(scores5)) {
    # Perform correlation test between the index column and the current column in scores5
    cor_result <- cor.test(index, scores5[[beta_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste("index", beta_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(index = character(),
                        p_value = numeric(),
                        rho = numeric(),
                        sum_of_squared_rank_differences = numeric(),
                        test_type = character(),
                        method = character(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract beta_diversity from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  beta_diversity <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(index = "index_m_h",
                        beta_diversity = beta_diversity,
                        p_value = cor_test$p.value,
                        rho = cor_test$estimate,
                        sum_of_squared_rank_differences = cor_test$statistic,
                        test_type = cor_test$alternative,
                        method = cor_test$method,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}

spearm_5m <- result_df
spearm_5m$Age_individual <- '5 months'

rm(index,cor_test_results,key,result_df,cor_result,cor_test,new_row,beta_column,beta_diversity,comparison)

#6 months
index <- meta6$index_m_h
cor_test_results <- list()

# Loop through each column in scores6 that is in col_needed
for (beta_column in col_needed) {
  # Check if the column exists in scores6
  if (beta_column %in% colnames(scores6)) {
    # Perform correlation test between the index column and the current column in scores6
    cor_result <- cor.test(index, scores6[[beta_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste("index", beta_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(index = character(),
                        p_value = numeric(),
                        rho = numeric(),
                        sum_of_squared_rank_differences = numeric(),
                        test_type = character(),
                        method = character(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract beta_diversity from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  beta_diversity <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(index = "index_m_h",
                        beta_diversity = beta_diversity,
                        p_value = cor_test$p.value,
                        rho = cor_test$estimate,
                        sum_of_squared_rank_differences = cor_test$statistic,
                        test_type = cor_test$alternative,
                        method = cor_test$method,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}

spearm_6m <- result_df
spearm_6m$Age_individual <- '6 months'

rm(index,cor_test_results,key,result_df,cor_result,cor_test,new_row,beta_column,beta_diversity,comparison)

#9 months
index <- meta9$index_m_h
cor_test_results <- list()

# Loop through each column in scores9 that is in col_needed
for (beta_column in col_needed) {
  # Check if the column exists in scores9
  if (beta_column %in% colnames(scores9)) {
    # Perform correlation test between the index column and the current column in scores9
    cor_result <- cor.test(index, scores9[[beta_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste("index", beta_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(index = character(),
                        p_value = numeric(),
                        rho = numeric(),
                        sum_of_squared_rank_differences = numeric(),
                        test_type = character(),
                        method = character(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract beta_diversity from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  beta_diversity <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(index = "index_m_h",
                        beta_diversity = beta_diversity,
                        p_value = cor_test$p.value,
                        rho = cor_test$estimate,
                        sum_of_squared_rank_differences = cor_test$statistic,
                        test_type = cor_test$alternative,
                        method = cor_test$method,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}

spearm_9m <- result_df
spearm_9m$Age_individual <- '9 months'

rm(index,cor_test_results,key,result_df,cor_result,cor_test,new_row,beta_column,beta_diversity,comparison)

#11 months
index <- meta11$index_m_h
cor_test_results <- list()

# Loop through each column in scores11 that is in col_needed
for (beta_column in col_needed) {
  # Check if the column exists in scores11
  if (beta_column %in% colnames(scores11)) {
    # Perform correlation test between the index column and the current column in scores11
    cor_result <- cor.test(index, scores11[[beta_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste("index", beta_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(index = character(),
                        p_value = numeric(),
                        rho = numeric(),
                        sum_of_squared_rank_differences = numeric(),
                        test_type = character(),
                        method = character(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract beta_diversity from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  beta_diversity <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(index = "index_m_h",
                        beta_diversity = beta_diversity,
                        p_value = cor_test$p.value,
                        rho = cor_test$estimate,
                        sum_of_squared_rank_differences = cor_test$statistic,
                        test_type = cor_test$alternative,
                        method = cor_test$method,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}

spearm_11m <- result_df
spearm_11m$Age_individual <- '11 months'

rm(index,cor_test_results,key,result_df,cor_result,cor_test,new_row,beta_column,beta_diversity,comparison)

#14 months
index <- meta14$index_m_h
cor_test_results <- list()

# Loop through each column in scores14 that is in col_needed
for (beta_column in col_needed) {
  # Check if the column exists in scores14
  if (beta_column %in% colnames(scores14)) {
    # Perform correlation test between the index column and the current column in scores14
    cor_result <- cor.test(index, scores14[[beta_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste("index", beta_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(index = character(),
                        p_value = numeric(),
                        rho = numeric(),
                        sum_of_squared_rank_differences = numeric(),
                        test_type = character(),
                        method = character(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract beta_diversity from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  beta_diversity <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(index = "index_m_h",
                        beta_diversity = beta_diversity,
                        p_value = cor_test$p.value,
                        rho = cor_test$estimate,
                        sum_of_squared_rank_differences = cor_test$statistic,
                        test_type = cor_test$alternative,
                        method = cor_test$method,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}

spearm_14m <- result_df
spearm_14m$Age_individual <- '14 months'

rm(index,cor_test_results,key,result_df,cor_result,cor_test,new_row,beta_column,beta_diversity,comparison)

# merging results m-h ====
spearm_5m$p_value_adj <- p.adjust(spearm_5m$p_value, method = "BH")
spearm_6m$p_value_adj <- p.adjust(spearm_6m$p_value, method = "BH")
spearm_9m$p_value_adj <- p.adjust(spearm_9m$p_value, method = "BH")
spearm_11m$p_value_adj <- p.adjust(spearm_11m$p_value, method = "BH")
spearm_14m$p_value_adj <- p.adjust(spearm_14m$p_value, method = "BH")

spearman_all_mh <- rbind(spearm_5m, spearm_6m,spearm_9m,spearm_11m,spearm_14m)

rm(spearm_5m, spearm_6m,spearm_9m,spearm_11m,spearm_14m)

spearman_all_mh <- spearman_all_mh %>% 
  mutate(sig_padj = if_else(
    condition = p_value_adj < 0.05, 
    true = "*", 
    false = "")
  )

rm(scores5,scores6,scores9,scores11,scores14,meta5,meta6,meta9,meta11,meta14)

# merging scores and index =====
meta <- full_join(scores_all,meta, by = "MMHP_SampleID")

# scatter plots m-h ====
ggplot(meta, aes(x=index_m_h, y=PC1)) + geom_point()+
  facet_wrap(~Age_individual)+
  geom_smooth(method = "lm", se = TRUE, color = "red")+
  ggtitle("Correlation between index m-h and PC1")

ggplot(meta, aes(x=index_m_h, y=PC2)) + geom_point()+
  facet_wrap(~Age_individual)+
  geom_smooth(method = "lm", se = TRUE, color = "red")+
  ggtitle("Correlation between index m-h and PC2")

# PCA plots =====
#without NA only 6-14 months, and age
tax_table(updated_phyloseq)[,"Species"] <- gsub("_", " ", tax_table(updated_phyloseq)[,"Species"])

plot_wih_legend <- updated_phyloseq %>% 
  ps_filter(Age_individual %in% c("5 months", "6 months", "9 months","11 months","14 months")) %>% 
  subset_samples(!is.na(index_m_h)) %>%  # Filter out NA values from column_A
  tax_transform("clr", rank = "Species") %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "index_m_h",auto_caption = NA,
           plot_taxa = 1:5, size = 0.8, 
           tax_lab_style = constraint_lab_style(size = 1.9, colour = "black",fontface = "italic"),
           axes = c(1, 2)) +
  labs(color = "Exposure index")+theme_bw()+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.position = "bottom",
        legend.key.height = unit(2, "mm"),
        legend.key.width  = unit(2, "mm"),
        legend.spacing.y  = unit(2, "mm"),   # vertical space between items
        legend.spacing.x  = unit(2, "mm"),
        legend.margin     = margin(0, 0, 0, 0, "mm"),
        legend.box.margin = margin(0, 0, 0, 0, "mm"))+
  ggtitle(NULL)

plot_wihout_legend <- updated_phyloseq %>% 
  ps_filter(Age_individual %in% c("5 months", "6 months", "9 months","11 months","14 months")) %>% 
  subset_samples(!is.na(index_m_h)) %>%
  tax_transform("clr", rank = "Species") %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(
    color = "index_m_h",
    auto_caption = NA,
    plot_taxa = 1:5,
    size = 0.8, 
    tax_vec_length = 1.5,      # Shortens arrow length (default ~1)
    tax_lab_length = 1.7,     #Positions labels closer to arrows (slightly longer 
    tax_lab_style = constraint_lab_style(size = 1.9, colour = "black", fontface = "italic", angle = 20),
    axes = c(1, 2)
  ) +
  ggtitle(NULL) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 4),
    axis.text.y  = element_text(size = 4),
    axis.title.x = element_text(size = 6), 
    axis.title.y = element_text(size = 6),
    legend.position = "none"   # hides the legend
  )

# saving figures for the paper ====
ggsave(
  filename   = "S3.3_PCA_TP_index_with_legend.svg",
  plot       = plot_wih_legend,
  width    = 85,
  height   = 90,
  units    = "mm",
  device   = svglite
)

ggsave(
  filename   = "S3.3_PCA_TP_index_without_legend.svg",
  plot       = plot_wihout_legend,
  width    = 85,
  height   = 90,
  units    = "mm",
  device   = svglite
)

