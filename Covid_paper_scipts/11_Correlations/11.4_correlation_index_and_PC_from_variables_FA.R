# LOAD LIBRARIES ====
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggplot2)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
index <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
scores <- read.delim("Output/Files/S3.4_PCA_index_factor_analysis_scores_from_FA.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data ====
index$Family_ID <- rownames(index)

scores_PC <- scores %>% select(Factor1,Factor2,Factor3,Factor4,Factor5)
scores_PC$Family_ID <- rownames(scores_PC)

scores_PC <- scores_PC[order(row.names(scores_PC)), ] # A -> Z
index <- index[order(row.names(index)), ] 

scat <- full_join(index,scores_PC, by = "Family_ID")

rownames(scores_PC) == rownames(index)

#removing families with NAs
remove_values <- c("1228","1225")
index <- subset(index, !Family_ID %in% remove_values)

index <- index$index_m_h 
scores_PC$Family_ID <- NULL

rm(scores,remove_values)

# Spearman correlations, no exact p-value ====
# Loop through each column in df2 and calculate Spearman correlation with X
cor_test_results <- list()

# Loop through each column in scores_PC
for (score_column in colnames(scores_PC)) {
  # Perform correlation test between the index column and the current column in scores_PC
  cor_result <- cor.test(index, scores_PC[[score_column]], method = "spearman", exact = FALSE)
  
  # Store the correlation test result along with column names
  cor_test_results[[paste("index", score_column, sep = "_vs_")]] <- cor_result
}

# Create an empty dataframe to store the results
result_df <- data.frame(index = character(),
                        pc = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract food item from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  pc <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(index = "index",
                        pc = pc,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}

rm(cor_result,cor_test,new_row,cor_test_results,pc,score_column,comparison,key)

# cleaning results ====
result_df$p_value_adj <- p.adjust(result_df$p_value, method = "BH")
result_df <- result_df %>% 
  mutate(sig_padj = if_else(
    condition = p_value_adj < 0.05, 
    true = "*", 
    false = "")
  )

pvalue <- round(result_df$p_value[1],3)
rho <- round(result_df$estimate[1],3)

# saving data ====
#write.table (result_df, "11.3_spearman_index_scores_correlation.txt", row.names=FALSE,sep = "\t")

# scatter plot ====
plot <- ggplot(scat, aes(x=index_m_h, y=Factor1)) + geom_point()+ labs(title = "")+
  ylab("PC1") + xlab("Exposure index") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15))+
  annotate("text", x = Inf, y = Inf, hjust = 8.1, vjust = 1.1, 
           label = paste("\nrho =", rho, "\np =", pvalue),   
           size = 4, fontface = "bold")+
  ggtitle("Spearman correlation")

plot

# saving figures for the paper ====
#ggsave("S11.3_index_spearman_PC1_correlation.pdf", plot, device = "pdf",width = 300, height = 200, dpi = 300, units = "mm")

#save both figures with dummy legend and without, later copy-paste legend instead of the dummy figure
ggsave(
  filename   = "S11.3_index_spearman_PC1_correlation.tiff",
  plot       = plot,
  width      = 170,      # full page width
  height     = 90,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

