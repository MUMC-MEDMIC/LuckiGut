# LOAD LIBRARIES ====
library(readxl) 
library(tidyverse)
library(phyloseq)
library(microViz)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ifdp_raw <- read.csv(file = "Input/6_combined_counts_IFDP.csv")
my_phyloseq <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# cleaning data ====
rownames(ifdp_raw) <- ifdp_raw$Fibers #make fibers into rownames
ifdp_raw <- ifdp_raw[, -1] #remove fiber column
ifdp_raw <- as.data.frame(t(ifdp_raw)) #make transposon
ifdp_raw$MMHP_SampleID <- rownames(ifdp_raw) #make rownames into a column
ifdp_raw$MMHP_SampleID <- gsub("[.-]", "_", ifdp_raw$MMHP_SampleID) #make the names of the sample the same format

#remaining some fibers to avoid problems in the future and spelling mistakes
names(ifdp_raw)[names(ifdp_raw) == "Beta-glucan"] <- "Betaglucan"
names(ifdp_raw)[names(ifdp_raw) == "Resistant starch"] <- "Resistant_starch"
names(ifdp_raw)[names(ifdp_raw) == "Xanthan "] <- "Xanthan"
names(ifdp_raw)[names(ifdp_raw) == "Galctomannan"] <- "Galactomannan"
names(ifdp_raw)[names(ifdp_raw) == "Carageenan "] <- "Carrageenan"
names(ifdp_raw)[names(ifdp_raw) == "Xloglucan"] <- "Xyloglucan"
ifdp_names <- setdiff(colnames(ifdp_raw), "MMHP_SampleID")

# clean phyloseq ====
my_phyloseq <- my_phyloseq %>% ps_join(ifdp_raw, by = "MMHP_SampleID") #adding ifdp to phyloseq
check <- as.data.frame(as.matrix(my_phyloseq@otu_table))

#taking metadata out
full_meta <- as.data.frame(as.matrix(my_phyloseq@sam_data))

#creating only fiber table for correlations
ifdp_names_extra <- c(ifdp_names, "Sample_aliquote","Age_individual","MMHP_SampleID","diet_class")
ifdp <- full_meta[, ifdp_names_extra, drop = FALSE]
convert_to_numeric <- function(df, columns) {
  df[columns] <- lapply(df[columns], function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
  })
  return(df)
}
ifdp <- convert_to_numeric(ifdp, ifdp_names)
rm(ifdp_raw,full_meta,ifdp_names_extra)
ifdp <- ifdp[order(row.names(ifdp)), ]

# arranging fiber data for correlations ====
# subset the fiber profiles per each TP
fiber_4m <- ifdp %>% filter(Age_individual == "4 months") %>% dplyr::select(-c("MMHP_SampleID","Sample_aliquote","Age_individual","diet_class"))
fiber_5m <- ifdp %>% filter(Age_individual == "5 months") %>% dplyr::select(-c("MMHP_SampleID","Sample_aliquote","Age_individual","diet_class"))
fiber_6m <- ifdp %>% filter(Age_individual == "6 months") %>% dplyr::select(-c("MMHP_SampleID","Sample_aliquote","Age_individual","diet_class"))
fiber_9m <- ifdp %>% filter(Age_individual == "9 months") %>% dplyr::select(-c("MMHP_SampleID","Sample_aliquote","Age_individual","diet_class"))
fiber_11m <- ifdp %>% filter(Age_individual == "11 months") %>% dplyr::select(-c("MMHP_SampleID","Sample_aliquote","Age_individual","diet_class"))
fiber_14m <- ifdp %>% filter(Age_individual == "14 months") %>% dplyr::select(-c("MMHP_SampleID","Sample_aliquote","Age_individual","diet_class"))

# filtering bacterial abundance > 0.01 ====
my_phyloseq <- phyloseq_validate(my_phyloseq, remove_undetected = TRUE)
ps_diet_comp <- my_phyloseq %>% tax_transform("compositional") %>% ps_get()

#filter per TP
ps_diet_4 <- ps_diet_comp %>% ps_filter(Age_individual == "4 months") 
ps_diet_5 <- ps_diet_comp %>% ps_filter(Age_individual == "5 months")
ps_diet_6 <- ps_diet_comp %>% ps_filter(Age_individual == "6 months") 
ps_diet_9 <- ps_diet_comp %>% ps_filter(Age_individual == "9 months") 
ps_diet_11 <- ps_diet_comp %>% ps_filter(Age_individual == "11 months") 
ps_diet_14 <- ps_diet_comp %>% ps_filter(Age_individual == "14 months")

#extract OTU table
otu_4 <- ps_diet_4 %>% otu_table() %>% as.data.frame()
otu_5 <- ps_diet_5 %>% otu_table() %>% as.data.frame()
otu_6 <- ps_diet_6 %>% otu_table() %>% as.data.frame()
otu_9 <- ps_diet_9 %>% otu_table() %>% as.data.frame()
otu_11 <- ps_diet_11 %>% otu_table() %>% as.data.frame()
otu_14 <- ps_diet_14 %>% otu_table() %>% as.data.frame()
rm(ps_diet_4,ps_diet_5,ps_diet_6,ps_diet_9,ps_diet_11,ps_diet_14,ps_diet_comp)

#calculate row mean per taxa
otu_4$mean <- rowMeans(otu_4)
otu_5$mean <- rowMeans(otu_5)
otu_6$mean <- rowMeans(otu_6)
otu_9$mean <- rowMeans(otu_9)
otu_11$mean <- rowMeans(otu_11)
otu_14$mean <- rowMeans(otu_14)

#select those bakut that are above 0.01
bakut_4 <- otu_4 %>% filter(mean >= 0.01) %>% rownames()
bakut_5 <- otu_5 %>% filter(mean >= 0.01) %>% rownames()
bakut_6 <- otu_6 %>% filter(mean >= 0.01) %>% rownames()
bakut_9 <- otu_9 %>% filter(mean >= 0.01) %>% rownames()
bakut_11 <- otu_11 %>% filter(mean >= 0.01) %>% rownames()
bakut_14 <- otu_14 %>% filter(mean >= 0.01) %>% rownames()
rm(otu_4,otu_5,otu_6,otu_9,otu_11,otu_14)

#combine into 1 vector
unik_baku <- c(bakut_4,bakut_5,bakut_6,bakut_9,bakut_11,bakut_14)
unik_baku <- sort(unique(unik_baku))
rm(bakut_4,bakut_5,bakut_6,bakut_9,bakut_11,bakut_14)

#remove everything below 1% from phyloseq
my_phyloseq_filtered <- my_phyloseq %>% tax_select(unik_baku)
rm(my_phyloseq)

# arranging fiber data for correlations ====
rank_viz <- "Species"

my_phyloseq_4m <- my_phyloseq_filtered %>% subset_samples(Age_individual == "4 months") %>% tax_transform("compositional", rank = rank_viz)
otu_4m <- data.frame(t(my_phyloseq_4m@otu_table))
otu_4m$sa_code <- rownames(otu_4m)
otu_4m <- subset(otu_4m, select = -c(sa_code))
otu_4m <- otu_4m[order(row.names(otu_4m)), ] #sorting row names in alphabetical order
rm(my_phyloseq_4m)

my_phyloseq_5m <- my_phyloseq_filtered %>% subset_samples(Age_individual == "5 months") %>% tax_transform("compositional", rank = rank_viz)
otu_5m <- data.frame(t(my_phyloseq_5m@otu_table))
otu_5m$sa_code <- rownames(otu_5m)
otu_5m <- subset(otu_5m, select = -c(sa_code))
otu_5m <- otu_5m[order(row.names(otu_5m)), ]
rm(my_phyloseq_5m)

my_phyloseq_6m <- my_phyloseq_filtered %>% subset_samples(Age_individual == "6 months") %>% tax_transform("compositional", rank = rank_viz)
otu_6m <- data.frame(t(my_phyloseq_6m@otu_table))
otu_6m$sa_code <- rownames(otu_6m)
otu_6m <- subset(otu_6m, select = -c(sa_code))
otu_6m <- otu_6m[order(row.names(otu_6m)), ]
rm(my_phyloseq_6m)

my_phyloseq_9m <- my_phyloseq_filtered %>% subset_samples(Age_individual == "9 months") %>% tax_transform("compositional", rank = rank_viz)
otu_9m <- data.frame(t(my_phyloseq_9m@otu_table))
otu_9m$sa_code <- rownames(otu_9m)
otu_9m <- subset(otu_9m, select = -c(sa_code))
otu_9m <- otu_9m[order(row.names(otu_9m)), ]
rm(my_phyloseq_9m)

my_phyloseq_11m <- my_phyloseq_filtered %>% subset_samples(Age_individual == "11 months") %>% tax_transform("compositional", rank = rank_viz)
otu_11m <- data.frame(t(my_phyloseq_11m@otu_table))
otu_11m$sa_code <- rownames(otu_11m)
otu_11m <- subset(otu_11m, select = -c(sa_code))
otu_11m <- otu_11m[order(row.names(otu_11m)), ]
rm(my_phyloseq_11m)

my_phyloseq_14m <- my_phyloseq_filtered %>% subset_samples(Age_individual == "14 months") %>% tax_transform("compositional", rank = rank_viz)
otu_14m <- data.frame(t(my_phyloseq_14m@otu_table))
otu_14m$sa_code <- rownames(otu_14m)
otu_14m <- subset(otu_14m, select = -c(sa_code))
otu_14m <- otu_14m[order(row.names(otu_14m)), ]
rm(my_phyloseq_14m)

# Calculate Spearman correlations, no exact p-value ====
#4 MONTHS
#checking names 
row1 <- rownames(fiber_4m)
row2 <- rownames(otu_4m)
setdiff(fiber,baku)

cor_test_results <- list() # List to store correlation test results

# Loop through each column of otu_4m
for (otu_column in colnames(otu_4m)) {
  # Loop through each column of fiber_4m
  for (fiber_column in colnames(fiber_4m)) {
    # Perform correlation test between the two columns
    cor_result <- cor.test(otu_4m[[otu_column]], fiber_4m[[fiber_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste(otu_column, fiber_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(bacteria = character(),
                        fiber = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract bacteria and fiber from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  bacteria <- comparison[1]
  fiber <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(bacteria = bacteria,
                        fiber = fiber,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}
result_df$ID <- rownames( result_df)
spearm_4m <- result_df
spearm_4m$Age_individual <- '4 months'
rm(otu_4m,fiber_4m,cor_test_results,result_df,cor_result,cor_test,new_row,bacteria, comparison,fiber,fiber_column,otu_column)

#5 MONTHS
cor_test_results <- list() # List to store correlation test results

# Loop through each column of otu_5m
for (otu_column in colnames(otu_5m)) {
  # Loop through each column of fiber_5m
  for (fiber_column in colnames(fiber_5m)) {
    # Perform correlation test between the two columns
    cor_result <- cor.test(otu_5m[[otu_column]], fiber_5m[[fiber_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste(otu_column, fiber_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(bacteria = character(),
                        fiber = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract bacteria and fiber from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  bacteria <- comparison[1]
  fiber <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(bacteria = bacteria,
                        fiber = fiber,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}
result_df$ID <- rownames( result_df)
spearm_5m <- result_df
spearm_5m$Age_individual <- '5 months'
rm(otu_5m,fiber_5m,cor_test_results,result_df,cor_result,cor_test,new_row,bacteria, comparison,fiber,fiber_column,otu_column)

#6 MONTHS
cor_test_results <- list() # List to store correlation test results

# Loop through each column of otu_6m
for (otu_column in colnames(otu_6m)) {
  # Loop through each column of fiber_6m
  for (fiber_column in colnames(fiber_6m)) {
    # Perform correlation test between the two columns
    cor_result <- cor.test(otu_6m[[otu_column]], fiber_6m[[fiber_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste(otu_column, fiber_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(bacteria = character(),
                        fiber = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract bacteria and fiber from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  bacteria <- comparison[1]
  fiber <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(bacteria = bacteria,
                        fiber = fiber,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}
result_df$ID <- rownames( result_df)
spearm_6m <- result_df
spearm_6m$Age_individual <- '6 months'
rm(otu_6m,fiber_6m,cor_test_results,result_df,cor_result,cor_test,new_row,bacteria, comparison,fiber,fiber_column,otu_column)

#9 MONTHS
cor_test_results <- list() # List to store correlation test results

# Loop through each column of otu_9m
for (otu_column in colnames(otu_9m)) {
  # Loop through each column of fiber_9m
  for (fiber_column in colnames(fiber_9m)) {
    # Perform correlation test between the two columns
    cor_result <- cor.test(otu_9m[[otu_column]], fiber_9m[[fiber_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste(otu_column, fiber_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(bacteria = character(),
                        fiber = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract bacteria and fiber from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  bacteria <- comparison[1]
  fiber <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(bacteria = bacteria,
                        fiber = fiber,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}
result_df$ID <- rownames( result_df)
spearm_9m <- result_df
spearm_9m$Age_individual <- '9 months'
rm(otu_9m,fiber_9m,cor_test_results,result_df,cor_result,cor_test,new_row,bacteria, comparison,fiber,fiber_column,otu_column)

#11 MONTHS
cor_test_results <- list() # List to store correlation test results

# Loop through each column of otu_11m
for (otu_column in colnames(otu_11m)) {
  # Loop through each column of fiber_11m
  for (fiber_column in colnames(fiber_11m)) {
    # Perform correlation test between the two columns
    cor_result <- cor.test(otu_11m[[otu_column]], fiber_11m[[fiber_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste(otu_column, fiber_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(bacteria = character(),
                        fiber = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract bacteria and fiber from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  bacteria <- comparison[1]
  fiber <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(bacteria = bacteria,
                        fiber = fiber,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}
result_df$ID <- rownames( result_df)
spearm_11m <- result_df
spearm_11m$Age_individual <- '11 months'
rm(otu_11m,fiber_11m,cor_test_results,result_df,cor_result,cor_test,new_row,bacteria, comparison,fiber,fiber_column,otu_column)

#14 MONTHS
cor_test_results <- list() # List to store correlation test results

# Loop through each column of otu_14m
for (otu_column in colnames(otu_14m)) {
  # Loop through each column of fiber_14m
  for (fiber_column in colnames(fiber_14m)) {
    # Perform correlation test between the two columns
    cor_result <- cor.test(otu_14m[[otu_column]], fiber_14m[[fiber_column]], method = "spearman", exact = FALSE)
    
    # Store the correlation test result along with column names
    cor_test_results[[paste(otu_column, fiber_column, sep = "_vs_")]] <- cor_result
  }
}

# Create an empty dataframe to store the results
result_df <- data.frame(bacteria = character(),
                        fiber = character(),
                        p_value = numeric(),
                        estimate = numeric(),
                        statistic = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each comparison result
for (key in names(cor_test_results)) {
  # Extract correlation test result
  cor_test <- cor_test_results[[key]]
  
  # Extract bacteria and fiber from the key
  comparison <- strsplit(key, "_vs_")[[1]]
  bacteria <- comparison[1]
  fiber <- comparison[2]
  
  # Create a new row for the result dataframe
  new_row <- data.frame(bacteria = bacteria,
                        fiber = fiber,
                        p_value = cor_test$p.value,
                        estimate = cor_test$estimate,
                        statistic = cor_test$statistic,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the result dataframe
  result_df <- rbind(result_df, new_row)
}
result_df$ID <- rownames( result_df)
spearm_14m <- result_df
spearm_14m$Age_individual <- '14 months'
rm(otu_14m,fiber_14m,cor_test_results,result_df,cor_result,cor_test,new_row,bacteria, comparison,fiber,fiber_column,otu_column)

# FDR correction and saving ====
#FDR adjustment
spearm_4m$p_value_adj <- p.adjust(spearm_4m$p_value, method = "BH")
spearm_5m$p_value_adj <- p.adjust(spearm_5m$p_value, method = "BH")
spearm_6m$p_value_adj <- p.adjust(spearm_6m$p_value, method = "BH")
spearm_9m$p_value_adj <- p.adjust(spearm_9m$p_value, method = "BH")
spearm_11m$p_value_adj <- p.adjust(spearm_11m$p_value, method = "BH")
spearm_14m$p_value_adj <- p.adjust(spearm_14m$p_value, method = "BH")

spearman_all <- rbind(spearm_4m,spearm_5m,spearm_6m,spearm_9m,spearm_11m,spearm_14m)

spearman_all <- spearman_all %>% 
  mutate(sig_padj = if_else(
    condition = p_value_adj < 0.05, 
    true = "*", 
    false = "")
  )

spearman_all_fixed <- spearman_all %>%
  mutate(bacteria = gsub("_", " ", bacteria))
spearman_all_fixed$ID <- NULL
names(spearman_all_fixed)[names(spearman_all_fixed) == "estimate"] <- "Correlation_coefficient_r"
names(spearman_all_fixed)[names(spearman_all_fixed) == "statistic"] <- "T_statistics"

# saving results ====
#write.table (spearman_all_fixed, "S4.1_baku_vs_fiber_per_TP_spearman_adjust.txt", row.names=FALSE,sep = "\t")

