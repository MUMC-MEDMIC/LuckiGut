# libraries ====
library(microViz)
library(tidyverse)
library(phyloseq)
library(FSA)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# import data ====
ifdp_raw <- read.csv(file = "Input/6_combined_counts_IFDP.csv")
baku_f_I <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# clean phyloseq data ====
keys_sample_name <- as.data.frame(as.matrix(baku_f_I@sam_data))
keys_sample_name <- keys_sample_name  %>%  select (MMHP_SampleID,Family_ID, Age_individual, diet_class,Sample_aliquote)

rank_viz <- "Species"
baku_f_I <- baku_f_I %>% tax_transform("compositional", rank = rank_viz)
reads <- as.data.frame(as.matrix(baku_f_I@otu_table))
rm(baku_f_I,rank_viz)

# read in fiber data ====
rownames(ifdp_raw) <- ifdp_raw$Fibers #make fibers into rownames
ifdp_raw <- ifdp_raw[, -1] #remove fiber column
ifdp_raw <- as.data.frame(t(ifdp_raw)) #make transposon
ifdp_raw$MMHP_SampleID <- rownames(ifdp_raw) #make rownames into a column
ifdp_raw$MMHP_SampleID <- gsub("[.-]", "_", ifdp_raw$MMHP_SampleID) #make the names of the sample the same format
#renaiming some fibers to avoid problems in the future and spelling mistakes
names(ifdp_raw)[names(ifdp_raw) == "Beta-glucan"] <- "Betaglucan"
names(ifdp_raw)[names(ifdp_raw) == "Resistant starch"] <- "Resistant_starch"
names(ifdp_raw)[names(ifdp_raw) == "Xanthan "] <- "Xanthan"
names(ifdp_raw)[names(ifdp_raw) == "Galctomannan"] <- "Galactomannan"
names(ifdp_raw)[names(ifdp_raw) == "Carageenan "] <- "Carrageenan"
names(ifdp_raw)[names(ifdp_raw) == "Xloglucan"] <- "Xyloglucan"

seq_codes <- keys_sample_name$MMHP_SampleID
ifdp_raw <- ifdp_raw[ifdp_raw$MMHP_SampleID %in% seq_codes, ] #selecting fiber data only for the samples that dietary analysis was done to
fibers <- colnames(ifdp_raw)
fibers <- setdiff(fibers, "MMHP_SampleID")

# add fiber reads to data ====
combined_df <- full_join(ifdp_raw,keys_sample_name, by = "MMHP_SampleID")
rm(keys_sample_name,ifdp_raw,seq_codes)

combined_df$Age_individual<- factor(combined_df$Age_individual, 
                                        levels = c("4 months",
                                                   "5 months",
                                                   "6 months",
                                                   "9 months",
                                                   "11 months",
                                                   "14 months"
                                        ))

# Kurska wallis test ====
# Loop over each bacterium name in the vector
# Initialize an empty data frame to store results
results_KW <- data.frame(bacteria = character(),
                         h_statistics = numeric(),#H statistic is asymptotically chi-square distributed, that is why it is chi-squared in R
                         df = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each bacterium name in the vector
for (fiber_name in fibers) {
  # Create the formula for the Kruskal-Wallis test dynamically
  formula <- as.formula(paste(fiber_name, "~ Age_individual"))
  
  # Perform the Kruskal-Wallis test
  result <- kruskal.test(formula, data = combined_df)
  
  # Add the results to the data frame
  results_KW <- rbind(results_KW, data.frame(fiber = fiber_name,
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
fiber_names_after_KW <- results_KW_sig$fiber 

rm(result,fiber_name,formula,results_KW_sig,fibers)

# Dunn test ====
combined_df$Age_individual <- gsub(" months", "", combined_df$Age_individual)
combined_df <- combined_df %>%
  mutate(Age_individual = recode(Age_individual,
                       `4` = 1,
                       `5` = 2,
                       `6` = 3,
                       `9` = 4,
                       `11` = 5,
                       `14` = 6))
combined_df$Age_individual <- as.numeric(combined_df$Age_individual)
combined_df$Age_individual <- as.factor(combined_df$Age_individual)

# Initialize an empty list to store results
dunn_results <- list()

# Loop over each bacterium name in the vector
for (fiber_name in fiber_names_after_KW) {
  # Create the formula for the Dunn's test dynamically
  formula <- as.formula(paste(fiber_name, "~ Age_individual"))
  
  # Perform Dunn's test with adjustment method "bh"
  result <- dunnTest(formula, data = combined_df, method = "bh")
  
  # Store the result in the list
  dunn_results[[fiber_name]] <- result
}

# Create an empty list to store the results
results_dunn <- list()

# Loop through each bacteria name in fiber_names_after_KW
for (fiber_name in fiber_names_after_KW) {
  # Replace "Bifidobacterium_longum" with the current bacteria name
  result_df <- dunn_results[[fiber_name]]$res
  
  # Add the result data frame to the list with the corresponding bacteria name
  results_dunn[[fiber_name]] <- result_df
}

# Merge the list of data frames into a single data frame
results_dunn <- bind_rows(results_dunn, .id = "fiber")

results_dunn <- results_dunn %>% 
  mutate(sig_p_adj = if_else(
    condition = P.adj < 0.05, 
    true = "*", 
    false = ""
  ))

# cleaning results ====
results_dunn <- results_dunn %>% mutate(Comparison = gsub("-", "_", Comparison))

results_dunn <- results_dunn %>% separate(Comparison, into = c("col1", "col2"), sep = "_")
results_dunn <- results_dunn %>%
  mutate(col1 = gsub(" ", "", col1),
         col2 = gsub(" ", "", col2))

results_dunn <- results_dunn %>%  mutate(col1 = recode(col1, 
                       `1` = "4 mon ", 
                       `2` = "5 mon ", 
                       `3` = "6 mon ", 
                       `4` = "9 mon ", 
                       `5` = "11 mon ", 
                       `6` = "14 mon "))

results_dunn <- results_dunn %>%  mutate(col2 = recode(col2, 
                                                       `1` = "- 4 mon", 
                                                       `2` = "- 5 mon", 
                                                       `3` = "- 6 mon", 
                                                       `4` = "- 9 mon", 
                                                       `5` = "- 11 mon", 
                                                       `6` = "- 14 mon"))

results_dunn <- results_dunn %>% unite(Comparison, col1, col2, sep = "")

# saving results ====
#write.table (results_KW, "S4.4_fibers_over_age_results_KW.txt", row.names=FALSE,sep = "\t")
#write.table (results_dunn, "S4.4_fibers_over_age_results_dunn.txt", row.names=FALSE,sep = "\t")




















