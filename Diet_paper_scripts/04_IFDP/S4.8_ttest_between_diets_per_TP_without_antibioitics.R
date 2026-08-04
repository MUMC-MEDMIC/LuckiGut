# 0.LOAD LIBRARIES ====
library(readxl) 
library(tidyverse)
library(phyloseq)
library(reshape2)
library(ComplexHeatmap)
library(plyr)
library(circlize) # for colorRamp2
require(lme4)
require(lmerTest)
library(microViz)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(svglite)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ifdp_raw <- read.csv(file = "Input/6_combined_counts_IFDP.csv")
my_phyloseq <- readRDS("Output/Files/S1.6_cleaned_phyloseq_object_MS_without_antibiotics.rds")

# cleaning data ====
rownames(ifdp_raw) <- ifdp_raw$Fibers #make fibers into rownames
ifdp_raw <- ifdp_raw[, -1] #remove fiber column
ifdp_raw <- as.data.frame(t(ifdp_raw)) #make transposon
ifdp_raw$MMHP_SampleID <- rownames(ifdp_raw) #make rownames into a column
ifdp_raw$MMHP_SampleID <- gsub("[.-]", "_", ifdp_raw$MMHP_SampleID) #make the names of the sample the same format

full_meta <- as.data.frame(as.matrix(my_phyloseq@sam_data))
full_meta <- full_meta[, c("MMHP_SampleID", "Sample_aliquote")] #select only sample names -> keys
seq_codes <- full_meta$MMHP_SampleID
ifdp <- ifdp_raw[ifdp_raw$MMHP_SampleID %in% seq_codes, ] #selecting fiber data only for the samples that dietary analysis was done to

my_phyloseq <- my_phyloseq %>% 
  ps_join(ifdp, by = "MMHP_SampleID")
fiber <- as.data.frame(samdat_tbl(my_phyloseq))
rm(my_phyloseq,ifdp_raw,full_meta,ifdp,seq_codes)

#renaiming some fibers to avoid problems in the future and spelling mistakes
names(fiber)[names(fiber) == "Beta-glucan"] <- "Betaglucan"
names(fiber)[names(fiber) == "Resistant starch"] <- "Resistant_starch"
names(fiber)[names(fiber) == "Xanthan "] <- "Xanthan"
names(fiber)[names(fiber) == "Galctomannan"] <- "Galactomannan"
names(fiber)[names(fiber) == "Carageenan "] <- "Carrageenan"
names(fiber)[names(fiber) == "Xloglucan"] <- "Xyloglucan"

sort(colnames(fiber))

col_to_keep <- c("Xanthan","Betaglucan","Cellulose","Xyloglucan","Dextran","Galactomannan","Galactoglucomannan",
                 "Glucomannan","Mannan","Carageenan","Galactan","Arabinogalactan","Arabinan","Rhamnogalacturonan",
                 "Pectin","Xylan","Arabinoxylan","Resistant_starch","Laminaran","Gellan","Inulin","Levan",
                 "Chitin","Alginate","diet_class","timepoint","Sample_aliquote")
fiber_filt <- fiber[, names(fiber) %in% col_to_keep]

#alphabeta ordering of fibers
sorted_col_names <- names(fiber_filt)[order(names(fiber_filt))] # Order column names alphabetically
fiber_filt <- fiber_filt[, sorted_col_names] # Reorder the columns
columns_to_move <- c("diet_class", "timepoint","Sample_aliquote") # Identify the columns to move to the end
fiber_columns <- setdiff(names(fiber_filt), columns_to_move) # Determine the remaining columns
new_order <- c(fiber_columns, columns_to_move) # Reorder the columns
fiber_filt <- fiber_filt[new_order]
rownames(fiber_filt) <- fiber_filt$Sample_aliquote
rm(col_to_keep,columns_to_move,new_order,sorted_col_names,fiber)

#creating varibales for the t-test
fiber_filt$CLASS_1 <- ifelse(fiber_filt$diet_class == "1", "1", "0") 
fiber_filt$CLASS_2 <- ifelse(fiber_filt$diet_class == "2", "1", "0")  
fiber_filt$CLASS_3 <- ifelse(fiber_filt$diet_class == "3", "1", "0")
custom_levels <- c("1", "0")
fiber_filt$CLASS_1 <- factor(fiber_filt$CLASS_1, levels = custom_levels)
fiber_filt$CLASS_2 <- factor(fiber_filt$CLASS_2, levels = custom_levels)
fiber_filt$CLASS_3 <- factor(fiber_filt$CLASS_3, levels = custom_levels)

str(fiber_filt)

# 3.Calculate statistics and t-test between diet classes at 4 months ====
selected_month <- fiber_filt %>% filter(timepoint == "4m")

# 3.1 bartlett test homogeneity of variances ====
result_bartlett <- data.frame(fiber = character(),
                              K_squared_bartlett = numeric(),
                              df_bartlett = integer(),
                              p_value_bartlett = numeric(),
                              stringsAsFactors = FALSE)
for (col_name in fiber_columns) {
  result <- bartlett.test(selected_month[[col_name]], selected_month$diet_class)
  result_bartlett <- rbind(result_bartlett, data.frame(fiber = col_name,
                                                       K_squared_bartlett = result$statistic,
                                                       df_bartlett = result$parameter,
                                                       p_value_bartlett = result$p.value))
}
rm(result,col_name)

# 3.2 shapiro test normality ====
# Create an empty data frame to store results
result_shapiro <- data.frame(fiber = character(),
                              W_statistic_shapiro = numeric(),
                              p_value_shapiro = numeric(),
                              stringsAsFactors = FALSE)

# Perform Shapiro-Wilk test for each row
for (col_name in fiber_columns) {
  # Perform Shapiro-Wilk test
  result <- shapiro.test(selected_month[[col_name]])
  
  # Store results in dataframe
  result_shapiro <- rbind(result_shapiro, data.frame(fiber = col_name,
                                                     W_value = result$statistic,
                                                     p_value_shapiro = result$p.value))
}

barn_shapiro <- full_join(result_shapiro,result_bartlett, by = "fiber")
barn_shapiro$Age <- "4 months"
rm(result,result_bartlett,result_shapiro,col_name)

# 3.3 t-test ====
# CLASS 1 vs other
# Preallocate the data frame with the maximum possible number of rows
class1 <- data.frame(fibers = character(), 
                        P_Value = numeric(), 
                        T_Statistic = numeric(), 
                        stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_1"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class1 <- rbind(class1, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class1$Adjusted_P_Value <- adjusted_p_values
class1$class <- "class1"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 2 vs other
# Preallocate the data frame with the maximum possible number of rows
class2 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_2"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class2 <- rbind(class2, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class2$Adjusted_P_Value <- adjusted_p_values
class2$class <- "class2"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 3 vs other
# Preallocate the data frame with the maximum possible number of rows
class3 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_3"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class3 <- rbind(class3, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class3$Adjusted_P_Value <- adjusted_p_values
class3$class <- "class3"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

#merging results together
result_t_test <- rbind(class1,class2,class3)
result_t_test$Age <- "4 months"
rm(class1,class2,class3,selected_month)

# 4.Calculate statistics and t-test between diet classes at 5 months ====
selected_month <- fiber_filt %>% filter(timepoint == "5m")

# 4.1 bartlett test homogeneity of variances ====
result_bartlett <- data.frame(fiber = character(),
                              K_squared_bartlett = numeric(),
                              df_bartlett = integer(),
                              p_value_bartlett = numeric(),
                              stringsAsFactors = FALSE)
for (col_name in fiber_columns) {
  result <- bartlett.test(selected_month[[col_name]], selected_month$diet_class)
  result_bartlett <- rbind(result_bartlett, data.frame(fiber = col_name,
                                                       K_squared_bartlett = result$statistic,
                                                       df_bartlett = result$parameter,
                                                       p_value_bartlett = result$p.value))
}
rm(result,col_name)

# 4.2 shapiro test normality ====
# Create an empty data frame to store results
result_shapiro <- data.frame(fiber = character(),
                             W_statistic_shapiro = numeric(),
                             p_value_shapiro = numeric(),
                             stringsAsFactors = FALSE)

# Perform Shapiro-Wilk test for each row
for (col_name in fiber_columns) {
  # Perform Shapiro-Wilk test
  result <- shapiro.test(selected_month[[col_name]])
  
  # Store results in dataframe
  result_shapiro <- rbind(result_shapiro, data.frame(fiber = col_name,
                                                     W_value = result$statistic,
                                                     p_value_shapiro = result$p.value))
}

barn_shapiro_semi <- full_join(result_shapiro,result_bartlett, by = "fiber")
barn_shapiro_semi$Age <- "5 months"
barn_shapiro <- rbind(barn_shapiro,barn_shapiro_semi)
rm(result,result_bartlett,result_shapiro,col_name,barn_shapiro_semi)

# 4.3 t-test ====
# CLASS 1 vs other
# Preallocate the data frame with the maximum possible number of rows
class1 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_1"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class1 <- rbind(class1, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class1$Adjusted_P_Value <- adjusted_p_values
class1$class <- "class1"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 2 vs other
# Preallocate the data frame with the maximum possible number of rows
class2 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_2"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class2 <- rbind(class2, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class2$Adjusted_P_Value <- adjusted_p_values
class2$class <- "class2"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 3 vs other
# Preallocate the data frame with the maximum possible number of rows
class3 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_3"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class3 <- rbind(class3, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class3$Adjusted_P_Value <- adjusted_p_values
class3$class <- "class3"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

#merging results together
result_t_test_semi <- rbind(class1,class2,class3)
result_t_test_semi$Age <- "5 months"
result_t_test <- rbind(result_t_test,result_t_test_semi)
rm(class1,class2,class3,selected_month,result_t_test_semi)

# 5.Calculate statistics and t-test between diet classes at 6 months ====
selected_month <- fiber_filt %>% filter(timepoint == "6m")

# 5.1 bartlett test homogeneity of variances ====
result_bartlett <- data.frame(fiber = character(),
                              K_squared_bartlett = numeric(),
                              df_bartlett = integer(),
                              p_value_bartlett = numeric(),
                              stringsAsFactors = FALSE)
for (col_name in fiber_columns) {
  result <- bartlett.test(selected_month[[col_name]], selected_month$diet_class)
  result_bartlett <- rbind(result_bartlett, data.frame(fiber = col_name,
                                                       K_squared_bartlett = result$statistic,
                                                       df_bartlett = result$parameter,
                                                       p_value_bartlett = result$p.value))
}
rm(result,col_name)

# 5.2 shapiro test normality tut ====
# Create an empty data frame to store results
result_shapiro <- data.frame(fiber = character(),
                             W_statistic_shapiro = numeric(),
                             p_value_shapiro = numeric(),
                             stringsAsFactors = FALSE)

# Perform Shapiro-Wilk test for each row
for (col_name in fiber_columns) {
  # Perform Shapiro-Wilk test
  result <- shapiro.test(selected_month[[col_name]])
  
  # Store results in dataframe
  result_shapiro <- rbind(result_shapiro, data.frame(fiber = col_name,
                                                     W_value = result$statistic,
                                                     p_value_shapiro = result$p.value))
}

barn_shapiro_semi <- full_join(result_shapiro,result_bartlett, by = "fiber")
barn_shapiro_semi$Age <- "6 months"
barn_shapiro <- rbind(barn_shapiro,barn_shapiro_semi)
rm(result,result_bartlett,result_shapiro,col_name,barn_shapiro_semi)

# 5.3 t-test ====
# CLASS 1 vs other
# Preallocate the data frame with the maximum possible number of rows
class1 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_1"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class1 <- rbind(class1, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class1$Adjusted_P_Value <- adjusted_p_values
class1$class <- "class1"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 2 vs other
# Preallocate the data frame with the maximum possible number of rows
class2 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_2"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class2 <- rbind(class2, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class2$Adjusted_P_Value <- adjusted_p_values
class2$class <- "class2"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 3 vs other
# Preallocate the data frame with the maximum possible number of rows
class3 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_3"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class3 <- rbind(class3, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class3$Adjusted_P_Value <- adjusted_p_values
class3$class <- "class3"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

#merging results together
result_t_test_semi <- rbind(class1,class2,class3)
result_t_test_semi$Age <- "6 months"
result_t_test <- rbind(result_t_test,result_t_test_semi)
rm(class1,class2,class3,selected_month,result_t_test_semi)

# 6.Calculate statistics and t-test between diet classes at 9 months ====
selected_month <- fiber_filt %>% filter(timepoint == "9m")

# 6.1 bartlett test homogeneity of variances ====
result_bartlett <- data.frame(fiber = character(),
                              K_squared_bartlett = numeric(),
                              df_bartlett = integer(),
                              p_value_bartlett = numeric(),
                              stringsAsFactors = FALSE)
for (col_name in fiber_columns) {
  result <- bartlett.test(selected_month[[col_name]], selected_month$diet_class)
  result_bartlett <- rbind(result_bartlett, data.frame(fiber = col_name,
                                                       K_squared_bartlett = result$statistic,
                                                       df_bartlett = result$parameter,
                                                       p_value_bartlett = result$p.value))
}
rm(result,col_name)

# 6.2 shapiro test normality tut ====
# Create an empty data frame to store results
result_shapiro <- data.frame(fiber = character(),
                             W_statistic_shapiro = numeric(),
                             p_value_shapiro = numeric(),
                             stringsAsFactors = FALSE)

# Perform Shapiro-Wilk test for each row
for (col_name in fiber_columns) {
  # Perform Shapiro-Wilk test
  result <- shapiro.test(selected_month[[col_name]])
  
  # Store results in dataframe
  result_shapiro <- rbind(result_shapiro, data.frame(fiber = col_name,
                                                     W_value = result$statistic,
                                                     p_value_shapiro = result$p.value))
}

barn_shapiro_semi <- full_join(result_shapiro,result_bartlett, by = "fiber")
barn_shapiro_semi$Age <- "9 months"
barn_shapiro <- rbind(barn_shapiro,barn_shapiro_semi)
rm(result,result_bartlett,result_shapiro,col_name,barn_shapiro_semi)

# 6.3 t-test ====
# CLASS 1 vs other
# Preallocate the data frame with the maximum possible number of rows
class1 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_1"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class1 <- rbind(class1, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class1$Adjusted_P_Value <- adjusted_p_values
class1$class <- "class1"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 2 vs other
# Preallocate the data frame with the maximum possible number of rows
class2 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_2"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class2 <- rbind(class2, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class2$Adjusted_P_Value <- adjusted_p_values
class2$class <- "class2"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 3 vs other
# Preallocate the data frame with the maximum possible number of rows
class3 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_3"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class3 <- rbind(class3, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class3$Adjusted_P_Value <- adjusted_p_values
class3$class <- "class3"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

#merging results together
result_t_test_semi <- rbind(class1,class2,class3)
result_t_test_semi$Age <- "9 months"
result_t_test <- rbind(result_t_test,result_t_test_semi)
rm(class1,class2,class3,selected_month,result_t_test_semi)

# 7.Calculate statistics and t-test between diet classes at 11 months ====
selected_month <- fiber_filt %>% filter(timepoint == "11m")

# 7.1 bartlett test homogeneity of variances ====
result_bartlett <- data.frame(fiber = character(),
                              K_squared_bartlett = numeric(),
                              df_bartlett = integer(),
                              p_value_bartlett = numeric(),
                              stringsAsFactors = FALSE)
for (col_name in fiber_columns) {
  result <- bartlett.test(selected_month[[col_name]], selected_month$diet_class)
  result_bartlett <- rbind(result_bartlett, data.frame(fiber = col_name,
                                                       K_squared_bartlett = result$statistic,
                                                       df_bartlett = result$parameter,
                                                       p_value_bartlett = result$p.value))
}
rm(result,col_name)

# 7.2 shapiro test normality tut ====
# Create an empty data frame to store results
result_shapiro <- data.frame(fiber = character(),
                             W_statistic_shapiro = numeric(),
                             p_value_shapiro = numeric(),
                             stringsAsFactors = FALSE)

# Perform Shapiro-Wilk test for each row
for (col_name in fiber_columns) {
  # Perform Shapiro-Wilk test
  result <- shapiro.test(selected_month[[col_name]])
  
  # Store results in dataframe
  result_shapiro <- rbind(result_shapiro, data.frame(fiber = col_name,
                                                     W_value = result$statistic,
                                                     p_value_shapiro = result$p.value))
}

barn_shapiro_semi <- full_join(result_shapiro,result_bartlett, by = "fiber")
barn_shapiro_semi$Age <- "11 months"
barn_shapiro <- rbind(barn_shapiro,barn_shapiro_semi)
rm(result,result_bartlett,result_shapiro,col_name,barn_shapiro_semi)

# 7.3 t-test ====
# CLASS 1 vs other
# Preallocate the data frame with the maximum possible number of rows
class1 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_1"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class1 <- rbind(class1, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class1$Adjusted_P_Value <- adjusted_p_values
class1$class <- "class1"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 2 vs other
# Preallocate the data frame with the maximum possible number of rows
class2 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_2"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class2 <- rbind(class2, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class2$Adjusted_P_Value <- adjusted_p_values
class2$class <- "class2"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 3 vs other
# Preallocate the data frame with the maximum possible number of rows
class3 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_3"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class3 <- rbind(class3, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class3$Adjusted_P_Value <- adjusted_p_values
class3$class <- "class3"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

#merging results together
result_t_test_semi <- rbind(class1,class2,class3)
result_t_test_semi$Age <- "11 months"
result_t_test <- rbind(result_t_test,result_t_test_semi)
rm(class1,class2,class3,selected_month,result_t_test_semi)

# 8.Calculate statistics and t-test between diet classes at 14 months ====
selected_month <- fiber_filt %>% filter(timepoint == "14m")

# 8.1 bartlett test homogeneity of variances ====
result_bartlett <- data.frame(fiber = character(),
                              K_squared_bartlett = numeric(),
                              df_bartlett = integer(),
                              p_value_bartlett = numeric(),
                              stringsAsFactors = FALSE)
for (col_name in fiber_columns) {
  result <- bartlett.test(selected_month[[col_name]], selected_month$diet_class)
  result_bartlett <- rbind(result_bartlett, data.frame(fiber = col_name,
                                                       K_squared_bartlett = result$statistic,
                                                       df_bartlett = result$parameter,
                                                       p_value_bartlett = result$p.value))
}
rm(result,col_name)

# 8.2 shapiro test normality tut ====
# Create an empty data frame to store results
result_shapiro <- data.frame(fiber = character(),
                             W_statistic_shapiro = numeric(),
                             p_value_shapiro = numeric(),
                             stringsAsFactors = FALSE)

# Perform Shapiro-Wilk test for each row
for (col_name in fiber_columns) {
  # Perform Shapiro-Wilk test
  result <- shapiro.test(selected_month[[col_name]])
  
  # Store results in dataframe
  result_shapiro <- rbind(result_shapiro, data.frame(fiber = col_name,
                                                     W_value = result$statistic,
                                                     p_value_shapiro = result$p.value))
}

barn_shapiro_semi <- full_join(result_shapiro,result_bartlett, by = "fiber")
barn_shapiro_semi$Age <- "14 months"
barn_shapiro <- rbind(barn_shapiro,barn_shapiro_semi)
rm(result,result_bartlett,result_shapiro,col_name,barn_shapiro_semi)

# 8.3 t-test ====
# CLASS 1 vs other
# Preallocate the data frame with the maximum possible number of rows
class1 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_1"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class1 <- rbind(class1, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class1$Adjusted_P_Value <- adjusted_p_values
class1$class <- "class1"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 2 vs other
# Preallocate the data frame with the maximum possible number of rows
class2 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_2"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class2 <- rbind(class2, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class2$Adjusted_P_Value <- adjusted_p_values
class2$class <- "class2"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

# CLASS 3 vs other
# Preallocate the data frame with the maximum possible number of rows
class3 <- data.frame(fibers = character(), 
                     P_Value = numeric(), 
                     T_Statistic = numeric(), 
                     stringsAsFactors = FALSE)

# Create an empty vector to store the p-values
p_values <- numeric()

# Loop through each column and perform t-test
for (col in fiber_columns) {
  # Create a formula for t-test
  formula <- as.formula(paste(col, "~ CLASS_3"))
  
  # Perform t-test
  result <- t.test(formula, data = selected_month, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value and t-statistic
  p_value <- result$p.value
  t_statistic <- result$statistic
  
  # Store p-value in the vector
  p_values <- c(p_values, p_value)
  
  # Store results in dataframe
  class3 <- rbind(class3, data.frame(fibers = col, P_Value = p_value, T_Statistic = t_statistic))
}

# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add adjusted p-values to the result dataframe
class3$Adjusted_P_Value <- adjusted_p_values
class3$class <- "class3"
rm(result,adjusted_p_values,col,formula,p_value,p_values,t_statistic)

#merging results together
result_t_test_semi <- rbind(class1,class2,class3)
result_t_test_semi$Age <- "14 months"
result_t_test <- rbind(result_t_test,result_t_test_semi)
rm(class1,class2,class3,selected_month,result_t_test_semi)

# 9. cleaning the combined df ====
result_t_test <- result_t_test %>% 
  mutate(sig_padj = if_else(
    condition = Adjusted_P_Value < 0.2 & Adjusted_P_Value >= 0.05, 
    true = "*", 
    false = if_else(Adjusted_P_Value < 0.05, "**", "")
  ))

#write.table (result_t_test, "S8_11_dietTP_vs_DF_t-test.txt", row.names=FALSE,sep = "\t")

result_t_test$P_Value <- NULL #removing p-values
result_t_test$Adjusted_P_Value <- NULL #removing p-values

#sig df
result_t_test_sig <- result_t_test %>% select(sig_padj, class,Age,fibers)
wide_sig <- pivot_wider(result_t_test_sig, names_from = fibers, values_from = sig_padj)
wide_sig$Age <- NULL
wide_sig$class <- NULL

#t-statistic df
result_t_test_stat <- result_t_test %>% select(T_Statistic, class,Age,fibers)
wide_stat <- pivot_wider(result_t_test_stat, names_from = fibers, values_from = T_Statistic)
wide_stat$class<- factor(wide_stat$class, 
                              levels = c("class1",
                                         "class2",
                                         "class3"
                              ))
wide_stat$Age<- factor(wide_stat$Age, 
                         levels = c("4 months",
                                    "5 months",
                                    "6 months",
                                    "9 months",
                                    "11 months",
                                    "14 months"
                         ))
wide_stat$Age <- gsub(" months", "", wide_stat$Age)
wide_stat <- wide_stat %>% mutate(rows = paste(class, Age, sep = "."))
wide_stat <- as.data.frame(wide_stat)
rownames(wide_stat) <- wide_stat$rows 
wide_stat$Age <- NULL
wide_stat$class <- NULL
wide_stat$rows <- NULL
wide_stat <- as.matrix(wide_stat)
column_labels <- gsub("_", " ", colnames(wide_stat))

rm(fiber_filt,result_t_test,result_t_test_stat,result_t_test_sig,custom_levels,fiber_columns,barn_shapiro)

# 10. heatmap, Panel 4A ====
my_palette <- colorRampPalette(c("#0072B2", "white", "#CD0000"))(100)

heatmap_with <- Heatmap(wide_stat, cluster_columns = T, cluster_rows = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(wide_sig[i,j], x, y, gp = gpar(fontsize = 6))
        },
        row_labels = gsub("c", "C", gsub("\\..*", "", rownames(wide_stat))),
        column_labels = column_labels,
        column_names_gp = gpar(fontsize = 6),
        column_title = NULL,
        row_names_gp = gpar(fontsize = 6),
        col = my_palette,
        row_split = c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3),rep(6,3)),
        row_title = NULL,
        row_names_side = "left",
        heatmap_legend_param = list(
          title = "T-statistic",
          at = c(-4, 0, 4),
          labels_gp = gpar(fontsize = 4),
          legend_height = unit(1, "cm"),
          legend_width = unit(0.5, "cm"),
          title_gp = gpar(fontsize = 6)),
        right_annotation = rowAnnotation(Age = anno_block(gp = gpar(fill = c("#BFD200","#AACC00", "#80B918",  "#55A630", "#2B9348" ,"#007F5F"), alpha = 0.7),
                                                          labels = c("4m","5m","6m","9m","11m","14m"),
                                                          labels_gp = gpar(col = "black", fontsize = 6),
                                                          show_name = T,
                                                          labels_rot = 0)))
heatmap_with 

# saving figure for the paper ====
# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 22.5
h_cm <- 17

tiff("S4.8_heatmap_without_antibiotics.tiff",
        width  = w_cm / 2.54,  # inches
        height = h_cm / 2.54,
        units  = "in",
        res    = 300)

draw(heatmap_with)  # important for ComplexHeatmap
dev.off()




