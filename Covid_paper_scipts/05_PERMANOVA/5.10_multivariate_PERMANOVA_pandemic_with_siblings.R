# Load packages ====
#statistics beta diversity, for infants PERMANOVA
#To test whether the groups are different with respect to centroid and dispersion, a PERMANOVA statistical test will be performed. 
#For this a multivariate extension of ANOVA will be used, as there are many OTU that will be used in the test. The extension is based 
#on distances between samples. The test compares distances of samples within the same group to distances of samples from different 
#groups. If the distance between samples from the different groups is much larger than samples from the same group, we conclude that 
#the groups are not equal.In order to test the significance of the result, a permutation test is used. Thus all samples are randomly 
#mixed over the groups and the test is repeated many times. If the ratio (between group distance / within group distance) is much 
#larger for the original data than for the permutations, we conclude there is a statistical significant difference.The test can be 
#applied in combination with any distance measure. Here we use the aitchison distance.
library(phyloseq)
library(microViz)
library(tidyverse)
library(vegan)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
permanova = read.delim("Output/Files/S5.2_PERMANOVA_results_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data ====
infant <- infant %>% ps_filter(older_siblings == "yes") 

final <- permanova %>% filter(last_run == "yes")

rm(permanova)

# Multivariable PERMANOVA for 1 week ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "1-2 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "1 week") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '1 week'
result$variables <- rownames(result)
week1_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 4 weeks ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "4 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "4 weeks") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '4 weeks'
result$variables <- rownames(result)
week4_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 8 weeks ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "8 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "8 weeks") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '8 weeks'
result$variables <- rownames(result)
week8_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 4 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "4 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "4 months") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '4 months'
result$variables <- rownames(result)
month4_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 5 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "5 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "5 months") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '5 months'
result$variables <- rownames(result)
month5_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 6 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "6 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "6 months") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '6 months'
result$variables <- rownames(result)
month6_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 9 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "9 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "9 months") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '9 months'
result$variables <- rownames(result)
month9_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 11 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "11 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "11 months") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '11 months'
result$variables <- rownames(result)
month11_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# Multivariable PERMANOVA for 14 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "14 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables <- final %>% filter(age == "14 months") %>% filter(variables != "Older siblings")
variables <- variables$wrok_var

# Run permanova once (no backward elimination)
permanova <- dis_ait_I %>%
  dist_permanova(
    seed = 1,
    variables = variables,
    n_processes = 1,
    n_perms = 9999
  )

# Format results (same as your loop body)
permanova_semi <- permanova@permanova
permanova_semi$param <- rownames(permanova_semi)

permanova_semi <- permanova_semi %>%
  dplyr::filter(!(param %in% c("Total", "Residual")))

# Save results (keeps your list-based workflow)
result <- permanova_semi
result$age <- '14 months'
result$variables <- rownames(result)
month14_permanova <- subset(result, select = -param)

rm(dis_ait_I, permanova, variables,permanova_semi)

# merging results ====
permanova_all <- rbind(week1_permanova,week4_permanova,week8_permanova,month4_permanova,month5_permanova,
                       month6_permanova,month9_permanova,month11_permanova,month14_permanova)
rm(week1_permanova,week4_permanova,week8_permanova,month4_permanova,month5_permanova,
      month6_permanova,month9_permanova,month11_permanova,month14_permanova,infant)

# saving data ====
#write.table (permanova_all, "S5.10_PERMANOVA_results_with_siblings_raw.txt", row.names=FALSE,sep = "\t")
