#NB! at 9 and 11 months the collinearity is absolute for bioinfo batch and covid
# 0. Load packages ====
#statistics beta diversity,for infants PERMANOVA
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

# 1. Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
permanova = read.delim("Output/Files/S5.1_PERMANOVA_results_raw.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
permanova_run1 <- permanova %>% filter(run == 1)
rm(permanova)

# 2.1 Multivariable PERMANOVA for 1 week ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "1-2 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "1 week")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 2.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '1 week'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
week1_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 3.1 Multivariable PERMANOVA for 4 weeks ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "4 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "4 weeks")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 3.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '4 weeks'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
week4_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 4.1 Multivariable PERMANOVA for 8 weeks ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "8 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "8 weeks")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 4.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '8 weeks'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
week8_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 5.1 Multivariable PERMANOVA for 4 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "4 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "4 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 5.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '4 months'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
month4_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 6.1 Multivariable PERMANOVA for 5 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "5 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "5 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 6.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '5 months'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
month5_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 7.1 Multivariable PERMANOVA for 6 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "6 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "6 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 7.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '6 months'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
month6_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 10.1 Multivariable PERMANOVA for 14 months ====
dis_ait_I <- infant %>% ps_filter(Age_individual == "14 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age, shipmentinfo_batch is identical to covid_270220
needed_tp <- permanova_run1 %>% filter(age == "14 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "bioinfo_batch")
rm(needed_tp)

#all wetlab variables "bioinfo_batch","box_Latvia","plate_Latvia","shipmentinfo_batch","QC_metagenomics_Latvia"
# not working "box_Latvia","plate_Latvia","shipmentinfo_batch"

#Running a loop, 0.2 p-value stop
stop_condition <- FALSE
results <- list()

while(!stop_condition) {
  tryCatch({
    # Run permanova
    permanova <- dis_ait_I %>%
      dist_permanova(
        seed = 1,
        variables = variables,
        n_processes = 1, n_perms = 9999
      )
    
    # Save results
    permanova_semi <- permanova@permanova
    permanova_semi$param <- rownames(permanova_semi)
    permanova_semi <- permanova_semi %>% filter(!(param %in% c("Total", "Residual")))
    results[[length(results) + 1]] <- permanova_semi
    
    # Check if highest p-value is less than or equal to 0.2
    if (max(permanova_semi$`Pr(>F)`) <= 0.2) {
      stop_condition <- TRUE
    } else {
      # Remove variable with highest p-value
      highst_pvalue <- which.max(permanova_semi$`Pr(>F)`)
      variables <- variables[-highst_pvalue]
    }
  })
}

# 10.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '14 months'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
month14_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 12. merging results ====
permanova_all <- rbind(week1_permanova,week4_permanova,week8_permanova,month4_permanova,month5_permanova,
                       month6_permanova,month14_permanova)
rm(week1_permanova,week4_permanova,week8_permanova,month4_permanova,month5_permanova,
      month6_permanova,month14_permanova,infant)

# saving data ====
#write.table (permanova_all, "S5.4_tech_PERMANOVA_results_raw.txt", row.names=FALSE,sep = "\t")




