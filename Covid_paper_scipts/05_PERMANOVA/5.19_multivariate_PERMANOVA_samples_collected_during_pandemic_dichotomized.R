# 0. Load packages ====
#applied in combination with any distance measure. Here we use the aitchison distance.
library(phyloseq)
library(microViz)
library(tidyverse)
library(vegan)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# 1. Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
time = read.delim("Output/Files/S10.10_time_spent_in_pandemic_dichotomized.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# selecting samples ====
infant_time <- ps_join(infant,time, by = "MMHP_SampleID")
infant_time <- subset_samples(infant_time, covid_270220 == "yes")

rm(infant,time)

# 2.1 Multivariable PERMANOVA for 1 week ====
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "1-2 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("atopy_father","atopy_mother","birthweight","breastfeeding_1_2w","days_group","delivery_place","delivery_type",
"formulafeeding_1_2w","furrypet_home_yn_6m","hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery",
'mat_weightgain',"older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "4 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("atopy_father","atopy_mother","birthweight","breastfeeding_4w","days_group",
             "delivery_place","delivery_type","formulafeeding_4w","furrypet_home_yn_6m",
             "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "8 weeks") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("atopy_father","atopy_mother","birthweight","breastfeeding_8w","days_group",
             "delivery_place","delivery_type","formulafeeding_8w","furrypet_home_yn_6m",
             "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "4 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("ageintrosolids","atopy_father","atopy_mother","birthweight","breastfeeding_4m","days_group",
             "daycare_6m","delivery_place","delivery_type","formulafeeding_4m","furrypet_home_yn_6m",
             "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "5 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("ageintrosolids","atopy_father","atopy_mother","birthweight","breastfeeding_5m","days_group",
             "daycare_6m","delivery_place","delivery_type","formulafeeding_5m","furrypet_home_yn_6m",
             "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "6 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("ageintrosolids","atopy_father","atopy_mother","birthweight","breastfeeding_6m","days_group",
             "daycare_6m","delivery_place","delivery_type","formulafeeding_6m","furrypet_home_yn_6m",
             "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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

# 8.1 Multivariable PERMANOVA for 9 months ====
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "9 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("ageintrosolids","atopy_father","atopy_mother","birthweight","breastfeeding_9m","days_group",
             "daycare_14m","delivery_place","delivery_type","formulafeeding_9m","furrypet_home_yn_14m",
             "hospital","inf_ab_9m","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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

# 8.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '9 months'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
month9_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 9.1 Multivariable PERMANOVA for 11 months ====
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "11 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
# inf_ab_11m exclude -> no variation

variables<-c("ageintrosolids","atopy_father","atopy_mother","birthweight","breastfeeding_11m","days_group",
             "daycare_14m","delivery_place","delivery_type","formulafeeding_11m","furrypet_home_yn_14m",
             "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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

# 9.2 cleaning results ====
#adding number of run 
for (i in seq_along(results)) {
  results[[i]]$run <- i
}
#merging together
permanova <- bind_rows(results, .id = "run")
#cleaning df
permanova$age <- '11 months'
permanova$variables <- rownames(permanova)
permanova$variables <- stringr::str_remove(permanova$variables, "\\.{2}.*")
permanova <- permanova %>% filter(!(variables %in% c("Total", "Residual")))
month11_permanova <- subset(permanova, select = -param)
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i)

# 10.1 Multivariable PERMANOVA for 14 months ====
dis_ait_I <- infant_time %>% ps_filter(Age_individual == "14 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
variables<-c("ageintrosolids","atopy_father","atopy_mother","birthweight","breastfeeding_14m","days_group",
             "daycare_14m","delivery_place","delivery_type","formulafeeding_14m","furrypet_home_yn_14m",
             "hospital","inf_ab_14m","mat_AB_use_final_preg","mat_AB_use_final_delivery","mat_weightgain",
             "older_siblings","pregn_weeks","sex","smoke_before_preg","twins.siblings_in_lucki")

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

# 11. merging results ====
permanova_all <- rbind(week1_permanova,week4_permanova,week8_permanova,month4_permanova,month5_permanova,
                       month6_permanova,month9_permanova,month11_permanova,month14_permanova)
rm(week1_permanova,week4_permanova,week8_permanova,month4_permanova,month5_permanova,
      month6_permanova,month9_permanova,month11_permanova,month14_permanova)

# saving data ====
#write.table (permanova_all, "S5.19_PERMANOVA_samples_time_spend_in_pandemic_continuous_results_dichotomized_raw.txt", row.names=FALSE,sep = "\t")
