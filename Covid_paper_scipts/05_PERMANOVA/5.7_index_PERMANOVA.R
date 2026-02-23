# Load packages ====
library(phyloseq)
library(microViz)
library(tidyverse)
library(vegan)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index = read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
permanova = read.delim("Output/Files/S5.1_PERMANOVA_results_raw.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

#merging with index 
index$Family_ID <- as.numeric(rownames(index))
infant_index <- ps_join(infant,index, by = "Family_ID")

#taking only thos sampel that have index data
baku_ei <- subset_samples(infant_index, !is.na(index_m_h)) #taking samples with index data
baku_ei <- baku_ei %>% ps_filter(covid_270220 == "yes") #taking samples after the covid

permanova_run1 <- permanova %>% filter(run == 1)

rm(infant,infant_index, index,permanova)

# Multivariable PERMANOVA for 5 months ====
dis_ait_I <- baku_ei %>% ps_filter(Age_individual == "5 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "5 months")
variables <- sort(unique(needed_tp$variables))
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

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

# cleaning results ====
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
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i,needed_tp)

# Multivariable PERMANOVA for 6 months ====
dis_ait_I <- baku_ei %>% ps_filter(Age_individual == "6 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "6 months")
variables <- sort(unique(needed_tp$variables))
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
variables <- variables[variables != "delivery_place"] #due to colinearity

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

# cleaning results ====
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
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i,needed_tp)

# Multivariable PERMANOVA for 9 months ====
dis_ait_I <- baku_ei %>% ps_filter(Age_individual == "9 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "9 months")
variables <- sort(unique(needed_tp$variables))
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

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

# cleaning results ====
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
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i,needed_tp)

# Multivariable PERMANOVA for 11 months ====
dis_ait_I <- baku_ei %>% ps_filter(Age_individual == "11 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "11 months")
variables <- sort(unique(needed_tp$variables))
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
variables <- variables[variables != "inf_ab_11m"] # no antibiotics in thsi group

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

# cleaning results ====
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
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i,needed_tp)

# Multivariable PERMANOVA for 14 months ====
dis_ait_I <- baku_ei %>% ps_filter(Age_individual == "14 months") %>% tax_agg("Species") %>% dist_calc("aitchison")

# variables for specific age
needed_tp <- permanova_run1 %>% filter(age == "14 months")
variables <- sort(unique(needed_tp$variables))
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

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

# cleaning results ====
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
rm(dis_ait_I, permanova, variables,highst_pvalue,permanova_semi,results,stop_condition, i,needed_tp)

# merging results ====
permanova_all <- rbind(month5_permanova,month6_permanova,month9_permanova,month11_permanova,month14_permanova)
rm(month5_permanova,month6_permanova,month9_permanova,month11_permanova,month14_permanova)

# saving table ====
#write.table (permanova_all, "S5.7_PERMANOVA_results_index_raw.txt", row.names=FALSE,sep = "\t")
