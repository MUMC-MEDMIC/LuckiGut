# Load packages ====
library(tidyverse)
library(phyloseq)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index = read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
permanova = read.delim("Output/Files/S5.1_PERMANOVA_results_raw.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
permanova_run1 <- permanova %>% filter(run == 1)
rm(permanova)

# cleaning meta and index ====
meta <- as.data.frame(as.matrix(infant@sam_data)) 
meta$mat_weightgain <- as.numeric(meta$mat_weightgain)
meta$pregn_weeks <- as.numeric(meta$pregn_weeks)
meta$ENS <- as.numeric(meta$ENS)
meta$shannon_index <- as.numeric(meta$shannon_index)
meta$birthweight <- as.numeric(meta$birthweight)
meta$observed_richness <- as.numeric(meta$observed_richness)
meta$ageintrosolids <- as.numeric(meta$ageintrosolids)

index$Family_ID <- rownames(index)
families <- index$Family_ID
meta <- meta %>% dplyr::filter(Family_ID %in% families)
meta_index <- full_join(index, meta, by = "Family_ID",multiple = "all")
meta_index$index_m_h <- as.numeric(meta_index$index_m_h)

meta_index <- meta_index %>% filter(covid_270220 == "yes") #only sample collected after covid started

rm(infant,families,index,meta)

#scaling data 
meta_index$birthweight <- scale(meta_index$birthweight)
meta_index$mat_weightgain <- scale(meta_index$mat_weightgain)
meta_index$ageintrosolids <- scale(meta_index$ageintrosolids)
meta_index$pregn_weeks <- scale(meta_index$pregn_weeks)
meta_index$index_m_h <- scale(meta_index$index_m_h)

# LM for 1 week ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "1 week")
variables <- sort(unique(needed_tp$variables))
variables <-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "1-2 weeks")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_1w <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_1w)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_1w)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_1w <- bind_rows(results, .id = "run")
results_1w$Age_individual <- "1 week"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 4 weeks ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "4 weeks")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "4 weeks")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_4w <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_4w)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_4w)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_4w <- bind_rows(results, .id = "run")
results_4w$Age_individual <- "4 weeks"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 8 weeks ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "8 weeks")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "8 weeks")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_8w <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_8w)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_8w)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_8w <- bind_rows(results, .id = "run")
results_8w$Age_individual <- "8 weeks"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 4 months ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "4 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "4 months")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_4m <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_4m)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_4m)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_4m <- bind_rows(results, .id = "run")
results_4m$Age_individual <- "4 months"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 5 months ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "5 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "5 months")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_5m <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_5m)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_5m)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_5m <- bind_rows(results, .id = "run")
results_5m$Age_individual <- "5 months"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 6 months ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "6 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
variables <- variables[variables != "delivery_place"] #this variable was not in final model in previous LR for alpha diversity
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "6 months")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_6m <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_6m)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_6m)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_6m <- bind_rows(results, .id = "run")
results_6m$Age_individual <- "6 months"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 9 months ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "9 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "9 months")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))
str(meta_age)

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_9m <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_9m)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_9m)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_9m <- bind_rows(results, .id = "run")
results_9m$Age_individual <- "9 months"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 11 months ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "11 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
variables <- variables[variables != "inf_ab_11m"] #removing antibiotics, as all "never"
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "11 months")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_11m <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_11m)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_11m)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_11m <- bind_rows(results, .id = "run")
results_11m$Age_individual <- "11 months"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# LM for 14 months ====
# variables for specific age from permanova
needed_tp <- permanova_run1 %>% filter(age == "14 months")
variables <- sort(unique(needed_tp$variables))
variables<-c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]
rm(needed_tp)

meta_age <- meta_index %>% filter(Age_individual == "14 months")
meta_age <- meta_age %>% select (all_of(c("ENS", variables)))

results <- list()

stop_condition <- FALSE

while (!stop_condition) {
  # Create the formula dynamically
  formula <- as.formula(paste("ENS ~", paste(variables, collapse = "+")))
  
  # Fit the model
  model_index_14m <- lm(formula, data = meta_age)
  
  # Get the summary of the model
  model_summary <- summary(model_index_14m)
  
  # Extract coefficients and p-values
  coefficients <- model_summary$coefficients
  
  # Calculate confidence intervals
  confidence_intervals <- confint(model_index_14m)
  
  # Create a data frame with coefficients, p-values, and confidence intervals
  p_values_df <- data.frame(
    Term = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std.Error = coefficients[, "Std. Error"],
    t.value = coefficients[, "t value"],
    p.value = coefficients[, "Pr(>|t|)"],
    CI.Lower = confidence_intervals[, 1],  # Lower bound of CI
    CI.Upper = confidence_intervals[, 2]   # Upper bound of CI
  )
  
  # Add results to the list
  results[[length(results) + 1]] <- p_values_df
  
  # Check if highest p-value is less than or equal to 0.2 (excluding intercept)
  if (max(p_values_df$p.value[-1]) <= 0.2) {
    stop_condition <- TRUE
  } else {
    # Find the variable with the highest p-value (excluding intercept)
    highest_pvalue_index <- which.max(p_values_df$p.value[-1])
    variable_to_remove <- variables[highest_pvalue_index]
    
    # Remove variable with highest p-value
    variables <- variables[-highest_pvalue_index]
  }
  
  # Stop if all variables have been removed
  if (length(variables) == 0) {
    stop_condition <- TRUE
  }
}

# Combine all results into a single dataframe
results_14m <- bind_rows(results, .id = "run")
results_14m$Age_individual <- "14 months"

rm(confidence_intervals,variables,meta_age,coefficients,model_summary,p_values_df,results,formula,highest_pvalue_index,stop_condition,variable_to_remove)

# merging together ====
all <- rbind(results_1w,results_4w,results_8w,results_4m,results_5m,results_6m,results_9m,results_11m,results_14m)
rm(results_1w,results_4w,results_8w,results_4m,results_5m,results_6m,results_9m,results_11m,results_14m)

# saving data ====
#write.table (all, "S6.9_LR_ENS_index_scaled_raw.txt", row.names=FALSE,sep = "\t")

# residuals ====
final_diagnostics_1w <- data.frame(
  fitted = fitted(model_index_1w),
  resid = residuals(model_index_1w),
  stud_resid = rstudent(model_index_1w),
  formula_used = paste(colnames(model_index_1w$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_1w$Age_individual <- "1 week"

final_diagnostics_4w <- data.frame(
  fitted = fitted(model_index_4w),
  resid = residuals(model_index_4w),
  stud_resid = rstudent(model_index_4w),
  formula_used = paste(colnames(model_index_4w$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_4w$Age_individual <- "4 weeks"

final_diagnostics_8w <- data.frame(
  fitted = fitted(model_index_8w),
  resid = residuals(model_index_8w),
  stud_resid = rstudent(model_index_8w),
  formula_used = paste(colnames(model_index_8w$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_8w$Age_individual <- "8 weeks"

final_diagnostics_4m <- data.frame(
  fitted = fitted(model_index_4m),
  resid = residuals(model_index_4m),
  stud_resid = rstudent(model_index_4m),
  formula_used = paste(colnames(model_index_4m$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_4m$Age_individual <- "4 months"

final_diagnostics_5m <- data.frame(
  fitted = fitted(model_index_5m),
  resid = residuals(model_index_5m),
  stud_resid = rstudent(model_index_5m),
  formula_used = paste(colnames(model_index_5m$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_5m$Age_individual <- "5 months"

final_diagnostics_6m <- data.frame(
  fitted = fitted(model_index_6m),
  resid = residuals(model_index_6m),
  stud_resid = rstudent(model_index_6m),
  formula_used = paste(colnames(model_index_6m$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_6m$Age_individual <- "6 months"

final_diagnostics_9m <- data.frame(
  fitted = fitted(model_index_9m),
  resid = residuals(model_index_9m),
  stud_resid = rstudent(model_index_9m),
  formula_used = paste(colnames(model_index_9m$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_9m$Age_individual <- "9 months"

final_diagnostics_11m <- data.frame(
  fitted = fitted(model_index_11m),
  resid = residuals(model_index_11m),
  stud_resid = rstudent(model_index_11m),
  formula_used = paste(colnames(model_index_11m$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_11m$Age_individual <- "11 months"

final_diagnostics_14m <- data.frame(
  fitted = fitted(model_index_14m),
  resid = residuals(model_index_14m),
  stud_resid = rstudent(model_index_14m),
  formula_used = paste(colnames(model_index_14m$model)[-1], collapse = " ,"))  # Shows final variables kept
final_diagnostics_14m$Age_individual <- "14 months"

rm(model_index_1w,model_index_4w,model_index_8w,model_index_4m,model_index_5m,model_index_6m,
   model_index_9m, model_index_11m, model_index_14m)

final_diagnostics <- rbind(final_diagnostics_1w,final_diagnostics_4w,final_diagnostics_8w,
                           final_diagnostics_4m,final_diagnostics_5m,final_diagnostics_6m,
                           final_diagnostics_9m,final_diagnostics_11m,final_diagnostics_14m)

rm(final_diagnostics_1w,final_diagnostics_4w,final_diagnostics_8w,
   final_diagnostics_4m,final_diagnostics_5m,final_diagnostics_6m,
   final_diagnostics_9m,final_diagnostics_11m,final_diagnostics_14m)

# saving data ====
#write.table (final_diagnostics, "S6.9_LR_ENS_index_scaled_residuals_raw.txt", row.names=FALSE,sep = "\t")




