# Load packages ====
library(tidyverse)
library(phyloseq)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
lr = read.delim("Output/Files/S6.6_LR_pandemic_rich_scaled_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning meta ====
infant <- infant %>% ps_filter(older_siblings == "yes") 
meta <- as.data.frame(as.matrix(infant@sam_data)) 
meta <- meta %>% select(!(older_siblings))

meta$mat_weightgain <- as.numeric(meta$mat_weightgain)
meta$pregn_weeks <- as.numeric(meta$pregn_weeks)
meta$ENS <- as.numeric(meta$ENS)
meta$shannon_index <- as.numeric(meta$shannon_index)
meta$birthweight <- as.numeric(meta$birthweight)
meta$observed_richness <- as.numeric(meta$observed_richness)
meta$ageintrosolids <- as.numeric(meta$ageintrosolids)

#scaling data 
meta$birthweight <- scale(meta$birthweight)
meta$mat_weightgain <- scale(meta$mat_weightgain)
meta$ageintrosolids <- scale(meta$ageintrosolids)
meta$pregn_weeks <- scale(meta$pregn_weeks)

rm(infant)

final_lr <- lr %>% filter(last_run == "yes")
final_lr$Term <- sub("yes$", "", final_lr$Term)
final_lr$Term <- sub("male$", "", final_lr$Term)
final_lr$Term <- sub("vaginal$", "", final_lr$Term)
final_lr$Term[final_lr$Term == "delivery_placehospital"] <- "delivery_place"
final_lr <- final_lr %>% filter(Term != "older_siblings") %>% filter(Term != "(Intercept)")

rm(lr)

# LM for 1 week ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "1 week")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "1-2 weeks") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_1w <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_1w)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_1w)

# Results table (same structure as in your loop)
results_1w <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_1w$Age_individual <- "1 week"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_1w,meta_age)

# LM for 4 weeks ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "4 weeks")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "4 weeks") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_4w <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_4w)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_4w)

# Results table (same structure as in your loop)
results_4w <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_4w$Age_individual <- "4 weeks"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_4w,meta_age)

# LM for 8 weeks ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "8 weeks")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "8 weeks") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_8w <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_8w)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_8w)

# Results table (same structure as in your loop)
results_8w <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_8w$Age_individual <- "8 weeks"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_8w,meta_age)

# LM for 4 months ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "4 months")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "4 months") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_4m <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_4m)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_4m)

# Results table (same structure as in your loop)
results_4m <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_4m$Age_individual <- "4 months"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_4m,meta_age)

# LM for 5 months ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "5 months")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "5 months") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_5m <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_5m)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_5m)

# Results table (same structure as in your loop)
results_5m <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_5m$Age_individual <- "5 months"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_5m,meta_age)

# LM for 6 months ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "6 months")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "6 months") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_6m <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_6m)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_6m)

# Results table (same structure as in your loop)
results_6m <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_6m$Age_individual <- "6 months"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_6m,meta_age)

# LM for 9 months ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "9 months")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "9 months") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_9m <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_9m)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_9m)

# Results table (same structure as in your loop)
results_9m <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_9m$Age_individual <- "9 months"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_9m,meta_age)

# LM for 11 months ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "11 months")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "11 months") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_11m <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_11m)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_11m)

# Results table (same structure as in your loop)
results_11m <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_11m$Age_individual <- "11 months"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_11m,meta_age)

# LM for 14 months ====
# variables for specific age from LR for all samples
variables <- final_lr %>% filter(Age_individual == "14 months")
variables <- variables$Term

meta_age <- meta %>% filter(Age_individual == "14 months") 
meta_age <- meta_age %>% select (all_of(c("observed_richness", variables)))

# Create the formula dynamically
formula <- as.formula(paste("observed_richness ~", paste(variables, collapse = "+")))

# Fit the model
model_index_14m <- lm(formula, data = meta_age)

# Extract coefficients, p-values, and confidence intervals
model_summary <- summary(model_index_14m)
coefficients <- model_summary$coefficients
confidence_intervals <- confint(model_index_14m)

# Results table (same structure as in your loop)
results_14m <- data.frame(
  Term     = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error= coefficients[, "Std. Error"],
  t.value  = coefficients[, "t value"],
  p.value  = coefficients[, "Pr(>|t|)"],
  CI.Lower = confidence_intervals[, 1],
  CI.Upper = confidence_intervals[, 2],
  row.names = NULL
)

results_14m$Age_individual <- "14 months"

rm(variables,coefficients,model_summary,formula,confidence_intervals,model_index_14m,meta_age)

# merging together ====
all <- rbind(results_1w,results_4w,results_8w,results_4m,results_5m,results_6m,results_9m,results_11m,results_14m)
rm(results_1w,results_4w,results_8w,results_4m,results_5m,results_6m,results_9m,results_11m,results_14m)

rm(meta,final_lr)

# saving data ====
#write.table (all, "S6.25_LR_pandemic_rich_scaled_raw_with_siblings.txt", row.names=FALSE,sep = "\t")
