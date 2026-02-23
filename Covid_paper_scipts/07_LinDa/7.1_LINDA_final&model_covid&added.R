# 0. Load in libraries ====
#Differential abundance analysis (DAA)
# more here: https://rdrr.io/github/zhouhj1994/LinDA/man/linda.html

library("LinDA")
library("phyloseq")
library("tidyverse")
library("ComplexHeatmap")

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# 1. Import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
final = read.delim("Output/Files/S5.2_PERMANOVA_final.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# 2. final model added covid 1 week ====
variables <- final %>% subset(age == "1 week") %>% pull(wrok_var)
variables <- c(variables, "covid_270220")

ind_l <- infant@sam_data$Age_individual == '1-2 weeks'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
meta_DAA$sex   <- factor(meta_DAA$sex, 
                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_father+atopy_mother+breastfeeding_1_2w+delivery_type+formulafeeding_1_2w+
                   mat_weightgain+older_siblings+pregn_weeks+sex+covid_270220', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_1w <- output %>%  mutate(Age_individual = "1 week")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 3. final model added covid 4 weeks ====
variables <- final %>% subset(age == "4 weeks") %>% pull(wrok_var)
variables <- c(variables, "covid_270220")

ind_l <- infant@sam_data$Age_individual == '4 weeks'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_father+atopy_mother+breastfeeding_4w+delivery_type+formulafeeding_4w+furrypet_home_yn_6m+
                   mat_weightgain+older_siblings+twins.siblings_in_lucki+covid_270220', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_4w <- output %>%  mutate(Age_individual = "4 weeks")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 4. final model added covid 8 weeks ====
variables <- final %>% subset(age == "8 weeks") %>% pull(wrok_var)
variables <- c(variables, "covid_270220")

ind_l <- infant@sam_data$Age_individual == '8 weeks'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
                                  levels = c("hospital","at home"))
meta_DAA$sex   <- factor(meta_DAA$sex, 
                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_father+atopy_mother+birthweight+breastfeeding_8w+delivery_place+delivery_type+
                   formulafeeding_8w+furrypet_home_yn_6m+older_siblings+pregn_weeks+sex+twins.siblings_in_lucki+covid_270220', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_8w <- output %>%  mutate(Age_individual = "8 weeks")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 5. final model added covid 4 months ====
variables <- final %>% subset(age == "4 months") %>% pull(wrok_var)
variables <- c(variables, "covid_270220")

ind_l <- infant@sam_data$Age_individual == '4 months'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ ageintrosolids+atopy_mother+birthweight+breastfeeding_4m+daycare_6m+delivery_type+mat_weightgain+
                   older_siblings+pregn_weeks+smoke_before_preg+twins.siblings_in_lucki+covid_270220', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_4m <- output %>%  mutate(Age_individual = "4 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 6. final model had covid already 5 months ====
variables <- final %>% subset(age == "5 months") %>% pull(wrok_var)
ind_l <- infant@sam_data$Age_individual == '5 months'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, variables, drop = FALSE]
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ ageintrosolids+atopy_mother+breastfeeding_5m+covid_270220+daycare_6m+delivery_type+
                   mat_weightgain+older_siblings', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_5m <- output %>%  mutate(Age_individual = "5 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 7. final model had covid already 6 months ====
variables <- final %>% subset(age == "6 months") %>% pull(wrok_var)
ind_l <- infant@sam_data$Age_individual == '6 months'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, variables, drop = FALSE]
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ ageintrosolids+atopy_father+atopy_mother+breastfeeding_6m+covid_270220+daycare_6m+
                   delivery_type+mat_weightgain+twins.siblings_in_lucki', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_6m <- output %>%  mutate(Age_individual = "6 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 8. final model had covid already 9 months ====
variables <- final %>% subset(age == "9 months") %>% pull(wrok_var)
ind_l <- infant@sam_data$Age_individual == '9 months'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, variables, drop = FALSE]
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
                                 levels = c("vaginal","caesarian section"))
#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_mother+breastfeeding_9m+covid_270220+daycare_14m+delivery_type+
                   formulafeeding_9m+older_siblings+smoke_before_preg', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_9m <- output %>%  mutate(Age_individual = "9 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 9. final model had covid already 11 months ====
variables <- final %>% subset(age == "11 months") %>% pull(wrok_var)
ind_l <- infant@sam_data$Age_individual == '11 months'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, variables, drop = FALSE]
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
#meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
#                                 levels = c("vaginal","caesarian section"))
meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_mother+birthweight+breastfeeding_11m+covid_270220+daycare_14m+
                   delivery_place+mat_AB_use_final_preg+older_siblings+twins.siblings_in_lucki', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_11m <- output %>%  mutate(Age_individual = "11 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 10. final model had covid already 14 months ====
variables <- final %>% subset(age == "14 months") %>% pull(wrok_var)
ind_l <- infant@sam_data$Age_individual == '14 months'
otu_DAA <- data.frame(infant@otu_table[, ind_l])
meta_DAA <- data.frame(infant@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, variables, drop = FALSE]
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)
#meta_DAA$delivery_type <- factor(meta_DAA$delivery_type, 
#                                 levels = c("vaginal","caesarian section"))

#meta_DAA$delivery_place <- factor(meta_DAA$delivery_place, 
#                                  levels = c("hospital","at home"))
#meta_DAA$sex   <- factor(meta_DAA$sex, 
#                         levels = c("male", "female"))
str(meta_DAA)
variables

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_father+covid_270220+daycare_14m+formulafeeding_14m+furrypet_home_yn_14m+
                   older_siblings+pregn_weeks+smoke_before_preg+twins.siblings_in_lucki', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_14m <- output %>%  mutate(Age_individual = "14 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# 11. combining all results together ====
final <- rbind(output_1w,output_4w,output_8w,output_4m,output_5m,output_6m,output_9m,output_11m,output_14m)
rm(output_1w,output_4w,output_8w,output_4m,output_5m,output_6m,output_9m,output_11m,output_14m)
#cleaning data
final$bacteria <- sub("\\.\\.\\..*", "", final$bacteria)

final$Age_individual<- factor(final$Age_individual, 
                              levels = c("1 week",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"
                              ))
final <- final %>% 
  mutate(sig_padj = if_else(
    condition = .$padj > 0.2, 
    true = "", 
    false = if_else(
      condition = .$padj <= 0.2&.$padj >= 0.05,
      true = "*",
      false = "**")))

final_sig_bakut_covid <- final %>% filter(padj <= 0.05) %>%  filter(variable == "covid_270220yes")
species_covidsig_21 <- sort(unique(final_sig_bakut_covid$bacteria))
species_all_216 <- sort(unique(final$bacteria))

# saving data ====
#write.table (final, "S7.1_LINDA_pandemic.txt", row.names=FALSE,sep = "\t")
