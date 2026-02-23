# Selecting only those TP, where covid was sig in PERMANOVA: 5, 6, 9, 11 and 14 months
# LOAD LIBRARIES ====
library("LinDA")
library("phyloseq")
library("tidyverse")
library("microViz")
library("ComplexHeatmap")
library("patchwork")

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
final = read.delim("Output/Files/S5.2_PERMANOVA_final.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
covid_bakut = read.delim("Output/Files/S7.2_LINDA_pandemic_covid_subset.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data ====
index$Family_ID <- rownames(index)
index$Family_ID <- as.numeric(index$Family_ID)

#merging with index 
infant_index <- ps_join(infant,index, by = "Family_ID")

# selecting only samples collected after covid
infant_covid <- infant_index %>% ps_filter(covid_270220 == "yes")
infant_covid <- subset_samples(infant_covid, !is.na(index_m_h)) #taking samples with index data

#check
plyr::count(infant_covid@sam_data$covid_270220)
plyr::count(infant_covid@sam_data$index_m_h)
rm(infant,infant_index,index)

#selecting bakut 
#picking up the bacteria of interest, they were either sig or below 0.2 of p
bakut <- unique(covid_bakut$bacteria)

# LINDA for 5 months ====
variables <- final %>% subset(age == "5 months") %>% pull(wrok_var)
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

ind_l <- infant_covid@sam_data$Age_individual == '5 months'
otu_DAA <- data.frame(infant_covid@otu_table[, ind_l])
otu_DAA <- otu_DAA[rownames(otu_DAA) %in% bakut, ]
meta_DAA <- data.frame(infant_covid@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)

meta_DAA$delivery_type<- factor(meta_DAA$delivery_type, 
                                levels = c("vaginal",
                                           "caesarian section"
                                ))
str(meta_DAA)

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ ageintrosolids+atopy_mother+breastfeeding_5m+daycare_6m+delivery_type+
                   mat_weightgain+older_siblings+index_m_h', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_5m <- output %>%  mutate(Age_individual = "5 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# LINDA for 6 months ====
variables <- final %>% subset(age == "6 months") %>% pull(wrok_var)
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

ind_l <- infant_covid@sam_data$Age_individual == '6 months'
otu_DAA <- data.frame(infant_covid@otu_table[, ind_l])
otu_DAA <- otu_DAA[rownames(otu_DAA) %in% bakut, ]
meta_DAA <- data.frame(infant_covid@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)

meta_DAA$delivery_type<- factor(meta_DAA$delivery_type, 
                                levels = c("vaginal",
                                           "caesarian section"
                                ))
str(meta_DAA)

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ ageintrosolids+atopy_father+atopy_mother+breastfeeding_6m+daycare_6m+delivery_type+
                   mat_weightgain+twins.siblings_in_lucki+index_m_h', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_6m <- output %>%  mutate(Age_individual = "6 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# LINDA for 9 months ====
variables <- final %>% subset(age == "9 months") %>% pull(wrok_var)
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

ind_l <- infant_covid@sam_data$Age_individual == '9 months'
otu_DAA <- data.frame(infant_covid@otu_table[, ind_l])
otu_DAA <- otu_DAA[rownames(otu_DAA) %in% bakut, ]
meta_DAA <- data.frame(infant_covid@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)

meta_DAA$delivery_type<- factor(meta_DAA$delivery_type, 
                                levels = c("vaginal",
                                           "caesarian section"
                                ))
str(meta_DAA)

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_mother+breastfeeding_9m+daycare_14m+delivery_type+formulafeeding_9m+older_siblings+
                   smoke_before_preg+index_m_h', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_9m <- output %>%  mutate(Age_individual = "9 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# LINDA for 11 months ====
variables <- final %>% subset(age == "11 months") %>% pull(wrok_var)
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

ind_l <- infant_covid@sam_data$Age_individual == '11 months'
otu_DAA <- data.frame(infant_covid@otu_table[, ind_l])
otu_DAA <- otu_DAA[rownames(otu_DAA) %in% bakut, ]
meta_DAA <- data.frame(infant_covid@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)

meta_DAA$delivery_place<- factor(meta_DAA$delivery_place, 
                                levels = c("hospital",
                                           "at home"
                                ))

str(meta_DAA)
plyr::count(meta_DAA$twins.siblings_in_lucki)

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_mother+birthweight+breastfeeding_11m+daycare_14m+delivery_place+
                   mat_AB_use_final_preg+older_siblings+twins.siblings_in_lucki+index_m_h', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_11m <- output %>%  mutate(Age_individual = "11 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# LINDA for 14 months ====
variables <- final %>% subset(age == "14 months") %>% pull(wrok_var)
variables <- c(variables, "index_m_h")
variables <- variables[variables != "covid_270220"]

ind_l <- infant_covid@sam_data$Age_individual == '14 months'
otu_DAA <- data.frame(infant_covid@otu_table[, ind_l])
otu_DAA <- otu_DAA[rownames(otu_DAA) %in% bakut, ]
meta_DAA <- data.frame(infant_covid@sam_data)
meta_DAA <- meta_DAA[ind_l, ]
meta_DAA <- meta_DAA[, c(variables), drop = FALSE] #adding covid variable to final core model
meta_DAA <- meta_DAA %>% mutate_if(is.numeric, as.numeric) %>% mutate_if(negate(is.numeric), as.character)

str(meta_DAA)
plyr::count(meta_DAA$twins.siblings_in_lucki)

formula_daa <- paste(variables, collapse = "+")
formula_daa <- paste0("'", paste("~", paste(variables, collapse = "+")), "'")
formula_daa

linda.obj <- linda(otu_DAA, meta_DAA, formula = '~ atopy_father+daycare_14m+formulafeeding_14m+furrypet_home_yn_14m+older_siblings+
                   pregn_weeks+smoke_before_preg+twins.siblings_in_lucki+index_m_h', alpha = 0.2,
                   prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#picking up data frames form linda's objects for significant adjusted p values
output <- bind_rows(linda.obj$output, .id = "variable")
output <- rownames_to_column(output, var = "bacteria")
output_14m <- output %>%  mutate(Age_individual = "14 months")

rm(variables,ind_l,otu_DAA,meta_DAA,linda.obj,output,formula_daa)

# combining all results together ====
final <- rbind(output_5m,output_6m,output_9m,output_11m,output_14m)
rm(output_5m,output_6m,output_9m,output_11m,output_14m)
#cleaning data
final$bacteria <- sub("\\.\\.\\..*", "", final$bacteria)

final$Age_individual<- factor(final$Age_individual, 
                              levels = c("5 months",
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

final <- final %>% 
  mutate(sig_p = if_else(
    condition = .$pvalue <= 0.05, 
    true = "*", 
    false = ""))

final_sig_bakut_covid <- final %>% filter(padj <= 0.05) %>%  filter(variable == "covid_270220yes")
species_covidsig_21 <- sort(unique(final_sig_bakut_covid$bacteria))
species_all_216 <- sort(unique(final$bacteria))


# saving data ====
#write.table (final, "S7.3_LINDA_final_sig_model_index_selected_bakut.txt", row.names=FALSE,sep = "\t")
