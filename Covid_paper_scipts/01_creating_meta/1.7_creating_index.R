# LOAD LIBRARIES ====
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggplot2)
library(GGally) #function "ggcorr"
library(ggcorrplot)
library(phyloseq)
library(microViz)
library(pheatmap)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# PART I. read in data ====
corona_raw <- read.delim("Output/Files/S1.2_corona_cleaned_data.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# creating a list of corona data ====
corona_list <- c("cor_home_in",
                 "cor_home_in_adultnr",
                 "cor_home_in_childnr",
                 "cor_home_out",
                 "cor_home_out_adultnr",
                 "cor_home_out_childnr",
                 "corm_work",
                 "corm_work_adultnr",
                 "corm_work_childnr",
                 "corm_work_hospital",
                 "corm_work_nothospit",
                 "corm_work_daycare",
                 "corm_work_primschool",
                 "corm_work_secschool",
                 "corm_work_bus",
                 "corm_work_otherocc",
                 "corm_publspace_in",
                 "corm_pspace_in_adunr",
                 "corm_pspace_chinr",
                 "corm_pspace_in_dist",
                 "corm_pspace_in_cap",
                 "corm_publicspace_out",
                 "corm_pspace_out_adnr",
                 "corm_pspace_out_chnr",
                 "corm_pspace_out_dist",
                 "corm_pspace_out_cap",
                 "corm_othhome_in",
                 "corm_othhome_in_adnr",
                 "corm_othhome_in_chnr",
                 "corm_othhome_in_dist",
                 "corm_othhome_in_cap",
                 "corm_othhome_out",
                 "corm_othhom_out_adnr",
                 "corm_othhom_out_chnr",
                 "corm_othhom_out_dist",
                 "corm_othhome_out_cap",
                 "corm_transp",
                 "corm_transp_adnr",
                 "corm_trans_chnr",
                 "corm_transp_dist",
                 "corm_transp_cap",
                 "corm_hands_water",
                 "corm_hands_watersoap",
                 "corm_hands_desinf",
                 "corf_work1",
                 "corf_work_adultnr1",
                 "corf_work_childnr1",
                 "corf_work_hospital1",
                 "corf_work_nothospit1",
                 "corf_work_daycare1",
                 "corf_work_primschoo1",
                 "corf_work_secschool1",
                 "corf_work_bus1",
                 "corf_work_otherocc1",
                 "corf_publspace_in1",
                 "corf_pspace_in_adun1",
                 "corf_pspace_chinr1",
                 "corf_pspace_in_dist1",
                 "corf_pspace_in_cap1",
                 "corf_publicspace_ou1",
                 "corf_pspace_out_adn1",
                 "corf_pspace_out_chn1",
                 "corf_pspace_out_dis1",
                 "corf_pspace_out_cap1",
                 "corf_othhome_in1",
                 "corf_othhome_in_adn1",
                 "corf_othhome_in_chn1",
                 "corf_othhome_in_dis1",
                 "corf_othhome_in_cap1",
                 "corf_othhome_out1",
                 "corf_othhom_out_adn1",
                 "corf_othhom_out_chn1",
                 "corf_othhom_out_dis1",
                 "corf_othhome_out_ca1",
                 "corf_transp1",
                 "corf_transp_adnr1",
                 "corf_trans_chnr1",
                 "corf_transp_dist1",
                 "corf_transp_cap1",
                 "corf_hands_water1",
                 "corf_hands_watersoa1",
                 "corf_hands_desinf1",
                 "corf_work_cap1",
                 "corf_work_dist_self1",
                 "corf_work_screen1",
                 "corm_work_cap",
                 "corm_work_dist_self",
                 "corm_work_screen",
                 "cor_hom_in_dist_self",
                 "cor_hom_in_dist_cha",
                 "cor_hom_in_dist_chic",
                 "cor_homout_dist_self",
                 "cor_homeout_dist_cha",
                 "cor_homeout_dist_chc")

# subsetting metadata for only covid variables (n=94) ====
corona <- corona_raw %>% select(all_of(corona_list),Family_ID) 

row.names(corona) <- corona$Family_ID
corona <- corona[, !names(corona) %in% "Family_ID"]

rm(corona_raw,corona_list)

# removing variables with only one level (n=86) ====
#creating function that counts number of levels
generate_levels_table <- function(df) {
  variable_names <- names(df)
  num_levels <- sapply(df, function(x) {
    num_unique <- length(unique(na.omit(x)))
    num_na <- sum(is.na(x))
    return(num_unique + ifelse(num_na > 0, 1, 0))
  })
  levels_table <- data.frame(Variable = variable_names, Levels = num_levels)
  return(levels_table)
}

levels_table <- generate_levels_table(corona)
level1 <- levels_table %>% filter(Levels == 1)
level2 <- levels_table %>% filter(Levels == 2)
list_level1 <- rownames(level1)
list_level2 <- rownames(level2)

#creating a function to count a levels for each of the variable
generate_count_table <- function(df, list_level) {
  count_table <- list()
  for(variable in list_level) {
    count_table[[variable]] <- table(factor(df[[variable]], levels = unique(df[[variable]]), exclude = NULL))
  }
  return(count_table)
}
generate_count_table(corona, list_level2)
# picking manually the variables with NA as a level
level_NA <- c("corf_trans_chnr1","corf_work_bus1","corf_work_primschoo1","corf_work_daycare1")

corona_levels <- corona %>% select(-all_of(c(list_level1,level_NA))) 
rm(corona,level1,level2,list_level1,list_level2,levels_table,level_NA)

# renaming variables ==== 
corona_fixed_col_names <- gsub("^cor_|cor", "", names(corona_levels))
names(corona_levels) <- corona_fixed_col_names
corona_fixed_col_names <- gsub("1", "", names(corona_levels))
names(corona_levels) <- corona_fixed_col_names
rm(corona_fixed_col_names)

# removing numeric data from df, except number of child and adults at work for parents (n=64) ====
numeric_columns <- c("f_othhom_out_adn",	"f_othhom_out_chn",	"f_othhome_in_adn",	"f_othhome_in_chn",	"f_pspace_chinr",	
                     "f_pspace_in_adun",	"f_pspace_out_adn",	"f_pspace_out_chn",
                     "home_in_adultnr",	"home_in_childnr",	"home_out_adultnr",	"home_out_childnr",	"m_othhom_out_adnr",	
                     "m_othhom_out_chnr",	"m_othhome_in_adnr",	"m_othhome_in_chnr",	"m_pspace_chinr",	"m_pspace_in_adunr",
                     "m_pspace_out_adnr","m_pspace_out_chnr","f_transp_adnr","m_transp_adnr")

corona_no_numb <- corona_levels %>% select(!all_of(numeric_columns)) 
corona_no_numb <- corona_no_numb[, order(names(corona_no_numb))]
rm(corona_levels,numeric_columns)

# creating essential jobs (n=58) ====
#essential jobs MOTHER
m_essen_job <- vector("character", length = nrow(corona_no_numb))
columns_to_run <- c("m_work_hospital", "m_work_nothospit","m_work_secschool","m_work_otherocc")
# Loop through each row
for (i in 1:nrow(corona_no_numb)) {
  # Initialize m_essen_job value as NA
  m_essen_job[i] <- NA
  
  # Loop through each column
  for (col in columns_to_run) {
    # Check if the value in the column is not NA
    if (!is.na(corona_no_numb[i, col])) {
      # Check if the value is "1"
      if (corona_no_numb[i, col] == "yes") {
        m_essen_job[i] <- "yes"
        break  # Exit the loop once "this" is assigned
      } else {
        next  # Move to the next column if value is "no"
      }
    }
  }
  
  # If all columns are "no" or NA, assign "that"
  if (is.na(m_essen_job[i])) {
    m_essen_job[i] <- "no"
  }
}

# Add m_essen_job as a new column to the dataframe
corona_no_numb$m_essen_job <- m_essen_job
rm(m_essen_job,col,i)

#essential jobs FATHER
f_essen_job <- vector("character", length = nrow(corona_no_numb))
columns_to_run2 <- c("f_work_hospital", "f_work_nothospit","f_work_secschool","f_work_otherocc")
for (i in 1:nrow(corona_no_numb)) {
  if (all(is.na(corona_no_numb[i, ]))) {
    X[i] <- NA
  } else {
    f_essen_job[i] <- NA
    
    for (col in columns_to_run2) {
      if (!is.na(corona_no_numb[i, col])) {
        if (corona_no_numb[i, col] == "yes") {
          f_essen_job[i] <- "yes"
          break  # Exit the loop once "this" is assigned
        } else {
          next  # Move to the next column if value is 0
        }
      }
    }
    
    if (is.na(f_essen_job[i])) {
      f_essen_job[i] <- "no"
    }
  }
}
corona_no_numb$f_essen_job <- f_essen_job
rm(f_essen_job,col,i)

used_var <- c(columns_to_run,columns_to_run2)
corona_no_numb <- corona_no_numb %>% select(!all_of(used_var))
rm(used_var,columns_to_run,columns_to_run2)

# recording all the variables to number ====
families <- rownames(corona_no_numb)
recoded <- corona_no_numb
#the higher is risk to get covid the higher is the number
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "no contact", 0, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "never", 5, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "mostly not", 4, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "sometimes", 3, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "often", 2, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "always", 1, .)))

recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "mother", 1, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "father", 2, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "mother, father", 3, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "mother, other", 4, .)))

recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "no contact", 0, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "contact on 1-2 days/week", 1, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "contact on 3-4 days/week", 2, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "contact on 5 or more days/week", 3, .)))

recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "(almost) never", 8, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "1-2 times/day", 7, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "3-4 times/day", 6, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "5-10 times/day", 5, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "11-15 times/day", 4, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "16-20 times/day", 3, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "21-25 times/day", 2, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == ">25 times/day", 1, .)))

recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "I'm not working at the moment", 0, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "I work from home (100%)", 1, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "On occasion I'm at work at the office (1-2 days/week)", 2, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "Often I'm at work at the office (3-4 days/week)", 3, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "Always I'm at work at the office (5 days/week)", 4, .)))

recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "no", 0, .)))
recoded <- recoded %>%  mutate(across(everything(), ~ ifelse(. == "yes", 1, .)))

recoded <- recoded %>% lapply(as.numeric) %>% as.data.frame()
row.names(recoded) <- families

# checking correlation (n=50) ====
# Compute the correlation matrix
correlation_matrix <- cor
ggcorr(recoded, label = TRUE, size = 1,hjust = 0.85) #default c("pairwise", "pearson"))
corr_plot <- ggcorr(recoded)
#to remove the variables that have above 0.5 correlation, fathers, as they have most of th NAs
correlations <- as.data.frame(corr_plot$data)
above_0.8 <- correlations %>% filter(label > 0.8) # no such numeric variables

ggplot(data = above_0.8, aes(x=x, y=y, fill=coefficient)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

variables_above_0.8 <- c("m_transp_cap","m_pspace_out_cap","f_transp_cap","f_othhome_out_ca","f_transp_dist",
                         "m_othhome_out_cap","m_transp_dist","m_work_screen")

recoded_low_correlation <- recoded %>% select(-all_of(variables_above_0.8))
recoded_low_correlation <- recoded[, !(names(recoded) %in% variables_above_0.8)]
rm(recoded,above_0.8,corr_plot,correlation_matrix,variables_above_0.8,correlations)

# creating new variables (n=50) na_no_numb,recoded_low_correlation,i,col,generate_count_table,generate_levels_table)

# combining old variables (n=48) into new ones (n=18) do not count Family_ID ====
#max value 18 -> exposure, 0 -> no exposure
#home
index_raw <- data.frame(matrix(ncol = 0, nrow = nrow(recoded_low_correlation))) #creating a new df

#checking what families have NAs
test <- recoded_low_correlation %>% select(!starts_with("f_"))
test <- test %>% filter(if_any(everything(), is.na)) %>% t() %>% as.data.frame() #only 2 families have NAs in none father questions
#only two questions have NAs -> will be treated as 0
rm(test)
#replacing NAs in non/father questions with 0 
recoded_low_correlation$hom_in_dist_chic[is.na(recoded_low_correlation$hom_in_dist_chic)] <- 0
recoded_low_correlation$m_pspace_in_dist[is.na(recoded_low_correlation$m_pspace_in_dist)] <- 0

index_raw$home_in_contacts <- rowSums(recoded_low_correlation[, c("home_in", "hom_in_dist_self", "hom_in_dist_cha","hom_in_dist_chic")], na.rm = FALSE) 
index_raw$home_out_contacts <- rowSums(recoded_low_correlation[, c("home_out", "homout_dist_self", "homeout_dist_cha","homeout_dist_chc")], na.rm = FALSE)
#public space
index_raw$m_pubspace_in_contacts <- rowSums(recoded_low_correlation[, c("m_publspace_in", "m_pspace_in_cap", "m_pspace_in_dist")], na.rm = FALSE) 
index_raw$f_pubspace_in_contacts <- rowSums(recoded_low_correlation[, c("f_publspace_in", "f_pspace_in_cap", "f_pspace_in_dist")], na.rm = FALSE)
index_raw$m_pubspace_out_contacts <- rowSums(recoded_low_correlation[, c("m_publicspace_out", "m_pspace_out_dist")], na.rm = FALSE) 
index_raw$f_pubspace_out_contacts <- rowSums(recoded_low_correlation[, c("f_publicspace_ou", "f_pspace_out_cap", "f_pspace_out_dis")], na.rm = FALSE) 

#other's home
index_raw$m_otherhome_in_contacts <- rowSums(recoded_low_correlation[, c("m_othhome_in", "m_othhome_in_cap", "m_othhome_in_dist")], na.rm = FALSE) 
index_raw$f_otherhome_in_contacts <- rowSums(recoded_low_correlation[, c("f_othhome_in", "f_othhome_in_cap", "f_othhome_in_dis")], na.rm = FALSE)
index_raw$m_otherhome_out_contacts <- rowSums(recoded_low_correlation[, c("m_othhome_out", "m_othhom_out_dist")], na.rm = FALSE) 
index_raw$f_otherhome_out_contacts <- rowSums(recoded_low_correlation[, c("f_othhome_out", "f_othhom_out_dis")], na.rm = FALSE) 
#transport, not used as only one variable with little variation in answered 
#index_raw$m_transport_contacts <- recoded_low_correlation$m_transp
#index_raw$f_transport_contacts <- recoded_low_correlation$f_transp
#hand washing
index_raw$f_low_hygiene <- rowSums(recoded_low_correlation[, c("f_hands_water", "f_hands_desinf","f_hands_watersoa")], na.rm = FALSE) 
index_raw$m_low_hygiene <- rowSums(recoded_low_correlation[, c("m_hands_water", "m_hands_desinf","m_hands_watersoap")], na.rm = FALSE) 
#protective gear at work
index_raw$f_protect_gear_work <- rowSums(recoded_low_correlation[, c("f_work_cap", "f_work_dist_self","f_work_screen")], na.rm = FALSE) 
index_raw$m_protect_gear_work <- rowSums(recoded_low_correlation[, c("m_work_cap", "m_work_dist_self")], na.rm = FALSE) 

index_raw$f_essen_job <- recoded_low_correlation$f_essen_job
index_raw$m_essen_job <- recoded_low_correlation$m_essen_job
index_raw$Family_ID <- rownames(recoded_low_correlation)

#max value 24 -> exposure, 0 -> no exposure
#mother
recoded_low_correlation <- recoded_low_correlation %>% mutate(m_work_adultnr_comb = ifelse(m_work_adultnr == 0, 0, 
                                                           ifelse(m_work_adultnr < 10, 5, 10)))
recoded_low_correlation <- recoded_low_correlation %>% mutate(m_work_childnr_comb = ifelse(m_work_childnr == 0, 0, 
                                                           ifelse(m_work_childnr < 10, 5, 10)))
index_raw$m_work_contact <- rowSums(recoded_low_correlation[, c("m_work", "m_work_adultnr_comb", "m_work_childnr_comb")], na.rm = FALSE)

#father
recoded_low_correlation <- recoded_low_correlation %>% mutate(f_work_adultnr_comb = ifelse(f_work_adultnr == 0, 0, 
                                                           ifelse(f_work_adultnr < 10, 5, 10)))
recoded_low_correlation <- recoded_low_correlation %>% mutate(f_work_childnr_comb = ifelse(f_work_childnr == 0, 0, 
                                                           ifelse(f_work_childnr < 10, 5, 10)))
index_raw$f_work_contact <- rowSums(recoded_low_correlation[, c("f_work", "f_work_adultnr_comb", "f_work_childnr_comb")], na.rm = FALSE)
rownames(index_raw) <- index_raw$Family_ID

# saving index ====
#write.table (index_raw, "S1.7_index_new_variables.txt", row.names=TRUE,sep = "\t")

# selecting variables that created index ====
var_sel_new_cleaned <- (c("home_in", "hom_in_dist_self", "hom_in_dist_cha","hom_in_dist_chic",
                      "home_out", "homout_dist_self", "homeout_dist_cha","homeout_dist_chc",
                      "m_publspace_in", "m_pspace_in_cap", "m_pspace_in_dist",
                      "f_publspace_in", "f_pspace_in_cap", "f_pspace_in_dist",
                      "m_publicspace_out", "m_pspace_out_dist",
                      "f_publicspace_ou", "f_pspace_out_cap", "f_pspace_out_dis",
                      "m_othhome_in", "m_othhome_in_cap", "m_othhome_in_dist",
                      "f_othhome_in", "f_othhome_in_cap", "f_othhome_in_dis",
                      "m_othhome_out", "m_othhom_out_dist",
                      "f_othhome_out", "f_othhom_out_dis",
                      "f_hands_water", "f_hands_desinf","f_hands_watersoa",
                      "m_hands_water", "m_hands_desinf","m_hands_watersoap",
                      "f_work_cap", "f_work_dist_self","f_work_screen",
                      "m_work_cap", "m_work_dist_self",
                      "f_essen_job",
                      "m_essen_job",
                      "m_work","m_work_adultnr_comb","m_work_childnr_comb", #combined columns
                      "f_work","f_work_adultnr_comb","f_work_childnr_comb"))

var_sel_old_cleaned <- (c("home_in", "hom_in_dist_self", "hom_in_dist_cha","hom_in_dist_chic",
                          "home_out", "homout_dist_self", "homeout_dist_cha","homeout_dist_chc",
                          "m_publspace_in", "m_pspace_in_cap", "m_pspace_in_dist",
                          "f_publspace_in", "f_pspace_in_cap", "f_pspace_in_dist",
                          "m_publicspace_out", "m_pspace_out_dist",
                          "f_publicspace_ou", "f_pspace_out_cap", "f_pspace_out_dis",
                          "m_othhome_in", "m_othhome_in_cap", "m_othhome_in_dist",
                          "f_othhome_in", "f_othhome_in_cap", "f_othhome_in_dis",
                          "m_othhome_out", "m_othhom_out_dist",
                          "f_othhome_out", "f_othhom_out_dis",
                          "f_hands_water", "f_hands_desinf","f_hands_watersoa",
                          "m_hands_water", "m_hands_desinf","m_hands_watersoap",
                          "f_work_cap", "f_work_dist_self","f_work_screen",
                          "m_work_cap", "m_work_dist_self",
                          "f_essen_job",
                          "m_essen_job",
                          "m_work","m_work_adultnr","m_work_childnr", #uncomb columns
                          "f_work","f_work_adultnr","f_work_childnr"))

corona_newvar_2 <- corona_no_numb %>% select(all_of(var_sel_old_cleaned))

# saving data ====
#write.table (corona_newvar_2, "S1.7_index_variables_selected&cleaned.txt", row.names=TRUE,sep = "\t")
rm(recoded_low_correlation,corona_newvar_2,var_sel_old_cleaned,var_sel_new_cleaned)

# translating to human language aka median ====
index_e <- index_raw
#map_exposure_max24 <- function(value) {
#  ifelse(is.na(value), NA,
#         ifelse(value == 0, "no exposure",
#                ifelse(value >= 1 & value <= 5, "little exposure",
#                       ifelse(value >= 6 & value <= 10, "moderate exposure",
#                              ifelse(value >= 11 & value <= 15, "high exposure",
#                                     ifelse(value >= 16 & value <= 24, "very high exposure", NA))))))
#}

# List of columns to process
list_of_columns_to_transform <- index_e %>% select(-f_essen_job, -m_essen_job) %>% colnames()
list_of_rows_to_transform <- index_e %>% select(-f_essen_job, -m_essen_job) %>% row.names()
df_to_transform <- index_e %>% select(-f_essen_job, -m_essen_job) 
list_of_columns_not_transform <- index_e %>% select(f_essen_job, m_essen_job) %>% colnames()
df_not_transform <- index_e %>% select(f_essen_job, m_essen_job)
#change all df to column list
column_list <- lapply(df_to_transform, function(x) list(x))
#remove NAs and 0 form the lists
column_list_no_NA_0 <- lapply(column_list, function(column) {
  column <- unlist(column)  # Convert column to atomic vector
  column[!is.na(column) & column != 0]
})
# Calculate median for each column in column_list_no_NA_0
medians <- sapply(column_list_no_NA_0, median)
rm(column_list_no_NA_0,column_list)
# Function to recode values based on median
recode_values <- function(column, median) {
  ifelse(is.na(column), NA, 
         ifelse(column == 0, 0, 
                ifelse(column < median, 1, 2)))
}
# Recode values in each column of df based on corresponding median
recode_df <- lapply(names(df_to_transform), function(col_name) {
  recode_values(df_to_transform[[col_name]], medians[col_name])
})
# Convert recoded list back to dataframe
recode_df <- as.data.frame(recode_df)

# Set column names of recoded dataframe
colnames(recode_df) <- list_of_columns_to_transform
row.names(recode_df) <- list_of_rows_to_transform
index_recoded <- cbind(recode_df, df_not_transform)

rm(df_not_transform,df_to_transform,index_raw,index_e,recode_df,list_of_columns_not_transform,list_of_columns_to_transform,
   list_of_rows_to_transform,medians)

# create final index f-m-h (n=18, fam = 32) ====
# Convert dataframe to long format for ggplot
df_long <- tidyr::pivot_longer(index_recoded, cols = everything())

# Plot histograms as a facet grid
ggplot(df_long, aes(x = value)) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(.~name, scales = "free") +
  labs(x = "Value", y = "Frequency", title = "Histograms of Columns")

#index_f_m_h
index_recoded$Family_ID <- NULL
index_recoded_no_NA <- na.omit(index_recoded)
index_f_m_h <- as.data.frame(rowSums(index_recoded_no_NA))
colnames(index_f_m_h)[colnames(index_f_m_h) == "rowSums(index_recoded_no_NA)"] <- "index_f_m_h"
hist(index_f_m_h$index_f_m_h,breaks = 10)
shapiro.test(index_f_m_h$index_f_m_h) # p-value = 0.7954 -> data is normally distributed
mean(index_f_m_h$index_f_m_h) #20.40625
sd(index_f_m_h$index_f_m_h) #5.802443
min(index_f_m_h$index_f_m_h) #9
max(index_f_m_h$index_f_m_h) #32

# saving fmh index ====
#write.table (index_f_m_h, "S1.7_index_fmh.txt", row.names=TRUE,sep = "\t")

rm(df_long,index_recoded_no_NA)

# create final index m-h (n=10, fam = 36) ====
index_sub <- index_recoded %>% select(-starts_with("f_"))
index_recoded_no_NA <- na.omit(index_sub)
index_m_h <- as.data.frame(rowSums(index_recoded_no_NA))
colnames(index_m_h)[colnames(index_m_h) == "rowSums(index_recoded_no_NA)"] <- "index_m_h"
hist(index_m_h$index_m_h,breaks = 10,
     main = "", 
     xlab = "Index values",
     cex.lab = 1,      # x and y axis label size
     cex.axis = 0.8)

shapiro.test(index_m_h$index_m_h) #p-value = 0.1706 -> data is normally distributed
mean(index_m_h$index_m_h) #11.36111
sd(index_m_h$index_m_h) #3.84078
min(index_m_h$index_m_h) #4
max(index_m_h$index_m_h) #18

# saving mh index ====
#write.table (index_m_h, "S1.7_index_mh.txt", row.names=TRUE,sep = "\t")

rm(index_recoded_no_NA,index_sub)

# additional info for the paper ====
summary(index_m_h$index_m_h)
shapiro.test(index_m_h$index_m_h)
