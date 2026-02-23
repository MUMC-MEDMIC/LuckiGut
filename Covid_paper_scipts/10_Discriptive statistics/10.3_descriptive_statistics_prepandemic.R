# load libraries ====
library(tidyverse)
library(phyloseq)
library(microViz)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
ps <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

#selecting family codes with covid samples 
ps_covid <- ps %>% ps_filter(covid_270220 == "yes")
covid_families <- unique(sort(ps_covid@sam_data$Family_ID))

ps_noncovid <- ps %>% ps_filter(!(Family_ID %in% covid_families)) 

# calculating number of samples and sets
meta <- as.data.frame(as.matrix(ps_noncovid@sam_data))
#subset for unique subject ID
unik = meta[!duplicated(meta$Family_ID),]
rm(ps_covid,ps,ps_noncovid,covid_families)

# selecting variables ====
milk <- c("breastfeeding_1_2w","breastfeeding_4w","breastfeeding_8w","breastfeeding_4m","breastfeeding_5m","breastfeeding_6m","breastfeeding_9m",
          "breastfeeding_11m","breastfeeding_14m")
formula <- c("formulafeeding_1_2w","formulafeeding_4w","formulafeeding_8w","formulafeeding_4m","formulafeeding_5m","formulafeeding_6m","formulafeeding_9m",
             "formulafeeding_11m","formulafeeding_14m")
anti <- c("inf_ab_1w","inf_ab_4w","inf_ab_8w","inf_ab_4m","inf_ab_5m","inf_ab_6m","inf_ab_9m","inf_ab_11m","inf_ab_14m")

other <- c("daycare_6m","daycare_14m","furrypet_home_yn_6m","furrypet_home_yn_14m")

common <- c("atopy_father","atopy_mother","covid_270220","delivery_place","delivery_type",
            "hospital","mat_AB_use_final_preg","mat_AB_use_final_delivery",
            "older_siblings","sex","smoke_before_preg","twins.siblings_in_lucki")

numeric <- c("ageintrosolids","birthweight","mat_weightgain","pregn_weeks","observed_richness","ENS")

all <- c(milk,anti,formula,other,common,numeric,"Age_individual", "Family_ID")

meta_perm <- unik %>% select(all_of(all))

rm(meta,milk,formula,other,common,all,unik)

# counting "antibiotic" questions ====
all_columns <- names(meta_perm)

antibio <- meta_perm %>% select(all_of(c(anti, "Age_individual")))
unique_answers <- unique(c(unlist(antibio[anti])))
unique_answers_noNA <- ifelse(is.na(unique_answers), "Missing", unique_answers)

# Initialize an empty data frame for the result
freq_table1 <- data.frame(matrix(0, nrow = length(unique_answers), ncol = length(anti)))
rownames(freq_table1) <- unique_answers_noNA
colnames(freq_table1) <- anti

# Fill the frequency table with counts
for (col in anti) {
  # Replace NA with the string "Missing" for counting
  meta_perm[[col]] <- ifelse(is.na(meta_perm[[col]]), "Missing", meta_perm[[col]])
  
  # Count the frequencies of answers
  counts <- plyr::count(meta_perm[[col]])
  # Populate the frequency table
  freq_table1[counts$x, col] <- counts$freq
}

remaining_columns <- setdiff(all_columns, anti)
check <- meta_perm %>% select(all_of(anti)) 
remaining_df <- meta_perm %>% select(-all_of(anti)) 

rm(col,anti,unique_answers,unique_answers_noNA,counts,check,antibio)

# counting "yes/no" questions ====
yes_columns <- remaining_columns[sapply(remaining_df, function(x) any(grepl("yes", x, ignore.case = TRUE), na.rm = TRUE))]

unique_answers <- unique(c(unlist(remaining_df[yes_columns])))
unique_answers_noNA <- ifelse(is.na(unique_answers), "Missing", unique_answers)

# Initialize an empty data frame for the result
freq_table2 <- data.frame(matrix(0, nrow = length(unique_answers), ncol = length(yes_columns)))
rownames(freq_table2) <- unique_answers_noNA
colnames(freq_table2) <- yes_columns

# Fill the frequency table with counts
for (col in yes_columns) {
  # Replace NA with the string "Missing" for counting
  remaining_df[[col]] <- ifelse(is.na(remaining_df[[col]]), "Missing", remaining_df[[col]])
  
  # Count the frequencies of answers
  counts <- plyr::count(remaining_df[[col]])
  # Populate the frequency table
  freq_table2[counts$x, col] <- counts$freq
}

remaining_columns <- setdiff(remaining_columns, yes_columns)
check <- remaining_df %>% select(all_of(yes_columns)) 
remaining_df <- remaining_df %>% select(-all_of(yes_columns)) 

rm(col,yes_columns,unique_answers,unique_answers_noNA,counts,check)

# counting "sex" questions ====
yes_columns <- remaining_columns[sapply(remaining_df, function(x) any(grepl("male", x, ignore.case = TRUE), na.rm = TRUE))]

unique_answers <- unique(c(unlist(remaining_df[yes_columns])))

# Initialize an empty data frame for the result
freq_table3 <- data.frame(matrix(0, nrow = length(unique_answers), ncol = length(yes_columns)))
rownames(freq_table3) <- unique_answers
colnames(freq_table3) <- yes_columns

# Fill the frequency table with counts
for (col in yes_columns) {
  # Replace NA with the string "Missing" for counting
  remaining_df[[col]] <- ifelse(is.na(remaining_df[[col]]), "Missing", remaining_df[[col]])
  
  # Count the frequencies of answers
  counts <- plyr::count(remaining_df[[col]])
  # Populate the frequency table
  freq_table3[counts$x, col] <- counts$freq
}

remaining_columns <- setdiff(remaining_columns, yes_columns)
check <- remaining_df %>% select(all_of(yes_columns)) 
remaining_df <- remaining_df %>% select(-all_of(yes_columns)) 

rm(col,yes_columns,unique_answers,counts,check)

# counting "delivery place" questions ====
yes_columns <- remaining_columns[sapply(remaining_df, function(x) any(grepl("home", x, ignore.case = TRUE), na.rm = TRUE))]

unique_answers <- unique(c(unlist(remaining_df[yes_columns])))
unique_answers_noNA <- ifelse(is.na(unique_answers), "Missing", unique_answers)

# Initialize an empty data frame for the result
freq_table4 <- data.frame(matrix(0, nrow = length(unique_answers), ncol = length(yes_columns)))
rownames(freq_table4) <- unique_answers_noNA
colnames(freq_table4) <- yes_columns

# Fill the frequency table with counts
for (col in yes_columns) {
  # Replace NA with the string "Missing" for counting
  remaining_df[[col]] <- ifelse(is.na(remaining_df[[col]]), "Missing", remaining_df[[col]])
  
  # Count the frequencies of answers
  counts <- plyr::count(remaining_df[[col]])
  # Populate the frequency table
  freq_table4[counts$x, col] <- counts$freq
}

remaining_columns <- setdiff(remaining_columns, yes_columns)
check <- remaining_df %>% select(all_of(yes_columns)) 
remaining_df <- remaining_df %>% select(-all_of(yes_columns)) 

rm(col,yes_columns,unique_answers,counts,check,unique_answers_noNA)

# counting "delivery type" questions ====
yes_columns <- remaining_columns[sapply(remaining_df, function(x) any(grepl("vaginal", x, ignore.case = TRUE), na.rm = TRUE))]

unique_answers <- unique(c(unlist(remaining_df[yes_columns])))
unique_answers_noNA <- ifelse(is.na(unique_answers), "Missing", unique_answers)

# Initialize an empty data frame for the result
freq_table5 <- data.frame(matrix(0, nrow = length(unique_answers), ncol = length(yes_columns)))
rownames(freq_table5) <- unique_answers_noNA
colnames(freq_table5) <- yes_columns

# Fill the frequency table with counts
for (col in yes_columns) {
  # Replace NA with the string "Missing" for counting
  remaining_df[[col]] <- ifelse(is.na(remaining_df[[col]]), "Missing", remaining_df[[col]])
  
  # Count the frequencies of answers
  counts <- plyr::count(remaining_df[[col]])
  # Populate the frequency table
  freq_table5[counts$x, col] <- counts$freq
}

remaining_columns <- setdiff(remaining_columns, yes_columns)
check <- remaining_df %>% select(all_of(yes_columns)) 
remaining_df <- remaining_df %>% select(-all_of(yes_columns)) 

rm(col,yes_columns,unique_answers,counts,check,unique_answers_noNA)

# counting numeric variables ====
remaining_df <- remaining_df %>% select(all_of(numeric))
remaining_df <- remaining_df %>% lapply(as.numeric)
#calculate shapiro test for all variables
shapiro <- remaining_df %>% lapply(function(x) shapiro.test(na.omit(x)))
result <- shapiro %>% sapply(`[`, c("statistic","p.value")) %>% t() %>% as.data.frame()
#idetify which one is a list
sapply(result, class)

result$statistic <- unlist(result$statistic)
result$p.value <- unlist(result$p.value)
#calculate standard deviation
sd <- remaining_df %>% lapply(function(x) sd(na.omit(x))) %>% as.data.frame() %>% t()
names(sd)[names(sd) == "V1"] <- "sd"
result <- cbind(result,sd)
#calculate variance
var <- remaining_df %>% lapply(function(x) var(na.omit(x))) %>% as.data.frame() %>% t()
names(var)[names(var) == "V1"] <- "variance"
result <- cbind(result,var)
#calculate the n
n <- remaining_df %>% lapply(function(x) length(na.omit(x))) %>% as.data.frame() %>% t()
names(n)[names(n) == "V1"] <- "n"
result <- cbind(result,n)
#calculate summary
num_variable <- remaining_df %>% lapply(function(x) summary(na.omit(x)))
num_variable <- num_variable %>% sapply(`[`, c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")) %>% t() %>% as.data.frame()
num_variable <- cbind(result,num_variable)
num_variable <- tibble::rownames_to_column(num_variable, "variable")

#calculate IQR
num_variable$IQR <- num_variable$`3rd Qu.` - num_variable$`1st Qu.`

rm(remaining_df,result,sd,shapiro,var,n,remaining_columns,numeric)

# calculating % ====
# calculating sum for 1 df
# Step 1: Calculate percentages for each column
df_percent1 <- freq_table1 / colSums(freq_table1) * 100
# Step 2: Add a row with the sum of the percentages (which should be 100)
sum_row <- colSums(df_percent1)
# Add the sum row to the data frame
df_percent1 <- rbind(df_percent1, sum_row)
#round up to 1 digit
df_percent1 <- round(df_percent1, 1)
rm(sum_row)

sum_row <- colSums(freq_table1)
freq_table1 <- rbind(freq_table1, sum_row)

# calculating sum for 2 df
df_percent2 <- freq_table2 / colSums(freq_table2) * 100
sum_row <- colSums(df_percent2)
df_percent2 <- rbind(df_percent2, sum_row)
df_percent2 <- round(df_percent2, 1)
rm(sum_row)
sum_row <- colSums(freq_table2)
freq_table2 <- rbind(freq_table2, sum_row)

# calculating sum for 3 df
df_percent3 <- freq_table3 / colSums(freq_table3) * 100
sum_row <- colSums(df_percent3)
df_percent3 <- rbind(df_percent3, sum_row)
df_percent3 <- round(df_percent3, 1)
rm(sum_row)
sum_row <- colSums(freq_table3)
freq_table3 <- rbind(freq_table3, sum_row)

# calculating sum for 4 df
df_percent4 <- freq_table4 / colSums(freq_table4) * 100
sum_row <- colSums(df_percent4)
df_percent4 <- rbind(df_percent4, sum_row)
df_percent4 <- round(df_percent4, 1)
rm(sum_row)
sum_row <- colSums(freq_table4)
freq_table4 <- rbind(freq_table4, sum_row)

# calculating sum for 5 df
df_percent5 <- freq_table5 / colSums(freq_table5) * 100
sum_row <- colSums(df_percent5)
df_percent5 <- rbind(df_percent5, sum_row)
df_percent5 <- round(df_percent5, 1)
rm(sum_row)
sum_row <- colSums(freq_table5)
freq_table5 <- rbind(freq_table5, sum_row)
rm(sum_row)

# combining freq and % ====
antibiotics <- freq_table1 %>% mutate(across(everything(), 
                                             ~ paste(.x, " (", df_percent1[[cur_column()]], "%)", sep = ""))) 
the_rest <- freq_table2 %>% mutate(across(everything(), 
                                          ~ paste(.x, " (", df_percent2[[cur_column()]], "%)", sep = "")))
sex <- freq_table3 %>% mutate(across(everything(), 
                                     ~ paste(.x, " (", df_percent3[[cur_column()]], "%)", sep = "")))
delivery_place <- freq_table4 %>% mutate(across(everything(), 
                                                ~ paste(.x, " (", df_percent4[[cur_column()]], "%)", sep = "")))
delivery_type <- freq_table5 %>% mutate(across(everything(), 
                                               ~ paste(.x, " (", df_percent5[[cur_column()]], "%)", sep = "")))
rm(freq_table1,freq_table2,freq_table3,freq_table4,freq_table5)
rm(df_percent1,df_percent2,df_percent3,df_percent4, df_percent5)

# writing results ====
#write.table (the_rest, "S10.2_descriptive_statistics_all.txt", row.names=TRUE,sep = "\t")
