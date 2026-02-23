# LOAD LIBRARIES ====
library(phyloseq)
library(tidyverse)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
index_var <- read.delim("Output/Files/S1.7_index_variables_selected&cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# recording variables to numeric ====
#the higher is risk to get covid the higher is the number
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "no contact", 0, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "never", 5, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "mostly not", 4, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "sometimes", 3, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "often", 2, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "always", 1, .)))

index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "mother", 1, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "father", 2, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "mother, father", 3, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "mother, other", 4, .)))

index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "no contact", 0, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "contact on 1-2 days/week", 1, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "contact on 3-4 days/week", 2, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "contact on 5 or more days/week", 3, .)))

index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "(almost) never", 8, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "1-2 times/day", 7, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "3-4 times/day", 6, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "5-10 times/day", 5, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "11-15 times/day", 4, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "16-20 times/day", 3, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "21-25 times/day", 2, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == ">25 times/day", 1, .)))

index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "I'm not working at the moment", 0, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "I work from home (100%)", 1, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "On occasion I'm at work at the office (1-2 days/week)", 2, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "Often I'm at work at the office (3-4 days/week)", 3, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "Always I'm at work at the office (5 days/week)", 4, .)))

index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "no", 0, .)))
index_var <- index_var %>%  mutate(across(everything(), ~ ifelse(. == "yes", 1, .)))

# principal components analysis ====
#removing fathers
index_var_no_f <- index_var %>% select(!starts_with("f_"))
#making data numeric 
index_var_no_f[] <- lapply(index_var_no_f, function(x) as.numeric(as.character(x)))

# Standardize the data (scale to mean=0 and sd=1)
normalized_df <- as.data.frame(scale(index_var_no_f))
# PCA removing NAs from the data
for_pca_nona <- normalized_df[complete.cases(normalized_df), ]
#performing PCA
data_pca_family <- prcomp(for_pca_nona)
summary(data_pca_family)

rm(index_var,index_var_no_f)

# collecting loadings and scores ====
loadings_pca <- as.data.frame(data_pca_family$rotation) #High loadings for a variable on a principal component suggest that the variable has a strong influence on that component.
scores_pca <- as.data.frame(data_pca_family$x) #Scores indicate the position of each sample in terms of the principal components.

# factor analysis ====
fa_fit <- factanal(for_pca_nona, factors = 5, rotation = "varimax", scores = "regression")

# collecting loadings and scores ====
loadings <- as.data.frame(unclass(fa_fit$loadings)) 
scores <- as.data.frame(fa_fit$scores) 

# saving data ====
#write.table (loadings_pca, "S3.4_PCA_index_factor_analysis_loadings_from_PCA.txt", row.names=TRUE,sep = "\t")
#write.table (scores_pca, "S3.4_PCA_index_factor_analysis_scores_from_PCA.txt", row.names=TRUE,sep = "\t")

#write.table (loadings, "S3.4_PCA_index_factor_analysis_loadings_from_FA.txt", row.names=TRUE,sep = "\t")
#write.table (scores, "S3.4_PCA_index_factor_analysis_scores_from_FA.txt", row.names=TRUE,sep = "\t")

# plotting PCA loadings ====
back_up_df <- loadings_pca
back_up_df$variables <- rownames(back_up_df)
#sorting only those variables that are above 0.1, or below -0.01
back_up_df_PC1 <- back_up_df[abs(back_up_df$PC1) > 0.1, ]
ggplot(data=back_up_df_PC1, aes(x = variables,y=PC1)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
list_PC1 <- row.names(back_up_df_PC1)

back_up_df_PC2 <- back_up_df[abs(back_up_df$PC2) > 0.1, ]
ggplot(data=back_up_df_PC2, aes(x = variables,y=PC2)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
list_PC2 <- row.names(back_up_df_PC2)

list_PC1_PC2 <-unique(c(list_PC1, list_PC2)) #list of variables, whose PC1 or PC2 is big enough
loadings_0.1 <- back_up_df[back_up_df$variables %in% list_PC1_PC2, ]
loadings_0.1_PCA <- loadings_0.1[, -which(names(loadings_0.1) %in% "variables")]

rm(list_PC1,list_PC1_PC2,list_PC2,normalized_df,back_up_df,back_up_df_PC1,back_up_df_PC2,loadings_0.1)

# plotting FA loadings ====
back_up_df <- loadings
back_up_df$variables <- rownames(back_up_df)
#sorting only those variables that are above 0.1, or below -0.01
back_up_df_PC1 <- back_up_df[abs(back_up_df$Factor1) > 0.1, ]
ggplot(data=back_up_df_PC1, aes(x = variables,y=Factor1)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
list_PC1 <- row.names(back_up_df_PC1)

back_up_df_PC2 <- back_up_df[abs(back_up_df$Factor2) > 0.1, ]
ggplot(data=back_up_df_PC2, aes(x = variables,y=Factor2)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
list_PC2 <- row.names(back_up_df_PC2)

list_PC1_PC2 <-unique(c(list_PC1, list_PC2)) #list of variables, whose PC1 or PC2 is big enough
loadings_0.1 <- back_up_df[back_up_df$variables %in% list_PC1_PC2, ]
loadings_0.1 <- loadings_0.1[, -which(names(loadings_0.1) %in% "variables")]

rm(list_PC1,list_PC1_PC2,list_PC2,normalized_df,back_up_df,back_up_df_PC1,back_up_df_PC2)

# PCA for PCA ====
plot(scores_pca[, 1], scores_pca[, 2], 
     xlab = "Principal Component 1", ylab = "Principal Component 2",
     main = "PCA of index variables")

# PCA for FA ====
plot(scores[, 1], scores[, 2], 
     xlab = "Principal Component 1", ylab = "Principal Component 2",
     main = "PCA of index variables")

# scree plot for PCA ====
# Create a data frame for plotting
explained_variance <- summary(data_pca_family)$importance[2,]
scree_data <- data.frame(
  Principal_Component = 1:length(explained_variance),
  Explained_Variance = explained_variance
)

# Create the scree plot using ggplot2
ggplot(scree_data, aes(x = Principal_Component, y = Explained_Variance)) +
  geom_bar(stat = "identity") +
  geom_line() +
  geom_point() +
  xlab("Principal Component") +
  ylab("Explained variance") +
  ggtitle("Scree plot for teh index")

rm(data_pca_family,for_pca_nona,loadings,loadings_0.1,scores,scree_data,explained_variance)

