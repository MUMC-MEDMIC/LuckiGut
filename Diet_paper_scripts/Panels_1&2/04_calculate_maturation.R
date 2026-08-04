library(randomForest)
library(tidyverse)
library(ggplot2)
library(vegan)
library(caret)





### STEPS (from paper PMID: )
# default parameters of randomForest regression: ntree=10000, mtry of p/3, where p is the number of input 97% ID OTUs
# regress the relative abundances of OTUs in the time series against chronologic age
# it identifies taxa that discriminate different periods of postnatal life

# input: rarefied OTU table at 2000 sequences per sample (here: just OTU table)
# 1. 100 iterations of RF algorithm to determine ranked lists of taxa in order ('feature importance')
# 2. rfcv function over 100 replicates to estimate the minimal number of top ranking age-discriminaroty taxa required for prediction (to choose a x-taxa model)
# -> training set: 12 healthy singletons (272 fecl samples)  (here: cross-validation on all samples)
# ->  nofurther parameter optimization
# -> output: sparse model of 24taxa (here: 44-taxa)
# 3. validation in other healthy children (13 singletons, 25 twins and triplets) (here in all samples from cross-validation)
# and applied to samples in children with SAM and MAM
# 4. smoothing spline function fit between microbiota age and chronologic age  in the validation set above 
# 5. Calculations of relative maturity and MAZ

# Read phyloseq object
phyloseq_obj <- readRDS(file = "generated_data/phyloseq_object.RDS")
  
# -> Make input data frame from OTU TABLE
  
tmp <- as.data.frame(t(phyloseq_obj@otu_table)) 
colnames(tmp) <- str_match(colnames(tmp), "s__(.*)")[,2] # only species names, skip if colnames correct
tmp <- tmp %>% #create a response variable 
  mutate(ChronoAge = as.integer(phyloseq_obj@sam_data$days_from_birth)) # columns with the date of the sample collection - chronological age


### 1. -> Split data into Train/Test set for model creation 
set.seed(100)
train <- sample(nrow(tmp), 0.66*nrow(tmp), replace = FALSE)
TrainSet <- tmp[train,]
TestSet <- tmp[-train,]

# RF algorithm  to determine ranked list of taxa ('feature importance')
rf.fit <- randomForest(ChronoAge ~ ., data = TrainSet, importance = TRUE, ntree = 1000, type = "regression")

# Plot and extract variables importance:
ImpData <- as.data.frame(importance(rf.fit)) %>%
  dplyr::arrange(desc(`%IncMSE`))
ImpData$Var.Names <- row.names(ImpData)

ImpData$Var.Names <- gsub("_"," ", ImpData$Var.Names)

# number of taxa in chosen model:
taxa_nr <- 44
# I ploted here 44 taxa (chosen model in 2.), you can choose how many top important taxa you want to plot 
ggplot(ImpData[1:taxa_nr,], aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  scale_x_discrete(limits = rev(ImpData[1:taxa_nr,]$Var.Names)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = "Species", y = "Increase in Mean Squared Error (MSE) (%)")


### 2. rfcv function over 100 replicates to estimate the minimal number of top ranking age-discriminaroty taxa required for prediction
# -> to later create the sparse model of minimal representation of top taxa 
rf.cv <- rfcv(trainx = TrainSet[,-ncol(TrainSet)], trainy = TrainSet[,ncol(TrainSet)], cv.fold = 100, step = 0.8)
error <- data.frame(n_spec = rf.cv$n.var, err = rf.cv$error.cv) 

# Save and plot the cross-validation error
ggplot(data = error, aes(x = err, y = n_spec)) +
  geom_point() +
  scale_y_reverse(breaks = c(326, 261, 209, 167, 134, 107, 85, 68, 55, 44, 35, 28, 22, 14, 4)) + # ! specify breaks 
  scale_x_continuous(position = "top") +
  geom_hline(yintercept=taxa_nr, linetype="dashed", color = "blue") + # choose your intercept to plot or skip (44 for chosen 44-taxa model)
  labs(x = "Cross-validation MSE",
       y = "Number of species") +
  theme_minimal() +
  theme(plot.background = element_rect(colour = "black", fill=NA, size=2))

ggsave("SP2_RF_taxa_importance_2.svg", width = 20, height = 20, unit = "cm")


#chosen model: 44 taxa

### 3. Before trainning/testing model, extract top x-taxa data 
sparse_model <- tmp[, colnames(tmp) %in% ImpData$Var.Names[1:taxa_nr]] %>% mutate(ChronoAge = as.integer(tmp$ChronoAge))

# -> train/test split with 100-k cross-validation
set.seed(100)
k_cv <- 100 # nr of cross-validations 
results <- list(length = 10)

for(cv in 1:k_cv) {
  # Create Train/Test split (66%/33%)
  folds <- caret::createFolds(phyloseq_big@sam_data$kindcode, k = 3)
  TrainSet <- sparse_model[-folds[[1]], ]
  TestSet <- sparse_model[folds[[1]], ]
  
  # Train model
  sparse_rf.fit <- randomForest(ChronoAge ~ ., data = TrainSet, ntree = 500, type = "regression")
  
  # Generate predictions (test model)
  sparse_rf.pred <- predict(sparse_rf.fit, TestSet[,-(taxa_nr+1)])
  
  # Save predictions
  tmp_pred <- rep(NA,nrow(tmp))
  tmp_pred[match(names(sparse_rf.pred), rownames(tmp))] <- sparse_rf.pred
  results[[cv]] <- tmp_pred
}

# -> Calculate mean of predictions for each sample:
pred_age <- rep(NA, nrow(tmp))
for(sample in 1:nrow(tmp)) {
  sample_vec <- rep(NA, k_cv)
  for(cv in 1:k_cv) {
    sample_vec[cv] <- results[[cv]][sample]
  }
  pred_age[sample] <- mean(sample_vec, na.rm = T)
}

# check if no NA in prediction results, if yes, increase the nr of cross validations of rerun it or drop samples:
sum(is.na(pred_age))



# calculate duration of breastfeeding
load("generated_data/timepoints_general_4m_to_14m.RData")
load("generated_data/clusters_solids.RData") 

# -> fill Na for breastfeeding at 4w 
clusters_solids$breastfeeding_4w[32] <- 0

# -> calculate
duration_dic <- c(4,5,6,9,11,14)
duration_breastfeeding <- vector("integer", length = nrow(phyloseq_obj@sam_data))
i <- 1
for(child in str_sub(phyloseq_obj@sam_data$kindcode, start = 1, end = 4)) {
  id <- match(child, clusters_solids$kindcode)
  duration <- 0
  # for tp 4m-14m:
  if(timepoints_general_4m_to_14m[[1]]$breastfeeding[id] == 0) {
    if(clusters_solids$breastfeeding_8w[id] == 1) {duration <- 2}
    else if(clusters_solids$breastfeeding_4w[id] == 1) {duration <- 1}
    else { duration <- 0 }
  }
  else {duration <- 0}
  # for tp 5-14m:
  for(t in 2:6) {
    if(timepoints_general_4m_to_14m[[t]]$breastfeeding[id] == 0) {break
    }
    else {duration <- duration_dic[t]}
  }
  duration_breastfeeding[i] <- duration
  i <- i+1
}



### 4. Visualize maturation prediction results and spline

# Create data frame with predictions and variables of interest 
microbial_age_df <- data.frame(kindcode = phyloseq_obj@sam_data$kindcode,
                               id = str_sub(phyloseq_obj@sam_data$kindcode, 1, 4),
                               diet_class = phyloseq_obj@sam_data$diet_class,
                               timepoint = factor(phyloseq_obj@sam_data$tp, levels = c("4m","5m","6m","9m","11m","14m")),
                               gender = phyloseq_obj@sam_data$bqq3_gender,
                               fur_pets = phyloseq_obj@sam_data$furrypet_indoors_6m,
                               siblings = as.integer(as.character(phyloseq_obj@sam_data$bqq45_siblings)),
                               birth_place = phyloseq_obj@sam_data$bqq13_delivery_place,
                               delivery_type = phyloseq_obj@sam_data$bqq14_delivery_type,
                               ageintrosolids = as.integer(as.character(phyloseq_obj@sam_data$ageintrosolids)),
                               birth_weigth = as.integer(as.character(phyloseq_obj@sam_data$bqq1_birthweight)),
                               duration_breastfeeding = duration_breastfeeding,
                               ChronoAge = as.integer(tmp$ChronoAge),
                               MicrobialAge = round(pred_age))
 
# Plot prediction results with loess regression line 
p1 <- ggplot(data = microbial_age_df, aes(x = ChronoAge, y = MicrobialAge)) +
  geom_point(aes(color = factor(diet_class))) + 
  geom_line(aes(x = ChronoAge, y = predict(loess(MicrobialAge~ ChronoAge, data = microbial_age_df))), color = "black", show.legend = F) +
  labs(x = "Chronological Age (days)",
       y = "Microbial Age (days)",
       title = bquote(bold("(A)")~"Microbial Age Predictions")) + 
  theme_bw() +
  guides(color=guide_legend(title="Dietary class")) +
  scale_x_continuous(limits = c(120, 440), breaks = c(124,155, 186,279, 341, 434)) + # adjust, if needed
  scale_y_continuous(limits = c(120, 440), breaks = c(124,155, 186,279, 341, 434)) +  # adjust, if needed
  geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed",   # adjust, if needed (for vertical lines where timepoints are)
             color = "grey", size=.5) +
  geom_label(data = data.frame(x = c(124,155, 186,279, 341, 434), y =  rep(430,6), label = c("4m","5m","6m","9m","11m","14m")),  # adjust, if needed (names of the timepoints on vertical lines )
            aes(x = x, y = y, label = label), color = "black", cex = 3) 

# Methods plot:
# ggplot(data = microbial_age_df, aes(x = ChronoAge, y = MicrobialAge)) +
#   geom_point(data = microbial_age_df %>% filter(ChronoAge == 214), color = "black") +
#   geom_line(aes(x = ChronoAge, y = predict(loess(MicrobialAge ~ ChronoAge, data = microbial_age_df))), color = "black", show.legend = F) +
#   geom_segment(x = 214, y = 165, xend = 214, yend = 207.4, linetype = "dashed", size = 0.15) +
#   geom_segment(x = 250, y = 165, xend = 214, yend = 165, arrow = arrow(length = unit(0.1, "cm"))) +  # arrow Ms
#   geom_label(aes(x = 250, y = 165,label = "Ms")) + 
#   geom_segment(x = 250, y = 210, xend = 214, yend = 210, arrow = arrow(length = unit(0.1, "cm"))) +  # arrow Mfit
#   geom_label(aes(x = 250, y = 210,label = "Mfit")) +  
#   geom_text(aes(x=350, y=225, label="Relative maturity = \nMs - Mfit"), colour="black", size = 5) +
#   labs(caption = "Ms: microbiota age of a sample\nMfit: microbiota age of interpolates spline fit")



## Calculate spline (baseline), median and sd (for MAZ calulcations ) 
microbial_age_df$spline <- predict(loess(MicrobialAge~ ChronoAge, data = microbial_age_df))
microbial_age_df <- microbial_age_df %>% 
  group_by(timepoint) %>% 
  mutate(med = median(MicrobialAge), 
         sd = sd(MicrobialAge),
         RelativeMaturity = MicrobialAge - spline) %>% # calculate Relative Maturity
  mutate(maz = (MicrobialAge - med) / sd) # calculate MAZ

# Save the microbial age dataframe 
write.csv(microbial_age_df, file = "generated_data/microbial_age_df.csv")

# -> Further visualize Relative Maturity and Maz, as boxplots for each timepoint

# Below with visualizations are inserted p values from statistical tests (pwc data) 
# It required to be adjusted for the x, xmin and xmax positions, as well as groups (I had problems with automatic command: add_xy_position())
require(ggpubr)
microbial_age_df <- microbial_age_df %>% mutate(timepoint = factor(timepoint, levels = c("4m","5m","6m","9m","11m","14m"), ordered = T))
pwc <- calcStats(microbial_age_df %>% dplyr::rename(month = timepoint), "RelativeMaturity", "diet_class")
pwc <- pwc %>% 
  mutate(y.position = get_y_position(formula = RelativeMaturity ~ diet_class, data = microbial_age_df %>% group_by(timepoint))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))
pwc$y.position[c(2 )] <- c(115)
p2 <- ggplot(data = microbial_age_df%>% mutate(diet_class = as.factor(diet_class)), 
         aes(x = timepoint, y = RelativeMaturity , color = diet_class )) +
  geom_boxplot(aes(color = diet_class), show.legend = F, outlier.shape=NA) + 
  geom_point(position=position_jitterdodge()) + 
  labs(x = "Time point",
       y = "Relative maturity",
       title = bquote(bold("(B)")~"Relative Microbiota Maturity")) +
  scale_color_discrete(name = "Dietary class") +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0)  
  

pwc <- calcStats(microbial_age_df %>% dplyr::rename(month = timepoint), "maz", "diet_class") # not plotted, as not stat difference 
p3 <- ggplot(data = microbial_age_df %>% mutate(diet_class = as.factor(diet_class)),
             aes(x =  timepoint, y = maz, color = diet_class)) +
  geom_boxplot(aes(color = diet_class), show.legend = F, outlier.shape = NA) + 
  geom_point(position=position_jitterdodge()) + 
  labs(x = "Time point",
       y = "MAZ",
       title = bquote(bold("(C)")~"MAZ")) +
  scale_color_discrete(name = "Dietary class") +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
  

# Generate final figure and save it
p1 + p2 + p3 + plot_layout(design = "
  0
  1
  2
", guides = "collect") & theme(legend.position = "bottom") 







calcStats <- function(data, variable, groups) {
  colnames(data)[which(colnames(data) == variable)] <- "variable"
  colnames(data)[which(colnames(data) == groups)] <- "groups"
  for(t in timepoints_names) {
    response <- (data %>% filter(month == t))[, "variable"]
    factor <- (data %>% filter(month == t))[, "groups"]
    if(bartlett.test(response, factor)$p.value >0.05) { # if >0.05 than groups have the same variance
      if(shapiro.test(response)$p.value >0.05) {#if >0.05 than  normal  (perform if bartlett.test >0.05
        tmp <- data %>% 
          filter(month == t) %>% 
          pairwise_t_test(variable ~ groups, p.adjust.method = "BH") %>%
          mutate(month = t, test = "ptt", statistic = "", .y. = variable) %>% select(-p.signif) 
        if(t == "4m") { # initialize pwc object
          pwc <- tmp
        } else {
          pwc <- rbind(pwc, tmp)
        }
        next
      }
    }
    tmp <- data %>% 
      filter(month == t) %>% 
      pairwise_wilcox_test(variable ~ groups, p.adjust.method = "BH", exact = F) %>%
      mutate(month = t, test = "pwt", .y. = variable) 
    
    if(t == "4m") { # initialize pwc object
      pwc <- tmp
    } else {
      pwc <- rbind(pwc, tmp)
    }
  }
  pwc <- pwc %>% mutate(month = factor(month, levels = timepoints_names))
  return(pwc)
}
