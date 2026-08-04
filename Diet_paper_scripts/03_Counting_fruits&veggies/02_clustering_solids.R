### LOAD LIBRARIES
require(logisticPCA) 
require(tidyverse) 
require(TSdist) # distance metrics
require(dtw) # dynamical time warping distance
require(superheat) # for heatmaps
require(clValid) 
require(reshape2) # for resaping data frames
require(factoextra) 
require(cluster) 
require(gridExtra) 
require(CGPfunctions)  # for 
require(pals)
require(dendextend) # for nicer dendograms
require(ggforce)
require(ggpubr)
require(cowplot) # for plotting multiple plots 
require(patchwork) # for combining plots
require(car) # for anova tests

### A. LOAD THE DIETARY DATA
load("generated_data/timepoints_general_4m_to_14m.RData") #list, with six elements (timepoints)
load("generated_data/other_variables.RData") #non-dietary variables 

### B. DATA PREPROCESSING:
# delete breast- and formulafeeding
# transform to matrices, for further steps 
for(i in 1:6) {
  timepoints_general_4m_to_14m[[i]] <- timepoints_general_4m_to_14m[[i]][,-c(1,2)] %>%
    mutate(across(everything(), ~as.numeric(as.character(.x))))
  timepoints_general_4m_to_14m[[i]] <- as.matrix(timepoints_general_4m_to_14m[[i]])
}
nr_children <- nrow(timepoints_general_4m_to_14m[[1]])
nr_fooditems <- ncol(timepoints_general_4m_to_14m[[1]])
nr_timepoints <- length(timepoints_general_4m_to_14m)


#################################### 
# 2. LOGISTIC PCA 
#################################### 
### A. RUN LOGISTIC PCA AND VISUALIZE IT (for 4m-14m)
timepoints_names <- c("4m", "5m", "6m", "9m", "11m", "14m")

# Important! You need to first define the value_s of parameter m for the logisticPCA

logpca_model <- list()
logpca_plot <- list()

for(i in 1:nr_timepoints) {
  logpca_model[[i]] <- logisticPCA(timepoints_general_4m_to_14m[[i]], k = 2, m = 3)
  colnames(logpca_model[[i]]$PCs) <- c("PC1", "PC2")
  logpca_model[[i]]$prop_deviance_expl <- round(logpca_model[[i]]$prop_deviance_expl, digits = 4)
   logpca_plot[[i]] <- plot(logpca_model[[i]], type = "scores", ) + 
    ggtitle(paste0(timepoints_names[i]," (",logpca_model[[i]]$prop_deviance_expl*100,"%)")) +
     scale_x_continuous(limits = c(-15,12)) + 
     scale_y_continuous(limits = c(-10,10)) + theme_minimal() + theme(legend.position = "none") 
    
}


logpca_plot[[1]] + logpca_plot[[2]] + logpca_plot[[3]] +logpca_plot[[4]] + logpca_plot[[5]] + logpca_plot[[6]] +
  plot_layout(design = "
  AB
  CD
  EF
", guides = "collect") 


### B. CREATE MATRICES FOR PC1, PC2 and one multivariate
timeseries_general_pc1 <- matrix(ncol = nr_timepoints, nrow =  nrow(logpca_model[[1]]$PCs))
rownames(timeseries_general_pc1) <- rownames(timepoints_general_4m_to_14m[[1]])
colnames(timeseries_general_pc1) <- timepoints_names
timeseries_general_pc2 <- matrix(ncol = nr_timepoints, nrow =  nrow(logpca_model[[1]]$PCs))
rownames(timeseries_general_pc2) <- rownames(timepoints_general_4m_to_14m[[1]])
colnames(timeseries_general_pc2) <- timepoints_names
for(i in 1:nr_timepoints) {
  timeseries_general_pc1[,i] <- logpca_model[[i]]$PCs[,1]
  timeseries_general_pc2[,i] <- logpca_model[[i]]$PCs[,2]
}

timeseries_general_m <- array(dim = c(nr_timepoints,2, nrow(logpca_model[[1]]$PCs))) #multidimensional
for(i in 1:nr_timepoints) {
  timeseries_general_m[i,1,] <- logpca_model[[i]]$PCs[,1]
  timeseries_general_m[i,2,] <- logpca_model[[i]]$PCs[,2]
}

### C.VISUALIZE LOADING PLOTS
#if ploted in the loop than stacked bar pplots (Fig P1C)
# heatmpas is 
tmp_plots <- list(length = 6)
for(i in 1:nr_timepoints) { 


  tmp <- as.data.frame(logpca_model[[i]]$U) %>% mutate(fooditem = colnames(timepoints_general_4m_to_14m[[1]]))

  ## A. arrange from positive to negative (if arrange like this, comment out B.)
  tmp <- tmp %>% arrange(desc(abs(V1)+abs(V2)))

  tmp_plots[[i]] <- ggplot(data = melt(tmp)) +
    geom_bar(aes(x = fooditem,
                 y = value, fill = variable),
             stat="identity", alpha = 0.8, position = "stack", width = 0.5) +
    theme_bw() +
    theme(axis.text.y = element_text(vjust = 0.5, hjust=1), axis.text.x = element_text(angle = 90)) +
    labs(y = "Loadings", x = "Food items", title = timepoints_names[i], fill = "") +
    scale_fill_manual(values = c("#f1c901", "#d46a92"), labels = c("PC1", "PC2")) +
    scale_y_continuous(breaks = c(-0.75,-0.5, -0.25,0, 0.25,0.5)) +
    coord_cartesian(ylim = c(-0.75,0.55)) +
    scale_x_discrete(limits = c("apple", "banana", "kiwi", "melon", "orange", "peach", "pear",  "strawberry", "tomato",
                                 "broccoli", "cauliflower", "carrots", "potato",
                                 "beans", "porr", "meat", "fish", "egg", "cheese","milk","yoghurt", "bread", "pasta", "rice",
                                 "soyprod", "pudding", "margarine", "butter"),
                      labels = c("apple", "banana", "kiwi", "melon", "orange", "peach", "pear",  "strawberry", "tomato",
                                 "broccoli", "cauliflower", "carrots", "potato",
                                 "beans", "porridge", "meat", "fish", "egg", "cheese","milk","yoghurt", "bread", "pasta", "rice",
                                 "soy products", "pudding", "margarine", "butter"))
  
  ##
  
  # B. arrange by food groups (if arrange like this, comment out A.)
  # 
  # tmp_plots[[i]] <- ggplot(data = melt(tmp), aes(x = fooditem,y = value, fill = variable)) +
  #   geom_bar(stat="identity", alpha = 0.8, position = "stack", width = 0.5) +
  #   theme_bw() +
  #   theme(axis.text.y = element_text(vjust = 0.5, hjust=1), axis.text.x = element_text(angle = 90)) +
  #   labs(y = "Loadings", x = "Food items", title = timepoints_names[i], fill = "") +
  #   scale_fill_manual(values = c("#f1c901", "#d46a92"), labels = c("PC1", "PC2")) +
  #   scale_y_continuous(breaks = c(-0.75, -0.5, -0.25,0, 0.25,0.5)) +
  #   coord_cartesian(ylim = c(-0.75,0.55)) +
  #   scale_x_discrete(limits = c("apple", "banana", "kiwi", "melon", "orange", "peach", "pear",  "strawberry", "tomato",
  #                               "broccoli", "cauliflower", "carrots", "potato",
  #                               "beans", "porr", "meat", "fish", "egg", "cheese","milk","yoghurt", "bread", "pasta", "rice",
  #                               "soyprod", "pudding", "margarine", "butter"),
  #                    labels = c("apple", "banana", "kiwi", "melon", "orange", "peach", "pear",  "strawberry", "tomato",
  #                               "broccoli", "cauliflower", "carrots", "potato",
  #                               "beans", "porridge", "meat", "fish", "egg", "cheese","milk","yoghurt", "bread", "pasta", "rice",
  #                               "soy products", "pudding", "margarine", "butter"))
}

tmp_plots[[1]] + tmp_plots[[2]] + tmp_plots[[3]] +
tmp_plots[[4]] + tmp_plots[[5]] + tmp_plots[[6]] + plot_layout(design = "
  AB
  CD
  EF
", guides = "collect")
 

# heatmap version 
# tmp <- tmp %>% mutate(timepoint = factor(timepoint, levels = timepoints_names, ordered =T),
#                       value = round(value, digits = 2))
# p1 <- coef_heatmap(df = tmp %>% filter(variable == "V1"), x = "fooditem", y = "timepoint", coef = "value", abs_coef_max = 0.5) + 
#   coord_flip() +
#   labs(y = "time point", title = "PC1", fill = "Loadings") +
#   theme_minimal() +
#   theme(axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank()
#   ) +
#   geom_text(aes(label = value))  
# 
# p2 <- coef_heatmap(df = tmp %>% filter(variable == "V2"), x = "fooditem", y = "timepoint", coef = "value", abs_coef_max = 0.5) + 
#   coord_flip() +
#   labs(y = "time point", title = "PC2", , fill = "Loadings") +
#   theme_minimal() +
#   theme(axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank()
#   ) +
#   geom_text(aes(label = value)) 
# 
# 
# 
# p1 + p2 + plot_layout(design = "
#   A
#   B
# ", guides = "collect")





### D. RESHAPE PC scores TO LONG FORMAT
df <- melt(timeseries_general_pc1)
df <- full_join(df, melt(timeseries_general_pc2), by=c("Var1","Var2"), suffix = c(".PC1", ".PC2")) 
df_reshaped <- df %>%
  dplyr::rename(child = Var1,
         month  = Var2,
         PC1 = value.PC1,
         PC2 = value.PC2) %>%
  mutate(month = factor(month),
         PC1 = round(PC1),
         PC2 = round(PC2),
         class = rep(0,nrow(df)))


#################################### 
# 3. HIERARCHICAL CLUSTERING 
#################################### 
### A. CALCULATE DISTANE MATRIX (multivaraite DTW)
distance <- matrix(ncol = nr_children, nrow = nr_children)
for(i in 1:nr_children) {
  for(j in 1:nr_children) {
    distance[i,j] <-  dtw( dist(timeseries_general_m[,,i],timeseries_general_m[,,j]), distance.only=T,window.type = "sakoechiba", window.size=2)$normalizedDistance 
  }
}
rownames(distance) <- rownames(timeseries_general_pc1)
colnames(distance) <- rownames(timeseries_general_pc1)

### B. PERFORM CLUSTERING & VISUALIZE DENDROGRAM (Fig P1A): 
clusters <- hclust(d = as.dist(distance), method = "complete")

# k: nr of clusters
plot(clusters %>% as.dendrogram %>%
       set("branches_k_color", k=3) %>% set("branches_lwd", 1.2) %>% set("labels_col", "white") %>% 
       set("leaves_pch", 19) %>% set("leaves_col", "black") %>% set("leaves_cex", .3),
     ylab = "Height", 
     main = "Clustering Dendrogram") 

### C. EXTRACT RESULTS
clusters_results_solids <- cutree(clusters, 3) # 3 clusters
df_reshaped$class <- rep(clusters_results_solids, nr_timepoints)


### D. ADDITIONAL VISUALIZATIONS 
# -> Time score plots colored by Dietary class (Fig P1B)
logpca_plot  <- list(length = nr_timepoints)
for(i in 1:nr_timepoints) {
  tmp_data <- as.data.frame(logpca_model[[i]]$PCs) %>% mutate(class = as.factor(clusters_results_solids))
  logpca_plot[[i]] <- ggplot(tmp_data, aes(PC1, PC2, color = class)) +
    geom_point(alpha = 0.7, position = "jitter") +
    labs(color = "Dietary class") + 
    ggtitle(paste0(timepoints_names[i]," (",logpca_model[[i]]$prop_deviance_expl*100,"%)")) +
    scale_x_continuous(limits = c(-15,12)) + 
    scale_y_continuous(limits = c(-10,10)) +
    theme_minimal()
}

logpca_plot[[1]] + logpca_plot[[2]] + logpca_plot[[3]] +logpca_plot[[4]] + logpca_plot[[5]] + logpca_plot[[6]] +
  plot_layout(design = "
  AB
  CD
  EF
", guides = "collect") 




# -> Animated one (for presentation)
####
# library(gganimate)
# library(gifski)
# tmp_data <- as.data.frame(logpca_model[[1]]$PCs) %>% mutate(time = 4, class = as.factor(clusters_results_solids))
# tmp <- c(4,5,6,9,11,14)
# for(i in 2:nr_timepoints) {
#   tmp_data <- rbind(tmp_data, (as.data.frame(logpca_model[[i]]$PCs) %>% mutate(time = tmp[i], class = as.factor(clusters_results_solids) )))
# }
# p <-ggplot(tmp_data, aes(PC1, PC2, color = class)) +
#   geom_point(alpha = 0.7, show.legend = FALSE) +
#   # Here comes the gganimate specific bits
#   labs(title = 'Time: {frame_time} m', x = 'PC1', y = 'PC2') +
#   transition_time(time) +
#   ease_aes('linear')
# animate(p, renderer = gifski_renderer())
# anim_save(filename = "thesis_plots/coloredPCAgif")
# ###








# # -> Time score plot, faceted by dietary class and colored for each timepoint 
# labels_.diet_class <- c("DIETARY CLASS 1", "DIETARY CLASS 2", "DIETARY CLASS 3")
# names(labels_.diet_class) <- c(1,2,3)
# 
# 
# ggplot(data = df_reshaped, aes(x=PC1, y=PC2, color = month)) +
#   geom_point(position = "jitter") +
#   facet_wrap(~class,  labeller = labeller(class = labels_.diet_class)) +
#   scale_color_manual("Time point", values = c("#EEE19E", "#4DAF4A",  "#00FFFF", "#0066FF" ,"#FF9900", "#450628")) + # color pallete or colors of timepoints 
#   geom_mark_ellipse(aes(color = month), expand = unit(0.5,"mm"), show.legend = FALSE) + 
#   theme_bw() +
#   theme(legend.position = "bottom",
#         strip.background = element_rect(
#           color="white", fill="white", size=1.5, linetype="solid"
#         ),
#         strip.text.x = element_text(face = "bold", size = 10)) +
#   geom_text(data = data.frame(label = c("n = 51", "n = 27", "n = 34"), # number of children in each cluster
#                               class   = c(1, 2, 3), 
#                               month = c("9m","9m","9m")),
#             mapping = aes(x=-10,y=-10, label = label), show.legend = FALSE) 




#################################### 
# 4. DESCRIBING / CHARACTERIZING THE CLUSTERS
#################################### 
### A. PRINT FOODITEMS COUNT PLOTS per timepoint, colored by cluster 
# Frist run function plotFooditemsCount (at the bottom of this script)
# Modify if needed, it's hard-coded for 3 clusters 
plotFooditemsCount(data = timepoints_general_4m_to_14m, groups = as.vector(clusters_results_solids))
# Modified for publication version below in the functions part: 

### B. CALCULATE THE TIME OF INTRODUCTION OF SOLID FOODS:
# set threshold (food considered introduces, when proportion of children consuming that food at earliest timestep >threshold 
threshold <- 0.5 
tmp_prop <- list(length = nr_timepoints)
# create list with fooditems proportion for each timepoint and cluster  (hardcoded for 3 clusters)
for(i in 1:nr_timepoints) { 
  indices <- as.vector(clusters_results_solids) # clustering h_mdtw
  tmp <- as.data.frame(timepoints_general_4m_to_14m[[i]][indices == 1, ]) 
  tmp_df <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 1),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(1,nr_fooditems))
  tmp <- as.data.frame(timepoints_general_4m_to_14m[[i]][indices == 2, ]) 
  tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 2),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(2,nr_fooditems))
  tmp_df <-rbind(tmp_df, tmp2)
  tmp <- as.data.frame(timepoints_general_4m_to_14m[[i]][indices == 3, ]) 
  tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 3),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(3,nr_fooditems))
  tmp_df <-rbind(tmp_df, tmp2)
  tmp_prop[[i]] <- tmp_df
}

# Calculate the time of introduction, from above proportions  
tmp_list <- list(length = 3)
for(gr in 1:3) { # calculate per cluster, when was each fooitem introduces (rule: introduecif prop>0.4)
  tmp_vector <- vector(mode = "character", length = nr_fooditems)
  for(fi in 1:nr_fooditems) { # per fooditem
    tmp_vector[fi] <- NA
    for(tp in 1:6) { # per timepoint 
      tmp2 <- tmp_prop[[tp]] %>% filter(fooditem == colnames(timepoints_general_4m_to_14m[[1]])[fi], class == gr)
      if (tmp2$prop >= threshold) {
        tmp_vector[fi] <- timepoints_names[tp]
        break
      }
    }
  }
  tmp_list[[gr]] <- tmp_vector
}
df_fooditems_introduction <- data.frame(fooditem = colnames(timepoints_general_4m_to_14m[[1]]), 
                                        cluster_1 = tmp_list[[1]], 
                                        cluster_2 = tmp_list[[2]],
                                        cluster_3 = tmp_list[[3]])

view(df_fooditems_introduction)
write.csv(df_fooditems_introduction, file = "06_timeofintro.csv")

### C. FOOD DIVERSITY ASSESMENT (Fig Supplementary):
# Food variety score (FVS): nr of singular food items (excl. breastfeeding / formulafeeding) max 28
# Dietary diversity score (DDS): nr of unique food groups (here: vegetables, fruits, legumes/nuts, meat, fish, eggs, diary, grains (bread, pasta, rice)) max 8 
# food allergen diversity (FAD): nr of main food allergens such as milk, egg, wheat, fish, soy, nuts max 6
fvs <- vector("integer")
dds <- vector("integer")
fad <- vector("integer")

# Calculate above metrics for each ifant at each timepoint 
for(tp in 1:nr_timepoints) {
  tmp_df <- as.data.frame(timepoints_general_4m_to_14m[[tp]]) 
  for(child in 1:nr_children) {
    fvs <- c(fvs, sum(tmp_df[child,]))
    tmp <-0
    if(1 %in% (tmp_df[child,] %>% select("cauliflower", "carrots", "broccoli", "potato"))) {tmp <- tmp +1} #vegetables
    if(1 %in% (tmp_df[child,] %>% select("banana", "pear", "apple", "melon", "peach", "kiwi", "orange", "strawberry", "tomato"))) {tmp <- tmp +1} #fruits
    if(1 %in% (tmp_df[child,] %>% select("beans", "porr"))) {tmp <- tmp +1} #legumes/nuts
    if(1 %in% (tmp_df[child,] %>% select("meat"))) {tmp <- tmp +1} #meat
    if(1 %in% (tmp_df[child,] %>% select("fish"))) {tmp <- tmp +1} #fish
    if(1 %in% (tmp_df[child,] %>% select("egg"))) {tmp <- tmp +1} #eggs
    if(1 %in% (tmp_df[child,] %>% select("milk", "yoghurt", "cheese"))) {tmp <- tmp +1} #diary
    if(1 %in% (tmp_df[child,] %>% select("bread", "rice", "pasta"))) {tmp <- tmp +1} #grains 
    dds <- c(dds, tmp)
    fad <- c(fad, sum(tmp_df[child,] %>% select("milk", "egg", "fish", "soyprod", "porr", "bread")))
  }
}

# Add metrics to reshaped data
df_reshaped <- df_reshaped %>% mutate(fvs = fvs, dds = dds, fad = fad)

# Visualize as boxplots with stats:
# FVS:
require(rstatix)
pwc <- calcStats(data = df_reshaped, variable = "fvs", groups = "class")
write.csv(pwc[,c(".y.","month","group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif")],
          file = "07_fvs_stats.csv")
pwc <- pwc %>% 
  mutate(y.position = get_y_position(formula = fvs ~ class, data = df_reshaped %>% group_by(month))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
        xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
              3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
        xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
              4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))

pwc$y.position[c(11, 12,13, 14,15, 18)] <- c(32, 30,27.76, 32, 30, 31)

p1 <- ggplot(data = df_reshaped, aes(x=month, y=fvs, fill=factor(class))) + 
  geom_boxplot() + ggtitle("Food variety score (FVS)") + 
  scale_y_continuous(limits = c(0,32.5), breaks = seq(0,28,4)) + 
  scale_fill_discrete(name="Dietary class") + 
  labs(x = "", y = "FVS", caption = "max = 28") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.caption = element_text(hjust = 0)) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = -2) 

# DDS:
pwc <- calcStats(data = df_reshaped, variable = "dds", groups = "class")
write.csv(pwc[,c(".y.","month","group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif")],
          file = "07_dds_stats.csv")
pwc <- pwc %>% 
  mutate(y.position = get_y_position(formula = dds ~ class, data = df_reshaped %>% group_by(month))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))

p2 <- ggplot(data = df_reshaped, aes(x=month, y=dds, fill=factor(class))) + 
  geom_boxplot() + ggtitle("Dietary diversity score (DDS)") +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8))  +
  scale_fill_discrete(name="Dietary class") + 
  labs(x = "", caption = "Groups: vegetables, fruits, legumes/nuts, meat, fish, eggs, diary, grains\nmax = 8",
       y = "DDS") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.caption = element_text(hjust = 0)) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0) 

# FAD
pwc <- calcStats(data = df_reshaped, variable = "fad", groups = "class")
write.csv(pwc[,c(".y.","month","group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif")],
          file = "07_fad_stats.csv")
pwc <- pwc %>% 
  mutate(y.position = get_y_position(formula = fad ~ class, data = df_reshaped %>% group_by(month))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))

p3 <- ggplot(data = df_reshaped, aes(x=month, y=fad, fill=factor(class))) + 
  geom_boxplot() + ggtitle("Food allergen diversity (FAD)") + 
  scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6)) + 
  scale_fill_discrete(name="Dietary class") + 
  labs(x = "Timepoint", y = "FAD", caption = "Allergens: milk, egg, wheat, fish, soy, nuts\nmax = 6") + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.caption = element_text(hjust = 0)) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0) 

combined <- p1 / p2 / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(combined)
# Additional significance stars added in graphical editor 


# Perform Statistical tests on it 
## -> differece between timepoints for each class
require(PMCMRplus)
for(col in c("fvs", "dds", "fad")) {
  print(col)
  for(cl in c(1,2,3)) {
    print(cl)
    response <- (df_reshaped %>% filter(class == cl))[,col]
    factor <- (df_reshaped %>% filter(class == cl))[,"month"]
    subject <- (df_reshaped %>% filter(class == cl))[,"child"]
    if(bartlett.test(response, factor)$p.value >0.05) { # if >0.05 than groups have the same variance
      if(shapiro.test(response)$p.value <0.05) {#if <0.05 than normal  (perform if bartlett.text >0.05
        print(aov(response ~ factor))
        print(pairwise.t.test(response, factor, paired = T, p.adjust.method = "BH"))
        next
      }
    } 
    print(friedman.test(response, factor, subject, p.adjust.method = "BH")) # perform if bartlett.test <0.05
    print(frdAllPairsConoverTest(response, factor, subject,p.adjust.method = "BH"))
  }
}





### D. LOOK AT THE EFFECT OF breast/formulafeeding 
# Caluclate proprotions of infant in a deitary cluster at each timepoint for breastfeeding and formulafeeding 
tmp_prop2 <- list(length = 6)
load("generated_data/timepoints.RData") # dietary_data with breastfeeding 
for(i in 1:6) { 
  indices <- as.vector(clusters_results_solids) # clustering h_mdtw
  tmp <- as.data.frame(timepoints[[i]][indices == 1, 1:2]) %>% mutate(breastfeeding = as.numeric(as.character(breastfeeding)),
                                                                      formulafeeding = as.numeric(as.character(formulafeeding)))
  tmp_df <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 1),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(1,2))
  
  tmp <- as.data.frame(timepoints[[i]][indices == 2, 1:2])  %>% mutate(breastfeeding = as.numeric(as.character(breastfeeding)),
                                                                       formulafeeding = as.numeric(as.character(formulafeeding)))
  tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 2),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(2,2))
  tmp_df <-rbind(tmp_df, tmp2)
  
  tmp <- as.data.frame(timepoints[[i]][indices == 3, 1:2]) %>% mutate(breastfeeding = as.numeric(as.character(breastfeeding)),
                                                                      formulafeeding = as.numeric(as.character(formulafeeding)))
  tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 3),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(3,2))
  tmp_df <-rbind(tmp_df, tmp2)
  tmp_prop2[[i]] <- tmp_df
}

# -> melt data into a data frame
df_fooditems_prop <- data.frame(rbind(tmp_prop2[[1]], tmp_prop2[[2]], tmp_prop2[[3]],
                                      tmp_prop2[[4]], tmp_prop2[[5]], tmp_prop2[[6]])) %>% 
  mutate(tp = c(rep(4, nr_timepoints), rep(5, nr_timepoints), rep(6, nr_timepoints), rep(9, nr_timepoints), rep(11, nr_timepoints), rep(14, nr_timepoints)))

# Perform linear regression with timepoint * dietary class interaction (only 6 samples, so no need for other models)
a <- lm(formula = prop ~ tp*class, data = df_fooditems_prop %>% filter(fooditem == "breastfeeding"))
a1 <- lm(formula = prop ~ tp*class, data = df_fooditems_prop %>% filter(fooditem == "formulafeeding"))
b <- data.frame(fooditem = c("breastfeeding","formulafeeding"),  # Create data frame to plot lines
                intercept = c(a$coefficients[[1]], a1$coefficients[[1]]),
                slope = c(a$coefficients[[4]],a$coefficients[[4]]))

# Plot scatter plots
ggplot(data = df_fooditems_prop, aes(x = tp, y = prop, color = factor(class))) +
  geom_point() + 
  facet_wrap(~fooditem) +
  geom_abline(data = b, aes(intercept = intercept, slope = slope), linetype = "dashed") + 
  geom_text(data = data.frame(fooditem = c("breastfeeding","formulafeeding"),
                              label = c("Anova, F = 1.284, p = 0.276", "Anova, F = 1.284, p = 0.276")),
            aes(x = 9, y = 0.01, label = label), color = "black", size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = c(4,5,6,9,11,14)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(2,"lines"),
        legend.position = "bottom",
        strip.background = element_rect(
          color="white", fill="white", size=1.5, linetype="solid"
        ),
        strip.text.x = element_text(face = "bold", size = 10)) +
  scale_color_discrete(name = "Dietary class") +
  labs(x = "Time", caption = "linear regression formula: prop ~ time*class") 



# -> Statistical tests performed on linear regression coeeficinets
Anova(a) # breastfeeding
Anova(a1) # formulafeeding




### E. INSPECT OTHER (NON-DIETARY) VARIABLES:
load("generated_data/other_variables.RData")

# First plot bar plots for categorical variables, colored by dietary class (Fig Supplementary)
col <- c(26:37, 12, 15:18) # nr of categorical columns in data
tmp_plot <- list(length = length(col))
for(i in 1:length(col)) {
  var <- col[i]
  cat(colnames(other_variables)[var])
  tmp <- as.matrix(table(other_variables$class, other_variables[,var])) /  as.vector(table(clusters_results_solids))
  print(tmp)
  #dev.new()
  tmp_plot[[i]] <- ggplot(data = as.data.frame(tmp), aes(x = Var2, y = Freq, fill = Var1)) +  
    geom_bar(position="dodge", stat="identity") +
    labs(x = "", y = "Proportion") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) + 
    ggtitle(colnames(other_variables)[var]) +
    scale_fill_discrete(name = "Dietary class")
}

# Combine plots together
tmp_plot[[1]] <- tmp_plot[[1]] + labs(title = "Breastfeeding at 4 months")
tmp_plot[[2]] <- tmp_plot[[2]] + labs(title = "Breastfeeding at 5 months", y ="")
tmp_plot[[3]] <- tmp_plot[[3]] + labs(title = "Breastfeeding at 6 months", y ="")
tmp_plot[[4]] <- tmp_plot[[4]] + labs(title = "Breastfeeding at 9 months", y ="")
tmp_plot[[5]] <- tmp_plot[[5]] + labs(title = "Breastfeeding at 11 months", y ="")
tmp_plot[[6]] <- tmp_plot[[6]] + labs(title = "Breastfeeding at 14 months", y ="")

tmp_plot[[7]] <- tmp_plot[[7]] + labs(title = "Formula feeding at 4 months")
tmp_plot[[8]] <- tmp_plot[[8]] + labs(title = "Formula feeding at 5 months", y ="")
tmp_plot[[9]] <- tmp_plot[[9]] + labs(title = "Formula feeding at 6 months", y ="")
tmp_plot[[10]] <- tmp_plot[[10]] + labs(title = "Formula feeding at 9 months", y ="")
tmp_plot[[11]] <- tmp_plot[[11]] + labs(title = "Formula feeding at 11 months", y ="")
tmp_plot[[12]] <- tmp_plot[[12]] + labs(title = "Formula feeding at 14 months", y ="")

tmp_plot[[13]] <-  tmp_plot[[13]] + scale_x_discrete(labels=c("At daycare","At home", "No")) + labs(title = "Furry pets")
tmp_plot[[14]] <-  tmp_plot[[14]] + scale_x_discrete(labels=c("Male","Female")) + labs(title = "Sex", y ="")
tmp_plot[[15]] <- tmp_plot[[15]] + scale_x_discrete(labels = c("Home", "Hospital")) + labs(title = "Delivery place", y ="")
tmp_plot[[16]] <- tmp_plot[[16]] + scale_x_discrete(labels = c("Vaginal", "C-section"))  + labs(title = "Delivery type", y ="")
tmp_plot[[17]] <-  tmp_plot[[17]] + labs(title = "Number of older siblings", y = "")


## Perform statistics on categorical variables, between dietary classes:
# -> Fisher's ExactTest on categorical variables
# might be error if some NA in data
for(i in 1:length(col)) {
  var <- col[i]
  # contigency table
  tmp_matrix = as.matrix(table(other_variables[,var], other_variables$class))
  test = chisq.test(tmp_matrix)
  print(colnames(other_variables)[var])
  print(test)
  if(test$p.value < 0.05) {
    #print(colnames(other_variables)[var])
  }
  tmp_plot[[i]] <- tmp_plot[[i]] + labs(caption = bquote(chi^2 == .(round(test$statistic, digits = 2)) * ", p =" ~ .(round(test$p.value, digits = 3))))
}

# Plot categorical variables
combined <- tmp_plot[[1]] + tmp_plot[[2]] + tmp_plot[[3]] + tmp_plot[[4]] + tmp_plot[[5]] + tmp_plot[[6]] + tmp_plot[[7]] + tmp_plot[[8]] + tmp_plot[[9]] + tmp_plot[[10]] + tmp_plot[[11]] + tmp_plot[[12]] & scale_x_discrete(labels=c("No","Yes"))
combined <- combined + tmp_plot[[13]] + tmp_plot[[14]] + tmp_plot[[15]] +tmp_plot[[16]] + tmp_plot[[17]] 

combined + plot_layout(design = "
  ABCDEF
  GHIJKL
  MNOPRR
", guides = "collect") & theme(legend.position = "bottom")

#ggsave("SP1_food_div_assesment.svg", width = 45, height = 20, units = "cm")
# '*Chi-square' added in graphical editor 

other_variables <- other_variables %>%
  mutate(ageintrosolids = as.numeric(as.character(ageintrosolids)),
         bqq1_birthweight = as.numeric(as.character(bqq1_birthweight)),
         bqq2_pregn_weeks = as.numeric(as.character(bqq2_pregn_weeks)),
         class = as.factor(class))

tmp_plot <- list(length = 3)
nr <- 1
for(i in c("ageintrosolids", "bqq1_birthweight", "bqq2_pregn_weeks")) {
  tmp <- other_variables %>% select("class", x = i)
  tmp_plot[[nr]] <- ggplot(data = tmp, aes(x = class, y = x, fill = as.factor(class))) + 
    geom_violin(width=1.4) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    labs(x = "Dietary class") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) + 
    theme(legend.position = "none") + 
    scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))
  nr <- nr + 1
}

pwc <- other_variables %>%
  pairwise_wilcox_test(ageintrosolids ~ class, exact = F) %>% 
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position()

tmp_plot[[1]] <- tmp_plot[[1]] + labs(y = "Age (weeks)", title = "Age solids introduced") + coord_cartesian(ylim = c(9,35)) + 
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T)
tmp_plot[[2]] <- tmp_plot[[2]] + labs(y = "Weight (g)", title = "Birth weight")
tmp_plot[[3]] <- tmp_plot[[3]] + labs(y = "Duration (weeks)", title = "Gestational age")


patchwork::wrap_plots(tmp_plot[[1]], tmp_plot[[2]], tmp_plot[[3]], 
          nrow = 1, ncol = 3)
ggsave("SP1_metadata2.svg", width = 25, height = 15, units = "cm")


# -> Perform statistical tests on numricial variables (ANOVA/KRUSKAL)
response <- as.numeric(as.character(other_variables$ageintrosolids)) #  or other numerical variable
factor <- other_variables$class
bartlett.test(response, factor)$p.value # if <0.05 than groups have different variance
shapiro.test(response)$p.value #if <0.05 than normal  (perform if bartlett.text >0.05)
aov(response ~ factor)
pairwise.t.test(response, factor)



kruskal.test(response ~ factor) # perform if bartlett.test <0.05
pairwise.wilcox.test(response, factor, exact=F) #performed for kruskal.test, post-hoc








### E. PLOT SLOP GRAPHS (Fig P1D)
# Create a data frame withproportion of infants in a dietary cluster consuming each food item
df_fooditems_prop <- list(length = 6)
for(i in 1:6) {
  df_fooditems_prop[[i]] <- data.frame(fooditem = colnames(timepoints_general_4m_to_14m[[1]]), 
                                       cluster_1 = tmp_prop[[i]] %>% filter(class == 1) %>% select(cluster_1 = prop), 
                                       cluster_2 = as.vector(tmp_prop[[i]] %>% filter(class == 2) %>% select(cluster_2 = prop)),
                                       cluster_3 = as.vector(tmp_prop[[i]] %>% filter(class == 3) %>% select(cluster_3 = prop)),
                                       month = timepoints_names[i])
}

df_fooditems_prop <- rbind(df_fooditems_prop[[1]], df_fooditems_prop[[2]], df_fooditems_prop[[3]],
                           df_fooditems_prop[[4]], df_fooditems_prop[[5]], df_fooditems_prop[[6]])
df_fooditems_prop <- df_fooditems_prop %>% mutate(month = factor(month, levels = timepoints_names))



## Slope graph per class (not used)
# choose class to plot as third argument (and change title)
newggslopegraph(df_fooditems_prop, month, cluster_3, fooditem, 
                LineThickness = 0.5) +
  labs(title="Proportion of fooditems in a dietary cluster", 
       subtitle="Dietary class 3") 

## Slope graph per fooditem (used)
# -> choose selection of fooditems to plot (or all)
selection <- c("apple", "banana", "beans", "broccoli", "butter",
               "carrots", "cauliflower", "kiwi", "margarine", "meat",
               "melon", "orange", "pasta", "pear", "porridge",
               "strawberry", "tomato", "yoghurt")
#selection <- c("bread", "cheese", "egg", "fish", "milk","peach","potato","pudding","rice","soyprod")
tmp <-  melt(df_fooditems_prop %>% filter(fooditem %in% selection)) 
tmp <- tmp %>%
  mutate(variable = factor(c(rep(1,nrow(tmp)/3), rep(2,nrow(tmp)/3), rep(3,nrow(tmp)/3))))

# -> Without selection
# tmp <- melt(df_fooditems_prop) %>% 
#   mutate(variable = factor(c(rep(1,168), rep(2,168), rep(3,168))))
#             

newggslopegraph(tmp, month, value, variable,
                Caption = NULL, SubTitle = NULL, Title = NULL)  +
  facet_wrap(~fooditem, nrow =2) +
  theme(
    strip.background = element_rect(
      color="black", fill="white", size=0.5, linetype="solid"
    ),
    strip.text.x = element_text(
      size = 10, color = "black", face = "bold"
    ),
    
  )

# Corrected in graphical editor
# porr -> porridge
# soyprod -> soyproducts


#################################### 
# 4. CLUSTERING OF DIET (WITHOUT BREAST/FORMULAFEEDING) : LOGISTIC PCA + H_CLUSTERING 
#################################### 
# A. LOGISTIC PCA (of fooditems)
logpca_model_diet <- list()
logpca_plot_diet <- list()
for(i in 1:6) {
  logpca_model_diet[[i]] <- logisticPCA(t(timepoints_general_4m_to_14m[[i]]), k = 2, m = 3)
  logpca_model_diet[[i]]$prop_deviance_expl <- round(logpca_model_diet[[i]]$prop_deviance_expl, digits = 4)
  colnames(logpca_model_diet[[i]]$PCs) <- c("PC1", "PC2")
  x<- as.data.frame(logpca_model_diet[[i]]$PCs)$PC1
  y <- as.data.frame(logpca_model_diet[[i]]$PCs)$PC2
  tmp <- data.frame(x = x, y = y)
  logpca_plot_diet[[i]] <- ggplot(data = tmp) + 
    geom_text(aes(x = x, y = y), label = colnames(timepoints_general_4m_to_14m[[i]])) +
    scale_x_continuous(limits = c(-30,30)) + 
    scale_y_continuous(limits = c(-22,20)) + theme_bw() + theme(legend.position = "none") +
  ggtitle(paste(timepoints_names[i]," (",logpca_model_diet[[i]]$prop_deviance_expl*100,"%)")) +
    labs(x = "PC1", y = "PC2")
}

# Plot andsavescore plots 
cowplot::plot_grid(logpca_plot_diet[[1]], logpca_plot_diet[[2]], logpca_plot_diet[[3]], logpca_plot_diet[[4]],logpca_plot_diet[[5]], logpca_plot_diet[[6]],
          nrow = 3, ncol = 2) 

ggsave(paste0("thesis_plots/11_diet_ordination.eps"), width = 20, height = 23, units = "cm")


# B. CREATE MULTIVARIATE MATRIX TO CALUCALTE DTW DISTANCE 
diet_m <- array(dim = c(6,2, nrow(logpca_model_diet[[1]]$PCs))) #multidimensional
for(i in 1:6) {
  diet_m[i,1,] <- logpca_model_diet[[i]]$PCs[,1]
  diet_m[i,2,] <- logpca_model_diet[[i]]$PCs[,2]
}

# 3. LONGITUDINAL CLUSTERING OF DIET (for all timepoints together)
# Calculate DTW distance
distance1 <- matrix(ncol = nr_fooditems, nrow = nr_fooditems)
for(i in 1:nr_fooditems) {
  for(j in 1:nr_fooditems) {
    distance1[i,j] <-  dtw( dist(diet_m[,,i],diet_m[,,j]), distance.only=T,window.type = "sakoechiba", window.size=2)$normalizedDistance 
  }
}
rownames(distance1) <- colnames(timepoints_general_4m_to_14m[[1]])
colnames(distance1) <- colnames(timepoints_general_4m_to_14m[[1]])

# hierarchical clustering, Ward's linkage
clusters <- hclust(d = as.dist(distance1), method = "ward.D")

# Plot and save dendrogram

plot(clusters %>% as.dendrogram %>%
      set("branches_lwd", 1.2) %>% set("labels_cex", 1) %>% 
       set("leaves_pch", 19) %>% set("leaves_col", "black") %>% set("leaves_cex", .5),
     ylab = "Height")




# D. CLUSTERING OF DIET per each timepoint separately:
# Calculate hamming distance

par(mfrow=c(2,3))
for(tp in 1:6) {
  distance1 <- matrix(ncol = nr_fooditems, nrow = nr_fooditems)
  for(i in 1:nr_fooditems) {
    for(j in 1:nr_fooditems) {
      distance1[i,j] <- sum(timepoints_general_4m_to_14m[[tp]][,i] != timepoints_general_4m_to_14m[[tp]][,j]) 
    }
  }
  rownames(distance1) <- colnames(timepoints_general_4m_to_14m[[1]])
  colnames(distance1) <- colnames(timepoints_general_4m_to_14m[[1]])
  clusters <- hclust(d = as.dist(distance1), method = "ward.D") # hierarchical clustering, Ward's linkage
  # Plot dendrogram 
  plot(clusters %>% as.dendrogram %>%
         set("branches_lwd", 1.2) %>% set("labels_cex", 1) %>% 
         set("leaves_pch", 19) %>% set("leaves_col", "black") %>% set("leaves_cex", .2),
       ylab = "Height", 
       main = timepoints_names[tp])
}
dev.off()


######################################
### 5. SAVE CLUSTERS SOLIDS
######################################

clusters_solids <- as.data.frame(other_variables) %>% 
  mutate(class = as.factor(clusters_results_solids))
save(clusters_solids, file = "generated_data/clusters_solids.RData")


# SAVE BREAST/FORMULA FEEDING vecors to use as sample data for children from phyloseq object 
# required: 03_metagenomics.R -> create phyloseq_big object
# required: 01_data.exploration.R -> timepoints list
# required: 02_clustering_solids.R -> clusters_solids data frame

# breastfeeding_sam <- ifelse(as.vector(sapply(phyloseq_big@sam_data$kindcode, function(x) {
#   kind <- as.numeric(str_sub(x, start = 1, end = 4))
#   tp <- as.numeric(str_sub(x, start = 8, end = 8)) - 3
#   timepoints[[tp]][which(clusters_solids$kindcode == kind), "breastfeeding"]
#   })) == 1,TRUE, FALSE)
# 
# formulafeeding_sam <- ifelse(as.vector(sapply(phyloseq_big@sam_data$kindcode, function(x) {
#   kind <- as.numeric(str_sub(x, start = 1, end = 4))
#   tp <- as.numeric(str_sub(x, start = 8, end = 8)) - 3
#   timepoints[[tp]][which(clusters_solids$kindcode == kind), "formulafeeding"]
# })) == 1, TRUE, FALSE)
# 
# saveRDS(breastfeeding_sam, file = "tmp_breastfeeding_sam.RDS")
# saveRDS(formulafeeding_sam, file = "tmp_formulafeeding_sam.RDS")





### FUNCTIONS
plotFooditemsCount <- function(data, groups) {
  stop <- length(data)
  tmp_plot <- list(length = stop)
  for(i in 1:stop) {
    indices <- groups # clustering 5/9
    tmp <- as.data.frame(data[[i]][indices == 1, ]) 
    tmp_df <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 1),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(1,ncol(data[[i]])))
    tmp <- as.data.frame(data[[i]][indices == 2, ]) 
    tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 2),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(2,ncol(data[[i]])))
    tmp_df <-rbind(tmp_df, tmp2)
    tmp <- as.data.frame(data[[i]][indices == 3, ]) 
    tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 3),digits=2)) %>% mutate(fooditem = colnames(tmp), class = rep(3,ncol(data[[i]])))
    tmp_df <-rbind(tmp_df, tmp2)
    
    tmp_plot[[i]] <- ggplot(data = tmp_df, aes(x = fooditem, y = prop, fill=as.factor(class))) + 
      geom_bar(stat="identity", position = "stack") + 
      geom_text(
        aes(label = ifelse(prop >0.1,prop,"")),
        size = 2,
        position = position_stack(vjust = 0.5),
        inherit.aes = TRUE) +
      ggtitle(timepoints_names[i]) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(limits = c(0, 3)) +
      labs(y = "Cumulative prop", x = "food items") + coord_flip() +
      theme_bw() + 
      scale_fill_discrete(name = "Dietary class")
  }
  if(stop == 3) {
    plot_grid(tmp_plot[[1]], tmp_plot[[2]], tmp_plot[[3]],  
              nrow = 1, ncol = 3)
  } else {
    combined <- tmp_plot[[1]] + tmp_plot[[2]] + tmp_plot[[3]] +tmp_plot[[4]] + tmp_plot[[5]] +tmp_plot[[6]]+ plot_layout(guides = "collect") & theme(legend.position = "bottom")
    print(combined)
    dev.off()
    print(combined)
  }
}


# plotFooditemsCount for Publication 
  #  - arrange/group food items according to food groups (maybe even colour the grouped food items in separate colours) 
  #-  test for significant differences in proportions between dietary clusters 
  # - and potentially only visualise the food items that significantly differ.

  
  # food groups:
  #vegetables: "cauliflower", "carrots", "broccoli", "potato"
  #fruits: "banana", "pear", "apple", "melon", "peach", "kiwi", "orange", "strawberry", "tomato"
  #legumes/nuts: "beans", "porr" 
  # meat: "meat"
  #fish: "fish"
  #eggs" "egg"
  #diary: "milk","yoghurt", "cheese"
  #grains: "bread", "rice", "pasta"
  #soyprod
  #pudding
  #margarine
  #butter

tmp_plot <- list(length = 6)
for(i in 1:6) {
  indices <-  as.vector(clusters_results_solids) 
  tmp <- as.data.frame(timepoints_general_4m_to_14m[[i]][indices == 1, ]) 
  tmp_df <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 1),digits=2), number = colSums(tmp)) %>% mutate(fooditem = colnames(tmp), class = rep(1,ncol(timepoints_general_4m_to_14m[[i]])))
  tmp <- as.data.frame(timepoints_general_4m_to_14m[[i]][indices == 2, ]) 
  tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 2),digits=2), number = colSums(tmp)) %>% mutate(fooditem = colnames(tmp), class = rep(2,ncol(timepoints_general_4m_to_14m[[i]])))
  tmp_df <-rbind(tmp_df, tmp2)
  tmp <- as.data.frame(timepoints_general_4m_to_14m[[i]][indices == 3, ]) 
  tmp2 <- tmp %>% summarise(prop = round(colSums(tmp)/sum(indices == 3),digits=2), number = colSums(tmp)) %>% mutate(fooditem = colnames(tmp), class = rep(3,ncol(timepoints_general_4m_to_14m[[i]])))
  tmp_df <-rbind(tmp_df, tmp2)
  tmp_df$prop <- tmp_df$prop + 0.001

  ## here statistical test between groups for significance
  fi_significance <- vector(length = 28)
  fi_stat <- vector(length = 28)
  for(j in 1:28) {
    tmp <- tmp_df %>% filter(fooditem == unique(tmp_df$fooditem)[j])
    res <- prop.test(x = tmp$number, n = c(51, 27, 34))
    fi_stat[j] <- res$statistic
    fi_significance[j] <- ifelse(is.na(res$p.value), 1, res$p.value)
  }
  # Print the statistics
  #cat(paste(round(fi_stat, 3),"\n"))
  #cat(paste(round(fi_significance, 3),"\n"))
  
  ###
  fi_order <- c("apple", "banana", "kiwi", "melon", "orange", "peach", "pear",  "strawberry", "tomato", 
                   "broccoli", "cauliflower", "carrots", "potato",
                   "beans", "porr", "meat", "fish", "egg", "cheese","milk","yoghurt", "bread", "pasta", "rice",
                   "soyprod", "pudding", "margarine", "butter")
  fi_colors <- c("orange","orange","orange","orange","orange","orange","orange","orange","orange",
                    "darkgreen","darkgreen","darkgreen","darkgreen",
                    "wheat2", "wheat2", "firebrick3", "slateblue3", "salmon", "steelblue2", "steelblue2","steelblue2",
                    "chocolate4","chocolate4","chocolate4", "black", "black", "black", "black" )

  # VERSION A: Only significantly different food items (comment the other one)
  # tmp_plot[[i]] <- ggplot(data = tmp_df %>% filter(fooditem %in% unique(tmp_df$fooditem)[fi_significance < 0.05]),
  #                         aes(x = factor(fooditem, level = fi_order), y = prop, fill=as.factor(class))) +
  #   geom_bar(stat="identity", position = "stack") +
  #   geom_text(
  #     aes(label = ifelse(prop >0.2,round(prop, digits = 2),"")),
  #     size = 2,
  #     position = position_stack(vjust = 0.5),
  #     inherit.aes = TRUE) +
  #   ggtitle(timepoints_names[i]) +
  #   labs(y = "Proportion in a dietary class", x = "Food items") +
  #   theme_bw() +
  #   scale_y_continuous(limits = c(0,3)) + 
  #   theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  #   coord_flip() +
  #   scale_fill_discrete(name = "Dietary class")
  
  # VERSION B: All food items (comment the other one)
  tmp_plot[[i]] <- ggplot(data = tmp_df,
                          aes(x = factor(fooditem, level = fi_order), y = prop, fill=as.factor(class))) +
    geom_bar(stat="identity", position = "stack") +
    geom_text(
      aes(label = ifelse(prop >0.2,round(prop, digits = 2),"")),
      size = 2,
      position = position_stack(vjust = 0.5),
      inherit.aes = TRUE) +
    ggtitle(timepoints_names[i]) +
    labs(y = "Proportion in a dietary class", x = "Food items") +
    theme_bw() +
    #coord_cartesian(y_lim = c(0,3)) +
    scale_y_continuous(limits = c(0,3.1)) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, colour = fi_colors)) +
    coord_flip() +
    scale_fill_discrete(name = "Dietary class")
}

tmp_plot[[1]] + tmp_plot[[2]] + tmp_plot[[3]] +tmp_plot[[4]] + tmp_plot[[5]] +tmp_plot[[6]] +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")


calcStats <- function(data, variable, groups) {
  colnames(data)[which(colnames(data) == variable)] <- "variable"
  colnames(data)[which(colnames(data) == groups)] <- "groups"
  for(t in timepoints_names) {
    response <- (data %>% filter(month == t))[, "variable"]
    factor <- (data %>% filter(month == t))[, "groups"]
    if(bartlett.test(response, factor)$p.value >0.05) { # if >0.05 than groups have the same variance
      if(shapiro.test(response)$p.value <0.05) {#if <0.05 than normal  (perform if bartlett.text >0.05
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

coef_heatmap(df = df, x = "variable", y = "groups", coef = "Estimate", abs_coef_max = 30) +
  expand_limits(y = c(0.5, 2)) + 
  coord_flip() +
  scale_x_discrete(limits = rev(df$variable)) +
  ylab("Coefficient") +
  theme_minimal() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
  ) +
  geom_text(aes(label = Estimate)) +
  geom_text(aes(y = 0.4, label = signif), vjust = 1)
dev.off()

coef_heatmap <- function(df, x, y,
                         abs_coef_max = 2,
                         min_size = 0.2,
                         max_size = 0.9,
                         coef = "Coefficient",
                         oob_fun = scales::oob_squish) {
  # type checks
  stopifnot(is.data.frame(df))
  stopifnot(is.character(x))
  stopifnot(is.character(y))
  stopifnot(is.numeric(abs_coef_max))
  stopifnot(is.numeric(min_size))
  stopifnot(is.numeric(max_size))
  stopifnot(is.character(coef))
  stopifnot(rlang::is_function(oob_fun))
  
  # scale absolute values of coefficients, for tile sizing scaling
  df[["abs_coef"]] <-
    df[[coef]] %>%
    abs() %>%
    oob_fun(range = c(0, abs_coef_max)) %>%
    scales::rescale(to = c(min_size, max_size))
  
  # make plot - will need customising later
  p <- df %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]])) +
    geom_tile(aes(fill = .data[[coef]], height = 0.5, width = 0.9)) +
    scale_fill_distiller(
      type = "div", oob = oob_fun,
      limits = c(-abs_coef_max, abs_coef_max)
    ) +
    theme_bw() +
    theme(panel.grid = element_line(linewidth = 0.1))
  
  return(p)
}

