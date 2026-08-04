library(microViz)
library(ComplexHeatmap)
library(circlize)

#Part 1====
#DMM species best model ====
#Dirichlet Multinomial Mixture (DMM) clustering, an unsupervised clustering method that uses Laplace approximation to 
#identify groups of communities (enterotypes) with similar composition
#Species table should be filtered by abundance and prevalence. Some rare species can screw the permutations 
#and lead to different cluster number

#Filtering only abundant and prevalent species, keep only taxa belonging to genera that have
#Filtering and keeping only those with a prevalence of 5% or more, and with relative abundance over 0.05% in all samples
phyloseq_big <- readRDS("generated_data/phyloseq_object.RDS")
phylo_species <- phyloseq_big %>%
  tax_filter(min_prevalence = 0.05, min_total_abundance  = 0.0001, tax_level = "Species")
species <- as.data.frame(phylo_species@otu_table)


#DMM app function creation
dmm_app<-function(df,only.clust=F, interactive=T){
  require(DirichletMultinomial)
  require(parallel)
  require(xtable)
  require(gtools)
  require(gdata)
  outpath<-invisible(tcltk::tk_choose.dir(default = getwd(), caption = "Select output directory"))
  dmm_otu_df<-as.matrix(df)
  colnames(dmm_otu_df) <- str_match(colnames(dmm_otu_df), "s__(.*)")[,2]
  fit<-mclapply(1:7, dmn, count=dmm_otu_df, verbose=T)
  lplc <- data.frame(Dir_Comp= seq(1:7),Value=sapply(fit, laplace))
  if(only.clust==F){
    p<-ggplot(lplc, aes(Dir_Comp, Value))+
      xlab("Number of Dirichlet Components")+
      ylab("Model Fit")+
      geom_line()+
      geom_point()
    print(p)
    invisible(readline(prompt="Press [enter] to continue\nPress [Esc] to interrupt\n"))
    best<-fit[[which.min(lplc$Value)]]
    p0 <- fitted(fit[[1]], scale=TRUE)
    p4 <- fitted(best, scale=TRUE)
    print(p4)
    colnames(p4) <- paste("m", 1:dim(p4)[2], sep="")
    meandiff <- colSums(abs(p4 - as.vector(p0)))
    diff <- rowSums(abs(p4 - as.vector(p0)))
    o.diff <- order(diff, decreasing=TRUE)
    c.diff <- cumsum(diff[o.diff]) / sum(diff)
    df_small <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 30)
    View(df_small)
    invisible(readline(prompt="Press [enter] to continue\nPress [Esc] to interrupt\n"))
    try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
    sav.ans<-invisible(readline(prompt="Save Table and plot? (y/n)"))
    if(startsWith(sav.ans,"y")){
      require(svglite)
      svglite(paste0(outpath,"/DMM_entero_bacteria_species_heatmap.svg"))
      try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
      dev.off()
      write.table(df_small, paste0(outpath,"/DMM_top30_driving_bacterial_species.txt"),quote = F, sep = "\t", col.names = T)
    }
    invisible(readline(prompt="Press [enter] to continue\nPress [Esc] to interrupt\n"))
    return(mixture(best, assign = TRUE))
  }
  else {
    best<-fit[[which.min(lplc$Value)]]
    p0 <- fitted(fit[[1]], scale=TRUE)
    p4 <- fitted(best, scale=TRUE)
    View(p4)
    colnames(p4) <- paste("m", 1:dim(p4)[2], sep="")
    meandiff <- colSums(abs(p4 - as.vector(p0)))
    diff <- rowSums(abs(p4 - as.vector(p0)))
    o.diff <- order(diff, decreasing=TRUE)
    c.diff <- cumsum(diff[o.diff]) / sum(diff)
    if(interactive==T){
      disp<-invisible(readline(prompt="Show Table and Heatmap?(y/n) "))
      if(startsWith(disp,"y",ignore.case = T)){
        df <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 10)
        View(df)
        try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
      }
      else {return(mixture(best, assign = TRUE))}
    }
    else {
      df <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 10)
      View(df)
      try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
      return(mixture(best, assign = TRUE))
    }
  }
}





#Creating DMM clusters
library(vegan)
library(ggplot2)
library(gdata)
library(gtools)
library(tidyverse)

OTU<-t(species)

#IMPORTANT: a window will pop, chose a folder to save data, otherwise the script won't run
#Plots 120-122
set.seed(190)
DMM_Cluster<-dmm_app(OTU,only.clust = F, interactive = T) 
#DMM_Cluster <- read.delim("DMM/22_DMM_cluster_bacteria_species_table.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))


# Change numbering of clusters:
tmp_dir <- c("1" = 4, "2" = 1, "3" = 2, "4" = 6, "5" = 5, "6" = 3)
tmp_DMM_Cluster <- tmp_dir[as.double(DMM_Cluster$DMM_Cluster)]
DMM_Cluster$DMM_Cluster <- tmp_DMM_Cluster



### My own code, to create nice DMM heatmap with mean composition of top 30 taxa and divercity indices for it
# data frame
# 1st: run above code to get 'best; object
library(plyr)
df_DMM <- rbind.fill(data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1) %>%
                                          tax_transform("compositional", rank = "Species"))@otu_table))), 
                 data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2) %>%
                                          tax_transform("compositional", rank = "Species"))@otu_table))),
                 data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3) %>%
                                          tax_transform("compositional", rank = "Species"))@otu_table))),
                 data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4) %>%
                                          tax_transform("compositional", rank = "Species"))@otu_table))),
                 data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5) %>%
                                          tax_transform("compositional", rank = "Species"))@otu_table))),
                 data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6) %>%
                                          tax_transform("compositional", rank = "Species"))@otu_table))))
df_DMM[is.na(df_DMM)] <- 0
colnames(df_DMM) <- gsub("_", " ", x = colnames(df_DMM))


df_small <- read.table("DMM/DMM_top30_driving_bacterial_species.txt")
rownames(df_small) <- gsub("_", " ", x = rownames(df_small))

Heatmap(t(df_DMM[,rownames(df_small)]), column_split = c(1:6), border = TRUE,
        col = colorRamp2(c(0, 0.8), c("#ffffff", "#005824")),
        cluster_rows = FALSE, cluster_columns = FALSE, 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(t(df_DMM[,rownames(df_small)])[i, j], digits = 3), x, y, gp = gpar(fontsize = 8))
        },
        column_names_gp = gpar(cex = 1), column_names_rot = 0,
        heatmap_legend_param = list(title = "Mean RA"),
        column_labels = paste0("n=",as.vector(table(DMM_Cluster))),
        top_annotation =  HeatmapAnnotation(bar = 1:6, 
                                            col = list(bar = c("1" = "#796248", "2" = "#a47852", "3" = "#cd9c58", "4" = "#e6cd98", "5" = "#f8edcf", "6" = "#d4ede8")),
                                            show_annotation_name = FALSE,
                                            show_legend = FALSE)) 


# -> Alpha diversity
df_dmm_diversity <- rbind(data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1))@otu_table)), dmm = 1),
                        data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2))@otu_table), index = "shannon"), 
                                   specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2))@otu_table)), dmm = 2),
                        data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3))@otu_table), index = "shannon"), 
                                   specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3))@otu_table)), dmm = 3),
                        data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4))@otu_table), index = "shannon"),
                                   specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4))@otu_table)), dmm = 4),
                        data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5))@otu_table), index = "shannon"),
                                   specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5))@otu_table)), dmm = 5),
                        data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6))@otu_table), index = "shannon"),
                                   specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6))@otu_table)), dmm = 6))

df_dmm_diversity <- df_dmm_diversity %>% mutate(dmm = factor(dmm, levels = c(1:6)))
p1 <- ggplot(data = df_dmm_diversity, aes(x=dmm, y=shannon, col = dmm)) + 
  geom_boxplot() +
  labs(x = "DMM cluster", y = "Shannon diversity") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = c(1:6)) +
  scale_color_manual(values = c("#796248","#a47852","#cd9c58", "#e6cd98","#f8edcf", "#d4ede8")) +
  theme(legend.position = "none")

p2 <-  ggplot(data = df_dmm_diversity, aes(x=dmm, y=specnr, col = dmm)) + 
  geom_boxplot() +
  labs(x = "DMM cluster", y = "Species number") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = c(1:6)) +
  scale_color_manual(values = c("#796248","#a47852","#cd9c58", "#e6cd98","#f8edcf", "#d4ede8")) +
  theme(legend.position = "none")

library(patchwork)
p1 + p2 + plot_layout(ncol = 2)

###

#creating MDSplot
library(viridis)
set.seed(190)
x<-metaMDS(OTU, distance = "bray")
points.x<-as.data.frame(x$points)
DMM_Cluster <- as.factor(DMM_Cluster)
#Plot_123
set.seed(190)
ggplot(points.x, aes(MDS1, MDS2, color=DMM_Cluster))+geom_point()+ggtitle("DMM_Clusters")+
  theme_bw()+scale_color_viridis(discrete = T, option = "viridis")

#save DMM_Cluster information per sample - change output folder accordingly
names(DMM_Cluster) <- phylo_species@sam_data$kindcode
write.table(as.data.frame(DMM_Cluster), "DMM/DMM_cluster_bacteria_species_table.txt", quote = F, row.names = T, col.names = T, sep = "\t")





#Part2====
library("ggplot2")
library("igraph")
library("slam")
library("scales")
library("dplyr")

#input====
#IMPORTANT: be sure you upload the correct cluster count
#create a file for transition. Should contain sample_name,family code = child code, age in weeks, and cluster form DMM
#dmm <- read.delim("DMM/DMM_cluster_bacteria_species_table.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
dmm <- DMM_Cluster
dmm <- tibble::rownames_to_column(dmm, "SampleID")
dmm$FamilyID <- as.integer(str_sub(dmm$SampleID, start = 1, end = 4))
dmm$DMM_Cluster <- as.character(dmm$DMM_Cluster)
dmm$Age_days_real <- as.integer(phylo_species@sam_data$days_from_birth)
dmm$diet_class <- as.character(phylo_species@sam_data$diet_class)
#FamilyID (child / family code) must be int
#Age_days_real (time points in days) must be  int
#DMM_Cluster (calculated wuth DMM function clusters) must be chr
for(class in 1:3) {
  md <- dmm %>% filter(diet_class == class)
  
  #modify====
  #Modify these values as you please, also further inthe script there are $column_name that should be modified to your own
  TIME_COLUMN     <- "Age_days_real"
  UNIQUE     <- "SampleID"
  TIME_POINTS     <- c(119, 153, 182, 273, 342, 425)
  TIME_FUDGE      <- 41            # +/- 1.5 months
  REQUIRE_ALL_TP  <- FALSE           # Drop subjects with missing time points
  CLUSTER_COLUMN  <- "DMM_Cluster"
  SUBJECT_COLUMN  <- "FamilyID"
  MINIMUM_PERCENT <- .004
  NODE_COLOR      <- "#072b5a"
  EDGE_COLOR      <- c("white","orange", "darkorange","darkorange", "darkred")
  EDGE_PCT_RANGE  <- c(0,1) #this is a frequency 0-> 100% (0.5 will be 50%, etc.)
  
  #error: Error in rgb(edgeColorRamp(edge_weights)/255) :color intensity NA, not in [0,1] -> adjust EDGE_PCT_RANGE between 0 and 0.9
  #error: something with parameters for graph -> restart R
  #error: Error in plot.new() : figure margins too large -> clean plots enviroment
  
  #double====
  #Double check for the user that they provided valid column names
  vars <- c(TIME_COLUMN, CLUSTER_COLUMN, SUBJECT_COLUMN ,UNIQUE)
  if (!all(vars %in% colnames(md))) {
    stop(sprintf("Column name not found: %s", paste(collapse=", ", setdiff(vars, colnames(md))))) }
  if (nrow(md) == 0) {
    stop("Metadata object 'md' doesn't have any rows.") }
  
  clusterNames <- unique(as.character(sort(md[[CLUSTER_COLUMN]])))
  timePtNames  <- unique(as.character(sort(TIME_POINTS)))
  stateNames   <- apply(expand.grid(clusterNames, timePtNames), 1L, paste, collapse="@")
  nStates      <- length(stateNames)
  
  #filter====
  #Filter out rows with missing data for time, subject, or cluster
  md <- md[!is.na(md[[TIME_COLUMN]]) & !is.na(md[[CLUSTER_COLUMN]])& !is.na(md[[UNIQUE]]) & !is.na(md[[SUBJECT_COLUMN]]), vars, drop=FALSE]
  
  #map====
  #Map samples to time points. For multiple, retain only closest
  # what i can also do: take timepoint form phyloseq and no mapping needed then (change to days)
  closestTP <- sapply(md[[TIME_COLUMN]], function (x) {
    x <- TIME_POINTS[which.min(abs(x - TIME_POINTS))]
  })
  residuals <- abs(closestTP - md[[TIME_COLUMN]])
  
  md[[TIME_COLUMN]] <- closestTP
  
  if (TIME_FUDGE <  1) inRange <- residuals <= TIME_FUDGE * closestTP
  if (TIME_FUDGE >= 1) inRange <- residuals <= TIME_FUDGE
  
  #IMPORTNANT! this step may remove samples that have values that differ a lot from given range
  md        <- md[inRange,,drop=FALSE]
  residuals <- residuals[inRange]
  
  md <- md[order(residuals),,drop=FALSE]
  #this steps will remove samples that are too close together on time scale from the same child/family/subject
  md <- md[!duplicated(paste(md[[TIME_COLUMN]], md[[SUBJECT_COLUMN]])),,drop=FALSE]
  
  
  #optionally====
  #Optionally require a subject to have samples from all the time points
  # if (REQUIRE_ALL_TP) {
  #   md <- plyr::ddply(md, SUBJECT_COLUMN, function (x) {
  #     if (nrow(x) == length(TIME_POINTS)) return (x)
  #     return (NULL)
  #   })
  # }
  #
  timePtCounts <- table(md[[TIME_COLUMN]])
  nodeCounts   <- unlist(as.list(table(paste(sep="@", md[[CLUSTER_COLUMN]], md[[TIME_COLUMN]]))))
  nodeCounts   <- setNames(nodeCounts / timePtCounts[sub("^.*@", "", names(nodeCounts))], names(nodeCounts))
  
  #count====
  #Count the number of subjects at each time point and/or cluster
  #table(md$Age_days_real, dnn = TIME_COLUMN)
  #table(md$DMM_Cluster,   dnn = CLUSTER_COLUMN)
  #table(md$DMM_Cluster, md$Age_days_real, dnn = c(CLUSTER_COLUMN, TIME_COLUMN))
  #table(md$DMM_Cluster, md$Age_days_real,  md$diet_class, dnn=c(CLUSTER_COLUMN, TIME_COLUMN, DIET_CLASS)) # dietary classes
  
  #assemble matrix====
  #Assemble a matrix to represent the number of each transition
  transitionMatrix <- matrix(0, nrow=nStates, ncol=nStates, dimnames=list(stateNames, stateNames))
  
  md[[CLUSTER_COLUMN]] <- paste(sep="@", md[[CLUSTER_COLUMN]], md[[TIME_COLUMN]])
  md[[TIME_COLUMN]]    <- as.numeric(factor(md[[TIME_COLUMN]]))
  
  plyr::d_ply(md[,vars], SUBJECT_COLUMN, function (x) {
    
    if (nrow(x) < 2) return (NULL)
    
    x <- x[order(x[[TIME_COLUMN]]),,drop=FALSE]
    
    for (i in seq_len(nrow(x) - 1)) {
      
      t1 <- x[i,   TIME_COLUMN]
      t2 <- x[i+1, TIME_COLUMN]
      if (t2 - t1 > 1) next
      
      c1 <- as.character(x[i,   CLUSTER_COLUMN])
      c2 <- as.character(x[i+1, CLUSTER_COLUMN])
      transitionMatrix[c1, c2] <<- transitionMatrix[c1, c2] + 1
    }
  })
  
  #rescale====
  #Rescale the transition matrix to percentages
  indices <- sort(rep(1:length(timePtNames), length(clusterNames)))
  for (timePt in 1:(length(timePtNames) - 1)) {
    i <- which(indices == timePt)
    transitionMatrix[i,] <- transitionMatrix[i,] / sum(transitionMatrix[i,])
  }
  
  #drop====
  #Drop cluster/timePt combos with zero observations (after filtering)
  y<-transitionMatrix
  no0<-y[y>0]
  quantile(no0,c(0,0.1,0.9,1))
  MINIMUM_PERCENT <- .0042
  y[which(y < MINIMUM_PERCENT)] <- 0
  #pheatmap::pheatmap(y,cluster_rows = F,cluster_cols = F)
  transitionMatrix[which(transitionMatrix < MINIMUM_PERCENT)] <- 0
  
  stateNames   <- stateNames[rowSums(transitionMatrix) | colSums(transitionMatrix)]
  clusterNames <- clusterNames[sapply(clusterNames, function (x) any(grep(sprintf("^%s@", x), stateNames)))]
  timePtNames  <- timePtNames[ sapply(timePtNames,  function (x) any(grep(sprintf("@%s$", x), stateNames)))]
  
  transitionMatrix <- transitionMatrix[stateNames, stateNames, drop=FALSE]
  
  
  #generate====
  #Generate the figure, axis labels, and color bar legend
  dev.new()
  g <- igraph::graph.adjacency(transitionMatrix, mode="upper", weighted=TRUE)
  layout <- matrix(unlist(strsplit(stateNames, "@", fixed=TRUE)), ncol=2, byrow=TRUE)
  layout <- matrix(ncol=2, c(
    as.numeric(factor(layout[,2], levels=timePtNames)), 
    as.numeric(factor(layout[,1], levels=rev(clusterNames))) ))
  vertex_weights <- scales::rescale(nodeCounts[names(V(g))]%>%as.numeric())
  edge_weights   <- scales::rescale(E(g)$weight, from=EDGE_PCT_RANGE) #bug is here
  edgeColorRamp <- colorRamp(EDGE_COLOR)
  node_colors <- c("1" = "#796248", "2" = "#a47852", "3" = "#cd9c58", "4" = "#e6cd98", "5" = "#f8edcf", "6" = "#d4ede8", "7" = "#9cd6ce")
  node_to_plot <- data.frame(dmm = as.vector(str_sub(names(nodeCounts),1,1)), node = as.vector(str_sub(names(nodeCounts),3,5)), col = as.vector(node_colors[str_sub(names(nodeCounts),1,1)]))
  plot.igraph(
    x                = g,
    xlim             = c(0, .8),
    ylim             = c(-1, 1),
    layout           = layout,
    vertex.label     = NA,
    vertex.size      = scales::rescale(sqrt(vertex_weights / 3.14)) * 25,
    #vertex.color     = NODE_COLOR,
    #vertex.color     = (node_to_plot %>% arrange(node))$col,
    #vertex.color = col_cpy,
    edge.width       = edge_weights * 35 + 1,#change 30 to other numbers to find perfect width of the edges=lines
    edge.color       = rgb(edgeColorRamp(edge_weights) / 255) 
  )
  # c("1" = "#796248", "2" = "#a47852", "3" = "#cd9c58", "4" = "#e6cd98", "5" = "#f8edcf", "6" = "#d4ede8", "7" = "#9cd6ce")
  xLabelPos = (seq_along(timePtNames)  - 1) *  2 / (length(timePtNames)  - 1) - 1
  yLabelPos = (seq_along(clusterNames) - 1) * -2 / (length(clusterNames) - 1) + 1
  
  text(-1.3, yLabelPos, clusterNames)
  text(xLabelPos, 1.5, c("4m","5m","6m","9m","11m","14m"))
  text(xLabelPos, 1.3, sprintf("n = %i", timePtCounts), cex=.6)
  
  with(
    cbreaks(EDGE_PCT_RANGE, pretty_breaks(10), percent),
    legend(
      x      = 1.5, 
      y      = 1, 
      fill   = rgb(edgeColorRamp(scales::rescale(breaks)) / 255), 
      legend = labels, 
      title  = "Transition\nFrequency", 
      bty    = "n" ))
  dev.copy2eps(device = 'postscript', file=paste0("thesis_plots_microbiome/16_Transitions_Over_Time_",class,".eps"), width = 15, ,height = 10, onefile=TRUE)
}

#save the model
#dev.copy2pdf(device = 'postscript', file="DMM/Transitions_Over_Time.pdf", height = 3.75, onefile=TRUE)













