# libraries ====
library("ggplot2")
library("igraph")
library("slam")
library("scales")
library("dplyr")
library("microViz")
library("stringr")

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
#IMPORTANT: be sure you upload the correct cluster count
#create a file for transition. Should contain sample_name,family code = child code, age in weeks, and cluster form DMM
dmm <- read.delim("Output/Files/S8.4_DMM_cluster_bacteria_species_table.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
phyloseq_big <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# preparing data ====
phylo_species <- phyloseq_big %>%
  tax_filter(min_prevalence = 0.05, min_total_abundance  = 0.0001, tax_level = "Species")

rm(phyloseq_big)

# preparing for plotting =====
dmm <- tibble::rownames_to_column(dmm, "SampleID")
dmm$FamilyID <- as.integer(str_sub(dmm$SampleID, start = 1, end = 4))
dmm$DMM_Cluster <- as.character(dmm$DMM_Cluster)
dmm$Age_days_real <- as.integer(phylo_species@sam_data$days_from_birth)
dmm$diet_class <- as.character(phylo_species@sam_data$diet_class)
#FamilyID (child / family code) must be int
#Age_days_real (time points in days) must be  int
#DMM_Cluster (calculated wuth DMM function clusters) must be chr

md1 <- dmm %>% filter(diet_class == "1")
md2 <- dmm %>% filter(diet_class == "2")
md3 <- dmm %>% filter(diet_class == "3")

# plotting for class 1 ====
# modify ====
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
  
# double ====
  #Double check for the user that they provided valid column names
  vars <- c(TIME_COLUMN, CLUSTER_COLUMN, SUBJECT_COLUMN ,UNIQUE)
  if (!all(vars %in% colnames(md1))) {
    stop(sprintf("Column name not found: %s", paste(collapse=", ", setdiff(vars, colnames(md1))))) }
  if (nrow(md1) == 0) {
    stop("Metadata object 'md1' doesn't have any rows.") }
  
  clusterNames <- unique(as.character(sort(md1[[CLUSTER_COLUMN]])))
  timePtNames  <- unique(as.character(sort(TIME_POINTS)))
  stateNames   <- apply(expand.grid(clusterNames, timePtNames), 1L, paste, collapse="@")
  nStates      <- length(stateNames)
  
# filter ====
  #Filter out rows with missing data for time, subject, or cluster
  md1 <- md1[!is.na(md1[[TIME_COLUMN]]) & !is.na(md1[[CLUSTER_COLUMN]])& !is.na(md1[[UNIQUE]]) & !is.na(md1[[SUBJECT_COLUMN]]), vars, drop=FALSE]
  
# map ====
  #Map samples to time points. For multiple, retain only closest
  # what i can also do: take timepoint form phyloseq and no mapping needed then (change to days)
  closestTP <- sapply(md1[[TIME_COLUMN]], function (x) {
    x <- TIME_POINTS[which.min(abs(x - TIME_POINTS))]
  })
  residuals <- abs(closestTP - md1[[TIME_COLUMN]])
  
  md1[[TIME_COLUMN]] <- closestTP
  
  if (TIME_FUDGE <  1) inRange <- residuals <= TIME_FUDGE * closestTP
  if (TIME_FUDGE >= 1) inRange <- residuals <= TIME_FUDGE
  
  #IMPORTNANT! this step may remove samples that have values that differ a lot from given range
  md1        <- md1[inRange,,drop=FALSE]
  residuals <- residuals[inRange]
  
  md1 <- md1[order(residuals),,drop=FALSE]
  #this steps will remove samples that are too close together on time scale fro the same child/family/subject
  md1 <- md1[!duplicated(paste(md1[[TIME_COLUMN]], md1[[SUBJECT_COLUMN]])),,drop=FALSE]
  
# optionally ====
  #Optionally require a subject to have samples from all the time points
  # if (REQUIRE_ALL_TP) {
  #   md1 <- plyr::ddply(md1, SUBJECT_COLUMN, function (x) {
  #     if (nrow(x) == length(TIME_POINTS)) return (x)
  #     return (NULL)
  #   })
  # }
  #
  timePtCounts <- table(md1[[TIME_COLUMN]])
  nodeCounts   <- unlist(as.list(table(paste(sep="@", md1[[CLUSTER_COLUMN]], md1[[TIME_COLUMN]]))))
  nodeCounts   <- setNames(nodeCounts / timePtCounts[sub("^.*@", "", names(nodeCounts))], names(nodeCounts))
  
# count ====
  #Count the number of subjects at each time point and/or cluster
  #table(md1$Age_days_real, dnn = TIME_COLUMN)
  #table(md1$DMM_Cluster,   dnn = CLUSTER_COLUMN)
  #table(md1$DMM_Cluster, md1$Age_days_real, dnn = c(CLUSTER_COLUMN, TIME_COLUMN))
  #table(md1$DMM_Cluster, md1$Age_days_real,  md1$diet_class, dnn=c(CLUSTER_COLUMN, TIME_COLUMN, DIET_CLASS)) # dietary classes
  
# assemble matrix ====
  #Assemble a matrix to represent the number of each transition
  transitionMatrix <- matrix(0, nrow=nStates, ncol=nStates, dimnames=list(stateNames, stateNames))
  
  md1[[CLUSTER_COLUMN]] <- paste(sep="@", md1[[CLUSTER_COLUMN]], md1[[TIME_COLUMN]])
  md1[[TIME_COLUMN]]    <- as.numeric(factor(md1[[TIME_COLUMN]]))
  
  plyr::d_ply(md1[,vars], SUBJECT_COLUMN, function (x) {
    
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

# rescale ====
  #Rescale the transition matrix to percentages
  indices <- sort(rep(1:length(timePtNames), length(clusterNames)))
  for (timePt in 1:(length(timePtNames) - 1)) {
    i <- which(indices == timePt)
    transitionMatrix[i,] <- transitionMatrix[i,] / sum(transitionMatrix[i,])
  }
  
# drop ====
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
  
# generate ====
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
    vertex.color     = (node_to_plot %>% arrange(node))$col,
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

dev.copy2pdf(device = 'postscript', file=paste0("S8.5_Transitions_Over_Time_class1",".pdf"), height = 3.75, onefile=TRUE)

# plotting for class 2 ====
# modify ====
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

# double ====
#Double check for the user that they provided valid column names
vars <- c(TIME_COLUMN, CLUSTER_COLUMN, SUBJECT_COLUMN ,UNIQUE)
if (!all(vars %in% colnames(md2))) {
  stop(sprintf("Column name not found: %s", paste(collapse=", ", setdiff(vars, colnames(md2))))) }
if (nrow(md2) == 0) {
  stop("Metadata object 'md2' doesn't have any rows.") }

clusterNames <- unique(as.character(sort(md2[[CLUSTER_COLUMN]])))
timePtNames  <- unique(as.character(sort(TIME_POINTS)))
stateNames   <- apply(expand.grid(clusterNames, timePtNames), 1L, paste, collapse="@")
nStates      <- length(stateNames)

# filter ====
#Filter out rows with missing data for time, subject, or cluster
md2 <- md2[!is.na(md2[[TIME_COLUMN]]) & !is.na(md2[[CLUSTER_COLUMN]])& !is.na(md2[[UNIQUE]]) & !is.na(md2[[SUBJECT_COLUMN]]), vars, drop=FALSE]

# map ====
#Map samples to time points. For multiple, retain only closest
# what i can also do: take timepoint form phyloseq and no mapping needed then (change to days)
closestTP <- sapply(md2[[TIME_COLUMN]], function (x) {
  x <- TIME_POINTS[which.min(abs(x - TIME_POINTS))]
})
residuals <- abs(closestTP - md2[[TIME_COLUMN]])

md2[[TIME_COLUMN]] <- closestTP

if (TIME_FUDGE <  1) inRange <- residuals <= TIME_FUDGE * closestTP
if (TIME_FUDGE >= 1) inRange <- residuals <= TIME_FUDGE

#IMPORTNANT! this step may remove samples that have values that differ a lot from given range
md2        <- md2[inRange,,drop=FALSE]
residuals <- residuals[inRange]

md2 <- md2[order(residuals),,drop=FALSE]
#this steps will remove samples that are too close together on time scale fro the same child/family/subject
md2 <- md2[!duplicated(paste(md2[[TIME_COLUMN]], md2[[SUBJECT_COLUMN]])),,drop=FALSE]

# optionally ====
#Optionally require a subject to have samples from all the time points
# if (REQUIRE_ALL_TP) {
#   md2 <- plyr::ddply(md2, SUBJECT_COLUMN, function (x) {
#     if (nrow(x) == length(TIME_POINTS)) return (x)
#     return (NULL)
#   })
# }
#
timePtCounts <- table(md2[[TIME_COLUMN]])
nodeCounts   <- unlist(as.list(table(paste(sep="@", md2[[CLUSTER_COLUMN]], md2[[TIME_COLUMN]]))))
nodeCounts   <- setNames(nodeCounts / timePtCounts[sub("^.*@", "", names(nodeCounts))], names(nodeCounts))

# count ====
#Count the number of subjects at each time point and/or cluster
#table(md2$Age_days_real, dnn = TIME_COLUMN)
#table(md2$DMM_Cluster,   dnn = CLUSTER_COLUMN)
#table(md2$DMM_Cluster, md2$Age_days_real, dnn = c(CLUSTER_COLUMN, TIME_COLUMN))
#table(md2$DMM_Cluster, md2$Age_days_real,  md2$diet_class, dnn=c(CLUSTER_COLUMN, TIME_COLUMN, DIET_CLASS)) # dietary classes

# assemble matrix ====
#Assemble a matrix to represent the number of each transition
transitionMatrix <- matrix(0, nrow=nStates, ncol=nStates, dimnames=list(stateNames, stateNames))

md2[[CLUSTER_COLUMN]] <- paste(sep="@", md2[[CLUSTER_COLUMN]], md2[[TIME_COLUMN]])
md2[[TIME_COLUMN]]    <- as.numeric(factor(md2[[TIME_COLUMN]]))

plyr::d_ply(md2[,vars], SUBJECT_COLUMN, function (x) {
  
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

# rescale ====
#Rescale the transition matrix to percentages
indices <- sort(rep(1:length(timePtNames), length(clusterNames)))
for (timePt in 1:(length(timePtNames) - 1)) {
  i <- which(indices == timePt)
  transitionMatrix[i,] <- transitionMatrix[i,] / sum(transitionMatrix[i,])
}

# drop ====
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

# generate ====
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
  vertex.color     = (node_to_plot %>% arrange(node))$col,
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

dev.copy2pdf(device = 'postscript', file=paste0("S8.5_Transitions_Over_Time_class2",".pdf"), height = 3.75, onefile=TRUE)

# plotting for class 3 ====
# modify ====
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

# double ====
#Double check for the user that they provided valid column names
vars <- c(TIME_COLUMN, CLUSTER_COLUMN, SUBJECT_COLUMN ,UNIQUE)
if (!all(vars %in% colnames(md3))) {
  stop(sprintf("Column name not found: %s", paste(collapse=", ", setdiff(vars, colnames(md3))))) }
if (nrow(md3) == 0) {
  stop("Metadata object 'md3' doesn't have any rows.") }

clusterNames <- unique(as.character(sort(md3[[CLUSTER_COLUMN]])))
timePtNames  <- unique(as.character(sort(TIME_POINTS)))
stateNames   <- apply(expand.grid(clusterNames, timePtNames), 1L, paste, collapse="@")
nStates      <- length(stateNames)

# filter ====
#Filter out rows with missing data for time, subject, or cluster
md3 <- md3[!is.na(md3[[TIME_COLUMN]]) & !is.na(md3[[CLUSTER_COLUMN]])& !is.na(md3[[UNIQUE]]) & !is.na(md3[[SUBJECT_COLUMN]]), vars, drop=FALSE]

# map ====
#Map samples to time points. For multiple, retain only closest
# what i can also do: take timepoint form phyloseq and no mapping needed then (change to days)
closestTP <- sapply(md3[[TIME_COLUMN]], function (x) {
  x <- TIME_POINTS[which.min(abs(x - TIME_POINTS))]
})
residuals <- abs(closestTP - md3[[TIME_COLUMN]])

md3[[TIME_COLUMN]] <- closestTP

if (TIME_FUDGE <  1) inRange <- residuals <= TIME_FUDGE * closestTP
if (TIME_FUDGE >= 1) inRange <- residuals <= TIME_FUDGE

#IMPORTNANT! this step may remove samples that have values that differ a lot from given range
md3        <- md3[inRange,,drop=FALSE]
residuals <- residuals[inRange]

md3 <- md3[order(residuals),,drop=FALSE]
#this steps will remove samples that are too close together on time scale fro the same child/family/subject
md3 <- md3[!duplicated(paste(md3[[TIME_COLUMN]], md3[[SUBJECT_COLUMN]])),,drop=FALSE]

# optionally ====
#Optionally require a subject to have samples from all the time points
# if (REQUIRE_ALL_TP) {
#   md3 <- plyr::ddply(md3, SUBJECT_COLUMN, function (x) {
#     if (nrow(x) == length(TIME_POINTS)) return (x)
#     return (NULL)
#   })
# }
#
timePtCounts <- table(md3[[TIME_COLUMN]])
nodeCounts   <- unlist(as.list(table(paste(sep="@", md3[[CLUSTER_COLUMN]], md3[[TIME_COLUMN]]))))
nodeCounts   <- setNames(nodeCounts / timePtCounts[sub("^.*@", "", names(nodeCounts))], names(nodeCounts))

# count ====
#Count the number of subjects at each time point and/or cluster
#table(md3$Age_days_real, dnn = TIME_COLUMN)
#table(md3$DMM_Cluster,   dnn = CLUSTER_COLUMN)
#table(md3$DMM_Cluster, md3$Age_days_real, dnn = c(CLUSTER_COLUMN, TIME_COLUMN))
#table(md3$DMM_Cluster, md3$Age_days_real,  md3$diet_class, dnn=c(CLUSTER_COLUMN, TIME_COLUMN, DIET_CLASS)) # dietary classes

# assemble matrix ====
#Assemble a matrix to represent the number of each transition
transitionMatrix <- matrix(0, nrow=nStates, ncol=nStates, dimnames=list(stateNames, stateNames))

md3[[CLUSTER_COLUMN]] <- paste(sep="@", md3[[CLUSTER_COLUMN]], md3[[TIME_COLUMN]])
md3[[TIME_COLUMN]]    <- as.numeric(factor(md3[[TIME_COLUMN]]))

plyr::d_ply(md3[,vars], SUBJECT_COLUMN, function (x) {
  
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

# rescale ====
#Rescale the transition matrix to percentages
indices <- sort(rep(1:length(timePtNames), length(clusterNames)))
for (timePt in 1:(length(timePtNames) - 1)) {
  i <- which(indices == timePt)
  transitionMatrix[i,] <- transitionMatrix[i,] / sum(transitionMatrix[i,])
}

# drop ====
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

# generate ====
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
  vertex.color     = (node_to_plot %>% arrange(node))$col,
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

dev.copy2pdf(device = 'postscript', file=paste0("S8.5_Transitions_Over_Time_class3",".pdf"), height = 3.75, onefile=TRUE)













