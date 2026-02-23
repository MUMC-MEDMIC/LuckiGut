# load libraries ====
# DMM transition ====
library("ggplot2")
library("igraph")
library("slam")
library("scales")
library("dplyr")
library("plyr")

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ==== 
cluster <- read.delim("Output/Files/S4.1_DMM_cluster_species_table.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
age <- read.delim("Input/17_age.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data ====
#optional step if rownames are the sample names
cluster <- tibble::rownames_to_column(cluster, "SampleID")
#file from meta with only sample codes and age in days
colnames(cluster)
#calcualte the real age in days between sample and DOB
#merge tables 
dmm <-left_join(cluster,age,by="SampleID")
dmm$DMM_Cluster <- as.character(dmm$DMM_Cluster)
dmm$Age_days_real <- as.numeric(dmm$Age_days_real)
#Family_ID (child / family code) must be int
#Age_days_real (time points in days) must be  int
#DMM_Cluster (calculated with DMM function clusters) must be chr

str(dmm)

#rearranging and recording the clusters
dmm$DMM_Cluster_recoded <- mapvalues(dmm$DMM_Cluster, from = c("1", "2", "3", "4", "5", "6","7"), to = c("6", "1", "5", "7", "4","2","3"))

# saving the renamed cluster ====
#write.table(dmm, "S4.2_DMM_cluster_species_table_renamed.txt", quote = F, row.names = T, col.names = T, sep = "\t")

# modify ====
md <- dmm
#Modify these values as you please, also further in the script there are $column_name that should be modified to your own
TIME_COLUMN     <- "Age_days_real"
UNIQUE     <- "SampleID"
TIME_POINTS     <- c(7,28,57,119,153, 182, 273, 342, 425)
TIME_FUDGE      <- 41            # +/- 1.5 months
REQUIRE_ALL_TP  <- FALSE           # Drop subjects with missing time points
CLUSTER_COLUMN  <- "DMM_Cluster_recoded"
SUBJECT_COLUMN  <- "Family_ID"
MINIMUM_PERCENT <- .004
NODE_COLOR      <- "#440154FF" #GL version: #072b5a
EDGE_COLOR      <- c("white","#E0D5E3ff","#C1AAC6ff","#A280AAff","#82568Dff","#632B71ff") #if you have 7 cluster you 
#have to have 6 colours in here, if 6 cluster 5 colours, etc.

EDGE_PCT_RANGE  <- c(0,1) #this is a frequency 0-> 100% (0.5 will be 50%, etc.)

#error: Error in rgb(edgeColorRamp(edge_weights)/255) :color intensity NA, not in [0,1] -> adjust EDGE_PCT_RANGE between 0 and 0.9
#error: something with parameters for graph -> restart R
#error: Error in plot.new() : figure margins too large -> clean plots enviroment

# double ====
#Double check for the user that they provided valid column names
vars <- c(TIME_COLUMN, CLUSTER_COLUMN, SUBJECT_COLUMN,UNIQUE)
if (!all(vars %in% colnames(md)))
  stop(sprintf("Column name not found: %s", paste(collapse=", ", setdiff(vars, colnames(md)))))
if (nrow(md) == 0)
  stop("Metadata object 'md' doesn't have any rows.")

clusterNames <- unique(as.character(sort(md[[CLUSTER_COLUMN]])))
timePtNames  <- unique(as.character(sort(TIME_POINTS)))
stateNames   <- apply(expand.grid(clusterNames, timePtNames), 1L, paste, collapse="@")
nStates      <- length(stateNames)

# filter ====
#Filter out rows with missing data for time, subject, or cluster
md <- md[!is.na(md[[TIME_COLUMN]]) & !is.na(md[[CLUSTER_COLUMN]])& !is.na(md[[UNIQUE]]) & !is.na(md[[SUBJECT_COLUMN]]), vars, drop=FALSE]

# map ====
#Map samples to time points. For multiple, retain only closest
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
#this steps will remove samples that are too close together on time scale fro the same child/family/subject
md <- md[!duplicated(paste(md[[TIME_COLUMN]], md[[SUBJECT_COLUMN]])),,drop=FALSE]

# optionally ====
#Optionally require a subject to have samples from all the time points
if (REQUIRE_ALL_TP) {
  md <- plyr::ddply(md, SUBJECT_COLUMN, function (x) {
    if (nrow(x) == length(TIME_POINTS)) return (x)
    return (NULL)
  })
}
timePtCounts <- table(md[[TIME_COLUMN]])
nodeCounts   <- unlist(as.list(table(paste(sep="@", md[[CLUSTER_COLUMN]], md[[TIME_COLUMN]]))))
nodeCounts   <- setNames(nodeCounts / timePtCounts[sub("^.*@", "", names(nodeCounts))], names(nodeCounts))

# count ====
#Count the number of subjects at each time point and/or cluster
table(md$Age_days_real, dnn = TIME_COLUMN)
table(md$DMM_Cluster,   dnn = CLUSTER_COLUMN)
table(md$DMM_Cluster, md$Age_days_real, dnn = c(CLUSTER_COLUMN, TIME_COLUMN))

# assemble matrix ====
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

# rescale ====
#Rescale the transition matrix to percentages
indices <- sort(rep(1:length(timePtNames), length(clusterNames)))
for (timePt in 1:(length(timePtNames) - 1)) {
  i <- which(indices == timePt)
  transitionMatrix[i,] <- transitionMatrix[i,] / sum(transitionMatrix[i,])
}

# drop ====
#Drop cluster/timePt combos with zero observations (after filtering)
y <-transitionMatrix
no0<-y[y>0]
quantile(no0,c(0,0.1,0.9,1))
MINIMUM_PERCENT <- .0042
y[which(y < MINIMUM_PERCENT)] <- 0
pheatmap::pheatmap(y,cluster_rows = F,cluster_cols = F)
transitionMatrix[which(transitionMatrix < MINIMUM_PERCENT)] <- 0

stateNames   <- stateNames[rowSums(transitionMatrix) | colSums(transitionMatrix)]
clusterNames <- clusterNames[sapply(clusterNames, function (x) any(grep(sprintf("^%s@", x), stateNames)))]
timePtNames  <- timePtNames[ sapply(timePtNames,  function (x) any(grep(sprintf("@%s$", x), stateNames)))]

transitionMatrix <- transitionMatrix[stateNames, stateNames, drop=FALSE]

# generate ====
pdf("S4.2_DMM_transition.pdf", width = w_inc, height = h_inc)

w_inc <- 17 / 2.54 #cm/inc
h_inc <- 9 / 2.54

# Set margins before plotting
old_par <- par(no.readonly = TRUE)  # Save current graphical parameters
par(mar = c(0, 0, 0, 0))  # Adjust these numbers to change margins

#Generate the figure, axis labels, and color bar legend
g <- igraph::graph.adjacency(transitionMatrix, mode="upper", weighted=TRUE)
layout <- matrix(unlist(strsplit(stateNames, "@", fixed=TRUE)), ncol=2, byrow=TRUE)
layout <- matrix(ncol=2, c(
  as.numeric(factor(layout[,2], levels=timePtNames)), 
  as.numeric(factor(layout[,1], levels=rev(clusterNames))) ))
vertex_weights <- scales::rescale(nodeCounts[names(V(g))]%>%as.numeric())
edge_weights   <- scales::rescale(E(g)$weight, from=EDGE_PCT_RANGE) #bug is here
edgeColorRamp <- colorRamp(EDGE_COLOR)
plot.igraph(
  x                = g,
  xlim             = c(-1, 2),#moves the graph in LEFT-RIGHT
  ylim             = c(-2, 1.5),#moves the graph in space
  layout           = layout,
  vertex.label     = NA,
  vertex.size      = scales::rescale(sqrt(vertex_weights / 3.14)) * 25,
  vertex.color     = NODE_COLOR,
  vertex.frame.width = 2,
  asp              = 1, ## Add this line to set a square aspect ratio, this helps a bit with blank spots in nodes, but distroys the overoll plot
  edge.width       = edge_weights * 35 + 1,#change 30 to other numbers to find perfect width of lines
  edge.color       = rgb(edgeColorRamp(edge_weights) / 255) 
)

xLabelPos = (seq_along(timePtNames)  - 1) *  2 / (length(timePtNames)  - 1) - 1
yLabelPos = (seq_along(clusterNames) - 1) * -2 / (length(clusterNames) - 1) + 1

#cex - changes size of "number of samples"
text(-1.3, yLabelPos, clusterNames,cex = 0.6) #change location of cluster numbers
text(xLabelPos, 1.3, timePtNames,cex = 0.6)#changes the location of the bar with "time point names"
#text(xLabelPos, 1.5, sprintf("n = %i", timePtCounts), cex=.7)#changes the location of the bar with "number of samples"
mtext("DMM clusters", side = 2, line = -8, at = 0, cex = 0.5) #in preview will be seen as overlap of x axes on th eplot, but if saved for publication looks ok
mtext("Age in days", side = 3, line = -1, at = 0, cex = 0.5)

#pretty_breaks - changes how much will be the transitional reverences broken 10%, 5% etc.
with(
  cbreaks(EDGE_PCT_RANGE, pretty_breaks(5), percent),
  legend(
    x      = 1.5, 
    y      = 1, #to get side box/legend box up
    fill   = rgb(edgeColorRamp(scales::rescale(breaks)) / 255), 
    legend = labels, 
    title  = "Transition\nFrequency", 
    cex = 0.5,
    bty    = "n" ))

par(old_par)
dev.off()






