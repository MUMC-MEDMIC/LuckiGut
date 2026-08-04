setwd("/Users/p70090895/Desktop/Master Thesis/Practical work")

### 0. LOAD PACKAGES 
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
library(tidyverse) 
library(ggplot2) 
library(readxl) 
library(vegan)
library(reshape2)
library(Skillings.Mack)
library(cowplot) #
library(phyloseq)
library(microbiome)
library(ComplexHeatmap)
library(microViz)
library(ggraph)
library(ggtext) 
library(DT) #
library(corncob) #
library(FSA)
library(ggforce)
library(patchwork)
library(RColorBrewer)
library(plyr)
library(glmmTMB)
library(circlize)
library(rstatix)
library(ggpubr)

### 1. LOAD AND READ DATA
load("generated_data/clusters_solids.RData") # results of clustering 
abundance <- read.table("input_files/MMHP-LCK-UM-780_abundance.tsv", header=TRUE)
read_counts <- read.table("input_files/MMHP-UM-LCK-780_reads_count.tsv", header=TRUE)
code <-  read_excel("input_files/SampleKey_Original_ID_16S_ID_MMHP_ID_JP14072022.xlsx")
code_copy <- code

dim(read_counts)

colnames(read_counts) <- str_remove(colnames(read_counts), regex(".txt.sort$", TRUE))
code_copy$MMHPSampleID <- gsub("-",".",code_copy$MMHPSampleID) 






# sample collection data  
sample_collection <- code <-  read_excel("input_files/sample_collection.xlsx")
sample_collection <- sample_collection %>% 
  mutate(dob_participant = as.POSIXct(dob_participant, format = "%d.%m.%Y"), # change format of date of birth and extract sample age (in days)
         dob_baby = as.POSIXct(dob_baby, format = "%d.%m.%Y"),
         Sample_date = as.POSIXct(Sample_date, format = "%d.%m.%Y")) %>%
  mutate(days = Sample_date - dob_baby)

### 2. DIVIDE INTO TIMEPOINTS
nr_children <- nrow(clusters_solids)
clusters_solids <- clusters_solids %>% dplyr::rename(diet_class = class)

# Check corresponding names and fill lists with readcounts 
timepoints_abundance <- list(length = 6) # for each timepoint 
timepoints_read_counts <- list(length = 6)
not_in_code <- vector("character")
not_in_abundance <- vector("character")
not_in_read_counts <- vector("character")
for(tp in 1:6) {
  timepoints_abundance[[tp]] <- matrix(ncol = nr_children, nrow = nrow(abundance)) 
  timepoints_read_counts[[tp]] <- matrix(ncol = nr_children, nrow = nrow(read_counts))
  for(child in 1:nr_children) {
    tmp_sampleID <- paste(clusters_solids$kindcode[child],".1.",tp+3,sep="")  #kind codes (ID) for this timepoint and child
    tmp_mmhp <- code_copy$MMHPSampleID[which(code_copy$OriginalSampleID == tmp_sampleID)] # corresponding mmhp ID
    # if sampleID -> mmhp ID is empty (not in code), than NA
    if(is_empty(tmp_mmhp)) { 
      not_in_code <- append(not_in_code, tmp_sampleID)
      timepoints_abundance[[tp]][,child] <- rep (NA, nrow(abundance))
      timepoints_read_counts[[tp]][,child] <- rep (NA, nrow(read_counts)) 
    }
    # if sampleID -> mmhp ID is NA (in code), than NA
    else if(is.na(tmp_mmhp)) { 
      not_in_code <- append(not_in_code, tmp_sampleID)
      timepoints_abundance[[tp]][,child] <- rep (NA, nrow(abundance)) 
      timepoints_read_counts[[tp]][,child] <- rep (NA, nrow(read_counts)) 
    }
    else { 
      # if sampleID -> mmhp ID present in input abundance table: than match
      if(tmp_mmhp %in% colnames(abundance)) { timepoints_abundance[[tp]][,child] <- abundance[,tmp_mmhp] }
      else { 
        timepoints_abundance[[tp]][,child] <- rep (NA, nrow(abundance)) 
        not_in_abundance <- append(not_in_abundance, tmp_sampleID)
      }
      # if sampleID -> mmhp ID present in input read_counts table: than match
      if(tmp_mmhp %in% colnames(read_counts)) { timepoints_read_counts[[tp]][,child] <- read_counts[,tmp_mmhp] }
      else { 
        timepoints_read_counts[[tp]][,child] <- rep (NA, nrow(read_counts)) 
        not_in_read_counts <- append(not_in_read_counts, tmp_sampleID)
        }
    }
  }
}
# Set col/row names
for(i in 1:6) {
  colnames(timepoints_read_counts[[i]]) <- clusters_solids$kindcode
  rownames(timepoints_read_counts[[i]]) <- rownames(read_counts)
  colnames(timepoints_abundance[[i]]) <- clusters_solids$kindcode
  rownames(timepoints_abundance[[i]]) <- rownames(abundance)
}
rm(tmp_mmhp)
rm(tmp_sampleID)


### 3. CALCULATE ALPHA DIVERSITY
# 3.1-> DIVERSITY (SHANNON INDEX), based on abundances:
df_diversity <- data.frame(kindcode = clusters_solids$kindcode,
                           diet_class = clusters_solids$diet_class,
                           tp_4m = vegan::diversity(t(timepoints_abundance[[1]]), index = "shannon"),
                           tp_5m = vegan::diversity(t(timepoints_abundance[[2]]), index = "shannon"),
                           tp_6m = vegan::diversity(t(timepoints_abundance[[3]]), index = "shannon"),
                           tp_9m = vegan::diversity(t(timepoints_abundance[[4]]), index = "shannon"),
                           tp_11m = vegan::diversity(t(timepoints_abundance[[5]]), index = "shannon"),
                           tp_14m = vegan::diversity(t(timepoints_abundance[[6]]), index = "shannon"))
 
# 3.2 -> MICROBIAL RICHNESS (Richness represents the number of species observed in each sample), based on read_counts: 
df_richness <- data.frame(kindcode = clusters_solids$kindcode,
                          diet_class = clusters_solids$diet_class,
                          tp_4m = vegan::specnumber(t(timepoints_read_counts[[1]])),
                          tp_5m = vegan::specnumber(t(timepoints_read_counts[[2]])),
                          tp_6m = vegan::specnumber(t(timepoints_read_counts[[3]])),
                          tp_9m = vegan::specnumber(t(timepoints_read_counts[[4]])),
                          tp_11m = vegan::specnumber(t(timepoints_read_counts[[5]])),
                          tp_14m = vegan::specnumber(t(timepoints_read_counts[[6]])))

# 3.3 -> MICROBIAL RICHNESS (Chao1), based on read counts : 
df_chao1 <- data.frame(kindcode = clusters_solids$kindcode,
                          diet_class = clusters_solids$diet_class,
                          tp_4m = vegan::estimateR(t(timepoints_read_counts[[1]]))[2,],
                          tp_5m = vegan::estimateR(t(timepoints_read_counts[[2]]))[2,],
                          tp_6m = vegan::estimateR(t(timepoints_read_counts[[3]]))[2,],
                          tp_9m = vegan::estimateR(t(timepoints_read_counts[[4]]))[2,],
                          tp_11m = vegan::estimateR(t(timepoints_read_counts[[5]]))[2,],
                          tp_14m = vegan::estimateR(t(timepoints_read_counts[[6]]))[2,])

# CREATE DATA FRAME for visualization
tmp_df <- melt(df_diversity) %>% 
  dplyr::rename(shannon = value,
                timepoint = variable) %>%
  mutate(kindcode = ifelse(timepoint == "tp_4m",paste0(kindcode, ".1.4"),
                           ifelse(timepoint == "tp_5m",paste0(kindcode, ".1.5"),
                                  ifelse(timepoint == "tp_6m",paste0(kindcode, ".1.6"),
                                         ifelse(timepoint == "tp_9m",paste0(kindcode, ".1.7"),
                                                ifelse(timepoint == "tp_11m",paste0(kindcode, ".1.8"), paste0(kindcode, ".1.9")))))),
         days_from_birth = sample_collection$days[match(kindcode, sample_collection$Sample_ID)],
         specnumber = (melt(df_richness))$value,
         chao1 = (melt(df_chao1))$value) 

# VISUZALIZE ALPHA DIVERSITY (and calculate statistics): 
# function calcStats: at the end of the script (run before proceeding here)
# -> boxplot
timepoints_names <- c("4m","5m","6m","9m","11m", "14m")
tmp_df <- tmp_df %>%
  mutate(month = str_remove(timepoint,"tp_"))

# Shannon index
pwc <- calcStats(data = tmp_df, variable = "shannon", groups = "diet_class")
write.csv(pwc[,c(".y.","month","group1","group2","n1","n2","p","p.adj","p.adj.signif")],
          file = "13_shannon_stats.csv")
pwc <- pwc %>% # manually added groups, x, xmin and xmax 
  mutate(y.position = get_y_position(formula = shannon ~ diet_class, data = tmp_df %>% group_by(month))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))
pwc$y.position[c(1,2)] <- c(2.7, 2.9)

p1 <- ggplot(data = tmp_df, aes(x=timepoint, y=shannon)) + 
  geom_boxplot(aes(fill = diet_class)) +
  labs(x = "Time point", y = "Shannon index", title =bquote(bold("(A)")~"Shannon index")) +
  theme_bw() +
  scale_fill_discrete(name="Dietary class") + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = timepoints_names) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0) 
  
# Specnumber:
pwc <- calcStats(data = tmp_df, variable = "specnumber", groups = "diet_class")
write.csv(pwc[,c(".y.","month","group1","group2","n1","n2","p","p.adj","p.adj.signif")],
          file = "13_specnumber_stats.csv")
pwc <- pwc %>%  # manually added groups, x, xmin and xmax 
  mutate(y.position = get_y_position(formula = specnumber ~ diet_class, data = tmp_df %>% group_by(month))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))
pwc$y.position[c(2,7,8, 14 )] <- c(70,75,85,100)

p2 <- ggplot(data = tmp_df, aes(x=timepoint, y=specnumber)) + 
  geom_boxplot(aes(fill=diet_class)) +
  labs(x = "Time point", y= "species number", title = bquote(bold("(B)")~"Species richness")) +
  theme_bw() +
  scale_fill_discrete(name="Dietary class") + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = timepoints_names) +
  scale_y_continuous(limits = c(0,120)) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0) 

#Chao1
pwc <- calcStats(data = tmp_df %>% filter(chao1 < 125), variable = "chao1", groups = "diet_class")
write.csv(pwc[,c(".y.","month","group1","group2","n1","n2","p","p.adj","p.adj.signif")],
          file = "13_chao_stats.csv")
pwc <- pwc %>% # manually added groups, x, xmin and xmax 
  mutate(y.position = get_y_position(formula = chao1 ~ diet_class, data = tmp_df %>% filter(chao1 < 125) %>% group_by(month))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))
pwc$y.position[c(2,7,8, 14 )] <- c(65,75,85,100) 

p3 <- ggplot(data = tmp_df %>% filter(chao1 < 125), aes(x=timepoint, y=chao1)) + 
  geom_boxplot(aes(fill=diet_class)) +
  labs(x = "Time point", y= "Chao1", title = bquote(bold("(C)")~"Species richness (Chao1)")) +
  theme_bw() +
  scale_fill_discrete(name="Dietary class") + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = timepoints_names) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0) 
    


# -> smoothing lines plot (loess regressions)
p11 <- ggplot(data = tmp_df, aes(x=days_from_birth, y=shannon, color=diet_class)) + 
  geom_smooth(level = FALSE, show.legend = FALSE) +
  xlab("Time (months)") + ylab("") +
  theme_bw() +
  scale_x_continuous(limits = c(120, 440), breaks = c(124,155, 186,279, 341, 434), labels = c(4,5,6,9,11,14)) +
  geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed", 
             color = "grey", size=.5) 

p21 <- ggplot(data = tmp_df, aes(x=days_from_birth, y=specnumber, color=diet_class)) + 
  geom_smooth(level = FALSE, show.legend = FALSE) +
  xlab("Time (months)") + ylab("") +
  theme_bw() +
  scale_x_continuous(limits = c(120, 440), breaks = c(124,155, 186,279, 341, 434), labels = c(4,5,6,9,11,14)) +
  geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed", 
             color = "grey", size=.5) 

p31 <- ggplot(data = tmp_df %>% filter(chao1 < 125), aes(x=days_from_birth, y=chao1, color=diet_class)) + 
  geom_smooth(level = FALSE, show.legend = FALSE) +
  xlab("Time (months)") + ylab("") +
  theme_bw() +
  scale_x_continuous(limits = c(120, 440), breaks = c(124,155, 186,279, 341, 434), labels = c(4,5,6,9,11,14)) +
  geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed", 
             color = "grey", size=.5) 

# Fig P2A
p1 + p11 + p2 + p21 + p3 + p31 + plot_layout(design = "
  001
  223
  445
", guides = "collect") & theme(legend.position = "bottom")



# CALCULATE STATISTICS
# for each metric:
#   for each class:
#    differences between time points (e.g. 4m-5m)

# ANOVA: normal distribution assumption
for(col in c("shannon", "specnumber")) {
  data <- tmp_df 
  if(col == "chao1") {data <- tmp_df %>% filter(chao1 < 125)}
  print(col)
  for(cl in c(1,2,3)) {
    print(cl)
    response <- (data %>% filter(diet_class == cl))[,col]
    factor <- (data %>% filter(diet_class == cl))[,"timepoint"]
    print(aov(response ~ factor))
    print(pairwise.t.test(response, factor, paired = T, p.adjust.method = "BH"))
  }
}



### 4. MAKE PHYLOSEQ OBJECT
# -> Taxonomy data
taxmat <- matrix(nrow = nrow(read_counts), ncol = 7)
rownames(taxmat) <- rownames(read_counts)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat[,1] <- str_match(rownames(read_counts), "k__(.*?)\\|p__")[,2] #domain
taxmat[,2] <- str_match(rownames(read_counts), "p__(.*?)\\|c__")[,2] #phylum
taxmat[,3] <- str_match(rownames(read_counts), "c__(.*?)\\|o__")[,2] #class
taxmat[,4] <- str_match(rownames(read_counts), "o__(.*?)\\|f__")[,2] #order
taxmat[,5] <- str_match(rownames(read_counts), "f__(.*?)\\|g__")[,2] #family
taxmat[,6] <- str_match(rownames(read_counts), "g__(.*?)\\|s__")[,2] #genus
taxmat[,7] <- str_match(rownames(read_counts), "s__(.*)")[,2] #species
TAX <- tax_table(taxmat)

# -> OTU TABLE
OTU <- cbind(as.matrix(timepoints_read_counts[[1]][ , colSums(is.na(timepoints_read_counts[[1]])) == 0]),
             as.matrix(timepoints_read_counts[[2]][ , colSums(is.na(timepoints_read_counts[[2]])) == 0]),
             as.matrix(timepoints_read_counts[[3]][ , colSums(is.na(timepoints_read_counts[[3]])) == 0]),
             as.matrix(timepoints_read_counts[[4]][ , colSums(is.na(timepoints_read_counts[[4]])) == 0]),
             as.matrix(timepoints_read_counts[[5]][ , colSums(is.na(timepoints_read_counts[[5]])) == 0]),
             as.matrix(timepoints_read_counts[[6]][ , colSums(is.na(timepoints_read_counts[[6]])) == 0]))
colnames(OTU) <- paste("sa",1:ncol(OTU), sep="")
OTU <- otu_table(OTU, taxa_are_rows = TRUE)

# -> Sample data
breastfeeding_sam <- readRDS("generated_data/tmp_breastfeeding_sam.RDS")
formulafeeding_sam <- readRDS("generated_data/tmp_formulafeeding_sam.RDS")  

samdat <- rbind(clusters_solids %>%  mutate(tp = rep("4m", nrow(clusters_solids))) %>% filter(kindcode %in% colnames(timepoints_read_counts[[1]][ , colSums(is.na(timepoints_read_counts[[1]])) == 0]) ) %>% mutate(kindcode = paste0(kindcode, ".1.4")),
                clusters_solids %>%  mutate(tp = rep("5m", nrow(clusters_solids))) %>% filter(kindcode %in% colnames(timepoints_read_counts[[2]][ , colSums(is.na(timepoints_read_counts[[2]])) == 0]) ) %>% mutate(kindcode = paste0(kindcode, ".1.5")),
                clusters_solids %>%  mutate(tp = rep("6m", nrow(clusters_solids))) %>% filter(kindcode %in% colnames(timepoints_read_counts[[3]][ , colSums(is.na(timepoints_read_counts[[3]])) == 0]) ) %>% mutate(kindcode = paste0(kindcode, ".1.6")),
                clusters_solids %>%  mutate(tp = rep("9m", nrow(clusters_solids))) %>% filter(kindcode %in% colnames(timepoints_read_counts[[4]][ , colSums(is.na(timepoints_read_counts[[4]])) == 0]) ) %>% mutate(kindcode = paste0(kindcode, ".1.7")),
                clusters_solids %>%  mutate(tp = rep("11m", nrow(clusters_solids))) %>% filter(kindcode %in% colnames(timepoints_read_counts[[5]][ , colSums(is.na(timepoints_read_counts[[5]])) == 0]) ) %>% mutate(kindcode = paste0(kindcode, ".1.8")),
                clusters_solids %>%  mutate(tp = rep("14m", nrow(clusters_solids))) %>% filter(kindcode %in% colnames(timepoints_read_counts[[6]][ , colSums(is.na(timepoints_read_counts[[6]])) == 0]) ) %>% mutate(kindcode = paste0(kindcode, ".1.9")))
samdat <- samdat %>% mutate(days_from_birth = as.numeric(sample_collection$days[match(samdat$kindcode, sample_collection$Sample_ID )]),
                            breastfeeding = breastfeeding_sam, 
                            formulafeeding = formulafeeding_sam)
row.names(samdat) <- paste("sa",1:ncol(OTU), sep="")
samdat <- sample_data(samdat)
sample_names(samdat) <- paste("sa",1:ncol(OTU), sep="")

# -> Phyloseq object 
phyloseq_big = phyloseq(OTU, samdat, TAX)
phyloseq_big <- phyloseq_validate(phyloseq_big, remove_undetected = TRUE)
phyloseq_big <- phyloseq_big %>%
  tax_filter(min_prevalence = 5) %>% 
  tax_fix() 

phyloseq_big_tmp <- phyloseq_big  %>%
  ps_mutate(
    CLASS_2 = ifelse(diet_class == 2, yes = 1, no = 0), # in reference to class 1
    CLASS_3 = ifelse(diet_class == 3, yes = 1, no = 0) # in reference to class 1
  ) 
#saveRDS(phyloseq_big_tmp, file = "generated_files/phyloseq_object.RDS")




### 5. ORDINATION
## 5.1 -> PCoA 
labels_.diet_class <- c("DIETARY CLASS 1", "DIETARY CLASS 2", "DIETARY CLASS 3")
names(labels_.diet_class) <- c(1,2,3)

# -> Original PCoA
phyloseq_big %>% 
  ps_mutate(tp = factor(tp, levels = c("4m", "5m", "6m", "9m", "11m", "14m"), ordered = T)) %>% # for visualization
  tax_transform("identity", rank = "Species") %>%
  dist_calc("aitchison") %>%  # aitchison distance
  ord_calc("PCoA") %>%
  ord_plot(color = "tp", size = 2,alpha=0.8) +
  geom_mark_ellipse(aes(color = tp), expand = unit(0.5,"mm"), show.legend = FALSE) + 
  facet_wrap(~diet_class, labeller = labeller(diet_class = labels_.diet_class))+
  scale_color_manual("Time point", values = c("#EEE19E", "#4DAF4A",  "#00FFFF", "#0066FF" ,"#FF9900", "#450628")) + # color pallete or colors of timepoints 
  scale_fill_manual("Time point", values = c("#EEE19E", "#4DAF4A",  "#00FFFF", "#0066FF" ,"#FF9900", "#450628")) + 
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = tp, y = tp), orientation = "y", show.legend = F) +
  ggside::geom_ysideboxplot(aes(fill = tp, x = tp), orientation = "x", show.legend = F) +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() +
  geom_text(data = data.frame(label = c("n = 184", "n = 95", "n = 110"), # number of samples in each cluster
                              diet_class   = c(1, 2, 3), 
                              tp = c("9m","9m","9m")),
            mapping = aes(x=-1.75,y=-4.2, label = label), show.legend = FALSE) +
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.text.x = element_text(face = "bold", size = 10)
  )

# -> Changed PCoA for article  (Fig P2B)
data <- phyloseq_big %>% 
  ps_mutate(tp = factor(tp, levels = c("4m", "5m", "6m", "9m", "11m", "14m"), ordered = T)) %>% # for visualization
  tax_transform("identity", rank = "Species") %>%
  dist_calc("aitchison") %>%  # aitchison distance
  ord_calc("PCoA") %>%
  ord_get()

data <- data.frame(MDS1 = data$CA$u[,1], diet_class = phyloseq_big@sam_data$diet_class,
                   tp = factor(phyloseq_big@sam_data$tp, levels = c("4m", "5m", "6m", "9m", "11m", "14m"), ordered = T))

pwc <- data %>%
  group_by(tp) %>%
  pairwise_t_test(MDS1~diet_class,p.adjust.method = "BH")

pwc <- pwc %>%  # manually added groups, x, xmin and xmax 
  mutate(y.position = get_y_position(formula = MDS1 ~ diet_class, data = data %>% group_by(tp))$y.position,
         groups = rep(list(V1 = c("1","2"), V2 = c("1","3"), V3 = c("2","3")),6),
         x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
         xmin = c(0.7333333, 0.7333333, 1.0000000, 1.7333333, 1.7333333, 2.0000000, 2.7333333, 2.7333333, 3.0000000, 3.7333333,
                  3.7333333, 4.0000000, 4.7333333, 4.7333333, 5.0000000, 5.7333333, 5.7333333, 6.0000000),
         xmax = c(1.000000, 1.266667, 1.266667, 2.000000, 2.266667, 2.266667, 3.000000, 3.266667, 3.266667, 4.000000, 4.266667,
                  4.266667, 5.000000, 5.266667, 5.266667, 6.000000, 6.266667, 6.266667))




ggplot(data = data, aes(x=tp, y=MDS1)) + 
  geom_boxplot(aes(fill=diet_class)) +
  labs(x = "Time point", y= "MDS1 [13.8%]") +
  theme_bw() +
  scale_fill_discrete(name="Dietary class") + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = timepoints_names) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0.01, hide.ns = T, bracket.nudge.y = 0) 






### In thesis: 
phyloseq_big %>% 
  ps_mutate(tp = factor(tp, levels = c("4m", "5m", "6m", "9m", "11m", "14m"), ordered = T)) %>% # for visualization
  tax_transform("identity", rank = "Species") %>%
  dist_calc("aitchison") %>%  # aitchison distance
  ord_calc("PCoA") %>%
  ord_plot(color = "tp", size = 2,alpha=0.8) +
  geom_mark_ellipse(aes(color = tp), expand = unit(0.5,"mm"), show.legend = FALSE) + 
  facet_wrap(~diet_class, labeller = labeller(diet_class = labels_.diet_class))+
  #scale_color_manual("Time point", values = c("#EEE19E", "#4DAF4A",  "#00FFFF", "#0066FF" ,"#FF9900", "#450628")) + # color pallete or colors of timepoints 
  #scale_fill_manual("Time point", values = c("#EEE19E", "#4DAF4A",  "#00FFFF", "#0066FF" ,"#FF9900", "#450628")) + 
  theme_bw() +
  labs(color = "Time point") +
  geom_text(data = data.frame(label = c("n = 184", "n = 95", "n = 110"), # number of samples in each cluster
                              diet_class   = c(1, 2, 3), 
                              tp = c("9m","9m","9m")),
            mapping = aes(x=-1.75,y=-4.2, label = label), show.legend = FALSE) 
###









# test statistical differeces between PC's

ord <- (phyloseq_big %>% 
  tax_transform("identity", rank = "Species") %>%
  dist_calc("aitchison") %>%  # aitchison / bray
  ord_calc("PCoA") %>% ord_get())$CA$u[,1:2] %>%
  as.data.frame() %>%
  mutate(diet_class = as.factor(phyloseq_big@sam_data$diet_class),
         timepoint = factor(phyloseq_big@sam_data$tp,levels = timepoints_names, ordered = T))

# Differences between classes at each timepoint:
ord %>%
  group_by(timepoint) %>%
  pairwise_t_test(MDS2~diet_class,p.adjust.method = "BH")

# Differences between timepoint for each class: (not paired, cause not all arguments have the same length)
pwc <- ord %>%
  group_by(diet_class) %>%
  pairwise_t_test(MDS2~timepoint, p.adjust.method = "BH")
print(pwc, n = 45)




## 5.2 -> RDA (with dietary classes as constrains) (Fig Supplementary)
plot_list <- list(length = 6)
for(i in 1:6) {
  plot_list[[i]] <-
    phyloseq_big_tmp %>%
    ps_filter(tp == timepoints_names[i]) %>%
    tax_transform("clr", rank = "Species") %>%
    ord_calc(
      constraints = c("CLASS_2", "CLASS_3"), # reference: CLASS_1
      method = "RDA", 
      scale_cc = FALSE 
    ) %>%
    ord_plot(
      colour = "diet_class", size = 2, alpha = 0.5,
      plot_taxa = 1:8,
      var_renamer = function(x) gsub("_", " ",x)) +
    theme_bw() +
    labs(color = "Dietary class", caption = "", title = timepoints_names[i])
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]] + plot_list[[6]] +
  plot_layout(design = "
  AB
  CD
  EF
", guides = "collect") & theme(legend.position = "bottom")



### 6. TEMPORAL CHANGE IN COMPOSITION - VISUALIZATION
microbial_age_df <- read.csv(file = "generated_data/microbial_age_df.csv") # for duration of breastfeeding

# Choose rank_viz c("Phylum", "Class", "Genus", "Species")
rank_viz <- "Species"


# 6.1. Relative abundance heatmap
tmp_df <- list(length = 3)
for(class in 1:3) { # for each diet_class create a data frame with mean relative abundances of taxa (rows: timepoints for each class, cols: taxa)
  tmp_df[[class]] <- rbind.fill(data.frame(t(rowMeans((phyloseq_big %>% ps_filter(diet_class == class, tp == "4m") %>%
                                                         tax_transform("compositional", rank = rank_viz))@otu_table))),
                                data.frame(t(rowMeans((phyloseq_big %>% ps_filter(diet_class == class, tp == "5m") %>%
                                                         tax_transform("compositional", rank = rank_viz))@otu_table))),
                                data.frame(t(rowMeans((phyloseq_big %>% ps_filter(diet_class == class, tp == "6m") %>%
                                                         tax_transform("compositional", rank = rank_viz))@otu_table))),
                                data.frame(t(rowMeans((phyloseq_big %>% ps_filter(diet_class == class, tp == "9m") %>%
                                                         tax_transform("compositional", rank = rank_viz))@otu_table))),
                                data.frame(t(rowMeans((phyloseq_big %>% ps_filter(diet_class == class, tp == "11m") %>%
                                                         tax_transform("compositional", rank = rank_viz))@otu_table))),
                                data.frame(t(rowMeans((phyloseq_big %>% ps_filter(diet_class == class, tp == "14m") %>%
                                                         tax_transform("compositional", rank = rank_viz))@otu_table))))
  tmp_df[is.na(tmp_df)] <- 0
}
tmp_df <- rbind.fill( tmp_df[[1]], tmp_df[[2]],  tmp_df[[3]]) # rows(class 1: tp1-6, class 2: tp1-6, class 3: tp-1-6)

# Choose most abundant taxa based on threshold (here > 1% at any time point)
taxa <-  colnames(tmp_df[,which(apply(tmp_df > 0.01, 2, any, na.rm=TRUE) == TRUE)])

# Create Heatmap
print(Heatmap(t(tmp_df[,taxa]), column_split = c(rep(1,6), rep(2,6), rep(3,6)),
              col = colorRamp2(c(0, 0.8), c("#ffffff", "#005824")),
              cluster_rows = T, cluster_columns = FALSE,
              border = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(t(tmp_df[,taxa])[i, j], digits = 2), x, y, gp = gpar(fontsize = 8))
              },
              row_labels = gsub("_", " ",rownames(t(tmp_df[,taxa]))),
              column_labels = rep(c("4","5","6","9","11","14"),3),
              column_names_gp = gpar(cex = 1), column_names_rot = 0,
              heatmap_legend_param = list(title = "Relative\nabundance")))

## 6.2. Abundance Line plots

# Merge relative abundances of all samples from each dietary class:
tmp_phylo <- phyloseq_big %>% ps_filter(diet_class == 1) %>%
  tax_transform("compositional", rank = rank_viz)
tmp_1 <- data.frame(t(tmp_phylo@otu_table)) %>%
  mutate(days = tmp_phylo@sam_data$days_from_birth,
         diet_class = tmp_phylo@sam_data$diet_class)

tmp_phylo <- phyloseq_big %>% ps_filter(diet_class == 2) %>%
  tax_transform("compositional", rank = rank_viz)
tmp_2 <- data.frame(t(tmp_phylo@otu_table)) %>%
  mutate(days = tmp_phylo@sam_data$days_from_birth,
         diet_class = tmp_phylo@sam_data$diet_class)

tmp_phylo <- phyloseq_big %>% ps_filter(diet_class == 3) %>%
  tax_transform("compositional", rank = rank_viz)
tmp_3 <- data.frame(t(tmp_phylo@otu_table)) %>%
  mutate(days = tmp_phylo@sam_data$days_from_birth,
         diet_class = tmp_phylo@sam_data$diet_class)

tmp_df2 <- rbind.fill(tmp_1, tmp_2, tmp_3) # bind all dietary classes
tmp_df2[is.na(tmp_df2)] <- 0
tmp_df2 <- melt(tmp_df2, id = c("days", "diet_class")) # change data to a long format for vizualization

# Visualize line plots (loess regression)
# Taxa chosen based on GLMM significance
if(rank_viz == "Genus") {
  variables <- c("Akkermansia", "Alistipes", "Bacteroides", "Blautia", "Eubacterium",
                 "Faecalibacterium", "Fusicatenibacter", "Klebsiella", "Lachnospiraceae_unclassified", "Parabacteroides", "Roseburia", "Veillonella")
  tmp_df2 <- tmp_df2 %>% filter(variable %in% variables)
} else if(rank_viz == "Species") {
  variables <- c("Akkermansia_muciniphila", "Alistipes_finegoldii", "Bacteroides_fragilis", "Bacteroides_ovatus", "Bacteroides_thetaiotaomicron",
                 "Bifidobacterium_adolescentis","Bifidobacterium_longum", "Eubacterium_rectale", "Faecalibacterium_prausnitzii",
                 "Fusicatenibacter_saccharivorans", "Klebsiella_pneumoniae", "Parabacteroides_distasonis", "Parabacteroides_merdae",
                 "Prevotella_copri", "Ruminococcus_gnavus", "Veillonella_parvula")
  tmp_df2 <- tmp_df2 %>% filter(variable %in% variables)
}

# Print line plots with relative abundance (fixed limits)
print(ggplot(data = tmp_df2, aes(x = days, y = value, color = diet_class)) +
        facet_wrap(~variable, scales = "free_y") +
        geom_point(alpha = 0.3) +
        geom_smooth(show.legend = T, se = F) +
        labs(x = "Time (months)", y = "Relative abundance", color = "Dietary class") +
        theme_minimal() +
        #coord_cartesian(ylim = c(0, 0.3), xlim = c(120,440)) + #For Genus
        coord_cartesian(ylim=c(0,0.1), xlim = c(120,440)) + # For Species
        scale_x_continuous(breaks = c(124,155, 186,279, 341, 434), labels = c(4,5,6,9,11,14)) +
        geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed",
                   color = "grey", size=.5) +
        theme(
          strip.background = element_rect(
            color="black", fill="white", size=0.5, linetype="solid"
          ),
          strip.text.x = element_text(
            size = 7, color = "black"
          )
        ))
#dev.off()

### For publication:differnt limits for each plot (Fig P3B)
if(rank_viz == "Species") {
  plots_list <- list(length = 16)
  limits_list <- list(c(0,0.1), c(0,0.075), c(0,0.1), c(0,0.05),
                      c(0,0.05), c(0,0.1), c(0,0.3), c(0,0.05),
                      c(0,0.25), c(0,0.05), c(0,0.05), c(0,0.05),
                      c(0,0.03), c(0,0.2), c(0,0.1), c(0,0.05))
  for(i in 1:16) {
    plots_list[[i]] <- ggplot(data = tmp_df2 %>% filter(variable == variables[i] ), aes(x = days, y = value, color = diet_class)) +
      geom_point(alpha = 0.3) +
      geom_smooth(show.legend = T, se = F) +
      labs(x = "", y = ifelse(i == 5,"Relative abundance",""), color = "Dietary class", title = gsub("_", " ", variables[i])) +
      theme_bw() +
      coord_cartesian(ylim=as.vector(limits_list[[i]]), xlim = c(120,440)) +
      scale_x_continuous(breaks = c(124,155, 186,279, 341, 434), labels = c(4,5,6,9,11,14)) +
      geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed",
                 color = "grey", size=.5) +
      theme(
        plot.title = element_text(
          color="black", face = "bold", hjust = 0.5, size = 10
        ),
        strip.text.x = element_text(
          size = 7, color = "black"
        )
      )
  }
  
  plots_list[[1]] + plots_list[[2]] + plots_list[[3]] + plots_list[[4]] +
    plots_list[[5]] + plots_list[[6]] + plots_list[[7]] + plots_list[[8]] +
    plots_list[[9]] + plots_list[[10]] + plots_list[[11]] + plots_list[[12]] +
    plots_list[[13]] + plots_list[[14]] + plots_list[[15]] + plots_list[[16]] +
    plot_layout(design = "
  ABCD
  EFGH
  IJKL
  MNOP
", guides = "collect") & theme(legend.position = "bottom")
  
  grid::grid.draw(grid::textGrob('Time (months)', x =.5, y=0.07))
  ggsave(file = "P3B_abundance_lines_species.svg", width = 40, height = 30, units = "cm")
  ########
}
# Fig Supplementary
if(rank_viz == "Genus") {
  plots_list <- list(length = 12)
  limits_list <- list(c(0,0.1), c(0,0.075), c(0,0.2), c(0,0.1),
                      c(0,0.05), c(0,0.3), c(0,0.05), c(0,0.1),
                      c(0,0.1), c(0,0.05), c(0,0.1), c(0,0.075))
  for(i in 1:12) {
    plots_list[[i]] <- ggplot(data = tmp_df2 %>% filter(variable == variables[i] ), aes(x = days, y = value, color = diet_class)) +
      geom_point(alpha = 0.3) +
      geom_smooth(show.legend = T, se = F) +
      labs(x = "", y = ifelse(i == 5,"Relative abundance",""), color = "Dietary class", title = gsub("_", " ", variables[i])) +
      theme_bw() +
      coord_cartesian(ylim=as.vector(limits_list[[i]]), xlim = c(120,440)) +
      scale_x_continuous(breaks = c(124,155, 186,279, 341, 434), labels = c(4,5,6,9,11,14)) +
      geom_vline(xintercept = c(124,155, 186,279, 341, 434), linetype="dashed",
                 color = "grey", size=.5) +
      theme(
        plot.title = element_text(
          color="black", face = "bold", hjust = 0.5, size = 10
        ),
        strip.text.x = element_text(
          size = 7, color = "black"
        )
      )
  }
  
  plots_list[[1]] + plots_list[[2]] + plots_list[[3]] + plots_list[[4]] +
    plots_list[[5]] + plots_list[[6]] + plots_list[[7]] + plots_list[[8]] +
    plots_list[[9]] + plots_list[[10]] + plots_list[[11]] + plots_list[[12]] +
    plot_layout(design = "
  ABCD
  EFGH
  IJKL
", guides = "collect") & theme(legend.position = "bottom")
  
  grid::grid.draw(grid::textGrob('Time (months)', x =.5, y=0.07))
  ggsave(file = "SP3_abundance_lines_genus.svg", width = 40, height = 30, units = "cm")
  ########
}



## C. Perform statistical modelling on composition (GLMM)
# Create a data frame with read counts for all samples and additional sample information
tmp <- as.data.frame(phyloseq_big %>%
                       tax_fix() %>%
                       tax_transform(trans = "identity", rank = rank_viz) %>%
                       otu_get(taxa = taxa, counts = TRUE)) %>%
  mutate(timepoint = factor(str_remove(phyloseq_big@sam_data$tp,"m"), levels = c("4","5","6","9","11","14")),
         id = str_sub(phyloseq_big@sam_data$kindcode, 1,4),
         CLASS_2 = phyloseq_big@sam_data$CLASS_2, #in reference to class 1
         CLASS_3 = phyloseq_big@sam_data$CLASS_3, #in reference to class 1
         delivery_home = ifelse(phyloseq_big@sam_data$bqq13_delivery_place == 1, yes = 1, no = 0), #in reference to hospital
         csection = ifelse(phyloseq_big@sam_data$bqq14_delivery_type == 4, yes = 1, no = 0),
         duration_breastfeeding = scale(as.numeric(as.character(microbial_age_df$duration_breastfeeding))),
         ageintrosolids = scale(as.numeric(as.character(phyloseq_big@sam_data$ageintrosolids))),
         birth_weight = scale(as.numeric(as.character(phyloseq_big@sam_data$bqq1_birthweight))),
         female = ifelse(phyloseq_big@sam_data$bqq3_gender == 2, yes = 1, no = 0), # in reference to male
         diet_class = factor(phyloseq_big@sam_data$diet_class))

# Create GLMMMs for each taxa
tmp_list <- list() #
i <- 1 #
for(ta in taxa) {
  print(paste(i,ta))
  tmp_ta <- tmp %>% dplyr::rename(taxa = ta)
  # Fit model:
  glmm <- glmmTMB(taxa ~ (CLASS_2 + CLASS_3)*timepoint + duration_breastfeeding + ageintrosolids + birth_weight + female + delivery_home + csection + (1 | id), data = tmp_ta, family = nbinom2)
  sm <- summary(glmm)
  if(is.na(sm$AICtab[1])) {
    print("Model convergence problem")
  }
  # Print results:
  print(round(sm$coefficients$cond[which(sm$coefficients$cond[,4] < 0.05), c(1,4)], digits=3))
  
  # Add results to a list for late vizualization of coefficients heatmap
  tmp_list[[i]] <- data.frame(round(sm$coefficients$cond,digits = 3),
                              taxa = ta,
                              sign = ifelse(sm$coefficients$cond[,4] < 0.05, ifelse(sm$coefficients$cond[,4] < 0.001, "***", "*"),""),
                              variable = c("Intercept", "Class 2","Class 3", "5m","6m","9m","11m","14m", "duration of breastfeeding",
                                           "ageintrosolids", "birth weigth", "female", "delivery at home", "C-section",
                                           "5m*Class 2", "6m*Class 2", "9m*Class 2","11m*Class 2", "14m*Class 2",
                                           "5m*Class 3","6m*Class 3","9m*Class 3","11m*Class 3","14m*Class 3"))
  i <- i + 1
}
# Extract results for chosen taxa (Genus coeff heatmap) (manually selected models without problems)
df <- rbind(tmp_list[[1]], tmp_list[[3]], tmp_list[[4]], tmp_list[[5]],
            tmp_list[[6]], tmp_list[[8]], tmp_list[[9]], tmp_list[[10]], tmp_list[[12]],
            tmp_list[[13]], tmp_list[[14]], tmp_list[[15]],
            tmp_list[[19]], tmp_list[[20]], tmp_list[[21]], tmp_list[[22]])
# Extract results for chosen taxa (Species coeff heatmap) (manualy selected models without problems)
df <- rbind(tmp_list[[1]], tmp_list[[2]], tmp_list[[3]], tmp_list[[4]], tmp_list[[5]],
            tmp_list[[6]], tmp_list[[7]], tmp_list[[8]], tmp_list[[9]], tmp_list[[11]], tmp_list[[12]],
            tmp_list[[13]], tmp_list[[14]], tmp_list[[15]], tmp_list[[16]], tmp_list[[17]],
            tmp_list[[18]],  tmp_list[[20]], tmp_list[[21]], tmp_list[[24]], tmp_list[[25]],
            tmp_list[[26]], tmp_list[[28]], tmp_list[[32]], tmp_list[[33]],
            tmp_list[[36]], tmp_list[[37]], tmp_list[[38]], tmp_list[[42]], tmp_list[[43]])
#nbinom2  19?, 23?,    27?     39? 40?  41?
# "Model convergence problem": 10, 22,29,  30, 34, 44
#
df$taxa <- gsub("_"," ", df$taxa)

# Visualize coefficients heatmap
# Function coef_heatmp at the end of the script !
df <- df %>% filter(variable != "Intercept") %>% # Filter out intercept
  dplyr::arrange(taxa)

if(rank_viz == "Species") {
  coef_heatmap(df = df, x = "variable", y = "taxa", coef = "Estimate", abs_coef_max = 30) +
    scale_x_discrete(limits = c("14m*Class 3",  "11m*Class 3",  "9m*Class 3",  "6m*Class 3",  "5m*Class 3",  "14m*Class 2",  "11m*Class 2",  "9m*Class 2",  "6m*Class 2",  "5m*Class 2","Class 3",  "Class 2",
                                "female", "duration of breastfeeding", "delivery at home", "C-section", "birth weigth", "ageintrosolids",
                                "14m",  "11m",  "9m",  "6m",  "5m"),
                     labels = c("14m*Class 3",  "11m*Class 3",  "9m*Class 3",  "6m*Class 3",  "5m*Class 3",  "14m*Class 2",  "11m*Class 2",  "9m*Class 2",  "6m*Class 2",  "5m*Class 2","Class 3",  "Class 2",
                                "Female", "Duration of breastfeeding",  "Delivery at home", "C-section", "Birth weigth",  "Age solids introduced",
                                "14m",  "11m",  "9m",  "6m",  "5m")) +
    scale_y_discrete(limits = unique(df$taxa),
                     position = "right") +
    coord_flip() +
    theme_minimal() +
    labs(x = "", y = "") +
    theme( panel.grid = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = -0.1)
    ) +
    geom_text(aes(label = sign), vjust = 0.7) +
    geom_vline(xintercept = 12.5, linetype = "solid") +
    geom_vline(xintercept = 18.5, linetype = "solid")  +
    theme(axis.text.x = element_text(face=c("bold","bold","plain", "plain", "bold", "bold","bold",
                                            "plain","plain","bold","plain","plain","plain","plain",
                                            "bold", "plain", "plain","plain","plain", "bold", "bold",
                                            "plain", "plain", "bold", "bold", "bold", "plain", "bold", "bold", "plain")))
  ggsave(filename = "P3A_GLMMspecies.svg", width = 30, height = 20, units = "cm")
}

if (rank_viz == "Genus") {
  coef_heatmap(df = df, x = "variable", y = "taxa", coef = "Estimate", abs_coef_max = 30) +
    scale_x_discrete(limits = c("14m*Class 3",  "11m*Class 3",  "9m*Class 3",  "6m*Class 3",  "5m*Class 3",  "14m*Class 2",  "11m*Class 2",  "9m*Class 2",  "6m*Class 2",  "5m*Class 2","Class 3",  "Class 2",
                                "female", "duration of breastfeeding", "delivery at home", "C-section", "birth weigth", "ageintrosolids",
                                "14m",  "11m",  "9m",  "6m",  "5m"),
                     labels = c("14m*Class 3",  "11m*Class 3",  "9m*Class 3",  "6m*Class 3",  "5m*Class 3",  "14m*Class 2",  "11m*Class 2",  "9m*Class 2",  "6m*Class 2",  "5m*Class 2","Class 3",  "Class 2",
                                "Female", "Duration of breastfeeding",  "Delivery at home", "C-section", "Birth weigth",  "Age solids introduced",
                                "14m",  "11m",  "9m",  "6m",  "5m")) +
    scale_y_discrete(limits = unique(df$taxa),
                     position = "right") +
    coord_flip() +
    theme_minimal() +
    labs(x = "", y = "") +
    theme( panel.grid = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = -0.1)
    ) +
    geom_text(aes(label = sign), vjust = 0.7) +
    geom_vline(xintercept = 12.5, linetype = "solid") +
    geom_vline(xintercept = 18.5, linetype = "solid")  +
    theme(axis.text.x = element_text(face=c("bold","bold","bold", "plain", "bold", "plain","bold","bold",
                                            "plain","plain","bold","bold","bold","plain", "bold", "bold")))
  ggsave(filename = "SP3_GLMMgenus.svg", width = 30, height = 20, units = "cm")
}
#dev.off()

write.csv(df %>% select(taxa, variable, Estimate, Std..Error, z.value, Pr...z.., sign),
          file = paste0("GLMM_",rank_viz, "_stats.csv"))





### 7. PERMANOVA
# What variables is the overall microbial community variation associated with?
# Assessing differences in community compositionwith permutational Multivariate Analysis of Variance, or perMANOVA.
# These tests are done on distances, meaning that they assess the differences between communities based on dissimilarity.

# -> For phyloseq object (longitudinal):
# permanova <- phyloseq_big_tmp %>%
#   ps_mutate(
#     duration_breastfeeding = microbial_age_df$duration_breastfeeding,
#     ageintrosolids = as.numeric(as.character(ageintrosolids)),
#     birth_weight = as.numeric(as.character(phyloseq_big_tmp@sam_data$bqq1_birthweight)),
#     pregn_weeks = as.numeric(as.character(bqq2_pregn_weeks))
#   ) %>%
#   tax_agg("Species") %>%
#   dist_calc("aitchison") %>%
#   dist_permanova(
#     variables = c("diet_class", "tp", "bqq3_gender", "bqq13_delivery_place", "bqq14_delivery_type", "duration of breastfeeding", "ageintrosolids",
#                   "birth_weight", "pregn_weeks", "fur_pets", "bqq45_siblings"),
#     n_perms = 9999,
#     seed = 12345, complete_cases = TRUE, verbose = "max"
#   )
# 
# 
# ggplot(data = data.frame(vars = c("dietary class", "time", "gender", "delivery place", "delivery type", "duration of breastfeeding", "ageintrosolids",
#                                   "birth weight", "pregnancy weeks", "furry pets", "siblings"), r2 = permanova@permanova$R2[1:11]),
#        aes(x = vars, y = r2)) +
#   geom_bar(stat = "identity") + 
#   coord_flip()



# -> PERMANOVA for each timepoint individually
permanova_list <- list(length = 6)
timepoints_names <- c("4m","5m","6m","9m","11m","14m")
for(i in 1:6) {
  print(timepoints_names[i])
  permanova_list[[i]] <- phyloseq_big %>%
          ps_mutate(
                    ageintrosolids = as.numeric(as.character(ageintrosolids)),
                    birth_weight = as.numeric(as.character(phyloseq_big@sam_data$bqq1_birthweight)),
                    pregn_weeks = as.numeric(as.character(bqq2_pregn_weeks)),
                    csection = ifelse(phyloseq_big@sam_data$bqq14_delivery_type == 4, yes = 1, no = 0)
          ) %>%
          ps_filter(tp == timepoints_names[i]) %>%
          tax_agg("Species") %>%
          dist_calc("aitchison") %>%
          dist_permanova(
            variables = c("diet_class", "bqq3_gender", "bqq13_delivery_place", "csection", "breastfeeding","formulafeeding" ,"ageintrosolids",
                          "birth_weight", "pregn_weeks", "furrypet_indoors_6m", "bqq45_siblings"),
            n_perms = 9999,
            seed = 12345, complete_cases = TRUE, verbose = "max"
          )
  print(permanova_list[[i]])
}

# Create a data frame with permanova results
df_permanova <- rbind(data.frame(vars = c("Dietary class", "Sex", "Delivery place", "Delivery type", "Breastfeeding", "Formulafeeding", "Age solids introduced",
                                          "Birth weight", "Gestational age", "Furry pets", "Older siblings"), r2 = permanova_list[[1]]@permanova$R2[1:11], tp = "4m", p = ifelse(permanova_list[[1]]@permanova$`Pr(>F)`[1:11] < 0.05,"*","" )),
                      data.frame(vars = c("Dietary class", "Sex", "Delivery place", "Delivery type", "Breastfeeding", "Formulafeeding", "Age solids introduced",
                                          "Birth weight", "Gestational age", "Furry pets", "Older siblings"), r2 = permanova_list[[2]]@permanova$R2[1:11], tp = "5m", p = ifelse(permanova_list[[2]]@permanova$`Pr(>F)`[1:11] < 0.05,"*","" )),
                      data.frame(vars = c("Dietary class", "Sex", "Delivery place", "Delivery type", "Breastfeeding", "Formulafeeding", "Age solids introduced",
                                          "Birth weight", "Gestational age", "Furry pets", "Older siblings"), r2 = permanova_list[[3]]@permanova$R2[1:11], tp = "6m", p = ifelse(permanova_list[[3]]@permanova$`Pr(>F)`[1:11] < 0.05,"*","" )),
                      data.frame(vars = c("Dietary class", "Sex", "Delivery place", "Delivery type", "Breastfeeding", "Formulafeeding", "Age solids introduced",
                                          "Birth weight", "Gestational age", "Furry pets", "Older siblings"), r2 = permanova_list[[4]]@permanova$R2[1:11], tp = "9m", p = ifelse(permanova_list[[4]]@permanova$`Pr(>F)`[1:11] < 0.05,"*","" )),
                      data.frame(vars = c("Dietary class", "Sex", "Delivery place", "Delivery type", "Breastfeeding", "Formulafeeding", "Age solids introduced",
                                          "Birth weight", "Gestational age", "Furry pets", "Older siblings"), r2 = permanova_list[[5]]@permanova$R2[1:11], tp = "11m", p = ifelse(permanova_list[[5]]@permanova$`Pr(>F)`[1:11] < 0.05,"*","" )),
                      data.frame(vars = c("Dietary class", "Sex", "Delivery place", "Delivery type", "Breastfeeding", "Formulafeeding", "Age solids introduced",
                                          "Birth weight", "Gestational age", "Furry pets", "Older siblings"), r2 = permanova_list[[6]]@permanova$R2[1:11], tp = "14m", p = ifelse(permanova_list[[6]]@permanova$`Pr(>F)`[1:11] < 0.05,"*","" )))
df_permanova$tp <- factor(df_permanova$tp, levels = timepoints_names, ordered = T)



# Visualize as bar plots
ggplot(data = df_permanova, aes(x = vars, y = r2)) +
  geom_bar(stat = "identity", fill = "darkgray") + 
  geom_text(aes(label=p), 
            size = 4,
            position = position_stack(vjust = 1.1),
            inherit.aes = TRUE) +
  scale_x_discrete(limits = c("Sex","Older siblings", "Gestational age", "Furry pets", "Formulafeeding", "Dietary class",
                              "Delivery type", "Delivery place", "Breastfeeding", "Birth weight", "Age solids introduced")) +
  facet_wrap(~tp, ) + 
  coord_flip() +
  theme_minimal() +
  labs(y = bquote(r^2), x = "", caption = "PERMANOVA on individual timepoints, species level")



# Save results
# write.csv(as.data.frame(permanova_list[[6]]@permanova[c(7,8,5,3,4,1,6,10,9,11,2),]) %>%
#             mutate(SumOfSqs = round(SumOfSqs, digits = 1),
#                    R2 = round(R2, digits = 3),
#                    `Pr(>F)` = round(`Pr(>F)`, digits = 3),
#                    F = round(F, digits = 3)),
#           file = "14m.csv")








### FUNCTIONS:

calcStats <- function(data, variable, groups) {
  colnames(data)[which(colnames(data) == variable)] <- "variable"
  colnames(data)[which(colnames(data) == groups)] <- "groups"
  for(t in timepoints_names) {
    response <- (data %>% filter(month == t))[, "variable"]
    factor <- (data %>% filter(month == t))[, "groups"]
    if(bartlett.test(response, factor)$p.value >0.05) { # if >0.05 than groups have the same variance
      if(shapiro.test(response)$p.value >0.05) {#if >0.05 than  normal  (perform if bartlett.test >0.05
        print("pairwise_t_test")
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
    print("pairwise_wilcox_test")
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
    geom_tile(aes(fill = .data[[coef]], height = .5, width = .9)) +
    scale_fill_distiller(
      type = "div", oob = oob_fun,
      limits = c(-abs_coef_max, abs_coef_max)
    ) +
    theme_bw() +
    theme(panel.grid = element_line(linewidth = 0.1))
  
  return(p)
}








