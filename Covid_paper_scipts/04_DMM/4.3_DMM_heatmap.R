# load libraries ====
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpubr)
library(tidyr)
library(phyloseq)
library(microViz)
library("ComplexHeatmap")

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ==== 
cluster <- read.delim("Output/Files/S4.2_DMM_cluster_species_table_renamed.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
bakut <- read.delim("Output/Files/S4.1_DMM_all_driving_bacterial_species.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data ====
cluster$MMHP_SampleID <- cluster$SampleID
cluster <- cluster %>% select(c(DMM_Cluster_recoded,MMHP_SampleID))

infant_rel <- infant %>% tax_transform(trans = "compositional")

#sanity checks 
#mat <- as(phyloseq::otu_table(infant_rel), "matrix")
#if (phyloseq::taxa_are_rows(infant_rel)) mat <- t(mat)

#range(mat, na.rm = TRUE)          # should be 0..1
#summary(rowSums(mat, na.rm=TRUE)) # should be ~1
#rm(mat)

bakut <- bakut %>% filter (Mean > 0.01)
bakut$bacteria <- rownames(bakut)
baku <- unique(bakut$bacteria)

# clean otu reads ====
otu_rel <- as.data.frame(as.matrix(infant_rel@otu_table))
otu_rel_t <- as.data.frame(as.matrix(t(otu_rel)))
otu_rel_t$MMHP_SampleID <- rownames(otu_rel_t)
otu_t_cluster <- full_join(cluster,otu_rel_t, by = "MMHP_SampleID")
otu_t_cluster$MMHP_SampleID <- NULL
rm(infant,infant_rel,cluster,otu_rel,otu_rel_t)

# calculating average ====
average <- otu_t_cluster %>%
  group_by(DMM_Cluster_recoded) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

average <- as.data.frame(as.matrix(average))
rownames(average) <- average$DMM_Cluster_recoded
average$DMM_Cluster_recoded <- NULL
average_t <- as.data.frame(as.matrix(t(average)))
average_t$bacteria <- rownames(average_t)

rm(average,otu_t_cluster)

# writing data ====
#write.table (average_t, "S4.3_DMM_all_driving_bacterial_species_clusters_renamed.txt", row.names=FALSE,sep = "\t")

# making long df ====
long_heatmap <- pivot_longer(data = average_t, cols = c("1","2","3","4","5","6","7"),names_to = "cluster", values_to = "rel_abundance")
long_heatmap_abun <- long_heatmap %>% filter(bacteria %in% baku)
long_heatmap_abun$rel_abundance <- as.numeric(long_heatmap_abun$rel_abundance)

long_heatmap_abun$bacteria <- gsub("_", " ", long_heatmap_abun$bacteria)

# heatmap without branches ====
plot <- ggplot(long_heatmap_abun, aes(cluster, bacteria, fill= rel_abundance)) + theme_bw()+
  geom_tile(color="black")+ggtitle("")+scale_fill_gradient(low="#EFE9F0ff", high="#632B71ff")+
  xlab("DMM cluster")+ylab("")+labs(fill = "Relative\nabundance")+theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15), 
    axis.title.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13))

plot

# preparing for branching ====
daa <- average_t %>% filter(bacteria %in% baku) %>% select(!(bacteria)) %>% as.matrix()
max_value <- round(max(daa),1)
min_value <- round(min(daa),1)

# heatmap with branches ====
my_palette <- colorRampPalette(c("#EFE9F0ff", "#632B71ff"))(100)
heatmap_with_legend <- Heatmap(matrix = daa,
                         cluster_columns = FALSE,
                         cluster_rows = TRUE,
                         row_labels = gsub("_", " ", rownames(daa)),
                         column_labels = gsub("_", " ", colnames(daa)),
                         column_names_rot = 0,
                         column_names_centered = TRUE, 
                         column_names_gp = gpar(fontsize = 9),
                         column_title = "DMM clusters",
                         column_title_side = "bottom", 
                         column_title_gp = gpar(fontsize = 6),
                         row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                         row_dend_width = unit(4, "cm"),  
                         col = my_palette,
                         heatmap_legend_param = list(
                           title = "Mean \nrelative \nabundance",
                           at = c(min_value,max_value),
                           labels_gp = gpar(fontsize = 6),
                           title_gp = gpar(fontsize = 6)))

heatmap_without_legend <- Heatmap(
  matrix = daa,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_labels = gsub("_", " ", colnames(daa)),
  column_title_gp = gpar(fontsize = 0.1),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 6, fontface = "italic"),
  row_labels = gsub("_", " ", rownames(daa)),
  row_dend_width = unit(1, "cm"),
  row_title_gp = gpar(fontsize = 0.1),
  col = my_palette,
  show_heatmap_legend = FALSE,
  column_title = ""
)

# saving figures for the paper ====
# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 11.4
h_cm <- 6

svglite("S4.3_DMM_heatmap_with_legend.svg",
        width  = w_cm / 2.54,
        height = h_cm / 2.54)

draw(heatmap_with_legend)  # important for ComplexHeatmap
dev.off()

# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 10
h_cm <- 6

svglite("S4.3_DMM_heatmap_without_legend.svg",
        width  = w_cm / 2.54,
        height = h_cm / 2.54)

draw(heatmap_without_legend)  # important for ComplexHeatmap
dev.off()









