# libraries ====
library(microViz)
library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(colorRamp2)
library(svglite)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in phyloseq data ====
baku_f_I <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
entero = read.delim("Input/10_DMM_cluster_per_sample_from_MS.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
baku = read.delim("Input/11_DMM_top30_driving_bacterial_species_MS.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# read in enterotype data ====
entero$Sample_aliquote <- rownames(entero)
entero <- entero %>% mutate(DMM_Cluster = recode(DMM_Cluster, "5" = 5, "3" = 2, "6" = 3, "4" = 6, "1" = 4, "2" = 1))
count_entero <- table(entero$DMM_Cluster)

baku_f_I <- ps_join(baku_f_I, entero, by = "Sample_aliquote")

# extracting mean reads per enterotype ====
df_DMM <- plyr::rbind.fill(data.frame(t(rowMeans((baku_f_I %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1) %>%
                                                    tax_transform("compositional", rank = "Species"))@otu_table))), 
                           data.frame(t(rowMeans((baku_f_I %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2) %>%
                                                    tax_transform("compositional", rank = "Species"))@otu_table))),
                           data.frame(t(rowMeans((baku_f_I %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3) %>%
                                                    tax_transform("compositional", rank = "Species"))@otu_table))),
                           data.frame(t(rowMeans((baku_f_I %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4) %>%
                                                    tax_transform("compositional", rank = "Species"))@otu_table))),
                           data.frame(t(rowMeans((baku_f_I %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5) %>%
                                                    tax_transform("compositional", rank = "Species"))@otu_table))),
                           data.frame(t(rowMeans((baku_f_I %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6) %>%
                                                    tax_transform("compositional", rank = "Species"))@otu_table))))
df_DMM[is.na(df_DMM)] <- 0

# selecting top 30 most abundant species ====
baku$bacteria <- rownames(baku)
bacteria_names <- unique(baku$bacteria)
bacteria_names30 <- bacteria_names[1:30]
df_DMM_30 <- df_DMM %>% select(all_of(bacteria_names30))
df_for_heatmap <- as.matrix(t(df_DMM_30))
table <- as.data.frame(as.matrix(t(df_DMM_30)))
table$bacteria <- row.names(table)

# heatmap ====
# Remove underscores from row names
row_names <- gsub("_", " ", rownames(df_for_heatmap))

ht <- Heatmap(df_for_heatmap,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(df_for_heatmap[i, j], digits = 3), x, y, gp = gpar(fontsize = 8))
        },
        column_names_gp = gpar(cex = 1),
        column_names_rot = 0,
        row_labels = row_names,
        col = colorRamp2(c(0, 0.4), c("#ffffff", "#005824")),
        column_labels = paste0("n=", as.vector(count_entero)),
        heatmap_legend_param = list(title = "Mean RA"),
        top_annotation = HeatmapAnnotation(
          bar = anno_simple(1:6, col = c("1" = "#796248", "2" = "#a47852", "3" = "#cd9c58", 
                                         "4" = "#e6cd98", "5" = "#f8edcf", "6" = "#d4ede8")),
          text = anno_text(1:6,rot = 0, gp = gpar(fontsize = 10, col = "black"), just = "right"),
          show_annotation_name = FALSE,
          show_legend = FALSE
        )
)

ht

# saving figure for the paper ====
tiff("S8.1_DMM_bacteria_per_enterotype.tiff", width = 11, height = 6, units = "in", res = 300, compression = "lzw")

draw(ht)

dev.off()


