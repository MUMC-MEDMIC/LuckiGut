# LOAD LIBRARIES ====
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(patchwork)
library(grid)
library(svglite)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
spearman_all = read.delim("Output/Files/S4.1_baku_vs_fiber_per_TP_spearman_adjust.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# clean data ====
spearm_4m <- spearman_all %>% filter (Age_individual == "4 months")
spearm_5m <- spearman_all %>% filter (Age_individual == "5 months")
spearm_6m <- spearman_all %>% filter (Age_individual == "6 months")
spearm_9m <- spearman_all %>% filter (Age_individual == "9 months")
spearm_11m <- spearman_all %>% filter (Age_individual == "11 months")
spearm_14m <- spearman_all %>% filter (Age_individual == "14 months")
rm(spearman_all)

# creating matrices for heatmaps ====
#separating p values and correlations into 2 df for heatmap
corr_4m <- dcast(spearm_4m, bacteria ~ fiber, value.var = "Correlation_coefficient_r")
rownames(corr_4m) <- corr_4m$bacteria
corr_4m$bacteria <- NULL 
pvalue_4m <- dcast(spearm_4m, bacteria ~ fiber, value.var = "p_value_adj") #reshaping data, baku will be row names, fiber, column names, filled with adjusted p-values
rownames(pvalue_4m) <- pvalue_4m$bacteria
pvalue_4m$bacteria <- NULL 

corr_5m <- dcast(spearm_5m, bacteria ~ fiber, value.var = "Correlation_coefficient_r")
rownames(corr_5m) <- corr_5m$bacteria
corr_5m$bacteria <- NULL 
pvalue_5m <- dcast(spearm_5m, bacteria ~ fiber, value.var = "p_value_adj")
rownames(pvalue_5m) <- pvalue_5m$bacteria
pvalue_5m$bacteria <- NULL 

corr_6m <- dcast(spearm_6m, bacteria ~ fiber, value.var = "Correlation_coefficient_r")
rownames(corr_6m) <- corr_6m$bacteria
corr_6m$bacteria <- NULL 
pvalue_6m <- dcast(spearm_6m, bacteria ~ fiber, value.var = "p_value_adj")
rownames(pvalue_6m) <- pvalue_6m$bacteria
pvalue_6m$bacteria <- NULL 

corr_9m <- dcast(spearm_9m, bacteria ~ fiber, value.var = "Correlation_coefficient_r")
rownames(corr_9m) <- corr_9m$bacteria
corr_9m$bacteria <- NULL 
pvalue_9m <- dcast(spearm_9m, bacteria ~ fiber, value.var = "p_value_adj")
rownames(pvalue_9m) <- pvalue_9m$bacteria
pvalue_9m$bacteria <- NULL 

corr_11m <- dcast(spearm_11m, bacteria ~ fiber, value.var = "Correlation_coefficient_r")
rownames(corr_11m) <- corr_11m$bacteria
corr_11m$bacteria <- NULL 
pvalue_11m <- dcast(spearm_11m, bacteria ~ fiber, value.var = "p_value_adj")
rownames(pvalue_11m) <- pvalue_11m$bacteria
pvalue_11m$bacteria <- NULL

corr_14m <- dcast(spearm_14m, bacteria ~ fiber, value.var = "Correlation_coefficient_r")
rownames(corr_14m) <- corr_14m$bacteria
corr_14m$bacteria <- NULL 
pvalue_14m <- dcast(spearm_14m, bacteria ~ fiber, value.var = "p_value_adj")
rownames(pvalue_14m) <- pvalue_14m$bacteria
pvalue_14m$bacteria <- NULL
rm(spearm_4m,spearm_5m,spearm_6m,spearm_9m,spearm_11m,spearm_14m)


#saving results
tp4 <- corr_4m
tp4$bacteria <- rownames(tp4)
tp4$age <- '4'
tp5 <- corr_5m
tp5$bacteria <- rownames(tp5)
tp5$age <- '5'
tp6 <- corr_6m
tp6$bacteria <- rownames(tp6)
tp6$age <- '6'
tp9 <- corr_9m
tp9$bacteria <- rownames(tp9)
tp9$age <- '9'
tp11 <- corr_11m
tp11$bacteria <- rownames(tp11)
tp11$age <- '11'
tp14 <- corr_14m
tp14$bacteria <- rownames(tp14)
tp14$age <- '14'
correlations <- rbind(tp4,tp5,tp6,tp9,tp11,tp14)

#write.table (correlations, "S3_18_baku_fiber_correlations.txt", row.names=TRUE,sep = "\t")
rm(tp4,tp5,tp6,tp9,tp11,tp14)

tp4 <- pvalue_4m
tp4$bacteria <- rownames(tp4)
tp4$age <- '4'
tp5 <- pvalue_5m
tp5$bacteria <- rownames(tp5)
tp5$age <- '5'
tp6 <- pvalue_6m
tp6$bacteria <- rownames(tp6)
tp6$age <- '6'
tp9 <- pvalue_9m
tp9$bacteria <- rownames(tp9)
tp9$age <- '9'
tp11 <- pvalue_11m
tp11$bacteria <- rownames(tp11)
tp11$age <- '11'
tp14 <- pvalue_14m
tp14$bacteria <- rownames(tp14)
tp14$age <- '14'
pvalues <- rbind(tp4,tp5,tp6,tp9,tp11,tp14)

#write.table (pvalues, "S3_19_baku_fiber_correlations_pvalues.txt", row.names=TRUE,sep = "\t")
rm(tp4,tp5,tp6,tp9,tp11,tp14)

# Visualize correlations as heatmap
corr_4m <- as.matrix(corr_4m)
corr_5m <- as.matrix(corr_5m)
corr_6m <- as.matrix(corr_6m)
corr_9m <- as.matrix(corr_9m)
corr_11m <- as.matrix(corr_11m)
corr_14m <- as.matrix(corr_14m)

#replacing adjusted p-values with *
pvalue_4m <- pvalue_4m %>% mutate(across(everything(), ~ ifelse(. <= 0.05, "*", ifelse(is.na(.), "", ""))))
pvalue_5m <- pvalue_5m %>% mutate(across(everything(), ~ ifelse(. <= 0.05, "*", ifelse(is.na(.), "", ""))))
pvalue_6m <- pvalue_6m %>% mutate(across(everything(), ~ ifelse(. <= 0.05, "*", ifelse(is.na(.), "", ""))))
pvalue_9m <- pvalue_9m %>% mutate(across(everything(), ~ ifelse(. <= 0.05, "*", ifelse(is.na(.), "", ""))))
pvalue_11m <- pvalue_11m %>% mutate(across(everything(), ~ ifelse(. <= 0.05, "*", ifelse(is.na(.), "", ""))))
pvalue_14m <- pvalue_14m %>% mutate(across(everything(), ~ ifelse(. <= 0.05, "*", ifelse(is.na(.), "", ""))))

# Supplementary figures 11 ====
my_palette <- colorRampPalette(c("#0072B2", "white", "#CD0000"))(100)
panel_4m <- Heatmap(corr_4m,cluster_columns = FALSE,cluster_rows = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(pvalue_4m[i, j],x, y, gp = gpar(fontsize = 8))
                    },
                    row_labels = gsub("_", " ", rownames(corr_4m)),
                    col = my_palette,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1)), column_title = "4 months") 

#TO MAKE CLUSTER cluster_columns = TRUE,cluster_rows = TRUE
panel_5m <- Heatmap(corr_5m,cluster_columns = FALSE,cluster_rows = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      value <- pvalue_5m[i, j]
                      if (!is.na(value)) {
                        grid.text(value, x, y, gp = gpar(fontsize = 8))
                      }
                    },
                    row_labels = gsub("_", " ", rownames(corr_5m)),
                    col = my_palette,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1)), column_title = "5 months")

panel_6m <- Heatmap(corr_6m,cluster_columns = FALSE,cluster_rows = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      value <- pvalue_6m[i, j]
                      if (!is.na(value)) {
                        grid.text(value, x, y, gp = gpar(fontsize = 8))
                      }
                    },
                    row_labels = gsub("_", " ", rownames(corr_6m)),
                    col = my_palette,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1)), column_title = "6 months")

panel_9m <- Heatmap(corr_9m,cluster_columns = FALSE,cluster_rows = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(pvalue_9m[i, j],x, y, gp = gpar(fontsize = 8))
                    },
                    row_labels = gsub("_", " ", rownames(corr_9m)),
                    col = my_palette,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1)), column_title = "9 months")

panel_11m <- Heatmap(corr_11m,cluster_columns = FALSE,cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(pvalue_11m[i, j],x, y, gp = gpar(fontsize = 8))
                     },
                     row_labels = gsub("_", " ", rownames(corr_11m)),
                     col = my_palette,
                     heatmap_legend_param = list(
                       title = "Correlation",
                       at = c(-1, 0, 1)), column_title = "11 months")

panel_14m <- Heatmap(corr_14m,cluster_columns = FALSE,cluster_rows = FALSE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(pvalue_14m[i, j],x, y, gp = gpar(fontsize = 8))
                     },
                     row_labels = gsub("_", " ", rownames(corr_14m)),
                     col = my_palette,
                     heatmap_legend_param = list(
                       title = "Correlation",
                       at = c(-1, 0, 1)), column_title = "14 months")

# all bacteria that had sig correlations with fibers
bacteria_fiber_abundant0.01 <- unique(c(rownames(corr_4m),rownames(corr_5m),rownames(corr_6m),rownames(corr_9m),rownames(corr_11m),rownames(corr_14m))) #n=31

panel_4m+panel_5m+panel_6m
panel_9m+panel_11m+panel_14m

# clustering panels ====
panel_4m_with <- Heatmap(corr_4m,cluster_columns = TRUE,cluster_rows = TRUE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(pvalue_4m[i, j],x, y, gp = gpar(fontsize = 4))
                    },
                    row_dend_width = unit(0.5, "cm"), 
                    row_labels = gsub("_", " ", rownames(corr_4m)),
                    row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                    column_names_gp = gpar(fontsize = 4),
                    column_dend_height = unit(0.5, "cm"),
                    col = my_palette,
                    show_heatmap_legend = TRUE,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1),
                      legend_height = unit(1.5, "cm"),
                      grid_width = unit(0.25, "cm"),
                      title_gp = gpar(fontsize = 5),
                      labels_gp = gpar(fontsize = 4)), column_title = "4 months", 
                    column_title_gp = gpar(fontsize = 6)) 

panel_4m_without <- Heatmap(corr_4m,cluster_columns = TRUE,cluster_rows = TRUE,
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           grid.text(pvalue_4m[i, j],x, y, gp = gpar(fontsize = 4))
                         },
                         row_dend_width = unit(0.5, "cm"), 
                         row_labels = gsub("_", " ", rownames(corr_4m)),
                         row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                         column_names_gp = gpar(fontsize = 4),
                         column_dend_height = unit(0.5, "cm"),
                         col = my_palette,
                         show_heatmap_legend = FALSE,
                         column_title = "4 months", 
                         column_title_gp = gpar(fontsize = 6)) 

panel_5m <- Heatmap(corr_5m,cluster_columns = TRUE,cluster_rows = TRUE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      value <- pvalue_5m[i, j]
                      if (!is.na(value)) {
                        grid.text(value, x, y, gp = gpar(fontsize = 4))
                      }
                    },
                    row_labels = gsub("_", " ", rownames(corr_5m)),
                    row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                    column_names_gp = gpar(fontsize = 4),
                    column_dend_height = unit(0.5, "cm"),
                    col = my_palette,
                    show_heatmap_legend = FALSE,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1),
                      title_gp = gpar(fontsize = 5),   # Title font size
                      labels_gp = gpar(fontsize = 5)),
                    column_title = "5 months", 
                    column_title_gp = gpar(fontsize = 6)) 

panel_6m <- Heatmap(corr_6m,cluster_columns = TRUE,cluster_rows = TRUE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      value <- pvalue_6m[i, j]
                      if (!is.na(value)) {
                        grid.text(value, x, y, gp = gpar(fontsize = 4))
                      }
                    },
                    row_labels = gsub("_", " ", rownames(corr_6m)),
                    row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                    column_names_gp = gpar(fontsize = 4),
                    column_dend_height = unit(0.5, "cm"),
                    col = my_palette,
                    show_heatmap_legend = FALSE,
                    heatmap_legend_param = list(
                      title = "Correlation",
                      at = c(-1, 0, 1)), column_title = "6 months", 
                    column_title_gp = gpar(fontsize = 6)) 

panel_9m <- Heatmap(corr_9m,cluster_columns = TRUE,cluster_rows = TRUE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(pvalue_9m[i, j],x, y, gp = gpar(fontsize = 4))
                    },
                    row_dend_width = unit(0.5, "cm"), 
                    row_labels = gsub("_", " ", rownames(corr_9m)),
                    row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                    column_names_gp = gpar(fontsize = 4),
                    column_dend_height = unit(0.5, "cm"),
                    col = my_palette,
                    show_heatmap_legend = FALSE,
                    column_title = "9 months", 
                    column_title_gp = gpar(fontsize = 6)) 

panel_11m <- Heatmap(corr_11m,cluster_columns = TRUE,cluster_rows = TRUE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(pvalue_11m[i, j],x, y, gp = gpar(fontsize = 4))
                     },
                     row_labels = gsub("_", " ", rownames(corr_11m)),
                     row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                     column_names_gp = gpar(fontsize = 4),
                     column_dend_height = unit(0.5, "cm"),
                     col = my_palette,
                     show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(
                       title = "Correlation",
                       at = c(-1, 0, 1),
                     title_gp = gpar(fontsize = 6),   # Title font size
                     labels_gp = gpar(fontsize = 6)),
                     column_title = "11 months", 
                     column_title_gp = gpar(fontsize = 6)) 

panel_14m <- Heatmap(corr_14m,cluster_columns = TRUE,cluster_rows = TRUE,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(pvalue_14m[i, j],x, y, gp = gpar(fontsize = 4))
                     },
                     row_labels = gsub("_", " ", rownames(corr_14m)),
                     row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                     column_names_gp = gpar(fontsize = 4),
                     column_dend_height = unit(0.5, "cm"),
                     col = my_palette,
                     show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(
                       title = "Correlation",
                       at = c(-1, 0, 1)), column_title = "14 months", 
                     column_title_gp = gpar(fontsize = 6)) 

# saving figure for the paper ====
w_cm <- 28/2
h_cm <- 17/2

svglite(
  "S4.2_456_months_IFDP_vs_baku_branching_with.svg",
  width  = w_cm / 2.54,   # svglite uses inches
  height = h_cm / 2.54
)

draw(panel_4m_with+panel_5m+panel_6m)  # important for ComplexHeatmap
dev.off()

svglite("S4.2_456_months_IFDP_vs_baku_branching_without.svg",
        width  = w_cm / 2.54,   # svglite uses inches
        height = h_cm / 2.54)

draw(panel_4m_without+panel_5m+panel_6m)  # important for ComplexHeatmap
dev.off()

svglite("S4.2_91114_months_IFDP_vs_baku_branching.svg",
        width  = w_cm / 2.54,   # svglite uses inches
        height = h_cm / 2.54)

draw(panel_9m+panel_11m+panel_14m)  # important for ComplexHeatmap
dev.off()


