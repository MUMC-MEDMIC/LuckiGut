# Load in libraries ====
library("tidyverse")
library("ComplexHeatmap")
library(ggplot2)
library(tidyverse)
library(svglite)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
linda = read.delim("Output/Files/S7.1_LINDA_pandemic.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning for barplot ====
linda$Age_individual<- factor(linda$Age_individual, 
                              levels = c("1 week",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"
                              ))

linda_sig <- linda %>% subset(padj <= 0.2)
linda_sig_21 <- linda %>% subset(padj <= 0.05) %>% filter(variable == "covid_270220yes")
linda_sig_covid <- linda_sig %>% subset(variable == "covid_270220yes")
sig_bakut_covid <- unique(linda_sig_covid$bacteria) #should be 44
linda_sig_21 <- unique(linda_sig_21$bacteria)

linda_covid <- filter(linda, bacteria %in% sig_bakut_covid)
linda_covid <- linda_covid %>% subset(variable == "covid_270220yes")

# saving data ====
#write.table (linda_sig_covid, "S7.2_LINDA_pandemic_covid_subset.txt", row.names=FALSE,sep = "\t")

rm(linda_sig_covid,linda_sig)

#covid significant only p < 0.050
ggplot(linda_covid, aes(Age_individual, str_replace_all(bacteria, "_", " "), fill= log2FoldChange)) + 
  geom_tile(color="black")+ggtitle("DAA of species for each TP, based on PERMANOVA final models, covid added, padj (<0.2)")+scale_fill_gradient2(low="#0072B2", high="#CD0000", mid = "white")+
  geom_text(aes(label = sig_padj), colour = "black",size = 3)+ labs(y = "Bacteria", x = "Infant's age")+theme_bw()

# preparing for branching ====
processed <- filter(linda, bacteria %in% sig_bakut_covid)
processed <- processed %>% subset(variable == "covid_270220yes")

cor_1w <- processed %>% filter(Age_individual == "1 week")
cor_4w <- processed %>% filter(Age_individual == "4 weeks")
cor_8w <- processed %>% filter(Age_individual == "8 weeks")
cor_4m <- processed %>% filter(Age_individual == "4 months")
cor_5m <- processed %>% filter(Age_individual == "5 months")
cor_6m <- processed %>% filter(Age_individual == "6 months")
cor_9m <- processed %>% filter(Age_individual == "9 months")
cor_11m <- processed %>% filter(Age_individual == "11 months")
cor_14m <- processed %>% filter(Age_individual == "14 months")

rownames(cor_1w) <- cor_1w$bacteria
rownames(cor_4w) <- cor_4w$bacteria
rownames(cor_8w) <- cor_8w$bacteria
rownames(cor_4m) <- cor_4m$bacteria
rownames(cor_5m) <- cor_5m$bacteria
rownames(cor_6m) <- cor_6m$bacteria
rownames(cor_9m) <- cor_9m$bacteria
rownames(cor_11m) <- cor_11m$bacteria
rownames(cor_14m) <- cor_14m$bacteria

p0 <- cor_1w %>% select(sig_padj)
p1 <- cor_4w %>% select(sig_padj)
p2 <- cor_8w %>% select(sig_padj)
p4 <- cor_4m %>% select(sig_padj)
p5 <- cor_5m %>% select(sig_padj)
p6 <- cor_6m %>% select(sig_padj)
p9 <- cor_9m %>% select(sig_padj)
p11 <- cor_11m %>% select(sig_padj)
p14 <- cor_14m %>% select(sig_padj)

p0 <- p0 %>% dplyr::rename(`1_week` = sig_padj)
p1 <- p1 %>% dplyr::rename(`4_weeks` = sig_padj)
p2 <- p2 %>% dplyr::rename(`8_weeks` = sig_padj)
p4 <- p4 %>% dplyr::rename(`4_months` = sig_padj)
p5 <- p5 %>% dplyr::rename(`5_months` = sig_padj)
p6 <- p6 %>% dplyr::rename(`6_months` = sig_padj)
p9 <- p9 %>% dplyr::rename(`9_months` = sig_padj)
p11 <- p11 %>% dplyr::rename(`11_months` = sig_padj)
p14 <- p14 %>% dplyr::rename(`14_months` = sig_padj)

pvalues <- cbind(p0,p1,p2,p4,p5,p6,p9,p11,p14)
pvalues <- as.matrix(pvalues)
rm(p0,p1,p2,p4,p5,p6,p9,p11,p14)

cor_1w <- cor_1w %>% select(log2FoldChange)
cor_4w <- cor_4w %>% select(log2FoldChange)
cor_8w <- cor_8w %>% select(log2FoldChange)
cor_4m <- cor_4m %>% select(log2FoldChange)
cor_5m <- cor_5m %>% select(log2FoldChange)
cor_6m <- cor_6m %>% select(log2FoldChange)
cor_9m <- cor_9m %>% select(log2FoldChange)
cor_11m <- cor_11m %>% select(log2FoldChange)
cor_14m <- cor_14m %>% select(log2FoldChange)

cor_1w <- cor_1w %>% dplyr::rename(`1_week` = log2FoldChange)
cor_4w <- cor_4w %>% dplyr::rename(`4_weeks` = log2FoldChange)
cor_8w <- cor_8w %>% dplyr::rename(`8_weeks` = log2FoldChange)
cor_4m <- cor_4m %>% dplyr::rename(`4_months` = log2FoldChange)
cor_5m <- cor_5m %>% dplyr::rename(`5_months` = log2FoldChange)
cor_6m <- cor_6m %>% dplyr::rename(`6_months` = log2FoldChange)
cor_9m <- cor_9m %>% dplyr::rename(`9_months` = log2FoldChange)
cor_11m <- cor_11m %>% dplyr::rename(`11_months` = log2FoldChange)
cor_14m <- cor_14m %>% dplyr::rename(`14_months` = log2FoldChange)

daa <- cbind(cor_1w,cor_4w,cor_8w,cor_4m,cor_5m,cor_6m,cor_9m,cor_11m,cor_14m)
daa <- as.matrix(daa)

rm(cor_1w,cor_4w,cor_8w,cor_4m,cor_5m,cor_6m,cor_9m,cor_11m,cor_14m)

# heatmap with branches ====
my_palette <- colorRampPalette(c("#0072B2", "white", "#CD0000"))(100)
plot_with_legend <- Heatmap(
  matrix = daa,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- pvalues[i, j]
    if (!is.na(value)) {
      # move text to the visual center of the cell
      grid.text(
        value,
        x = x,
        y = y - height * 0.4,                      # base center
        gp = gpar(fontsize = 8),
        just = "centre"              # ensure centered both horizontally & vertically
      )
    }
  },
  row_labels = gsub("_", " ", rownames(daa)),
  row_names_gp = gpar(fontsize = 4.5, fontface = "italic"),
  column_labels = gsub("_", " ", colnames(pvalues)),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 7),
  col = my_palette,
  heatmap_legend_param = list(
    title = expression(log[2]*"FC"),
    at = c(-8, 0, 8),
    labels_gp = gpar(fontsize = 4),
    legend_height = unit(1, "cm"),
    legend_width = unit(0.5, "cm"),
    title_gp = gpar(fontsize = 6)
  ),
  column_title = ""
)

plot_with_legend

# saving figures for the paper ====
# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 17
h_cm <- 6.7

svglite("S7.2_Linda_heatmap_with_legend.svg",
        width  = w_cm / 2.54,
        height = h_cm / 2.54)

draw(plot_with_legend)  # important for ComplexHeatmap
dev.off()

plot_without_legend <- Heatmap(
  matrix = daa,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_heatmap_legend = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- pvalues[i, j]
    if (!is.na(value)) {
      grid.text(
        value,
        x = x,
        y = y - height * 0.4,  # decrease or increase 0.15 to adjust
        gp = gpar(fontsize = 8),
        just = "centre"         # ensure centered around (x, y)
      )
    }
  },
  row_labels = gsub("_", " ", rownames(daa)),
  row_names_gp = gpar(fontsize = 4.5, fontface = "italic"),
  column_labels = gsub("_", " ", colnames(pvalues)),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 7),
  col = my_palette,
  column_title = ""
)
plot_without_legend

# saving figures for the paper ====

# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 17
h_cm <- 6.7

svglite("S7.2_Linda_heatmap_without_legend.svg",
        width  = w_cm / 2.54,
        height = h_cm / 2.54)

draw(plot_without_legend)  # important for ComplexHeatmap
dev.off()

# volcano plot ====
ggplot(final, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.2)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed") + 
  theme_minimal()
colnames(final)
plyr::count(final$bacteria)

# MA Plot (Mean vs. Fold Change Plot) ====
ggplot(final, aes(x = baseMean, y = log2FoldChange, color = padj < 0.2)) +
  geom_point() +
  scale_x_log10() +
  theme_minimal()

# Bar Plot of Differential Abundant Taxa ==== 
vari <- c("#F5B8B6","#003856","#A1B4A5")

processed_poster <- processed_poster %>%
  mutate(bacteria = gsub("_", " ", bacteria))

processed_poster %>% 
  filter(padj < 0.05) %>% 
  ggplot(aes(x = reorder(bacteria, log2FoldChange), y = log2FoldChange, fill = Age_individual)) +
  geom_col() + 
  facet_grid(~Age_individual) + 
  coord_flip() + 
  scale_fill_manual(values = vari) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +  # Add horizontal line at y = 0
  theme_minimal()+
  labs(x = "", y = "Microbial Abundance Change")  +
  labs(fill = " ")

# preparing for heatmap poster ====
poster_cov <- linda %>% subset(variable == "covid_270220yes") #selecting only covid
poster_cov_age <- poster_cov %>% filter(Age_individual %in% c("5 months","6 months","9 months","11 months", "14 months"))
poster_cov_age_padj <- poster_cov_age %>% filter(padj < 0.05)

sig_bakut_poster <- unique(poster_cov_age_padj$bacteria) #should be 20

poster <- filter(linda, bacteria %in% sig_bakut_poster)
poster <- poster %>% subset(variable == "covid_270220yes")
poster <- poster %>% filter(Age_individual %in% c("5 months","6 months","9 months","11 months", "14 months"))

rm(poster_cov,poster_cov_age, poster_cov_age_padj)

cor_5m <- poster %>% filter(Age_individual == "5 months")
cor_6m <- poster %>% filter(Age_individual == "6 months")
cor_9m <- poster %>% filter(Age_individual == "9 months")
cor_11m <- poster %>% filter(Age_individual == "11 months")
cor_14m <- poster %>% filter(Age_individual == "14 months")

rownames(cor_5m) <- cor_5m$bacteria
rownames(cor_6m) <- cor_6m$bacteria
rownames(cor_9m) <- cor_9m$bacteria
rownames(cor_11m) <- cor_11m$bacteria
rownames(cor_14m) <- cor_14m$bacteria

p5 <- cor_5m %>% select(sig_padj)
p6 <- cor_6m %>% select(sig_padj)
p9 <- cor_9m %>% select(sig_padj)
p11 <- cor_11m %>% select(sig_padj)
p14 <- cor_14m %>% select(sig_padj)

p5 <- p5 %>% dplyr::rename(`5_months` = sig_padj)
p6 <- p6 %>% dplyr::rename(`6_months` = sig_padj)
p9 <- p9 %>% dplyr::rename(`9_months` = sig_padj)
p11 <- p11 %>% dplyr::rename(`11_months` = sig_padj)
p14 <- p14 %>% dplyr::rename(`14_months` = sig_padj)

pvalues <- cbind(p5,p6,p9,p11,p14)
pvalues <- as.matrix(pvalues)
rm(p5,p6,p9,p11,p14)

cor_5m <- cor_5m %>% select(log2FoldChange)
cor_6m <- cor_6m %>% select(log2FoldChange)
cor_9m <- cor_9m %>% select(log2FoldChange)
cor_11m <- cor_11m %>% select(log2FoldChange)
cor_14m <- cor_14m %>% select(log2FoldChange)

cor_5m <- cor_5m %>% dplyr::rename(`5_months` = log2FoldChange)
cor_6m <- cor_6m %>% dplyr::rename(`6_months` = log2FoldChange)
cor_9m <- cor_9m %>% dplyr::rename(`9_months` = log2FoldChange)
cor_11m <- cor_11m %>% dplyr::rename(`11_months` = log2FoldChange)
cor_14m <- cor_14m %>% dplyr::rename(`14_months` = log2FoldChange)

daa <- cbind(cor_5m,cor_6m,cor_9m,cor_11m,cor_14m)
daa <- as.matrix(daa)

rm(cor_5m,cor_6m,cor_9m,cor_11m,cor_14m)

# Heatmap poster ====
library(circlize)
my_palette <- colorRamp2(c(min(daa, na.rm = TRUE), 0, max(daa, na.rm = TRUE)), 
                         c("#003664", "white", "#962C2C")) 

Heatmap(matrix = daa,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- pvalues[i, j]
          if (!is.na(value)) {
            grid.text(value, x, y, gp = gpar(fontsize = 8))
          }
        },
        row_labels = gsub("_", " ", rownames(daa)),
        column_labels = gsub("_", " ", colnames(pvalues)),
        column_names_rot = 0,  # This sets column names to horizontal
        col = my_palette,
        heatmap_legend_param = list(
          title = "Log \nchange",
          at = c(-7, 0, 8) #check values from daa table to set these values
        ),
        column_title = "Abundace of bacteria in pandemic samples")

