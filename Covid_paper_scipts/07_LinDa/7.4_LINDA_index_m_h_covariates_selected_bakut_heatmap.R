# LOAD LIBRARIES ====
library(ggplot2)
library(tidyverse)
library("microViz")
library("ComplexHeatmap")
library("patchwork")

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
linda = read.delim("Output/Files/S7.3_LINDA_final_sig_model_index_selected_bakut.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# barplot preparing barplot ====
#sig bakut
linda_sig <- linda %>% subset(sig_padj != "")
linda_sig_covid <- linda_sig %>% subset(variable == "index_m_h")

sig_bakut_covid_padj <- unique(linda_sig_covid$bacteria)

linda_sig_padj <- filter(linda, bacteria %in% sig_bakut_covid_padj)
linda_sig_padj <- filter(linda_sig_padj, variable == "index_m_h")

ggplot(linda_sig_padj, aes(Age_individual, bacteria, fill= log2FoldChange)) + theme_bw()+
  geom_tile(color="black")+ggtitle("Index+covariates per TP sig p.adj bakut duo 0.2")+scale_fill_gradient2(low="#0047AB", high="#610023", mid = "white")+
  geom_text(aes(label = sig_padj), colour = "black",size = 3)

#all bakut 
linda_index <- linda %>% subset(variable == "index_m_h")
ggplot(linda_index, aes(Age_individual, bacteria, fill= log2FoldChange)) + theme_bw()+
  geom_tile(color="black")+ggtitle("Index+covariates per TP sig p.adj bakut duo 0.2")+scale_fill_gradient2(low="#0047AB", high="#610023", mid = "white")+
  geom_text(aes(label = sig_padj), colour = "black",size = 3)

# preparing for barplot with grouping ====
#sort correlations per TP
processed <- linda %>% subset(variable == "index_m_h")
cor_5m <- processed %>% filter(Age_individual == "5 months")
cor_6m <- processed %>% filter(Age_individual == "6 months")
cor_9m <- processed %>% filter(Age_individual == "9 months")
cor_11m <- processed %>% filter(Age_individual == "11 months")
cor_14m <- processed %>% filter(Age_individual == "14 months")

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

cor_5m <- cor_5m %>% dplyr::rename(`6_months` = log2FoldChange)
cor_6m <- cor_6m %>% dplyr::rename(`6_months` = log2FoldChange)
cor_9m <- cor_9m %>% dplyr::rename(`9_months` = log2FoldChange)
cor_11m <- cor_11m %>% dplyr::rename(`11_months` = log2FoldChange)
cor_14m <- cor_14m %>% dplyr::rename(`14_months` = log2FoldChange)

daa <- cbind(cor_5m,cor_6m,cor_9m,cor_11m,cor_14m)
daa <- as.matrix(daa)
rm(cor_5m,cor_6m,cor_9m,cor_11m,cor_14m)

# barplot with branches ====
my_palette <- colorRampPalette(c("#0039A6", "white", "#D52B1E"))(100)

heatmap_with_legend <- Heatmap(
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
  column_labels = gsub("_", " ", colnames(pvalues)),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 4.5, fontface = "italic"),
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

# saving figures for the paper ====
# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 8.5
h_cm <- 9

svglite("S7.4_Linda_index_heatmap_with_legend.svg",
        width  = w_cm / 2.54,
        height = h_cm / 2.54)

draw(heatmap_with_legend)  # important for ComplexHeatmap
dev.off()

heatmap_without_legend <- Heatmap(
  matrix = daa,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- pvalues[i, j]
    if (!is.na(value)) {
      grid.text(
        value,
        x = x,
        y = y - height * 0.4,
        gp = gpar(fontsize = 8),
        just = "centre"
      )
    }
  },
  row_labels = gsub("_", " ", rownames(daa)),
  column_labels = gsub("_", " ", colnames(pvalues)),
  column_names_rot = 45,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 4.5, fontface = "italic"),
  col = my_palette,
  show_heatmap_legend = FALSE,
  column_title = ""
)

# saving figures for the paper ====
# svg device sizes are in INCHES -> convert cm to inches
w_cm <- 8.5
h_cm <- 9

svglite("S7.4_Linda_index_heatmap_without_legend.svg",
        width  = w_cm / 2.54,
        height = h_cm / 2.54)

draw(heatmap_without_legend)  # important for ComplexHeatmap
dev.off()

# checking abundance for sig species ====
sig_baku <- linda_sig_covid$bacteria
otu <- as.data.frame(as.matrix(baku_f_I_covid@otu_table))
otu$bacteria <- rownames(otu)
otu <- otu %>% filter(bacteria %in% sig_baku)
otu <- as.data.frame(as.matrix(t(otu)))
otu$MMHP_SampleID <- rownames(otu)

meta <- as.data.frame(as.matrix(baku_f_I_covid@sam_data))
meta_otu <- full_join(otu,meta, by = "MMHP_SampleID")
meta_otu <- subset(meta_otu, !grepl("bacteria", MMHP_SampleID, ignore.case = TRUE))

age <- c("5 months","6 months", "9 months", "11 months", "14 months")
meta_otu <- meta_otu %>% filter(Age_individual %in% age)
meta_otu <- meta_otu[!is.na(meta_otu$index_m_h), ]
meta_otu$Gordonibacter_pamelaeae <- as.numeric(meta_otu$Gordonibacter_pamelaeae)
meta_otu$Blautia_obeum <- as.numeric(meta_otu$Blautia_obeum)
meta_otu$Prevotella_timonensis <- as.numeric(meta_otu$Prevotella_timonensis)
meta_otu$index_m_h <- as.numeric(meta_otu$index_m_h)

meta_otu$Age_individual<- factor(meta_otu$Age_individual, 
                                 levels = c("5 months",
                                            "6 months",
                                            "9 months",
                                            "11 months",
                                            "14 months"
                                 ))

rm(otu,meta,sig_baku)

# scatterplots ====
plot1 <- meta_otu %>% filter(Age_individual == "9 months") %>% 
  ggplot(aes(x=index_m_h, y=Gordonibacter_pamelaeae)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Gordonibacter pamelaeae at 9 months") +
  labs(x = "", 
       y = "Abundance") +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 14))

plot2 <- meta_otu %>% filter(Age_individual == "14 months") %>% 
  ggplot(aes(x=index_m_h, y=Gordonibacter_pamelaeae)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Gordonibacter pamelaeae at 14 months") +
  labs(x = "", 
       y = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 14))

plot3 <- meta_otu %>% filter(Age_individual == "9 months") %>% 
  ggplot(aes(x=index_m_h, y=Blautia_obeum)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Blautia obeum at 9 months") +
  labs(x = "Exposure index", 
       y = "Abundance") +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 14))

plot4 <- meta_otu %>% filter(Age_individual == "9 months") %>% 
  ggplot(aes(x=index_m_h, y=Prevotella_timonensis)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Prevotella timonensis at 9 months") +
  labs(x = "Exposure index", 
       y = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 14))

comb <- (plot1+plot2)
comb
