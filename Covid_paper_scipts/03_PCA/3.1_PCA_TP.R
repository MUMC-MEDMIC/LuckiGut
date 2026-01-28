# load libraries ====
#PCA / PCoA with vectors, beta diversity
#Centered Log-Ratio (clr) Transformation is used 
#transforms your data by taking the log of the ratio between observed frequencies x and their geometric mean G(x),because
#ord_calc method was used aka PCA, the taxa loading vectors can be drawn. PCA combines taxa abundances into new dimensions. The first axes display the greatest variation in your microbial data
#vectors contribute the most to the differences
library(viridis)
library(vegan)
library(phyloseq)
library(microViz)
library(ggplot2)
library(stringr)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in the data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# PCA ====
tax_table(infant)[,"Species"] <- gsub("_", " ", tax_table(infant)[,"Species"])

plot1 <- infant %>% 
  tax_transform("clr", rank = "Species") %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age_individual", 
           plot_taxa = 1:5, 
           tax_lab_style = tax_lab_style(
             size = 5, 
             colour = "black",
             fontface = "italic"),
           auto_caption = NA)+
  scale_color_viridis(discrete = TRUE)+
  labs(color = "Infants age")+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.position = "left")

plot1

# saving figures for the paper ====
ggsave("S3.1_PCA_TP.tiff", plot1, width = 400, height = 300, dpi = 300, units = "mm")
#ggsave("S3.1_PCA_TP.pdf", plot1, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")








