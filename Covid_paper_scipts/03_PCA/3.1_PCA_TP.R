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
library(patchwork)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in the data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

infant@sam_data$Age_individual <- gsub("weeks", "wks", infant@sam_data$Age_individual)
infant@sam_data$Age_individual <- gsub("months", "mos", infant@sam_data$Age_individual)
infant@sam_data$Age_individual <- gsub("1-2", "1", infant@sam_data$Age_individual)

infant@sam_data$Age_individual<- factor(infant@sam_data$Age_individual, levels = c("1 wks", "4 wks", "8 wks",
                                                                                   "4 mos", "5 mos", "6 mos", "9 mos", "11 mos", "14 mos"))
# PCA ====
tax_table(infant)[,"Species"] <- gsub("_", " ", tax_table(infant)[,"Species"])

plot_with_legend <- infant %>% 
  tax_transform("clr", rank = "Species") %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age_individual", 
           plot_taxa = 1:5, 
           size = 0.8,
           stroke = 0,
           tax_vec_style_sel = vec_tax_sel(linewidth = 0.2),
           tax_lab_style = tax_lab_style(
             size = 1.9, 
             colour = "black",
             fontface = "italic",
             angle = 10),
           auto_caption = NA)+
  scale_color_viridis(discrete = TRUE)+
  labs(color = "Infants \nage")+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "mm"),
        legend.key.width  = unit(0.5, "mm"),
        legend.spacing.y  = unit(2.1, "mm"),   # vertical space between items
        legend.spacing.x  = unit(0.1, "mm"),
        legend.margin     = margin(0, 0, 0, 0, "mm"),
        legend.box.margin = margin(0, 0, 0, 0, "mm"))

plot_without_legend <- infant %>% 
  tax_transform("clr", rank = "Species") %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age_individual", 
           plot_taxa = 1:5, 
           size = 0.8,
           stroke = 0,
           tax_vec_length = 1.2,      # Shortens arrow length (default ~1)
           tax_lab_length = 1.5,   
           tax_vec_style_sel = vec_tax_sel(linewidth = 0.2),
           tax_lab_style = tax_lab_style(
             size = 1.9, 
             colour = "black",
             fontface = "italic",
             angle = 10),
           auto_caption = NA)+
  scale_color_viridis(discrete = TRUE)+
  labs(color = "Infants \nage")+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2))

# saving figures for the paper ====
ggsave(
  filename   = "S3.1_PCA_TP_with_legend.svg",
  plot       = plot_with_legend,
  width    = 109,
  height   = 60,
  units    = "mm",
  device   = svglite
)

ggsave(
  filename   = "S3.1_PCA_TP_without_legend.svg",
  plot       = plot_without_legend,
  width    = 90,
  height   = 55,
  units    = "mm",
  device   = svglite
)











