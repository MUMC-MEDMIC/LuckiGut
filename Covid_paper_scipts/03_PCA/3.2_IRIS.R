# load libraries====
library(microViz)
library(viridis)
library(stringr)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# iris plot ====
clr_pca <- infant %>%
  tax_transform("clr") %>% ord_calc(method = "PCA")

iris <- clr_pca %>%
  ord_plot_iris(
    tax_level = "Species", n_taxa = 12,
    tax_lab_style = tax_lab_style(type = "text", size = 2.5, fontface = "bold.italic"),
    anno_colour = "Age_individual",
    taxon_renamer = function(x) str_replace_all(stringr::str_remove_all(x, "^G: | [ae]t rel."), "_", " ")) +
  scale_color_viridis(discrete = TRUE)+
  guides(colour = "none")+
  theme(legend.text = element_text(size = 13))

iris

# saving figures for the paper ====
ggsave("S3.2_iris_TP.tiff", iris, width = 400, height = 300, dpi = 300, units = "mm")
#ggsave("S3.2_iris_TP.pdf", iris, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")

