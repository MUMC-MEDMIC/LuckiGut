# load libraries====
library(microViz)
library(viridis)
library(stringr)
library(svglite)
library(ggplot2)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# iris plot ====
clr_pca <- infant %>%
  tax_transform("clr") %>% ord_calc(method = "PCA")

iris_with_legend <- clr_pca %>%
  ord_plot_iris(
    tax_level = "Species", n_taxa = 12,
    tax_lab_style = tax_lab_style(type = "text", size = 2.5, fontface = "bold.italic"),
    anno_colour = "Age_individual",
    taxon_renamer = function(x)
      str_replace_all(stringr::str_remove_all(x, "^G: | [ae]t rel."), "_", " ")
  ) +
  scale_color_viridis(discrete = TRUE, name = NULL) +   # remove legend title
  theme(
    legend.position = "right",
    legend.title = element_blank(),                     # extra-safe
    legend.text  = element_text(size = 4, face = "italic"),
    legend.key.height = unit(2, "mm"),
    legend.key.width  = unit(2, "mm"),
    legend.spacing.x  = unit(1, "mm"),
    legend.spacing.y  = unit(1, "mm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1)))

iris_without_legend <- clr_pca %>%
  ord_plot_iris(
    tax_level = "Species", n_taxa = 12,
    tax_lab_style = tax_lab_style(type = "text", size = 2.5, fontface = "bold.italic"),
    anno_colour = "Age_individual",
    taxon_renamer = function(x)
      str_replace_all(stringr::str_remove_all(x, "^G: | [ae]t rel."), "_", " ")
  ) +
  scale_color_viridis(discrete = TRUE, name = NULL) +   # remove legend title
  theme(
    legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 1)))

# saving figures for the paper ====
ggsave(
  filename   = "S3.2_iris_TP_with_legend.svg",
  plot       = iris_with_legend,
  width    = 100,
  height   = 60,
  units    = "mm",
  device   = svglite
)

ggsave(
  filename   = "S3.2_iris_TP_without_legend.svg",
  plot       = iris_without_legend,
  width    = 54.67,
  height   = 60,
  units    = "mm",
  device   = svglite
)

