# LOAD LIBRARIES ====
library(tidyverse)
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)

# read in data ====
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

#read in phyloseq
my_phyloseq <-  readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# clean age ==== 

# calculate aitchinson distance ====
aitchison_dists <- my_phyloseq %>%  tax_transform("identity") %>% dist_calc("aitchison")

# PCoA ====
pcoa_ord <- ord_calc(aitchison_dists, method = "PCoA")
ord_plot(pcoa_ord, colour = "Age_individual",shape = "diet_class")+scale_color_viridis(discrete = TRUE)+
  labs(
    colour = "Age individual",   # Adjust legend title for colour
    shape = "Diet class"         # Adjust legend title for shape
  )

#change % manually here
plot <- ord_plot(pcoa_ord, colour = "Age_individual", shape = "diet_class") +
  scale_color_viridis(discrete = TRUE) +
  labs(
    x = "PCo1 [13.8%]",  # Change x-axis label
    y = "PCo2 [4.7%]",  # Change y-axis label
    colour = "Age individual",   # Adjust legend title for colour
    shape = "Diet class"         # Adjust legend title for shape
  )

# saving figure for the paper ====
ggsave(
  filename   = "S6_PCoA_diet_age.tiff",
  plot       = plot,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)
