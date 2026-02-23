# load libraries ====
library(microViz)
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(viridis)
library(vegan)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
dmm <- read.delim("Output/Files/S4.2_DMM_cluster_species_table_renamed.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
  
# preparing data for PCA ====
dmm <- dmm %>% select(all_of(c("DMM_Cluster_recoded", "SampleID")))
dmm <- dmm %>% rename(MMHP_SampleID = SampleID)
dmm$DMM_Cluster_recoded <- factor(dmm$DMM_Cluster_recoded)

infant_merged <- ps_join(infant,dmm, by = "MMHP_SampleID")

# David's PCA ====
infant_merged %>%
  tax_transform("clr", rank = "Species") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "DMM_Cluster_recoded") +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "DMM clusters") +
  theme_bw()

# David's PCoA ====
infant_merged %>%
  tax_transform("identity", rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "DMM_Cluster_recoded") +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "DMM clusters") +
  theme_bw()

# PCA ====
set.seed(190) 

plot <- infant_merged %>% 
  tax_transform("clr", rank = "Species") %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "DMM_Cluster_recoded",auto_caption = NA,size = 0.5)+scale_color_viridis(discrete = TRUE)+
  labs(color = "DMM clusters")+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11), 
        axis.title.y = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  stat_ellipse(aes(group = DMM_Cluster_recoded, color = DMM_Cluster_recoded), 
               type = "t",    # T-distribution (95% confidence)
               level = 0.95,  # 95% confidence ellipses
               linetype = 2,  # Dashed lines
               size = 0.5)

plot

# saving figures for the paper ====
#ggsave("S4.4_DMM_PCA.pdf", plot, device = "pdf",width = 500, height = 300, dpi = 300, units = "mm")

ggsave(
  filename   = "S4.4_DMM_PCA.tiff",
  plot       = plot,
  width      = 160,      # full page width
  height     = 60,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

