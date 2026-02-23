# load libraries ====
library(tidyverse)
library(ggplot2)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
index = read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data =====
index$Family_ID <- rownames(index)

# histogram ====
plot <- index %>% ggplot(aes(x = index_m_h)) + 
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_bw() +
  xlab("Index values") + ylab("Frequency") +
  coord_cartesian(xlim = c(3, 20)) +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

ggsave(
  filename   = "S10.9_index_mh_histogram.tiff",
  plot       = plot,
  width      = 85,      # full page width
  height     = 60,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

