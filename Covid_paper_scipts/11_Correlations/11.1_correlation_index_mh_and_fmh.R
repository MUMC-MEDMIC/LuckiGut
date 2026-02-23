# LOAD LIBRARIES ====
library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
index_f_m_h <- read.delim("Output/Files/S1.7_index_fmh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
index_m_h <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# correlating index fmh and mh ====
index_m_h$Family_ID <- rownames(index_m_h)
index_f_m_h$Family_ID <- rownames(index_f_m_h)
test <- full_join(index_m_h,index_f_m_h, by = "Family_ID")

correlation <- cor.test(test$index_f_m_h, test$index_m_h, method = "spearman", exact = FALSE) # p - 1.376e-15, rho 0.9403346 
rho <- round(correlation$estimate,3)
p_val <- correlation$p.value  

# plotting correlation ====
plot <- ggplot(test, aes(x=index_m_h, y=index_f_m_h)) + geom_point()+ labs(title = "")+
  ylab("Exposure index with fathers") + xlab("Exposure index without fathers") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15))+
  annotate("text", x = Inf, y = Inf, hjust = 4.1, vjust = 1.1, 
           label = paste("\nrho =", rho, "\np =", "< 0.001"),   
           size = 4, fontface = "bold")+
  ggtitle("Spearman correlation")

plot

# saving figures for the paper ====
#ggsave("S11.1_index_spearman_correlation.pdf", plot, device = "pdf",width = 300, height = 200, dpi = 300, units = "mm")

#save both figures with dummy legend and without, later copy-paste legend instead of the dummy figure
ggsave(
  filename   = "S11.1_index_spearman_correlation.tiff",
  plot       = plot,
  width      = 170,      # full page width
  height     = 90,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

