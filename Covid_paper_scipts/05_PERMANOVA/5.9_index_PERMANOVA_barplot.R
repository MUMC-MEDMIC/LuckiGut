# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S5.8_index_PERMANOVA_results_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# clean data ====
model$age<- factor(model$age, levels = c("5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))

model$variables <- factor(model$variables, 
                          levels = c("Complementary foods age",
                                     "Breastfeeding",
                                     "Formula feeding",
                                     "Probiotics",
                                     "Daycare",
                                     "Exposure index",
                                     "Pets pregnancy",
                                     "Older siblings",
                                     "Pets",
                                     "Siblings in cohort",
                                     "Antibiotics",
                                     "Collected during pandemic",
                                     "Mother's atopy",
                                     "Father's atopy",
                                     "Antepartum antibiotics",
                                     "Intrapartum antibiotics",
                                     "Maternal weight gain",
                                     "Prenatal smoking",
                                     "Birth weight",
                                     "Delivery type",
                                     "Delivery place",
                                     "Hospitalisation",
                                     "Gestational age",
                                     "Sex"))
model$run <- as.numeric(model$run)

# preparing for final model ====
final <- model %>% filter(last_run == "yes")

# Barplot final model ====
plot <- ggplot(data=final, aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("")+facet_grid(~age)+
  geom_text(aes(label = sig_padj), colour = "black",size = 5) + ggtitle("")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Infant's health" = "#CC79A7",
                               "Perinatal"="darkgreen"))+theme_bw()+
  labs(fill = "Exposure \ncategories")+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 11),
        legend.position = "none")

# saving figures for the paper ====
ggsave(
  filename   = "S5.9_index_PERMANOVA_final.tiff",
  plot       = plot,
  width      = 225,      # full page width
  height     = 160,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)





