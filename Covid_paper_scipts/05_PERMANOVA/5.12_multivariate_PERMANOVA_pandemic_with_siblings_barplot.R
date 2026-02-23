# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S5.11_PERMANOVA_results_cleaned_with_siblings.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# Cleaning data for barplot ====
model$age<- factor(model$age, levels = c("1 week",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))

model$variables[model$variables == "Collected during pandemic"] <- "Pandemic"
model$variables[model$variables == "Complementary foods age"] <- "Complementary age"

model$variables <- factor(model$variables, 
                          levels = c("Complementary age",
                                     "Breastfeeding",
                                     "Formula feeding",
                                     "Probiotics",
                                     "Daycare",
                                     "Pets pregnancy",
                                     "Older siblings",
                                     "Pets",
                                     "Siblings in cohort",
                                     "Antibiotics",
                                     "Pandemic",
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

# Barplot final model ====
plot <- ggplot(data=model, aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("")+facet_grid(~age)+
  geom_text(aes(label = sig_padj), colour = "black",size = 4) + ggtitle("")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Perinatal"="darkgreen"))+theme_bw()+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        legend.position = "none")

# saving figures for the paper ====
ggsave(
  filename   = "S5.12_PERMANOVA_final_with_siblings.tiff",
  plot       = plot,
  width      = 170,      # full page width
  height     = 90,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)
