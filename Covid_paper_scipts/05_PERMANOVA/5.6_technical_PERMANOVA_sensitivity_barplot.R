# Load packages ====
library(tidyverse)
library(ggplot2)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S5.5_tech_PERMANOVA_results_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning for the barplot ====
model$age<- factor(model$age, levels = c("1 week",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "14 months"))

model$variables <- factor(model$variables, 
                          levels = c("Complementary foods age",
                                     "Breastfeeding",
                                     "Formula feeding",
                                     "Probiotics",
                                     "Daycare",
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
                                     "Sex",
                                     "Sequencing batch"))
model$run <- as.numeric(model$run)

# preparing for final model ====
final <- model %>% filter(last_run == "yes")

# Barplot final model ====
plot <- ggplot(data=final, aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("")+facet_grid(~age)+
  geom_text(aes(label = sig_padj), colour = "black",size = 3) + ggtitle("")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Perinatal"="darkgreen",
                               "Technical variable"="#b9b5b5"))+theme_bw()+
  labs(fill = "Exposure \ncategories")+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        strip.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),        # Shrink key boxes
        legend.key.height = unit(0.2, "cm"),     # Make thinner
        legend.key.width = unit(0.2, "cm"),      # Make narrower
        legend.spacing.y = unit(0.1, "cm"))

plot

# saving figures for the paper ====
ggsave(
  filename   = "S5.6_tech_PERMANOVA_final.tiff",
  plot       = plot,
  width      = 170,      # full page width
  height     = 90,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)


