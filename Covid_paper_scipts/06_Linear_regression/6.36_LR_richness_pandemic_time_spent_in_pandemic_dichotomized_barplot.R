# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S6.35_LR_pandemic_rich_scaled_time_spent_in_pandemic_dichotomized_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# Cleaning data for barplot ====
model$Age_individual<- factor(model$Age_individual, levels = c("1 week",
                                                               "4 weeks",
                                                               "8 weeks",
                                                               "4 months",
                                                               "5 months",
                                                               "6 months",
                                                               "9 months",
                                                               "11 months",
                                                               "14 months"))

model$variables <- factor(model$variables, 
                          levels = c("Complementary foods age (wks)",
                                     "Breastfeeding (ref: no)",
                                     "Formula feeding (ref: no)",
                                     "Probiotics (ref: no)",
                                     "Daycare (ref: no)",
                                     "Pets pregnancy (ref: no)",
                                     "Older siblings (ref: no)",
                                     "Pets (ref: no)",
                                     "Siblings in cohort (ref: no)",
                                     "Antibiotics (ref: no)",
                                     "Time spent in pandemic (ref: below median)",
                                     "Mother's atopy (ref: no)",
                                     "Father's atopy (ref: no)",
                                     "Antepartum antibiotics (ref: no)",
                                     "Intrapartum antibiotics (ref: no)",
                                     "Maternal weight gain (kg)",
                                     "Prenatal smoking (ref: no)",
                                     "Birth weight (g)",
                                     "Vaginal delivery (ref: C-section)",
                                     "Delivery: hospital (ref: home)",
                                     "Hospitalisation (ref: no)",
                                     "Gestational age (wks)",
                                     "Male (ref. female)"))

model$run <- as.numeric(model$run)

# preparing for final model ====
final <- model %>% filter(last_run == "yes") %>% filter(Term != "(Intercept)")

# Barplot final model ====
plot <- ggplot(data=final, aes(x=variables, y=Estimate,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+ 
  ylab("Expected change in ENS") + xlab ("Variable")+facet_grid(~Age_individual)+
  geom_text(aes(label = sig_padj), colour = "black",size = 2) + ggtitle("")+ xlab("")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Infant's health" = "#CC79A7",
                               "Perinatal"="darkgreen"))+theme_bw()+
  labs(fill = "Exposure \ncategories")+
  theme(axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 7), 
        axis.title.y = element_text(size = 1),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        strip.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),        # Shrink key boxes
        legend.key.height = unit(0.2, "cm"),     # Make thinner
        legend.key.width = unit(0.2, "cm"),      # Make narrower
        legend.spacing.y = unit(0.1, "cm"),
        panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
        panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))      # Reduce vertical spacing)
plot 

# saving figures for the paper ====
ggsave(
  filename   = "S6.36_LR_richness_scaled_dichotomized.tiff",
  plot       = plot,
  width      = 170,      # full page width
  height     = 100,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)






