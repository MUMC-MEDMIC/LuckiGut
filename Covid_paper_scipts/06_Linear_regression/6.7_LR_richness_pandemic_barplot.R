# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S6.6_LR_pandemic_rich_scaled_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

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
                                     "Collected pandemic (ref: before)",
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
plot_without_legend <- ggplot(data=final, aes(x=variables, y=Estimate, fill=Exposure_categories)) +
  geom_bar(stat = "identity") + coord_flip() + 
  ylab("Expected change in observed richness") +  
  facet_grid(~Age_individual) +
  geom_text(aes(label = sig_padj), colour = "black", size = 3,position = position_fill(vjust = 0),hjust = -0.1) + 
  ggtitle(NULL) + xlab(NULL) +
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" = "red3",
                               "Parental" = "#0072B2",
                               "Perinatal" = "darkgreen")) +
  theme_bw() +
  labs(fill = "Exposure \\ncategories") +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7), 
        axis.title.y = element_text(size = 1),
        legend.position = "none",  # Add this line
        strip.text = element_text(size = 7),
        panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
        panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)),
        panel.spacing = unit(0.15, "cm")
  )     # Reduce vertical spacing)

plot_without_legend

# saving figures for the paper ====
ggsave(
  filename   = "S6.7_LR_pandemic_rich_scaled_without_legend.svg",
  plot       = plot_without_legend,
  width    = 170,
  height   = 60,
  units    = "mm",
  device   = svglite
)





