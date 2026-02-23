# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model_without <- read.delim("Output/Files/S6.22_LR_pandemic_rich_scaled_cleaned_without_outliers.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
model_with <- read.delim("Output/Files/S6.6_LR_pandemic_rich_scaled_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# Cleaning data for barplot ====
model_with <- model_with %>% filter(Age_individual == "4 months" | Age_individual == "5 months")
model_with$outliers <- "with outliers"

model_without <- model_without %>% filter(Age_individual == "4 months" | Age_individual == "5 months")
model_without$outliers <- "without outliers"

model_with <- model_with %>% mutate(Age_individual_new = paste(Age_individual, outliers, sep = "_"))
model_without <- model_without %>% mutate(Age_individual_new = paste(Age_individual, outliers, sep = "_"))

model <- rbind(model_with,model_without)

plyr::count(model$Age_individual_new)
model$Age_individual_new<- factor(model$Age_individual_new, levels = c("4 months_with outliers",
                                                                       "4 months_without outliers",
                                                                       "5 months_with outliers",
                                                                       "5 months_without outliers"))

model$variables <- factor(model$variables, 
                          levels = c("Compelmentary foods age (wks)",
                                     "Breastfeeding (ref: no)",
                                     "Formulafeeding (ref: no)",
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
plot <- ggplot(data=final, aes(x=variables, y=Estimate,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+ 
  ylab("Expected change (coefficient)") + xlab ("Variable")+facet_grid(~Age_individual_new)+
  geom_text(aes(label = sig_padj), colour = "black",size = 4) + ggtitle("Observed richness")+ xlab("")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Perinatal"="darkgreen"))+theme_bw()+
  labs(fill = "Exposure \ncategories")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15))

plot

# saving figures for the paper ====
ggsave("S6.23_LR_pandemic_rich_scaled_without_outliers.tiff", plot, width = 400, height = 300, dpi = 300, units = "mm")
#ggsave("S6.23_LR_pandemic_rich_scaled_without_outliers.pdf", plot, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")







