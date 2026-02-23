# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S5.2_PERMANOVA_results_cleaned.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

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
model$run <- as.numeric(model$run)

# plotting down flow of PERMANOVA ====
# 1 week
model %>% filter(age == "1 week") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 1 week")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()
#4 weeks
model %>% filter(age == "4 weeks") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 4 weeks")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#8 weeks
model %>% filter(age == "8 weeks") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 8 weeks")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#4 months
model %>% filter(age == "4 months") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 4 months")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#5 months
model %>% filter(age == "5 months") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 5 months")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#6 months
model %>% filter(age == "6 months") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 6 months")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#9 months
model %>% filter(age == "9 months") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 9 months")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#11 months
model %>% filter(age == "11 months") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 11 months")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

#14 months
model %>% filter(age == "14 months") %>%  ggplot(aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+ facet_grid(~run)+
  geom_text(aes(label = sig_adj), colour = "black",size = 5) + ggtitle("PERMANOVA downflow for 14 months")+
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="darkgreen",
                               "Infant's health"="#CD0000",
                               "Parental"="#0074B6",
                               "Perinatal"="#CC79A7"))+theme_bw()

# preparing for final model ====
final <- model %>% filter(last_run == "yes")

# saving data ====
#write.table (final, "S5.3_PERMANOVA_final.txt", row.names=FALSE,sep = "\t")

# Barplot final model ====
plot_with_legend <- ggplot(data=final, aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") +facet_grid(~age)+ggtitle(NULL) + xlab(NULL) +
  geom_text(aes(label = sig_padj), colour = "black",size = 3) +
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Perinatal"="darkgreen"))+theme_bw()+
  labs(fill = "Exposure \ncategories")+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        strip.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),        # Shrink key boxes
        legend.key.height = unit(0.2, "cm"),     # Make thinner
        legend.key.width = unit(0.2, "cm"),      # Make narrower
        legend.spacing.y = unit(0.1, "cm"),
        panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
        panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))      # Reduce vertical spacing)

plot_without_legend <- ggplot(data=final, aes(x=variables, y=R2,fill=Exposure_categories)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") +facet_grid(~age)+ggtitle(NULL) + xlab(NULL) +
  geom_text(aes(label = sig_padj), colour = "black",size = 3) +
  scale_fill_manual(values = c("Dietary" = "orange",
                               "Environmental" ="red3",
                               "Parental"="#0072B2",
                               "Perinatal"="darkgreen"))+theme_bw()+
  labs(fill = "Exposure \ncategories")+
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        strip.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),        # Shrink key boxes
        legend.key.height = unit(0.2, "cm"),     # Make thinner
        legend.key.width = unit(0.2, "cm"),      # Make narrower
        legend.spacing.y = unit(0.1, "cm"),
        panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
        panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)),
        panel.spacing = unit(0.15, "cm"),
        legend.position = "none")      # Reduce vertical spacing)

# saving figures for the paper ====
ggsave(
  filename   = "S5.3_mult_PERMANOVA_pandemic_final_with_legend.svg",
  plot       = plot_with_legend,
  width    = 170,
  height   = 60,
  units    = "mm",
  device   = svglite
)

ggsave(
  filename   = "S5.3_mult_PERMANOVA_pandemic_final_without_legend.svg",
  plot       = plot_without_legend,
  width    = 162,
  height   = 60,
  units    = "mm",
  device   = svglite
)

# barplot for poster only covid =====
final_narrow <- final %>% filter(!(age %in% c("1 week", "4 weeks", "8 weeks", "4 months"))) %>% 
  filter(wrok_var == "covid_270220")

final_narrow %>% 
  ggplot(aes(x=variables, y=R2,fill=age)) +
  geom_bar(stat = "identity")+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab (" ")+facet_grid(~age)+
  geom_text(aes(label = sig_padj), colour = "black",size = 5) + ggtitle("Permanova. Effect size of pandemic")+
  theme_bw()+scale_fill_manual(values = c("5 months" = "#003664",
                                          "6 months" ="#2975A7",
                                          "9 months"="#3EB2BD",
                                          "11 months"="#DAE0E4",
                                          "14 months"="#CCC1B8"))+
  theme(legend.position = "none")

# barplot for poster last TP ====
final_lastTP <- final %>% filter(!(age %in% c("1 week", "4 weeks", "8 weeks", "4 months"))) 
ggplot(data=final_lastTP, aes(x=variables, y=R2,fill=Theme)) +
  geom_bar(stat = "identity")+coord_flip()+scale_y_continuous(labels=scales::percent)+ 
  ylab("% explained (R^2)") + xlab ("Variable")+facet_grid(~age)+
  geom_text(aes(label = sig_padj), colour = "black",size = 5) + ggtitle("Core final models")+
  scale_fill_manual(values = c("Diet" = "orange",
                               "Environmental exposures" ="darkgreen",
                               "Infant health"="#CD0000",
                               "Parental health"="#0072B2",
                               "Perinatal covariates"="#CC79A7"))+theme_bw()

# circular barplot ====
plot1<- final %>% filter(age == "1 week") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("1-2 weeks") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot2<- final %>% filter(age == "4 weeks") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("4 weeks") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot3<- final %>% filter(age == "8 weeks") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("8 weeks") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot4<- final %>% filter(age == "4 months") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("4 months") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot5<- final %>% filter(age == "5 months") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("5 months") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot6<- final %>% filter(age == "6 months") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("6 months") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot7<- final %>% filter(age == "9 months") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("9 months") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot8<- final %>% filter(age == "11 months") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("11 months") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

plot9<- final %>% filter(age == "14 months") %>% ggplot(aes(x = variables, y = R2, fill = Exposure_categories)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = sig_padj), colour = "black", size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Dietary" = "orange",
    "Environmental" = "red3",
    "Parental" = "#0072B2",
    "Perinatal" = "darkgreen"
  )) +
  theme_minimal() +
  ggtitle("14 months") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 8)
  )+coord_radial(inner.radius = 0.01)

(plot1 + plot2 + plot3 + 
    plot4 + plot5 + plot6 +
    plot7 + plot8 + plot9) + 
  plot_layout(ncol = 3, nrow = 3)

plots <- list(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9)
plots <- list(plot1,plot4,plot7,plot2,plot5,plot8,plot3,plot6,plot9)
# Remove margins for all plots to minimize space
plots <- lapply(plots, function(p) p + theme(plot.margin = unit(c(0.1,0.9,0.1,1.3), "cm")))
# Arrange in 3x3 grid
all <- wrap_plots(plots, ncol = 3, nrow = 3)

# saving circular plot ====
ggsave(
  filename   = "S5.3_mult_PERMANOVA_pandemic_final_circle.tiff",
  plot       = all,
  width      = 170,      # full page width
  height     = 90,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

# saving figures for the paper ====
#ggsave("S5.3_mult_PERMANOVA_pandemic_final_circle.pdf", all, device = "pdf",width = 600, height = 300, dpi = 300, units = "mm")









