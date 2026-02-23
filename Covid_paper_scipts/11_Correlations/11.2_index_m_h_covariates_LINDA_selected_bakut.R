# LOAD LIBRARIES ====
library("phyloseq")
library("tidyverse")
library("microViz")
library("patchwork")

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
linda_index = read.delim("Output/Files/S7.3_LINDA_final_sig_model_index_selected_bakut.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# checking abundance for sig species ====
linda_sig_covid <- linda_index %>% filter(variable == "index_m_h") %>% filter(padj <= 0.2)
sig_baku <- linda_sig_covid$bacteria
rm(linda_index,linda_sig_covid)

otu <- as.data.frame(as.matrix(infant@otu_table))
otu$bacteria <- rownames(otu)
otu <- otu %>% filter(bacteria %in% sig_baku)
otu <- as.data.frame(as.matrix(t(otu)))
otu$MMHP_SampleID <- rownames(otu)

index$Family_ID <- rownames(index)
meta <- as.data.frame(as.matrix(infant@sam_data))
meta_index <- full_join(index,meta, by = "Family_ID",multiple = "all")
meta_index <- meta_index %>% filter(!is.na(index_m_h))

meta_otu <- full_join(otu,meta_index, by = "MMHP_SampleID")
meta_otu <- subset(meta_otu, !grepl("bacteria", MMHP_SampleID, ignore.case = TRUE))

rm(index,meta,meta_index,otu,infant)

age <- c("5 months", "6 months", "9 months", "11 months", "14 months")
meta_otu <- meta_otu %>% filter(Age_individual %in% age)
meta_otu$Gordonibacter_pamelaeae <- as.numeric(meta_otu$Gordonibacter_pamelaeae)
meta_otu$Blautia_obeum <- as.numeric(meta_otu$Blautia_obeum)
meta_otu$Prevotella_timonensis <- as.numeric(meta_otu$Prevotella_timonensis)
meta_otu$index_m_h <- as.numeric(meta_otu$index_m_h)

rm(sig_baku,age)

# scatter plots bacteria per TP ====
meta_otu %>%
  ggplot(aes(x=index_m_h, y=Gordonibacter_pamelaeae)) + 
  geom_point() +
  facet_wrap(~Age_individual) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Correlation between Exposure index and Gordonibacter pamelaeae abundance") +
  labs(x = "Exposure index", 
       y = "Abundance") +
  theme_minimal()

meta_otu %>%
  ggplot(aes(x=index_m_h, y=Blautia_obeum)) + 
  geom_point() +
  facet_wrap(~Age_individual) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Correlation between Exposure index and Blautia obeum abundance") +
  labs(x = "Exposure index", 
       y = "Abundance") +
  theme_minimal()

meta_otu %>%
  ggplot(aes(x=index_m_h, y=Prevotella_timonensis)) + 
  geom_point() +
  facet_wrap(~Age_individual) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggtitle("Correlation between Exposure index and Prevotella timonensis abundance") +
  labs(x = "Exposure index", 
       y = "Abundance") +
  theme_minimal()

# scatter plots baku per TP ====
# extracting rho and p values ====
# Filter data and fit lm model
data_9m <- meta_otu %>% filter(Age_individual == "9 months")
model9 <- lm(Gordonibacter_pamelaeae ~ index_m_h, data = data_9m)

data_14m <- meta_otu %>% filter(Age_individual == "14 months")
model14 <- lm(Gordonibacter_pamelaeae ~ index_m_h, data = data_14m)

# Extract stats
r_squared9 <- signif(summary(model9)$r.squared, 3) #R-squared / coefficient  
rho9 <- signif(cor(data_9m$index_m_h, data_9m$Gordonibacter_pamelaeae, use = "complete.obs"), 3) #Pearson's correlation coefficient 
p_val9 <- round(summary(model9)$coefficients[2, 4], 3)  #p-value for the slope , forcing to round up to 3

r_squared14 <- round(summary(model14)$r.squared * 1000, 0) / 1000 # forcing to round up to 3
rho14 <- signif(cor(data_14m$index_m_h, data_14m$Gordonibacter_pamelaeae, use = "complete.obs"), 3) #Pearson's correlation coefficient 
p_val14 <- signif(summary(model14)$coefficients[2, 4], 3) #p-value for the slope 

# scatterplots ====
plot9 <- data_9m %>% 
  ggplot(aes(x = index_m_h, y = Gordonibacter_pamelaeae)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, 
           label = paste("R² =", r_squared9, "\nrho =", rho9, "\np =", p_val9),   
           size = 3, fontface = "bold") +
  ggtitle(expression(italic("Gordonibacter pamelaeae")~"at 9 months")) + 
  labs(x = "", y = "Abundance") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 1), 
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 9))+
  scale_y_continuous(labels = scales::label_scientific(digits = 2))

plot9

plot14 <- data_14m %>% 
  ggplot(aes(x = index_m_h, y = Gordonibacter_pamelaeae)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, 
           label = paste("R² =", r_squared14, "\nrho =", rho14, "\np =", p_val14),
           size = 3, fontface = "bold") +
  ggtitle(expression(italic("Gordonibacter pamelaeae")~"at 14 months")) + 
  labs(x = "", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 1), 
        axis.title.y = element_text(size = 1),
        plot.title = element_text(size = 9))+
  scale_y_continuous(labels = scales::label_scientific(digits = 2))

plot14

sci1 <- function(x) formatC(x, format = "e", digits = 1)
plot9  <- plot9  + scale_y_continuous(labels = sci1)
plot14 <- plot14 + scale_y_continuous(labels = sci1)

comb <- (plot9+plot14)
comb

# saving figures for the paper ====
#ggsave("S11.2_Scatter_plot_baku_TP_index.pdf", comb, device = pdf,width = 300, height = 200, dpi = 300, units = "mm")

ggsave(
  filename   = "S11.2_Scatter_plot_baku_TP_index.tiff",
  plot       = comb,
  width      = 170,      # full page width
  height     = 85,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

