# libraries ====
library(tidyverse)
library(microViz)
library(phyloseq)
library(viridis)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

index$Family_ID <- as.numeric(rownames(index))

#merging with index 
baku_f_I_index <- ps_join(infant,index, by = "Family_ID")
baku_f_I_rel <- baku_f_I_index %>% tax_transform(trans = "compositional")

# filtering Gordoni bacter ====
gordoni <- subset_taxa(baku_f_I_rel, (Species == "Gordonibacter_pamelaeae"))
gordoni_filt <- sample_data(gordoni)[, c("MMHP_SampleID", "Age_individual", "index_m_h","covid_270220")]

gordoni_filt <- as.data.frame(as.matrix(gordoni_filt))
gordoni_filt <- gordoni_filt %>% mutate(comp = ifelse(covid_270220 == "no", "prepandemic", NA))
mean <- mean(as.numeric(gordoni_filt$index_m_h), na.rm = TRUE)
gordoni_filt <- gordoni_filt %>%
  mutate(comp = ifelse(is.na(comp), 
                       ifelse(index_m_h < mean, "below index mean", "above index mean"), 
                       comp))
meta <- gordoni_filt %>% mutate(comp = ifelse(is.na(comp), "pandemic", comp))

plyr::count(meta$comp)

otu <- as.data.frame(as.matrix(t(gordoni@otu_table)))
otu$MMHP_SampleID <- rownames(otu)

all <- full_join(otu,meta, by = "MMHP_SampleID")

rm(baku_f_I,baku_f_I_rel,baku_f_I_rel,gordoni,gordoni_filt,index,mean)

all$Age_individual<- factor(all$Age_individual, levels = c("1-2 weeks",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))

#making new groups
subggroup <- all %>% filter(comp != "prepandemic")
subggroup$comp <- "pandemic"

all <- all %>% filter(comp != "pandemic")

all_with_doubles <- rbind(all,subggroup)

all_with_doubles$comp <- factor(all_with_doubles$comp, levels = c("prepandemic",
                                                           "pandemic",
                                                           "above index mean",
                                                           "below index mean"))

#all_with_doubles$Gordonibacter_pamelaeae_log10 <- log10(all_with_doubles$Gordonibacter_pamelaeae+1)
all_with_doubles$Gordonibacter_pamelaeae_sqrt <- sqrt(all_with_doubles$Gordonibacter_pamelaeae+0.001)

# boxplot all =====
all_with_doubles %>%
  ggplot( aes(x=comp, y=Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Scaled Gordonibacter pamelaeae reads per TP") +
  xlab("")+facet_wrap(~Age_individual)

# boxplot per TP ====
plot1 <- all_with_doubles %>% filter(Age_individual == "1-2 weeks") %>% 
ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +  ggtitle("Gordonibacter pamelaeae abundace per group") +
  xlab("")+facet_wrap(~Age_individual) +ylab("")
  
plot2 <- all_with_doubles %>% filter(Age_individual == "4 weeks") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+facet_wrap(~Age_individual)+ylab("")

plot3 <- all_with_doubles %>% filter(Age_individual == "8 weeks") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+facet_wrap(~Age_individual)+ylab("")

plot4 <- all_with_doubles %>% filter(Age_individual == "4 months") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+facet_wrap(~Age_individual)+
  ylab("Square root of relative abundance")

plot5 <- all_with_doubles %>% filter(Age_individual == "5 months") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+facet_wrap(~Age_individual)+ylab("")

plot6 <- all_with_doubles %>% filter(Age_individual == "6 months") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+facet_wrap(~Age_individual)+ylab("")+ylim(0.03,0.05)

plot7 <- all_with_doubles %>% filter(Age_individual == "9 months") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 1), 
    axis.title.y = element_text(size = 10)) +
  xlab("")+facet_wrap(~Age_individual)+ylab("Square root of relative abundance")

plot8 <- all_with_doubles %>% filter(Age_individual == "11 months") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+facet_wrap(~Age_individual)+ylab("")+ylim(0.03,0.05)

plot9 <- all_with_doubles %>% filter(Age_individual == "14 months") %>% 
  ggplot(aes(x = comp, y = Gordonibacter_pamelaeae_sqrt, fill=comp)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 1), 
    axis.title.y = element_text(size = 13)) +
  xlab("")+facet_wrap(~Age_individual)+ylab("")+ylim(0.03,0.05)

combined_plot <- (plot1 + plot2 + plot3) / (plot4 + plot5 + plot6)/ (plot7 + plot8 + plot9)
combined_plot

# plot only for sig TP ====
combined_plot2 <- (plot7 + plot9)
combined_plot2

# saving figures for the paper ====
#ggsave("S8.4_gordoni_boxplot.pdf", combined_plot2, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")

ggsave(
  filename   = "S8.4_gordoni_boxplot.tiff",
  plot       = combined_plot2,
  width      = 170,      # full page width
  height     = 85,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

