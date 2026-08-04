# libraries ====
library(tidyverse)
library(microViz)
library(viridis)
library(broom)
library(FSA)
library(ggpubr) #for manual p values
library(patchwork)
library(svglite)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ps <- ps_diet <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
baku <- read.delim("Output/Files/S5.3_cleaned_reads_bakut.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
kegg_codes <- read.delim("Output/Files/S5.5_dunn_kegg_between_diets_per_TP.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
scfa_meta <- read.delim("Output/Files/S5.1_cleaned_scfa_meta.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# turning to long format ====
baku_long <- baku %>%
  pivot_longer(
    cols = starts_with("MMHP"), # Columns to pivot
    names_to = "MMHP_SampleID",      # Name for the new variable column
    values_to = "read_count"         # Name for the new value column
  )

#removing bacteria with 0 reads 
#baku_long <- baku_long %>% filter(read_count > 0)

# adding scfa meta to baku list ====
bakut_long_scfa <- right_join(scfa_meta, baku_long, by = "KEGG_code",multiple = "all")

# taking kegg codes out, that were sig ====
kegg_codes <- kegg_codes %>% filter(p_value_adj < 0.05)
kegg <- kegg_codes
kegg_codes <- unique(sort(kegg_codes$KEGG_code))

bakut_long_scfa_sig  <- bakut_long_scfa %>% filter(KEGG_code %in% kegg_codes)

rm(baku, baku_long,bakut_long_scfa,scfa_meta,kegg_codes)

kegg <- kegg %>% filter(Age_individual == "4 months") %>% filter(p_value_adj < 0.05) %>% filter(scfa == "butyrate")
kegg_4mon_sig <- sort(unique(kegg$KEGG_code))

# adding meta from phyloseq ====
meta <- as.data.frame(as.matrix(ps@sam_data))
meta <- meta %>% select(MMHP_SampleID,starts_with("diet"),Age_individual,ageintrosolids,sex)

bakut_long_scfa_sig_meta <- inner_join(meta, bakut_long_scfa_sig, by = "MMHP_SampleID",multiple = "all")

rm(ps, meta,bakut_long_scfa_sig)

check <- length(unique(bakut_long_scfa_sig_meta$MMHP_SampleID))

# counting top 10 abundant species and preparing data ====
bakut_long_scfa_sig_meta_4mon<- bakut_long_scfa_sig_meta %>% filter(Age_individual == "4 months")%>% filter(KEGG_code %in% kegg_4mon_sig)
top_species <- bakut_long_scfa_sig_meta_4mon %>%
  group_by(species) %>%
  summarise(total_reads = sum(read_count), .groups = "drop") %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = 10) %>%
  pull(species)  %>% sort()

new_df <- bakut_long_scfa_sig_meta_4mon %>%
  mutate(species_grouped = ifelse(species %in% top_species, species, "Other")) %>%
  group_by(KEGG_code,diet_class,Age_individual, species_grouped,scfa) %>% 
  summarise(read_count = sum(read_count), .groups = "drop") %>%
  arrange(KEGG_code,diet_class,Age_individual, scfa, desc(read_count))

new_df$diet_class <- gsub(" ", "", new_df$diet_class)
new_df$diet_class <- as.character(new_df$diet_class)
new_df$diet_class <- factor(new_df$diet_class, 
                          levels = c("1",
                                     "2",
                                     "3"))

# changing colours ====
plyr::count(new_df$species_grouped)
species_colors <- c(
  "Acidaminococcus intestini" =  "#8B9D77", 
  "Bacteroides caccae" = "#DBD3D1",
  "Bacteroides dorei" = "#84869B",
  "Bacteroides faecis" = "#A56266",
  "Bacteroides thetaiotaomicron" = "#34686C",
  "Bacteroides vulgatus" = "#B3B2A5",
  "Enterococcus faecalis" =  "#C5D3E3",
  "Escherichia coli" = "#003856",
  "Flavonifractor plautii" =  "#F5B8B6",
  "Other" = "#849E8A", 
  "unclassified" = "#D3B8C5"
)

new_df$species_color <- species_colors[new_df$species_grouped]

plyr::count(new_df$Age_individual)

df3 <- new_df %>% filter(scfa == "butyrate")

# butyrate 4 months only =====
#taking sig at 4 months only kegg codes
df3_4mon <- df3 %>% filter(Age_individual == "4 months") %>% filter(KEGG_code %in% kegg_4mon_sig)
df3_4mon$species_grouped <- gsub("_", " ", df3_4mon$species_grouped)

diet_but_1 <- df3_4mon %>% filter(KEGG_code == "K00074")
diet_but_3 <- df3_4mon %>% filter(KEGG_code == "K00929")
diet_but_5 <- df3_4mon %>% filter(KEGG_code == "K14534")
diet_but_6 <- df3_4mon %>% filter(KEGG_code == "K17865")
diet_but_7 <- df3_4mon %>% filter(KEGG_code == "K18014")

dunn_but_kegg_1 <- kegg %>% filter(KEGG_code == "K00074")
dunn_but_kegg_3 <- kegg %>% filter(KEGG_code == "K00929")
dunn_but_kegg_5 <- kegg %>% filter(KEGG_code == "K14534")
dunn_but_kegg_6 <- kegg %>% filter(KEGG_code == "K17865")
dunn_but_kegg_7 <- kegg %>% filter(KEGG_code == "K18014")

# plotting boxplots =====
plot1_with <- ggplot(diet_but_1, aes(x = diet_class, y = read_count)) +
  geom_bar(
    aes(fill = species_grouped),
    position = "fill",
    stat = "identity"
  ) +
  facet_grid(~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL) +
  theme_minimal() +
  scale_fill_manual(values = species_colors) +
  labs(fill = "Top 4 species") +
  coord_cartesian(ylim = c(0, 1.05)) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 1))) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 5, face = "italic"),
    legend.title = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(b = 1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(3, "mm"))
plot1_with <- plot1_with + stat_pvalue_manual(dunn_but_kegg_1,y.position = 1.02, step.increase = 0.00005,tip.length = 0,
                     label = "sig_adj",hide.ns = TRUE)

plot3_with <- ggplot(diet_but_3, aes(x = diet_class, y = read_count)) +
  geom_bar(
    aes(fill = species_grouped),
    position = "fill",
    stat = "identity"
  ) +
  facet_grid(~ KEGG_code) +  # Separate plots for SCFA and TP
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL) +
  theme_minimal() +
  scale_fill_manual(values = species_colors) +
  labs(fill = "Top 10 species") +  # in reality it is 8, but merged with plot 1 legend
  coord_cartesian(ylim = c(0, 1.10)) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 1))) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 5, face = "italic"),
    legend.title = element_text(size = 6),
    strip.text = element_text(size = 7, margin = margin(b = 1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(3, "mm"))
plot3_with <- plot3_with + stat_pvalue_manual(dunn_but_kegg_3,y.position = 1.02, step.increase = 0.0000025,tip.length = 0,
                           label = "sig_adj",hide.ns = TRUE)

plot5 <- ggplot(diet_but_5, aes(x = diet_class, y = read_count)) +
  geom_bar(aes(fill = species_grouped),position = "fill", stat = "identity") +
  facet_grid( ~ KEGG_code) +  # Separate plots for SCFA and TP
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL) + theme_minimal()+
  scale_fill_manual(values = species_colors) +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot5<- plot5 + stat_pvalue_manual(dunn_but_kegg_5,y.position = 1.02, step.increase = 0.000001,tip.length = 0,
                           label = "sig_adj",hide.ns = TRUE)

plot6 <- ggplot(diet_but_6, aes(x = diet_class, y = read_count)) +
  geom_bar(aes(fill = species_grouped),position = "fill", stat = "identity") +
  facet_grid( ~ KEGG_code) +  # Separate plots for SCFA and TP
  ggtitle("Butyrate producing species per KEGG code at 4 months") +
  xlab("Diet class") +
  ylab("Relative abundance") + theme_minimal()+
  scale_fill_manual(values = species_colors) +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(b = 1)), #grid size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot6<- plot6 + stat_pvalue_manual(dunn_but_kegg_6,y.position = 1.02, step.increase = 0.000001,tip.length = 0,
                           label = "sig_adj",hide.ns = TRUE)

plot7 <- ggplot(diet_but_7, aes(x = diet_class, y = read_count)) +
  geom_bar(aes(fill = species_grouped),position = "fill", stat = "identity") +
  facet_grid( ~ KEGG_code) +  # Separate plots for SCFA and TP
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab("Relative abundance") + theme_minimal()+
  scale_fill_manual(values = species_colors)+
  coord_cartesian(ylim = c(0, 1.05)) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot7<- plot7 + stat_pvalue_manual(dunn_but_kegg_7,y.position = 1.02, step.increase = 0.000001,tip.length = 0,
                           label = "sig_adj",hide.ns = TRUE)


combined_plot_with <- (plot6 + plot5 + plot1_with) / (plot7 + plot3_with)+
  plot_layout() & 
  theme(plot.margin = margin(1, 0, 2, 0)) #top, right, bottom, left
combined_plot_with

plot1_without <- ggplot(diet_but_1, aes(x = diet_class, y = read_count)) +
  geom_bar(aes(fill = species_grouped),position = "fill", stat = "identity") +
  facet_grid( ~ KEGG_code) +  # Separate plots for SCFA and TP
  ggtitle(NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))+
  xlab("Diet class") +
  ylab(NULL) + theme_minimal()+
  scale_fill_manual(values = species_colors) + 
  labs(fill = "Top 4 species")+
  coord_cartesian(ylim = c(0, 1.05)) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.position = "none",
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(3, "mm"))+
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 1)))
plot1_without <- plot1_without + stat_pvalue_manual(dunn_but_kegg_1,y.position = 1.02, step.increase = 0.00005,tip.length = 0,
                                                    label = "sig_adj",hide.ns = TRUE)

plot3_without <- ggplot(diet_but_3, aes(x = diet_class, y = read_count)) +
  geom_bar(aes(fill = species_grouped),position = "fill", stat = "identity") +
  facet_grid( ~ KEGG_code) +  # Separate plots for SCFA and TP
  ggtitle(NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))+
  xlab("Diet class") +
  ylab(NULL) + theme_minimal()+
  scale_fill_manual(values = species_colors) + 
  labs(fill = "Top 8 species")+
  coord_cartesian(ylim = c(0, 1.10)) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.position = "none",
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(3, "mm"))+
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 1)))
plot3_without <- plot3_without + stat_pvalue_manual(dunn_but_kegg_3,y.position = 1.02, step.increase = 0.0000025,tip.length = 0,
                                                    label = "sig_adj",hide.ns = TRUE)

combined_plot_without <- (plot6 + plot5 + plot1_without) / (plot7 + plot3_without+plot_spacer())+
  plot_layout() & 
  theme(plot.margin = margin(1, 0, 2, 0)) #top, right, bottom, left
combined_plot_without

# saving figure for the paper ====
ggsave(
  filename   = "S5.6_scfa_production_per_baku_with.svg",
  plot       = combined_plot_with,
  width    = 80,
  height   = 90,
  units    = "mm",
  device   = svglite
)

ggsave(
  filename   = "S5.6_scfa_production_per_baku_without.svg",
  plot       = combined_plot_without,
  width    = 80,
  height   = 90,
  units    = "mm",
  device   = svglite
)

