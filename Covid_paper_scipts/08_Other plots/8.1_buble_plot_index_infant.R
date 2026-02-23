# load libraries ====
library(tidyverse)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index = read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# selecting families with index ====
# calculating number of samples and sets
meta <- as.data.frame(as.matrix(infant@sam_data))
rm(infant)

index$Family_ID <- rownames(index)
families <- sort(unique(index$Family_ID))
index <- meta %>% dplyr::filter(Family_ID %in% families)
plyr::count(index$Age_individual)

rm(families)

# selecting variables ==== 
index <- index %>% select(all_of(c("Age_individual","covid_270220","Family_ID")))

# incorporating missing TP ====
samples <- index %>% group_by(Family_ID) %>% count(Age_individual)
age <- unique(index$Age_individual)

df_completed <- index %>%
  distinct(Family_ID) %>%
  crossing(Age_individual = age) %>%
  left_join(index, by = c("Family_ID", "Age_individual")) %>% 
  mutate(across(everything(), ~replace_na(., "missing/pandemic")))

df_completed <- df_completed %>%
  mutate(Infant_ID = dense_rank(Family_ID))

df_completed$bubble_size <- 1

df_completed$Age_individual<- factor(df_completed$Age_individual, levels = c(
                                         "1-2 weeks",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))
# cleaning meta for plot ====
df_completed <- df_completed %>%
  mutate(covid_270220 = recode(covid_270220,  
                       "no" = "prepandemic",
                       "yes" = "pandemic"))
df_completed$Infant_ID <- as.character(df_completed$Infant_ID)
df_completed$Infant_ID <- factor(df_completed$Infant_ID, levels = as.character(1:36))

# filtering out samples that are in the freezer, but were not sequenced ====
df_completed_filt <- df_completed %>%
  filter(!(Family_ID %in% c("1228", "1226", "1221") &
             Age_individual %in% c("14 months", "11 months", "9 months", "6 months")))

# checking the grey families in the beginning ====
grey <- df_completed_filt %>% filter(Infant_ID %in% c(1,2,3,4, 5, 6, 7,8, 14,18,20))
grey <- sort(unique(grey$Family_ID))

meta <- meta %>% filter(Family_ID %in% grey)
dob <- unique(meta$DOB_infant)

ref <- as.Date("2020-02-27")

babies_born_before_covid <- dob[dob < ref]

rm(dob,ref)

grey <- meta %>% filter(Family_ID %in% grey) %>% filter(DOB_infant %in% babies_born_before_covid) 
grey <- sort(unique(grey$Family_ID))
grey <- df_completed_filt %>% filter(Family_ID %in% grey)
babies_born_before_covid <- as.numeric(sort(unique(grey$Infant_ID)))

rm(meta,grey)

# marking those who were born before covid 
df_completed_filt <- df_completed_filt %>%
  mutate(
    covid_270220 = if_else(Infant_ID == "4" & Age_individual == "1-2 weeks",
                           "missing/prepandemic",   
                           covid_270220))
df_completed_filt <- df_completed_filt %>%
  mutate(
    covid_270220 = if_else(Infant_ID == "6" & Age_individual == "1-2 weeks",
                           "missing/prepandemic",   
                           covid_270220))

plyr::count(df_completed_filt$covid_270220)


# bubble plot ====
custom_colors <- c(
  "missing/pandemic" = "grey80",
  "prepandemic" = "#0039A6",
  "pandemic" = "#D52B1E",
  "missing/prepandemic" = "#669999"
)

plot1 <- df_completed_filt %>% ggplot(aes(x=Age_individual, y=Infant_ID, size = bubble_size,color = covid_270220)) +
  geom_line(aes(group = Infant_ID), color = "black", linewidth = 0.5) +
  geom_point(aes(fill = covid_270220), 
             shape = 21, 
             color = "black", 
             stroke = 0.5, 
             size = 4)+ylab("") + xlab ("")+
  theme_minimal()+
  theme(
    # legend.position = "none",
    axis.text.x = element_text(size = 9),        # hides labels
    axis.title.x = element_blank(),       # hides axis title
    axis.ticks.x = element_blank(),       # hides ticks
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, margin = margin(b = 5)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_fill_manual(values = custom_colors,name = "Sample status" )+
  scale_y_discrete(limits = rev) +
  guides(size = "none")

plot1

# saving figures for the paper ====
#ggsave("S8.1_Bubble_plot_index.pdf", plot1, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")

ggsave(
  filename   = "S8.1_Bubble_plot_index.tiff",
  plot       = plot1,
  width      = 225,      # full page width
  height     = 160,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

# saving data ====
#write.table (df_completed_filt, "S8.1_recoded_infnats_for_bubble_plot.txt", row.names=FALSE,sep = "\t")

