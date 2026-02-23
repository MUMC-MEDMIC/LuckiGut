# load libraries ====
library(tidyverse)
library(phyloseq)
library(scales)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index = read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
timeline <- read.delim("Input/18_covid_restriction_timeline.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
bubble_codes <- read.delim("Output/Files/S8.1_recoded_infnats_for_bubble_plot.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# selecting families with index ====
# calculating number of samples and sets
meta <- as.data.frame(as.matrix(infant@sam_data))
rm(infant)

index$Family_ID <- rownames(index)
families <- sort(unique(index$Family_ID))
index <- meta %>% dplyr::filter(Family_ID %in% families)
plyr::count(index$Age_individual)

rm(families,meta)

# selecting dob and last sample date ====
last_sample_date <- index %>%
  arrange(Family_ID, Sample_date_filled) %>%   # ensure dates are in time order
  group_by(Family_ID) %>%                     # one group per subject
  slice_tail(n = 1)                        # keep last row in each group

last_sample_date <- last_sample_date %>% select(Family_ID, Sample_date_filled)
last_sample_date <- last_sample_date %>% rename(last_date = Sample_date_filled)

# creating keys ====
bubble_codes <- bubble_codes %>% select(Family_ID,Infant_ID) %>% unique()
bubble_codes$Family_ID <- as.character(bubble_codes$Family_ID)

# selecting DOB for timeline ====
dob <- index %>% select(Family_ID, DOB_infant) %>% unique()

# cleaning timeline ====
timeline$Start.Date <- as.Date(timeline$Start.Date, format = "%d/%m/%Y")
timeline$covid_timeline_dates <- timeline$Start.Date

timeline <- timeline %>% rename(common_date = Start.Date)

# merging together ====
all <- full_join(bubble_codes,last_sample_date, by = "Family_ID")
all <- full_join(all,dob, by = "Family_ID")

all <- all %>% rename(first_date = DOB_infant)

all$first_date <- as.Date(all$first_date)
all$last_date <- as.Date(all$last_date)
all$Infant_ID <- factor(all$Infant_ID)

all$common_date <- all$first_date

all_with_timeline <- full_join(all,timeline, by = "common_date")

rm(dob,last_sample_date, bubble_codes,index)

# selecting dates for filtering ====
oldest_sample_date <- min(all$last_date, na.rm = TRUE)
youngest_sample_date <- max(all$last_date, na.rm = TRUE)

# preparing for the plot ====
relevant_timeline <- timeline %>% filter(covid_timeline_dates < youngest_sample_date) 

# gantt plot ====
vlines <- tibble::tibble(date = as.Date(c("2020-02-20")), #"2020-03-12"
                        label = c("Pandemic reaches the NL")) # "First measures"

windows <- tibble::tibble(
  start = as.Date(c("2020-03-12", "2020-12-14","2021-12-19")),
  end   = as.Date(c("2020-05-11","2021-02-08","2022-01-15")),
  window = c("Lockdown1", "Lockdown2", "Lockdown3"))

# first measures in place 2020-03-12
# last measures in place 2023-03-10
# last measures relevant for the cohort 2022-04-11

plot1 <- all %>%
  ggplot(aes(y = Infant_ID)) +
  geom_rect(data = windows,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = window),
    inherit.aes = FALSE,
    alpha = 0.8) +
  scale_fill_manual(values = c("Lockdown1" = "darkgreen", 
                               "Lockdown2" = "red3", 
                               "Lockdown3" = "orange")) +
  geom_segment(aes(x = first_date, xend = last_date, yend = Infant_ID),
    linewidth = 3, lineend = "butt", colour = "grey80") +
  # Add vertical line
  geom_vline(data = vlines,
    aes(xintercept = date),
    linetype = "dashed",
    color = "black",
    linewidth = 1) +
  geom_text(
    data = vlines,
    aes(x = date+10, y = 2.5, label = label),
    vjust = 1.5, hjust = 0,
    inherit.aes = FALSE,
    size = 6) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90,size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10)) +
  scale_y_discrete(limits = rev)

# saving figures for the paper ====
#ggsave("S8.5_timeline_pandemic_plot.pdf", plot1, device = "pdf",width = 500, height = 300, dpi = 300, units = "mm")

ggsave(
  filename   = "S8.5_timeline_pandemic_plot.tiff",
  plot       = plot1,
  width      = 225,      # full page width
  height     = 160,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)
