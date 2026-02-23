# LOAD LIBRARIES ====
library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
index_raw <- read.delim("Output/Files/S1.7_index_new_variables.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# preparing for barplot ====
index_raw <- index_raw %>% select(-Family_ID)

df_long_1 <- tidyr::pivot_longer(index_raw, cols = everything())

df_long_1 <- df_long_1 %>% mutate(
  parent = case_when(
    grepl("^m_", name) ~ "mother",
    grepl("^f_", name) ~ "father",
    TRUE              ~ "family"
  )
)

df_long_1$name[df_long_1$name == "home_in_contacts"] <- "Contact at home"
df_long_1$name[df_long_1$name == "home_out_contacts"] <- "Contact outside"
df_long_1$name[df_long_1$name == "f_essen_job"] <- "Essential job"
df_long_1$name[df_long_1$name == "m_essen_job"] <- "Essential job"
df_long_1$name[df_long_1$name == "f_low_hygiene"] <- "Low hygiene"
df_long_1$name[df_long_1$name == "m_low_hygiene"] <- "Low hygiene"
df_long_1$name[df_long_1$name == "f_work_contact"] <- "Contact at work"
df_long_1$name[df_long_1$name == "m_work_contact"] <- "Contact at work"
df_long_1$name[df_long_1$name == "f_protect_gear_work"] <- "PG at work"
df_long_1$name[df_long_1$name == "m_protect_gear_work"] <- "PG at work"
df_long_1$name[df_long_1$name == "f_otherhome_in_contacts"] <- "Contact at OH"
df_long_1$name[df_long_1$name == "m_otherhome_in_contacts"] <- "Contact at OH"
df_long_1$name[df_long_1$name == "f_pubspace_in_contacts"] <- "Contact inside PS"
df_long_1$name[df_long_1$name == "m_pubspace_in_contacts"] <- "Contact inside PS"
df_long_1$name[df_long_1$name == "f_pubspace_out_contacts"] <- "Contact outside PS"
df_long_1$name[df_long_1$name == "m_pubspace_out_contacts"] <- "Contact outside PS"
df_long_1$name[df_long_1$name == "f_otherhome_out_contacts"] <- "Contact outside OH"
df_long_1$name[df_long_1$name == "m_otherhome_out_contacts"] <- "Contact outside OH"

my_cols <- c(
  "family" = "orange",
  "father" = "red3",
  "mother" = "darkgreen")

# family barplots ====
home <- df_long_1 %>% filter(name == "Contact at home") %>% ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

outhome <- df_long_1 %>% filter(name == "Contact outside") %>% ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

dummy <- df_long_1 %>% ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 5, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  labs(fill = "Question")+theme(
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

dummy2 <- df_long_1 %>% ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 5, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  labs(fill = "Question")+theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9),
    strip.text = element_text(size = 9))

# father's plots ====
workf <- df_long_1 %>% filter(name == "Contact at work") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

essenf <- df_long_1 %>% filter(name == "Essential job") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

hyginf <- df_long_1 %>% filter(name == "Low hygiene") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_x_continuous(breaks = scales::breaks_width(2)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

gearf <- df_long_1 %>% filter(name == "PG at work") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

otherhomef <- df_long_1 %>% filter(name == "Contact at OH") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

indoorf <- df_long_1 %>% filter(name == "Contact inside PS") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

outpubf <- df_long_1 %>% filter(name == "Contact outside PS") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

otheroutf <- df_long_1 %>% filter(name == "Contact outside OH") %>% filter(parent == "father") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

# mother's plots ====
workm <- df_long_1 %>% filter(name == "Contact at work") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

essenm <- df_long_1 %>% filter(name == "Essential job") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_x_continuous(breaks = scales::breaks_width(1)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

hyginm <- df_long_1 %>% filter(name == "Low hygiene") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_x_continuous(breaks = scales::breaks_width(2)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

gearm <- df_long_1 %>% filter(name == "PG at work") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks = scales::breaks_width(5))+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

otherhomem <- df_long_1 %>% filter(name == "Contact at OH") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

indoorm <- df_long_1 %>% filter(name == "Contact inside PS") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_x_continuous(breaks = scales::breaks_width(2))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

outpubm <- df_long_1 %>% filter(name == "Contact outside PS") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

otheroutm <- df_long_1 %>% filter(name == "Contact outside OH") %>% filter(parent == "mother") %>% 
  ggplot(aes(x = value, fill = parent)) +
  stat_bin(binwidth = 1, geom = "bar", colour = "white")+theme_bw()+
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_fill_manual(values = my_cols) +
  xlab("") + #ylab("")+
  facet_wrap(.~name, scales = "free")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 9))

# plotting together ====
all <- wrap_plots(
  hyginf, essenf, workf, gearf,
  hyginm, essenm, workm, gearm,
  otheroutf, otherhomef, indoorf, outpubf,
  otheroutm, otherhomem, indoorm, outpubm,
  home, outhome, plot_spacer(),dummy,
  ncol = 4, nrow = 5)

all

all2 <- wrap_plots(
  hyginf, essenf, workf, gearf,
  hyginm, essenm, workm, gearm,
  otheroutf, otherhomef, indoorf, outpubf,
  otheroutm, otherhomem, indoorm, outpubm,
  home, outhome, plot_spacer(),dummy2,
  ncol = 4, nrow = 5)

all2

# saving figures for the paper ====
#save both figures with dummy legend and without, later copy-paste legend instead of the dummy figure
ggsave(
  filename   = "S8.2_index_histogram_new_variable2.tiff",
  plot       = all,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

ggsave(
  filename   = "S8.2_index_histogram_new_variable.tiff",
  plot       = all2,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)

#ggsave("S8.2_index_histogram_new_variable.pdf", all, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")



























