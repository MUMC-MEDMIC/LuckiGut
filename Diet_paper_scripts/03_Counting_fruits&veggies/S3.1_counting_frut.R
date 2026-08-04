# libraries ====
library(microViz)
library(tidyverse)
library(phyloseq)
library(broom)
library(FSA)
library(ggpubr) #for manual p values
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
baku_f_I <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
load("Output/Files/S1.2_timepoints_general_4m_to_14m_famcodes.RData") #list, with six elements (timepoints)

# merging food into one df ====
mon4 <- timepoints_general_4m_to_14m[[1]]
mon5 <- timepoints_general_4m_to_14m[[2]]
mon6 <- timepoints_general_4m_to_14m[[3]]
mon9 <- timepoints_general_4m_to_14m[[4]]
mon11 <- timepoints_general_4m_to_14m[[5]]
mon14 <- timepoints_general_4m_to_14m[[6]]

rownames(mon4) <- mon4$Family_ID
rownames(mon5) <- mon5$Family_ID
rownames(mon6) <- mon6$Family_ID
rownames(mon9) <- mon9$Family_ID
rownames(mon11) <- mon11$Family_ID
rownames(mon14) <- mon14$Family_ID

mon4$Family_ID <- NULL
mon5$Family_ID <- NULL
mon6$Family_ID <- NULL
mon9$Family_ID <- NULL
mon11$Family_ID <- NULL
mon14$Family_ID <- NULL

rm(timepoints_general_4m_to_14m)

# counting how many fruits there are ====
# making all values numeric
mon4 <- mon4 %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
mon5 <- mon5 %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
mon6 <- mon6 %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
mon9 <- mon9 %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
mon11 <- mon11 %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
mon14 <- mon14 %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

# counting fruits
fruits <- c("apple", "banana", "kiwi", "melon", "orange", "peach", "pear", "strawberry")

mon4 <- mon4 %>%   mutate(fruits = rowSums(select(mon4, all_of(fruits)), na.rm = TRUE))
mon5 <- mon5 %>%   mutate(fruits = rowSums(select(mon5, all_of(fruits)), na.rm = TRUE))
mon6 <- mon6 %>%   mutate(fruits = rowSums(select(mon6, all_of(fruits)), na.rm = TRUE))
mon9 <- mon9 %>%   mutate(fruits = rowSums(select(mon9, all_of(fruits)), na.rm = TRUE))
mon11 <- mon11 %>%   mutate(fruits = rowSums(select(mon11, all_of(fruits)), na.rm = TRUE))
mon14 <- mon14 %>%   mutate(fruits = rowSums(select(mon14, all_of(fruits)), na.rm = TRUE))

mon4 <- mon4 %>% select(fruits)
mon5 <- mon5 %>% select(fruits)
mon6 <- mon6 %>% select(fruits)
mon9 <- mon9 %>% select(fruits)
mon11 <- mon11 %>% select(fruits)
mon14 <- mon14 %>% select(fruits)

#taking rownames out
mon4$Family_ID <- rownames(mon4)
mon5$Family_ID <- rownames(mon5)
mon6$Family_ID <- rownames(mon6)
mon9$Family_ID <- rownames(mon9)
mon11$Family_ID <- rownames(mon11)
mon14$Family_ID <- rownames(mon14)

# adding age
mon4$Age_individual <- "4 months"
mon5$Age_individual <- "5 months"
mon6$Age_individual <- "6 months"
mon9$Age_individual <- "9 months"
mon11$Age_individual <- "11 months"
mon14$Age_individual <- "14 months"

food <- rbind(mon4,mon5,mon6,mon9,mon11,mon14)

food <- food %>%  mutate(merged = paste(Age_individual, Family_ID, sep = "@"))
food$Family_ID <- NULL
food$Age_individual <- NULL

rm(mon4,mon5,mon6,mon9,mon11,mon14)

# merging with other metadata ====
meta <- as.data.frame(as.matrix(baku_f_I@sam_data))
meta <- meta %>%  mutate(merged = paste(Age_individual, Family_ID, sep = "@"))

mega <- full_join(food, meta, by = "merged")
mega <- mega %>% filter(!is.na(Sample_aliquote)) 

mega$Age_individual<- factor(mega$Age_individual, 
                             levels = c("4 months",
                                        "5 months",
                                        "6 months",
                                        "9 months",
                                        "11 months",
                                        "14 months"
                             ))

mega$diet_class <- as.factor(mega$diet_class)

rm(baku_f_I,food,meta,fruits)

# Kruskal Wallis test, as n count is small ====
kruskal <- mega %>%
  group_by(Age_individual) %>% 
  summarise(
    kruskal_test = list(kruskal.test(fruits ~ diet_class)),  # Apply Kruskal-Wallis test
    .groups = 'drop'
  ) %>%
  mutate(test_results = map(kruskal_test, tidy)) %>%
  unnest(test_results)

#add significance
kruskal <- kruskal %>% 
  mutate(sig = if_else(
    condition = .$p.value > 0.05, 
    true = NA, 
    false = if_else(
      condition = .$p.value < 0.001,
      true = "***",
      false = if_else(
        condition = .$p.value < 0.05&.$p.value > 0.01,
        true = "*",
        false = "**"))))

kruskal <- data.frame(
  h_statistic = kruskal$statistic,
  parameter = kruskal$parameter,
  p_value = kruskal$p.value,
  Age_individual = kruskal$Age_individual,
  method = kruskal$method,
  sig = kruskal$sig
)

#write.table (kruskal, "S3.1_kruskal_count_fruits.txt", row.names=TRUE,sep = "\t")

# Dunn test ====
dunn <- mega  %>% 
  group_by(Age_individual) %>%
  summarise(
    dunn_test = list(
      dunnTest(fruits ~ diet_class, method = "bh")),
    .groups = "drop"
  ) %>%
  mutate(pairwise_results = map(dunn_test, ~.$res)) %>%
  unnest(pairwise_results)

dunn <- dunn %>% 
  mutate(sig = if_else(
    condition = .$P.unadj > 0.05, 
    true = NA, 
    false = if_else(
      condition = .$P.unadj < 0.001,
      true = "***",
      false = if_else(
        condition = .$P.unadj < 0.05&.$P.unadj > 0.01,
        true = "*",
        false = "**"))))

dunn <- dunn %>% 
  mutate(sig_adj = if_else(
    condition = .$P.adj > 0.05, 
    true = NA, 
    false = if_else(
      condition = .$P.adj < 0.001,
      true = "***",
      false = if_else(
        condition = .$P.adj < 0.05&.$P.adj > 0.01,
        true = "*",
        false = "**"))))


dunn <- data.frame(
  Comparison = dunn$Comparison,
  Z = dunn$Z,
  p_value = dunn$P.unadj,
  p_value_adj = dunn$P.adj,
  Age_individual = dunn$Age_individual,
  sig_adj = dunn$sig_adj,
  sig = dunn$sig
)

dunn <- dunn %>% separate(Comparison, into = c("group1", "group2"), sep = "-")
dunn$group1 <- gsub(" ", "", dunn$group1)
dunn$group2 <- gsub(" ", "", dunn$group2)
dunn$group1 <- as.factor(dunn$group1)
dunn$group2 <- as.factor(dunn$group2)

#write.table (dunn, "S3.1_dunn_count_fruits.txt", row.names=TRUE,sep = "\t")

# preparing for boxplot ====
mon4 <- mega %>% filter(Age_individual == "4 months")
mon5 <- mega %>% filter(Age_individual == "5 months")
mon6 <- mega %>% filter(Age_individual == "6 months")
mon9 <- mega %>% filter(Age_individual == "9 months")
mon11 <- mega %>% filter(Age_individual == "11 months")
mon14 <- mega %>% filter(Age_individual == "14 months")

dunn4 <- dunn %>% filter(Age_individual == "4 months")
dunn5 <- dunn %>% filter(Age_individual == "5 months")
dunn6 <- dunn %>% filter(Age_individual == "6 months")
dunn9 <- dunn %>% filter(Age_individual == "9 months")
dunn11 <- dunn %>% filter(Age_individual == "11 months")
dunn14 <- dunn %>% filter(Age_individual == "14 months")

# boxplot ====
ggplot(mega, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class)) + 
  facet_grid( ~ Age_individual) +
  ggtitle("") +
  xlab("Diet class") +
  ylab("Fruit count")  + theme_minimal()+
  theme(legend.position = "none")

# boxplot with significance ====
plot1 <- ggplot(mon4, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 1.2, alpha = 0.7,
    color = "black", shape = 16
  )+
    ggtitle("") +
  facet_grid( ~ Age_individual) +
  xlab("") +
  ylab("Count of fruits") + theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 10)),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.4),
        panel.grid.minor = element_blank())
plot1 <- plot1 + stat_pvalue_manual(dunn4,y.position = 5.5, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot2 <- ggplot(mon5, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 1.2, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle("") +
  facet_grid( ~ Age_individual) +
  xlab("") +
  ylab("") + theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 10)),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.4),
        panel.grid.minor = element_blank())
plot2 <- plot2 + stat_pvalue_manual(dunn5,y.position = 5.5, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot3 <- ggplot(mon6, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 1.2, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle("") +
  facet_grid( ~ Age_individual) +
  xlab("") +
  ylab("") + theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 10)),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.4),
        panel.grid.minor = element_blank())
plot3 <- plot3 + stat_pvalue_manual(dunn6,y.position = 6.5, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot4 <- ggplot(mon9, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 1.2, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle("") +
  facet_grid( ~ Age_individual) +
  xlab("") +
  ylab("Count of fruits") + theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 10)),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.4),
        panel.grid.minor = element_blank())
plot4 <- plot4 + stat_pvalue_manual(dunn9,y.position = 8.5, step.increase = 0.05,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot5 <- ggplot(mon11, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 1.2, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle("") +
  facet_grid( ~ Age_individual) +
  xlab("Diet class") +
  ylab("") + theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 10)),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.4),
        panel.grid.minor = element_blank())
plot5 <- plot5 + stat_pvalue_manual(dunn11,y.position = 8.5, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot6 <- ggplot(mon14, aes(x = diet_class, y = fruits)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 1.2, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle("") +
  facet_grid( ~ Age_individual) +
  xlab("") +
  ylab("") + theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 7, margin = margin(b = 10)),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.4),
        panel.grid.minor = element_blank())
plot6 <- plot6 + stat_pvalue_manual(dunn14,y.position = 8.5, step.increase = 0.05,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

combined_plot <- (plot1 + plot2 + plot3) / (plot4 + plot5 + plot6)
combined_plot

# saving figure for the paper ====
ggsave(
  filename   = "S3.1_fruit_count.tiff",
  plot       = combined_plot,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)
