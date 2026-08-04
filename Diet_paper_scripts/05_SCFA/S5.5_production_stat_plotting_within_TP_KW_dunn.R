# libraries ====
library(tidyverse)
library(microViz)
library(viridis)
library(broom)
library(ggpubr) #for manual p values
library(FSA)
library(patchwork)
library(svglite)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
scfa <- read.delim("Output/Files/S5.4_long_scfa_diet_paper.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
check <- unique(scfa$KEGG_code) #should be 49 with doubles
scfa$diet_class <- gsub(" ", "", scfa$diet_class)
keys <- scfa %>% select(KEGG_code,scfa)

# count kegg code sper scfa ====
count <- scfa %>% select(scfa, KEGG_code)
count <- count %>% distinct(scfa, KEGG_code, .keep_all = TRUE)
plyr::count(count$scfa)
plyr::count(count$KEGG_code)
tp <- count %>% group_by(KEGG_code) %>% tally() %>% filter(n >= 2) 

# prepare data ====
scfa$diet_class <- as.character(scfa$diet_class)
scfa$Age_individual<- factor(scfa$Age_individual, 
                             levels = c("4 months",
                                        "5 months",
                                        "6 months",
                                        "9 months",
                                        "11 months",
                                        "14 months"))

scfa$diet_class <- factor(scfa$diet_class, 
                          levels = c("1",
                                     "2",
                                     "3"))

diet <- scfa %>% na.omit(diet_class)
length(unique(diet$MMHP_SampleID)) #should be 389
rm(scfa)

diet$log_read_counts <- log10(diet$read_count + 1)  
diet$sq_read_counts <- sqrt(diet$read_count+0.001)
double_kegg_codes <- c("K00140","K01034","K01035","K01899","K01900","K01902","K01903","K01913") # all double codes

# normally distributed data? ====
diet1 <- diet %>% filter(diet_class == "1") 
ggplot(diet1, aes(x=read_count)) + 
  geom_histogram(color="black", fill="white")+ facet_wrap(~Age_individual)+
  ggtitle("SCFA related production gene distribution per TP for diet 1")
ggplot(diet1, aes(sample = read_count)) + 
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~Age_individual) +
  ggtitle("SCFA related production gene Q-Q Plot per TP for diet 1")

diet2 <- diet %>% filter(diet_class == "2") 
hist(diet2$read_count)

diet3 <- diet %>% filter(diet_class == "3") 
hist(diet3$read_count)

rm(diet1,diet2,diet3)

# 3.kruskal for KEGG codes at each TP between diet type ====
kruskal_3_keggs <- diet %>% filter(!(scfa == "acetate" & KEGG_code %in% double_kegg_codes)) %>%
  group_by(Age_individual, KEGG_code) %>% 
  summarise(
    kruskal_test = list(kruskal.test(read_count ~ diet_class)),  # Apply Kruskal-Wallis test
    .groups = 'drop'
  ) %>%
  mutate(test_results = map(kruskal_test, tidy)) %>%
  unnest(test_results)

#add significance
kruskal_3_keggs <- kruskal_3_keggs %>% 
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

#if this step does not work, restart R
kruskal_3_keggs <- data.frame(
  h_statistic = kruskal_3_keggs$statistic,
  parameter = kruskal_3_keggs$parameter,
  p_value = kruskal_3_keggs$p.value,
  Age_individual = kruskal_3_keggs$Age_individual,
  KEGG_code = kruskal_3_keggs$KEGG_code,
  method = kruskal_3_keggs$method,
  sig = kruskal_3_keggs$sig
)
check <- unique(kruskal_3_keggs$KEGG_code) #should be 49

# 3.Dunn ====
# Initialize a list to store failed KEGG codes
failed_kegg <- list()

# Perform Dunn test and capture failures
dunn_3_kegg <- diet  %>% filter(!(scfa == "acetate" & KEGG_code %in% double_kegg_codes)) %>%
  group_by(Age_individual, KEGG_code) %>%
  summarise(
    dunn_test = list(tryCatch(
      dunnTest(read_count ~ diet_class, method = "bh"),
      error = function(e) {
        # Capture failed KEGG codes
        failed_kegg <<- append(failed_kegg, unique(KEGG_code))
        NULL  # Return NULL for failed tests
      }
    )),
    .groups = "drop"
  ) %>%
  mutate(pairwise_results = map(dunn_test, ~.$res)) %>%
  unnest(pairwise_results)

#add significance
dunn_3_kegg <- dunn_3_kegg %>% 
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

dunn_3_kegg <- dunn_3_kegg %>% 
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

dunn_3_kegg <- data.frame(
  Comparison = dunn_3_kegg$Comparison,
  Z = dunn_3_kegg$Z,
  p_value = dunn_3_kegg$P.unadj,
  p_value_adj = dunn_3_kegg$P.adj,
  Age_individual = dunn_3_kegg$Age_individual,
  KEGG_code = dunn_3_kegg$KEGG_code,
  sig_adj = dunn_3_kegg$sig_adj,
  sig = dunn_3_kegg$sig
)

dunn_3_kegg <- dunn_3_kegg %>% separate(Comparison, into = c("group1", "group2"), sep = "-")
dunn_3_kegg$group1 <- gsub(" ", "", dunn_3_kegg$group1)
dunn_3_kegg$group2 <- gsub(" ", "", dunn_3_kegg$group2)
check <- unique(dunn_3_kegg$KEGG_code) #should be 48

failed_kegg <- unique(unlist(failed_kegg))
rm(failed_kegg)

# adding scfa to kegg codes ====
#separating codes into double and non double
keys_unik <- keys %>% distinct(KEGG_code, scfa, .keep_all = TRUE)
keys_doubl <- keys_unik %>% filter(KEGG_code %in% double_kegg_codes)
keys_no_doubl <- keys_unik %>% filter(!KEGG_code %in% double_kegg_codes)
rm(keys_unik,keys)

#filtering kruskal into double and no double codes
krus_doubl <- kruskal_3_keggs %>% filter(KEGG_code %in% double_kegg_codes)
krus_no_doubl <- kruskal_3_keggs %>% filter(!KEGG_code %in% double_kegg_codes)

#filtering dunn into double and no double codes
dunn_doubl <- dunn_3_kegg %>% filter(KEGG_code %in% double_kegg_codes)
dunn_no_doubl <- dunn_3_kegg %>% filter(!KEGG_code %in% double_kegg_codes)

#adding scfa to non double codes
kru_final <- full_join(keys_no_doubl,krus_no_doubl, by = "KEGG_code",multiple = "all")
dun_final <- full_join(keys_no_doubl,dunn_no_doubl, by = "KEGG_code",multiple = "all")
dun_final <- dun_final[!is.na(dun_final$Z), ]
rm(keys_no_doubl,dunn_no_doubl,krus_no_doubl,kruskal_3_keggs,dunn_3_kegg)

#adding scfa to double codes
kru <- full_join(keys_doubl,krus_doubl, by = "KEGG_code",multiple = "all")
dun <- full_join(keys_doubl,dunn_doubl, by = "KEGG_code",multiple = "all")
rm(dunn_doubl,keys_doubl,krus_doubl)

#merging tables
kru_final <- rbind(kru_final,kru)
dun_final <- rbind(dun_final,dun)

check <- unique(kru_final$KEGG_code) #should be 49
check <- unique(dun_final$KEGG_code) #should be 48

rm(kru,dun)

#write.table (kru_final, "S5.5_kruskal_kegg_between_diets_per_TP.txt", row.names=FALSE,sep = "\t")
#write.table (dun_final, "S5.5_dunn_kegg_between_diets_per_TP.txt", row.names=FALSE,sep = "\t")

# plotting 4 months butyrate ====
#selecting only butyrate
rm(kru_final,count,tp,check,double_kegg_codes)
diet_but <- diet %>% filter(scfa == "butyrate")
diet_but <- diet_but %>% filter(Age_individual == "4 months")

but_dun <- dun_final %>% filter(scfa == "butyrate")
but4_dunn <- but_dun %>% filter(Age_individual == "4 months")

plyr::count(but4_dunn$KEGG_code)

diet_but_1 <- diet_but %>% filter(KEGG_code == "K00074")
diet_but_3 <- diet_but %>% filter(KEGG_code == "K00929")
diet_but_5 <- diet_but %>% filter(KEGG_code == "K14534")
diet_but_6 <- diet_but %>% filter(KEGG_code == "K17865")
diet_but_7 <- diet_but %>% filter(KEGG_code == "K18014")

dunn_but_kegg_1 <- but4_dunn %>% filter(KEGG_code == "K00074")
dunn_but_kegg_3 <- but4_dunn %>% filter(KEGG_code == "K00929")
dunn_but_kegg_5 <- but4_dunn %>% filter(KEGG_code == "K14534")
dunn_but_kegg_6 <- but4_dunn %>% filter(KEGG_code == "K17865")
dunn_but_kegg_7 <- but4_dunn %>% filter(KEGG_code == "K18014")

plot1 <- ggplot(diet_but_1, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.15, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle(NULL) +
  facet_grid( ~ KEGG_code) +
  xlab("Diet class") +
  ylab(NULL) + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot1 <- plot1 + stat_pvalue_manual(dunn_but_kegg_1,y.position = 4, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)
plot1

plot3 <- ggplot(diet_but_3, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.15, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL) +theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot3 <- plot3 + stat_pvalue_manual(dunn_but_kegg_3,y.position = 4.5, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot5 <- ggplot(diet_but_5, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.15, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL)  + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot5 <- plot5 + stat_pvalue_manual(dunn_but_kegg_5,y.position = 4, step.increase = 0.0001,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot6 <- ggplot(diet_but_6, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.15, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle("Abundance of KEGG genes for butyrate production at 4 months") +
  xlab("Diet class") +
  ylab("log10(read counts)")  + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        plot.title = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot6 <- plot6 + stat_pvalue_manual(dunn_but_kegg_6,y.position = 4, step.increase = 0.0001,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

plot7 <- ggplot(diet_but_7, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.15, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab("log10(read counts)")  + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot7 <- plot7 + stat_pvalue_manual(dunn_but_kegg_7,y.position = 4, step.increase = 0.001,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)


combined_plot <- (plot6 + plot5 + plot1) / (plot7 + plot3) +
  plot_layout() & 
  theme(plot.margin = margin(1, 0, 2, 0)) #top, right, bottom, left
combined_plot

# plotting 4 months butyrate presentation ====
#selecting only butyrate
plot1 <- ggplot(diet_but_1, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.5, alpha = 0.7,
    color = "black", shape = 16
  )+
  ggtitle(NULL) +
  facet_grid( ~ KEGG_code) +
  xlab("Diet class") +
  ylab(NULL) + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 12, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot1 <- plot1 + stat_pvalue_manual(dunn_but_kegg_1,y.position = 4, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE,size = 10)
plot1

plot3 <- ggplot(diet_but_3, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.5, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL) +theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 12, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot3 <- plot3 + stat_pvalue_manual(dunn_but_kegg_3,y.position = 4.5, step.increase = 0.1,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE,size = 10)

plot5 <- ggplot(diet_but_5, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.5, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab(NULL)  + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 12, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot5 <- plot5 + stat_pvalue_manual(dunn_but_kegg_5,y.position = 4, step.increase = 0.0001,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE,size = 10)

plot6 <- ggplot(diet_but_6, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.5, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle("Abundance of KEGG genes for butyrate production at 4 months") +
  xlab("Diet class") +
  ylab("log10(read counts)")  + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 12, margin = margin(b = 1)),
        plot.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot6 <- plot6 + stat_pvalue_manual(dunn_but_kegg_6,y.position = 4, step.increase = 0.0001,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE,size = 10)

plot7 <- ggplot(diet_but_7, aes(x = diet_class, y = log_read_counts)) +
  geom_boxplot(aes(fill = diet_class), outlier.shape = NA) + 
  geom_point(
    position = position_jitter(width = 0.15, height = 0),           # jitter to avoid overlap
    size = 0.5, alpha = 0.7,
    color = "black", shape = 16
  )+
  facet_grid( ~ KEGG_code) +
  ggtitle(NULL) +
  xlab("Diet class") +
  ylab("log10(read counts)")  + theme_minimal()+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 17), 
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 12, margin = margin(b = 1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot7 <- plot7 + stat_pvalue_manual(dunn_but_kegg_7,y.position = 4, step.increase = 0.001,tip.length = 0,
                                    label = "sig_adj",hide.ns = TRUE)

combined_plot <- (plot6 + plot5 + plot1) / (plot7 + plot3) +
  plot_layout() & 
  theme(plot.margin = margin(1, 0, 2, 0)) #top, right, bottom, left
combined_plot

# saving figure for the paper ====
ggsave(
  filename   = "S5.5_scfa_production.svg",
  plot       = combined_plot,
  width    = 80,
  height   = 90,
  units    = "mm",
  device   = svglite
)

