# load in packages ====
library(tidyverse) 
library(phyloseq)
library(microViz)
library(plyr)
library(reshape2)
library(dplyr)
library(readxl)
library(LinDA)
library(ComplexHeatmap)
library(patchwork)
library(gridExtra)
library(svglite)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# loading in MS phyloseq ====
ps_diet <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

#cleaning data
ps_diet <- ps_diet  %>%
  ps_mutate(
    timepoint_number = factor(str_remove(ps_diet@sam_data$timepoint,"m"), levels = c("4","5","6","9","11","14")),
  )

# filtering abundance > 0.01 ====
ps_diet <- phyloseq_validate(ps_diet, remove_undetected = TRUE)
ps_diet_comp <- ps_diet %>% tax_transform("compositional") %>% ps_get()

#filter per TP
ps_diet_4 <- ps_diet_comp %>% ps_filter(timepoint == "4m") 
ps_diet_5 <- ps_diet_comp %>% ps_filter(timepoint == "5m")
ps_diet_6 <- ps_diet_comp %>% ps_filter(timepoint == "6m") 
ps_diet_9 <- ps_diet_comp %>% ps_filter(timepoint == "9m") 
ps_diet_11 <- ps_diet_comp %>% ps_filter(timepoint == "11m") 
ps_diet_14 <- ps_diet_comp %>% ps_filter(timepoint == "14m")

#extract OTU table
otu_4 <- ps_diet_4 %>% otu_table() %>% as.data.frame()
otu_5 <- ps_diet_5 %>% otu_table() %>% as.data.frame()
otu_6 <- ps_diet_6 %>% otu_table() %>% as.data.frame()
otu_9 <- ps_diet_9 %>% otu_table() %>% as.data.frame()
otu_11 <- ps_diet_11 %>% otu_table() %>% as.data.frame()
otu_14 <- ps_diet_14 %>% otu_table() %>% as.data.frame()
rm(ps_diet_4,ps_diet_5,ps_diet_6,ps_diet_9,ps_diet_11,ps_diet_14,ps_diet_comp)

#calculate row mean per taxa
otu_4$mean <- rowMeans(otu_4)
otu_5$mean <- rowMeans(otu_5)
otu_6$mean <- rowMeans(otu_6)
otu_9$mean <- rowMeans(otu_9)
otu_11$mean <- rowMeans(otu_11)
otu_14$mean <- rowMeans(otu_14)

#select those bakut that are above 0.01
bakut_4 <- otu_4 %>% filter(mean >= 0.01) %>% rownames()
bakut_5 <- otu_5 %>% filter(mean >= 0.01) %>% rownames()
bakut_6 <- otu_6 %>% filter(mean >= 0.01) %>% rownames()
bakut_9 <- otu_9 %>% filter(mean >= 0.01) %>% rownames()
bakut_11 <- otu_11 %>% filter(mean >= 0.01) %>% rownames()
bakut_14 <- otu_14 %>% filter(mean >= 0.01) %>% rownames()
rm(otu_4,otu_5,otu_6,otu_9,otu_11,otu_14)

#combine into 1 vector
unik_baku <- c(bakut_4,bakut_5,bakut_6,bakut_9,bakut_11,bakut_14)
unik_baku <- unique(unik_baku)
rm(bakut_4,bakut_5,bakut_6,bakut_9,bakut_11,bakut_14)

#remove everything below 1% from phyloseq
ps_diet_sel <- ps_diet %>% tax_select(unik_baku)
rm(unik_baku)

# A.Class 1 to others ====
# A.1 LINDA preparations  ====
num_samples <- nsamples(ps_diet_sel)
ind_l <- rep(TRUE, num_samples)
otu_DAA <- as.data.frame(ps_diet_sel@otu_table[, ind_l])
meta_DAA <- cbind.data.frame(Diet = factor(ps_diet_sel@sam_data$diet_class[ind_l]),
                             Timepoint = factor(ps_diet_sel@sam_data$timepoint_number[ind_l]),
                             Solids = as.numeric(ps_diet_sel@sam_data$ageintrosolids[ind_l]),
                             Birth_weigth = as.numeric(ps_diet_sel@sam_data$birth_weigth[ind_l]),
                             Breastfeeding = as.numeric(ps_diet_sel@sam_data$duration_breastfeeding[ind_l]),
                             Place = factor(ps_diet_sel@sam_data$delivery_place_word[ind_l]),
                             Delivery = factor(ps_diet_sel@sam_data$delivery_type_bi_word[ind_l]),
                             Kindcode = factor(ps_diet_sel@sam_data$Sample_aliquote[ind_l]),
                             Sex = factor(ps_diet_sel@sam_data$sex_word[ind_l]))

meta_DAA$Timepoint<- factor(meta_DAA$Timepoint, 
                           levels = c("4",
                                      "5",
                                      "6",
                                      "9",
                                      "11",
                                      "14"
                           ))

meta_DAA$Delivery<- factor(meta_DAA$Delivery, 
                           levels = c("vaginal",
                                      "caesarian section"
                           ))

meta_DAA$Sex<- factor(meta_DAA$Sex, 
                      levels = c("male",
                                 "female"
                      ))

meta_DAA$Place<- factor(meta_DAA$Place, 
                        levels = c("hospital",
                                   "at home"
                        ))

meta_DAA <- separate(meta_DAA, Kindcode, into = c("ID", "musor1","musor2"), sep = "\\.")

# A.2 running LINDA ====
linda.obj_inter <- linda(otu_DAA, meta_DAA, formula = '~Diet*Timepoint+Solids+Birth_weigth+Breastfeeding+Place+Delivery+Sex+(1|ID)', alpha = 0.05,
                         prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#output
output_inter <- bind_rows(linda.obj_inter$output, .id = "variable")
output_inter <- rownames_to_column(output_inter, var = "bacteria")
output_inter$bacteria <- gsub(".*s__", "", output_inter$bacteria)
output_inter$bacteria <- sub("\\.\\.\\..*", "", output_inter$bacteria)
output_inter$bacteria <- gsub("_", " ", output_inter$bacteria)

output_round_inter <- output_inter %>% 
  mutate(sig_padj_02 = if_else(
    condition = padj < 0.2 & padj >= 0.05, 
    true = "*", 
    false = if_else(padj < 0.05, "**", "")
  ))

output_round_inter <- output_round_inter %>% 
  mutate(sig_padj_01 = if_else(
    condition = padj < 0.1 & padj >= 0.05, 
    true = "*", 
    false = if_else(padj < 0.05, "**", "")
  ))

output_round_inter <- output_round_inter %>% 
  mutate(sig_padj_005 = if_else(padj < 0.05, "*", ""))

output_round_inter <- output_round_inter %>% 
  mutate(sig_pvalue = if_else(
    condition = pvalue <= 0.05, 
    true = "*", 
    false = "")
  )

output_round_inter <- output_round_inter %>% 
  mutate(sig_pvalue = if_else(
    condition = pvalue <= 0.05, 
    true = "*", 
    false = "")
  )

output_round_inter <- output_round_inter %>% 
  mutate(sig_padj = if_else(
    condition = padj <= 0.05, 
    true = "*", 
    false = "")
  )

# Tighten the threshold & show robustness
output_round_inter <- output_round_inter %>% 
  mutate(q_value = if_else(
    condition = padj < 0.2 & padj >= 0.1, 
    true = "Suggestive", 
    false = if_else(condition = padj < 0.1 & padj > 0.05, 
                    true = "Moderate", 
                    false = if_else(condition = padj < 0.05, 
                                    true = "Strong", 
                                    false = ""))))

# A.3 cleaning for plot ====
plyr::count(output_round_inter$variable)

output_round_inter$variable = gsub("Birth_weigth", "Birth weight (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Breastfeeding", "Duration of breastfeeding (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Deliverycaesarian section", "C-section (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet2:Timepoint5", "5m*Class2 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet2:Timepoint6", "6m*Class2 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet2:Timepoint9", "9m*Class2 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet2:Timepoint11", "11m*Class2 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet2:Timepoint14", "14m*Class2 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint5", "5m*Class3 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint6", "6m*Class3 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint9", "9m*Class3 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint11", "11m*Class3 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint14", "14m*Class3 (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet2", "Class 2 (ref.1)", output_round_inter$variable)#must be after interaction term
output_round_inter$variable = gsub("Diet3", "Class 3 (ref.1)", output_round_inter$variable)#must be after interaction term
output_round_inter$variable = gsub("Placeat home", "Delivery at home (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Sexfemale", "Female (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Solids", "Age solids introduced (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint5", "5m (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint6", "6m (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint9", "9m (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint11", "11m (ref.1)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint14", "14m (ref.1)", output_round_inter$variable)

plyr::count(output_round_inter$variable)

output_round_inter$variable<- factor(output_round_inter$variable, 
                                     levels = c("5m (ref.1)",
                                                "6m (ref.1)",
                                                "9m (ref.1)",
                                                "11m (ref.1)",
                                                "14m (ref.1)",
                                                "Age solids introduced (ref.1)",
                                                "Birth weight (ref.1)",
                                                "C-section (ref.1)",
                                                "Delivery at home (ref.1)",
                                                "Duration of breastfeeding (ref.1)",
                                                "Female (ref.1)",
                                                "Class 2 (ref.1)",
                                                "Class 3 (ref.1)",
                                                "5m*Class2 (ref.1)",
                                                "6m*Class2 (ref.1)",
                                                "9m*Class2 (ref.1)",
                                                "11m*Class2 (ref.1)",
                                                "14m*Class2 (ref.1)",
                                                "5m*Class3 (ref.1)",
                                                "6m*Class3 (ref.1)",
                                                "9m*Class3 (ref.1)",
                                                "11m*Class3 (ref.1)",
                                                "14m*Class3 (ref.1)"
                                     ))

plyr::count(output_round_inter$variable)

result_1to_other <- output_round_inter

rm(otu_DAA,output_inter,linda.obj_inter,output_round_inter,meta_DAA,num_samples,ind_l)

# B.Class 2 to others ====
# B.1 LINDA preparations ====
num_samples <- nsamples(ps_diet_sel)
ind_l <- rep(TRUE, num_samples)
otu_DAA <- as.data.frame(ps_diet_sel@otu_table[, ind_l])
meta_DAA <- cbind.data.frame(Diet = factor(ps_diet_sel@sam_data$diet_class[ind_l]),
                             Timepoint = factor(ps_diet_sel@sam_data$timepoint_number[ind_l]),
                             Solids = as.numeric(ps_diet_sel@sam_data$ageintrosolids[ind_l]),
                             Birth_weigth = as.numeric(ps_diet_sel@sam_data$birth_weigth[ind_l]),
                             Breastfeeding = as.numeric(ps_diet_sel@sam_data$duration_breastfeeding[ind_l]),
                             Place = factor(ps_diet_sel@sam_data$delivery_place_word[ind_l]),
                             Delivery = factor(ps_diet_sel@sam_data$delivery_type_bi_word[ind_l]),
                             Kindcode = factor(ps_diet_sel@sam_data$Sample_aliquote[ind_l]),
                             Sex = factor(ps_diet_sel@sam_data$sex_word[ind_l]))

meta_DAA$Timepoint<- factor(meta_DAA$Timepoint, 
                            levels = c("4",
                                       "5",
                                       "6",
                                       "9",
                                       "11",
                                       "14"
                            ))

meta_DAA$Delivery<- factor(meta_DAA$Delivery, 
                           levels = c("vaginal",
                                      "caesarian section"
                           ))

meta_DAA$Sex<- factor(meta_DAA$Sex, 
                      levels = c("male",
                                 "female"
                      ))

meta_DAA$Place<- factor(meta_DAA$Place, 
                        levels = c("hospital",
                                   "at home"
                        ))

meta_DAA <- separate(meta_DAA, Kindcode, into = c("ID", "musor1","musor2"), sep = "\\.")

meta_DAA$Diet<- factor(meta_DAA$Diet, 
                       levels = c("2",
                                  "3",
                                  "1"
                       ))
# B.2 running LINDA ====
linda.obj_inter <- linda(otu_DAA, meta_DAA, formula = '~Diet*Timepoint+Solids+Birth_weigth+Breastfeeding+Place+Delivery+Sex+(1|ID)', alpha = 0.05,
                         prev.cut = 0, lib.cut = 1, winsor.quan = NULL, corr.cut = 0.1, p.adj.method = 'BH',type = "count")

#output
output_inter <- bind_rows(linda.obj_inter$output, .id = "variable")
output_inter <- rownames_to_column(output_inter, var = "bacteria")
output_inter$bacteria <- gsub(".*s__", "", output_inter$bacteria)
output_inter$bacteria <- sub("\\.\\.\\..*", "", output_inter$bacteria)
output_inter$bacteria <- gsub("_", " ", output_inter$bacteria)

output_round_inter <- output_inter %>% 
  mutate(sig_padj_02 = if_else(
    condition = padj < 0.2 & padj >= 0.05, 
    true = "*", 
    false = if_else(padj < 0.05, "**", "")
  ))

output_round_inter <- output_round_inter %>% 
  mutate(sig_padj_01 = if_else(
    condition = padj < 0.1 & padj >= 0.05, 
    true = "*", 
    false = if_else(padj < 0.05, "**", "")
  ))

output_round_inter <- output_round_inter %>% 
  mutate(sig_padj_005 = if_else(padj < 0.05, "*", ""))

output_round_inter <- output_round_inter %>% 
  mutate(sig_pvalue = if_else(
    condition = pvalue <= 0.05, 
    true = "*", 
    false = "")
  )

output_round_inter <- output_round_inter %>% 
  mutate(sig_padj = if_else(
    condition = padj <= 0.05, 
    true = "*", 
    false = "")
  )

# Tighten the threshold & show robustness
output_round_inter <- output_round_inter %>% 
  mutate(q_value = if_else(
    condition = padj < 0.2 & padj >= 0.1, 
    true = "Suggestive", 
    false = if_else(condition = padj < 0.1 & padj > 0.05, 
                    true = "Moderate", 
                    false = if_else(condition = padj < 0.05, 
                                    true = "Strong", 
                                    false = ""))))

# B.3 cleaning for plot ====
plyr::count(output_round_inter$variable)

output_round_inter$variable = gsub("Birth_weigth", "Birth weight (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Breastfeeding", "Duration of breastfeeding (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Deliverycaesarian section", "C-section (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet1:Timepoint5", "5m*Class1 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet1:Timepoint6", "6m*Class1 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet1:Timepoint9", "9m*Class1 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet1:Timepoint11", "11m*Class1 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet1:Timepoint14", "14m*Class1 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint5", "5m*Class3 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint6", "6m*Class3 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint9", "9m*Class3 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint11", "11m*Class3 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet3:Timepoint14", "14m*Class3 (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Diet1", "Class 1 (ref.2)", output_round_inter$variable) #must be after interaction term
output_round_inter$variable = gsub("Diet3", "Class 3 (ref.2)", output_round_inter$variable) #must be after interaction term
output_round_inter$variable = gsub("Placeat home", "Delivery at home (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Sexfemale", "Female (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Solids", "Age solids introduced (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint5", "5m (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint6", "6m (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint9", "9m (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint11", "11m (ref.2)", output_round_inter$variable)
output_round_inter$variable = gsub("Timepoint14", "14m (ref.2)", output_round_inter$variable)

plyr::count(output_round_inter$variable)

output_round_inter$variable<- factor(output_round_inter$variable, 
                                     levels = c("5m (ref.2)",
                                                "6m (ref.2)",
                                                "9m (ref.2)",
                                                "11m (ref.2)",
                                                "14m (ref.2)",
                                                "Age solids introduced (ref.2)",
                                                "Birth weight (ref.2)",
                                                "C-section (ref.2)",
                                                "Delivery at home (ref.2)",
                                                "Duration of breastfeeding (ref.2)",
                                                "Female (ref.2)",
                                                "Class 1 (ref.2)",
                                                "Class 3 (ref.2)",
                                                "5m*Class1 (ref.2)",
                                                "6m*Class1 (ref.2)",
                                                "9m*Class1 (ref.2)",
                                                "11m*Class1 (ref.2)",
                                                "14m*Class1 (ref.2)",
                                                "5m*Class3 (ref.2)",
                                                "6m*Class3 (ref.2)",
                                                "9m*Class3 (ref.2)",
                                                "11m*Class3 (ref.2)",
                                                "14m*Class3 (ref.2)"
                                     ))

rm(otu_DAA,output_inter,linda.obj_inter)

result_2to_other <- output_round_inter

# merging and saving results ====
#uniting for the paper all results
all <- rbind(result_1to_other,result_2to_other)
#write.table (all, "S2.1_LINDA_results_all.txt", row.names=FALSE,sep = "\t")

# 6.cleaning for plot ====
#cleaning LINDA 1 to others
strains_sig <- result_1to_other[grepl("\\*|\\*\\*", result_1to_other$sig_pvalue), ]
strains_vector <- strains_sig$bacteria #selecting only significant strains
strains_vector <- strains_vector %>% unique() %>% sort()
result_1to_other_subset <- result_1to_other[result_1to_other$bacteria %in% strains_vector, ]

#cleaning LINDA 2 to others
result_2to_other_subset <- result_2to_other[result_2to_other$bacteria %in% strains_vector, ]

all <- rbind(result_1to_other_subset,result_2to_other_subset)
all <- filter(all, !(variable %in% c("Age solids introduced (ref.2)",
                                     "Birth weight (ref.2)",
                                     "C-section (ref.2)",
                                     "Delivery at home (ref.2)",
                                     "Duration of breastfeeding (ref.2)",
                                     "Female (ref.2)",
                                     "Class 1 (ref.2)",
                                     "5m*Class1 (ref.2)",
                                     "6m*Class1 (ref.2)",
                                     "9m*Class1 (ref.2)",
                                     "11m*Class1 (ref.2)",
                                     "14m*Class1 (ref.2)",
                                     "5m (ref.2)",
                                     "6m (ref.2)",
                                     "9m (ref.2)",
                                     "11m (ref.2)",
                                     "14m (ref.2)")))

all$variable<- factor(all$variable, 
                      levels = c("5m (ref.1)",
                                 "6m (ref.1)",
                                 "9m (ref.1)",
                                 "11m (ref.1)",
                                 "14m (ref.1)",
                                 "5m (ref.2)",
                                 "6m (ref.2)",
                                 "9m (ref.2)",
                                 "11m (ref.2)",
                                 "14m (ref.2)",
                                 "Age solids introduced (ref.1)",
                                 "Birth weight (ref.1)",
                                 "C-section (ref.1)",
                                 "Delivery at home (ref.1)",
                                 "Duration of breastfeeding (ref.1)",
                                 "Female (ref.1)",
                                 "Class 2 (ref.1)",
                                 "Class 3 (ref.1)",
                                 "Class 3 (ref.2)",
                                 "5m*Class2 (ref.1)",
                                 "6m*Class2 (ref.1)",
                                 "9m*Class2 (ref.1)",
                                 "11m*Class2 (ref.1)",
                                 "14m*Class2 (ref.1)",
                                 "5m*Class3 (ref.1)",
                                 "6m*Class3 (ref.1)",
                                 "9m*Class3 (ref.1)",
                                 "11m*Class3 (ref.1)",
                                 "14m*Class3 (ref.1)",
                                 "5m*Class3 (ref.2)",
                                 "6m*Class3 (ref.2)",
                                 "9m*Class3 (ref.2)",
                                 "11m*Class3 (ref.2)",
                                 "14m*Class3 (ref.2)"
                      ))

all <- all %>%
  mutate(
    q_stars = case_when(
      q_value == "Strong" ~ "***",
      q_value == "Moderate" ~ "**",
      q_value == "Suggestive" ~ "*",
      TRUE ~ ""                       # blank for "empty"
    ),
    q_value = factor(q_value, levels = c("Strong","Moderate","Suggestive"))
  )

# 7.Panel 3 A ====
panel_3A_with <- all %>%
  ggplot(aes(x = reorder(variable, desc(variable)), y = bacteria )) +
  geom_tile(aes(fill = log2FoldChange, height = .5, width = .9)) +
  scale_y_discrete(limits = sort(unique(all$bacteria)), position = "right") +
  coord_flip() +
  theme_minimal() +
  scale_fill_gradient2(
    low = "#0072B2", high = "#CD0000", mid = "white", midpoint = 0,
    limits = c(-13, 20),
    guide = guide_colorbar(barheight = 4, barwidth = 0.3)
  ) +
  # 2) Draw the stars and map color to create a legend
  geom_text(aes(label = q_stars, color = q_value), vjust = 0.7, size = 2, na.rm = TRUE, show.legend = TRUE) +
  # 3) Manual legend for the stars
  scale_color_manual(
    name = "Significance",
    values = c(Strong = "black", Moderate = "black", Suggestive = "black"),
    labels = c("Strong", "Moderate", "Suggestive"),
    na.translate = FALSE
  ) +
  guides(
    fill  = guide_colorbar(order = 1, title = "Log2 fold \nchange \nin abundance"),
    color = guide_legend(order = 2, override.aes = list(label = c("***","**","*"), size = 4))
  ) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = -0.03, size = 7),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title= element_text(size = 6),
    legend.box  = "vertical" # stack colorbar and stars legend
  ) +
  geom_vline(xintercept = c(29.5,24.5,18.5,15.5,10.5,5.5), linetype = "solid") +
  labs(y = "", x = "")

panel_3A_with

panel_3A_without <- all %>%
  ggplot(aes(x = reorder(variable, desc(variable)), y = bacteria )) +
  geom_tile(aes(fill = log2FoldChange, height = .5, width = .9)) +
  scale_y_discrete(limits = sort(unique(all$bacteria)), position = "right") +
  coord_flip() +
  theme_minimal() +
  scale_fill_gradient2(
    low = "#0072B2", high = "#CD0000", mid = "white", midpoint = 0,
    limits = c(-13, 20),
    guide = "none"          # remove fill legend
  ) +
  geom_text(aes(label = q_stars, color = q_value),
            vjust = 0.7, size = 2, na.rm = TRUE, show.legend = FALSE) +  # no color legend
  scale_color_manual(
    values = c(Strong = "black", Moderate = "black", Suggestive = "black"),
    na.translate = FALSE
  ) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = -0.03, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "none",   # make sure no legends are drawn
    plot.margin = margin(t = 5.5, r = 35, b = 5.5, l = 5.5, unit = "pt") ) +
  geom_vline(xintercept = c(29.5,24.5,18.5,15.5,10.5,5.5), linetype = "solid") +
  labs(y = "", x = "")

panel_3A_without

# saving figure for the paper ====
ggsave(
  filename   = "S2.1_linda_heatmap_with.svg",
  plot       = panel_3A_with ,
  width    = 170,
  height   = 180,
  units    = "mm",
  device   = svglite
)

ggsave(
  filename   = "S2.1_linda_heatmap_without.svg",
  plot       = panel_3A_without ,
  width    = 170,
  height   = 160,
  units    = "mm",
  device   = svglite
)

