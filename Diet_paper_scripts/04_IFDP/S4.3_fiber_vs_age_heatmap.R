# LOAD LIBRARIES ====
library(readxl) 
library(tidyverse)
library(phyloseq)
library(microViz)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)# for colorRamp2

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ifdp_raw <- read.csv(file = "Input/6_combined_counts_IFDP.csv")
my_phyloseq <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# cleaning data ====
#read and clean fiber data
rownames(ifdp_raw) <- ifdp_raw$Fibers #make fibers into rownames
ifdp_raw <- ifdp_raw[, -1] #remove fiber column
ifdp_raw <- as.data.frame(t(ifdp_raw)) #make transposon
ifdp_raw$MMHP_SampleID <- rownames(ifdp_raw) #make rownames into a column
ifdp_raw$MMHP_SampleID <- gsub("[.-]", "_", ifdp_raw$MMHP_SampleID) #make the names of the sample the same format

#renaiming some fibers to avoid problems in the future and spelling mistakes
names(ifdp_raw)[names(ifdp_raw) == "Beta-glucan"] <- "Betaglucan"
names(ifdp_raw)[names(ifdp_raw) == "Resistant starch"] <- "Resistant_starch"
names(ifdp_raw)[names(ifdp_raw) == "Xanthan "] <- "Xanthan"
names(ifdp_raw)[names(ifdp_raw) == "Galctomannan"] <- "Galactomannan"
names(ifdp_raw)[names(ifdp_raw) == "Carageenan "] <- "Carrageenan"
names(ifdp_raw)[names(ifdp_raw) == "Xloglucan"] <- "Xyloglucan"
ifdp_names <- setdiff(colnames(ifdp_raw), "MMHP_SampleID")

# clean phyloseq ====
my_phyloseq <- my_phyloseq %>% ps_join(ifdp_raw, by = "MMHP_SampleID") #adding ifdp to phyloseq
check <- as.data.frame(as.matrix(my_phyloseq@otu_table))

#taking metadata out
full_meta <- as.data.frame(as.matrix(my_phyloseq@sam_data))

#creating only fiber table for correlations
ifdp_names_extra <- c(ifdp_names, "Sample_aliquote","Age_individual","MMHP_SampleID","diet_class")
ifdp <- full_meta[, ifdp_names_extra, drop = FALSE]
convert_to_numeric <- function(df, columns) {
  df[columns] <- lapply(df[columns], function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
  })
  return(df)
}
ifdp <- convert_to_numeric(ifdp, ifdp_names)
rm(ifdp_raw,full_meta,ifdp_names_extra)

# filtering abundance > 0.01 ====
my_phyloseq <- phyloseq_validate(my_phyloseq, remove_undetected = TRUE)
ps_diet_comp <- my_phyloseq %>% tax_transform("compositional") %>% ps_get()

#filter per TP
ps_diet_4 <- ps_diet_comp %>% ps_filter(Age_individual == "4 months") 
ps_diet_5 <- ps_diet_comp %>% ps_filter(Age_individual == "5 months")
ps_diet_6 <- ps_diet_comp %>% ps_filter(Age_individual == "6 months") 
ps_diet_9 <- ps_diet_comp %>% ps_filter(Age_individual == "9 months") 
ps_diet_11 <- ps_diet_comp %>% ps_filter(Age_individual == "11 months") 
ps_diet_14 <- ps_diet_comp %>% ps_filter(Age_individual == "14 months")

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
my_phyloseq_filtered <- my_phyloseq %>% tax_select(unik_baku)
rm(unik_baku,my_phyloseq)

# cleaning data for the heatmap ====
# Create sample data (for each column, a metadata)
sample_data <- as.data.frame(as.matrix(sample_data(my_phyloseq_filtered)))
sample_data <- sample_data[, !colnames(sample_data) %in% ifdp_names]
col_to_keep <- c("MMHP_SampleID","Sample_aliquote","Age_individual", "diet_class")
sample_data <- sample_data[, colnames(sample_data) %in% col_to_keep]
sample_data$kindcode_copy <- sample_data$Sample_aliquote
sample_data <- separate(sample_data, kindcode_copy, into = c("id", "musor", "musor1"), sep = "\\.")
sample_data$musor<- NULL
sample_data$musor1<- NULL
sample_data$Age_individual <- gsub("months", "", sample_data$Age_individual)

#creating table for supplementary fiber per diet per TP
mean_df <- ifdp[, colnames(ifdp) %in% ifdp_names]
mean_df <- as.data.frame(t(mean_df))

df <- as.data.frame(t(mean_df)) %>% mutate(timepoint = sample_data$Age_individual, 
                                           diet_class = sample_data$diet_class) %>%
  group_by(timepoint, diet_class) %>%
  summarize_all(mean) %>%
  arrange(diet_class)

# Specify columns to exclude
exclude_columns <- c("timepoint", "diet_class")
# Sort the remaining columns alphabetically
sorted_columns <- sort(ifdp_names)
# Reorder the columns: first the excluded ones, then the sorted ones
df_sorted <- df[, c(exclude_columns, sorted_columns)]

#write.table (df_sorted, "S4_12_mean_DF_DC.txt", row.names=FALSE,sep = "\t")
rm(exclude_columns,sorted_columns)

# Heatmap of fibers, Supplementary Figure 13 ====
plot <- Heatmap(t(df[,3:26]), column_split = c(rep(1,6), rep(2,6), rep(3,6)),
        col = colorRamp2(c(0, 0.1), c("#ffffff", "#005824")),
        border = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(t(df[,3:26])[i, j], digits = 3), x, y, gp = gpar(fontsize = 8))
        },
        cluster_columns = FALSE, column_labels = rep(c("4","5","6","9","11","14"),3),
        column_names_gp = gpar(cex = 1), column_names_rot = 0,
        heatmap_legend_param = list(title = "Mean DF \ndegradation \ncapabilities"))

# saving figure for the paper ====
w_cm <- 22.5
h_cm <- 17

tiff("S4.3_Mean_df_degradation_capacity_per_TP.tiff",
     width  = w_cm / 2.54,  # inches
     height = h_cm / 2.54,
     units  = "in",
     res    = 300)

draw(plot)  # important for ComplexHeatmap
dev.off()
