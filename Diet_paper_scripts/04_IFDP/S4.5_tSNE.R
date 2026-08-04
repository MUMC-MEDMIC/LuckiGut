# LOAD LIBRARIES ====
library(ggplot2)
library(tidyverse)
library(Rtsne)
library(phyloseq)
library(microViz)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ifdp_raw <- read.csv(file = "Input/6_combined_counts_IFDP.csv")
my_phyloseq <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# clean fiber data ====
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

# cleaning in phyloseq ====
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

# data cleaning ====
ifdp_clean <- ifdp %>% select(-c("Age_individual","diet_class","MMHP_SampleID","Sample_aliquote"))
tsne <- Rtsne(as.data.frame(ifdp_clean),dims=3, perplexity=50)

# t-SNE, Supplementary Figure 12  ====
set.seed(100)

plot <- ggplot(data = as.data.frame(tsne$Y) %>% 
         mutate(diet_class = ifdp$diet_class,
                timepoint = ifdp$Age_individual),
       aes(x = V1, y = V2)) + 
  geom_point(aes(color = timepoint, shape = diet_class)) +
  #geom_mark_ellipse(aes(color = timepoint), expand = unit(0.5,"mm"), alpha = 0.2)  +
  scale_color_manual("Time point", labels = c("4m","5m","6m","9m","11m","14m"),
                     values = c("#BFD200","#AACC00", "#80B918",  "#55A630", "#2B9348" ,"#007F5F")) +
  scale_shape_manual("Dietary class", values = c(19,4,17)) +
  labs(x = "t-SNE1", y = "t-SNE 2") +
  theme_bw()

# saving figure for the paper ====
ggsave(
  filename   = "S4.5_tSNE.tiff",
  plot       = plot,
  width      = 225,      # full page width
  height     = 170,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)


