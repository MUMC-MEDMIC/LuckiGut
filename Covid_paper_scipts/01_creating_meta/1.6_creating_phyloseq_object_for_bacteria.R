# load libraries ====
#how to convert metagenomic file to phyloseq tables, https://mvuko.github.io/meta_phyloseq/
#meta should have the same amount of obs. as otu has variables!

library(tidyverse)
library(phyloseq)
library(microViz)
library(vegan)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
otu_raw = read.delim("Output/Files/S1.5_OTU_bacteria.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
meta_raw = read.delim("Output/Files/S1.4_full_meta.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
read_count = read.delim("Output/Files/S1.5_read_counts.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# OTU ====
#check structure
str(otu_raw) #OTU ID's have to be the rownames
#change the names of the rows in otu from whole string to only species names to make it easy to analyse
otu <- tibble::rownames_to_column(otu_raw, "taxa")
otu[c('k','p','c','o','f','g','x')]<- str_split_fixed(otu$taxa, ";", n = 7)
otu[c('z','s')]<- str_split_fixed(otu$x, "__", n = 2)
rownames(otu) <- otu$s
otu = subset(otu, select = -c(k,p,c,o,f,g,z,x,s,taxa))

# TAXA ====
taxa_raw <- otu_raw %>% tibble::rownames_to_column (var = "SampleID") %>% select(SampleID) %>% 
  separate(SampleID, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),";")  
#to clean taxa from garbage
taxa <- data.frame(row.names = row.names(taxa_raw),
                         Kingdom = str_replace(taxa_raw[,1], "k__",""),
                         Phylum = str_replace(taxa_raw[,2], "p__",""),
                         Class = str_replace(taxa_raw[,3], "c__",""),
                         Order = str_replace(taxa_raw[,4], "o__",""),
                         Family = str_replace(taxa_raw[,5], "f__",""),
                         Genus = str_replace(taxa_raw[,6], "g__",""),
                         Species = str_replace(taxa_raw[,7], "s__",""),
                         stringsAsFactors = FALSE)
#check structure
str(taxa) 
#the output is a data frame of characters, and we need taxa to be recognized as factors
taxa<- taxa %>% mutate_if(is.character, as.factor)
#change the names of the rows in taxa to only species names to make it easy to analyse
rownames(taxa) <- taxa$Species
#write taxa file down
#write.table (taxa, "S4_105_taxa_bacteria.txt", row.names=TRUE,sep = "\t")

# META ====
#adding read counts to meta
meta_raw_counts = full_join(meta_raw,read_count, by = "MMHP_SampleID")
#the output is a data frame of characters, and we need taxa to be recognized as factors
meta<- meta_raw_counts %>% mutate_if(is.character, as.factor)
#correct order or type of variable
meta$Age_individual<- factor(meta$Age_individual, levels = c("1-2 weeks", "4 weeks", "8 weeks","4 months", "5 months", "6 months", "9 months", "11 months", "14 months", "mother", "father"))
sapply(meta,class)#check what class are the variables

# adding ENS, shannon, simpson, inversimpson to the data to the data ====
#calculate shannon index
ens_sha <- otu %>% t() %>% diversity(index = "shannon")%>% as.data.frame()
#calculates the exponential value of shannon index
names(ens_sha)[names(ens_sha) == "."] <- "shannon_index"
ens_sha <- ens_sha %>% mutate(exp(shannon_index))
names(ens_sha)[names(ens_sha) == "exp(shannon_index)"] <- "ENS"
#make the row names into column
ens_sha <- tibble::rownames_to_column(ens_sha, "MMHP_SampleID")
#calculate simpson
simpson <- otu %>% t() %>% diversity(index = "simpson")%>% as.data.frame()
names(simpson)[names(simpson) == "."] <- "simpson_index"
simpson <- tibble::rownames_to_column(simpson, "MMHP_SampleID")#make the row names into column
#calculate invers-simpson
invsimpson <- otu %>% t() %>% diversity(index = "invsimpson")%>% as.data.frame()
names(invsimpson)[names(invsimpson) == "."] <- "invsimpson_index"
invsimpson <- tibble::rownames_to_column(invsimpson, "MMHP_SampleID")#make the row names into column
#spieces observed richness
richness <- otu %>% t() %>% specnumber()%>% as.data.frame()
names(richness)[names(richness) == "."] <- "observed_richness"
richness <- tibble::rownames_to_column(richness, "MMHP_SampleID")#make the row names into column
#merge indexes 
indexes <- inner_join(ens_sha, simpson, by = "MMHP_SampleID")
indices <- inner_join(richness,invsimpson, by = "MMHP_SampleID")
ind_final <- inner_join(indexes,indices, by = "MMHP_SampleID")
#merge meat with calculated ENS table
meta <- inner_join(ind_final, meta, by = "MMHP_SampleID")
rm(ens_sha,indexes,indices,invsimpson,simpson,ind_final,richness)
#make sample names row names
rownames(meta)<- meta$MMHP_SampleID

# creating phyloseq object ====
#change the names of the rows in taxa and otu from whole string to only species names to make it easy to analyse
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)
#transform data to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(meta)
rm(taxa_raw,otu_raw,meta_raw)
#and put them in one object
phylo_object_bacteria<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)
#check if ok
sample_sums(phylo_object_bacteria) 
sample_names(phylo_object_bacteria) 
rank_names(phylo_object_bacteria)        
rm(phylo_OTU,phylo_TAX,otu_mat,phylo_samples,tax_mat)
rm(taxa,meta,otu)
#adding observed richness to the data, chao1 had the same value sas observed_richness
phylo_object_bacteria <- ps_calc_richness(
  phylo_object_bacteria,
  "Species",
  index = "observed",
  detection = 0
)

# save phyloseq object ====
#saveRDS(phylo_object_bacteria, 'S1.6_phyloseq_object_bacteria.rds')

# filtering species ====
# filtering with 5% prevalence
baku_filter_species <- phylo_object_bacteria %>% tax_filter(min_prevalence = 0.05, tax_level = "Species")
#saveRDS(baku_filter_species, 'S1.6_phyloseq_object_bacteria_filtered.rds')

#check that species are different
#raw <- as.data.frame(phylo_object_bacteria@otu_table) #728 obs
#prevalence <- as.data.frame(baku_filter_species@otu_table) #216 obs
#rm(raw, prevalence)

# filtering by individual ====
baku_f_I.M <- baku_filter_species %>% ps_filter(Individual != "father")
#check that number of samples is different
#all <- as.data.frame(baku_filter_species@sam_data) #930 obs
#m_i <- as.data.frame(baku_f_I.M@sam_data) #924 obs
#rm(all, m_i)
#saveRDS(baku_f_I.M, 'S1.6_phyloseq_object_bacteria_filtered_mother_infant.rds')

baku_f_I <- baku_f_I.M %>% ps_filter(Individual != "mother")
#check that number of samples is different
#m_i <- as.data.frame(baku_f_I.M@sam_data) #924 obs
#only_i <- as.data.frame(baku_f_I@sam_data) #808 obs
#rm(m_i, only_i)

#saveRDS(baku_f_I, 'S1.6_phyloseq_object_bacteria_filtered_infant.rds')
rm(baku_f_I.M,baku_filter_species,phylo_object_bacteria)









