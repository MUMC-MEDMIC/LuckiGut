# libraries ====
library(microViz)
library(vegan)
library(ggplot2)
library(gdata)
library(gtools)
library(DirichletMultinomial)
library(plyr)
library(ComplexHeatmap)
library(colorRamp2)
library(patchwork)
library(viridis)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in phyloseq data ====
phyloseq_big <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

# DMM species best model ====
#Dirichlet Multinomial Mixture (DMM) clustering, an unsupervised clustering method that uses Laplace approximation to 
#identify groups of communities (enterotypes) with similar composition
#Species table should be filtered by abundance and prevalence. Some rare species can screw the permutations 
#and lead to different cluster number
#Filtering only abundant and prevalent species, keep only taxa belonging to genera that have
#Filtering and keeping only those with a prevalence of 5% or more, and with relative abundance over 0.05% in all samples
phylo_species <- phyloseq_big %>%
  tax_filter(min_prevalence = 0.05, min_total_abundance  = 0.0001, tax_level = "Species")
species <- as.data.frame(phylo_species@otu_table)


# DMM app function creation ====
dmm_app<-function(df,only.clust=F, interactive=T){
  require(DirichletMultinomial)
  require(parallel)
  require(xtable)
  require(gtools)
  require(gdata)
  outpath<-invisible(tcltk::tk_choose.dir(default = getwd(), caption = "Select output directory"))
  dmm_otu_df<-as.matrix(df)
  #colnames(dmm_otu_df) <- str_match(colnames(dmm_otu_df), "s__(.*)")[,2]
  fit<-mclapply(1:7, dmn, count=dmm_otu_df, verbose=T)
  lplc <- data.frame(Dir_Comp= seq(1:7),Value=sapply(fit, laplace))
  if(only.clust==F){
    p<-ggplot(lplc, aes(Dir_Comp, Value))+
      xlab("Number of Dirichlet Clusters")+
      ylab("Model Fit")+
      geom_line()+
      geom_point()
    print(p)
    invisible(readline(prompt="Press [enter] to continue\nPress [Esc] to interrupt\n"))
    best<-fit[[which.min(lplc$Value)]]
    p0 <- fitted(fit[[1]], scale=TRUE)
    p4 <- fitted(best, scale=TRUE)
    print(p4)
    colnames(p4) <- paste("cluster", 1:dim(p4)[2], sep="")
    meandiff <- colSums(abs(p4 - as.vector(p0)))
    diff <- rowSums(abs(p4 - as.vector(p0)))
    o.diff <- order(diff, decreasing=TRUE)
    c.diff <- cumsum(diff[o.diff]) / sum(diff)
    df_small <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 30)
    all_baku <- cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff)
    View(df_small)
    invisible(readline(prompt="Press [enter] to continue\nPress [Esc] to interrupt\n"))
    try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
    sav.ans<-invisible(readline(prompt="Save Table and plot? (y/n)"))
    if(startsWith(sav.ans,"y")){
      require(svglite)
      svglite(paste0(outpath,"/DMM_entero_bacteria_species_heatmap.svg"))
      try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
      dev.off()
      write.table(df_small, paste0(outpath,"/S8.4_DMM_top30_driving_bacterial_species.txt"),quote = F, sep = "\t", col.names = T)
      write.table(all_baku, paste0(outpath,"/S8.4_DMM_all_driving_bacterial_species.txt"),quote = F, sep = "\t", col.names = T)
    }
    invisible(readline(prompt="Press [enter] to continue\nPress [Esc] to interrupt\n"))
    return(mixture(best, assign = TRUE))
  }
  else {
    best<-fit[[which.min(lplc$Value)]]
    p0 <- fitted(fit[[1]], scale=TRUE)
    p4 <- fitted(best, scale=TRUE)
    View(p4)
    colnames(p4) <- paste("m", 1:dim(p4)[2], sep="")
    meandiff <- colSums(abs(p4 - as.vector(p0)))
    diff <- rowSums(abs(p4 - as.vector(p0)))
    o.diff <- order(diff, decreasing=TRUE)
    c.diff <- cumsum(diff[o.diff]) / sum(diff)
    if(interactive==T){
      disp<-invisible(readline(prompt="Show Table and Heatmap?(y/n) "))
      if(startsWith(disp,"y",ignore.case = T)){
        df <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 10)
        View(df)
        try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
      }
      else {return(mixture(best, assign = TRUE))}
    }
    else {
      df <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 10)
      View(df)
      try(heatmapdmn(dmm_otu_df, fit[[1]], best, 30),silent = T)
      return(mixture(best, assign = TRUE))
    }
  }
}

# Creating DMM clusters ====
OTU<-t(species)

#IMPORTANT: a window will pop, chose a folder to save data, otherwise the script won't run
#Plots 120-122
set.seed(190)
DMM_Cluster<-dmm_app(OTU,only.clust = F, interactive = T) 

# Change numbering of clusters:
tmp_dir <- c("1" = 5, "2" = 6, "3" = 4, "4" = 1, "5" = 2, "6" = 3)
tmp_DMM_Cluster <- tmp_dir[DMM_Cluster]
names(tmp_DMM_Cluster) <- names(DMM_Cluster)
DMM_Cluster <- tmp_DMM_Cluster

# heatmap ====
### My own code, to create nice DMM heatmap with mean composition of top 30 taxa and divercity indices for it
# data frame
# 1st: run above code to get 'best; object
df_DMM <- rbind.fill(data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1) %>%
                                              tax_transform("compositional", rank = "Species"))@otu_table))), 
                     data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2) %>%
                                              tax_transform("compositional", rank = "Species"))@otu_table))),
                     data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3) %>%
                                              tax_transform("compositional", rank = "Species"))@otu_table))),
                     data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4) %>%
                                              tax_transform("compositional", rank = "Species"))@otu_table))),
                     data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5) %>%
                                              tax_transform("compositional", rank = "Species"))@otu_table))),
                     data.frame(t(rowMeans((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6) %>%
                                              tax_transform("compositional", rank = "Species"))@otu_table))))
df_DMM[is.na(df_DMM)] <- 0

df_small <- read.table("Output/Files/S8.4_DMM_top30_driving_bacterial_species.txt")
cairo_ps("H:/Penders_lab/Theoretical_work/Paper_2/S8.4_DMM_heatmap.eps", height = 10, width = 12)
Heatmap(t(df_DMM[,rownames(df_small)]), column_split = c(1:6), border = TRUE,
        col = colorRamp2(c(0, 0.4), c("#ffffff", "#005824")),
        cluster_rows = FALSE, cluster_columns = FALSE, 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(t(df_DMM[,rownames(df_small)])[i, j], digits = 3), x, y, gp = gpar(fontsize = 8))
        },
        column_names_gp = gpar(cex = 1), column_names_rot = 0,
        heatmap_legend_param = list(title = "Prev"),
        column_labels = paste0("n=",as.vector(table(DMM_Cluster))),
        top_annotation =  HeatmapAnnotation(bar = 1:6, 
                                            col = list(bar = c("1" = "#796248", "2" = "#a47852", "3" = "#cd9c58", "4" = "#e6cd98", "5" = "#f8edcf", "6" = "#d4ede8")),
                                            show_annotation_name = FALSE,
                                            show_legend = FALSE)) 
dev.off()

# Alpha diversity ====
df_dmm_diversity <- rbind(data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 1))@otu_table)), dmm = 1),
                          data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 2))@otu_table)), dmm = 2),
                          data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3))@otu_table), index = "shannon"), 
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 3))@otu_table)), dmm = 3),
                          data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4))@otu_table), index = "shannon"),
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 4))@otu_table)), dmm = 4),
                          data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5))@otu_table), index = "shannon"),
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 5))@otu_table)), dmm = 5),
                          data.frame(shannon = vegan::diversity(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6))@otu_table), index = "shannon"),
                                     specnr = vegan::specnumber(t((phyloseq_big %>% ps_mutate(dmm = DMM_Cluster) %>% ps_filter(dmm == 6))@otu_table)), dmm = 6))

df_dmm_diversity <- df_dmm_diversity %>% mutate(dmm = factor(dmm, levels = c(1:6)))
p1 <- ggplot(data = df_dmm_diversity, aes(x=dmm, y=shannon, col = dmm)) + 
  geom_boxplot() +
  labs(x = "DMM cluster", y = "Shannon diversity") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = c(1:6)) +
  scale_color_manual(values = c("#796248","#a47852","#cd9c58", "#e6cd98","#f8edcf", "#d4ede8")) +
  theme(legend.position = "none")

p2 <-  ggplot(data = df_dmm_diversity, aes(x=dmm, y=specnr, col = dmm)) + 
  geom_boxplot() +
  labs(x = "DMM cluster", y = "Species number") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_discrete(labels = c(1:6)) +
  scale_color_manual(values = c("#796248","#a47852","#cd9c58", "#e6cd98","#f8edcf", "#d4ede8")) +
  theme(legend.position = "none")

p1 + p2 + plot_layout(ncol = 2)

#creating MDSplot ====
set.seed(190)
x<-metaMDS(OTU, distance = "bray")
points.x<-as.data.frame(x$points)
DMM_Cluster <- as.factor(DMM_Cluster)
set.seed(190)
ggplot(points.x, aes(MDS1, MDS2, color=DMM_Cluster))+geom_point()+ggtitle("DMM_Clusters")+
  theme_bw()+scale_color_viridis(discrete = T, option = "viridis")

#save DMM_Cluster information per sample - change output folder accordingly
names(DMM_Cluster) <- phylo_species@sam_data$Sample_aliquote
#write.table(as.data.frame(DMM_Cluster), "S8.4_DMM_cluster_bacteria_species_table.txt", quote = F, row.names = T, col.names = T, sep = "\t")


