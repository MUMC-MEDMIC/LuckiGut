# load libraries ====
library(microViz)
library(vegan)
library(ggplot2)
library(gdata)
library(gtools)
library(viridis)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# DMM species best model ====
#Dirichlet Multinomial Mixture (DMM) clustering, an unsupervised clustering method that uses Laplace approximation to 
#identify groups of communities (enterotypes) with similar composition
#Species table should be filtered by abundance and prevalence. Some rare species can screw the permutations 
#and lead to different cluster number

#Filtering only abundant and prevalent species, keep only taxa belonging to genera that have done in previous scripts
#Filtering and keeping only those with a prevalence of 5% or more, and with relative abundance over 0.05% in all samples done in previous scripts
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
species <- as.data.frame(infant@otu_table)

#DMM app fucntion creation, 
dmm_app<-function(df,only.clust=F, interactive=T){
  require(DirichletMultinomial)
  require(parallel)
  require(xtable)
  require(gtools)
  require(gdata)
  outpath<-invisible(tcltk::tk_choose.dir(default = getwd(), caption = "Select output directory"))
  dmm_otu_df<-as.matrix(df)
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
    colnames(p4) <- paste("cluster", 1:dim(p4)[2], sep="") #was m before not cluster
    meandiff <- colSums(abs(p4 - as.vector(p0)))
    diff <- rowSums(abs(p4 - as.vector(p0)))
    o.diff <- order(diff, decreasing=TRUE)
    c.diff <- cumsum(diff[o.diff]) / sum(diff)
    df_small <- head(cbind(Mean=p0[o.diff], p4[o.diff,], diff=diff[o.diff], c.diff), 10)
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
      write.table(df_small, paste0(outpath,"/S4.1_DMM_top10_driving_bacterial_species.txt"),quote = F, sep = "\t", col.names = T)
      write.table(all_baku, paste0(outpath,"/S4.1_DMM_all_driving_bacterial_species.txt"),quote = F, sep = "\t", col.names = T)
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

#Creating DMM clusters
OTU<-t(species)

#IMPORTANT: a window will pop, chose a folder to save data, otherwise the script won't run
#save each plot as pdf manually also!
set.seed(400)
DMM_Cluster<-dmm_app(OTU,only.clust = F, interactive = T) 

# saving data ====
#write.table(as.data.frame(DMM_Cluster), "S4.1_DMM_cluster_species_table.txt", quote = F, row.names = T, col.names = T, sep = "\t")









