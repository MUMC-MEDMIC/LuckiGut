# load libraries ====
library(tidyverse)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# selecting meta ====
meta <- as.data.frame(as.matrix(infant@sam_data))
meta <- meta %>% dplyr::distinct(Family_ID, .keep_all = TRUE) 

meta <- meta %>% select(all_of(c("delivery_type","older_siblings"))) 
meta <- meta %>% mutate(vaginal = ifelse(delivery_type == "vaginal", "yes", "no"))

rm(infant)

# creating 2x2 table ====
#rows = older siblings, cols = vaginal birth
tab <- table(OlderSiblings = meta$older_siblings,
             VaginalBirth  = meta$vaginal)

tab2 <- tab[c("yes","no"), c("yes","no")]

tab2

# Compute OR (and CI) fisher ====
#odds of vaginal birth (yes) in infants with older siblings (yes) compared to infants without older siblings (no)
ft <- fisher.test(tab2, alternative = "two.sided")  # gives OR + exact CI
ft$estimate              # odds ratio
ft$conf.int              # 95% CI
ft$p.value

# chi-square test ====
chisq.test(tab2, correct = FALSE)   # Pearson chi-square
chisq.test(tab2, correct = TRUE)    # with Yates correction (more conservative)

a <- tab2["yes","yes"]
b <- tab2["yes","no"]
c <- tab2["no","yes"]
d <- tab2["no","no"]

OR <- (a*d)/(b*c)

# log(OR) SE and Wald 95% CI
se_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
CI <- exp(log(OR) + c(-1, 1) * 1.96 * se_logOR)

OR
CI                        
  
# if any cell is 0 : a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5                      
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        