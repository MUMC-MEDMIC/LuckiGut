# libraries ====
library(ggpubr) #for manual p values
library(tidyverse)
library(MASS)
library(matrixcalc)
library(viridis)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import and clean data ====
distance <- read.delim("Output/Files/S9.1_aitchinson_distance_between_TP_beta_diversity.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning ====
distance <- distance %>% mutate(pairs_new = recode(pairs, 
                                                   "1-2 weeks vs 4 weeks" = "1_4w", 
                                                   "4 weeks vs 8 weeks" = "4_8w",
                                                   "8 weeks vs 4 months" = "8_4m",
                                                   "4 months vs 5 months" = "4_5m",
                                                   "5 months vs 6 months" = "5_6m",
                                                   "6 months vs 9 months" = "6_9m",
                                                   "9 months vs 11 months" = "9_11m",
                                                   "11 months vs 14 months" = "11_14m"))

mack <- distance %>% dplyr::select(all_of(c("pairs_new","Family_ID1","value")))

distance$pairs<- factor(distance$pairs, 
                        levels = c("1-2 weeks vs 4 weeks",
                                   "4 weeks vs 8 weeks",
                                   "8 weeks vs 4 months",
                                   "4 months vs 5 months",
                                   "5 months vs 6 months",
                                   "6 months vs 9 months",
                                   "9 months vs 11 months",
                                   "11 months vs 14 months"))

# Ensure factors are properly set
str(mack)
mack <- mack %>%
  mutate(
    Family_ID1 = as.factor(Family_ID1),
    pairs_new = as.factor(pairs_new))

# creating matrix ====
# Reshape to matrix: rows = treatments, columns = families
mack_matrix <- mack %>%
  pivot_wider(
    id_cols = pairs_new,
    names_from = Family_ID1,
    values_from = value
  ) %>%
  column_to_rownames("pairs_new") %>%
  as.matrix()

# removing block (ak Family with only one observation)
sparse_blocks_names <- colnames(mack_matrix)[which(colSums(!is.na(mack_matrix)) == 1)]
mack_matrix_clean <- mack_matrix[, !(colnames(mack_matrix) %in% sparse_blocks_names)]

rm(mack,mack_matrix,sparse_blocks_names)

# skilling mack test formula ====
Ski.Mack.2<-function (y, groups = NULL, blocks = NULL, simulate.p.value = FALSE, 
                      B = 10000) 
{
  require(MASS)
  require(matrixcalc)
  specify_decimal <- function(x, k) {formatC(x,format = "e", digits = k)}
  options(digits = 10)
  if (!is.matrix(y)) {
    blocks <- factor(blocks)
    k <- nlevels(blocks)
    groups <- factor(groups)
    t <- nlevels(groups)
    d <- matrix(unlist(split(y, blocks)), ncol = k, byrow = FALSE)
    if (any(is.na(groups)) || any(is.na(blocks))) {
      stop("NA's are not allowed in groups or blocks")
    }
    if (any(diff(c(length(y), length(groups), length(blocks))))) {
      stop("y, groups and blocks must have the same length")
    }
    for (i in 1:k) {
      if (sum(is.na(d[, i])) >= (t - 1)) {
        stop("Block#", i, " has only one observation. Please remove this block")
      }
    }
  }
  if (is.matrix(y)) {
    d <- as.matrix(y)
    k <- dim(d)[2]
    t <- dim(d)[1]
  }
  for (i in 1:k) {
    if (sum(is.na(d[, i])) >= (t - 1)) {
      stop("Block#", i, " has only one observation. Please remove this block")
    }
  }
  y.rank.NA <- matrix(NA, nrow = t, ncol = k)
  y.rank <- matrix(NA, nrow = t, ncol = k)
  for (i in 1:k) {
    y.rank.NA[, i] <- rank(d[, i], ties.method = c("average"), 
                           na.last = "keep")
    y.rank[, i] <- ifelse(is.na(y.rank.NA[, i]) == TRUE, 
                          (sum(!is.na(y.rank.NA[, i])) + 1)/2, y.rank.NA[, 
                                                                         i])
  }
  A = matrix(NA, nrow = 1, ncol = t)
  for (j in 1:t) {
    Ai = matrix(NA, nrow = 1, ncol = k)
    for (i in 1:k) {
      Ai[1, i] <- sqrt(12/(sum(!is.na(y.rank.NA[, i])) + 
                             1)) * (y.rank[j, i] - (sum(!is.na(y.rank.NA[, 
                                                                         i])) + 1)/2)
    }
    A[1, j] <- sum(Ai)
  }
  transform.matrix <- function(mat) {
    num.of.rows <- nrow(mat)
    num.of.cols <- ncol(mat)
    output <- matrix(rep(0, num.of.rows * num.of.rows), c(num.of.rows, 
                                                          num.of.rows))
    for (i in seq(1, num.of.rows - 1)) for (j in seq(i + 
                                                     1, num.of.rows)) {
      output[i, j] <- -sum(mat[i, ] * mat[j, ])
      output[j, i] <- output[i, j]
    }
    for (i in 1:num.of.rows) {
      output[i, i] <- (-1) * sum(output[, i])
    }
    output
  }
  CovMat <- transform.matrix(!is.na(y.rank.NA))
  Cov.Inv <- ginv(CovMat)
  rank.cov <- matrix.rank(Cov.Inv)
  T <- A %*% Cov.Inv %*% t(A)
  SM.test <- function(y) {
    if (is.matrix(y)) {
      d <- as.matrix(y)
      k <- dim(d)[2]
      t <- dim(d)[1]
    }
    y.rank.NA <- matrix(NA, nrow = t, ncol = k)
    y.rank <- matrix(NA, nrow = t, ncol = k)
    for (i in 1:k) {
      y.rank.NA[, i] <- rank(d[, i], ties.method = c("average"), 
                             na.last = "keep")
      y.rank[, i] <- ifelse(is.na(y.rank.NA[, i]) == TRUE, 
                            (sum(!is.na(y.rank.NA[, i])) + 1)/2, y.rank.NA[, 
                                                                           i])
    }
    A = matrix(NA, nrow = 1, ncol = t)
    for (j in 1:t) {
      Ai = matrix(NA, nrow = 1, ncol = k)
      for (i in 1:k) {
        Ai[1, i] <- sqrt(12/(sum(!is.na(y.rank.NA[, i])) + 
                               1)) * (y.rank[j, i] - (sum(!is.na(y.rank.NA[, 
                                                                           i])) + 1)/2)
      }
      A[1, j] <- sum(Ai)
    }
    transform.matrix <- function(mat) {
      num.of.rows <- nrow(mat)
      num.of.cols <- ncol(mat)
      output <- matrix(rep(0, num.of.rows * num.of.rows), 
                       c(num.of.rows, num.of.rows))
      for (i in seq(1, num.of.rows - 1)) for (j in seq(i + 
                                                       1, num.of.rows)) {
        output[i, j] <- -sum(mat[i, ] * mat[j, ])
        output[j, i] <- output[i, j]
      }
      for (i in 1:num.of.rows) {
        output[i, i] <- (-1) * sum(output[, i])
      }
      output
    }
    CovMat <- transform.matrix(!is.na(y.rank.NA))
    Cov.Inv <- ginv(CovMat)
    T.SM <- A %*% Cov.Inv %*% t(A)
  }
  if (simulate.p.value == TRUE) {
    logic.d <- is.na(d)
    num.missing <- matrix(NA, nrow = k, ncol = 1)
    for (i in 1:k) {
      num.missing[i, 1] <- sum(logic.d[, i])
    }
    remove <- function(x, logic) {
      x <- x[logic == FALSE]
      return(x)
    }
    simulated.SM <- matrix(NA, nrow = B, ncol = 1)
    for (b.num in 1:B) {
      y.sim <- matrix(NA, nrow = t, ncol = k)
      for (i in 1:k) {
        if (num.missing[i, 1] == 0) {
          count <- 0
          rand.int <- sample(y.rank[, i], replace = FALSE)
          for (j in 1:t) {
            count = count + 1
            y.sim[j, i] <- rand.int[count]
          }
        }
        if (num.missing[i, 1] > 0) {
          count <- 0
          rand.int <- sample(remove(y.rank[, i], logic.d[, 
                                                         i]), replace = FALSE)
          for (j in 1:t) {
            if (logic.d[j, i] == TRUE) {
              y.sim[j, i] <- NA
            }
            if (logic.d[j, i] == FALSE) {
              count = count + 1
              y.sim[j, i] <- rand.int[count]
            }
          }
        }
      }
      simulated.SM[b.num] <- SM.test(y.sim)
    }
    sim.pval <- round(sum(simulated.SM >= as.numeric(T))/B, 
                      digits = 7)
  }
  if (simulate.p.value == TRUE) {
    cat("\n")
    cat(c("Skillings-Mack Statistic = "), specify_decimal(T, 
                                                          6), c(", p-value = "), specify_decimal(pchisq(T, 
                                                                                                        df = rank.cov, lower.tail = FALSE), 6), "\n")
    cat(c("Note: the p-value is based on the chi-squared distribution with d.f. = "), 
        rank.cov, "\n")
    cat(c("Based on B = "), specify_decimal(B, 0), c(", Simulated p-value = "), 
        specify_decimal(sim.pval, 6), "\n")
    cat("\n")
  }
  if (simulate.p.value == FALSE) {
    cat("\n")
    cat(c("Skillings-Mack Statistic = "), specify_decimal(T, 
                                                          2), c(", p-value = "), specify_decimal(pchisq(T, 
                                                                                                        df = rank.cov, lower.tail = FALSE), 2), "\n")
    cat(c("Note: the p-value is based on the chi-squared distribution with d.f. = "), 
        rank.cov, "\n")
    cat("\n")
  }
  return(list( Statistics = as.numeric(specify_decimal(T,2)),
               p.value=as.numeric(specify_decimal(pchisq(T,df = rank.cov, lower.tail = FALSE), 2)),
               Nblocks = k, Ntreatments = t, rawdata = d, rankdata = y.rank, 
               varCovarMatrix = CovMat, adjustedSum = A))
}

# running the test ====
Ski.Mack.2(y = mack_matrix_clean, simulate.p.value = TRUE)

# wilcox post hoc ==== 
plyr::count(mack$pairs_new)
wil_1 <- distance %>% filter(pairs_new == "1_4w" | pairs_new == "4_8w") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))
wil_2 <- distance %>% filter(pairs_new == "4_8w" | pairs_new == "8_4m") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))
wil_3 <- distance %>% filter(pairs_new == "8_4m" | pairs_new == "4_5m") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))
wil_4 <- distance %>% filter(pairs_new == "4_5m" | pairs_new == "5_6m") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))
wil_5 <- distance %>% filter(pairs_new == "5_6m" | pairs_new == "6_9m") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))
wil_6 <- distance %>% filter(pairs_new == "6_9m" | pairs_new == "9_11m") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))
wil_7 <- distance %>% filter(pairs_new == "9_11m" | pairs_new == "11_14m") %>% subset(duplicated(Family_ID1) | duplicated(Family_ID1, fromLast=TRUE))

#paired: a logical value specifying that we want to compute a paired Wilcoxon test
test_1 <- wilcox.test(value ~ pairs_new, data = wil_1, alternative="two.sided",paired=T)
stat.test <- data.frame(
  p_value = test_1$p.value,
  Z = qnorm(test_1$p.value / 2),
  group1 = '1-2 weeks vs 4 weeks',
  group2 = '4 weeks vs 8 weeks',
  W_value = test_1$statistic,
  alternative = test_1$alternative,
  method = test_1$method
)

test_2 <- wilcox.test(value ~ pairs_new, data = wil_2, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_2$p.value,
  Z = qnorm(test_2$p.value / 2),
  group1 = '4 weeks vs 8 weeks',
  group2 = '8 weeks vs 4 months',
  W_value = test_2$statistic,
  alternative = test_2$alternative,
  method = test_2$method
)
stat.test <- rbind(stat.test, s) 

test_3 <- wilcox.test(value ~ pairs_new, data = wil_3, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_3$p.value,
  Z = qnorm(test_3$p.value / 2),
  group1 = '8 weeks vs 4 months',
  group2 = '4 months vs 5 months',
  W_value = test_3$statistic,
  alternative = test_3$alternative,
  method = test_3$method
)
stat.test <- rbind(stat.test, s) 

test_4 <- wilcox.test(value ~ pairs_new, data = wil_4, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_4$p.value,
  Z = qnorm(test_4$p.value / 2),
  group1 = '4 months vs 5 months',
  group2 = '5 months vs 6 months',
  W_value = test_4$statistic,
  alternative = test_4$alternative,
  method = test_4$method
)
stat.test <- rbind(stat.test, s) 

test_5 <- wilcox.test(value ~ pairs_new, data = wil_5, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_5$p.value,
  Z = qnorm(test_5$p.value / 2),
  group1 = '5 months vs 6 months',
  group2 = '6 months vs 9 months',
  W_value = test_5$statistic,
  alternative = test_5$alternative,
  method = test_5$method
)
stat.test <- rbind(stat.test, s) 

test_6 <- wilcox.test(value ~ pairs_new, data = wil_6, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_6$p.value,
  Z = qnorm(test_6$p.value / 2),
  group1 = '6 months vs 9 months',
  group2 = '9 months vs 11 months',
  W_value = test_6$statistic,
  alternative = test_6$alternative,
  method = test_6$method
)
stat.test <- rbind(stat.test, s) 

test_7 <- wilcox.test(value ~ pairs_new, data = wil_7, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_7$p.value,
  Z = qnorm(test_7$p.value / 2),
  group1 = '9 months vs 11 months',
  group2 = '11 months vs 14 months',
  W_value = test_7$statistic,
  alternative = test_7$alternative,
  method = test_7$method
)
stat.test <- rbind(stat.test, s) 

rm(test_1,test_2,test_3,test_4,test_5,test_6,test_7, s)
rm(wil_1,wil_2,wil_3,wil_4,wil_5,wil_6,wil_7)

#renaming columns in p-value df
wilcox = stat.test
names(wilcox)[names(wilcox) == "Freq"] <- "p_value" 
wilcox$p_value_test_adj <- p.adjust(wilcox$p_value, method = "BH")
rm(stat.test)

#add significance
wilcox <- wilcox %>% 
  mutate(sig_p_adj = if_else(
    condition = .$p_value_test_adj > 0.05, 
    true = NA, 
    false = if_else(
      condition = .$p_value_test_adj < 0.001,
      true = "***",
      false = if_else(
        condition = .$p_value_test_adj < 0.05&.$p_value_test_adj > 0.01,
        true = "*",
        false = "**"))))

#write.table (wilcox, "S23_wilcox_beta_diversity.txt", row.names=FALSE,sep = "\t")

# shortening the names for the figure ====
distance$pairs <- str_replace_all(distance$pairs, c(
  "months" = "mos",
  "weeks"  = "wks"
))

wilcox$group1 <- str_replace_all(wilcox$group1, c(
  "months" = "mos",
  "weeks"  = "wks"
))

wilcox$group2 <- str_replace_all(wilcox$group2, c(
  "months" = "mos",
  "weeks"  = "wks"
))

distance$pairs<- factor(distance$pairs, 
                        levels = c("1-2 wks vs 4 wks",
                                   "4 wks vs 8 wks",
                                   "8 wks vs 4 mos",
                                   "4 mos vs 5 mos",
                                   "5 mos vs 6 mos",
                                   "6 mos vs 9 mos",
                                   "9 mos vs 11 mos",
                                   "11 mos vs 14 mos"))

# Create the boxplot with p values ====
boxplot <- ggplot(distance, aes(x = pairs, y = value, color = pairs)) +
  geom_boxplot() +theme_bw()+
  labs(x = "", y = "Aitchison distance", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(0, 90) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 1), 
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 1),
        legend.title = element_text(size = 1))+
  scale_color_viridis(discrete = TRUE)

#adding p-values
boxplot <- boxplot + stat_pvalue_manual(wilcox,
                                        y.position = 80, step.increase = 0.0000001,label = "sig_p_adj",hide.ns = TRUE,
                                        size = 5, tip.length = 0)
boxplot

# saving figures for the paper ====
#ggsave("S9.2_aitchinson_distance_skilling_wilcox_boxplot.pdf", boxplot, device = "pdf",width = 400, height = 300, dpi = 300, units = "mm")

ggsave(
  filename   = "S9.2_aitchinson_distance_skilling_wilcox_boxplot.tiff",
  plot       = boxplot,
  width      = 84,      # full page width
  height     = 60,       # max = 225 mm
  units      = "mm",
  dpi        = 300,
  compression = "lzw",
  device     = "tiff"
)



