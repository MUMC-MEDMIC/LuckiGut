# load libraries ====
library(microViz)
library(ggplot2)
library(tidyverse)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# load in data ====
model = read.delim("Output/Files/S6.17_LR_pandemic_ENS_scaled_residuals_raw_without_outliers.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# prepare data for the plot =====
model <- model %>% filter(Age_individual == "4 months" | Age_individual == "5 months")
model$Age_individual <- factor(model$Age_individual , levels = c("4 months",
                                                                           "5 months"))
# QQ plots ====
model %>%
  ggplot(aes(sample = stud_resid)) +
  stat_qq(color = "black") +
  stat_qq_line(color = "red", linewidth = 1) +
  facet_wrap(~Age_individual, scales = "free") +
  labs(
    title = "Pandemic: Q-Q Plots. ENS studentized residuals",
    x = "Theoretical Quantiles", 
    y = "Studentized Residuals"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Residuals vs fitted ====
model %>%
  ggplot(aes(x = fitted, y = stud_resid)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, col = "red", lwd = 1) +
  facet_wrap(~Age_individual, scales = "free") +
  labs(title = "Pandemic: Studentized residuals vs fitted. ENS") +
  theme_minimal()

# selecting outliers ====
outlier4 <- model %>% filter(Age_individual == "4 months")
i_max <- which.max(outlier4$stud_resid)
outlier4_max_rowname <- rownames(outlier4)[i_max]

outlier5 <- model %>% filter(Age_individual == "5 months")
i_max <- which.max(outlier5$stud_resid)
outlier5_max_rowname <- rownames(outlier5)[i_max]

rm(outlier4,outlier5,i_max)
