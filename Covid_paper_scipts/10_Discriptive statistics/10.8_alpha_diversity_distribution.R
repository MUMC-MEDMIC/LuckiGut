# load libraries ====
library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyverse)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in the data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# taking meta out ====
meta <- as.data.frame(as.matrix(infant@sam_data))

meta$Age_individual <- factor(meta$Age_individual , levels = c("1-2 weeks",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))

rm(infant)

# Shapiro-Wilk normality tests ====
alpha <- c("ENS", "shannon_index","observed_richness")

shapiro_res <- meta %>%
  select(Age_individual, all_of(alpha)) %>%
  pivot_longer(cols = all_of(alpha),
               names_to = "variable",
               values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  group_by(Age_individual, variable) %>%
  summarise(
    statistic = {
      x <- value[!is.na(value)]
      if (length(x) < 3) NA_real_ else shapiro.test(x)$statistic
    },
    p.value = {
      x <- value[!is.na(value)]
      if (length(x) < 3) NA_real_ else shapiro.test(x)$p.value
    },
    .groups = "drop"
  )
# optional: reorder columns
shapiro_res <- shapiro_res[, c("Age_individual", "variable", "statistic", "p.value")]

# Q-Q plots ====
meta %>% ggplot(aes(sample = as.numeric(ENS))) +
  stat_qq() +
  stat_qq_line(col = "blue") +facet_wrap(~Age_individual)+
  labs(title = "Q-Q plot for ENS", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

meta %>% ggplot(aes(sample = as.numeric(observed_richness))) +
  stat_qq() +
  stat_qq_line(col = "blue") +facet_wrap(~Age_individual)+
  labs(title = "Q-Q plot for Richness", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

# histogram ====
meta %>% ggplot(aes(x = as.numeric(ENS))) +
  geom_histogram(binwidth = 0.5) +facet_wrap(~Age_individual)+
  labs(x = "ENS values", y = "Frequency", title = "Histogram of ENS")

meta %>% ggplot(aes(x = as.numeric(observed_richness))) +
  geom_histogram(binwidth = 0.5) +facet_wrap(~Age_individual)+
  labs(x = "Observed richness values", y = "Frequency", title = "Histogram of Observed richness")

# Q-Q plot of residuals (studentized) for all TP ====
#fitted model objects 
m_ens <- lm(ENS ~ Age_individual, data = meta)
m_rich <- lm(observed_richness ~ Age_individual, data = meta)

# THEN plot residuals as qq
qqnorm(rstudent(m_ens), main = "Q-Q plot: studentized residuals for ENS")
qqline(rstudent(m_ens), col = "red")

qqnorm(rstudent(m_rich), main = "Q-Q plot: studentized residuals for richness")
qqline(rstudent(m_rich), col = "red")

#plot residuals vs fitted
plot(fitted(m_ens), rstudent(m_ens),
     xlab = "Fitted values",
     ylab = "Studentized residuals",
     main = "Residuals vs fitted for ENS")
abline(h = 0, col = "red")

plot(fitted(m_rich), rstudent(m_rich),
     xlab = "Fitted values",
     ylab = "Studentized residuals",
     main = "Residuals vs fitted for ENS")
abline(h = 0, col = "red")

# Q-Q plot of residuals (studentized) per TP ====
#Fit Intercept-Only Models Per Timepoint
models_ens <- meta %>%
  group_by(Age_individual) %>%
  do(model = lm(ENS ~ 1, data = .)) %>%  # Intercept-only model per group
  deframe()  # Convert to named list

#Extract Residuals Per Timepoint
ens_resids <- map_dfr(names(models_ens), ~ {
  model <- models_ens[[.x]]  
  data.frame(
    Age_individual = .x,      
    fitted         = fitted(model),
    resid          = residuals(model),
    stud_resid     = rstudent(model)
  )
})

ens_resids$Age_individual <- factor(ens_resids$Age_individual , levels = c("1-2 weeks",
                                                               "4 weeks",
                                                               "8 weeks",
                                                               "4 months",
                                                               "5 months",
                                                               "6 months",
                                                               "9 months",
                                                               "11 months",
                                                               "14 months"))

#Q-Q Plots Per Timepoint
ens_resids %>%
  ggplot(aes(sample = stud_resid)) +
  stat_qq(color = "black") +
  stat_qq_line(color = "red", linewidth = 1) +
  facet_wrap(~Age_individual, scales = "free") +
  labs(
    title = "ENS Studentized Residuals: Q-Q Plots by Timepoint",
    x = "Theoretical Quantiles", 
    y = "Studentized Residuals"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Q-Q Plots Per Timepoint
ens_resids %>%
  ggplot(aes(x = fitted, y = stud_resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, col = "red", lwd = 1) +
  facet_wrap(~Age_individual, scales = "free") +
  labs(title = "ENS Studentized Residuals vs Fitted by Timepoint",
       x = "Fitted Values", y = "Studentized Residuals") +
  theme_minimal()

ens_resids %>%
  ggplot(aes(x = fitted, y = stud_resid)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, col = "red", lwd = 1) +
  facet_wrap(~Age_individual, scales = "free") +
  labs(title = "ENS Studentized Residuals vs Fitted by Timepoint") +
  theme_minimal()











