wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/")
setwd(wd)

library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stats)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(viridis)
library(car)
library(multcomp)

#Read in data
fit <- read_csv("./fitness_data/Fitness_data_simplified.csv")
samples <- read_csv("../02_degs/assembled_transcripts_final/sample.info.all.modified.csv")

#Clean up the clone names (some spelling errors)
fit$Clone[fit$Clone == 'SWP4'] <- 'CBSWP4'
fit$Clone[fit$Clone == 'NWP7'] <- 'CBNWP7'
fit$Clone[fit$Clone == 'ARLCO10'] <- 'ARLC10'
fit$Clone[fit$Clone == 'LPI'] <- 'LP1'

#Remove CBSWP4
fit <- fit[fit$Clone != 'CBSWP4', ]

#Convert Clone to Genotype Clone
fit$Clone[fit$Clone == 'ARLC10'] <- 'C1'
fit$Clone[fit$Clone == 'BING2'] <- 'C2'
fit$Clone[fit$Clone == 'CBNWP7'] <- 'C3'
fit$Clone[fit$Clone == 'CF1'] <- 'C4'
fit$Clone[fit$Clone == 'LF8'] <- 'C4'
fit$Clone[fit$Clone == 'LP1'] <- 'C5'
fit$Clone[fit$Clone == 'PWA10'] <- 'C5'
fit$Clone[fit$Clone == 'PWA7'] <- 'C6'

fit$Temperature <- as.factor(fit$Temperature)


#MC Method
##########################################

#Calculated the total cumulative offspring for each individual
cumulative_offspring <- fit %>%
  group_by(Unique_ID, Clone, Temperature) %>%
  summarise(Total_Babies = sum(Babies, na.rm = TRUE)) %>%
  ungroup()


# Calculate mean reproduction for each group (Clone and Temperature)
fit_residuals <- cumulative_offspring %>%
  group_by(Clone, Temperature) %>%
  mutate(mean_reproduction = mean(Total_Babies, na.rm = TRUE)) %>%
  ungroup()

# Ccalculate the residual for each individual
fit_residuals <- fit_residuals %>%
  mutate(residual_reproduction = Total_Babies - mean_reproduction)

#Test for normality
shapiro_test <- shapiro.test(fit_residuals$residual_reproduction)
print(shapiro_test)

leveneTest(residual_reproduction ~ Clone * Temperature, data = fit_residuals)

bartlett.test(residual_reproduction ~ Temperature, data = fit_residuals)

qq_plot <- ggplot(fit_residuals, aes(sample = residual_reproduction)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Q-Q Plot of Residuals")

print(qq_plot)

#Perform ANOVA

anova_results <- aov(residual_reproduction ~ Clone * Temperature, data = fit_residuals)
summary(anova_results)

#Now, performing an anova on residuals seems weird because the means are essentially 0.
#I think this is why the results of the ANOVA are 1


#Method of Hartman et al.
##########################################
###For all data
shapiro_test_all <- shapiro.test(cumulative_offspring$Total_Babies)
print(shapiro_test_all)
#not normal

bartlett_all <- bartlett.test(Total_Babies ~ Temperature, data = cumulative_offspring)
print(bartlett_all)
#variances are barely non-sig

anova_all <- aov(Total_Babies ~ Clone * Temperature, data = cumulative_offspring)
summary(anova_all)

kruskal_all <- kruskal.test(Total_Babies ~ Clone , data = cumulative_offspring)
summary(kruskal_all)



###For each clone
c1_data <- cumulative_offspring %>%
  filter(Clone == 'C1')

shapiro_test_c1 <- shapiro.test(c1_data$Total_Babies)
print(shapiro_test_c1)
#normal

bartlett_c1 <- bartlett.test(Total_Babies ~ Temperature, data = c1_data)
print(bartlett_c1)
#one group with high variance

#anova_c1 <- aov(Total_Babies ~ Temperature, data = c1_data)
#summary(anova_c1)

#dunnett_test_c1 <- glht(anova_c1, 
#                     linfct = mcp(Temperature = "Dunnett")) 
#summary(dunnett_test_c1)

wilcox.test(Total_Babies ~ Temperature, data = c1_data, alternative = "two.sided", exact = FALSE)

###
c2_data <- cumulative_offspring %>%
  filter(Clone == 'C2')

shapiro_test_c2 <- shapiro.test(c2_data$Total_Babies)
print(shapiro_test_c2)
#normal

bartlett_c2 <- bartlett.test(Total_Babies ~ Temperature, data = c2_data)
print(bartlett_c2)
#variances likely equal

anova_c2 <- aov(Total_Babies ~ Temperature, data = c2_data)
summary(anova_c2)

dunnett_test_c2 <- glht(anova_c2, 
                        linfct = mcp(Temperature = "Dunnett")) 
summary(dunnett_test_c2)

###
c3_data <- cumulative_offspring %>%
  filter(Clone == 'C3')

shapiro_test_c3 <- shapiro.test(c3_data$Total_Babies)
print(shapiro_test_c3)
#normal

bartlett_c3 <- bartlett.test(Total_Babies ~ Temperature, data = c3_data)
print(bartlett_c3)
#variances likely equal

anova_c3 <- aov(Total_Babies ~ Temperature, data = c3_data)
summary(anova_c3)

dunnett_test_c3 <- glht(anova_c3, 
                        linfct = mcp(Temperature = "Dunnett")) 
summary(dunnett_test_c3)

###
c4_data <- cumulative_offspring %>%
  filter(Clone == 'C4')

shapiro_test_c4 <- shapiro.test(c4_data$Total_Babies)
print(shapiro_test_c4)
#normal

bartlett_c4 <- bartlett.test(Total_Babies ~ Temperature, data = c4_data)
print(bartlett_c4)
#variances likely equal


anova_c4 <- aov(Total_Babies ~ Temperature, data = c4_data)
summary(anova_c4)

dunnett_test_c4 <- glht(anova_c4, 
                        linfct = mcp(Temperature = "Dunnett")) 
summary(dunnett_test_c4)

###
c5_data <- cumulative_offspring %>%
  filter(Clone == 'C5')

shapiro_test_c5 <- shapiro.test(c5_data$Total_Babies)
print(shapiro_test_c5)
#normal

bartlett_c5 <- bartlett.test(Total_Babies ~ Temperature, data = c5_data)
print(bartlett_c5)
#variances likely equal


anova_c5 <- aov(Total_Babies ~ Temperature, data = c5_data)
summary(anova_c5)

dunnett_test_c5 <- glht(anova_c5, 
                        linfct = mcp(Temperature = "Dunnett")) 
summary(dunnett_test_c5)

###
c6_data <- cumulative_offspring %>%
  filter(Clone == 'C6')

shapiro_test_c6 <- shapiro.test(c6_data$Total_Babies)
print(shapiro_test_c6)
#normal

bartlett_c6 <- bartlett.test(Total_Babies ~ Temperature, data = c6_data)
print(bartlett_c6)
#variances likely equal

anova_c6 <- aov(Total_Babies ~ Temperature, data = c6_data)
summary(anova_c6)

dunnett_test_c6 <- glht(anova_c6, 
                        linfct = mcp(Temperature = "Dunnett")) 
summary(dunnett_test_c6)

###########################################

c1_data <- fit_residuals %>%
  filter(Clone == 'C1')

mean(c1_data$residual_reproduction)

c1_data_20 <- c1_data %>%
  filter(Temperature == 20)

c1_data_25 <- c1_data %>%
  filter(Temperature == 25)

mean(c1_data_20$residual_reproduction)
mean(c1_data_25$residual_reproduction)
median(c1_data_20$residual_reproduction)
median(c1_data_25$residual_reproduction)


anova_c1 <- aov(residual_reproduction ~ Temperature, data = c1_data)
summary(anova_c1)
t_test_results <- t.test(residual_reproduction ~ Temperature, data = c1_data)
print(t_test_results)

###
c2_data <- fit_residuals %>%
  filter(Clone == 'C2')

mean(c2_data$residual_reproduction)

c2_data_20 <- c2_data %>%
  filter(Temperature == 20)

c2_data_25 <- c2_data %>%
  filter(Temperature == 25)

mean(c2_data_20$residual_reproduction)
mean(c2_data_25$residual_reproduction)
median(c2_data_20$residual_reproduction)
median(c2_data_25$residual_reproduction)


anova_c2 <- aov(residual_reproduction ~ Temperature, data = c2_data)
summary(anova_c2)
t_test_results_c2 <- t.test(residual_reproduction ~ Temperature, data = c2_data)
print(t_test_results_c2)

###
c3_data <- fit_residuals %>%
  filter(Clone == 'C3')

mean(c3_data$residual_reproduction)

c3_data_20 <- c3_data %>%
  filter(Temperature == 20)

c3_data_25 <- c3_data %>%
  filter(Temperature == 25)

mean(c3_data_20$residual_reproduction)
mean(c3_data_25$residual_reproduction)
median(c3_data_20$residual_reproduction)
median(c3_data_25$residual_reproduction)


anova_c3 <- aov(residual_reproduction ~ Temperature, data = c3_data)
summary(anova_c3)
t_test_results_c3 <- t.test(residual_reproduction ~ Temperature, data = c3_data)
print(t_test_results_c3)


c1 <- cumulative_offspring %>%
  filter(Clone == "C1")

anova_results_correct <- aov(Total_Babies ~ Temperature, data = c1)
summary(anova_results_correct)
residuals <- resid(anova_results)
plot(fit_residuals$residual_reproduction, residuals, xlab="Measurement", ylab="Residuals")

