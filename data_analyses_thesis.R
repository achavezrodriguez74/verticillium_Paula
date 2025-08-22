# Call packages ----

library(olsrr)
library(multcomp)
library(multcompView)
library(ggpubr)
library(ARTool)
library(emmeans)
library(rcompanion)
library(scales)
library(ggsignif)
library(rstatix)  # https://github.com/kassambara/rstatix
library(tibble)
library(dplyr)
# Conover-Iman
library(DescTools)
library(tidyverse)
library(readxl)
library(bayesplot)
library(brms)

# PCR graphs ----

# Call data
PCR_data = read_excel("datasets/All data combined.xlsx", sheet = "Thesis")
PCR_raw  = PCR_data

# Transform to Log2
PCR_data = PCR_data %>%  mutate(across(c(Expression), function(x) log2(x)))

# Delete WT
PCR_data = PCR_data %>% filter(Treatment != "WT")

# Plot 

# SXM ----
PCR_data_SXM   = PCR_data %>% filter(Stage == "SXM")
PCR_data_SXM.1 = PCR_data_SXM %>% group_by(Gen,Treatment,Order) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_PCR_SXM = ggplot(PCR_data_SXM.1, aes(x=Gen, y=avg, fill=as.factor(Order))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#d0d1e6", "#045a8d", "#fee391", "#e34a33",
                                          "#74c476"), labels = c("ΔVel1", "oeVel1", "ΔVel2", 
                                                                 "oeVel2", "ΔVel2_IDD")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_PCR_SXM

pdf("Figures/Figure_PCR_SXM.pdf",
    width=12,height=12*3/5)
print(Figure_PCR_SXM)
dev.off()

# PDM ----
PCR_data_PDM   = PCR_data %>% filter(Stage == "PDM")
PCR_data_PDM.1 = PCR_data_PDM %>% group_by(Gen,Treatment,Order) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_PCR_PDM = ggplot(PCR_data_PDM.1, aes(x=Gen, y=avg, fill=as.factor(Order))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#d0d1e6", "#045a8d", "#fee391", "#e34a33",
                                          "#74c476"), labels = c("ΔVel1", "oeVel1", "ΔVel2", 
                                                                 "oeVel2", "ΔVel2_IDD")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_PCR_PDM

pdf("Figures/Figure_PCR_PDM.pdf",
    width=12,height=12*3/5)
print(Figure_PCR_PDM)
dev.off()

# PLATE ----
PCR_data_PLATE   = PCR_data %>% filter(Stage == "PLATE")
PCR_data_PLATE.1 = PCR_data_PLATE %>% group_by(Gen,Treatment,Order) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_PCR_PLATE = ggplot(PCR_data_PLATE.1, aes(x=Gen, y=avg, fill=as.factor(Order))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#d0d1e6", "#045a8d", "#fee391", "#e34a33",
                                          "#74c476"), labels = c("ΔVel1", "oeVel1", "ΔVel2", 
                                                                 "oeVel2", "ΔVel2_IDD")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_PCR_PLATE

pdf("Figures/Figure_PCR_PLATE.pdf",
    width=12,height=12*3/5)
print(Figure_PCR_PLATE)
dev.off()

# Bayesian Methods for Group Comparison ----
library(brms)
library(bayestestR)

# SXM ----
PCR_raw_SXM   = PCR_data %>% filter(Stage == "SXM")
PCR_raw_SXM$Treatment = as.factor(PCR_raw_SXM$Treatment)

# LLM1 - ΔVel1
model.LLM1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_SXM[1:6,])

summary(model.LLM1.ΔVel1)
posterior_summary(model.LLM1.ΔVel1)
describe_posterior(model.LLM1.ΔVel1)
rope(model.LLM1.ΔVel1)
posterior = as_draws_df(model.LLM1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.00125

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 7.64e-03 | [-0.49, 0.51] | 51.45% | [-0.20, 0.20] |    75.47% | 1.003 | 1664.00
TreatmentΔVel1 |     3.71 | [ 2.97, 4.45] | 99.88% | [-0.20, 0.20] |        0% | 1.003 | 1584.00

# LLM1 - oeVel1
model.LLM1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(1:3,7:9),])

summary(model.LLM1.oeVel1)
posterior_summary(model.LLM1.oeVel1)
describe_posterior(model.LLM1.oeVel1)
rope(model.LLM1.oeVel1)
posterior.1 = as_draws_df(model.LLM1.oeVel1)
1- mean(posterior.1$b_TreatmentWT > 0) # 0.99275

# Summary of Posterior Distribution 

Parameter   | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
  (Intercept) |   0.22 | [ 0.12,  0.32] | 99.65% | [-0.01, 0.01] |        0% | 1.001 | 1347.00
TreatmentWT |  -0.22 | [-0.36, -0.07] | 99.28% | [-0.01, 0.01] |        0% | 1.000 | 1612.00

# LLM1 - ΔVel2
model.LLM1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(1:3,10:12),])

summary(model.LLM1.ΔVel2)
posterior_summary(model.LLM1.ΔVel2)
describe_posterior(model.LLM1.ΔVel2)
rope(model.LLM1.ΔVel2)
posterior.2 = as_draws_df(model.LLM1.ΔVel2)
1- mean(posterior.2$b_TreatmentΔVel2 > 0) # 0.01175

# Summary of Posterior Distribution 

Parameter      |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | -3.45e-03 | [-0.90, 0.81] | 50.42% | [-0.10, 0.10] |    24.08% | 1.002 | 1736.00
TreatmentΔVel2 |      1.68 | [ 0.39, 2.95] | 98.83% | [-0.10, 0.10] |        0% | 1.002 | 1741.00

# LLM1 - oeVel2
model.LLM1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(1:3,13:15),])

summary(model.LLM1.oeVel2)
posterior_summary(model.LLM1.oeVel2)
describe_posterior(model.LLM1.oeVel2)
rope(model.LLM1.oeVel2)
posterior.2 = as_draws_df(model.LLM1.oeVel2)
1- mean(posterior.2$b_TreatmentWT > 0) # 0.10025

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.27 | [-0.61, 0.09] | 95.17% | [-0.02, 0.02] |     1.68% | 1.001 | 1413.00
TreatmentWT |   0.27 | [-0.24, 0.79] | 89.98% | [-0.02, 0.02] |     2.84% | 1.000 | 1430.00

# LLM1 - ΔVel2_IDD
model.LLM1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(1:3,16:18),])

summary(model.LLM1.ΔVel2_IDD)
posterior_summary(model.LLM1.ΔVel2_IDD)
describe_posterior(model.LLM1.ΔVel2_IDD)
rope(model.LLM1.ΔVel2_IDD)
posterior.3 = as_draws_df(model.LLM1.ΔVel2_IDD)
1- mean(posterior.3$b_TreatmentΔVel2_IDD > 0) # 0.88775

# Summary of Posterior Distribution 

Parameter          |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
  (Intercept)        | 2.36e-04 | [-0.70, 0.73] | 50.05% | [-0.04, 0.04] |    13.68% | 1.005 | 1221.00
TreatmentΔVel2_IDD |    -0.48 | [-1.54, 0.61] | 88.78% | [-0.04, 0.04] |     3.58% | 1.006 |  675.00

# AML1 - ΔVel1
model.AML1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_SXM[19:24,])

summary(model.AML1.ΔVel1)
posterior_summary(model.AML1.ΔVel1)
describe_posterior(model.AML1.ΔVel1)
rope(model.AML1.ΔVel1)
posterior.1a = as_draws_df(model.AML1.ΔVel1)
1- mean(posterior.1a$b_TreatmentΔVel1 > 0) # 0.15575

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
  (Intercept)    |  -0.02 | [-1.41, 1.48] | 51.30% | [-0.08, 0.08] |    12.79% | 1.002 | 1626.00
TreatmentΔVel1 |   0.86 | [-1.25, 2.82] | 84.42% | [-0.08, 0.08] |     4.63% | 1.001 | 1463.00

# AML1 - oeVel1
model.AML1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(19:21,25:27),])

summary(model.AML1.oeVel1)
posterior_summary(model.AML1.oeVel1)
describe_posterior(model.AML1.oeVel1)
rope(model.AML1.oeVel1)
posterior.2b = as_draws_df(model.AML1.oeVel1)
mean(posterior.2b$b_TreatmentWT > 0) # 0.019

# Summary of Posterior Distribution 

Parameter   | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
  (Intercept) |   2.05 | [ 0.84,  3.27] | 99.85% | [-0.13, 0.13] |        0% | 1.000 | 1654.00
TreatmentWT |  -2.05 | [-3.75, -0.36] | 98.60% | [-0.13, 0.13] |        0% | 1.000 | 1771.00

# AML1 - ΔVel2
model.AML1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(19:21,28:30),])

summary(model.AML1.ΔVel2)
posterior_summary(model.AML1.ΔVel2)
describe_posterior(model.AML1.ΔVel2)
rope(model.AML1.ΔVel2)
posterior.3b = as_draws_df(model.AML1.ΔVel2)
1- mean(posterior.3b$b_TreatmentΔVel2 > 0) # 0.971

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 5.24e-03 | [-0.89, 0.98] | 50.78% | [-0.08, 0.08] |    19.00% | 1.001 | 1645.00
TreatmentΔVel2 |    -1.23 | [-2.66, 0.07] | 97.10% | [-0.08, 0.08] |     1.05% | 1.000 | 1587.00

# AML1 - oeVel2
model.AML1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(19:21,31:33),])

summary(model.AML1.oeVel2)
posterior_summary(model.AML1.oeVel2)
describe_posterior(model.AML1.oeVel2)
rope(model.AML1.oeVel2)
posterior.4b = as_draws_df(model.AML1.oeVel2)
1- mean(posterior.4b$b_TreatmentWT > 0) # 0.931

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.97 | [-0.06, 1.96] | 96.97% | [-0.07, 0.07] |     1.11% | 1.003 | 1614.00
TreatmentWT |  -0.98 | [-2.43, 0.53] | 93.10% | [-0.07, 0.07] |     2.42% | 1.002 | 1645.00

# AML1 - ΔVel2_IDD
model.AML1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(19:21,34:36),])

summary(model.AML1.ΔVel2_IDD)
posterior_summary(model.AML1.ΔVel2_IDD)
describe_posterior(model.AML1.ΔVel2_IDD)
rope(model.AML1.ΔVel2_IDD)
posterior.5b = as_draws_df(model.AML1.ΔVel2_IDD)
1- mean(posterior.5b$b_TreatmentΔVel2_IDD > 0) # 0.9075

# Summary of Posterior Distribution 

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
  (Intercept)        |   0.02 | [-1.54, 1.72] | 51.00% | [-0.11, 0.11] |    13.97% | 1.001 | 2151.00
TreatmentΔVel2_IDD |  -1.29 | [-3.67, 0.87] | 90.75% | [-0.11, 0.11] |     3.84% | 1.000 | 2663.00

# NML1 - ΔVel1
model.NML1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_SXM[37:42,])

summary(model.NML1.ΔVel1)
posterior_summary(model.NML1.ΔVel1)
describe_posterior(model.NML1.ΔVel1)
rope(model.NML1.ΔVel1)
posterior.1b = as_draws_df(model.NML1.ΔVel1)
1- mean(posterior.1b$b_TreatmentΔVel1 > 0) # 0.002

# Summary of Posterior Distribution 

Parameter      |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | -1.96e-04 | [-0.28, 0.30] | 50.08% | [-0.05, 0.05] |    37.84% | 1.002 | 1334.00
TreatmentΔVel1 |      0.84 | [ 0.39, 1.25] | 99.50% | [-0.05, 0.05] |        0% | 1.002 | 1484.00

# NML1 - oeVel1
model.NML1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(37:39,43:45),])

summary(model.NML1.oeVel1)
posterior_summary(model.NML1.oeVel1)
describe_posterior(model.NML1.oeVel1)
rope(model.NML1.oeVel1)
posterior.2b = as_draws_df(model.NML1.oeVel1)
1- mean(posterior.2b$b_TreatmentWT > 0) # 0.40

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.03 | [-0.30, 0.23] | 62.62% | [-0.01, 0.01] |     9.29% | 1.002 | 1452.00
TreatmentWT |   0.03 | [-0.33, 0.40] | 59.33% | [-0.01, 0.01] |     6.66% | 1.001 | 1452.00

# NML1 - ΔVel2
model.NML1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(37:39,46:48),])

summary(model.NML1.ΔVel2)
posterior_summary(model.NML1.ΔVel2)
describe_posterior(model.NML1.ΔVel2)
rope(model.NML1.ΔVel2)
posterior.3b = as_draws_df(model.NML1.ΔVel2)
1- mean(posterior.3b$b_TreatmentΔVel2 > 0) # 0.1575

# Summary of Posterior Distribution 

Parameter      |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | -1.88e-03 | [-0.80, 0.76] | 50.38% | [-0.04, 0.04] |    12.34% | 1.000 | 1211.00
TreatmentΔVel2 |      0.44 | [-0.72, 1.54] | 84.25% | [-0.04, 0.04] |     4.05% | 1.000 | 1234.00

# NML1 - oeVel2
model.NML1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(37:39,49:51),])

summary(model.NML1.oeVel2)
posterior_summary(model.NML1.oeVel2)
describe_posterior(model.NML1.oeVel2)
rope(model.NML1.oeVel2)
posterior.4b = as_draws_df(model.NML1.oeVel2)
mean(posterior.4b$b_TreatmentWT > 0) # 0.32

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.14 | [-0.43, 0.67] | 74.20% | [-0.03, 0.03] |     8.05% | 1.003 | 1712.00
TreatmentWT |  -0.14 | [-0.95, 0.67] | 68.00% | [-0.03, 0.03] |     6.37% | 1.002 | 1950.00

# NML1 - ΔVel2_IDD
model.NML1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(37:39,52:54),])

summary(model.NML1.ΔVel2_IDD)
posterior_summary(model.NML1.ΔVel2_IDD)
describe_posterior(model.NML1.ΔVel2_IDD)
rope(model.NML1.ΔVel2_IDD)
posterior.5b = as_draws_df(model.NML1.ΔVel2_IDD)
mean(posterior.5b$b_TreatmentΔVel2_IDD > 0) # 0.9075

# Summary of Posterior Distribution 

Parameter          |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------------
  (Intercept)        | -1.13e-03 | [-0.29, 0.31] | 50.50% | [-0.02, 0.02] |    13.92% | 1.001 | 1491.00
TreatmentΔVel2_IDD |     -0.23 | [-0.69, 0.21] | 91.38% | [-0.02, 0.02] |     2.74% | 1.003 | 1377.00

# PDM ----
PCR_raw_PDM   = PCR_data %>% filter(Stage == "PDM")
PCR_raw_PDM$Treatment = as.factor(PCR_raw_PDM$Treatment)

# LLM1 - ΔVel1
model.LLM1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PDM[1:6,])

summary(model.LLM1.ΔVel1)
posterior_summary(model.LLM1.ΔVel1)
describe_posterior(model.LLM1.ΔVel1)
rope(model.LLM1.ΔVel1)
posterior = as_draws_df(model.LLM1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.002

# Summary of Posterior Distribution 

Parameter      |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | -1.65e-03 | [-0.57, 0.54] | 50.38% | [-0.11, 0.11] |    42.18% | 1.003 | 1329.00
TreatmentΔVel1 |      1.87 | [ 1.14, 2.67] | 99.80% | [-0.11, 0.11] |        0% | 1.002 | 1155.00

# LLM1 - oeVel1
model.LLM1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(1:3,7:9),])

summary(model.LLM1.oeVel1)
posterior_summary(model.LLM1.oeVel1)
describe_posterior(model.LLM1.oeVel1)
rope(model.LLM1.oeVel1)
posterior.1 = as_draws_df(model.LLM1.oeVel1)
1- mean(posterior.1$b_TreatmentWT > 0) # 0.02

# Summary of Posterior Distribution 

Parameter   | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
  (Intercept) |  -0.66 | [-1.14, -0.21] | 99.22% | [-0.04, 0.04] |        0% | 1.000 | 2541.00
TreatmentWT |   0.66 | [ 0.02,  1.37] | 97.75% | [-0.04, 0.04] |     0.47% | 1.001 | 2344.00

# LLM1 - ΔVel2
model.LLM1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(1:3,10:12),])

summary(model.LLM1.ΔVel2)
posterior_summary(model.LLM1.ΔVel2)
describe_posterior(model.LLM1.ΔVel2)
rope(model.LLM1.ΔVel2)
posterior.2 = as_draws_df(model.LLM1.ΔVel2)
1- mean(posterior.2$b_TreatmentΔVel2 > 0) # 0.025

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
  (Intercept)    |  -0.01 | [-0.65, 0.60] | 52.00% | [-0.06, 0.06] |    20.50% | 1.001 | 1957.00
TreatmentΔVel2 |   0.87 | [ 0.00, 1.90] | 97.50% | [-0.06, 0.06] |     0.18% | 1.000 | 1914.00

# LLM1 - oeVel2
model.LLM1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(1:3,13:15),])

summary(model.LLM1.oeVel2)
posterior_summary(model.LLM1.oeVel2)
describe_posterior(model.LLM1.oeVel2)
rope(model.LLM1.oeVel2)
posterior.3 = as_draws_df(model.LLM1.oeVel2)
1- mean(posterior.3$b_TreatmentWT > 0) # 0.43425

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.04 | [-0.69, 0.71] | 55.62% | [-0.03, 0.03] |    10.61% | 1.001 | 1316.00
TreatmentWT |   0.05 | [-0.93, 1.10] | 56.57% | [-0.03, 0.03] |     7.74% | 1.001 | 1800.00

# LLM1 - ΔVel2_IDD
model.LLM1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(1:3,16:18),])

summary(model.LLM1.ΔVel2_IDD)
posterior_summary(model.LLM1.ΔVel2_IDD)
describe_posterior(model.LLM1.ΔVel2_IDD)
rope(model.LLM1.ΔVel2_IDD)
posterior.4 = as_draws_df(model.LLM1.ΔVel2_IDD)
1- mean(posterior.4$b_TreatmentΔVel2_IDD > 0) # 0.43575

# Summary of Posterior Distribution 

Parameter          |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
  (Intercept)        | 3.15e-04 | [-0.36, 0.37] | 50.12% | [-0.02, 0.02] |    10.89% | 1.001 | 1866.00
TreatmentΔVel2_IDD |     0.03 | [-0.47, 0.54] | 56.43% | [-0.02, 0.02] |     7.50% | 1.003 | 1789.00

# AML1 - ΔVel1
model.AML1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PDM[19:24,])

summary(model.AML1.ΔVel1)
posterior_summary(model.AML1.ΔVel1)
describe_posterior(model.AML1.ΔVel1)
rope(model.AML1.ΔVel1)
posterior.a = as_draws_df(model.AML1.ΔVel1)
1- mean(posterior.a$b_TreatmentΔVel1 > 0) # 0.0085

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
  (Intercept)    |  -0.03 | [-1.46, 1.33] | 51.68% | [-0.16, 0.16] |    24.76% | 1.001 | 2137.00
TreatmentΔVel1 |   2.72 | [ 0.84, 4.71] | 99.15% | [-0.16, 0.16] |        0% | 1.000 | 2187.00

# AML1 - oeVel1
model.AML1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(19:21,25:27),])

summary(model.AML1.oeVel1)
posterior_summary(model.AML1.oeVel1)
describe_posterior(model.AML1.oeVel1)
rope(model.AML1.oeVel1)
posterior.1a = as_draws_df(model.AML1.oeVel1)
mean(posterior.1a$b_TreatmentWT > 0) # 0.2945

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.34 | [-0.81, 1.39] | 77.65% | [-0.06, 0.06] |     8.03% | 1.001 | 1811.00
TreatmentWT |  -0.33 | [-1.92, 1.19] | 70.55% | [-0.06, 0.06] |     6.42% | 1.002 | 1794.00

# AML1 - ΔVel2
model.AML1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(19:21,28:30),])

summary(model.AML1.ΔVel2)
posterior_summary(model.AML1.ΔVel2)
describe_posterior(model.AML1.ΔVel2)
rope(model.AML1.ΔVel2)
posterior.2a = as_draws_df(model.AML1.ΔVel2)
1- mean(posterior.2a$b_TreatmentΔVel2 > 0) # 0.0005

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 1.45e-03 | [-0.47, 0.50] | 50.25% | [-0.11, 0.11] |    49.37% | 1.003 | 2244.00
TreatmentΔVel2 |     2.03 | [ 1.35, 2.68] | 99.95% | [-0.11, 0.11] |        0% | 1.002 | 2267.00

# AML1 - oeVel2
model.AML1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(19:21,31:33),])

summary(model.AML1.oeVel2)
posterior_summary(model.AML1.oeVel2)
describe_posterior(model.AML1.oeVel2)
rope(model.AML1.oeVel2)
posterior.3a = as_draws_df(model.AML1.oeVel2)
1- mean(posterior.3a$b_TreatmentWT > 0) # 0.25225

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.60 | [-2.21, 1.13] | 80.58% | [-0.09, 0.09] |     7.08% | 1.001 | 2008.00
TreatmentWT |   0.62 | [-1.79, 3.18] | 74.78% | [-0.09, 0.09] |     6.39% | 1.000 | 2274.00

# AML1 - ΔVel2_IDD
model.AML1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(19:21,34:36),])

summary(model.AML1.ΔVel2_IDD)
posterior_summary(model.AML1.ΔVel2_IDD)
describe_posterior(model.AML1.ΔVel2_IDD)
rope(model.AML1.ΔVel2_IDD)
posterior.4a = as_draws_df(model.AML1.ΔVel2_IDD)
mean(posterior.4a$b_TreatmentΔVel2_IDD > 0) # 0.04225

# Summary of Posterior Distribution 

Parameter          |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
  (Intercept)        | 7.87e-03 | [-1.02, 1.14] | 50.68% | [-0.09, 0.09] |    17.68% | 1.001 | 1938.00
TreatmentΔVel2_IDD |    -1.24 | [-2.89, 0.25] | 95.78% | [-0.09, 0.09] |     1.42% | 1.001 | 1967.00

# NML1 - ΔVel1
model.NML1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PDM[37:42,])

summary(model.NML1.ΔVel1)
posterior_summary(model.NML1.ΔVel1)
describe_posterior(model.NML1.ΔVel1)
rope(model.NML1.ΔVel1)
posterior.1b = as_draws_df(model.NML1.ΔVel1)
1- mean(posterior.1b$b_TreatmentΔVel1 > 0) # 0.083

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
  (Intercept)    |  -0.01 | [-1.25, 1.16] | 50.75% | [-0.09, 0.09] |    14.63% | 1.001 | 1881.00
TreatmentΔVel1 |   1.09 | [-0.67, 2.77] | 91.70% | [-0.09, 0.09] |     2.47% | 1.000 | 2121.00

# NML1 - oeVel1
model.NML1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(37:39,43:45),])

summary(model.NML1.oeVel1)
posterior_summary(model.NML1.oeVel1)
describe_posterior(model.NML1.oeVel1)
rope(model.NML1.oeVel1)
posterior.2b = as_draws_df(model.NML1.oeVel1)
mean(posterior.2b$b_TreatmentWT > 0) # 0.09825

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.98 | [-0.23, 2.29] | 95.30% | [-0.08, 0.08] |     2.11% | 1.001 | 1881.00
TreatmentWT |  -1.02 | [-2.79, 0.74] | 90.18% | [-0.08, 0.08] |     2.89% | 1.000 | 1903.00

# NML1 - ΔVel2
model.NML1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(37:39,46:48),])

summary(model.NML1.ΔVel2)
posterior_summary(model.NML1.ΔVel2)
describe_posterior(model.NML1.ΔVel2)
rope(model.NML1.ΔVel2)
posterior.3b = as_draws_df(model.NML1.ΔVel2)
1- mean(posterior.3b$b_TreatmentΔVel2 > 0) # 0.0215

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 1.81e-03 | [-1.13, 1.03] | 50.30% | [-0.10, 0.10] |    22.63% | 1.001 | 1942.00
TreatmentΔVel2 |     1.65 | [ 0.07, 3.13] | 97.85% | [-0.10, 0.10] |     0.21% | 1.001 | 1827.00

# NML1 - oeVel2
model.NML1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(37:39,49:51),])

summary(model.NML1.oeVel2)
posterior_summary(model.NML1.oeVel2)
describe_posterior(model.NML1.oeVel2)
rope(model.NML1.oeVel2)
posterior.4b = as_draws_df(model.NML1.oeVel2)
mean(posterior.4b$b_TreatmentWT > 0) # 0.1145

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.81 | [-0.29, 1.94] | 94.62% | [-0.07, 0.07] |     2.05% | 1.000 | 1628.00
TreatmentWT |  -0.83 | [-2.43, 0.73] | 88.55% | [-0.07, 0.07] |     3.47% | 1.002 | 1340.00

# NML1 - ΔVel2_IDD
model.NML1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(37:39,52:54),])

summary(model.NML1.ΔVel2_IDD)
posterior_summary(model.NML1.ΔVel2_IDD)
describe_posterior(model.NML1.ΔVel2_IDD)
rope(model.NML1.ΔVel2_IDD)
posterior.5b = as_draws_df(model.NML1.ΔVel2_IDD)
mean(posterior.5b$b_TreatmentΔVel2_IDD > 0) # 0.06475

# Summary of Posterior Distribution 

Parameter          |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------------
  (Intercept)        | -3.54e-03 | [-1.02, 1.00] | 50.48% | [-0.08, 0.08] |    16.37% | 1.001 | 2097.00
TreatmentΔVel2_IDD |     -1.03 | [-2.54, 0.54] | 93.53% | [-0.08, 0.08] |     1.84% | 1.003 | 1527.00

# PLATE ----
PCR_raw_PLATE   = PCR_data %>% filter(Stage == "PLATE")
PCR_raw_PLATE$Treatment = as.factor(PCR_raw_PLATE$Treatment)

# LLM1 - ΔVel1
model.LLM1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[1:6,])

summary(model.LLM1.ΔVel1)
posterior_summary(model.LLM1.ΔVel1)
describe_posterior(model.LLM1.ΔVel1)
rope(model.LLM1.ΔVel1)
posterior = as_draws_df(model.LLM1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.000

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 2.95e-03 | [-0.67, 0.78] | 50.48% | [-0.23, 0.23] |    64.24% | 1.001 | 1840.00
TreatmentΔVel1 |     4.18 | [ 3.10, 5.15] |   100% | [-0.23, 0.23] |        0% | 1.000 | 1884.00

# LLM1 - oeVel1
model.LLM1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(1:3,7:9),])

summary(model.LLM1.oeVel1)
posterior_summary(model.LLM1.oeVel1)
describe_posterior(model.LLM1.oeVel1)
rope(model.LLM1.oeVel1)
posterior.1 = as_draws_df(model.LLM1.oeVel1)
mean(posterior.1$b_TreatmentWT > 0) # 0.4005

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.02 | [-0.16, 0.19] | 64.50% | [-0.01, 0.01] |     9.08% | 1.001 | 1745.00
TreatmentWT |  -0.02 | [-0.27, 0.24] | 59.95% | [-0.01, 0.01] |     6.26% | 1.003 | 1701.00

# LLM1 - ΔVel2
model.LLM1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(1:3,10:12),])

summary(model.LLM1.ΔVel2)
posterior_summary(model.LLM1.ΔVel2)
describe_posterior(model.LLM1.ΔVel2)
rope(model.LLM1.ΔVel2)
posterior.2 = as_draws_df(model.LLM1.ΔVel2)
1- mean(posterior.2$b_TreatmentΔVel2 > 0) # 0.000

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 2.58e-03 | [-0.20, 0.23] | 51.12% | [-0.14, 0.14] |    92.37% | 1.001 | 1604.00
TreatmentΔVel2 |     2.58 | [ 2.29, 2.87] |   100% | [-0.14, 0.14] |        0% | 1.001 | 1654.00

# LLM1 - oeVel2
model.LLM1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(1:3,13:15),])

summary(model.LLM1.oeVel2)
posterior_summary(model.LLM1.oeVel2)
describe_posterior(model.LLM1.oeVel2)
rope(model.LLM1.oeVel2)
posterior.3 = as_draws_df(model.LLM1.oeVel2)
1- mean(posterior.3$b_TreatmentWT > 0) # 0.08175

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.37 | [-0.79, 0.08] | 95.75% | [-0.03, 0.03] |     1.29% | 1.001 | 1133.00
TreatmentWT |   0.37 | [-0.27, 1.00] | 91.83% | [-0.03, 0.03] |     2.53% | 1.003 |  975.00

# LLM1 - ΔVel2_IDD
model.LLM1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(1:3,16:18),])

summary(model.LLM1.ΔVel2_IDD)
posterior_summary(model.LLM1.ΔVel2_IDD)
describe_posterior(model.LLM1.ΔVel2_IDD)
rope(model.LLM1.ΔVel2_IDD)
posterior.4 = as_draws_df(model.LLM1.ΔVel2_IDD)
mean(posterior.4$b_TreatmentΔVel2_IDD > 0) # 0.239

# Summary of Posterior Distribution 

Parameter          |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
  (Intercept)        | 6.25e-04 | [-0.47, 0.42] | 50.20% | [-0.02, 0.02] |    12.34% | 1.000 | 1558.00
TreatmentΔVel2_IDD |    -0.16 | [-0.81, 0.50] | 76.10% | [-0.02, 0.02] |     5.42% | 1.001 | 1375.00

# AML1 - ΔVel1
model.AML1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[19:24,])

summary(model.AML1.ΔVel1)
posterior_summary(model.AML1.ΔVel1)
describe_posterior(model.AML1.ΔVel1)
rope(model.AML1.ΔVel1)
posterior.a = as_draws_df(model.AML1.ΔVel1)
mean(posterior.a$b_TreatmentΔVel1 > 0) # 0.00425

# Summary of Posterior Distribution 

Parameter      |   Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | 4.20e-03 | [-0.36,  0.36] | 51.32% | [-0.06, 0.06] |    34.68% | 1.003 | 1576.00
TreatmentΔVel1 |    -0.97 | [-1.44, -0.48] | 99.58% | [-0.06, 0.06] |        0% | 1.004 | 1379.00

# AML1 - oeVel1
model.AML1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(19:21,25:27),])

summary(model.AML1.oeVel1)
posterior_summary(model.AML1.oeVel1)
describe_posterior(model.AML1.oeVel1)
rope(model.AML1.oeVel1)
posterior.1a = as_draws_df(model.AML1.oeVel1)
1 - mean(posterior.1a$b_TreatmentWT > 0) # 0.138

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.15 | [-0.37, 0.10] | 92.25% | [-0.01, 0.01] |     2.92% | 1.000 | 1462.00
TreatmentWT |   0.15 | [-0.22, 0.49] | 86.20% | [-0.01, 0.01] |     3.89% | 1.000 | 1303.00

# AML1 - ΔVel2
model.AML1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(19:21,28:30),])

summary(model.AML1.ΔVel2)
posterior_summary(model.AML1.ΔVel2)
describe_posterior(model.AML1.ΔVel2)
rope(model.AML1.ΔVel2)
posterior.2a = as_draws_df(model.AML1.ΔVel2)
mean(posterior.2a$b_TreatmentΔVel2 > 0) # 0.00275

# Summary of Posterior Distribution 

Parameter      |   Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | 2.03e-03 | [-0.34,  0.35] | 50.50% | [-0.06, 0.06] |    37.45% | 1.000 | 2120.00
TreatmentΔVel2 |    -1.05 | [-1.55, -0.53] | 99.72% | [-0.06, 0.06] |        0% | 1.001 | 2148.00

# AML1 - oeVel2
model.AML1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(19:21,31:33),])

summary(model.AML1.oeVel2)
posterior_summary(model.AML1.oeVel2)
describe_posterior(model.AML1.oeVel2)
rope(model.AML1.oeVel2)
posterior.3a = as_draws_df(model.AML1.oeVel2)
1- mean(posterior.3a$b_TreatmentWT > 0) # 0.29925

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.16 | [-0.70, 0.39] | 76.35% | [-0.03, 0.03] |     6.87% | 1.001 | 1666.00
TreatmentWT |   0.16 | [-0.66, 0.97] | 70.08% | [-0.03, 0.03] |     6.29% | 1.001 | 1564.00

# AML1 - ΔVel2_IDD
model.AML1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(19:21,34:36),])

summary(model.AML1.ΔVel2_IDD)
posterior_summary(model.AML1.ΔVel2_IDD)
describe_posterior(model.AML1.ΔVel2_IDD)
rope(model.AML1.ΔVel2_IDD)
posterior.4a = as_draws_df(model.AML1.ΔVel2_IDD)
mean(posterior.4a$b_TreatmentΔVel2_IDD > 0) # 0.3555

# Summary of Posterior Distribution 

Parameter          |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------------
  (Intercept)        | -2.79e-03 | [-1.13, 1.13] | 50.35% | [-0.05, 0.05] |    10.42% | 1.003 | 1271.00
TreatmentΔVel2_IDD |     -0.20 | [-1.90, 1.40] | 64.45% | [-0.05, 0.05] |     6.66% | 1.006 | 1131.00

# NML1 - ΔVel1
model.NML1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[37:42,])

summary(model.NML1.ΔVel1)
posterior_summary(model.NML1.ΔVel1)
describe_posterior(model.NML1.ΔVel1)
rope(model.NML1.ΔVel1)
posterior.1b = as_draws_df(model.NML1.ΔVel1)
1- mean(posterior.1b$b_TreatmentΔVel1 > 0) # 0.0065

# Summary of Posterior Distribution 

Parameter      |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
  (Intercept)    | 1.20e-03 | [-0.56, 0.52] | 50.22% | [-0.07, 0.07] |    31.13% | 1.001 | 1332.00
TreatmentΔVel1 |     1.23 | [ 0.52, 1.99] | 99.35% | [-0.07, 0.07] |        0% | 1.001 | 1718.00

# NML1 - oeVel1
model.NML1.oeVel1    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(37:39,43:45),])

summary(model.NML1.oeVel1)
posterior_summary(model.NML1.oeVel1)
describe_posterior(model.NML1.oeVel1)
rope(model.NML1.oeVel1)
posterior.2b = as_draws_df(model.NML1.oeVel1)
mean(posterior.2b$b_TreatmentWT > 0) # 0.1015

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |   0.23 | [-0.08, 0.52] | 94.83% | [-0.02, 0.02] |     2.24% | 1.000 | 1731.00
TreatmentWT |  -0.23 | [-0.63, 0.21] | 89.85% | [-0.02, 0.02] |     3.05% | 1.000 | 1626.00

# NML1 - ΔVel2
model.NML1.ΔVel2    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(37:39,46:48),])

summary(model.NML1.ΔVel2)
posterior_summary(model.NML1.ΔVel2)
describe_posterior(model.NML1.ΔVel2)
rope(model.NML1.ΔVel2)
posterior.3b = as_draws_df(model.NML1.ΔVel2)
1- mean(posterior.3b$b_TreatmentΔVel2 > 0) # 0.26775

# Summary of Posterior Distribution 

Parameter      |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
  (Intercept)    | -1.20e-03 | [-0.42, 0.46] | 50.20% | [-0.02, 0.02] |    10.66% | 1.001 | 1273.00
TreatmentΔVel2 |      0.13 | [-0.53, 0.71] | 73.22% | [-0.02, 0.02] |     6.45% | 1.001 | 1598.00

# NML1 - oeVel2
model.NML1.oeVel2    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(37:39,49:51),])

summary(model.NML1.oeVel2)
posterior_summary(model.NML1.oeVel2)
describe_posterior(model.NML1.oeVel2)
rope(model.NML1.oeVel2)
posterior.4b = as_draws_df(model.NML1.oeVel2)
1 - mean(posterior.4b$b_TreatmentWT > 0) # 0.07975

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
  (Intercept) |  -0.20 | [-0.45, 0.05] | 95.80% | [-0.02, 0.02] |     1.37% | 1.004 | 1430.00
TreatmentWT |   0.20 | [-0.13, 0.57] | 92.03% | [-0.02, 0.02] |     2.45% | 1.003 | 1270.00

# NML1 - ΔVel2_IDD
model.NML1.ΔVel2_IDD    = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(37:39,52:54),])

summary(model.NML1.ΔVel2_IDD)
posterior_summary(model.NML1.ΔVel2_IDD)
describe_posterior(model.NML1.ΔVel2_IDD)
rope(model.NML1.ΔVel2_IDD)
posterior.5b = as_draws_df(model.NML1.ΔVel2_IDD)
mean(posterior.5b$b_TreatmentΔVel2_IDD > 0) # 0.18225

# Summary of Posterior Distribution 

Parameter          |   Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
  (Intercept)        | 2.53e-03 | [-0.65, 0.73] | 50.42% | [-0.04, 0.04] |    12.92% | 1.001 | 1448.00
TreatmentΔVel2_IDD |    -0.33 | [-1.35, 0.62] | 81.77% | [-0.04, 0.04] |     4.89% | 1.001 | 1817.00

# Microesclerotia ----

# LLM1 ----

plates   = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/plates_LLM1.txt",dec=".")
model    = aov(Normalized ~ Biological.Rep, data = plates)

# Assumptions
# 2. Normality
norm = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
variance = ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Biological.Rep, data = plates)

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

# Testing a non-parametric model
stat.test = aov(Normalized ~ Biological.Rep, data = plates) %>%
  tukey_hsd()

write.csv(stat.test, "statistical_results/LLM1_microsclearotia.csv")

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_melanization = ggboxplot(plates, x = "Biological.Rep", y = "Normalized",
                                ylab = "Normalized melanization", xlab = "")  #+ 
#  stat_pvalue_manual(stat.test, label = "p.adj.signif", 
#                     y.position = c(1.7, 1.8, 1.9, 1.7,2.0,1.9))
figure_melanization

pdf("Figures/figure_melanization_LLM1.pdf",
    width=7,height=7*3/5)
print(figure_melanization)
dev.off()

# Figure assumptions

assumptions_LLM1_microsclerotia = ggarrange(norm, variance,
                                            labels = c("A", "B"),
                                            ncol = 2, nrow = 1)
assumptions_LLM1_microsclerotia

pdf("Figures/assumptions_LLM1_microsclerotia.pdf",
    width=7,height=7*3/5)
print(assumptions_LLM1_microsclerotia)
dev.off()

# AML1 dark ----

plates   = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/plates_AML1.txt",dec=".")
model    = aov(Normalized ~ Treatment, data = plates)

# Assumptions
# 2. Normality
norm_dark = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
variance_dark = ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Treatment, data = plates)

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

# Testing a non-parametric model
stat.test = aov(Normalized ~ Treatment, data = plates) %>%
  tukey_hsd()

write.csv(stat.test, "statistical_results/AML1_microsclearotia_dark.csv")

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_melanization = ggboxplot(plates, x = "Treatment", y = "Normalized",
                                ylab = "Normalized melanization", xlab = "")  

figure_melanization

pdf("Figures/figure_melanization_dark_AML1.pdf",
    width=6,height=6*3/5)
print(figure_melanization)
dev.off()

# AML1 light ----

plates_light = read_excel("datasets/Quantification AML1 light.xlsx",
                          sheet = "ESTADISTICAA")
model        = aov(Normalized ~ Biological.Rep, data = plates_light)

# Assumptions
# 2. Normality
norm_light = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
variance_light = ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Biological.Rep, data = plates_light)

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

# Testing a non-parametric model
stat.test = aov(Normalized ~ Treatment, data = plates_light) %>%
  tukey_hsd()

write.csv(stat.test, "statistical_results/AML1_microsclearotia_light.csv")

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_melanization_light = ggboxplot(plates_light, x = "Treatment", y = "Normalized",
                                      ylab = "Normalized melanization", xlab = "") #+ 
#stat_pvalue_manual(stat.test, label = "p.adj.signif", 
#y.position = c(1.7, 1.8, 1.9))
figure_melanization_light

pdf("Figures/figure_melanization_light_AML1.pdf",
    width=6,height=6*3/5)
print(figure_melanization_light)
dev.off()

# Figure assumptions

assumptions_AML1_microsclerotia = ggarrange(norm_dark, variance_dark,
                                            norm_light,variance_light,
                                            labels = c("A", "B",
                                                       "C", "D"),
                                            ncol = 2, nrow = 2)
assumptions_AML1_microsclerotia

pdf("Figures/assumptions_AML1_microsclerotia.pdf",
    width=6,height=6*3/5)
print(assumptions_AML1_microsclerotia)
dev.off()

# NML1 ----

nml1_stress   = read_excel("datasets/New quantification.xlsx",
                           sheet = "Sheet1")
# Barplot

nml1_stress_1 = nml1_stress %>% group_by(Treatment,Plot) %>%
  summarize(Normalized_mean = mean(Normalized),
            Normalized_sd   = sd(Normalized))

NML1_microsclerotia = ggplot(nml1_stress_1, aes(x=Plot, y=Normalized_mean, fill=Plot)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=Normalized_mean-Normalized_sd, ymax=Normalized_mean+Normalized_sd), width=.2,
                position=position_dodge(.9)) + scale_fill_manual(values=c('#525252','#969696',"#d9d9d9")) + 
  theme_classic()
NML1_microsclerotia

pdf("Figures/NML1_microsclerotia.pdf",
    width=5,height=5*3/5)
print(NML1_microsclerotia)
dev.off()

# AML1 glucosa dark ----

plant_glucose = read_excel("datasets/Glucose AML1 quantification.xlsx", 
                           sheet = "Statistics")
model         = aov(Normalized ~ Treatment, data = plant_glucose)

# Assumptions
# 2. Normality
norm_AML1_glu = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
var_AML1_glu = ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Treatment, data = plant_glucose)

# Figure assumptions

assumptions_AML1_glu_microsclerotia = ggarrange(norm_AML1_glu, var_AML1_glu,
                                                labels = c("A", "B"),
                                                ncol = 2, nrow = 1)
assumptions_AML1_glu_microsclerotia

pdf("Figures/assumptions_AML1_glu_microsclerotia.pdf",
    width=6,height=6*3/5)
print(assumptions_AML1_glu_microsclerotia)
dev.off()

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

# Testing a non-parametric model
stat.test = aov(Normalized ~ Treatment, data = plant_glucose) %>%
  tukey_hsd()
write.csv(stat.test, "statistical_results/stat.test_AML1_microsclerotia_glu.csv")

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_melanization_glucose = ggboxplot(plant_glucose, x = "Treatment", y = "Normalized",
                                        ylab = "Normalized melanization", xlab = "") #+ 
#  stat_pvalue_manual(stat.test, label = "p.adj.signif", 
#                     y.position = c(1.7, 1.8, 1.9, 2.1, 2.3))
figure_melanization_glucose

pdf("Figures/figure_AML1_glucose.pdf",
    width=6,height=6*3/5)
print(figure_melanization_glucose)
dev.off()

# Radius ----

# AML1 Glucose ----

AML1_colony_size = read_excel("datasets/AML1 colony size in glucose.xlsx",
                              sheet = "statistics")
model            = aov(Radio ~ Condition*Treatment, data = AML1_colony_size)

# Assumptions 
# 2. Normality
norm_AML1_glu = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
var_AML1_glu  = ols_plot_resid_fit(model)
# bartlett.test(Radio ~ Condition:Treatment, data = AML1_colony_size)

# Figure assumptions

assumptions_AML1_glu_microsclerotia = ggarrange(norm_AML1_glu, var_AML1_glu,
                                            labels = c("A", "B"),
                                            ncol = 2, nrow = 1)
assumptions_AML1_glu_microsclerotia

pdf("Figures/assumptions_AML1_glu_radius.pdf",
    width=6,height=6*3/5)
print(assumptions_AML1_glu_microsclerotia)
dev.off()

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

stat.test = aov(Radio ~ Condition*Treatment, data = AML1_colony_size) %>%
  tukey_hsd()

write.csv(stat.test, "statistical_results/stat.test_AML1_radius_glu.csv")

# Barplot

AML1_colony_size_1 = AML1_colony_size %>% group_by(Condition,Treatment,Plots) %>%
  summarize(radio_mean = mean(Radio),
            radio_sd   = sd(Radio))

radius_AML1 = ggplot(AML1_colony_size_1, aes(x=Plots, y=radio_mean, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=radio_mean-radio_sd, ymax=radio_mean+radio_sd), width=.2,
                position=position_dodge(.9)) + scale_fill_manual(values=c('#525252','lightgray')) + 
  theme_classic()
radius_AML1

pdf("Figures/radius_AML1.pdf",
    width=7,height=7*3/5)
print(radius_AML1)
dev.off()

# AML1 Sucrose ----

AML1_colony_sucrose = read_excel("datasets/AML1 colony size in sucrose.xlsx",
                              sheet = "statistics")
model = brm(
  Radius ~ Treatment * Condition,
  data = AML1_colony_sucrose,
  family = Gamma(link = "log"),  # or try gaussian if Gamma feels too sensitive
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(gamma(2, 0.1), class = "shape")  # adjust based on your data
  ),
  chains = 4, iter = 4000, warmup = 500
)

summary(model)
posterior_summary(model)

# Post-hoc comparisons

emm = emmeans(model, ~ Condition * Treatment)
pairs(emm)

samples = as.data.frame(posterior_samples(emm))

# Plot distributions for selected groups
ggplot(samples, aes(x = Treatment, y = emmean, color = Condition)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD),
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(title = "Posterior Means with 95% Credible Intervals",
       y = "Estimated Radius", x = "Treatment") +
  theme_minimal()

# Barplot

AML1_colony_size_1 = AML1_colony_sucrose %>% group_by(Condition,Treatment,Plot) %>%
  summarize(radio_mean = mean(Radius),
            radio_sd   = sd(Radius))

radius_AML1 = ggplot(AML1_colony_size_1, aes(x=Plot, y=radio_mean, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=radio_mean-radio_sd, ymax=radio_mean+radio_sd), width=.2,
                position=position_dodge(.9)) + scale_fill_manual(values=c('#525252','lightgray')) + 
  theme_classic()
radius_AML1

pdf("Figures/radius_AML1_sucrose.pdf",
    width=5,height=5*3/5)
print(radius_AML1)
dev.off()

# Sporas ----

# LLM1 ----

# Call data
spores_LLM1 = read_excel("datasets/Spores three main genes.xlsx",
                         sheet = "LLM1 spores")

# Statistical model
model  = aov(Spores ~ Treatment, data = spores_LLM1)
summary(model)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model)
# 1. Homogeneity of variances
ols_plot_resid_fit(model)
bartlett.test(Spores ~ Treatment, data = spores_LLM1)

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

stat.test = aov(Spores ~ Plot, data = spores_LLM1) %>%
  tukey_hsd()

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_spores_LLM1 = ggbarplot(spores_LLM1, x = "Plot", y = "Spores", 
                               ylab = "Normalized conidia formation", xlab = "", add = "mean_se", fill = "#525252") + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                     y.position = c(1.1, 1.2, 1.3, 1.4)) + 
  theme_classic()
figure_spores_LLM1

pdf("Figures/figure_spores_LLM1.pdf",
    width=6,height=6*3/5)
print(figure_spores_LLM1)
dev.off()

# AML1 ----

# Call data
spores_AML1 = read_excel("datasets/Spores three main genes.xlsx",
                         sheet = "AML1 spores")

# Statistical model
model  = aov(Spores ~ Treatment, data = spores_AML1)
summary(model)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model)
# 1. Homogeneity of variances
ols_plot_resid_fit(model)
bartlett.test(Spores ~ Treatment, data = spores_AML1)

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

stat.test = aov(Spores ~ Plot, data = spores_AML1) %>%
  tukey_hsd()

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_spores_AML1 = ggbarplot(spores_AML1, x = "Plot", y = "Spores", 
                               ylab = "Normalized conidia formation", xlab = "", add = "mean_se", fill = "#525252") + 
  theme_classic()
figure_spores_AML1

pdf("Figures/figure_spores_AML1.pdf",
    width=6,height=6*3/5)
print(figure_spores_AML1)
dev.off()

# NML1 ----

# Call data
spores_NML1 = read_excel("datasets/Spores three main genes.xlsx",
                         sheet = "NML1 spores")

# Statistical model
model  = aov((Spores) ~ Treatment, data = spores_NML1)
summary(model)

# Assumptions (not normal and variances no homogenous)
# 2. Normality
ols_plot_resid_qq(model)
# 1. Homogeneity of variances
ols_plot_resid_fit(model)
bartlett.test((Spores) ~ Treatment, data = spores_NML1)

# Non-parametric alternative (Kruskal-Wallis Test)

kruskal.test((Spores) ~ Treatment, data = spores_NML1)

# Pairwise interactions

stat.test = pairwise.wilcox.test(spores_NML1$Spores, spores_NML1$Treatment,
                                 p.adjust.method = "none")

# Figure

figure_spores_NML1 = ggbarplot(spores_NML1, x = "Plot", y = "Spores", 
                               ylab = "Normalized conidia formation", xlab = "", add = "mean_se", fill = "#525252") + 
  theme_classic()
figure_spores_NML1

pdf("Figures/figure_spores_NML1.pdf",
    width=6,height=6*3/5)
print(figure_spores_NML1)
dev.off()
