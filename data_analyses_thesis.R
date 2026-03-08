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
library(bayestestR)

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

# Preliminary PCR results ----
PCR_preliminary = read_excel("datasets/OE summary Real Time.xlsx", 
                             sheet = "OE")
# Transform to Log2
PCR_preliminary = PCR_preliminary %>%  mutate(across(c(Expression), function(x) log2(x)))

# Delete WT
PCR_preliminary   = PCR_preliminary %>% filter(Treatment != "WT")
PCR_preliminary.1 = PCR_preliminary %>% group_by(Stage) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))
# Plot
Figure_PCR_preliminary = ggplot(PCR_preliminary.1, aes(x=Stage, y=avg, fill=as.factor(Stage))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#e34a33", "#045a8d", "#74c476"), labels = c("PDM", "PLATE", "SXM")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_PCR_preliminary

pdf("Figures/Figure_PCR_preliminary.pdf",
    width=12,height=12*3/5)
print(Figure_PCR_preliminary)
dev.off()

# Bayesian Methods for Group Comparison ----

# SXM ----
PCR_raw_SXM   = PCR_raw %>% filter(Stage == "SXM")
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
norm_LLM1_conidia = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
var_LLM1_conidia  = ols_plot_resid_fit(model)
bartlett.test(Spores ~ Treatment, data = spores_LLM1)

# Figure assumptions

assumptions_LLM1_conidia = ggarrange(norm_LLM1_conidia, var_LLM1_conidia,
                                                labels = c("A", "B"),
                                                ncol = 2, nrow = 1)
assumptions_LLM1_conidia

pdf("Figures/assumptions_LLM1_conidia.pdf",
    width=6,height=6*3/5)
print(assumptions_LLM1_conidia)
dev.off()

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

stat.test = aov(Spores ~ Plot, data = spores_LLM1) %>%
  tukey_hsd()

write.csv(stat.test, "statistical_results/stat.test_LLM1_conidia.csv")

# stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_spores_LLM1 = ggbarplot(spores_LLM1, x = "Plot", y = "Spores", 
                               ylab = "Normalized conidia formation", xlab = "", add = "mean_se", fill = "#525252") + 
  #  stat_pvalue_manual(stat.test, label = "p.adj.signif", 
  #                     y.position = c(1.1, 1.2, 1.3, 1.4)) + 
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
norm_AML1_conidia = ols_plot_resid_qq(model)
# 1. Homogeneity of variances
var_AML1_conidia  = ols_plot_resid_fit(model)
bartlett.test(Spores ~ Treatment, data = spores_AML1)

# Figure assumptions

assumptions_AML1_conidia = ggarrange(norm_AML1_conidia, var_AML1_conidia,
                                     labels = c("A", "B"),
                                     ncol = 2, nrow = 1)
assumptions_AML1_conidia

pdf("Figures/assumptions_AML1_conidia.pdf",
    width=6,height=6*3/5)
print(assumptions_AML1_conidia)
dev.off()

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure

stat.test = aov(Spores ~ Plot, data = spores_AML1) %>%
  tukey_hsd()

write.csv(stat.test, "statistical_results/stat.test_AML1_conidia.csv")

# stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

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

# Plant health ----

# LLM1 ----
LLM1_plant_main = read_excel("datasets/ALL DATA STATISTICS_plant.xlsx",
                             sheet = "LLM1_stats_main")
# Aligned ranks anova 
# https://rcompanion.org/handbook/F_16.html

LLM1_plant_main$Plot = factor(LLM1_plant_main$Plot,
                              levels=unique(LLM1_plant_main$Plot))

model = art(as.numeric(Rank) ~ Plot,data = LLM1_plant_main)
anova(model)

# Post-hoc comparisons 

model.lm = artlm(model, "Plot")
marginal = emmeans(model.lm,~ Plot)
test     = pairs(marginal,adjust = "tukey")

# Reorganizing data for plotting

LLM1_plant_main_1  = LLM1_plant_main %>% count(Plot,Rank)
LLM1_plant_main_2  = LLM1_plant_main %>% count(Plot)

levels                = unique(LLM1_plant_main_1$Plot)
LLM1_plant_main_plot  = c()

for (i in levels){
  temp_1 = LLM1_plant_main_1  %>% filter(Plot == i)
  temp_2 = LLM1_plant_main_2 %>% filter(Plot == i)
  temp_3 = temp_1 %>% mutate(Percentage = temp_1$n/temp_2$n)
  LLM1_plant_main_plot = rbind(LLM1_plant_main_plot,temp_3)
}

# Plot

LLM1_plant_main_plot$Rank = factor(LLM1_plant_main_plot$Rank,
                                   levels=unique(LLM1_plant_main_plot$Rank))
LLM1_plant_main_plot$Plot = factor(LLM1_plant_main_plot$Plot,
                                   levels=unique(LLM1_plant_main_plot$Plot))

figure_plant_LLM1 = ggplot(LLM1_plant_main_plot, aes(Plot, Percentage, fill = (Rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant_LLM1

# Save

pdf("Figures/figure_plant_LLM1.pdf",
    width=8,height=8*3/5)
print(figure_plant_LLM1)
dev.off()

# AML1 ----
AML1_plant_main = read_excel("datasets/ALL DATA STATISTICS_plant.xlsx",
                             sheet = "AML1_stats")
# Aligned ranks anova 
# https://rcompanion.org/handbook/F_16.html

AML1_plant_main$Plot = factor(AML1_plant_main$Plot,
                              levels=unique(AML1_plant_main$Plot))

model = art(as.numeric(Rank) ~ Plot,data = AML1_plant_main)
anova(model)

# Post-hoc comparisons 

model.lm = artlm(model, "Plot")
marginal = emmeans(model.lm,~ Plot)
test     = pairs(marginal,adjust = "tukey")

# Reorganizing data for plotting

AML1_plant_main_1  = AML1_plant_main %>% count(Plot,Rank)
AML1_plant_main_2  = AML1_plant_main %>% count(Plot)

levels                = unique(AML1_plant_main_1$Plot)
AML1_plant_main_plot  = c()

for (i in levels){
  temp_1 = AML1_plant_main_1  %>% filter(Plot == i)
  temp_2 = AML1_plant_main_2 %>% filter(Plot == i)
  temp_3 = temp_1 %>% mutate(Percentage = temp_1$n/temp_2$n)
  AML1_plant_main_plot = rbind(AML1_plant_main_plot,temp_3)
}

# Plot

AML1_plant_main_plot$Rank = factor(AML1_plant_main_plot$Rank,
                                   levels=unique(AML1_plant_main_plot$Rank))
AML1_plant_main_plot$Plot = factor(AML1_plant_main_plot$Plot,
                                   levels=unique(AML1_plant_main_plot$Plot))

figure_plant_AML1 = ggplot(AML1_plant_main_plot, aes(Plot, Percentage, fill = (Rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant_AML1

# Save

pdf("Figures/figure_plant_AML1.pdf",
    width=8,height=8*3/5)
print(figure_plant_AML1)
dev.off()

# NML1 ----
NML1_plant_main = read_excel("datasets/ALL DATA STATISTICS_plant.xlsx",
                             sheet = "NML1_stats")
# Aligned ranks anova 
# https://rcompanion.org/handbook/F_16.html

NML1_plant_main$Plot = factor(NML1_plant_main$Plot,
                              levels=unique(NML1_plant_main$Plot))

model = art(as.numeric(Rank) ~ Plot,data = NML1_plant_main)
anova(model)

# Post-hoc comparisons 

model.lm = artlm(model, "Plot")
marginal = emmeans(model.lm,~ Plot)
test     = pairs(marginal,adjust = "tukey")

# Reorganizing data for plotting

NML1_plant_main_1  = NML1_plant_main %>% count(Plot,Rank)
NML1_plant_main_2  = NML1_plant_main %>% count(Plot)

levels                = unique(NML1_plant_main_1$Plot)
NML1_plant_main_plot  = c()

for (i in levels){
  temp_1 = NML1_plant_main_1  %>% filter(Plot == i)
  temp_2 = NML1_plant_main_2 %>% filter(Plot == i)
  temp_3 = temp_1 %>% mutate(Percentage = temp_1$n/temp_2$n)
  NML1_plant_main_plot = rbind(NML1_plant_main_plot,temp_3)
}

# Plot

NML1_plant_main_plot$Rank = factor(NML1_plant_main_plot$Rank,
                                   levels=unique(NML1_plant_main_plot$Rank))
NML1_plant_main_plot$Plot = factor(NML1_plant_main_plot$Plot,
                                   levels=unique(NML1_plant_main_plot$Plot))

figure_plant_NML1 = ggplot(NML1_plant_main_plot, aes(Plot, Percentage, fill = (Rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant_NML1

# Save

pdf("Figures/figure_plant_NML1.pdf",
    width=8,height=8*3/5)
print(figure_plant_NML1)
dev.off()

# LLM1 SI ----
LLM1_plant_main = read_excel("datasets/ALL DATA STATISTICS_plant.xlsx",
                             sheet = "LLM1_stats_SI")
# Aligned ranks anova 
# https://rcompanion.org/handbook/F_16.html

LLM1_plant_main$Plot = factor(LLM1_plant_main$Plot,
                              levels=unique(LLM1_plant_main$Plot))

model = art(as.numeric(Rank) ~ Plot,data = LLM1_plant_main)
anova(model)

# Post-hoc comparisons 

model.lm = artlm(model, "Plot")
marginal = emmeans(model.lm,~ Plot)
test     = pairs(marginal,adjust = "tukey")

# Reorganizing data for plotting

LLM1_plant_main_1  = LLM1_plant_main %>% count(Plot,Rank)
LLM1_plant_main_2  = LLM1_plant_main %>% count(Plot)

levels                = unique(LLM1_plant_main_1$Plot)
LLM1_plant_main_plot  = c()

for (i in levels){
  temp_1 = LLM1_plant_main_1  %>% filter(Plot == i)
  temp_2 = LLM1_plant_main_2 %>% filter(Plot == i)
  temp_3 = temp_1 %>% mutate(Percentage = temp_1$n/temp_2$n)
  LLM1_plant_main_plot = rbind(LLM1_plant_main_plot,temp_3)
}

# Plot

LLM1_plant_main_plot$Rank = factor(LLM1_plant_main_plot$Rank,
                                   levels=unique(LLM1_plant_main_plot$Rank))
LLM1_plant_main_plot$Plot = factor(LLM1_plant_main_plot$Plot,
                                   levels=unique(LLM1_plant_main_plot$Plot))

figure_plant_LLM1_sup = ggplot(LLM1_plant_main_plot, aes(Plot, Percentage, fill = (Rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant_LLM1_sup

# Save

pdf("Figures/figure_plant_LLM1_sup.pdf",
    width=6,height=6*3/5)
print(figure_plant_LLM1_sup)
dev.off()

# Other genes ----
OTHER_plant_main = read_excel("datasets/ALL DATA STATISTICS_plant.xlsx",
                              sheet = "OTHER_stats")
# Aligned ranks anova 
# https://rcompanion.org/handbook/F_16.html

OTHER_plant_main$Plot = factor(OTHER_plant_main$Plot,
                               levels=unique(OTHER_plant_main$Plot))

model = art(as.numeric(Rank) ~ Plot,data = OTHER_plant_main)
anova(model)

# Post-hoc comparisons 

model.lm = artlm(model, "Plot")
marginal = emmeans(model.lm,~ Plot)
test     = pairs(marginal,adjust = "tukey")

# Reorganizing data for plotting

OTHER_plant_main_1  = OTHER_plant_main %>% count(Plot,Rank)
OTHER_plant_main_2  = OTHER_plant_main %>% count(Plot)

levels                = unique(OTHER_plant_main_1$Plot)
OTHER_plant_main_plot  = c()

for (i in levels){
  temp_1 = OTHER_plant_main_1  %>% filter(Plot == i)
  temp_2 = OTHER_plant_main_2 %>% filter(Plot == i)
  temp_3 = temp_1 %>% mutate(Percentage = temp_1$n/temp_2$n)
  OTHER_plant_main_plot = rbind(OTHER_plant_main_plot,temp_3)
}

# Plot

OTHER_plant_main_plot$Rank = factor(OTHER_plant_main_plot$Rank,
                                    levels=unique(OTHER_plant_main_plot$Rank))
OTHER_plant_main_plot$Plot = factor(OTHER_plant_main_plot$Plot,
                                    levels=unique(OTHER_plant_main_plot$Plot))

figure_plant_OTHER_sup = ggplot(OTHER_plant_main_plot, aes(Plot, Percentage, fill = (Rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant_OTHER_sup

# Save

pdf("Figures/figure_plant_OTHER_sup.pdf",
    width=8,height=8*3/5)
print(figure_plant_OTHER_sup)
dev.off()

# PCR Complex ----

# Call data
vel_expression     = read_excel("datasets/Summary Velvet expression in LLM1.xlsx",
                            sheet = "statistics")
vel_expression_raw = vel_expression

# Transform to Log2
vel_expression = vel_expression %>%  mutate(across(c(Expression), 
                                                   function(x) log2(x)))

# Delete WT
vel_expression = vel_expression %>% filter(Treatment != "WT")

# Vel1 ----

vel_expression_vel1      = vel_expression %>% filter(Protein == "VEL1")

vel_expression_vel1_plot = vel_expression_vel1 %>% group_by(Stage,Treatment,Plot) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_expression_vel1 = ggplot(vel_expression_vel1_plot, aes(x=Stage, y=avg, fill=as.factor(Plot))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#e5f5e0", "#a1d99b", "#238b45", "#084594"), 
                    labels = c("ΔLLM1", "LLM1-C", "oeLLM1","ΔVEL2")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_expression_vel1

pdf("Figures/Figure_expression_vel1.pdf",
    width=6,height=6*3/5)
print(Figure_expression_vel1)
dev.off()

# Vel2 ----

vel_expression_vel2      = vel_expression %>% filter(Protein == "VEL2")

vel_expression_vel2_plot = vel_expression_vel2 %>% group_by(Stage,Treatment,Plot) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_expression_vel2 = ggplot(vel_expression_vel2_plot, aes(x=Stage, y=avg, fill=as.factor(Plot))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#e5f5e0", "#a1d99b", "#238b45", "#fff7bc","#7bccc4","#084594"), 
                    labels = c("ΔLLM1", "LLM1-C", "oeLLM1","ΔVEL1ΔLLM1","ΔVEL1oeLLM1","ΔVEL1")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_expression_vel2

pdf("Figures/Figure_expression_vel2.pdf",
    width=6,height=6*3/5)
print(Figure_expression_vel2)
dev.off()

# Bayesian Methods for Group Comparison ----

# VEL1 ----

# SXM ----

vel1_expression_raw_SXM  = vel_expression_raw %>% filter(Stage == "SXM" & Protein == "VEL1")
vel1_expression_raw_SXM$Treatment    = as.factor(vel1_expression_raw_SXM$Treatment)

# WT - ΔLLM1 ----

model.vel1.SXM.ΔLLM1     = brm(Expression ~ Treatment, data = vel1_expression_raw_SXM[1:6,])

posterior_summary(model.vel1.SXM.ΔLLM1)
describe_posterior(model.vel1.SXM.ΔLLM1)
rope(model.vel1.SXM.ΔLLM1)
posterior = as_draws_df(model.vel1.SXM.ΔLLM1)
1- mean(posterior$b_TreatmentΔLLM1 > 0) # 0.09375

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.79, 1.20] |   100% | [-0.01, 0.01] |        0% | 1.000 | 1552.00
TreatmentΔLLM1 |   0.15 | [-0.14, 0.42] | 90.62% | [-0.01, 0.01] |     2.71% | 1.000 | 1855.00

# WT - LLM1-C ----

model.vel1.SXM.LLM1.C     = brm(Expression ~ Treatment, data = vel1_expression_raw_SXM[c(1:3,7:9),])

posterior_summary(model.vel1.SXM.LLM1.C)
describe_posterior(model.vel1.SXM.LLM1.C)
rope(model.vel1.SXM.LLM1.C)
posterior = as_draws_df(model.vel1.SXM.LLM1.C)
1- mean(posterior$b_TreatmentWT > 0) # 0.086

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   0.89 | [ 0.75, 1.03] |   100% | [-0.01, 0.01] |        0% | 1.000 | 1989.00
TreatmentWT |   0.11 | [-0.08, 0.32] | 91.40% | [-0.01, 0.01] |     2.63% | 1.000 | 1651.00

# WT - oeLLM1 ----

model.vel1.SXM.oeLLM1     = brm(Expression ~ Treatment, data = vel1_expression_raw_SXM[c(1:3,10:12),])

posterior_summary(model.vel1.SXM.oeLLM1)
describe_posterior(model.vel1.SXM.oeLLM1)
rope(model.vel1.SXM.oeLLM1)
posterior = as_draws_df(model.vel1.SXM.oeLLM1)
mean(posterior$b_TreatmentWT > 0) # 0.0335

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   1.24 | [ 1.06, 1.45] |   100% | [-0.02, 0.02] |        0% | 1.000 | 2149.00
TreatmentWT |  -0.25 | [-0.54, 0.04] | 96.65% | [-0.02, 0.02] |     1.13% | 1.001 | 2328.00

# WT - ΔVEL2 ----

model.vel1.SXM.ΔVEL2     = brm(Expression ~ Treatment, data = vel1_expression_raw_SXM[c(1:3,13:15),])

posterior_summary(model.vel1.SXM.ΔVEL2)
describe_posterior(model.vel1.SXM.ΔVEL2)
rope(model.vel1.SXM.ΔVEL2)
posterior = as_draws_df(model.vel1.SXM.ΔVEL2)
1-mean(posterior$b_TreatmentΔVEL2  > 0) # 0.24325

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.72, 1.32] | 99.92% | [-0.01, 0.01] |        0% | 1.000 | 1883.00
TreatmentΔVEL2 |   0.11 | [-0.33, 0.53] | 75.67% | [-0.01, 0.01] |     5.82% | 1.002 | 1753.00

# PDM ----

vel1_expression_raw_PDM  = vel_expression_raw %>% filter(Stage == "PDM" & Protein == "VEL1")
vel1_expression_raw_PDM$Treatment    = as.factor(vel1_expression_raw_PDM$Treatment)

# WT - ΔLLM1 ----

model.vel1.PDM.ΔLLM1     = brm(Expression ~ Treatment, data = vel1_expression_raw_PDM[c(1:8),],
                               iter = 4000)

posterior_summary(model.vel1.PDM.ΔLLM1)
describe_posterior(model.vel1.PDM.ΔLLM1)
rope(model.vel1.PDM.ΔLLM1)
posterior = as_draws_df(model.vel1.PDM.ΔLLM1)
1-mean(posterior$b_TreatmentΔLLM1  > 0) # 0.0395

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.38, 1.59] | 99.60% | [-0.06, 0.06] |        0% | 1.000 | 2450.00
TreatmentΔLLM1 |   0.73 | [-0.18, 1.60] | 96.05% | [-0.06, 0.06] |     1.76% | 1.000 | 2143.00

# WT - LLM1-C ----

model.vel1.PDM.LLM1.C     = brm(Expression ~ Treatment, data = vel1_expression_raw_PDM[c(1:4,9:12),])

posterior_summary(model.vel1.PDM.LLM1.C)
describe_posterior(model.vel1.PDM.LLM1.C)
rope(model.vel1.PDM.LLM1.C)
posterior = as_draws_df(model.vel1.PDM.LLM1.C)
mean(posterior$b_TreatmentWT  > 0) # 0.25525

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   1.25 | [ 0.64, 1.84] | 99.85% | [-0.04, 0.04] |        0% | 1.000 | 2365.00
TreatmentWT |  -0.24 | [-1.09, 0.61] | 74.48% | [-0.04, 0.04] |     7.26% | 1.001 | 2199.00

# WT - oeLLM1 ----

model.vel1.PDM.oeLLM1     = brm(Expression ~ Treatment, data = vel1_expression_raw_PDM[c(1:4,13:16),],
                                iter = 4000)

posterior_summary(model.vel1.PDM.oeLLM1)
describe_posterior(model.vel1.PDM.oeLLM1)
rope(model.vel1.PDM.oeLLM1)
posterior = as_draws_df(model.vel1.PDM.oeLLM1)
mean(posterior$b_TreatmentWT  > 0) # 0.00

# Summary of Posterior Distribution 

Parameter   | Median |         95% CI |   pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------
(Intercept) |   2.64 | [ 2.47,  2.82] | 100% | [-0.09, 0.09] |        0% | 1.000 | 4151.00
TreatmentWT |  -1.64 | [-1.88, -1.40] | 100% | [-0.09, 0.09] |        0% | 1.000 | 4444.00

# WT - ΔVEL2 ----

model.vel1.PDM.ΔVEL2     = brm(Expression ~ Treatment, data = vel1_expression_raw_PDM[c(1:4,17:19),],
                               iter = 4000)

posterior_summary(model.vel1.PDM.ΔVEL2)
describe_posterior(model.vel1.PDM.ΔVEL2)
rope(model.vel1.PDM.ΔVEL2)
posterior = as_draws_df(model.vel1.PDM.ΔVEL2)
mean(posterior$b_TreatmentΔVEL2   > 0) # 0.180125

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.78, 1.20] |   100% | [-0.01, 0.01] |        0% | 1.000 | 4226.00
TreatmentΔVEL2 |  -0.12 | [-0.45, 0.21] | 81.99% | [-0.01, 0.01] |     5.33% | 1.000 | 3857.00

# PLATE ----

vel1_expression_raw_PLATE  = vel_expression_raw %>% filter(Stage == "PLATE" & Protein == "VEL1")
vel1_expression_raw_PLATE$Treatment    = as.factor(vel1_expression_raw_PLATE$Treatment)

# WT - ΔLLM1 ----

model.vel1.PLATE.ΔLLM1     = brm(Expression ~ Treatment, data = vel1_expression_raw_PLATE[c(1:8),],
                                 iter = 4000)

posterior_summary(model.vel1.PLATE.ΔLLM1)
describe_posterior(model.vel1.PLATE.ΔLLM1)
rope(model.vel1.PLATE.ΔLLM1)
posterior = as_draws_df(model.vel1.PLATE.ΔLLM1)
1-mean(posterior$b_TreatmentΔLLM1   > 0) # 0.05

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.82, 1.18] |   100% | [-0.02, 0.02] |        0% | 1.000 | 4313.00
TreatmentΔLLM1 |   0.20 | [-0.06, 0.47] | 95.00% | [-0.02, 0.02] |     2.14% | 1.000 | 4145.00

# WT - LLM1-C ----

model.vel1.PLATE.LLM1.C     = brm(Expression ~ Treatment, data = vel1_expression_raw_PLATE[c(1:4,9:12),],
                                  iter = 4000)

posterior_summary(model.vel1.PLATE.LLM1.C)
describe_posterior(model.vel1.PLATE.LLM1.C)
rope(model.vel1.PLATE.LLM1.C)
posterior = as_draws_df(model.vel1.PLATE.LLM1.C)
1-mean(posterior$b_TreatmentWT   > 0) # 0.09225

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   0.92 | [ 0.83, 1.01] |   100% | [-0.01, 0.01] |        0% | 1.001 | 3842.00
TreatmentWT |   0.08 | [-0.05, 0.22] | 90.77% | [-0.01, 0.01] |     3.50% | 1.001 | 4255.00

# WT - oeLLM1 ----

model.vel1.PLATE.oeLLM1     = brm(Expression ~ Treatment, data = vel1_expression_raw_PLATE[c(1:4,13:16),],
                                  iter = 4000)

posterior_summary(model.vel1.PLATE.oeLLM1)
describe_posterior(model.vel1.PLATE.oeLLM1)
rope(model.vel1.PLATE.oeLLM1)
posterior = as_draws_df(model.vel1.PLATE.oeLLM1)
mean(posterior$b_TreatmentWT   > 0) # 0.171125

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   1.06 | [ 0.95, 1.16] |   100% | [-0.01, 0.01] |        0% | 1.000 | 4404.00
TreatmentWT |  -0.06 | [-0.20, 0.08] | 82.89% | [-0.01, 0.01] |     5.57% | 1.000 | 4020.00

# WT - ΔVEL2 ----

model.vel1.PLATE.ΔVEL2     = brm(Expression ~ Treatment, data = vel1_expression_raw_PLATE[c(1:4,17:19),],
                                 iter = 4000)

posterior_summary(model.vel1.PLATE.ΔVEL2)
describe_posterior(model.vel1.PLATE.ΔVEL2)
rope(model.vel1.PLATE.ΔVEL2)
posterior = as_draws_df(model.vel1.PLATE.ΔVEL2)
mean(posterior$b_TreatmentΔVEL2   > 0) # 0.00325

# Summary of Posterior Distribution 

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.91,  1.09] |   100% | [-0.02, 0.02] |        0% | 1.001 | 3917.00
TreatmentΔVEL2 |  -0.28 | [-0.43, -0.14] | 99.67% | [-0.02, 0.02] |        0% | 1.001 | 3290.00

# VEL2 ----

# SXM ----

VEL2_expression_raw_SXM  = vel_expression_raw %>% filter(Stage == "SXM" & Protein == "VEL2")
VEL2_expression_raw_SXM$Treatment    = as.factor(VEL2_expression_raw_SXM$Treatment)

# WT - ΔLLM1 ----

model.VEL2.SXM.ΔLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_SXM[c(1:6),],
                               iter = 4000)

posterior_summary(model.VEL2.SXM.ΔLLM1)
describe_posterior(model.VEL2.SXM.ΔLLM1)
rope(model.VEL2.SXM.ΔLLM1)
posterior = as_draws_df(model.VEL2.SXM.ΔLLM1)
1-mean(posterior$b_TreatmentΔLLM1   > 0) # 0.1955

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.67, 1.32] | 99.87% | [-0.02, 0.02] |        0% | 1.001 | 4222.00
TreatmentΔLLM1 |   0.15 | [-0.29, 0.62] | 80.45% | [-0.02, 0.02] |     5.05% | 1.001 | 4019.00

# WT - LLM1-C ----

model.VEL2.SXM.LLM1.C     = brm(Expression ~ Treatment, data = VEL2_expression_raw_SXM[c(1:3,7:9),],
                                iter = 4000)

posterior_summary(model.VEL2.SXM.LLM1.C)
describe_posterior(model.VEL2.SXM.LLM1.C)
rope(model.VEL2.SXM.LLM1.C)
posterior = as_draws_df(model.VEL2.SXM.LLM1.C)
mean(posterior$b_TreatmentWT   > 0) # 0.464875

# Summary of Posterior Distribution 

Parameter   |    Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept) |      1.01 | [ 0.80, 1.21] |   100% | [-0.01, 0.01] |        0% | 1.000 | 3809.00
TreatmentWT | -8.77e-03 | [-0.30, 0.28] | 53.51% | [-0.01, 0.01] |     7.42% | 1.001 | 3118.00

# WT - oeLLM1 ----

model.VEL2.SXM.oeLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_SXM[c(1:3,10:12),],
                                iter = 4000)

posterior_summary(model.VEL2.SXM.oeLLM1)
describe_posterior(model.VEL2.SXM.oeLLM1)
rope(model.VEL2.SXM.oeLLM1)
posterior = as_draws_df(model.VEL2.SXM.oeLLM1)
mean(posterior$b_TreatmentWT   > 0) # 0.218375

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   1.27 | [ 0.59, 1.90] | 99.38% | [-0.03, 0.03] |        0% | 1.000 | 3820.00
TreatmentWT |  -0.27 | [-1.20, 0.70] | 78.16% | [-0.03, 0.03] |     5.24% | 1.000 | 3611.00

# WT - ΔVEL1ΔLLM1 ----

model.VEL2.SXM.ΔVEL1ΔLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_SXM[c(1:3,13:15),],
                                    iter = 4000)

posterior_summary(model.VEL2.SXM.ΔVEL1ΔLLM1)
describe_posterior(model.VEL2.SXM.ΔVEL1ΔLLM1)
rope(model.VEL2.SXM.ΔVEL1ΔLLM1)
posterior = as_draws_df(model.VEL2.SXM.ΔVEL1ΔLLM1)
1-mean(posterior$b_TreatmentΔVEL1ΔLLM1   > 0) # 0.1505

# Summary of Posterior Distribution 

Parameter           | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------------
(Intercept)         |   0.99 | [ 0.08, 1.97] | 98.06% | [-0.05, 0.05] |        0% | 1.001 | 3689.00
TreatmentΔVEL1ΔLLM1 |   0.57 | [-0.81, 1.95] | 84.95% | [-0.05, 0.05] |     4.26% | 1.001 | 3175.00

# WT - ΔVEL1oeLLM1 ----

model.VEL2.SXM.ΔVEL1oeLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_SXM[c(1:3,16:18),],
                                     iter = 4000)

posterior_summary(model.VEL2.SXM.ΔVEL1oeLLM1)
describe_posterior(model.VEL2.SXM.ΔVEL1oeLLM1)
rope(model.VEL2.SXM.ΔVEL1oeLLM1)
posterior = as_draws_df(model.VEL2.SXM.ΔVEL1oeLLM1)
1-mean(posterior$b_TreatmentΔVEL1oeLLM1   > 0) # 0.3805

# Summary of Posterior Distribution 

Parameter            | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
(Intercept)          |   0.99 | [ 0.31, 1.67] | 99.16% | [-0.03, 0.03] |        0% | 1.001 | 4407.00
TreatmentΔVEL1oeLLM1 |   0.10 | [-0.82, 1.07] | 61.95% | [-0.03, 0.03] |     6.71% | 1.001 | 4294.00

# WT - ΔVEL1 ----

model.VEL2.SXM.ΔVEL1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_SXM[c(1:3,19:21),],
                               iter = 4000)

posterior_summary(model.VEL2.SXM.ΔVEL1)
describe_posterior(model.VEL2.SXM.ΔVEL1)
rope(model.VEL2.SXM.ΔVEL1)
posterior = as_draws_df(model.VEL2.SXM.ΔVEL1)
1-mean(posterior$b_TreatmentΔVEL1   > 0) # 0.074375

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.17, 1.88] | 98.45% | [-0.06, 0.06] |        0% | 1.001 | 4022.00
TreatmentΔVEL1 |   0.77 | [-0.47, 1.90] | 92.56% | [-0.06, 0.06] |     2.08% | 1.001 | 3452.00

# PDM ----

VEL2_expression_raw_PDM  = vel_expression_raw %>% filter(Stage == "PDM" & Protein == "VEL2")
VEL2_expression_raw_PDM$Treatment    = as.factor(VEL2_expression_raw_PDM$Treatment)

# WT - ΔLLM1 ----

model.VEL2.PDM.ΔLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PDM[c(1:8),],
                               iter = 4000)

posterior_summary(model.VEL2.PDM.ΔLLM1)
describe_posterior(model.VEL2.PDM.ΔLLM1)
rope(model.VEL2.PDM.ΔLLM1)
posterior = as_draws_df(model.VEL2.PDM.ΔLLM1)
1-mean(posterior$b_TreatmentΔLLM1   > 0) # 0.02225

# Summary of Posterior Distribution 

Parameter      | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [0.62, 1.40] | 99.95% | [-0.04, 0.04] |        0% | 1.000 | 4032.00
TreatmentΔLLM1 |   0.59 | [0.03, 1.13] | 97.78% | [-0.04, 0.04] |     0.13% | 1.001 | 3991.00

# WT - LLM1-C ----

model.VEL2.PDM.LLM1.C     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PDM[c(1:4,9:12),],
                                iter = 4000)

posterior_summary(model.VEL2.PDM.LLM1.C)
describe_posterior(model.VEL2.PDM.LLM1.C)
rope(model.VEL2.PDM.LLM1.C)
posterior = as_draws_df(model.VEL2.PDM.LLM1.C)
mean(posterior$b_TreatmentWT   > 0) # 0.008

# Summary of Posterior Distribution 

Parameter   | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
(Intercept) |   1.42 | [ 1.20,  1.64] |   100% | [-0.03, 0.03] |        0% | 1.000 | 4074.00
TreatmentWT |  -0.42 | [-0.72, -0.12] | 99.20% | [-0.03, 0.03] |        0% | 1.000 | 4421.00

# WT - oeLLM1 ----

model.VEL2.PDM.oeLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PDM[c(1:4,13:16),],
                                iter = 4000)

posterior_summary(model.VEL2.PDM.oeLLM1)
describe_posterior(model.VEL2.PDM.oeLLM1)
rope(model.VEL2.PDM.oeLLM1)
posterior = as_draws_df(model.VEL2.PDM.oeLLM1)
mean(posterior$b_TreatmentWT   > 0) # 0.00075

# Summary of Posterior Distribution 

Parameter   | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
(Intercept) |   2.26 | [ 1.90,  2.60] |   100% | [-0.07, 0.07] |        0% | 1.001 | 3965.00
TreatmentWT |  -1.26 | [-1.74, -0.75] | 99.92% | [-0.07, 0.07] |        0% | 1.001 | 3730.00

# WT - ΔVEL1ΔLLM1 ----

model.VEL2.PDM.ΔVEL1ΔLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PDM[c(1:4,17:19),],
                                    iter = 4000)

posterior_summary(model.VEL2.PDM.ΔVEL1ΔLLM1)
describe_posterior(model.VEL2.PDM.ΔVEL1ΔLLM1)
rope(model.VEL2.PDM.ΔVEL1ΔLLM1)
posterior = as_draws_df(model.VEL2.PDM.ΔVEL1ΔLLM1)
mean(posterior$b_TreatmentΔVEL1ΔLLM1   > 0) # 0.13825

# Summary of Posterior Distribution 

Parameter           | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------------
(Intercept)         |   1.00 | [ 0.56, 1.44] | 99.95% | [-0.03, 0.03] |        0% | 1.000 | 4067.00
TreatmentΔVEL1ΔLLM1 |  -0.31 | [-0.97, 0.37] | 86.17% | [-0.03, 0.03] |     4.39% | 1.002 | 3393.00

# WT - ΔVEL1oeLLM1 ----

model.VEL2.PDM.ΔVEL1oeLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PDM[c(1:4,20:22),],
                                     iter = 5000)

posterior_summary(model.VEL2.PDM.ΔVEL1oeLLM1)
describe_posterior(model.VEL2.PDM.ΔVEL1oeLLM1)
rope(model.VEL2.PDM.ΔVEL1oeLLM1)
posterior = as_draws_df(model.VEL2.PDM.ΔVEL1oeLLM1)
1-mean(posterior$b_TreatmentΔVEL1oeLLM1   > 0) # 0.322

# Summary of Posterior Distribution 

Parameter            | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
(Intercept)          |   1.00 | [ 0.05, 1.99] | 97.84% | [-0.06, 0.06] |     0.18% | 1.000 | 4106.00
TreatmentΔVEL1oeLLM1 |   0.29 | [-1.23, 1.79] | 67.80% | [-0.06, 0.06] |     7.66% | 1.000 | 3594.00

# WT - ΔVEL1 ----

model.VEL2.PDM.ΔVEL1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PDM[c(1:4,23:25),],
                               iter = 4000)

posterior_summary(model.VEL2.PDM.ΔVEL1)
describe_posterior(model.VEL2.PDM.ΔVEL1)
rope(model.VEL2.PDM.ΔVEL1)
posterior = as_draws_df(model.VEL2.PDM.ΔVEL1)
mean(posterior$b_TreatmentΔVEL1   > 0) # 0.00325

# Summary of Posterior Distribution 

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.87,  1.11] |   100% | [-0.02, 0.02] |        0% | 1.000 | 3828.00
TreatmentΔVEL1 |  -0.35 | [-0.52, -0.16] | 99.67% | [-0.02, 0.02] |        0% | 1.000 | 4390.00

# PLATE ----

VEL2_expression_raw_PLATE  = vel_expression_raw %>% filter(Stage == "PLATE" & Protein == "VEL2")
VEL2_expression_raw_PLATE$Treatment    = as.factor(VEL2_expression_raw_PLATE$Treatment)

# WT - ΔLLM1 ----

model.VEL2.PLATE.ΔLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PLATE[c(1:8),],
                                 iter = 4000)

posterior_summary(model.VEL2.PLATE.ΔLLM1)
describe_posterior(model.VEL2.PLATE.ΔLLM1)
rope(model.VEL2.PLATE.ΔLLM1)
posterior = as_draws_df(model.VEL2.PLATE.ΔLLM1)
1-mean(posterior$b_TreatmentΔLLM1   > 0) # 0.1145

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.45, 1.51] | 99.61% | [-0.04, 0.04] |        0% | 1.001 | 4005.00
TreatmentΔLLM1 |   0.40 | [-0.33, 1.18] | 88.55% | [-0.04, 0.04] |     4.37% | 1.000 | 4197.00

# WT - LLM1-C ----

model.VEL2.PLATE.LLM1.C     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PLATE[c(1:4,9:12),],
                                  iter = 4000)

posterior_summary(model.VEL2.PLATE.LLM1.C)
describe_posterior(model.VEL2.PLATE.LLM1.C)
rope(model.VEL2.PLATE.LLM1.C)
posterior = as_draws_df(model.VEL2.PLATE.LLM1.C)
mean(posterior$b_TreatmentWT   > 0) # 0.04425

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   1.05 | [ 1.01, 1.09] |   100% | [ 0.00, 0.00] |        0% | 1.000 | 4385.00
TreatmentWT |  -0.05 | [-0.11, 0.01] | 95.58% | [ 0.00, 0.00] |     1.87% | 1.000 | 3971.00

# WT - oeLLM1 ----

model.VEL2.PLATE.oeLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PLATE[c(1:4,13:16),],
                                  iter = 4000)

posterior_summary(model.VEL2.PLATE.oeLLM1)
describe_posterior(model.VEL2.PLATE.oeLLM1)
rope(model.VEL2.PLATE.oeLLM1)
posterior = as_draws_df(model.VEL2.PLATE.oeLLM1)
1-mean(posterior$b_TreatmentWT   > 0) # 0.436125

# Summary of Posterior Distribution 

Parameter   | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------
(Intercept) |   0.98 | [ 0.76, 1.20] | 99.96% | [-0.01, 0.01] |        0% | 1.001 | 3795.00
TreatmentWT |   0.02 | [-0.29, 0.33] | 56.39% | [-0.01, 0.01] |     9.12% | 1.000 | 4751.00

# WT - ΔVEL1ΔLLM1 ----

model.VEL2.PLATE.ΔVEL1ΔLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PLATE[c(1:4,17:19),],
                                      iter = 4000)

posterior_summary(model.VEL2.PLATE.ΔVEL1ΔLLM1)
describe_posterior(model.VEL2.PLATE.ΔVEL1ΔLLM1)
rope(model.VEL2.PLATE.ΔVEL1ΔLLM1)
posterior = as_draws_df(model.VEL2.PLATE.ΔVEL1ΔLLM1)
1-mean(posterior$b_TreatmentΔVEL1ΔLLM1   > 0) # 0.012

# Summary of Posterior Distribution 

Parameter           | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)         |   1.00 | [0.85, 1.16] |   100% | [-0.02, 0.02] |        0% | 1.000 | 3385.00
TreatmentΔVEL1ΔLLM1 |   0.32 | [0.08, 0.55] | 98.80% | [-0.02, 0.02] |        0% | 1.000 | 3301.00

# WT - ΔVEL1oeLLM1 ----

model.VEL2.PLATE.ΔVEL1oeLLM1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PLATE[c(1:4,20:22),],
                                       iter = 4000)

posterior_summary(model.VEL2.PLATE.ΔVEL1oeLLM1)
describe_posterior(model.VEL2.PLATE.ΔVEL1oeLLM1)
rope(model.VEL2.PLATE.ΔVEL1oeLLM1)
posterior = as_draws_df(model.VEL2.PLATE.ΔVEL1oeLLM1)
mean(posterior$b_TreatmentΔVEL1oeLLM1   > 0) # 0.4315

# Summary of Posterior Distribution 

Parameter            | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------------
(Intercept)          |   1.00 | [ 0.50, 1.54] | 99.78% | [-0.03, 0.03] |        0% | 1.000 | 4071.00
TreatmentΔVEL1oeLLM1 |  -0.05 | [-0.81, 0.68] | 56.85% | [-0.03, 0.03] |     8.21% | 1.000 | 3688.00

# WT - ΔVEL1 ----

model.VEL2.PLATE.ΔVEL1     = brm(Expression ~ Treatment, data = VEL2_expression_raw_PLATE[c(1:4,23:25),],
                                 iter = 6000)

posterior_summary(model.VEL2.PLATE.ΔVEL1)
describe_posterior(model.VEL2.PLATE.ΔVEL1)
rope(model.VEL2.PLATE.ΔVEL1)
posterior = as_draws_df(model.VEL2.PLATE.ΔVEL1)
1-mean(posterior$b_TreatmentΔVEL1   > 0) # 0.1448333

# Summary of Posterior Distribution 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   0.99 | [ 0.11, 1.84] | 98.23% | [-0.06, 0.06] |        0% | 1.000 | 5891.00
TreatmentΔVEL1 |   0.60 | [-0.72, 1.94] | 85.52% | [-0.06, 0.06] |     4.60% | 1.000 | 5795.00

# PCR Light interaction ----

# Call data
light_expression     = read_excel("datasets/WC and FRQ - my genes.xlsx",
                                  sheet = "statistics")
light_expression_raw = light_expression

# Transform to Log2
light_expression = light_expression %>%  mutate(across(c(Expression), 
                                                       function(x) log2(x)))

# Delete WT
light_expression = light_expression %>% filter(Treatment != "WT")

# SXM ----

light_expression.SXM   = light_expression %>% filter(Stage == "SXM")

light_expression.SXM_plot = light_expression.SXM %>% group_by(Genes,Treatment,Plot) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_light_expression.SXM = ggplot(light_expression.SXM_plot, aes(x=Genes, y=avg, fill=as.factor(Plot))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#fec44f", "#d95f0e"), 
                    labels = c("ΔFRQ", "ΔWC1")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_light_expression.SXM

pdf("Figures/Figure_light_expression.SXM.pdf",
    width=6,height=6*3/5)
print(Figure_light_expression.SXM)
dev.off()

# PDM ----

light_expression.PDM   = light_expression %>% filter(Stage == "PDM")

light_expression.PDM_plot = light_expression.PDM %>% group_by(Genes,Treatment,Plot) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_light_expression.PDM = ggplot(light_expression.PDM_plot, aes(x=Genes, y=avg, fill=as.factor(Plot))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#fec44f", "#d95f0e"), 
                    labels = c("ΔFRQ", "ΔWC1")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_light_expression.PDM

pdf("Figures/Figure_light_expression.PDM.pdf",
    width=6,height=6*3/5)
print(Figure_light_expression.PDM)
dev.off()

# PLATE ----

light_expression.PLATE   = light_expression %>% filter(Stage == "PLATE")

light_expression.PLATE_plot = light_expression.PLATE %>% group_by(Genes,Treatment,Plot) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_light_expression.PLATE = ggplot(light_expression.PLATE_plot, aes(x=Genes, y=avg, fill=as.factor(Plot))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#fec44f", "#d95f0e"), 
                    labels = c("ΔFRQ", "ΔWC1")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_light_expression.PLATE

pdf("Figures/Figure_light_expression.PLATE.pdf",
    width=6,height=6*3/5)
print(Figure_light_expression.PLATE)
dev.off()

# Bayesian Methods for Group Comparison ----

# SXM ----

# LLM1 ----

light_expression_raw_SXM  = light_expression_raw %>% filter(Stage == "SXM" & Genes == "LLM1")
light_expression_raw_SXM$Treatment    = as.factor(light_expression_raw_SXM$Treatment)

# WT - ΔFRQ ----

model.light.SXM.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_SXM[1:6,],
                               iter = 6000)

posterior_summary(model.light.SXM.ΔFRQ)
describe_posterior(model.light.SXM.ΔFRQ)
rope(model.light.SXM.ΔFRQ)
posterior = as_draws_df(model.light.SXM.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.4145833

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.66, 1.33] | 99.92% | [-0.01, 0.01] |        0% | 1.000 | 5147.00
TreatmentΔFRQ |  -0.04 | [-0.52, 0.45] | 58.54% | [-0.01, 0.01] |     6.91% | 1.000 | 4883.00

# WT - ΔWC1 ----

model.light.SXM.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_SXM[c(1:3,7:9),],
                               iter = 6000)

posterior_summary(model.light.SXM.ΔWC1)
describe_posterior(model.light.SXM.ΔWC1)
rope(model.light.SXM.ΔWC1)
posterior = as_draws_df(model.light.SXM.ΔWC1)
mean(posterior$b_TreatmentΔWC1 > 0) # 0.1526667

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.83, 1.17] |   100% | [-0.01, 0.01] |        0% | 1.002 | 4827.00
TreatmentΔWC1 |  -0.10 | [-0.34, 0.14] | 84.73% | [-0.01, 0.01] |     3.96% | 1.001 | 5050.00

# AML1 ----

light_expression_raw_SXM  = light_expression_raw %>% filter(Stage == "SXM" & Genes == "AML1")
light_expression_raw_SXM$Treatment    = as.factor(light_expression_raw_SXM$Treatment)

# WT - ΔFRQ ----

model.light.SXM.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_SXM[c(1:6),],
                               iter = 6000)

posterior_summary(model.light.SXM.ΔFRQ)
describe_posterior(model.light.SXM.ΔFRQ)
rope(model.light.SXM.ΔFRQ)
posterior = as_draws_df(model.light.SXM.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.005583333

Summary of Posterior Distribution 

Parameter     | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.78,  1.21] |   100% | [-0.03, 0.03] |        0% | 1.001 | 5566.00
TreatmentΔFRQ |  -0.54 | [-0.85, -0.23] | 99.44% | [-0.03, 0.03] |        0% | 1.002 | 4878.00

# WT - ΔWC1 ----

model.light.SXM.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_SXM[c(1:3,7:9),],
                               iter = 6000)

posterior_summary(model.light.SXM.ΔWC1)
describe_posterior(model.light.SXM.ΔWC1)
rope(model.light.SXM.ΔWC1)
posterior = as_draws_df(model.light.SXM.ΔWC1)
mean(posterior$b_TreatmentΔWC1 > 0) # 0.001416667

Summary of Posterior Distribution 

Parameter     | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.90,  1.10] |   100% | [-0.02, 0.02] |        0% | 1.000 | 5331.00
TreatmentΔWC1 |  -0.37 | [-0.51, -0.23] | 99.86% | [-0.02, 0.02] |        0% | 1.000 | 6318.00

# NML1 ----

light_expression_raw_SXM  = light_expression_raw %>% filter(Stage == "SXM" & Genes == "NML1")
light_expression_raw_SXM$Treatment    = as.factor(light_expression_raw_SXM$Treatment)

# WT - ΔFRQ ----

model.light.SXM.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_SXM[c(1:6),],
                               iter = 6000)

posterior_summary(model.light.SXM.ΔFRQ)
describe_posterior(model.light.SXM.ΔFRQ)
rope(model.light.SXM.ΔFRQ)
posterior = as_draws_df(model.light.SXM.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.2806667

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.59, 1.43] | 99.88% | [-0.02, 0.02] |        0% | 1.002 | 4810.00
TreatmentΔFRQ |  -0.13 | [-0.70, 0.46] | 71.93% | [-0.02, 0.02] |     5.75% | 1.001 | 4457.00

# WT - ΔWC1 ----

model.light.SXM.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_SXM[c(1:3,7:9),],
                               iter = 6000)

posterior_summary(model.light.SXM.ΔWC1)
describe_posterior(model.light.SXM.ΔWC1)
rope(model.light.SXM.ΔWC1)
posterior = as_draws_df(model.light.SXM.ΔWC1)
1-mean(posterior$b_TreatmentΔWC1 > 0) # 0.2058333

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.69, 1.31] | 99.92% | [-0.02, 0.02] |        0% | 1.000 | 4802.00
TreatmentΔWC1 |   0.13 | [-0.31, 0.59] | 79.42% | [-0.02, 0.02] |     5.13% | 1.000 | 4459.00

# PDM ----

# LLM1 ----

light_expression_raw_PDM  = light_expression_raw %>% filter(Stage == "PDM" & Genes == "LLM1")
light_expression_raw_PDM$Treatment    = as.factor(light_expression_raw_PDM$Treatment)

# WT - ΔFRQ ----

model.light.PDM.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_PDM[c(1:6),],
                               iter = 6000)

posterior_summary(model.light.PDM.ΔFRQ)
describe_posterior(model.light.PDM.ΔFRQ)
rope(model.light.PDM.ΔFRQ)
posterior = as_draws_df(model.light.PDM.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.1805833

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.73, 1.26] | 99.98% | [-0.01, 0.01] |        0% | 1.000 | 4941.00
TreatmentΔFRQ |  -0.13 | [-0.52, 0.26] | 81.94% | [-0.01, 0.01] |     4.35% | 1.000 | 5382.00

# WT - ΔWC1 ----

model.light.PDM.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_PDM[c(1:3,7:9),],
                               iter = 6000)

posterior_summary(model.light.PDM.ΔWC1)
describe_posterior(model.light.PDM.ΔWC1)
rope(model.light.PDM.ΔWC1)
posterior = as_draws_df(model.light.PDM.ΔWC1)
mean(posterior$b_TreatmentΔWC1 > 0) # 0.2265

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.66, 1.34] | 99.92% | [-0.02, 0.02] |        0% | 1.000 | 5344.00
TreatmentΔWC1 |  -0.13 | [-0.62, 0.38] | 77.35% | [-0.02, 0.02] |     4.82% | 1.001 | 5251.00

# AML1 ----

light_expression_raw_PDM  = light_expression_raw %>% filter(Stage == "PDM" & Genes == "AML1")
light_expression_raw_PDM$Treatment    = as.factor(light_expression_raw_PDM$Treatment)

# WT - ΔFRQ ----

model.light.PDM.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_PDM[c(1:6),],
                               iter = 6000)

posterior_summary(model.light.PDM.ΔFRQ)
describe_posterior(model.light.PDM.ΔFRQ)
rope(model.light.PDM.ΔFRQ)
posterior = as_draws_df(model.light.PDM.ΔFRQ)
1-mean(posterior$b_TreatmentΔFRQ > 0) # 0.05175

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   0.99 | [-0.49, 2.44] | 92.95% | [-0.12, 0.12] |     3.60% | 1.000 | 6103.00
TreatmentΔFRQ |   1.61 | [-0.50, 3.67] | 94.83% | [-0.12, 0.12] |     1.86% | 1.000 | 6002.00

# WT - ΔWC1 ----

model.light.PDM.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_PDM[c(1:3,7:9),],
                               iter = 6000)

posterior_summary(model.light.PDM.ΔWC1)
describe_posterior(model.light.PDM.ΔWC1)
rope(model.light.PDM.ΔWC1)
posterior = as_draws_df(model.light.PDM.ΔWC1)
1-mean(posterior$b_TreatmentΔWC1 > 0) # 0.005416667

Summary of Posterior Distribution 

Parameter     | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [0.54, 1.43] | 99.74% | [-0.07, 0.07] |        0% | 1.000 | 4896.00
TreatmentΔWC1 |   1.15 | [0.51, 1.81] | 99.46% | [-0.07, 0.07] |        0% | 1.001 | 4834.00

# NML1 ----

light_expression_raw_PDM  = light_expression_raw %>% filter(Stage == "PDM" & Genes == "NML1")
light_expression_raw_PDM$Treatment    = as.factor(light_expression_raw_PDM$Treatment)

# WT - ΔFRQ ----

model.light.PDM.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_PDM[c(1:6),],
                               iter = 6000)

posterior_summary(model.light.PDM.ΔFRQ)
describe_posterior(model.light.PDM.ΔFRQ)
rope(model.light.PDM.ΔFRQ)
posterior = as_draws_df(model.light.PDM.ΔFRQ)
1-mean(posterior$b_TreatmentΔFRQ > 0) # 8.333333e-05

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   0.99 | [-0.29, 2.29] | 95.32% | [-0.42, 0.42] |    11.12% | 1.000 | 5163.00
TreatmentΔFRQ |   7.69 | [ 6.02, 9.61] | 99.99% | [-0.42, 0.42] |        0% | 1.000 | 6059.00

# WT - ΔWC1 ----

model.light.PDM.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_PDM[c(1:3,7:9),],
                               iter = 6000)

posterior_summary(model.light.PDM.ΔWC1)
describe_posterior(model.light.PDM.ΔWC1)
rope(model.light.PDM.ΔWC1)
posterior = as_draws_df(model.light.PDM.ΔWC1)
1-mean(posterior$b_TreatmentΔWC1 > 0) # 0.01016667

Summary of Posterior Distribution 

Parameter     | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)   |   0.91 | [-2.87,  4.64] | 71.83% | [-0.47, 0.47] |    19.62% | 1.000 | 5887.00
TreatmentΔWC1 |   7.57 | [ 1.84, 13.04] | 98.98% | [-0.47, 0.47] |        0% | 1.001 | 5027.00

# PLATE ----

# LLM1 ----

light_expression_raw_PLATE  = light_expression_raw %>% filter(Stage == "PLATE" & Genes == "LLM1")
light_expression_raw_PLATE$Treatment    = as.factor(light_expression_raw_PLATE$Treatment)

# WT - ΔFRQ ----

model.light.PLATE.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_PLATE[c(1:6),],
                                 iter = 6000)

posterior_summary(model.light.PLATE.ΔFRQ)
describe_posterior(model.light.PLATE.ΔFRQ)
rope(model.light.PLATE.ΔFRQ)
posterior = as_draws_df(model.light.PLATE.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.006083333

Summary of Posterior Distribution 

Parameter     | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.82,  1.17] | 99.97% | [-0.03, 0.03] |        0% | 1.001 | 3962.00
TreatmentΔFRQ |  -0.44 | [-0.67, -0.21] | 99.39% | [-0.03, 0.03] |        0% | 1.000 | 4112.00

# WT - ΔWC1 ----

model.light.PLATE.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_PLATE[c(1:3,7:9),],
                                 iter = 6000)

posterior_summary(model.light.PLATE.ΔWC1)
describe_posterior(model.light.PLATE.ΔWC1)
rope(model.light.PLATE.ΔWC1)
posterior = as_draws_df(model.light.PLATE.ΔWC1)
mean(posterior$b_TreatmentΔWC1 > 0) # 0.027

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.81, 1.18] |   100% | [-0.02, 0.02] |        0% | 1.000 | 5339.00
TreatmentΔWC1 |  -0.26 | [-0.51, 0.01] | 97.30% | [-0.02, 0.02] |     0.70% | 1.000 | 4865.00

# AML1 ----

light_expression_raw_PLATE  = light_expression_raw %>% filter(Stage == "PLATE" & Genes == "AML1")
light_expression_raw_PLATE$Treatment    = as.factor(light_expression_raw_PLATE$Treatment)

# WT - ΔFRQ ----

model.light.PLATE.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_PLATE[c(1:6),],
                                 iter = 6000)

posterior_summary(model.light.PLATE.ΔFRQ)
describe_posterior(model.light.PLATE.ΔFRQ)
rope(model.light.PLATE.ΔFRQ)
posterior = as_draws_df(model.light.PLATE.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.3410833

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.84, 1.16] |   100% | [-0.01, 0.01] |        0% | 1.000 | 5400.00
TreatmentΔFRQ |  -0.03 | [-0.27, 0.18] | 65.89% | [-0.01, 0.01] |     6.44% | 1.001 | 4800.00

# WT - ΔWC1 ----

model.light.PLATE.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_PLATE[c(1:3,7:9),],
                                 iter = 6000)

posterior_summary(model.light.PLATE.ΔWC1)
describe_posterior(model.light.PLATE.ΔWC1)
rope(model.light.PLATE.ΔWC1)
posterior = as_draws_df(model.light.PLATE.ΔWC1)
1-mean(posterior$b_TreatmentΔWC1 > 0) # 0.0165

Summary of Posterior Distribution 

Parameter     | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [0.79, 1.23] | 99.98% | [-0.02, 0.02] |        0% | 1.001 | 5110.00
TreatmentΔWC1 |   0.39 | [0.07, 0.69] | 98.35% | [-0.02, 0.02] |        0% | 1.001 | 4799.00

# NML1 ----

light_expression_raw_PLATE  = light_expression_raw %>% filter(Stage == "PLATE" & Genes == "NML1")
light_expression_raw_PLATE$Treatment    = as.factor(light_expression_raw_PLATE$Treatment)

# WT - ΔFRQ ----

model.light.PLATE.ΔFRQ     = brm(Expression ~ Treatment, data = light_expression_raw_PLATE[c(1:6),],
                                 iter = 6000)

posterior_summary(model.light.PLATE.ΔFRQ)
describe_posterior(model.light.PLATE.ΔFRQ)
rope(model.light.PLATE.ΔFRQ)
posterior = as_draws_df(model.light.PLATE.ΔFRQ)
mean(posterior$b_TreatmentΔFRQ > 0) # 0.47525

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.53, 1.43] | 99.71% | [-0.02, 0.02] |        0% | 1.001 | 3421.00
TreatmentΔFRQ |  -0.01 | [-0.63, 0.64] | 52.48% | [-0.02, 0.02] |     6.86% | 1.000 | 3541.00

# WT - ΔWC1 ----

model.light.PLATE.ΔWC1     = brm(Expression ~ Treatment, data = light_expression_raw_PLATE[c(1:3,7:9),],
                                 iter = 6000)

posterior_summary(model.light.PLATE.ΔWC1)
describe_posterior(model.light.PLATE.ΔWC1)
rope(model.light.PLATE.ΔWC1)
posterior = as_draws_df(model.light.PLATE.ΔWC1)
1-mean(posterior$b_TreatmentΔWC1 > 0) # 0.1115833

Summary of Posterior Distribution 

Parameter     | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)   |   1.00 | [ 0.47, 1.52] | 99.56% | [-0.03, 0.03] |        0% | 1.001 | 5435.00
TreatmentΔWC1 |   0.35 | [-0.40, 1.08] | 88.84% | [-0.03, 0.03] |     3.20% | 1.000 | 5257.00

# PCR graphs other genes ----

# Call data
PCR_others = read_excel("datasets/All data combined.xlsx", sheet = "Thesis_2")
PCR_raw    = PCR_others

# Transform to Log2
PCR_others = PCR_others %>%  mutate(across(c(Expression), function(x) log2(x)))

# Delete WT
PCR_others = PCR_others %>% filter(Treatment != "WT")

# Plot 

# SXM ----
PCR_others_SXM   = PCR_others %>% filter(Stage == "SXM")
PCR_others_SXM.1 = PCR_others_SXM %>% group_by(Gen,Treatment,Order) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_others_SXM = ggplot(PCR_others_SXM.1, aes(x=Gen, y=avg, fill=as.factor(Order))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#d0d1e6", "#045a8d", "#fee391", "#e34a33",
                                          "#74c476"), labels = c("ΔVel1", "oeVel1", "ΔVel2", 
                                                                 "oeVel2", "ΔVel2_IDD")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_others_SXM

pdf("Figures/Figure_others_SXM.pdf",
    width=12,height=12*3/5)
print(Figure_others_SXM)
dev.off()

# PDM ----
PCR_others_PDM   = PCR_others %>% filter(Stage == "PDM")
PCR_others_PDM.1 = PCR_others_PDM %>% group_by(Gen,Treatment,Order) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_others_PDM = ggplot(PCR_others_PDM.1, aes(x=Gen, y=avg, fill=as.factor(Order))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#d0d1e6", "#045a8d", "#fee391", "#e34a33",
                                          "#74c476"), labels = c("ΔVel1", "oeVel1", "ΔVel2", 
                                                                 "oeVel2", "ΔVel2_IDD")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_others_PDM

pdf("Figures/Figure_others_PDM.pdf",
    width=12,height=12*3/5)
print(Figure_others_PDM)
dev.off()

# PLATE ----
PCR_others_PLATE   = PCR_others %>% filter(Stage == "PLATE")
PCR_others_PLATE.1 = PCR_others_PLATE %>% group_by(Gen,Treatment,Order) %>% 
  summarize(avg = mean(Expression), n = n(), 
            sd = sd(Expression), se = sd/sqrt(n))

Figure_others_PLATE = ggplot(PCR_others_PLATE.1, aes(x=Gen, y=avg, fill=as.factor(Order))) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(name = "", values = c("#d0d1e6", "#045a8d", "#fee391", "#e34a33",
                                          "#74c476"), labels = c("ΔVel1", "oeVel1", "ΔVel2", 
                                                                 "oeVel2", "ΔVel2_IDD")) + 
  xlab("") + ylab("Expression") + theme_classic()
Figure_others_PLATE

pdf("Figures/Figure_others_PLATE.pdf",
    width=12,height=12*3/5)
print(Figure_others_PLATE)
dev.off()

# Bayesian Methods for Group Comparison ----

# SXM ----
PCR_raw_SXM   = PCR_raw %>% filter(Stage == "SXM")
PCR_raw_SXM$Treatment = as.factor(PCR_raw_SXM$Treatment)

# MET1 - ΔVel1
model.MET1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_SXM[1:6,],
                           iter = 6000)
summary(model.MET1.ΔVel1)
posterior_summary(model.MET1.ΔVel1)
describe_posterior(model.MET1.ΔVel1)
rope(model.MET1.ΔVel1)
posterior = as_draws_df(model.MET1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.01141667

Parameter      | Median |          95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
------------------------------------------------------------------------------------------------
(Intercept)    |  -3.62 | [-28.24, 17.91] | 63.92% | [-2.99, 2.99] |    23.10% | 1.000 | 7174.00
TreatmentΔVel1 |  43.77 | [  8.06, 79.04] | 98.86% | [-2.99, 2.99] |        0% | 1.000 | 7293.00

# MET1 - ΔVel2
model.MET1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(1:3,7:9),],
                           iter = 6000)
summary(model.MET1.ΔVel2)
posterior_summary(model.MET1.ΔVel2)
describe_posterior(model.MET1.ΔVel2)
rope(model.MET1.ΔVel2)
posterior = as_draws_df(model.MET1.ΔVel2)
1- mean(posterior$b_TreatmentΔVel2 > 0) # 

Parameter      | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [0.37, 1.66] | 99.24% | [-0.12, 0.12] |        0% | 1.000 | 5692.00
TreatmentΔVel2 |   2.14 | [1.17, 3.05] | 99.75% | [-0.12, 0.12] |        0% | 1.001 | 5539.00

# MET1 - ΔVel2_IDD
model.MET1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(1:3,10:12),],
                               iter = 6000)
summary(model.MET1.ΔVel2_IDD)
posterior_summary(model.MET1.ΔVel2_IDD)
describe_posterior(model.MET1.ΔVel2_IDD)
rope(model.MET1.ΔVel2_IDD)
posterior = as_draws_df(model.MET1.ΔVel2_IDD)
1- mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.1205833

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [-0.05, 1.97] | 97.22% | [-0.06, 0.06] |     0.84% | 1.001 | 5410.00
TreatmentΔVel2_IDD |   0.72 | [-0.70, 2.18] | 87.94% | [-0.06, 0.06] |     3.64% | 1.001 | 5715.00

# CTZ1 - ΔVel1
model.CTZ1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(13:18),],
                           iter = 6000)
summary(model.CTZ1.ΔVel1)
posterior_summary(model.CTZ1.ΔVel1)
describe_posterior(model.CTZ1.ΔVel1)
rope(model.CTZ1.ΔVel1)
posterior = as_draws_df(model.CTZ1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.00

Parameter      | Median |         95% CI |   pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.95,  1.06] | 100% | [-0.04, 0.04] |        0% | 1.000 | 6744.00
TreatmentΔVel1 |  -0.80 | [-0.88, -0.73] | 100% | [-0.04, 0.04] |        0% | 1.000 | 6487.00

# CTZ1 - ΔVel2
model.CTZ1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(13:15,19:21),],
                           iter = 6000)
summary(model.CTZ1.ΔVel2)
posterior_summary(model.CTZ1.ΔVel2)
describe_posterior(model.CTZ1.ΔVel2)
rope(model.CTZ1.ΔVel2)
posterior = as_draws_df(model.CTZ1.ΔVel2)
1- mean(posterior$b_TreatmentΔVel2 > 0) # 0.11125

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [-0.01, 2.04] | 97.42% | [-0.07, 0.07] |     0.63% | 1.000 | 6261.00
TreatmentΔVel2 |   0.79 | [-0.79, 2.22] | 88.88% | [-0.07, 0.07] |     3.48% | 1.000 | 5756.00

# CTZ1 - ΔVel2_IDD
model.CTZ1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(13:15,22:24),],
                               iter = 6000)
summary(model.CTZ1.ΔVel2_IDD)
posterior_summary(model.CTZ1.ΔVel2_IDD)
describe_posterior(model.CTZ1.ΔVel2_IDD)
rope(model.CTZ1.ΔVel2_IDD)
posterior = as_draws_df(model.CTZ1.ΔVel2_IDD)
1- mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.1949167

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   0.99 | [-0.13, 1.97] | 96.43% | [-0.06, 0.06] |     1.00% | 1.000 | 5659.00
TreatmentΔVel2_IDD |   0.50 | [-0.94, 2.05] | 80.51% | [-0.06, 0.06] |     5.00% | 1.000 | 4194.00

# YAL1 - ΔVel1
model.YAL1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(25:30),],
                           iter = 6000)
summary(model.YAL1.ΔVel1)
posterior_summary(model.YAL1.ΔVel1)
describe_posterior(model.YAL1.ΔVel1)
rope(model.YAL1.ΔVel1)
posterior = as_draws_df(model.YAL1.ΔVel1)
mean(posterior$b_TreatmentΔVel1 > 0) # 0.04258333

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.71, 1.27] | 99.85% | [-0.02, 0.02] |        0% | 1.000 | 3811.00
TreatmentΔVel1 |  -0.32 | [-0.70, 0.11] | 95.74% | [-0.02, 0.02] |     1.14% | 1.000 | 4797.00

# YAL1 - ΔVel2
model.YAL1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(25:27,31:33),],
                           iter = 8000)
summary(model.YAL1.ΔVel2)
posterior_summary(model.YAL1.ΔVel2)
describe_posterior(model.YAL1.ΔVel2)
rope(model.YAL1.ΔVel2)
posterior = as_draws_df(model.YAL1.ΔVel2)
mean(posterior$b_TreatmentΔVel2 > 0) # 0.010125

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.81,  1.20] | 99.98% | [-0.02, 0.02] |        0% | 1.001 | 6052.00
TreatmentΔVel2 |  -0.40 | [-0.69, -0.12] | 98.99% | [-0.02, 0.02] |        0% | 1.001 | 4939.00

# YAL1 - ΔVel2_IDD
model.YAL1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_SXM[c(25:27,34:36),],
                               iter = 6000)
summary(model.YAL1.ΔVel2_IDD)
posterior_summary(model.YAL1.ΔVel2_IDD)
describe_posterior(model.YAL1.ΔVel2_IDD)
rope(model.YAL1.ΔVel2_IDD)
posterior = as_draws_df(model.YAL1.ΔVel2_IDD)
mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.08741667

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [ 0.64, 1.34] | 99.96% | [-0.02, 0.02] |        0% | 1.001 | 5334.00
TreatmentΔVel2_IDD |  -0.30 | [-0.80, 0.20] | 91.26% | [-0.02, 0.02] |     2.65% | 1.001 | 5824.00

# PDM ----
PCR_raw_PDM   = PCR_raw %>% filter(Stage == "PDM")
PCR_raw_PDM$Treatment = as.factor(PCR_raw_PDM$Treatment)

# MET1 - ΔVel1
model.MET1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PDM[1:6,],
                           iter = 6000)
summary(model.MET1.ΔVel1)
posterior_summary(model.MET1.ΔVel1)
describe_posterior(model.MET1.ΔVel1)
rope(model.MET1.ΔVel1)
posterior = as_draws_df(model.MET1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.01691667

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   0.97 | [-0.58, 2.44] | 91.32% | [-0.16, 0.16] |     5.78% | 1.000 | 6330.00
TreatmentΔVel1 |   2.46 | [ 0.27, 4.64] | 98.31% | [-0.16, 0.16] |        0% | 1.000 | 6323.00

# MET1 - ΔVel2
model.MET1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(1:3,7:9),],
                           iter = 6000)
summary(model.MET1.ΔVel2)
posterior_summary(model.MET1.ΔVel2)
describe_posterior(model.MET1.ΔVel2)
rope(model.MET1.ΔVel2)
posterior = as_draws_df(model.MET1.ΔVel2)
1- mean(posterior$b_TreatmentΔVel2 > 0) # 0.1528333 

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.43, 1.58] | 99.57% | [-0.03, 0.03] |        0% | 1.001 | 4646.00
TreatmentΔVel2 |   0.33 | [-0.51, 1.15] | 84.72% | [-0.03, 0.03] |     4.07% | 1.001 | 4680.00

# MET1 - ΔVel2_IDD
model.MET1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(1:3,10:12),],
                               iter = 6000)
summary(model.MET1.ΔVel2_IDD)
posterior_summary(model.MET1.ΔVel2_IDD)
describe_posterior(model.MET1.ΔVel2_IDD)
rope(model.MET1.ΔVel2_IDD)
posterior = as_draws_df(model.MET1.ΔVel2_IDD)
1- mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.001833333

Parameter          | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [0.70, 1.29] | 99.98% | [-0.06, 0.06] |        0% | 1.000 | 5297.00
TreatmentΔVel2_IDD |   1.09 | [0.66, 1.53] | 99.82% | [-0.06, 0.06] |        0% | 1.001 | 4646.00

# CTZ1 - ΔVel1
model.CTZ1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(13:18),],
                           iter = 6000)
summary(model.CTZ1.ΔVel1)
posterior_summary(model.CTZ1.ΔVel1)
describe_posterior(model.CTZ1.ΔVel1)
rope(model.CTZ1.ΔVel1)
posterior = as_draws_df(model.CTZ1.ΔVel1)
mean(posterior$b_TreatmentΔVel1 > 0) # 0.0006666667

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.88,  1.12] |   100% | [-0.04, 0.04] |        0% | 1.001 | 5184.00
TreatmentΔVel1 |  -0.66 | [-0.83, -0.49] | 99.93% | [-0.04, 0.04] |        0% | 1.001 | 4989.00

# CTZ1 - ΔVel2
model.CTZ1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(13:15,19:21),],
                           iter = 6000)
summary(model.CTZ1.ΔVel2)
posterior_summary(model.CTZ1.ΔVel2)
describe_posterior(model.CTZ1.ΔVel2)
rope(model.CTZ1.ΔVel2)
posterior = as_draws_df(model.CTZ1.ΔVel2)
1- mean(posterior$b_TreatmentΔVel2 > 0) # 0.00001

Parameter      | Median |         95% CI |   pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.96,  1.03] | 100% | [-0.04, 0.04] |        0% | 1.000 | 5756.00
TreatmentΔVel2 |  -0.66 | [-0.71, -0.61] | 100% | [-0.04, 0.04] |        0% | 1.001 | 5571.00

# CTZ1 - ΔVel2_IDD
model.CTZ1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(13:15,22:24),],
                               iter = 6000)
summary(model.CTZ1.ΔVel2_IDD)
posterior_summary(model.CTZ1.ΔVel2_IDD)
describe_posterior(model.CTZ1.ΔVel2_IDD)
rope(model.CTZ1.ΔVel2_IDD)
posterior = as_draws_df(model.CTZ1.ΔVel2_IDD)
mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.1081667

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [ 0.29, 1.69] | 99.02% | [-0.04, 0.04] |        0% | 1.000 | 4425.00
TreatmentΔVel2_IDD |  -0.48 | [-1.46, 0.50] | 89.18% | [-0.04, 0.04] |     3.01% | 1.000 | 4667.00

# YAL1 - ΔVel1
model.YAL1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(25:30),],
                           iter = 6000)
summary(model.YAL1.ΔVel1)
posterior_summary(model.YAL1.ΔVel1)
describe_posterior(model.YAL1.ΔVel1)
rope(model.YAL1.ΔVel1)
posterior = as_draws_df(model.YAL1.ΔVel1)
1-mean(posterior$b_TreatmentΔVel1 > 0) # 0.088

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.65, 1.36] | 99.86% | [-0.02, 0.02] |        0% | 1.000 | 5497.00
TreatmentΔVel1 |   0.29 | [-0.22, 0.80] | 91.20% | [-0.02, 0.02] |     2.34% | 1.000 | 4506.00

# YAL1 - ΔVel2
model.YAL1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(25:27,31:33),],
                           iter = 6000)
summary(model.YAL1.ΔVel2)
posterior_summary(model.YAL1.ΔVel2)
describe_posterior(model.YAL1.ΔVel2)
rope(model.YAL1.ΔVel2)
posterior = as_draws_df(model.YAL1.ΔVel2)
1-mean(posterior$b_TreatmentΔVel2 > 0) # 0.01783333

Parameter      | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [0.39, 1.60] | 99.43% | [-0.06, 0.06] |        0% | 1.000 | 5497.00
TreatmentΔVel2 |   0.94 | [0.09, 1.84] | 98.22% | [-0.06, 0.06] |        0% | 1.001 | 5206.00

# YAL1 - ΔVel2_IDD
model.YAL1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_PDM[c(25:27,34:36),],
                               iter = 6000)
summary(model.YAL1.ΔVel2_IDD)
posterior_summary(model.YAL1.ΔVel2_IDD)
describe_posterior(model.YAL1.ΔVel2_IDD)
rope(model.YAL1.ΔVel2_IDD)
posterior = as_draws_df(model.YAL1.ΔVel2_IDD)
mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.08741667

Parameter          | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [ 0.84,  1.14] | 99.99% | [-0.04, 0.04] |        0% | 1.000 | 4930.00
TreatmentΔVel2_IDD |  -0.74 | [-0.95, -0.50] | 99.99% | [-0.04, 0.04] |        0% | 1.000 | 5574.00

# PLATE ----
PCR_raw_PLATE   = PCR_raw %>% filter(Stage == "PLATE")
PCR_raw_PLATE$Treatment = as.factor(PCR_raw_PLATE$Treatment)

# MET1 - ΔVel1
model.MET1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[1:6,],
                           iter = 6000)
summary(model.MET1.ΔVel1)
posterior_summary(model.MET1.ΔVel1)
describe_posterior(model.MET1.ΔVel1)
rope(model.MET1.ΔVel1)
posterior = as_draws_df(model.MET1.ΔVel1)
1- mean(posterior$b_TreatmentΔVel1 > 0) # 0.001

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   0.97 | [-2.75,  4.64] | 74.36% | [-0.73, 0.73] |    31.39% | 1.001 | 5490.00
TreatmentΔVel1 |  12.95 | [ 7.82, 18.18] | 99.90% | [-0.73, 0.73] |        0% | 1.000 | 5593.00

# MET1 - ΔVel2
model.MET1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(1:3,7:9),],
                           iter = 6000)
summary(model.MET1.ΔVel2)
posterior_summary(model.MET1.ΔVel2)
describe_posterior(model.MET1.ΔVel2)
rope(model.MET1.ΔVel2)
posterior = as_draws_df(model.MET1.ΔVel2)
mean(posterior$b_TreatmentΔVel2 > 0) # 0.0105 

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.74,  1.26] | 99.97% | [-0.03, 0.03] |        0% | 1.000 | 4478.00
TreatmentΔVel2 |  -0.51 | [-0.90, -0.13] | 98.95% | [-0.03, 0.03] |        0% | 1.000 | 5381.00

# MET1 - ΔVel2_IDD
model.MET1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(1:3,10:12),],
                               iter = 6000)
summary(model.MET1.ΔVel2_IDD)
posterior_summary(model.MET1.ΔVel2_IDD)
describe_posterior(model.MET1.ΔVel2_IDD)
rope(model.MET1.ΔVel2_IDD)
posterior = as_draws_df(model.MET1.ΔVel2_IDD)
1- mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.05191667

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   0.99 | [ 0.34, 1.58] | 99.10% | [-0.04, 0.04] |        0% | 1.000 | 5493.00
TreatmentΔVel2_IDD |   0.64 | [-0.26, 1.59] | 94.81% | [-0.04, 0.04] |     1.53% | 1.000 | 3459.00

# CTZ1 - ΔVel1
model.CTZ1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(13:18),],
                           iter = 6000)
summary(model.CTZ1.ΔVel1)
posterior_summary(model.CTZ1.ΔVel1)
describe_posterior(model.CTZ1.ΔVel1)
rope(model.CTZ1.ΔVel1)
posterior = as_draws_df(model.CTZ1.ΔVel1)
mean(posterior$b_TreatmentΔVel1 > 0) # 0.0000

Parameter      | Median |         95% CI |   pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.98,  1.02] | 100% | [-0.02, 0.02] |        0% | 1.001 | 4534.00
TreatmentΔVel1 |  -0.38 | [-0.40, -0.35] | 100% | [-0.02, 0.02] |        0% | 1.001 | 4788.00

# CTZ1 - ΔVel2
model.CTZ1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(13:15,19:21),],
                           iter = 6000)
summary(model.CTZ1.ΔVel2)
posterior_summary(model.CTZ1.ΔVel2)
describe_posterior(model.CTZ1.ΔVel2)
rope(model.CTZ1.ΔVel2)
posterior = as_draws_df(model.CTZ1.ΔVel2)
mean(posterior$b_TreatmentΔVel2 > 0) # 0.001666667

Parameter      | Median |         95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
-----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.90,  1.10] |   100% | [-0.02, 0.02] |        0% | 1.000 | 5402.00
TreatmentΔVel2 |  -0.35 | [-0.49, -0.22] | 99.83% | [-0.02, 0.02] |        0% | 1.000 | 6466.00

# CTZ1 - ΔVel2_IDD
model.CTZ1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(13:15,22:24),],
                               iter = 9000)
summary(model.CTZ1.ΔVel2_IDD)
posterior_summary(model.CTZ1.ΔVel2_IDD)
describe_posterior(model.CTZ1.ΔVel2_IDD)
rope(model.CTZ1.ΔVel2_IDD)
posterior = as_draws_df(model.CTZ1.ΔVel2_IDD)
1-mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.1992778

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [ 0.57, 1.45] | 99.79% | [-0.02, 0.02] |        0% | 1.000 | 8344.00
TreatmentΔVel2_IDD |   0.19 | [-0.42, 0.80] | 80.07% | [-0.02, 0.02] |     4.86% | 1.000 | 9293.00

# YAL1 - ΔVel1
model.YAL1.ΔVel1     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(25:30),],
                           iter = 6000)
summary(model.YAL1.ΔVel1)
posterior_summary(model.YAL1.ΔVel1)
describe_posterior(model.YAL1.ΔVel1)
rope(model.YAL1.ΔVel1)
posterior = as_draws_df(model.YAL1.ΔVel1)
1-mean(posterior$b_TreatmentΔVel1 > 0) # 0.02975

Parameter      | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
----------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [ 0.40, 1.59] | 99.38% | [-0.05, 0.05] |        0% | 1.001 | 4799.00
TreatmentΔVel1 |   0.82 | [-0.06, 1.68] | 97.02% | [-0.05, 0.05] |     1.13% | 1.000 | 4926.00

# YAL1 - ΔVel2
model.YAL1.ΔVel2     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(25:27,31:33),],
                           iter = 6000)
summary(model.YAL1.ΔVel2)
posterior_summary(model.YAL1.ΔVel2)
describe_posterior(model.YAL1.ΔVel2)
rope(model.YAL1.ΔVel2)
posterior = as_draws_df(model.YAL1.ΔVel2)
1-mean(posterior$b_TreatmentΔVel2 > 0) # 0.01083333

Parameter      | Median |       95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
---------------------------------------------------------------------------------------------
(Intercept)    |   1.00 | [0.51, 1.48] | 99.61% | [-0.05, 0.05] |        0% | 1.002 | 3399.00
TreatmentΔVel2 |   0.88 | [0.22, 1.54] | 98.92% | [-0.05, 0.05] |        0% | 1.001 | 4767.00

# YAL1 - ΔVel2_IDD
model.YAL1.ΔVel2_IDD     = brm(Expression ~ Treatment, data = PCR_raw_PLATE[c(25:27,34:36),],
                               iter = 6000)
summary(model.YAL1.ΔVel2_IDD)
posterior_summary(model.YAL1.ΔVel2_IDD)
describe_posterior(model.YAL1.ΔVel2_IDD)
rope(model.YAL1.ΔVel2_IDD)
posterior = as_draws_df(model.YAL1.ΔVel2_IDD)
mean(posterior$b_TreatmentΔVel2_IDD > 0) # 0.3263333

Parameter          | Median |        95% CI |     pd |          ROPE | % in ROPE |  Rhat |     ESS
--------------------------------------------------------------------------------------------------
(Intercept)        |   1.00 | [ 0.55, 1.47] | 99.81% | [-0.02, 0.02] |        0% | 1.000 | 3235.00
TreatmentΔVel2_IDD |  -0.11 | [-0.77, 0.51] | 67.37% | [-0.02, 0.02] |     6.45% | 1.000 | 3704.00