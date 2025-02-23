setwd("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado")

# calling packages ####
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)

# calling the data ####
df         = read.delim("gene_normalized.txt",dec=".")
# One-sample t-test ####
# Chr1g15480 (p-value = 0.01456) ####
Chr1g15480 = df %>% filter(Code == "Chr1g15480")
shapiro.test(Chr1g15480$dVEL2)
res        = t.test(Chr1g15480$dVEL2, mu = 1)
res
# Chr4g08340 (p-value = 0.02605) ####
Chr4g08340 = df %>% filter(Code == "Chr4g08340")
shapiro.test(Chr4g08340$dVEL2)
res        = t.test(Chr4g08340$dVEL2, mu = 1)
res
# Chr4g09760  (p-value = 0.1482) ####
Chr4g09760 = df %>% filter(Code == "Chr4g09760")
shapiro.test(Chr4g09760$dVEL2)
res        = t.test(Chr4g09760$dVEL2, mu = 1)
res
# Chr1g17500  ####
Chr1g17500 = df %>% filter(Code == "Chr1g17500")
shapiro.test(Chr1g17500$dVEL2)
res        = t.test(Chr1g17500$dVEL2, mu = 1)
res
# Chr2g09080  ####
Chr2g09080 = df %>% filter(Code == "Chr2g09080")
shapiro.test(Chr2g09080$dVEL2)
res        = t.test(Chr2g09080$dVEL2, mu = 1)
res

df_1       = df %>% filter(Code == "Chr1g15480"|Code == "Chr4g08340"|Code == "Chr4g09760")
ggplot(data=df_1, aes(x=Code, y=dVEL2/3)) +
  geom_bar(stat="identity") + theme_minimal()

# 

# calling the data ####
df_08340   = read.delim("08340.txt",dec=".")
res <- wilcox.test(df_08340$WT,df_08340$c08340,exact = FALSE)
res

df_llm1    = read.delim("llm1.txt",dec=".")
res.1 <- wilcox.test(df_llm1$WT,df_llm1$cLlm1_1,exact = FALSE)
res.1
res.2 <- wilcox.test(df_llm1$WT,df_llm1$cLlm1_2,exact = FALSE)
res.2
res.3 <- wilcox.test(df_llm1$WT,df_llm1$Comp_Llm1,exact = FALSE)
res.3

# New analysis ----

OE_plant   = read.delim("OE_plant_experiment.txt",dec=".")
WT_d08340  = wilcox.test(OE_plant$WT,OE_plant$d_08340,exact = FALSE)
WT_d08340  # (p-value = 0.4954)
WT_oe08340 = wilcox.test(OE_plant$WT,OE_plant$oe08340,exact = FALSE)
WT_oe08340 # (p-value = 0.8827)
WT_oeLLM1  = wilcox.test(OE_plant$WT,OE_plant$oeLLM1,exact = FALSE)
WT_oeLLM1  # (p-value = 0.008983)

# New analysis (01/14/2024) ----

OE_plant   = read.delim("OE_plant_experiment_14_04_2024.txt",dec=".")
WT_oeLM1   = wilcox.test(OE_plant$WT,OE_plant$oeLlm1,exact = FALSE)
WT_oeLM1  # (p-value = 0.4037)
WT_dMet1   = wilcox.test(OE_plant$WT,OE_plant$d_Met1,exact = FALSE)
WT_dMet1 # (p-value = 0.93)

# New analysis (29/05/2024)
OE_plant   = read.delim("Plant_29_05_24.txt",dec=".")
OE_comb    = read.delim("Plant_29_05_24_comb.txt",dec=".")
WT_oeLM1   = wilcox.test(OE_plant$WT,OE_plant$oeLLM1.1a,exact = FALSE)
WT_oeLM1  # (p-value = 0.006009)
WT_oeLM2   = wilcox.test(OE_plant$WT,OE_plant$oeLLM1.2a,exact = FALSE)
WT_oeLM2  # (p-value = 0.03027)
WT_MET1   = wilcox.test(OE_plant$WT,OE_plant$ΔMET1,exact = FALSE)
WT_MET1  # (p-value = 0.533)
WT_comb   = wilcox.test(OE_plant$WT,OE_comb$oeLLM1,exact = FALSE)
WT_comb  # (p-value = 0.004743) # https://stats.stackexchange.com/questions/197235/wilcoxon-rank-test-where-sample-sizes-are-very-different

# New analysis (30/05/2024)
OE_plant   = read.delim("Plant_30_05_24.txt",dec=".")
WT_oeLM1   = wilcox.test(OE_plant$WT,OE_plant$oeLLM1,exact = FALSE)
WT_oeLM1

# Analysis (07/08/2024)
OE_plant   = read.delim("plant_07_08_24.txt",dec=".")
WT_oeLM1   = wilcox.test(OE_plant$WT,OE_plant$oeLLM1,exact = FALSE)
WT_oeLM1
WT_dLLM1   = wilcox.test(OE_plant$WT,OE_plant$dLLM1,exact = FALSE)
WT_dLLM1
oeLM1_dLLM1   = wilcox.test(OE_plant$oeLLM1,OE_plant$dLLM1,exact = FALSE)
oeLM1_dLLM1

# Plates
library(olsrr)
library(multcomp)
library(multcompView)
library(ggpubr)

plates   = read.delim("platos.txt",dec=".")
model    = aov(Normalized ~ variable, data = plates)

# Assumptions

par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(model, which = 3)
# 2. Normality
plot(model, which = 2)

# Assumptions 2
# 2. Normality
ols_plot_resid_qq(model)
# 1. Homogeneity of variances
ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Biological.Rep, data = plates)

# Tukey
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# t-test
t.test(plates$Normalized[1:6] , plates$Normalized[25:30],
       alternative = "two.sided", var.equal = FALSE)

# Figure

my_comparisons = list( c("WT", "ΔLLM1"), c("WT", "LLM1-NoGFP"), 
                       c("WT", "LLM1+GFP-C"),  c("WT", "oeLLM1"),
                       c("ΔLLM1", "oeLLM1"))

ggboxplot(plates, x = "Biological.Rep", y = "Normalized",
          ylab = "Normalized melanization", xlab = "") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")




