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

# LLM1----

plates   = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/plates.txt",dec=".")
model    = aov(Normalized ~ Biological.Rep, data = plates)

# Assumptions----
# 2. Normality
ols_plot_resid_qq(model)
# 1. Homogeneity of variances
ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Biological.Rep, data = plates)

# Tukey----
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Figure----

# Testing a non-parametric model
stat.test = aov(Normalized ~ Biological.Rep, data = plates) %>%
  tukey_hsd()

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_melanization = ggboxplot(plates, x = "Biological.Rep", y = "Normalized",
          ylab = "Normalized melanization", xlab = "") + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                     y.position = c(1.7, 1.8, 1.9, 1.7,2.0,1.9))
figure_melanization

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_melanization.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_melanization)
dev.off()

# Plants total ----

plants.total = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/new_data_plant.txt",dec=".")

# reorganizing the data

# healthy == 1; weak.symptoms == 2; strong.symptoms == 3; 
# very.strong.symptoms == 4

MOCK      = c(rep(1, each = 177), rep(2, each = 3))
sym1      = rep("MOCK", each = length(MOCK))

WT        = c(rep(1, each = 31),rep(2, each = 66),
         rep(3, each = 65),rep(4, each = 18))
sym2      = rep("WT", each = length(WT))

ΔLLM1     = c(rep(1, each = 59), rep(2, each = 54),
          rep(3, each = 32),rep(4, each = 5))
sym3      = rep("ΔLLM1", each = length(ΔLLM1))

comp.LLM1 = c(rep(1, each = 11), rep(2, each = 14),
              rep(3, each = 32),rep(4, each = 3))
sym4      = rep("comp.LLM1", each = length(comp.LLM1))

LLM1.GFP  = c(rep(1, each = 10), rep(2, each = 7),
              rep(3, each = 10),rep(4, each = 3))
sym5      = rep("LLM1.GFP", each = length(LLM1.GFP))

oeLLM1    = c(rep(1, each = 59), rep(2, each = 68),
              rep(3, each = 44),rep(4, each = 9))
sym6      = rep("oeLLM1", each = length(oeLLM1))

treatments = c(sym1,sym2,sym3,sym4,sym5,sym6)
rank       = c(MOCK,WT,ΔLLM1,comp.LLM1,LLM1.GFP,oeLLM1)

total     = as.data.frame(cbind(treatments,rank))

# Aligned ranks anova ----
# https://rcompanion.org/handbook/F_16.html

total$treatments = factor(total$treatments,
                      levels=unique(total$treatments))

model = art(as.numeric(rank) ~ treatments,data = total)
anova(model)

# Post-hoc comparisons ----

model.lm = artlm(model, "treatments")
marginal = emmeans(model.lm,~ treatments)
test     = pairs(marginal,adjust = "tukey")

# Figure plant ----

df = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/df_LLM1_plant.txt")

df$rank = factor(df$rank,levels=unique(df$rank))
df$treatments = factor(df$treatments,levels=unique(df$treatments))

figure_plant = ggplot(df, aes(treatments, Percent, fill = (rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_plant.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_plant)
dev.off()

# Plants without first treatment ----

plants.total.1 = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/new_data_plant_without.txt",dec=".")

# reorganizing the data

# healthy == 1; weak.symptoms == 2; strong.symptoms == 3; 
# very.strong.symptoms == 4

MOCK      = c(rep(1, each = 148), rep(2, each = 2))
sym1      = rep("MOCK", each = length(MOCK))

WT        = c(rep(1, each = 22),rep(2, each = 52),
              rep(3, each = 58),rep(4, each = 18))
sym2      = rep("WT", each = length(WT))

ΔLLM1     = c(rep(1, each = 25), rep(2, each = 28),
              rep(3, each = 32),rep(4, each = 5))
sym3      = rep("ΔLLM1", each = length(ΔLLM1))

comp.LLM1 = c(rep(1, each = 11), rep(2, each = 14),
              rep(3, each = 32),rep(4, each = 3))
sym4      = rep("comp.LLM1", each = length(comp.LLM1))

LLM1.GFP  = c(rep(1, each = 10), rep(2, each = 7),
              rep(3, each = 10),rep(4, each = 3))
sym5      = rep("LLM1.GFP", each = length(LLM1.GFP))

oeLLM1    = c(rep(1, each = 31), rep(2, each = 45),
              rep(3, each = 35),rep(4, each = 9))
sym6      = rep("oeLLM1", each = length(oeLLM1))

treatments = c(sym1,sym2,sym3,sym4,sym5,sym6)
rank       = c(MOCK,WT,ΔLLM1,comp.LLM1,LLM1.GFP,oeLLM1)

total.1     = as.data.frame(cbind(treatments,rank))

# Aligned ranks anova ----

total.1$treatments = factor(total.1$treatments,
                          levels=unique(total.1$treatments))

model.1 = art(as.numeric(rank) ~ treatments,data = total.1)
anova(model.1)

# Post-hoc comparisons ----

model.lm1 = artlm(model.1, "treatments")
marginal.1 = emmeans(model.lm1,~ treatments)
pairs(marginal.1,adjust = "tukey")

# Sporulation ----

sporas = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/sporulation.txt")
model  = aov(Normalized ~ Biological.Rep, data = sporas)

# Assumptions----
# 2. Normality
ols_plot_resid_qq(model)
# 1. Homogeneity of variances
ols_plot_resid_fit(model)
bartlett.test(Normalized ~ Biological.Rep, data = sporas)

# Tukey----
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# Testing a non-parametric model
stat.test = aov(Normalized ~ Biological.Rep, data = sporas) %>%
  tukey_hsd()

stat.test = stat.test %>% filter(!(p.adj.signif == "ns"))

figure_spores = ggbarplot(sporas, x = "Biological.Rep", y = "Normalized", 
          ylab = "Normalized conidia formation", xlab = "", add = "mean_se") + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", 
    y.position = c(1.1, 1.2, 1.3, 1.4))
figure_spores

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_sporulation.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_spores)
dev.off()

# Real time PCR LLM1 ----

realtimePCR     = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/realtimePCR.txt")

# SXM (normal not homogeneous variances) ----
realtimePCR.SXM   = realtimePCR %>% filter(Medium == "SXM")
model.SXM         = aov((Values) ~ Treatments, data = realtimePCR.SXM)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.SXM)
e = resid(model.SXM)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.SXM)
bartlett.test((Values) ~ Treatments, data = realtimePCR.SXM)

# PDM (normal not homogeneous variances) ----
realtimePCR.PDM   = realtimePCR %>% filter(Medium == "PDM")
model.PDM         = aov((Values) ~ Treatments, data = realtimePCR.PDM)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.PDM)
e = resid(model.PDM)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.PDM)
bartlett.test(Values ~ Treatments, data = realtimePCR.PDM)

# PLATE (not normal not homogeneous variances) ----
realtimePCR.PLATE = realtimePCR %>% filter(Medium == "PLATE")
model.PLATE       = aov((Values) ~ Treatments, data = realtimePCR.PLATE)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.PLATE)
e = resid(model.PLATE)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.PLATE)
bartlett.test(Values ~ Treatments, data = realtimePCR.PLATE)

#library(dunn.test)
#dunn.test(realtimePCR.SXM$Values, realtimePCR.SXM$Treatments, kw=TRUE,
#          method="none")

# Plotting

stat.test = ConoverTest(realtimePCR.SXM$Values, realtimePCR.SXM$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
# stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("Vel1","Vel2","WT","Vel2","WT","WT")
group2    = c("IDD","IDD","IDD","Vel1","Vel1","Vel2")
stat.test.SXM = cbind(stat.test,group1,group2)
stat.test.SXM = stat.test.SXM %>% replace(is.na(.), "***")

figure_pcr.SXM = ggbarplot(realtimePCR.SXM, x = "Treatments", y = "Values", 
                       ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.SXM), label = "p.adj.signif",
                     y.position = c(18, 23, 25, 29, 32, 27))

figure_pcr.SXM

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.SXM.LLM1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.SXM)
dev.off()

# PDM ----
#PDM.pairwise = ConoverTest(realtimePCR.PDM$Values, realtimePCR.PDM$Treatments,method="holm")
#PDM.pairwise

#dunn.test(realtimePCR.PDM$Values, realtimePCR.PDM$Treatments, kw=TRUE,
#          method="none")

# Plotting

stat.test = ConoverTest(realtimePCR.PDM$Values, realtimePCR.PDM$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("Vel1","Vel2","WT","WT")
group2    = c("IDD","IDD","Vel1","Vel2")
stat.test.PDM = cbind(stat.test,group1,group2)

figure_pcr.PDM = ggbarplot(realtimePCR.PDM, x = "Treatments", y = "Values", 
                       ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.PDM), label = "p.adj.signif",
                     y.position = c(8, 12, 16, 20))

figure_pcr.PDM

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.PDM.LLM1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.PDM)
dev.off()

# PLATE ----
#PLATE.pairwise = ConoverTest(realtimePCR.PLATE$Values, realtimePCR.PLATE$Treatments,method="holm")
#PLATE.pairwise

#dunn.test(realtimePCR.PLATE$Values, realtimePCR.PLATE$Treatments, kw=TRUE,
#          method="none",list=TRUE)

# Plotting

stat.test = ConoverTest(realtimePCR.PLATE$Values, realtimePCR.PLATE$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("Vel1","Vel2","WT","WT")
group2    = c("IDD","IDD","Vel1","Vel2")
stat.test.PLATE = cbind(stat.test,group1,group2)

figure_pcr.PLATE = ggbarplot(realtimePCR.PLATE, x = "Treatments", y = "Values", 
                       ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.PLATE), label = "p.adj.signif",
                     y.position = c(26, 12, 29, 32))

figure_pcr.PLATE

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.PLATE.LLM1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.PLATE)
dev.off()

# AML ----

# Plants total ----

total_AML = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/new_data_plant_AML.txt")

# Aligned ranks anova ----
# https://rcompanion.org/handbook/F_16.html

total_AML$treatments = factor(total_AML$treatments,
                              levels=unique(total_AML$treatments))

model = art(as.numeric(rank) ~ treatments,data = total_AML)
anova(model)

# Post-hoc comparisons ----

model.lm = artlm(model, "treatments")
marginal = emmeans(model.lm,~ treatments)
test     = pairs(marginal,adjust = "tukey")

# Figure plant ----

df = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/df_AML_plant.txt")

df$rank = factor(df$rank,levels=unique(df$rank))
df$treatments = factor(df$treatments,levels=unique(df$treatments))

figure_plant_AML = ggplot(df, aes(treatments, Percent, fill = (rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("")
figure_plant_AML

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_plant_AML.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_plant_AML)
dev.off()

# Real time PCR AML1 ----
# Not normally and not homogenous variances

realtimePCR = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/realtimePCR_AML1.txt")

# SXM (not normal not homogeneous variances) ----
realtimePCR.SXM   = realtimePCR %>% filter(Medium == "SXM")
model.SXM         = aov((Values) ~ Treatments, data = realtimePCR.SXM)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.SXM)
e = resid(model.SXM)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.SXM)
bartlett.test((Values) ~ Treatments, data = realtimePCR.SXM)

# Plotting

stat.test = ConoverTest(realtimePCR.SXM$Values, realtimePCR.SXM$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
#stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
#group1    = c("Vel2","WT","Vel2","WT","WT")
#group2    = c("IDD","IDD","Vel1","Vel1","Vel2")
#stat.test.SXM = cbind(stat.test,group1,group2)
stat.test.SXM = stat.test

figure_pcr.SXM.AML1 = ggbarplot(realtimePCR.SXM, x = "Treatments", y = "Values", 
                                ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") # + 
#  stat_pvalue_manual((stat.test.SXM), label = "p.adj.signif",
#                     y.position = c(10, 18, 22, 26,30))

figure_pcr.SXM.AML1

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.SXM.AML1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.SXM.AML1)
dev.off()

# PDM (normal not homogeneous variances) ----
realtimePCR.PDM   = realtimePCR %>% filter(Medium == "PDM")
model.PDM         = aov((Values) ~ Treatments, data = realtimePCR.PDM)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.PDM)
e = resid(model.PDM)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.PDM)
bartlett.test(Values ~ Treatments, data = realtimePCR.PDM)

# Plotting

stat.test = ConoverTest(realtimePCR.PDM$Values, realtimePCR.PDM$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("ΔVel1","ΔVel2","ΔVel1","ΔVel2")
group2    = c("IDD","IDD","WT","WT")
stat.test.PDM = cbind(stat.test,group1,group2)

figure_pcr.PDM.AML1 = ggbarplot(realtimePCR.PDM, x = "Treatments", y = "Values", 
                                ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.PDM), label = "p.adj.signif",
                     y.position = c(12, 18, 15, 8))

figure_pcr.PDM.AML1

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.PDM.AML1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.PDM.AML1)
dev.off()

# PLATE (normal not homogeneous variances) ----
realtimePCR.PLATE = realtimePCR %>% filter(Medium == "PLATE")
model.PLATE       = aov((Values) ~ Treatments, data = realtimePCR.PLATE)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.PLATE)
e = resid(model.PLATE)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.PLATE)
bartlett.test(Values ~ Treatments, data = realtimePCR.PLATE)

# Plotting

stat.test = ConoverTest(realtimePCR.PLATE$Values, realtimePCR.PLATE$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("ΔVel2")
group2    = c("WT")
stat.test.PLATE = cbind(stat.test,group1,group2)

figure_pcr.PLATE.AML1 = ggbarplot(realtimePCR.PLATE, x = "Treatments", y = "Values", 
                                  ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.PLATE), label = "p.adj.signif",
                     y.position = c(1.5))

figure_pcr.PLATE.AML1

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.PLATE.AML1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.PLATE.AML1)
dev.off()

# Real time PCR NML1 ----
# Not normally and not homogenous variances

realtimePCR = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/realtimePCR_NML1.txt")

# SXM (not normal not homogeneous variances) ----
realtimePCR.SXM   = realtimePCR %>% filter(Medium == "SXM")
model.SXM         = aov((Values) ~ Treatments, data = realtimePCR.SXM)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.SXM)
e = resid(model.SXM)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.SXM)
bartlett.test((Values) ~ Treatments, data = realtimePCR.SXM)

# Plotting

stat.test = ConoverTest(realtimePCR.SXM$Values, realtimePCR.SXM$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("ΔVel1","ΔVel2","ΔVel1","ΔVel2")
group2    = c("IDD","IDD","WT","WT")
stat.test.SXM = cbind(stat.test,group1,group2)

figure_pcr.SXM.NML1 = ggbarplot(realtimePCR.SXM, x = "Treatments", y = "Values", 
                                ylab = "Relative Normalized Expression", xlab = "", add = "mean_se")  + 
  stat_pvalue_manual((stat.test.SXM), label = "p.adj.signif",
                     y.position = c(4, 6, 8, 4))

figure_pcr.SXM.NML1

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.SXM.NML1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.SXM.NML1)
dev.off()

# PDM (normal not homogeneous variances) ----
realtimePCR.PDM   = realtimePCR %>% filter(Medium == "PDM")
model.PDM         = aov((Values) ~ Treatments, data = realtimePCR.PDM)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.PDM)
e = resid(model.PDM)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.PDM)
bartlett.test(Values ~ Treatments, data = realtimePCR.PDM)

# Plotting

stat.test = ConoverTest(realtimePCR.PDM$Values, realtimePCR.PDM$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("ΔVel1","ΔVel2")
group2    = c("IDD","IDD")
stat.test.PDM = cbind(stat.test,group1,group2)

figure_pcr.PDM.NML1 = ggbarplot(realtimePCR.PDM, x = "Treatments", y = "Values", 
                                ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.PDM), label = "p.adj.signif",
                     y.position = c(5, 7))

figure_pcr.PDM.NML1

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.PDM.NML1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.PDM.NML1)
dev.off()

# PLATE (normal not homogeneous variances) ----
realtimePCR.PLATE = realtimePCR %>% filter(Medium == "PLATE")
model.PLATE       = aov((Values) ~ Treatments, data = realtimePCR.PLATE)

# Assumptions
# 2. Normality
ols_plot_resid_qq(model.PLATE)
e = resid(model.PLATE)
shapiro.test(e)
# 1. Homogeneity of variances
ols_plot_resid_fit(model.PLATE)
bartlett.test(Values ~ Treatments, data = realtimePCR.PLATE)

# Plotting

stat.test = ConoverTest(realtimePCR.PLATE$Values, realtimePCR.PLATE$Treatments,method="holm")
stat.test = as.data.frame(stat.test[1])
stat.test$Treatments = row.names(stat.test)
stat.test = stat.test[, c("Treatments", "mean.rank.diff", "pval")]
stat.test = stat.test %>% mutate(p.adj.signif = case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                                          pval <= 0.001 & pval > 0.0001 ~ "***",
                                                          pval <= 0.0001 & pval > 0.00001 ~ "****",
                                                          pval > 0.05 ~ "NS"))
stat.test = stat.test %>% filter(!(p.adj.signif == "NS"))
group1    = c("ΔVel1")
group2    = c("IDD")
stat.test.PLATE = cbind(stat.test,group1,group2)

figure_pcr.PLATE.NML1 = ggbarplot(realtimePCR.PLATE, x = "Treatments", y = "Values", 
                                  ylab = "Relative Normalized Expression", xlab = "", add = "mean_se") + 
  stat_pvalue_manual((stat.test.PLATE), label = "p.adj.signif",
                     y.position = c(4))

figure_pcr.PLATE.NML1

png("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/Figures/figure_pcr.PLATE.NML1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(figure_pcr.PLATE.NML1)
dev.off()


