# Melanization----

# Plates----
library(olsrr)
library(multcomp)
library(multcompView)
library(ggpubr)
library(ARTool)
library(emmeans)
library(rcompanion)
library(scales)
library(ggsignif)

plates   = read.delim("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/Project_8 (Reviews)/Family/Pere/doctorado/datasets/plates.txt",dec=".")
model    = aov(Normalized ~ variable, data = plates)

# Assumptions----
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

# Tukey----
summary(model)
TUKEY = TukeyHSD(model, conf.level=.95)
TUKEY

# t-test
t.test(plates$Normalized[1:6] , plates$Normalized[25:30],
       alternative = "two.sided", var.equal = FALSE)

# Figure----

my_comparisons = list( c("WT", "ΔLLM1"), c("WT", "LLM1-NoGFP"), 
                       c("WT", "LLM1+GFP-C"),  c("WT", "oeLLM1"),
                       c("ΔLLM1", "oeLLM1"))

figure_melanization = ggboxplot(plates, x = "Biological.Rep", y = "Normalized",
          ylab = "Normalized melanization", xlab = "") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
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

figure_plant = ggplot(total, aes(treatments, as.numeric(rank), fill = (rank))) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) + 
  scale_fill_manual(name = "", values = c("#ff7f00", "#fdc086", "#ffff99",
                                          "#7fc97f"),
                    breaks=c('4', '3', '2', '1'),
                    labels = c("very strong","strong","weak","healthy")) + 
  ylab("# Plants [%]") + xlab("") + 
  geom_signif(comparisons = list(c("WT", "ΔLLM1")), map_signif_level = TRUE,
              y_position = 1.1,tip_length = 0.01) + 
  geom_signif(comparisons = list(c("WT", "oeLLM1")), map_signif_level = TRUE,
                                            y_position = 1.05, tip_length = 0.01) + 
  geom_signif(comparisons = list(c("ΔLLM1", "oeLLM1")), map_signif_level = TRUE,
              y_position = 1.0, tip_length = 0.01) + 
  geom_signif(comparisons = list(c("WT", "comp.LLM1")), map_signif_level = TRUE,
              y_position = 0.90, tip_length = 0.01) + 
  geom_signif(comparisons = list(c("WT", "LLM1.GFP")), map_signif_level = TRUE,
              y_position = 0.95, tip_length = 0.01)
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


