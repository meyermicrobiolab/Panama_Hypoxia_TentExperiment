library(DHARMa)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(gridExtra)
library(dplyr)
library(car)
library(logistf)

remove(list = ls())

# bleaching ---------------------------------------------------------------------
df <- read.csv(file = "survival_tent.csv")
str(df)

# Full Tent (Dark Blue), Partial Tent (Green) and Open Plot (Orange).
# open= green
# partial= orange
# FT= blue

# ordered: full, partial, open
df$treatment <- ordered(df$treatment, levels = c("hypoxic", "tented_control", "control"))

treat.labs <- c("fully closed", "partially closed", "open")
names(treat.labs) <- c("hypoxic","tented_control", "control")
treat_col <- c("darkblue", "#E69F00","#009E73")

spp.labs <- c("A. tenuifolia", "S. siderea")             
names(spp.labs) <- c("A_tenuifolia", "S_siderea")

day.labs <- c("Before (0h)","After (24h)", "Recovery (72h)")
names(day.labs) <- c("Day_0", "Day_2", "Day_5")
df$bleaching<-df$bleaching*100

df.summary2 <- df %>% 
  group_by(day, treatment, spp) %>%
  summarise(
    sd = sd(bleaching),
    n = n(), 
    se = sd / sqrt(n),
    bleaching = mean(bleaching))
df.summary2

p_bleach <-  ggplot(df, aes(day, bleaching, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.4)) + 
  geom_line(aes(group = treatment), data = df.summary2) +
  geom_point(aes(group = treatment), size=3,data = df.summary2) +
  geom_errorbar(aes(ymin = bleaching-se, ymax = bleaching + se),data = df.summary2, width = 0.1)+
  theme_bw(base_size = 14) +
  scale_x_discrete(labels =day.labs) +
  scale_color_manual(values = treat_col, labels = treat.labs, name = "Treatment") +
  labs(x = "Time",
       y =  'Bleaching (%)') +
  theme(legend.position="bottom") +
  theme(strip.text = element_text(face = "italic"))+
  facet_grid(spp ~ ., scales = "fixed", 
             labeller = labeller(spp = spp.labs))
p_bleach

# live tissue:
df.summary3 <- df %>% 
  group_by(day, treatment, spp) %>%
  summarise(
    sd = sd(Live_tissue),
    n = n(), 
    se = sd / sqrt(n),
    Live_tissue = mean(Live_tissue)
  )
df.summary3


## stats ---------------------------------------------------------------------
df$bleaching<-df$bleaching/100
df_aten <- df[which(df$spp == "A_tenuifolia"),]
df_sid <- df[which(df$spp != 'A_tenuifolia'),]
str(df_sid)

# aten - no response - no test
# ssid — Firth penalized logistic regression used due to complete separation
# (control and tented_control never bleached; standard GLM fails to converge)

df_sid$treatment <- factor(df_sid$treatment, ordered = FALSE) # bleaching process for understanding the results
m_bleach_ssid_firth <- logistf(bleaching ~ treatment + day, data = df_sid) 
summary(m_bleach_ssid_firth)

# Bleaching in S. siderea was significantly associated with treatment and 
# day (Firth penalized logistic regression; LR χ²₄ = 17.70, p = 0.001). Hypoxic 
# corals had substantially higher odds of bleaching than tented controls or open 
# controls (OR = 0.027 for controls vs. hypoxic, 95% CI [0.0002, 0.40], p = 0.005).
# Bleaching was significantly elevated on Day 2 relative to Day 0 (OR = 43.6, p = 0.005) but did not
# differ significantly on Day 5 (OR = 3.6, p = 0.43), suggesting partial recovery."


# Pairwise day comparisons, S. siderea, fully closed (hypoxic) treatment only:
# Bleaching only occurred in the hypoxic treatment, so day contrasts are only
# meaningful there. Tukey on aov() is not appropriate for a binary response, so
# each pair of days is tested with a Firth logistic model, releveling the
# reference day for each comparison, then p-values are Holm-adjusted.

df_sid_hyp <- droplevels(subset(df_sid, treatment == "hypoxic"))
df_sid_hyp$day <- factor(df_sid_hyp$day)
table(df_sid_hyp$day, df_sid_hyp$bleaching) # n = 6 corals per day

day_pairs <- combn(levels(df_sid_hyp$day), 2, simplify = FALSE)
day_contrasts <- lapply(day_pairs, function(p) {
  dd <- df_sid_hyp
  # p[1] is the reference, so the coefficient for p[2] is log-odds of p[2] vs p[1]
  dd$day <- relevel(dd$day, ref = p[1])
  m <- logistf(bleaching ~ day, data = dd)
  nm <- paste0("day", p[2])
  data.frame(contrast = paste(p[2], "-", p[1]),
             log_odds = unname(m$coefficients[nm]),
             lower95  = unname(m$ci.lower[nm]),   # profile penalized likelihood CI
             upper95  = unname(m$ci.upper[nm]),
             p        = unname(m$prob[nm]))
})
day_contrasts <- do.call(rbind, day_contrasts)

# Holm adjustment across the 3 day contrasts
day_contrasts$p_holm <- p.adjust(day_contrasts$p, method = "holm")
day_contrasts

# Day_2 > Day_0 (Holm p = 0.012); Day_5 vs Day_0 n.s. (p = 0.44);
# Day_5 vs Day_2 borderline (Holm p = 0.055).
# Day 0 had no bleached corals, so Firth's penalty drives the effect sizes:
# interpret direction and p-values, not the magnitudes of the log-odds.


# mortality ---------------------------------------------------------------
p_mort <-  ggplot(df, aes(day, Live_tissue, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.4)) + 
  geom_line(aes(group = treatment), data = df.summary3) +
  geom_point(aes(group = treatment), size=3,data = df.summary3) +
  geom_errorbar(aes(ymin = Live_tissue-sd, ymax = Live_tissue + sd), data = df.summary3, width = 0.2)+
  theme_bw(base_size = 14) +
  scale_x_discrete(labels =day.labs) +
  scale_color_manual(values = treat_col, labels = treat.labs, name = "Treatment") +
  labs(x = "Time",
       y =  'Live tissue (%)') +
  theme(legend.position="bottom") +
  theme(strip.text = element_text(face = "italic"))+
  facet_grid(spp ~ ., scales = "fixed", 
             labeller = labeller(spp = spp.labs)) 
p_mort

## stats ---------------------------------------------------------------------
# aten only, as sid had no response/change
# did any tissue loss occur? recode as binary:
# interesting question is really how tissue loss varies across days within the hypoxic group:
# Live_tissue is a percentage, so any loss means < 100 (< 1 would only flag fully dead colonies)

tissue_loss <- as.integer(df_aten$Live_tissue < 100)
logistf(tissue_loss ~ day, data = df_aten, subset = treatment == "hypoxic")
m_mort_aten_firth <- logistf(tissue_loss ~ day, data = df_aten, subset = treatment == "hypoxic")
summary(m_mort_aten_firth)
exp(cbind(OR = coef(m_mort_aten_firth),
          lower = m_mort_aten_firth$ci.lower,
          upper = m_mort_aten_firth$ci.upper))

# Pairwise day comparisons, A. tenuifolia, fully closed (hypoxic) only
# Same approach as the S. siderea bleaching contrasts above: each pair of days
# is tested with a Firth logistic model (releveling the reference day), and
# p-values are Holm-adjusted across the 3 contrasts.

df_aten_hyp <- droplevels(subset(df_aten, treatment == "hypoxic"))
df_aten_hyp$tissue_loss <- as.integer(df_aten_hyp$Live_tissue < 100) # any loss
df_aten_hyp$day <- factor(df_aten_hyp$day)
table(df_aten_hyp$day, df_aten_hyp$tissue_loss) # no loss on Day 0; all corals lost tissue by Day 2

aten_day_pairs <- combn(levels(df_aten_hyp$day), 2, simplify = FALSE)
aten_day_contrasts <- lapply(aten_day_pairs, function(p) {
  dd <- df_aten_hyp
  # p[1] is the reference, so the coefficient for p[2] is log-odds of p[2] vs p[1]
  dd$day <- relevel(dd$day, ref = p[1])
  m <- logistf(tissue_loss ~ day, data = dd)
  nm <- paste0("day", p[2])
  data.frame(contrast = paste(p[2], "-", p[1]),
             log_odds = unname(m$coefficients[nm]),
             lower95  = unname(m$ci.lower[nm]),   # profile penalized likelihood CI
             upper95  = unname(m$ci.upper[nm]),
             p        = unname(m$prob[nm]))
})
aten_day_contrasts <- do.call(rbind, aten_day_contrasts)

# Holm adjustment across the 3 day contrasts
aten_day_contrasts$p_holm <- p.adjust(aten_day_contrasts$p, method = "holm")
aten_day_contrasts

# Significance letters for Fig. 1 -------------------------------------------
# Compact letters from the Holm-adjusted Firth contrasts above. Days sharing a
# letter are not significantly different (Holm p >= 0.05). Letters compare days
# within the fully closed (hypoxic) treatment of a single species only, so they
# are added to one panel per plot: S. siderea (bleaching), A. tenuifolia (tissue loss).
letters_from <- function(tab) {
  # multcompLetters() needs a named vector of p-values with names like "A-B"
  pv <- setNames(tab$p_holm, gsub(" - ", "-", tab$contrast))
  multcompView::multcompLetters(pv, threshold = 0.05)$Letters
}

# y = 112 puts letters in a row above the data (max 100%) so they don't overlap points
lab_bleach <- df.summary2 |>
  ungroup() |>
  filter(treatment == "hypoxic", spp == "S_siderea") |>
  mutate(label = letters_from(day_contrasts)[day], y = 112)
lab_bleach

lab_mort <- df.summary3 |>
  ungroup() |>
  filter(treatment == "hypoxic", spp == "A_tenuifolia") |>
  mutate(label = letters_from(aten_day_contrasts)[day], y = 112)
lab_mort

# treatment/day/spp columns in the label data place and color the letters correctly
p_bleach_lab <- p_bleach +
  geom_text(aes(y = y, label = label), data = lab_bleach,
            show.legend = FALSE, size = 5, fontface = "bold")
p_mort_lab <- p_mort +
  geom_text(aes(y = y, label = label), data = lab_mort,
            show.legend = FALSE, size = 5, fontface = "bold")

# Fig. 1 ------------------------------------------------------------------
Figure1 <- ggarrange(p_bleach_lab, p_mort_lab, 
                     ncol= 2, nrow = 1,
                     labels = c("A)","B)"),
                     common.legend = TRUE, 
                     legend = "bottom")

Figure1

ggsave(Figure1, filename = "Figure1.pdf", width=10, height=7, dpi=600, bg = "transparent")

# PAM -------------------------------------------------------------------
df_pam <- read.csv(file = "PAM.csv")
str(df_pam)
df_pam$treatment <- ordered(df_pam$treatment, levels = c("hypoxic", "tented_control", "control"))

p_pam <- ggboxplot(df_pam, x = "treatment", y = "FoFm",
                   fill = "treatment",
                   add = "jitter") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = treat_col) +
  scale_x_discrete(labels = treat.labs)+
  theme(legend.position="none") +
  labs(x = "treatment", y = "Photosynthetic yield (Fo:Fm)") +
  facet_grid(spp ~ .) +  
  stat_compare_means(method = "anova", label.y = .75)+ # Add global p-value
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "hypoxic")
p_pam

# nice figure:

#ordered: 

df.summary4 <- df_pam %>% 
  group_by(treatment, spp) %>%
  summarise(
    sd = sd(FoFm),
    n = n(), 
    se = sd / sqrt(n),
    FoFm = mean(FoFm))
df.summary4

p_pam <-  ggplot(df_pam, aes(treatment, FoFm, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) + 
  geom_line(aes(group = treatment), data = df.summary4) +
  geom_point(aes(group = treatment), size=3, data = df.summary4) +
  geom_errorbar(aes(ymin = FoFm-sd, ymax = FoFm + sd), data = df.summary4, width = 0.2)+
  theme_bw(base_size = 14) +
  scale_x_discrete(labels =treat.labs) +
  scale_color_manual(values = treat_col, labels = treat.labs, name = "Treatment") +
  labs(x = "Treatment", y = "Photosynthetic yield (Fo:Fm)") +
  theme(strip.text = element_text(face = "italic")) +
  theme(legend.position="bottom") +
  facet_grid(spp ~ ., scales = "fixed", 
             labeller = labeller(spp = spp.labs)) 
p_pam

## stats --------------------------------------------------------------------
df_aten_pam <- df_pam[ which(df_pam$spp == "A_tenuifolia"),]
df_sid_pam <- df_pam[which(df_pam$spp != 'A_tenuifolia'),]

#interpret results with tukey to make comparisons on plot
anova_aten <-aov(FoFm ~ treatment, data = df_aten_pam, na.action=na.omit)
summary(anova_aten) 
# tukey2 <- TukeyHSD(anova_aten)
# print(tukey2)
# tukey.cld <- multcompLetters4(anova_aten, tukey2)
# print(tukey.cld)

# 1. Homogeneity of variances
plot(anova_aten, 1)
leveneTest(FoFm ~ treatment, data = df_aten_pam) # homogeneity of variance is sig. bc of the corals that are dead/dying here

# what about other tests that don't care about homogeneity?
oneway.test(FoFm ~ treatment, data = df_aten_pam)
pairwise.t.test(df_aten_pam$FoFm, df_aten_pam$treatment,
                p.adjust.method = "BH", pool.sd = FALSE)
# both other options show significance, so there is consistency

# 2. Normality
plot(anova_aten, 2) # eh.. ok-ish, outliers are there

# Extract the residuals
aov_residuals <- residuals(object = anova_aten)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

# Note that a non-parametric alternative to one-way ANOVA is Kruskal-Wallis rank sum test
# which can be used when ANOVA assumptions are not met.

kruskal.test(FoFm ~ treatment, data = df_aten_pam)

# it's got a sig. p value, and as a non-parametric alternative to ANOVA (rank sum test, it can be used when
# the normality assumptions aren't met, again we are good. 
# 0.05764

data_summary <- group_by(df_aten_pam, treatment) %>%
  summarise(mean=mean(FoFm), sd=sd(FoFm), n = n(),  
            se = sd / sqrt(n)) %>%
  arrange(desc(mean))

# adding the compact letter display to the table with means and sd
cld1 <- as.data.frame.list(tukey.cld$treatment)
data_summary$Tukey <- cld1$Letters
print(data_summary)

# ssids -  no response
kruskal.test(FoFm ~ treatment, data = df_sid_pam)

# Zoox --------------------------------------------------------------------
df_zoox <- read.csv(file = "SMR.csv")
df_zoox$treatment <- ordered(df_zoox$treatment, levels = c("hypoxic", "tented_control", "control"))

df.summary5 <- df_zoox %>% 
  group_by(treatment, spp) %>%
  summarise(
    sd = sd(zoox),
    n = n(), 
    se = sd / sqrt(n),
    zoox = mean(zoox))
df.summary5

p_zoox <-  ggplot(df_zoox, aes(treatment, zoox, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.01)) + 
  geom_line(aes(group = treatment), data = df.summary5) +
  #ylim(0,6)+
  geom_point(aes(group = treatment), size=3, data = df.summary5) +
  geom_errorbar(aes(ymin = zoox-sd, ymax = zoox + sd), data = df.summary5, width = 0.2)+
  theme_bw(base_size = 14) +
  scale_x_discrete(labels =treat.labs) +
  scale_color_manual(values = treat_col, labels = treat.labs, name = "Treatment") +
  labs(x = "Treatment",   
       y = expression(paste("Symbiodinium density (cells/cm"^2, ")"))) +
  theme(strip.text = element_text(face = "italic")) +
  theme(legend.position="bottom") +
  facet_grid(spp ~ ., scales = "free", 
             labeller = labeller(spp = spp.labs)) 
p_zoox

## stats -------------------------------------------------------------------
df_aten_zoo <- df_zoox[ which(df_zoox$spp == "A_tenuifolia"),]

anova_aten <-aov(zoox ~ treatment, data = df_aten_zoo, na.action=na.omit)
summary(anova_aten) 

# 1. Homogeneity of variances
plot(anova_aten, 1)
leveneTest(zoox ~ treatment, data = df_aten_zoo) # homogeneity of variance is OK

# 2. Normality
plot(anova_aten, 2) # eh.. ok-ish, outliers are there

# Extract the residuals
aov_residuals <- residuals(object = anova_aten)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals) # good

df_sid_zoo <- df_zoox[which(df_zoox$spp != 'A_tenuifolia'),]
anova_aten <-aov(zoox ~ treatment, data = df_sid_zoo, na.action=na.omit)
summary(anova_aten) # nada, again, assumptions good?

# 1. Homogeneity of variances
plot(anova_aten, 1)
leveneTest(zoox ~ treatment, data = df_sid_zoo) # homogeneity of variance is OK

# 2. Normality
plot(anova_aten, 2) # eh.. ok-ish, outliers are there

# Extract the residuals
aov_residuals <- residuals(object = anova_aten)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals) # good


# Fig. 2 ------------------------------------------------------------------

Figure2 <- ggarrange(p_zoox, p_pam, 
                     ncol= 2, nrow = 1,
                     labels = c("A)","B)"),
                     common.legend = TRUE, 
                     legend = "bottom")

Figure2

ggsave(Figure2, filename = "Figure2.pdf", width=10, height=7, dpi=600, bg = "transparent")

# SMR --------------------------------------------------------------------
df_smr <- read.csv(file = "SMR.csv")
str(df_smr)
df_smr$treatment <- ordered(df_smr$treatment, levels = c("hypoxic", "tented_control", "control"))


df.summary5 <- df_smr %>% 
  group_by(treatment, spp) %>%
  summarise(
    sd = sd(SMR_umolO2_hr_cm2),
    n = n(), 
    se = sd / sqrt(n),
    SMR_umolO2_hr_cm2 = mean(SMR_umolO2_hr_cm2))
df.summary5

p_smr <-  ggplot(df_smr, aes(treatment, SMR_umolO2_hr_cm2, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.001)) + 
  geom_point(aes(group = treatment), size=3, data = df.summary5) +
  geom_errorbar(aes(ymin = SMR_umolO2_hr_cm2-se, ymax = SMR_umolO2_hr_cm2 + se), data = df.summary5, width = 0.2)+
  theme_bw(base_size = 14) +
  scale_x_discrete(labels =treat.labs) +
  scale_color_manual(values = treat_col, labels = treat.labs, name = "Treatment") +
  labs(x = "Treatment",   
       y = expression(paste("Oxygen consumption (O" [2], " h" ^-1, "cm"^2, ")"))) +
  theme(legend.position="bottom") +
  theme(strip.text = element_text(face = "italic")) +
  facet_grid(. ~ spp, scales = "free", 
             labeller = labeller(spp = spp.labs)) 
p_smr

df.summary6 <- df_smr %>% 
  group_by(treatment, spp) %>%
  summarise(
    sd = sd(Abs_mgO2_min),
    n = n(), 
    se = sd / sqrt(n),
    Abs_mgO2_min = mean(Abs_mgO2_min))
df.summary6

p_smr_abs <-  ggplot(df_smr, aes(treatment, Abs_mgO2_min, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.001)) + 
  geom_point(aes(group = treatment), size=3, data = df.summary6) +
  geom_errorbar(aes(ymin = Abs_mgO2_min-se, ymax = Abs_mgO2_min + se), 
                data = df.summary6, width = 0.2)+
  theme_bw(base_size = 14) +
  scale_x_discrete(labels =treat.labs) +
  scale_color_manual(values = treat_col, labels = treat.labs, name = "Treatment") +
  labs(x = "Treatment",   
       y = expression(paste("Oxygen consumption (O" [2], " h" ^-1, ")"))) +
  theme(legend.position="bottom") +
  theme(strip.text = element_text(face = "italic")) +
  facet_grid(. ~ spp, scales = "free", 
             labeller = labeller(spp = spp.labs)) 
p_smr_abs

ggsave(p_smr, filename = "Figure3.pdf", width=8, height=5, 
       dpi=600, bg = "transparent")

## stats -------------------------------------------------------------------
# aten
df_aten_smr <- df_smr[ which(df_smr$spp == "A_tenuifolia"),]

anova_aten <-aov(SMR_umolO2_hr_cm2 ~ treatment, data = df_aten_smr, na.action=na.omit)
summary(anova_aten) 

# 1. Homogeneity of variances
plot(anova_aten, 1)
leveneTest(SMR_umolO2_hr_cm2 ~ treatment, data = df_aten_smr) # homogeneity of variance is OK

# 2. Normality
plot(anova_aten, 2) # eh.. ok-ish, outliers are there

# Extract the residuals
aov_residuals <- residuals(object = anova_aten)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals) # good

df_sid_smr <- df_smr[which(df_smr$spp != 'A_tenuifolia'),]
anova_sid <-aov(SMR_umolO2_hr_cm2 ~ treatment, data = df_sid_smr, na.action=na.omit)
summary(anova_sid) # nada, again, assumptions good?

# 1. Homogeneity of variances
plot(anova_sid, 1)
leveneTest(SMR_umolO2_hr_cm2 ~ treatment, data = df_sid_smr) # homogeneity of variance is OK

# 2. Normality
plot(anova_sid, 2) # eh.. ok-ish, outliers are there

# Extract the residuals
aov_residuals <- residuals(object = anova_sid)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals) # good

