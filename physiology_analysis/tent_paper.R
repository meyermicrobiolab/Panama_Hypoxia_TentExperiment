library(DHARMa)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(gridExtra)
library(dplyr)
library(car)
remove(list = ls())

# bleaching ---------------------------------------------------------------------
df <- read.csv(file = "survival_tent.csv")
str(df)

# ordered: full, partial, open
df$treatment <- ordered(df$treatment, levels = c("hypoxic", "tented_control", "control"))

treat.labs <- c("fully closed", "partially closed", "open")
names(treat.labs) <- c("hypoxic","tented_control", "control")
treat_col <- c("#D55E00", "#0072B2","#56B4E9") # from colorblind friendly palette

spp.labs <- c("Agaricia tenuifolia", "Siderastrea siderea")             
names(spp.labs) <- c("A_tenuifolia", "S_siderea")

day.labs <- c("Before (0h)","During (24h)", "After (72h)")
names(day.labs) <- c("Day_0", "Day_2", "Day_5")
df$bleaching<-df$bleaching*100

df.summary2 <- df %>% 
  group_by(day, treatment, spp) %>%
  summarise(
    sd = sd(bleaching),
    n = n(), 
    se = sd / sqrt(n),
    bleaching = mean(bleaching)
  )
df.summary2

p_bleach <-  ggplot(df, aes(day, bleaching, color = treatment)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.4)) + 
  geom_line(aes(group = treatment), data = df.summary2) +
  geom_point(aes(group = treatment), size=3,data = df.summary2) +
  geom_errorbar(aes(ymin = bleaching-se, ymax = bleaching+se),data = df.summary2, width = 0.1)+
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
df_aten <- df[ which(df$spp == "A_tenuifolia"),]
df_sid <- df[which(df$spp != 'A_tenuifolia'),]
str(df_sid)

# aten
m_bleach_aten <- glm(data = df_aten, bleaching ~ treatment + day, 
                     na.action = na.exclude, 
                     family=binomial)
# check model assumptions, use Dharma becasue it can peak into these relationships with binary data:
summary(m_bleach_aten)
simulationOutput <- simulateResiduals(fittedModel = m_bleach_aten, plot = T)

# ssid
m_bleach_ssid <- glm(data = df_sid, bleaching ~ treatment + day, 
                     na.action = na.exclude, 
                     family=binomial)
# check model assumptions, use Dharma becasue it can peak into these relationships with binary data:
summary(m_bleach_ssid)
simulationOutput <- simulateResiduals(fittedModel = m_bleach_ssid, plot = T)

anova(m_bleach_ssid, test = "Chi") # likelihood tests to confirm statistically significant effects
#interpret results with tukey to make comparisons on plot
m <-aov(bleaching ~ treatment + day, data = df_sid, na.action=na.omit)
summary(m) # agrees with GLM results
tukey2 <- TukeyHSD(m)
print(tukey2)
tukey.cld <- multcompLetters4(m, tukey2)
print(tukey.cld)

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
df_aten$Live_tissue<-df_aten$Live_tissue/100
df_sid$Live_tissue<-df_sid$Live_tissue/100

# aten
m_mort_aten <- glm(data = df_aten, Live_tissue ~ treatment + day, 
                  family=quasibinomial,  na.action = na.exclude)
# check model assumptions, use Dharma becasue it can peak into these relationships with binary data:
summary(m_mort_aten)
simulationOutput <- simulateResiduals(fittedModel = m_mort_aten, plot = T)
anova(m_mort_aten, test = "Chi") # likelihood tests to confirm statistically significant effects

#interpret results with tukey to make comparisons on plot
m <-aov(Live_tissue ~ treatment * day, data = df_aten, na.action=na.omit)
summary(m) # agrees with GLM results
tukey2 <- TukeyHSD(m)
print(tukey2)
tukey.cld <- multcompLetters4(m, tukey2)
print(tukey.cld)

# ssid
m_mort_ssid <- glm(data = df_sid, Live_tissue ~ treatment + day, 
                     na.action = na.exclude, 
                     family=binomial)
# check model assumptions, use Dharma becasue it can peak into these relationships with binary data:
summary(m_mort_ssid)
simulationOutput <- simulateResiduals(fittedModel = m_mort_ssid, plot = T)
anova(m_mort_ssid, test = "Chi") # likelihood tests to confirm statistically significant effects
# nada, no response. thats what we see in the experiment. that is good.


# Fig. 1 ------------------------------------------------------------------
Figure1 <- ggarrange(p_bleach, p_mort, 
                     ncol= 2, nrow = 1,
                     labels = c("A)","B)"),
                     common.legend = TRUE, 
                     legend = "none")

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
       y = expression(paste("Symbiodiniaceae density (cells/cm"^2, ")"))) +
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
                     legend = "none")

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
  theme(legend.position="none") +
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

