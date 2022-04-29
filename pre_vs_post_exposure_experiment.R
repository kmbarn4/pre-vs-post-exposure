#Data analysis for:
#Bd post-exposure prophylaxis experiment ("headstart")

### Survival analysis

# Load survival data
headstart_survival <- read.csv("/Users/kmbarnett/Desktop/2022 winter/headstart data/clean/headstart_survival.csv")
head(headstart_survival)

library("KMsurv")
library("survival")

# Cox proportional hazards model
survival_model <- coxph(Surv(mortality_day,mortality_status)~treatment, data=headstart_survival)
summary(survival_model)

## Bd Load Data
headstart_GE <- read.csv("/Users/kmbarnett/Desktop/2022 winter/headstart data/clean/headstart_Bd_loads_final.csv")
head(headstart_GE)

##Bd prevalence model
library(glmmTMB)
prev_model = glmmTMB(Infected ~ Treatment, family="binomial", data = headstart_GE)
null_prev_model = glmmTMB(Infected ~ 1, family="binomial", data = headstart_GE)
anova(prev_model, null_prev_model)

#upload prevalence data
prev <- read.csv("/Users/kmbarnett/Desktop/2022 winter/headstart data/data/prev_confidence_intervals.csv", header = T)
head(prev)

#plot prevalence
library(ggplot2)
positions <- c("control", "pre-exposure", "post-exposure")
ggplot(prev) + 
  geom_bar( aes(x=treatment, y=prevalence), stat="identity", colour="black", fill="skyblue", alpha=0.7) + 
  #scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, by = 0.25)) +
  scale_x_discrete(limits = positions) +
  geom_errorbar( aes(x=treatment, ymin=lower_CI, ymax=upper_CI), width=0.4, colour="black", alpha=0.9, size=1.3) + theme_classic()

##Bd infection intensity model
library(glmmTMB)
infection_intensity_model = glmmTMB(round(ALL.Bd.GE) ~ Treatment, family="nbinom2", ziformula = ~1, data = headstart_GE)
summary(infection_intensity_model)

##post-hoc pairwise comparisons of infection intensity
#bonferroni correction for multiple testing p = 0.017 

#comparing pre-exposure and post-exposure bd infection intensities
pre_vs_post_model = glmmTMB(round(ALL.Bd.GE) ~ Treatment, family="nbinom2", ziformula = ~1, data = subset(headstart_GE, Treatment != "control"))
summary(pre_vs_post_model)

#comparing control vs. pre-exposure bd infection intensities
asw_vs_pre_model = glmmTMB(round(ALL.Bd.GE) ~ Treatment, family="nbinom2", ziformula = ~1, data = subset(headstart_GE, Treatment != "post-exposure"))
summary(asw_vs_pre_model)

#comparing control vs. post-exposure bd infection intensities
asw_vs_post_model = glmmTMB(round(ALL.Bd.GE) ~ Treatment, family="nbinom2", ziformula = ~1, data = subset(headstart_GE, Treatment != "pre-exposure"))
summary(asw_vs_post_model)

#extract mean Bd intensity from model using emmeans function (response scale)
#this doesn't really work since it gives negative confidence intervals
library(emmeans)
emm = emmeans(infection_intensity_model, ~Treatment, type = "response")
emm

#extract mean Bd intensity from model using emmeans function (log scale)
library(emmeans)
emm_log = emmeans(infection_intensity_model, ~Treatment)
emm_log

#emmeans csv
log_means <- read.csv("/Users/kmbarnett/Desktop/2022 winter/headstart data/clean/log_emmeans_inf_intensity.csv")
head(log_means)

#plot the means
library(ggplot2)
positions <- c("control", "pre-exposure", "post-exposure")
ggplot(log_means) + 
  labs(y = "ln(Bd GE)") +
  geom_bar( aes(x=treatment, y=mean), stat="identity", colour="black", fill="darkorange", alpha=0.7) + scale_x_discrete(limits = positions) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  geom_errorbar( aes(x=treatment, ymin=lower_CI, ymax=upper_CI), width=0.4, colour="black", alpha=0.9, size=1.3) + theme_classic()

