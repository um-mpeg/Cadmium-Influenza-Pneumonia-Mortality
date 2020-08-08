##################################################################################
# This R code can be used to reproduce the paper 
# "Environmental Cadmium, Mortality from Influenza and Pneumonia in U.S. Adults"
# by Sung Kyun Park and Howard Hu
# Under review at Environ Health Perspect
# R code is written by Sung Kyun Park, sungkyun@umich.edu ###
##################################################################################

##################################################################
# Kaplan-Meier curves with IPW: Figure 1
# Data from NH-III and NH99-06
##################################################################
library(survival)
library(epiDisplay)
library(Hmisc)
library(ggplot2)
library(tidyverse)
#### load the survey package
library(survey)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\envie\\Dropbox (University of Michigan)\\My Documents\\Cadmium\\lung\\EHP_revision_analysis")

load("cd45.nh3.rda")
names(cd45)

attach(cd45)
summary(cd45[which(m_inf==1),]$PERMTH_EXM)

#install.packages("survival")
#install.packages("ggplot2") # for plots
#install.packages("survminer") # for plots
library("survival")
library("ggplot2")
library("survminer")
library(grid)

## check survtime
summ(PERMTH_EXM)
summ(age.die)

### Create weight for low vs high Cd
table(cd.cr2, cd.cr2n)
cd.wt <- glm(cd.cr2n ~ age+factor(sex)+factor(raceth)
             +factor(educ)+factor(htn_bp)+chol+bmi+factor(phase), 
             family = binomial(), data = cd45)
summary(cd.wt)

p.cd.obs <- ifelse(cd.cr2n == 0, 1 - predict(cd.wt, type = "response"),
                   predict(cd.wt, type = "response"))

## pred<-predict(cd.wt, type = "response")
## length(pred)

cd45$w <- 1/p.cd.obs
summary(cd45$w)
sd(cd45$w)
summ(cd45$w)
sum(cd45$w)
tapply(cd45$w, cd.cr2, sum)

## stabilized w
pd.cd<-predict(cd.wt, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(cd.cr2n ~1, family = binomial(), data=cd45)
summary(numer.fit)
pn.cd <- predict(numer.fit, type = "response")

cd45$sw <- ifelse(cd.cr2n == 0, ((1-pn.cd)/(1-pd.cd)),
                       (pn.cd/pd.cd))
summary(cd45$sw)
tapply(cd45$sw, cd.cr2, sum)

### weighted survfit
### using age as the time scale
summary(age.die/12)
tapply(age.die/12, m_inf, max)

fit1 <- survfit(Surv(HSAITMOR/12, age.die/12, m_inf) ~ cd.cr2, weights=sw, data=cd45)
summary(fit1)

km.nh3.age<-ggsurvplot(fit1, data = cd45, xlab="Age", xlim=c(45,105),
                   conf.int=F, linetype = "strata",
                   legend.labs = c("Low", "High"),
                   break.time.by=10, ggtheme = theme_bw(), censor=F, fun="cumhaz",
                   title="A. NHANES-III: Creatinine-corrected urinary cadmium", risk.table = TRUE)
km.nh3.age

########## NH 1999-2006 ##########

load("C:\\Users\\envie\\Dropbox (University of Michigan)\\My Documents\\Cadmium\\lung\\EHP_revision_analysis\\cd45.99a.rda")
names(cd45.99)

### Create weight for low vs high Cd
cd.wt.99 <- glm(bcd2n ~ age+factor(sex)+factor(raceth1)
             +factor(educ)+factor(htn_bp)+chol+bmi+factor(cycle), 
             family = binomial(), data = cd45.99)
summary(cd.wt.99)

p.cd.obs.99 <- ifelse(cd45.99$bcd2n == 0, 1 - predict(cd.wt.99, type = "response"),
                   predict(cd.wt.99, type = "response"))

cd45.99$w <- 1/p.cd.obs.99
summary(cd45.99$w)
sd(cd45.99$w)
summ(cd45.99$w)

## stabilized w
pd.cd.99<-predict(cd.wt.99, type = "response")

# estimation of numerator of ip weights
numer.fit.99 <- glm(bcd2n ~1, family = binomial(), data=cd45.99)
summary(numer.fit.99)
pn.cd.99 <- predict(numer.fit.99, type = "response")

cd45.99$sw <- ifelse(cd45.99$bcd2n == 0, ((1-pn.cd.99)/(1-pd.cd.99)),
                  (pn.cd.99/pd.cd.99))
summary(cd45.99$sw)
tapply(cd45.99$sw, cd45.99$bcd2, sum)

### weighted survfit
### using age as the time scale

summ(cd45.99$ridagemn/12)
summ(cd45.99$age.die/12)
tapply(cd45.99$age.die/12, cd45.99$m_inf, max)

fit2 <- survfit(Surv(ridagemn/12, age.die/12, m_inf) ~ bcd2, weights=sw, data=cd45.99)
summary(fit2)

km.nh99.age<-ggsurvplot(fit2, data = cd45.99, xlab="Age", xlim=c(45,105),
                       conf.int=F, linetype = "strata",
                       legend.labs = c("Low", "High"),
                       break.time.by=10, ggtheme = theme_bw(), censor=F, fun="cumhaz",
                       title="B. NHANES 1999-2006: Blood cadmium", risk.table = TRUE)
km.nh99.age

