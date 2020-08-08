##################################################################################
# This R code can be used to reproduce the analysis of NH-III data in the paper 
# "Environmental Cadmium, Mortality from Influenza and Pneumonia in U.S. Adults"
# by Sung Kyun Park and Howard Hu
# Under review at Environ Health Perspect
# R code is written by Sung Kyun Park, sungkyun@umich.edu ###
##################################################################################

library(haven)
library(survival)
library(epiDisplay)
library(Hmisc)
library(ggplot2)
library(survey)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\envie\\Dropbox (University of Michigan)\\My Documents\\Cadmium\\lung\\EHP_revision_analysis")

load("cd45.nh3.rda")
names(cd45)
summary(cd45)
summ(cd45)
table(phase)

attach(cd45)
bpdsn45<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=cd45, nest=T)

### incident rate
summ(PERMTH_EXM)
summary(PERMTH_EXM)
summary(PERMTH_EXM/12)
sum(PERMTH_EXM)
sum(wt_mh)
svymean(~PERMTH_EXM, bpdsn45)
confint(svymean(~PERMTH_EXM, bpdsn45))
sum(m_inf)
1000*(sum(m_inf)/sum(PERMTH_EXM))*12
1000*(sum(m_inf*wt_mh)/sum(PERMTH_EXM*wt_mh))*12
1000*(svytotal(~m_inf, bpdsn45)/svytotal(~PERMTH_EXM, bpdsn45))*12

ci.poisson(sum(m_inf), sum(PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(m_inf*wt_mh), sum(PERMTH_EXM*wt_mh)/(1000*12), alpha=.05)

######################################################################
### TABLE 1 Characteristics ###
######################################################################

summary(ucd)
### DL=0.01, how many?
lod<-ifelse(ucd==0.01, 1, 0)
tab1(lod)
# 239/7173 (3.3%)

svymean(~ucd, bpdsn45)
svymean(~cd.cr, bpdsn45)
exp(svymean(~log(ucd), bpdsn45))
exp(confint(svymean(~log(ucd), bpdsn45)))
exp(svymean(~log(cd.cr), bpdsn45))
exp(confint(svymean(~log(cd.cr), bpdsn45)))

exp(svymean(~log(ucreat), bpdsn45))
exp(confint(svymean(~log(ucreat), bpdsn45)))

exp(svymean(~log(scot), bpdsn45, na.rm=T))
exp(confint(svymean(~log(scot), bpdsn45, na.rm=T)))
summary(scot)


summ(cd45)
svymean(~age+bmi+chol, bpdsn45)
confint(svymean(~age+bmi+chol, bpdsn45))

svytable(~sex, Ntotal=100,design=bpdsn45)
svytable(~raceth, Ntotal=100,design=bpdsn45)
svytable(~educ, Ntotal=100,design=bpdsn45)
svytable(~smk, Ntotal=100,design=bpdsn45)
svytable(~htn_bp, Ntotal=100,design=bpdsn45)
svytable(~pir1, Ntotal=100,design=bpdsn45)

svytable(~sex+m_inf, Ntotal=100, design=bpdsn45)
tabpct(sex, m_inf)

## check weighted distribution
#svyby(~ucd, ~smk, bpdsn45, svyquantile, quantiles=c(0.25,0.5,0.75))
svyquantile(~ucd, bpdsn45, c(.25,.5,.75))
svyquantile(~cd.cr, bpdsn45, c(.25,.5,.75))
svyquantile(~cd.cr, bpdsn45, c(1/3,2/3))
quantile(cd.cr, c(1/3,2/3))

######################################################################
### TABLE 2 Survival Analysis ###
######################################################################

#########################################
### use age.die instead of PERMTH_EXM ###
#########################################

## ALL SUBJECTS

#model 1
svycox.cd.cr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                       , bpdsn45)
summary(svycox.cd.cr)

#model 2
svycox.cd.cr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                       +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr)

#model 3 
svycox.cd.cr.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+I(2-sex)+factor(raceth)+factor(phase)
                             +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.lscot)

#### Never smokers
nh3.never.45<-subset(cd45, smk==1)

bpdsn.never.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=nh3.never.45, nest=T)

#Model 1
svycox.cd.cr.never<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                             , bpdsn.never.45)
summary(svycox.cd.cr.never)

#Model 2
svycox.cd.cr.never<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                             +factor(sex)+factor(htn_bp)+chol+bmi, bpdsn.never.45)
summary(svycox.cd.cr.never)

### Model 3: additional adjustment for scot
svycox.cd.cr.never.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                                   +factor(sex)+factor(htn_bp)+chol+bmi+log(scot), bpdsn.never.45)
summary(svycox.cd.cr.never.lscot)


#### Ever smokers
nh3.ever.45<-subset(cd45, smk!=1)

bpdsn.ever.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=nh3.ever.45, nest=T)

#Model 1
svycox.cd.cr.ever<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                             , bpdsn.ever.45)
summary(svycox.cd.cr.ever)

#Model 2
svycox.cd.cr.ever<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                             +factor(sex)+factor(htn_bp)+chol+bmi, bpdsn.ever.45)
summary(svycox.cd.cr.ever)

### Model 3: additional adjustment for scot
svycox.cd.cr.ever.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                                   +factor(sex)+factor(htn_bp)+chol+bmi+log(scot), bpdsn.ever.45)
summary(svycox.cd.cr.ever.lscot)

#################################
### Tertiles: For Table S3 ######
#################################

tabpct(m_inf, cd.cr3)

#model 1
svycox.cd.cr3<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd.cr3+age+factor(sex)+factor(raceth)+factor(phase)
                        , bpdsn45)
summary(svycox.cd.cr3)

svycox.cd.cr3.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~as.numeric(cd.cr3)+age+factor(sex)+factor(raceth)+factor(phase)
                              , bpdsn45)
summary(svycox.cd.cr3.trend)

#model 2
svycox.cd.cr3<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd.cr3+age+factor(sex)+factor(raceth)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr3)

svycox.cd.cr3.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~as.numeric(cd.cr3)+age+factor(sex)+factor(raceth)+factor(phase)
                              +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr3.trend)

#model 3 (adj for smk and serum cotinine)
svycox.cd.cr3.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd.cr3+age+factor(sex)+factor(raceth)+factor(phase)
                              +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr3.lscot)

svycox.cd.cr3.lscot.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~as.numeric(cd.cr3)+age+factor(sex)+factor(raceth)+factor(phase)
                                    +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr3.lscot.trend)

##########################################
### Sensitivity analysis: For Table S4 ###
##########################################

### additional adjustment for pir
svycox.cd.cr.pir1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                            +pir+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.pir1)

## additional adjustment for iron
svycox.cd.cr.iron<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                            +factor(sex)+factor(htn_bp)+chol+bmi+scot+iron, bpdsn45)
summary(svycox.cd.cr.iron)

svycox.cd.cr.tib<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                           +factor(sex)+factor(htn_bp)+chol+bmi+scot+tib, bpdsn45)
summary(svycox.cd.cr.tib)

svycox.cd.cr.pct<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                           +factor(sex)+factor(htn_bp)+chol+bmi+scot+pct, bpdsn45)
summary(svycox.cd.cr.pct)

svycox.cd.cr.fer<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                           +factor(sex)+factor(htn_bp)+chol+bmi+scot+fer, bpdsn45)
summary(svycox.cd.cr.fer)

### Additional adjustment for packyrs
summary(packyrs)
svycox.cd.cr.pack1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                             +factor(smk)+log(scot)+packyrs+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.pack1)

#########################################
### EXCLUDING TOP-CODED AGE: Table S5 ###
#########################################

### how many age 90 and older?

cd90<-subset(cd45, age<90)
summary(cd90)
summ(cd90)
## 69 subj excluded

cd90over<-subset(cd45, age>=90)

tabpct(cd45$m_inf, cd45$never, percent="col")
tabpct(cd90$m_inf, cd90$never, percent="col")
tabpct(cd90over$m_inf, cd90over$never, percent="col")

sum(cd45$m_inf) #141
sum(cd90$m_inf) #137

ci.poisson(sum(cd90$m_inf), sum(cd90$PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(cd90$m_inf*cd90$wt_mh), sum(cd90$PERMTH_EXM*cd90$wt_mh)/(1000*12), alpha=.05)

bpdsn90<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=cd90, nest=T)

#model 1
svycox.cd.cr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                       , bpdsn90)
summary(svycox.cd.cr)

#model 2
svycox.cd.cr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                       +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn90)
summary(svycox.cd.cr)

#model 3 
svycox.cd.cr.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+I(2-sex)+factor(raceth)+factor(phase)
                             +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn90)
summary(svycox.cd.cr.lscot)

#### Never smokers
nh3.never.90<-subset(cd90, smk==1)

bpdsn.never.90<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=nh3.never.90, nest=T)

#Model 1
svycox.cd.cr.never<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                           , bpdsn.never.90)
summary(svycox.cd.cr.never)

#Model 2
svycox.cd.cr.never<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                             +factor(sex)+factor(htn_bp)+chol+bmi, bpdsn.never.90)
summary(svycox.cd.cr.never)

### Model 3: additional adjustment for scot
svycox.cd.cr.never.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                                   +factor(sex)+factor(htn_bp)+chol+bmi+log(scot), bpdsn.never.90)
summary(svycox.cd.cr.never.lscot)

#### Ever smokers
nh3.ever.90<-subset(cd90, smk!=1)

bpdsn.ever.90<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=nh3.ever.90, nest=T)

#Model 1
svycox.cd.cr.ever<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
                             , bpdsn.ever.90)
summary(svycox.cd.cr.ever)

#Model 2
svycox.cd.cr.ever<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                             +factor(sex)+factor(htn_bp)+chol+bmi, bpdsn.ever.90)
summary(svycox.cd.cr.ever)

### Model 3: additional adjustment for scot
svycox.cd.cr.ever.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+factor(educ)+factor(raceth)+factor(phase)
                                   +factor(sex)+factor(htn_bp)+chol+bmi+log(scot), bpdsn.ever.90)
summary(svycox.cd.cr.ever.lscot)


###################################################
### Table S1: cov-adj standardization + cov adj ###
###################################################

#model 1
svycox.cd.ccr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.ccr/0.75)+age+factor(sex)+factor(raceth)+factor(phase)
                       , bpdsn45)
summary(svycox.cd.ccr)

#model 2
svycox.cd.ccr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.ccr/0.75)+age+factor(sex)+factor(raceth)+factor(phase)
                       +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.ccr)

#model 3 
svycox.cd.ccr.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.ccr/0.75)+age+I(2-sex)+factor(raceth)+factor(phase)
                             +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.ccr.lscot)

### Tertiles
#model 1
svycox.cd3<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd3+age+factor(sex)+factor(raceth)+factor(phase)
                        , bpdsn45)
summary(svycox.cd3)

svycox.cd3.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~as.numeric(cd3)+age+factor(sex)+factor(raceth)+factor(phase)
                              , bpdsn45)
summary(svycox.cd3.trend)

#model 2
svycox.cd3<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd3+age+factor(sex)+factor(raceth)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd3)

svycox.cd3.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~as.numeric(cd3)+age+factor(sex)+factor(raceth)+factor(phase)
                              +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd3.trend)

#model 3 (adj for smk and serum cotinine)
svycox.cd3.lscot<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd3+age+factor(sex)+factor(raceth)+factor(phase)
                              +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd3.lscot)

svycox.cd3.lscot.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~as.numeric(cd3)+age+factor(sex)+factor(raceth)+factor(phase)
                                    +factor(smk)+log(scot)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd3.lscot.trend)


#######################################################
### Table S2: Interaction b/w Cd and smoking status ###
#######################################################

tabpct(lcdns, m_inf, percent="row")
tabpct(hcdns, m_inf, percent="row")
tabpct(lcdes, m_inf, percent="row")
tabpct(hcdes, m_inf, percent="row")

svycox.cd.smk<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~hcdns+lcdes+hcdes+age+I(2-sex)+factor(raceth)+factor(phase)
                         +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.smk)

svycox.cd.smk1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~cd.cr2*I(1-never)+age+I(2-sex)+factor(raceth)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.smk1)

##################################################################
# Parametric estimation of survival curves via hazards model
##################################################################

# creation of person-month data

###############################################################################
#install.packages("splitstackshape")
library("splitstackshape")

### need to increase memory
memory.size()
memory.limit(20000)

cd45.surv <- expandRows(cd45, "PERMTH_EXM", drop=F) 
cd45.surv$time <- sequence(rle(cd45.surv$SEQN)$lengths)-1
cd45.surv$age.mon <- cd45.surv$time+cd45.surv$HSAITMOR 
cd45.surv$event <- ifelse(cd45.surv$age.mon==cd45.surv$age.die-1 & 
                            cd45.surv$m_inf==1, 1, 0)
cd45.surv$age.mon.sq <- cd45.surv$age.mon^2

summ(cd45.surv$age.mon)

save(cd45.surv, file="cd45.forlogistic.rda")

## DO NOT RUN
###############################################################################
load("cd45.forlogistic.rda")
names(cd45.surv)

# fit of parametric hazards model

## using svyglm
bpdsn45.surv<-svydesign(id=~psu, strata=~strata, weights=~wt_mh, data=cd45.surv, nest=T)

svyglm.cd.cr <- svyglm(event~I(cd.cr/0.86)+I(2-sex)+factor(raceth)+factor(phase)
                       +factor(educ)+factor(htn_bp)+chol+bmi
                       +age.mon + age.mon.sq, family=quasibinomial(), bpdsn45.surv)
summary(svyglm.cd.cr)
logistic.display(svyglm.cd.cr)

#model 2 (for comparison, this is in Table 2)
svycox.cd.cr<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)+I(2-sex)+factor(raceth)+factor(phase)+factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr)

### Interaction between Cd and smoking ###
#using svycox
#svycox.cd.smk<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~hcdns+lcdes+hcdes+age+I(2-sex)+factor(raceth)+factor(phase)
#                         +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)

svyglm.cd.smk <- svyglm(event~hcdns+lcdes+hcdes+age+I(2-sex)+factor(raceth)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi
                        +age.mon + age.mon.sq, family=quasibinomial(), bpdsn45.surv)
summary(svyglm.cd.smk)
logistic.display(svyglm.cd.smk)
summary(svyglm.cd.smk)$cov.unscaled

# using interaction
table(never)

svyglm.cd.smk1 <- svyglm(event~cd.cr2*I(1-never)+age+I(2-sex)+factor(raceth)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi
                        +age.mon + age.mon.sq, family=quasibinomial(), bpdsn45.surv)
summary(svyglm.cd.smk1)
logistic.display(svyglm.cd.smk1)



####################################################################
### Table S6. Competing Risk: Subdistribution and Cause-specific ###
####################################################################

### Load the library 'cmprsk'
library(cmprsk)

##################################################################################
##########  DO NOT RUN, ALREADY RUN, MOVE TO THE END OF ##########
table(cd45$m_total)
table(cd45$m_total, cd45$m_inf)

cd45$c_inf<-ifelse(cd45$m_inf==1, 1, ifelse(cd45$m_inf==0&cd45$m_total==1, 2, 0))
table(cd45$c_inf)
# c_inf is the influenza death indicator:
# 1: influenza/pneumonia death.
# 2: Non-influenza death.
# 0: Censored observation: alive at end of follow-up.

#save(cd45, file="cd45.nh3.rda")

cd45$c_other<-ifelse(cd45$c_inf==2, 1, 0)
# Create variable denoting occurrence of non-influenza death.
#################################################################################

names(cd45)
table(c_inf)

######################################################################
# Subdistribution hazard models and cause-specific hazard models.
######################################################################

cov.mat<-cd45[,c("cd.cr","age","raceth","sex","educ","htn_bp","chol","bmi")]
cov.mat$cd.cr<-cov.mat$cd.cr/0.86

crr.1 <- crr(PERMTH_EXM,c_inf,cov.mat,failcode=1,cencode=0)
# Subdistribution hazard model for influenza death.

crr.2 <- crr(PERMTH_EXM,c_inf,cov.mat,failcode=2,cencode=0)
# Subdistribution hazard model for non-influenza death.

cox.1 <- coxph(Surv(HSAITMOR, age.die, m_inf) ~ I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
               +factor(educ)+factor(htn_bp)+chol+bmi, data=cd45)
# Cause-specific hazard for influenza death.
cox.2 <- coxph(Surv(HSAITMOR, age.die, c_other) ~ I(cd.cr/0.86)+factor(sex)+factor(raceth)+factor(phase)
               +factor(educ)+factor(htn_bp)+chol+bmi, data=cd45)
# Cause-specific hazard for non-influenza death.
summary(crr.1)
summary(crr.2)
summary(cox.1)
summary(cox.2)


########################################################################
##### Figure 2: Effect modification by Age, Sex and Race/ethnicity #####
##### Also for Table S8 #####
########################################################################

#using model 2
tapply(age, agecat, summary)
## instead of age 65, use age 90 as the cut-off

cd45$age90<-ifelse(cd45$age>=90, 1, 0)
tab1(cd45$age90)

svycox.cd.cr.age<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(agecat)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age)
svycox.cd.cr.age1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*I(3-agecat)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age1)

svycox.cd.cr.age<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age90+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age)
svycox.cd.cr.age1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*I(1-age90)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age1)

#use age65.90
cd45$age65.90<-ifelse(cd45$age<65, 1, ifelse(cd45$age>=65&cd45$age<90, 2, 3))
tab1(cd45$age65.90)

cd45$age65<-ifelse(cd45$age65.90==2, 1, 0)
table(cd45$age65, cd45$age65.90)

cd45$age45<-ifelse(cd45$age65.90==1, 1, 0)
table(cd45$age45, cd45$age65.90)

svycox.cd.cr.age45<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(age65.90)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age45)
svycox.cd.cr.age65<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age45+I(cd.cr/0.86)*age90+factor(sex)+factor(raceth)+factor(phase)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age65)
svycox.cd.cr.age90<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age45+I(cd.cr/0.86)*age65+factor(sex)+factor(raceth)+factor(phase)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age90)
svycox.cd.cr.age.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age65.90+factor(sex)+factor(raceth)+factor(phase)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age.trend)


svycox.cd.cr.sex<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(sex)+factor(raceth)+factor(phase)
                       +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.sex)
svycox.cd.cr.sex1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*I(2-sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.sex1)

svycox.cd.cr.race<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(raceth)+factor(sex)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.race)
svycox.cd.nhb<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*nhw+I(cd.cr/0.86)*mex+I(cd.cr/0.86)*oth+factor(sex)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.nhb)

#using model 3
tapply(age, agecat, summary)

svycox.cd.cr.age<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(agecat)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age)
svycox.cd.cr.age1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*I(3-agecat)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age1)

svycox.cd.cr.age45<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(age65.90)+factor(sex)+factor(raceth)+factor(phase)
                             +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age45)
svycox.cd.cr.age65<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age45+I(cd.cr/0.86)*age90+factor(sex)+factor(raceth)+factor(phase)
                             +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age65)
svycox.cd.cr.age90<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age45+I(cd.cr/0.86)*age65+factor(sex)+factor(raceth)+factor(phase)
                             +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age90)
svycox.cd.cr.age.trend<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*age65.90+factor(sex)+factor(raceth)+factor(phase)
                                 +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age.trend)

svycox.cd.cr.sex<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.sex)
svycox.cd.cr.sex1<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*I(2-sex)+factor(raceth)+factor(phase)
                            +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.sex1)

svycox.cd.cr.race<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*factor(raceth)+factor(sex)+factor(phase)
                            +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.race)
svycox.cd.nhb<-svycoxph(Surv(HSAITMOR, age.die, m_inf)~I(cd.cr/0.86)*nhw+I(cd.cr/0.86)*mex+I(cd.cr/0.86)*oth+factor(sex)+factor(phase)
                        +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.nhb)

