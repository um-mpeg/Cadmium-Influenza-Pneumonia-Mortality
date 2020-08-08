##################################################################################
# This R code can be used to reproduce the analysis of NH99-06 data in the paper 
# "Environmental Cadmium, Mortality from Influenza and Pneumonia in U.S. Adults"
# by Sung Kyun Park and Howard Hu
# Under review at Environ Health Perspect
# R code is written by Sung Kyun Park, sungkyun@umich.edu ###
##################################################################################

library(haven)
library(survival)
library(epiDisplay)
library(Hmisc)
#### load the survey package
library(survey)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\envie\\Dropbox (University of Michigan)\\My Documents\\Cadmium\\lung\\EHP_revision_analysis")

load("cd45.99.rda")
names(cd45)
summary(cd45)
summary(ucd)

attach(cd45)
### survey
bpdsn45<-svydesign(id=~psu, strata=~strata, weights=~wtmec8yr, data=cd45, nest=T)

### incident rate
summ(permth_exm)
summary(permth_exm)
summary(permth_exm/12)
sum(permth_exm)
sum(wtmec8yr)
svymean(~permth_exm, bpdsn45)
confint(svymean(~permth_exm, bpdsn45))
sum(m_inf)
sum(m_inf[is.na(ucd3)])
1000*(sum(m_inf)/sum(permth_exm))*12
1000*(sum(m_inf*wtmec8yr)/sum(permth_exm*wtmec8yr))*12
1000*(svytotal(~m_inf, bpdsn45)/svytotal(~permth_exm, bpdsn45))*12

ci.poisson(sum(m_inf), sum(permth_exm)/(1000*12), alpha=.05)
ci.poisson(sum(m_inf*wtmec8yr), sum(permth_exm*wtmec8yr)/(1000*12), alpha=.05)

######################################################################
### TABLE 1 Characteristics ###
######################################################################
summary(bcd)
summary(bcd[cycle==1])
summary(bcd[cycle==2])
summary(bcd[cycle==3])
summary(bcd[cycle==4])

svymean(~bcd, bpdsn45)
exp(svymean(~log(bcd), bpdsn45))
exp(confint(svymean(~log(bcd), bpdsn45)))

exp(svymean(~log(scot), bpdsn45, na.rm=T))
exp(confint(svymean(~log(scot), bpdsn45, na.rm=T)))
summary(scot)

summ(cd45)
svymean(~age+bmi+chol, bpdsn45)
confint(svymean(~age+bmi+chol, bpdsn45))

svytable(~sex, Ntotal=100,design=bpdsn45)
svytable(~raceth1, Ntotal=100,design=bpdsn45)
svytable(~educ, Ntotal=100,design=bpdsn45)
svytable(~smk, Ntotal=100,design=bpdsn45)
svytable(~htn_bp, Ntotal=100,design=bpdsn45)
svytable(~pir1, Ntotal=100,design=bpdsn45)

quantile(bcd, 0.8)-quantile(bcd, 0.2)

######################################################################
### TABLE 2 Survival Analysis ###
######################################################################

#########################################
### use age.die instead of permth_exm ###
#########################################

## ALL SUBJECTS

### model 1
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     , bpdsn45)
summary(svycox.bcd)

### model 2
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.bcd)

### model 3
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+I(2-sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(smk)+log(scot)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.bcd)

#### Never smokers
nh.never<-subset(cd45, smk==0)
bpdsn.never<-svydesign(id=~psu, strata=~strata, weights=~wtmec8yr, data=nh.never, nest=T)

### model 1
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     , bpdsn.never)
summary(svycox.bcd)

### model 2
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn.never)
summary(svycox.bcd)

### model 3
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+log(scot)+factor(htn_bp)+chol+bmi, bpdsn.never)
summary(svycox.bcd)

#### Ever smokers
nh.ever<-subset(cd45, smk!=0)
bpdsn.ever<-svydesign(id=~psu, strata=~strata, weights=~wtmec8yr, data=nh.ever, nest=T)

### model 1
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     , bpdsn.ever)
summary(svycox.bcd)

### model 2
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn.ever)
summary(svycox.bcd)

### model 3
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+log(scot)+factor(htn_bp)+chol+bmi, bpdsn.ever)
summary(svycox.bcd)

##############################
### Tertiles: For Table S3 ###
##############################

tabpct(m_inf, bcd3)

#model 1
svycox.bcd3<-svycoxph(Surv(ridagemn, age.die, m_inf)~bcd3+factor(sex)+factor(raceth1)+factor(cycle)
                      , bpdsn45)
summary(svycox.bcd3)

svycox.bcd3.trend<-svycoxph(Surv(ridagemn, age.die, m_inf)~as.numeric(bcd3)+factor(sex)+factor(raceth1)+factor(cycle)
                            , bpdsn45)
summary(svycox.bcd3.trend)

### model 2
svycox.bcd3<-svycoxph(Surv(ridagemn, age.die, m_inf)~bcd3+factor(sex)+factor(raceth1)+factor(cycle)
                      +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.bcd3)

svycox.bcd3.trend<-svycoxph(Surv(ridagemn, age.die, m_inf)~as.numeric(bcd3)+factor(sex)+factor(raceth1)+factor(cycle)
                            +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.bcd3.trend)

### model 3
svycox.bcd3<-svycoxph(Surv(ridagemn, age.die, m_inf)~bcd3+factor(sex)+factor(raceth1)+factor(cycle)
                      +factor(educ)+factor(smk)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.bcd3)

svycox.bcd3.trend<-svycoxph(Surv(ridagemn, age.die, m_inf)~as.numeric(bcd3)+factor(sex)+factor(raceth1)+factor(cycle)
                            +factor(educ)+factor(smk)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.bcd3.trend)

##########################################
### Sensitivity analysis: For Table S4 ###
##########################################

### additional adjustment for pir
svycox.bcd.pir<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                         +factor(educ)+factor(htn_bp)+chol+bmi+pir, bpdsn45)
summary(svycox.bcd.pir)

### Additional adjustment for packyrs
summary(packyrs)
svycox.bcd.pack<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                          +factor(educ)+factor(smk)+log(scot)+factor(htn_bp)+chol+bmi+packyrs, bpdsn45)
summary(svycox.bcd.pack)

#########################################
### EXCLUDING TOP-CODED AGE: Table S5 ###
#########################################

### how many age 85 and older?
summ(cd45)

cd85<-subset(cd45, age<85)
summ(cd85)
## 370 subj excluded

cd85over<-subset(cd45, age>=85)

tabpct(cd45$m_inf, cd45$never, percent="col")
tabpct(cd85$m_inf, cd85$never, percent="col")
tabpct(cd85over$m_inf, cd85over$never, percent="col")


sum(cd45$m_inf) #56
sum(cd85$m_inf) #42

ci.poisson(sum(cd85$m_inf), sum(cd85$permth_exm)/(1000*12), alpha=.05)
ci.poisson(sum(cd85$m_inf*cd85$wtmec8yr), sum(cd85$permth_exm*cd85$wtmec8yr)/(1000*12), alpha=.05)

bpdsn85<-svydesign(id=~psu, strata=~strata, weights=~wtmec8yr, data=cd85, nest=T)

## ALL SUBJECTS

### model 1
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     , bpdsn85)
summary(svycox.bcd)

### model 2
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn85)
summary(svycox.bcd)

### model 3
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(smk)+log(scot)+factor(htn_bp)+chol+bmi, bpdsn85)
summary(svycox.bcd)


#### Never smokers
nh.never.85<-subset(cd85, smk==0)
bpdsn.never.85<-svydesign(id=~psu, strata=~strata, weights=~wtmec8yr, data=nh.never.85, nest=T)

### model 1
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     , bpdsn.never.85)
summary(svycox.bcd)

### model 2
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn.never.85)
summary(svycox.bcd)

### model 3
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+log(scot)+factor(htn_bp)+chol+bmi, bpdsn.never.85)
summary(svycox.bcd)

#### Ever smokers
nh.ever.85<-subset(cd85, smk!=0)
bpdsn.ever.85<-svydesign(id=~psu, strata=~strata, weights=~wtmec8yr, data=nh.ever.85, nest=T)

### model 1
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     , bpdsn.ever.85)
summary(svycox.bcd)

### model 2
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn.ever.85)
summary(svycox.bcd)

### model 3
svycox.bcd<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
                     +factor(educ)+log(scot)+factor(htn_bp)+chol+bmi, bpdsn.ever.85)
summary(svycox.bcd)


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

#save(cd45, file="cd45.99.rda")

cd45$c_other<-ifelse(cd45$c_inf==2, 1, 0)
# Create variable denoting occurrence of non-influenza death.
#################################################################################

names(cd45)
tab1(c_inf)

######################################################################
# Subdistribution hazard models and cause-specific hazard models.
######################################################################

cov.mat<-cd45[,c("bcd","age","raceth1","sex","educ","htn_bp","chol","bmi")]
cov.mat$bcd<-cov.mat$bcd/0.5
summ(cd45$bcd)
summ(cov.mat$bcd)

crr.1 <- crr(permth_exm,c_inf,cov.mat,failcode=1,cencode=0)
# Subdistribution hazard model for influenza death.

crr.2 <- crr(permth_exm,c_inf,cov.mat,failcode=2,cencode=0)
# Subdistribution hazard model for non-influenza death.

cox.1 <- coxph(Surv(ridagemn, age.die, m_inf) ~ I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
               +factor(educ)+factor(htn_bp)+chol+bmi, data=cd45)
# Cause-specific hazard for influenza death.
cox.2 <- coxph(Surv(ridagemn, age.die, c_other) ~ I(bcd/0.5)+factor(sex)+factor(raceth1)+factor(cycle)
               +factor(educ)+factor(htn_bp)+chol+bmi, data=cd45)
# Cause-specific hazard for non-influenza death.
summary(crr.1)
summary(crr.2)
summary(cox.1)
summary(cox.2)

###################################################################
##### Figure 2: Effect modification by Sex and Race/ethnicity #####
###################################################################

#using model 2
tapply(age, agecat, summary)
## instead of age 65, use age 85 as the cut-off

cd45$age85<-ifelse(cd45$age>=85, 1, 0)
tab1(cd45$age85)

#use age65.85
cd45$age65.85<-ifelse(cd45$age<65, 1, ifelse(cd45$age>=65&cd45$age<85, 2, 3))
tab1(cd45$age65.85)

cd45$age65<-ifelse(cd45$age65.85==2, 1, 0)
table(cd45$age65, cd45$age65.85)

cd45$age45<-ifelse(cd45$age65.85==1, 1, 0)
table(cd45$age45, cd45$age65.85)

### Model 2
svycox.cd.age<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(agecat)+factor(sex)+factor(raceth1)+factor(cycle)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.age)
svycox.cd.age1<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*I(3-agecat)+factor(sex)+factor(raceth1)+factor(cycle)
                         +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.age1)


svycox.cd.age<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age85+factor(sex)+factor(raceth1)+factor(cycle)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.age)
svycox.cd.age1<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*I(1-age85)+factor(sex)+factor(raceth1)+factor(cycle)
                         +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.age1)

#use age65.85
svycox.cd.cr.age45<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(age65.85)+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age45)
svycox.cd.cr.age65<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age45+I(bcd/0.5)*age85+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age65)
svycox.cd.cr.age85<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age45+I(bcd/0.5)*age65+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age85)
svycox.cd.cr.age.trend<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age65.85+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.cr.age.trend)

svycox.cd.sex<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(sex)+factor(raceth1)+factor(cycle)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.sex)
svycox.cd.sex1<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*I(2-sex)+factor(raceth1)+factor(cycle)
                         +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.sex1)

svycox.cd.race<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(raceth1)+factor(sex)+factor(cycle)
                         +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.race)

svycox.cd.nhb<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*nhw+I(bcd/0.5)*mex+I(bcd/0.5)*oth+factor(sex)+factor(cycle)
                        +factor(educ)+factor(htn_bp)+chol+bmi, bpdsn45)
summary(svycox.cd.nhb)

names(cd45)

### Model 3
svycox.cd.age<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(agecat)+factor(sex)+factor(raceth1)+factor(cycle)
                        +factor(smk)+factor(educ)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.cd.age)
svycox.cd.age1<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*I(3-agecat)+factor(sex)+factor(raceth1)+factor(cycle)
                         +factor(smk)+factor(educ)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.cd.age1)

#use age65.85
svycox.cd.cr.age45<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(age65.85)+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age45)
svycox.cd.cr.age65<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age45+I(bcd/0.5)*age85+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age65)
svycox.cd.cr.age85<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age45+I(bcd/0.5)*age65+factor(sex)+factor(raceth1)+factor(cycle)
                             +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age85)
svycox.cd.cr.age.trend<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*age65.85+factor(sex)+factor(raceth1)+factor(cycle)
                                 +factor(educ)+factor(htn_bp)+chol+bmi+factor(smk)+log(scot), bpdsn45)
summary(svycox.cd.cr.age.trend)


svycox.cd.sex<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(sex)+factor(raceth1)+factor(cycle)
                        +factor(smk)+factor(educ)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.cd.sex)
svycox.cd.sex1<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*I(2-sex)+factor(raceth1)+factor(cycle)
                         +factor(smk)+factor(educ)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.cd.sex1)

svycox.cd.race<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*factor(raceth1)+factor(sex)+factor(cycle)
                         +factor(smk)+factor(educ)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.cd.race)

svycox.cd.nhb<-svycoxph(Surv(ridagemn, age.die, m_inf)~I(bcd/0.5)*nhw+I(bcd/0.5)*mex+I(bcd/0.5)*oth+factor(sex)+factor(cycle)
                        +factor(smk)+factor(educ)+factor(htn_bp)+chol+bmi+log(scot), bpdsn45)
summary(svycox.cd.nhb)


