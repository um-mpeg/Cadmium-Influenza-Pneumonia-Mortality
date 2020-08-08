##################################################################################
# This R code can be used to reproduce the paper 
# "Environmental Cadmium, Mortality from Influenza and Pneumonia in U.S. Adults"
# by Sung Kyun Park and Howard Hu
# Under review at Environ Health Perspect
# R code is written by Sung Kyun Park, sungkyun@umich.edu ###
##################################################################################

###################################################################
##### Figure 2: Effect modification by Sex and Race/ethnicity #####
###################################################################

library(ggplot2)
library(grid)

### NHANES-III

modif<-c("Age","Age","Age","Sex","Sex","Race/Ethnicity","Race/Ethnicity")
group <- c("45-64","65-89","90+","Male","Female","NHW","NHB")
HR  <- c(1.1232,1.1842,1.6257,1.1357,1.1784,1.2422,1.2208) 
lower <- c(1.0060,1.0696,0.8863,1.0455,1.0271,1.0147,1.1380)
upper <- c(1.254,1.311,2.982,1.2336,1.352,1.521,1.310)


class(modif)
class(HR)

modif1<-factor(modif, levels=unique(modif))
group1<-factor(group, levels=unique(group))

RR_data <- data.frame(modif1, group1, HR, lower, upper)

p2 = ggplot(data=RR_data,
            aes(x = group1,y = HR, ymin = lower, ymax = upper))+
  geom_pointrange(aes(col=modif1))+
  geom_hline(aes(),yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% CI)")+
  geom_errorbar(aes(ymin=lower, ymax=upper, col=modif1),width=0.2,cex=1)+ 
  facet_wrap(~modif1,strip.position="top",nrow=1,scales = "free_x") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        legend.position="none")+
  scale_y_log10(breaks=c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,3.0))+guides(fill=FALSE)+
  ggtitle("NHANES-III")
p2


### Continuous NHANES

group99 <- c("45-64","65-84","85+","Male","Female","NHW","NHB")
group99.1<-factor(group99, levels=unique(group99))
HR2  <- c(0.9586,1.2732,1.2669,1.1478,1.1320,1.1342,1.3843) 
lower2 <- c(0.59402,1.0605,0.6603,0.8877,0.8901,0.9411,1.0308)
upper2 <- c(1.5469,1.529,2.431,1.484,1.439,1.3670,1.8590)

RR_data1 <- data.frame(modif1, group99.1, HR2, lower2, upper2)

p1 = ggplot(data=RR_data1,
            aes(x = group99.1,y = HR2, ymin = lower2, ymax = upper2))+
  geom_pointrange(aes(col=modif1))+
  geom_hline(aes(),yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% CI)")+
  geom_errorbar(aes(ymin=lower2, ymax=upper2, col=modif1),width=0.2,cex=1)+ 
  facet_wrap(~modif1,strip.position="top",nrow=1,scales = "free_x") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        legend.position="none")+
  scale_y_log10(breaks=c(0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5))+guides(fill=FALSE)+
  ggtitle("NHANES 1999-2006")
p1

grid.newpage()

## Next, use the pushViewport() function to define the various frames (viewports)
## in the grid graphic system

pushViewport(viewport(layout=grid.layout(2,1)))

## and then use the print function to print the objects into the viewport.

print(p2, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(p1, vp=viewport(layout.pos.row=2, layout.pos.col=1))


### meta-analysis
library(rmeta)
cycle<-c("NHANES-III", "NHANES 1999-2006")

# Model 2
coef<-c(0.140605,0.130216)
se<-c(0.044990,0.089577)

model<-meta.summaries(coef, se, names=cycle, method=c("fixed","random"), logscale=T)
#method="random" estimates and adds a heterogeneity variance
summary(model)
plot(model)

# Model 3
coef<-c(0.088469,0.133981)
se<-c(0.075128,0.100181)
model<-meta.summaries(coef, se, names=cycle, method=c("fixed","random"), logscale=T)
#method="random" estimates and adds a heterogeneity variance
summary(model)

# Never smokers
coef<-c(0.2357824,0.5372925)
se<-c(0.0637788,0.3015130)
model<-meta.summaries(coef, se, names=cycle, method=c("fixed","random"), logscale=T)
#method="random" estimates and adds a heterogeneity variance
summary(model)
