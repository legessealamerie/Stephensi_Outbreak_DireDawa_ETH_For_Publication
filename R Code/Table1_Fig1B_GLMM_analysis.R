############### GLMM Analysis for urban malaria Dire Dawa Ethiopia ##########
# required library
library(ggplot2)
library(tidyr)
library(dplyr)
library(MASS)
library(sjstats)
library(lme4)
library(glmulti)
library(broom)
library(tidyverse)
## importing  data 
dd=read.csv("dd_cc_individual_data.csv",header=T,na.strings="")
### to see variable names 
colnames(dd)
#creating a case/control ID: 1/0
dd$case=NA
dd$case[dd$casexp=="Case"]=1
dd$case[dd$casexp=="Control"]=0
#creating a malaria positive ID: 1/0
dd$mal=NA
dd$mal[dd$malpos=="Positive"]=1
dd$mal[dd$malpos=="Negative"]=0
dd$qpcr=NA
dd$qpcr[dd$pfqpos=="Positive"]=1
dd$qpcr[dd$pfqpos=="Negative"]=0
## tabulation of outcome variables 
table(dd$mal)
table(dd$qpcr)
#removing missing malpos data
dd2=dd[-which(dd$malpos=="NA"),]
dd$SITE=NA
dd$SITE[dd$site=="DD City"]=1
dd$SITE[dd$site=="DDU"]=0
### re-coding age category
dd %>% mutate(agecat=recode(agecat,"Under 5"="0","5 - 15 Years"="1","Above 15 Years"="2"))

##############################################################################
#############Table 1################
# Univariate Analysis
#############################################################################
# Stephensi presence
steph=glmer(mal~(1|idhh/case)+as.factor(stephpos),data=dd,family=binomial,
            nAGQ = 0)
## confidence interval
se <- sqrt(diag(vcov(steph))) # standard error 
tab <- cbind(Est = fixef(steph), LL = fixef(steph) - 1.96 * se, 
             UL = fixef(steph) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## larvae presence
larv=glmer(mal~(1|idhh/case)+as.factor(larvaepos),data=dd,family=binomial,nAGQ = 0)
summary(larv)
se <- sqrt(diag(vcov(larv)))
tab <- cbind(Est = fixef(larv), LL = fixef(larv) - 1.96 * se, UL = fixef(larv) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
##' study site  @city and @University
s=glmer(mal~(1|idhh/case)+site,data=dd,family=binomial,nAGQ = 1)
summary(s)
se <- sqrt(diag(vcov(s)))
tab <- cbind(Est = fixef(s), LL = fixef(s) - 1.96 * se, UL = fixef(s) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## site and sex interaction 
s_sex=glmer(mal~(1|idhh/case)+site*sex,data=dd,family=binomial,nAGQ = 0)
## summary result
summary(s_sex)
se <- sqrt(diag(vcov(s_sex)))
tab <- cbind(Est = fixef(s_sex), LL = fixef(s_sex) - 1.96 * se, UL = fixef(s_sex) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## site and age 
s_age=glmer(mal~(1|idhh/case)+site*agecat,data=dd,family=binomial,nAGQ = 0)
summary(s_age)
se <- sqrt(diag(vcov(s_age)))
tab <- cbind(Est = fixef(s_age), LL = fixef(s_age) - 1.96 * se, UL = fixef(s_age) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)

## sex
sex=glmer(mal~(1|idhh/case/SITE)+as.factor(sex),data=dd,family=binomial,nAGQ = 0)
summary(sex)
se <- sqrt(diag(vcov(sex)))
tab <- cbind(Est = fixef(sex), LL = fixef(sex) - 1.96 * se, UL = fixef(sex) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## Age 
age= glmer(mal~(1|idhh/case) +agecat,data=dd,family=binomial,nAGQ = 0)
summary(age)
se <- sqrt(diag(vcov(age)))
tab <- cbind(Est = fixef(age), LL = fixef(age) - 1.96 * se, UL = fixef(age) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## adult presence
adult= glmer(mal~(1|idhh/case) +adultpos_2,data=dd,family=binomial,nAGQ = 0)
summary(adult)
se <- sqrt(diag(vcov(adult)))
tab <- cbind(Est = fixef(adult), LL = fixef(adult) - 1.96 * se, UL = fixef(adult) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## water body presence in the neighborhood
water= glmer(mal~(1|idhh/case) +as.factor(waterbody),data=dd,family=binomial,nAGQ = 0)
summary(water)
se <- sqrt(diag(vcov(water)))
tab <- cbind(Est = fixef(water), LL = fixef(water) - 1.96 * se, UL = fixef(water) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
## use of spray 
spr= glmer(mal~(1|idhh/case) +as.factor(spray),data=dd,
           family=binomial,nAGQ = 0)
summary(spr)
se <- sqrt(diag(vcov(spr)))
tab <- cbind(Est = fixef(spr), LL = fixef(spr) - 1.96 * se,
             UL = fixef(spr)  + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
##############################################################################
## Multivariate Analysis 
## individual and household level factors 
ind_hh=glmer(mal~(1|idhh/case) +site+ as.factor(stephpos)+site*sex+site*agecat  +as.factor(waterbody)+as.factor(spray),data=dd,family=binomial,nAGQ = 1)
summary(ind_hh)
performance::icc(ind_hh)
se <- sqrt(diag(vcov(ind_hh)))
tab <- cbind(OR = fixef(ind_hh), LL = fixef(ind_hh) - 1.96 * se, UL = fixef(ind_hh) + 1.96 * se)
## odds ratios with 95% CI
print(exp(tab), digits=3)
##############################################################################

### forest plot 
#Required library
library(devtools)
devtools::install_github("strengejacke/sjPlot")
require(sjPlot)
set_theme(base = theme_classic())
plot_model(ind_hh,rm.terms=c("as.factor(waterbody)Yes","as.factor(agecat)Above 15 Years"),show.values = T,title = "Forest plot the full model")
plot_model(ind_hh,rm.terms=c("as.factor(waterbody)Yes","as.factor(agecat)Above 15 Years"),title = "Forest plot the full model")

## performance of the model 
performance::r2(ind_hh)
###############################################################################################
#Extracting the effect of random effects (so the background OR not explained by the fixed effects)
randoms=ranef(ind_hh,condVar=TRUE)[[1]]
variances <- as.numeric(attr(randoms, "postVar"))
res=data.frame(WES = rownames(randoms), mean_effect = randoms$`(Intercept)`+sum(coef(summary(ind_hh))[,"Estimate"]))
res$lower <- res$mean_effect - 2* sqrt(variances)
res$upper <- res$mean_effect + 2* sqrt(variances)
res$mean_effect <- exp(res$mean_effect)
res$lower <- exp(res$lower)
res$upper <- exp(res$upper)
res$WES <- reorder(res$WES, res$mean_effect, mean)
require(ggplot2)
ggplot(data=res,aes(x=WES,y=mean_effect))+geom_point(col="blue")+geom_errorbar(width=1,aes(ymin=lower,ymax=upper),
                                                                               col="black")+labs(title="Malaria Positiviy",y="Odds Ratio",x=NULL)

res$Case=sub('*:DD City','',as.character(res$WES))
res$Case=sub('*:DDU','',as.character(res$Case))
res$Case=substr(res$Case,nchar(res$Case),nchar(res$Case))

#ODDS between cases and controls (as random effects)
apply(res[res$Case==1,2:4],2,quantile,prob=0.5)
apply(res[res$Case==0,2:4],2,quantile,prob=0.5)

#ODDS between sites (as random effects)
apply(res[1:157,2:4],2,quantile,prob=0.5)#City
apply(res[158:290,2:4],2,quantile,prob=0.5)#DDU

#Finally odds city & cases VS city & controls (as random effects)
x=which(res$Case[1:157]==1)
apply(res[x,2:4],2,quantile,prob=0.5)
x=which(res$Case[1:157]==0)
apply(res[x,2:4],2,quantile,prob=0.5)

#Finally odds DDU & cases VS DDU & controls (as random effects)
x=which(res$Case[158:290]==1)
apply(res[c(158:290)[x],2:4],2,quantile,prob=0.5)
x=which(res$Case[158:290]==0)
apply(res[c(158:290)[x],2:4],2,quantile,prob=0.5)
###############################################################################
######## Figure 1B ########
## Taking qpcr as outcome variable
## City and University together taking only index family and control family
### frequency table of case category 
table(dd$casecat)
### importing only index and control family
case_cat<-read.csv("Case_category_data.csv",header=T,na.strings = "")

#creating a malaria positive ID microscopy/RDT: 1/0
case_cat$mal=NA
case_cat$mal[case_cat$malpos=="Positive"]=1
case_cat$mal[case_cat$malpos=="Negative"]=0
#creating a qPCR positive  ID : 1/0
case_cat$qpcr=NA
case_cat$qpcr[case_cat$pfqpos=="Positive"]=1
case_cat$qpcr[case_cat$pfqpos=="Negative"]=0
## tabulation of outcome variables 
table(case_cat$mal)
table(case_cat$qpcr)
### qPCR result as an outcome
all_qpcr=glmer(qpcr~(1|idhh)+factor(casecat),data=case_cat,family=binomial,nAGQ = 1)
summary(all_qpcr)
##microscopy/rdt result as an outcome 
all_mic_rdt=glmer(mal~(1|idhh)+factor(casecat),data=case_cat,family=binomial,nAGQ = 1)
summary(all_mic_rdt)
###############################################################################
table(case_cat$site)
## recoding site 
case_cat$SITE=NA
case_cat$SITE[case_cat$site=="DD City"]=1
case_cat$SITE[case_cat$site=="DDU"]=0

# filter only city 
city<-filter(case_cat,SITE==1)
# qPCR result as an outcome
city_qpcr=glmer(qpcr~(1|idhh)+factor(casecat),data=city,family=binomial,nAGQ = 1)
summary(city_qpcr)
##microscopy/rdt result as an outcome 
city_mic_rdt=glmer(mal~(1|idhh)+factor(casecat),data=city,family=binomial,nAGQ = 1)
summary(city_mic_rdt)
# University site 
univ<-filter(case_cat,SITE==0) 
univ_qpcr=glmer(qpcr~(1|idhh)+factor(casecat),data=univ,family=binomial,nAGQ = 0)
summary(univ_qpcr)
##microscopy/rdt result as an outcome 
univ_mic_rdt=glmer(mal~(1|idhh)+factor(casecat),data=univ,family=binomial,nAGQ = 1)
summary(city_mic_rdt)

##############################################################################