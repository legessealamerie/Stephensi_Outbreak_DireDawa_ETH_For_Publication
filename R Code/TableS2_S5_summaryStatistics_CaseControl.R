
###################### Summary Statistics Table S2#############################
# required library
library(tidyr)
library(dplyr)
library(tidyverse)
dd=read.csv("dd_cc_individual_data.csv",header=T,na.strings="")
dd_hh=read.csv("HH_level_entomology_data.csv",header=T,na.strings="")
colnames(dd)
#creating a case/control ID: 1/0
dd$case=NA
dd$case[dd$casexp=="Case"]=1
dd$case[dd$casexp=="Control"]=0
#creating a malaria positive ID: 1/0
dd$mal=NA
dd$mal[dd$malpos=="Positive"]=1
dd$mal[dd$malpos=="Negative"]=0
table(dd$case)
table(dd$mal)
#removing missing malpos data
dd2=dd[-which(dd$malpos=="NA"),]
dd$SITE=NA
dd$SITE[dd$site=="DD City"]=1
dd$SITE[dd$site=="DDU"]=0
table(dd$agecat)
dd$agecat<-recode(dd$agecat,'Under 5'="0","5 - 15 Years"="1","Above 15 Years"="2")
dd$agecat <- factor(dd$agecat,
                    levels = c(0,1, 2),
                    labels = c("Under 5", "5 - 15 Years", "Above 15 Years"))
### malaria prevalence in two site frequency table 
xtabs(~site+mal+casexp,dd)

# tabulation
xtabs( ~ site+casexp, dd)
###############################################################################
# summary statistics individual level for categorical variable  
#We use a for loop such that for all variables in “vars”, do as follows.
vars <- c("sex", "mal", "travel","irs","repellent","spray")

for(var in vars){
  
  #create a table for one variable at a time among cases and control
  freq_table <- table(dd[,var],dd$site,dd$casexp)
  #print variable name to as a header in the output
  print(var)
  # table of proportions
  print(prop.table(freq_table))
  #rename column names as per your output
  colnames(freq_table) <- c("n", "%")
  #print table with (n & %)
  print(freq_table)
  # repeats the loop till all variables in “vars” are complete
}
### summary statistics for numeric variable ########
summary(dd$age)
###############################################################################
########## summary statistics household level ################## 
# Household level distribution case and control
xtabs(~site+casexp,dd_hh)
## summary statistics household level 
colnames(dd_hh)
#We use a for loop such that for all variables in “vars”, do as follows.
vars <- c("stephlarv","stephadult","stephpos","livestockpres","eave",
          "waterbodytype","waterbody","irs")
for(var in vars){
  
  #create a table for one variable at a time among cases and control
  freq_table_hh <- table(dd_hh[,var],dd_hh$site,dd_hh$casexp)
  #print variable name to as a header in the output
  print(var)
  # table of proportions
  print(prop.table(freq_table_hh))
  #rename column names as per your output
  colnames(freq_table_hh) <- c("n", "%")
  #print table with (n & %)
  print(freq_table_hh)
  # repeats the loop till all variables in “vars” are complete
}
## To change the variables in to numeric 
dd_hh$dist_river<-as.numeric(dd_hh$dist_river)
dd_hh$dist_artificial_cont<-as.numeric(dd_hh$dist_artificial_cont)
### to summarize the variables using mean
dd_hh %>%
  group_by(site,casexp) %>%
  summarise_at(c("dist_river", "dist_artificial_cont"), mean, na.rm = TRUE)

###############################################################################

#### Case control analysis to identify risk factors Table S5
###############################################################################
library(plyr)
library(survival)
library(epiR)
## tabulation of case and site
table(dd$casexp, dd$site)
## water body presence in the neighborhood
waterbody <- table(dd$waterbody, dd$casexp)
#Apply epi.2by2 function to table.
cc_waterbody <- epi.2by2(waterbody, method = "case.control")
## spray use 
spray <- table(dd$spray, dd$casexp)
#Apply epi.2by2 function to table.
cc_spray<- epi.2by2(spray, method = "case.control")
## larvea presence
larvea <- table(dd$larvaepos, dd$casexp)
#Apply epi.2by2 function to table.
cc_larvea<- epi.2by2(larvea, method = "case.control")
## An. stephensi adult presence
adult <- table(dd$adultpos_2, dd$casexp)
#Apply epi.2by2 function to table.
cc_adult<- epi.2by2(adult, method = "case.control")
## An. stephensi larvae/adult presence
steppos <- table(dd$stephpos,dd$casexp)
#Apply epi.2by2 function to table.
cc_stephpos <- epi.2by2(steppos, method = "case.control")
##LLIN use
IRS <- table(dd$irs,dd$casexp)
#Apply epi.2by2 function to table.
cc_IRS <- epi.2by2(IRS, method = "case.control")
##Travel History
travel <- table(dd$travel,dd$casexp)
#Apply epi.2by2 function to table.
cc_travel <- epi.2by2(travel, method = "case.control")
##Eave opened
eave <- table(dd$eave,dd$casexp)
#Apply epi.2by2 function to table.
cc_eave <- epi.2by2(eave, method = "case.control")
#Livestock presence 
LS <- table(dd$livestockpres,dd$casexp)
#Apply epi.2by2 function to table.
cc_LS <- epi.2by2(LS, method = "case.control")
#Distance from Artifical container  
ds <- table(dd$distance_Artif_con_cat,dd$casexp)
#Apply epi.2by2 function to table.
cc_ds <- epi.2by2(ds, method = "case.control")
## sex 
sex <- table(dd$sex,dd$casexp)
#Apply epi.2by2 function to table.
cc_sex <- epi.2by2(sex, method = "case.control")
###############################################################################
