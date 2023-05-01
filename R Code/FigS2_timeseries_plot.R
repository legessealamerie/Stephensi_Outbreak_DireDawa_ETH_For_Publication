
### Malaria  Trend Analysis using  Dire Dawa  secondary data collected from 
## 34 Health centers.
## required library
library(trend)
### importing data
ds=read.csv("Secondary_data_dd.csv",header = T)  
colnames(ds)
## create time series data 
dd_ts <- ts(ds$positive_case,start=c(2019,1),frequency=12)
my.limit=c(2019,2020,2021,2022)
my.limit=as.numeric(my.limit)
### time series plot from 2019-2022
plot(dd_ts,ylab="Number of case", xlab="year",xlim=c(2019,2023))
## to display the time series data 
dd_ts
### trend Analysis to check monotonic increasing or decreasing trend in malaria
## Seasonal Mann-Kendall trend test
trend_test <- smk.test(dd_ts)
trend_test
summary(trend_test)
###################################