library(rgdal)
library(raster)
library(ggplot2)
library(spatstat)
library(plotrix)
library(fields)
library(leaflet)
library(sf)
library(dplyr)
library(tidyr)
library(maptools)
library(RColorBrewer)
library(lattice)
library(geoR)
library(plotrix)
library(car)
library(sp)

# Moran's I and spatial dependencies
library(spdep) # Spatial Dependence: Weighting Schemes, Statistics and Models
library(ape) # Analyses of Phylogenetics and Evolution
library(pgirmess) # Data Analysis in Ecology

library(sf)

#### read in data and plot to check coordinates correct ####
dd5 <- read.csv("Data/dd3.individual.final.classifications.csv")

dd5$longitude <- as.integer(dd5$longitude)
dd5$latitude <- as.integer(dd5$latitude)

dd5.full <- dd5 %>% drop_na(longitude)

#transform coordinates to  utm 37N
wkt1 <- sf::st_crs(32637)[[2]]  
sps.37N <- SpatialPoints(dd5.full[, c('longitude', 'latitude')], sp::CRS(wkt1)) #transform warns, this can be ignored
sps.84 <- spTransform(sps.37N, CRS("+proj=longlat +datum=WGS84"))
dd5.full[, c("long", "lat")] <- coordinates(sps.84)

#subsite to include city points only
dd5.full.goro <- subset(dd5.full, site == "DD City")

#### GLOBAL MORANS I ####
#testing for spatial autocorrelation in household malaria prevalence at point-level 
# logistic transformation (log odds) to make more normal
dd5.full.goro$log_odds <- logit(dd5.full.goro$malpos.hh) 

#aggregate to hh level
dd5.clus <- dd5.full.goro %>%
  group_by(idhh) %>%
  filter(row_number()==1)

dd5.clus <- subset(dd5.clus, select = c("malpos.hh", "long", "lat", "log_odds"))

xy=cbind(dd5.clus$long, dd5.clus$lat)

coords<-coordinates(xy) # set spatial coordinates to create a spatial object
IDs<-row.names(as.data.frame(coords))

# choose a distance d such that pairs of points with distances less than d are neighbors and those further apart are not. 
Neigh_nb<-knn2nb(knearneigh(coords, k=1, longlat = TRUE), row.names=IDs)     # using the "spdep" package
# assigns at least one neighbor to each and calculates the distances between
dsts<-unlist(nbdists(Neigh_nb,coords)) # returns the distance between nearest neighbors for each point
summary(dsts)

max_1nn<-max(dsts)
max_1nn # maximum distance to provide at least one neighbor to each point

quantile(dsts, prob=c(.25,.5,.75), type=1)


# Create different neighbor structures based upon distance
Neigh_kd1<-dnearneigh(coords,d1=0, d2=0.02, row.names=IDs)   # neighbors within maximum distance


#assign weights; 
weights<-nb2listw(Neigh_kd1, style="W", zero.policy = TRUE)   # row standardized binary weights, using minimum distance for one neighbor
print.listw(weights, zero.policy = TRUE)     
options(scipen=999)

moran.test(dd5.clus$log_odds , listw=weights, na.action = na.omit, zero.policy = TRUE)  #using row standardised weights

#simulate random permutations of prevalence  - then compare to observed statistic
bperm<-moran.mc(dd5.clus$log_odds , listw=weights,nsim=999, na.action = na.omit, zero.policy = TRUE)
bperm

#### LOCAL MORANS I ####
# First calculate the local Moran's I around each point based on the spatial weights object (binary based on at least one neighbor)
I <-localmoran(dd5.clus$log_odds, weights, na.action = na.omit, zero.policy = TRUE) # "spdep" package

# Print 'LISA' for each point
Coef<-printCoefmat(data.frame(I[IDs,], row.names=row.names(coords),
                              check.names=FALSE))

# Plot the spatial data against its spatially lagged values (the weighted mean of its neighbors)                         
par(mfrow=c(1,1))

nci<-moran.plot(dd5.clus$log_odds, listw=weights, 
                xlab="Log prevalence", ylab="Spatially lagged log prev", labels=T, pch=16, col="grey", zero.policy = TRUE)
text(c(3,3, -5,-5),c(0.9, -1.9,0.9,-1.9), c("High-High", "High-Low", "Low-High", "Low-Low"), cex=0.8)


#code clusters as high high etc 
x<-dd5.clus$log_odds
lhx<-cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=T)

# Map points that are local outliers in the plot
infl<-nci$is_inf==T # find which points are statistically significant outliers
sum(infl==T)    #14 true (12% - more than would expect by chance)

#take significant clusters 
wx<-stats::lag(weights,dd5.clus$log_odds, zero.policy = TRUE)
lhwx<-cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=T)
lhlh<-interaction(lhx,lhwx,infl,drop=T)

names<-rep("none", length(lhlh))
names[lhlh=="L.L.TRUE"]<-"LL"
names[lhlh=="H.L.TRUE"]<-"HL"
names[lhlh=="L.H.TRUE"]<-"LH"
names[lhlh=="H.H.TRUE"]<-"HH"

# map to show local clusters
dd_localM<-as.data.frame(cbind(xy,names))
colnames(dd_localM)<-c("long", "lat", "names")
dd_localM[c("long", "lat")] <- lapply( dd_localM[c("long", "lat")], function(x) as.numeric(as.character(x)) )
factpal <- colorFactor(c( "red","coral4","coral","green","lightgrey"), names)
leaflet(dd_localM) %>% addTiles() %>% addCircleMarkers(~long, ~lat, fillOpacity=1,
                                                       color= ~factpal(names), radius=4, stroke=TRUE, weight=1) %>% 
  addLegend(pal = factpal, values = ~names, title="Class")

table(dd_localM$names)

### prevalence and clustering plot #### 
clus <- subset(dd_localM, names %in% c("HH", "LL", "HL"))
clus <- cbind(clus, malpos.hh=NA)


prev <- subset(dd5.clus, select = c("long", "lat", "malpos.hh"))
prev <- cbind(prev, names=NA)
prev

map <- rbind(clus, prev)


pal = colorNumeric("Reds", prev$malpos.hh)
leaflet(prev) %>%addTiles() %>% addCircleMarkers(~long, ~lat, fillOpacity=1,
                                                 fillColor= ~pal(malpos.hh), radius = 4, stroke=TRUE, 
                                                 color = "black", weight=1) %>% 
  addLegend("bottomright", pal = pal, values = ~malpos.hh, 
            na.label = "", title = "Household prevalence") %>% 
  addMarkers()

write.csv(map, "Data/clusters results.csv")
