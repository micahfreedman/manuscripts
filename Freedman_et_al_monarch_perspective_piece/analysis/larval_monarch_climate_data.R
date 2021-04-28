##### Extracting monarch occurrence data for larvae, eggs, and pupae using to iNatTools package

#

library(remotes)
library(iNatTools)

#monarchs_pupae <- iNat(taxon_id = 48662, term_id = 1, term_value_id = 4)
#started at 1:18 pm
#Sys.time()
#finished at 1:20
#nrow(monarchs_pupae) #returns 3484 records
#monarchs_pupae$observed_on

#cannot handle all larval records simultaneously because there are >10,000, so split across months
#monarchs_larvae_june <- iNat(taxon_id = 48662, term_id = 1, term_value_id = 6, month = 6)
#2711 larval records for June

#monarchs_larvae_july <- iNat(taxon_id = 48662, term_id = 1, term_value_id = 6, month = 7)
#5083 larval records for July

#monarchs_larvae_august <- iNat(taxon_id = 48662, term_id = 1, term_value_id = 6, month = 8)
#7999 records for August

#now for eggs
#monarchs_eggs <- iNat(taxon_id = 48662, term_id = 1, term_value_id = 7)
#1219 records of eggs

##  Now group altogether
#immature_monarchs <- rbind(monarchs_pupae, monarchs_larvae_june, monarchs_larvae_july, monarchs_larvae_august, monarchs_eggs)


#only keep columns for quality_grade, observed_on, captive, location

#monarch_records <- immature_monarchs[,colnames(immature_monarchs) %in% c('quality_grade', 'observed_on', 'captive', 'location')]

#write.csv(monarch_records, file = '~/Desktop/immature_monarchs_inat.csv')

monarch_records <- read.csv(file = '~/Desktop/immature_monarchs_inat.csv')

monarch_records$latitude <- as.numeric(gsub(",.*$", "", monarch_records$location))
monarch_records$longitude <- as.numeric(sub('.*,', '', monarch_records$location))

monarch_records <- monarch_records[monarch_records$quality_grade=='research' & monarch_records$captive==FALSE,]

monarch_records$year <- as.numeric(substr(monarch_records$observed_on, 1, 4))
hist(monarch_records$year, breaks = 50) #can see exponential increase in observations

monarch_records$month <- as.numeric(trimws(sub("^[^-]*-([^-]*)-.*$", "\\1", monarch_records$observed_on)))
hist(monarch_records$month) #as expected, vast majority of records are in June, July, August

monarch_records$day <- as.numeric(sub('.*\\-', '', monarch_records$observed_on))
hist(monarch_records$day) #as expected, no real pattern based on day of the month

monarch_records_summer <- monarch_records[monarch_records$month %in% c(6,7,8),]
#down to 17629 records for summer

#####################3

library(rworldmap)
library(raster)
library(sp) #manipulationg spatial data
library(ggplot2)
library(devtools) #needed to download prism from github
library(reshape2) ##melting dataframes
library(dplyr) #data wrangling
library(raster) ##working with raster data
library(rgdal)
library(prism) ##prism data access
library(RANN) #includes functions for nearest neighbor mathcing
library(geosphere)
library(maps)

north_am <- getMap(resolution = 'low')
plot(north_am, xlim=range(-130,-35), ylim = range(-15,50))
points(monarch_records_summer$longitude,monarch_records_summer$latitude,col='blue',cex=1,pch=20) #includes all records regardless of location

setwd('~/Documents/grad_school/Manuscripts/Perspective piece/analysis/')

r <- raster::getData("worldclim",var="bio",res=10) #takes a while to access data (file is 123.3 MB), and requires an internet connection

r2 <- r[[c(1:19)]]

names(r2) <- c("BIO1 = Annual Mean Temperature", "BIO2 = Mean Diurnal Range", "BIO3 = Isothermality",
               "BIO4 = Temperature Seasonality", "BIO5 = Max Temperature of Warmest Month", "BIO6 = Min Temperature of Coldest Month", "BIO7 = Temperature Annual Range (BIO5-BIO6)", "BIO8 = Mean Temperature of Wettest Quarter", "BIO9 = Mean Temperature of Driest Quarter", "BIO10 = Mean Temperature of Warmest Quarter", "BIO11 = Mean Temperature of Coldest Quarter", "BIO12 = Annual Precipitation", "BIO13 = Precipitation of Wettest Month", "BIO14 = Precipitation of Driest Month", "BIO15 = Precipitation Seasonality (Coefficient of Variation)", "BIO16 = Precipitation of Wettest Quarter", "BIO17 = Precipitation of Driest Quarter", "BIO18 = Precipitation of Warmest Quarter", "BIO19 = Precipitation of Coldest Quarter") #rename columns to give their full names/descriptions

monarch_records_summer <- na.omit(monarch_records_summer)

coords <- data.frame(x=monarch_records_summer$longitude,y=monarch_records_summer$latitude)
points <- SpatialPoints(coords, proj4string = r@crs)
values <- as.data.frame(extract(r2,points)) #takes a little while; this matches coordinates to their corresponding Bioclim data

monarch_records_summer <- cbind(monarch_records_summer, values)


monarch_records_summer$east.west <- ifelse(monarch_records_summer$long < -105, 'west', 'east')
monarch_records_summer[which(monarch_records_summer$long > -113.5 & monarch_records_summer$long < -105 & monarch_records_summer$lat >45 & monarch_records_summer$lat < 60),]$east.west <- 'east'
monarch_records_summer[which(monarch_records_summer$long > -114 & monarch_records_summer$long < -112 & monarch_records_summer$lat<47 & monarch_records_summer$lat > 45),]$east.west <- 'west'
monarch_records_summer[which(monarch_records_summer$long > -106 & monarch_records_summer$long < -104 & monarch_records_summer$lat<41 & monarch_records_summer$lat >39 ),]$east.west <- 'east'

monarch_records_summer$resident <- ifelse(monarch_records_summer$lat < 31, 'resident', 'non-resident')
monarch_records_summer[monarch_records_summer$lat < 34 & monarch_records_summer$lat > 32 & monarch_records_summer$long < -116 & monarch_records_summer$long > -119 ,]$resident <- 'resident'
monarch_records_summer[monarch_records_summer$lat > 32 & monarch_records_summer$lat < 35 & 
              monarch_records_summer$long < -75 & monarch_records_summer$long > -80.8,]$resident <- 'resident'
monarch_records_summer[monarch_records_summer$lat < 31 & monarch_records_summer$lat > 28 & monarch_records_summer$long < -97 & monarch_records_summer$long > -105,]$resident <- 'non-resident'

ggplot(monarch_records_summer[monarch_records_summer$resident!='resident',], aes(x = BIO18...Precipitation.of.Warmest.Quarter, col = east.west))+
  geom_density(aes(fill = east.west), adjust = 4, alpha = 0.2, size = 1)+
  theme_bw()+
  scale_color_manual(values = c('yellow','orange'))+
  scale_fill_manual(values = c('yellow','orange'))+
  geom_vline(xintercept = median(monarch_records_summer[monarch_records_summer$east.west=='west',]$BIO18...Precipitation.of.Warmest.Quarter, na.rm=T), col = 'black', lty = 2)+
  geom_vline(xintercept = median(monarch_records_summer[monarch_records_summer$east.west=='east',]$BIO18...Precipitation.of.Warmest.Quarter, na.rm=T), col = 'black', lty = 2)+
  geom_hline(yintercept=0, colour="white", size=1)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  xlab('BIO18 - Precipitation of Warmest Quarter')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16), legend.box.background = element_rect(colour = "black"))

(figs3 = ggplot(monarch_records_summer[monarch_records_summer$resident!='resident',], aes(x = BIO10...Mean.Temperature.of.Warmest.Quarter, col = east.west))+
  geom_density(aes(fill = east.west), adjust = 4, alpha = 0.2, size = 1)+
  theme_bw()+
  scale_color_manual(values = c('yellow','orange'))+
  scale_fill_manual(values = c('yellow','orange'))+
  geom_hline(yintercept=0, colour="white", size=1)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  xlab('BIO10 = Mean Temperature of Warmest Quarter')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16), legend.box.background = element_rect(colour = "black")))

ggsave(figs3, file = '~/Desktop/figs3.pdf', height = 6, width = 6)

##################3

prism_monarchs <- monarch_records_summer[monarch_records_summer$resident=='non-resident',]

aggregate(BIO10...Mean.Temperature.of.Warmest.Quarter ~ east.west, FUN = median, data = prism_monarchs) #no real difference when looking at mean temperature of warmest quarter, although this also takes into account night-time low temperatures

nrow(prism_monarchs)

plot(north_am, xlim=range(-130,-35), ylim = range(-15,50))
points(prism_monarchs$longitude,prism_monarchs$latitude,col=c('blue','orange')[factor(prism_monarchs$east.west)],cex=1,pch=20)

prism_monarchs$country <- map.where(database = 'world', prism_monarchs$longitude, prism_monarchs$latitude) #weirdly includes additional information beyound just country codes

prism_monarchs$country <- gsub(":.*$", "", prism_monarchs$country)

prism_monarchs <- na.omit(prism_monarchs)

prism_monarchs <- prism_monarchs[prism_monarchs$month %in% c(7,8),]

canadian.records <- prism_monarchs[prism_monarchs$country=='Canada',]

setwd('./PRISM_tmax_30yr_normal_4kmM2_all_bil/')

get_prism_normals(type = 'tmax', resolution = '4km', mon = 7:8, keepZip = TRUE) #get 4km resolution data for July and August
prism_archive_ls() #first two entries are the ones containing the 30 year normals for July and August, respectively

pd_image(ls_prism_data()[1,1]) ###this shows 30 year normals for July
pd_image(ls_prism_data()[2,1]) #same for August

prism.monarchs.coords <- prism_monarchs[,c(7,6)]
names(prism.monarchs.coords) <- c('long','lat')

mystack <- ls_prism_data()[c(1,2),1] %>%  prism_stack(.)

mycrs <- mystack@crs@projargs

coordinates(prism.monarchs.coords) <- c('long','lat')

proj4string(prism.monarchs.coords) <- CRS(mycrs)

data <- data.frame(coordinates(prism.monarchs.coords), extract(mystack, prism.monarchs.coords))

names(data)[3:4] <- c('July_Normals_Tmax', 'August_Normals_Tmax')

###### now load all Canadian data

ontario <- read.csv('~/Downloads/climate-normals-ontario.csv')
quebec <- read.csv('~/Downloads/climate-normals-quebec.csv')
alberta <- read.csv('~/Downloads/climate-normals-alberta.csv')
saskatchewan <- read.csv('~/Downloads/climate-normals-saskatchewan.csv')
bc <- read.csv('~/Downloads/climate-normals-BC.csv')
manitoba <- read.csv('~/Downloads/climate-normals-manitoba.csv')
new_brunswick <- read.csv('~/Downloads/climate-normals-new_brunswick.csv')
newfoundland <- read.csv('~/Downloads/climate-normals-newfoundland.csv')
nova_scotia <- read.csv('~/Downloads/climate-normals-nova_scotia.csv')
prince_edward <- read.csv('~/Downloads/climate-normals-prince_edward.csv')

canada <- rbind(ontario,quebec,alberta,saskatchewan,bc,manitoba,new_brunswick,newfoundland,nova_scotia,prince_edward)

nrow(canada) #over 560,000 lines of data

#restrict to only the necessary elements: tmax in July and August, latitude, longitude (maybe also keep station ID info)

names(canada)

columns.keep <- c('x','y','STATION_NAME','PROVINCE_CODE','E_NORMAL_ELEMENT_NAME','MONTH','VALUE')

canada.keep <- canada[,columns.keep]

levels(canada.keep$E_NORMAL_ELEMENT_NAME) #want to save entries for 'Mean daily max temperature deg C' and 'Mean daily min temperature deg C', as well as months 7 and 8 (july and august)

canada.keep <- canada.keep[canada.keep$E_NORMAL_ELEMENT_NAME %in% c("Mean daily max temperature deg C","Mean daily min temperature deg C") & canada.keep$MONTH %in% c(7,8),]

nrow(canada.keep) #restricts data to 2512 entries; corresponds to data for 628 stations (out of 686 total)

head(canada.keep)

names(canada.keep)[c(1,2)] <- c('longitude','latitude')

prism_monarchs_canada <- prism_monarchs[prism_monarchs$country=='Canada',]

HavMat <- distm(prism_monarchs[,c('longitude','latitude')], canada.keep[,c('longitude','latitude')], fun=distHaversine)

head(HavMat)

prism_monarchs_canada$new.var <- prism_monarchs_canada$STATION_NAME[max.col(HavMat)]


d <- pointDistance(prism_monarchs_canada[,7:6], canada.keep[,1:2], lonlat=TRUE, allpairs=T) 
i <- apply(d, 1, which.min)

prism_monarchs_canada$STATION_NAME <- canada.keep$STATION_NAME[i]

use.canada <- merge(prism_monarchs_canada, canada.keep[,c(3,4,5,7)], by = 'STATION_NAME')

library(plyr)

use.canada$E_NORMAL_ELEMENT_NAME <- revalue(use.canada$E_NORMAL_ELEMENT_NAME, c("Mean daily max temperature deg C"="Tmax", "Mean daily min temperature deg C"="Tmin"))

head(use.canada)

use.canada$measure <- as.factor(paste(use.canada$month, use.canada$E_NORMAL_ELEMENT_NAME, sep = "_"))

use.canada$id <- paste(use.canada$observed_on, use.canada$latitude)

dst <- dcast(use.canada, id ~ E_NORMAL_ELEMENT_NAME, value.var = 'VALUE', fun.aggregate = mean) #takes the mean July and August temperature for each location, reports as Tmax, but still gives four records for each locatoin; aggregate again

dst <- merge(dst, use.canada, by = 'id')

dst <- aggregate(Tmax ~ id + latitude + longitude, data = dst, FUN = mean) 
dst$long <- dst$longitude
dst$lat <- dst$latitude

data$Tmax <- rowMeans(data[,c(3,4)], na.rm = T)
data = na.omit(data)
names(data)[1:2] <- c('longitude','latitude')
data = data[,c(2,1,5)]
dst = dst[,c(2,3,4)]

fin.dat <- rbind(data, dst)

fin.dat$east.west <- ifelse(fin.dat$long < -105, 'west', 'east')
fin.dat[which(fin.dat$long > -113.5 & fin.dat$long < -105 & fin.dat$lat >45 & fin.dat$lat < 60),]$east.west <- 'east'
fin.dat[which(fin.dat$long > -114 & fin.dat$long < -112 & fin.dat$lat<47 & fin.dat$lat > 45),]$east.west <- 'west'
fin.dat[which(fin.dat$long > -106 & fin.dat$long < -104 & fin.dat$lat<41 & fin.dat$lat >39 ),]$east.west <- 'east'


fin.dat$resident <- ifelse(fin.dat$lat < 31, 'resident', 'non-resident')
fin.dat[fin.dat$lat < 34 & fin.dat$lat > 32 & fin.dat$long < -116.5 & fin.dat$long > -119 ,]$resident <- 'resident'
fin.dat[fin.dat$lat > 32 & fin.dat$lat < 35 & 
          fin.dat$long < -75 & fin.dat$long > -80.8,]$resident <- 'resident'
fin.dat[fin.dat$lat < 31 & fin.dat$lat > 28 & fin.dat$long < -97 & fin.dat$long > -105,]$resident <- 'non-resident'

fin.dat[fin.dat$resident=='resident',]$east.west <- 'resident'

table(fin.dat$east.west) #looks right

fin.dat.plot <- droplevels(fin.dat[!fin.dat$east.west=='resident',])

fin.dat.plot$east.west <- as.factor(fin.dat.plot$east.west)

fin.dat.plot$east.west <- revalue(fin.dat.plot$east.west, c('east' = 'eastern', 'west' = 'western'))

ggplot(fin.dat.plot, aes(x = Tmax, col = east.west))+
  geom_density(aes(fill = east.west), adjust = 4, alpha = 0.2, size = 1)+
  theme_bw()+
  scale_color_manual(values = c('yellow','orange'))+
  scale_fill_manual(values = c('yellow','orange'))+
  geom_vline(xintercept = median(fin.dat.plot[fin.dat.plot$east.west=='western',]$Tmax, na.rm=T), col = 'black', lty = 2)+
  geom_vline(xintercept = median(fin.dat.plot[fin.dat.plot$east.west=='eastern',]$Tmax, na.rm=T), col = 'black', lty = 2)+
  geom_hline(yintercept=0, colour="white", size=1)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  xlab('Maximum Daily Temperature (July - August) (°C)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16), legend.box.background = element_rect(colour = "black"))

aggregate(Tmax ~ east.west, FUN = median, fin.dat.plot) #as with adults, immature stages also experience warmer summertime high temperatures in the West

plot(north_am, xlim=range(-130,-60), ylim = range(20,40))
points(fin.dat.plot$longitude,fin.dat.plot$latitude,col=fin.dat.plot$east.west,cex=1,pch=20)

library(sf)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

(figs1 <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = fin.dat.plot, aes(x = longitude, y = latitude, col = east.west)) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE)+
  scale_color_manual(values = c('yellow','orange'))+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  theme(panel.background = element_rect(fill = 'aliceblue'))+
  theme(legend.position = 'none')+
  theme(axis.title = element_blank(), axis.text = element_blank()))

ggsave(figs1, file = '~/Desktop/figs1.pdf', height = 8, width = 8)
