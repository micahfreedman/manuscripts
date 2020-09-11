#### load gbif data 

start_time <- Sys.time()
gbif.monarchs <- read.csv('~/Downloads/monarchs_gbif.csv') #be patient: file is 92 Mb and contains over 212,000 records; takes about 21 seconds to load
end_time <- Sys.time()
start_time - end_time

names(gbif.monarchs)

cols.save <- c('gbifID','countryCode','decimalLatitude','decimalLongitude','coordinatePrecision','elevation','eventDate','day','month','year')

gbif <- gbif.monarchs[,cols.save]

#choose criteria for which records to keep

#specify only records from June, July, August (meteorological summer)

gbif <- gbif[gbif$month %in% c(6,7,8),] #resticts dataset to 28,869 observations

#use United States and Canada specimens, since this is the core of the summer breeding range

gbif <- gbif[gbif$countryCode %in% c('US','CA','MX'),] #contains 27132 observations

#restrict to specimens that have decimal coordinates available

gbif <- gbif[!is.na(gbif$decimalLatitude),] #almost all were georeferenced, so now have 26,759 observations

#plot observations

library(rworldmap)

north_am <- getMap(resolution = 'low')

plot(north_am, xlim=range(-130,-35), ylim = range(-15,50))

points(gbif$decimalLongitude,gbif$decimalLatitude,col='blue',cex=1,pch=20)

#make the same map, but separate points according to whether they came from pre- or post-1960

gbif$decade <- as.factor(with(gbif, ifelse(year < 1980, 'Pre-1980','Post-1980')))

plot(north_am, xlim=range(-125,-65), ylim = range(30,40))

points(gbif$decimalLongitude,gbif$decimalLatitude,col=as.factor(levels(gbif$decade)),cex=0.5,pch=20)
legend('bottomleft', legend = c('Pre-1980','Post-1980'), col = c('black','red'), pch = 20)

###now pull climate data from bioclim

library(raster)
library(sp)

r <- getData("worldclim",var="bio",res=2.5)

r2 <- r[[c(1:19)]]

names(r2) <- c("BIO1 = Annual Mean Temperature", "BIO2 = Mean Diurnal Range", "BIO3 = Isothermality",
               "BIO4 = Temperature Seasonality", "BIO5 = Max Temperature of Warmest Month", "BIO6 = Min Temperature of Coldest Month", "BIO7 = Temperature Annual Range (BIO5-BIO6)", "BIO8 = Mean Temperature of Wettest Quarter", "BIO9 = Mean Temperature of Driest Quarter", "BIO10 = Mean Temperature of Warmest Quarter", "BIO11 = Mean Temperature of Coldest Quarter", "BIO12 = Annual Precipitation", "BIO13 = Precipitation of Wettest Month", "BIO14 = Precipitation of Driest Month", "BIO15 = Precipitation Seasonality (Coefficient of Variation)", "BIO16 = Precipitation of Wettest Quarter", "BIO17 = Precipitation of Driest Quarter", "BIO18 = Precipitation of Warmest Quarter", "BIO19 = Precipitation of Coldest Quarter")

lats <- gbif$decimalLatitude
lons <- gbif$decimalLongitude

coords <- data.frame(x=lons,y=lats)


start_time <- Sys.time()
points <- SpatialPoints(coords, proj4string = r@crs)
end_time

values <- extract(r2,points)

(values<-as.data.frame(values))

(bioclim.vars <- cbind(gbif, values))

#tentatively split map in west and east

bioclim.vars$eastern.vs.western <- ifelse(bioclim.vars$decimalLongitude < -107.5, 'western', 'eastern')

plot(north_am, xlim=range(-125,-65), ylim = range(30,50))

plot.bioclim.vars <- bioclim.vars[bioclim.vars$decimalLatitude > 25,]

ew.cols <- c('blue','orange')

points(plot.bioclim.vars$decimalLongitude,plot.bioclim.vars$decimalLatitude,col=ew.cols[as.factor((plot.bioclim.vars$eastern.vs.western))],cex=0.5,pch=20)

##### also may want to exclude putatively non-migratory butterflies from south Florida and southern California

##########

library(ggplot2)

ggplot(bioclim.vars, aes(x = BIO1...Annual.Mean.Temperature, col = eastern.vs.western))+
  stat_density(aes(fill = eastern.vs.western), adjust = 2.5, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='western',]$BIO1...Annual.Mean.Temperature, na.rm = T), col = 'orange', lty = 2)+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='eastern',]$BIO1...Annual.Mean.Temperature, na.rm = T), col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))

ggplot(bioclim.vars, aes(x = BIO10...Mean.Temperature.of.Warmest.Quarter, col = eastern.vs.western))+
  stat_density(aes(fill = eastern.vs.western), adjust = 2.5, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='western',]$BIO10...Mean.Temperature.of.Warmest.Quarter, na.rm = T), col = 'orange', lty = 2)+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='eastern',]$BIO10...Mean.Temperature.of.Warmest.Quarter, na.rm = T), col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))

aggregate(BIO1...Annual.Mean.Temperature ~ eastern.vs.western, mean, data = bioclim.vars) #confirms that western monarchs are captured in areas with slightly higher mean annual temperature

aggregate(BIO2...Mean.Diurnal.Range ~ eastern.vs.western, mean, data = bioclim.vars) #greater mean diurnal range for western samples

aggregate(BIO3...Isothermality ~ eastern.vs.western, mean, data = bioclim.vars) #much greater isothermality for western samples

aggregate(BIO4...Temperature.Seasonality ~ eastern.vs.western, mean, data = bioclim.vars) #east north american has higher temperature seasonality

aggregate(BIO5...Max.Temperature.of.Warmest.Month ~ eastern.vs.western, mean, data = bioclim.vars) #nearly identical for eastern and western: maybe western isn't actually that much hotter?

ggplot(bioclim.vars, aes(x = BIO5...Max.Temperature.of.Warmest.Month, col = eastern.vs.western))+
  stat_density(aes(fill = eastern.vs.western), adjust = 2.5, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='western',]$BIO5...Max.Temperature.of.Warmest.Month, na.rm = T), col = 'orange', lty = 2)+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='eastern',]$BIO5...Max.Temperature.of.Warmest.Month, na.rm = T), col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))


aggregate(BIO10...Mean.Temperature.of.Warmest.Quarter ~ eastern.vs.western, mean, data = bioclim.vars) #eastern monarchs actually come from places that have warmer summers overall

aggregate(BIO15...Precipitation.Seasonality..Coefficient.of.Variation. ~ eastern.vs.western, mean, data = bioclim.vars) #not surprsingly, coefficient of variation in precip is much more variable for western samples

ggplot(bioclim.vars, aes(x = BIO15...Precipitation.Seasonality..Coefficient.of.Variation., col = eastern.vs.western))+
  stat_density(aes(fill = eastern.vs.western), adjust = 2.5, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='western',]$BIO15...Precipitation.Seasonality..Coefficient.of.Variation., na.rm = T), col = 'orange', lty = 2)+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='eastern',]$BIO15...Precipitation.Seasonality..Coefficient.of.Variation., na.rm = T), col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))

ggplot(bioclim.vars, aes(x = BIO2...Mean.Diurnal.Range, col = eastern.vs.western))+
  stat_density(aes(fill = eastern.vs.western), adjust = 2.5, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='western',]$BIO2...Mean.Diurnal.Range, na.rm = T), col = 'orange', lty = 2)+
  geom_vline(xintercept = mean(bioclim.vars[bioclim.vars$eastern.vs.western=='eastern',]$BIO2...Mean.Diurnal.Range, na.rm = T), col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))




#### alternative approach: see if it's possible to get daily temperature data for particular locations

#library("rnoaa")
#library("lawn")

#m.noaa <- gbif

#my.id <- 'hXlboUhvfeLHNskwddRzxfoPbvOBGoVZ'

#ncdc(datasetid = 'GHCND', stationid = 'GHCND:USW00014895', startdate = '2013-10-01',
    # enddate = '2013-10-01', token = my.id) #returns output, so the general function does seem to work

#now need to match coordinates to their nearest NOAA station

#station_data <- ghcnd_stations() #unclear why, but this only works sometimes and not others....usually returns a curl_fetch_memory error, with no interpretable reason for why; interface is not stable, avoid for any serious uses

#drop all entries that don't have lat/long:

#m.noaa <- m.noaa[!is.na(m.noaa$decimalLatitude),]#drops down to 27745 entries

#m.noaa <- m.noaa[!m.noaa$decimalLatitude==0,] #drop all 0 entries, now down to 27711

#test.locs <- m.noaa[,c(1,3:4)] 

#colnames(test.locs)[c(1,2,3)] <- c('id','latitude','longitude') #CRITICAL: must change the name of the first column to 'ID' in order for function to work

#nn <- meteo_nearby_stations(lat_lon_df = test.locs[1:500,], station_data = station_data,
                      #limit = 1, var = c("TMAX","TMIN"),
                      #year_min = 1871, year_max = 2019) #takes about 15 minutes to run, but does work; however, seems glitchy and doesn't actually match (stations it returns are very obviously wrong); seems like too unstable of a function to use

#lat_lon_df <- data.frame(id = c("sydney", "brisbane"),
                        # latitude = c(-33.8675, -27.4710),
                         #longitude = c(151.2070, 153.0234))
#nearby_stations <-  meteo_nearby_stations(lat_lon_df = lat_lon_df,
                                        #  station_data = station_data, radius = 10) #this does seem to work as well (example taken straight from documentation)

#miami <- data.frame(id = "miami", latitude = 25.7617, longitude = -80.1918)

# Get all stations within 50 kilometers
#meteo_nearby_stations(lat_lon_df = miami, station_data = ghcnd_stations(),
                      #radius = 50, var = c("PRCP", "TMAX"),
                      #year_min = 1992, year_max = 1992)

## try using nearest neighbor function to pinpoint nearest station and its ID

library(dplyr)
library(sf)
library(RANN)

#station_data_matching <- station_data$data[,c(1,4,5,7,9)] #issue is that NOAA website is horrendous and only allows 25 station IDs to be accessed at a time; completely useless

#closest <- nn2(station_data_matching[,c(2,5)], test.locs[,2:3], k = 1, searchtype = "radius", radius = 0.5) 

#station_data[closest$nn.idx[1],] #return information for closest station for row 1 of observations
#test.locs[1,] #confirms that this is in fact very close

#test.locs$station_ID <- station_data[closest[1]$nn.idx,]
#test.locs$distance_to_station <- closest[2]$nn.dists

#test.locs <- test.locs[test.locs$distance_to_station < 0.2,] #restrict to observations that are within 0.2 decimal degrees (~20 km) of the nearest recording station

#nrow(test.locs) #leaves 26,505 records

#hist(test.locs$distance_to_station) #most observations come from within 0.01 decimal degrees, which is <1 km away and should be very good



######## try working with prism data

library(devtools) #needed to download prism from github
library(reshape2) ##melting dataframes
library(dplyr) #data wrangling
library(raster) ##working with raster data
library(sp) ##manipulationg spatial data
library(rgdal)
library(prism) ##prism data access

options(prism.path = '~/Downloads/PRISM_tmax_30yr_normal_4kmM2_all_bil/')

setwd('~/Downloads/PRISM_tmax_30yr_normal_4kmM2_all_bil/')

get_prism_monthlys(type = 'tmax', years=1990:2018, mon = 1:12, keepZip = TRUE)
ls_prism_data(name=TRUE)

prism_image(ls_prism_data()[441,1]) #note: in order for this function to work, must choose the unzipped _bil files, e.g. PRISM_tmax_stable_4kmM2_199001_bil

get_prism_normals(type = 'tmax', resolution = '4km', mon = 7:8, keepZip = TRUE)
ls_prism_data()

(july.image <- prism_image(ls_prism_data()[43,1])) ###this shows 30 year normals for July
prism_image(ls_prism_data()[51,1]) #same for August

RS.july <- prism_stack(ls_prism_data()[43,1])
proj4string(RS.july)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

df.july <- data.frame(rasterToPoints(RS.july))
july.df <- melt(df.july, c("x", "y"))
names(july.df)[c(1:2, 4)] <- c("longitude", "latitude", "july_value")
nrow(july.df) #even with 4 km data, produces 481631 occcurrences

RS.august <- prism_stack(ls_prism_data()[51,1])
proj4string(RS.august)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

df.august <- data.frame(rasterToPoints(RS.august))
august.df <- melt(df.august, c("x", "y"))
names(august.df)[c(1:2, 4)] <- c("longitude", "latitude", "august_value")

#now cbind these two dataframes

summer.df <- cbind(july.df, august.df[,4])
names(summer.df)[5] <- 'august_value'

###restrict gbif records to only US samples from July, and from 1980-2018 (period of data collection)

prism.monarchs.summer <- gbif[gbif$countryCode=='US' & gbif$year %in% c(1980:2018) & gbif$month %in% 7:8,] 
table(prism.monarchs.summer$month) #leaves 18209 records, with 6307 from July and 11902 from August

prism.monarchs.summer <- prism.monarchs.summer[!prism.monarchs.summer$decimalLatitude == 0,] #drop one value that has 0,0 for coordinates

#do nearest neighbor matching for these two dataframes

head(prism.monarchs.summer)

closest.prism <- nn2(summer.df[,c(1,2)], prism.monarchs.summer[,c(4,3)], k = 1, searchtype = 'radius', radius = 0.5)
  
names(closest.prism)

### alternative: 

prism.monarchs.coords <- prism.monarchs.summer[,c(4,3)]
prism.monarchs.coords$id <- seq(1:nrow(prism.monarchs.coords))
names(prism.monarchs.coords) <- c('long','lat')


mystack <- ls_prism_data()[c(43,51),1] %>%  prism_stack(.)
mycrs <- mystack@crs@projargs

coordinates(prism.monarchs.coords) <- c('long','lat')

proj4string(prism.monarchs.coords) <- CRS(mycrs)

data <- data.frame(coordinates(prism.monarchs.coords), extract(mystack, prism.monarchs.coords))

names(data)[3:4] <- c('July_Normals_Tmax', 'August_Normals_Tmax')

data$east.west <- ifelse(data$long < -105 & data$lat < 42.5, 'western', ifelse(
  data$long < -106 & data$lat < 47, 'western', ifelse(data$long < -111 & data$lat < 60, 'western', 'eastern')))

data$ew2 <- ifelse(data$long < -102, 'western', 'eastern')

#plot data

plot(north_am, xlim=range(-125,-65), ylim = range(30,50))

points(data$lon,data$lat,col=ew.cols[as.factor((data$ew2))],cex=0.5,pch=20)
                           
data$summer.normals <- rowMeans(data[,3:4], na.rm = T)

ggplot(data, aes(x = summer.normals, col = ew2))+
  geom_density(aes(fill = ew2), adjust = 2.5, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = median(data[data$ew2=='western',]$summer.normals, na.rm=T), col = 'orange', lty = 2)+
  geom_vline(xintercept = median(data[data$ew2=='eastern',]$summer.normals, na.rm=T), 
             col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))

### eliminate Flordia from the eastern samples and LA and San Diego from western samples

data2 <- data[!data$lat < 30,]

plot(north_am, xlim=range(-125,-65), ylim = range(30,50))

points(data2$lon,data2$lat,col=ew.cols[as.factor((data2$ew2))],cex=0.5,pch=20)

ggplot(data2, aes(x = summer.normals, col = ew2))+
  geom_density(aes(fill = ew2), adjust = 4, alpha = 0.4)+
  theme_bw()+
  scale_color_manual(values = c('blue','orange'))+
  scale_fill_manual(values = c('blue','orange'))+
  geom_vline(xintercept = median(data2[data2$ew2=='western',]$summer.normals, na.rm=T), col = 'orange', lty = 2)+
  geom_vline(xintercept = median(data2[data2$ew2=='eastern',]$summer.normals, na.rm=T), 
             col = 'blue', lty = 2)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  xlab('Maximum Daily Temperature (July - August)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))

table(data2$east.west) #more than 10x as many eastern as western samples






###### load all canadian data

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

ggplot(canada.keep[canada.keep$E_NORMAL_ELEMENT_NAME=='Mean daily max temperature deg C' & canada.keep$MONTH==7,], aes(x = VALUE))+
  geom_density()+
  facet_wrap(~PROVINCE_CODE, scales = 'free')

### need to match canadian records to their nearest station

gbif.canada <- gbif[gbif$countryCode=='CA' & gbif$month %in% c(7,8),]

nrow(gbif.canada) #2879 canadian records from july and august

gbif.canada <- gbif.canada[!gbif.canada$decimalLatitude == 0,] #drop all entries where coordinates are not known

nrow(gbif.canada) #2855 records, with a mean collection latitude of 45.13

names(gbif.canada)[3:4] <- c('latitude','longitude')

library(geosphere)

#mat <- distm(gbif.canada[,c('longitude','latitude')], canada.keep[,c('longitude','latitude')], fun=distVincentyEllipsoid) #distVincentyEllipsoid is the most computationally intensive option here and takes a while to run, about 30 minutes; consider redoing with distHaversine function instead

HavMat <- distm(gbif.canada[,c('longitude','latitude')], canada.keep[,c('longitude','latitude')], fun=distHaversine) #this only takes about 10 seconds to run

head(HavMat)

gbif.canada$STATION_NAME <- canada.keep$STATION_NAME[max.col(-HavMat)]

use.canada <- merge(gbif.canada, canada.keep[,c(3,4,5,7)], by = 'STATION_NAME')

use.canada$E_NORMAL_ELEMENT_NAME <- droplevels(use.canada$E_NORMAL_ELEMENT_NAME)#drop unused climate variable levels

#rename factor levels

library(plyr)

use.canada$E_NORMAL_ELEMENT_NAME <- revalue(use.canada$E_NORMAL_ELEMENT_NAME, c("Mean daily max temperature deg C"="Tmax", "Mean daily min temperature deg C"="Tmin"))

head(use.canada)

use.canada$measure <- as.factor(paste(use.canada$month, use.canada$E_NORMAL_ELEMENT_NAME, sep = "_"))

dst <- dcast(use.canada, gbifID ~ E_NORMAL_ELEMENT_NAME, value.var = 'VALUE', fun.aggregate = mean) #takes the mean July and August temperature for each location, reports as Tmax, but still gives four records for each locatoin; aggregate again
dst <- merge(dst, use.canada[,c(2,4,5)])

dst <- aggregate(Tmax ~ gbifID + latitude + longitude, data = dst, FUN = mean) 

head(data) #revisit original dataframe with climate data

data$east.west <- NULL
data$ew2 <- NULL

dst$long <- dst$longitude
dst$lat <- dst$latitude
dst$July_Normals_Tmax <- NA
dst$August_Normals_Tmax <- NA
dst$summer.normals <- dst$Tmax
dst$gbifID <- NULL
dst$Tmax <- NULL
dst$longitude <- NULL
dst$latitude <- NULL

fin.dat <- rbind(data, dst)

head(fin.dat)

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

(fig1b <- ggplot(fin.dat.plot, aes(x = summer.normals, col = east.west))+
  geom_density(aes(fill = east.west), adjust = 4, alpha = 0.2, size = 1)+
  theme_bw()+
  scale_color_manual(values = c('yellow','orange'))+
  scale_fill_manual(values = c('yellow','orange'))+
  geom_vline(xintercept = median(fin.dat.plot[fin.dat.plot$east.west=='western',]$summer.normals, na.rm=T), col = 'black', lty = 2)+
  geom_vline(xintercept = median(fin.dat.plot[fin.dat.plot$east.west=='eastern',]$summer.normals, na.rm=T), col = 'black', lty = 2)+
    geom_hline(yintercept=0, colour="white", size=1)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  xlab('Maximum Daily Temperature (July - August) (°C)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16), legend.box.background = element_rect(colour = "black")))
  
ggsave(file = '~/Desktop/fig1b.tiff', fig1b, width = 6, height = 6, dpi = 200)


#####



gbif.summer <- gbif[gbif$countryCode %in% c('US','CA') & gbif$month %in% c(7,8),]

head(gbif.summer)

library(ggmap)

register_google(key = "AIzaSyB5lTW12MLtkc40A3ZPGs5qaVfUcmIsS-M") 

nMap <- get_map("Topeka, Kansas",zoom=4,maptype="satellite",source="google")

ggmap(nMap)+
  theme_bw()+
  geom_point(data = gbif.summer, aes (x = decimalLongitude, y = decimalLatitude),col = 'white',size = 0.5)+
  theme(legend.position = 'none')+
  geom_vline(xintercept = seq(-70, -120, by = -10), col = 'white', size = 0.3)+
  geom_hline(yintercept = seq(20, 60, by = 10), col = 'white', size = 0.3)+
  geom_vline(xintercept = seq(-70, -120, by = -2.5), col = 'white', size = 0.1)+
  geom_hline(yintercept = seq(20, 60, by = 2.5), col = 'white', size = 0.1)


names(gbif.summer)[3:4] <- c('lat','long')

ggmap(nMap)+
  theme_bw()+
  geom_point(data = gbif.summer, aes (x = long, y = lat, col = east.west),size = 0.5)+
  theme(legend.position = 'none')+
  geom_vline(xintercept = seq(-70, -120, by = -10), col = 'white', size = 0.3)+
  geom_hline(yintercept = seq(10, 60, by = 10), col = 'white', size = 0.3)+
  geom_vline(xintercept = seq(-70, -120, by = -2.5), col = 'white', size = 0.1)+
  geom_hline(yintercept = seq(10, 60, by = 2.5), col = 'white', size = 0.1)+
  geom_vline(xintercept = seq(-70, -120, by = -0.5), col = 'white', size = 0.05)+
  geom_hline(yintercept = seq(10, 60, by = 0.5), col = 'white', size = 0.05)+
  scale_color_manual(values = c('yellow','orange'))

gbif.summer$east.west <- ifelse(gbif.summer$long < -105, 'west', 'east')
gbif.summer[which(gbif.summer$long > -113.5 & gbif.summer$long < -105 & gbif.summer$lat >45 & gbif.summer$lat < 60),]$east.west <- 'east'
gbif.summer[which(gbif.summer$long > -114 & gbif.summer$long < -112 & gbif.summer$lat<47 & gbif.summer$lat > 45),]$east.west <- 'west'
gbif.summer[which(gbif.summer$long > -106 & gbif.summer$long < -104 & gbif.summer$lat<41 & gbif.summer$lat >39 ),]$east.west <- 'east'

#### get map of where observations are in january and use this to get rid of non-migratory populations

gbif.jan <- gbif.monarchs[gbif.monarchs$countryCode=='US' & gbif.monarchs$month %in% c(1,2),]
nrow(gbif.jan) #leaves 1106 observations

gbif.jan<-gbif.jan[!is.na(gbif.jan$decimalLatitude),]

gbif.jan$resident <- ifelse(gbif.jan$decimalLatitude < 31, 'resident', 'non-resident')
gbif.jan[gbif.jan$decimalLatitude < 34.5 & gbif.jan$decimalLatitude > 32 & gbif.jan$decimalLongitude < -116.5 & gbif.jan$decimalLongitude > -119 ,]$resident <- 'resident'
gbif.jan[gbif.jan$decimalLatitude > 32 & gbif.jan$decimalLatitude < 35 & gbif.jan$decimalLongitude > -81 & gbif.jan$decimalLongitude < 82.5,]$resident <- 'resident'

ggmap(nMap)+
  theme_bw()+
  geom_point(data = gbif.jan, aes (x = decimalLongitude, y = decimalLatitude, 
                                   col = resident),size = 0.5)+
  theme(legend.position = 'none')+
  geom_vline(xintercept = seq(-70, -120, by = -10), col = 'white', size = 0.3)+
  geom_hline(yintercept = seq(10, 60, by = 10), col = 'white', size = 0.3)+
  geom_vline(xintercept = seq(-70, -120, by = -2.5), col = 'white', size = 0.1)+
  geom_hline(yintercept = seq(10, 60, by = 2.5), col = 'white', size = 0.1)+
  geom_vline(xintercept = seq(-70, -120, by = -0.5), col = 'white', size = 0.05)+
  geom_hline(yintercept = seq(10, 60, by = 0.5), col = 'white', size = 0.05)+
  scale_color_manual(values = c('green','white'))

###now recreate original summer map, but this time using the same filter to specify monarchs from year-round breeding areas

gbif.summer$resident <- ifelse(gbif.summer$lat < 31, 'resident', 'non-resident')
gbif.summer[gbif.summer$lat < 34 & gbif.summer$lat > 32 & gbif.summer$long < -116 & gbif.summer$long > -119 ,]$resident <- 'resident'
gbif.summer[gbif.summer$lat > 32 & gbif.summer$lat < 35 & 
              gbif.summer$long < -75 & gbif.summer$long > -80.8,]$resident <- 'resident'
gbif.summer[gbif.summer$lat < 31 & gbif.summer$lat > 28 & gbif.summer$long < -97 & gbif.summer$long > -105,]$resident <- 'non-resident'

gbif.summer$east.west <- as.factor(gbif.summer$east.west)

levels(gbif.summer$east.west) <- c(levels(gbif.summer$east.west),"resident")

gbif.summer[gbif.summer$resident=='resident',]$east.west <- 'resident'

ggmap(nMap)+
  theme_bw()+
  geom_point(data = gbif.summer, aes (x = long, y = lat, col = east.west),size = 0.5)+
  theme(legend.position = 'none')+
  #geom_vline(xintercept = seq(-70, -120, by = -10), col = 'white', size = 0.3)+
  #geom_hline(yintercept = seq(10, 60, by = 10), col = 'white', size = 0.3)+
  #geom_vline(xintercept = seq(-70, -120, by = -2.5), col = 'white', size = 0.1)+
  #geom_hline(yintercept = seq(10, 60, by = 2.5), col = 'white', size = 0.1)+
  #geom_vline(xintercept = seq(-70, -120, by = -0.5), col = 'white', size = 0.05)+
  #geom_hline(yintercept = seq(10, 60, by = 0.5), col = 'white', size = 0.05)+
  scale_color_manual(values = c('yellow','orange', 'white'))

nMapz3 <- get_map("Topeka, Kansas",zoom=3,maptype="satellite",source="google")

(zz <- ggmap(nMapz3)+
  theme_bw()+
  geom_point(data = gbif.summer, aes (x = long, y = lat, col = east.west),size = 0.3)+
  theme(legend.position = 'none')+
  #geom_vline(xintercept = seq(-30, -160, by = -10), col = 'white', size = 0.3)+
  #geom_hline(yintercept = seq(-20, 80, by = 10), col = 'white', size = 0.3)+
  scale_color_manual(values = c('yellow','orange', 'white')))
  #geom_vline(xintercept = seq(-30, -160, by = -2.5), col = 'white', size = 0.1)+
  #geom_hline(yintercept = seq(-20, 80, by = 2.5), col = 'white', size = 0.1)+
  #geom_vline(xintercept = seq(-30, -160, by = -0.5), col = 'white', size = 0.05)+
  #geom_hline(yintercept = seq(-20, 80, by = 0.5), col = 'white', size = 0.05))

pdf(file = '~/Desktop/summer_monarch_map.pdf', width = 12, height = 12)

zz

dev.off()
  

#ggsave(file="~/Desktop/zz.tiff", plot=zz, width=6, height=6, dpi = 600)


ggmap(nMap)+
  theme_bw()+
  geom_point(data = gbif.summer, aes (x = long, y = lat, col = resident),size = 0.5)+
  theme(legend.position = 'none')+
  geom_vline(xintercept = seq(-70, -120, by = -10), col = 'white', size = 0.3)+
  geom_hline(yintercept = seq(10, 60, by = 10), col = 'white', size = 0.3)+
  geom_vline(xintercept = seq(-70, -120, by = -2.5), col = 'white', size = 0.1)+
  geom_hline(yintercept = seq(10, 60, by = 2.5), col = 'white', size = 0.1)+
  geom_vline(xintercept = seq(-70, -120, by = -0.5), col = 'white', size = 0.05)+
  geom_hline(yintercept = seq(10, 60, by = 0.5), col = 'white', size = 0.05)+
  scale_color_manual(values = c('green','white'))






head(gbif.summer)

r <- getData("worldclim",var="bio",res=2.5)

r2 <- r[[c(1:19)]]

names(r2) <- c("BIO1 = Annual Mean Temperature", "BIO2 = Mean Diurnal Range", "BIO3 = Isothermality",
               "BIO4 = Temperature Seasonality", "BIO5 = Max Temperature of Warmest Month", "BIO6 = Min Temperature of Coldest Month", "BIO7 = Temperature Annual Range (BIO5-BIO6)", "BIO8 = Mean Temperature of Wettest Quarter", "BIO9 = Mean Temperature of Driest Quarter", "BIO10 = Mean Temperature of Warmest Quarter", "BIO11 = Mean Temperature of Coldest Quarter", "BIO12 = Annual Precipitation", "BIO13 = Precipitation of Wettest Month", "BIO14 = Precipitation of Driest Month", "BIO15 = Precipitation Seasonality (Coefficient of Variation)", "BIO16 = Precipitation of Wettest Quarter", "BIO17 = Precipitation of Driest Quarter", "BIO18 = Precipitation of Warmest Quarter", "BIO19 = Precipitation of Coldest Quarter")

lats <- gbif.summer$lat
lons <- gbif.summer$long

coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r2,points)

(values<-as.data.frame(values))

(bioclim.vars <- cbind(gbif.summer, values))

bioclim.vars$east.west <- revalue(bioclim.vars$east.west, c('east' = 'eastern', 'west' = 'western'))

Fig1c <- ggplot(bioclim.vars[bioclim.vars$east.west!='resident',], aes(x = BIO18...Precipitation.of.Warmest.Quarter, col = east.west))+
  geom_density(aes(fill = east.west), adjust = 4, alpha = 0.2, size = 1)+
  theme_bw()+
  scale_color_manual(values = c('yellow','orange'))+
  scale_fill_manual(values = c('yellow','orange'))+
  geom_vline(xintercept = median(bioclim.vars[bioclim.vars$east.west=='western',]$BIO18...Precipitation.of.Warmest.Quarter, na.rm=T), col = 'black', lty = 2)+
  geom_vline(xintercept = median(bioclim.vars[bioclim.vars$east.west=='eastern',]$BIO18...Precipitation.of.Warmest.Quarter, na.rm=T), col = 'black', lty = 2)+
  geom_hline(yintercept=0, colour="white", size=1)+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))+
  xlab('BIO18 - Precipitation of Warmest Quarter')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16), legend.box.background = element_rect(colour = "black"))

Fig1c

ggsave(file = '~/Desktop/fig1c.tiff', Fig1c, width = 6, height = 6, dpi = 200)

