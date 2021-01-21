###### Script for getting BioClim data for shrubs sampled in island/mainland comparisons

shrub.bearings = read.csv(file="~/Downloads/shrubs_coordinates.csv",head=T) #load in plant coordinates

head(shrub.bearings) #need to split into separate lat/long columns; first remove the \xa1 that results from degree symbol

shrub.bearings$Position <- gsub('\xa1', '', shrub.bearings$Position) #removes that particular string from all entries; now need to split this text string at the space

split.co <- strsplit(shrub.bearings$Position, " ")

for(i in 1:length(split.co)){
  shrub.bearings$lat[i] <- split.co[[i]][1]
} #loop through elements of split.co to extract latitudes

for(i in 1:length(split.co)){
  shrub.bearings$lon[i] <- split.co[[i]][2]
} #same for longitudes

shrub.bearings$lat <- as.numeric(gsub('N', '', shrub.bearings$lat)) #drop N and convert to numeric
shrub.bearings$lon <- as.numeric(gsub('W', '', shrub.bearings$lon))*-1 #drop W, convert to numeric, multiply by -1 to imply longitude is west

head(shrub.bearings)

plant.spp <- (str_extract(shrub.bearings$Name, "[aA-zZ]+")) #extract text from ID
plant.spp <- revalue(plant.spp, c("CERC"="Cercocarpus","CMEG"="Ceanothus","DEND"="Dendromecon","HET"="Heteromeles","PRUN"="Prunus","QPAC"="Quercus","QBER"="Quercus")) #rename species so that they match

plant.num <- (str_extract(shrub.bearings$Name, "[0-9]+")) #extract numbers from ID
plant.ID <- paste(plant.spp, plant.num, sep = "_")

shrub.bearings$plant.ID <- plant.ID
shrub.bearings$species <- plant.spp

leaf.data.spatial <- merge(shrub.bearings, raw.dat, by = "plant.ID")

df1=aggregate(lat ~ plant.ID + Site, data = leaf.data.spatial, FUN = mean)
df2=aggregate(lon ~ plant.ID + Site, data = leaf.data.spatial, FUN = mean)

plant.mapping <- merge(df1,df2, by = c("plant.ID","Site"))
spp2 <- (str_extract(plant.mapping$plant.ID, "[aA-zZ]+"))
plant.mapping$species <- gsub('_', '', spp2)

#####
#now need to convert into a spatially explicit dataframe
#####

plant.mapping.sp <- SpatialPointsDataFrame(plant.mapping[,c(3:4)], plant.mapping[,-c(3:4)])

crs.geo<-CRS("+init=epsg:4326")
proj4string(plant.mapping.sp) <- crs.geo #this is doing the same things as the "coordinates" function and is just combining lat and lon into a readable single entity for spatial mapping objects

is.projected(plant.mapping.sp) #points now are indeed projections


##### first, add in bioclim variables

library(raster)
library(sp)

r <- getData("worldclim",var="bio",res=2.5)

r2 <- r[[c(1:19)]]

names(r2) <- c("BIO1 = Annual Mean Temperature", "BIO2 = Mean Diurnal Range", "BIO3 = Isothermality",
               "BIO4 = Temperature Seasonality", "BIO5 = Max Temperature of Warmest Month", "BIO6 = Min Temperature of Coldest Month", "BIO7 = Temperature Annual Range (BIO5-BIO6)", "BIO8 = Mean Temperature of Wettest Quarter", "BIO9 = Mean Temperature of Driest Quarter", "BIO10 = Mean Temperature of Warmest Quarter", "BIO11 = Mean Temperature of Coldest Quarter", "BIO12 = Annual Precipitation", "BIO13 = Precipitation of Wettest Month", "BIO14 = Precipitation of Driest Month", "BIO15 = Precipitation Seasonality (Coefficient of Variation)", "BIO16 = Precipitation of Wettest Quarter", "BIO17 = Precipitation of Driest Quarter", "BIO18 = Precipitation of Warmest Quarter", "BIO19 = Precipitation of Coldest Quarter")

#specify latitudes and longitudes for each sampled point to extract climate data

lats <- plant.mapping$lat
lons <- plant.mapping$lon

coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r2,points)
(values<-as.data.frame(values))

(bioclim.vars <- cbind(plant.mapping, values))

#### get prinicipal components of climatic data

drop.columns.pca.clim <- c("plant.ID", "Site", "lat", "lon", "species") #select columns to exlude from analysis
bioclim.pca <- bioclim.vars[,!(names(bioclim.vars) %in% drop.columns.pca.clim)] #drop columns with the aforementioned names

(fit.clim <- princomp(bioclim.pca))

pc.sum <- Reduce("+",fit.clim[[1]]) #get sum of PC scores

fit.clim[[1]][1] / pc.sum #first PC explains 83.6% of variance
fit.clim[[1]][2] / pc.sum #second PC explains 13.3% of variance

library(factoextra)

fviz_pca_var(fit.clim, col.var = "black")

bioclim.vars$PC1 <- fit.clim$scores[,1] #extract first principal component and append to original df
bioclim.vars$PC2 <- fit.clim$scores[,2] #same for second PC

#add island/mainland column

islands <- c("Santa Rosa", "Santa Cruz", "Catalina")

bioclim.vars$IM <- ifelse(bioclim.vars$Site %in% islands, "Island", "Mainland")

ggplot(bioclim.vars, aes(x = PC1, y = PC2*-1, col = Site))+
  geom_jitter(size = 5)+
  theme_dark()+
  xlab('PC1 (83.6%) - Temperature Seasonality')+
  ylab('PC2 (13.3%) - Annual Precipitation')+
  scale_color_manual(values = c('yellow','red','gold1', 'tomato1','khaki1','pink'))+
  theme(legend.position = 'none')+
  theme(axis.text = element_blank())+
  annotate('text',x=-300,y=-175,label='Catalina',color = 'yellow',size=5)+
  annotate('text',x=-625,y=90,label='Santa Cruz',color = 'gold1',size=5, angle = 60)+
  annotate('text',x=-725,y=-25,label='Santa Rosa',color = 'khaki1',size=5)+
  annotate('text',x=25,y=75,label='Gaviota',color = 'red',size=5, angle = 60)+
  annotate('text',x=500,y=50,label='Santa Monicas',color = 'tomato1',size=5)+
  annotate('text',x=900,y=0,label='Stunt Ranch',color = 'pink',size=5)

leaf.data <- merge(bioclim.vars[,c(1,3:4,25:26)], raw.dat, by = "plant.ID")