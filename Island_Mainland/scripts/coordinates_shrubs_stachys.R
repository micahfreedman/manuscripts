#### Script for creating mapping data for Channel Islands / mainland comparisons

setwd('~/Documents/GitHub/Sites/manuscripts/Island_Mainland/')

#load relevant libraries (some might not actually be needed?)

library(rworldmap)
library(ggmap)
library(maps)
library(mapdata)
library(stringr)
library(rworldxtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgdal)
library(ggsn)
library(cowplot)
library(ggspatial)

### first do zoomed out map of California 

#register_google(key = "api key here-M")  #one-time API key entry, required for using ggmap

sat_california <- get_map(location = 'California', zoom = 6, maptype = 'hybrid')

#pdf(file = './figures/Fig1a1.pdf', height = 6, width = 6)  
ggmap(sat_california)+
  annotate("rect", xmin = -121, xmax = -117.5, ymin = 32.6, ymax = 35, 
           col = 'red', fill = 'white', alpha = 0.3)+ #note: bounding box is just approximate
  theme(axis.title = element_blank(), axis.text = element_blank())
#dev.off()

world <- ne_countries(scale='large',returnclass = 'sf')
states <- ne_download(scale = 'medium', type = "states", returnclass = 'sf')
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = 'sf') #include a layer for major lakes


#pdf(file = './figures/Fig1a1.pdf', height = 4, width = 6)
ggplot(data = world)+
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = states, fill = 'antiquewhite1')+
  geom_sf(data = lakes, fill = "lightcyan", color = "gray40")+
  theme(panel.grid.major = element_line(colour = 'white', linetype = "dashed",size = 0.5), 
        panel.background = element_rect(fill = "lightcyan"),panel.border = element_rect(fill = NA))+
  coord_sf(xlim = c(-130, -110), ylim = c(25, 45)) +
  annotate("rect", xmin = -121, xmax = -117, ymin = 32.5, ymax = 35, col = 'red', fill = 'white', alpha = 0.2)
#dev.off()  

sites<-read.csv("./data_files/site_coordinates.csv",row.names=NULL,header = T) #read in coordinates for site-level locations

sites$IM <- ifelse(sites$Site %in% c('Catalina','Santa_Cruz','Santa_Rosa'), 'island', 
                   ifelse(sites$Site %in% c('Rancho_Santa_Ana','Santa_Barbara'), 'common garden', 'mainland'))

#pdf(file = './figures/Fig1a2.pdf', height = 4, width = 6)
ggplot(data = world) +
  geom_sf(fill = "antiquewhite1") +
  geom_point(data = sites, aes(x = -long, y = lat, col = IM), size = 4, alpha = 0.5)+
  scale_color_manual(values = c('black','blue','red'))+
    coord_sf(xlim = c(-121, -117), ylim = c(32.5, 35)) +
    theme(panel.grid.major = element_line(colour = gray(0.9), linetype = "dashed",size = 0.5), 
          panel.background = element_rect(fill = "lightcyan"),panel.border = element_rect(fill = NA))+
  annotation_scale(line_width = 0.5)+
  theme(legend.position = 'none', axis.title = element_blank())+
  annotate("text", x= -118.9, y= 33.3, label = "Catalina", cex = 4, col = 'blue')+
  annotate("text", x= -119.3, y= 33.85, label = "Santa Cruz", cex = 4, col = 'blue')+
  annotate("text", x = -120.5, y = 34.2, label = "Santa Rosa", cex = 4, col = 'blue')+
  annotate("text", x = -120.2, y = 34.7, label = "Gaviota", cex = 4, col = 'red')+
  annotate("segment", x = -119.2, xend = -118.81, y = 34.4, yend = 34.075, col = 'red')+
  annotate("text", x = -118.7, y = 34.5, label = "Santa Monica Mtns.", cex = 4, col = 'red')+
  annotate("segment", x = -118.65, xend = -118.3, y = 34.092, yend = 34.2, col = 'red')+
  annotate("text", x = -117.9, y = 34.3, label = "Stunt Ranch", cex = 4, col = 'red')+
  annotate("text", x = -119.2, y = 34.8, label = "Santa Barbara\nBotanic Garden", cex = 4, col = 'black')+
  annotate("segment", x = -119.7, xend = -119.5, y = 34.45, yend = 34.65, col = 'black')+
  annotate("text", x = -117.5, y = 33.9, label = "Rancho Santa Ana\nBotanic Garden", cex = 4, col = 'black')
#dev.off()

####################

############ site maps ###########

#### now create separate maps for each individual site

shrub.bearings <- read.csv(file="./data_files/shrubs_coordinates.csv",head=T) #load in plant coordinates
leaf.data <- read.csv(file="./data_files/chaparral_leaf_morphology.csv") #also need to load in main dataset to get info on island ID for each plant

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

leaf.data$plant.ID <- paste(str_extract(leaf.data$Species, "[^_]+"), leaf.data$Plant, sep = "_")

leaf.data.spatial <- merge(shrub.bearings, leaf.data, by = "plant.ID")

df1=aggregate(lat ~ plant.ID + Site, data = leaf.data.spatial, FUN = mean)
df2=aggregate(lon ~ plant.ID + Site, data = leaf.data.spatial, FUN = mean)

plant.mapping <- merge(df1,df2, by = c("plant.ID","Site"))
spp2 <- (str_extract(plant.mapping$plant.ID, "[aA-zZ]+"))
plant.mapping$species <- gsub('_', '', spp2)


#get map of all plants from Catalina

catalina.map <- get_map(location = c(-118.45,33.4), maptype = "satellite", source = "google", zoom = 11) #just manually enter points based on mean location of sampled plants

table(with(plant.mapping[plant.mapping$Site == "Catalina",], species)) #all five species present

#pdf('./figures/FigS1a.pdf', height = 6, width = 6)
(catalina.final <- ggmap(catalina.map)+
    geom_jitter(data = plant.mapping[plant.mapping$Site == "Catalina",], aes(x = lon, y = lat), 
                pch = 21, size = 2, col = 'black', fill = 'yellow')+
    theme(axis.title = element_blank(), axis.text = element_text(size = 8))+
    theme(legend.position = "none")+
    ggtitle('Catalina')+
    theme(plot.title = element_text(hjust = 0.5, size = 18))+
  scalebar(x.min = -118.7, x.max = -118.5,
                 y.min = 33.28, y.max = 33.35, 
                 dist = 3, dist_unit = "km", model = 'WGS84', 
                 transform = TRUE, box.fill = c('orange','white'), height = 0.05, 
           st.dist = 0.2, st.color = 'white'))
#dev.off()
  
#now get maps for other locations

#average lat and lon for each site for selecting coordinates
aggregate(lat ~ Site, data = plant.mapping, FUN = mean)
aggregate(lon ~ Site, data = plant.mapping, FUN = mean)

### Stunt Ranch

stunt.map <- get_map(location = c(-118.655,34.0925), maptype = "satellite", source = "google", zoom = 15)

#pdf('./figures/FigS1e.pdf', height = 6, width = 6)
(stunt.final <- ggmap(stunt.map)+
    geom_jitter(data = plant.mapping[plant.mapping$Site == "Stunt_Ranch",], aes(x = lon, y = lat), 
                fill = 'yellow', pch = 21, size = 2)+
    theme(axis.title = element_blank())+
    theme(axis.text = element_text(size = 8))+
    theme(legend.position = 'none')+
    ggtitle('Stunt Ranch')+
    theme(plot.title = element_text(hjust = 0.5, size = 18))+
    scalebar(x.min = -118.665, x.max = -118.65,
             y.min = 34.085, y.max = 34.09, 
             dist = 0.5, dist_unit = "km", model = 'WGS84', 
             transform = TRUE, box.fill = c('orange','white'), height = 0.05, 
             st.dist = 0.2, st.color = 'white'))
#dev.off()

### Gaviota

gav.map <- get_map(location = c(-120.208,34.5), maptype = "satellite", source = "google", zoom = 15)

#pdf('./figures/FigS1d.pdf', height = 6, width = 6)
(gav.final <- ggmap(gav.map)+
    geom_jitter(data = plant.mapping[plant.mapping$Site == "Gaviota",], aes(x = lon, y = lat), 
                fill = 'yellow', pch = 21, size = 2)+
    theme(axis.title = element_blank())+
    theme(axis.text = element_text(size = 8))+
    theme(legend.position = 'none')+
    ggtitle('Gaviota')+
    theme(plot.title = element_text(hjust = 0.5, size = 18))+
    scalebar(x.min = -120.225, x.max = -120.205,
             y.min = 34.492, y.max = 34.495, 
             dist = 0.5, dist_unit = "km", model = 'WGS84', 
             transform = TRUE, box.fill = c('orange','white'), height = 0.08, 
             st.dist = 0.2, st.color = 'white'))
dev.off()

### Santa Cruz Island

cruz.map <- get_map(location = c(-119.73,34), maptype = "satellite", source = "google", zoom = 11)

#pdf('./figures/FigS1b.pdf', height = 6, width = 6)
(cruz.final <- ggmap(cruz.map)+
    geom_jitter(data = plant.mapping[plant.mapping$Site == "Santa_Cruz",], aes(x = lon, y = lat), 
                fill = 'yellow', pch = 21, size = 2)+
    theme(axis.title = element_blank())+
    theme(axis.text = element_text(size = 8))+
    theme(legend.position = 'none')+
    ggtitle('Santa Cruz')+
    theme(plot.title = element_text(hjust = 0.5, size = 18))+
    scalebar(x.min = -119.9, x.max = -119.8,
             y.min = 33.85, y.max = 33.9, 
             dist = 3, dist_unit = "km", model = 'WGS84', 
             transform = TRUE, box.fill = c('orange','white'), height = 0.08, 
             st.dist = 0.2, st.color = 'white'))
#dev.off()

### Santa Monica Mountains

smm.map <- get_map(location = c(-118.8,34.075), maptype = "satellite", source = "google", zoom = 14)

#pdf('./figures/FigS1f.pdf', height = 6, width = 6)
(smm.final <- ggmap(smm.map)+
    geom_jitter(data = plant.mapping[plant.mapping$Site == "Santa_Monicas",], aes(x = lon, y = lat), 
                fill = 'yellow', pch = 21, size = 2)+
    theme(axis.title = element_blank())+
    theme(axis.text = element_text(size = 8))+
    theme(legend.position = 'none')+
    ggtitle('Santa Monica Mtns.')+
    theme(plot.title = element_text(hjust = 0.5, size = 18))+
    scalebar(x.min = -118.82, x.max = -118.8,
             y.min = 34.06, y.max = 34.07, 
             dist = 0.5, dist_unit = "km", model = 'WGS84', 
             transform = TRUE, box.fill = c('orange','white'), height = 0.08, 
             st.dist = 0.2, st.color = 'white'))
#dev.off()

### Santa Rosa

rosa.map <- get_map(location = c(-120.07,33.985), maptype = "satellite", source = "google", zoom = 12)

#pdf('./figures/FigS1c.pdf', height = 6, width = 6)
(rosa.final <- ggmap(rosa.map)+
    geom_jitter(data = plant.mapping[plant.mapping$Site == "Santa_Rosa",], aes(x = lon, y = lat), 
                fill = 'yellow', pch = 21, size = 2)+
    theme(axis.title = element_blank())+
    theme(axis.text = element_text(size = 8))+
    theme(legend.position = 'none')+
    ggtitle('Santa Rosa')+
    theme(plot.title = element_text(hjust = 0.5, size = 18))+
    scalebar(x.min = -120.15, x.max = -120.1,
             y.min = 33.92, y.max = 33.95, 
             dist = 2, dist_unit = "km", model = 'WGS84', 
             transform = TRUE, box.fill = c('orange','white'), height = 0.08, 
             st.dist = 0.2, st.color = 'white'))
#dev.off()

### attempt to plot all together instead of exporting separately

pdf(file = "./figures/Supplemental_Figures/Figure.S1.pdf", height = 12, width = 16)

(site_maps <- plot_grid(catalina.final, cruz.final, rosa.final, 
          gav.final, stunt.final, smm.final, 
          nrow = 2, align = 'h')) #doesn't render quite right
dev.off()


###### Finally, get map for Stachys sampling locations

stachys_coordinates <- read.csv('./data_files/stachys_coordinates.csv')

pdf('./figures/FigSB_map.pdf', height = 6, width = 6)
ggplot(data = world) +
  geom_sf(fill = "antiquewhite1") +
  geom_point(data = stachys_coordinates, aes(x = -long, y = lat, fill = Site), size = 4, pch = 21)+
  scale_fill_manual(values = c('orange','darkorange4','black','blue','red','darkblue','gold'))+
  coord_sf(xlim = c(-121, -118), ylim = c(33, 35)) +
  theme(panel.grid.major = element_line(colour = gray(0.9), linetype = "dashed",size = 0.5), 
        panel.background = element_rect(fill = "lightcyan"),panel.border = element_rect(fill = NA))+
  annotation_scale(line_width = 0.5)+
  theme(legend.position = 'none', axis.title = element_blank())+
  annotate("text", x= -119.3, y= 33.85, label = "Santa Cruz", cex = 4, col = 'blue')+
  annotate("text", x = -120.5, y = 33.8, label = "Santa Rosa", cex = 4, col = 'darkblue')+
  annotate("text", x = -120.3, y = 34.7, label = "Gaviota", cex = 4, col = 'darkorange4')+
  annotate("text", x = -118.6, y = 34.3, label = "Santa Monica Mtns.", cex = 4, col = 'red')+
  annotate("text", x = -120.2, y = 34.35, label = "El Capitan", cex = 4, col = 'orange')+
  annotate("text", x = -118.8, y = 33.95, label = "Zuma", cex = 4, col = 'gold')+
  annotate("text", x = -119.2, y = 34.66, label = "Santa Barbara\nBotanic Garden", cex = 4, col = 'black')
dev.off()
  



ggplot(data = usa) +
  geom_sf(fill = "antiquewhite1") +
  geom_point(data = stachys_coordinates, aes(x = -long, y = lat, fill = Site), size = 4, pch = 21)+
  scale_fill_manual(values = c('orange','darkorange4','black','blue','red','darkblue','gold'))+
  coord_sf(xlim = c(-121, -118), ylim = c(33, 35)) +
  theme(panel.grid.major = element_line(colour = gray(0.9), linetype = "dashed",size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),panel.border = element_rect(fill = NA))+
  theme(legend.position = 'none', axis.title = element_blank())+
  annotate("text", x= -119.3, y= 33.85, label = "Santa Cruz", cex = 4, col = 'blue')+
  annotate("text", x = -120.5, y = 33.8, label = "Santa Rosa", cex = 4, col = 'darkblue')+
  annotate("text", x = -120.3, y = 34.7, label = "Gaviota", cex = 4, col = 'darkorange4')+
  annotate("text", x = -118.6, y = 34.3, label = "Santa Monica Mtns.", cex = 4, col = 'red')+
  annotate("text", x = -120.2, y = 34.35, label = "El Capitan", cex = 4, col = 'orange')+
  annotate("text", x = -118.8, y = 33.95, label = "Zuma", cex = 4, col = 'gold')+
  annotate("text", x = -119.2, y = 34.66, label = "Santa Barbara\nBotanic Garden", cex = 4, col = 'black')
dev.off()

##########################################