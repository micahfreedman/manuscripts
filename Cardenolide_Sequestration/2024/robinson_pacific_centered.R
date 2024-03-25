library(rnaturalearth)
library(ggplot2)
library(sf)
library(terra)

sf_use_s2(FALSE)

world_map = ne_countries(scale = "medium", returnclass = "sf", type = "countries") #background map 

lakes <- ne_download(scale = "small", type = "lakes", category = "physical", returnclass = 'sf') #include a layer for major lakes

robinson = st_crs("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", )

world_map_transform <- st_transform(world_map, robinson) #project world map boundaries
lakes_transform <- st_transform(lakes, robinson) #project lakes object

########

crsrobin <- "+proj=robin +lon_0=-180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world_robinson2 <- world_map %>%
  st_break_antimeridian(lon_0 = 180) %>% # insert this before transformation
  # (& ignore the warning)
  st_transform(crs = crsrobin)

ggplot() +
  geom_sf(data = world_robinson2, fill = "lemonchiffon", color = "gray40")+
  geom_sf(data = lakes_transform, fill = "lightcyan", color = "gray40")+
  theme(panel.background = element_rect(fill='lightcyan'))+
  coord_sf(xlim = c(-8000000, 12000000), ylim = c(-4800000, 5700000))

##########

crs_mollweide <- "+proj=moll +lon_0=-180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world_mollweide <- world_map %>%
  st_break_antimeridian(lon_0 = 180) %>% # insert this before transformation
  # (& ignore the warning)
  st_transform(crs = crs_mollweide)

st_bbox(world_mollweide)

ggplot() +
  geom_sf(data = world_mollweide, fill = "lemonchiffon", color = "gray40")+
  geom_sf(data = lakes_transform, fill = "lightcyan", color = "gray40")+
  theme(panel.background = element_rect(fill='lightcyan'))+
  coord_sf(xlim = c(-6000000, 12000000), ylim = c(-4800000, 6000000))

##########
