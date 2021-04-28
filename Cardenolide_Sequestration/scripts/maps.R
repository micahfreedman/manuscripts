### maps

library('ggmap')
library('sf')
library("rnaturalearth")
library("rnaturalearthdata")
library('rgeos')
library('rgdal')
library('maps')

world <- ne_countries(scale= 'medium',returnclass = 'sf') #scale = 110 or scale = 'medium' seem to work well; scale = 10 takes a long time

ggplot(data = world)+
  geom_sf()+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  theme(panel.background = element_rect(fill = 'aliceblue'))

plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

#pdf('./figures/Figx1a.pdf', height = 12, width = 24)

plot.map("world", fill=TRUE, col="cornsilk", bg="aliceblue", center = 200, ylim=c(-60, 90), mar=c(0,0,0,0))

#dev.off()

########

#now get map of the Marianas

library(rnaturalearthhires)

new_world <- ne_countries(scale= 10,returnclass = 'sf')

marianas <- subset(new_world, admin == c("Guam", "Northern Mariana Islands"))



(map.marianas <- ggplot(data = marianas) +
  geom_sf(fill = "antiquewhite1")+
  ylim(c(13,15.5))+
  xlim(c(144.4,146))+
  theme(axis.title = element_blank())+
  theme(panel.grid.major = element_line(colour = gray(0.9), linetype = "dashed",size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),panel.border = element_rect(fill = NA))+
  annotate("text", x= 145.2, y= 13.5, label = "Guam", cex = 4)+
  annotate("text", x= 145.5, y= 14.2, label = "Rota", cex = 4)+
  annotate("text", x= 145.6, y= 14.8, label = "Tinian", cex = 4)+
  annotate("text", x= 145.5, y= 15.2, label = "Saipan", cex = 4)+
  annotate("text", x= 144.8, y= 15.5, label = "Mariana Islands", cex = 6))

ggsave(map.marianas, file = './figures/Fig5xa.pdf', height = 8, width = 5)
