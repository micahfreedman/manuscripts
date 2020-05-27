#### Monarch global wing morphology project

library(ggplot2)
library(lme4)
library(lmerTest)
library(plyr)
library(ggmap)
library(mapproj)
library(car)
library(emmeans)

## read in most up-to-date data

setwd('~/Documents/GitHub/Sites/manuscripts/Freedman_et_al_monarch_global_wing_morphology/data_and_analysis/')

wings = read.csv(file="./wings_04.25.20.csv",header=TRUE)

## reformat dates into separate columns for day, month, and year of collection

date.split <- t(sapply(wings$Collection_Date, function(x) substring(x, first=c(1,5,7), last=c(4,6,8))))
wings <- cbind(wings, (date.split))

names(wings)[names(wings) == c('1','2','3')] <- c('year','month','day')

#change year back into a continuous numeric variable
wings$year = as.numeric(levels(wings$year))[wings$year]
#do the same for day
wings$day = as.numeric(levels(wings$day))[wings$day]
#replace 0s with NA and get rid of any values greater than 31
wings$day[wings$day == 0] = NA
wings$day[wings$day >= 32] = NA

wings$month <- revalue(wings$month, c('00' = NA))
wings$month <- month.abb[wings$month]
wings$month
names(wings)

levels(wings$Collection)

wings$Region <- revalue(wings$Region, c("North_America" = "North America", "Central_America"="Central America", "South_America"="South America", "Pacific_Islands"="Pacific Islands")) #rename levels for regions of interest

## Add average columns

wings$MLength = rowMeans(cbind(wings$LLength, wings$RLength),na.rm=T)
wings$MWidth = rowMeans(cbind(wings$LWidth, wings$RWidth),na.rm=T)
wings$MArea = rowMeans(cbind(wings$LArea, wings$RArea),na.rm=T)
wings$MPerimeter = rowMeans(cbind(wings$LPerimeter,wings$RPerimeter),na.rm=T)
wings$aspectratio = wings$MLength/wings$MWidth
wings$roundness = (4*pi*wings$MArea)/(wings$MPerimeter)^2

names(wings)

### Add geocode column ###

for (i in 1:nrow(wings))
  wings$geocode_location_full[i]<-paste(wings$Country.Archipelago[i], wings$Island.State[i], wings$Site.City[i], sep=" ")

##Geocoding step

wings$geocode_location_full

#register_google(key = "mykey") #do not need to run here; this is just to demonstrate how coordinates were generated in the first place

#coordinates <- geocode(wings$geocode_location_full)
#wings$lat <- coordinates[[2]]
#wings$lon <- coordinates[[1]] #these values are already included in the loaded dataset

#####

#omit reared specimens and specimens from David_M collection

wings <- wings[wings$Wild.Caught. =='yes',]

wings <- wings[wings$Collection != 'David_M',] #this particular collection was photographed in a way that made individual monarchs look abnormally large (10-15% larger than normal); omit for this reason

#excluse sole gyandromorph from data

wings = subset(wings, wings$Sex != "gyandromorph")

####add column for julian dates

date<-as.Date(as.character(wings$Collection_Date), "%Y%m%d")

wings$julian.date <- as.numeric(format(date, "%j"))

###############################

### Now create map; begin by excluding ambiguous cases

#### create filter to exclude non-migratory North American specimens 

winter.months <- c('Dec','Jan','Feb')

#first, filter any North American specimen collected below 30N, with the exception of overwintering monarchs

ambiguous1 <- with(wings, wings[which(Region == 'North America' & lat < 30),]) #180 records

#filter any U.S. records outside of CA collected in Dec., Jan., Feb.

ambiguous2 <- with(wings, wings[which(Country.Archipelago == 'United_States' & Island.State != 'California' & month %in% winter.months),]) #only records come from Florida, which is good; includes 22 records

#filter any southern California records from LA basin and San Diego

ambiguous3 <- with(wings, wings[which(Island.State %in% c('California','Arizona') & lat < 34.1),]) #216 records

#filter any records from Cuba and the Bahamas from October - March

ambiguous4 <- with(wings, wings[which(Country.Archipelago %in% c('Cuba','Bahamas') & month %in% c('Oct','Nov',winter.months,'Mar')),]) #19 records

#filter any records from Mexican locations considered to be Central America between October - March

ambiguous5 <- with(wings, wings[which(Region == 'Central America' & Country.Archipelago == 'Mexico' & month %in% c('Oct','Nov',winter.months,'Mar')),]) #75 records

ambiguous <- rbind(ambiguous1,ambiguous2,ambiguous3,ambiguous4,ambiguous5)

nrow(ambiguous) #512 total observations

### now, reclassify these locations as "ambiguous" under Region; for mapping, show these locations in white; also classify British transient specimens as their own region, shown in black

levels(wings$Region) <- c(levels(wings$Region), 'Ambiguous','British Isles','Pacific Ambiguous')

wings[which(wings$Index %in% ambiguous$Index),]$Region <- 'Ambiguous'
wings[which(wings$Country.Archipelago == 'British_Isles'),]$Region <- 'British Isles'
wings[which(wings$Country.Archipelago %in% c('Australia','New_Zealand') & wings$lat < -28),]$Region <- 'Pacific Ambiguous'
wings[which(wings$Island.State == 'Western_Australia'),]$Region <- 'Pacific Islands'

table(wings$Region)

### prepare map of global monarch collections

library('sf')
library("rnaturalearth")
library("rnaturalearthdata")
library('rgeos')
library('rgdal')

world <- ne_countries(scale= 'medium',returnclass = 'sf') #scale = 110 or scale = 'medium' seem to work well; scale = 10 takes a long time

(fig1.raw <- ggplot(data = world)+
  geom_sf()+
  theme_bw()+
  geom_point(data = wings[!(wings$Region %in% c('Asia',"Indian_Ocean")),],
             aes(x = lon, y = lat, col = Region), size = 0.8)+
  theme(legend.position = 'none')+
  theme(axis.title = element_blank())+
  scale_color_manual(values = c('purple', 'orange','gold','green4','blue','brown1','lawngreen','black','cyan'))+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  theme(panel.background = element_rect(fill = 'aliceblue')))

## actual Figure 1, including table of specimen counts, made in Keynote

####

# group Canary Islands and Madeira with Morocco, based on proximity

wings[which(wings$Country.Archipelago %in% c('Madeira','Canary_Islands')),]$Country.Archipelago = 'Morocco'

#########################


################### generate principle components for size based on wing length, wing area, width width; first, drop all rows with NA entries for these measures

wings.pca.size <- na.omit(wings[,c("Index",'MLength','MWidth','MArea')])

fit.size = princomp(wings.pca.size[,c('MLength','MWidth','MArea')])

summary(fit.size)
names(fit.size)
head(fit.size$scores)

#first two principle components account for 95.6% of the total variance for the 3 wing size related variables

PC1 = fit.size$scores[,1]
wings.pca.size$Size.PC1 = PC1

#merge size PCA with original wings dataframe

wings = merge(wings, wings.pca.size, by = intersect(names(wings), names(wings.pca.size)))

wings$year.scaled = scale(wings$year)

##### plot of size principle component by region

wings.plot <- wings

#first set North America as the "reference level" to compare all others to
wings.plot$Region = factor(wings.plot$Region, levels=c("North America","Central America","South America","Caribbean","Pacific Islands","Atlantic"))
wings.plot$Region <- droplevels(wings.plot$Region)

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

ggplot(wings.plot[!is.na(wings.plot$Region),], aes(x=Region,y=Size.PC1))+
  geom_point(aes(col = Region), position = position_jitterdodge(3.2), size = 1)+
  geom_boxplot(outlier.size = 0, outlier.color = 'white', alpha = 0, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 1)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank())+
  scale_color_manual(name="Region",values=c("limegreen","gold","brown1","orange","blue",'purple'))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Size PC1")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16))

######

#### Run model to compare size PC between regions, including change through time

model1.size = lmer(Size.PC1 ~ Region*year.scaled + Sex + (1|Collection) + (1|Country.Archipelago) + smoothing + image_type  + (1|observer) + scale(abs(lat)), data = wings.plot) #here, run a model that uses only butterflies from the regions of interest and omits individuals whose migratory status is ambiguous

summary(model1.size)
ranef(model1.size)

####

model_all.size = lmer(Size.PC1 ~ Region*year.scaled + Sex + (1|Collection) + (1|Country.Archipelago) + smoothing + image_type  + (1|observer) + scale(abs(lat)), data = wings) #for this model, include ambiguous specimens

summary(model_all.size)

plot.trends.size <- data.frame(emtrends(model1.size, "Region", var="year.scaled", data = wings.plot)) #pbkrtest.limit = 5000)) #get marginal means using the model 

plot.trends.size$t <- plot.trends.size$year.scaled.trend / plot.trends.size$SE
plot.trends.size$p <- 2* pt(abs(plot.trends.size$t), df = plot.trends.size$df, lower.tail = FALSE)

plot.trends.size$Region <- factor(plot.trends.size$Region, c('Atlantic','Pacific Islands','Caribbean','South America','Central America','North America'))

ggplot(plot.trends.size[!is.na(plot.trends.size$Region),], aes(x = year.scaled.trend, y = Region, col = Region))+
  theme_classic(base_size = 16)+
  geom_point(size = 3)+
  geom_errorbarh(aes(xmax = asymp.UCL, xmin = asymp.LCL),
                 size = 1, height = 0.2)+
  scale_color_manual(values = c('purple','blue','orange','brown1','gold','limegreen'))+
  theme(legend.position = 'none')+
  xlab('Slope Estimate: Wing Size')+
  geom_vline(xintercept = 0, lty = 2)+ 
  theme(axis.title.y = element_blank()) #looks fine, shows clear difference in trajectory between recently established non-migratory populations and other populations

plot.trends.size.all <- data.frame(emtrends(model_all.size, "Region", var="year.scaled", data = wings)) #can edit this line of code to add argument pbkrtest.limit = 6000, and it will then return degrees of freedom estimates; takes about 15 minutes to run doing this and is not important for interpretation here

ggplot(plot.trends.size.all[plot.trends.size.all$Region != 'British Isles',], aes(x = year.scaled.trend, y = Region, col = Region))+
  theme_classic(base_size = 20)+
  geom_point(size = 3)+
  geom_errorbarh(aes(xmax = asymp.UCL, xmin = asymp.LCL),
                 size = 1, height = 0.2)+
  theme(legend.position = 'none')+
  xlab('Slope Estimate: Wing Size')+
  geom_vline(xintercept = 0, lty = 2)+ 
  theme(axis.title.y = element_blank()) #reinforces the point that Atlantic and Pacific truly are outliers; only negative slopes even when ambiguous cases from Pacific and North America are separately included

############

#now do the same thing for shape principal component

################### generate principle components for size based on wing length, wing area, width width; first, drop all rows with NA entries for these measures

wings.pca.shape <- na.omit(wings[,c("Index",'roundness','aspectratio')])

fit.shape = princomp(wings.pca.shape[,c('roundness','aspectratio')])

summary(fit.shape)

PC1 = fit.shape$scores[,1]
wings.pca.shape$Shape.PC1 = PC1*-10 #here multiply this value by 10 just so the values are on a comparable scale to those for size PC1

###
#merge shape PCA with original wings dataframe

wings = merge(wings, wings.pca.shape, by = intersect(names(wings), names(wings.pca.shape)))

##### plot of size principle component by region

wings.plot <- wings

#first set North America as the "reference level" to compare all others to
wings.plot$Region = factor(wings.plot$Region, levels=c("North America","Central America","South America","Caribbean","Pacific Islands","Atlantic"))
wings.plot$Region <- droplevels(wings.plot$Region)

ggplot(wings.plot[!is.na(wings.plot$Region) & wings.plot$Shape.PC1 <4,], aes(x=Region,y=Shape.PC1))+
  geom_point(aes(col = Region), position = position_jitterdodge(3.2), size = 1)+
  geom_boxplot(outlier.size = 0, outlier.color = 'white', alpha = 0.6,outlier.alpha=0,size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank(), )+
  scale_color_manual(name="Region",values=c("limegreen","gold","brown1","orange","blue",'purple'))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Shape PC1")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16)) #lots of elongated outliers among Pacific populations

model1.shape = lmer(Shape.PC1 ~ Region*year.scaled + Sex + (1|Country.Archipelago) + scale(abs(lat)) + smoothing + image_type + (1|observer) + (smoothing|Collection), data = wings.plot)
summary(model1.shape)



model_all.shape = lmer(Shape.PC1 ~ Region*year.scaled + Sex + (1|Country.Archipelago) + scale(abs(lat)) + smoothing + image_type + (1|observer) + (smoothing|Collection), data = wings) #same as before; run second model that includes ambiguous migratory status individuals
summary(model_all.shape) 

wing.shape.slopes = data.frame(emtrends(model_all.shape, "Region", var="year.scaled", data = wings))

wing.shape.slopes$Region <- factor(wing.shape.slopes$Region, c('Atlantic','Pacific Islands','Caribbean','South America','Central America','North America'))

ggplot(wing.shape.slopes[!is.na(wing.shape.slopes$Region),], aes(x = year.scaled.trend, y = Region, col = Region))+
  theme_classic(base_size = 16)+
  geom_point(size = 3)+
  geom_errorbarh(aes(xmax = asymp.UCL, xmin = asymp.LCL),
                 size = 1, height = 0.2)+
  scale_color_manual(values = c('purple','blue','orange','brown1','gold','limegreen'))+
  theme(legend.position = 'none')+
  xlab('Slope Estimate: Wing Shape')+
  geom_vline(xintercept = 0, lty = 2)+
  theme(axis.title.y = element_blank())


wing.shape.slopes$t <- wing.shape.slopes$year.scaled.trend / wing.shape.slopes$SE
wing.shape.slopes$p <- 2* pt(abs(wing.shape.slopes$t), df = wing.shape.slopes$df, lower.tail = FALSE)

wing.shape.slopes

###now make figure 3, showing all observations along size and shape axes

wings.plot$Region = factor(wings.plot$Region, levels=c("North America","Central America","South America","Caribbean","Pacific Islands","Atlantic"))
wings.plot$Region <- droplevels(wings.plot$Region)

ggplot(wings.plot[!is.na(wings.plot$Region) & wings.plot$Shape.PC1 <8,], aes(x=Size.PC1,y=Shape.PC1))+
  geom_jitter(aes(col = Region), size = 0.5, alpha = 0.3)+
  theme_classic()+
  stat_ellipse(aes(col = Region), size = 1, lty = 1, level = 0.95)+
  theme(axis.text = element_text(size = 16))+
  theme(axis.title = element_text(size = 16))+
  theme(legend.title=element_blank(), )+
  scale_color_manual(name="Region",values=c("limegreen","gold","brown1","orange","blue",'purple'))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  xlab('Size PC1')+
  ylab("Shape PC1")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  theme(legend.justification = 'top')+
  ylim(c(-2,3)) #note: for purposes of display, dataset does not include extreme shape outliers from the Pacific

#############################

## plot trend through time for Pacific, Atlantic, and North America

#################

North_America <- subset(wings, wings$Region == 'North America')
Caribbean <- subset(wings, wings$Region == 'Caribbean')

#subset out only Marianas, Samoa, Hawaii, New Caledonia, Australia, PNG

presentation1 <- wings[wings$Region=='Pacific Islands' & wings$Country.Archipelago %in% c('Hawaii','Samoa','New_Caledonia','New_Guinea','Australia', 'Marianas', 'Fiji'),]

presentation1$Country.Archipelago <- revalue(presentation1$Country.Archipelago, c("New_Caledonia" = "New Caledonia", 'New_Guinea' = 'New Guinea'))

#add column for year of establishment for Pacific populations
presentation1$establishment.year <- with(presentation1, ifelse(Country.Archipelago == 'Australia', 1871, ifelse(Country.Archipelago == 'Fiji', 1877, ifelse(Country.Archipelago == 'Hawaii', 1841, ifelse(Country.Archipelago == 'Marianas', 1887, ifelse(Country.Archipelago == 'New Caledonia', 1868, ifelse(Country.Archipelago == 'New Guinea', 1875, ifelse(Country.Archipelago == 'Samoa', 1863, NA))))))))

(fig4a <- ggplot(presentation1[!is.na(presentation1$Country.Archipelago),], aes(x=year,y=Size.PC1))+
  geom_point(col = 'blue')+
  theme_classic()+
  facet_grid(~Country.Archipelago)+
  stat_smooth(method = 'lm', fill = 'blue',col="black",alpha=0.25,fullrange=F)+
  xlab("Year")+
  ylab("Size PC1")+
  theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size=16))+
  geom_hline(yintercept = mean(North_America$Size.PC1, na.rm=T), lty = 2)+
  geom_hline(yintercept = mean(Caribbean$Size.PC1, na.rm=T), lty =2)+
  scale_x_continuous(limits = c(1840,2025))+
  geom_vline(aes(xintercept = establishment.year), col = 'red', lty = 2, size = 0.5)+
  ggtitle('Pacific Islands')+
  theme(plot.title = element_text(hjust = 0.5, size = 20)))

### now for North America

ggplot(North_America, aes(x = year, y = Size.PC1))+
  geom_point(col = 'limegreen')+
  theme_classic()+
  geom_smooth(method='lm', col = 'black')+
  ylab('Size PC1')+
  xlab("Year")+
  theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size=16))+
  geom_hline(yintercept = mean(North_America$Size.PC1, na.rm=T), lty = 2)+
  geom_hline(yintercept = mean(Caribbean$Size.PC1, na.rm=T), lty =2)+
  scale_x_continuous(limits = c(1850,2025))+
  ylim(c(-4,4))+
  ggtitle('North America')+
  theme(plot.title = element_text(hjust = 0.5, size = 20))

## and the Atlantic

ggplot(wings[wings$Region=='Atlantic',], aes(x = year, y = Size.PC1))+
  geom_point(col = 'purple')+
  theme_classic()+
  geom_smooth(method='lm', col = 'black')+
  ylab('Size PC1')+
  xlab("Year")+
  theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size=16))+
  geom_hline(yintercept = mean(North_America$Size.PC1, na.rm=T), lty = 2)+
  geom_hline(yintercept = mean(Caribbean$Size.PC1, na.rm=T), lty =2)+
  scale_x_continuous(limits = c(1850,2025))+
  ylim(c(-4,4))+
  ggtitle('Atlantic')+
  theme(plot.title = element_text(hjust = 0.5, size = 20))

##add central america for comparison

ggplot(wings[wings$Region=='Central America',], aes(x = year, y = Size.PC1))+
  geom_point(col = 'gold')+
  theme_classic()+
  geom_smooth(method='lm', col = 'black')+
  ylab('Size PC1')+
  xlab("Year")+
  theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size=16))+
  geom_hline(yintercept = mean(North_America$Size.PC1, na.rm=T), lty = 2)+
  geom_hline(yintercept = mean(Caribbean$Size.PC1, na.rm=T), lty =2)+
  scale_x_continuous(limits = c(1850,2025))+
  ylim(c(-4,4))+
  ggtitle('Central America')+
  theme(plot.title = element_text(hjust = 0.5, size = 20))

################
################################################################################################################################################################################################################################################################################################################################

#### common garden data

adults <- read.csv(file = './adult_morphology.csv')

library(multcomp)

#make size PC, as done above

adults <- adults[!is.na(adults$RArea),]

adults$MArea <- rowMeans(cbind(adults$LArea, adults$RArea), na.rm = T)
adults$MLength <- rowMeans(cbind(adults$LLength,adults$RLength), na.rm = T)
adults$MWidth <- rowMeans(cbind(adults$LWidth,adults$RWidth), na.rm=T)
adults$MPer <- rowMeans(cbind(adults$LPer,adults$RPer), na.rm=T)
adults$roundness <- (4*pi*adults$MArea)/(adults$MPer)^2
adults$aspect_ratio <- adults$MLength / adults$MWidth

ad.col.keep.size <- c('MArea','MLength','MWidth')

adults.size.keep <- adults[,ad.col.keep.size]

fit.size = princomp(adults.size.keep, cor = T)
summary(fit.size) #not surprisingly first PC accounts for 96.4% of variance
adults$Size.PC1 = fit.size$scores[,1]

adults$Mon.Pop <- factor(adults$Mon.Pop, c('ENA','CA','HI','GU','AU','PR'))
adults$Mon.Pop <- revalue(adults$Mon.Pop, c('ENA' = 'Eastern North America','CA' = 'Western North America','HI' = 'Hawaii','GU' = 'Guam', 'AU' = 'Australia', 'PR' = 'Puerto Rico'))

ggplot(adults, aes(x = Mon.Pop, y = Size.PC1))+
  geom_point(aes(col = Mon.Pop), position = position_jitterdodge(3.2), size = 1.7)+
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2, size = 5)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank() )+
  scale_color_viridis_d()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Size PC1")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16), legend.position = 'bottom')

#alternate color for same figure

adults$pa <- revalue(adults$Mon.Pop, c('Eastern North America' = 'ENA', 'Western North America' = 'WNA', 'Hawaii' = 'HI', 'Guam' = 'GU', 'Australia' = 'AU', 'Puerto Rico' = 'PR'))

ggplot(adults, aes(x = pa, y = Size.PC1))+
  geom_point(aes(col = pa), position = position_jitterdodge(3.2), size = 1.7)+
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2, size = 5)+
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 60, size =16, hjust = 1))+
  scale_color_manual(values = c('limegreen','limegreen','blue','blue','blue','orange'))+
  ylab("Size PC1")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.position = 'none')+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})

ad.col.keep.shape <- c('aspect_ratio','roundness') 

adults.shape.keep <- adults[,ad.col.keep.shape]
fit.shape = princomp(adults.shape.keep, cor = T)
summary(fit.shape) #first pc explains 78% of shape variance (lower than in wild caught butterflies)

adults$Shape.PC1 <- fit.shape$scores[,1]

(figs2a <- ggplot(adults, aes(x = Mon.Pop, y = Shape.PC1))+
  geom_point(aes(col = Mon.Pop), position = position_jitterdodge(3.2), size = 1.7)+
  geom_boxplot(outlier.size = 0, outlier.color = 'white', alpha = 0.5, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2, size = 5)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank() )+
  scale_color_viridis_d()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Shape PC1")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16), legend.position = 'bottom'))

#pdf(file = './Second submission/figures/FigS2a.pdf', width = 8, height = 8)
#figs2a
#dev.off()

dates <- as.Date(adults$Date, format='%m/%d/%Y')

start_date = as.Date('2018-05-15', tz="UTC") #define earliest date on which caterpillars were added to plants (May 15)

adults$exp.days <- as.numeric(difftime(dates, start_date , units = "days")) + 730486


#test for size differences among populations
loc.ad.adult.size.model1 <- lmer(Size.PC1 ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group) + GH + scale(OE) + Sex + Year, data = adults)

summary(loc.ad.adult.size.model1)
Anova(loc.ad.adult.size.model1) #monarch population explains far more variation than host plant ID
data.frame(emmeans(loc.ad.adult.size.model1, 'Mon.Pop'))


#now estimate broad sense heritability of wing size
wing.size.heritability <- lmer(Size.PC1 ~ (1|maternal_family), adults)
as.data.frame(VarCorr(wing.size.heritability))$vcov[1] / sum(as.data.frame(VarCorr(wing.size.heritability))$vcov) #broad sense heritability for wing size is 0.289


summary(glht(loc.ad.adult.size.model1, linfct = mcp(Mon.Pop = "Tukey"))) #besides Guam, all populations are smaller than the migratory North American populations
summary(glht(loc.ad.adult.size.model1, linfct = mcp(Species = "Tukey"))) #no significant pairwise differences among milkweed hosts

#now for forewing shape
loc.ad.adult.shape.model1 <- lmer(Shape.PC1 ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family) + GH + Year + scale(OE) + Sex + Year, data = adults)

Anova(loc.ad.adult.shape.model1) #again, monarch population has a strong impact on wing shape; surprisingly, milkweed species explains a greater proportion of variation in wing shape than wing size; no difference between sexes in wing shape
data.frame(emmeans(loc.ad.adult.shape.model1, 'Mon.Pop'))

summary(glht(loc.ad.adult.shape.model1, linfct = mcp(Mon.Pop = "Tukey"))) #Puerto Rico significantly less elongated than all other populations; no other pairwise differences
summary(glht(loc.ad.adult.shape.model1, linfct = mcp(Species = "Tukey"))) #only one pairwise diff among milkweed species (GOPH more elongated than AINC)

#now estimate broad sense heritability of wing shape
wing.shape.heritability <- lmer(Shape.PC1 ~ (1|maternal_family), adults)
as.data.frame(VarCorr(wing.shape.heritability))$vcov[1] / sum(as.data.frame(VarCorr(wing.shape.heritability))$vcov) #broad sense heritability for wing size is 0.367

#eclosion mass

loc.ad.adult.mass.model1 <- lmer(emergence_weight ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family) + GH + Year + scale(OE) + Sex + Year, data = adults)

Anova(loc.ad.adult.mass.model1)

summary(glht(loc.ad.adult.mass.model1, linfct = mcp(Mon.Pop = "Tukey")))
summary(glht(loc.ad.adult.mass.model1, linfct = mcp(Species = "Tukey")))

(figs2b <- ggplot(adults, aes(x = Mon.Pop, y = emergence_weight))+
  geom_point(aes(col = Mon.Pop), position = position_jitterdodge(3.2), size = 1.7)+
  geom_boxplot(outlier.size = 0, outlier.color = 'white', alpha = 0.5, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2, size = 5)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank() )+
  scale_color_viridis_d()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Eclosion Mass (g)")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16), legend.position = 'bottom'))

#pdf(file = './Second submission/figures/FigS2b.pdf', width = 8, height = 8)
#figs2b
#dev.off()

#now estimate broad sense heritability of eclosion mass
eclosion.mass.heritability <- lmer(emergence_weight ~ (1|maternal_family), adults)
as.data.frame(VarCorr(eclosion.mass.heritability))$vcov[1] / sum(as.data.frame(VarCorr(eclosion.mass.heritability))$vcov) #broad sense heritability for eclosion mass is 0.245

#now do wing loading, using wet mass
adults$wing.loading <- adults$emergence_weight / adults$MArea

loc.ad.adult.wl.model1 <- lmer(wing.loading ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family) + GH + Year + scale(OE) + Sex + Year, data = adults)

Anova(loc.ad.adult.wl.model1) #again, monarch population explains lot of variation; monarch sex is highly significant as well

summary(glht(loc.ad.adult.wl.model1, linfct = mcp(Mon.Pop = "Tukey")))
summary(glht(loc.ad.adult.wl.model1, linfct = mcp(Species = "Tukey")))


(figs2c <- ggplot(adults, aes(x = Mon.Pop, y = wing.loading))+
  geom_point(aes(col = Mon.Pop), position = position_jitterdodge(3.2), size = 1.7)+
  geom_boxplot(outlier.size = 0, outlier.color = 'white', alpha = 0.5, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2, size = 5)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank() )+
  scale_color_viridis_d()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Wing Loading (Total Mass (mg) / Forewing Area (cm^2))")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16), legend.position = 'bottom'))

#pdf(file = './Second submission/figures/FigS2c.pdf', width = 8, height = 8)
#figs2c
#dev.off()

#now estimate broad sense heritability of wing loading
wl.heritability <- lmer(wing.loading ~ (1|maternal_family), adults)
as.data.frame(VarCorr(wl.heritability))$vcov[1] / sum(as.data.frame(VarCorr(wl.heritability))$vcov) #broad sense heritability for wing loading is 0.231

#now load thorax and abdomen mass data

thor.ab <- read.csv(file = './thorax_abdomen_mass.csv')

#merge with the rest of the data

tsx <- merge(thor.ab, adults, by = 'Cat.ID')

tsx$drymass <- tsx$thorax + tsx$abdomen

#first check correlation between wet and dry mass

plot(drymass ~ emergence_weight, tsx) #clearly very tightly correlated, so probably not worth adding another plot explicitly for dry mass

tsx$ab.thor.ratio <- tsx$thorax / tsx$abdomen

loc.ad.thor.ab <- lmer(ab.thor.ratio ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family) + GH + scale(OE) + Sex.y, data = tsx[tsx$fresh == 'yes',])

Anova(loc.ad.thor.ab) #species differ, more modest differences among sex and monarch populations; strong effect of OE, whereby it decreases thorax mass more strongly than abdominal mass (makes sense in the context of transmission dynamics -- should favor reproduction and vertical transmission but not dispersal)

summary(glht(loc.ad.thor.ab, linfct = mcp(Mon.Pop = "Tukey"))) #no significant pairwise differences among populations
summary(glht(loc.ad.thor.ab, linfct = mcp(Species = "Tukey")))

plot.tsx <- tsx[tsx$fresh=='yes',] #restrict to only individuals that were collected on the same day that they eclosed for plotting (same data was used for analysis above)

(figs2d=ggplot(plot.tsx[!is.na(plot.tsx$Mon.Pop),], aes(x = Mon.Pop, y = ab.thor.ratio))+
  geom_point(aes(col = Mon.Pop), position = position_jitterdodge(3.2), size = 1.7)+
  geom_boxplot(outlier.size = 0, outlier.color = 'white', alpha = 0.5, size = 0.7)+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 2, size = 5)+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank())+
  theme(legend.title=element_blank() )+
  scale_color_viridis_d()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylab("Thorax Mass / Abdomen Mass")+
  theme(axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=16), legend.position = 'bottom'))

#pdf(file = './Second submission/figures/FigS2d.pdf', width = 8, height = 8)
#figs2d
#dev.off()

#and finally estimate broad sense heritability of ab/thor ratio
abthor.heritability <- lmer(ab.thor.ratio ~ (1|maternal_family), tsx)
as.data.frame(VarCorr(abthor.heritability))$vcov[1] / sum(as.data.frame(VarCorr(abthor.heritability))$vcov) #broad sense heritability for abdomen to thorax ratio is 0.094

#####

#now do MANOVA that incorporates multiple morphological traits (for now, stick to wing size, wing shape, body mass)

adults_manova_model1 <- manova(cbind(Size.PC1, Shape.PC1, emergence_weight) ~ Species + Mon.Pop + sym.allo + (OE) + Sex + Error(maternal_family), data = adults)

summary(adults_manova_model1)

adults_manova_model2 <- manova(cbind(Size.PC1, Shape.PC1, emergence_weight, ab.thor.ratio) ~ Species + Mon.Pop + sym.allo + (OE) + Sex.x + Error(maternal_family), data = tsx)

summary(adults_manova_model2)

#### create dedicated column for wing loading

adults$wing.loading <- adults$emergence_weight / adults$MArea

############################################
#correlations among traits
############################################
names(adults)
names(thor.ab)

corr.all.adult <- merge(adults, thor.ab[thor.ab$fresh=='yes',], by = 'Cat.ID', all = T)

names(corr.all.adult)

corr.all.adult$thor.ab.ratio <- corr.all.adult$thorax/corr.all.adult$abdomen
corr.all.adult$dry.mass <- corr.all.adult$thorax + corr.all.adult$abdomen

names.keep <- c('Cat.ID','Year','Species','Mon.Pop','Size.PC1','Shape.PC1','wing.loading','Sex.y','thorax','abdomen','emergence_weight')
names.keep2 <- c('Size.PC1','Shape.PC1','wing.loading','thorax','abdomen','emergence_weight','dry.mass','thor.ab.ratio')
full.names.cor <- c('MLength','MWidth','MArea','roundness','aspect_ratio','Size.PC1','Shape.PC1','emergence_weight','dry.mass','wing.loading','thorax','abdomen','thor.ab.ratio')

emb <- na.omit(corr.all.adult[,names.keep2])
names(emb) <- c('Size PC1', 'Shape PC1', 'Wing Loading','Thorax Mass', 'Abdomen Mass','Total Mass', 'Total Mass (Dry)', 'Thorax to Abdomen Ratio')
emb <- emb[,c(1,2,6,7,3,4,5,8)]

M <- cor(emb)

library(corrplot)

corrplot(M, method = 'circle')

library(PerformanceAnalytics)

chart.Correlation(emb)

dmc <- na.omit(corr.all.adult[,full.names.cor])

names(dmc) <- c('Forewing Length', 'Forewing Width', 'Forewing Area','Forewing Roundness','Aspect Ratio','Size PC1','Shape PC1','Total Mass','Total Mass (Dry)','Wing Loading','Thorax Mass','Abdomen Mass','Thorax to Abdomen Ratio')

chart.Correlation(dmc)

N = cor(dmc)

figs3a <- corrplot(N, method = 'circle')

#pdf(file = './Second submission/figures/figs3a.pdf', height = 5, width = 5)
#figs3a
#dev.off()

###############

#figure with early records; here, use specimens from the first 50 years post-establishment

earlies <- wings[c(which(wings$Country.Archipelago=='Marianas' & wings$year < 1938),which(wings$Country.Archipelago == 'Australia' & wings$Region=='Pacific Islands' & wings$year < 1922),which(wings$Country.Archipelago == 'Samoa' & wings$year < 1919),which(wings$Country.Archipelago == 'New_Guinea' & wings$year < 1926),which(wings$Country.Archipelago == 'Hawaii' & wings$year < 1901),which(wings$Country.Archipelago == 'Fiji' & wings$year < 1917),which(wings$Country.Archipelago =='New_Caledonia' & wings$year < 1919),which(wings$Country.Archipelago =='Morocco' & wings$year < 1937),which(wings$Country.Archipelago == 'French_Polynesia' & wings$year < 1933),which(wings$Country.Archipelago == 'Indonesia' & wings$year < 1921),which(wings$Country.Archipelago == 'Tonga' & wings$year < 1914),which(wings$Country.Archipelago == 'Vanuatu' & wings$year < 1919),which(wings$Country.Archipelago == 'Solomon_Islands' & wings$year < 1937),which(wings$Country.Archipelago == 'Taiwan' & wings$year < 1940),which(wings$Country.Archipelago == 'Marquesas' & wings$year < 1934),which(wings$Region == 'North America' & wings$year < 1921),which(wings$Country.Archipelago == 'British_Isles'),which(wings$Country.Archipelago =='Madeira' & wings$year < 1900), which(wings$Country.Archipelago=='Azores' & wings$year < 1910)),] #notably, nothing to include for New Zealand, Spain, Morocco, Portugal, Cape Verde, Micronesia

nrow(earlies)

#add UK as its own separate region

levels(earlies$Region) <- c(levels(earlies$Region),'British Isles')

earlies[which(earlies$Country.Archipelago == 'British_Isles'),]$Region <- 'British Isles'

early.means.size <- aggregate(Size.PC1 ~ Region, earlies, mean)
names(early.means.size)[2] <- 'mean.size'
early.means.shape <- aggregate(Shape.PC1 ~ Region, earlies, mean)
names(early.means.shape)[2] <- 'mean.shape'
early.se.size <- aggregate(Size.PC1 ~ Region, earlies, function(x) sd(x)/sqrt(length(x)))
names(early.se.size)[2] <- 'size.se'
early.se.shape <- aggregate(Shape.PC1 ~ Region, earlies, function(x) sd(x)/sqrt(length(x)))
names(early.se.shape)[2] <- 'shape.se'

early.size <- merge(early.means.size, early.se.size)
early.shape <- merge(early.means.shape, early.se.shape)

plot.early <- merge(early.size, early.shape, by  = 'Region')
plot.early$Region <- factor(plot.early$Region, levels = c('North America','Atlantic','Pacific Islands','British Isles'))

(fig2a <- ggplot(plot.early[which(plot.early$Region %in% c('North America','Pacific Islands','Atlantic','British Isles')),], aes(x = mean.size, y = mean.shape, group = Region, col = Region))+
    theme_bw()+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = mean.shape - shape.se, ymax = mean.shape + shape.se), size = 0.5, width = 0.02)+
    geom_errorbarh(aes(xmin = mean.size - size.se, xmax = mean.size + size.se), size = 0.5, height = 0.02)+
    scale_color_manual(name="Region",values=
                         c("limegreen","purple",'blue','black'))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    ylab("Shape PC1")+
    xlab('Size PC1')+
    theme(axis.title = element_text(size=16))+
    theme(legend.text = element_text(size=16))+
    theme(axis.text = element_text(size = 16))+
    theme(legend.title = element_blank(), legend.position = c(0.79,0.15))+
    theme(legend.background = element_rect(colour = 'black', fill = 'white'))+
    geom_hline(yintercept = 0, lty = 2)+
    geom_vline(xintercept = 0, lty = 2))

#pdf('./Second submission/figures/Figure_range.expansion.pdf', width = 5.5, height = 5)

#fig2a

#dev.off()

### run a quick model to test for differences in wing size and shape

#set North America as the reference level

earlies$Region <- relevel(earlies$Region, ref = 'North America')
earlies$expansion <- ifelse(earlies$Region =='North America','ancestral','expansion')

summary(lm(Size.PC1 ~ Region + Sex + year , earlies))
summary(lm(Shape.PC1 ~ Region + Sex + year, earlies)) #this works for quick analysis, but better approach might be to use parameters estimated from global models to account for various things like collection ID

######### finally test for differences in coefficient of variation between populations

library(cvequality)

with(na.omit(wings[!wings$Region %in% c('Asia','Indian_Ocean'),c('MArea','Region')]), asymptotic_test(x = MArea, y = droplevels(Region))) #confirms that CVs do differ between locations; can try pairwise comparisons

library(cvcqv)

#get CV for size and shape for various regions; start with North America

cv_versatile(
  adults[adults$Mon.Pop=='Eastern North America',]$MArea, 
  na.rm = TRUE, 
  digits = 3, 
  method = "basic", 
  correction = TRUE, 
  alpha = 0.05
)

cv_versatile(
  wings[wings$Region=='North America',]$MArea, 
  na.rm = TRUE, 
  digits = 3, 
  method = "basic", 
  correction = TRUE, 
  alpha = 0.05
)


wild_coef.vars.area <- list()
for(i in 1:length(levels(wings.plot$Region))){
  wild_coef.vars.area[i] <- cv_versatile(
    wings.plot[which(wings.plot$Region == 
                       levels(wings.plot$Region)[i]),]$MArea, 
    na.rm = TRUE, 
    digits = 3, 
    method = "basic", 
    correction = TRUE, 
    alpha = 0.05
  )[2]
}
names(wild_coef.vars.area) <- levels(wings.plot$Region)
wild_coef.vars.area<- ldply(wild_coef.vars.area)

wild_coef.vars.area$.id <- factor(wild_coef.vars.area$.id, levels = c('North America','Central America','South America','Caribbean','Pacific Islands','Atlantic'))
ggplot(wild_coef.vars.area, aes(x = .id, y = est, col = .id))+
  theme_bw()+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2, size = 1)+
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
  scale_color_manual(values = c('limegreen','gold','brown1','orange','blue','purple'))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), axis.title.x = element_blank())+
  ylab('Coefficient of Variation - Forewing Area')


wild_coef.vars.shape <- list()
for(i in 1:length(levels(wings.plot$Region))){
  wild_coef.vars.shape[i] <- cv_versatile(
    wings.plot[which(wings.plot$Region == 
                       levels(wings.plot$Region)[i]),]$aspectratio, 
    na.rm = TRUE, 
    digits = 3, 
    method = "basic", 
    correction = TRUE, 
    alpha = 0.05
  )[2]
}
names(wild_coef.vars.shape) <- levels(wings.plot$Region)
(wild_coef.vars.shape <- ldply(wild_coef.vars.shape)) #perhaps confusing: pacific and atlantic populations show increased CV for wing shape relative to all other populations in wild caught data
wild_coef.vars.shape$.id <- factor(wild_coef.vars.shape$.id, levels = c('North America','Central America','South America','Caribbean','Pacific Islands','Atlantic'))
ggplot(wild_coef.vars.shape, aes(x = .id, y = est, col = .id))+
  theme_bw()+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2, size = 1)+
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
  scale_color_manual(values = c('limegreen','gold','brown1','orange','blue','purple'))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), axis.title.x = element_blank())+
  ylab('Coefficient of Variation - Aspect Ratio') #note: this figure includes all data; figure shown in supplement excluded extreme shape outliers that were more than 4 SD away from the mean for all locations


com.gar_coef.vars.area.mp <- list()
for(i in 1:length(levels(adults$Mon.Pop))){
  com.gar_coef.vars.area.mp[i] <- cv_versatile(
    adults[which(adults$Mon.Pop == 
                       levels(adults$Mon.Pop)[i]),]$MArea, 
    na.rm = TRUE, 
    digits = 3, 
    method = "basic", 
    correction = TRUE, 
    alpha = 0.05
  )[2]
}
names(com.gar_coef.vars.area.mp) <- levels(adults$Mon.Pop)
(com.gar_coef.vars.area.mp <- ldply(com.gar_coef.vars.area.mp)) #not strong evidence for differences in CV among each individual population; try pooling by region instead
ggplot(com.gar_coef.vars.area.mp, aes(x = .id, y = est, col = .id))+
  geom_point()+
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2)


adults$region <- as.factor(ifelse(adults$Mon.Pop %in% c('Western North America','Eastern North America'), 'North America', ifelse(adults$Mon.Pop == 'Puerto Rico', 'Caribbean', 'Pacific Islands')))
com.gar_coef.vars.area.reg <- list()
for(i in 1:length(levels(adults$region))){
  com.gar_coef.vars.area.reg[i] <- cv_versatile(
    adults[which(adults$region == 
                   levels(adults$region)[i]),]$MArea, 
    na.rm = TRUE, 
    digits = 3, 
    method = "basic", 
    correction = TRUE, 
    alpha = 0.05
  )[2]
}
names(com.gar_coef.vars.area.reg) <- levels(adults$region)
(com.gar_coef.vars.area.reg<- ldply(com.gar_coef.vars.area.reg)) #corroborates previous analysis, though all are overlapping and not sifnificantly different
#CV of wing area for common garden reared monarchs
com.gar_coef.vars.area.reg$.id <- factor(com.gar_coef.vars.area.reg$.id, c('North America','Caribbean','Pacific Islands'))
ggplot(com.gar_coef.vars.area.reg, aes(x = .id, y = est, col = .id))+
  theme_bw()+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2, size = 1)+
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
  scale_color_manual(values = c('limegreen','orange','blue'))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), axis.title.x = element_blank())+
  ylab('Coefficient of Variation - Forewing Area')


com.gar_coef.vars.shape.reg <- list()
for(i in 1:length(levels(adults$region))){
  com.gar_coef.vars.shape.reg[i] <- cv_versatile(
    adults[which(adults$region == 
                   levels(adults$region)[i]),]$aspect_ratio, 
    na.rm = TRUE, 
    digits = 3, 
    method = "basic", 
    correction = TRUE, 
    alpha = 0.05
  )[2]
}
names(com.gar_coef.vars.shape.reg) <- levels(adults$region)
(com.gar_coef.vars.shape.reg<- ldply(com.gar_coef.vars.shape.reg)) #again all are overlapping
com.gar_coef.vars.shape.reg$.id <- factor(com.gar_coef.vars.shape.reg$.id, c('North America','Caribbean','Pacific Islands'))
ggplot(com.gar_coef.vars.shape.reg, aes(x = .id, y = est, col = .id))+
  theme_bw()+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.2, size = 1)+
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
  scale_color_manual(values = c('limegreen','orange','blue'))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1), axis.title.x = element_blank())+
  ylab('Coefficient of Variation - Aspect Ratio')
  


  


cv_versatile(
  wings[wings$Region=='Caribbean',]$MArea, 
  na.rm = TRUE, 
  digits = 3, 
  method = "basic", 
  correction = TRUE, 
  alpha = 0.05
)[2]








##### north american monarchs: latitude effect when only looking between march and june?

sec.analy <- lmer(Size.PC1 ~ lat + Sex + year + (1|Collection), North_America)
summary(sec.analy) #overall get an effect of latitude; now restrict to only eastern NA

sec.analy2 <- lmer(Size.PC1 ~ scale(lat) + Sex + year + (1|Collection), North_America[North_America$lon > -110,])
summary(sec.analy2) #still get an effect when looking at only eastern monarchs; now restirct to only eastern monarchs collected between March - July

sec.analy3 <- lmer(Size.PC1 ~ lat + Sex + year + (1|Collection), North_America[North_America$month %in% c('Mar','Apr','May','Jun','Jul') & North_America$lon > -110,])
summary(sec.analy3) #still get an effect of latitude here; pattern seems to be real, although excluding July (and June + July) erases the effect; the analysis above includes 278 records, while excluding July leaves only 128


