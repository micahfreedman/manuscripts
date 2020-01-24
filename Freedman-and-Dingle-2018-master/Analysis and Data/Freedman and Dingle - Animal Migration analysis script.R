library(ggplot2)
library(lme4)
library(plyr)
library(rworldmap)
library(MuMIn)
library(car)
library(lsmeans)
library(multcomp)
library(sjPlot)
library(AICcmodavg)
library(stringr)
library(PerformanceAnalytics)
library(lmerTest)
library(svglite)
library(gridExtra)
library(cowplot)

options(contrasts=c("contr.sum", "contr.poly"))

wings <- read.csv(file="~/Downloads/North_American_wings.csv",header=TRUE)

head(wings)
names(wings)

####Split dates into day, month, and year, then rename columns accordingly and convert to numeric values

date.split <- t(sapply(wings$Collection_Date, function(x) substring(x, first=c(1,5,7), last=c(4,6,8))))

wings <- cbind(wings, date.split)

names(wings)[names(wings) == '1'] <- 'year'
names(wings)[names(wings) == '2'] <- 'month'
names(wings)[names(wings) == '3'] <- 'day'

wings$year <- as.numeric(levels(wings$year))[wings$year]
wings$day <- as.numeric(levels(wings$day))[wings$day]

####Convert days with entry of 0 to 15 (so as to not discard data...)

wings$day[wings$day == 0] <- 15


####Replace months with factor levels and arrange chronologically

months <- c("January","February","March","April","May","June","July","August",
            "September","October","November","December")
wings$month <- months[ wings$month ]
wings$month <- as.factor(wings$month)
wings$month = factor(wings$month, levels = c("January","February","March","April","May","June","July","August","September","October","November","December"))


####Also add column for Julian day of collection

date<-as.Date(as.character(wings$Collection_Date), "%Y%m%d")

wings$julian.date <- as.numeric(format(date, "%j"))

####Now add another column to capture that Julian dates in this dataset cannot be treated as a linear factor.  Instead, treat July 1st as the most "summer-like" date in the dataset, which will have the minimum value for the cyclic Julian day measure.  Every other day will be expressed as abs(x) away from this date. July 1st correponds to Julian Day 183 in a normal year, which is almost exactly halfway through the year.


wings$julian.date.continuous <- abs(wings$julian.date - 183)
hist(wings$julian.date,breaks=50)

####Histogram of frequency of observations across the full year. Peaks in abundance correspond to single collection efforts, such as those from the Yang et al. 2016 paper. Can also peak in abundance of observations around late summer.

hist(wings$julian.date.continuous,breaks = 50)

####"Folded"" version of the previous histogram, such that the x-axis corresponds to number of days away from July 1st

####Add one more column to account for cyclic nature of seasonality in Julian date measure: scale according to cosine wave, with maximum value at June 21 (corresponding to longest days of the year), and minimum value at December 22 (shortest days of the year). Here, the expectation is that butterflies collecting during summer periods might be expected to be smaller, all other things considered.

wings$daylength.cycle <- -cos(pi*(wings$julian.date + 10)/182)

(FigS1b = ggplot(wings, aes( x = julian.date, daylength.cycle))+
  geom_jitter(pch = 21, cex = 1)+
  theme_bw()+
  ylab("Daylength Index")+
  xlab("Julian Date")+
  theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 14,angle=75, vjust=0.5), axis.text.y = element_text(size =14))+
  scale_x_continuous(minor_breaks = 20, breaks = seq(1,365,15))+
  geom_vline(xintercept=172, col = 'red', size = 1, lty = 2)+
  geom_vline(xintercept = c(80,264), col = 'orange', size = 1, lty = 2)+
  geom_vline(xintercept=355, col = 'purple', size = 1, lty = 2)+
  annotate('text', x = 10, y = 0.9, label = '(b)', size = 8))

(FigS1a = ggplot(wings, aes( x = julian.date))+
  geom_histogram(binwidth = 9)+
  theme_bw()+
  ylab("Count of Monarch Observations")+
  xlab("Julian Date")+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14))+
  annotate('text', x = 2, y = 135, label = '(a)', size = 8))


FigS1 <- grid.arrange(FigS1a,FigS1b, nrow = 2)


ggsave(file="~/Documents/grad_school/Manuscripts/Animal Migration/Figures/FigS1.pdf", plot=FigS1, width=8, height=8)


####This plot shows all of the observations throughout the calendar year (x-axis), and the shape captures the nature of how one might expect for seasonality to progress (i.e. steepest changes in seasonality occur at equixones, corresponding roughly to days 90 and 270)

###Add columns for average wing values

wings$MLength <- rowMeans(cbind(wings$LLength, wings$RLength),na.rm=T)
wings$MWidth <- rowMeans(cbind(wings$LWidth, wings$RWidth),na.rm=T)
wings$MArea <- rowMeans(cbind(wings$LArea, wings$RArea),na.rm=T)*100
wings$MPerimeter <- rowMeans(cbind(wings$LPerimeter,wings$RPerimeter),na.rm=T)
wings$aspectratio <- wings$MLength/wings$MWidth
wings$roundness <- (4*pi*wings$MArea)/(wings$MPerimeter)^2


####Drop entries for butterflies that were not wild-caught

wings <- wings[wings$Wild.Caught. == "yes",]


####Drop the sole gyandromorph from the data

wings <- wings[wings$Sex != "gyandromorph",]

####Drop entries coded as North America but more than 1 degree of latitude south from overwintering locations


wings <- wings[wings$lat > 19.5,]

nrow(wings)



####Dataframe with 1825 entries

####Add a column for eastern vs. western North America; set everything west of the 110th parallel as western, and everything east of this line eastern

wings$East.versus.West <- ifelse(wings$lon < -110 & wings$lat, "West", "East")

####Generate map of North America, with all points labelled. Florida non-migrants are labelled in aquamarine and will be excluded from the main analysis of North American migrants.

#Figure 2

cols = c('blue','forestgreen')
north.am.map <- getMap(resolution = 'low')
plot(north.am.map, xlim=range(-130,-68), ylim = range(18,50))
points(wings$lon, wings$lat, col = cols[as.factor(wings$East.versus.West)], cex = 1, pch = 20)
points(wings[wings$Island.State=="Florida" & wings$lat < 28,]$lon, wings[wings$Island.State=="Florida" & wings$lat < 28,]$lat, col = 'aquamarine', cex = 1, pch = 20)
abline(v=-110, lty = 2)
rect(xleft = -84, ybottom = 23.5, xright = -79, ytop = 28, lty = 2, border = 'red')


####Exclude south Florida non-migrants from dataset


wings.drop <- which(wings$Island.State == "Florida" & wings$lat < 28)

wings <- wings[-wings.drop,]

nrow(wings)


#### Restricts dataset to 1773 observations

#### Build a linear model to investigate the effects of the following variables on wing area:

####First, dabble with a few test models to determine which random effects should be included

wing.area.model.test1 <- lmer( MArea ~ scale(year) + Sex + scale(daylength.cycle) + East.versus.West + scale(lat) + (scale(lat)|East.versus.West) + (scale(lon)|East.versus.West) + scale(lon) + overwintering. + (1|Collection) + (1|observer) + (1|image_type) + (1|Island.State), data = wings)

####Compare to a reduced model without latitude and longitude nested in east vs. west

wing.area.model.test2 <- lmer( MArea ~ scale(year) + Sex + scale(daylength.cycle) + East.versus.West + scale(lat) + scale(lon) + overwintering. + (1|Collection) + (1|observer) + (1|image_type) + (1|Island.State), data = wings)

AIC(wing.area.model.test1,wing.area.model.test2)

####Simpler model 2 is preferred. Try another model, this time with no effect of state

wing.area.model.test3 <- lmer( MArea ~ scale(year) + Sex + East.versus.West + scale(lat) + scale(lon) + overwintering. + (1|Collection) + (1|observer) + (1|image_type) + scale(daylength.cycle) , data = wings)

AIC(wing.area.model.test1,wing.area.model.test2,wing.area.model.test3)

####Model 3 is slightly preferred. 

####Try another model without longitude, as this is partly covered by the east versus west term

wing.area.model.test4 <- lmer( MArea ~ scale(year) + Sex + East.versus.West + scale(lat) + overwintering. + (1|Collection) + (1|observer) + (1|image_type) + scale(daylength.cycle) + (1|Island.State), data = wings)

AIC(wing.area.model.test1,wing.area.model.test2,wing.area.model.test3,wing.area.model.test4)

####Best performing model is 4, although it's close.

####Try another model without state of collection as a random effect, and an even simpler model that does not include east vs. west or state at all. Also remove image type and observer to see if they matter.

wing.area.model.test5 <- lmer( MArea ~ scale(year) + Sex + East.versus.West + scale(lat) + overwintering. + (1|Collection) + (1|observer) + (1|image_type) + scale(daylength.cycle), data = wings)

wing.area.model.test6 <- lmer( MArea ~ scale(year) + Sex + scale(lat) + overwintering. + (1|Collection) + (1|observer) + (1|image_type) + scale(daylength.cycle), data = wings)

wing.area.model.test7 <- lmer( MArea ~ scale(year) + Sex + East.versus.West + scale(lat) + overwintering. + (1|Collection) + scale(daylength.cycle) + (1|Island.State), data = wings)


AIC(wing.area.model.test1,wing.area.model.test2,wing.area.model.test3,wing.area.model.test4,wing.area.model.test5,wing.area.model.test6,wing.area.model.test7) #model7 is pretty clearly the best fitting model here, which includes all of the suspected fixed effects

####Start with model that includes interactions between sex and year, east vs. west and year, and sex with east vs. west. Begin with full 3 way interaction between these factors, plus an interaction between overwintering status and east vs. west and overwintering status * sex

model1.wing.area <- lmer(MArea ~ scale(year)*Sex*East.versus.West + overwintering.*Sex + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

#reduce 3-way interaction in two way interactions

model2.wing.area <- lmer(MArea ~ scale(year)*Sex + Sex*East.versus.West + scale(year)*East.versus.West + overwintering.*Sex + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model1.wing.area,model2.wing.area) #model2 is preferred

#drop year*sex interaction

model3.wing.area <- lmer(MArea ~ Sex*East.versus.West + scale(year)*East.versus.West + overwintering.*Sex + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model1.wing.area,model2.wing.area,model3.wing.area) #3 is strongly preferred

#drop overwintering * sex interaction

model4.wing.area <- lmer(MArea ~ Sex*East.versus.West + scale(year)*East.versus.West + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model1.wing.area,model2.wing.area,model3.wing.area,model4.wing.area) #model 4 is best

#next drop interaction between east vs. west and year

model5.wing.area <- lmer(MArea ~ Sex*East.versus.West + scale(year) + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model4.wing.area, model5.wing.area) #model5 is best

#next drop interaction between east vs. west and sex

model6.wing.area <- lmer(MArea ~ Sex + scale(year) + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model5.wing.area,model6.wing.area) #model6 is best

#try another model without an interaction between east vs. west and overwintering

model7.wing.area <- lmer(MArea ~ Sex + scale(year) + overwintering. + East.versus.West + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model6.wing.area,model7.wing.area) #model6 is still better

#try another model without east versus west altogether

model8.wing.area <- lmer(MArea ~ Sex + scale(year) + overwintering. + scale(lat) + scale(daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

AIC(model6.wing.area,model7.wing.area,model8.wing.area) #this time, model 8 is preferred and doesn't include east versus west at all, suggesting that it probably isn't an important source of variation in the data

#try a model that leaves in east versus west and interaction term but drops Julian date term

model9.wing.area <- lmer(MArea ~ Sex + scale(year) + overwintering.*East.versus.West + scale(lat) + (1|Collection) + (1|Island.State), data = wings)

AIC(model6.wing.area, model7.wing.area, model8.wing.area, model9.wing.area) #do not drop Julian date

####Model 8 still seems to fit the data best. However, because we are specifically interested in the hypothesis that there is a difference between eastern versus western butterflies, use model 6, which included an interaction between eastern vs. western and overwintering status, such that overwintering eastern butterflies are larger than overwintering western butterflies.

####Graphical representation and summary of main effects

plot_model(model6.wing.area)

summary(model6.wing.area)

####Overwintering butterflies are signficantly larger than those collected in summer breeding areas. Butterflies from higher latitudes are larger, even when taking into account Mexican overwintering locations that are quite far south. Butterflies have become signficantly larger through time. Males are signficantly larger than females. East versus west is not included in this model, because

####Plot random effects from model

sjp.lmer(model6.wing.area)

####Only real thing to note here is that the random effects associated with state are fairly small and all seem to be "shrunken" toward the overall mean; biggest outlier seems to be California, where butterflies are small.  However, among collections, there is one signficant outlier (David M), which likely reflects some error in his photographing process. His images are wayyyy bigger than would be expected given all of the other data.

model6.not.scaled <- lmer(MArea ~ Sex + (year) + overwintering.*East.versus.West + (lat) + (daylength.cycle) + (1|Collection) + (1|Island.State), data = wings)

summary(model6.not.scaled)

##Table 1

####ANOVA table summarizing contributions of main effects in wing area analysis

Anova(model6.wing.area, type = 3)

###Generate figures to graphically depict some of the effects described in the model

##Figure 5a - Change through time in wing area, using raw data

(Fig5a <- ggplot(wings, aes( x = year, y = MArea))+
  geom_jitter(col = 'blue')+
  geom_smooth(method = 'lm', col = 'black', fill = 'grey')+
  theme_bw()+
  ylim(250,1250)+
  labs(x = "Year of Collection", y = bquote('Wing Area ('*mm^2*')'))+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)))


####Alternate Figure 1b, using a histogram instead of scatterplot

####first need to split up data in time chunks

wings$time.chunk <- ifelse(wings$year < 1900, "1875-1900",
                           ifelse(wings$year < 1925 & wings$year > 1900, "1900-1925",
                                  ifelse(wings$year < 1950 & wings$year > 1925, "1925-1950",
                                         ifelse(wings$year < 1975 & wings$year > 1950, "1950-1975",
                                                ifelse(wings$year < 2000 & wings$year > 1975, "1975-2000", "2000-Present")))))

ggplot(wings[!is.na(wings$time.chunk),], aes(x = time.chunk, y = MArea))+
  geom_boxplot()+
  theme_bw()+
  labs(x = "Year of Collection", y = bquote('Wing Area ('*cm^2*')'))+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14, angle = 60, hjust=1))


####Shows same general pattern: with the exception of 1875-1900, each successive time period has seen an increase in size

###Figure showing latitudinal pattern; use residuals after accounting for collection ID and excluding overwintering individuals here


lat.figure <- wings[wings$overwintering. == "no",]

(Fig4a <- ggplot(lat.figure, aes( x = lat, y = MArea))+
  geom_point(col = 'blue')+
  geom_smooth(method = 'lm', col = 'black', fill = 'grey')+
  theme_bw()+
  labs(x = "Latitude of Collection", y = bquote('Wing Area ('*mm^2*')'))+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)))
    
(Fig4b <- ggplot(lat.figure, aes( x = lat, y = roundness))+
  geom_point(col = 'blue')+
  geom_smooth(method = 'lm', col = 'black', fill = 'grey')+
  theme_bw()+
  labs(x = "Latitude of Collection", y = 'Wing Roundness')+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)))

plot_grid(Fig4a,Fig4b, ncol = 2, labels = c('(a)','(b)'), label_y = 0.95)


####Actually does not come out strongly when plotting raw data (in part due to effect of collection, where Mexican overwintering butterflies from David M collection distort southerly wing areas). However, after excluding overwintering individuals, pattern looks pretty clear.

####interaction between overwintering status and east vs. west





ggplot (wings[!is.na(wings$overwintering.),], aes( x = East.versus.West, y = MArea, fill = overwintering.))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c('orange','purple'), name = "Overwintering\nStatus")+
  labs(y = bquote('Wing Area ('*mm^2*')'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(size=18), axis.title = element_text(size = 22), legend.title = element_text(size = 18), legend.text = element_text(size = 14))+theme(legend.position = c(0.75,0.9))+theme(legend.title.align = 0.5, legend.box.background = element_rect(colour = "black"))+
  annotate('text', x = 0.81, y = 820, label = 'n=705', size = 4)+
  annotate('text', x = 1.19, y = 975, label = 'n=65', size = 4)+
  annotate('text', x = 1.81, y = 800, label = 'n=576', size = 4)+
  annotate('text', x = 2.19, y = 905, label = 'n=221', size = 4)

#65 overwintering in eastern, 705 non
#221 overwintering in western, 576 non


##Now conduct a separate analysis for wing roundness; go through same model selection procedure


model1.wing.roundness <- lmer(roundness ~ scale(year)*Sex*East.versus.West + overwintering.*Sex + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

#reduce 3-way interaction in two way interactions

model2.wing.roundness <- lmer(roundness ~ scale(year)*Sex + Sex*East.versus.West + scale(year)*East.versus.West + overwintering.*Sex + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model1.wing.roudness,model2.wing.roudness) #model2 is preferred

#drop year*sex interaction

model3.wing.roundness <- lmer(roundness ~ Sex*East.versus.West + scale(year)*East.versus.West + overwintering.*Sex + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model1.wing.roudness,model2.wing.roudness,model3.wing.roudness) #3 is strongly preferred

#drop overwintering * sex interaction

model4.wing.roundness <- lmer(roundness ~ Sex*East.versus.West + scale(year)*East.versus.West + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model1.wing.roundness,model2.wing.roundness,model3.wing.roundness,model4.wing.roundness) #model 4 is best

#next drop interaction between east vs. west and year

model5.wing.roundness <- lmer(roundness ~ Sex*East.versus.West + scale(year) + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model4.wing.roundness, model5.wing.roundness) #model5 is best

#next drop interaction between east vs. west and sex

model6.wing.roundness <- lmer(roundness ~ Sex + scale(year) + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model5.wing.roundness,model6.wing.roundness) #model6 is best

#try another model without an interaction between east vs. west and overwintering

model7.wing.roundness <- lmer(roundness ~ Sex + scale(year) + overwintering. + East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model6.wing.roundness,model7.wing.roundness) #model7 is better

#try another model without east versus west altogether

model8.wing.roundness <- lmer(roundness ~ Sex + scale(year) + overwintering. + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model6.wing.roundness,model7.wing.roundness,model8.wing.roundness) #this time, model 8 is preferred and doesn't include east versus west at all, suggesting that it probably isn't an important source of variation in the data

#try a  simple model, with no east vs. west or overwintering 

model9.wing.roundness <- lmer(roundness ~ Sex + scale(year) + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model8.wing.roundness,model9.wing.roundness)

#try very simple model with only latitude, sex, and random effects

model10.wing.roundness <- lmer(roundness ~ Sex + scale(lat) + (smoothing|Collection) + (1|Island.State), data = wings)

AIC(model9.wing.roundness,model10.wing.roundness)

#models get better by dropping nearly all predictors, but then can't estimate effect of these predictors at all

####For this analysis, use the same best-fitting model from above / don't bother doing whole model selection procedure; however, do include smoothing effect, since this is likely to affect roundness measurements



model6.wing.roundness <- lmer(roundness ~ Sex + scale(year) + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State), data = wings)

summary(model6.wing.roundness)

plot_model(model6.wing.roundness)

sjp.lmer(model6.wing.roundness)


####Based on results above, seems like only significant predictor of wing roudness is latitude, and in the expected direction (i.e. butterflies at higher latitudes are more elongate). In contrast to wing area, there is no significant effect through time, with sex, or at overwintering sites. There is a modest effect of the Julian date of collection, and in the direction that might be expected: butterflies collected closer to July 1st indeed have marginally rounder wings than those collected at other times of year.

##### ANCOVA analysis for statement about correlation between roundness and wing area: "Male forewings are also slightly more elongated, perhaps because of the inherent correlation between wing size and elongation"

model6b.wing.roundness <- lmer(roundness ~ Sex + scale(year) + overwintering.*East.versus.West + scale(lat) + scale(daylength.cycle) + (smoothing|Collection) + (1|Island.State) + MArea, data = wings)

summary(model6b.wing.roundness)

model6.wing.roundness.simple <- lmer(roundness ~ Sex + overwintering. + scale(lat) + (smoothing|Collection) + (1|Island.State) + MArea, data = wings)

summary(model6.wing.roundness.simple)

AIC(model6b.wing.roundness,model6.wing.roundness.simple)

###Table 1, Part II

####ANOVA table summarizing wing roudness results

Anova(model6.wing.roundness, type = 3)

## Now load the Flockhart data for comparison.

####Data here are the raw data from the Flockhart et al. (2017) paper in Movement Ecology, available as a supplement online

flockhart.data <- read.csv('~/Downloads/40462_2017_98_MOESM2_ESM.csv',header = T)

head(flockhart.data)


(Fig5b <- ggplot(flockhart.data, aes(x = WinterYear, y = WingArea_mm))+
  geom_point(col = 'goldenrod')+
  geom_smooth(method = 'lm', col = 'black', fill = 'grey')+
  theme_bw()+
  ylim(250,1250)+
  labs(x = "Collection Year", y = bquote('Wing Area ('*mm^2*')'))+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)))

(Fig5c <- ggplot(wings[wings$year>=1974,], aes(x = year, y = MArea))+
    geom_point(col = 'red')+
    geom_smooth(method = 'lm', col = 'black', fill = 'grey')+
    theme_bw()+
    ylim(250,1250)+
    labs(x = "Collection Year", y = bquote('Wing Area ('*mm^2*')'))+
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)))

plot_grid(Fig5a,Fig5b,Fig5c, nrow=1, labels = c('(a)','(b)','(c)'))



####Also indicates positive trend through time, although somewhat subtle. Check basic linear model to see if wing area has increased through time, after accounting for migration distance


flockhart.area.model1 <- lm(WingArea_mm ~ scale(WinterYear) + scale(dist.centroid.km) + Sex, flockhart.data)

summary(flockhart.area.model1)

####Significant effect of year through time, such that wings have increased by about 5.3mm^2 per year since 1974. This corresponds to an overall increase of about 0.2cm^2 over the 40 year time window.

####Compare this to the same set of dates from the whole dataset.  Select subset of dates from 1974-present and do not include overwintering individuals, and then run the same best-fitting model as before

wings.contemporary <- wings[wings$year >= 1974  & wings$overwintering. == "no",]

nrow(wings.contemporary)

wings.contemporary.area <- lmer( MArea ~ scale(year) + Sex + East.versus.West + scale(lat) + (1|Collection) + scale(daylength.cycle) + (1|Island.State), data = wings.contemporary)

summary(wings.contemporary.area)

####Get similar result when comparing just the subset of data from 1974-present in the breeding areas. This further corroborates the result of an increase in size through time.

####Compare the same models using wing length as the response variable

flockhart.length.model1 <- lm(WingLength_mm ~ scale(WinterYear) + scale(dist.centroid.km) + Sex, flockhart.data)

wings.contemporary.length <- lmer( MLength ~ scale(year) + Sex + East.versus.West + scale(lat) + (1|Collection) + scale(daylength.cycle) + (1|Island.State), data = wings.contemporary)

summary(flockhart.length.model1)
summary(wings.contemporary.length)

####Again, both models pull out a significant increase through time. However, with these data, the pattern seems to be more pronounced in the Flockhart overwintering data than for the data from non-overwintering sites. 

##Explore patterns of how host plant identity may influence wing morphology

####Here, we are interested in the hypothesis that changes through time in wing morphology could reflect changes through time in the assemblage of available milkweed plants.  Specifically, species on which monarchs have been shown to grow large (e.g. Aslcepias syriaca) may have become relatively more common over the sampling window.  This could generate an increase in size through time that has little or nothing to do with natural selection favoring large butterflies.

####Load data from 2017 adaptation experiment that reared caterpillars from multiple maternal families on a range of host plants. This dataset comprises 232 adult butterflies from four populations (eastern North America, western North America, Hawaii, Australia) reared on five milkweed species (A. syriaca, A. fascicularis, A. curassavica, G. physocarpus, G. fruticosus).

adaptation.experiment.2017 <- read.csv(file = '~/Documents/grad_school/pacific_monarchs/local_adaptation/assignments_plant.csv')

head(adaptation.experiment.2017)

adaptation.experiment.2017$Sex[adaptation.experiment.2017$Sex=='']=NA
adaptation.experiment.2017$Sex = droplevels(adaptation.experiment.2017$Sex)

####Create column for mean wing area

adaptation.experiment.2017$MArea <- rowMeans(cbind(adaptation.experiment.2017$Larea,adaptation.experiment.2017$Rarea),na.rm = T)


####Extract butterfly population ID from monarch column

adaptation.experiment.2017$Monarch.population <- as.factor((str_extract(adaptation.experiment.2017$Monarch, "[aA-zZ]+"))) #extract text from ID

####Combine two Gomphocarpus species into one, since they are so functionally similar. Then select only observations where adults actually eclosed and were measurable.

adaptation.experiment.2017$Species <- revalue(adaptation.experiment.2017$Species, c("GOPH"="Gomphocarpus spp.", "GOFR"="Gomphocarpus spp.", "ASFA"="A. fascicularis", "ASYR"="A. syriaca", "ASCU"="A. curassavica"))

adaptation.experiment.2017 <- adaptation.experiment.2017[!is.na(adaptation.experiment.2017$Sex),]


###Figure 3

####Plot of raw data: butterfly wing area as a function of host plant ID

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

ggplot(adaptation.experiment.2017, aes(x = Species, y = MArea*100))+
  geom_boxplot(fill='forestgreen')+
  theme_bw()+
  labs(y = bquote('Wing Area ('*mm^2*')'))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14, face = 'italic', angle = 60, hjust = 1), axis.text.y = element_text(size = 14))+
  theme(axis.title = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 1)+
  theme(legend.position = 'none')

adaptation.experiment.2017$Monarch.population <- revalue(adaptation.experiment.2017$Monarch.population, c("AU"="Australia", "CA"="California", "HI"="Hawaii","NA"="E. North Amer."))

ggplot(adaptation.experiment.2017, aes(x = Monarch.population, y = MArea, fill = Monarch.population))+
  geom_boxplot()+
  theme_bw()+
  labs(y = bquote('Wing Area ('*cm^2*')'))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14, face = 'italic', angle = 60, hjust = 1))+
  theme(axis.title = element_text(size = 16))+
  stat_summary(fun.data = give.n, geom = "text", vjust = 1)+
  scale_fill_manual(values = c('red','orange','blue','forestgreen'))+
  theme(legend.position = 'none')

####Seems to indicate that butterflies across populations grow larger on A. syriaca. Now run a linear model to test for the effects of host species.

host.plant.area.model1 <- lmer(MArea ~ Monarch.population + Species + (1|Monarch) + Sex + (1|GH) + (1|Group), adaptation.experiment.2017)

summary(host.plant.area.model1)

Anova(host.plant.area.model1)

####Pulls out a highly significant overall effect of species.  By contrast, the impacts of butterfly sex and even monarch population are relatively modest, suggesting that wing area does indeed have a very strong environmental component.

####Pairwise comparisons between host plants and monarch populations

species.comparison <- glht(host.plant.area.model1, mcp(Species="Tukey"))

summary(species.comparison)

pop.comparison <- glht(host.plant.area.model1, mcp(Monarch.population="Tukey"))

summary(pop.comparison)

####Asclepias syriaca produces butterflies that are signficant larger than butterflies reared on A. curassavica and A. fascicularis, and marginally larger than butterflies reared on Gomphocarpus spp.

###Check to see if host species also affects wing roundness

adaptation.experiment.2017$MPerimeter <- rowMeans(cbind(adaptation.experiment.2017$Lperimeter,adaptation.experiment.2017$Rperimeter),na.rm = T)

adaptation.experiment.2017$roundness <- (4*pi*adaptation.experiment.2017$MArea)/(adaptation.experiment.2017$MPerimeter)^2

host.plant.roundness.model1 <- lmer(roundness ~ Species + (1|Monarch) + Sex + (1|GH) + Monarch.population + (1|Group), adaptation.experiment.2017)

summary(host.plant.roundness.model1)

Anova(host.plant.roundness.model1)

levels(adaptation.experiment.2017$Monarch)

####In contrast to wing area, host species does not explain a signficant portion of variation in wing roundness. The only significant predictor here was butterfly sex, with females slightly rounder than females (at least when reared under "summer" conditions in a common greenhouse environment).

###Figure and tables for supplemental information

###Table S1:

####Summary of specimens by collection, including date ranges; first, condense some collections (treated separately because of differnet imaging techniques, but group together for table)

levels(wings$Collection)

wings.table <- wings

wings.table$Collection <- revalue(wings.table$Collection, c("AMNH" = "American Museum of Natural History", "Bishop_Museum" = "Bishop Museum of Hawaii", "Bohart" = "Bohart Museum (UC Davis)", "Cal_Academy" = "California Academy of Sciences", "California_Miscellaneous" = "Miscellaneous", "CUIC" = "Cornell University Insect Collection", "David_M" = "Mexican Overwintering - David M.", "Essig" = "Essig Museum (UC Berkeley)", "Graham_Owen_PC" = "Miscellaneous", "Harvard" = "Harvard Museum of Comparative Zoology", "LA_expt_2017" = "Miscellaneous", "LACMNH" = "Los Angeles County Museum of Natural History", "LHY" = "Yang et al. (2016), Ecography", "Li_et_al_2016" = "Li et al. (2016), Animal Migration", "McGuire_Center_1" = "McGuire Center (University of Florida)", "McGuire_Center_2" = "McGuire Center (University of Florida)", "Mexico" = "Miscellaneous", "MZ" = "Personal Collection: Myron Zalucki", "Peabody_Museum" = "Peabody Museum (Yale University)", "Peabody_Museum_X1" = "Peabody Museum (Yale University)", "Riverside" = "UC Riverside Insect Collection", "Zalucki_Paine_Malcolm" = "Miscellaneous"))

collection.table <- as.data.frame(with(wings.table, table(Collection)))

min <- vector(length = length(levels(wings.table$Collection)))
max <- vector(length = length(levels(wings.table$Collection)))


for(i in 1:length(levels(wings.table$Collection))){
  min[i] <- min(wings.table[wings.table$Collection == levels(wings.table$Collection)[i],]$year,na.rm=T)
}
for(i in 1:length(levels(wings.table$Collection))){
  max[i] <- max(wings.table[wings.table$Collection == levels(wings.table$Collection)[i],]$year,na.rm=T)
}

(collection.table$range <- paste(min, max, sep = "-"))

(colnames(collection.table)[2:3] <- c("Sample Count", "Collection Date Range"))

collection.table

write.csv(collection.table, file = "~/Downloads/Freedman_Dingle_TableS1.csv")

###Table S2: Complete specimen repository

####Here, just use raw data file

###Table S3: AIC model comparison table

(wing.area.model.comparison <- aictab(list(model1.wing.area,model2.wing.area,model3.wing.area,model4.wing.area,model5.wing.area,model6.wing.area,model7.wing.area), modnames = c('model1.wing.area','model2.wing.area','model3.wing.area','model4.wing.area','model5.wing.area','model6.wing.area','model7.wing.area')))

write.csv(wing.area.model.comparison, file="~/Downloads/model_comparison.csv")

####Brief examination of correlation between wing length and area, for figure S2

#generate a correlation plot between all measurements

corr.vec <- c("MLength","MWidth", "MArea", "MPerimeter", "aspectratio", "roundness")

wings.correlation <- wings[, names(wings) %in% corr.vec]

chart.Correlation(wings.correlation)


####Explore relationship between body size and wing area from host plant experiment


ggplot(adaptation.experiment.2017, aes(x = Weight, y = MArea*100))+
  geom_point(aes(col = Species))+
  geom_smooth(method = 'lm', col = 'black', lty = 2)+
  scale_color_manual(values = c('red','orange','blue','purple'))+
  theme_bw()+
  labs(x = "Wet Body Mass (mg)", y = bquote('Wing Area ('*mm^2*')'))+
  theme(legend.position = c(0.75,0.15))+
  theme(legend.box.background = element_rect(colour = "black", size = 1))+
  theme(legend.title.align=0.5)+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14))+
  theme(legend.title = element_text(size =16), legend.text = element_text(size = 14, face = 'italic'))
  
####Add figure for wing measurements

smoothing <- read.csv('~/Downloads/wing_smoothing.csv')

ggplot(smoothing, aes(x= -X, y =-Y))+
  geom_polygon(fill = 'grey50', color = 'orange', size =1, lty = 3)+
  coord_fixed()+
  theme_blank()+
  theme(axis.text = element_blank(), axis.title = element_blank())
  



