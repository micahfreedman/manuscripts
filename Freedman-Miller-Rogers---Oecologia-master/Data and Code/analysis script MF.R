############################################################################
## Script for analysis of Micah's invasive ant and honeydew-producing insect project
## Authors: Micah Freedman & Haldre Rogers
## last updated 02Jan2018
#############################################################################

#read in data

antdata <- read.table("~/Dropbox/Invasive Ants and Scales/data_and_scripts/ant_hpi.txt", header=T, na.string=c("na", "NA", ""))

#load libraries
library(ggplot2)
library(AICcmodavg)
library(lme4)

########## Data exploration and cleanup ##############

str(antdata) #make sure everything is in right format. Leaf number is a factor, change to numeric
antdata$leaf_number<-as.numeric(antdata$leaf_number)
antdata$ants_present_leaves<-as.factor(antdata$ants_present_leaves)
antdata$hpi_status<-as.factor(antdata$hpi_present) #make new factor called hpi_status for figures
levels(antdata$hpi_status)[levels(antdata$hpi_status)=="1"] <- "HPI Present"
levels(antdata$hpi_status)[levels(antdata$hpi_status)=="0"] <- "HPI Absent"
levels(antdata$site) #check sites for mistakes
levels(antdata$species) #check species for mistakes 

#Replace factor level labels for sites, tree species, islands

site.names <- factor(c("Anao","Forbi","Gravel Pit", "Isang", "LADT", "MTR1", "NBlas", "Palii", "Racetrack", "Ritidian", "SBlas")) # vector of new site names

tree.spp <- factor(c("Aglaia","Cynometra","Eugenia","Ficus","Meiogyne","Leucaena","Macaranga","Morinda","Ochrosia opp.", "Ochrosia mar.","Premna","Psychotria","Triphasia")) #vector of new species names (uses some outdated taxonomy, e.g. Guamia; can update if needed...)

island.names <- factor(c("Guam","Rota","Saipan"))

levels(antdata$site) <- site.names #sub in new site names
levels(antdata$species) <- tree.spp #sub in new species names
levels(antdata$island) <- island.names #sub in new island names

#add column for whether tree species are native / non-native

is.native <- ifelse(antdata$species %in% c('Aglaia','Cynometra',
                                          'Eugenia','Ficus','Meiogyne','Macaranga',
                                          'Ochrosia opp.','Ochrosia mar.','Premna','Psychotria'),"native","non-native")
antdata <- cbind(antdata,is.native) 

##### Quick exploration of data

with(antdata, table(hpi_present, ants_present_leaves)) #almost never have hpi present without ants

with(antdata, table(island, ants_present_leaves, hpi_present)) #ants are frequently present without hpi's on saipan, and situations with neither ants nor hpi's are common across all three islands, but most situations of hpi presence occur on Guam (42 out of 59 instances)
with(antdata, table(species, ants_present_leaves, hpi_present))
with(antdata, table(is.native, hpi_present,ants_present_leaves))

ggplot(antdata, aes(x=as.factor(ants_present_leaves)))+
  geom_bar()+
  facet_grid(island~hpi_status, scales="free_y")+
  xlab("Ant Presence")+
  ylab("# of Observations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Full Dataset")

#Looks like majority of "HPI present" observations come from Guam; also appears to be a strong association between presence of ants and HPIs, regardless of island

#check to see how well represented species are between islands and sites
with(antdata, table(island,species)) #a few instances where no individuals of a particular species were sampled from a particular island; besides Ficus, all species have >7 overall observations

##### Stats #####

#center and scale leaf number (continuous predictor)

#use hpi presence as response variable; first try models with all interaction terms and no site-level term

antdata$leaf_number <- scale(antdata$leaf_number)

model1a <- glm(hpi_present ~ ants_present_leaves*species*island + leaf_number, family = binomial, data=antdata)
summary(model1a) #model runs, but error estimates are overinflated
model2a <- glm(hpi_present ~ ants_present_leaves*species + island + leaf_number, family = binomial, data = antdata)
summary(model2a) #same issue; try with no interaction terms
model3a <- glm(hpi_present ~ ants_present_leaves + species + island + leaf_number, family = binomial, data = antdata) #Ficus seems to be causing trouble; remove from dataset and try full models again

#drop Ficus as a species from dataset (only 7 observations)

antdata <- subset(antdata, antdata$species != "Ficus")
antdata$species <- droplevels(antdata$species) #drop factor level corresponding to Ficus

#try models again with full set of interaction terms and site as a random effect

require(optimx)

model2a <- glmer(hpi_present ~ ants_present_leaves*species*island + leaf_number + (1|site), family = binomial, data=antdata, control = glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))) #ants_present is extremely strong predictor, but model fails to converge even with alternative optimizer

m2a.init <- getME(model2a,c("theta","fixef")) #try adjusting theta
model2a.2 <- update(model2a, start = m2a.init) #still fails to converge

model2b <- update(model2a, form = ~ ants_present_leaves+species*island + leaf_number + (1|site)) #model with only interaction between species and island 

m2b.init <- getME(model2b,c("theta","fixef")) #again try with adjusted theta
model2b.2 <- update(model2b, start = m2b.init) #still fails to converge

model2c <- update(model2a, form = ~ ants_present_leaves*species+island + leaf_number + (1|site))

m2c.init <- getME(model2c,c("theta","fixef"))
model2c.2 <- update(model2c, start = m2c.init) #still fails to converge

model2d <- update(model2a, form = ~ ants_present_leaves*island + species + leaf_number + (1|site)) #this one actually does converge

model2e <- update(model2a, form = ~ island + ants_present_leaves + species + leaf_number + (1|site)) #likewise converges

AIC(model2a,model2b,model2c,model2d,model2e) #best model seems to be 2e, which has all covariates but no interaction terms

require(AICcmodavg)

(model2.comp = aictab(list(model2a,model2b,model2c,model2d,model2e), modnames=c("model2a","model2b","model2c","model2d","model2e"))) #create AIC comparison table

write.table(model2.comp, "~/Desktop/model1.csv", sep=",") #write output of AICc table to .csv file for figure S3

#indeed, model2e receives 87% of the model weight

summary(model2e)

require(sjPlot)

plot_model(model2e) #plot fixed effects
sjp.glmer(model2e) #plot random (site level) effects

#generate confidence intervals for pairwise differences between factors for island and species

require(multcomp)

(island.2e <- glht(model2e, linfct = mcp(island = "Tukey")))
confint(island.2e) #confidence intervals for Guam-Saipan comparison do not overlap 0, suggesting signficantly fewer HPIs on Saipan than on Guam
summary(island.2e)

(spp.2e <- glht(model2e, linfct = mcp(species = "Tukey")))
confint(spp.2e) 
summary(spp.2e) #lone significant pairwise species-level difference is Cynometra < Macaranga

ants.2e <- glht(model2e, linfct = mcp(ants_present_leaves = "Tukey"))
confint(ants.2e) 
summary(ants.2e)

#get odds ratios for fixed effects and their confidence intervals

exp(confint(model2e,parm="beta_",method="Wald"))

p2<-sjp.glmer(model2e, type = "eff", facet.grid = FALSE, 
        show.ci = TRUE, prnt.plot = FALSE) #generate plots of all effects

p2$plot.list[1] #show plot for predicted HPI probablities by island, with 95% CIs

p2$data[1:3,] #extract confidence intervals produced

###### test same models with ant presence as the response variable and HPI presence as a predictor

model5a <- glmer(ants_present_leaves ~ hpi_status*species*island + leaf_number + (1|site), family = binomial, data=antdata, control = glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))) #ants_present is extremely strong predictor. Model fails to converge though, likely because of too many predictors for the relatively small number of data points
model5b <- update(model5a, form = ~ hpi_status+species*island + leaf_number + (1|site)) #model also does not converge
model5c <- update(model5a, form = ~ hpi_status*species+island + leaf_number + (1|site)) #likewise
model5d <- update(model5a, form = ~ hpi_status*island+species + leaf_number +(1|site)) #model does converge
model5e <- update(model5a, form = ~ hpi_status + island + species + leaf_number + (1|site))

AIC(model5a,model5b,model5c,model5d,model5e) #again, model with no interaction terms performs best

(model5.comp = aictab(list(model5a,model5b,model5c,model5d,model5e), modnames=c("model5a","model5b","model5c","model5d","model5e"))) #create AIC comparison table for this set of models

write.table(model5.comp, "~/Desktop/model2.csv", sep=",") #write output of AICc table to .csv file for figure S3

#again, model 5e gets 82% of the model weight

summary(model5e)

#generate confidence intervals 

island.5e <- glht(model5e, linfct = mcp(island = "Tukey"))
confint(island.5e) #no signficant pairwise differences between islands, although Guam-Rota is close
summary(island.5e)

species.5e <- glht(model5e, linfct = mcp(species = "Tukey"))
confint(species.5e) #no significant pairwise differences between species
summary(species.5e)

ants.5e <- glht(model5e, linfct = mcp(hpi_status = "Tukey"))
confint(ants.5e) 
summary(ants.5e)

p5<-sjp.glmer(model5e, type = "eff", facet.grid = FALSE, 
              show.ci = TRUE, prnt.plot = FALSE)

p5$plot.list[2]

p5$data[3:5,]

#### use counts of ants ascending/descending trunk as a response variable, with HPI status of sampled leaves, tree species, island as predictors and a Poisson distribution

hist(antdata$sum_ants,breaks=40)

ggplot(antdata, aes(x = sum_ants, col = island))+
  geom_density() #shows that all are somewhat zero-inflated but that this is least pronounced on Guam, which also has all of the higher outliers

#initial attempts show that Poisson models both are highly overdispersed, so try alternate distributions...

#try zero inflated poisson mixture model

require(pscl)

zin.m8 <- zeroinfl(sum_ants ~ island + species + hpi_status, data = antdata) #can use a zero-inflated model to separately model zero and non-zero responses, and then use Poisson distribution for the remaining count data; however, this approach does not support mixed effects models, so site cannot be included as a random intercept term in this model

summary(zin.m8) #take-home message from zero inflated model: compared to Guam, Rota has much higher proportion of zero values, though Saipan-Guam contrast for zero-nonzero is not actually different; among trees that do have ants present on trunks, Guam has significantly more than either Rota or Saipan

require(car)
require(MASS)

library(glmmADMB) #try glmmADMB package for zero-inflated model that does inlcude random effects

model8.zinfl2 <- glmmadmb(sum_ants ~ hpi_status + island + species + (1|site), data = antdata, zeroInflation = TRUE, family="poisson") #takes about 15 seconds to run

summary(model8.zinfl2) #ant sum significantly lower on Rota, modestly lower on Saipan

island.8z <- glht(model8.zinfl2, linfct = mcp(island = "Tukey"))
confint(island.8z)
summary(island.8z) #According to this model, Guam has significantly more ants on the trunks than Rota, but not Saipan

species.8z <- glht(model8.zinfl2, linfct = mcp(species = "Tukey"))
confint(species.8z)
summary(species.8z) #suggests a decent number of species-level pairwise differences in ant abundance

#try separate set of models using glm.nb. Here, models have issues with site being included as a random effect and run indefinitely. Instead, try models that do not inlclude site.

model8a = glm.nb(sum_ants ~ island*species*hpi_status ,data = antdata) 
model8b = update(model8a, form = ~ island*hpi_status + species)
model8c = update(model8a, form = ~ island*species + hpi_status)
model8d = update(model8a, form = ~ hpi_status*species + island)
model8e = update(model8a, form = ~ species+hpi_status+island)

(model8.comp = aictab(list(model8a,model8b,model8c,model8d,model8e), modnames=c("model8a","model8b","model8c","model8d","model8e"))) #does not run because aictab not yet designed for glm.nb objects

bbmle::AICtab(model8a,model8b,model8c,model8d,model8e) #model 8b is preferred slightly

summary(model8b)
summary(model8e) #give very similar answers; for sake of consistency use model8e, since it has the same form as previous models without interaction terms between covariates

AIC(model8e, zin.m8, model8.zinfl2) #compare AIC scores across model classes (maybe not advisable since they use different likelihood functions); model8e (standard negative binomial) is clearly favored

#try including site as a random effect in model8e to see if it can run

model8f <- glmer.nb(sum_ants ~ island + species + hpi_status + (1|site), data = antdata) #takes about 30 seconds to run, but does not converge

m8f.init <- getME(model8f,c("theta","fixef"))
model8f.2 <- update(model8f, start = m8f.init, control=glmerControl(optCtrl=list(maxfun=20000))) #converges quickly

summary(model8f.2) #output nearly identical to model8f; 

AIC(model8e, model8f.2) #given the choice of these two models, model including site as a random effect is preferred, so stick with it

model8k <- glmer.nb(sum_ants ~ island*hpi_status + species + (1|site), data = antdata) #compare other model that performed well, but this time include site
m8k.init <- getME(model8k,c("theta","fixef"))
model8k.2 <- update(model8k, start = m8k.init, control=glmerControl(optCtrl=list(maxfun=20000))) #still doesn't quite converge

AIC(model8f.2,model8k.2) #with addition of site term, now model with no interactions is favored after all

AICc(model8f.2)
AICc(model8k.2)

#get confidence intervals

island.8f2 <- glht(model8f.2, linfct = mcp(island = "Tukey"))
confint(island.8f2) #signficantly more ants on Guam than either Rota or Saipan, but no difference between Saipan and Rota
summary(island.8f2)

species.8f2 <-glht(model8f.2, linfct = mcp(species = "Tukey"))
confint(species.8f2) #no signficant pairwise differences between islands
summary(species.8f2)

ants.8f2 <- glht(model8f.2, linfct = mcp(hpi_status = "Tukey"))
confint(ants.8f2) 
summary(ants.8f2)

p8<-sjp.glm(model8f.2, type = "eff", facet.grid = FALSE, 
              show.ci = TRUE, prnt.plot = FALSE)

p8$plot.list[1]

p8$data[1:3,]

#############

#Model output for figures

#############

### Fig 1a - Probability of HPIs across islands

antdata$mod2e <- fitted(model2e) #create new column with fitted values for model 2e

fig1a.values <- p2$data[1:3,]

ggplot(antdata, aes(x = island, y = mod2e))+
  geom_point(size=0.5, col = 'white')+
  ylab("Probability of HPIs")+
  stat_summary(fun.y = mean, geom = "point", aes(col = island), size=5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar",aes(col=island),size=1.2, width = 0.2)+
  theme_classic()+
  theme_classic()+
  scale_color_manual(values = c('orange','blue','cyan'))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none") #while this figure works, doesn't make full sense, because the means and confidence intervals are derived after the logit tranform, so the values aren't appropriate

(Fig1a<-ggplot(fig1a.values, aes(x = labels, y = y, col = labels))+
  geom_point(size=4)+
  ylab("Probability of HPIs")+
  theme_classic()+
  theme_classic()+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  scale_color_manual(values = c('orange','blue','cyan'))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")+
  coord_cartesian(ylim = c(0,1))+
  annotate("text", x = 1, y = 0.7, label = 'A', size = 6)+
  annotate("text", x = 2, y = 0.6, label = 'AB', size = 6)+
  annotate('text', x = 3, y = 0.35, label = 'B', size = 6))+
  border(color = "black", size = 0.8, linetype = NULL)
  #use this figure instead

## Figure 1b - probability of ants across islands

fig1b.values <- p5$data[3:5,]  

(Fig1b<-ggplot(fig1b.values, aes(x = labels, y = y, col=labels))+
  geom_point(size=4)+
  ylab("Probability of Ants")+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  theme_classic()+
  scale_color_manual(values = c('orange','blue','cyan'))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")+
  coord_cartesian(ylim = c(0,1))+
  annotate("text", x = 1, y = 0.95, label = 'A', size = 6)+
  annotate("text", x = 2, y = 0.75, label = 'A', size = 6)+
  annotate('text', x = 3, y = 0.9, label = 'A', size = 6))+
  border(color = "black", size = 0.8, linetype = NULL)
  
fig1c.values <- p8$data[1:3,]

## Figure 1c - Counts of ants/minute on tree trunks; use glmer.nb model incorporating site as a random effect

(Fig1c <- ggplot(fig1c.values, aes(x = labels, y = y, col = labels))+
  geom_point(size=4)+
  ylab("Ants/Minute")+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  theme_classic()+
  scale_color_manual(values = c('orange','blue','cyan'))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")+
  annotate("text", x = 1, y = 11, label = 'A', size = 6)+
  annotate("text", x = 2, y = 4.5, label = 'B', size = 6)+
  annotate('text', x = 3, y = 5, label = 'B', size = 6)+
  coord_cartesian(ylim = c(0,11)))

antdata$mod.zin8 <- fitted(zin.m8) #model output and associated figure for first zero-inflated poisson model (using zinb function)

ggplot(antdata, aes(x = island, y = mod.zin8))+
  geom_point(size=0.5, col = 'white')+
  ylab("Ants/Minute")+
  stat_summary(fun.y = mean, geom = "point", aes(col = island), size=5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar",aes(col=island),size=1.2, width = 0.2)+
  theme_classic()+
  theme_classic()+
  scale_color_manual(values = c('orange','blue','cyan'))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=20))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none") #looks similar to previous figure

antdata$mod8z <- fitted(model8.zinfl2) #model output for second zero-inflated poisson model; can see that the estimate for Guam has notably shrunk towards 0, since this class of model disregards many of the zeros in the data (which are likely true zeros) as process error

ggplot(antdata, aes(x = island, y = mod8z))+
  geom_point(size=0.5, col = 'white')+
  ylab("Ants/Minute")+
  stat_summary(fun.y = mean, geom = "point", aes(col = island), size=5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar",aes(col=island),size=1.2, width = 0.2)+
  theme_classic()+
  theme_classic()+
  scale_color_manual(values = c('orange','blue','cyan'))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=20))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")

#### assemble Figure 1

require(ggpubr)

ggarrange(Fig1a, Fig1b, Fig1c,
          labels = c("(a)", "(b)", "(c)"),
          ncol = 3, nrow = 1)

### Fig S1 - model fitted predictions for each species; HPI probability

fig.s1 <- p2$data[6:17,]

(s1a <- ggplot(fig.s1, aes(x = labels, y = y))+
  geom_point(size=4)+
  ylab("Probability of HPIs")+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  theme_classic()+
  xlab('Tree Species')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12,face = "italic"))+
  theme(axis.title.x = element_text(size=16, vjust = -1.5))+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")+
  theme(strip.text.x = element_text(size = 16)))

##with ant presence as response

fig.s2 <- p5$data[6:17,]

(s2a <- ggplot(fig.s2, aes(x = labels, y = y))+
  geom_point(size=4)+
  ylab("Probability of Ants")+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  theme_classic()+
  xlab('Tree Species')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12,face = "italic"))+
  theme(axis.title.x = element_text(size=16, vjust = -1.5))+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")+
  theme(strip.text.x = element_text(size = 16)))

#### species-level differences, with ant abundance as response

fig.s3 <- p8$data[4:15,]

(s3a <- ggplot(fig.s3, aes(x = labels, y = y))+
  geom_point(size=4)+
  ylab("Ant Abundance")+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  theme_classic()+
  xlab('Tree Species')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12,face = "italic"))+
  theme(axis.title.x = element_text(size=16, vjust = -1.5))+
  theme(axis.title.y=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(legend.position="none")+
  theme(strip.text.x = element_text(size = 16)))

#bind all three together (same x axis on all of them)

ggarrange(s1a, s2a, s3a,
          labels = c("(a)", "(b)", "(c)"),
          ncol = 1, nrow = 3) #looks weird 

##use facet wrap instead

fig.s1$panel <- "Probability of HPIs"
fig.s2$panel <- "Probability of Ants"
fig.s3$panel <- "Ant Abundance"

plotteroo <- rbind(fig.s1,fig.s2,fig.s3)

ggplot(transform(plotteroo,panel=factor(panel,levels=c('Probability of HPIs','Probability of Ants','Ant Abundance'))), aes(x = labels, y = y))+
  geom_point(size=4)+
  facet_wrap(~panel, scales = 'free_y', ncol = 1)+
  geom_errorbar(aes(x=labels, ymax=conf.low, ymin=conf.high),width=0.2, size=1.2)+
  theme_classic()+
  xlab('Tree Species')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=12,face = "italic"))+
  theme(axis.title.x = element_text(size=16, vjust = -1.5))+
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  theme(strip.text.x = element_text(size = 16))

#######################

####  ant community statistics and figures

source("~/Dropbox/Invasive Ants and Scales/data_and_scripts/biostats.R")

require(vegan)

ants <- read.table(file="~/Dropbox/Invasive Ants and Scales/data_and_scripts/matrix_all_ants.txt",header = T)

#the above is a dataframe that includes all 15 species of ants that were recorded on leaves; for this particular analysis, I'm interested in abundances (as opposed to just presence/absence of individual taxa) and so will be using Bray-Curtis dissimilarity

sppUse<-as.data.frame(ants[,(4:ncol(ants))])
#create dataframe of only abundances

scree <- nmds.scree(x=sppUse, distance="bray", k=10, autotransform=F, trymax=15)
#generates a scree plot; generally here you look for an inflection point at a given k number of inflection points, and then use k+1 ordination axes; no distinct "elbow" here, but starts leveling off around 5, so maybe use k=6

ants.bc.mat <- vegdist(decostand(sppUse, 'total')) #generate distance matrix for stats

(full_model = adonis(decostand(ants.bc.mat, 'total') ~ island + site + species, data = ants)) #run model to determine which factors explain variance in distance matrix
(nested = adonis(decostand(ants.bc.mat, 'total') ~ island/site + species, data = ants)) #alternative model with site nested within island; results very similar

library(RVAideMemoire)

pairwise.perm.manova(ants.bc.mat,ants$island,nperm=1000, p.method = 'bonferroni') #Guam-Rota and Guam-Saipan comparisons are significantly different, but not Saipan-Rota comparison

sppNMDS1<-metaMDS(comm=sppUse, distance="bray", k=2, autotransform=F, trymax=500)
#for first try, use only two ordination axes and simulate over 500 random starts
sppNMDS1

stressplot(sppNMDS1)

sppNMDS2<-metaMDS(comm=sppUse, distance="bray", k=6, autotransform=F, trymax=2500)
#for comparison, use six ordination axes and simulate over 2500 starts; still does not converge
sppNMDS2

stressplot(sppNMDS2)

data.scores <- as.data.frame(scores(sppNMDS1))
data.scores$species <- ants$species
data.scores$island <- ants$island
data.scores$site = ants$site
head(data.scores)

require(plyr)

find_hull <- function(data.scores) data.scores[chull(data.scores$NMDS1, data.scores$NMDS2), ]
hulls.island <- ddply(data.scores, "island", find_hull)
hulls.site = ddply(data.scores, "site", find_hull)

plot.island <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, colour=island), fill=island) +
  geom_jitter(data=data.scores,aes(x=NMDS1,y=NMDS2),size=4) + 
  geom_polygon(aes(fill = island),data = hulls.island, alpha = 0.1, size =1,show.legend = FALSE) +
  labs(x = "NMDS1", y = "NMDS2")+
  theme_classic()+
  theme(plot.title = element_text(size=18))+
  scale_color_manual("Island\n",labels = c("Guam", "Rota", "Saipan"), values = c("orange", "blue", "cyan"))+
  scale_fill_manual(values = c('orange','blue','cyan'))+
  ylab('NMDS Axis 2')+
  xlab('NMDS Axis 1')+
  theme(axis.title.x = element_text(size =12),axis.title.y = element_text(size=12))+
  theme(legend.title = element_text(size =14), legend.text = element_text(size=12))
plot.island #figure for NMDS using only two ordination axes

#generate figure for k=6 ordination

data.scores2 <- as.data.frame(scores(sppNMDS2))
data.scores2$species <- ants$species
data.scores2$island <- ants$island
data.scores2$site = ants$site
head(data.scores2)

find_hull <- function(data.scores2) data.scores2[chull(data.scores2$NMDS1, data.scores2$NMDS2), ]
hulls.island <- ddply(data.scores2, "island", find_hull)
hulls.site = ddply(data.scores2, "site", find_hull)

plot.island2 <- ggplot(data = data.scores2, aes(x = NMDS1, y = NMDS2, colour=island), fill=island) +
  geom_jitter(data=data.scores2,aes(x=NMDS1,y=NMDS2),size=4) + 
  geom_polygon(aes(fill = island),data = hulls.island, alpha = 0.1, size =1,show.legend = FALSE) +
  labs(x = "NMDS1", y = "NMDS2")+
  theme_classic()+
  theme(plot.title = element_text(size=18))+
  scale_color_manual("Island\n",labels = c("Guam", "Rota", "Saipan"), values = c("orange", "blue", "cyan"))+
  scale_fill_manual(values = c('orange','blue','cyan'))+
  ylab('NMDS Axis 2')+
  xlab('NMDS Axis 1')+
  theme(axis.text = element_text(size =16),axis.title = element_text(size=16))+
  theme(legend.title = element_text(size =16), legend.text = element_text(size=16))
plot.island2 #plot for ordination using all six axes (suggests Rota is basically a subset of Saipan, with Guam largely not overlapping either Rota or Saipan)

plot.island2 + theme(legend.position = c(0.1,0.1)) + theme(legend.title = element_blank()) + ylab("Ant Community NMDS Axis 2") + xlab("Ant Community NMDS Axis 1") #reposition legend and remove title

##### based on request for reviewer 1, generate ant community matrix omitting T. albipes and see how islands differ

ants_no.techno <- ants[,-4]

sppUse3<-as.data.frame(ants_no.techno[,(4:ncol(ants_no.techno))])

#drop rows with all 0s

sppUse3 <- sppUse3[-which(rowSums(sppUse3)==0),]

#create dataframe of only abundances

scree3 <- nmds.scree(x=sppUse3, distance="bray", k=5, autotransform=F, trymax=10)

##cannot do analysis without T. albipes, because without it we only have 79 rows with data present, compared to 137 rows for the whole dataset

##################

anosim(sppUse, ants$island, permutations = 1000, distance = "bray")

#anosim stands for analysis of similarity and is the commonly used metric for determining whether clusters of interest can be considered significanlty different.  Here, there is indeed significant separation based on island.

anosim(sppUse, ants$site, permutations = 1000, distance = "bray")

#also get signficant difference among sites, but this is less compelling since these should really be nested within island; anosim is also sensitive to including too many comparisons, since this increases the chance that any two clusters will be different; so, it's probably better just to stick to the island-wide comparison, which is pretty clear

anosim(sppUse, ants$site:ants$island, permutations = 1000, distance = "bray")
#nested site within island doesn't change result at all

anosim(sppUse, ants$species, permutations = 1000, distance = "bray")
#tree species ID does not explain a signficant proportion of the variance

##### multivariate statistics for assessing importance of individual ant species in determining community composition

library(mvabund)

ants.mat <- data.matrix(sppUse)

ants.mv <- manyglm(ants.mat ~ island + site + species, data = ants, family = 'negative.binomial') 

anova(ants.mv) #warning: takes about 3 minutes to run; suggests that island, followed by site, followed by species are all important determinants of the distance matrix (order concordant with results suggested by anosim)

ants.mv #shows strongest determinants of ant community composition: 1. T. albipes (by far) 2. T. melanocephalum 3. A. gracillipes






######### phylogeny ###########

require(ape)
require(phytools)
marianas = read.newick(file="~/Dropbox/Invasive Ants and Scales/data_and_scripts/guam.phylo.txt")

marianas[[3]][11] <- "Ochrosia_oppositifolia"
marianas[[3]][13] <- "Meiogyne_cylindrocarpa"

m2 = collapse.singles(marianas)
plot(m2)
axisPhylo(1, las = 2)

############
#retrospective power analysis, as per the request of reviewer 2
############

pwr.2p2n.test(h = 0.47, n1 = 103, sig.level = 0.05, power = 0.8) #based on 103 sampled branches from Guam, significance level of 0.05, and the log response ratio given in Mooney et al. 2010

#suggests we would need at least 54 samples from Rota to reliably detect a difference




