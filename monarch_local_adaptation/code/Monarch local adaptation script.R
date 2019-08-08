#### Statistical analysis for local adaptation experiment

#load all libraries

library(lme4)
library(car)
library(sjPlot)
library(ggplot2)
library(lmerTest)
library(car)
library(plyr)
library(multcomp)
library(sjstats)
library(cvequality)
library(MuMIn)
library(emmeans)

#load data

setwd('~/Documents/grad_school/Manuscripts/Monarch local adaptation/')

lad <- read.csv("./Data/cats.temp.csv")

lad$Mon.Pop <- factor(lad$Mon.Pop, levels(lad$Mon.Pop)[c(4,5,2,3,1,6)])

head(lad)
str(lad)

#convert year, usage, plant.id into factors

lad[,c('Year','Usage','Plant.ID')] <- lapply(lad[,c('Year','Usage','Plant.ID')], factor)

#change all values coded as GOFR to GOPH

lad$Species <- factor(gsub('GOFR', 'GOPH', lad$Species))

#add column for sympatric/allopatric status

lad$sym.allo <- with(lad, ifelse(Species %in% c('AINC','ASYR') & Mon.Pop == 'ENA', 'sympatric',
                          ifelse(Species %in% c('ASPEC','ASFA') & Mon.Pop == 'CA', 'sympatric',
                          ifelse(Species == 'GOPH' & Mon.Pop %in% c('AU','HI'), 'sympatric',
                          ifelse(Species == 'ASCU' & Mon.Pop %in% c('GU','PR'), 'sympatric', 
                          'allopatric')))))

#add column for monarch maternal family, since some IDs are repeated between years

lad$Mon.ID <- with(lad, paste(Mon.Pop, ID, Year, sep = '-'))

#add another column that describes the number of larvae placed on each plant; first create a column that captures each individual monarch family x plant genotype combination

lad$reps <- with(lad, paste(Species, Plant.ID, Mon.ID, sep = "-"))
length(unique(lad$reps)) #828 unique combinations

larvae <- aggregate(Caterpillar.. ~ reps, FUN = max, lad) #how many caterpillars within each of these replicates?
colnames(larvae)[2] <- 'Larvae'

lad <- merge(lad, larvae, by = "reps")
lad$Larvae <- as.numeric(lad$Larvae)

lad$Larvae <- ifelse(lad$Larvae > 5, 5, lad$Larvae) #replace the three entries with >5 larvae as 5, just for the sake of treating larvae as a continuous variable

dates <- as.Date(lad$Date, format='%m/%d/%Y')

start_date = as.Date('2018-05-15', tz="UTC") #define earliest date on which caterpillars were added to plants (May 15)

lad$exp.days <- as.numeric(difftime(dates, start_date , units = "days"))

###

#create quick figure to show distribution of when caterpillars were reared across years

mon.pop.cols <- c('blue','cornflowerblue','purple','forestgreen','orange','hotpink')

lad$Mon.Pop <- factor(lad$Mon.Pop, c('AU','HI','GU','ENA','CA','PR'))

FigS1a <- ggplot(lad, aes(x = exp.days))+
  geom_histogram(bins = 40, fill = 'black')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))+
  xlab('Days Since May 15th')+
  ylab('Number of Larvae Measured')

FigS1a

ggsave(file="./Figures/Figure_S1a.tiff", plot=FigS1a, width=7, height=4)

FigS1b <- ggplot(lad, aes(x = exp.days, fill = Mon.Pop))+
  geom_histogram(bins = 40)+
  facet_wrap(~Mon.Pop*Year, nrow = 2)+
  theme_bw()+
  scale_fill_manual(values = mon.pop.cols, guide = F)+
  theme(strip.text = element_text(size = 16), axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))+
  xlab('Days Since May 15th')+
  ylab('Number of Larvae Measured')

FigS1b

ggsave(file="./Figures/Figure_S1b.tiff", plot=FigS1b, width=7, height=4)

###

###

lad$is.migratory <- ifelse(lad$Mon.Pop %in% c('ENA','CA'), 'migratory', 'resident')

#begin statistical analyses, starting with day 8 larval mass

hist(lad$Weight, breaks=100) #long right hand tail, with bimodal peaks that probably correspond to discrete larval instars

hist(log(lad$Weight), breaks = 100) #does produce a generally more normal looking distribution

#first try complicated model with 4-way covariate interaction

model1.cat.growth <- lmer((Weight) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage +  Species*GH*scale(exp.days)*Year, data = lad)

hist(resid(model1.cat.growth), breaks = 100) #looks quite good, although ranef() shows that individual plant IDs are being apportioned a huge amount of the variance that probably should be attributed to either monarch families or populations (e.g. very large fifth instar caterpillars from some pops should not just completely be a byproduct of the plant that they were reared on)

r.squaredGLMM(model1.cat.growth)

ranef(model1.cat.growth)

#now do same model but with no covariate interactions

model2.cat.growth <- lmer((Weight) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad)

r.squaredGLMM(model2.cat.growth)

#now try model 1, but this time with larval mass log transformed

model3.cat.growth <- lmer((log(Weight)) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage +  Species*GH*scale(exp.days)*Year, data = lad)

r.squaredGLMM(model3.cat.growth)

#model 2 but with mass log transformed

model4.cat.growth <- lmer((log(Weight)) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad)

summary(model4.cat.growth)

r.squaredGLMM(model4.cat.growth)

hist(resid(model4.cat.growth), breaks = 100) #also looks good, though perhaps with a slightly longer left hand tail

#model5.cat.growth <- glmer((Weight) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad, family = 'Gamma') #warning: takes about two minutes to run

#hist(resid(model5.cat.growth), breaks = 100)

Anova(model4.cat.growth)

pop.comparison <- glht(model4.cat.growth, mcp(Mon.Pop="Tukey"))
summary(pop.comparison)

milkweed.comparison <- glht(model4.cat.growth, mcp(Species="Tukey"))
summary(milkweed.comparison)

day8.spp <- as.data.frame(emmeans(model4.cat.growth, "Species")) #note: need to exponentiate to transform back to raw values
day8.mon <- as.data.frame(emmeans(model4.cat.growth, "Mon.Pop"))
sym.allo.contrast <- as.data.frame(emmeans(model4.cat.growth, 'sym.allo'))
day8.gh <- as.data.frame(emmeans(model4.cat.growth, "GH"))
day8.yr <- as.data.frame(emmeans(model4.cat.growth, "Year"))
day8.use <- as.data.frame(emmeans(model4.cat.growth, "Usage"))
day8.days <- as.data.frame(emmeans(model4.cat.growth, "exp.days"))


lad_model4 <- lad[-which(is.na(lad$Weight)),]

lad_model4$output <- exp(fitted(model4.cat.growth))

#First, plot raw data on caterpillar mass on day 8; then plot model fitted values showing contrast between sympatric and allopatric

agg.pop <- aggregate(Weight ~ Mon.ID + Species + Mon.Pop + sym.allo, data = lad, FUN = mean) #aggregate raw data so that each maternal family x host plant is a unique data point

agg.pop$Mon.Pop <- factor(agg.pop$Mon.Pop, levels(agg.pop$Mon.Pop)[c(4,5,2,3,1,6)])

Fig4a1 <- ggplot(agg.pop[agg.pop$Weight < 1.4,], aes( x = Mon.Pop, y = Weight, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Day 8 Mass (g)')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  facet_wrap(~Species)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='gold'))+
  annotate('rect', ymin = 0, ymax = 1.25, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 0, ymax = 1.25, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)

ggplot(agg.pop[agg.pop$Species=='GOPH',], aes( x = Mon.Pop, y = Weight, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Day 8 Mass (g)')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle('GOPH')+
  theme(plot.title = element_text(hjust = 0.5, face = 'italic', size = 16))

ggplot(agg.pop[agg.pop$Species %in% c('ASPEC','AINC','ASYR','ASFA'),], 
       aes( x = Mon.Pop, y = Weight, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Day 8 Mass (g)')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  facet_wrap(~Species)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='gold'))

goph.model <- lmer(log(Weight) ~ Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad[lad$Species=='GOPH',])

nam.plant.model <- lmer(log(Weight) ~ Species + is.migratory + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad[lad$Species %in% c('ASYR','ASPEC','AINC','ASFA'),])

Anova(nam.plant.model)

na.anc.der <- as.data.frame(emmeans(nam.plant.model, 'is.migratory'))

na.anc.der$is.migratory <- revalue(na.anc.der$is.migratory, c('resident' = 'derived', 'migratory' = 'ancestral'))

ggplot(na.anc.der, aes( x = is.migratory, y = exp(emmean), col = is.migratory))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
  scale_color_manual(values = c('blue','forestgreen'), guide = F)+
  theme_bw()+
  ylab('Day 8 Mass (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())

Anova(goph.model)

goph.la <- as.data.frame(emmeans(goph.model, 'sym.allo'))

'ggplot(goph.la, aes( x = sym.allo, y = exp(emmean), col = sym.allo))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  ylab('Day 8 Mass (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())'

Fig4a1

ggsave(file="./Figures/Fig4a1.tiff", plot=Fig4a1, width=6, height=7)

#quick alternate figure for talk

model.ancestral1 <- lmer((log(Weight)) ~ Species + is.migratory + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad)

summary(model.ancestral1)

Anova(model.ancestral1)

mig.res.em <- as.data.frame(emmeans(model.ancestral1, 'is.migratory'))
mig.res.em$is.migratory <- revalue(mig.res.em$is.migratory, c('migratory' = 'ancestral', 'resident' = 'derived'))

Fig4a3 <- ggplot(mig.res.em, aes( x = is.migratory, y = exp(emmean), col = is.migratory))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
  scale_color_manual(values = c('blue','forestgreen'), guide = F)+
  theme_bw()+
  ylab('Day 8 Mass (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())

ggsave(file="./Figures/Fig4a3.tiff", plot=Fig4a3, width=3, height=7)


#now plot marginal means summarizing sympatric/allopatric contrast, with error bars correspond to 95% confidence intervals

Fig4a2 <- ggplot(sym.allo.contrast, aes( x = sym.allo, y = exp(emmean), col = sym.allo))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  ylab('Day 8 Mass (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())  

Fig4a2

ggsave(file="./Figures/Fig4a2.tiff", plot=Fig4a2, width=3, height=7)

############

#next run models for survival to day 8

lad$survival = ifelse(is.na(lad$Weight) == TRUE, 0, 1) #add column for survival

#use same model structure as model 4

model1.survival <- glmer(survival ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad, family = 'binomial') #warning: takes about two minutes to run

summary(model1.survival)

Anova(model1.survival)

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
} #define function to convert estimated marginal means back into survival percentages

day8.surv.spp <- as.data.frame(emmeans(model1.survival, "Species"))
day8.surv.pop <- as.data.frame(emmeans(model1.survival, "Mon.Pop"))
day8.surv.sym <- as.data.frame(emmeans(model1.survival, "sym.allo"))
day8.surv.use <- as.data.frame(emmeans(model1.survival, "Usage"))
day8.surv.gh <- as.data.frame(emmeans(model1.survival, "GH"))
day8.surv.yr <- as.data.frame(emmeans(model1.survival, "Year"))

lad_model_survival <- lad #create copy of dataframe

lad_model_survival$output <- fitted(model1.survival) #add fitted survival results back to original dataframe

fitted.survival <- aggregate(output ~ Mon.ID + Species + Mon.Pop + sym.allo + Plant.ID, data = lad_model_survival, FUN = mean)
fitted.survival2 <- aggregate(output ~ Species + Mon.Pop + sym.allo, fitted.survival, mean)
fitted.survival2.sd <- aggregate(output ~ Species + Mon.Pop + sym.allo, data = fitted.survival, FUN = sd)
names(fitted.survival2.sd)[4] <- 'SD'

plot.survival <- merge(fitted.survival2,fitted.survival2.sd)

plot.survival$Mon.Pop <- factor(plot.survival$Mon.Pop, c('ENA','CA','HI','GU','AU','PR'))

Fig4b1 <- ggplot(plot.survival, aes( x = Mon.Pop, y = output, col = sym.allo))+
  geom_point(size = 3)+
  theme_bw()+
  scale_color_manual(values = c('black','red'), guide = F)+
  xlab('Monarch Population')+
  ylab('Proportion Survival to Day 8')+
  geom_errorbar(aes(ymin = output - SD, ymax = output + SD), width = 0.3, size = 1)+
  facet_wrap(~Species)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='gold'))+
  annotate('rect', ymin = 0.3, ymax = 1, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 0.3, ymax = 1, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)

Fig4b1

ggsave(file="./Figures/Fig4b1.tiff", plot=Fig4b1, width=6, height=7)

#now plot sym/allo contrast for survival

plot.sym.allo.surv <- as.data.frame(emmeans(model1.survival, "sym.allo"))

Fig4b2 <- ggplot(plot.sym.allo.surv, aes( x = sym.allo, y = logit2prob(emmean), col = sym.allo))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = logit2prob(asymp.UCL), ymin = logit2prob(asymp.LCL)), width = 0.2, size=1)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  ylab('Proportion Survival to Day 8')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank()) 

Fig4b2

ggsave(file="./Figures/Fig4b2.tiff", plot=Fig4b2, width=3, height=7)

#now generate ancestral/derived contrast

lad$ancestral <- ifelse(lad$Mon.Pop %in% c('ENA','CA'), 'ancestral', 'derived')

modelx1.survival <- glmer(survival ~ Species + ancestral + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad, family = 'binomial') #warning: takes about two minutes to run

summary(modelx1.survival)

Anova(modelx1.survival)

day8.surv.anc <- as.data.frame(emmeans(modelx1.survival, "ancestral"))

plot.anc.surv <- as.data.frame(emmeans(modelx1.survival, "ancestral"))

Fig4b3 <- ggplot(plot.anc.surv, aes( x = ancestral, y = logit2prob(emmean), col = ancestral))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = logit2prob(asymp.UCL), ymin = logit2prob(asymp.LCL)), width = 0.2, size=1)+
  scale_color_manual(values = c('blue','forestgreen'), guide = F)+
  theme_bw()+
  ylab('Proportion Survival to Day 8')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank()) 

Fig4b3

ggsave(file="./Figures/Fig4b3.tiff", plot=Fig4b3, width=3, height=7)

#################################

#now load in adult data

#2017 data

adults.2017 <- read.csv(file = '~/Downloads/2017_adults.csv')
str(adults.2017) #ID is indeed a factor
adults.2017$Group <- as.factor(adults.2017$Group)
#change ID into a variable called Cat.ID to merge with 2018 data
colnames(adults.2017)[1] <- 'Cat.ID'

adults.2018 <- read.csv(file = '~/Downloads/Adults 2018.csv')
str(adults.2018)


#now load caterpillar ID info from other sheet

all.data <- read.csv(file='~/Downloads/experiment_master.csv')
str(all.data)
colnames(all.data)
colnames(adults.2018)
colnames(adults.2017)

##restrict 2017 and 2018 adult data to just wing measurements and Cat.ID

keep.columns <- c('Cat.ID','LLength','LWidth','LArea','LPer','RLength','RWidth','RArea','RPer','HW.Mass','RHWArea')

keep.2017 <- adults.2017[,colnames(adults.2017) %in% keep.columns]
keep.2017$HW.Mass <- NA
keep.2018 <- adults.2018[,colnames(adults.2018) %in% keep.columns]
adult.data <- rbind(keep.2017,keep.2018)

####merge all based on Cat.ID

full <- merge(all.data, adult.data, by = 'Cat.ID', all.x = T)

str(full)

#drop values with no OE entry

adults.lad <- full[-which(is.na(full$OE)),] #leaves behind 1079 records

adults.lad$days_to_pupation <- as.numeric(as.Date(adults.lad$Pupation, format='%m/%d/%y')) - as.numeric(as.Date(adults.lad$Date, format='%m/%d/%y')) + 8

adults.lad$days_to_eclosion <- as.numeric(as.Date(adults.lad$Eclosion, format='%m/%d/%y')) - as.numeric(as.Date(adults.lad$Date, format='%m/%d/%y')) + 8

hist(adults.lad$days_to_eclosion)

hist(adults.lad$days_to_eclosion,breaks = 30)

ggplot(adults.lad, aes( x = days_to_pupation, y = days_to_eclosion))+
  geom_point()+
  geom_smooth(method = 'lm')

head(adults.lad)

#create column for monarch ID;  first drop entries with 'NA' for monarch population

adults.lad <- adults.lad[-which(is.na(adults.lad$Mon.Pop)),]

adults.lad$maternal_family <- with(adults.lad, paste(Mon.Pop, ID, Year, sep = "_"))

length(unique(adults.lad$maternal_family)) #82 unique maternal families

#convert all usage > 2 into 2

adults.lad$Usage <- ifelse(adults.lad$Usage > 2, 2, adults.lad$Usage)

# make usage and year factors

adults.lad$Year <- as.factor(adults.lad$Year)
adults.lad$Usage <- as.factor(adults.lad$Usage)

#create column for adult weight at emergence

adults.lad$emergence_weight <- adults.lad$Total.weight - adults.lad$Env.weight

hist(adults.lad$emergence_weight)

adults.lad$emergence_weight <- ifelse(adults.lad$emergence_weight > 1, NA, adults.lad$emergence_weight)

adults.lad <- adults.lad[adults.lad$emergence_weight < 1,] #drop the one entry with emergence weight greater than one (this observation also has some other wonky data)

hist(adults.lad$emergence_weight)

#create new column to split butterflies into "uninfected" versus "heavily infected"

adults.lad$infection_status <- ifelse(adults.lad$OE >= 2, "infected", "uninfected")

table(adults.lad$infection_status) #overall 345 infected vs. 720 uninfected

table(adults.lad$infection_status, adults.lad$Mon.Pop, adults.lad$Species)

#add another column for binomial variable of infection status

adults.lad$infection_status_binomial <- with(adults.lad, ifelse(infection_status=='infected', 1, 0))

#combine GOFR into GOPH

adults.lad$Species <- revalue(adults.lad$Species, c('GOFR'='GOPH'))

#create new column for local vs. foreign tolerance mechanism

adults.lad$sym.allo <- with(adults.lad, ifelse(Species %in% c('ASYR','AINC') & Mon.Pop == 'ENA', 'sympatric', ifelse(Species %in% c('ASFA','ASPEC') & Mon.Pop == 'CA', 'sympatric', ifelse(Species == 'ASCU' & Mon.Pop %in% c('PR','GU'), 'sympatric', ifelse(Species == 'GOPH' & Mon.Pop %in% c('AU','HI'), 'sympatric', 'allopatric')))))

###

ad.dates <- as.Date(adults.lad$Date, format='%m/%d/%y')

adults.lad$exp.days <- as.numeric(difftime(ad.dates, start_date , units = "days"))

#now run a model looking at days to eclosion and using the same parameters as all of the previous models

hist(adults.lad$days_to_eclosion) #restrict analysis to individuals that pupated within 35 days (this drops 14 observations out of about 1050)

adults.lad$Mon.Pop <- factor(adults.lad$Mon.Pop, c('AU','HI','GU','ENA','CA','PR'))

eclosion.model1 <- lmer(days_to_eclosion ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad[adults.lad$days_to_eclosion < 35,])

summary(eclosion.model1) #including exp.days generates a non-convergence error message, but omitting it does not meaningfully change any estimates

Anova(eclosion.model1)

ecl.spp <- as.data.frame(emmeans(eclosion.model1, 'Species'))
ecl.pop <- as.data.frame(emmeans(eclosion.model1, 'Mon.Pop'))
ecl.sym <- as.data.frame(emmeans(eclosion.model1, 'sym.allo'))
ecl.use <- as.data.frame(emmeans(eclosion.model1, 'Usage'))
ecl.gh <- as.data.frame(emmeans(eclosion.model1, 'GH'))
ecl.yr <- as.data.frame(emmeans(eclosion.model1, 'Year'))
ecl.sex <- as.data.frame(emmeans(eclosion.model1, 'Sex'))


#plot the raw data

agg.eclosion <- aggregate(days_to_eclosion ~ Species + Mon.Pop + sym.allo + maternal_family + Plant.ID, data = adults.lad[adults.lad$days_to_eclosion <35,], FUN = mean)

agg.eclosion$Mon.Pop <- factor(agg.eclosion$Mon.Pop, c('ENA','CA','HI','GU','AU','PR'))

Fig4c1 <- ggplot(agg.eclosion, aes( x = Mon.Pop, y = days_to_eclosion, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Days to Eclosion')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  facet_wrap(~Species)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='gold'))+
  scale_y_reverse()+
  annotate('rect', ymin = 15, ymax = 35, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 15, ymax = 35, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)

Fig4c1

ggsave(file="./Figures/Fig4c1.tiff", plot=Fig4c1, width=6, height=7)

#now make a figure of model output for days to eclosion

plot.agg.ecl.sym <- as.data.frame(emmeans(eclosion.model1, 'sym.allo'))

Fig4c2 <- ggplot(plot.agg.ecl.sym, aes( x = sym.allo, y = emmean, col = sym.allo))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  ylab('Days to Eclosion')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())+
  scale_y_reverse(labels = scales::number_format(accuracy = 0.2))

Fig4c2

ggsave(file="./Figures/Fig4c2.tiff", plot=Fig4c2, width=3, height=7)

#

adults.lad$ancestral <- ifelse(adults.lad$Mon.Pop %in% c('ENA','CA'), 'ancestral', 'derived')

eclosion.model.x1 <- lmer(days_to_eclosion ~ Species + ancestral + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad[adults.lad$days_to_eclosion < 35,])

Anova(eclosion.model.x1)

plot.ecl.anc <- as.data.frame(emmeans(eclosion.model.x1, 'ancestral'))

Fig4c3 <- ggplot(plot.ecl.anc, aes( x = ancestral, y = emmean, col = ancestral))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  scale_color_manual(values = c('blue','forestgreen'), guide = F)+
  theme_bw()+
  ylab('Days to Eclosion')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())+
  scale_y_reverse()

Fig4c3

ggsave(file="./Figures/Fig4c3.tiff", plot=Fig4c3, width=3, height=7)

###

#finally create figure for mass at eclosion; first build model

adult.mass.model1 <- lmer(emergence_weight ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad) #gives warning about singularity; however, simplifying model by dropping non-significant terms does not meaningfully change estimates of other parameters, so stick with this model for consistency; isSingular also returns FALSE, so should be alright

summary(adult.mass.model1)
Anova(adult.mass.model1)

em.pop <- as.data.frame(emmeans(adult.mass.model1, 'Mon.Pop'))
em.spp <- as.data.frame(emmeans(adult.mass.model1, 'Species'))
em.sym <- as.data.frame(emmeans(adult.mass.model1, 'sym.allo'))
em.use <- as.data.frame(emmeans(adult.mass.model1, 'Usage'))
em.gh <- as.data.frame(emmeans(adult.mass.model1, 'GH'))
em.yr <- as.data.frame(emmeans(adult.mass.model1, 'Year'))
em.sex <- as.data.frame(emmeans(adult.mass.model1, 'Sex'))


#now aggregate raw data for figure

adult.mass.agg <- aggregate(emergence_weight ~ Species + Mon.Pop + sym.allo + maternal_family + Plant.ID, mean, data = adults.lad)

adult.mass.agg$Mon.Pop <- factor(adult.mass.agg$Mon.Pop, c('ENA','CA','HI','GU','AU','PR'))

Fig4d1 <- ggplot(adult.mass.agg, aes( x = Mon.Pop, y = emergence_weight, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Adult Eclosion Mass (g)')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  facet_wrap(~Species)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='gold'))+
  annotate('rect', ymin = 0.1, ymax = 1, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 0.1, ymax = 1, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)

Fig4d1

ggsave(file="./Figures/Fig4d1.tiff", plot=Fig4d1, width=6, height=7)

#now highlight sym/allo contrast; first put fitted values with original dataframe

plot.agg.ecl.sym <- as.data.frame(emmeans(adult.mass.model1, 'sym.allo'))

Fig4d2 <- ggplot(plot.agg.ecl.sym, aes( x = sym.allo, y = emmean, col = sym.allo))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  ylab('Adult Eclosion Mass (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())

Fig4d2

ggsave(file="./Figures/Fig4d2.tiff", plot=Fig4d2, width=3, height=7)

adult.mass.model.x1 <- lmer(emergence_weight ~ Species + ancestral + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad)

Anova(adult.mass.model.x1)

plot.eclm.anc <- as.data.frame(emmeans(adult.mass.model.x1, 'ancestral'))

Fig4d3 <- ggplot(plot.eclm.anc, aes( x = ancestral, y = emmean, col = ancestral))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  scale_color_manual(values = c('blue','forestgreen'), guide = F)+
  theme_bw()+
  ylab('Adult Eclosion Mass (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())

Fig4d3

ggsave(file="./Figures/Fig4d3.tiff", plot=Fig4d3, width=3, height=7)


###### Don't include as a primary figure, but quickly check on wing area as a response variable

adults.lad$MArea <- rowMeans(cbind(adults.lad$LArea, adults.lad$RArea), na.rm = T)
adults.lad$Aspect_Ratio <- rowMeans(cbind(adults.lad$LLength,adults.lad$RLength),na.rm=T) / rowMeans(cbind(adults.lad$RWidth,adults.lad$LWidth),na.rm = T)

adult.area.model1 <- lmer(MArea ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + Year + scale(exp.days), data = adults.lad)

summary(adult.area.model1)
Anova(adult.area.model1) #Wing area is explained almost entirely by population of origin, with minimal effects based on host plant. There is also no signficant sympatric/allopatric contrast, which makes sense given that eclosion mass also did not have this contrast and wing area vs. body mass are almost perfect correlated.

wing.spp <- as.data.frame(emmeans(adult.area.model1, 'Species'))
wing.pop <- as.data.frame(emmeans(adult.area.model1, 'Mon.Pop'))
wing.sym <- as.data.frame(emmeans(adult.area.model1, 'sym.allo'))
wing.use <- as.data.frame(emmeans(adult.area.model1, 'Usage'))
wing.gh <- as.data.frame(emmeans(adult.area.model1, 'GH'))
wing.yr <- as.data.frame(emmeans(adult.area.model1, 'Year'))
wing.sex <- as.data.frame(emmeans(adult.area.model1, 'Sex'))

adult.aspectratio.model1 <- lmer(Aspect_Ratio ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + OE + Sex, data = adults.lad)

summary(adult.aspectratio.model1) #very unexpectedly, monarchs tend to have less elongated wings when they are reared on their sympatric host plants; however, the effect is still modest relative to inherent differences among populations. By contrast, host plant doesn't seem to affect wing shape very much.
Anova(adult.aspectratio.model1) 

summary(glht(adult.area.model1, mcp(Mon.Pop="Tukey"))) #PR smaller than everything besides AU; ENA larger than everything besides GU and CA
summary(glht(adult.area.model1, mcp(Species="Tukey"))) #no significant pairwise differences among species

summary(glht(adult.aspectratio.model1, mcp(Mon.Pop="Tukey"))) #PR rounder than everything; HI more elongated than Guam
summary(glht(adult.aspectratio.model1, mcp(Species="Tukey"))) #no significant pairwise differences among milkweed species





###### Now generate figures for coefficient of variation in performance

cv.model <- aggregate(Weight ~ Species + Mon.Pop + sym.allo + Usage + Mon.ID + Pop + Plant.ID, FUN = cv, data = lad)

names(cv.model)[8] <- 'CV'

cv.weight.model1 <- lmer(CV ~ Species + Mon.Pop + sym.allo + Usage + (1|Mon.ID) + Usage + (1|Pop/Plant.ID), data = cv.model)

Anova(cv.weight.model1) #when only non-zero caterpillar weights are included, sym/allo contrast is not significant

#####

lad.zeros <- lad #create new DF to use for second test of CV
lad.zeros[is.na(lad.zeros$Weight),]$Weight <- 0 #cannot use this in aggregate with CV as the function, since it cannot be calculated for entries with mean of 0

zer.agg.mass <- aggregate(Weight ~ Species + Mon.Pop + Plant.ID + sym.allo + GH + Year + Usage + Mon.ID + Larvae + Pop, mean, data = lad.zeros)

repl.nas <- which(aggregate(Weight ~ Species + Mon.Pop + Plant.ID + sym.allo + GH + Year + Usage + Mon.ID + Larvae + Pop, mean, data = lad.zeros)$Weight == 0)

zer.cv.model <- lad.zeros[!lad.zeros$Plant.ID %in% zer.agg.mass[repl.nas,]$Plant.ID,]


zer.cv.model.agg <- aggregate(Weight ~ Species + Mon.Pop + Plant.ID + sym.allo + GH + Year + Usage + Mon.ID + Larvae + Pop, FUN = cv, data = zer.cv.model)

names(zer.cv.model.agg)[11] <- 'CV'

cv.weight.model2 <- lmer(CV ~ Species + Mon.Pop + sym.allo + (1|Mon.ID) + GH + Year + Usage + Larvae + (1|Pop/Plant.ID), data = zer.cv.model.agg)

summary(cv.weight.model2)

Anova(cv.weight.model2) #when zeros are included, returns a significant sym/allo contrast for caterpillars reared on sympatric versus allopatric host plants

emmeans(cv.weight.model2, 'Species')
emmeans(cv.weight.model2, 'Mon.Pop')
emmeans(cv.weight.model2, 'sym.allo')


MuMIn::r.squaredGLMM(cv.weight.model2) #fixed effects explain 7.6% of variation, random effects explain 31%

hist(resid(cv.weight.model2),breaks= 100)

summary(glht(cv.weight.model2, mcp(Mon.Pop="Tukey"))) #no difference in CV between monarch populations

summary(glht(cv.weight.model2, mcp(Species="Tukey"))) #ASYR has the most variable performance, ASPEC also quite variable; AINC least variable

### plot coefficients of variation across milkweed species and monarch populations; first as a boxplot

ggplot(zer.cv.model.agg, aes( x = Mon.Pop, y = CV, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Coefficient of Variation (Larval Mass)')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  facet_wrap(~Species)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='gold'))

plot.cv.species <- as.data.frame(emmeans(cv.weight.model2, 'Species'))
plot.cv.species$Species <- factor(plot.cv.species$Species, c('ASPEC','ASYR','ASFA','AINC','ASCU','GOPH'))

fig5a <- ggplot(plot.cv.species, aes( x = Species, y = emmean))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  theme_bw()+
  ylab('Coefficient of Variation (Larval Mass)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16))+
  annotate('text', x = 1, y = 1.2, label = 'BC')+
  annotate('text', x = 2, y = 1.15, label = 'C')+
  annotate('text', x = 3, y = 1.1, label = 'ABC')+
  annotate('text', x = 4, y = 0.85, label = 'A')+
  annotate('text', x = 5, y = 0.97, label = 'AB')+
  annotate('text', x = 6, y = 1, label = 'ABC')

fig5a

ggsave(file="./Figures/Fig5a.tiff", plot=fig5a, width=6, height=4.5)

plot.cv.monpop <- as.data.frame(emmeans(cv.weight.model2, 'Mon.Pop'))
plot.cv.monpop$Mon.Pop <- factor(plot.cv.monpop$Mon.Pop, c('ENA','CA','HI','GU','AU','PR'))

fig5b <- ggplot(plot.cv.monpop, aes( x = Mon.Pop, y = emmean))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  theme_bw()+
  ylab('Coefficient of Variation (Larval Mass)')+
  xlab('Monarch Population')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16))+
  annotate('text', x = 1, y = 1.08, label = 'A')+
  annotate('text', x = 2, y = 1.08, label = 'A')+
  annotate('text', x = 3, y = 1.2, label = 'A')+
  annotate('text', x = 4, y = 0.96, label = 'A')+
  annotate('text', x = 5, y = 0.96, label = 'A')+
  annotate('text', x = 6, y = 1.12, label = 'A')

fig5b

ggsave(file="./Figures/Fig5b.tiff", plot=fig5b, width=6, height=4.5)

plot.cv.sym.allo <- as.data.frame(emmeans(cv.weight.model2, 'sym.allo'))

fig5c <- ggplot(plot.cv.sym.allo, aes( x = sym.allo, y = emmean, col = sym.allo))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  ylab('Coefficient of Variation (Larval Mass)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
        axis.title.x = element_blank())

fig5c

ggsave(file="./Figures/Fig5c.tiff", plot=fig5c, width=3, height=9)



##### try cv analysis, but this time without aggregating by plant ID

zer.cv.model.agg2 <- aggregate(Weight ~ Species + Mon.Pop + sym.allo + GH + Usage + Mon.ID + Pop, FUN = cv, data = zer.cv.model) #gives 688 replicates of maternal family * host species
names(zer.cv.model.agg2)[8] <- 'CV'

cv.weight.model3 <- lmer(CV ~ Species + Mon.Pop + sym.allo + (1|Mon.ID) + GH  + Usage + (1|Pop), data = zer.cv.model.agg2)

summary(cv.weight.model3)
Anova(cv.weight.model3) #gives generally similar results to the previous model that used plant genotype x maternal family as the unit of replication
r.squaredGLMM(cv.weight.model3)

#finally try very coarse analysis that includes only species, population, and sym/allo

simp <- aggregate(Weight ~ Species + Mon.Pop + sym.allo, cv, data = lad)

Anova(lm(Weight ~ Species + Mon.Pop + sym.allo, simp)) #washes out a lot of the variation; still maintains finding that variation on AINC is lower than all other species, but sympatric/allopatric contrast disappears

##### secondary analysis: what happens when North American populations and species are lumped together into one big sympatric complex?

lad$sym.allo2 <- ifelse(lad$Mon.Pop %in% c('ENA','CA') & lad$Species %in% c('ASPEC','ASYR','AINC','ASFA'), 'sympatric', ifelse(lad$Mon.Pop %in% c('HI','AU') & lad$Species == 'GOPH', 'sympatric', ifelse(lad$Mon.Pop %in% c('GU','PR') & lad$Species == 'ASCU', 'sympatric', 'allopatric')))

#now redo same analysis as before, but this time with new sym.allo term

model_x.cat.growth <- lmer((log(Weight)) ~ Species + Mon.Pop + sym.allo2 + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad)

summary(model_x.cat.growth)

Anova(model_x.cat.growth) #still recapitulates local adaptation result, although magnitude of sympatric/allopatric contrast is weaker

### alternative approach: analyze data with ENA and CA only


model_y.cat.growth <- lmer((log(Weight)) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad[lad$Mon.Pop %in% c('ENA','CA'),])

summary(model_y.cat.growth)

Anova(model_y.cat.growth) #no significant within-North America effect of local adaptation, although the effect is in the same direction that would be expected, with sym>allo

### correlation between adult mass and wingspan

summary(lm(MArea ~ emergence_weight, data = adults.lad))

ggplot(adults.lad, aes ( x = emergence_weight, y = MArea))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = 'lm', col = 'red')+
  ylab(bquote('Forewing area ('*cm^2*')'))+
  xlab('Mass at emergence (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

summary(lm(MArea ~ emergence_weight, data = adults.lad[adults.lad$Mon.Pop=='PR',]))
       
ggplot(na.omit(adults.lad), aes ( x = emergence_weight, y = MArea))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = 'lm', col = 'red')+
  ylab(bquote('Forewing area ('*cm^2*')'))+
  xlab('Mass at emergence (g)')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))+
  facet_wrap(~Mon.Pop)+
  theme(strip.text = element_text(size = 16))


##### predictions slides

library(ggsignif)

pred2 <- na.omit(read.csv(file = '~/Documents/grad_school/Manuscripts/Monarch local adaptation/temp files/pred2.csv'))

pred2$Population <- as.factor(pred2$Population)
pred2$Host <- as.factor(pred2$Host)
pred2$Population
pred2$sym.allo <- ifelse(pred2$Population == 1 & pred2$Host == 1, 'sympatric',
                         ifelse(pred2$Population==2 & pred2$Host == 2,'sympatric',
                                ifelse(pred2$Population == 3 & pred2$Host ==3, 'sympatric', 'allopatric')))

ggplot(pred2, aes(x = Population, y = Performance, fill = Host, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0, size = 0.8, position = position_dodge(0.85))+
  theme_bw()+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))+
  scale_fill_manual(values = c('gold','lightgreen','cornflowerblue','thistle1'),
                    name = 'Host Plant')+
  scale_color_manual(values = c('black','red'), name = '')+
  geom_vline(xintercept = c(1.5,2.5), lty = 2)+
  ylim(c(0,8.5))+
  xlab('Herbivore Population')+
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))+
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 0.5, ymax = 1.7, fill = 'grey', col = 'black')+
  annotate("text")
  

ggplot(pred2, aes(x = Population, y = Performance, fill = Host, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0, size = 0.8, position = position_dodge(0.85))+
  theme_bw()+
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))+
  scale_fill_manual(values = c('gold','lightgreen','cornflowerblue','purple'),
                    name = 'Host Plant')+
  scale_color_manual(values = c('black','red'), name = '')+
  geom_vline(xintercept = c(1.5,2.5), lty = 2)+
  ylim(c(0,8.5))+
  xlab('Herbivore Population')+
  geom_signif(y_position=c(7.5, 6.8, 5, 8.2), 
              xmin = c(0.55,1.55,2.55, 0.55), xmax = c(1.45,2.45,3.45, 3.45),
              annotation = c('Pop 1 mean = 5.6', 
                             'Pop 2 mean = 3.5', 
                             'Pop 3 mean = 2.8',
                             'Prediction 2'), 
              color = 'black', tip_length = 0.03)+
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20))+
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 0.2, ymax = 1.4, fill = 'grey', col = 'black')+
  annotate('text', x = 1, y = 1.2, label = 'Prediction 3', size = 4)+
  annotate('text', x = 1, y = 0.8, label = 'allopatric CV = 0.44', size = 3)+
  annotate('text', x = 1, y = 0.5, label = 'sympatric CV = 0.19', color = 'red', size = 3)+
  annotate('rect', xmin = 0.55, xmax = 1.45, ymin = 1.8, ymax = 3, fill = 'grey', col = 'black')+
  annotate('text', x = 1, y = 2.8, label = 'Prediction 1', size = 4)+
  annotate('text', x = 1, y = 2.4, label = 'allopatric mean = 3.64', size = 3)+
  annotate('text', x = 1, y = 2.1, label = 'sympatric mean = 5.06', color = 'red', size = 3)


aggregate(Performance ~ Population, FUN = mean, pred2)

aggregate(Performance ~ sym.allo + Population, FUN = cv, pred2)
aggregate(Performance ~ Population, FUN = cv, pred2)
aggregate(Performance ~ sym.allo, FUN = cv, pred2)


########

#basic crossing reaction norms plot

basic.la <- read.csv('~/Downloads/demos.csv')

ggplot(basic.la[basic.la$Scenario==1,], 
       aes(x = Host.Plant, y = Performance, group = Population, col = Population))+
  geom_point()+
  theme_bw()+
  geom_line()+
  theme(strip.text = element_blank())+
  scale_color_manual(values = c('blue','forestgreen'))+
  ylim(c(0,1))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16))+
  xlab('Host Plant')

ggplot(basic.la, aes(x = Host.Plant, y = Performance, group = Population, col = Population))+
  geom_point()+
  theme_bw()+
  geom_line()+
  facet_wrap(~Scenario)+
  theme(strip.text = element_blank())+
  scale_color_manual(values = c('blue','forestgreen'))+
  ylim(c(0,1))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16))+
  xlab("Host Plant")

ggplot(basic.la[basic.la$Scenario %in% c(1,3),], 
       aes(x = Host.Plant, y = Performance, group = Population, col = Population))+
  geom_point()+
  theme_bw()+
  geom_line()+
  facet_wrap(~Scenario, nrow = 1)+
  theme(strip.text = element_blank())+
  scale_color_manual(values = c('blue','forestgreen'))+
  ylim(c(0,1))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16))+
  xlab("Host Plant")




sharon.test <- lmer(emergence_weight ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad[adults.lad$Mon.Pop %in% c('AU','HI','GU','PR'),])

  


###### other things to check -- see if equalizing the number of maternal families changes anything about the variance analysis

table(zer.cv.model.agg$Mon.ID, zer.cv.model.agg$Mon.Pop) #minimum number of families is 12 for ENA and HI

coefs.cv.model <- vector(length = 1000)
t.cv.model <- vector(length = 100)

for(i in 1:1000){

CA<-zer.cv.model.agg[zer.cv.model.agg$Mon.Pop=='CA',]
CA.keep <- sample(levels(factor(CA$Mon.ID)), size = 11, replace = F)

ENA<-zer.cv.model.agg[zer.cv.model.agg$Mon.Pop=='ENA',]
ENA.keep <- sample(levels(factor(ENA$Mon.ID)), size = 11, replace = F)

HI<-zer.cv.model.agg[zer.cv.model.agg$Mon.Pop=='HI',]
HI.keep <- sample(levels(factor(HI$Mon.ID)), size = 11, replace = F)

AU<-zer.cv.model.agg[zer.cv.model.agg$Mon.Pop=='AU',]
AU.keep <- sample(levels(factor(AU$Mon.ID)), size = 11, replace = F)

GU<-zer.cv.model.agg[zer.cv.model.agg$Mon.Pop=='GU',]
GU.keep <- sample(levels(factor(GU$Mon.ID)), size = 11, replace = F)

PR<-zer.cv.model.agg[zer.cv.model.agg$Mon.Pop=='PR',]
PR.keep <- sample(levels(factor(PR$Mon.ID)), size = 11, replace = F)

cv.rand.samp <- zer.cv.model.agg[zer.cv.model.agg$Mon.ID %in% c(CA.keep,ENA.keep,GU.keep,AU.keep,HI.keep,PR.keep),]

#now redo analysis with random sample

cv.weight.model.x <- lmer(CV ~ Species + Mon.Pop + sym.allo + (1|Mon.ID) + GH + Year + Usage + Larvae + (1|Pop/Plant.ID), data = cv.rand.samp)

coefs.cv.model[i] <- summary(cv.weight.model.x)$coefficients[12]
#coefs.cv.model[i] <- summary(cv.weight.model.x)$coefficients[12,4]

}

mean(coefs.cv.model)
mean(t.cv.model)

pt(mean(coefs.cv.model), 598.3, lower.tail = TRUE) #average t statistic across iterations is -1.65; corresponds to a p value of 0.049 assuming a two-tailed distribution

hist(coefs.cv.model, breaks = 50, xlab = 'Estimated coefficient', main = NULL)
abline(v = -8.34, col = 'red', lw = 2)
abline(v = (-8.34 + 3.918), col = 'red', lw = 1, lty =2)
abline(v = (-8.34 - 3.918), col = 'red', lw = 1, lty =2)
abline(v = mean(coefs.cv.model), col = 'blue', lw = 2)
abline(v = mean(coefs.cv.model)+sd(coefs.cv.model)*1.96, col = 'blue', lw = 1, lty = 2)
abline(v = mean(coefs.cv.model)-sd(coefs.cv.model)*1.96, col = 'blue', lw = 1, lty = 2)




summary(cv.weight.model2) #actual value from full model is -8.34

#add in latex data here

#load in latex data

latex <- read.csv('~/Downloads/Monarch local adaptation 2018 - Latex.csv')

latex.merge <- merge(lad, latex, by = 'Plant.ID')

latex.merge$Latex.conc <- latex.merge$Latex.conc - min(latex.merge$Latex.conc, na.rm=T)

latex.merge$Dmass <- latex.merge$Dmass - min(latex.merge$Dmass, na.rm=T)

latex.merge$Species.x <- revalue(latex.merge$Species.x, c('GOFR' = 'GOPH'))


latex.merge.agg <- aggregate(Weight ~ reps + Species.x + Mon.Pop + Latex.conc, mean, data =latex.merge)

lamagg <- aggregate(Latex.conc ~ reps + Species.x, mean, data =latex.merge)
lamagg$Species.x <- factor(lamagg$Species.x, c('ASPEC','ASYR','ASFA','AINC','ASCU','GOPH'))

latex.merge.agg <- latex.merge.agg[!is.na(latex.merge.agg$Species.x),]
latex.merge.agg$Species.x <- droplevels(latex.merge.agg$Species.x)

mcb <- c('darkolivegreen2','purple','orange','darkorange4','darkseagreen4','blue')
hdh <- c('darkorange4','darkseagreen4','orange','darkolivegreen2','purple','blue')

ggplot(lamagg, aes(x = Species.x, y = Latex.conc, fill = Species.x))+
  geom_boxplot(outlier.shape = NA, alpha = 0.1)+
  theme_bw()+
  geom_point(aes(color = Species.x), position = position_jitterdodge(1.5), pch = 20, size = 3)+
  scale_color_manual(values = hdh)+
  scale_fill_manual(values = hdh)+
  ylab('Latex (mg/g)')+
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab('Milkweed Species')

ggplot(na.omit(latex.merge[latex.merge$Latex.conc < 25,]), aes(x = (Latex.conc), y = Weight))+
  geom_point(aes(col = Species.x))+
  geom_smooth(method = 'lm', formula = y ~ poly(x, 2), col = 'black', lty = 4 )+
  scale_color_manual(values = mcb, name = 'Species')+
  theme_classic()+theme(legend.position = c(0.8,0.7),legend.text = element_text(size =14),
                        legend.title = element_blank())+
  xlab('Latex Production (mg/0.1g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16))+
  ylab('Caterpillar Mass (g)')


ggplot(latex.merge, aes(x = Latex.conc, y = Weight))+
  geom_point(aes(col = Species.x))+
  facet_wrap(~Species.x, scales = 'free')+
  geom_smooth(method = 'lm', col = 'black')+
  scale_color_manual(values = mcb)+
  theme_bw()+theme(legend.position = 'none')+
  xlab('Latex Concentration (mg/0.1g)')+
  ylab('Day 8 Caterpillar Mass (g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16))

ggplot(latex.merge, aes(x = Latex.conc, y = Weight))+
  geom_point(aes(col = Species.x))+
  facet_wrap(~Species.x, scales = 'free')+
  geom_smooth(method = 'lm', col = 'black', size = 0.5, alpha = 0.2)+
  scale_color_manual(values = mcb)+
  theme_bw()+theme(legend.position = 'none')+
  xlab('Latex (mg/g)')+
  ylab('Day 8 Caterpillar Mass (g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16))

ggplot(latex.merge, aes(x = Latex.conc, y = Weight))+
  geom_point(aes(col = Species.x))+
  facet_wrap(~Species.x*Mon.Pop, scales = 'free')+
  geom_smooth(method = 'lm') #takes about 20 seconds to render

ggplot(latex.merge, aes(x = Latex.conc, y = Weight))+
  geom_point(aes(col = Species2))+
  facet_wrap(~Species2 + Mon.Pop, scales = 'free')+
  geom_smooth(method = 'lm')

yamsnax <- aggregate(Latex.conc ~ Species.x + Plant.ID, mean, data = latex.merge)

ggplot(yamsnax, aes(x = Species.x, y = Latex.conc, fill = Species.x))+
  geom_boxplot(outlier.shape = NA, alpha = 0.2)+theme_bw()+
  geom_point(position = position_jitter(width = 0.2), aes(col = Species.x), size = 0.8)+
  scale_color_manual(values = mcb)+
  scale_fill_manual(values = mcb)+
  theme(legend.position = 'none')+
  xlab('Species')+
  ylab('Standardized latex production (mg latex / g plant tissue)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16), axis.text.x = element_text(angle = 60, hjust = 1))

scarm <- aggregate(Dmass ~  Plant.ID + Species.x, mean, data = latex.merge)

ggplot(scarm, aes(x = Species.x, y = Dmass, fill = Species.x))+
  geom_boxplot(outlier.shape = NA, alpha = 0.2)+theme_bw()+
  geom_point(position = position_jitter(width = 0.2), aes(col = Species.x), size = 0.8)+
  scale_color_manual(values = mcb)+
  scale_fill_manual(values = mcb)+
  theme(legend.position = 'none')+
  ylab('Absolute latex production per leaf (mg)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16), axis.text.x = element_text(angle = 60, hjust = 1), axis.title.x = element_blank())

model.latex <- lmer(Weight ~ Species.x + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + Dmass + GH, data = latex.merge)

summary(model.latex)

#add in cardenolide data here

cardenolides <- read.csv('~/Documents/grad_school/Manuscripts/Monarch sequestration/concatenated_cardenolides.csv')

cardenolides$Year <- as.factor(cardenolides$Year)

masses <- read.csv('~/Documents/grad_school/Manuscripts/Monarch sequestration/masses.csv')

head(cardenolides)

cardenolides.standardized <- cbind(cardenolides[,1:7], cardenolides[,8:ncol(cardenolides)] / cardenolides$Digitoxin)

standards <- cardenolides.standardized[cardenolides.standardized$Tissue=='Standard',]
standard.means <- colSums(standards[,8:ncol(cardenolides)]) / nrow(standards) #generate average values for standards, to be subtracted from the rest of the dataframe

noop <- sweep(cardenolides.standardized[,8:ncol(cardenolides.standardized)], 2, standard.means)

noop[noop<0] = 0

cards <- cbind(cardenolides[,1:7],noop)

nrow(cards) #424 observations

nrow(cards[cards$Tissue=='Wing',]) #241 observations for wings

colnames(masses)[1] <- 'Sample'

ham <- merge(cards, masses, all.x = TRUE) #merge, but include rows where mass is NA

masses$Species <- revalue(masses$Species, c('GOFR' = 'GOPH')) #need to revalue so that GOFR is now GOPH

(average.masses <- aggregate(Mass ~ Species + Tissue, masses, mean)) #use these only for leaf tissues; for wings, use predicted regression values from correlation between hind wing area and mass

ham$spec.tiss <- paste(ham$Species, ham$Tissue, sep = "_") #create matching column
average.masses$spec.tiss <- paste(average.masses$Species, average.masses$Tissue, sep = "_") #matching column

replacements <- merge(ham[is.na(ham$Mass),], average.masses, all.x = TRUE) 

#now sub in actual values for replacements

replacements[replacements$spec.tiss==average.masses$spec.tiss[2],]$Mass <- average.masses$Mass[2]
replacements[replacements$spec.tiss==average.masses$spec.tiss[3],]$Mass <- average.masses$Mass[3]
replacements[replacements$spec.tiss==average.masses$spec.tiss[5],]$Mass <- average.masses$Mass[5]
replacements[replacements$spec.tiss==average.masses$spec.tiss[6],]$Mass <- average.masses$Mass[6]
replacements[replacements$spec.tiss==average.masses$spec.tiss[7],]$Mass <- average.masses$Mass[7]
replacements[replacements$spec.tiss==average.masses$spec.tiss[8],]$Mass <- average.masses$Mass[8]
replacements[replacements$spec.tiss==average.masses$spec.tiss[9],]$Mass <- average.masses$Mass[9]
replacements[replacements$spec.tiss==average.masses$spec.tiss[11],]$Mass <- average.masses$Mass[11]
replacements[replacements$spec.tiss==average.masses$spec.tiss[12],]$Mass <- average.masses$Mass[12]


ham <- ham[!is.na(ham$Mass),]

ham <- rbind(ham, replacements)

ham$Digitoxin <- NULL
ham$spec.tiss <- NULL

ham <- ham[,c(1:5,78,8:77)]

names(ham)

ham$cumulative <- rowSums(ham[,7:ncol(ham)])

ham$concentration <- (0.15 * ham$cumulative) / (ham$Mass/1000 * 0.8)

aggregate(concentration ~ Tissue + Species, mean, data = ham) #sequestration ratio for AINC is about 2.9, ASCU is about 1.6, ASFA is about 2.1, ASPEC is about 9.1, ASYR is about 15.9, GOPH is about 0.8


nrow(ham[ham$Tissue=='Leaf',]) #has 163 leaf tissue samples

leaves <- ham[ham$Tissue=='Leaf',]

colnames(leaves)[1] <- 'Plant.ID'

plot.card.perf <- merge(leaves[,c(1,78)], lad, by = 'Plant.ID')

scaleFUN <- function(x) sprintf("%.0f", x)

ggplot(plot.card.perf, aes(x = log(concentration*1000), y = Weight, col = Species))+
  geom_point()+
  facet_wrap(~Species, scales = 'free_x')+
  geom_smooth(method = 'lm')+
  scale_color_manual(values = mc2)+
  theme_bw()+theme(legend.position = 'none')+
  xlab('log(Cardenolides (ug/g))')+
  ylab('Day 8 Caterpillar Mass (g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16))+
  scale_y_continuous(labels=scaleFUN)

card.agg <- aggregate(Weight ~ Plant.ID + concentration + Species, plot.card.perf, mean)

ggplot(card.agg, aes(x = log(concentration*1000), y = Weight))+
  geom_point(aes(col = Species))+
  scale_color_manual(values = mc2)+
  theme_bw()+theme(legend.position = 'none')+
  xlab('log(Cardenolides (ug/g))')+
  ylab('Day 8 Caterpillar Mass (g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16))+
  geom_smooth(method = 'lm', col = 'black')

#general question: what happens if all non-surviving caterpillars are assigned a weight of 0?

lad2 <- lad
lad2$Weight[is.na(lad2$Weight)] <- 0

#now run same model as before

modelx4.cat.growth <- lmer((log(Weight+0.001)) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad2)

summary(modelx4.cat.growth) #again gives a significant sym/allo effect, although result is conditional upon log transforming the weight measure, which requires adding 0.001 so as to not take the log of 0 values; residuals are not even close to being normally distributed though, suggesting this isn't an appropriate model


#### is there LA within north america?

lad.americas <- subset(lad, lad$Mon.Pop %in% c('ENA','CA'))

americas.only <- lmer(log(Weight+0.001) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad.americas)

summary(americas.only) #modestly higher weight in sympatric combinations

amer.sym.allo.contrast <- as.data.frame(emmeans(americas.only, 'sym.allo'))

americas.only.survival <- glmer(survival ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad.americas, family = 'binomial')

summary(americas.only.survival) #survival not higher on sympatric hosts for North America

