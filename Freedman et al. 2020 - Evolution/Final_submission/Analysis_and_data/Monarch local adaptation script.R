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
library(lemon)

#load data

setwd('~/Documents/grad_school/Manuscripts/Monarch local adaptation/Final_submission/Analysis_and_data/')

lad <- read.csv("./larval_data.csv")

head(lad)
str(lad)

lad$Mon.Pop <- factor(lad$Mon.Pop, levels(lad$Mon.Pop)[c(4,5,2,3,1,6)]) #rearrange factor levels for monarch population

#convert year, usage, plant.id into factors

lad[,c('Year','Usage','Plant.ID')] <- lapply(lad[,c('Year','Usage','Plant.ID')], factor)

#change all values coded as GOFR to GOPH (previously Australian Gomphocarpus was treated as G. fruticosus, not G. physocarpus)

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

table(lad$Larvae) #vast majority of replicates come from plants that had 5 larvae added initially

############################################

dates <- as.Date(lad$Date, format='%m/%d/%Y') ## create vector of dates used in experiment

start_date = as.Date('2018-05-15', tz="UTC") #define earliest date on which caterpillars were added to plants (May 15)

lad$exp.days <- as.numeric(difftime(dates, start_date , units = "days")) ###create new variable (exp.days) that defines how many days after the start of the experiment a particular replicate was used. This term is also a stand-in for plant age and the photoperiod that larvae were experiencing during their development.

###

#create figure to show distribution of when caterpillars were reared across years

mon.pop.cols <- c('blue','cornflowerblue','purple','forestgreen','orange','hotpink')

lad$Mon.Pop <- factor(lad$Mon.Pop, c('AU','HI','GU','ENA','CA','PR'))

FigS1a <- ggplot(lad, aes(x = exp.days))+
  geom_histogram(bins = 40, fill = 'black')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))+
  xlab('Days Since May 15th')+
  ylab('Number of Larvae Measured')

FigS1a ### the majority of rearing occurred in the period 60-150 days after the onset of the experiment (May 15th), which corresponds to July-September

#ggsave(file="./Figures/Figure_S1a.tiff", plot=FigS1a, width=7, height=4)

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

#ggsave(file="./Figures/Figure_S1b.tiff", plot=FigS1b, width=7, height=4)

###############

lad$ancestral_pop <- ifelse(lad$Mon.Pop %in% c('ENA','CA'), 'ancestral', 'derived') ##define whether monarch populations are ancestral or derived

###############

## Begin statistical analyses

###############

#start with day eight larval mass

hist(lad$Weight, breaks=100) #long right hand tail, with bimodal peaks that probably correspond to discrete larval instars

hist(log(lad$Weight), breaks = 100) #log transformation produces a distribution that is less skewed, though still not totally normally distributed

### Build linear model for subsequent analyses. Previous analyses using more complicated model structure with two- and three-way interactions between fixed effects seemed to be overfit, so use simpler model that includes only the following terms:

#1. Fixed effect of milkweed species ('Species')
#2. Fixed effect of monarch population ('Mon.Pop')
#3. Fixed effect of sympatric/allopatric status ('sym.allo')
#4. Random effect of plant ID nested inside milkweed seed family (1|Pop/Plant.ID)
#5. Random effect of monarch maternal family (1|Mon.ID)
#6. Random effect of position within greenhouse (1|Group)
#7. Covariate for whether a plant was being used for the first or second time ('Usage')
#8. Covariate for which greenhouse caterpillars were reared in ('GH')
#9. Covariate for the time at which rearing took place ('exp.days'), described above; this term is centered and scaled since it is continuous
#10. Covariate for the year in which rearing took place (2017 versus 2018)

model4.cat.growth <- lmer((log(Weight)) ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad) 

summary(model4.cat.growth)

r.squaredGLMM(model4.cat.growth) #fixed effects alone explain about 30% of variance; adding random factors explains another 36% or so

hist(resid(model4.cat.growth), breaks = 100) # histogram of model residuals; looks good, though perhaps with a slightly long left hand tail

Anova(model4.cat.growth) #can use default type II ANOVA since there are no interaction terms
#This is the output reported in the main text for larval growth on day 8, and also in the supplementary materials

pop.comparison <- glht(model4.cat.growth, mcp(Mon.Pop="Tukey"))
summary(pop.comparison) ## tests for significant pairwise performance differences among monarch populations, using TukeyHSD correction for multiple comparisons

summary(glht(model4.cat.growth, mcp(Species="Tukey"))) ##does the same for milkweed species

## get estimated marginal means for all factor levels

(day8.spp <- as.data.frame(emmeans(model4.cat.growth, "Species"))) #note: need to exponentiate to transform back to raw values
day8.mon <- as.data.frame(emmeans(model4.cat.growth, "Mon.Pop"))
sym.allo.contrast <- as.data.frame(emmeans(model4.cat.growth, 'sym.allo'))
day8.gh <- as.data.frame(emmeans(model4.cat.growth, "GH"))
day8.yr <- as.data.frame(emmeans(model4.cat.growth, "Year"))
day8.use <- as.data.frame(emmeans(model4.cat.growth, "Usage"))
day8.days <- as.data.frame(emmeans(model4.cat.growth, "exp.days"))

### Prepare figure showing day 8 larval mass across host plant x monarch population combinations -- use raw data instead of model output

agg.pop <- aggregate(Weight ~ Mon.ID + Species + Mon.Pop + sym.allo, data = lad, FUN = mean) #aggregate raw data so that each maternal family x host plant is a unique data point

mon.order <- c('CA','ENA','HI','GU','AU','PR') #set up vector of monarch populations in order, to use after operations that re-arrange factor levels alphabetically
spp.order <- c('ASPEC','ASYR','ASFA','AINC','ASCU','GOPH') #same thing for milkweed species

agg.pop$Mon.Pop <- factor(agg.pop$Mon.Pop, mon.order)
agg.pop$Species <- factor(agg.pop$Species, spp.order)

(Fig4a1 <- ggplot(agg.pop[agg.pop$Weight < 1.4,], aes( x = Mon.Pop, y = Weight, col = sym.allo))+
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
    theme(strip.background =element_rect(fill='cornsilk'))+
    annotate('rect', ymin = 0, ymax = 1.25, xmin = 0.5, xmax = 2.4, 
             col = 'blue', fill = NA, size = 0.7)+
    annotate('rect', ymin = 0, ymax = 1.25, xmin = 2.5, xmax = 6.5, 
             col = 'forestgreen', fill = NA, size = 0.7))

Fig4a1 <- Fig4a1 + facet_rep_wrap(~Species, repeat.tick.labels = 'bottom') #replot data with factor levels below each panel

Fig4a1

### Now create figure with marginal means of performance across milkweed species

day8.spp$Species <- factor(day8.spp$Species, spp.order) #rearrange factor levels

summary(glht(model4.cat.growth, linfct = mcp(Species = "Tukey"))) #check for significant pairwise differences among species

mcb <- c('darkorange4','darkseagreen4','orange','darkolivegreen2','purple','blue') #establish vector of milkweed colors, to use throughout 

(sup.fig.emm.cat.growth <- ggplot(day8.spp, aes( x = Species, y = exp(emmean), col = Species))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Day 8 Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Milkweed Species')+
    ylim(c(0,0.5))+
    scale_color_manual(values = mcb)+
    annotate('text', x = 1, y = 0.18, label = 'BC', size = 6)+
    annotate('text', x = 2, y = 0.23, label = 'B', size = 6)+
    annotate('text', x = 3, y = 0.35, label = 'AB', size = 6)+
    annotate('text', x = 4, y = 0.45, label = 'A', size = 6)+
    annotate('text', x = 5, y = 0.48, label = 'A', size = 6)+
    annotate('text', x = 6, y = 0.35, label = 'AB', size = 6))
##This figure shows up in one of the supplements

### Now do the same for marginal means of monarch populations

day8.mon$Mon.Pop <- factor(day8.mon$Mon.Pop, mon.order)

(sup.fig.emm.cat.growth2 <- ggplot(day8.mon, aes( x = Mon.Pop, y = exp(emmean)))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Day 8 Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1))+
    xlab('Monarch Population')+
    annotate('text', x = 1, y = 0.4, label = 'A', size = 6)+
    annotate('text', x = 2, y = 0.47, label = 'A', size = 6)+
    annotate('text', x = 3, y = 0.24, label = 'BC', size = 6)+
    annotate('text', x = 4, y = 0.38, label = 'AB', size = 6)+
    annotate('text', x = 5, y = 0.41, label = 'AB', size = 6)+
    annotate('text', x = 6, y = 0.22, label = 'C', size = 6))
##This figure is also in one of the supplements; letters are based on Tukey test comparisons 

################################

#now plot marginal means summarizing sympatric/allopatric contrast, with error bars correspond to 95% confidence intervals

(Fig4a2 <- ggplot(sym.allo.contrast, aes( x = sym.allo, y = exp(emmean), col = sym.allo))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
    scale_color_manual(values = c('black','red'), guide = F)+
    theme_bw()+
    ylab('Day 8 Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    annotate('text', x = 1, y = 0.28, label = 't = 3.968\np < 0.001', size = 5)+
    annotate('text', x = 1.5, y = 0.25, label = '***', size = 8))

Fig4a2

#ggsave(file="./Figures/Fig4a2.tiff", plot=Fig4a2, width=3, height=7)

##################################

#to generate last panel in figure, need to run separate model that replaces mon.pop term with ancestral/derived population status

model.ancestral1 <- lmer((log(Weight)) ~ Species + ancestral_pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad)

summary(model.ancestral1) #returns very similar output to previous model, and includes strong ancestral / derived contrast term

anc.em <- as.data.frame(emmeans(model.ancestral1, 'ancestral_pop'))

(Fig4a3 <- ggplot(anc.em, aes( x = ancestral_pop, y = exp(emmean), col = ancestral_pop))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = exp(upper.CL), ymin = exp(lower.CL)), width = 0.2, size=1)+
    scale_color_manual(values = c('blue','forestgreen'), guide = F)+
    theme_bw()+
    ylab('Day 8 Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    annotate('text', x = 2, y = 0.35, label = 't = 3.233\np = 0.002', size = 5)+
    annotate('text', x = 1.5, y = 0.25, label = '**', size = 8))

#ggsave(file="./Figures/Fig4a3.tiff", plot=Fig4a3, width=3, height=7)

##########################################

## based on reviewer comment, also run model that includes interaction between sym/allo and ancestral/derived status of host plants

model.cat.growth.rev2 <- lmer((log(Weight)) ~ Species + Mon.Pop + ancestral_pop*sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad) #this is the same model as used previous analyses, but now with an interaction term with sym.allo

summary(model.cat.growth.rev2) #no evidence for significant interaction between ancestral/derived population status and sympatric/allopatric status, as might be expected if derived populations are especially specialized on their host plants

## can also run alternative model to test for whether local adaptation effect is driven by especially reduced performance of derived populations on ancestral hosts

#first set up contrast for whether host plants are ancestral or derived

lad$anc_der_host <- ifelse(lad$Species %in% c('ASCU','GOPH'), 'derived_host', 'ancestral_host')

model.cat.growth.xx <- lmer((log(Weight)) ~ Species + Mon.Pop + ancestral_pop*anc_der_host + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad) #same as previous models, with with interaction between ancestral/derived status of both monarch populations and host plant species; here an interaction might indicate that derived populations struggle on ancestral hosts, or vice versa

summary(model.cat.growth.xx) ## get modest interaction term, but in the opposite direction to that predicted: derived populations have a DISADVANTAGE overall on derived hosts; part of this effect could be the terrible performance of Puerto Rico on GOPH, although some of that should also be captured in the sym.allo term

#########################################

#next run models for survival to day 8

#########################################

lad$survival = ifelse(is.na(lad$Weight) == TRUE, 0, 1) #add column for survival

#use same model structure as model 4

model1.survival <- glmer(survival ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad, family = 'binomial') #heads up: takes about two minutes to run

summary(model1.survival) #still get significant sym.allo term

Anova(model1.survival) #this is the output reported in the results section

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
} #define function to convert estimated marginal means back into survival percentages

## get marginal means for all survival terms

day8.surv.spp <- as.data.frame(emmeans(model1.survival, "Species")) ##note: all of these calculations give a warning message since original model had some difficult estimating variance-covariance matrix (although still produced sensible/interpretable output that seems to make sense); data tables returned here show Inf for the df term
day8.surv.pop <- as.data.frame(emmeans(model1.survival, "Mon.Pop"))
day8.surv.sym <- as.data.frame(emmeans(model1.survival, "sym.allo"))
day8.surv.use <- as.data.frame(emmeans(model1.survival, "Usage"))
day8.surv.gh <- as.data.frame(emmeans(model1.survival, "GH"))
day8.surv.yr <- as.data.frame(emmeans(model1.survival, "Year"))

day8.surv.spp$Species <- factor(day8.surv.spp$Species, spp.order)

summary(glht(model1.survival, linfct = mcp(Species = "Tukey")))

(sup.fig.emm.surv.spp <- ggplot(day8.surv.spp, aes( x = Species, y = logit2prob(emmean), 
                                                    col = Species))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = logit2prob(asymp.LCL), 
                      ymin = logit2prob(asymp.UCL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Proportion Survival to Day 8')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Milkweed Species')+
    ylim(c(0.55,0.95))+
    scale_color_manual(values = mcb)+
    annotate('text', x = 1, y = 0.8, label = 'C', size = 6)+
    annotate('text', x = 2, y = 0.82, label = 'BC', size = 6)+
    annotate('text', x = 3, y = 0.87, label = 'ABC', size = 6)+
    annotate('text', x = 4, y = 0.93, label = 'A', size = 6)+
    annotate('text', x = 5, y = 0.88, label = 'AB', size = 6)+
    annotate('text', x = 6, y = 0.85, label = 'ABC', size = 6))
## Supplemental figure showing survival differences across milkweed species

#ggsave(file="./Figures/sup.fig.emm.surv.spp.tiff", plot = sup.fig.emm.surv.spp, width=4, height=7)

summary(glht(model1.survival, linfct = mcp(Mon.Pop = "Tukey")))

day8.surv.pop$Mon.Pop <- factor(day8.surv.pop$Mon.Pop, mon.order)

(sup.fig.emm.surv.pop <- ggplot(day8.surv.pop, aes( x = Mon.Pop, y = logit2prob(emmean)))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = logit2prob(asymp.LCL), 
                      ymin = logit2prob(asymp.UCL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Proportion Survival to Day 8')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Monarch Population')+
    ylim(c(0.55,0.95))+
    annotate('text', x = 1, y = 0.84, label = 'A', size = 6)+
    annotate('text', x = 2, y = 0.84, label = 'A', size = 6)+
    annotate('text', x = 3, y = 0.8, label = 'A', size = 6)+
    annotate('text', x = 4, y = 0.92, label = 'A', size = 6)+
    annotate('text', x = 5, y = 0.92, label = 'A', size = 6)+
    annotate('text', x = 6, y = 0.91, label = 'A', size = 6))
### Supplemental figure showing that survival does not differ significantly in any pairwise comparisons across monarch populations

#ggsave(file="./Figures/sup.fig.emm.surv.pop.tiff", plot = sup.fig.emm.surv.pop, width=4, height=7)

### Create Fig4b1

agg.pop.surv <- aggregate(survival ~ Mon.ID + Species + Mon.Pop + sym.allo, data = lad, FUN = mean) #aggregate raw data so that each maternal family x host plant is a unique data point

agg.pop.surv$Mon.Pop <- factor(agg.pop.surv$Mon.Pop, mon.order)
agg.pop.surv$Species <- factor(agg.pop.surv$Species, spp.order)

Fig4b1 <- ggplot(agg.pop.surv, aes( x = Mon.Pop, y = survival, col = sym.allo))+
  geom_boxplot(outlier.alpha = 0)+
  scale_color_manual(values = c('black','red'), guide = F)+
  theme_bw()+
  xlab('Monarch Population')+
  ylab('Proportion Survival to Day 8')+
  geom_point(aes(fill = sym.allo), size = 2, alpha = 0.4, pch = 21, col = 'black',
             position = position_jitterdodge())+
  scale_fill_manual(values = c('black','red'), guide = F)+
  facet_wrap(~Species)+
  ylim(c(0,1.1))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size =16),
        strip.text = element_text(size = 16, face = 'italic'), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.background =element_rect(fill='cornsilk'))+
  annotate('rect', ymin = 0, ymax = 1.05, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 0, ymax = 1.05, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)
## As per suggestions of reviewers, show survival data in boxplot form for the sake of consistency. 

(Fig4b1 <- Fig4b1 + facet_rep_wrap(~Species, repeat.tick.labels = 'bottom'))

#ggsave(file="./Figures/Fig4b1.tiff", plot=Fig4b1, width=6, height=7)

### Make Fig4b2, showing contrast between allopatric and sympatric for survival

plot.sym.allo.surv <- as.data.frame(emmeans(model1.survival, "sym.allo"))

(Fig4b2 <- ggplot(plot.sym.allo.surv, aes( x = sym.allo, y = logit2prob(emmean), col = sym.allo))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = logit2prob(asymp.UCL), ymin = logit2prob(asymp.LCL)), width = 0.2, size=1)+
    scale_color_manual(values = c('black','red'), guide = F)+
    theme_bw()+
    ylab('Proportion Survival to Day 8')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    annotate('text', x = 1, y = 0.82, label = 't = 1.995\np = 0.046', size = 5)+
    annotate('text', x = 1.5, y = 0.77, label = '*', size = 8))

#ggsave(file="./Figures/Fig4b2.tiff", plot=Fig4b2, width=3, height=7)

#now generate ancestral/derived contrast

lad$ancestral <- ifelse(lad$Mon.Pop %in% c('ENA','CA'), 'ancestral', 'derived')

#### And finally, do the same for Fig4b3 -- first need to run a new model with ancestral/derived status of populations, as before

modelx1.survival <- glmer(survival ~ Species + ancestral_pop + sym.allo + (1|Pop/Plant.ID) + (1|Mon.ID)  + (1|Group)  + Usage + GH + scale(exp.days) + Year, data = lad, family = 'binomial') #heads up: takes about two minutes to run

summary(modelx1.survival) #no significant contrast here (not surprising given relatively weak monarch population effect in overall model)

Anova(modelx1.survival)

plot.anc.surv <- as.data.frame(emmeans(modelx1.survival, "ancestral_pop"))

(Fig4b3 <- ggplot(plot.anc.surv, aes( x = ancestral_pop, y = logit2prob(emmean), col = ancestral_pop))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = logit2prob(asymp.UCL), ymin = logit2prob(asymp.LCL)), width = 0.2, size=1)+
    scale_color_manual(values = c('blue','forestgreen'), guide = F)+
    theme_bw()+
    ylab('Proportion Survival to Day 8')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    annotate('text', x = 1, y = 0.81, label = 't = 1.189\np = 0.234', size = 5)) #if anything, derived populations have modestly higher survival across hosts, despite the fact that they seem to grow more slowly

#ggsave(file="./Figures/Fig4b3.tiff", plot=Fig4b3, width=3, height=7)

###########################################################

## Now load in adult data; this includes information on things like days to pupation, days to eclosion, OE infection status, sex of resulting butterflies, and all adult morphological information

###########################################################

#2017 data

adults.2017 <- read.csv(file = './2017_adults.csv')
str(adults.2017) #ID is indeed a factor
adults.2017$Group <- as.factor(adults.2017$Group)
#change ID into a variable called Cat.ID to merge with 2018 data
colnames(adults.2017)[1] <- 'Cat.ID'

adults.2018 <- read.csv(file = './Adults 2018.csv')
str(adults.2018)

#now load caterpillar ID info from other sheet

all.data <- read.csv(file='./experiment_master.csv')
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

full <- merge(all.data, adult.data, by = 'Cat.ID')

table(full$Cat.ID) ##looks correct

#drop values with no OE entry

adults.lad <- full[-which(is.na(full$OE)),] #leaves behind 1079 records

adults.lad$days_to_pupation <- as.numeric(as.Date(adults.lad$Pupation, format='%m/%d/%y')) - as.numeric(as.Date(adults.lad$Date, format='%m/%d/%y')) + 8
#create new variable for days to pupation: need to add eight because date variable refers to day on which larvae were measured, which was eight days after they were set up
adults.lad$days_to_eclosion <- as.numeric(as.Date(adults.lad$Eclosion, format='%m/%d/%y')) - as.numeric(as.Date(adults.lad$Date, format='%m/%d/%y')) + 8
#as with days to pupation, days to eclosion gets an 8 added to reflect the first eight days prior to measurement

hist(adults.lad$days_to_eclosion,breaks = 30) #normally distributed, but with a few extremely long outliers -- these were almost all butterflies that were severely diseased or malformed; anything beyond about 35 days can be safely discarded (some breeders apparently use 26-28 days as a hard cutoff)

head(adults.lad) #shows dataframe used for all subsequent analyses of performance in pupa and adults

ggplot(adults.lad, aes( x = days_to_pupation, y = days_to_eclosion))+
  geom_point()+
  geom_smooth(method = 'lm') #not surprisingly, get pretty much perfect correspondence between putation date and eclosion date, although there's a bit of spread (days as a pupa ranged from as little as 7 to as much as 10, probably reflecting some 2-day old pupa that were recorded as day 1)

#create column for monarch ID;  first drop entries with 'NA' for monarch population, since these were typically individuals found in the greenhouse without any identifying information or that weren't recorded when plants were set up initially

adults.lad <- adults.lad[-which(is.na(adults.lad$Mon.Pop)),]

adults.lad$maternal_family <- with(adults.lad, paste(Mon.Pop, ID, Year, sep = "_"))

length(unique(adults.lad$maternal_family)) #82 unique maternal families

table(adults.lad$maternal_family) #show numbers of adult butterflies that resulted from each individual maternal family

#convert all usage > 2 into 2 (this only applies to two plants that were used 3 times)

adults.lad$Usage <- ifelse(adults.lad$Usage > 2, 2, adults.lad$Usage)

# make usage and year factors

adults.lad$Year <- as.factor(adults.lad$Year)
adults.lad$Usage <- as.factor(adults.lad$Usage)

#create column for adult weight at emergence

adults.lad$emergence_weight <- adults.lad$Total.weight - adults.lad$Env.weight

adults.lad <- adults.lad[adults.lad$emergence_weight < 1,] #drop the one entry with emergence weight greater than one (this observation also has some other wonky data)

hist(adults.lad$emergence_weight) #looks nice and normally distributed

###### add new column to break up adults into OE infection categories; only consider monarchs with OE status of greater than 1 to be infected enough to warrant inclusion in "infected" category; anything lower probably indicates incidental transfer (this is anything less than 10 counted OE spores)

adults.lad$infection_status_binomial <- with(adults.lad, ifelse(OE >= 2, 1, 0))

adults.lad$Species <- revalue(adults.lad$Species, c('GOFR'='GOPH')) #combine GOFR into GOPH, as before

#create sympatric/allopatric column for adult data

adults.lad$sym.allo <- with(adults.lad, ifelse(Species %in% c('ASYR','AINC') & Mon.Pop == 'ENA', 'sympatric', ifelse(Species %in% c('ASFA','ASPEC') & Mon.Pop == 'CA', 'sympatric', ifelse(Species == 'ASCU' & Mon.Pop %in% c('PR','GU'), 'sympatric', ifelse(Species == 'GOPH' & Mon.Pop %in% c('AU','HI'), 'sympatric', 'allopatric')))))

###

ad.dates <- as.Date(adults.lad$Date, format='%m/%d/%y') #as before, create vector of dates

adults.lad$exp.days <- as.numeric(difftime(ad.dates, start_date , units = "days")) ## create new variable for experimental days (as before)

##########################################

## Analyses for adult performance metrics

##########################################

eclosion.model1 <- lmer(days_to_eclosion ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad[adults.lad$days_to_eclosion < 35,]) #restrict analysis to only adults that eclosed within 35 days of being added; this only excludes 7 observations, all of which were diseased in some way

summary(eclosion.model1)

Anova(eclosion.model1) #these are the values reported in the results section

### generate marginal means across all factors

ecl.spp <- as.data.frame(emmeans(eclosion.model1, 'Species'))
ecl.pop <- as.data.frame(emmeans(eclosion.model1, 'Mon.Pop'))
ecl.sym <- as.data.frame(emmeans(eclosion.model1, 'sym.allo'))
ecl.use <- as.data.frame(emmeans(eclosion.model1, 'Usage'))
ecl.gh <- as.data.frame(emmeans(eclosion.model1, 'GH'))
ecl.yr <- as.data.frame(emmeans(eclosion.model1, 'Year'))
ecl.sex <- as.data.frame(emmeans(eclosion.model1, 'Sex'))


#plot the raw data

agg.eclosion <- aggregate(days_to_eclosion ~ Species + Mon.Pop + sym.allo + maternal_family + Plant.ID, data = adults.lad[adults.lad$days_to_eclosion <35,], FUN = mean)

agg.eclosion$Mon.Pop <- factor(agg.eclosion$Mon.Pop, mon.order)
agg.eclosion$Species <- factor(agg.eclosion$Species, spp.order)

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
  theme(strip.background =element_rect(fill='cornsilk'))+
  scale_y_reverse()+
  annotate('rect', ymin = 15, ymax = 35, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 15, ymax = 35, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)

Fig4c1 <- Fig4c1 + facet_rep_wrap(~Species, repeat.tick.labels = 'bottom')

Fig4c1

#now make a figure of model output for days to eclosion focused on sym.allo contrast

plot.agg.ecl.sym <- as.data.frame(emmeans(eclosion.model1, 'sym.allo'))

(Fig4c2 <- ggplot(plot.agg.ecl.sym, aes( x = sym.allo, y = emmean, col = sym.allo))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
    scale_color_manual(values = c('black','red'), guide = F)+
    theme_bw()+
    ylab('Days to Eclosion')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    scale_y_reverse(labels = scales::number_format(accuracy = 0.2))+
    annotate('text', x = 1, y = 22.1, label = 't = 2.192\np = 0.028', size = 5)+
    annotate('text', x = 1.5, y = 22.4, label = '*', size = 8))

#ggsave(file="./Figures/Fig4c2.tiff", plot=Fig4c2, width=3, height=7)#ggsave(file="./Figures/Fig4c1.tiff", plot=Fig4c1, width=6, height=7)

### add a new column for whether monarch populations are ancestral

adults.lad$ancestral_pop <- ifelse(adults.lad$Mon.Pop %in% c('ENA','CA'), 'ancestral', 'derived')

### run model with this term instead of monarch population to get Fig4c3

eclosion.model.x1 <- lmer(days_to_eclosion ~ Species + ancestral_pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad[adults.lad$days_to_eclosion < 35,])

Anova(eclosion.model.x1)

plot.ecl.anc <- as.data.frame(emmeans(eclosion.model.x1, 'ancestral_pop'))

(Fig4c3 <- ggplot(plot.ecl.anc, aes( x = ancestral_pop, y = emmean, col = ancestral_pop))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
    scale_color_manual(values = c('blue','forestgreen'), guide = F)+
    theme_bw()+
    ylab('Days to Eclosion')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    scale_y_reverse()+
    annotate('text', x = 1, y = 22.8, label = 't = 2.213\np = 0.001', size = 5)+
    annotate('text', x = 1.5, y = 22.4, label = '*', size = 8))

#ggsave(file="./Figures/Fig4c3.tiff", plot=Fig4c3, width=3, height=7)

##Make supplemental figure based on marginal means for time eclosion across milkweed species

ecl.spp$Species <- factor(ecl.spp$Species, spp.order)

summary(glht(eclosion.model1, linfct = mcp(Species = "Tukey")))

(sup.fig.emm.ecl.spp <- ggplot(ecl.spp, aes( x = Species, y = (emmean), col = Species))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = (lower.CL), ymin = (upper.CL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Days to Eclosion')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Milkweed Species')+
    ylim(c(20,25.5))+
    scale_y_reverse()+
    scale_color_manual(values = mcb)+
    annotate('text', x = 1, y = 23, label = 'C', size = 6)+
    annotate('text', x = 2, y = 22, label = 'BC', size = 6)+
    annotate('text', x = 3, y = 21, label = 'AB', size = 6)+
    annotate('text', x = 4, y = 20.5, label = 'A', size = 6)+
    annotate('text', x = 5, y = 20, label = 'A', size = 6)+
    annotate('text', x = 6, y = 20.5, label = 'A', size = 6))

#ggsave(file="./Figures/sup.fig.emm.ecl.spp.tiff", plot = sup.fig.emm.ecl.spp, width=4, height=7)

##and finally do the same for monarch population

summary(glht(eclosion.model1, linfct = mcp(Mon.Pop = "Tukey")))

ecl.pop$Mon.Pop <-factor(ecl.pop$Mon.Pop, mon.order)

(sup.fig.emm.ecl.pop <- ggplot(ecl.pop, aes( x = Mon.Pop, y = (emmean)))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = (lower.CL), ymin = (upper.CL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Days to Eclosion')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Monarch Population')+
    ylim(c(20,25.5))+
    scale_y_reverse()+
    annotate('text', x = 1, y = 21.2, label = 'AB', size = 6)+
    annotate('text', x = 2, y = 20.5, label = 'BC', size = 6)+
    annotate('text', x = 3, y = 22, label = 'AB', size = 6)+
    annotate('text', x = 4, y = 21, label = 'A', size = 6)+
    annotate('text', x = 5, y = 20.8, label = 'AB', size = 6)+
    annotate('text', x = 6, y = 22.2, label = 'C', size = 6))

#ggsave(file="./Figures/sup.fig.emm.ecl.pop.tiff", plot = sup.fig.emm.ecl.pop, width=4, height=7)

#####################################

#Now create figures for mass at eclosion

adult.mass.model1 <- lmer(emergence_weight ~ Species + Mon.Pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad) #gives warning about singularity; however, simplifying model by dropping non-significant terms does not meaningfully change estimates of other parameters, so stick with this model for consistency
isSingular(adult.mass.model1) #isSingular also returns FALSE, so should be alright

summary(adult.mass.model1)
Anova(adult.mass.model1) #this is the output reported in the results

##get marginal means across factor levels

em.pop <- as.data.frame(emmeans(adult.mass.model1, 'Mon.Pop'))
em.spp <- as.data.frame(emmeans(adult.mass.model1, 'Species'))
em.sym <- as.data.frame(emmeans(adult.mass.model1, 'sym.allo'))
em.use <- as.data.frame(emmeans(adult.mass.model1, 'Usage'))
em.gh <- as.data.frame(emmeans(adult.mass.model1, 'GH'))
em.yr <- as.data.frame(emmeans(adult.mass.model1, 'Year'))
em.sex <- as.data.frame(emmeans(adult.mass.model1, 'Sex'))

#now aggregate raw data for figure

adult.mass.agg <- aggregate(emergence_weight ~ Species + Mon.Pop + sym.allo + maternal_family + Plant.ID, mean, data = adults.lad)

adult.mass.agg$Mon.Pop <- factor(adult.mass.agg$Mon.Pop, mon.order)
adult.mass.agg$Species <- factor(adult.mass.agg$Species, spp.order)

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
  theme(strip.background =element_rect(fill='cornsilk'))+
  annotate('rect', ymin = 0.1, ymax = 1, xmin = 0.5, xmax = 2.4, 
           col = 'blue', fill = NA, size = 0.7)+
  annotate('rect', ymin = 0.1, ymax = 1, xmin = 2.5, xmax = 6.5, 
           col = 'forestgreen', fill = NA, size = 0.7)

(Fig4d1 <- Fig4d1 + facet_rep_wrap(~Species, repeat.tick.labels = 'bottom'))

#ggsave(file="./Figures/Fig4d1.tiff", plot=Fig4d1, width=6, height=7)

### now highlight sym/allo contrast; first put fitted values with original dataframe

plot.agg.ecl.sym <- as.data.frame(emmeans(adult.mass.model1, 'sym.allo'))

(Fig4d2 <- ggplot(plot.agg.ecl.sym, aes( x = sym.allo, y = emmean, col = sym.allo))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
    scale_color_manual(values = c('black','red'), guide = F)+
    theme_bw()+
    ylab('Adult Eclosion Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    annotate('text', x = 1, y = 0.595, label = 't = 0.027\np > 0.9', size = 5))+
  ylim(c(0.54,0.6))

#ggsave(file="./Figures/Fig4d2.tiff", plot=Fig4d2, width=3, height=7)

### now build model with ancestral population status

adult.mass.model.x1 <- lmer(emergence_weight ~ Species + ancestral_pop + sym.allo + (1|Pop/Plant.ID) + (1|maternal_family)  + (1|Group)  + Usage + GH + Year + scale(OE) + Sex + scale(exp.days), data = adults.lad)

Anova(adult.mass.model.x1)

plot.eclm.anc <- as.data.frame(emmeans(adult.mass.model.x1, 'ancestral_pop'))

(Fig4d3 <- ggplot(plot.eclm.anc, aes( x = ancestral_pop, y = emmean, col = ancestral_pop))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.2, size=1)+
    scale_color_manual(values = c('blue','forestgreen'), guide = F)+
    theme_bw()+
    ylab('Adult Eclosion Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.title.x = element_blank())+
    annotate('text', x = 2, y = 0.595, label = 't = 1.668\np = 0.095', size = 5))

#ggsave(file="./Figures/Fig4d3.tiff", plot=Fig4d3, width=3, height=7)

em.spp$Species <- factor(em.spp$Species, spp.order)

summary(glht(adult.mass.model1, linfct = mcp(Species = "Tukey")))

(sup.fig.emm.em.spp <- ggplot(em.spp, aes( x = Species, y = (emmean), col = Species))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = (lower.CL), ymin = (upper.CL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Eclosion Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Milkweed Species')+
    scale_color_manual(values = mcb)+
    ylim(c(0.5,0.64))+
    annotate('text', x = 1, y = 0.6, label = 'A', size = 6)+
    annotate('text', x = 2, y = 0.632, label = 'A', size = 6)+
    annotate('text', x = 3, y = 0.61, label = 'A', size = 6)+
    annotate('text', x = 4, y = 0.61, label = 'A', size = 6)+
    annotate('text', x = 5, y = 0.605, label = 'A', size = 6)+
    annotate('text', x = 6, y = 0.61, label = 'A', size = 6)) ##Supplemental figure with marginal means across milkweed species

#ggsave(file="./Figures/sup.fig.emm.em.spp.tiff", plot = sup.fig.emm.em.spp, width=4, height=7)

em.pop$Mon.Pop <- factor(em.pop$Mon.Pop, mon.order)

summary(glht(adult.mass.model1, linfct = mcp(Mon.Pop = "Tukey")))

(sup.fig.emm.em.pop <- ggplot(em.pop, aes( x = Mon.Pop, y = (emmean)))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymax = (lower.CL), ymin = (upper.CL)), width = 0.2, size=1)+
    theme_bw()+
    ylab('Eclosion Mass (g)')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size =16),
          axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'none')+
    xlab('Monarch Population')+
    ylim(c(0.48,0.64))+
    annotate('text', x = 1, y = 0.62, label = 'A', size = 6)+
    annotate('text', x = 2, y = 0.63, label = 'A', size = 6)+
    annotate('text', x = 3, y = 0.62, label = 'AB', size = 6)+
    annotate('text', x = 4, y = 0.62, label = 'AB', size = 6)+
    annotate('text', x = 5, y = 0.605, label = 'AB', size = 6)+
    annotate('text', x = 6, y = 0.58, label = 'B', size = 6)) ## supplemental figure with marginal means across monarch populations

#ggsave(file="./Figures/sup.fig.emm.em.pop.tiff", plot = sup.fig.emm.em.pop, width=4, height=7)

#######################################################

## Now do analyses based on coefficient of variation in performance

#######################################################

cv.model <- aggregate(Weight ~ Species + Mon.Pop + sym.allo + Usage + Mon.ID + Pop + Plant.ID, FUN = cv, data = lad)

names(cv.model)[8] <- 'CV'

cv.weight.model1 <- lmer(CV ~ Species + Mon.Pop + sym.allo + Usage + (1|Mon.ID) + Usage + (1|Pop/Plant.ID), data = cv.model)

summary(cv.weight.model1) #when only non-zero caterpillar weights are included, sym/allo contrast is not significant, although direction of effect is still consistent with lower CV in sympatric combinations

#### now redo analysis with 0s included 

lad.zeros <- lad #create new DF to use for second test of CV
lad.zeros[is.na(lad.zeros$Weight),]$Weight <- 0 #cannot use this in aggregate with CV as the function, since it cannot be calculated for entries with mean of 0

zer.agg.mass <- aggregate(Weight ~ Species + Mon.Pop + Plant.ID + sym.allo + GH + Year + Usage + Mon.ID + Larvae + Pop, mean, data = lad.zeros)

repl.nas <- which(aggregate(Weight ~ Species + Mon.Pop + Plant.ID + sym.allo + GH + Year + Usage + Mon.ID + Larvae + Pop, mean, data = lad.zeros)$Weight == 0)

zer.cv.model <- lad.zeros[!lad.zeros$Plant.ID %in% zer.agg.mass[repl.nas,]$Plant.ID,]


zer.cv.model.agg <- aggregate(Weight ~ Species + Mon.Pop + Plant.ID + sym.allo + GH + Year + Usage + Mon.ID + Larvae + Pop, FUN = cv, data = zer.cv.model)

names(zer.cv.model.agg)[11] <- 'CV'

cv.weight.model2 <- lmer(CV ~ Species + Mon.Pop + sym.allo + (1|Mon.ID) + GH + Year + Usage + Larvae + (1|Pop/Plant.ID), data = zer.cv.model.agg)

summary(cv.weight.model2) #when zeros are included, returns a significant sym/allo contrast for caterpillars reared on sympatric versus allopatric host plants; CV is significantly lower in sympatric compared to allopatric combinations

Anova(cv.weight.model2) 

##marginal means for model of CV

(emmeans(cv.weight.model2, 'Species'))
(emmeans(cv.weight.model2, 'Mon.Pop'))
emmeans(cv.weight.model2, 'sym.allo')

MuMIn::r.squaredGLMM(cv.weight.model2) #fixed effects explain very little of the overall variation in CV (only 7.6%); random effects add more explanatory power

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

#ggsave(file="./Figures/Fig5a.tiff", plot=fig5a, width=6, height=4.5)

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

#ggsave(file="./Figures/Fig5b.tiff", plot=fig5b, width=6, height=4.5)

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

#ggsave(file="./Figures/Fig5c.tiff", plot=fig5c, width=3, height=9)

###### other things to check -- see if equalizing the number of maternal families changes anything about the variance analysis; this corresponds to analyses describes in one of the appendices

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
} ### heads up: this will take about 10-15 minutes with a sample of 1,000

mean(coefs.cv.model)
mean(t.cv.model)

pt(mean(coefs.cv.model), 598.3, lower.tail = TRUE) #average t statistic across iterations is -1.65; corresponds to a p value of 0.047 assuming a two-tailed distribution

summary(cv.weight.model2) #actual value from full model is -8.34

hist(coefs.cv.model, breaks = 50, xlab = 'Estimated coefficient', main = NULL)
abline(v = -8.34, col = 'red', lw = 2)
abline(v = (-8.34 + 3.918), col = 'red', lw = 1, lty =2)
abline(v = (-8.34 - 3.918), col = 'red', lw = 1, lty =2)
abline(v = mean(coefs.cv.model), col = 'blue', lw = 2)
abline(v = mean(coefs.cv.model)+sd(coefs.cv.model)*1.96, col = 'blue', lw = 1, lty = 2)
abline(v = mean(coefs.cv.model)-sd(coefs.cv.model)*1.96, col = 'blue', lw = 1, lty = 2)

###############################################################

#### Add in data on latex production and cardenolides in milkweeds

###############################################################

latex <- read.csv('./Monarch local adaptation 2018 - Latex.csv')

head(latex) #contains information on latex production at the level of each individual plant; reports a couple of metrics:
#1. Disc mass = mass of leaf discs used in cardenolide measurements
#2. Leaf mass = mass of excised leaf used when measuring latex production
#3. Paper.Start.Mass = mass of filter paper used for collecting latex, prior to latex addition
#4. Dmass = change in mass after latex added
#5. Latex.conc = Dmass/Leaf.mass * 100

latex.merge <- merge(lad, latex, by = 'Plant.ID') ## combine data on latex production with original larval performance dataframe

latex.merge$Latex.conc <- latex.merge$Latex.conc - min(latex.merge$Latex.conc, na.rm=T) #redefine latex.conc measurement

latex.merge$Dmass <- latex.merge$Dmass - min(latex.merge$Dmass, na.rm=T) #redefine dmass measurement

latex.merge$Species.x <- revalue(latex.merge$Species.x, c('GOFR' = 'GOPH')) #reclassify GOFR as GOPH

#####


lamagg <- aggregate(Latex.conc ~ Species.x + Plant.ID, mean, data = latex.merge) ### get latex concentration across all plants
lamagg$Species.x <- factor(lamagg$Species.x, spp.order) #reorder species after aggregating

hdh <- c('darkorange4','darkseagreen4','orange','darkolivegreen2','purple','blue')

### plot latex production by host species, same as what is shown in Fig3c

(fig3c <- ggplot(lamagg, aes(x = Species.x, y = Latex.conc, fill = Species.x))+
    geom_boxplot(outlier.shape = NA, alpha = 0.2)+theme_bw()+
    geom_point(position = position_jitter(width = 0.2), aes(col = Species.x), size = 1.2)+
    scale_color_manual(values = mcb)+
    scale_fill_manual(values = mcb)+
    theme(legend.position = 'none')+
    xlab('Milkweed Species')+
    ylab('Latex (mg latex / g plant tissue)')+
    theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16), axis.text.x = element_text(angle = 60, hjust = 1)))

#ggsave(file="./Figures/Fig3c.tiff", plot=fig3c, width=4, height=7)

##check for correlation between larval mass and latex production, pooled across all data

ggplot(na.omit(latex.merge[latex.merge$Latex.conc < 25,]), aes(x = (Latex.conc), y = Weight))+
  geom_point(aes(col = Species.x))+
  geom_smooth(method = 'lm', formula = y ~ poly(x, 2), col = 'black', lty = 4 )+
  scale_color_manual(values = mcb, name = 'Species')+
  theme_classic()+theme(legend.position = c(0.8,0.7),legend.text = element_text(size =14),
                        legend.title = element_blank())+
  xlab('Latex Production (mg/0.1g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16))+
  ylab('Caterpillar Mass (g)')

latex.merge$Species.x <- factor(latex.merge$Species.x, spp.order) #put species back in order for this dataframe

### generate fig3d, which shows larval performance as a function of latex production, faceted across different milkweed species

(fig3d <- ggplot(latex.merge, aes(x = Latex.conc, y = Weight))+
    geom_point(aes(col = Species.x))+
    facet_wrap(~Species.x, scales = 'free')+
    geom_smooth(method = 'lm', col = 'black')+
    scale_color_manual(values = mcb)+
    theme_bw()+theme(legend.position = 'none')+
    xlab('Latex Concentration (mg/0.1g)')+
    ylab('Day 8 Caterpillar Mass (g)')+
    theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16)))

#ggsave(file="./Figures/Fig3d.tiff", plot=fig3d, width=8, height=7)

###################################################

## Finally add cardenolide data

###################################################

cardenolides <- read.csv('./concatenated_cardenolides.csv') 

head(cardenolides) #shows all cardenolide data, with each peak assigned its own column; includes data on cardenolide concentrations in monarch wing tissue, which is the subject of a separate manuscript; this analysis focuses only on leaf concentrations

cardenolides$Year <- as.factor(cardenolides$Year)

masses <- read.csv('./masses.csv') ### load data on dry mass of leaf discs

cardenolides <- cardenolides[cardenolides$Tissue!='Wing',] #restrict analysis hereafter to only leaf tissues

cardenolides.standardized <- cbind(cardenolides[,1:7], cardenolides[,8:ncol(cardenolides)] / cardenolides$Digitoxin) ### standardize cardenolide data by dividing each sample by its corresponding digitoxin concentration (this is an internal standard that was included to allow for comparisons across samples)

standards <- cardenolides.standardized[cardenolides.standardized$Tissue=='Standard',]
standard.means <- colSums(standards[,8:ncol(cardenolides)]) / nrow(standards) #generate average values for standards, to be subtracted from the rest of the dataframe

drops <- sweep(cardenolides.standardized[,8:ncol(cardenolides.standardized)], 2, standard.means) #subtract standard means from each corresponding row

drops[drops<0] = 0 #replace all negative concentration values with 0

cards <- cbind(cardenolides[,1:7],drops)

table(cards$Tissue) #229 observations, 200 of which are leaf tissue

#### Now calculate concentrations based on leaf mass for each corresponding sample

colnames(masses)[1] <- 'Sample' 

ham <- merge(cards, masses, all.x = TRUE) #merge, but include rows where mass is NA

masses$Species <- revalue(masses$Species, c('GOFR' = 'GOPH')) #need to revalue so that GOFR is now GOPH

(average.masses <- aggregate(Mass ~ Species + Tissue, masses, mean)) #use these only for leaf tissues; for wings, use predicted regression values from correlation between hind wing area and mass

ham$spec.tiss <- paste(ham$Species, ham$Tissue, sep = "_") #create matching column
average.masses$spec.tiss <- paste(average.masses$Species, average.masses$Tissue, sep = "_") #matching column

replacements <- merge(ham[is.na(ham$Mass),], average.masses, all.x = TRUE) 

#now sub in actual values for replacements; logic here is to take samples that did not have leaf tissue mass and just interpolate their mass using the species-specific average calculated from other samples

replacements[replacements$spec.tiss==average.masses$spec.tiss[2],]$Mass <- average.masses$Mass[2]
replacements[replacements$spec.tiss==average.masses$spec.tiss[3],]$Mass <- average.masses$Mass[3]
replacements[replacements$spec.tiss==average.masses$spec.tiss[5],]$Mass <- average.masses$Mass[5]
replacements[replacements$spec.tiss==average.masses$spec.tiss[6],]$Mass <- average.masses$Mass[6]
replacements[replacements$spec.tiss==average.masses$spec.tiss[7],]$Mass <- average.masses$Mass[7]

ham <- ham[!is.na(ham$Mass),] #drop samples with no mass

ham <- rbind(ham, replacements) #add the replacement samples back to the original dataframe

ham$Digitoxin <- NULL #drop the standard colum
ham$spec.tiss <- NULL

ham <- ham[,c(1:5,79,8:77)] #restrict the dataframe to only the relevant columns

names(ham)

ham$cumulative <- rowSums(ham[,7:ncol(ham)]) #generate a cumulative measure based on mass-corrected concentrations

ham$concentration <- (0.15 * ham$cumulative) / (ham$Mass/1000 * 0.8) #convert into final measurement; 0.15 value corresponds to the known concentration of digitoxin initially injected; 0.8 value comes from the fact that during extraction procedure, only 800 uL of 1 mL supernatant was retained for analysis after leaf tissue was ground; values are divided by 1000 for unit conversion

aggregate(concentration ~ Tissue + Species, mean, data = ham) #sequestration ratio for AINC is about 2.9, ASCU is about 1.6, ASFA is about 2.1, ASPEC is about 9.1, ASYR is about 15.9, GOPH is about 0.8

leaves <- ham[ham$Tissue=='Leaf',] #already done, but reiterate that we only want to analyze leaf tissue

colnames(leaves)[1] <- 'Plant.ID' #change name of first column to allow for merging with original dataframe

plot.card.perf <- merge(leaves[,c(1,78)], lad, by = 'Plant.ID')

scaleFUN <- function(x) sprintf("%.0f", x)

plot.card.perf$Species <- factor(plot.card.perf$Species, spp.order) #rearrange species to plot

ggplot(plot.card.perf, aes(x = log(concentration*1000), y = Weight, col = Species))+
  geom_point()+
  facet_wrap(~Species, scales = 'free_x')+
  geom_smooth(method = 'lm')+
  scale_color_manual(values = mcb)+
  theme_bw()+theme(legend.position = 'none')+
  xlab('ln(Cardenolides (ug/g))')+
  ylab('Day 8 Caterpillar Mass (g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16))+
  scale_y_continuous(labels=scaleFUN) ## no evidence for consistent relationship between performance and cardenolide concentration within species

card.agg <- aggregate(Weight ~ Plant.ID + concentration + Species, plot.card.perf, mean) ##look across all species together

ggplot(card.agg, aes(x = log(concentration*1000), y = Weight))+
  geom_point(aes(col = Species))+
  scale_color_manual(values = mcb)+
  theme_bw()+theme(legend.position = 'none')+
  xlab('ln(Cardenolides (ug/g))')+
  ylab('Day 8 Caterpillar Mass (g)')+
  theme(axis.text = element_text(size=14),axis.title=element_text(size=16),strip.text = element_text(size=16))+
  geom_smooth(method = 'lm', col = 'black') ## if anything, performance is better on higher cardenolide species (but this may reflect negative correlation between latex production and cardenolide concentration)

#### Generate Fig3a based on cardenolide concentrations across species

leaves$Species <- factor(leaves$Species, c('ASPEC','ASYR','ASFA','AINC','ASCU','GOPH'))

(fig3a <- ggplot(na.omit(leaves[leaves$Tissue=='Leaf' & leaves$concentration >0.001,]), 
                 aes(x = Species , y = log(concentration*1000), fill = Species))+
    geom_boxplot(outlier.shape = NA, alpha = 0.1)+
    theme_bw()+
    geom_point(aes(color = Species), position = position_jitterdodge(1.5), pch = 20, size = 3)+
    scale_color_manual(values = hdh)+
    scale_fill_manual(values = hdh)+
    ylab('ln(Cardenolide Concentration (ug/g))')+
    theme(legend.position = 'none')+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    xlab('Milkweed Species'))

#### Do multivariate analysis of cardenolide profiles across species for Fig3c

library(vegan)
library(plotly)

ham <- ham[-which(ham$cumulative == 0),] #### first drop all rows that have 0 as the cumulative concentration (this corresponds to something like 3 samples that did not have any detectable peaks, all from A. incarnata and A. fascicularis)

leaves.compoundNMDS <- metaMDS(comm=ham[ham$Tissue=='Leaf',7:76], distance="bray", k=10, autotransform=F, trymax=100) ### run MDS, using bray-curtis dissimilarities to take into account differences in both compound presence-absence and abundance; use k = 10 axes and run across 100 random starts

stressplot(leaves.compoundNMDS) #looks good

leaves.data.scores <- as.data.frame(cbind(ham[ham$Tissue=='Leaf',1:6],scores(leaves.compoundNMDS))) #retain only numeric scores

leaves.data.scores <- leaves.data.scores[-which(leaves.data.scores$Species=='ASYR' & leaves.data.scores$NMDS1 < 0),] #drop one outlier that was probably mislabeled (ASYR sample that groups with GOPH)

leaf.axis1 <- aggregate(NMDS1 ~ Species, mean, data = leaves.data.scores) #calculate means for each species along MDS1
leaf.axis2 <- aggregate(NMDS2 ~ Species, mean, data = leaves.data.scores) #same for MDS2

leaf.axis1.sd <- aggregate(NMDS1 ~ Species, sd, data = leaves.data.scores) #now get SD for MDS1
names(leaf.axis1.sd)[2] <- 'SD1'
leaf.axis2.sd <- aggregate(NMDS2 ~ Species, sd, data = leaves.data.scores) #same for MDS2
names(leaf.axis2.sd)[2] <- 'SD2'

df_list <- list(leaf.axis1, leaf.axis2, leaf.axis1.sd,leaf.axis2.sd)
leaf.species <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)

leaf.species$Species <- factor(leaf.species$Species, c('ASFA','ASPEC','AINC','ASYR','ASCU','GOPH')) #rearrange spp.

mc2 <- c('darkolivegreen2','darkviolet','orange','darkorange4','darkseagreen4','blue') #define color vector

(fig3b <- ggplot(leaf.species, aes(x = NMDS1, y = NMDS2, group = Species, col = Species))+
    geom_point(size = 4, pch = 17)+
    geom_errorbar(aes(ymin = NMDS2 - SD2, ymax = NMDS2 + SD2))+
    geom_errorbarh(aes(xmin = NMDS1 - SD1, xmax = NMDS1 + SD1))+
    theme_bw()+
    scale_color_manual(values = mc2)+
    geom_point(data = leaves.data.scores, aes(x = NMDS1, y = NMDS2), alpha = 0.5)+
    theme(axis.text = element_text(size = 16), legend.text = element_text(size = 16),
          legend.title = element_blank(), axis.title = element_text(size =16))) ###Plot NMDS axes for milkweed cardenolide concentrations

#ggsave(file="./Figures/Fig3b.tiff", plot=fig3b, width=8, height=7)

dev.off()

### End of analysis


