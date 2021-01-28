##### Analysis of leaf morphology for island/mainland shrub comparisons

library(plyr)
library(lme4)
library(emmeans)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(merTools)

setwd("~/Documents/GitHub/Sites/manuscripts/Island_Mainland/")

leaf.data <- read.csv(file = "./data_files/chaparral_leaf_morphology.csv", header = T)

head(leaf.data) #dataframe with 5664 rows, each corresponding to a measurement from a single leaf

#########
#Description of columns
#########

#plant.ID = unique identifier for plants
#IM = categorical variable for whether sampled plant came from island or mainland
#Site = categorical variable defining 6 sites
#Exclsure = yes/no variable for whether sampled plants came from inside or outside of deer exclosure fences; relevant only for Catalina Island
#Species = categorical variable for species
#Plant = index for individual plants
#Position = categorical variable for position within plant canopy; lower = within deer browsing zone, upper = above deer browsing zone
#Aspect = caterogical variable defining which direction a particular plant is facing, if on a slope
#Elevation
#Diameter = measurements of stem diameter, in some cases multiple stems
#1st_year = categorical variable for whether leaf tissue is fully expanded
#Leaf_ID = grouping variable for leaves within branches
#Leaf_length = response variable for length along leaf's longest axis (excluding petiole)
#Leaf_Area
#Leaf_Area_corrected = leaf area after accounting for any tissue lost to herbivores
#Area_no_spines = leaf area after interpolating between marginal spine apices

#Combine Dendromecon into a single species, and replace species IDs with just genus

leaf.data$Species <- revalue(leaf.data$Species, c("Ceanothus_megacarpus" = "Ceanothus", "Cercocarpus_betuloides" = "Cercocarpus", "Dendromecon_harfordii" = "Dendromecon", "Dendromecon_rigida"= "Dendromecon", "Heteromeles_arbutifolia" = "Heteromeles", "Prunus_ilicifolia" = "Prunus"))

#rename sites as well

leaf.data$Site <- revalue(leaf.data$Site, c("Santa_Cruz" = "Santa Cruz", "Santa_Monicas"="Santa Monicas", "Santa_Rosa" = "Santa Rosa", "Stunt_Ranch" = "Stunt Ranch"))

#create new indexing variable for each plant

for( i in 1:nrow(leaf.data))
  leaf.data$plant.ID[i] = paste(leaf.data$Species[i], leaf.data$Plant[i], sep="_")

length(unique(leaf.data$plant.ID)) #339 unique plants sampled

#Convert aspect into a continuous variable; start with north facing aspect; assign due north a value of 3, due south a value of -3, NNW = 2, SSW = -2, etc. A more formal treatment would involve treating this as a circular variable, which would involve converting aspects into degrees (0-360) and then including separate sine and cosine model terms. This should work as a coarse simplification, especially since aspects were not recorded precisely in the first place. Here, the expectation is that northward aspects should be associated with increased leaf area and decreased specific leaf area.

leaf.data$Aspect.N = with(leaf.data, ifelse(Aspect == "N", 3,
        ifelse(Aspect == "NNW" | Aspect == "NNE", 2,
        ifelse(Aspect == "WNW" | Aspect == "ENE", 1,
        ifelse(Aspect == "E" | Aspect == "W", 0,
        ifelse(Aspect == "ESE" | Aspect == "WSW", -1,
        ifelse(Aspect == "SSE" | Aspect == "SSW", -2, -3)))))))

#East facing aspect

leaf.data$Aspect.E = with(leaf.data, ifelse(Aspect == "E", 3,
        ifelse(Aspect == "ENE" | Aspect == "ESE", 2,
        ifelse(Aspect == "NNE" | Aspect == "SSE", 1,
        ifelse(Aspect == "N" | Aspect == "S", 0,
        ifelse(Aspect == "NNW" | Aspect == "SSW", -1,
        ifelse(Aspect == "WSW" | Aspect == "WNW", -2, -3)))))))

### a few entries did not have an aspect recorded; for these, treat them as 0,0

leaf.data$Aspect.N[is.na(leaf.data$Aspect.N)] <- 0
leaf.data$Aspect.E[is.na(leaf.data$Aspect.E)] <- 0

#combine stem diameter measurements into a single measurement of basal area

leaf.data$Stem_Area = with(leaf.data, (Diameter1_.mm.*pi/100) + (Diameter2_.mm.*pi/100))

#add column for true leaf area (leaf area - internal damage)

leaf.data$True_area = (leaf.data$Leaf_Area_no.petiole - leaf.data$Internal_area_correction)

#add column for proportion of leaf lost to herbivory

leaf.data$Herbivory = 1 - (1 - (((leaf.data$Leaf_Area_corrected - leaf.data$True_area)/leaf.data$Leaf_Area_corrected) + (leaf.data$Internal_area_correction/leaf.data$True_area)))

hist(leaf.data$Herbivory, breaks = 25)

#add column for marginal spinescence

leaf.data$spines = ifelse(leaf.data$Herbivory > 0.1, NA, (leaf.data$Area_no_spines - leaf.data$Leaf_Area_corrected) / leaf.data$Leaf_Area_corrected * 100)

#create new indexing variable for each branch

for( i in 1:nrow(leaf.data))
  leaf.data$branch.ID[i] = paste(leaf.data$Species[i], leaf.data$Plant[i], leaf.data$Position[i], sep="_")

length(unique(leaf.data$branch.ID)) #585 unique branches sampled

#reorder sites and remove underscores

leaf.data$Site = factor(leaf.data$Site, c("Santa Cruz","Santa Rosa","Catalina","Santa Monicas","Stunt Ranch","Gaviota",'RSA', "SBBG"))

#filter out first year leaves; these are leaves had not yet fully expanded, which are expected to be morphologically different than mature leaves

table(leaf.data$X1st_year.) #will get rid of 809 observations

leaf.data <- leaf.data[leaf.data$X1st_year. == 0,]

leaf.data$Source.IM <- ifelse(leaf.data$Source %in% c('Santa Cruz','Santa Rosa','Catalina'), 'Island','Mainland')

####################################
#### Actual statistical analyses ###
####################################

## center and scale continuous predictors

head(leaf.data)

leaf.data$Stem_Area_scaled = scale(leaf.data$Stem_Area)
leaf.data$Elevation_scaled = scale(leaf.data$Elevation_.m.)

#for the time being, omit observations from the botanical gardens, and only focus on samples collected from their original location

insitu.leaf.data <- leaf.data[leaf.data$Site %in% c('Catalina','Gaviota','Santa Cruz','Santa Monicas','Santa Rosa','Stunt Ranch'),] #down to 4096 observations

###First question: is leaf area different between islands and mainland? what about between positions within the canopy?

### General comment on error variance structure:
#All species are represented at all sites (with a few minor exceptions, e.g. no Cercocarpus on Santa Rosa), so this should be a crossed and nested random effect; notation can be either (1|site) + (1|Species) or (1|Site:Species); plant ID is a nested random effect within species, and branch ID is a nested random effect within plant ID

model_leaf_area = lmer(True_area ~ IM + (IM|Species) + Position + (1|Site/plant.ID/branch.ID) + scale(Aspect.N) + scale(Aspect.E) + Elevation_scaled + Stem_Area_scaled, data = insitu.leaf.data) #this model includes the elevation at which plants were sampled and their basal stem area as covariates; elevation doesn't make much sense to include, since it's not really evenly sampled across sites (e.g. all samples from a site typically came from a fairly narrow elevational band); basal stem area is a proxy for the age of the plant, which could potentially be important, especially for island taxa with a history of recent herbivore removal. 

summary(model_leaf_area)

model_leaf_area_simple <- lmer(True_area ~ IM + (IM|Species) + Position + scale(Aspect.N) + scale(Aspect.E) + (1|Site/plant.ID/branch.ID), data = insitu.leaf.data) 

## Note: the first model gives a notice about the boundary fit being singular. This relates to the estimation of the random slopes for the site terms, which all end up being 0. If site is instead specified as a fixed effect, this warning goes away. However, this model structure (without plant.ID nested within site) reduces the available degrees of freedom for estimating other parameters. It's less important for this model and more so for other models using smaller subsets of data. The other possible issue is that there are only two levels (or sometimes just one) for branch ID within each plant.ID, though this is just a feature of the data structure, and removing branch ID (not supported) doesn't resolve the singularity warning. 

summary(model_leaf_area_simple) #Calculated degrees of freedom for island/mainland contrast should reflect that analysis is comparing island/mainland status across 5 species x 6 sites, which does appear to be the case. Analysis suggests that, overall, mainland species have smaller leaves, although the difference is not significant. Leaves from the upper portion of the canopy are also somewhat smaller. 

ranef(model_leaf_area_simple)$Species #shows species-level intercepts, compared to the overall model intercept of 7.809 cm^2. Also gives the island/mainland adjustment contrast specific to that species, which is in addition to the island/mainland overall contrast of -3.492. So, for example, mainland Prunus has a mean value of 7.809 - 3.409 + 6.372 - 5.062 = 5.710 By contrast, an island Prunus has a mean value of 7.809 + 6.372 = 14.181. The differences between other species are smaller, but the overall effect of mainland species having smaller leaves is consistent.

## get model estimated values for each combination of species*island/mainland

newdata <- with(insitu.leaf.data, expand.grid(Species=unique(Species), IM=unique(IM), Aspect.N = mean(Aspect.N), Aspect.E = mean(Aspect.E), Position = unique(Position)))

newdata$predicted.values <- predict(model_leaf_area_simple, newdata, re.form = ~(IM|Species))
(ndag <- aggregate(predicted.values ~ IM + Species, newdata, mean))

randomSims <- REsim(model_leaf_area_simple, n.sims = 500)

plotREsim(randomSims)

### General point: for global models that include all species, the inclusion of the (IM|Species) term means that models are comparing the island/mainland effect across each individual species, with the within-species estimates calculated as the pooled best-fit averages from all of the individual plants sampled from each island. This is why the degrees of freedom for the island/mainland term is comparatively low: fundamentally the model is estimating five within species slopes, and using these to get a pooled global island/mainland slope. Alternatively, and for the sake of clarity, it is also instructive to run species-specific models, which will not include the (IM|Species) term and will instead calculate d.f. from the number of sampled plants across island versus mainland. These models won't affect the within-species effect sizes, but the corresponding t and p values will be affected based on the d.f. calculation.

for(i in 1:length(unique(insitu.leaf.data$Species))){
  spx <- subset(insitu.leaf.data, insitu.leaf.data$Species == unique(insitu.leaf.data$Species)[i])
  model_spx <- lmer(True_area ~ IM + Position + (1|Site/plant.ID/branch.ID), data = spx)
  print(summary(model_spx)[10])}
#For leaf area, t and p values for island/mainland effect are below:
#Ceanothus IM: t = -3.853, p = 0.0125*
#Cercocarpus IM: t = -3.119, p = 0.0712
#Dendromecon IM: t = -2.588, p = 0.0787
#Heteromeles IM: t = -1.822, p = 0.1344
#Prunus IM: t = -4.126, p = 0.0151*

## alternative model formulation that does not treat species as a random intercept

model_leaf_area_simple_x <- lmer(True_area ~ IM*Species + Position + scale(Aspect.E) + scale(Aspect.N) + (1|Site/plant.ID/branch.ID), data = insitu.leaf.data)

summary(model_leaf_area_simple_x) #does not return an overall significant island/mainland effect and has much different d.f. calculation for the island/mainland term (this model treats each individual plant as a replicate, rather than the species to which they belong); does enable estimation of marginal means based on the species*IM term

emm.leaf.area <- emmeans(model_leaf_area_simple, specs = pairwise ~ IM, type = "response",
                         pbkrtest.limit = 4090) #takes about a minute to run
plot.emm.leaf.area <- data.frame(emm.leaf.area[[1]])

ggplot(plot.emm.leaf.area, aes(x = IM, y = emmean, col = IM))+
  theme_bw(base_size = 16)+
  geom_point(size = 4, pch = 22, aes(fill = IM), col = 'black')+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 1)+
  scale_color_manual(values = c('blue','red'))+
  scale_fill_manual(values = c('blue','red'))+
  facet_wrap(~Species, nrow = 1)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  theme(legend.position = 'bottom', legend.title = element_blank()) #estimated marginal means and confidence intervals for each species*IM combination; these values appear similar to what is shown in the boxplots of raw data. Can see that lower estimate overlaps 0 for Ceanothus, which highlights the issue of using this model format (allows for negative predicted values).

#### Make Figure 2a, with boxplots of leaf area 

#first, get branch level means for plotting

(leaf.area.agg = aggregate( True_area ~ Species + IM + Site + plant.ID, data = insitu.leaf.data, FUN = mean))

#pdf('./figures/Fig2a.pdf', height = 8, width = 8)
ggplot(leaf.area.agg, aes(x=Site, y=True_area))+
  geom_boxplot(aes(fill=IM), alpha = 0.3, outlier.color = 'white')+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21)+
  facet_wrap(~Species, nrow = 2, scales = 'free')+
  theme_bw()+
  labs(y = expression ('Leaf Area'~(cm^2)))+
  scale_fill_manual(values = c('blue','red'))+
  scale_color_manual(values = c('blue','red'))+
  theme(axis.text.x = element_text(size=12, angle=60,hjust=1), axis.text.y = element_text(size =12))+
  theme(legend.title=element_blank(),legend.text=element_text(size=14), 
        axis.title.y = element_text(size = 16))+
  theme(axis.title.x = element_blank())+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  theme(legend.position = c(0.8,0.2))
#dev.off()

# now, for comparison, get values from common garden shrubs

common.gardens <- leaf.data[leaf.data$Site %in% c('RSA','SBBG'),] #759 observations

with(aggregate(True_area ~ plant.ID + Species + Source, mean, data = common.gardens)
, table(Species, Source)) #shows the numbers of plants measured from common gardens, separated by their original place of collection. Noteworthy: no plants from Santa Rosa included in common garden sampling design, and lots of mainland locations that weren't part of original sampling scheme. 

with(aggregate(True_area ~ plant.ID + Species + Source.IM, mean, data = common.gardens)
     , table(Species, Source.IM)) 

model_leaf_area_cg <- lmer(True_area ~ Source.IM + (Source.IM|Species) + Position + Site + (1|plant.ID/branch.ID) + scale(Aspect.E) + scale(Aspect.N), data = common.gardens) #note: this model does not include a term for the plant's provenance, since there are 13 levels (many of which only have 1 or 2 observations); instead, provenances are pooled into island versus mainland (the Source.IM term).

emmeans(model_leaf_area_cg, specs = pairwise ~ Source.IM, type = "response")

summary(model_leaf_area_cg) #Mainland plants have modestly smaller leaves, although the difference is not significant
ranef(model_leaf_area_cg)[[3]]

agg.com.gar.leaf.area <- do.call(data.frame, aggregate(True_area ~ Source.IM + Species,  function(x) c(mean = mean(x), sd = sd(x)), data = aggregate(True_area ~ Species + Source.IM + Position, mean, data = common.gardens)))

ggplot(agg.com.gar.leaf.area, aes(x = Source.IM, y = True_area.mean))+
  geom_point()+
  geom_errorbar(aes(ymin = True_area.mean - True_area.sd, 
                    ymax = True_area.mean + True_area.sd), width = 0.2)+
  facet_wrap(~Species, nrow = 1)

###As above, run regressions for each species separately. Here, no site term is included, since all common garden samples came from one of two locations (RSA or SBBG). Only Prunus ends up being significantly different.

for(i in 1:length(unique(common.gardens$Species))){
  spx <- subset(common.gardens, common.gardens$Species == unique(common.gardens$Species)[i])
  model_spx <- lmer(True_area ~ Source.IM + Position + (1|plant.ID/branch.ID), data = spx)
  print(summary(model_spx)[10])}

#Common garden IM contrasts:
#Dendromecon: t = -0.016, p = 0.9879
#Heteromeles: t = -0.227, p = 0.8265
#Prunus: t = -4.880, p = 0.0011**
#Ceanothus: t = -1.547, p = 0.1827
#Cercocarpus: t = 0.3835, p = 0.7146

emm.area.cg <- emmeans(model_leaf_area_cg, specs = pairwise ~ Source.IM, type = "response") #pools together all species from the common garden setup
emm.area.cg <- data.frame(emm.area.cg[[1]])

#plot for visualization of pooled common garden data, but don't include as a figure
ggplot(emm.area.cg, aes(x = Source.IM, y = emmean, col = Source.IM))+
  theme_bw(base_size = 14)+
  geom_point(size = 4, pch = 21, aes(fill = Source.IM), col = 'black')+
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, size = 1)+
  scale_color_manual(values = c('blue','red'))+
  scale_fill_manual(values = c('blue','red'))+
  labs(y = expression ('Leaf Area'~(cm^2)))+
  theme(legend.position = 'none', axis.title.x = element_blank())+
  theme(strip.text = element_text(size=14, face = 'italic')) #averaged across levels of species

############################################

### Next, make figure for specific leaf area

#For this, first need to load in separate file that contains masses of pooled leaves

masses <- read.csv(file = './data_files/shrub_leaf_masses.csv', head = T)

for( i in 1:nrow(masses)){
  masses$SLA.ID[i] <- paste(masses$Species[i], masses$Number[i], masses$Position[i], sep = "_")
}

#create new column that combines plant.ID and canopy position for matching with SLA data

for( i in 1:nrow(leaf.data))
  leaf.data$SLA.ID[i] = paste(leaf.data$plant.ID[i], leaf.data$Position[i], sep="_")

leaf.data$SLA.ID #looks right

sla.dataframe <- merge(masses[,4:5], leaf.data, by = "SLA.ID")

sla.dataframe$SLA <- sla.dataframe$True_area / sla.dataframe$Mass #create new column for SLA (leaf area divided by tissue mass)

hist(sla.dataframe$SLA) #consider removing values above 25

#Use the same model structure as above, except here there are not multiple measurements within branches, so there are fewer overall observations (but the same number of replicates). 

sla.model <- aggregate(SLA ~ plant.ID + Position + IM + Site + Species + Aspect.N + Aspect.E + Elevation_scaled + Stem_Area_scaled, sla.dataframe[sla.dataframe$SLA < 20,], mean) #need to use aggregate before passing  to model, since for SLA, mass measurements were at the level of branches (sets of leaves) rather than individual leaves

model_specific_leaf_area = lmer(SLA ~ IM + (IM|Species) + Position + (1|Site/plant.ID) + scale(Aspect.N) + scale(Aspect.E), data = sla.model)

summary(model_specific_leaf_area) #Conclusions: no significant difference between island and mainland, though mainland plants tend to have somewhat thicker leaves. Leaves that come from the upper canopy, which are less susceptible to self-shading and have more direct sun exposure, are significantly thicker. Leaves from north facing and west facing slopes tend to be thinner. 

#get marginal means to compare IM 
(emm.sla <- emmeans(model_specific_leaf_area, specs = pairwise ~ IM, type = "response"))

sla.agg <- aggregate(SLA ~ Species + IM + Site + plant.ID, data = sla.dataframe, mean) #now aggregate again, but this time ignore canopy position for the sake of plotting

#pdf('./figures/Fig2b.pdf', height = 8, width = 8)
ggplot(sla.agg[sla.agg$SLA < 20,], aes(x=Site, y=SLA))+
  geom_boxplot(aes(fill=IM), alpha = 0.3, outlier.color = 'white')+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21)+
  facet_wrap(~Species, nrow = 2, scales = 'free')+
  theme_bw()+
  labs(y = expression ('Specific Leaf Area'~(cm^2/g)))+
  scale_fill_manual(values = c('blue','red'))+
  scale_color_manual(values = c('blue','red'))+
  theme(axis.text.x = element_text(size=12, angle=60,hjust=1), axis.text.y = element_text(size =12))+
  theme(legend.title=element_blank(),legend.text=element_text(size=14), 
        axis.title.y = element_text(size = 16))+
  theme(axis.title.x = element_blank())+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  theme(legend.position = c(0.8,0.2))
#dev.off()

## Species-level differences in SLA (none are significantly different)

for(i in 1:length(unique(sla.model$Species))){
  spx <- subset(sla.model, sla.model$Species == unique(sla.model$Species)[i])
  model_spx <- lmer(SLA ~ IM + Position + (1|Site/plant.ID), data = spx)
  print(summary(model_spx)[10])}

##island/mainland contrasts:
#Prunus: t = -1.796, p = 0.1421
#Cercocarpus: t = -0.890, p = 0.3763
#Ceanothus: t = -1.988, p = 0.1046
#Dendromecon: t = 0.542, p = 0.5908
#Heteromeles: t = -0.208, p = 0.8446

## no comparable measurements of SLA for common garden samples; samples from RSA were never massed

#### Next measurement for Figure 2c: marginal spinescence

leaf.data$spines

#here, restrict the analysis to only Prunus and Heteromeles, since they are the only two with meaningful spines that would deter vertebrate herbivores. 

spine.data <- leaf.data[leaf.data$Species %in% c('Prunus','Heteromeles'),]

hist(spine.data$spines)

#use the same model structure as before

#first define dataframe: omit SBBG and RSA as sites, and run na.omit to drop all NA values for spinescence (these result from instances where there was >10% herbivory)

spine.data.model <- spine.data[!spine.data$Site %in% c('RSA','SBBG'),]

model_leaf_spines = lmer(spines ~ IM*Position + (IM|Species) + (1|Site/plant.ID/branch.ID) + scale(Aspect.N) + scale(Aspect.E), data = spine.data.model) #

summary(model_leaf_spines) #Conclusion: mainland plants are significantly more spiny. Leaves from the upper canopy are less spiny (spinescence heteroblasty), but this pattern is more pronounced on the mainland.

spines.agg <- aggregate(spines ~ Species + Site + IM + plant.ID, data = spine.data.model, mean)

spines.agg

#pdf('./figures/Fig2c.pdf', height = 6, width = 5)
ggplot(spines.agg, aes(x=Site, y=spines))+
  geom_boxplot(aes(fill=IM), outlier.color = 'white', alpha = 0.3)+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21)+
  facet_wrap(~Species, nrow = 1, scales = 'free')+
  theme_bw()+
  ylab('Spinescence')+
  scale_fill_manual(values = c('blue','red'))+
  scale_color_manual(values = c('blue','red'))+
  theme(axis.text.x = element_text(size=12, angle=60,hjust=1), axis.text.y = element_text(size =12))+
  theme(legend.title=element_blank(),legend.text=element_text(size=14), 
        axis.title.y = element_text(size = 16))+
  theme(axis.title.x = element_blank())+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  theme(legend.position = 'bottom')
#dev.off()

#get marginal means for spinescence, averaged over levels of position (this will be shown separately). In this case, definitely looks like both heteromeles and prunus should be significantly differ between island and mainland.

(emm.spines.IM <- emmeans(model_leaf_spines, specs = pairwise ~ IM, type = "response")) #gives ridiculously large confidence intervals for mainland plants

### generate separate figure to include on the side showing comparison of upper vs. lower spinescence in island vs. mainland

emm.spines <- emmeans(model_leaf_spines, specs = pairwise ~ IM*Position, type = "response")
emm.spines <- data.frame(emm.spines[[2]])

#pdf('./figures/Fig2d.pdf', height = 6, width = 4)
ggplot(emm.spines, aes(x = Position, y = estimate, group = IM, col = IM))+
  theme_bw(base_size = 16)+
  geom_point(size = 4, pch = 22, aes(fill = IM), col = 'black')+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.2, size = 1)+
  scale_color_manual(values = c('blue','red'))+
  scale_fill_manual(values = c('blue','red'))+
  xlab('Canopy Position')+
  ylab('Spinescence')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(strip.text = element_text(size=16, face = 'italic')) #captures the pattern whereby the difference between upper and lower canopy position is more pronounced on the mainland, reflective of spinescence heteroblasty
#dev.off()

### Do species-level comparisons 

summary(lmer(spines ~ IM*Position + (1|Site/plant.ID/branch.ID), data = spine.data.model[spine.data.model$Species=='Heteromeles',]))
summary(lmer(spines ~ IM*Position + (1|Site/plant.ID/branch.ID), data = spine.data.model[spine.data.model$Species=='Prunus',]))

#Spinescence island/mainland contrasts:
#Heteromeles: t = 3.063, p = 0.0341
#Prunus: t = 8.008, p = 0.0002

## spinescence of common garden samples

spines.common.garden <- spine.data[spine.data$Site %in% c('SBBG','RSA'),] #leaves 314 samples, 24 of which are NAs
spines.common.garden$Source.IM <- ifelse(spines.common.garden$Source %in% c('Santa Rosa','Santa Cruz','Catalina'), 'Island', 'Mainland')

with(aggregate(spines ~ plant.ID + Species + Source.IM, mean, data = spines.common.garden)
     , table(Species, Source.IM)) #table showing number of plants sampled for spinescence at common gardens for both Prunus and Heteromeles

model_spines_cg <- lmer(spines ~ Source.IM*Position + (Source.IM|Species) + (1|plant.ID/branch.ID) + scale(Aspect.E) + scale(Aspect.N), data = spines.common.garden[!is.na(spines.common.garden$spines),]) #

summary(model_spines_cg)

(emm.spines.cg <- emmeans(model_spines_cg, specs = pairwise ~ Source.IM, type = "response"))
emm.spines.cg <- data.frame(emm.spines.cg[[1]])

#pdf('./figures/Fig3b.pdf', height = 3, width = 4)
ggplot(emm.spines.cg, aes(x = Source.IM, y = emmean, col = Source.IM))+
  theme_bw(base_size = 14)+
  geom_point(size = 4, pch = 21, aes(fill = Source.IM), col = 'black')+
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, size = 1)+
  scale_color_manual(values = c('blue','red'))+
  scale_fill_manual(values = c('blue','red'))+
  ylab('Spinescence')+
  theme(legend.position = 'none', axis.title.x = element_blank())+
  theme(strip.text = element_text(size=14, face = 'italic'))
#dev.off()

## again do each species separately for common garden samples

summary(lmer(spines ~ Source.IM*Position + (1|plant.ID/branch.ID) + Aspect.E + Aspect.N, data = spines.common.garden[spines.common.garden$Species=="Heteromeles",])) 
summary(lmer(spines ~ Source.IM*Position + (1|plant.ID/branch.ID) + Aspect.E + Aspect.N, data = spines.common.garden[spines.common.garden$Species=="Prunus",]))

#Common garden spinescence comparisons:
#Heteromeles: t = 1.952, p = 0.0949
#Prunus: t = 8.597, p < 0.0001

###### Finally, create figure for cyanogenic glycosides

cyanide <- read.csv(file = './data_files/cyanide_measurements.csv', header = T)

head(cyanide)

cyanide$IM <- ifelse(cyanide$Site %in% c('Catalina','Santa Cruz','Santa Rosa'), 'Island', 'Mainland')

# drop Cercocarpus samples for analysis, since they had too little to be useful

cyanide <- na.omit(cyanide[cyanide$Species != 'Cercocarpus',]) #drops 25 measurements, leaving 194 independent measures

#also drop samples collected from Botanic Gardens for now

insitu_cyanide <- cyanide[!cyanide$Source %in% c('SBBG','RSA'),] #drops another 29 measurements, leaving 165 measures

insitu_cyanide$PlantID <- paste(insitu_cyanide$Species, insitu_cyanide$PlantID, sep = "_")

##plot raw data, for the time ignoring tissue age and only focusing on island/mainland contrast

plot.cyanide <- aggregate(Tissue.Concentration ~ Species + Site + IM + PlantID, mean, data = insitu_cyanide)

plot.cyanide$Site <- factor(plot.cyanide$Site, levels = c('Santa Cruz','Santa Rosa','Catalina','Santa Monicas','Stunt Ranch','Gaviota'))

#pdf('./figures/Fig2e.pdf', height = 6, width = 5)
ggplot(plot.cyanide, aes(x=Site, y=Tissue.Concentration))+
  geom_boxplot(aes(fill=IM), outlier.color = 'white', alpha = 0.3)+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21)+
  facet_wrap(~Species, nrow = 1, scales = 'free')+
  theme_bw()+
  ylab('Cyanide Concentration (mg HCN / g)')+
  scale_fill_manual(values = c('blue','red'))+
  scale_color_manual(values = c('blue','red'))+
  theme(axis.text.x = element_text(size=12, angle=60,hjust=1), axis.text.y = element_text(size =12))+
  theme(legend.title=element_blank(),legend.text=element_text(size=14), 
        axis.title.y = element_text(size = 16))+
  theme(axis.title.x = element_blank())+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  theme(legend.position = 'bottom')
#dev.off()

## statistical model, same structure as before, except again treat species as a fixed effect this time, since there are only two levels. Also don't need to include branch ID, since only a single leaf was sampled from each branch.

model_cyanide = lmer(Tissue.Concentration ~ IM*Age + (IM|Species) + (1|Site/PlantID), data = insitu_cyanide)

summary(model_cyanide) #No significant difference between island and mainland, but younger tissue consistently has higher concentrations

emm.cyanide <- emmeans(model_cyanide, specs = pairwise ~ IM*Age, type = "response")
emm.cyanide <- data.frame(emm.cyanide[[1]])

#pdf('./figures/Fig2f.pdf', height = 6, width = 4)
ggplot(emm.cyanide, aes(x = Age, y = emmean, group = IM, col = IM))+
  theme_bw(base_size = 16)+
  geom_point(size = 4, pch = 22, aes(fill = IM), col = 'black')+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, size = 1)+
  scale_color_manual(values = c('blue','red'))+
  scale_fill_manual(values = c('blue','red'))+
  xlab('Tissue Age')+
  ylab('Cyanide Concentration (mg HCN / g)')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(strip.text = element_text(size=16, face = 'italic'))
#dev.off()

### now get species level cyanide comparisons

summary(lmer(Tissue.Concentration ~ IM*Age + (1|Site/PlantID), data = insitu_cyanide[insitu_cyanide$Species=='Heteromeles',]))
summary(lmer(Tissue.Concentration ~ IM*Age + (1|Site/PlantID), data = insitu_cyanide[insitu_cyanide$Species=='Prunus',]))

#Cyanide comparisons for island/mainland contrasts:
#Heteromeles: t = 2.430, p = 0.0553
#Prunus: t = 3.353, p = 0.0148

## Cyanide comparison for common garden plants

#overall model

summary(lmer(Tissue.Concentration ~ IM*Age + (IM|Species) + (1|PlantID), data = cyanide[cyanide$Source %in% c('SBBG','RSA'),]))
#no overall consistent difference between island and mainland

#species level comparisons of common garden plants

summary(lmer(Tissue.Concentration ~ IM + Age + (1|PlantID), data = cyanide[cyanide$Species =='Heteromeles' & cyanide$Source %in% c('SBBG','RSA'),]))
summary(lmer(Tissue.Concentration ~ IM + Age + (1|PlantID), data = cyanide[cyanide$Species =='Prunus' & cyanide$Source %in% c('SBBG','RSA'),]))
#Heteromeles: t = 1.316, p = 0.2129
#Prunus: t = 2.389, p = 0.0359*

###### Finally, for Catalina only, compare trait data inside versus outside of exclosures

catalina <- subset(leaf.data, leaf.data$Site == 'Catalina')

with(aggregate(True_area ~ plant.ID + Species + Exclosure, mean, data = catalina)
     , table(Species, Exclosure)) #Prunus almost all outside the exclosures, but other species spread pretty evenly

#build model to see if leaf area was higher inside exclosures

summary(lmer(True_area ~ Exclosure + (Exclosure|Species) + Position + Aspect.E + Aspect.N + (1|plant.ID/branch.ID), data = catalina)) #df looks correct. No significant difference in leaf area inside versus outside of exclosures, though there is a slight positive effect of being inside an exclosure for leaf area.

ranef((lmer(True_area ~ Exclosure + (Exclosure|Species) + Position + Aspect.E + Aspect.N + (1|plant.ID/branch.ID), data = catalina))) #only big increase is for Prunus, though sample sizes inside exclosures are too small (n = 1) to be meaningful

#next test for SLA differences

summary(lmer(SLA ~ Exclosure + (Exclosure|Species) + Position + Aspect.E + Aspect.N + (1|plant.ID), data = sla.dataframe[sla.dataframe$Site=='Catalina',])) #upper leaves have lower SLA (as expected), but again, no significant differences inside versus outside exclosures

#now for spinescence

summary(lmer(spines ~ Exclosure + (Exclosure|Species) + Position + Aspect.E + Aspect.N + (1|plant.ID/branch.ID), data = catalina[catalina$Species %in% c('Heteromeles','Prunus'),])) #lower spines on north-facing plants and branches from the upper canopy, but again no difference between inside versus outside exclosures

#lastly test CNglcs

cng.catalina <- cyanide[cyanide$Site=='Catalina' & cyanide$Source == 'Catalina',] #leaves 36 observations from 19 plants (11 Heteromeles, 8 Prunus), but only six of these (5 Heteromeles, 1 Prunus) are from inside of exclosures

cng.catalina$plant.ID <- paste(cng.catalina$Species, cng.catalina$PlantID, sep = "_")

cng.catalina.model = aggregate(Tissue.Concentration ~ IM + Age + Species + plant.ID + Exclosure, mean, data = merge(cng.catalina, catalina[,colnames(catalina) %in% c('Exclosure','plant.ID')], by = 'plant.ID', all = F))

summary(lmer(Tissue.Concentration ~ Exclosure + Age + (Exclosure|Species) + (1|plant.ID), data = cng.catalina.model)) #no difference inside versus outside of exclosures, but young tissue has higher concentrations, as expected


########

##Herbivory data

hist(leaf.data$Herbivory) #strongly zero-inflated, probably makes sense to start using a negative binomial (?) distribution as a first pass

leaf.data$Herbivory.binom <- ifelse(leaf.data$Herbivory > 0, 1, 0 )

Herbivory.model <- glmer(Herbivory.binom ~ IM + Aspect.N + Aspect.E + Position + (IM|Species) + (1|Site/plant.ID/branch.ID), family = 'binomial', leaf.data)

summary(Herbivory.model) #no indication that island plants are any more likely to have evidence of herbivory; only significant model term suggests that leaves on westward facing aspects are more likely to experience herbivory

#####

#get effect size measurements for all traits

#start with leaf area

library(effectsize)

cohens_d(leaf.area.agg$True_area, leaf.area.agg$IM) #overall effect size for insularity: 0.60, with 95% confidence intervals between 0.36-0.83; don't actually report this value, though, since it averages all species and pools their variances in a way that doesn't account for small differences in sample sizes across island vs. mainland sites

output <- matrix(ncol = 3, nrow = 5)
for(i in 1:length(levels(leaf.area.agg$Species))){
  species <- levels(leaf.area.agg$Species)[i]
  cohd <- cohens_d(leaf.area.agg[leaf.area.agg$Species==species,]$True_area, leaf.area.agg[leaf.area.agg$Species==species,]$IM)
  output[i,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}

plot.cohd.1 <- data.frame('species' = levels(leaf.area.agg$Species), 'cohens_D' = output[,1],
                        'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'IS')

#Ceanothus: 1.58, 0.97-2.18
#Cercocarpus: 0.99, 0.39-1.57
#Dendromecon: 1.37, 0.75-1.98
#Heteromeles: 0.60, 0.10-1.10
#Prunus: 1.86, 1.24-2.47

#now for common garden

leaf.area.cg.agg <- aggregate(True_area ~ Source.IM + Species + plant.ID, data = common.gardens, mean)

cohens_d(True_area ~ Source.IM, data = leaf.area.cg.agg)
#get nearly identical effect size, but with wider confidence intervals: 0.59, (-0.05) - (1.22)

output <- matrix(ncol = 3, nrow = 5)
for(i in 1:length(levels(leaf.area.cg.agg$Species))){
  species <- levels(leaf.area.cg.agg$Species)[i]
  cohd <- cohens_d(True_area ~ Source.IM, data = leaf.area.cg.agg[leaf.area.cg.agg$Species==species,])
  output[i,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}

plot.cohd.2 <- data.frame('species' = levels(leaf.area.agg$Species), 'cohens_D' = output[,1],
                          'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'CG')

#Ceanothus: 1.12, (-0.57) - (1.12)
#Cercocarpus: -0.27, (-1.65, 1.14)
#Dendromecon: -0.03, (-1.73, 1.67)
#Heteromeles: 0.21, (-1.12, 1.52)
#Prunus: 3.09, (1.12, 4.99)

plot.cohd <- rbind(plot.cohd.1, plot.cohd.2)

#pdf('./figures/test.Fig.3a.pdf', width = 9, height = 4)
ggplot(plot.cohd, aes(x = location, y = cohens_D, col = location))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1, size = 1)+
  geom_point(size = 4, pch = 21, col = 'black', aes(fill = location))+
  theme_bw(base_size = 16)+
  facet_wrap(~species, nrow = 1)+
  theme(legend.position = 'none', axis.title.x = element_blank())+
  ylab('Insularity Effect Size (Leaf Area)')+
  scale_color_manual(values = c('purple','dodgerblue'))+
  scale_fill_manual(values = c('purple','dodgerblue'))+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  ylim(c(-2,5.5))
#dev.off()

#### now do effect size for SLA

cohens_d(SLA ~ IM, data = sla.agg)
#effect size is 0.23, with CIs of 0.00-0.47

calc_cohens_d <- function(dataset, response, comparison) {
  for(i in 1:length(levels(dataset$Species))){
    species <- dataset$Species[i]
    cohenD <- cohens_d(response ~ comparison, data = dataset)
    print(cohenD)
  }
}

calc_cohens_d(dataset = sla.agg, response = SLA, comparison = IM)

output <- matrix(ncol = 3, nrow = 5)
for(i in 1:length(levels(sla.agg$Species))){
  species <- levels(sla.agg$Species)[i]
  cohd <- cohens_d(SLA ~ IM, data = sla.agg[sla.agg$Species==species,])
  output[i,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}
#Ceanothus: 0.90, 0.34-1.46
#Cercocarpus: 0.53, -0.04-1.09
#Dendromecon: 0.00, -0.55 - 0.57
#Heteromeles: 0.11, -0.38 - 0.61
#Prunus: 0.65, 0.12 - 1.17

plot.cohd.sla <- data.frame('species' = levels(sla.agg$Species), 'cohens_D' = output[,1],
                            'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'IS')

#pdf('./figures/test.Fig.3b.pdf', width = 9, height = 4)
ggplot(plot.cohd.sla, aes(x = location, y = cohens_D))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1, size = 1, col = 'purple')+
  geom_point(size = 4, pch = 21, col = 'black', fill = 'purple')+
  theme_bw(base_size = 16)+
  facet_wrap(~species, nrow = 1)+
  theme(legend.position = 'none', axis.title.x = element_blank())+
  ylab('Insularity Effect Size (SLA)')+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  ylim(c(-1,2))
#dev.off()

#no comparable measurements for common garden plants for SLA

#### Now for spinescence

cohens_d(spines ~ IM, data = spines.agg)
#Much stronger overall effect: -1.73, (-2.19) - (-1.27)
output <- matrix(ncol = 3, nrow = 2)
for(i in 4:5){
  species <- levels(spines.agg$Species)[i]
  cohd <- cohens_d(spines ~ IM, data = spines.agg[spines.agg$Species==species,])
  output[i-3,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}

plot.cohd.1.spines <- data.frame('species' = levels(spines.agg$Species)[4:5], 'cohens_D' = output[,1],
                                 'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'IS')

#Heteromeles: -1.61, (-2.21) - (-1.01)
#Prunus: -3.40, (-4.34) - (-2.44)

spines.common.garden.agg <- aggregate(spines ~ Source.IM + Species + plant.ID, data = spines.common.garden, mean)

output <- matrix(ncol = 3, nrow = 2)
for(i in 4:5){
  species <- levels(spines.common.garden.agg$Species)[i]
  cohd <- cohens_d(spines ~ Source.IM, data = 
                     spines.common.garden.agg[spines.common.garden.agg$Species==species,])
  output[i-3,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}

plot.cohd.2.spines <- data.frame('species' = levels(spines.agg$Species)[4:5], 'cohens_D' = output[,1],
                                 'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'CG')

#Heteromeles: -0.23, (-1.54) - (1.10)
#Prunus: -4.36, (-6.77) - (-1.90)

plot.cohd.spines <- rbind(plot.cohd.1.spines, plot.cohd.2.spines)

#pdf('./figures/test.Fig.3c.pdf', width = 4, height = 4)
ggplot(plot.cohd.spines, aes(x = location, y = cohens_D, col = location))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1, size = 1)+
  geom_point(size = 4, pch = 21, col = 'black', aes(fill = location))+
  theme_bw(base_size = 16)+
  facet_wrap(~species, nrow = 1)+
  theme(legend.position = 'none', axis.title.x = element_blank())+
  ylab('Insularity Effect Size (Spinescence)')+
  scale_color_manual(values = c('purple','dodgerblue'))+
  scale_fill_manual(values = c('purple','dodgerblue'))+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  ylim(c(-9,3))
#dev.off()

#and finally, do the same for cyanogenic glycosides

cyanide.in.situ <- cyanide[!cyanide$Source %in% c('SBBG','RSA'),]

cyanide.agg <- aggregate(Tissue.Concentration ~ IM + Species + PlantID, mean, data = cyanide.in.situ) #leaves observations from 87 individual plants

cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg)
#substantially negative: -0.89, (-1.33) - (-0.43)

cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg[cyanide.agg$Species=='Heteromeles',])
#for Heteromeles: -1.05, (-1.70) - (-0.38)
cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg[cyanide.agg$Species=='Prunus',])
#for Prunus: -1.34, (-1.34) - (-0.67)

cyanide.agg.cg <- aggregate(Tissue.Concentration ~ PlantID + Species + IM, data = cyanide[cyanide$Source %in% c('SBBG','RSA'),], mean) #21 observations

cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg.cg[cyanide.agg.cg$Species=='Heteromeles',])
#still get pretty pronounced effect size for Heteromeles: -0.98, (-2.30) - (0.40)
cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg.cg[cyanide.agg.cg$Species=='Prunus',])
#and even stronger for Prunus: -1.77, (-3.17) - (-0.30)

output <- matrix(ncol = 3, nrow = 2)
for(i in 2:3){
  species <- levels(cyanide.agg$Species)[i]
  cohd <- cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg[cyanide.agg$Species==species,])
  output[i-1,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}

plot.cohd.1.cng <- data.frame('species' = levels(cyanide.agg$Species)[2:3], 'cohens_D' = output[,1],
                                 'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'IS')

output <- matrix(ncol = 3, nrow = 2)
for(i in 2:3){
  species <- levels(cyanide.agg.cg$Species)[i]
  cohd <- cohens_d(Tissue.Concentration ~ IM, data = cyanide.agg.cg[cyanide.agg.cg$Species==species,])
  output[i-1,] <- c(cohd[[1]], cohd[[3]], cohd[[4]])
}

plot.cohd.2.cng <- data.frame('species' = levels(cyanide.agg$Species)[2:3], 'cohens_D' = output[,1],
                              'lower_CI' = output[,2], 'upper_CI' = output[,3], 'location' = 'CG')

plot.cohd.cng <- rbind(plot.cohd.1.cng, plot.cohd.2.cng)

#pdf('./figures/test.Fig.3d.pdf', width = 4, height = 4)
ggplot(plot.cohd.cng, aes(x = location, y = cohens_D, col = location))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1, size = 1)+
  geom_point(size = 4, pch = 21, col = 'black', aes(fill = location))+
  theme_bw(base_size = 16)+
  facet_wrap(~species, nrow = 1)+
  theme(legend.position = 'none', axis.title.x = element_blank())+
  ylab('Insularity Effect Size (CNglc Conc.)')+
  scale_color_manual(values = c('purple','dodgerblue'))+
  scale_fill_manual(values = c('purple','dodgerblue'))+
  theme(strip.text = element_text(size=16, face = 'italic'))+
  ylim(c(-4,1))
#dev.off()
