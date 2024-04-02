##### Analysis of leaf morphology for island/mainland shrub comparisons

#load libraries

library(plyr)
library(lme4)
library(emmeans)
library(ggplot2)
library(lmerTest)
library(ggeffects)

leaf.data <- read.csv(file = "./data_files/Shrubs/Morphology/chaparral_leaf_morphology.csv", header = T)

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

leaf.data$Species <- as.factor(leaf.data$Species)

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

model_leaf_area_full = lmer(True_area ~ IM + (IM|Species) + Position + (1|Site/plant.ID/branch.ID) + scale(Aspect.N) + scale(Aspect.E) + Elevation_scaled + Stem_Area_scaled, data = insitu.leaf.data) #this model includes the elevation at which plants were sampled and their basal stem area as covariates; elevation doesn't make much sense to include, since it's not really evenly sampled across sites (e.g. all samples from a site typically came from a fairly narrow elevational band); basal stem area is a proxy for the age of the plant, which could potentially be important, especially for island taxa with a history of recent herbivore removal. 

summary(model_leaf_area_full) 

model_leaf_area_simple <- lmer(log(True_area) ~ IM*Species + scale(Aspect.N) + scale(Aspect.E) + Position + (1|Site/plant.ID/branch.ID), data = insitu.leaf.data)

## Note: the first model gives a notice about the boundary fit being singular. This relates to the estimation of the random slopes for the site terms, which all end up being 0. If site is instead specified as a fixed effect, this warning goes away. However, this model structure (without plant.ID nested within site) reduces the available degrees of freedom for estimating other parameters. It's less important for this model and more so for other models using smaller subsets of data. The other possible issue is that there are only two levels (or sometimes just one) for branch ID within each plant.ID, though this is just a feature of the data structure, and removing branch ID (not supported) doesn't resolve the singularity warning. 

summary(model_leaf_area_simple) #note that the island/mainland contrast reported here is not informative on its own, since it is based on the reference species level (Ceanothus)

emm.leaf.area.spp <- emmeans(model_leaf_area_simple, specs = ~ IM + Species, type = "response", pbkrtest.limit = 4090) #get estimated marginal means for each species

emm.leaf.area.spp.contrasts <- contrast(object = emm.leaf.area.spp, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.leaf.area.spp <- data.frame(emmeans(model_leaf_area_simple, specs = pairwise ~ IM + Species, type = "response", pbkrtest.limit = 4090)[[1]]) #generate data frame for plotting

(leaf_area_model_output <- data.frame(predict_response(model = model_leaf_area_simple, terms = c("IM", "Species"), type = "random", back_transform = T, interval = "confidence"))) #similar to values produced by emmeans, but for sake of simplicity, do not use these values

(fig2a <- ggplot(data = plot.emm.leaf.area.spp, aes(x = IM, y = response, col = IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  labs(y = expression ('Leaf Area'~(cm^2)))+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))

#get branch level means for plotting

(leaf.area.agg = aggregate( True_area ~ Species + IM + Site + plant.ID, data = insitu.leaf.data, FUN = mean))

#pdf('./figures/Fig2a.pdf', height = 4, width = 10)
fig2a+
  geom_point(data = leaf.area.agg, aes(x = IM, y = True_area, fill = IM), size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red'))
#dev.off()

#### Plot raw data for supplemental figure

#pdf('./figures/Supplemental_Figures/FigS7A.pdf', height = 8, width = 8)
ggplot(leaf.area.agg, aes(x=Site, y=True_area))+
  geom_boxplot(aes(fill=IM), alpha = 0.3, outlier.color = 'white')+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21, size = 3, alpha = 0.5)+
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

model_leaf_area_common_garden <- lmer(log(True_area) ~ Source.IM*Species + Position + Site + (1|plant.ID/branch.ID) + scale(Aspect.E) + scale(Aspect.N), data = common.gardens) #note: this model does not include a term for the plant's provenance, since there are 13 levels (many of which only have 1 or 2 observations); instead, provenances are pooled into island versus mainland (the Source.IM term).

summary(model_leaf_area_common_garden)

(leaf_area_cg_model_output <- data.frame(predict_response(model = model_leaf_area_common_garden, terms = c("Source.IM", "Species"), type = "random", back_transform = T, interval = "confidence"))) 

emm.leaf.area.spp.common.garden <- emmeans(model_leaf_area_common_garden, specs = ~ Source.IM + Species, type = "response") #get estimated marginal means for each species

emm.leaf.area.spp.common.garden.contrasts <- contrast(object = emm.leaf.area.spp.common.garden, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.leaf.area.spp.common.garden <- data.frame(emmeans(model_leaf_area_common_garden, specs = pairwise ~ Source.IM + Species, type = "response")[[1]]) #generate data frame for plotting

(figS2a <- ggplot(data = plot.emm.leaf.area.spp.common.garden, aes(x = Source.IM, y = response, col = Source.IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = Source.IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  labs(y = expression ('Leaf Area'~(cm^2)))+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))

cg.leaf.data <- leaf.data[leaf.data$Site %in% c('SBBG','RSA'),]

(leaf.area.cg.agg = aggregate( True_area ~ Species + Source.IM + Site + plant.ID, data = cg.leaf.data, FUN = mean))

#pdf('./figures/Supplemental_Figures/FigS8A.pdf', height = 4, width = 10)
figS2a+
  geom_point(data = leaf.area.cg.agg, aes(x = Source.IM, y = True_area, fill = Source.IM), 
             size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red'))
#dev.off()

############################################

### Next, make figure for specific leaf area

#For this, first need to load in separate file that contains masses of pooled leaves

masses <- read.csv(file = './data_files/Shrubs/Morphology/shrub_leaf_masses.csv', head = T)

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

model_specific_leaf_area = lmer(SLA ~ IM*Species + Position + (1|Site/plant.ID) + scale(Aspect.N) + scale(Aspect.E), data = sla.model)

summary(model_specific_leaf_area) #Conclusions: somewhat significant difference between island and mainland. Leaves that come from the upper canopy, which are less susceptible to self-shading and have more direct sun exposure, are significantly thicker. Leaves from north facing and west facing slopes tend to be thinner. 

anova(model_specific_leaf_area)

emm.sla.spp <- emmeans(model_specific_leaf_area, specs = ~ IM + Species, type = "response") #get estimated marginal means for each species

emm.sla.spp.contrasts <- contrast(object = emm.sla.spp, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.sla.spp <- data.frame(emmeans(model_specific_leaf_area, specs = pairwise ~ IM + Species, type = "response")[[1]]) #generate data frame for plotting

(fig2b <- ggplot(data = plot.emm.sla.spp, aes(x = IM, y = emmean, col = IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  labs(y = expression ('SLA'~(cm^2/g)))+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))

sla.agg <- aggregate(SLA ~ Species + IM + Site + plant.ID, data = sla.dataframe, mean) #now aggregate again, but this time ignore canopy position for the sake of plotting

#pdf('./figures/Fig2b.pdf', height = 4, width = 10)
fig2b+
  geom_point(data = sla.agg[sla.agg$SLA < 20,], aes(x = IM, y = SLA, fill = IM), size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red')) #note that two high outlier values for Cercocarpus are truncated
#dev.off()  

#plot raw data

#pdf('./figures/Supplemental_Figures/FigS7B.pdf', height = 8, width = 8)
ggplot(sla.agg[sla.agg$SLA < 20,], aes(x=Site, y=SLA))+
  geom_boxplot(aes(fill=IM), alpha = 0.3, outlier.color = 'white')+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21, size = 3, alpha = 0.5)+
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

#### Next measurement for Figure 2c: marginal spinescence

leaf.data$spines

#here, restrict the analysis to only Prunus and Heteromeles, since they are the only two with meaningful spines that would deter vertebrate herbivores. 

spine.data <- leaf.data[leaf.data$Species %in% c('Prunus','Heteromeles'),]

hist(spine.data$spines)

#use the same model structure as before

#first define dataframe: omit SBBG and RSA as sites, and run na.omit to drop all NA values for spinescence (these result from instances where there was >10% herbivory)

spine.data.model <- spine.data[!spine.data$Site %in% c('RSA','SBBG'),]

model_leaf_spines = lmer(spines ~ IM*Position + IM*Species + (1|Site/plant.ID/branch.ID) + scale(Aspect.N) + scale(Aspect.E), data = spine.data.model) #

summary(model_leaf_spines) #Conclusion: mainland plants are not significantly more spiny overall. Leaves from the upper canopy are less spiny (spinescence heteroblasty), but this pattern is more pronounced on the mainland.

emm.spines.spp <- emmeans(model_leaf_spines, specs = ~ IM + Species, type = "response") #get estimated marginal means for each species

emm.spines.spp.contrast <- contrast(object = emm.spines.spp, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.spines <- data.frame(emmeans(model_leaf_spines, specs = pairwise ~ IM + Species, type = "response")[[1]]) #generate data frame for plotting

(fig2c <- ggplot(data = plot.emm.spines, aes(x = IM, y = emmean, col = IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  ylab('Spinescence (%)')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))

#plot raw data

spines.agg <- aggregate(spines ~ Species + Site + IM + plant.ID, data = spine.data.model, mean)

#pdf('./figures/Fig2c.pdf', height = 4, width = 4)
fig2c+
  geom_point(data = spines.agg, aes(x = IM, y = spines, fill = IM), size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red'))
#dev.off()

#pdf('./figures/Supplemental_Figures/FigS7C.pdf', height = 6, width = 5)
ggplot(spines.agg, aes(x=Site, y=spines))+
  geom_boxplot(aes(fill=IM), outlier.color = 'white', alpha = 0.3)+
  geom_point(aes(fill=IM), position = position_jitterdodge(0.5), pch = 21, size = 3, alpha = 0.5)+
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

#and finally get marginal means for canopy position for spinescence

emm.spines.position <- data.frame(emmeans(model_leaf_spines, specs = ~ IM + Species + Position, type = "response"))

#pdf('./figures/Supplemental_Figures/FigS9.pdf', height = 5, width = 6)
ggplot(data = emm.spines.position, aes(x = Position, y = emmean, col = IM, group = IM))+
  theme_bw(base_size = 14)+
  geom_point(position = position_dodge(0.1), size = 4)+
  geom_line(position = position_dodge(0.1))+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.1))+
  facet_wrap(~Species, scales = 'free_y')+
  ylab('Spinescence (%)')+
  xlab('Canopy Position')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  scale_color_manual(values = c('blue','red'))+
  theme(strip.text = element_text(face = "italic", size = 16))
#dev.off()

## spinescence of common garden samples

spines.common.garden <- spine.data[spine.data$Site %in% c('SBBG','RSA'),] #leaves 314 samples, 24 of which are NAs
spines.common.garden$Source.IM <- ifelse(spines.common.garden$Source %in% c('Santa Rosa','Santa Cruz','Catalina'), 'Island', 'Mainland')

with(aggregate(spines ~ plant.ID + Species + Source.IM, mean, data = spines.common.garden)
     , table(Species, Source.IM)) #table showing number of plants sampled for spinescence at common gardens for both Prunus and Heteromeles

model_spines_common_garden <- lmer(log(spines + 0.001) ~ Source.IM*Position + Source.IM*Species + (1|plant.ID/branch.ID) + scale(Aspect.E) + scale(Aspect.N), data = spines.common.garden[!is.na(spines.common.garden$spines),]) #

summary(model_spines_common_garden)

emm.spines.spp.common.garden <- emmeans(model_spines_common_garden, specs = ~ Source.IM + Species, type = "response") #get estimated marginal means for each species

emm.spines.spp.common.garden.contrast <- contrast(object = emm.spines.spp.common.garden, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.spines.common.garden <- data.frame(emmeans(model_spines_common_garden, specs = pairwise ~ Source.IM + Species, type = "response")[[1]]) #generate data frame for plotting

#plot common garden difference in leaf spines

(figs8b <- ggplot(data = plot.emm.spines.common.garden, aes(x = Source.IM, y = response, col = Source.IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = Source.IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  ylab('Spinescence (%)')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))

(leaf.spines.cg = aggregate( spines ~ Species + Source.IM + Site + plant.ID, data = cg.leaf.data, FUN = mean))

#pdf('./figures/Supplemental_Figures/FigS8B.pdf', height = 4, width = 4)

figs8b+
  geom_point(data = leaf.spines.cg[leaf.spines.cg$Species %in% c('Prunus','Heteromeles'),], 
             aes(x = Source.IM, y = spines, fill = Source.IM), 
             size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red'))

#dev.off()


###### Finally, create figure for cyanogenic glycosides

cyanide <- read.csv(file = './data_files/Shrubs/Cyanide/cyanide_measurements.csv', header = T)

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

#pdf('./figures/Supplemental_Figures/FigS7D.pdf', height = 6, width = 5)
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

model_cyanide = lmer(log(Tissue.Concentration) ~ IM*Species + IM*Age + (1|Site/PlantID), data = insitu_cyanide)

summary(model_cyanide) #No significant difference between island and mainland, but younger tissue consistently has higher concentrations

emm.cyanide.spp <- emmeans(model_cyanide, specs = ~ IM + Species, type = "response") #get estimated marginal means for each species

emm.cyanide.spp.contrasts <- contrast(object = emm.cyanide.spp, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.cyanide <- data.frame(emmeans(model_cyanide, specs = pairwise ~ IM + Species, type = "response")[[1]]) #generate data frame for plotting

(fig2d <- ggplot(data = plot.emm.cyanide, aes(x = IM, y = response, col = IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  ylab('CNglcs (mg HCN/g)')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))
#dev.off()

#pdf('./figures/Fig2d.pdf', height = 4, width = 4)
fig2d+
  geom_point(data = plot.cyanide, aes(x = IM, y = Tissue.Concentration, fill = IM), size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red'))
#dev.off()

#and finally get marginal means for canopy position for spinescence

emm.cyanide.age <- data.frame(emmeans(model_cyanide, specs = ~ IM + Species + Age, type = "response"))

emm.cyanide.age$Age <- factor(emm.cyanide.age$Age, levels = c('young','old')) #reorder so that young tissue is plotted first

#pdf('./figures/Supplemental_Figures/FigS10.pdf', height = 5, width = 6)
ggplot(data = emm.cyanide.age, aes(x = Age, y = response, col = IM, group = IM))+
  theme_bw(base_size = 14)+
  geom_point(position = position_dodge(0.1), size = 4)+
  geom_line(position = position_dodge(0.1))+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.1))+
  facet_wrap(~Species, scales = 'free_y')+
  ylab('CNglc content (mg HCN/g)')+
  xlab('Leaf Age')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  scale_color_manual(values = c('blue','red'))+
  theme(strip.text = element_text(face = "italic", size = 16))
#dev.off()

## Cyanide comparison for common garden plants

#overall model

model_cyanide_common_garden <- lmer(log(Tissue.Concentration) ~ IM*Age + IM*Species + (1|PlantID), data = cyanide[cyanide$Source %in% c('SBBG','RSA'),])

summary(model_cyanide_common_garden)
#no significant difference between island and mainland

emm.cyanide.spp.common.garden <- emmeans(model_cyanide_common_garden, specs = ~ IM + Species, type = "response") #get estimated marginal means for each species

emm.cyanide.spp.common.garden.contrasts <- contrast(object = emm.cyanide.spp.common.garden, method = "pairwise", by = "Species") #generates all species-level contrasts of interest

plot.emm.cyanide.common.garden <- data.frame(emmeans(model_cyanide_common_garden, specs = pairwise ~ IM + Species, type = "response")[[1]]) #generate data frame for plotting

(figs8c <- ggplot(data = plot.emm.cyanide.common.garden, aes(x = IM, y = response, col = IM))+
  geom_point(size = 4)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = IM), width = 0.2, size = 1)+
  facet_wrap(~Species, scales = "free_y", nrow = 1)+
  scale_color_manual(values = c('blue','red'))+
  ylab('Estimated CNglcs (mg HCN/g)')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  theme(strip.text = element_text(face = "italic", size = 16)))

ggplot(insitu_cyanide, aes(x = IM, y = Tissue.Concentration, col = Age))+
  geom_boxplot()

cg_cyanide <- cyanide[cyanide$Source %in% c('SBBG','RSA'),] #leaves 29 measurements

plot.cyanide.cg <- aggregate(Tissue.Concentration ~ Species + Site + IM + PlantID, mean, data = cg_cyanide)

#pdf('./figures/Supplemental_Figures/FigS8C.pdf', height = 4, width = 4)

figs8c+
  geom_point(data = plot.cyanide.cg, aes(x = IM, y = Tissue.Concentration, fill = IM), 
             size = 2, alpha = 0.2, pch = 21,
             position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('blue','red'))
  
#dev.off()

##############################

###### Finally, for Catalina only, compare trait data inside versus outside of exclosures

catalina <- subset(leaf.data, leaf.data$Site == 'Catalina')

with(aggregate(True_area ~ plant.ID + Species + Exclosure, mean, data = catalina)
     , table(Species, Exclosure)) #Prunus almost all outside the exclosures, but other species spread pretty evenly

#build model to see if leaf area was higher inside exclosures

catalina_model_leaf_area <- (lmer(log(True_area) ~ Exclosure*Species + Position + Aspect.E + Aspect.N + (1|plant.ID/branch.ID), data = catalina)) #df looks correct. No significant difference in leaf area inside versus outside of exclosures, though there is a slight positive effect of being inside an exclosure for leaf area.

summary(catalina_model_leaf_area) #leaf area is modestly larger inside exclosures

#next test for SLA differences

summary(lmer(SLA ~ Exclosure*Species + Position + Aspect.E + Aspect.N + (1|plant.ID), data = sla.dataframe[sla.dataframe$Site=='Catalina',])) #upper leaves have lower SLA (as expected), but again, no significant differences inside versus outside exclosures

#now for spinescence

summary(lmer(spines ~ Exclosure*Species + Position + Aspect.E + Aspect.N + (1|plant.ID/branch.ID), data = catalina[catalina$Species %in% c('Heteromeles','Prunus'),])) #lower spines on north-facing plants and branches from the upper canopy, but again no difference between inside versus outside exclosures

#lastly test CNglcs

cng.catalina <- cyanide[cyanide$Site=='Catalina' & cyanide$Source == 'Catalina',] #leaves 36 observations from 19 plants (11 Heteromeles, 8 Prunus), but only six of these (5 Heteromeles, 1 Prunus) are from inside of exclosures

cng.catalina$plant.ID <- paste(cng.catalina$Species, cng.catalina$PlantID, sep = "_")

cng.catalina.model = aggregate(Tissue.Concentration ~ IM + Age + Species + plant.ID + Exclosure, mean, data = merge(cng.catalina, catalina[,colnames(catalina) %in% c('Exclosure','plant.ID')], by = 'plant.ID', all = F))

summary(lmer(Tissue.Concentration ~ Exclosure*Species + Age*Species + (1|plant.ID), data = cng.catalina.model)) #no difference inside versus outside of exclosures, but young tissue has higher concentrations, as expected


########

##Herbivory data

hist(leaf.data$Herbivory) #strongly zero-inflated, probably makes sense to start using a negative binomial (?) distribution as a first pass

leaf.data$Herbivory.binom <- ifelse(leaf.data$Herbivory > 0, 1, 0 )

Herbivory.model <- glmer(Herbivory.binom ~ IM*Species + Aspect.N + Aspect.E + Position + (1|Site/plant.ID/branch.ID), family = 'binomial', leaf.data)

summary(Herbivory.model) #no indication that island plants are any more likely to have evidence of herbivory; only significant model term suggests that leaves on westward facing aspects are more likely to experience herbivory, plus some species-level differences in levels of herbivory

anova(Herbivory.model)

#####

#get effect size measurements for all traits for comparison

#first load Bowen and Van Vuren data

bvv97 <- data.frame(readxl::read_xlsx('./data_files/Shrubs/Bowen and Van Vuren Effect Sizes.xlsx'))

#start with leaf area

library(effectsize)

d1 <- data.frame(cohens_d(leaf.area.agg$True_area, leaf.area.agg$IM)) #in situ leaf area

d2 <- data.frame(cohens_d(leaf.area.cg.agg$True_area, leaf.area.cg.agg$Source.IM)) #common garden leaf area

#then SLA

d3 <- data.frame(cohens_d(sla.agg$SLA, sla.agg$IM)) #in situ SLA, no common garden estimates available

#then spinescence

spines.agg.cg <- aggregate(spines ~ Species + Site + Source.IM + plant.ID, data = spine.data[spine.data$Site %in% c('SBBG','RSA'),], mean)

d4 <- data.frame(cohens_d(spines.agg$spines, spines.agg$IM)) #in situ spinescence

d5 <- data.frame(cohens_d(spines.agg.cg$spines, spines.agg.cg$Source.IM)) #in situ spinescence

#and finally CNglcs

cyanide.agg <- aggregate(Tissue.Concentration ~ Species + Site + IM + PlantID, data = insitu_cyanide, mean)

d6 <- data.frame(cohens_d(cyanide.agg$Tissue.Concentration, cyanide.agg$IM)) #in situ cyanide

cyanide.agg.cg <- aggregate(Tissue.Concentration ~ Species + Site + IM + PlantID, data = cyanide[cyanide$Source %in% c('SBBG','RSA'),], mean)

d7 <- data.frame(cohens_d(cyanide.agg.cg$Tissue.Concentration, cyanide.agg.cg$IM)) #common garden cyanide

eff_size_df1 <- rbind(d1, d2, d3, d4, d5, d6, d7)

eff_size_df1$trait <- c('Leaf Area', 'Leaf Area', 'SLA', 'Spinescence', 'Spinescence', 'CNglc Content', 'CNglc Content')
eff_size_df1$comparison <- c('This study (in situ)', 'This study (botanic garden)', 'This study (in situ)', 'This study (in situ)', 'This study (botanic garden)', 'This study (in situ)', 'This study (botanic garden)')

eff_size_df1

eff_size_df <- data.frame(rbind(eff_size_df1,
                      c(1.27109515, 'SD', 0.556537944, 1.985652355, 'Leaf Area', 'Bowen and Van Vuren (1997)'), 
                      c(-2.14438919, 'SD', -2.327746765, -1.961031616, 'Spinescence', 'Bowen and Van Vuren (1997)')))
#note: values inputted above are the mean effect sizes (from each of the five measured species) from BVV and their corresponding SDs

eff_size_df$Cohens_d <- as.numeric(eff_size_df$Cohens_d)
eff_size_df$CI_low <- as.numeric(eff_size_df$CI_low)
eff_size_df$CI_high <- as.numeric(eff_size_df$CI_high)

eff_size_df$trait <- as.factor(eff_size_df$trait)
eff_size_df$trait <- factor(eff_size_df$trait, levels = c('Leaf Area', 'SLA', 'Spinescence', 'CNglc Content'))

eff_size_df$comparison <- factor(eff_size_df$comparison, levels = c('Bowen and Van Vuren (1997)', 'This study (in situ)', 'This study (botanic garden)'))

#pdf('./figures/Fig4.pdf', height = 5, width = 6)
ggplot(eff_size_df, aes(x = comparison, y = Cohens_d, col = comparison))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, size = 1)+
  facet_wrap(~trait, nrow = 1)+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'bottom',
        legend.title = element_blank(), legend.box = 'vertical', axis.ticks = element_blank())+
  geom_hline(yintercept = 0, lty = 2, size = 0.5)+
  ylab("Standardized Insularity Effect (Cohen's D)")+
  scale_color_manual(values = c('darkblue','turquoise','darkgrey'))+
  theme(legend.direction = 'vertical')
#dev.off()

#now do the same but with t-statistics, and report separately for each species; start with in situ samples

t1 <- data.frame(emm.leaf.area.spp.contrasts)
t1$Trait <- 'Leaf Area'
t1$Comparison <- 'This study (in situ)'

t2 <- data.frame(emm.sla.spp.contrasts)
t2$Trait <- 'SLA'
t2$Comparison <- 'This study (in situ)'

t3 <- data.frame(emm.spines.spp.contrast)
t3$Trait <- 'Spinescence'
t3$Comparison <- 'This study (in situ)'

t4 <- data.frame(emm.cyanide.spp.contrasts)
t4$Trait <- 'CNglc Content'
t4$Comparison <- 'This study (in situ)'

#and then do common garden samples
t5 <- data.frame(emm.leaf.area.spp.common.garden.contrasts)
t5$Trait <- 'Leaf Area'
t5$Comparison <- 'This study (common garden)'

t6 <- data.frame(emm.spines.spp.common.garden.contrast)
t6$Trait <- 'Spinescence'
t6$Comparison <- 'This study (common garden)'

t7 <- data.frame(emm.cyanide.spp.common.garden.contrasts)
t7$Trait <- 'CNglc Content'
t7$Comparison <- 'This study (common garden)'

#now combine

names(bvv97)[3] <- 'Species'
names(bvv97)[4] <- 't.ratio'

bvv97$Comparison <- 'Bowen and Van Vuren (1997)'

columns.keep <- c('Trait','Species','t.ratio','Comparison')

t_stat_combined <- rbind(t1[,names(t1) %in% columns.keep], 
                         t2[,names(t2) %in% columns.keep],
                         t3[,names(t3) %in% columns.keep],
                         t4[,names(t4) %in% columns.keep],
                         t5[,names(t5) %in% columns.keep],
                         t6[,names(t6) %in% columns.keep],
                         t7[,names(t7) %in% columns.keep],
                         bvv97[,names(bvv97) %in% columns.keep])

t_stat_combined$t.ratio <- as.numeric(t_stat_combined$t.ratio)

t_stat_combined <- t_stat_combined[!t_stat_combined$Trait %in% c('Phenols','Tannins'),]
t_stat_combined <- t_stat_combined[!(t_stat_combined$Species %in% c('Ceanothus','Dendromecon','Cercocarpus') & t_stat_combined$Trait == 'Spinescence'),]

t_stat_combined$Trait <- factor(t_stat_combined$Trait, levels = c('Leaf Area', 'SLA', 'Spinescence', 'CNglc Content'))

t_stat_combined$Comparison <- factor(t_stat_combined$Comparison, levels = c('Bowen and Van Vuren (1997)', 'This study (in situ)', 'This study (common garden)'))

#pdf('./figures/FigSxx1.pdf', height = 7, width = 7)
ggplot(t_stat_combined, aes(x = Comparison, y = t.ratio, fill = Comparison))+
  geom_hline(yintercept = 0, lty = 2, size = 0.5)+
  geom_hline(yintercept = 1.96, lty = 1, size = 0.1, color = 'red')+
  geom_hline(yintercept = -1.96, lty = 1, size = 0.1, color = 'red')+
  geom_point(size = 4, pch = 21, alpha = 0.8)+
  facet_wrap(. ~ Trait*Species, scales = 'free_x', ncol = 5)+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'bottom',
        legend.title = element_blank(), legend.box = 'vertical', axis.ticks = element_blank())+
  theme(legend.direction = 'vertical')+
  ylab("Insularity Effect (t statistic)")+
  scale_fill_manual(values = c('darkblue','turquoise','darkgrey'))+
  theme(strip.text.x = element_text(size = 10, face = "italic"))+
  ylim(c(-10.5,10.5))
#dev.off()

