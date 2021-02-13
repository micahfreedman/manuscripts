#### Script for analyzing cardenolide data

library(ggplot2)
library(dplyr)
library(vegan)
library(plotly)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(multcomp)
library(stringr)

setwd('~/Documents/GitHub/Sites/manuscripts/Cardenolide_Sequestration/')

raw_data <- read.csv('./data/concatenated_cardenolides.csv') #read in data. *note that negative values have already been replaced with 0s in this dataframe

str(raw_data)

#Description of columns: 
#Year = year in which samples were run on the HPLC (2018 vs. 2019)
#Plate = identity of the 96-well plate that sample was included in. This shouldn't matter unless there are systematic issues with loss of column integrity over the course of sampling.
#Injection = index value for the order in which samples were injected
#Sample = ID of the tissue being analyzed 
#Species = self-explanatory
#Tissue = leaf or wing; also includes entries that say standard, which correspond to runs where only the internal standard was injected
#Position = position within a given plate. Not actually imoportant.
#Digitoxin = integrated peak area of the internal standard, digitoxin
#X0.487 denotes the first peak, which had an RT of 0.487 minutes. All of the subsequent columns use the same naming system.


#First operation: subtract away background that is present in internal standard only samples. Take the mean value for each column from the varioius standard samples, and then subtract this from the rest of the dataframe. Replace any negative values with 0. Then, finally divide each row by its corresponding standard concentration.

subtract_baseline <- function(dfx) {
  standard.means <- colSums(dfx[dfx$Tissue == 'Standard', 9:ncol(dfx)], na.rm = T) / 
    nrow(dfx[dfx$Tissue=='Standard',])
  dfx2 <- sweep(dfx[,9:ncol(dfx)], 2, standard.means)
  dfx2[dfx2 < 0] <- 0
  dfx3 <- cbind(dfx[,1:8], dfx2)
}

cardenolides <- subtract_baseline(raw_data)

#Next operation: need to merge this dataframe with a separate dataframe containing the corresponding mass of the tissue sampled. Then, use this to adjust and get a concentration.

masses <- read.csv('./data/tissue_mass.csv') #read in data corresponding to the mass of the tissue sampled

str(masses) #merge based on "Sample" column, which has the same indexing in the cardenolides dataframe

cardenolides <- merge(masses[,c(1,4)], cardenolides, by = 'Sample')

#Next operation: get a cumulative cardenolide value, adjusted for tissue mass and internal standard concentration, for each sample. Here, the conversion is based on the IS concentration of 15 mg digitoxin / 100 mL. Then, for each sample, 0.8 mL out of the original 1 mL of extract was vacuum evaporated, so need a factor of 0.8 in there as well. Lastly, divide by the tissue mass/1000 to get cardenolides (mg/g), and add a column for cumulative concentration.  

standardize_cards <- function(x){
  std_vals <- (cardenolides[,10:ncol(cardenolides)] / cardenolides$Digitoxin) * 0.8 * 0.15 / 
    (cardenolides$Mass/1000)
  std_vals$concentration <- rowSums(std_vals)
  full_data <- cbind(cardenolides[,1:9], std_vals)
}

cardenolides <- standardize_cards(cardenolides) #These values correspond to cardenolides (per unit of 0.15 mg digitoxin / g tissue)

cardenolides <- cardenolides[-which(is.na(cardenolides$concentration)),] #drop the one row where the concentration of internal standard wasn't recorded, thereby producing NAs across rows

#################

#Figure out sequestration ratio: divide the average sequestered concentration for wings by the average leaf concentration. Aggregate based on species and tissue, and get both mean and standard errors.

sequestration_ratio <- function(x){
  wing_means <- do.call(data.frame, aggregate(concentration ~ Species + Tissue, 
                                              function(y) c(mean(y), sd(y) / sqrt(length(y))), data = x))
  names(wing_means)[3:4] <- c('mean','SE')
  ratio <- wing_means[wing_means$Tissue == 'Wing',]$mean / wing_means[wing_means$Tissue == 'Leaf',]$mean
  ratio.table <- data.frame(cbind(levels(factor(x$Species)), ratio))
  names(ratio.table)[1] <- 'Species'
  output <- list(wing_means, ratio.table)
  print(output)
}

sequestration_ratio(cardenolides)

######

#Plot raw concentrations

#first, rearrange milkweed species according to phylogenetic position

cardenolides$Species <- factor(cardenolides$Species, c('GOPH','ASCU','AINC','ASFA','ASYR','ASPEC'))

#make vector of colors for plotting
ascl.colors <- c('blue','purple','coral','gold','dodgerblue','turquoise')

pdf('./figures/Figx2.pdf', height = 8, width = 5)
ggplot(cardenolides, aes(x = Tissue, y = concentration, fill = Species))+
  geom_boxplot(outlier.color = 'white', alpha = 0.3)+
  geom_point(position = position_jitterdodge(1), pch = 21)+
  scale_fill_manual(values = ascl.colors)+
  facet_wrap(~Species, ncol = 2, scales = 'free')+
  theme_light(base_size = 16)+
  theme(legend.position = 'none')+
  ylab('Cardenolide Concentration (mg/g)')
dev.off()

#######

#Multivariate analysis of disparity in leaf and wing samples

#First, analyze concurrently. Need to drop all rows with cumulative values of 0. Also need to specify only numeric columns, not including the final concentration column

cards.multivar <- function(x){
  df1 <- x[-which(x$concentration==0),]
  df2 <- df1[,10:ncol(df1)-1]
}

cards.multivariate <- cards.multivar(cardenolides)

cards_nmds <- metaMDS(comm = cards.multivariate, distance = 'bray', k = 10, trymax = 20)

stressplot(cards_nmds)

cards.data.scores <- cbind(cardenolides[-which(cardenolides$concentration==0),1:9], scores(cards_nmds))

pdf('./figures/Figx3.pdf', height = 8, width = 8)
ggplot(cards.data.scores, aes(x = NMDS1, y = NMDS2, col = Species, shape = Tissue))+
  geom_point(size = 2)+
  scale_shape_manual(values = c(15,19))+
  scale_color_manual(values = ascl.colors)+
  theme_light(base_size = 16)+
  theme(legend.position = c(0.2, 0.77), legend.background = element_rect(fill = 'grey95', color = 'black'), 
        legend.title.align = 0.5, legend.title = element_text(face = 'bold'))+
  ylim(-0.8,1.2)+
  guides(col = guide_legend(ncol = 2, order = 1))
dev.off()

plot_ly(cards.data.scores[cards.data.scores$Tissue=='Wing',], x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Species, colors = ascl.colors) #3d version of the same figure, but only showing points for wings. Might have issues trying to run with R versions < 4.0.0

##### 

#Now generate a polarity index for each sample. Here, the idea is to determine whether the sequestered cardenolides from some species are disproportionately composed of polar compounds (expected to be the case for ASYR and ASPEC)

polarity.index <- function(x){
  vec.RT <- as.numeric(sub('.', '', names(x)[10:79]))
  polSum <- rowSums(data.frame(mapply(`*`,cardenolides[10:79],vec.RT)))
  pol.ind <- 1 - (polSum / max(polSum))
  print(pol.ind)
}

cardenolides$polarity.index <- polarity.index(cardenolides)

PI.plot <- do.call(data.frame, aggregate(polarity.index ~ Species, cardenolides[cardenolides$Tissue=='Wing',], function(y) c(mean(y), sd(y))))
names(PI.plot)[2:3] <- c('Polarity_Index', 'SD')

pdf('./figures/Figx3c.pdf', height = 8, width = 3)
ggplot(PI.plot, aes(x = Species, y = Polarity_Index, col = Species, fill = Species))+
  geom_point(size = 4, pch = 21)+
  geom_errorbar(aes(ymin = Polarity_Index - SD, ymax = Polarity_Index + SD), width = 0.2, size = 1)+
  theme_light(base_size = 16)+
  ylab('Polarity Index - Wing Tissue')+
  scale_color_manual(values = ascl.colors)+
  scale_fill_manual(values = ascl.colors)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

###################

## Now, combine cardenolide data with experimental information for statistical analyses

adult.info <- read.csv('./data/adult_data.csv')

head(adult.info)

names(adult.info)[1] <- 'Sample'

wing.cardenolides <- merge(adult.info[,names(adult.info) %in% c('Sample','Usage','Plant.ID','Pop','GH','Group','Mon.Pop','maternal_family','sym.allo','Sex')], cardenolides[cardenolides$Tissue == 'Wing',], by = 'Sample')

table(wing.cardenolides$maternal_family) #440 total samples from 73 maternal families

names(wing.cardenolides) #Description of columns:
#Sample: self-explanatory
#Usage: refers to whether a particular plant was being used for the first or second time
#Plant.ID: refers to the identity of the plant that was used for rearing
#Pop: provenance of the plant used for rearing
#GH: refers to either 86W or Core, the two greenhouses where rearing occurred
#Group: position of plants on greenhouse benches
#Mon.Pop: the population of origin for monarch samples
#Sex: male or female
#maternal_family: maternal family of origin for monarch samples
#sym.allo: refers to whether a given monarch was reared on a sympatric or allopatric host plant

#plot of raw data to explore patterns across species and monarch populations

#pdf('./figures/FigSx1.pdf', height = 8, width = 7)
ggplot(wing.cardenolides, aes(x = Mon.Pop, y = concentration))+
  geom_boxplot(aes(fill = Species), outlier.color = 'white', alpha = 0.3)+
  geom_point(aes(fill = Species), position = position_jitter(0.2), pch =21)+
  facet_wrap(~Species, scales = 'free')+
  scale_fill_manual(values = ascl.colors)+
  theme_light(base_size = 14)+
  xlab('Monarch Population')+
  ylab('Cardenolide Concentration (mg/g)')+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))
#dev.off()

#A couple of things to note based on above plot. AINC and ASFA can largely be ignored because the values are so consistently low across all populations. For the other species, the monarch population from Guam has fairly low sequestraiton on many of the species, including its sympatric species (ASCU). By far the most interesting outlier is the Puerto Rican population, which has the highest mean sequestration on ASCU and GOPH (despite no history with it!) and by far the lowest on ASYR and ASPEC (both involve sequestration of predominantly polar compounds). Also interesting to note that variation on GOPH seems to be substantially lower than across other species: perhaps this reflects something about sequestration on GOPH being primarily passive.

#First test: use model structure to directly test for local adaptation. Use the same model structure as in Evolution paper. In all models, exclude AINC and ASFA from the input data, since they have such low levels of cardenolide to begin with.

local.adaptation.seq.model <- lmer(concentration ~ Mon.Pop + Species + sym.allo + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])

summary(local.adaptation.seq.model) #intercept is based on Australian monarchs on ASCU. All species are lower than ASCU (checks out), and Puerto Rico and Guam both sequester less than Australia, according to this model. However, there is no effect of local adaptation: in fact, sequestration on sympatric host plants is actually modestly lower (this is probably due to Guam and ASCU being a sympatric combo)

#For the next model, don't include a sympatric/allopatric term, but do include an interaction between monarch population and milkweed species, since there did seem to be 

sequestration.model <- lmer(concentration ~  Mon.Pop*Species + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])

summary(sequestration.model) #Once again ASCU and AU are the reference levels. ASCU higher than all other species. The strongest interaction term is for Puerto Rico on ASYR, which has a strongly negative estimate. Males sequester modestly less than females (I think this is a general pattern that has been reported before).

Anova(sequestration.model, type = 3) #does indeed indicate a highly significant monarch population x species effect, although that by itself is not indicative of local adaptation

emmeans(sequestration.model, specs = ~Mon.Pop*Species, data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),]) #estimated marginal means for each species x monarch population combination

summary(glht(sequestration.model, linfct = mcp(Mon.Pop = "Tukey"))) #gives an error message related to interaction term. As an alternative, can create a new data column that coresponds to the combination of mon.pop x species, and then use this as a predictor.

wing.cardenolides$combination <- paste(wing.cardenolides$Mon.Pop, wing.cardenolides$Species, sep = ' x ')

sequestration.model.x <- lmer(concentration ~ combination + Mon.Pop + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])

summary(glht(sequestration.model.x, linfct = mcp(Mon.Pop = "Tukey"))) #Can't get estimates because model includes too many levels for the combination term (n = 36) for the model to be fully specified.

emm.mon.pop <- as.data.frame(emmeans(sequestration.model, specs = ~Mon.Pop, data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])) #marginal means for monarch population only

ggplot(emm.mon.pop, aes(x = Mon.Pop, y = emmean))+
  geom_point()+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3)

#### Now, for the sake of comparison, create a reaction norm plot showing Puerto Rico and eastern North America with ASCU and ASYR. Build a function that enables pairwise reaction norms comparisons for any combination of two species and two populations.

rxn.norm.plot <- function(df1, pops, spp){
  df.use <- df1[df1$Species %in% spp & df1$Mon.Pop %in% pops,]
  agg.data <- do.call(data.frame, aggregate(concentration ~ Mon.Pop + Species, df.use, 
                        function(y) c(mean(y), sd(y) / sqrt(length(y)))))
  names(agg.data)[3:4] <- c('Concentration', 'SE')
  ggplot(agg.data, aes(x = Species, y = Concentration, group = Mon.Pop, col = Mon.Pop))+
    geom_point(size = 3)+
    geom_line(size = 1)+
    geom_errorbar(aes(ymin = Concentration - SE, ymax = Concentration + SE), width = 0.1, size = 1)+
    theme_classic(base_size = 20)+
    scale_color_manual(values = c('darkgreen','orange'))+
    ylab('Wing Cardenolide Concentration (mg/g)')+
    theme(legend.position = c(0.85, 0.76), legend.title = element_blank())
}

pdf('./figures/Figx4.pdf', height = 6, width = 4)
rxn.norm.plot(wing.cardenolides, pops = c('ENA','PR'), spp = c('ASYR','ASCU'))
dev.off()

## Related to the above figure, conduct an analysis to see whether there is a potential tradeoff in sequestration ability across temperate versus tropical species (here, this is a standin for sequestration mode, as the temperate species require more active sequestration)

wing.cardenolides$polar.spp <- ifelse(wing.cardenolides$Species %in% c('ASCU','GOPH'), 'tropical', 'temperate')

trop.model.seq <- lmer(concentration ~ polar.spp*Mon.Pop + (1|Species/Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('AINC','ASFA'),])

summary(trop.model.seq) #Puerto Rico struggles mightily on temperate species, not surprisingly. Guam is actually modestly better on temperate species than would be suggested, but this effect would go away after correcting for multiple comparisons.
Anova(trop.model.seq, type = 3) #again, reiterates that populations vary greatly in their ability to sequester depending on whether the host plant is tropical or temperate in origin

######## Q for marianas comps: do ASCU genotypes from Rota really have lower cardenolides than from Guam

cardenolides[cardenolides$Species=='ASCU' & cardenolides$Tissue=='Leaf',]$Sample

merge(cardenolides[cardenolides$Species=='ASCU' & cardenolides$Tissue=='Leaf',], adult.info[adult.info$Species=='ASCU', names(adult.info) %in% c('Sample','Pop')], by = 'Sample')

adult.info[adult.info$Species=='ASCU',names(adult.info) %in% c('Sample','Pop')]

ascu.leaf.cards <- cardenolides[cardenolides$Species=='ASCU' & cardenolides$Tissue=='Leaf',]
ascu.leaf.cards$Sample
ascu.leaf.cards$Pop <- ifelse(ascu.leaf.cards$Sample %in% c('1004','1005','1014','1017','1018','1020','1024','1025','1027','1028','1030','1034','1037','1042','1059','1062','1057','172','186','195','203','216','218','224','L170','L175','L188','L190','L199','L206','L215'), 'Guam', 'Rota')

summary(lm(concentration ~ Pop, data = ascu.leaf.cards)) #not significantly higher, but def in the same direction as field samples

ggplot(ascu.leaf.cards, aes(x = Pop, y = concentration))+
  geom_boxplot()



####### Unfortunately can not include the 5 monarchs collected in Guam in 2018, since they were analyzed in Rachel's lab on a different column


