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

#make vector of colors for plotting
ascl.colors <- c('coral','purple','gold','turquoise','dodgerblue','blue')

pdf('./figures/Figx2.pdf', height = 8, width = 5)
ggplot(cardenolides, aes(x = Tissue, y = concentration, fill = Species))+
  geom_boxplot(outlier.color = 'white', alpha = 0.3)+
  geom_point(position = position_jitterdodge(0.5), pch = 21)+
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

pdf('./figures/Figx3.pdf', height = 8, width = 10)
ggplot(cards.data.scores, aes(x = NMDS1, y = NMDS2, col = Species))+
  geom_point(aes(shape = Tissue), size = 2)+
  scale_shape_manual(values = c(15,19))+
  scale_color_manual(values = ascl.colors)+
  theme_light(base_size = 16)
dev.off()

plot_ly(cards.data.scores[cards.data.scores$Tissue=='Wing',], x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Species, colors = ascl.colors) #3d version of the same figure, but only showing points for wings. Might have issues trying to run with R versions < 4.0.0

##### 

#Now generate a polarity index for each sample. Here, the idea is to determine whether the sequestered cardenolides from some species are disproportionately composed of polar compounds (expected to be the case for ASYR and ASPEC)



###################

## Now, combine cardenolide data with experimental information for statistical analyses

adult.info <- read.csv('./data/adult_data.csv')

head(adult.info)

names(adult.info)[1] <- 'Sample'

wing.cardenolides <- merge(adult.info[,names(adult.info) %in% c('Sample','Usage','Plant.ID','Pop','GH','Group','Mon.Pop','maternal_family','sym.allo')], cardenolides[cardenolides$Tissue == 'Wing',], by = 'Sample')

table(wing.cardenolides$maternal_family) #440 total samples from 73 maternal families

names(wing.cardenolides) #Description of columns:
#Sample: self-explanatory
#Usage: refers to whether a particular plant was being used for the first or second time
#Plant.ID: refers to the identity of the plant that was used for rearing
#Pop: provenance of the plant used for rearing
#GH: refers to either 86W or Core, the two greenhouses where rearing occurred
#Group: position of plants on greenhouse benches
#Mon.Pop: the population of origin for monarch samples
#maternal_family: maternal family of origin for monarch samples
#sym.allo: refers to whether a given monarch was reared on a sympatric or allopatric host plant

sequestration.model <- lmer(concentration ~  Mon.Pop*Species + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])

summary(sequestration.model)

Anova(sequestration.model, type = 3) #does indeed indicate a highly significant monarch population x species effect, although that by itself is not indicative of local adaptation

emmeans(sequestration.model, specs = ~Mon.Pop*Species, data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])

summary(glht(sequestration.model, linfct = mcp(Mon.Pop = "Tukey")))



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

