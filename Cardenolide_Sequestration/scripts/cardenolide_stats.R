######################################################################################
####################### Script for analyzing cardenolide data ########################
######################################################################################

#Note that analyses presented here do not necessarily correspond to the order in which they appear within the manuscript. To locate code used to generate a particular figure or table, use CTRL-F (find) to search it by name -- the naming convention is Figure1, Figure2, FigureS1, TableS7, etc. Code for some figures may be located in other scripts.

#load libraries

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
library(corrplot)
library(equatiomatic)

#setwd() - specify root directory here; this should be the file downloaded from Dryad

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

######################################################################################
######################################################################################
######################################################################################

#First operation: subtract away background that is present in internal standard only samples. Take the mean value for each column from the varioius standard samples, and then subtract this from the rest of the dataframe. Replace any negative values with 0. Then, finally divide each row by its corresponding standard concentration.

subtract_baseline <- function(dfx) {
  standard.means <- colSums(dfx[dfx$Tissue == 'Standard', 9:ncol(dfx)], na.rm = T) / 
    nrow(dfx[dfx$Tissue=='Standard',])
  dfx2 <- sweep(dfx[,9:ncol(dfx)], 2, standard.means)
  dfx2[dfx2 < 0] <- 0
  dfx3 <- cbind(dfx[,1:8], dfx2)
}

cardenolides <- subtract_baseline(raw_data)

######################################################################################
######################################################################################
######################################################################################

#Next operation: need to merge this dataframe with a separate dataframe containing the corresponding mass of the tissue sampled. Then, use this to adjust and get a concentration.

masses <- read.csv('./data/tissue_mass.csv') #read in data corresponding to the mass of the tissue sampled

str(masses) #merge based on "Sample" column, which has the same indexing in the cardenolides dataframe

cardenolides <- merge(masses[,c(1,4)], cardenolides, by = 'Sample')

######################################################################################
######################################################################################
######################################################################################

#Next operation: get a cumulative cardenolide value, adjusted for tissue mass and internal standard concentration, for each sample. Here, the conversion is based on the IS concentration of 15 mg digitoxin / 100 mL. Then, for each sample, 0.8 mL out of the original 1 mL of extract was vacuum evaporated, so need a factor of 0.8 in there as well. Lastly, divide by the tissue mass/1000 to get cardenolides (mg/g), and add a column for cumulative concentration.  

standardize_cards <- function(x){
  std_vals <- (cardenolides[,10:ncol(cardenolides)] / cardenolides$Digitoxin) * 0.8 * 0.15 / 
    (cardenolides$Mass/1000)
  std_vals$concentration <- rowSums(std_vals)
  full_data <- cbind(cardenolides[,1:9], std_vals)
}

cardenolides <- standardize_cards(cardenolides) #These values correspond to cardenolides (per unit of 0.15 mg digitoxin / g tissue)

cardenolides <- cardenolides[-which(is.na(cardenolides$concentration)),] #drop the one row where the concentration of internal standard wasn't recorded, thereby producing NAs across rows

######################################################################################
######################################################################################
######################################################################################

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

sequestration_ratio(cardenolides) #Note: these are "aggregated" values and reflect mean wing divided by mean plant concentrations

######################################################################################
######################################################################################
######################################################################################

#Plot raw concentrations

#first, rearrange milkweed species according to phylogenetic position (inferred based on Fishbein et al. papers)

cardenolides$Species <- factor(cardenolides$Species, c('GOPH','ASCU','AINC','ASFA','ASYR','ASPEC'))

#make vector of colors for plotting
ascl.colors <- c('blue','purple','coral','gold','dodgerblue','turquoise')

#as per reviewer request, use full species names instead of abbreviations

cardenolides$Species <- recode_factor(cardenolides$Species, GOPH = 'G. physocarpus',
                               ASCU = 'A. curassavica',
                               AINC = 'A. incarnata',
                               ASFA = 'A. fascicularis',
                               ASYR = 'A. syriaca',
                               ASPEC = 'A. speciosa')

#Make Figure 3 - within-species comparisons of leaf and wing concentrations

#pdf('./figures/Figure3.pdf', height = 6, width = 6)
ggplot(cardenolides, aes(x = Tissue, y = concentration, fill = Species))+
  geom_boxplot(outlier.color = 'white', alpha = 0.3)+
  geom_point(position = position_jitterdodge(1), pch = 21)+
  scale_fill_manual(values = ascl.colors)+
  facet_wrap(~Species, ncol = 3, scales = 'free')+
  theme_light(base_size = 16)+
  theme(legend.position = 'none')+
  ylab('Cardenolide Concentration (mg/g)')+
  theme(strip.text = element_text(face = 'italic', size = 14))
#dev.off()

#separately plot only species used in primary analysis (use this plot for presentation at Evolution meeting). Here, A. incarnata and A. fascicularis are omitted because of their consistently low levels of sequestered cardenolides.
ggplot(cardenolides[!cardenolides$Species %in% c('A. incarnata','A. fascicularis'),], aes(x = Tissue, y = concentration, fill = Species))+
  geom_boxplot(outlier.color = 'white', alpha = 0.3)+
  geom_point(position = position_jitterdodge(1), pch = 21, size = 3)+
  scale_fill_manual(values = ascl.colors[c(1,2,5,6)])+
  facet_wrap(~Species, nrow = 1)+
  theme_light(base_size = 20)+
  theme(legend.position = 'none')+
  ylab('Cardenolide Concentration (mg/g)')+
  theme(strip.text = element_text(face = 'italic', size = 20))

######################################################################################
######################################################################################
######################################################################################

#Multivariate analysis of disparity in leaf and wing samples

#First, analyze concurrently. Need to drop all rows with cumulative values of 0. Also need to specify only numeric columns, not including the final concentration column

cards.multivar <- function(x){
  df1 <- x[-which(x$concentration==0),]
  df2 <- df1[,10:ncol(df1)-1]
}

cards.multivariate <- cards.multivar(cardenolides)

cards_nmds <- metaMDS(comm = cards.multivariate, distance = 'bray', k = 10, trymax = 20) #run with k = 10 dimensions and 20 random starts. General rule of thumb is that stress values should be below 0.1. This takes about 90 seconds to run on my computer.

stressplot(cards_nmds)

cards.data.scores <- cbind(cardenolides[-which(cardenolides$concentration==0),1:9], cards_nmds$points) #merge sample data with MDS scores

#drop one GOPH wing sample that appears to be an indexing error (groups with AINC and ASFA)

cards.data.scores <- cards.data.scores[!(cards.data.scores$Species=='GOPH' & cards.data.scores$MDS1 > 0.3),]

#Make figure 2B - NMDS plot of all leaf and wing samples shown concurrently

#pdf('./figures/Figure2B.pdf', height = 8, width = 8)
ggplot(cards.data.scores, aes(x = MDS1, y = MDS2, col = Species, shape = Tissue))+
  geom_point(size = 2)+
  scale_shape_manual(values = c(21,15))+
  scale_color_manual(values = ascl.colors)+
  theme_light(base_size = 16)+
  theme(legend.position = c(0.25, 0.82), legend.background = element_rect(fill = 'grey95', color = 'black'), 
        legend.title.align = 0.5, legend.title = element_text(face = 'bold'))+
  theme(legend.text = element_text(face = 'italic'))+
  ylim(-0.8,1.2)+
  guides(col = guide_legend(ncol = 2, order = 1))
#dev.off()

plot_ly(cards.data.scores[cards.data.scores$Tissue=='Wing',], x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~Species, colors = ascl.colors, type = 'scatter3d', mode = 'markers') #3d version of the same figure, but only showing points for wings. Might have issues trying to run with R versions < 4.0.0

######################################################################################
######################################################################################
######################################################################################

#Within each milkweed species, test for whether leaf and wing cardenolide composition is different
cards.permanova <- cardenolides[-which(rowSums(cardenolides[,10:79]) == 0),] #select only rows with non-zero cumulative concentrations

lw.permanova <- function(x){
  output <- list()
    for(i in 1:length(unique(x$Species))){
    sppX <- x[x$Species == levels(x$Species)[i],]
    cards.dist.mat <- vegdist(sppX[,10:79], method = 'bray')
    output[[i]] <- adonis2(cards.dist.mat ~ Tissue, data = sppX, permutations = 1000)
    print(levels(x$Species)[i])
    print(output[[i]])}
}

lw.permanova(cards.permanova) #print output of six MANOVAs, arranged based on species identity; order is same as that produced using levels(cards.permanova$Species), i.e. GOPH, ASCU, AINC, ASFA, ASYR, ASPEC

#Include these results as part of Table S4

#GOPH: F = 103.46
#ASCU: F = 57.19
#AINC: F = 15.07
#ASYR: F = 26.96
#ASFA: F = 20.11
#ASPEC: F = 9.46

######################################################################################
######################################################################################
######################################################################################

#Now generate a polarity index for each sample. Here, the idea is to determine whether the sequestered cardenolides from some species are disproportionately composed of polar compounds (expected to be the case for A. syriaca and A. speciosa)

#drop rows where sum is 0 from main dataframe
cardenolides.polarity <- cardenolides[-which(rowSums(cardenolides[,10:79])<0.1),]

polarity.index <- function(dfx){
  cards.numeric <- dfx[,10:79] #can't include samples with sum 0 rows
  rel.abund <- sweep(cards.numeric, 1, rowSums(cards.numeric), "/")
  vec.RT <- as.numeric(sub('.', '', names(dfx)[10:79]))
  polSum <- rowSums(data.frame(mapply(`*`,rel.abund,vec.RT)))
  pol.ind <- 1 - (polSum / max(polSum))
  print(pol.ind)
}

cardenolides.polarity$polarity <- polarity.index(cardenolides.polarity)

PI.plot <- do.call(data.frame, aggregate(polarity ~ Species, cardenolides.polarity[cardenolides.polarity$Tissue=='Wing',], function(y) c(mean(y), sd(y))))
names(PI.plot)[2:3] <- c('polarity', 'SD') #generate means and standard deviations for polarity index values across milkweed species for wing tissue samples

#Make Figure S2-A

#pdf('./figures/FigureS2A.pdf', height = 8, width = 3)
ggplot(PI.plot, aes(x = Species, y = polarity, col = Species, fill = Species))+
  geom_point(size = 4, pch = 21)+
  geom_errorbar(aes(ymin = polarity - SD, ymax = polarity + SD), width = 0.2, size = 1)+
  theme_light(base_size = 16)+
  ylab('Polarity Index - Wing Tissue')+
  scale_color_manual(values = ascl.colors)+
  scale_fill_manual(values = ascl.colors)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, face = 'italic'))
#dev.off()

######################################################################################
######################################################################################
######################################################################################

## Now, combine cardenolide data with experimental information for statistical analyses

adult.info <- read.csv('./data/adult_data.csv')

head(adult.info)

names(adult.info)[1] <- 'Sample'

#First, get info for leaf samples

cardenolides$Plant.ID <- cardenolides$Sample

leaves <- merge(cardenolides[cardenolides$Tissue=='Leaf',], unique(adult.info[,names(adult.info) %in% c('Plant.ID','Pop')]), by = 'Plant.ID')

#As per reviewer request, generate means and standard errors of leaf cardenolide content for each milkweed population

milkweed.pops <- do.call(data.frame, aggregate(concentration ~ Pop + Species, data = leaves, function(x) c(length(x), mean(x), sd(x)/sqrt(length(x)))))
names(milkweed.pops)[3:5] <- c('N', 'Mean', 'SE')

#Include as part of contents of Table S1
milkweed.pops #entries with NA for the standard error had sample sizes of 1

#######
#Now get info for wing samples

wing.cardenolides <- merge(adult.info[,names(adult.info) %in% c('Sample','Usage','Pop','GH','Group','Mon.Pop','maternal_family','sym.allo','Sex','Plant.ID')], cardenolides[cardenolides$Tissue == 'Wing', 1:ncol(cardenolides)-1], by = 'Sample') #omit the plant.ID from the cardenolides dataframe to avoid issues with duplicate indexing

monarch.wings.byPop <- do.call(data.frame, aggregate(concentration ~ Pop + Species, data = wing.cardenolides, function(x) c(length(x), mean(x), sd(x)/sqrt(length(x)))))
names(monarch.wings.byPop)[3:5] <- c('N', 'Mean', 'SE')

#Remaining contents of Table S1
monarch.wings.byPop

ps1 <- rbind(milkweed.pops, monarch.wings.byPop) #combine datasets for plotting, add new column for tissue type
ps1$tissue <- c(rep('leaf', nrow(milkweed.pops)), rep('wing', nrow(monarch.wings.byPop)))

#Plot to show relationship between individual milkweed genotypes and the level of sequestered wing cardenolides from monarchs on those genotypes. Include as Figure S11.

#pdf('./figures/FigureS10.pdf', height = 10, width = 10)
ggplot(ps1, aes(x = Pop, y = Mean, col = Species, shape = tissue))+
  theme_bw(base_size = 16)+
  geom_point()+
  geom_errorbar(aes(ymax = Mean + SE, ymin = Mean - SE), width = 0.2)+
  facet_wrap(~Species, scales = 'free')+
  scale_color_manual(values = ascl.colors)+
  ylab('Cardenolide Concentration')+
  xlab('Milkweed Population')+
  theme(axis.text.x = element_text(angle = 90, vjust = 1))
#dev.off()

######################################################################################
######################################################################################
######################################################################################

table(wing.cardenolides$maternal_family)

nrow(wing.cardenolides)

length(unique(wing.cardenolides$maternal_family)) #440 total samples from 73 maternal families

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

#### For Table S2, get number of replicates for each maternal family across each species

write.csv(table(wing.cardenolides$maternal_family, wing.cardenolides$Species), file = './figures/TableS2.csv')

#as per reviewer request, rename populations so that they aren't just abbreviations

wing.cardenolides$Mon.Pop <- recode_factor(wing.cardenolides$Mon.Pop, AU = 'Australia',
                                      CA = 'W. N. America',
                                      ENA = 'E. N. America',
                                      GU = 'Guam',
                                      HI = 'Hawaii',
                                      PR = 'Puerto Rico')

wing.cardenolides$Mon.Pop <- factor(wing.cardenolides$Mon.Pop, levels = c('E. N. America', 'W. N. America', 'Hawaii', 'Guam', 'Australia', 'Puerto Rico')) #rearrange order of populations so that they reflect patterns of relatedness

#plot of raw data to explore patterns across species and monarch populations
ggplot(wing.cardenolides, aes(x = Mon.Pop, y = concentration))+
  geom_boxplot(aes(fill = Species), outlier.color = 'white', alpha = 0.3)+
  geom_point(aes(fill = Species), position = position_jitter(0.25), pch =21)+
  facet_wrap(~Species, scales = 'free')+
  scale_fill_manual(values = ascl.colors)+
  theme_light(base_size = 16)+
  xlab('Monarch Population')+
  ylab('Wing Cardenolide Concentration (mg/g)')+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))

#A couple of things to note based on above plot. AINC and ASFA can largely be ignored because the values are so consistently low across all populations. For the other species, the monarch population from Guam has fairly low sequestraiton on many of the species, including its sympatric species (ASCU). By far the most interesting outlier is the Puerto Rican population, which has the highest mean sequestration on ASCU and GOPH (despite no history with it!) and by far the lowest on ASYR and ASPEC (both involve sequestration of predominantly polar compounds). Also interesting to note that variation on GOPH seems to be substantially lower than across other species: perhaps this reflects something about sequestration on GOPH being primarily passive.

#First test: use model structure to directly test for local adaptation. Use the same model structure as in Freedman et al. (2020) Evolution paper. In all models, exclude AINC and ASFA from the input data, since they have such low levels of cardenolide to begin with. This is the second model form described in the manuscript.

local.adaptation.seq.model <- lmer(concentration ~  Mon.Pop + Species + sym.allo + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),])

summary(local.adaptation.seq.model) #intercept is based on Australian monarchs on GOPH. All species are lower than ASCU (checks out), and Guam sequesters less than Australia, according to this model. However, there is no effect of local adaptation, with only a very modest positive intercept for sympatric status.
ranef(local.adaptation.seq.model) 

Anova(local.adaptation.seq.model, type = 3) #These are the results presented in Table S8

#As per reviewer suggestion, present equation for the above model in non-R syntax

equatiomatic::extract_eq(local.adaptation.seq.model) #Note that output here was then fed into a LaTeX rendering software (https://quicklatex.com) and truncated because of the large number of model terms

#For the next model, don't include a sympatric/allopatric term, but do include an interaction between monarch population and milkweed species, since there did seem to be one. This is the first model form described in the manuscript.

sequestration.model <- lmer(concentration ~  Mon.Pop*Species + Sex + (1|Pop/Plant.ID) + (1|maternal_family) , data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),])

summary(sequestration.model) #Once again ASCU and AU are the reference levels. ASCU higher than all other species. The strongest interaction term is for Puerto Rico on ASYR, which has a strongly negative estimate. Males sequester modestly less than females (I think this is a general pattern that has been reported before).

ranef(sequestration.model)
hist(ranef(sequestration.model)[[2]][[1]], breaks = 9) #distribution of random intercepts for maternal families

Anova(sequestration.model, type = 3) #These are the results presented in Table S8

#as per suggestion of reviewers, express this model in formal linear model syntax. Once again, output was truncated for brevity.

equatiomatic::extract_eq(sequestration.model)

######################################################################################
######################################################################################
######################################################################################

#As per suggestions of reviewers, reformulate model so that it includes plant population explicitly nested within milkweed species, and milkweed species and monarch population as random effects themselves

sequestration.model2 <- lmer(concentration ~ Sex + (1|Mon.Pop/maternal_family) + (1|Species/Pop/Plant.ID), data = wing.cardenolides[!wing.cardenolides %in% c('A. fascicularis','A. incarnata'),]) 

summary(sequestration.model2)

coef(sequestration.model2)
#gives warning about singularity: zero variance among maternal families within this model. In practice, this results because of the very large number of observation levels (e.g. 355 for Plant.ID nested within Pop nested within Species) compared to the relatively small number of total observations (439)

sequestration.model3 <- lmer(concentration ~ Sex + (1|Mon.Pop) + (1|Species/Pop/Plant.ID), data = wing.cardenolides[!wing.cardenolides %in% c('A. fascicularis','A. incarnata'),]) #this model runs without problems, although does not include the term for maternal family and also does not accommodate an interaction between monarch population and milkweed species

sequestration.model4 <- lmer(concentration ~ Sex + (Mon.Pop|Species) + (1|maternal_family) + (1|Pop/Plant.ID), data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),]) #also runs, but cannot accommodate the full range of random effects (plant population, maternal family, plant ID) without singularity errors

summary(sequestration.model4)

hh = merTools::REsim(sequestration.model4, n.sims = 10000) #runs without throwing an error message, but note that variance estimates are basically 0 (reflecting the warning about singularity that the initial model call prodoced) for many of the milkweed*monarch.pop combinations

tail(hh, n = 24) #note: overall model intercept was 5.56

#as per reviewer suggestion, also test sym.allo model where A. curassavica is treated as sympatric for ENA and CA

wing.cardenolides$sym.allo2 <- with(wing.cardenolides, ifelse(Mon.Pop %in% c('E. N. America','W. N. America','Guam','Puerto Rico') & Species=='A. curassavica','sympatric', ifelse(Mon.Pop=='E. N. America' & Species %in% c('A. syriaca','A. incarnata'), 'sympatric', ifelse(Mon.Pop=='W. N. America' & Species %in% c('A. speciosa','A. fascicularis'), 'sympatric', ifelse(Mon.Pop %in% c('Hawaii','Australia') & Species=='G. physocarpus', 'sympatric', 'allopatric')))))

local.adaptation.seq.model2 <- lmer(concentration ~ Mon.Pop + Species + sym.allo2 + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),])

summary(local.adaptation.seq.model2) #does not change inferences about sym.allo term

coef(sequestration.model4)
#compared to model that treated milkweed species and monarch population as fixed effects, the values are generally similar, but the random effects structure tends to reduce variation within each milkweed species due to shrinkage

######################################################################################
######################################################################################
######################################################################################

#### For Figure S9, show variation among maternal families in sequestration amounts. Use random intercepts.

randoms<-ranef(sequestration.model, condVar = TRUE)
mat.fams <- attr(ranef(sequestration.model, condVar = TRUE)[[2]], "postVar") #maternal families
mw.pops <- attr(ranef(sequestration.model, condVar = TRUE)[[3]], "postVar") #milkweed populations

rand.interc1<-randoms$maternal_family
rand.interc2<-randoms$Pop

df1<-data.frame(Intercepts=randoms$maternal_family[,1],
               sd.interc=2*sqrt(mat.fams[,,1:length(mat.fams)]),
               lev.names=rownames(rand.interc1))
df2<-data.frame(Intercepts=randoms$Pop[,1],
                sd.interc=2*sqrt(mw.pops[,,1:length(mw.pops)]),
                lev.names=rownames(rand.interc2))

df1$lev.names<-factor(df1$lev.names,levels=df1$lev.names[order(df1$Intercepts)])
df1$pop <- sub("\\_.*", "", df1$lev.names)

df1 <- df1[order(df1$pop),]

df2$lev.names<-factor(df2$lev.names,levels=df2$lev.names[order(df2$Intercepts)])
df2$species <- ifelse(df2$lev.names %in% c('1','2','3','4','5','6','10','11','N1','N2','D','S'), 'GOPH', ifelse(df2$lev.names %in% c('C1','C13','C15','C2','C21','C39','C46','ITH1','ITH2'), 'ASYR', ifelse(df2$lev.names %in% c('G','GS','R'), 'ASCU', 'ASPEC')))


#pdf('./figures/FigureS9.pdf', height = 12, width = 8)
ggplot(df1,aes(lev.names,Intercepts, col = pop))+
  geom_hline(yintercept=0) +geom_errorbar(aes(ymin=Intercepts-sd.interc, 
       ymax=Intercepts+sd.interc), width=0,color="black") + geom_point(size = 1)+
  theme_bw() + xlab("Monarch Population") + ylab("")+
  theme(axis.text.x=element_text(size=rel(1.2)),
        axis.title.x=element_text(size = 12),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_blank())+ 
  scale_color_manual(values = c('cyan','limegreen','forestgreen','purple','blue','orange'))+
  coord_flip()+
  ylab('Random Intercept')+
  facet_wrap(~pop, scales = 'free_y', ncol = 1)+
  theme(legend.position = 'NONE')
#dev.off()

#Figure S10 - show random intercepts for each plant population within each species

#pdf('./figures/FigureS10.pdf', height = 12, width = 8)
ggplot(df2, aes(lev.names,Intercepts, col = species))+
  geom_hline(yintercept=0) +geom_errorbar(aes(ymin=Intercepts-sd.interc, 
       ymax=Intercepts+sd.interc), width=0,color="black") + geom_point(size = 1)+
  theme_bw() + xlab("Milkweed Population") + ylab("")+
  theme(axis.text.x=element_text(size=rel(1.2)),
        axis.title.x=element_text(size = 12),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_blank())+ 
  scale_color_manual(values = c('purple','turquoise','dodgerblue','blue'))+
  coord_flip()+
  ylab('Random Intercept')+
  facet_wrap(~species, scales = 'free_y', ncol = 1)+
  theme(legend.position = 'NONE')
#dev.off()

######################################################################################
######################################################################################
######################################################################################

#Create plots showing outputs (estimated marginal means) from statistical models

emm.combos <- as.data.frame(emmeans(sequestration.model, specs = ~Mon.Pop*Species, data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),])) #estimated marginal means for each species x monarch population combination

write.csv(emm.combos, file = './figures/TableS11.csv', row.names = FALSE) #write the output of the above call to Table S11

######

#Use model output to create Figure 5A, showing marginal means of overall sequestration across species for each population

emm.populations <- as.data.frame(emmeans(sequestration.model, specs = ~Mon.Pop, data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),]))

emm.populations

#pdf('./figures/Figure5A.pdf', height = 8, width = 6)
ggplot(emm.populations, aes(x = Mon.Pop, y = emmean))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, size = 1)+
  theme_light(base_size = 24)+
  ylab('Estimated Marginal Mean - \nCardenolide Concentration (mg/g)')+
  xlab('Monarch Population')+
  ggtitle(label = 'Overall sequestration')+
  theme(legend.title = element_blank(), legend.position = 'bottom')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  annotate("text", x = 1, y = 8.2, label = "A", size = 6)+
  annotate("text", x = 2, y = 8.5, label = "A", size = 6)+
  annotate("text", x = 3, y = 8.7, label = "A", size = 6)+
  annotate("text", x = 4, y = 7.5, label = "A", size = 6)+
  annotate("text", x = 5, y = 9, label = "A", size = 6)+
  annotate("text", x = 6, y = 7.7, label = "A", size = 6)
#dev.off()

#plot values for ASCU from the above

emm.combos$sym.allo <- ifelse(emm.combos$Mon.Pop %in% c('Puerto Rico','Guam') & 
                                emm.combos$Species == 'A. curassavica', 'sympatric',
                              ifelse(emm.combos$Mon.Pop %in% c('Hawaii','Australia') & 
                                       emm.combos$Species=='G. physocarpus', 'sympatric',
                                     ifelse(emm.combos$Mon.Pop=='E. N. America' & emm.combos$Species == 'A. syriaca', 'sympatric', ifelse(emm.combos$Mon.Pop=='W. N. America' & emm.combos$Species=='A. speciosa', 'sympatric', 'allopatric'))))

#Create Figure 5B, this time focusing only on marginal means from A. curassavica. Here the goal is to highlight that Guam is lower than all other populations on this host.

#pdf('./figures/Figure5B.pdf', height = 8, width = 6)
ggplot(emm.combos[emm.combos$Species=='A. curassavica',], aes(x = Mon.Pop, y = emmean))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, size = 1)+
  theme_light(base_size = 24)+
  ylab('Estimated Marginal Mean - \nCardenolide Concentration (mg/g)')+
  xlab('Monarch Population')+
  ggtitle(label = 'A. curassavica')+
  #scale_color_manual(values = c('black','red'))+
  theme(plot.title = element_text(hjust = 0.5, face = 'italic'))+
  annotate("text", x = 5, y = 18, label = "A", size = 6)+
  annotate("text", x = 2, y = 16, label = "AB", size = 6)+
  annotate("text", x = 1, y = 16, label = "A", size = 6)+
  annotate("text", x = 4, y = 11.5, label = "B", size = 6)+
  annotate("text", x = 3, y = 14, label = "AB", size = 6)+
  annotate("text", x = 6, y = 18, label = "A", size = 6)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))
#dev.off()

#for purposes of display, truncate values for lower confidence limits at 0

emm.combos$lower.CL <- ifelse(emm.combos$lower.CL < 0, 0, emm.combos$lower.CL)

#Previously Figure 4A: marginal means across each host x population combination, color coded based on sympatric/allopatric status

#pdf('./figures/Figure4_old.pdf', height = 12, width = 6)
ggplot(emm.combos, aes(x = Mon.Pop, y = emmean, col = sym.allo))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, size = 1)+
  theme_light(base_size = 16)+
  ylab('Estimated Marginal Mean - \nCardenolide Concentration (mg/g)')+
  xlab('Monarch Population')+
  scale_color_manual(values = c('black','red'))+
  facet_wrap(~Species, scales = 'free', nrow = 2)+
  theme(legend.title = element_blank(), legend.position = 'bottom')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  theme(strip.text = element_text(face = 'italic'))
#dev.off()

#Add updated Figure 4a, which shows raw data using boxplots for each combination.

#pdf('./figures/Figure4a.pdf', height = 8, width = 5)
ggplot(wing.cardenolides[!(wing.cardenolides$Species %in% c('A. incarnata','A. fascicularis')),], 
       aes(x = Mon.Pop, y = concentration))+
  geom_boxplot(aes(col = sym.allo), outlier.color = 'white')+
  geom_point(aes(fill = sym.allo), position = position_jitter(0.25), size = 2, pch = 21, alpha = 0.5)+
  theme_light(base_size = 16)+
  scale_color_manual(values = c('black','red'))+
  scale_fill_manual(values = c('black','red'))+
  facet_wrap(~Species, scales = 'free', nrow = 2)+
  theme(legend.title = element_blank(), legend.position = 'bottom')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  theme(strip.text = element_text(face = 'italic'))+
  xlab('Monarch Population')+
  ylab('Hindwing Cardenolide Concentration (mg/g)')
#dev.off()

#Add updated Figure 4B, which shows marginal means for the sympatric versus allopatric contrast

emm.sym.allo <- as.data.frame(emmeans(local.adaptation.seq.model, specs = ~sym.allo, data = wing.cardenolides[!wing.cardenolides$Species %in% c('A. fascicularis','A. incarnata'),]))

#pdf('./figures/Figure4b.pdf', height = 8, width = 2.5)
ggplot(emm.sym.allo, aes(x = sym.allo, y = emmean, col = sym.allo))+
  geom_point(size = 6)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, size = 1)+
  theme_light(base_size = 16)+
  ylab('Estimated Marginal Mean - \nCardenolide Concentration (mg/g)')+
  theme(legend.title = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  scale_color_manual(values = c('black','red'))
#dev.off()  

#plot for presentation at Evolution meeting, with no sympatric/allopatric contrast shown
ggplot(emm.combos, aes(x = Mon.Pop, y = emmean))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, size = 1)+
  theme_light(base_size = 20)+
  ylab('Estimated Marginal Mean - \nCardenolide Concentration (mg/g)')+
  xlab('Monarch Population')+
  facet_wrap(~Species, nrow = 1)+
  theme(legend.title = element_blank(), legend.position = 'bottom')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))+
  theme(strip.text = element_text(size = 20, face = 'italic'))

summary(glht(sequestration.model, linfct = mcp(Mon.Pop = ""))) #Explicitly test for population-level differences in sequestration. Gives a warning message related to interaction term. As an alternative, can create a new data column that corresponds to the combination of mon.pop x species, and then use this as a predictor.

wing.cardenolides$combination <- paste(wing.cardenolides$Mon.Pop, wing.cardenolides$Species, sep = ' x ')

sequestration.model.combo <- lmer(concentration ~ combination + Mon.Pop + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = wing.cardenolides[!wing.cardenolides$Species %in% c('ASFA','AINC'),])

summary(glht(sequestration.model.combo, linfct = mcp(Mon.Pop = "Tukey"))) #Can't get estimates because model includes too many levels for the combination term (n = 36) for the model to be fully specified.

######################################################################################
######################################################################################
######################################################################################

#### Now, for the sake of comparison, create a reaction norm plot showing Puerto Rico and eastern North America with ASCU and ASYR. Build a function that enables pairwise reaction norms comparisons for any combination of two species and two populations.

rxn.norm.plot <- function(df1, pops, spp){
  df.use <- df1[df1$Species %in% spp & df1$Mon.Pop %in% pops,]
  agg.data <- do.call(data.frame, aggregate(concentration ~ Mon.Pop + Species, df.use, 
                        function(y) c(mean(y), sd(y) / sqrt(length(y)))))
  names(agg.data)[3:4] <- c('concentration', 'SE')
  ggplot(agg.data, aes(x = Species, y = concentration, group = Mon.Pop, fill = Mon.Pop))+
    geom_errorbar(aes(ymin = concentration - SE, ymax = concentration + SE, col = Mon.Pop), 
                  width = 0.1, size = 1.5)+
    geom_point(size = 5, pch = 21)+
    geom_line(size = 2, aes(col = Mon.Pop))+
    theme_light(base_size = 20)+
    scale_fill_manual(values = c('darkgreen','orange'))+
    scale_color_manual(values = c('darkgreen','orange'))+
    ylab('Wing Cardenolide Concentration (mg/g)')+
    theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), 
          legend.background = element_rect(fill = 'white', colour = 'black'))
}

(rxn.norm1 <- rxn.norm.plot(wing.cardenolides, pops = c('E. N. America','Puerto Rico'), spp = c('A. syriaca','A. curassavica')))

#To above figure, add family-specific reaction norms

family_averages <- aggregate(concentration ~ maternal_family + Mon.Pop + Species, wing.cardenolides[wing.cardenolides$Mon.Pop %in% c('E. N. America','Puerto Rico') & wing.cardenolides$Species %in% c('A. syriaca','A. curassavica'),], mean)

#pdf('./figures/Fig4C.pdf', height = 8, width = 6)
rxn.norm1 +
  geom_point(data = family_averages, aes(col = Mon.Pop), alpha = 0.2, size = 3)+
  scale_color_manual(values = c('darkgreen','orange'))+
  geom_line(data = family_averages, aes(group = maternal_family, col = Mon.Pop), alpha = 0.3, size = 1)+
  theme(axis.text.x = element_text(face = 'italic'))
#dev.off()

#statistical summary of above plot

rxn.norm.model <- lmer(concentration ~ Species*Mon.Pop + Sex + (1|maternal_family) + (1|Plant.ID), data = wing.cardenolides[wing.cardenolides$Species %in% c('A. syriaca','A. curassavica') & wing.cardenolides$Mon.Pop %in% c('Puerto Rico','E. N. America'),])

summary(rxn.norm.model)

######################################################################################
######################################################################################
######################################################################################

## Test whether the polarity of sequestered cardenolides is different for Puerto Rican monarchs

head(wing.cardenolides)

d2 <- merge(wing.cardenolides[,1:9], cardenolides.polarity, by = 'Sample')

d2$pri <- ifelse(d2$Mon.Pop=='Puerto Rico', 'Puerto\nRico', 'Other\nPopulations')

aggregate(polarity ~ pri + Species, data = d2[d2$Tissue=='Wing',], mean)
aggregate(polarity ~ pri + Species, data = d2[d2$Tissue=='Wing',], sd)

polarity.model <- lmer(polarity ~ pri + Species + (1|maternal_family), data = d2[d2$Species %in% c('A. syriaca','A. speciosa'),])

summary(polarity.model)

plot.polarity <- data.frame(emmeans(polarity.model, "pri", "Species"))

#For Figure S2-B, make plot showing lower polarity of sequestered compounds for Puerto Rican population on A. syriaca and A. speciosa. Use the same color scheme as Figure S2-A.

#pdf('./figures/FigS2b.pdf', height = 8, width = 3)
ggplot(plot.polarity, aes(x = pri, y = emmean, col = Species))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), width = 0.3, size= 0.7)+
  facet_wrap(~Species)+
  theme_light()+
  ylab('Polarity Index - Wing Tissue')+
  theme(axis.title.x = element_blank(), axis.text = element_text(size = 14), 
        axis.title.y = element_text(size = 16), axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(strip.text = element_text(size = 14), legend.position = 'none')+
  scale_color_manual(values = c('dodgerblue','turquoise'))
#dev.off()

######################################################################################
######################################################################################
######################################################################################

## Plot correlation between sequestration from GOPH and ASCU across populations

df3 <- aggregate(concentration ~   Species + Mon.Pop + maternal_family, mean, data = wing.cardenolides[wing.cardenolides$Species =='A. curassavica',])
names(df3)[4] <- 'ASCU_conc'

df4 <- aggregate(concentration ~  Species + Mon.Pop + maternal_family, mean, data = wing.cardenolides[wing.cardenolides$Species =='G. physocarpus',])
names(df4)[4] <- 'GOPH_conc'

plot.card.cor <- merge(df3[,c(3,4)], df4, by = 'maternal_family')

ggplot(plot.card.cor, aes(x = ASCU_conc, y = GOPH_conc))+
  geom_point(aes(col = Mon.Pop))

df5 <- aggregate(concentration ~ Species + Mon.Pop + maternal_family, mean, data = wing.cardenolides[wing.cardenolides$Species%in%c('A. curassavica','G. physocarpus'),])

ggplot(df5, aes(x = Species, y = concentration, col = Mon.Pop))+
  geom_point()+
  geom_line(aes(group = maternal_family))+
  scale_color_manual(values = c('cyan','limegreen','forestgreen','purple','blue','orange'))+
  facet_wrap(~Mon.Pop)+
  theme_light() #does seem to be some variation across maternal families in rank order of sequestration

######################################################################################
######################################################################################
######################################################################################

######## For Marianas comparisons: do ASCU genotypes from Rota have lower cardenolides than from Guam, or was this only the case in wild-collected plant tissue. This information goes into Supplemental Appendix 1 but is no longer part of the main manuscript.

cardenolides[cardenolides$Species=='A. curassavica' & cardenolides$Tissue=='Leaf',]$Sample

merge(cardenolides[cardenolides$Species=='A. curassavica' & cardenolides$Tissue=='Leaf',], adult.info[adult.info$Species=='ASCU', names(adult.info) %in% c('Sample','Pop')], by = 'Sample')

adult.info[adult.info$Species=='A. curassavica',names(adult.info) %in% c('Sample','Pop')]

ascu.leaf.cards <- cardenolides[cardenolides$Species=='A. curassavica' & cardenolides$Tissue=='Leaf',]
ascu.leaf.cards$Sample
ascu.leaf.cards$Pop <- ifelse(ascu.leaf.cards$Sample %in% c('1004','1005','1014','1017','1018','1020','1024','1025','1027','1028','1030','1034','1037','1042','1059','1062','1057','172','186','195','203','216','218','224','L170','L175','L188','L190','L199','L206','L215'), 'Guam', 'Rota')

guam.rota.ascu <- lm(concentration ~ Pop, data = ascu.leaf.cards) #not significantly higher, but def in the same direction as field samples
summary(guam.rota.ascu)

emmeans(guam.rota.ascu, specs = 'Pop')

ggplot(ascu.leaf.cards, aes(x = Pop, y = concentration))+
  geom_boxplot()

######################################################################################
######################################################################################
######################################################################################

#Create plots for what was originally Figure 5D: compare sequestration by monarchs from Guam on A. curassavica to monarchs from all other populations

#Unfortunately can not include the 5 monarchs collected in Guam in 2018 for comparison with wild caught monarchs from Guam in 2015, since they were analyzed in Rachel's lab on a different column

#### Add a column for adjusted cardenolide concentration (average wing concentration for Guam and other populations divided by global mean leaf concentration for A. curassavica), to allow for direct comparison with wild-caught butterflies

spp.leaf.means <- aggregate(concentration ~ Species, mean, data = cardenolides[cardenolides$Tissue=='Leaf',])

ASCU <- wing.cardenolides[wing.cardenolides$Species=='A. curassavica',]

ASCU$adjusted_values <- ASCU$concentration / spp.leaf.means[spp.leaf.means$Species=='A. curassavica',][[2]]

ASCU$guam <- ifelse(ASCU$Mon.Pop=='GU', 'Guam', 'Others')

aggregate(adjusted_values ~ guam, mean, data = ASCU)

#pdf('./figures/Figure5D_old.pdf',height = 7, width = 3)
ggplot(ASCU, aes(x = guam, y = adjusted_values, col = guam))+
  geom_boxplot(outlier.color = 'white', aes(fill = guam), alpha = 0.1)+
  geom_point(position = position_jitterdodge(0.25))+
  theme_light(base_size = 16)+
  theme(axis.title.x = element_blank(), legend.position = 'none')+
  ylab('Adjusted Cardenolide Concentration')+
  scale_color_manual(values = c('mediumorchid','orange'))+
  scale_fill_manual(values = c('mediumorchid','orange'))+
  ggtitle('Reared Monarchs -\nHost Plant Adjusted')+
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

######################################################################################
######################################################################################
######################################################################################

#For supplemental tables, show which individual compounds were most prevalent across each species

ascu.comp.aves <- colMeans(wing.cardenolides[wing.cardenolides$Species=='A. curassavica',19:88])

sort(ascu.comp.aves, decreasing = TRUE) #shows that frugoside is most abundant compound in wings of A. curassavica reared monarchs

######

#Based on standards that Rachel ran sent from Anurag, we know the following compounds. 
#5.93 is frugoside
#7.40 is calactin
#6.60 is calotropin
#1.12 is aspecioside (should be in ASYR and APSEC)

sum(sort(ascu.comp.aves, decreasing = TRUE)[c(1,4,5)]) / sum(ascu.comp.aves) #Frugoside, calactin, and calotropin comprise 51% of the total cardenolides sequestered from ASCU. Frugoside is the most abundant compound, with about 24% of the total amount. 

#do the same for other species

goph.comp.aves <- colMeans(wing.cardenolides[wing.cardenolides$Species=='G. physocarpus',19:88])

sort(goph.comp.aves, decreasing = TRUE) #For GOPH, frugoside is also the most abundant compound. Calactin is the fourth most abundant compound, and calotropin is sixth.
sum(sort(goph.comp.aves, decreasing = TRUE)[c(1,4,6)]) / sum(goph.comp.aves) #these compounds comprise 46% of the total sequestered from Gomphocarpus

asyr.comp.aves <- colMeans(wing.cardenolides[wing.cardenolides$Species=='A. syriaca',19:88])
sort(asyr.comp.aves, decreasing = TRUE) #For ASYR, 1.06 (likely to be aspecioside, which registered at 1.12) is by far the most abundant compound.

sort(asyr.comp.aves, decreasing = TRUE)[1] / sum(asyr.comp.aves) #Aspecioside comprises 48.3% of the total amount sequestered from ASYR

aspec.comp.aves <- colMeans(wing.cardenolides[wing.cardenolides$Species=='A. speciosa',19:88])
sort(aspec.comp.aves, decreasing = TRUE)

sort(aspec.comp.aves, decreasing = TRUE)[1] / sum(aspec.comp.aves) #Aspecioside comprises 44.0% of the total amount sequestered from ASPEC

######################################################################################
######################################################################################
######################################################################################

#Test for correlation between development time and sequestration

adult.info$Date <- as.Date(adult.info$Date, "%m/%d/%Y")

adult.info$development_time <- as.numeric(difftime(as.Date(adult.info$Eclosion, "%m/%d/%Y"),adult.info$Date, units = 'days') + 8)

df6 <- merge(adult.info[,c(1,20,31)], wing.cardenolides, by = 'Sample')

#pdf('./Figures/FigureS8.pdf',height = 8, width = 8)
ggplot(na.omit(df6[df6$development_time > 15 & df6$development_time < 30,]), 
       aes(x = development_time, y = concentration))+
  geom_point()+
  facet_wrap(~Species, scales = 'free_y')+
  stat_smooth(method = 'lm', formula = y ~ poly(x, 1), col = 'black')+
  theme_light(base_size = (16))+
  ylab('Cardenolide Concentration (mg/g)')+
  xlab('Days to Eclosion')+
  theme(legend.title = element_blank())
#dev.off()

summary(lmer(concentration ~ scale(development_time) + Sex + Species + Mon.Pop + (1|Plant.ID) + (1|maternal_family), data = df6)) #don't get same pattern as Agrawal paper, where sequestration was correlated with development rate (higher sequestration in slower developing caterpillars), at least for a model with a single slope estimate for development time

model.dev.time <- (lmer(concentration ~  Sex + (scale(development_time)|Species) + Mon.Pop + (1|Plant.ID) + (1|maternal_family), data = df6))
coef(model.dev.time) #when each species is allowed to have a varying slope for this relationship, then there is indeed a positive relationship between development time and sequestration on A. curassavica, consistent with the Agrawal et al. (2021) study. However, it is the only species for which this pattern is apparent, and it is still a relatively modest relationship (an increase of 1 SD of development time leads to a 2.6% increase in sequestered cardenolides)

######################################################################################
######################################################################################
######################################################################################

# Explicitly test for local adaptation while taking into account development time. Create a new metric (sequestered cardenolides divided by development time) and use this as a response variable.

df6$cards.adjusted <- df6$concentration/df6$development_time

df7 <- df6[df6$development_time > 15,] #exclude errant sample with a development time less than 15 days (probably a date/data entry error)

model.dev.time.adjusted <- lmer(cards.adjusted ~ Mon.Pop + Species + sym.allo + Sex + (1|Pop/Plant.ID) + (1|maternal_family), data = df7[!df7$Species %in% c('A. fascicularis','A. incarnata'),])

summary(model.dev.time.adjusted) #No change in inference, as suggested by the lack of a correlation between sequestration and development time shown in plots previous. Seems pretty clear that there is no overall sympatric sequestration advantage.

#Redo reaction norm figure

aggregate(cards.adjusted ~ Mon.Pop + Species, data = df6, mean)

df6$pr.vs.all <- ifelse(df6$Mon.Pop == 'PR', 'PR', 'other')

aggregate(cards.adjusted ~ pr.vs.all + Species, df6, mean)

adj.means <- aggregate(cards.adjusted ~ Mon.Pop + Species, df6[!df6$Species %in% c('A. fascicularis','A. incarnata'),], mean)
adj.sd <- aggregate(cards.adjusted ~ Mon.Pop + Species, df6[!df6$Species %in% c('A. fascicularis','A. incarnata'),], sd)
names(adj.sd)[3] <- 'SD'
plot.adj <- merge(adj.means, adj.sd)

(0.313 - 0.242) / 0.242  #29.3% higher adjusted sequestration on GOPH
(0.705 - 0.542) / 0.542 #30.0% higher adjusted sequestration of ASCU

#pdf('./figures/FigureS7a.pdf', height = 6, width = 6)
ggplot(plot.adj, aes(x = Mon.Pop, y = cards.adjusted))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = cards.adjusted-SD, ymax = cards.adjusted+SD), width = 0.3, size = 1)+
  theme_light(base_size = 16)+
  ylab('Sequestered Cardenolides (mg/g)\nper Day of Development')+
  xlab('Monarch Population')+
  facet_wrap(~Species, scales = 'free')+
  theme(legend.title = element_blank(), legend.position = 'bottom')+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))
#dev.off()

#Incorporating development time eliminates some of the modest sequestration advantage of Puerto Rican monarchs on ASCU, but they still have an advantage of GOPH. 

######################################################################################
######################################################################################
######################################################################################

#Within each species, test for correlations in production of major compounds

#pdf('./figures/FigureS3.pdf', height = 8, width = 8)
#par(mfrow = c(2, 2))

ASCU <- wing.cardenolides[wing.cardenolides$Species=='A. curassavica',]
ascu.M <- cor(ASCU[,names(sort(ascu.comp.aves, decreasing = TRUE)[1:6])])
corrplot(ascu.M, method = 'circle', title = 'A. curassavica', type = "upper", col=colorRampPalette(c("blue","white","orange"))(200), outline = 'black', diag = FALSE, mar = c(0,0,2,0))
GOPH <- wing.cardenolides[wing.cardenolides$Species=='G. physocarpus',]
goph.M <- cor(GOPH[,names(sort(goph.comp.aves, decreasing = TRUE)[1:6])])
corrplot(goph.M, method = 'circle', title = 'G. physocarpus', type = "upper", col=colorRampPalette(c("blue","white","orange"))(200), outline = 'black', diag = FALSE, mar = c(0,0,2,0))
ASYR <- wing.cardenolides[wing.cardenolides$Species=='A. syriaca',]
asyr.M <- cor(ASYR[,names(sort(asyr.comp.aves, decreasing = TRUE)[1:5])])
corrplot(asyr.M, method = 'circle', title = 'A. syriaca', type = "upper", col=colorRampPalette(c("blue","white","orange"))(200), outline = 'black', diag = FALSE, mar = c(0,0,2,0))
ASPEC <- wing.cardenolides[wing.cardenolides$Species=='A. speciosa',]
aspec.M <- cor(ASPEC[,names(sort(aspec.comp.aves, decreasing = TRUE)[1:5])])
corrplot(aspec.M, method = 'circle', title = 'A. speciosa', type = "upper", col=colorRampPalette(c("blue","white","orange"))(200), outline = 'black', diag = FALSE, mar = c(0,0,2,0))

#dev.off()
#par(mfrow = c(1, 1))

######################################################################################
######################################################################################
######################################################################################

#For Table S10, calculate coefficient of variation in sequestration from each species

aggregate(concentration ~ Species, function(x) sd(x)/mean(x), data = wing.cardenolides)

#Coefficient of variation is lowest on the two high cardenolide species, ASCU and GOPH.

######################################################################################
######################################################################################
######################################################################################

#Create figure showing multivariate disparity in sequestered cardenolides in relation to monarch population

#Start with NMDS plot for ASYR separated by population

ascl.nam.cards_nmds <- metaMDS(comm = wing.cardenolides[wing.cardenolides$Species == 'A. syriaca', 19:88], distance = 'bray', k = 3, trymax = 20)

nam.data.scores <- cbind(wing.cardenolides[wing.cardenolides$Species == 'A. syriaca',1:18], ascl.nam.cards_nmds$points)

ggplot(nam.data.scores, aes(x = MDS1, y = MDS2, col = Mon.Pop))+
  geom_point()+
  stat_ellipse(aes(col = Mon.Pop), level = 0.95)

#Alternative approach: take overall NMDS scores from more comprehensive overall analysis, then facet by species

wing.cardenolides.nmds <- wing.cardenolides[wing.cardenolides$Sample != '1003.1' & wing.cardenolides$Sample != '418.1' & wing.cardenolides$Sample != '822.1',] #first, omit three samples that have wonky NMDS scores (huge outliers that kind of make wonder whether they were either misassigned or contained only trace amounts; note that these were also omitted from earlier NMDS analyses)

ascl.cards_nmds <- metaMDS(comm = wing.cardenolides.nmds[wing.cardenolides.nmds$Species %in% c('A. syriaca','A. speciosa','G. physocarpus','A. curassavica'), 19:88], distance = 'bray', k = 5, trymax = 200) #takes about 5 minutes to run with 500 iterations and k = 10

ascl.wing.data.scores <- cbind(wing.cardenolides.nmds[wing.cardenolides.nmds$Species %in% c('A. syriaca','A. speciosa','G. physocarpus','A. curassavica'),1:18], ascl.cards_nmds$points)

#pdf('./figures/FigureS4.pdf', height = 8, width = 8)
ggplot(ascl.wing.data.scores, aes(x = MDS1, y = MDS2, col = Mon.Pop, fill = Mon.Pop))+
  geom_point()+
  facet_wrap(~Species, scales = 'free')+
  stat_ellipse(aes(col = Mon.Pop), level = 0.95, geom = 'polygon', alpha = 0.05)+
  theme_light(base_size = 16)+
  scale_color_manual(values = c('cyan','limegreen','forestgreen','purple','blue','orange'))+
  scale_fill_manual(values = c('cyan','limegreen','forestgreen','purple','blue','orange'))+
  theme(legend.position = 'bottom', legend.title = element_blank())
#dev.off()

#Conduct MANOVA to formally assess variance in profile of sequestered cardenolides that is attributable to milkweed species, monarch population, their interaction, and butterfly sex. Note that these results are reported in Table S6

ascl.cards.dist.mat <- vegdist(wing.cardenolides.nmds[wing.cardenolides.nmds$Species %in% c('A. syriaca','A. speciosa','G. physocarpus','A. curassavica'),19:88], method = 'bray')

adonis2(ascl.cards.dist.mat ~ Mon.Pop*Species + Sex, data = wing.cardenolides.nmds[wing.cardenolides.nmds$Species %in% c('A. syriaca','A. speciosa','G. physocarpus','A. curassavica'),], permutations = 1000) #Table S6

######################################################################################
######################################################################################
######################################################################################

#Test for correlation between leaf and wing cardenolides, using only matched samples (i.e. those where a particular wing sample can be associated with the particular host plant from which it came)

wing.cardenolides$sample_match <- wing.cardenolides$Plant.ID

leaf.cardenolides <- cardenolides[cardenolides$Tissue=='Leaf',]

leaf.cardenolides$sample_match <- readr::parse_number(leaf.cardenolides$Sample)

lw.cor <- merge(wing.cardenolides[,c('Mon.Pop','Sex','maternal_family','Species','concentration','sample_match')], leaf.cardenolides[,c('concentration','sample_match')], by = 'sample_match')

names(lw.cor)[6:7] <- c('wing_conc', 'leaf_conc')

lw.cor <- lw.cor[-which(lw.cor$Species=='A. syriaca' & lw.cor$leaf_conc > 3),] #drop one ASYR outlier whose leaf concentration is 10x higher than all other samples

#Create Figure S5, which shows correlation between paired leaf and wing cardenolide samples across species. 

#pdf('./figures/FigureS5.pdf', height = 6, width = 6)
ggplot(lw.cor[!lw.cor$Species %in% c('AINC','ASFA'),], 
       aes(x = leaf_conc, y = wing_conc))+
  theme_light(base_size = 16)+
  geom_point()+
  facet_wrap(~Species, scales = 'free')+
  stat_smooth(method = 'lm')+
  xlab('Leaf Cardenolide Concentration (mg/g)')+
  ylab('Wing Cardenolide Concentration (mg/g)')
#dev.off()

#Run a model that explicitly uses leaf concentration as a predictor variable 

lw.cor.model1 <- lmer(wing_conc ~ leaf_conc + Species*Mon.Pop + Sex + (1|maternal_family) + (1|sample_match), data = lw.cor[lw.cor$Species %in% c('A. curassavica','A. syriaca','A. speciosa','G. physocarpus'),])

summary(lw.cor.model1) #only get a weak positive association between leaf concentration and wing concentration

#Now test a model that does not include leaf concentration as a predictor variable for comparison.

lw.cor.model2 <- lmer(wing_conc ~ Species*Mon.Pop + Sex + (1|maternal_family) + (1|sample_match), data = lw.cor[lw.cor$Species %in% c('A. curassavica','A. syriaca','A. speciosa','G. physocarpus'),])

AIC(lw.cor.model1, lw.cor.model2) #model without leaf chemistry is preferred. Present these data in Table S7.

##### Get individual-level sequestration ratio values (corresponding to slopes of lines shown in Figure S1)

lw.cor$ratio <- lw.cor$wing_conc / lw.cor$leaf_conc

lw.cor$ratio <- ifelse(lw.cor$ratio == "Inf", lw.cor$ratio == NA, lw.cor$ratio) #replace all "inf" values (resulting from leaf concentrations of 0) with NA

aggregate(ratio ~ Species, mean, data = lw.cor)
aggregate(ratio ~ Species, length, data = lw.cor)
aggregate(ratio ~ Species, function(x) sd(x)/sqrt(length(x)), data = lw.cor)

#Make Figure S1. First create the dataframe containing the values that will be displayed.

Concentration <- c(lw.cor$leaf_conc, lw.cor$wing_conc)
Tissue <- c(rep('Leaf', nrow(lw.cor)), rep('Wing', nrow(lw.cor)))
Species <- rep(lw.cor$Species, 2)
Sample <- rep(lw.cor$sample_match, 2)

plot_lw.cor <- data.frame(Sample, Tissue, Species, Concentration)

#pdf('./figures/FigureS1.pdf', height = 6, width = 6)
ggplot(plot_lw.cor, aes(x = Tissue, y = Concentration, fill = Species))+
  geom_boxplot(outlier.color = 'white', alpha = 0.1)+
  geom_point(pch = 21)+
  geom_line(aes(group = Sample, col = Species))+
  scale_fill_manual(values = ascl.colors)+
  scale_color_manual(values = ascl.colors)+
  facet_wrap(~Species, ncol = 3, scales = 'free')+
  theme_light(base_size = 16)+
  theme(legend.position = 'none')+
  ylab('Cardenolide Concentration (mg/g)')+
  theme(strip.text = element_text(face = 'italic', size = 14))
#dev.off()
