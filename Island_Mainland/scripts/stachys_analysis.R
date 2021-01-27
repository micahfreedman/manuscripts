##### Analysis of data pertaining to Stachys common garden

library(ggplot2)
library(vegan)
library(plotly)
library(plyr)
library(emmeans)

setwd('~/Documents/GitHub/Sites/manuscripts/Island_Mainland/')

stachys_sbbg <- read.csv('./data_files/Field_Setup_SBBG.csv')

head(stachys_sbbg)

stachys_sbbg$plant.ID <- paste(stachys_sbbg$Population, stachys_sbbg$Genotype, sep = "_")

with(stachys_sbbg[stachys_sbbg$Year==2017 & stachys_sbbg$Garden=='SBBG',], table(Population))
with(stachys_sbbg[stachys_sbbg$Year==2017 & stachys_sbbg$Garden=='SBBG',], table(plant.ID))

###################

stachys_leaf_vocs <- read.csv('./data_files/stachys_leaf_vocs.csv')

head(stachys_leaf_vocs)
names(stachys_leaf_vocs)
stachys_leaf_vocs$SourcePop <- revalue(stachys_leaf_vocs$SourcePop, c('El_Capitan' = 'El Capitan', 'Santa_Cruz' = 'Santa Cruz', 'Santa_Monica_Mtns' = "Santa Monicas", 'Santa_Rosa' = 'Santa Rosa'))
stachys_leaf_vocs$SourcePop <- factor(stachys_leaf_vocs$SourcePop, c('Gaviota','El Capitan','Santa Monicas','Zuma','Santa Cruz','Santa Rosa'))
names(stachys_leaf_vocs)[2] <- 'Population'

stachys_vocs_compound_list <- read.csv('./data_files/stachys_compound_list.csv')
head(stachys_vocs_compound_list)

#first, standardize samples based on recorded value of internal standard

stachys_standardized <- cbind(stachys_leaf_vocs[,1:4], stachys_leaf_vocs[,5:ncol(stachys_leaf_vocs)] / stachys_leaf_vocs$RT_15.864_internal_standard_tetralin)

stachys_standardized$RT_15.864_internal_standard_tetralin <- NULL #drop the internal standard column now that it's no longer needed

stachys_standardized$RT_38.141_phytol <- NULL #also drop phytol, which is not really a volatile compound

### make nMDS plot of compounds

stachys_NMDS<-metaMDS(comm=stachys_standardized[stachys_standardized$Site=='SBBG',5:ncol(stachys_standardized)], distance="bray", k=8, autotransform=F, trymax=1000) #does not always converge, but sometimes does with as few as 100 iterations

stressplot(stachys_NMDS)

nmds.data.scores <- as.data.frame(cbind((stachys_standardized[stachys_standardized$Site=='SBBG',c(1:4)]),scores(stachys_NMDS)))

pdf('./figures/Fig4b.pdf', height = 4, width = 6)
ggplot(nmds.data.scores, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(fill = Population), size = 4, pch = 21)+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))
dev.off()

#test with only mainland populations

stachys_mainland <- stachys_standardized[!stachys_standardized$Population %in% c('Santa Cruz', 'Santa Rosa'),]

stachys_mainland_NMDS<-metaMDS(comm=stachys_mainland[,5:ncol(stachys_mainland)], distance="bray", k=8, autotransform=F, trymax=2500) 

nmds.data.scores2 <- as.data.frame(cbind((stachys_mainland[,c(1:4)]),scores(stachys_mainland_NMDS)))

ggplot(nmds.data.scores2, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(fill = Population), size = 4, pch = 21)+
  theme_bw(base_size = 16)+
  scale_fill_manual(values = c('darkorange4','orange', 'red','gold')) #no obvious separation based on population of origin for mainland samples; instead, all seem to be composed of basically the same subset of compounds
dev.off()

## make plot of cumulative emission levels

stachys_standardized$cumulative <- rowSums(stachys_standardized[,5:ncol(stachys_standardized)])

ggplot(stachys_standardized, aes(x = Population, y = (cumulative)))+
  theme_bw(base_size = 16)+
  geom_boxplot(aes(fill = Population), alpha = 0.3)+
  geom_point(aes(col = Population), position = position_jitterdodge(0.5))+
  scale_fill_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  scale_color_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  ylab('ln(Cumulative Volatile Content)')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(legend.position = 'NULL') #gets the point across, but probably not worth including as a main figure

## split up dataset into compound classes, then make table of their production across populations

aliphatics <- cbind(stachys_standardized[,c(1:4)],stachys_standardized[,c(names(stachys_standardized) %in% stachys_vocs_compound_list[stachys_vocs_compound_list$Compound_Class=='aliphatic',]$Header_Compound)]) 
aromatics <- cbind(stachys_standardized[,c(1:4)],stachys_standardized[,c(names(stachys_standardized) %in% stachys_vocs_compound_list[stachys_vocs_compound_list$Compound_Class=='aromatic',]$Header_Compound)]) 
monoterpenes <- cbind(stachys_standardized[,c(1:4)],stachys_standardized[,c(names(stachys_standardized) %in% stachys_vocs_compound_list[stachys_vocs_compound_list$Compound_Class=='monoterpene',]$Header_Compound)]) 
sesquiterpenes <- cbind(stachys_standardized[,c(1:4)],stachys_standardized[,c(names(stachys_standardized) %in% stachys_vocs_compound_list[stachys_vocs_compound_list$Compound_Class=='sesquiterpene',]$Header_Compound)]) 
unknowns <- cbind(stachys_standardized[,c(1:4)],stachys_standardized[,c(names(stachys_standardized) %in% stachys_vocs_compound_list[stachys_vocs_compound_list$Compound_Class=='unknown',]$Header_Compound)]) 

aliphatics$cumulative <- rowSums(aliphatics[,5:ncol(aliphatics)])
aromatics$cumulative <- rowSums(aromatics[,5:ncol(aromatics)])
monoterpenes$cumulative <- rowSums(monoterpenes[,5:ncol(monoterpenes)])
sesquiterpenes$cumulative <- rowSums(sesquiterpenes[,5:ncol(sesquiterpenes)])
unknowns$cumulative <- rowSums(unknowns[,5:ncol(unknowns)])

aggregate(cumulative ~ Population, function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))), data = aliphatics)
aggregate(cumulative ~ Population, function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))), data = aromatics)
aggregate(cumulative ~ Population, function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))), data = monoterpenes)
aggregate(cumulative ~ Population, function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))), data = sesquiterpenes)
aggregate(cumulative ~ Population, function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))), data = unknowns)

########

#now get raw chromatograms to show for demonstration

dem.indexes <- stachys_standardized[stachys_standardized$cumulative %in% aggregate(cumulative ~ Population, stachys_standardized, max)$cumulative,]$PlantID #specify the samples from each population with the highest overall production of leaf VOCs

elcap <- read.csv('./data_files/stachys_chromatograms/EC81.D/tic_front.csv', skip = 2, col.names = c('RT','value'))
elcap$Population <- 'El Capitan'
gav <- read.csv('./data_files/stachys_chromatograms/GAV102.D/tic_front.csv', skip = 2, col.names = c('RT','value'))
gav$Population <- 'Gaviota'
sci <- read.csv('./data_files/stachys_chromatograms/SCI16.D/tic_front.csv', skip = 2, col.names = c('RT','value'))
sci$Population <- 'Santa Cruz'
smm <- read.csv('./data_files/stachys_chromatograms/SMM108.D/tic_front.csv', skip = 2, col.names = c('RT','value'))
smm$Population <- 'Santa Monicas'
sri <- read.csv('./data_files/stachys_chromatograms/SRI12.D/tic_front.csv', skip = 2, col.names = c('RT','value'))
sri$Population <- 'Santa Rosa'
zum <- read.csv('./data_files/stachys_chromatograms/ZUM_104.D/tic_front.csv', skip = 2, col.names = c('RT','value'))
zum$Population <- 'Zuma'

chromatograms.stachys <- rbind(elcap, gav, sci, smm, sri, zum)
chromatograms.stachys$Population <- factor(chromatograms.stachys$Population, c('Gaviota','El Capitan','Santa Monicas','Zuma','Santa Cruz','Santa Rosa'))

pdf('./figures/Fig4a.pdf', height = 9, width = 6)
ggplot(chromatograms.stachys[chromatograms.stachys$RT < 40,], aes(x = RT, y = value))+
  theme_dark(base_size = 10)+
  geom_line(aes(col = Population), size = 0.3)+
  facet_wrap(~Population, nrow = 6, strip.position = 'right')+
  scale_color_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  theme(legend.position = 'none')+
  ylab('Signal Intensity')+
  xlab('Retention Time (Minutes)')+
  theme(strip.text = element_text(size = 13))
dev.off()

#### final comparison: SCI genotypes grown at SBBG versus on SCI

sci_NMDS<-metaMDS(comm=stachys_standardized[stachys_standardized$Population=='Santa Cruz',5:ncol(stachys_standardized)], distance="bray", k=3, autotransform=F, trymax=250) 

sci.nmds.data.scores <- as.data.frame(cbind((stachys_standardized[stachys_standardized$Population=='Santa Cruz',1:4]),scores(sci_NMDS)))

pdf("./figures/Supplemental_Figures/Figure.Sx.pdf", height = 4, width = 4)
ggplot(sci.nmds.data.scores, aes(x = NMDS1, y = NMDS2, shape = Site, fill = Site))+
  geom_point(size = 4)+
  scale_fill_manual(values = c('blue','cyan'))+
  scale_shape_manual(values = c(21,22))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'bottom', legend.title = element_blank()) #no differentiation between SCI genotpyes grown on Santa Cruz versus on the mainland
dev.off()

####

stachys.dist.mat<-vegdist(stachys_standardized[5:ncol(stachys_standardized)], method='bray')

adonis2(stachys.dist.mat ~ SourceIM, data = stachys_standardized, permutations = 999, method = 'bray')

########################################################

#Analysis of Stachys biomass across years

stachys_biomass <- read.csv('./data_files/Field_Setup_SBBG.csv')

head(stachys_biomass)

stachys_biomass$IM <- ifelse(stachys_biomass$Population %in% c('Santa_Rosa','Santa_Cruz'), 'Island','Mainland')

stachys_biomass$Population = factor(stachys_biomass$Population, levels=c("El_Capitan","Gaviota","Santa_Monicas","Zuma","Santa_Cruz","Santa_Rosa"))

stachys_biomass$Population <- revalue(stachys_biomass$Population, c('El_Capitan' = 'El Capitan', 'Santa_Monicas' = 'Santa Monicas', 'Santa_Cruz' = 'Santa Cruz', 'Santa_Rosa' = 'Santa Rosa'))

stachys_biomass$total_mass <- with(stachys_biomass, (Cumulative_Mass - Bag_Mass) + (Inside_Bag_Mass - Inside_Bag_Only_Mass))

sbbg_stachys_biomass <- stachys_biomass[stachys_biomass$Garden == 'SBBG',] #filter out biomass records from the SCI garden

pdf(file = './figures/Fig6a.pdf', height = 6, width = 5)
ggplot(sbbg_stachys_biomass, aes(x = Population, y = total_mass))+
  geom_boxplot(aes(fill = Population), alpha = 0.3, outlier.color = 'white')+
  geom_point(aes(fill = Population), col = 'black', pch = 21, position = position_jitterdodge(1))+
  facet_wrap(~Year)+
  scale_fill_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  scale_color_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  theme_bw(base_size = 16)+theme(axis.title.x = element_blank(),
                                 axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(legend.position = 'none')+
  ylab('Aboveground Dry Biomass (g)')
dev.off()

####

stbu_model <- lmer(total_mass ~ IM + (1|Population/Genotype) + Year + (1|Row) + (1|Column), sbbg_stachys_biomass)

summary(stbu_model) #Strong year effects (not surprising given supplemental watering in 2016), and decent island/mainland effect, such that island plants generally had higher biomass 

emmeans(stbu_model, specs = pairwise ~ IM + Year, type = "response")

######

#now also analyze SLA

stachys.sla <- read.csv('./data_files/sla_sbbg.csv')
                                   
head(stachys.sla)
names(stachys.sla)[1] <- 'Index'

stachys.sla <- merge(stachys.sla, stachys_sbbg, by = 'Index')
stachys.sla$Population <- sub("_", " ", stachys.sla$Population)
stachys.sla$Population <- factor(stachys.sla$Population, c('Gaviota','El Capitan','Santa Monicas','Zuma','Santa Cruz','Santa Rosa'))

pdf(file = './figures/Fig6b.pdf', height = 6, width = 3)
ggplot(stachys.sla[stachys.sla$Garden == 'SBBG' & stachys.sla$Year==2017,], 
       aes(x = Population, y = SLA))+
  geom_boxplot(aes(fill = Population), alpha = 0.3, outlier.colour = 'white')+
  geom_point(aes(fill = Population), col = 'black', pch = 21, position = position_jitterdodge(1))+
  scale_fill_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  scale_color_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  theme_bw(base_size = 16)+theme(axis.title.x = element_blank(),
                                 axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(legend.position = 'none')+
  labs(y = expression ('Specific Leaf Area'~(cm^2/g)))+
  ylim(c(1.25,4))
dev.off()

stachys.sla$IM <- ifelse(stachys.sla$Population %in% c('Santa Cruz', 'Santa Rosa'), 'Island',"Mainland")

model_stachys_sla <- lmer(SLA ~ IM + (1|Population/Genotype) + (1|Row) + (1|Column) + leaves, data = stachys.sla[stachys.sla$Garden == 'SBBG' & stachys.sla$Year==2017,]) #includes # of leaves as a covariate (vast majority of samples had 6 leaves sampled, but a few had 8)

summary(model_stachys_sla) #island plants have significantly higher SLA
#now plot SLA for Santa Cruz genotypes from the island versus the mainland

ggplot(stachys.sla[stachys.sla$Year==2017 & stachys.sla$Population=='Santa Cruz',],
       aes(x = Garden, y = SLA, col = Garden, fill = Garden))+
  geom_boxplot(alpha = 0.2, outlier.color = 'white')+
  geom_point()+
  geom_line(aes(group = plant.ID), col = 'black', lty = 2)+
  scale_color_manual(values = c('blue','dodgerblue'))+
  scale_fill_manual(values = c('blue','dodgerblue'))

## plants grown on SCI had substantially lower SLA, probably due to receiving less water and being in generally much drier soil


