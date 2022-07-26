######################################################################################
######## Script for processing data from Guam and Rota samples collected in 2015 #####
######################################################################################

#This script was used to produce Figures 5B and 5C as well Figure S1

library(ggplot2)
library(emmeans)
library(vegan)

#setwd()

marianas.wing.cardenolides.full <- read.csv(file = 'marianas_wing_cardenolides.csv')

head(marianas.wing.cardenolides)

#start by omitting Saipan, which only had two observations (both associated with Calotropis)
marianas.wing.cardenolides <- marianas.wing.cardenolides.full[marianas.wing.cardenolides.full$Site!='Saipan',]

#run basic linear model comparing total sequestration amounts between Guam and Rota, controlling for sex and the mass of the wing tissue 

guam.rota.model1 <- lm(Total ~ Site + sex + wingmass, data = marianas.wing.cardenolides)

summary(guam.rota.model1) #No significant difference between Guam and Rota. Males have modestly less than females

emmeans(guam.rota.model1, specs = 'Site')

#make boxplot showing relationship

#pdf('../figures/Figure5B.pdf',height = 7, width = 3)
ggplot(marianas.wing.cardenolides, aes(x = Site, y = Total, col = Site))+
  geom_boxplot(outlier.color = 'white', aes(fill = Site), alpha = 0.1)+
  geom_point(position = position_jitterdodge(0.25))+
  theme_light(base_size = 16)+
  theme(axis.title.x = element_blank(), legend.position = 'none')+
  ylab('Wing Cardenolide Concentration (mg/g)')+
  scale_color_manual(values = c('mediumorchid','turquoise4'))+
  scale_fill_manual(values = c('mediumorchid','turquoise4'))+
  ggtitle('Wild Monarchs -\nNon-Adjusted')+
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#compare leaf cardenolides from naturally-occuring plants

marianas.leaf.cardenolides <- read.csv(file = 'marianas_leaf_cardenolides.csv')

head(marianas.leaf.cardenolides)

marianas.leaf.cardenolides <- marianas.leaf.cardenolides[marianas.leaf.cardenolides$species=='Asclepias_curassavica',] #restrict comparison to only leaf samples from A. curassavica

guam.rota.leaf.model <- lm(Total ~ location, data = marianas.leaf.cardenolides)

summary(guam.rota.leaf.model)

emmeans(guam.rota.leaf.model, specs = 'location')

ggplot(marianas.leaf.cardenolides, aes(x = location, y = Total, col = location))+
  geom_boxplot(outlier.color = 'white', aes(fill = location), alpha = 0.1)+
  geom_point(position = position_jitterdodge(0.25))+
  theme_light(base_size = 16)+
  theme(axis.title.x = element_blank(), legend.position = 'none')+
  ylab('Leaf Cardenolide Concentration (mg/g)')+
  scale_color_manual(values = c('mediumorchid','turquoise4'))+
  scale_fill_manual(values = c('mediumorchid','turquoise4')) #modestly lower cardenolide concentrations on Rota compared to Guam when comparing wild-caught plants

#now use these average values from Guam and Rota ASCU to correct for sequestration levels in monarchs

marianas.wing.cardenolides$adjusted_values <- with(marianas.wing.cardenolides, ifelse(Site == 'Guam', Total/3.96, Total/2.72))

aggregate(adjusted_values ~ Site, marianas.wing.cardenolides, mean)

guam.rota.model2 <- lm(adjusted_values ~  Site + wingmass + sex, data = marianas.wing.cardenolides)

summary(guam.rota.model2) #now, the difference between Guam and Rota is significant

#pdf('../figures/Figure5C.pdf',height = 7, width = 3)
ggplot(marianas.wing.cardenolides, aes(x = Site, y = adjusted_values, col = Site))+
  geom_boxplot(outlier.color = 'white', aes(fill = Site), alpha = 0.1)+
  geom_point(position = position_jitterdodge(0.25))+
  theme_light(base_size = 16)+
  theme(axis.title.x = element_blank(), legend.position = 'none')+
  ylab('Adjusted Cardenolide Concentration')+
  scale_color_manual(values = c('mediumorchid','turquoise4'))+
  scale_fill_manual(values = c('mediumorchid','turquoise4'))+
  ggtitle('Wild Monarchs -\nHost Plant Adjusted')+
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

######################################################################################
######################################################################################
######################################################################################

# Figure S1 - NMDS plot of cardenolide composition across islands (here, including Saipan)

nmds.guam <- marianas.wing.cardenolides.full[,4:20]

nmds1 <- metaMDS(comm = nmds.guam, distance = 'bray', k = 3, trymax = 500) 

stressplot(nmds1)

plot.nmds.guam <- cbind(marianas.wing.cardenolides.full[,1:3], scores(nmds1))

#pdf('../figures/FigureS1.pdf',height = 6, width = 6)
ggplot(plot.nmds.guam, aes(x = NMDS1, y = NMDS2, fill = Site))+
  geom_point(size = 4, pch = 21)+
  theme_light(base_size = 16)+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  scale_fill_manual(values = c('mediumorchid','turquoise','darkgreen'))
#dev.off()