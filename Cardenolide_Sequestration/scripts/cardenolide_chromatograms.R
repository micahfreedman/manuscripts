### Plotting raw data from cardenolide readout

library(dplyr)
library(ggplot2)

setwd("~/Documents/GitHub/Sites/manuscripts/Cardenolide_Sequestration/data/chromatograms/")

all_chromatograms <- list() #empty list to read in data

for(i in seq_along(list.files())){
  all_chromatograms[[i]] <- read.table(file = list.files()[i], skip = 44, header = F)
} #loop through directory; each file has metadata in the header that spans the first 44 lines

names(all_chromatograms) <- sub("\\_.*", "", list.files()) #name each dataframe according to its file name (chooses all characters prior to the first underscore)

sample_id <- names(all_chromatograms)

all_chromatograms <- mapply(cbind, all_chromatograms, sample = sample_id, SIMPLIFY = F) #adds a sample ID column to each dataframe 

df1 <- do.call(rbind, all_chromatograms) #converts into a single dataframe

names(df1) <- c('RT', 'Step', 'Value', 'sample') #rename columns

df1$Value <- as.numeric(ifelse(df1$Value < 0, 0, df1$Value)) #replace all signal intensity values that are negative with 0s

df1$species <- with(df1, ifelse(sample %in% c('1014.2','1020','1030.1','1039.3','1051.1'), 'ASCU', 
                          ifelse(sample %in% c('406.2','433'), 'GOPH', 
                          ifelse(sample %in% c('536','547.1'), 'AINC', 
                          ifelse(sample %in% c('702','749.1'), 'ASPEC', 
                          ifelse(sample %in% c('801.1P','805.4','807.1AU','844'), 'ASYR', 'ASFA'))))))

ggplot(df1, aes(x = RT, y = Value, col = species))+
  geom_line()+
  facet_wrap(~sample, scales = 'free') #quick plot to see all data. Note that the internal standard is just after 10 minutes, so anything after is not important. Also, large peak immediately at the start of each chromatogram was not a cardenolide, so this can be omitted.

ggplot(df1[df1$RT > 0.7 & df1$RT < 12 & df1$sample %in% c('1020','1051.1'),], 
       aes(x = RT, y = Value, col = sample))+
  geom_line()+
  ylab('Signal Intensity (mAU)') #example of a leaf and wing sample together; baseline is slightly higher for leaf tissue, but can see correspondence between some peaks (e.g. 2.3 mins, 5.8 mins, and the IS)

############

#actual figure to display: choose representative wing chromatograms for each milkweed species and plot together

display.samples <- c('406.2','547.1','749.1','805.4','940.1','1051.1') 

df1$species <- factor(df1$species, c('GOPH','ASCU','AINC','ASFA','ASYR','ASPEC'))
ascl.colors <- c('blue','purple','coral','gold','dodgerblue','turquoise')

#pdf('../../figures/Figx.pdf', height = 10, width = 4)
ggplot(df1[df1$RT > 0.7 & df1$RT < 11 & df1$sample %in% display.samples,],
       aes(x = RT, y = Value, col = species))+
  geom_line(size = 0.5)+
  geom_area(aes(fill = species), alpha = 0.3)+
  theme_light(base_size = 18)+
  facet_wrap(~species, ncol = 1, strip.position = 'right')+
  ylab('Signal Intensity (mAU)')+
  xlab('Retention Time (Minutes)')+
  scale_fill_manual(values = ascl.colors)+
  scale_color_manual(values = ascl.colors)+
  theme(legend.position = 'none')+
  ylim(c(0,150))
#dev.off()

#######

#Plotting single wing samples for display

ggplot(df1[df1$RT > 0.7 & df1$RT < 12 & df1$sample =='805.4',], 
       aes(x = RT, y = Value))+
  geom_line(size = 0.5)+
  geom_area(aes(fill = species), alpha = 0.3)+
  theme_light(base_size = 18)+
  ylab('Signal Intensity (mAU)')+
  xlab('Retention Time (Minutes)')+
  theme(legend.position = 'none')+
  ylim(c(0,150))

#######

#Contrast Puerto Rico and 805.4 sample on ASYR

df2 <- df1[df1$sample %in% c('807.1AU','801.1P'),]

df2$display.value <- ifelse(df2$sample=='807.1AU', df2$Value, df2$Value*-1)

df2$pop <- ifelse(df2$sample=='801.1P', "PR", "Australia")

pdf('~/Documents/GitHub/Sites/manuscripts/Cardenolide_Sequestration/figures/chromatogram_asyr_comp.pdf', height = 5, width = 7)
ggplot(df2[df2$RT > 0.7 & df2$RT < 12,], 
       aes(x = RT, y = display.value))+
  geom_line(size = 0.5, aes(col = pop))+
  geom_area(aes(fill = pop), alpha = 0.2)+
  theme_dark(base_size = 14)+
  ylab('Signal Intensity')+
  xlab('Retention Time (Minutes)')+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylim(c(-75,75))+
  scale_fill_manual(values = c('cyan','orange'))+scale_color_manual(values = c('cyan','orange'))+
  theme(axis.text.y = element_blank())+
  annotate("text", x = 1.1, y = -35, label = 'aspecioside', angle = 90, color = 'white')+
  annotate("text", x = 10, y = -35, label = 'digitoxin\n(value truncated)',angle = 90, color = 'white')
dev.off()

#Now do the same on ASCU (1039.3 = PR, 1030.1 = ENA)

df3 <- df1[df1$sample %in% c('1051.1','1039.3'),]

df3$display.value <- ifelse(df3$sample=='1051.1', df3$Value, df3$Value*-1)

ggplot(df3[df3$RT > 0.7 & df3$RT < 12,], 
       aes(x = RT, y = display.value))+
  geom_line(size = 0.5, aes(col = sample))+
  geom_area(aes(fill = sample), alpha = 0.1)+
  theme_light(base_size = 18)+
  ylab('Signal Intensity')+
  xlab('Retention Time (Minutes)')+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c('orange','forestgreen'))+scale_color_manual(values = c('orange','forestgreen'))+
  theme(axis.text.y = element_blank())+
  ylim(c(-150,150))

########################################################


