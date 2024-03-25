library(ggplot2)
library(emmeans)
library(lme4)
library(scales)

anet_stbu = data.frame(readxl::read_excel(path = '~/Documents/GitHub/Sites/manuscripts/Island_Mainland/data_files/Anet-stbu.xlsx'))

head(anet_stbu)
table(anet_stbu$Provenance)

anet_stbu$Provenance <- as.factor(anet_stbu$Provenance)

anet_stbu$Provenance <- factor(anet_stbu$Provenance, labels = c("Mainland","Island"))

#simplest model possible, with gas exchange as a function of provenance (island vs. mainland)

model_anet1 <- lm(GasEx ~  Provenance, anet_stbu)

model_anet2 <- lm(GasEx ~  Provenance, anet_stbu[anet_stbu$ID != 107,]) #omit the one point with an exceptionally high SLA

summary(model_anet1)

summary(model_anet2)

################################################################################################

sla_sbbg = read.csv('~/Documents/GitHub/Sites/manuscripts/Island_Mainland/data_files/sla_sbbg.csv')

head(sla_sbbg)

comparison = merge(anet_stbu, sla_sbbg)

comparison$Site <- factor(comparison$Site, c('Gaviota','El Capitan','Santa Monicas','Zuma','Santa Cruz','Santa Rosa'))

#pdf('./figures/Fig5D.pdf', height = 6, width = 3.5)
ggplot(comparison, aes(x = SLA, y = GasEx))+
  theme_bw(base_size = 14)+
  geom_point(aes(fill = Site), size = 4, pch = 21)+
  geom_smooth(method = 'lm', col = 'black', alpha = 0.2, size = 0.5, lty = 2)+
  ylab(bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')'))+
  labs(x = expression ('SLA'~(cm^2/g)))+
  scale_fill_manual(values = c('darkorange4','orange', 'red','gold', 'blue','darkblue'))+
  scale_y_continuous(limit=c(0,11),oob = squish)+
  guides(fill = guide_legend(nrow = 3))+
  theme(legend.position = 'bottom', legend.title = element_blank())
#dev.off()  

model1_gasEx <- (lm(GasEx ~ SLA*Provenance, comparison)) #no evidence for interaction between provenance and SLA

summary(model1_gasEx)

model2_gasEx <- (lm(GasEx ~ SLA + Provenance, comparison)) #negative relationship between SLA and gas exchange; interestingly, after accounting for SLA differences, island plants have higher rates of gas exchange

summary(model2_gasEx)

summary(lm(GasEx ~ SLA + Provenance, comparison[comparison$SLA < 3.5,])) #after omitting the exceptional SLA point, the relationship between gas exchange and SLA goes away

emm.anet <- data.frame(emmeans(model2_gasEx, specs = ~ Provenance, type = "response"))

#pdf('./figures/Fig5C.pdf', height = 6, width = 2)
ggplot(emm.anet, aes(x = Provenance, y = emmean, col = Provenance))+
  theme_bw(base_size = 14)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.5)+
  geom_point(data = anet_stbu, aes(x = Provenance, y = GasEx, col = Provenance),
             position = position_jitterdodge(0.5), alpha = 0.2)+
  theme(legend.position = 'none')+
  scale_color_manual(values = c('red','blue'))+
  ylab(bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')'))
#dev.off()