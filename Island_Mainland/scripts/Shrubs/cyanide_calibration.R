library(ggplot2)

calibration <- read.csv(file = './data_files/Shrubs/Cyanide/cyanide_calibration.csv')

summary(lm(abs ~ conc, calibration)) #returns formula that can be used for relating absorbance and concentration

summary(lm(abs ~ conc, calibration))[[9]] #extract R^2 for correlation between concentration and absorbance

#pdf('./figures/Supplemental_Figures/FigS13.pdf', height = 6, width = 6)
ggplot(calibration, aes(x = conc, y = abs))+
  theme_bw(base_size = 14)+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', se = F)+
  ylab('Absorbance Value')+
  xlab('Potassium Cyanide Concentration (mg/L)')
#dev.off()
