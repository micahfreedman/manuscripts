### Western Monarch Population Viability Analysis
### Schultz, Brown, Pelton and Crone
### Submitted to Biological Conservation, March 16, 2017

# code for running analysis de novo

##set working directory
##setwd("/Users/...")


library(msm)
library(car)
library(MARSS)
library(mgcv)
library(ggplot2)
library(ggrepel)

WMTC = read.csv("~/Downloads/1-s2.0-S0006320717304809-mmc1.csv")
dat=t(WMTC) #Transpose since MARSS needs time ACROSS columns
years = dat[1,]      
n = nrow(dat)-1
dat = dat[2:nrow(dat),]
nyrs = length(years) # years in the data set
n
# Create figure with # sites counted over time
count.nums = !is.na(dat)
plot(years, apply(count.nums, MARGIN = 2, FUN = sum), xlab = "year", ylab = "# sites counted", main = "", type = "o")

# Create figure with histogram of # years with counts vs # of sites
hist(apply(count.nums, MARGIN = 1, FUN = sum), main = "", breaks = 5:25, xlab = "# years with counts", ylab = "# sites")

#estimate parameters
Z.model = factor(rep("1,", n))
R.model = "diagonal and unequal" 

#note - on our computers, this takes about 20 minutes to run
kem1 = MARSS(dat, model=list(Z=Z.model, R=R.model))

#note - on our computers, this takes about 2 days to run - but CIs not needed to generate mu, sigma-squared or run model checks
#CIs = MARSSparamCIs(kem1,method = "parametric",nboot = 100)

# plot of estimated population sizes, with confidence limits
kem1.UI = kem1$states + 1.96*kem1$states.se
kem1.LI = kem1$states - 1.96*kem1$states.se

#Figure with abundance index
plot(years, exp((kem1$states)), ylim = (c(0, max(exp(kem1$states)))), xlab = "year", ylab = "abundance index")
polygon(c(years, rev(years)), (c(exp(kem1.UI), rev(exp(kem1.LI)))), col = "grey90", border = NA)
points(years, exp(kem1$states), type = "o", pch = 19, lwd = 2)

# delta method to estimate standard error of u
deltamethod(g = ~ exp((x18+x19+x20+x21+x22+x23+x24+x25+x26+x27)/10)/exp((x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17)/17), mean = c(kem1$states[20:36], kem1$states[1:10]), cov = diag(c(kem1$states.se[20:36]^2, kem1$states.se[1:10]^2))) # standard error of the estimated ratio of geometric mean abundance in 1980's vs. 2000's = 6.8
# so ratio is 28.7 +/- SE 0f 6.8, approx CI = 15-42; so if there are 300,000 now, there were at least 4.5 million in 1980's


#scaled abundance; year with greatest # of sum of counts in recent years is 2015, which was 293,435. 
#kem1$states in 2015 = 4.03. So, exp(4.03)*5216 = ~293,435; use 5216 as the scaling factor
scale_factor = 5216/1000000  # use "millions of butterflies" on the y-axis
plot(years, scale_factor*exp((kem1$states)), ylim = (c(0, scale_factor*1.8*max(exp(kem1$states)))), xlab = "year", ylab = "abundance index (millions of butterflies)")
polygon(c(years, rev(years)), (c(scale_factor*exp(kem1.UI), rev(scale_factor*exp(kem1.LI)))), col = "grey90", border = NA)
points(years, scale_factor*exp(kem1$states), type = "o", pch = 19, lwd = 2)
mean(scale_factor*exp((kem1$states))[21:36])

# scaled abundance on a log-base 10 scale
plot(years, log(scale_factor*exp((kem1$states)),10), ylim = log((c(0.03, scale_factor*1.8*max(exp(kem1$states)))),10), xlab = "year", ylab = "abundance index (millions of butterflies)", yaxt = "n", cex= 1.9, col = "white")
polygon(c(years, rev(years)), log((c(scale_factor*exp(kem1.UI), rev(scale_factor*exp(kem1.LI)))),10), col = 'grey90', border = NA)
points(years, log(scale_factor*exp(kem1$states),10), type = "o", pch = 19, lwd = 2, cex = 1.3, col = 'orange')
axis(side = 2, at = log(c(0.1,0.3,1,3,10),10), labels = c(0.1,0.3,1,3,10),  cex = 1.9)


# calculating annual growth rates (log-lambda) and standard errors of log-lambda using the delta method:
lambdas = (kem1$states[2:36])-(kem1$states[1:35])
lam.se = array()
for(i in 1:length(lambdas))
{
  lam.se[i] = deltamethod(~x1-x2, mean = kem1$states[i:(i+1)], cov = matrix(c(kem1$states.se[i]^2,0,0,kem1$states.se[i+1]^2), nrow = 2, byrow = T))
  
}
cbind(lambdas, lam.se)
lam.ui = lambdas + 1.96*lam.se
lam.li = lambdas - 1.96*lam.se

# plot of annual growth rates
lam.years = 1981:2015 # year of start of annual growth rate
plot(lam.years, lambdas, type = "o", xlab = "year", ylab = expression(paste("growth rate, ln( ", lambda, " )")), ylim = c(-3,3), col = "white", cex = 1.9)
polygon(c(1981:2015, rev(1981:2015)), (c(lam.ui, rev(lam.li))), col = "grey90", border = NA)
points(1981:2015, lambdas, type = "o", pch = 19, lwd = 2, cex= 1.3)
points(c(1981,2015), c(0,0), type = "l", lty = "dotted")


### model checks
# covariates for model checks
lam.years = 1981:2015 # year of start of annual growth rate
decades = c(rep("e", 10), rep("n",10), rep("t", 15)) # categorical variable for 80's, 90's and 21st century (2000-2016)

## tests for change in mean log-lambda over time
q1 = lm(lambdas ~ lam.years) # linear model
  Anova(q1)
q1g = gam(lambdas ~ s(lam.years)) # generalized additive model - this is the most general test and the one included in the manuscript
  summary(q1g)
q1f = lm(lambdas ~ decades) # model with decade as categorical variable, for comparison with gam
  Anova(q1f)

# plot of GAM model for publication
plot(q1g, xlab = "Year", rug = F, ylab = expression(paste("smoothed trend in growth rate, ln(  ", lambda, ")")), shade = T, shift = coef(q1g)[1], cex = 1.9) # about as boring a line as you can get

## tests for change in magnitude of residuals [i.e., increase or decrease in variance]
lam.resid = lambdas - mean(lambdas)# residual growth rate - we used this because it had the nicest pattern of residuals and esitmated parameters, i.e., approx. normal residuals and no negative estimates of the parameter (compare to squared residuals)

q2 = lm((lam.resid)^2 ~ lam.years)# linear model
  Anova(q2)
  
q2g = gam((lam.resid)^2 ~ s(lam.years)) # generalized additive model - this is the most general test and the one included in the manuscript
  summary(q2g) 
q2f = lm((lam.resid)^2 ~ decades) # model with decade as categorical variable, for comparison with gam
  Anova(q2f)
q2f.means = lm((lam.resid)^2 ~ -1+decades) # model with decade as categorical variable, for comparison with gam
  summary(q2f.means)
  
  
# plot of GAM model for publication
plot(q2g, ylim = c(-0.6, 2), xlab = "Year", rug = F, ylab = "smoothed trend in environmental stochasticity", shade = T, shift = coef(q2g)[1], cex = 1.9) # some evidence for more extreme values in the 1990's than other years.  No evidence for differences in variance between the 1980's and 2000's [ML curve for each period falls within the CI's for the other period]

# approx. estimates of SD in log-lambda
mean(predict(q2g)[1:8]) # 1980's
mean(predict(q2g)[9:18]) # 1990's
mean(predict(q2g)[19:35]) # 2000's
mean(predict(q2g))  # average value of squared deviations from gam
(kem1$par$Q) # average value of SD lambda in MARSS model

# comparison of residuals of |resid| with squared deviations, and the log of |resid|
q2a = gam(abs(lam.resid)~ s(lam.years))
q2b = gam(log(abs(lam.resid)) ~ s(lam.years))

hist(resid(q2a)) #residuals from absoulte value model
hist(resid(q2g)) # residuals from the squared-deviations model
hist(resid(q2b)) # residuals from the log-transformed |residuals| model

## autocorrelations
# test for significant autocorrelation of growth rates
acf(lambdas)
acf(lambdas, plot = F) # no significant lags in plot
cor.test(lambdas[1:34], lambdas[2:35])
cor.test(lambdas[1:33], lambdas[3:35]) 
cor.test(lambdas[1:32], lambdas[4:35])
cor.test(lambdas[1:31], lambdas[5:35]) 
cor.test(lambdas[1:30], lambdas[6:35]) 


# test for significant autocorrelation of abundance
acf(t(kem1$states)) # positive autocorrelation, but no negative values at all
cor.test(t(kem1$states)[1:34], t(kem1$states)[2:35])
cor.test(t(kem1$states)[1:33], t(kem1$states)[3:35])
cor.test(t(kem1$states)[1:32], t(kem1$states)[4:35])
cor.test(t(kem1$states)[1:31], t(kem1$states)[5:35]) 
cor.test(t(kem1$states)[1:30], t(kem1$states)[6:35])

#save file
#save.image(file="WM_PVA.RData")


##### create output .csv for comparison with eastern data

western_abundance <- data.frame("year" = years, "abundance" = as.vector(kem1$states), "abundance_scaled" = as.vector(scale_factor*exp(kem1$states)),"LI" = as.vector(scale_factor*exp(kem1.LI)), "UI" = scale_factor*exp(as.vector(kem1.UI)), "source" = rep("Schutlz", 36))

write.csv(western_abundance, file = '~/Desktop/god_never_again.csv', row.names = F)

### append data from 2017-2019 to this

wab_recent <- data.frame("year" = c(2017,2018,2019,2020), "abundance" = rep(NA, 4), "abundance_scaled" = c(0.192, 0.027, 0.029, 0.002),"LI" = rep(NA,4), "UI" = rep(NA,4), "source" = rep("Xerces", 4))

western_abundance <- rbind(western_abundance, wab_recent)

f1_format <- theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'))

(Fig1d <- ggplot(western_abundance, aes(x = year, y = abundance_scaled))+
  theme_dark()+
  f1_format+
  geom_point(aes(shape = source), size = 4, col = 'orange')+
  geom_line(col = 'orange')+
  geom_ribbon(aes(ymin=LI, ymax= UI), alpha = 0.3, fill = 'orange')+
  scale_y_log10(breaks = c(30, 10, 3, 1, 0.3, 0.1, 0.03))+
  ylab('Scaled Abundance (Millions of Butterflies)')+
  scale_shape_manual(values = c(16,18))+
  theme(legend.position = 'none')+
  xlim(c(1980,2019))+
  xlab('Year')+
  geom_vline(xintercept = 2016.5, lty = 2)+
  scale_x_continuous(n.breaks = 10))

ggsave(file="~/Desktop/test.pdf", plot=Fig1d, width=6, height=5)

#### read in eastern North American data

eastern <- read.csv('~/Documents/GitHub/Sites/manuscripts/Freedman_et_al_monarch_perspective_piece/analysis/eastern_overwintering.csv')

(Fig1c <- ggplot(eastern, aes(x = year, y = hectares))+
  theme_dark()+
  f1_format+
  geom_point(size = 4, col = 'yellow')+
  geom_line(col = 'yellow')+
  ylab('Abundance (Hectares Occupied)')+
  theme(legend.position = 'none')+
  xlab('Year')+
  stat_smooth(method = 'lm', col = 'black', fill = 'yellow', alpha = 0.25, 
              formula = y ~ splines::bs(x, 5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'), limits = c(0,20), 
                     oob = scales::squish)+
  scale_x_continuous(breaks = seq(1980, 2020, by = 5), lim = c(1980,2020)))

ggsave(file="~/Desktop/test2.pdf", plot=Fig1c, width=6, height=5)

## comparison between eastern and western overwintering numbers

plot.both <- merge(eastern, western_abundance, by =  'year')

summary(lm(abundance_scaled ~ acreage, data = plot.both))

(Fig1b <- ggplot(plot.both, aes(x = (hectares), y = abundance_scaled))+
  theme_classic()+
  f1_format+
  scale_y_log10(breaks = c(30, 10, 3, 1, 0.3, 0.1, 0.03))+
  geom_point(col = 'black', size = 2)+
  geom_text_repel(aes(label = year), size = 3)+
  ylab('Western Abundance (Millions of Butterflies)')+
  xlab('Eastern Abundance (Hectares Occupied)')+
  annotate("text", x = 17, y = 9.5, label = 'r = 0.030\n')+
  annotate("text", x = 17, y = 7.5, label = 'p = 0.374')+
  annotate("rect", xmin = 14.5, xmax = 19.5, ymin = 5.5, ymax = 15, alpha = 0.1, col = 'black'))

ggsave(file = '~/Desktop/test3.pdf', width = 6, height = 5)

#### now make same figure, but with one year time lag to see if there is any correlation (i.e. does eastern population size in year 'n' predict western population size in year 'n+1'?)



ggplot(plot.both[2:nrow(plot.both),], aes(x = acreage[1:nrow(plot.both)-1], y = abundance_scaled[2:nrow(plot.both)]))+
  theme_classic()+
  f1_format+
  scale_y_log10(breaks = c(30, 10, 3, 1, 0.3, 0.1, 0.03))+
  geom_point(col = 'black', size = 2)+
  stat_smooth(method = 'lm', col = 'black', fill = 'gray', 
              alpha = 0.4, lty = 2, size = 0.5)+
  geom_text_repel(aes(label = year), size = 3)+
  ylab("Western Abundance (year n+1)")+
  xlab("Eastern Abundance (year n)")

summary(lm(abundance_scaled[2:nrow(plot.both)] ~ acreage[1:nrow(plot.both)-1], data = plot.both[2:nrow(plot.both),])) #here, actually do get a significant predictive effect of previous year's eastern abundance on western abundance

pb97 <- plot.both[!plot.both$year %in% c(1995,1996),]

ggplot(pb97[2:nrow(pb97),], aes(x = acreage[1:nrow(pb97)-1], y = abundance_scaled[2:nrow(pb97)]))+
  theme_classic()+
  f1_format+
  scale_y_log10(breaks = c(30, 10, 3, 1, 0.3, 0.1, 0.03))+
  geom_point(col = 'black', size = 2)+
  stat_smooth(method = 'lm', col = 'black', fill = 'gray', 
              alpha = 0.4, lty = 2, size = 0.5)+
  geom_text_repel(aes(label = year), size = 3)+
  ylab("Western Abundance (year n+1)")+
  xlab("Eastern Abundance (year n)")

summary(lm(abundance_scaled[2:nrow(pb97)] ~ acreage[1:nrow(pb97)-1], data = pb97[2:nrow(pb97),])) #effect goes away when '96 and '97 excluded
