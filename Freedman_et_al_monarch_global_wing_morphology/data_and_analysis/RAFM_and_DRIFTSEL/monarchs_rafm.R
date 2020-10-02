library(snpR); library(ggplot2); library(reshape2)
##########################
#get data
setwd("~/Documents/GitHub/Sites/manuscripts/Freedman_et_al_monarch_global_wing_morphology/data_and_analysis/RAFM_and_DRIFTSEL/")
raw_genos <- readRDS("flt_snps_clean.RDS")

#get into snpRdata format

genos <- raw_genos[,-c(1:3)]
snp_meta <- raw_genos[,1:3]
sample_meta <- data.frame(pop = substr(colnames(raw_genos)[-c(1:3)], 1, 3))
formatted_monarch_snps <- import.snpR.data(genos, snp.meta = snp_meta, sample.meta = sample_meta, mDat = "NN")


##########################
#get the coancestry matrix
source("./RAFM.r")

#now get into RAFM format

Rgenos <- format_snps(formatted_monarch_snps, "rafm", facets = "pop", n_samp = 1000) #takes a few minutes to run with 1000 snps

table(Rgenos$pop)

Rgenos.keep <- Rgenos[Rgenos$pop %in% c(1,12,3,4,7),]

#get the coancestries
RAFMall <- do.all(Rgenos.keep, 1000, 50, 25)  #specifies to use the 1000 snps, run analysis for 1000 iterations, with a burnin of 50 and sampling every 25 iterations
#need to make sure plot window is large enough, or you get an error about plot.new() figure margins being too small, and everything is for naught


######################################################################

#above code was used to generate ancestry-coancestry matrix, whose information is saved as theta; load this below

thetas <- read.csv('~/Documents/grad_school/Manuscripts/Global wing morphology/RAFM/RAFM_1000SNPs_10000_theta.csv') ### use the pre-generated values here

#get a matrix of the coancestry medians within and between each pop
comat <- matrix(NA, 5, 5)
for(i in 1:5){
  for(j in 1:5){
    comat[i,j] <- summary(thetas[i,j,])[3]
  }
}

comat[lower.tri(comat)] <- NA

comat <- as.data.frame(comat)
colnames(comat) <- 1:5
comat$pop <- c('Eastern\nNorth\nAmerica','Guam','Hawaii','Queensland','Western\nNorth\nAmerica')

mcomat <- melt(comat, id.vars = "pop")

mcomat$variable <- revalue(as.character(mcomat$variable), c('1' = "Eastern\nNorth\nAmerica", '2' = 'Guam', '3' = 'Hawaii', '4' = 'Queensland', '5' = 'Western\nNorth\nAmerica'))

#plot
ggplot(mcomat, aes(pop, variable, fill = value)) + 
  geom_tile(color = 'black') + 
  theme_minimal() +
  scale_fill_gradient2(low = "white", high = "red", na.value = 'powderblue') + geom_text(aes(label = round(value, 3)))+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(size = 12), axis.title = element_blank())+
  theme(legend.position = 'none')
