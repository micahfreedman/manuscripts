####### DRIFTSEL analysis 

###### load in all data and related scripts

setwd('~/Documents/GitHub/Sites/manuscripts/Freedman_et_al_monarch_global_wing_morphology/data_and_analysis/RAFM_and_DRIFTSEL/DRIFTSEL/input_data/')

source('../driftsel.r')

monarch.pedigree <- read.csv('./complete_monarchs_pedigree_driftsel.csv')

theta <- read.csv('../../RAFM_1000SNPs_10000_theta.csv')

#convert this into a list of 5x5 matrices summarizing theta

monarch.covars <- read.csv('./complete_monarch_covariates_driftsel.csv')

monarch.traits <- read.csv('./complete_monarch_traits_driftsel.csv')

monarchs.samp <- MH(theta, data.matrix(monarch.pedigree2), data.matrix(monarch.covars2), data.matrix(monarch.traits2)[,c(2:3,6:7)], 100000, 10000, 500, alt = T) #wing size, wing shape, abdomen mass, thorax mass

neut.test(monarchs.samp$pop.ef, monarchs.samp$G, monarchs.samp$theta, silent = F) #signature of selection value is only 0.369, suggesting that selection by itself is not especially strong relative to drift

monarchs.samp.size <- MH(theta, data.matrix(monarch.pedigree2), data.matrix(monarch.covars2), data.matrix(monarch.traits2)[,c(1,3)], 5000, 500, 100, alt = T)

neut.test(monarchs.samp.size$pop.ef, monarchs.samp.size$G, monarchs.samp.size$theta, silent = F) #univariate signal of selection on wing size: 0.431

##

monarchs.samp.shape <- MH(theta, data.matrix(monarch.pedigree2), data.matrix(monarch.covars2), data.matrix(monarch.traits2)[,c(1,4)], 5000, 500, 100, alt = T)

neut.test(monarchs.samp.shape$pop.ef, monarchs.samp.shape$G, monarchs.samp.shape$theta, silent = F) #univariate signal of selection on wing shape: 0.374


monarchs.samp.thorax <- MH(theta, data.matrix(monarch.pedigree2), data.matrix(monarch.covars2), data.matrix(monarch.traits2)[,c(1,6)], 5000, 500, 100, alt = T)

neut.test(monarchs.samp.thorax$pop.ef, monarchs.samp.thorax$G, monarchs.samp.shape$theta, silent = F) #univariate signal of selection on thorax weight: 0.416

monarchs.samp.abdomen <- MH(theta, data.matrix(monarch.pedigree2), data.matrix(monarch.covars2), data.matrix(monarch.traits2)[,c(1,7)], 5000, 500, 100, alt = T)

neut.test(monarchs.samp.abdomen$pop.ef, monarchs.samp.abdomen$G, monarchs.samp.abdomen$theta, silent = F) #univariate signal of selection on adbomen weight: 0.322

monarchs.samp.wl <- MH(theta, data.matrix(monarch.pedigree2), data.matrix(monarch.covars2), data.matrix(monarch.traits2)[,c(1,5)], 5000, 500, 100, alt = T)

neut.test(monarchs.samp.wl$pop.ef, monarchs.samp.wl$G, monarchs.samp.wl$theta, silent = F) #univariate signal of selection on wing loading: 0.456 (this is the single strongest signature of selection, which maybe makes sense)
