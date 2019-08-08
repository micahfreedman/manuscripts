library(ape)

milkweed.phylo <- read.nexus('~/Downloads/milkweeds-fasta/trees_milkweed_100000')

str(milkweed.phylo)

consensus.milkweed <- consensus(milkweed.phylo[20000:100000], p = 1)

consensus.milkweed$tip.label[6] <- 'Asclepias incarnata'

is.binary.tree(consensus.milkweed) #has a polytomy

resolved.milkweed.phy <- multi2di(consensus.milkweed)

resolved.milkweed.phy$tip.label[1] <- 'GOPH'
resolved.milkweed.phy$tip.label[2] <- 'ASPEC'
resolved.milkweed.phy$tip.label[3] <- 'ASYR'
resolved.milkweed.phy$tip.label[4] <- 'ASCU'
resolved.milkweed.phy$tip.label[5] <- 'ASFA'
resolved.milkweed.phy$tip.label[6] <- 'AINC'

plot(resolved.milkweed.phy, direction = 'downwards', cex = 2)

plot(resolved.milkweed.phy, cex = 2, tip.col = c('blue','darkorange4','darkseagreen4','purple','orange','darkolivegreen2'))



