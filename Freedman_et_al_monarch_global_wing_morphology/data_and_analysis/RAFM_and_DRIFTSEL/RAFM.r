########################################################################################
########################################################################################
#
# ---- RAFM 1.2 source code --------------------------
#
########################################################################################
########################################################################################

adjust.this <-
function(rates, props, w, a, dirichlet=TRUE){
	adjustment = exp(a*(rates - 0.44))
	if( dirichlet ){
		output = props / adjustment
		output[which(output>10^3)] = 10^3
		output[which(output<10^(-3))] = 10^(-3)
		}
	else{ # only Gaussian distribution in this case
		output = props*adjustment
		}
	return(output)
	}

AFM <-
function(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1){

	# Allele counts
	dat <- compress.data(dat)
	
	# Initial values
	epsilon <- 10^(-eps)	
	nloci <- length(dat)
	npop <- nrow(dat[[1]])
	pAnc <- dat
	z <- dat
	logalpha <- rnorm(npop, 1, sqrt(2))
	for( j in 1:nloci ){
		nal <- ncol(dat[[j]])
		pAnc[[j]] <- rdirtrunc(rep(1,nal), epsilon)
		for( i in 1:npop ){
			z[[j]][i,] <- rdirtrunc(pAnc[[j]]*exp(logalpha[i]), epsilon)
			}
		}
	kap <- matrix(0.2/(npop-1), ncol=npop, nrow=npop)
	diag(kap) <- 0.8
		
	# Initial proposals
	propAnc <- rep(10, nloci)
	propz <- matrix(100, nrow=npop, ncol=nloci)
	propkap <- rep(100, npop)
	propalpha <- rep(0.1, npop)
	iterno <- 0 # needed for adaptation only
	
	# Prior parameters
	prioralpha <- pria
	priorkap <- matrix(0.2/(npop-1), npop, npop)
	diag(priorkap) <- 0.8
	priorkap <- prik*priorkap
	priorAnc <- pAnc
	for( j in 1:nloci ){
		priorAnc[[j]] <- rep(priq, length(priorAnc[[j]]))
		}	

	# Initial likelihood
	clike1 <- rep(0, nloci)
	clike2 <- matrix(0, ncol=nloci, nrow=npop)
	clike3 <- rep(0, npop)
	clike4 <- rep(0, npop)
	clike5 <- matrix(0, ncol=nloci, nrow=npop)
	for( j in 1:nloci ){
		clike1[j] <- l1(pAnc[[j]], priorAnc[[j]], epsilon)
		for( i in 1:npop ){
			clike2[i,j] <- l2(z[[j]][i,], logalpha[i], pAnc[[j]], epsilon)
			clike5[i,j] <- l5(kap[i,], z[[j]], dat[[j]][i,])
			}
		}
	for( i in 1:npop ){
		clike3[i] <- l3(logalpha[i], prioralpha)
		clike4[i] <- l4(kap[i,], priorkap[i,], i, epsilon)
		}

	# Accept ratios
	acAnc <- rep(0, nloci)
	acz <- matrix(0, nrow=npop, ncol=nloci)
	ackap <- rep(0, npop)
	acalpha <- rep(0, npop)
	
	# Output variables
	output <- c(kap, logalpha)

	# Monte Carlo loop
	for( i in 1:nMC ){
		if(i %% 1000 == 0){cat("Progress:", i, "\n")}
		# Status
			if( round(i/100)==i/100 ){
				print(paste("iter", i), quote=F)
				}
			iterno <- iterno + 1
			adjust <- (i <= burnin)

		# Updating ancestral frequencies
		for( j in 1:nloci ){
			oldlike <- clike1[j] + sum(clike2[,j])
			oldp <- pAnc[[j]]
			newp <- rdirtrunc(propAnc[j]*oldp, epsilon)
			nl1 <- l1(newp, priorAnc[[j]], epsilon)
			nl2 <- rep(0, npop)
			for( i in 1:npop ){
				nl2[i] <- l2(z[[j]][i,], logalpha[i], newp, epsilon)
				}
			there <- ddirtrunc(newp, propAnc[j]*oldp, epsilon, log=T)
			here <- ddirtrunc(oldp, propAnc[j]*newp, epsilon, log=T)
			accept <- nl1 + sum(nl2) + here - oldlike - there
			if( is.nan(accept) ){ accept <- -Inf }
			if( log(runif(1,0,1)) < accept ){
				pAnc[[j]] <- newp
				clike1[j] <- nl1
				clike2[,j] <- nl2
				acAnc[j] <- acAnc[j] + 1
				}
			}

		# Updating kappa
		for( i in 1:npop ){
			oldkap <- kap[i,]
			newkap <- rdirtrunc(oldkap*propkap[i], epsilon)
			oldlike <- clike4[i] + sum(clike5[i,])
			nl4 <- l4(newkap, priorkap[i,], i, epsilon)
			nl5 <- rep(0, nloci)
			for( j in 1:nloci ){
				nl5[j] <- l5(newkap, z[[j]], dat[[j]][i,])
				}
			here <- ddirtrunc(oldkap, propkap[i]*newkap, epsilon)
			there <- ddirtrunc(newkap, propkap[i]*oldkap, epsilon)
			accept <- nl4 + sum(nl5) + here - oldlike - there
			if( is.nan(accept) ){ accept <- -Inf }
			if( log(runif(1,0,1)) < accept ){
				kap[i,] <- newkap
				clike4[i] <- nl4
				clike5[i,] <- nl5
				ackap[i] <- ackap[i] + 1
				}
			}		

		# Updating alpha
		for( i in 1:npop ){
			olda <- logalpha[i]
			newa <- rnorm(1, olda, propalpha[i])
			oldlike <- clike3[i] + sum(clike2[i,])
			nl3 <- l3(newa, prioralpha)
			nl2 <- rep(0, nloci)
			for( j in 1:nloci ){
				nl2[j] <- l2(z[[j]][i,], newa, pAnc[[j]], epsilon)
				}
			accept <- nl3 + sum(nl2) - oldlike # symmetric
			if( is.nan(accept) ){
				accept <- -Inf 
				}
			if( log(runif(1,0,1)) < accept ){
				logalpha[i] <- newa
				clike2[i,] <- nl2
				clike3[i] <- nl3
				acalpha[i] <- acalpha[i] + 1
				}
			}
			
		# Updating frequencies in lineages
		for( j in 1:nloci ){
			for( i in 1:npop ){
				oldlike <- clike2[i,j] + sum(clike5[,j])
				oldz <- z[[j]][i,]
				newz <- rdirtrunc(propz[i,j]*oldz, epsilon)
				nl2 <- l2(newz, logalpha[i], pAnc[[j]], epsilon)
				nl5 <- rep(0, npop)
				Z <- z[[j]]
				Z[i,] <- newz
				for( k in 1:npop ){
					nl5[k] <- l5(kap[k,], Z, dat[[j]][k,])
					}
				here <- ddirtrunc(oldz, propz[i,j]*newz, epsilon)
				there <- ddirtrunc(newz, propz[i,j]*oldz, epsilon)
				accept <- nl2 + sum(nl5) + here - oldlike - there
				if( is.nan(accept) ){ accept <- -Inf }
				if( log(runif(1,0,1)) < accept ){
					z[[j]][i,] <- newz
					clike2[i,j] <- nl2
					clike5[,j] <- nl5
					acz[i,j] <- acz[i,j] + 1
					}
				}
			}
			output <- rbind(output, c(kap, logalpha)) # output variable
			
		# Adjusting proposals
		if( adjust ){
			w <- 1 - 0.1*exp(-i/1000)
			a <- 0.1*exp(-i/1000)
			iterno <- w*iterno
			ackap <- w*ackap
			acAnc <- w*acAnc
			acz <- w*acz
			acalpha <- w*acalpha
			ratekap <- ackap / iterno
			rateAnc <- acAnc / iterno
			ratez <- acz / iterno
			ratealpha <- acalpha / iterno
			propkap <- adjust.this(ratekap, propkap, w, a)
			propAnc <- adjust.this(rateAnc, propAnc, w, a)
			propz <- adjust.this(ratez, propz, w, a)
			propalpha <- adjust.this(ratealpha, propalpha, w, a, dirichlet=FALSE)
			}	
			} # MC loop closes

	# Burnin & thinning
	output <- output[(burnin+1):nMC,]
	imax <- floor(nrow(output) / thinning)
	totake <- thinning*1:imax
	output <- output[totake,]

	# Data out
	nmc_ <- nrow(output)
	nc <- ncol(output)
	kapm <- array(NA, dim=c(npop, npop, nmc_))
	alpham <- array(NA, dim=c(npop, nmc_))
	for( i in 1:nmc_ ){
		kapm[,,i] <- matrix(output[i,1:(npop^2)], ncol=npop, nrow=npop)
		logalpha <- output[i,(npop^2+1):nc]
		alpham[,i] <- logalpha
		}
	print("posteriors written", quote=F)
	return(list(kapm,alpham))
	}

CDFTruncatedBeta <-
function(x, alpha, beta, a, b){
	F = (pbeta(x, alpha, beta) - pbeta(a, alpha, beta)) / (pbeta(b, alpha, beta) - pbeta(a, alpha, beta))
	return(F)
	}

compress.data <-
function(rawDat){
	pops = rawDat[,1]
	DNA = rawDat[,2:ncol(rawDat)]
	dna = matrix(NA, ncol=ncol(DNA)/2, nrow=2*nrow(DNA))
	for( i in 1:(ncol(DNA)/2) ){
		dna[,i] = c(DNA[,2*i], DNA[,-1+2*i])
		}
	pops = rep(pops, 2)
	uniqs = unique(pops)
	nloci = ncol(dna)
	dat = list(1:nloci)
	for( j in 1:nloci ){
		dnaj = dna[,j]
		nas = which(is.na(dnaj))
		if( length(nas) > 0 ){
			dnaj = dnaj[-nas]
			popsj = pops[-nas]
			}
		else{ popsj = pops }
		genotypes = c(unique(dnaj), "unobs")
		datj = matrix(0, nrow=length(uniqs), ncol=length(genotypes))
		for( p in 1:length(uniqs) ){
			thispop = dnaj[which(popsj==uniqs[p])]
			datj[p, ] = sapply(genotypes, function(x) length(which(thispop==x)))
			}
		dat[[j]] = datj
		}
	return(dat)
	}

ddirtrunc <-
function(x, alpha, epsilon, log=TRUE){
	n = length(alpha)
	b = 1
	at = n*epsilon
	bt = n
	alphat = sum(alpha)
	xi = max(c( epsilon, 1-bt+b ))
	eta = min(c( b, 1-at+epsilon ))
	p = logPDFTruncatedBeta(x[n-1], alpha[n-1], alphat-alpha[n-1], xi, eta, epsilon)
	if( n > 2 ){
		for( k in (n-2):1 ){
			xi = max(c( epsilon/(1 - sum(x[(k+1):(n-1)])) , 1 - (bt - b*(n-k)) / (1-sum(x[(k+1):(n-1)])) ))
			eta = min(c( b/(1 - sum(x[(k+1):(n-1)])) , 1 - (at - epsilon*(n-k)) / (1-sum(x[(k+1):(n-1)])) ))
			tmp = x[k] / (1-sum(x[(k+1):(n-1)]))
			p = p + logPDFTruncatedBeta(tmp, alpha[k], alphat - sum(alpha[k:(n-1)]), xi, eta, epsilon)
			}
		for( k in (n-2):1 ){
			p = p - log(1 - sum(x[(k+1):(n-1)]))
			}
		}	
	if( !log ){ p = exp(p) }
	return(p)
	}

# The workhorse
#	
do.all <-
function(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1){
	primary <- AFM(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1)
	kap <- primary[[1]]
	alpha <- primary[[2]]
	secondary <- gen.theta(kap, alpha)
	theta <- secondary[[1]]
	fst <- secondary[[2]]
	alpha <- exp(alpha)
	npop = dim(theta)[1]
	print("Drawing trace of coancestry matrix theta, assess convergence.", quote=F)
	if( npop > 4 ){ 
		npop = 4
		print("N.B. Omitting populations 5+ from graph.", quote=F)
		}
    #dev.new() #may need to update depending on OS
	par(mfcol=c(npop,npop))
	for( i in 1:npop ){
		for( j in 1:npop ){
			headl = paste("theta[ ",i," , ",j," ]", sep="")
			plot( theta[i,j,], type='l', main=headl, xlab="", ylab="" )
			}
		}
	out = as.list(1:4)
	names(out) = c("theta", "fst", "kappa", "alpha")
	out[[1]] = theta
	out[[2]] = fst
	out[[3]] = kap
	out[[4]] = alpha
	return(out)
	}	

gen.fst <-
function(theta){
	npop = ncol(theta)
	offdiag = (sum(theta) - sum(diag(theta))) / (npop^2 - npop)
	fst = (mean(diag(theta)) - offdiag) / (1 - offdiag)
	return(fst)
	}

gen.theta <-
function(kapm, alpham){

	# Eq. 12, Karhunen & Ovaskainen
	alpham = exp(alpham) # away from log-normal scale
	fstm = rep(NA, ncol(alpham))
	thetam = array(NA, dim=dim(kapm))
	npop = nrow(alpham)
	for( i in 1:length(fstm) ){
		kap = kapm[,,i]
		alpha = alpham[,i]
		theta = matrix(NA, nrow=npop, ncol=npop)
		for( j in 1:npop ){
			for( k in 1:j ){
				theta[k,j] = sum( kap[k,]*kap[j,] / (alpha + 1) )
				theta[j,k] = theta[k,j]
				}
			}
		thetam[,,i] = theta
		fstm[i] = gen.fst(theta)
		}

	return(list(thetam, fstm))
	}

heat.theta <-
function(thetapost, mean=TRUE, verbose=TRUE){
	mod = as.matrix(0*thetapost[,,1])
	central = lower = upper = mod
	npop = dim(thetapost)[1]
	for( i in 1:npop ){
		for( j in 1:npop ){
			X = thetapost[i,j,]
			lower[i,j] = lower[j,i] = quantile(X, 0.025)
			upper[i,j] = upper[j,i] = quantile(X, 0.975)
			if( mean ){ 
				central[i,j] = central[j,i] = mean(X)
				}
			else{ central[i,j] = central[j,i] = median(X) }
			}
		}
	central = 0.5*(central + t(central))
	lower = 0.5*(lower + t(lower))
	upper = 0.5*(upper + t(upper))
	out = as.list(1:3)
	out[[1]] = central
	out[[2]] = lower
	out[[3]] = upper
	names(out) = c("mean", "lower", "upper")
	if( !mean ){
		names(out)[1] = "median"
		}
	heatmap(central, symm=T, margins=c(5,5), xlab="subpopulations, note permutation", ylab="", main="Coancestry matrix")
	if( verbose ){
		return(out)
		}
	}

l1 <-
function(p_, priorAnc_, epsilon){
	return(ddirtrunc(p_, priorAnc_, epsilon, log=TRUE))
	}

l2 <-
function(z_, logalpha_, p_, epsilon){
	return(ddirtrunc(z_, exp(logalpha_)*p_, epsilon, log=TRUE))
	}

l3 <-
function(logalpha_, prioralpha_){
	return(dnorm(logalpha_, prioralpha_[1], sqrt(prioralpha_[2]), log=TRUE))
	}

l4 <-
function(kap_, priorkap_, i_, epsilon){
	if( max(kap_)==kap_[i_] ){
		f = ddirtrunc(kap_, priorkap_, epsilon, log=TRUE)
		}
	else{ f = -Inf }
	return(f)	
	}

l5 <-
function(kap_, z_, dat){
	K = matrix(kap_, nrow=nrow(z_), ncol=ncol(z_))
	p = apply(K*z_, 2, sum)
	f = dmultinom(dat, size=sum(dat), prob=p, log=TRUE)
 	return(f)
	}

logBeta <-
function(a){
	return( sum(lgamma(a)) - lgamma(sum(a)) )
	}

logPDFDirichlet <-
function(x, a){
	return( sum((a-1)*log(x)) - logBeta(a) )
	}

logPDFTruncatedBeta <-
function(x, alpha, beta, a, b, eps){
	if( x < a || x > b ){
		f = -Inf
		}
	else{	
		f = (alpha-1)*log(x) + (beta-1)*log(1-x) - logBeta(c(alpha,beta))
		normalizing = log(pbeta(b, alpha, beta) - pbeta(a, alpha, beta))
		if( normalizing == -Inf ){ # no place for x to move between a and b
			f = 0
			}
		else{
			f = f - normalizing
			}
		}
	return(f)
	}

RandomDirichlet <-
function(a){
	eps = 10^(-7)
	x = rgamma(length(a), a)
	if( max(x) < 10^(-100) ){
		x = rep(0, length(a))
		x[sample(1:length(x))] = 1
		}
	x = x / sum(x)
	x[which(x<eps)] = eps
	x[which(x>1-eps)] = 1 - eps
	x = x / sum(x)
	return(x)
	}

RandomTruncatedBeta <-
function(alpha, beta, a, b){
	quantil = pbeta(a, alpha, beta) + runif(1,0,1)*(pbeta(b, alpha, beta) - pbeta(a, alpha, beta))
	x = qbeta(quantil, alpha, beta)
	return(x)
	}

rdirtrunc <-
function(alpha, epsilon){
	n = length(alpha)
	at = n*epsilon
	bt = n
	b = 1
	alphat = sum(alpha)
	x = rep(0,n)
	xi = max(c(epsilon, 1 - bt + b))
	eta = min(c(b, 1 - at + epsilon))
	x[n-1] = RandomTruncatedBeta(alpha[n-1], alphat-alpha[n-1], xi, eta)
	if( n > 2 ){
		for( k in (n-2):1 ){
			xi = max(c( epsilon/(1-sum(x[(k+1):(n-1)])) , 1 - (bt - b*(n-k))/(1 - sum(x[(k+1):(n-1)])) ))
			eta = min(c( b/(1-sum(x[(k+1):(n-1)])) , 1 - (at - epsilon*(n-k))/(1 - sum(x[(k+1):(n-1)])) ))
			tmp = RandomTruncatedBeta(alpha[k], alphat - sum(alpha[k:(n-1)]), xi, eta)
			x[k] = tmp*(1 - sum(x[(k+1):(n-1)]))
			}
		}
	x[n] = 1 - sum(x)
	x[which(x<epsilon)] = epsilon
	x[which(x>b)] = b
	return(x)
	}

viz.theta <-
function(thetapost, distance=T, center=F, main=NA){

	# Multidimensional scaling
	npop = dim(thetapost)[1]
	theta = D = matrix(NA, npop, npop)
	for( i in 1:npop ){
		for( j in 1:i ){
			theta[i,j] = theta[j,i] = mean(thetapost[i,j,])
			} }
	for( i in 1:npop ){
		for( j in 1:i ){
			D2 = 2*(theta[i,i] + theta[j,j] - 2*theta[i,j])
			D[i,j] = D[j,i] = sqrt(D2)
			} }
	D = rbind(D, sqrt(2*diag(theta)))
	D = cbind(D, c(sqrt(2*diag(theta)), 0))
	D = sqrt(2/pi)*D # scaling SD -> E absolute value
	fit = cmdscale(D, k=2)
	if( center ){ fit[nrow(fit),] = c(0,0) }
		
	# Plotting
	if( is.na(main) ){
		main = "Expected drift distances (units of ancestral SD)"
		if( !is.na(distance) ){
			if( !distance ){
				main = "Population-level coancestry"
				}
			}
		}	
	eps = max(abs(fit))*0.1
	M = max(abs(fit)) + eps
	cx1 = 1 - 1*is.na(distance)
	plot(fit[,1], fit[,2], ylab="MDS 2", main=main, pch=16, xlim=c(-M,M), ylim=c(-M,M), axes=T, xlab="MDS 1", cex=cx1)
	eps = max(abs(fit))*0.1
	colvec = 1:9
	colvec[7] = "orange"
	colvec[9] = "chartreuse"
	eps = (1 - 1*is.na(distance))*eps
	text(fit[1:npop,1]+eps, fit[1:npop,2], 1:npop, col=colvec)
	text(fit[npop+1,1]+eps, fit[npop+1,2], "A")
		
	# Distance information
	if( !is.na(distance) ){
		for( i in 1:(npop+1) ){
			for( j in 1:(npop+1) ){
				lines(c(fit[i,1],fit[j,1]), c(fit[i,2], fit[j,2]))
				}
			}	
		eps = 0.8*eps
		if( distance ){ dd = round(D,2) }
		else{
			self = diag(theta)
			dd = rbind(theta,self)
			dd = cbind(dd,c(self,0))
			dd = round(dd,2)
			}
		for( i in 1:nrow(fit) ){
			if( i>1 ){
				for( j in 1:(i-1) ){
					text(0.5*(fit[i,1]+fit[j,1]), 0.5*(fit[i,2]+fit[j,2])+eps, dd[i,j], cex=0.75)			
					}
				}
			}
		}
	}
	
viz.theta.exact <-
function(thetapost, distance=T, center=F, main=NA){

	# Multidimensional scaling
	nmc = dim(thetapost)[3]
	npop = dim(thetapost)[1]
	fit = matrix(0, ncol=2, nrow=npop+1)
	for( k in 1:nmc ){
		npop = dim(thetapost)[1]
		theta = D = matrix(NA, npop, npop)
		for( i in 1:npop ){
			for( j in 1:i ){
				theta[i,j] = theta[j,i] = thetapost[i,j,k]
				} }
		for( i in 1:npop ){
			for( j in 1:i ){
				D2 = 2*(theta[i,i] + theta[j,j] - 2*theta[i,j])
				D[i,j] = D[j,i] = sqrt(D2)
				} }
		D = rbind(D, sqrt(2*diag(theta)))
		D = cbind(D, c(sqrt(2*diag(theta)), 0))
		D = sqrt(2/pi)*D # scaling SD -> E absolute value
		fitk = cmdscale(D, k=2)
		if( center ){ fitk[nrow(fitk),] = c(0,0) }
		fit = fit + fitk / nmc # gradually towards mean
		}
		
	# Plotting
	if( is.na(main) ){
		main = "Expected drift distances (units of ancestral SD)"
		if( !is.na(distance) ){
			if( !distance ){
				main = "Population-level coancestry"
				}
			}
		}	
	eps = max(abs(fit))*0.1
	M = max(abs(fit)) + eps
	cx1 = 1 - 1*is.na(distance)
	plot(fit[,1], fit[,2], ylab="MDS 2", main=main, pch=16, xlim=c(-M,M), ylim=c(-M,M), axes=T, xlab="MDS 1", cex=cx1)
	eps = max(abs(fit))*0.1
	colvec = 1:9
	colvec[7] = "orange"
	colvec[9] = "chartreuse"
	eps = (1 - 1*is.na(distance))*eps
	text(fit[1:npop,1]+eps, fit[1:npop,2], 1:npop, col=colvec)
	text(fit[npop+1,1]+eps, fit[npop+1,2], "A")
		
	# Distance information
	if( !is.na(distance) ){
		for( i in 1:(npop+1) ){
			for( j in 1:(npop+1) ){
				lines(c(fit[i,1],fit[j,1]), c(fit[i,2], fit[j,2]))
				}
			}	
		eps = 0.8*eps
		if( distance ){ dd = round(D,2) }
		else{
			self = diag(theta)
			dd = rbind(theta,self)
			dd = cbind(dd,c(self,0))
			dd = round(dd,2)
			}
		for( i in 1:nrow(fit) ){
			if( i>1 ){
				for( j in 1:(i-1) ){
					text(0.5*(fit[i,1]+fit[j,1]), 0.5*(fit[i,2]+fit[j,2])+eps, dd[i,j], cex=0.75)			
					}
				}
			}
		}
	}
