########################################################################################
########################################################################################
#
# ---- driftsel 2.1.2 source code --------------------------
#
########################################################################################
########################################################################################

calc.blocks <-
function(blocks){
	blocks = as.matrix(blocks[,-1], nrow=nrow(blocks)) # IDs off
	nb = ncol(blocks)
	dimmat = nrow(blocks)
	output = array(0, dim=c(dimmat, dimmat, nb))
	for( i in 1:nb ){
		VB = matrix(0, ncol=dimmat, nrow=dimmat)
		levs = unique(blocks[,i])
		for( k in levs ){
			targ = which(blocks[,i]==k)
			VB[targ,targ] = 1
			}
		output[,,i] = VB
		}
	return(output)
	}

calc.loads <-
function(ped, THP){
	n = nrow(ped)
	npop = ncol(THP)
	wpop = matrix(0, nrow=n, ncol=npop)
	
	# Founders
	founders = which(!is.na(ped[,4]))
	for( i in founders ){
		wpop[i,ped[i,4]] = wpop[i,ped[i,4]] + 1
		}
	
	# Offspring
	maxf = max(founders)
	for( i in (maxf+1):n ){
		s = ped[i,2] # id codes
		d = ped[i,3]
		s = min(which(ped[,1]==s)) # row numbers
		d = min(which(ped[,1]==d))
		wpop[i,] = 0.5*(wpop[s,] + wpop[d,]) # 'pure loading'		
		}
		
	# Out	
	return(wpop)
	}

calc.M.alt <-
function(wpop, THP, ped){

	# Population-level effects
	n = nrow(wpop)
	npop = ncol(THP)
	Cov_pp = theta = matrix(0, n, n)
	for( i in 1:n ){
		for( j in 1:i ){
			Cov_pp[i,j] = t(wpop[i,]) %*% THP %*% wpop[j,]
			Cov_pp[j,i] = Cov_pp[i,j]
			}
		}
		
	# Total coancestry
	founders = which(!is.na(ped[,4]))
	for( i in founders ){
		for( j in 1:i ){
			if( i==j ){ 
				theta[i,j] = 0.5 + 0.5*THP[ped[i,4],ped[i,4]] 
				}
			else{ theta[i,j] = THP[ped[i,4],ped[j,4]] }
			theta[j,i] = theta[i,j]
			}
		}
	offs = which(is.na(ped[,4]))
	for( i in offs ){
		s = ped[i,2] # id codes
		d = ped[i,3]
		s = min(which(ped[,1]==s)) # row numbers
		d = min(which(ped[,1]==d))
		for( j in 1:i ){
			if( j<i ){
				theta[i,j] = 0.5*(theta[j,s] + theta[j,d])
				}
			else{ # self
				theta[i,i] = 0.5 + 0.5*theta[s,d]
				}
			theta[j,i] = theta[i,j]
			}
		}
	
	# Residual covariance
	M = theta - Cov_pp
	return(M)
	}

calc.M <-
function(pars, THP, THB, censped){
	n = nrow(THB)
	dif = matrix(0, nrow=n, ncol=n)
	for( i in 1:n ){
		for( j in 1:i ){	
			x = 0
			if( (censped[i]+censped[j])==0 ){
				if( THB[i,j]>=0.25 ){ # full-sibs or self
					S = pars[i,3]
					D = pars[i,4]
					x = 0.125*( THP[S,S] + THP[D,D] )
					if( i==j ){ x = 2*x }
					}
				if( THB[i,j]==0.125 ){ # half-sibs
					if( pars[i,3]==pars[j,3] ){
						S = pars[i,3]
						x = 0.125*THP[S,S]
						}
					if( pars[i,4]==pars[j,4] ){
						D = pars[i,4]
						x = 0.125*THP[D,D]
						}
					}
				}
			dif[i,j] = dif[j,i] = x
			}
		}
	M = THB - dif	
	return(M)
	}

calc.W.alt <-
function(wpop, covars, traits, THP){

	# Dimensions
	nt = ncol(traits)
	np = ncol(THP)
	nx = ncol(covars)
	n = nrow(covars)

	# Parental contributions
	indic = wpop
		
	# Censoring in covariates
	censed = which( apply(covars,1,function(x)length(which(is.na(x)))) > 0 )
	if( length(censed) > 0 ){
		nclass = apply(covars[-censed,,drop=F],2,function(x)length(unique(x)))	
		binaries = which(nclass<=2)
		covars[censed,-binaries] = 0
		covars[censed,1] = 1 # censored observations go to corner class
		for( j in binaries ){	
			covars[censed,j] = mean(covars[,j], na.rm=T) # can be either
			}
		}	

	# covars and indic are constants accross traits,
	# however they need to be changed into design matrices
	X = matrix(0, nrow=nt*n, ncol=nt*nx)
	Z = matrix(0, nrow=nt*n, ncol=nt*np)
	for( i in 1:nt ){
		lo = 1 + n*(i-1)
		hi = n*i
		lox = 1 + nx*(i-1)
		hix = nx*i
		loz = 1 + np*(i-1)
		hiz = np*i
		X[lo:hi,lox:hix] = covars
		Z[lo:hi,loz:hiz] = indic
		}
		
	W = cbind(X, Z)
	return(W)
	}

calc.V <-
function(g, THP, G, priors){ # changes in the default application
	if( is.na(priors$fixed[[2]][1,1]) ){
		ng = length(g)
		n = ng - nrow(G)*nrow(THP)
		V = matrix(0, ncol=ng, nrow=ng)
		VF = 100*diag(n)
		Va = 2*G %x% THP
		V[1:n,1:n] = VF
		V[(n+1):ng,(n+1):ng] = Va # i.e. Va is a function of (G, THP)
		}
	else{ V = priors$fixed[[2]] }	
	return(V)
	}

calc.W <-
function(pars, covars, traits, THP){

	# Dimensions
	nt = ncol(traits)
	np = ncol(THP)
	nx = ncol(covars)
	n = nrow(covars)

	# Parental contributions
	indic = matrix(0, ncol=np, nrow=n)
	for( i in 1:n ){
		if( !is.na(pars[i,3]) ){
			indic[i,pars[i,3]] = indic[i,pars[i,3]] + 0.5  
			}
		if( !is.na(pars[i,4]) ){
			indic[i,pars[i,4]] = indic[i,pars[i,4]] + 0.5
			}
		}
		
	# Censoring in covariates
	censed = which( apply(covars,1,function(x)length(which(is.na(x)))) > 0 )
	if( length(censed) > 0 ){
		nclass = apply(covars[-censed,,drop=F],2,function(x)length(unique(x)))	
		binaries = which(nclass<=2)
		covars[censed,-binaries] = 0
		covars[censed,1] = 1 # censored observations go to corner class
		for( j in binaries ){	
			covars[censed,j] = mean(covars[,j], na.rm=T) # can be either
			}
		}
			
	# covars and indic are constants accross traits,
	# however they need to be changed into design matrices
	X = matrix(0, nrow=nt*n, ncol=nt*nx)
	Z = matrix(0, nrow=nt*n, ncol=nt*np)
	for( i in 1:nt ){
		lo = 1 + n*(i-1)
		hi = n*i
		lox = 1 + nx*(i-1)
		hix = nx*i
		loz = 1 + np*(i-1)
		hiz = np*i
		X[lo:hi,lox:hix] = covars
		Z[lo:hi,loz:hiz] = indic
		}
		
	W = cbind(X, Z)
	return(W)
	}

check.check <-
function(X, limes){
	solve(X/limes)
	solve(X*limes)
	return(FALSE) # this replaces TRUE which indicates problems
	}

check.order <-
function(traits, binary){
  if( !is.na(binary) ){
    if( length(binary) < ncol(traits) ){
      guide = 1:ncol(traits)
      nonbin = guide[-binary]
      if( min(binary) < max(nonbin) ){
        print("Order of traits changed in output for computational convenience,", quote=F)
        print("present order:", quote=F) 
        print(guide[c(nonbin,binary)])
        print("N.B. G matrix, pop.ef, etc.", quote=F)
        }
      }
    }
  }

const.sig <-
function(E, VR, G, M, blocks, Blist, bbit, sparse){
	Sig = (2*G %x% M) + (E %x% VR) # so far so good
	if( bbit ){
		for( i in 1:(dim(blocks)[3]) ){
			VB = blocks[,,i]
			if( sparse ){
				VB = as.matrix.csr(VB)
				}
			Sig = Sig + (Blist[,,i] %x% VB)
			}
		}
	return(Sig)
	}

convert.2 <-
function(dat, dd, nb){
	nn = nrow(dat)
	out2 = list(1:nb)
	for( k in 1:nb ){
		outthis = array(NA, dim=c(dd, dd, nn))
		targ = ((k-1)*dd^2) + 1:(dd^2)
		datthis = dat[,targ,drop=F]
		for( i in 1:nn ){
			outthis[,,i] = matrix(datthis[i,], nrow=dd, ncol=dd)
			}
		out2[[k]] = outthis
		}
	return(out2)
	}

convert <-
function(dat, dimy, npop, dimx){
	out = list(1:5)
	n = nrow(dat)
	out[[1]] = array(NA, dim=c(dimx, dimy, n))
	out[[2]] = array(NA, dim=c(npop, dimy, n))
	out[[3]] = array(NA, dim=c(dimy, dimy, n))
	out[[4]] = array(NA, dim=c(dimy, dimy, n))
	out[[5]] = array(NA, dim=c(npop, npop, n))
	for( i in 1:n ){
		dati = dat[i,]
		out[[1]][,,i] = matrix(dati[1:(dimx*dimy)], nrow=dimx)
		dati = dati[-(1:(dimx*dimy))]
		out[[2]][,,i] = matrix(dati[1:(npop*dimy)], nrow=npop)
		dati = dati[-(1:(npop*dimy))]
		out[[3]][,,i] = matrix(dati[1:(dimy*dimy)], nrow=dimy)
		dati = dati[-(1:(dimy*dimy))]
		out[[4]][,,i] = matrix(dati[1:(dimy*dimy)], nrow=dimy)
		dati = dati[-(1:(dimy*dimy))]
		out[[5]][,,i] = matrix(dati[1:(npop*npop)], nrow=npop)		
		}
	names(out) = c("fixed.ef","pop.ef","G","E","theta")
	return(out)
	}

diwish_ <-
function(x, v, S, log=F){
	f = diwish(x, v, S)
	if( log ){
		f = log(f)
		}
	f = max(c(f, -1.0e280))
	return(f)
	}

dwish_ <-
function(x, v, S, log=F){ # non-normalized!
	S = as.matrix(S)
	x = as.matrix(x)
	S_ = pseudoinverse(S) # to avoid numerical singularity
	p = ncol(x)
	f = rep(NA,3)
	f[1] = 0.5*(v-p-1)*log(det(x))
	f[2] = 0.5*sum(diag(S_%*%x))
	f[3] = 0.5*v*log(det(S))
	f[which(f>1.0e280)] = 1.0e280
	f[which(f<(-1.0e280))] = -1.0e280
	f = f[1] - f[2] - f[3]
	if( !log ){
		f = exp(f)
		}
	return(f)
	}

ellipsis <-
function(mu, Sig, prob.mass, lwd, col){
	radius <- sqrt(qchisq(prob.mass, 2))
	theta <- 2*pi*(0:360) / 360
	unit.circle <- cbind(cos(theta), sin(theta))
	Q <- base::chol(Sig, pivot=TRUE)
	order <- order(attr(Q, "pivot"))
	ellipse <- t(mu + radius * t(unit.circle %*% Q[, order]))
	lines(ellipse[,1], ellipse[,2], lwd=lwd, col=col)
	}

H.test <-
function(popefpost, Gpost, THpost, env, silent=T){

	# Stupid scaling again
	nx = dim(popefpost)[1] * dim(popefpost)[2]
	npop = dim(THpost)[1]
	nmc = dim(THpost)[3]
  	env = env[,-1] # removing ids
	env = into.dist(env)

	# Posterior of Mantel statistic for phenotypic distance
	mobs = mrand = rep(NA, nmc)
	for( i in 1:nmc ){
		G = Gpost[,,i] # G matrix
		theta = THpost[,,i] # theta	
		C2 = matrix(mvrnorm( 1, rep(0, nx), 2*G %x% theta ), nrow=npop)
		D2 = into.dist(popefpost[,,i]) # observed distance matrix		
		C2 = into.dist(C2) # distance matrix from neutral model
		mobs[i] = mantel.stat(D2, env) # observed product moment
		mrand[i] = mantel.stat(C2, env) # randomized product moment
		}
	
	# Fancy figures	
	if( !silent ){
  	par(mfcol=c(1,2))
  	plot(mobs, type='l', ylab="observed product moment", xlab="iteration", main="Convergence check")
  	plot(mrand, type='l', ylab="product moment from neutral model", xlab="iteration", main="Convergence check")
  	}

 	out = length(which(mobs>mrand)) / length(mobs)
	return(out)
	}

ind.theta <-
function(pars, censped){
	n = nrow(pars)
	th = matrix(0, nrow=n, ncol=n)
	for( i in 1:n ){
		for( j in 1:i ){
			x = 0
			if( i==j ){ x = 0.5 }
			else{
				if( (censped[i]+censped[j])==0 ){
					x = 0.125*( 1*(pars[i,1]==pars[j,1]) + 1*(pars[i,2]==pars[j,2]) )
					}
				}	
			th[i,j] = th[j,i] = x
			}
		}
	return(th)
	}

into.dist <-
function(covars){
	if( is.vector(covars) ){ covars = as.matrix(covars) }
	npop = nrow(covars)
	D = matrix(NA, npop, npop)
	for( j in 1:npop ){
		for( k in 1:j ){
			D[j,k] = sqrt(sum( (covars[k,] - covars[j,])^2 ))
			D[k,j] = D[j,k]
			}
		}	
	return(D)
	}

likelihood <-
function(y, g, W, cSig, sparse){
	mu = W %*% g
	f = mvdnorm.chol(y, mu, cSig, sparse, log=T)
	return(f)
	}

mantel.stat <-
function(A, B){
	D = A*B
	out = 0
	for( i in 1:nrow(D) ){
		for( j in 1:i ){
			out = out + D[i,j] } }
	return(out) }

MH <-
function(poster, ped, covars, traits, nmc, burnin, thin, blocks=NA, priors=NA, tmp=NA, binary=NA, alt=F){

	# Dependencies
	library(corpcor)
	library(MCMCpack)
	library(SparseM)

	# Automate conversion
	ped <- as.matrix(ped)
	covars <- as.matrix(covars)
	traits <- as.matrix(traits)
	blocks <- as.matrix(blocks)
	
	# Phenotypic data
	X <- as.matrix(covars[,-1], ncol=ncol(covars)-1)
	X <- cbind(1, X) # 1 for grand mean
	y <- c(traits[,-1])
	N <- length(y)
	cens <- which(is.na(y))
	if( length(cens)>0 ){ y <- y[-cens] }
	if( N>=200 ){ sparse <- T } # sparse matrix tricks
	else{ sparse <- F }
	
	# Coancestry matrix & the other 'Kroenecker matrix'
	nop <- dim(poster)[1]
	if( max(poster)==0 ){ # Bayes spirit
    		for( i in 1:(dim(poster)[3]) ){
      			poster[,,i] <- 10^(-4)*diag(nop)
      			}
    		}
	w <- sample(1:(dim(poster)[3]), 1)
	ww <- w
	THP <- as.matrix(poster[,,w])
	if( alt ){ # one generation
		ped <- ped[,-1]
		censped <- apply(ped, 1, function(x)length(which(is.na(x))))		
		THB <- ind.theta(ped, censped)
		M <- calc.M(ped, THP, THB, censped)
		}
	else{
		wpop <- calc.loads(ped, THP)
		M <- calc.M.alt(wpop, THP, ped)
		}
	VR <- diag(ncol(M))
	if( sparse ){
		M <- as.matrix.csr(M)
		VR <- as.matrix.csr(VR) 
		}

	# Fundamental covariance matrices
	traits <- as.matrix(traits[,-1], ncol=ncol(traits)-1)
	ntr <- ncol(traits)
	G <- rwish(ntr, diag(ntr)/ntr) # rwish_ doesn't work yet
	E <- rwish(ntr, diag(ntr)/ntr)

	# Treating random/block effects
	bbit <- (length(blocks) > 1)
	if( bbit ){
		nb <- ncol(blocks) - 1
		blocks <- calc.blocks(blocks) # covariates into design matrices
		Blist <- array(NA, dim=c(ntr, ntr, nb))
		for( i in 1:nb ){
			Blist[,,i] <- rwish(ntr, diag(ntr)/ntr)
			}
		}
	else{ Blist <- NA }
	
	# Covariance matrix
	Sig <- const.sig(E, VR, G, M, blocks, Blist, bbit, sparse)
	if( length(cens)>0 ){ Sig <- Sig[-cens,-cens] }
	if( sparse ){ 
		Sig <- as.matrix.csr(Sig)
		if( is.na(tmp) ){ cSig <- SparseM::chol(Sig) }
		else{ cSig <- SparseM::chol(Sig, tmpmax=tmp) }
		}
	else{ cSig <- base::chol(Sig) }
	
	# Default priors
	if( length(priors)==1 ){
		priors <- as.list(1:4)
		names(priors) <- c("fixed","G","E","random")
		ng <- ncol(traits)*(ncol(THP)+ncol(X))
		priors[[1]] <- list(rep(0,ng), matrix(NA,ncol=1,nrow=1)) # for g
		priors[[2]] <- list(ntr+1, diag(ntr)/(ntr+1)) # for G
		priors[[3]] <- list(ntr+1, diag(ntr)/(ntr+1)) # for E
		if( bbit ){
			priors[[4]] <- as.list(1:nb)
			for( i in 1:nb ){ # for blocks stuff
				priors[[4]][[i]] <- list(ntr+1, diag(ntr)/(ntr+1))
				}
			}
		}
	else{
		priors <- priors
		names(priors) <- c("fixed","G","E","random")
		priors$fixed[[2]] <- as.matrix(priors$fixed[[2]])
		}
		
	# Design matrices
	g <- rnorm(ncol(traits)*(ncol(THP)+ncol(X)), 0, sqrt(max(poster)))
	V <- calc.V(g, THP, G, priors)
	if( alt ){ W <- calc.W(ped, X, traits, THP) }
	else{ W <- calc.W.alt(wpop, X, traits, THP) }
	if( length(cens)>0 ){ W <- W[-cens,] }			

	# Treating non-normal variables
	if( !is.na(binary) ){
		y_ <- y
		y <- update.y(y_, NA, g, W, Sig, binary, ntr, first=T)
		}
	
	# Densities
	like <- likelihood(y, g, W, cSig, sparse)
	priG <- pri.GG(G, priors)
	priE <- pri.E(E, priors)
	prig <- pri.g(g, THP, G, sparse, priors)
	if( bbit ){
		priBlocks <- sapply(1:nb, function(x) pri.B(Blist[,,x], x, priors))
		}
		
	output <- c(g,G,E,THP)
	if( bbit ){
		out2 <- c(Blist, recursive=T)
		}

	# Adjustment parameters
	dltG <- dltE <- 2*ncol(G)
	dltb <- 1
	iterno <- Eup <- Gup <- bup <- 0
	if( bbit ){ 
		dltB <- rep(2*ncol(G), nb) 
		Bup <- rep(0, nb) # not to be confused with bup
		}
	else{ dltB <- 10 }
	props <- rep( 1/(dim(poster)[3]), dim(poster)[3] )
	

	# A MASSIVE MCMC LOOP starts -
	for( i in 1:nmc ){

		# Updating g, nicely packed
		V <- calc.V(g, THP, G, priors)
		try(g <- update.g(W, cSig, V, y, sparse, priors), silent=T) # numerical singularity causes fail
		like <- likelihood(y, g, W, cSig, sparse)
		prig <- pri.g(g, THP, G, sparse, priors)

		# Updating latent liabilities, if any
		if( !is.na(binary) ){ 
			y <- update.y(y_, y, g, W, Sig, binary, ntr, first=F) 
			like <- likelihood(y, g, W, cSig, sparse)
			}

		# Updating G, random walk
		newG <- G
		try(newG <- rwish_(dltG, G/dltG), silent=T)
		here <- dwish_(G, dltG, newG/dltG, log=T)
		there <- dwish_(newG, dltG, G/dltG, log=T)
		newpriG <- pri.GG(newG, priors)	
		newprig <- pri.g(g, THP, newG, sparse, priors)
		newSig <- const.sig(E, VR, newG, M, blocks, Blist, bbit, sparse)
		if( length(cens)>0 ){ newSig <- newSig[-cens,-cens] }
		if( sparse ){ 
			newSig <- as.matrix.csr(newSig)
			if( is.na(tmp) ){ newcSig <- SparseM::chol(newSig) }
			else{ newcSig <- SparseM::chol(newSig, tmpmax=tmp) }
			}
		else{ newcSig <- base::chol(newSig) }
		newlike <- likelihood(y, g, W, newcSig, sparse)
		accept <- newlike + newpriG + newprig + here - like - priG - prig - there
		a <- log(runif(1,0,1))
		if( !is.na(accept) ){
		if( a < accept ){
			G <- newG
			Sig <- newSig
			cSig <- newcSig
			like <- newlike
			priG <- newpriG
			prig <- newprig
			Gup <- Gup + 1
			}}

		# Updating E, random walk
		newE <- E
		try(newE <- rwish_(dltE, E/dltE), silent=T)
		here <- dwish_(E, dltE, newE/dltE, log=T)
		there <- dwish_(newE, dltE, E/dltE, log=T)
		newpriE <- pri.E(newE, priors)
		newSig <- const.sig(newE, VR, G, M, blocks, Blist, bbit, sparse)
		if( length(cens)>0 ){ newSig <- newSig[-cens,-cens] }
		if( sparse ){ 
			newSig <- as.matrix.csr(newSig)
			if( is.na(tmp) ){ newcSig <- SparseM::chol(newSig) }
			else{ newcSig <- SparseM::chol(newSig, tmpmax=tmp) }
			}
		else{ newcSig <- base::chol(newSig) }
		newlike <- likelihood(y, g, W, newcSig, sparse)
		accept <- newlike + newpriE + here - like - priE - there
		a <- log(runif(1,0,1))
		if( !is.na(accept) ){
		if( a < accept ){
			E <- newE
			Sig <- newSig
			cSig <- newcSig
			like <- newlike
			priE <- newpriE
			Eup <- Eup + 1
			}}
			
		# Updating theta, independence sampler
		if( runif(1,0,1) < 0.5 ){
			k <- sample(1:(dim(poster)[3]), 1, prob=props)
			}
		else{ k <- sample(1:(dim(poster)[3]), 1) }
		newTHP <- as.matrix(poster[,,k])
		if( alt ){ newM <- calc.M(ped, newTHP, THB, censped) }
		else{ newM <- calc.M.alt(wpop, newTHP, ped) }
		if( sparse ){ newM <- as.matrix.csr(newM) }
		newSig <- const.sig(E, VR, G, newM, blocks, Blist, bbit, sparse)
		if( length(cens)>0 ){ newSig <- newSig[-cens,-cens] }
		if( sparse ){ 
			newSig <- as.matrix.csr(newSig)
			if( is.na(tmp) ){ newcSig <- SparseM::chol(newSig) }
			else{ newcSig <- SparseM::chol(newSig, tmpmax=tmp) }
			}
		else{ newcSig <- base::chol(newSig) }
		newprig <- pri.g(g, newTHP, G, sparse, priors)
		newlike <- likelihood(y, g, W, newcSig, sparse)
		here <- log(props[w] + 0.5)
		there <- log(props[k] + 0.5)
		accept <- newlike + newprig - like - prig + here - there	
		a <- log(runif(1,0,1))
		if( !is.na(accept) ){
		if( a < accept ){	
			THP <- newTHP
			M <- newM
			Sig <- newSig
			cSig <- newcSig
			prig <- newprig
			like <- newlike
			w <- k
			}}
			
		# Updating block/random effects, random walk
		if( bbit ){
			for( k in 1:nb ){
				newB <- Blist
				try(newB[,,k] <- rwish_(dltB[k], newB[,,k,drop=F]/dltB[k]), silent=T)
				here <- dwish_(Blist[,,k], dltB[k], newB[,,k]/dltB[k], log=T)
				there <- dwish_(newB[,,k], dltB[k], Blist[,,k]/dltB[k], log=T)
				newpri <- pri.B(newB[,,k], k, priors)
				newSig <- const.sig(E, VR, G, M, blocks, newB, bbit, sparse)
				if( length(cens)>0 ){ newSig <- newSig[-cens,-cens] }
				if( sparse ){ 
					newSig <- as.matrix.csr(newSig)
					if( is.na(tmp) ){ newcSig <- SparseM::chol(newSig) }
					else{ newcSig <- SparseM::chol(newSig, tmpmax=tmp) }
					}
				else{ newcSig <- base::chol(newSig) }
				newlike <- likelihood(y, g, W, newcSig, sparse)
				accept <- newlike + newpri + here - like - priBlocks[k] - there
				a <- log(runif(1,0,1))
				if( !is.na(accept) ){
				if( a < accept ){
					Blist <- newB
					Sig <- newSig
					cSig <- newcSig
					like <- newlike
					priBlocks[k] <- newpri
					Bup[k] <- Bup[k] + 1
					}}
				}		
			}

		# Adapting proposals of G, E, etc.
		iterno <- iterno + 1
		if( i < burnin ){
			wgh <- 1 - 0.1*exp(-i/1000)
			a <- 0.1*exp(-i/1000)
			iterno <- wgh*iterno	
			Gup <- wgh*Gup
			Eup <- wgh*Eup
			rateG <- Gup / iterno
			rateE <- Eup / iterno
			dltG_ <- exp(a*(0.23 - rateG))*dltG
			dltE_ <- exp(a*(0.23 - rateE))*dltE
			if( dltG_ > ntr ){
			if( dltG_ < 200 ){
				dltG <- dltG_ } }
			if( dltE_ > ntr ){
			if( dltE_ < 200 ){
				dltE <- dltE_ } }
			if( bbit ){
				for( i in 1:length(dltB) ){
					Bup[i] <- wgh*Bup[i]
					rateB <- Bup[i] / iterno
					dltB_ <- exp(a*(0.23 - rateB))*dltB[i]
					  if( dltB_ > ntr ){
							if( dltB_ < 200 ){
								dltB[i] <- dltB_ }}
					}
				}
			}
			
		# Adapting proposal of theta
		ww <- c(ww, w)
		if( i==burnin ){ props <- study(ww, dim(poster)[3], maxpar=30) }
			
		# Info for user
		if( round(i/100)==i/100 ){ print(paste("iter",i), quote=F) }
		output <- rbind(output, c(g,G,E,THP))
		if( bbit ){
			out2 <- rbind(out2, c(Blist, recursive=T))
			}
		} # MC LOOP CLOSES

		
	# Thinning and burn-in
	output <- output[(burnin+1):nrow(output),]
	imax <- floor(nrow(output) / thin)
	totake <- thin*1:imax
	if( length(totake>1) ){
		output <- output[totake,]
		}	
	if( bbit ){
		out2 <- out2[(burnin+1):nrow(out2),,drop=F]
		if( length(totake>1) ){
			out2 <- out2[totake,,drop=F]
			}
		}
		
	# Conversion
	out <- convert(output, ncol(traits), ncol(THP), ncol(X))
	if( bbit ){ 
		out2 <- convert.2(out2, ncol(G), nb)
		out <- list(out, out2)
		names(out) <- c("main.result","block.ef")
		}
	
	return(out)	
	}

mvdnorm.chol <-
function(x, mu, cSig, sparse, log=F){
	dif = x - mu
	if( sparse ){ 
		SS = t(dif) %*% SparseM::backsolve(cSig, dif)
		}
	else{ SS = t(dif) %*% xx(cSig, dif) }
	SS = c(as.matrix(SS))
	if( sparse ){
		f1 = -(cSig@log.det)
		}
	else{ f1 = - sum(log(diag(cSig))) }
	f1 = min(c( f1, 1.0e280 ))
	SS = max(c( SS, -1.0e280 ))
	f1 = max(c( f1, -1.0e280 ))
	f = f1 - 0.5*SS
	if( !log ){
		f = exp(f)
		}
	return(f)
	}

mvdnorm <-
function(x, mu, Sig, log=F, invert=T, sparse=F){
	if( sparse ){
		cSig = SparseM::chol(Sig)
		f = mvdnorm.chol(x, mu, cSig, TRUE, log=T)
		}
	else{
		Sig_ = solve(Sig)
		SS = t(x - mu) %*% Sig_ %*% (x - mu)
		f1 = -0.5*log(det(Sig))
		}
	f1 = min(c( f1, 1.0e280 ))
	SS = max(c( SS, -1.0e280 ))
	f1 = max(c( f1, -1.0e280 ))
	f = f1 - 0.5*SS	
	if( !log ){
		f = exp(f)
		}
	return(f)
	}

neut.test <-
function(popefpost, Gpost, THpost, silent=T, G.off=F, th.off=F, main=NA){

	# Turning off-diagonals off if needed
	if( th.off ){
		npop = dim(THpost)[1]
		if( npop>1 ){
			for( i in 2:npop ){
				for( j in 1:(i-1) ){
					THpost[i,j,] = THpost[j,i,] = 0 }}}}
	if( G.off ){
		ntrait = dim(Gpost)[1]
		if( ntrait > 1 ){
			for( i in 2:ntrait ){
				for( j in 1:(i-1) ){
					Gpost[i,j,] = Gpost[j,i,] = 0 }}}}

	# Mahalanobis!
	apost = popefpost
	nmc = dim(apost)[3]
	D = rep(NA, nmc)
	for( i in 1:nmc ){
		a = c(apost[,,i])
		G = Gpost[,,i]
		THP = THpost[,,i]
		Sig = 2*G %x% THP
		Sig_ = solve(Sig)
		mu = rep(0, length(a))
		D2 = t(mu - a) %*% Sig_ %*% (mu - a)
		D[i] = D2
		}
	
	# Posterior probabilities
	dfr = (dim(Gpost)[1])*(dim(apost)[1])
	cdf = pchisq(D, dfr)
	out = mean(cdf)
	if( !silent ){
		if( is.na(main) ){ main = "Posterior trace" }
		plot(cdf, xlab="iteration", ylab="signal of selection (S)", main=main, type='l', lwd=2, ylim=c(0,1))
		lines(c(-length(cdf),2*length(cdf)), rep(0.2,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.8,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.5,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.05,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.95,2), lty='dashed')		
		}
	return(out)
	}

Nh <-
function(invSig, mu_Sig, sparse){
	C = base::chol(invSig)
	u = rnorm(length(mu_Sig))
	gg = solve(t(C), mu_Sig)
	fi = solve(C, gg)
	u_ = solve(C, u)
	x = fi + u_
	return(x)
	}

optim.beta <-
function(gridx, gridy, maxpar=6){
	eps = (gridx[2] - gridx[1]) / 2
	grdim = maxpar*10
	arange = brange = (1:grdim) / 10
	errormat = matrix(NA, ncol=grdim, nrow=grdim)
	for( i in 1:grdim ){
		for( j in 1:grdim ){
			a = arange[i]
			b = brange[j]
			f = sapply( gridx, function(x) (pbeta(x+eps, a, b) - pbeta(x-eps, a, b)) )
			f = f / sum(f)
			errormat[i,j] = sum((f - gridy)^2)
			}
		}
	num = which(errormat==min(errormat))
	ai = as.integer(grdim*(num/grdim - floor(num/grdim)))
	if( ai==0 ){ ai = grdim }
	bi = ceiling(num/grdim)
	a = arange[ai]
	b = brange[bi]
	f = sapply( gridx, function(x) (pbeta(x+eps, a, b) - pbeta(x-eps, a, b)) )
	f = f / sum(f)
	out = c(a, b)
	return(out)
	}

pri.B <-
function(B, k, priors){
	prithis = priors$random[[k]]
	f = dwish_(B, prithis[[1]], prithis[[2]], log=T)
	return(f)
	}

pri.E <-
function(E, priors){
	f = dwish_(E, priors$E[[1]], priors$E[[2]], log=T)
	return(f)
	}

pri.g <-
function(g, THP, G, spar, priors){
	V = calc.V(g, THP, G, priors)
	f <- NA # unavoidable psd fix
	try(f <- mvdnorm(g, priors$fixed[[1]], V, log=T, sparse=F), silent=T)
	return(f)
	}

pri.GG <-
function(G, priors){
	f = dwish_(G, priors$G[[1]], priors$G[[2]], log=T)
	return(f)
	}

reord <-
function(traits, binary){
	if( is.na(binary) ){
		out = traits
		}
	else{
		if( length(binary) < (ncol(traits) - 1) ){
			id = traits[,1]
			traits = traits[,-1]
			guide = 1:ncol(traits)
			nonbin = guide[-binary]
			out = traits[,c(nonbin, binary)]
			out = cbind(id, out)
			}
		else{
			out = traits
			}
		}    
   return(out)
   }

riwish_ <-
function(v, S){
	sing = TRUE
	while( sing ){
		X = riwish(v, S)
		try(sing <- check.check(X) , silent=T)
		}
	return(X)
	}

rwish_ <-
function(v, S){
	limes = 200 # hard-coded maximal degree of freedom
	sing = TRUE
	while( sing ){
		X = rwish(v, S)
		try(sing <- check.check(X, limes) , silent=T)
		}
	X = 0.5*(X + t(X))
	return(X)
	}

S.test <-
function(popefpost, Gpost, THpost, silent=T, G.off=F, th.off=F, main=NA){

	# Turning off-diagonals off if needed
	if( th.off ){
		npop = dim(THpost)[1]
		if( npop>1 ){
			for( i in 2:npop ){
				for( j in 1:(i-1) ){
					THpost[i,j,] = THpost[j,i,] = 0 }}}}
	if( G.off ){
		ntrait = dim(Gpost)[1]
		if( ntrait > 1 ){
			for( i in 2:ntrait ){
				for( j in 1:(i-1) ){
					Gpost[i,j,] = Gpost[j,i,] = 0 }}}}

	# Mahalanobis!
	apost = popefpost
	nmc = dim(apost)[3]
	D = rep(NA, nmc)
	for( i in 1:nmc ){
		a = c(apost[,,i])
		G = Gpost[,,i]
		THP = THpost[,,i]
		Sig = 2*G %x% THP
		Sig_ = solve(Sig)
		mu = rep(0, length(a))
		D2 = t(mu - a) %*% Sig_ %*% (mu - a)
		D[i] = D2
		}
	
	# Posterior probabilities
	dfr = (dim(Gpost)[1])*(dim(apost)[1])
	cdf = pchisq(D, dfr)
	out = mean(cdf)
	if( !silent ){
		if( is.na(main) ){ main = "Posterior trace" }
		plot(cdf, xlab="iteration", ylab="signal of selection (S)", main=main, type='l', lwd=2, ylim=c(0,1))
		lines(c(-length(cdf),2*length(cdf)), rep(0.2,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.8,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.5,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.05,2), lty='dashed')
		lines(c(-length(cdf),2*length(cdf)), rep(0.95,2), lty='dashed')		
		}
	return(out)
	}

scale.pop.ef <-
function(posterior){
	nmc = dim(posterior)[3]
	phenodist = as.list(1:nmc)
	npop = dim(posterior)[1]
	ntrait = dim(posterior)[2]
	nc = npop*ntrait
	phenodist = as.list( rep(NA,nmc) )
	for( i in 1:nmc ){
		D = matrix(NA, npop, npop)
		pheno = posterior[,,i]
		phenodist[[i]] = into.dist(pheno)
		}
	return(phenodist)
	}

study.2 <-
function(whistory, maxk, maxpar=6){
	whistory = round(whistory/maxk, 1)
	gridx = (1:9)/10
	gridy = sapply(gridx, function(x) length(which(whistory==x)))
	gridy = gridy / sum(gridy)
	parvec = optim.beta(gridx, gridy, maxpar=maxpar)
	gridprop = (1:maxk) / (maxk+1)
	props = sapply(gridprop, function(x) dbeta(x, parvec[1], parvec[2]))
	props = props / sum(props)
	return(props)
	}

study <-
function(whistory, maxk, maxpar=30){
	y = whistory / (maxk+1)
	arange = brange = (10:(10*maxpar)) / 10
	grdim = length(arange)
	likemat = matrix(NA, ncol=grdim, nrow=grdim)
	for( i in 1:grdim ){
		for( j in 1:grdim ){
			a = arange[i]
			b = brange[j]
			likemat[i,j] = sum(dbeta(y, a, b, log=T))
			}
		}
	num = which(likemat==max(likemat))
	ai = as.integer(grdim*(num/grdim - floor(num/grdim)))
	if( ai==0 ){ ai = grdim }
	bi = ceiling(num/grdim)
	a = arange[ai]
	b = brange[bi]
	gridprop = (1:maxk) / (maxk+1)
	props = sapply(gridprop, function(x) dbeta(x, a, b))
	props = props / sum(props)
	return(props)
	}

trunc1 <-
function(U, mux, sdx, signx){
	gridx = (1:9999) / 10^4
	if( signx==-1 ){
		up = pnorm(0, mux, sdx)
		gridx = gridx*up
		}
	if( signx==1 ){
		lo = 1 - pnorm(0, mux, sdx)
		gridx = 1 - lo*gridx
		}
	vals = qnorm(gridx, mux, sdx)
	if( max(vals)==-Inf ){ vals = rep(-1.0e-07, 10^4) }
	else{ vals[which(vals==-Inf)] = min(vals[which(vals>-Inf)]) }
	if( min(vals)==Inf ){ vals = rep(1.0e-07, 10^4) }
	else{ vals[which(vals==Inf)] = max(vals[which(vals<Inf)]) }
	U = ceiling(9999*U)
	out = vals[U]
	return(out)
	}

truncrnorm <-
function(n, mux, sdx, signx){ # inefficient, calculates bins n times; yet redundant
	out = sapply(1:n, function(x) trunc1(runif(1,0,1), mux, sdx, signx) )
	return(out)
	}

update.g <-
function(W, cSig, V, y, sparse, priors){		
	V_ = solve(V)
	if( sparse ){
		gam_ = (t(W) %*% SparseM::backsolve(cSig,W)) + V_
		gam_ = as.matrix(gam_)
		mu_sig = t(W) %*% SparseM::backsolve(cSig,matrix(y,ncol=1))
		mu_sig = mu_sig + V_ %*% priors$fixed[[1]]
		mu_sig = as.matrix(mu_sig)
		}
	else{
		gam_ = (t(W) %*% xx(cSig,W)) + V_
		mu_sig = t(W) %*% xx(cSig,matrix(y,ncol=1))
		mu_sig = mu_sig + V_ %*% priors$fixed[[1]]
		}
	g = Nh(gam_, mu_sig, sparse)
	return(g)
	}

update.y <-
function(y_, y, g, W, Sig, binary, ntr, first=FALSE){

	# Basics
	mu <- W %*% g
	nin <- length(y) / ntr
	bins <- NA
	for( b in binary ){
		binb <- (1:nin) + (b-1)*nin
		bins <- c(bins, binb)
		}
	bins <- bins[-1]
	nonbins <- (1:length(y))[-bins]

	# First iteration
	if( first ){
		y <- y_
		for( i in bins ){
			if( y_[i]==1 ){ y[i] <- 0.01 }
			if( y_[i]==0 ){ y[i] <- -0.01 }
			}
		}
	
	# Other iterations
	else{
		y_[which(y_==0)] <- -1 # sign for normal truncation
		y <- y
		V <- solve(Sig)
		V <- as.matrix(V)
		Sig <- as.matrix(Sig)
		for( j in bins ){
			B <- -V[j,-j] / V[j,j]
			mu.j <- mu[j] + B %*% (y[-j] - mu[-j])
			s2.j <- Sig[j,j] - B %*% Sig[-j,j]
			if( s2.j <= 0 ){ s2.j <- 1.0e-04 }
			y[j] <- truncrnorm(1, mu.j, sqrt(s2.j), y_[j])
			}
		}		
	return(y)
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

viz.traits <-
function(fixedpost, popefpost, Gpost, THpost, traits, siz=0.5, main=NA){

    # Posterior means
    fixedpost = fixedpost[1,traits,] # ancestral mean
    popefpost = popefpost[,traits,]
    Gpost = Gpost[traits,traits,]
    popef = matrix(NA, nrow=dim(popefpost)[1], ncol=dim(popefpost)[2])
    G = matrix(NA, nrow=dim(Gpost)[1], ncol=dim(Gpost)[2])
    TH = matrix(NA, nrow=dim(THpost)[1], ncol=dim(THpost)[2])
    mu = rep(NA, dim(Gpost)[1])
    for( i in 1:length(traits) ){
		mu[i] = mean(fixedpost[i,])
		for( j in 1:length(traits) ){ G[i,j] = mean(Gpost[i,j,]) }
		for( j in 1:dim(popefpost)[1] ){ popef[j,i] = mean(popefpost[j,i,]) } # this way
        }
    for( i in 1:nrow(TH) ){
		for( j in 1:ncol(TH) ){ TH[i,j] = mean(THpost[i,j,]) }
		}
                
    # Drawing
	npop = ncol(TH)
	ntr = length(traits)
	npop = nrow(popef)
	colvec = 1:9
	colvec[7] = "orange"
	colvec[9] = "chartreuse"		
	i = 1
	j = 2
	ei = mu[i]
	ej = mu[j]
	Gthis = matrix(c(G[i,i],G[i,j],G[j,i],G[j,j]), ncol=2)
	xlab = paste("trait", traits[i])
	ylab = paste("trait", traits[j])
	M = max(c( sqrt(2*G[i,i]*diag(TH)) , sqrt(2*G[j,j]*diag(TH)), popef ))
	plot(ei, ej, pch=16, cex=0, xlim=c(ei-M,ei+M), ylim=c(ej-M,ej+M), xlab=xlab, ylab=ylab, main=main)
	for( k in 1:npop ){
		mu = c(ei,ej)
		Sigma = 2*TH[k,k]*Gthis
		if( k > 9 ){ k = 1 + k%%8 }
		ellipsis(mu, Sigma, siz, 1, colvec[k])
		}
	text(ei+popef[,i], ej+popef[,j], 1:npop, cex=1, col=colvec)
	text(ei, ej, "A")
	}

xx <-
function(cA, B){
	nc = ncol(B)
	product = matrix(NA, nrow=nrow(cA), ncol=ncol(B))
	for( j in 1:nc ){
		yep = solve(t(cA), B[,j])
		product[,j] = solve(cA, yep)
		}
	return(product)
	}
