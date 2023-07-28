# Functions to be repeatedly used

library(mvtnorm)

# Fitness (multivariate Gaussian fitness function w/o covariance)
# Calculate the overall fitness given a set of traits' values (d, a vector) and SDs of their respective fitness functions (a, also a vector) 
# Overall fitness calculated as function of Euclidean distance from the global optimum; each trait's fitness function's SD becomes a scaling coefficient (remains mathematically equivalent to original definition)
fitness <- function(d,a){
	# Remove traits that do not affect fitness
	ds=d[which(a>0)]
	as=a[which(a>0)]
	ds=ds/as # Rescale the phenotypes; when the fitness function is narrow, the distance from optimum is rescaled to be larger
	D=sqrt(sum(ds^2)) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits
	w=dnorm(D,mean=0,sd=1)/dnorm(0,mean=0,sd=1)
	return(w)
}

# Fitness (multivariate Gaussian fitness function w/ covariance)
# Parameters: multivariate phenotype (a vector), optimal phenotype (a vector), matrix of selection gradient (a matrix)
# Zero eigenvalues of beta correspond to dimensions with no selection
fitness <- function(d,beta){
	dnew=as.numeric(d%*%eigen(beta)$vectors) # Convert to value of orthogonal traits
	ds=dnew[which(eigen(beta)$values>0)] # Remove dimensions under no selection
	ds=ds/eigen(beta)$values[which(eigen(beta)$values>0)] # Rescale the phenotypes
	D=sqrt(sum(ds^2)) # Overall fitness calculated as geometric mean of fitness values calculated from orthogonal traits
	w=dnorm(D,mean=0,sd=1)/dnorm(0,mean=0,sd=1)
	return(w)
}
# Sanity check: whether two fitness calculation functions give the same result when there is no covariance in beta; need to name two functions differently before testing
fitness.uv(c(0.1,0.1),c(1,1))
fitness.mv(c(0.1,0.1),rbind(c(1,0),c(0,1)))
# Sanity check: whether covariance in beta affects fitness
fitness(c(0.1,0.1),rbind(c(1,0),c(0,1)))
fitness(c(0.1,0.1),rbind(c(1,0.5),c(0.5,1)))
fitness(c(0.1,0.1),rbind(c(1,-0.5),c(-0.5,1)))

# Fixation probability
# Calculate fixation probability of a mutation given ancestral fitness, mutant fitness, and effective population size (Ne)
fix.prob <- function(wa,wm,Ne){
	if(wa>0){
		s=(wm/wa)-1 # Coefficient of selection
		if(s==0){
			p=1/(2*Ne)
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
		}
	}else{ # The ancestral fitness is close to 0 (close enough to be recognized as zero by R)
		if(wm==0){ # Both ancestral and mutant fitness are close to 0
			p=1/(2*Ne) # The mutation is considered neutral
		}else{ # Ancestral fitness is close to 0 while mutant fitness isn't
			p=1 # The mutation is considered strongly beneficial
		}
	}
	return(p)
}
# Test/sanity check
wa=fitness(c(1e-4,1e-4),rbind(c(1,0),c(0,1)))
wm=fitness(c(1e-3,1e-3),rbind(c(1,0),c(0,1)))
fix.prob(wa,wm,1e3)
fix.prob(wa,wm,1e4)

# Generate mutants to obtain a distribution of fitness effect and fixation probability (WT background)
ma <- function(par,n){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	sig.np=c(sig1,sig2)
	mp=rbind(c(siga^2,r*siga*sigb),c(r*siga*sigb,sigb^2))
	
	out=matrix(0,nrow=n,ncol=2)
	type=sample(1:3,n,prob=c(U1,U2,Up)/(U1+U2+Up),replace=TRUE) # Type of each mutation
	for(i in 1:n){
		if(type[i]<3){ # The mutation is not pleiotropic
			out[i,type[i]]=rnorm(1,mean=0,sd=sig.np[type[i]])
		}else{
			out[i,]=as.numeric(rmvnorm(1,sigma=mp))
		}
	}
	return(out)
}

# Simulate a single lineage and write phenotype thru time (no correlational selection)
# Parameters: mutational parameters (vector of length 8), optimal phenotype (vector of length 2), SDs of fitness functions (vector of length 2), Ne, time of simulation
# a=c(0,0) for neutral evolution, opt!=c(0,0) for directional selection
sim <- function(par,opt,a,Ne,T){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	sig.np=c(sig1,sig2)
	mp=rbind(c(siga^2,r*siga*sigb),c(r*siga*sigb,sigb^2))

	phe=matrix(0,nrow=2,ncol=(T+1))
	for(t in 2:(T+1)){
		phe[,t]=phe[,(t-1)]
		nmut=rpois(1,lambda=(U1+U2+Up)) # Total number of mutations that would occur in this time step
		if(nmut>=1){
			type=sample(1:3,nmut,prob=c(U1,U2,Up)/(U1+U2+Up),replace=TRUE) # Type of each mutation
			for(i in 1:nmut){
				phe_mutant=phe[,t]
				if(type[i]<3){ # The mutation is not pleiotropic
					effect=rnorm(1,mean=0,sd=sig.np[type[i]])
					phe_mutant[type[i]]=phe_mutant[type[i]]+effect # Calculate mutant phenotype
				}else{ # The mutation is pleiotropic
					effect=rmvnorm(1,sigma=mp)
					phe_mutant=phe_mutant+as.numeric(effect) # Calculate mutant phenotype
				}
				wa=fitness(phe[,t]-opt,a);wm=fitness(phe_mutant-opt,a)
				fp=fix.prob(wa,wm,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					phe[,t]=phe_mutant # If the mutation is fixed, add its effect onto the population mean before the next mutation is examined
				}
			}
		}
	}
	return(phe)
}

# Simulate a single lineage and write phenotype thru time (w/ correlational selection)
# Parameters: mutational parameters (vector of length 8), optimal phenotype (vector of length 2), matrix of selection gradient (2x2 matrix), Ne, time of simulation
# Zero eigenvalues of beta correspond to dimensions with no selection
sim <- function(par,opt,beta,Ne,T){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	sig.np=c(sig1,sig2)
	mp=rbind(c(siga^2,r*siga*sigb),c(r*siga*sigb,sigb^2))

	phe=matrix(0,nrow=2,ncol=(T+1))
	for(t in 2:(T+1)){
		phe[,t]=phe[,(t-1)]
		nmut=rpois(1,lambda=(U1+U2+Up)) # Total number of mutations that would occur in this time step
		if(nmut>=1){
			type=sample(1:3,nmut,prob=c(U1,U2,Up)/(U1+U2+Up),replace=TRUE) # Type of each mutation
			for(i in 1:nmut){
				phe_mutant=phe[,t]
				if(type[i]<3){ # The mutation is not pleiotropic
					effect=rnorm(1,mean=0,sd=sig.np[type[i]])
					phe_mutant[type[i]]=phe_mutant[type[i]]+effect # Calculate mutant phenotype
				}else{ # The mutation is pleiotropic
					effect=rmvnorm(1,sigma=mp)
					phe_mutant=phe_mutant+as.numeric(effect) # Calculate mutant phenotype
				}
				wa=fitness(phe[,t]-opt,beta);wm=fitness(phe_mutant-opt,beta)
				fp=fix.prob(wa,wm,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					phe[,t]=phe_mutant # If the mutation is fixed, add its effect onto the population mean before the next mutation is examined
				}
			}
		}
	}
	return(phe)
}
# Test
par=c(0,.1,0,.1,2,0.9,.1,.1)
opt=c(10,-10)
beta=rbind(c(1,0),c(0,1))
Ne=1e3
T=1e3
x=sim(par,opt,beta,Ne,T)
plot(1:(T+1),x[1,])

# Time taken to adapt
# Parameters: fitness values through time (a vector) and Ne
t.adapt <- function(w,Ne){
	for(t in 2:length(w)){
		s=1-w[t]
		if(s<(1/(2*Ne))){
			ta=t
			break
		}
	}
	return(t)
}

# Simulate multiple traits (no pleiotropy)
# Parameters: Effect size SD (a vector), per-trait mutation rate (assumed to be constaint), optimal phenotype, selection matrix, Ne, time
sim.mv.npl <- function(sig2,U,opt,beta,Ne,T){
	ntrait=length(sig2)
	phe=matrix(0,nrow=ntrait,ncol=(T+1))
	for(t in 2:(T+1)){
		phe[,t]=phe[,(t-1)]
		nmut=rpois(1,lambda=ntrait*U) # Total number of mutations that would occur in this time step
		if(nmut>=1){
			type=sample(1:ntrait,nmut,replace=TRUE) # Type of each mutation
			for(i in 1:nmut){
				phe_mutant=phe[,t]
				effect=rnorm(1,mean=0,sd=sqrt(sig2[type]))
				phe_mutant[type[i]]=phe_mutant[type[i]]+effect # Calculate mutant phenotype
				wa=fitness(phe[,t]-opt,beta);wm=fitness(phe_mutant-opt,beta)
				fp=fix.prob(wa,wm,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					phe[,t]=phe_mutant # If the mutation is fixed, add its effect onto the population mean before the next mutation is examined
				}
			}
		}
	}
	return(phe)
}
# Test
ntrait=10
sig2=exp((1-ntrait):0)
U=1
opt=rep(0,ntrait)
beta=matrix(0,nrow=ntrait,ncol=ntrait);diag(beta)=10
Ne=1e3;T=1e3
x=sim(sig2,U,opt,beta,Ne,T)

# Simulate multiple traits (universal pleiotropy)
# Parameters: Effect size SD (a vector), per-trait mutation rate (assumed to be constaint), optimal phenotype, selection matrix, Ne, time
sim.mv.upl <- function(sig2,U,opt,beta,Ne,T){
	ntrait=length(sig2)
	phe=matrix(0,nrow=ntrait,ncol=(T+1))
	for(t in 2:(T+1)){
		phe[,t]=phe[,(t-1)]
		nmut=rpois(1,lambda=U) # Total number of mutations that would occur in this time step
		if(nmut>=1){
			for(i in 1:nmut){
				phe_mutant=phe[,t]
				effect=rnorm(ntrait,mean=rep(0,ntrait),sd=sqrt(sig2))
				phe_mutant=phe_mutant+effect # Calculate mutant phenotype
				wa=fitness(phe[,t]-opt,beta);wm=fitness(phe_mutant-opt,beta)
				fp=fix.prob(wa,wm,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					phe[,t]=phe_mutant # If the mutation is fixed, add its effect onto the population mean before the next mutation is examined
				}
			}
		}
	}
	return(phe)
}

