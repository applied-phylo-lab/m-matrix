# Function(s) to simulate evolution

library(mvtnorm)

# Fitness (multivariate Gaussian fitness function)
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

# Fixation probability
# Calculate fixation probability of a mutation given ancestral phenotype (x1), mutant phenotype (x2), SD of fitness function (a), and effective population size (Ne)
# x1, x2, and a are all numbers when a single trait is considered; they are vectors of the same length if multiple traits are considered
fix.prob <- function(x1,x2,a,Ne){
	wa=fitness(x1,a) # Ancestral fitness
	wm=fitness(x2,a) # Mutant fitness
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

# Simulate evolution
# Parameters: mutational parameters (vector of length 8), optimal phenotype (vector of length 2), SDs of fitness functions (vector of length 2), Ne, number of replicate lineages, time of simulation
# a=c(0,0) for neutral evolution, opt!=c(0,0) for directional selection
sim <- function(par,opt,a,Ne,Nrep,T){
	# Read and interpret the mutational parameters
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	sig.np=c(sig1,sig2)
	mp=rbind(c(siga^2,r*siga*sigb),c(r*siga*sigb,sigb^2))

	phe.end=matrix(0,nrow=Nrep,ncol=2) # End-point phenotypes of all replicate lineages
	phe.all=list() # All phenotypes through time
	mean.all=matrix(0,nrow=2,ncol=(T+1)) # Mean phenotype through time
	var.all=matrix(0,nrow=2,ncol=(T+1)) # Variance through time
	for(n in 1:Nrep){
		phe=matrix(0,nrow=2,ncol=(T+1)) # Phenotypes through time for this replicate lineage
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
					fp=fix.prob(phe[,t]-opt,phe_mutant-opt,a,Ne)
					if.fix=rbinom(n=1,size=1,prob=fp)
					if(if.fix==1){
						phe[,t]=phe_mutant # If the mutation is fixed, add its effect onto the population mean before the next mutation is examined
					}
				}
			}
		}
		phe.end[n,]=phe[,(T+1)]
		phe.all[[n]]=phe
	}
	# Calculate mean and variance through time
	for(t in 2:(T+1)){
		d=matrix(0,nrow=Nrep,ncol=2)
		for(n in 1:Nrep){
			d[n,]=phe.all[[n]][,t]
		}
		mean.all[1,t]=mean(d[,1]);mean.all[2,t]=mean(d[,2])
		var.all[1,t]=var(d[,1]);var.all[2,t]=var(d[,2])
	}
	# Put data matrices into a list
	out=list(phe.end,mean.all,var.all)
	names(out)=c("phe_end","mean.all","var.all") # Rename elements in the list
	return(out)
}

# Calculate time taken before the optimum is reached
# Input: mean phenotype through time (2x(T+1) matrix), optimum (vector of length 2), SDs of fitness functions (vector of length 2), Ne
t.adapt <- function(dm,opt,a,Ne){
	tmin=0
	if(a[1]!=0&a[2]!=0){
		for(t in 2:(T+1)){
			w=fitness(dm[,t]-opt,a)
			s=1-w
			if(s<(1/(2*Ne))){
				tmin[1]=t
				break
			}
		}
	}
	return(tmin)
}

# Test run
# Use parameter sets that give clear predictions
# Parameters used for all tests
Ne=1e2
Nrep=100
T=1000

# Parameter set 1: no pleiotropy (Up=0), stabilizing selection
par=c(1,0.1,1,0.1,0,0,0.1,0.1)
opt=c(0,0)
a=c(1,1)
out1=sim(par,opt,a,Ne,Nrep,T)
var(out1[[1]][,1]) # Check variance
var(out1[[1]][,1])/((par[1]*(par[2]^2)+par[5]*(par[7]^2))*T) # Check ratio of variance and expected variance

# Parameter set 2: strong pleiotropy (U1=U2=0, Up>0, r=0), stabilizing selection
par=c(0,0.1,0,0.1,1,0,0.1,0.1)
opt=c(0,0)
a=c(1,1)
out2=sim(par,opt,a,Ne,Nrep,T)
var(out2[[1]][,1]) # Check variance
var(out2[[1]][,1])/((par[1]*(par[2]^2)+par[5]*(par[7]^2))*T) # Check ratio of variance and expected variance

# Parameter set 3: no pleiotropy (Up=0), trait 2 under stabilizing selection
par=c(1,0.1,1,0.1,0,0,0.1,0.1)
opt=c(0,0)
a=c(0,0.1) # Set strong stabilizing selection on trait 2 to make effect of constraint easy to see
out3=sim(par,opt,a,Ne,Nrep,T)
var(out3[[1]][,1]) # Check variance
var(out3[[1]][,1])/((par[1]*(par[2]^2)+par[5]*(par[7]^2))*T) # Check ratio of variance and expected variance

# Parameter set 4: strong pleiotropy (U1=U2=0, Up>0, r=0), trait 2 under stabilizing selection
par=c(0,0.1,0,0.1,1,0,0.1,0.1)
opt=c(0,0)
a=c(0,0.1) # Set strong stabilizing selection on trait 2 to make effect of constraint easy to see
out4=sim(par,opt,a,Ne,Nrep,T)
var(out4[[1]][,1]) # Check variance
var(out4[[1]][,1])/((par[1]*(par[2]^2)+par[5]*(par[7]^2))*T) # Check ratio of variance and expected variance

# Parameter set 5: no pleiotropy, univariate directional selection
par=c(1,0.1,1,0.1,0,0,0.1,0.1)
opt=c(2,0) # Trait 1 under directional selection
a=c(0.1,0.01) # Set strong stabilizing selection on trait 2 to make effect of constraint easy to see
Ne=1e2
Nrep=100
T=1000
out5=sim(par,opt,a,Ne,Nrep,T)
mean(out5[[1]][,1]) # Check end-point mean phenotype
t.adapt(out5[[2]],opt,a,Ne) # Check time taken to adapt

# Parameter set 6: strong pleiotropy, univariate directional selection
par=c(0,0.1,0,0.1,1,0,0.1,0.1)
opt=c(2,0) # Trait 1 under directional selection
a=c(0.1,0.01) # Set strong stabilizing selection on trait 2 to make effect of constraint easy to see
out6=sim(par,opt,a,Ne,Nrep,T)
mean(out6[[1]][,1]) # Check end-point mean phenotype
t.adapt(out6[[2]],opt,a,Ne) # Check time taken to adapt

