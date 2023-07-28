setwd("/Users/daohanji/Desktop/M_Matrix/")
library(mvtnorm)

# Simulate more than 2 traits
# Only consider two extreme cases: no pleiotropy vs universal pleiotropy
# Mutation variance proportion to effect size variance

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

# Function for simulations
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

ntrait=10 # Number of traits
sig2=exp((1-ntrait):0) # Mutational variances (effect size variances)
U=1 # Per-trait mutation rate (set constant for all traits)
Ne=1e3

# Stabilizing selection
# No pleiotropy
opt=rep(0,ntrait) # Optimal phenotype
beta=matrix(0,nrow=ntrait,ncol=ntrait);diag(beta)=10 # Selection matrix
T=5e4
Nrep=500
phe.all=matrix(0,nrow=ntrait*Nrep,ncol=(T+1+2))
for(n in 1:Nrep){
	row.start=(n-1)*ntrait+1
	row.end=n*ntrait
	phe.all[row.start:row.end,1]=n
	phe.all[row.start:row.end,2]=1:ntrait
	phe=sim.mv.npl(sig2,U,opt,beta,Ne,T)
	phe.all[row.start:row.end,3:(ncol(phe.all))]=phe
}
write.table(phe.all,file="phe_npl.txt",sep="\t")
# Extract end-point phenotype and rearrange
phe.end=matrix(0,nrow=Nrep,ncol=ntrait)
for(n in 1:Nrep){
	sub=phe.all[which(phe.all[,1]==n),]
	phe.end[n,]=sub[,ncol(sub)]
}
write.table(phe.end,file="phe_end_npl.txt",sep="\t")
rmat=cov(phe.end)
write.table(rmat,file="rmat_npl.txt",sep="\t")

# Stabilizing selection
# Universal pleiotropy
opt=rep(0,ntrait) # Optimal phenotype
beta=matrix(0,nrow=ntrait,ncol=ntrait);diag(beta)=10 # Selection matrix
T=5e4
Nrep=500
phe.all=matrix(0,nrow=ntrait*Nrep,ncol=(T+1+2))
for(n in 1:Nrep){
	row.start=(n-1)*ntrait+1
	row.end=n*ntrait
	phe.all[row.start:row.end,1]=n
	phe.all[row.start:row.end,2]=1:ntrait
	phe=sim.mv.upl(sig2,U,opt,beta,Ne,T)
	phe.all[row.start:row.end,3:(ncol(phe.all))]=phe
}
write.table(phe.all,file="phe_upl.txt",sep="\t")
# Extract end-point phenotype and rearrange
phe.end=matrix(0,nrow=Nrep,ncol=ntrait)
for(n in 1:Nrep){
	sub=phe.all[which(phe.all[,1]==n),]
	phe.end[n,]=sub[,ncol(sub)]
}
write.table(phe.end,file="phe_end_upl.txt",sep="\t")
rmat=cov(phe.end)
write.table(rmat,file="rmat_upl.txt",sep="\t")


