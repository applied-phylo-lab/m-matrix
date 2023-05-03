# Simulate evolution with mutational parameters giving rise to zero mutational covariance

setwd("/Users/daohanji/Desktop/M_Matrix/")
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
# Input: mean phenotype through time (2x(T+1) matrix), variance through time (2x(T+1) matrix), optimum (vector of length 2), SDs of fitness functions (vector of length 2), number of replicate lineages
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

# Alter Up, siga, and sigb such that mutational covariance is unchanged
# Adjust deleteriousness of pleiotropic mutations
# Arguments: the original parameter values (par), change to Up
# Components of par: U1, sig1, U2, sig2, Up, r, siga, sigb
par.convert.type1 <- function(par,c){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	par.new=par
	par.new[5]=par.new[5]*c
	par.new[7]=par.new[7]/sqrt(c)
	par.new[8]=par.new[8]/sqrt(c)
	return(par.new)
}

# Alter Up and r such that mutational covariance is unchanged, then alter of U1*sig1^2 and U2*sig2^2 to keep variances unchanged
# Adjust contribution of pleiotropic mutations to each trait's mutational variance
# Arguments: the original parameter values (par), change to r (r becomes c*r), and what other parameter values to change (v, a vector; would be c(1,3) if only mutation rates are being adjusted)
# Components of par: U1, sig1, U2, sig2, Up, r, siga, sigb
par.convert.type2 <- function(par,c,v){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	if(c>(1/abs(par[6]))){
		return("Error")
	}else{
		par.new=par
		par.new[6]=par.new[6]*c # Rescale r;
		par.new[5]=par.new[5]/c # Rescale Up
		c1=(U1*(sig1^2)+(1-(1/c))*Up*(siga^2))/(U1*(sig1^2))
		c2=(U2*(sig2^2)+(1-(1/c))*Up*(sigb^2))/(U2*(sig2^2))
		if(c1>0&c2>0){
			par.new[v[1]]=par.new[v[1]]*c1
			par.new[v[2]]=par.new[v[2]]*c2
			return(par.new)
		}else{
			return("Error")
		}
	}
}

par.default=c(1,0.1,1,0.1,1,0,0.1,0.1) # Starting parameter set
par.all=par.default # Data matrix to store all parameter sets

# Conversion type 1
# All rescaling coefficients to consider
c.all=10^((-5:5)/2)
for(c in 1:length(c.all)){
	par.new=par.convert.type1(par.default,c.all[c])
	if(length(which(par.new==par.default))!=8){
		par.all=rbind(par.all,par.new)
	}
}

# Conversion type 2
if(par.default[6]!=0){
	c.all=((1:9)/10)/par.default[6] # Let r take 0.1, 0.2,..., 0.9 (applicable only if it does not start from 0)
	for(c in 1:length(c.all)){
		par.new=par.convert.type2(par.default,c.all[c],c(1,3))
		if((length(par.new)==8)&(length(which(par.new==par.default))!=8)){
			par.all=rbind(par.all,par.new)
		}
	}
}else{
	c.all=5/(1:10) # Let Up be multiplied by 0.2, ..., 2
	for(c in 1:length(c.all)){
		par.new=par.convert.type2(par.default,c.all[c],c(1,3))
		if((length(par.new)==8)&(length(which(par.new==par.default))!=8)){
			par.all=rbind(par.all,par.new)
		}
	}
}

Ne=1e3
Nrep=200
T=1e4

# Columns: mean of trait 1; mean of trait 2; variance of trait 1; variance of trait 2; covariance; correlation; time taken to adapt (applicable only if there is directional selection)
out.stabilizing=matrix(0,nrow=nrow(par.all),ncol=6) # Scenario 1: Both traits under stabilizing selection
#out.qn=matrix(0,nrow=nrow(par.all),ncol=6) # Scenario 2: One trait under stabilizing selection, the other is neutral
out.dir=matrix(0,nrow=nrow(par.all),ncol=7) # Scenario 3: One trait is under stabilizing selection, the other is under directional selection
#out.dir.biv1=matrix(0,nrow=nrow(par.all),ncol=7) # Scenario 4.1: Both traits are under stabilizing selection (selection aligned w/ mutational covariance)
#out.dir.biv2=matrix(0,nrow=nrow(par.all),ncol=7) # Scenario 4.2: Both traits are under stabilizing selection (selection misaligned w/ mutational covariance)
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	# Scenario 1: Both traits under stabilizing selection
	sim.out=sim(par,c(0,0),c(1,1),Ne,Nrep,T)
	out.stabilizing[c,1]=sim.out[[2]][1,(T+1)]
	out.stabilizing[c,2]=sim.out[[2]][2,(T+1)]
	out.stabilizing[c,3]=sim.out[[3]][1,(T+1)]
	out.stabilizing[c,4]=sim.out[[3]][2,(T+1)]
	out.stabilizing[c,5]=cor(sim.out[[1]][,1],sim.out[[1]][,2])*sqrt(sim.out[[3]][1,(T+1)]*sim.out[[3]][2,(T+1)])
	out.stabilizing[c,6]=cor(sim.out[[1]][,1],sim.out[[1]][,2])

	# Scenario 2
	#sim.out=sim(par,c(0,0),c(1,0),Ne,Nrep,T)
	#out.qn[c,1]=sim.out[[2]][1,(T+1)]
	#out.qn[c,2]=sim.out[[2]][2,(T+1)]
	#out.qn[c,3]=sim.out[[3]][1,(T+1)]
	#out.qn[c,4]=sim.out[[3]][2,(T+1)]
	#out.qn[c,5]=cor(sim.out[[1]][,1],sim.out[[1]][,2])*sqrt(sim.out[[3]][1,(T+1)]*sim.out[[3]][2,(T+1)])
	#out.qn[c,6]=cor(sim.out[[1]][,1],sim.out[[1]][,2])

	# Scenario 3
	sim.out=sim(par,c(0,2),c(1,1),Ne,Nrep,T)
	out.dir[c,1]=sim.out[[2]][1,(T+1)]
	out.dir[c,2]=sim.out[[2]][2,(T+1)]
	out.dir[c,3]=sim.out[[3]][1,(T+1)]
	out.dir[c,4]=sim.out[[3]][2,(T+1)]
	out.dir[c,5]=cor(sim.out[[1]][,1],sim.out[[1]][,2])*sqrt(sim.out[[3]][1,(T+1)]*sim.out[[3]][2,(T+1)])
	out.dir[c,6]=cor(sim.out[[1]][,1],sim.out[[1]][,2])
	out.dir[c,7]=t.adapt(sim.out[[2]],c(0,2),c(1,1),Ne)

	# Scenario 4.1
	#sim.out=sim(par,c(2,2),c(1,1),Ne,Nrep,T)
	#out.dir.biv1[c,1]=sim.out[[2]][1,(T+1)]
	#out.dir.biv1[c,2]=sim.out[[2]][2,(T+1)]
	#out.dir.biv1[c,3]=sim.out[[3]][1,(T+1)]
	#out.dir.biv1[c,4]=sim.out[[3]][2,(T+1)]
	#out.dir.biv1[c,5]=cor(sim.out[[1]][,1],sim.out[[1]][,2])*sqrt(sim.out[[3]][1,(T+1)]*sim.out[[3]][2,(T+1)])
	#out.dir.biv1[c,6]=cor(sim.out[[1]][,1],sim.out[[1]][,2])
	#out.dir.biv1[c,7]=t.adapt(sim.out[[2]],c(2,2),c(1,1),Ne)
	
	# Scenario 4.2
	#sim.out=sim(par,c(2,-2),c(1,1),Ne,Nrep,T)
	#out.dir.biv2[c,1]=sim.out[[2]][1,(T+1)]
	#out.dir.biv2[c,2]=sim.out[[2]][2,(T+1)]
	#out.dir.biv2[c,3]=sim.out[[3]][1,(T+1)]
	#out.dir.biv2[c,4]=sim.out[[3]][2,(T+1)]
	#out.dir.biv2[c,5]=cor(sim.out[[1]][,1],sim.out[[1]][,2])*sqrt(sim.out[[3]][1,(T+1)]*sim.out[[3]][2,(T+1)])
	#out.dir.biv2[c,6]=cor(sim.out[[1]][,1],sim.out[[1]][,2])
	#out.dir.biv2[c,7]=t.adapt(sim.out[[2]],c(2,-2),c(1,1),Ne)
}

write.table(data.frame(par.all,out.stabilizing),file="sim_out_stabilizing_zero_cor.txt",sep="\t")
write.table(data.frame(par.all,out.dir),file="sim_out_directional_zero_cor.txt",sep="\t")

