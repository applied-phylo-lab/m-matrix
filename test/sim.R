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

# Simulate a single lineage and write phenotype thru time
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
				fp=fix.prob(phe[,t]-opt,phe_mutant-opt,a,Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					phe[,t]=phe_mutant # If the mutation is fixed, add its effect onto the population mean before the next mutation is examined
				}
			}
		}
	}
	return(phe)
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

par.default=c(1,0.1,1,0.1,1,0.5,0.1,0.1) # Starting parameter set
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

rownames(par.all)=1:nrow(par.all)
write.table(par.all,file="par_all.txt",sep="\t")

Ne=1e3

# Univariate directional selection
Nrep=100
T=5e3
opt=c(2,0)
width=c(1,1)
out.all=matrix(0,nrow=Nrep*nrow(par.all)*2,ncol=(T+1+10))
no=rep(0,nrow(out.all))
colnames(out.all)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","rep","trait",0:T)
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row.start=row=(c-1)*2*Nrep+1
	row.end=row=c*2*Nrep
	no[row.start:row.end]=c
	for(i in row.start:row.end){
		out.all[i,1:8]=par
	}
	for(n in 1:Nrep){
		row=(c-1)*2*Nrep+2*n-1
		out.all[row:(row+1),9]=n
		out.all[row:(row+1),10]=c(1,2)
		phe=sim(par,opt,width,Ne,T)
		out.all[row:(row+1),11:(ncol(out.all))]=phe
	}
}
out.all=data.frame(no,out.all)
write.table(out.all,file="out_dir.txt",sep="\t")

# Mulivariate directional selection
# Selection gradient misaligned with mutational covariance
Nrep=100
T=5e3
opt=c(2,-2)
width=c(1,1)
out.all=matrix(0,nrow=Nrep*nrow(par.all)*2,ncol=(T+1+10))
no=rep(0,nrow(out.all))
colnames(out.all)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","rep","trait",0:T)
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row.start=row=(c-1)*2*Nrep+1
	row.end=row=c*2*Nrep
	no[row.start:row.end]=c
	for(i in row.start:row.end){
		out.all[i,1:8]=par
	}
	for(n in 1:Nrep){
		row=(c-1)*2*Nrep+2*n-1
		out.all[row:(row+1),9]=n
		out.all[row:(row+1),10]=c(1,2)
		phe=sim(par,opt,width,Ne,T)
		out.all[row:(row+1),11:(ncol(out.all))]=phe
	}
}
out.all=data.frame(no,out.all)
write.table(out.all,file="out_dir2.txt",sep="\t")

# Stabilizing selection
Nrep=100
T=5e4
opt=c(0,0)
width=c(5,5)
out.all=matrix(0,nrow=Nrep*nrow(par.all)*2,ncol=(T+1+10))
no=rep(0,nrow(out.all))
colnames(out.all)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","rep","trait",0:T)
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row.start=row=(c-1)*2*Nrep+1
	row.end=row=c*2*Nrep
	no[row.start:row.end]=c
	for(i in row.start:row.end){
		out.all[i,1:8]=par
	}
	for(n in 1:Nrep){
		row=(c-1)*2*Nrep+2*n-1
		out.all[row:(row+1),9]=n
		out.all[row:(row+1),10]=c(1,2)
		phe=sim(par,opt,width,Ne,T)
		out.all[row:(row+1),11:(ncol(out.all))]=phe
	}
}
out.all=data.frame(no,out.all)
out.end=data.frame(out.all[,1:11],out.all[,ncol(out.all)])
write.table(out.all,file="out_stab.txt",sep="\t")
write.table(out.end,file="out_stab_end.txt",sep="\t")

