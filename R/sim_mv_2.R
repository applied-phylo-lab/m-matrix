setwd("/Users/daohanji/Desktop/M_Matrix/revision")
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


# Simulate multiple traits (within-module universal pleiotropy)
# Parameters: Effect size SD (a vector), per-trait mutation rate (assumed to be constaint), optimal phenotype, selection matrix, Ne, time, module information (list of vectors)
sim.mv.pl <- function(sig2,U,opt,beta,Ne,T,module){
	ntrait=length(sig2)
	nmodule=length(module)
	module.size=lengths(module)
	phe=matrix(0,nrow=ntrait,ncol=(T+1))
	for(t in 2:(T+1)){
		phe[,t]=phe[,(t-1)]
		nmut=rpois(1,lambda=nmodule*U) # Total number of mutations that would occur in this time step
		if(nmut>=1){
			type=sample(1:nmodule,nmut,replace=TRUE,prob=module.size/sum(module.size)) # Type of each mutation
			for(i in 1:nmut){
				phe_mutant=phe[,t]
				effect=rnorm(length(module[[type[i]]]),mean=0,sd=sqrt(sig2[module[[type[i]]]]))
				phe_mutant[module[[type[i]]]]=phe_mutant[module[[type[i]]]]+effect # Calculate mutant phenotype
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
# Universal pleiotropy within module, no pleiotropy between modules
#module=list(c(1,3,5,7,9),c(2,4,6,8,10))
module=list(1:5,6:10)
#module=list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10))
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
	phe=sim.mv.pl(sig2,U,opt,beta,Ne,T,module)
	phe.all[row.start:row.end,3:(ncol(phe.all))]=phe
}
write.table(phe.all,file="phe_pl_2m.txt",sep="\t")
# Extract end-point phenotype and rearrange
phe.end=matrix(0,nrow=Nrep,ncol=ntrait)
for(n in 1:Nrep){
	sub=phe.all[which(phe.all[,1]==n),]
	phe.end[n,]=sub[,ncol(sub)]
}
write.table(phe.end,file="phe_end_pl_2m.txt",sep="\t")
rmat=cov(phe.end)
write.table(rmat,file="rmat_pl_2m.txt",sep="\t")

# Multivariate
ntrait=10 # Number of traits
sig2=exp((1-ntrait):0) # Mutational variances (effect size variances)
mmat=matrix(0,nrow=ntrait,ncol=ntrait);diag(mmat)=sig2

# 2 modules
setwd("./out_2m")
rmat<-read.table("rmat_pl_2m.txt",sep="\t") # Read R matrix
rmat=data.matrix(rmat) # convert to matrix format
cor(log(diag(mmat)),log(diag(rmat)))
lm(log(diag(rmat))~log(diag(mmat)))
# Plot
d=data.frame(c(rep("A",nrow(d)/2),rep("B",nrow(d)/2)),log(diag(mmat)),log(diag(rmat)));colnames(d)=c("mod","vm","vr")
g<-ggplot(d,aes(x=vm,y=vr,colour=mod))
g=g+geom_point(show.legend=FALSE)
g=g+scale_color_manual(values=c("A"="orange","B"="purple"))
g=g+geom_smooth(method="lm",colour="black")
g=g+geom_smooth(method="lm",data=subset(d,mod=="A"),colour="orange",fill="yellow")
g=g+geom_smooth(method="lm",data=subset(d,mod=="B"),colour="purple",fill="blue")
g=g+theme_classic()
g=g+xlab(expression(paste("ln(M variance)")))+ylab(expression(paste("ln(R variance)")))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_pl_2m.pdf",plot=g,width=5.5,height=5)
cor.test(d[which(d$mod=="A"),1],d[which(d$mod=="A"),2]);lm(d[which(d$mod=="A"),2]~d[which(d$mod=="A"),1])
cor.test(d[which(d$mod=="B"),1],d[which(d$mod=="B"),2]);lm(d[which(d$mod=="B"),2]~d[which(d$mod=="B"),1])
setwd("..")

setwd("./out_2m_alt")
rmat<-read.table("rmat_pl_2m.txt",sep="\t") # Read R matrix
rmat=data.matrix(rmat) # convert to matrix format
cor(log(diag(mmat)),log(diag(rmat)))
lm(log(diag(rmat))~log(diag(mmat)))
# Plot
mod=c("A","B","A","B","A","B","A","B","A","B")
d=data.frame(mod,log(diag(mmat)),log(diag(rmat)));colnames(d)=c("mod","vm","vr")
g<-ggplot(d,aes(x=vm,y=vr,colour=mod))
g=g+geom_point(show.legend=FALSE)
g=g+scale_color_manual(values=c("A"="orange","B"="purple"))
g=g+geom_smooth(method="lm",colour="black")
g=g+geom_smooth(method="lm",data=subset(d,mod=="A"),colour="orange",fill="yellow")
g=g+geom_smooth(method="lm",data=subset(d,mod=="B"),colour="purple",fill="blue")
g=g+theme_classic()
g=g+xlab(expression(paste("ln(M variance)")))+ylab(expression(paste("ln(R variance)")))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_pl_2m.pdf",plot=g,width=5.5,height=5)
cor.test(d[which(d$mod=="A"),2],d[which(d$mod=="A"),3]);lm(d[which(d$mod=="A"),3]~d[which(d$mod=="A"),2])
cor.test(d[which(d$mod=="B"),2],d[which(d$mod=="B"),3]);lm(d[which(d$mod=="B"),3]~d[which(d$mod=="B"),2])
setwd("..")

# 5 modules
setwd("./out_5m")
rmat<-read.table("rmat_pl_5m.txt",sep="\t") # Read R matrix
rmat=data.matrix(rmat) # convert to matrix format
cor(log(diag(mmat)),log(diag(rmat)))
lm(log(diag(rmat))~log(diag(mmat)))
# Plot
d=data.frame(log(diag(mmat)),log(diag(rmat)));colnames(d)=c("vm","vr")
g<-ggplot(d,aes(x=vm,y=vr))
g=g+geom_point()+geom_smooth(method="lm")
g=g+theme_classic()
g=g+xlab(expression(paste("ln(M variance)")))+ylab(expression(paste("ln(R variance)")))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_pl_5m.pdf",plot=g,width=5.5,height=5)


