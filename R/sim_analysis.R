# Analyze output of simulations

library(ggplot2)

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

# Effective population size used throughout the study
Ne=1e3

# Zero mutational covariance
setwd("./zero_cor")
par.all<-read.table("par_all_zero_cor.txt",sep="\t")
frac=(par.all[,5]*(par.all[,7]^2))/(par.all[,1]*(par.all[,2]^2)+par.all[,5]*(par.all[,7]^2)) # Fraction of mutational variance contributed by pleiotropic mutations

# Directional selection
# Parameters should be the same as those used for the simulations (see sim_zero_cor.R)
opt=c(10,0)
beta=rbind(c(1,0),c(0,1))
flist=c("out_dir_zero_cor.txt","out_dir_zero_cor_set2.txt")
for(i in 1:length(flist)){
	di<-read.table(flist[i],sep="\t")
	if(i==1){
		d=di
	}else{
		d=rbind(d,di)
	}
}
#d<-read.table("out_dir_zero_cor.txt",sep="\t")
# Calculate fitness thru time for each lineage, and calculate mean and variance of fitness at each time point
w.all=matrix(0,nrow=nrow(d)/2,ncol=ncol(d)-11) # Data matrix for all fitness values
out=matrix(0,nrow=nrow(par.all),ncol=2) # Data matrix for mean of integrated fitness (i.e., mean over time) and its variance
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row=which(d[,1]==c) # Extract rows corresponding to the mutational parameter combination under consideration
	dsub=d[row,]
	dsub1=dsub[which(dsub$trait==1),];dsub1=dsub1[,12:ncol(dsub1)] # Extract values of trait 1 thru time
	dsub2=dsub[which(dsub$trait==2),];dsub2=dsub2[,12:ncol(dsub2)] # Extract values of trait 2 thru time

	w=matrix(0,nrow=nrow(dsub1),ncol=ncol(dsub1)) # Data matrix to store all fitness values
	for(n in 1:nrow(dsub1)){
		# Prepare to write the fitness values
		row.write=(c-1)*nrow(dsub1)+n
		for(t in 1:ncol(dsub1)){
			phe=c(dsub1[n,t],dsub2[n,t])
			w[n,t]=fitness(phe-opt,beta)
			w.all[row.write,t]=w[n,t] # Write into the matrix that stores all fitness values
		}
	}
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
write.table(w.all,file="fitness_dir_zero_cor.txt",sep="\t")
# Calculate integrated fitness from the fitness data file
Nrep=500
w.all<-read.table("fitness_dir_zero_cor.txt",sep="\t")
out=matrix(0,nrow=nrow(par.all),ncol=2)
for(c in 1:nrow(par.all)){
	row.start=(c-1)*Nrep+1
	row.end=c*Nrep
	w=w.all[row.start:row.end,]
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
out=data.frame(par.all,out[,1],sqrt(out[,2])) # Generate the data frame for plot
colnames(out)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","wim","wisd") # Rename columns
# Make plots
# Series 1 (effect size distribution of pleiotropic mutations vary)
sub1=out[1:9,];sub2=rbind(out[1,],out[10:19,])
#plot(log10(sub1[,5]),sub1[,9]) # Quick check using basic plot functions
g<-ggplot(sub1,aes(x=log10(Up),y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir_series1_zero_cor.pdf",plot=g,width=5.5,height=5)
# Series 2 (contribution of pleiotropic mutations to mutational variances vary)
g<-ggplot(sub2,aes(x=Up,y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab(expression(U[P]))+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))

ggsave("plot_dir_series2_zero_cor.pdf",plot=g,width=5.5,height=5)

# Stabilizing selection
# End-point variance
flist=c("out_stab_end_zero_cor.txt","out_stab_end_zero_cor_set2.txt") # Combine a list of files
for(i in 1:length(flist)){
	di<-read.table(flist[i],sep="\t")
	if(i==1){
		d=di
	}else{
		d=rbind(d,di)
	}
}
#d<-read.table("out_stab_end_zero_cor.txt",sep="\t")
out=matrix(0,nrow=nrow(par.all),ncol=4)
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row=which(d[,1]==c)
	dsub=d[row,]
	dsub1=dsub[which(dsub$trait==1),]
	dsub2=dsub[which(dsub$trait==2),]
	out[c,1]=var(dsub1[,12])
	out[c,2]=var(dsub2[,12])
	out[c,3]=cor(dsub1[,12],dsub2[,12])
	out[c,4]=cor.test(dsub1[,12],dsub2[,12])$p.value
}
# Calculate SDs of the above estimates
error=matrix(0,nrow=nrow(out),ncol=3)
ss=nrow(d)/(2*nrow(par.all)) # Sample size
error[,1]=sqrt(2*(out[,1]^2)/(ss-1))
error[,2]=sqrt(2*(out[,2]^2)/(ss-1))
error[,3]=(1-out[,3]^2)/(ss-1)
out=data.frame(par.all,out,error)
colnames(out)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","var1","var2","rcor","rcor.pv","var1.sd","var2.sd","rcor.sd")
sub1=out[1:9,];sub2=rbind(out[1,],out[10:19,])
g<-ggplot()+geom_line(data=sub1,aes(x=log10(Up),y=var1),colour="purple")+geom_line(data=sub1,aes(x=log10(Up),y=var2),colour="orange")
g=g+geom_pointrange(data=sub1,aes(x=log10(Up),y=var1,ymin=var1-var1.sd,ymax=var1+var1.sd),colour="purple")+geom_pointrange(data=sub1,aes(x=log10(Up),y=var2,ymin=var2-var2.sd,ymax=var2+var2.sd),colour="orange")
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Variance")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
g=g+ylim(c(0,0.01))
ggsave("plot_stab_series1_zero_cor.pdf",plot=g,width=5.5,height=5)
g<-ggplot()+geom_line(data=sub2,aes(x=Up,y=var1),colour="purple")+geom_line(data=sub2,aes(x=Up,y=var2),colour="orange")
g=g+geom_pointrange(data=sub2,aes(x=Up,y=var1,ymin=var1-var1.sd,ymax=var1+var1.sd),colour="purple")+geom_pointrange(data=sub2,aes(x=Up,y=var2,ymin=var2-var2.sd,ymax=var2+var2.sd),colour="orange")
g=g+theme_classic()
g=g+xlab(expression(U[P]))+ylab("Variance")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
g=g+ylim(c(0,0.01))
ggsave("plot_stab_series2_zero_cor.pdf",plot=g,width=5.5,height=5)

# Variance thru time (check whether the curve is smooth)
d<-read.table("out_stab_zero_cor.txt",sep="\t")
# Pcik one parameter combination as example
c=10
par=par.all[c,]
row=which(d[,1]==c)
dsub=d[row,]
dsub1=dsub[which(dsub[,11]==1),];dsub1=dsub1[,12:ncol(dsub1)]
dsub2=dsub[which(dsub[,11]==2),];dsub2=dsub2[,12:ncol(dsub2)]
v=rep(0,ncol(dsub1))
for(t in 1:ncol(dsub1)){
	v[t]=var(dsub1[,t])
}
plot(1:length(v),v)

setwd("..")

# Non-zero mutational covariance
setwd("./non_zero_cor")
par.all<-read.table("par_all.txt",sep="\t")
frac=(par.all[,5]*(par.all[,7]^2))/(par.all[,1]*(par.all[,2]^2)+par.all[,5]*(par.all[,7]^2)) # Fraction of mutational variance contributed by pleiotropic mutations

# Directional selection on one trait
# Parameters should be the same as those used for the simulations (see sim.R)
opt=c(10,0)
beta=rbind(c(1,0),c(0,1))
flist=c("out_dir.txt","out_dir_set2.txt") # Combine a list of files
for(i in 1:length(flist)){
	di<-read.table(flist[i],sep="\t")
	if(i==1){
		d=di
	}else{
		d=rbind(d,di)
	}
}
#d<-read.table("out_dir.txt",sep="\t")
# Calculate fitness thru time for each lineage, and calculate mean and variance of fitness at each time point
w.all=matrix(0,nrow=nrow(d)/2,ncol=ncol(d)-11) # Data matrix for all fitness values
out=matrix(0,nrow=nrow(par.all),ncol=2) # Data matrix for mean of integrated fitness (i.e., mean over time) and its variance
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row=which(d[,1]==c) # Extract rows corresponding to the mutational parameter combination under consideration
	dsub=d[row,]
	dsub1=dsub[which(dsub$trait==1),];dsub1=dsub1[,12:ncol(dsub1)] # Extract values of trait 1 thru time
	dsub2=dsub[which(dsub$trait==2),];dsub2=dsub2[,12:ncol(dsub2)] # Extract values of trait 2 thru time

	w=matrix(0,nrow=nrow(dsub1),ncol=ncol(dsub1)) # Data matrix to store all fitness values
	for(n in 1:nrow(dsub1)){
		# Prepare to write the fitness values
		row.write=(c-1)*nrow(dsub1)+n
		for(t in 1:ncol(dsub1)){
			phe=c(dsub1[n,t],dsub2[n,t])
			w[n,t]=fitness(phe-opt,beta)
			w.all[row.write,t]=w[n,t] # Write into the matrix that stores all fitness values
		}
	}
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
write.table(w.all,file="fitness_dir.txt",sep="\t")
Nrep=500
w.all<-read.table("fitness_dir.txt",sep="\t")
out=matrix(0,nrow=nrow(par.all),ncol=2)
for(c in 1:nrow(par.all)){
	row.start=(c-1)*Nrep+1
	row.end=c*Nrep
	w=w.all[row.start:row.end,]
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
out=data.frame(par.all,out[,1],sqrt(out[,2])) # Generate the data frame for plot
colnames(out)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","wim","wisd") # Rename columns
# Make plots
# Series 1 (effect size distribution of pleiotropic mutations vary)
sub1=out[1:9,];sub2=rbind(out[1,],out[10:15,])
#plot(log10(sub1[,5]),sub1[,9]) # Quick check using basic plot functions
g<-ggplot(sub1,aes(x=log10(Up),y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir_series1.pdf",plot=g,width=5.5,height=5)
# Series 2 (contribution of pleiotropic mutations to mutational variances vary)
g<-ggplot(sub2,aes(x=r,y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab("Correlation coefficient")+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir_series2.pdf",plot=g,width=5.5,height=5)


# Directional selection on two traits, w/o correlational selection
# Parameters should be the same as those used for the simulations (see sim.R)
opt=c(10,-10)
beta=rbind(c(1,0),c(0,1))
flist=c("out_dir1.txt","out_dir1_set2.txt") # Combine a list of files
for(i in 1:length(flist)){
	di<-read.table(flist[i],sep="\t")
	if(i==1){
		d=di
	}else{
		d=rbind(d,di)
	}
}
#d<-read.table("out_dir1.txt",sep="\t")
# Calculate fitness thru time for each lineage, and calculate mean and variance of fitness at each time point
w.all=matrix(0,nrow=nrow(d)/2,ncol=ncol(d)-11) # Data matrix for all fitness values
out=matrix(0,nrow=nrow(par.all),ncol=2) # Data matrix for mean of integrated fitness (i.e., mean over time) and its variance
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row=which(d[,1]==c) # Extract rows corresponding to the mutational parameter combination under consideration
	dsub=d[row,]
	dsub1=dsub[which(dsub$trait==1),];dsub1=dsub1[,12:ncol(dsub1)] # Extract values of trait 1 thru time
	dsub2=dsub[which(dsub$trait==2),];dsub2=dsub2[,12:ncol(dsub2)] # Extract values of trait 2 thru time

	w=matrix(0,nrow=nrow(dsub1),ncol=ncol(dsub1)) # Data matrix to store all fitness values
	for(n in 1:nrow(dsub1)){
		# Prepare to write the fitness values
		row.write=(c-1)*nrow(dsub1)+n
		for(t in 1:ncol(dsub1)){
			phe=c(dsub1[n,t],dsub2[n,t])
			w[n,t]=fitness(phe-opt,beta)
			w.all[row.write,t]=w[n,t] # Write into the matrix that stores all fitness values
		}
	}
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
write.table(w.all,file="fitness_dir1.txt",sep="\t")
Nrep=500
w.all<-read.table("fitness_dir1.txt",sep="\t")
out=matrix(0,nrow=nrow(par.all),ncol=2)
for(c in 1:nrow(par.all)){
	row.start=(c-1)*Nrep+1
	row.end=c*Nrep
	w=w.all[row.start:row.end,]
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
out=data.frame(par.all,out[,1],sqrt(out[,2])) # Generate the data frame for plot
colnames(out)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","wim","wisd") # Rename columns
# Make plots
# Series 1 (effect size distribution of pleiotropic mutations vary)
sub1=out[1:9,];sub2=rbind(out[1,],out[10:15,])
#plot(log10(sub1[,5]),sub1[,9]) # Quick check using basic plot functions
g<-ggplot(sub1,aes(x=log10(Up),y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir1_series1.pdf",plot=g,width=5.5,height=5)
# Series 2 (contribution of pleiotropic mutations to mutational variances vary)
g<-ggplot(sub2,aes(x=r,y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab("Correlation coefficient")+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir1_series2.pdf",plot=g,width=5.5,height=5)

# Directional selection on two traits, w/ correlational selection
# Parameters should be the same as those used for the simulations (see sim.R)
opt=c(10,-10)
beta=rbind(c(1,-0.9),c(-0.9,1))
flist=c("out_dir2.txt","out_dir2_set2.txt") # Combine a list of files
for(i in 1:length(flist)){
	di<-read.table(flist[i],sep="\t")
	if(i==1){
		d=di
	}else{
		d=rbind(d,di)
	}
}
#d<-read.table("out_dir2.txt",sep="\t")
# Calculate fitness thru time for each lineage, and calculate mean and variance of fitness at each time point
w.all=matrix(0,nrow=nrow(d)/2,ncol=ncol(d)-11) # Data matrix for all fitness values
out=matrix(0,nrow=nrow(par.all),ncol=2) # Data matrix for mean of integrated fitness (i.e., mean over time) and its variance
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row=which(d[,1]==c) # Extract rows corresponding to the mutational parameter combination under consideration
	dsub=d[row,]
	dsub1=dsub[which(dsub$trait==1),];dsub1=dsub1[,12:ncol(dsub1)] # Extract values of trait 1 thru time
	dsub2=dsub[which(dsub$trait==2),];dsub2=dsub2[,12:ncol(dsub2)] # Extract values of trait 2 thru time

	w=matrix(0,nrow=nrow(dsub1),ncol=ncol(dsub1)) # Data matrix to store all fitness values
	for(n in 1:nrow(dsub1)){
		# Prepare to write the fitness values
		row.write=(c-1)*nrow(dsub1)+n
		for(t in 1:ncol(dsub1)){
			phe=c(dsub1[n,t],dsub2[n,t])
			w[n,t]=fitness(phe-opt,beta)
			w.all[row.write,t]=w[n,t] # Write into the matrix that stores all fitness values
		}
	}
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
write.table(w.all,file="fitness_dir2.txt",sep="\t")
Nrep=500
w.all<-read.table("fitness_dir2.txt",sep="\t")
out=matrix(0,nrow=nrow(par.all),ncol=2)
for(c in 1:nrow(par.all)){
	row.start=(c-1)*Nrep+1
	row.end=c*Nrep
	w=w.all[row.start:row.end,]
	w.int=rowMeans(w)
	out[c,1]=mean(w.int)
	out[c,2]=var(w.int)
}
out=data.frame(par.all,out[,1],sqrt(out[,2])) # Generate the data frame for plot
colnames(out)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","wim","wisd") # Rename columns
# Make plots
# Series 1 (effect size distribution of pleiotropic mutations vary)
sub1=out[1:9,];sub2=rbind(out[1,],out[10:15,])
#plot(log10(sub1[,5]),sub1[,9]) # Quick check using basic plot functions
g<-ggplot(sub1,aes(x=log10(Up),y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir2_series1.pdf",plot=g,width=5.5,height=5)
# Series 2 (contribution of pleiotropic mutations to mutational variances vary)
g<-ggplot(sub2,aes(x=r,y=wim))
g=g+geom_line()+geom_pointrange(aes(ymin=wim-wisd,ymax=wim+wisd))
g=g+theme_classic()
g=g+xlab("Correlation coefficient")+ylab("Mean integrated fitness")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
ggsave("plot_dir2_series2.pdf",plot=g,width=5.5,height=5)

# Stabilizing selection
# End-point variance
flist=c("out_stab1_end.txt","out_stab1_end_set2.txt","out_stab1_end_set3.txt") # Combine a list of files
for(i in 1:length(flist)){
	di<-read.table(flist[i],sep="\t")
	if(i==1){
		d=di
	}else{
		d=rbind(d,di)
	}
}
#d<-read.table("out_stab_end_zero_cor.txt",sep="\t")
out=matrix(0,nrow=nrow(par.all),ncol=4)
for(c in 1:nrow(par.all)){
	par=par.all[c,]
	row=which(d[,1]==c)
	dsub=d[row,]
	dsub1=dsub[which(dsub$trait==1),]
	dsub2=dsub[which(dsub$trait==2),]
	out[c,1]=var(dsub1[,12])
	out[c,2]=var(dsub2[,12])
	out[c,3]=cor(dsub1[,12],dsub2[,12])
	out[c,4]=cor.test(dsub1[,12],dsub2[,12])$p.value
}
# Calculate SDs of the above estimates
error=matrix(0,nrow=nrow(out),ncol=3)
ss=nrow(d)/(2*nrow(par.all)) # Sample size
error[,1]=sqrt(2*(out[,1]^2)/(ss-1))
error[,2]=sqrt(2*(out[,2]^2)/(ss-1))
error[,3]=(1-out[,3]^2)/(ss-1)
out=data.frame(par.all,out,error)
colnames(out)=c("U1","sig1","U2","sig2","Up","r","siga","sigb","var1","var2","rcor","rcor.pv","var1.sd","var2.sd","rcor.sd")
sub1=out[1:9,];sub2=rbind(out[1,],out[10:15,])
g<-ggplot()+geom_line(data=sub1,aes(x=log10(Up),y=var1),colour="purple")+geom_line(data=sub1,aes(x=log10(Up),y=var2),colour="orange")
g=g+geom_pointrange(data=sub1,aes(x=log10(Up),y=var1,ymin=var1-var1.sd,ymax=var1+var1.sd),colour="purple")+geom_pointrange(data=sub1,aes(x=log10(Up),y=var2,ymin=var2-var2.sd,ymax=var2+var2.sd),colour="orange")
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Variance")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
g=g+ylim(c(0,0.01))
ggsave("plot_stab_series1.pdf",plot=g,width=5.5,height=5)
g<-ggplot()+geom_line(data=sub2,aes(x=r,y=var1),colour="purple")+geom_line(data=sub2,aes(x=r,y=var2),colour="orange")
g=g+geom_pointrange(data=sub2,aes(x=r,y=var1,ymin=var1-var1.sd,ymax=var1+var1.sd),colour="purple")+geom_pointrange(data=sub2,aes(x=r,y=var2,ymin=var2-var2.sd,ymax=var2+var2.sd),colour="orange")
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],U[P])))+ylab("Variance")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
g=g+ylim(c(0,0.01))
ggsave("plot_stab_series2.pdf",plot=g,width=5.5,height=5)

setwd("..")

# Multivariate
setwd("./mv")
Nrep=100
ntrait=10 # Number of traits
sig2=exp((1-ntrait):0) # Mutational variances (effect size variances)
mmat=matrix(0,nrow=ntrait,ncol=ntrait);diag(mmat)=sig2

# No pleiotropy
rmat<-read.table("rmat_npl.txt",sep="\t") # Read R matrix
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
ggsave("plot_npl.pdf",plot=g,width=5.5,height=5)

# Universal pleiotropy
rmat<-read.table("rmat_upl.txt",sep="\t") # Read R matrix
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
ggsave("plot_upl.pdf",plot=g,width=5.5,height=5)

# Calculate variance through time and phylogenetic signal
d<-read.table("phe_npl.txt",sep="\t")
var.all=matrix(0,nrow=ntrait,ncol=ncol(d))
for(i in 1:ntrait){
	dsub=d[which(d[,2]==i),]
	for(t in 3:ncol(dsub)){
		var.all[i,(t-2)]=var(dsub[,t])
	}
}
write.table(var.all,file="var_all_npl.txt",sep="\t")
cor.all=rep(0,ntrait)
for(i in 1:ntrait){
	cor.all[i]=cor(1:ncol(var.all),var.all[i,])
}
write.table(cor.all,file="physig_npl.txt",sep="\t")

d<-read.table("phe_upl.txt",sep="\t")
var.all=matrix(0,nrow=ntrait,ncol=ncol(d))
for(i in 1:ntrait){
	dsub=d[which(d[,2]==i),]
	for(t in 3:ncol(dsub)){
		var.all[i,(t-2)]=var(dsub[,t])
	}
}
write.table(var.all,file="var_all_upl.txt",sep="\t")
cor.all=rep(0,ntrait)
for(i in 1:ntrait){
	cor.all[i]=cor(1:ncol(var.all),var.all[i,])
}
write.table(cor.all,file="physig_upl.txt",sep="\t")






