# Analysis of results from simulations where all traits are under stabilizing selection



library(ggplot2)
library(resample) # Function to use: colVars

npl.all=c(0,10,20,30,40,50) # Number of pleiotropic loci
npop=50 # Number of populations per simulation
ntrait=10 # Number of traits
setwd(paste("./",ntrait,"t",sep="")) # Go to directory corresponding to the focal number of traits

# End-point result
# Columns of each data file: trait, population mean, within-population variance
out=matrix(0,nrow=length(npl.all),ncol=3);out[,1]=npl.all # Summary data file; columns: number of pleiotropic loci, R variance (mean across traits), mean G variance (averaged across populations, then across traits)
out_vr=matrix(0,nrow=length(npl.all),ncol=ntrait) # Data matrix for R variances
out_vg=matrix(0,nrow=length(npl.all),ncol=ntrait) # Data matrix for mean G variances
sd_vr=matrix(0,nrow=length(npl.all),ncol=ntrait) # Sampling SD of R variance
se_vg=matrix(0,nrow=length(npl.all),ncol=ntrait) # SE of mean G variance
for(npl in npl.all){
  fn=paste("sim_out_end_",ntrait,"t_",npl,"pl.txt",sep="") # Data file to read
  d<-read.table(fn,sep="\t")
  # Arrange population means and variances into tables
  tab_mean=matrix(0,nrow=npop,ncol=ntrait)
  tab_var=matrix(0,nrow=npop,ncol=ntrait)
  for(i in 1:npop){
    dsub=d[(ntrait*(i-1)+1):(ntrait*i),]
    tab_mean[i,]=dsub[,2]
    tab_var[i,]=dsub[,3]
  }
  # Write R-matrix
  rmat=cov(tab_mean)
  fn_rmat=paste("rmat_",ntrait,"t_",npl,"pl.txt",sep="")
  write.table(data.frame(rmat),file=fn_rmat,sep="\t")
  # Write information into summary data matrices
  row=which(out[,1]==npl)
  out_vr[row,]=diag(rmat)
  out_vg[row,]=colMeans(tab_var)
  sd_vr[row,]=sqrt(2*(diag(rmat)^2)/(npop-1))
  se_vg[row,]=sqrt(colVars(tab_var)/npop)
  out[row,2]=mean(diag(rmat))
  out[row,3]=mean(colMeans(tab_var))
}
# Write output (trait-by-trait)
fn=paste("sum_all_stab_",ntrait,"t.txt",sep="")
write.table(data.frame(out_vr,out_vg,sd_vr,se_vg),file=fn,sep="\t")
# Write output (averaged across traits) & generate data frame used for plotting
out=data.frame(out,matrix(0,nrow=length(npl.all),ncol=2))
for(i in 1:nrow(out)){
  out[i,4]=mean(sd_vr[i,])
  out[i,5]=mean(se_vg[i,])
}
colnames(out)=c("npl","vr","vg","sd_vr","se_vg") # Name columns before writing file
fn=paste("sum_stab_",ntrait,"t.txt",sep="")
write.table(out,file=fn,sep="\t")

# Temporal dynamics
# Process data files and remove the starting lines
vr_temp_all=list() # R variance over time (a list of matrices; each matrix for a value of npl)
vg_temp_all=list() # Mean G variance over time (a list of matrices; each matrix for a value of npl)
for(i in 1:length(npl.all)){
  fn=paste("sim_out_all_",ntrait,"t_",npl.all[i],".txt",sep="") # Data file to read
  d<-readLines(fn)
  d=d[(length(d)-100*npop*ntrait+1):length(d)]
  fn_new=paste("sim_out_all_",ntrait,"t_",npl.all[i],"_processed.txt",sep="")
  writeLines(d,fn_new)
  d<-read.table(fn_new,sep="\t") # Re-read rearranged file
  vr_temp=matrix(0,nrow=ntrait,ncol=100) # R variance over time
  vg_temp=matrix(0,nrow=ntrait,ncol=100) # Mean G variance over time
  for(t in (1:100)*100){
    dt=d[which(d[,1]==t),]
    for(j in 1:ntrait){
      dt_sub=dt[which(dt[,3]==j),]
      vr_temp[j,t/100]=var(dt_sub[,4])
      vg_temp[j,t/100]=mean(dt_sub[,5])
    }
  }
  vr_temp_all[[i]]=vr_temp
  vg_temp_all[[i]]=vg_temp
}
# Arrange into a data frame for ggplot (columns: npl, time, R variance, G variance)
out_temp_sum=matrix(0,nrow=length(npl.all)*100,ncol=4)
for(i in 1:length(npl.all)){
  # Rows of the output data matrix to fill
  row_start=100*(i-1)+1
  row_end=100*i
  
  out_temp_sum[row_start:row_end,1]=npl.all[i]
  out_temp_sum[row_start:row_end,2]=100*(1:100)
  out_temp_sum[row_start:row_end,3]=colMeans(vr_temp_all[[i]])
  out_temp_sum[row_start:row_end,4]=colMeans(vg_temp_all[[i]])
}
# Write data file
colnames(out_temp_sum)=c("npl","t","vr","vg") # Name columns before writing file
fn=paste("sum_temp_stab_",ntrait,"t.txt",sep="")
write.table(out_temp_sum,file=fn,sep="\t")
# Generate plots
d<-read.table(fn,sep="\t")
d$npl=factor(d$npl,levels=npl.all,ordered=TRUE) # Factorize
# R variance over time
g<-ggplot(d,aes(x=t,y=log10(vr),group=npl,col=npl))
g=g+geom_point(size=2)+geom_line(lwd=1)
g=g+theme_classic()
g=g+xlab("")+ylab("")
# G variance over time
h<-ggplot(d,aes(x=t,y=log10(vg),group=npl,col=npl))
h=h+geom_point(size=2)+geom_line(lwd=1)
h=h+theme_classic()
h=h+xlab("")+ylab("")
# Save plots
fn1=paste("plot_temp_vr_stab_",ntrait,"t.pdf",sep="")
fn2=paste("plot_temp_vg_stab_",ntrait,"t.pdf",sep="")
ggsave(fn1,plot=g,width=7,height=5)
ggsave(fn2,plot=h,width=7,height=5)

setwd("..")

# Plot end-point results
# Results from different trait numbers to be shown in one plot, so data files need to be combined
ntrait.all=c(2,5,10)
d=NULL
for(i in 1:length(ntrait.all)){
  setwd(paste("./",ntrait.all[i],"t",sep="")) # Go to the corresponding directory
  fn=paste("sum_stab_",ntrait.all[i],"t.txt",sep="") # Find the file
  di<-read.table(fn,sep="\t")
  di=data.frame(rep(ntrait.all[i],nrow(di)),di) # Add a column of trait number
  if(i==1){
    d=di
  }else{
    d=rbind(d,di)
  }
  setwd("..") # Back to last directory
}
colnames(d)[1]="ntrait" # Rename the first column
write.table(d,file="sum_stab.txt",sep="\t") # Write data file

d<-read.table("sum_stab.txt",sep="\t") # (Re-)read data file
d$ntrait=factor(d$ntrait,levels=ntrait.all,ordered=TRUE)
# Plot G variance
g<-ggplot(d,aes(x=npl,y=log10(vg),group=ntrait,col=ntrait))
g=g+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=log10(vg-se_vg),ymax=log10(vg+se_vg)),width=1.5)
g=g+scale_color_manual(values=c("purple","darkgrey", "orange"))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")
# Plot R variance
h<-ggplot(d,aes(x=npl,y=log10(vr),group=ntrait,col=ntrait))
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=log10(vr-sd_vr),ymax=log10(vr+sd_vr)),width=1.5)
h=h+scale_color_manual(values=c("purple","darkgrey", "orange"))
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")
# Save plots
ggsave("plot_stab_vg.pdf",plot=g,width=7,height=5)
ggsave("plot_stab_vr.pdf",plot=h,width=7,height=5)



