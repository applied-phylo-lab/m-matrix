# Analysis of results from simulations where one trait is under directional selection while all other traits are under stabilizing selection
# Plot mean of trait 1 (z1)


library(ggplot2)
library(resample) # Function to use: colVars

npl.all=c(0,10,20,30,40,50) # Number of pleiotropic loci
npop=50 # Number of populations per simulation
ntrait=10 # Number of traits
setwd(paste("./",ntrait,"t",sep="")) # Go to directory corresponding to the focal number of traits

# End-point result
out=matrix(0,nrow=length(npl.all),ncol=3);out[,1]=npl.all # Summary data file (columns: npl, cross-population mean z1, SE of cross-population mean z1)
for(npl in npl.all){
  fn=paste("sim_dir_out_end_",ntrait,"t_",npl,"pl.txt",sep="") # Data file to read
  d<-read.table(fn,sep="\t")
  # Arrange population means and variances into tables
  tab_mean=matrix(0,nrow=npop,ncol=ntrait) # Data matrix for means (each row for a population, each column for a trait)
  #tab_var=matrix(0,nrow=npop,ncol=ntrait) # Data matrix for G variance (each row for a population, each column for a trait)
  for(i in 1:npop){
    dsub=d[(ntrait*(i-1)+1):(ntrait*i),]
    tab_mean[i,]=dsub[,2]
    #tab_var[i,]=dsub[,3]
  }
  # Write R-matrix
  #rmat=cov(tab_mean)
  #fn_rmat=paste("rmat_dir_",ntrait,"t_",npl,"pl.txt",sep="")
  #write.table(data.frame(rmat),file=fn_rmat,sep="\t")
  # Write information into summary data matrices
  row=which(out[,1]==npl) # Find the row to write data into
  out[row,2]=mean(tab_mean[,1])
  out[row,3]=sd(tab_mean[,1])/sqrt(npop)
}
# Write output, generating data file used for plotting
fn=paste("sum_dir_z1_",ntrait,"t.txt",sep="")
write.table(out,file=fn,sep="\t")

# Temporal dynamics
# Process data files and remove the starting lines
z_temp=matrix(0,nrow=length(npl.all),ncol=100) # Data matrix for mean z1 over time
for(i in 1:length(npl.all)){
  fn=paste("sim_out_all_",ntrait,"t_",npl.all[i],".txt",sep="") # Data file to read
  d<-readLines(fn)
  d=d[(length(d)-100*npop*ntrait+1):length(d)]
  fn_new=paste("sim_out_all_",ntrait,"t_",npl.all[i],"_processed.txt",sep="")
  writeLines(d,fn_new)
  d<-read.table(fn_new,sep="\t") # Re-read rearranged file
  dsub=d[which(d[,3]==1),] # Extract rows corresponding to z1
  for(t in (1:100)*100){
    dt=dsub[which(dsub[,1]==t),]
    z_temp[i,t/100]=mean(dt[,4])
  }
}
# Arrange into a data frame for ggplot (columns: npl, time, mean z1)
out_temp_sum=matrix(0,nrow=length(npl.all)*100,ncol=3)
for(i in 1:length(npl.all)){
  # Rows of the output data matrix to fill
  row_start=100*(i-1)+1
  row_end=100*i
  out_temp_sum[row_start:row_end,1]=npl.all[i]
  out_temp_sum[row_start:row_end,2]=100*(1:100)
  out_temp_sum[row_start:row_end,3]=z_temp[i,]
}
# Write data file
colnames(out_temp_sum)=c("npl","t","z1") # Name columns before writing file
fn=paste("sum_temp_dir_",ntrait,"t.txt",sep="")
write.table(out_temp_sum,file=fn,sep="\t")
# Generate plots
d<-read.table(fn,sep="\t")
d$npl=factor(d$npl,levels=npl.all,ordered=TRUE) # Factorize
# Mean z1 variance over time
g<-ggplot(d,aes(x=t,y=z1,group=npl,col=npl))
g=g+geom_point(size=2)+geom_line(lwd=1)
g=g+theme_classic()
g=g+xlab("")+ylab("")
# Save plots
fn1=paste("plot_temp_z1_dir_",ntrait,"t.pdf",sep="")
ggsave(fn1,plot=g,width=7,height=5)

setwd("..")

# Plot end-point results
# Results from different trait numbers to be shown in one plot, so data files need to be combined
ntrait.all=c(2,5,10)
d=NULL
for(i in 1:length(ntrait.all)){
  setwd(paste("./",ntrait.all[i],"t",sep="")) # Go to the corresponding directory
  fn=paste("sum_dir_z1_",ntrait.all[i],"t.txt",sep="") # Find the file
  di<-read.table(fn,sep="\t")
  di=data.frame(rep(ntrait.all[i],nrow(di)),di) # Add a column of trait number
  if(i==1){
    d=di
  }else{
    d=rbind(d,di)
  }
  setwd("..") # Back to last directory
}
colnames(d)=c("ntrait","npl","z1","se") # Rename the first column
write.table(d,file="sum_dir.txt",sep="\t") # Write data file

d<-read.table("sum_dir.txt",sep="\t") # (Re-)read data file
d$ntrait=factor(d$ntrait,levels=ntrait.all,ordered=TRUE)
# Plot mean z1
g<-ggplot(d,aes(x=npl,y=z1,group=ntrait,col=ntrait))
g=g+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=z1-se,ymax=z1+se),width=1.5)
g=g+scale_color_manual(values=c("black", "orange","purple"))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")
# Save plots
ggsave("plot_dir_z1.pdf",plot=g,width=7,height=5)



