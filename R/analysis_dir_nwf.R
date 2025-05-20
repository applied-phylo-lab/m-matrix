# Analysis of results from simulations where one trait is under directional selection while all other traits are under stabilizing selection (non-WF simulations)
# Plot mean of trait 1 (z1)


library(ggplot2)
library(resample) # Function to use: colVars

npl.all=c(0,10,20,30,40,50) # Number of pleiotropic loci
npop=50 # Number of populations per simulation
ntrait=2 # Number of traits
setwd(paste("./",ntrait,"t",sep="")) # Go to directory corresponding to the focal number of traits

# End-point result
# Note: these files have population no. as the first column
out=matrix(0,nrow=length(npl.all),ncol=5);out[,1]=npl.all # Summary data file (columns: npl, cross-population mean z1, SE of cross-population mean z1, mean N, SE of mean N)
for(npl in npl.all){
  fn=paste("sim_dir_nwf_out_end_",ntrait,"t_",npl,"pl.txt",sep="") # Data file to read
  d<-read.table(fn,sep="\t")
  # Arrange population means and variances into tables
  tab_mean=matrix(0,nrow=npop,ncol=ntrait) # Data matrix for means (each row for a population, each column for a trait)
  #tab_var=matrix(0,nrow=npop,ncol=ntrait) # Data matrix for G variance (each row for a population, each column for a trait)
  tab_N=rep(0,npop)
  for(i in 1:npop){
    dsub=d[(ntrait*(i-1)+1):(ntrait*i),]
    tab_mean[i,]=dsub[,4]
    #tab_var[i,]=dsub[,5]
    tab_N[i]=dsub[1,3]
  }
  # Write R-matrix
  #rmat=cov(tab_mean)
  #fn_rmat=paste("rmat_dir_",ntrait,"t_",npl,"pl.txt",sep="")
  #write.table(data.frame(rmat),file=fn_rmat,sep="\t")
  # Write information into summary data matrices
  row=which(out[,1]==npl) # Find the row to write data into
  out[row,2]=mean(tab_mean[,1])
  out[row,3]=sd(tab_mean[,1])/sqrt(npop)
  out[row,4]=mean(tab_N)
  out[row,5]=sd(tab_N)/sqrt(npop)
}
# Write output, generating data file used for plotting
fn=paste("sum_dir_nwf_z1_",ntrait,"t.txt",sep="")
write.table(out,file=fn,sep="\t")

# Temporal dynamics
# Process data files and remove the starting lines
z_temp=matrix(0,nrow=length(npl.all),ncol=100) # Data matrix for mean z1 over time
N_temp=matrix(0,nrow=length(npl.all),ncol=100) # Data matrix for mean N over time
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
    z_temp[i,t/100]=mean(dt[,5])
    N_temp[i,t/100]=mean(dt[,4])
  }
}
# Arrange into a data frame for ggplot (columns: npl, time, mean z1, mean N)
out_temp_sum=matrix(0,nrow=length(npl.all)*100,ncol=4)
for(i in 1:length(npl.all)){
  # Rows of the output data matrix to fill
  row_start=100*(i-1)+1
  row_end=100*i
  out_temp_sum[row_start:row_end,1]=npl.all[i]
  out_temp_sum[row_start:row_end,2]=100*(1:100)
  out_temp_sum[row_start:row_end,3]=z_temp[i,]
  out_temp_sum[row_start:row_end,4]=N_temp[i,]
}
# Write data file
colnames(out_temp_sum)=c("npl","t","z1","N") # Name columns before writing file
fn=paste("sum_temp_dir_nwf_",ntrait,"t.txt",sep="")
write.table(out_temp_sum,file=fn,sep="\t")
# Generate plots
d<-read.table(fn,sep="\t")
d$npl=factor(d$npl,levels=npl.all,ordered=TRUE) # Factorize
# Mean z1 variance over time
g<-ggplot(d,aes(x=t,y=z1,group=npl,col=npl))
g=g+geom_point(size=2)+geom_line(lwd=1)
g=g+theme_classic()
g=g+xlab("")+ylab("")
# Mean N over time
h<-ggplot(d,aes(x=t,y=N,group=npl,col=npl))
h=h+geom_point(size=2)+geom_line(lwd=1)
h=h+theme_classic()
h=h+scale_y_continuous(breaks=c(0,200,400,600,800,1000),limits=c(0,1100))
h=h+xlab("")+ylab("")
# Save plots
fn1=paste("plot_temp_z1_dir_nwf_",ntrait,"t.pdf",sep="")
ggsave(fn1,plot=g,width=7,height=5)
fn2=paste("plot_temp_N_dir_nwf_",ntrait,"t.pdf",sep="")
ggsave(fn2,plot=h,width=7,height=5)

setwd("..")

# Plot end-point results
# Results from different trait numbers to be shown in one plot, so data files need to be combined
ntrait.all=c(2,5,10)
d=NULL
for(i in 1:length(ntrait.all)){
  setwd(paste("./",ntrait.all[i],"t",sep="")) # Go to the corresponding directory
  fn=paste("sum_dir_nwf_z1_",ntrait.all[i],"t.txt",sep="") # Find the file
  di<-read.table(fn,sep="\t")
  di=data.frame(rep(ntrait.all[i],nrow(di)),di) # Add a column of trait number
  if(i==1){
    d=di
  }else{
    d=rbind(d,di)
  }
  setwd("..") # Back to last directory
}
colnames(d)=c("ntrait","npl","z1","se_z1","N","se_N") # Rename the first column
write.table(d,file="sum_dir_nwf.txt",sep="\t") # Write data file

d<-read.table("sum_dir_nwf.txt",sep="\t") # (Re-)read data file
d$ntrait=factor(d$ntrait,levels=ntrait.all,ordered=TRUE)
# Plot mean z1
g<-ggplot(d,aes(x=npl,y=z1,group=ntrait,col=ntrait))
g=g+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=z1-se_z1,ymax=z1+se_z1),width=1.5)
g=g+scale_color_manual(values=c("black", "orange","purple"))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")
# Plot N
h<-ggplot(d,aes(x=npl,y=N,group=ntrait,col=ntrait))
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=N-se_N,ymax=N+se_N),width=1.5)
h=h+scale_color_manual(values=c("black", "orange","purple"))
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+scale_y_continuous(breaks=c(0,200,400,600,800,1000),limits=c(0,1050))
h=h+xlab("")+ylab("")
# Save plots
ggsave("plot_dir_nwf_z1.pdf",plot=g,width=7,height=5)
ggsave("plot_dir_nwf_N.pdf",plot=h,width=7,height=5)



