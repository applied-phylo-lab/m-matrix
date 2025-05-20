# Rate of beneficial mutations and the fraction of mutations that are beneficial

library(ggplot2)

fitness<-function(d,width){
	sub=which(width>0)
	D=sqrt(sum((d[sub]/width[sub])^2))
	w=dnorm(D,mean=0,sd=1)/dnorm(0,mean=0,sd=1)
	return(w)
}

ntrait.all=c(2,5,10)
L.per.trait=50
L.all=c(0,10,20,30,40,50)
u=1e-6;nmut=1e6

out=matrix(0,nrow=length(ntrait.all)*length(L.all),ncol=6)
for(i in 1:length(ntrait.all)){
	ntrait=ntrait.all[i]
	sig=rep(1,ntrait)
	opt=c(20,rep(0,ntrait-1));width=c(10,rep(1,ntrait-1))
	z0=rep(0,ntrait)
	w0=fitness(z0-opt,width)

	row=(i-1)*length(L.all)
	for(j in 1:length(L.all)){
		L=L.all[j]
		L.total=L+(L.per.trait-L)*ntrait
		nb=0
		for(i in 1:nmut){
			type=sample(c("z1","other","pl"),1,prob=c(L.per.trait-L,(L.per.trait-L)*(ntrait-1),L)/L.total)
			if(type=="z1"){
				nb=nb+1
			}else{
				if(type=="pl"){
					z=z0+rnorm(ntrait,mean=0,sd=sig)
					w=fitness(z-opt,width)
					if(w>w0){
						nb=nb+1
					}
				}
			}
		}
		fb=nb/nmut
		fb.sd=sqrt(fb*(1-fb)/nmut)
		rb=fb*L.total*u
		rb.sd=fb.sd*L.total*u

		row=row+1
		out[row,1]=ntrait
		out[row,2]=L
		out[row,3]=fb
		out[row,4]=fb.sd
		out[row,5]=rb
		out[row,6]=rb.sd
	}
}

d=data.frame(out)
colnames(d)=c("ntrait","L_upl","fb","fb.sd","rb","rb.sd")
fb.max=d$fb+d$fb.sd;fb.min=d$fb-d$fb.sd
rb.max=d$rb+d$rb.sd;rb.min=d$rb-d$rb.sd
d=data.frame(d,fb.max,fb.min,rb.min,rb.min)
d$ntrait=factor(d$ntrait,levels=c(2,5,10),ordered=TRUE)

g<-ggplot(d,aes(x=L_upl,y=rb,colour=ntrait))
g=g+geom_line(linewidth=2)+geom_point(size=3)+geom_errorbar(aes(ymin=rb.min,ymax=rb.max))
g=g+theme_classic()+xlab("")+ylab("")
g=g+scale_color_manual(values=c("gray15","orange","purple"))
ggsave(g,file="rate_beneficial.pdf")

g<-ggplot(d,aes(x=L_upl,y=fb,colour=ntrait))
g=g+geom_line(linewidth=2)+geom_point(size=3)+geom_errorbar(aes(ymin=fb.min,ymax=fb.max))
g=g+theme_classic()+xlab("")+ylab("")
g=g+scale_color_manual(values=c("gray15","orange","purple"))
ggsave(g,file="frac_beneficial.pdf")

