# Obtain an identical M matrix by manipulate a second parameter given the first parameter value change

# Calculate the analytical solution of the observed mutational correlation
m_as <- function(par){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	m1=U1*rbind(c(sig1^2,0),c(0,0))
	m2=U2*rbind(c(0,0),c(0,sig2^2))
	mp=Up*rbind(c(siga^2,r*siga*sigb),c(r*siga*sigb,sigb^2))
	m=m1+m2+mp
	return(m)
}

# Conversion type 1
# Change Up and r such that mutational covariance is unchanged, then change values of U1, sig1, U2, and/or sig2 to keep variances unchanged
# Arguments: the original parameter values (par), change to r (r becomes c*r), and what other parameter values to change (v, a vector)
# Components of par: U1, sig1, U2, sig2, Up, r, siga, sigb
par.convert <- function(par,c,v){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	if(c>(1/abs(par[6]))){
		return("Error")
	}else{
		par.new=par
		par.new[6]=par.new[6]*c # Rescale r;
		par.new[5]=par.new[5]/c # Rescale Up
		c1=(U1*(sig1^2)+(1-(1/c))*Up*(siga^2))/(U1*(sig1^2))
		c2=(U2*(sig2^2)+(1-(1/c))*Up*(sigb^2))/(U2*(sig2^2))
		# Rescale U1 and U2
		par.new[v[1]]=par.new[v[1]]*c1
		par.new[v[2]]=par.new[v[2]]*c2
		return(par.new)
	}
}

# Conversion type 2
# Change Up, siga, and sigb such that mutational covariance is unchanged, then change values of U1, sig1, U2, and/or sig2 to keep variances unchanged
# Arguments: the original parameter values (par), change to Up (Up becomes Up*r). (No 3rd argument because no other parameters need to change)
# Components of par: U1, sig1, U2, sig2, Up, r, siga, sigb
par.convert <- function(par,c){
	U1=par[1];sig1=par[2];U2=par[3];sig2=par[4];Up=par[5];r=par[6];siga=par[7];sigb=par[8]
	par.new=par
	par.new[5]=par.new[5]*c # Rescale Up
	par.new[7]=par.new[7]/sqrt(c) # Rescale siga
	par.new[8]=par.new[8]/sqrt(c) # Rescale sigb
	return(par.new)
}


