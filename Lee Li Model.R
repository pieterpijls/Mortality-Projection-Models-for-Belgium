################################################
### Question 3
################################################
# set working directory
setwd("C:/PIETER/Pieter/KUL/ADVANCED LIFE/Assignment 2")

# Install packages
library(reshape)
library(ggplot2)
library(gridExtra)
library(cowplot)
# define some parameters
Country=c("BE")
Sex=c("Male", "Female")
xv=c(0:101) #we do not take 110+ into account because there are almost no datapoints
yv=c(1960:2015)

nage = length(xv)
nyear = length(yv)
rname = years = yv[1]:yv[length(yv)]
cname = 0:(length(xv)-1)

### read in data from HMD ###

Deaths_BE = read.table("Deaths_1x1_BE.txt",header=TRUE,dec=".",sep="",skip=2)
Exposures_BE = read.table("Exposures_1x1_BE.txt",header=TRUE,dec=".",sep="",skip=2)

# put data into matrix format
DeathMat_BE_M = matrix(data=0, ncol=nage, nrow=nyear, dimnames=list(c(rname),c(cname)) )
DeathMat_BE_F = matrix(data=0, ncol=nage, nrow=nyear, dimnames=list(c(rname),c(cname)) )

ExpMat_BE_M = matrix(data=0, ncol=nage, nrow=nyear, dimnames=list(c(rname),c(cname)) )
ExpMat_BE_F = matrix(data=0, ncol=nage, nrow=nyear, dimnames=list(c(rname),c(cname)) )

x = Deaths_BE[,-5]
xMmatrix = x[,-3]
xFmatrix = x[,-4]
for(i in 1:nyear){
  for(j in 1:nage){
    DeathMat_BE_M[i,j] = as.numeric(as.character(subset(xMmatrix, Year==yv[i] & Age==xv[j])$Male))
    DeathMat_BE_F[i,j] = as.numeric(as.character(subset(xFmatrix, Year==yv[i] & Age==xv[j])$Female))
  }
}

x = Exposures_BE[,-5]
xMmatrix = x[,-3]
xFmatrix = x[,-4]
for(i in 1:nyear){
  for(j in 1:nage){
    ExpMat_BE_M[i,j] = as.numeric(as.character(subset(xMmatrix, Year==yv[i] & Age==xv[j])$Male))
    ExpMat_BE_F[i,j] = as.numeric(as.character(subset(xFmatrix, Year==yv[i] & Age==xv[j])$Female))
  }
}
# Avoid problem of taking logarithm of zero
DeathMat_BE_M[DeathMat_BE_M==0]=0.01
DeathMat_BE_F[DeathMat_BE_F==0]=0.01

ExpMat_BE_M[ExpMat_BE_M==0]=0.01
ExpMat_BE_F[ExpMat_BE_F==0]=0.01

dtx_M <- DeathMat_BE_M
etx_M <- ExpMat_BE_M
dtx_F <- DeathMat_BE_F
etx_F <- ExpMat_BE_F
dtx_ALL <- DeathMat_BE_M + DeathMat_BE_F
etx_ALL <- ExpMat_BE_M + ExpMat_BE_F
# Construct weight matrix (1 = data available, 0 = data not available)

wa  = matrix(data=1,ncol=dim(etx_M)[2],nrow=dim(etx_M)[1]);
I = dim(wa)[1]
J = dim(wa)[2]
for(i in 1:I){
  for(j in 1:J){
    if(is.na(dtx_M[i,j])){
      wa[i,j] = 0
    }
  }
}
dimnames(wa) = dimnames(etx_M)

# Fit Multipopulation model Li-Lee (M1-M1) 
#
# function FitLiLee_BIS(): 
##		log m(t,x) = beta1(x) + beta2(x)*kappa2(t) + Poisson error
# OUTPUT:
# common factors: ...$SVD$A.x, ...$NR$A.x (zero if Ax=FALSE), ...$SVD$B.x, ...$SVD$K.t or ...$NR$B.x, ...$NR$K.t
# country specific parameters: for example ...$SVD$a.xi$BE, ...$SVD$b.xi$BE, ...$SVD$k.ti$BE	
# country specific mortality:	for example ...$SVD$m.txi$BE

## STEP1: GET COMMON FACTORS B.x AND K.t FROM LEE-CARTER FIT
# fill in parameters function 
int = xv*0  
constraints = "ALL"  
exclude.coh = FALSE

dtx <- dtx_ALL
etx <- etx_ALL
  
  mtx	= dtx/etx      			# matrix of death rates
  qtx	= 1-exp(-mtx)   			# matrix of mortality rates
  
  n = length(xv)				# number of ages
  m = length(yv)				# number of years
  
  cy= (yv[1]-xv[n]):(yv[m]-xv[1])  	# cohort approximate years of birth
  
  # Initialise parameter vectors
  beta1v=int
  beta2v=(1:n)*0
  beta3v=(1:n)*0			# dummy vector, this will stay at 0
  kappa2v=(1:m)*0
  gamma3v=(1:(n+m-1))*0		# dummy vector, this will stay at 0
  ia=array((1:m),c(m,n))	# matrix of year indexes, i, for the data
  ja=t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
  ya=ia-ja		 		# matrix of year of birth indexes for the data
  imj=(1-n):(m-1)			# the range of values taken by i-j
  lg=n+m-1		 		# number of different values taken by i-j
  ca=ya+yv[1]-xv[1]		# matrix of years of birth
  
  
  ww=cy*0+1	 # this is a vector of 1's and 0's with
  # a 0 if the cohort is completely excluded
  #   for(k in 1:lg)
  #   {
  #     ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
  #   }
  
  # Stage 0
  # Gives initial estimates for beta1(x), beta2(x) and kappa2(t)
  mx=mean(xv)
  for(j in 1:n)
  {
    if(all(int==0)){
      beta1v[j]=sum(log(mtx[,j])*wa[,j])/sum(wa[,j])
    }
    beta2v[j]=1/n
  }
  kappa2v=(m:1)-(m+1)/2

  llmaxM2B=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
    #   b1,b3,k2,g3 are given
    #   solve for b2
    if(min(dv=0)) {
      print(b1)
      print(b2)
      print(b3)
      print(k2)
      print(g3)
      print(dv)
      print(ev)}
    b21=b2
    b20=b21-1
    thetat=k2*ev*exp(b1+b3*g3)
    s1=sum(dv*k2*wv)
    
    while(abs(b21-b20) > 1e-10)
    {
      b20=b21; 
      f0=sum((exp(b20*k2)*thetat)*wv)-s1;
      df0=sum((exp(b20*k2)*k2*thetat)*wv)
      b21=b20-f0/df0
    }
    b21
  }
  
  llmaxM2D=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
    #   b1,b2,b3,g3 are given
    #   solve for k2
    k21=k2
    k20=k21-1
    thetat=b2*ev*exp(b1+b3*g3)
    s1=sum(dv*b2*wv)
    while(abs(k21-k20) > 1e-10)
    {
      k20=k21
      f0=sum((exp(k20*b2)*thetat)*wv)-s1
      df0=sum((exp(k20*b2)*b2*thetat)*wv)
      k21=k20-f0/df0
    }
    k21
  }  
  
  
  
    
  # Stage 1: iterate
  l0=-1000000
  l1=-999999
  iteration=0
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001
  while(abs(l1-l0) > 1e-20)
  {
    iteration=iteration+1
    
    l0=l1
    
    #     mhat=mtx*0
    #     for(i in 1:m)
    #     {
    #       mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    #     }
    #     epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    #     l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    #     
    #    if(iteration==1) 
    # Stage 1B optimise over the beta2(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]	# exposure
      beta2v[j]=llmaxM2B(beta1v[j],beta2v[j],beta3v[j],
                         kappa2v,gamma3v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    #cat(l1,"->")
    
    # Stage 1D optimise over the kappa2(t)
    for(i in 1:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]	# exposure
      kappa2v[i]=llmaxM2D(beta1v,beta2v,beta3v,
                          kappa2v[i],gamma3v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    #cat(l1,"->")	
    
    # Now apply the constraints
    if(constraints=="ALL"){
      fac21=mean(kappa2v)
      fac22=sum(beta2v)
      kappa2v=fac22*(kappa2v-fac21)    # sum kappa2=0
      beta2v=beta2v/fac22             # sum beta2=1
      beta1v=beta1v+beta2v*fac22*fac21 # => adjust beta1 
    }
    else{
      fac22=sum(beta2v)
      kappa2v=fac22*kappa2v           # => adjust kappa2
      beta2v=beta2v/fac22  		  # sum beta2=1
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    #cat(l1,"->")
    
    
    # Stage 1A optimise over the beta1(x)
    if(all(int==0)){
      for(j in 1:n)
      {		 		 
        # cycle through the range of years
        wv=1	    # can be set to a vector of weights
        # to e.g. exclude duff years
        wv=wa[,j]
        s1=sum(wv*dtx[,j])
        s2=sum(wv*etx[,j]*exp(beta2v[j]*kappa2v+beta3v[j]*gamma3v[(n+1-j):(n+m-j)]))
        beta1v[j]=log(s1)-log(s2)
      }
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    #cat(l1,"->")
    
    
  }		 # end while loop
  
  # calculate number of parameters and deduct the number of constraints
  # also count number of parameters in 'beta1v'!

  
    npar_tot=length(beta1v)+length(beta2v)+length(kappa2v)-2

  
  # Calculate the BIC
  BIC=-2*l1+log(sum(wa))*npar_tot
  
  A.x = beta1v
  B.x = beta2v		# note: sum B.x = 1
  K.t = kappa2v		# note: sum K.t = 0 


  
## STEP2: GET POPULATION SPECIFIC PARAMETERS (NR LEE-CARTER) +  POPULATION SPECIFIC DEATH RATES
  alpha.x	= c()	# initialize list alphas
  beta.x	= c()	# initialize list betas
  kappa.t	= c()	# initialize list kappas
  m.txi	= c() 	# initialize list death rates
  
  #fit701(xv, yvSPEC, (Data$UNI[[c]]$etx) * exp(rep(A.x, each=length(yvSPEC), times=1) + K.t %o% B.x) , Data$UNI[[c]]$dtx , Data$UNI[[c]]$wa ,xv*0,"ALL",exclude.coh) # adjust exposure!

  #fit701 = function(xv,yv,etx,dtx,wa,int,constraints,exclude.coh){
    yvSPEC = yv
    etx <- etx_M*exp(rep(A.x, each=length(yvSPEC), times=1) + K.t %o% B.x)
    dtx <- dtx_M
    int <- xv*0
    constraints = "ALL"
    exclude.coh = FALSE
    
    mtx	= dtx/etx      			# matrix of death rates
    qtx	= 1-exp(-mtx)   			# matrix of mortality rates
    
    n = length(xv)				# number of ages
    m = length(yv)				# number of years
    
    cy= (yv[1]-xv[n]):(yv[m]-xv[1])  	# cohort approximate years of birth
    
    # Initialise parameter vectors
    beta1v=int
    beta2v=(1:n)*0
    beta3v=(1:n)*0			# dummy vector, this will stay at 0
    kappa2v=(1:m)*0
    gamma3v=(1:(n+m-1))*0		# dummy vector, this will stay at 0
    ia=array((1:m),c(m,n))	# matrix of year indexes, i, for the data
    ja=t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
    ya=ia-ja		 		# matrix of year of birth indexes for the data
    imj=(1-n):(m-1)			# the range of values taken by i-j
    lg=n+m-1		 		# number of different values taken by i-j
    ca=ya+yv[1]-xv[1]		# matrix of years of birth
    
    
    ww=cy*0+1	 # this is a vector of 1's and 0's with
    # a 0 if the cohort is completely excluded
    #   for(k in 1:lg)
    #   {
    #     ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
    #   }
    
    # Stage 0
    # Gives initial estimates for beta1(x), beta2(x) and kappa2(t)
    mx=mean(xv)
    for(j in 1:n)
    {
      if(all(int==0)){
        beta1v[j]=sum(log(mtx[,j])*wa[,j])/sum(wa[,j])
      }
      beta2v[j]=1/n
    }
    kappa2v=(m:1)-(m+1)/2
    
    # Stage 1: iterate
    l0=-1000000
    l1=-999999
    iteration=0
    # l1 is the latest estimate of the log-likelihood
    # l0 is the previous estimate
    # we continue to iterate if the improvement in log-likelihood
    # exceeds 0.0001
    while(abs(l1-l0) > 1e-20)
    {
      iteration=iteration+1
      
      l0=l1
      
      #     mhat=mtx*0
      #     for(i in 1:m)
      #     {
      #       mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      #     }
      #     epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      #     l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #     
      #    if(iteration==1) 
      # Stage 1B optimise over the beta2(x)
      for(j in 1:n)
      {		 		 
        # cycle through the range of years
        dv=dtx[,j]	# actual deaths
        ev=etx[,j]	# exposure
        beta2v[j]=llmaxM2B(beta1v[j],beta2v[j],beta3v[j],
                           kappa2v,gamma3v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")
      
      # Stage 1D optimise over the kappa2(t)
      for(i in 1:m)
      {		 		 
        # cycle through the range of years
        dv=dtx[i,]	# actual deaths
        ev=etx[i,]	# exposure
        kappa2v[i]=llmaxM2D(beta1v,beta2v,beta3v,
                            kappa2v[i],gamma3v[(n+i-1):i],dv,ev,wv=wa[i,])
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")	
      
      # Now apply the constraints
      if(constraints=="ALL"){
        fac21=mean(kappa2v)
        fac22=sum(beta2v)
        kappa2v=fac22*(kappa2v-fac21)    # sum kappa2=0
        beta2v=beta2v/fac22             # sum beta2=1
        beta1v=beta1v+beta2v*fac22*fac21 # => adjust beta1 
      }
      else{
        fac22=sum(beta2v)
        kappa2v=fac22*kappa2v           # => adjust kappa2
        beta2v=beta2v/fac22  		  # sum beta2=1
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")
      
      
      # Stage 1A optimise over the beta1(x)
      if(all(int==0)){
        for(j in 1:n)
        {		 		 
          # cycle through the range of years
          wv=1	    # can be set to a vector of weights
          # to e.g. exclude duff years
          wv=wa[,j]
          s1=sum(wv*dtx[,j])
          s2=sum(wv*etx[,j]*exp(beta2v[j]*kappa2v+beta3v[j]*gamma3v[(n+1-j):(n+m-j)]))
          beta1v[j]=log(s1)-log(s2)
        }
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")
      
      
    }		 # end while loop
    
    # calculate number of parameters and deduct the number of constraints
    # also count number of parameters in 'beta1v'!
    npar=length(beta1v)+length(beta2v)+length(kappa2v)-2
    
    
    # Calculate the BIC
    BIC=-2*l1+log(sum(wa))*npar
    
    alpha.x_M	= beta1v
    beta.x_M 	= beta2v		# note: sum = 1
    kappa.t_M 	= kappa2v	 		# note: sum = 0
    ll_M <- l1
    npar_M = npar + npar_tot
    BIC_M = -2*ll_M+log(sum(wa))*npar_M
    residuals_M = epsilon

    # death rates Males
    m.tx_M = t(exp(rep(A.x, each=1, times=length(yvSPEC)) + rep(alpha.x_M, each=1, times=length(yvSPEC)) + B.x%o%K.t + beta.x_M%o%kappa.t_M))

    
    #fit701 = function(xv,yv,etx,dtx,wa,int,constraints,exclude.coh){
    yvSPEC = yv
    etx <- etx_F*exp(rep(A.x, each=length(yvSPEC), times=1) + K.t %o% B.x)
    dtx <- dtx_F
    int <- xv*0
    constraints = "ALL"
    exclude.coh = FALSE
    
    mtx	= dtx/etx      			# matrix of death rates
    qtx	= 1-exp(-mtx)   			# matrix of mortality rates
    
    n = length(xv)				# number of ages
    m = length(yv)				# number of years
    
    cy= (yv[1]-xv[n]):(yv[m]-xv[1])  	# cohort approximate years of birth
    
    # Initialise parameter vectors
    beta1v=int
    beta2v=(1:n)*0
    beta3v=(1:n)*0			# dummy vector, this will stay at 0
    kappa2v=(1:m)*0
    gamma3v=(1:(n+m-1))*0		# dummy vector, this will stay at 0
    ia=array((1:m),c(m,n))	# matrix of year indexes, i, for the data
    ja=t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
    ya=ia-ja		 		# matrix of year of birth indexes for the data
    imj=(1-n):(m-1)			# the range of values taken by i-j
    lg=n+m-1		 		# number of different values taken by i-j
    ca=ya+yv[1]-xv[1]		# matrix of years of birth
    
    
    ww=cy*0+1	 # this is a vector of 1's and 0's with
    # a 0 if the cohort is completely excluded
    #   for(k in 1:lg)
    #   {
    #     ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
    #   }
    
    # Stage 0
    # Gives initial estimates for beta1(x), beta2(x) and kappa2(t)
    mx=mean(xv)
    for(j in 1:n)
    {
      if(all(int==0)){
        beta1v[j]=sum(log(mtx[,j])*wa[,j])/sum(wa[,j])
      }
      beta2v[j]=1/n
    }
    kappa2v=(m:1)-(m+1)/2
    
    # Stage 1: iterate
    l0=-1000000
    l1=-999999
    iteration=0
    # l1 is the latest estimate of the log-likelihood
    # l0 is the previous estimate
    # we continue to iterate if the improvement in log-likelihood
    # exceeds 0.0001
    while(abs(l1-l0) > 1e-20)
    {
      iteration=iteration+1
      
      l0=l1
      
      #     mhat=mtx*0
      #     for(i in 1:m)
      #     {
      #       mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      #     }
      #     epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      #     l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #     
      #    if(iteration==1) 
      # Stage 1B optimise over the beta2(x)
      for(j in 1:n)
      {		 		 
        # cycle through the range of years
        dv=dtx[,j]	# actual deaths
        ev=etx[,j]	# exposure
        beta2v[j]=llmaxM2B(beta1v[j],beta2v[j],beta3v[j],
                           kappa2v,gamma3v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")
      
      # Stage 1D optimise over the kappa2(t)
      for(i in 1:m)
      {		 		 
        # cycle through the range of years
        dv=dtx[i,]	# actual deaths
        ev=etx[i,]	# exposure
        kappa2v[i]=llmaxM2D(beta1v,beta2v,beta3v,
                            kappa2v[i],gamma3v[(n+i-1):i],dv,ev,wv=wa[i,])
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")	
      
      # Now apply the constraints
      if(constraints=="ALL"){
        fac21=mean(kappa2v)
        fac22=sum(beta2v)
        kappa2v=fac22*(kappa2v-fac21)    # sum kappa2=0
        beta2v=beta2v/fac22             # sum beta2=1
        beta1v=beta1v+beta2v*fac22*fac21 # => adjust beta1 
      }
      else{
        fac22=sum(beta2v)
        kappa2v=fac22*kappa2v           # => adjust kappa2
        beta2v=beta2v/fac22  		  # sum beta2=1
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")
      
      
      # Stage 1A optimise over the beta1(x)
      if(all(int==0)){
        for(j in 1:n)
        {		 		 
          # cycle through the range of years
          wv=1	    # can be set to a vector of weights
          # to e.g. exclude duff years
          wv=wa[,j]
          s1=sum(wv*dtx[,j])
          s2=sum(wv*etx[,j]*exp(beta2v[j]*kappa2v+beta3v[j]*gamma3v[(n+1-j):(n+m-j)]))
          beta1v[j]=log(s1)-log(s2)
        }
      }
      
      mhat=mtx*0
      for(i in 1:m)
      {
        mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
      }
      epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
      l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
      #cat(l1,"->")
      
      
    }		 # end while loop
    
    # calculate number of parameters and deduct the number of constraints
    # also count number of parameters in 'beta1v'!
    
    npar=length(beta1v)+length(beta2v)+length(kappa2v)-2
    
    
    # Calculate the BIC
    BIC=-2*l1+log(sum(wa))*npar
    
    alpha.x_F	= beta1v
    beta.x_F 	= beta2v		# note: sum = 1
    kappa.t_F 	= kappa2v	 		# note: sum = 0
    ll_F <- l1
    npar_F = npar + npar_tot
    BIC_F = -2*ll_F+log(sum(wa))*npar_F
    residuals_F = epsilon
    
    # death rates Females
    m.tx_F = t(exp(rep(A.x, each=1, times=length(yvSPEC)) + rep(alpha.x_F, each=1, times=length(yvSPEC)) + B.x%o%K.t + beta.x_F%o%kappa.t_F))
    
        
###### Plot fitted parameters
    par(mfrow=c(3,3))
    plot(A.x,type="l",xlab="Age",ylab="ALL",main=expression(A[x]))
    plot(B.x,type="l",xlab="Age",ylab=" ",main=expression(B[x]))
    plot(yv, K.t,type="l",xlab="Year",ylab=" ",main=expression(K[t]))
    plot(alpha.x_M,type="l",xlab="Age",ylab="Males",main=expression(alpha[x]^M))
    plot(beta.x_M,type="l",xlab="Age",ylab=" ",main=expression(beta[x]^M))
    plot(yv, kappa.t_M,type="l",xlab="Year",ylab=" ",main=expression(kappa[t]^M))
    plot(alpha.x_F,type="l",xlab="Age",ylab="Females",main=expression(alpha[x]^F))
    plot(beta.x_F,type="l",xlab="Age",ylab=" ",main=expression(beta[x]^F))
    plot(yv, kappa.t_F,type="l",xlab="Year",ylab=" ",main=expression(kappa[t]^F))
    
    #residuals heatmap
    library('fields')
    library('RColorBrewer')
    cols = rev(brewer.pal(11,'RdYlBu'))
    cols = c(cols,rep(cols[11],5))
    rf <- colorRampPalette(cols)   # make colors
    r <- rf(64)
    yv_fig <- c(yv[1],yv[6],yv[11],yv[16],yv[21],yv[26],yv[31],yv[36],yv[41],yv[46],yv[51],yv[56])
    xv_fig <- c(xv[1],xv[6],xv[11],xv[16],xv[21],xv[26],xv[31],xv[36],xv[41],xv[46],xv[51],xv[56],
                xv[61],xv[66],xv[71],xv[76],xv[81],xv[86],xv[91],xv[96],xv[101])
    
    par(fig=c(0,0.45,0,1))
    image(t(residuals_M),xlab="Age",ylab="Year",col=r,main="Male",axes=F)
    mtext(text=yv_fig, side=2, line=0.3, at=seq(0,1,5/length(yv)), cex=0.8,las=1)
    mtext(text=xv_fig, side=1, line=0.3, at=seq(0,1,5/length(xv)), las=2, cex=0.8)
    par(fig=c(0.1,0.5,0,1),new=T)
    image.plot(t(residuals_M),legend.only=T)
    par(fig=c(0.52,0.95,0,1),new=T)
    image(t(residuals_F),xlab="Age",ylab="Year",col=r,main="Female",axes=F)
    mtext(text=yv_fig, side=2, line=0.3, at=seq(0,1,5/length(yv)), cex=0.8,las=1)
    mtext(text=xv_fig, side=1, line=0.3, at=seq(0,1,5/length(xv)), las=2, cex=0.8)
    par(fig=c(0.6,1,0,1),new=T)
    image.plot(t(residuals_F),legend.only=T)
    
    #######################""
    # Male 
    dataplot1 <- melt(residuals_M)
    colnames(dataplot1) <- c("Year","Age","Residual")
    
    cols <- rev(rainbow(9)[-9])
    hex <- c("#FF0000", "#FFA500", "#FFFF00", "#008000", "#9999ff", "#000066")

    
    p1 <- ggplot(data = dataplot1, aes(x = dataplot1$Age, y = dataplot1$Year))  + geom_tile(aes(fill = dataplot1$Residual)) + scale_fill_gradient2() 
    p1 <- p1  + scale_fill_gradientn(colours = r, limits=c(-5,15), name="Value")
    p1 <- p1 + expand_limits(x = 0, y = 1960) + ylim(1960,2015) + xlim(0,101)

    p1 <- p1 + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    p1 <- p1 + xlab("Age") + ylab("Year") + ggtitle("Male")+theme(plot.title = element_text(lineheight=.8, face="bold"))+theme(plot.title = element_text(hjust = 0.5))
    p1
    
    
    # Female
    dataplot2 <- melt(residuals_F)
    colnames(dataplot2) <- c("Year","Age","Residual")
    
    cols <- rev(rainbow(6)[-6])
    hex <- c("#FF0000", "#FFA500", "#FFFF00", "#008000", "#9999ff", "#000066")

    
    p2 <- ggplot(data = dataplot2, aes(x = dataplot2$Age, y = dataplot2$Year))  + geom_tile(aes(fill = dataplot2$Residual)) + scale_fill_gradient2() 
    p2 <- p2  + scale_fill_gradientn(colours = r, limits=c(-5,15), name="Value")
    p2 <- p2 + expand_limits(x = 0, y = 1960) + ylim(1960,2015) + xlim(0,101)

    p2 <- p2 + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    p2 <- p2 + xlab("Age") + ylab("Year") + ggtitle("Female") +theme(plot.title = element_text(lineheight=.8, face="bold"))+theme(plot.title = element_text(hjust = 0.5))
    p2
    
    plotq3 <- grid.arrange(p1, p2, nrow = 1)
    ggsave("plotq3.png", plot = plotq3, width = 12,height = 5)

    
    ###############"""
    
    
    

    fit_M <- list(A.x = A.x, B.x = B.x, K.t = K.t, a.x = alpha.x_M, b.x = beta.x_M, k.t = kappa.t_M, m.tx = m.tx_M, ll = ll_M, npar = npar_M, BIC = BIC_M, residuals = residuals_M)
    fit_F <- list(A.x = A.x, B.x = B.x, K.t = K.t, a.x = alpha.x_F, b.x = beta.x_F, k.t = kappa.t_F, m.tx = m.tx_F, ll = ll_F, npar = npar_F, BIC = BIC_F, residuals = residuals_F)
    
#################################    
#### 	12. PROJECT PARAMETERS
#################################    
#### project for males
    
    #ProjectParameters.SUR.MultiPopulation_BIS=function(fit, nAhead, nSim, CountrySPEC, Methods, ArimaSpec, constraint){
    #fill in inputs
    fit <- fit_M
    nAhead <- 40
    nSim <- 1000
    Methods=c("NR")	
    
    sur.SVD.k.ti	= list("country" = "AR1.0")
    names(sur.SVD.k.ti) = CountrySPEC
    sur.SVD.K.t		= list("RWD")
    sur.NR.k.ti		= list("country" = "AR1.0")
    names(sur.NR.k.ti) = CountrySPEC
    sur.NR.K.t		= list("RWD")
    
    SURSeries = list("SVD" = list("k.ti" = sur.SVD.k.ti, "K.t" = sur.SVD.K.t), "NR" = list("k.ti" = sur.NR.k.ti , "K.t" = sur.NR.K.t))
    ArimaSpec <- SURSeries
    constraint = F
    
      projList_M=fit		# initialize list
      
        
        fit.i=fit	# output method i
        
        
        # step 1. collect k.ti, K.t and their specifications + store in list 'kappaTS' and 'kappaM'
        
        kappaTS = list(K.t = fit.i$K.t, k.t = fit.i$k.t)
        kappaM  = list(K.t = ArimaSpec[[2]]$K.t[[1]], k.t = ArimaSpec[[2]]$k.ti[[1]])
        
        
        # step 2. project kappa TS
        
        #kappaSIM = ProjectSUR.1.model.MultiPop_BIS(TS.list=kappaTS, M.list=kappaM, nSim, nAhead, constraint)
        #ProjectSUR.1.model.MultiPop_BIS=function(TS.list, M.list, nSim, nAhead, constraint){
        
        TS.list=kappaTS
        M.list=kappaM
        
          TS_n = length(TS.list)		# number of kappa time series
          nY   = length(TS.list[[1]])	# number of years in calibration period
          
          ### step 1. get input (string) systemfit
          
          formula.list=list()
          
          CreateFormulaString_BIS=function(TS, Name, index){
            
            if(Name=="RWD"){
              ts_n 	  = length(TS)
              ts_name = deparse(substitute(TS.list))
              
              ts_tailstring = paste(ts_name,"[[", index, "]][", 2, ":", ts_n, "]",sep="")                   
              ts_headstring = paste(ts_name, "[[", index, "]][", 1, ":", ts_n-1, "]",sep="")   
              formulastring = paste(ts_tailstring, "-", ts_headstring, " ~ 1",sep="" )
            }
            
            if(Name=="AR1.0"){
              ts_n 	  = length(TS)
              ts_name = deparse(substitute(TS.list))
              
              ts_tailstring = paste(ts_name,"[[", index, "]][", 2, ":", ts_n, "]",sep="")                   
              ts_headstring = paste(ts_name, "[[", index, "]][", 1, ":", ts_n-1, "]",sep="")   
              formulastring = paste(ts_tailstring, " ~ -1 + ", ts_headstring,sep="")
            }
            
            if(Name=="AR1.1"){
              ts_n 	  = length(TS)
              ts_name = deparse(substitute(TS.list))
              
              ts_tailstring = paste(ts_name,"[[", index, "]][", 2, ":", ts_n, "]",sep="")                   
              ts_headstring = paste(ts_name, "[[", index, "]][", 1, ":", ts_n-1, "]",sep="")   
              formulastring = paste(ts_tailstring, " ~ 1 + ", ts_headstring,sep="")
            }
            
            return(formulastring)
            
          }
          
          
          for(i in 1:TS_n){
            formula 	 = CreateFormulaString_BIS(TS=TS.list[[i]],Name=M.list[[i]], index=i)
            formula.list = append(formula.list, list(as.formula(formula)))
          }
          
          ### step 2. perform SUR regression
          
          # first fit unconstrained problem
          regression_results = systemfit(formula=formula.list, method="SUR", maxiter = 10000, methodResidCov="noDfCor")

          ### step 3. simulate
          Simulate_BIS=function(TS, Name, Coeff, Errors){
            
            if(Name=="RWD"){
              ts_sim = filter(Coeff + Errors, filter=c(1), init=tail(TS,1), method="recursive")
            }
            
            if(Name=="AR1.0"){
              ts_sim = filter(Errors, filter=c(Coeff[1]), init=tail(TS,1), method="recursive")
            }
            
            if(Name=="AR1.1"){
              ts_sim = filter(Coeff[1] + Errors, filter=c(Coeff[2]), init=tail(TS,1), method="recursive")
            }
            
            return(ts_sim)
            
          }
          
          output = TS.list	# initialize list of length TS_n
          
          for(i in 1:TS_n){ 
            if(nSim!=0){
              output[[i]]  	= array(dim=c(nY + nAhead, nSim))
            }
            if(nSim==0){
              output[[i]]  	= array(dim=c(nY + nAhead, 1))
            }
            output[[i]][1:nY,]= TS.list[[i]]
          }
          
          if(nAhead!=0){	
            if(nSim!=0){
              for(i in 1:nSim){
                errors = t( mvrnorm(n=nAhead, mu=rep(0, TS_n), Sigma=regression_results$residCov) )
                
                for(j in 1:TS_n){
                  output[[j]][(nY+1):(nY+nAhead),i] = Simulate_BIS(TS.list[[j]], M.list[[j]], regression_results$eq[[j]]$coefficients, errors[j,])
                }
                
              }
            }
            if(nSim==0){
              errors = matrix(0, nrow=TS_n, ncol=nAhead)
              
              for(j in 1:TS_n){
                output[[j]][(nY+1):(nY+nAhead),1] = Simulate_BIS(TS.list[[j]], M.list[[j]], regression_results$eq[[j]]$coefficients, errors[j,])
              }
              
              
            }
          }
          
          
        kappaSIM <- list(K.t=output[[1]], k.t=output[[2]], drift = unname(regression_results$eq[[1]]$coefficients),
                      a = unname(regression_results$eq[[2]]$coefficients), C = unname(regression_results$residCov))
          
        #}
        
      
        # step 3. get results
        
        K.t.proj 	= kappaSIM$K.t
        k.t.proj 	= list(kappaSIM$k.t)
        
        names(k.t.proj) = "Males"
        
        projList_M[[1]] = list(K.t=K.t.proj, k.t=k.t.proj, drift = kappaSIM$drift, a = kappaSIM$a, C = kappaSIM$C)
 
        
        
#### project for females
        
        #ProjectParameters.SUR.MultiPopulation_BIS=function(fit, nAhead, nSim, CountrySPEC, Methods, ArimaSpec, constraint){
        #fill in inputs
        fit <- fit_F
        
        
        projList_F=fit		# initialize list
        
        
        fit.i=fit	# output method i
        
        
        # step 1. collect k.ti, K.t and their specifications + store in list 'kappaTS' and 'kappaM'
        
        kappaTS = list(K.t = fit.i$K.t, k.t = fit.i$k.t)
        kappaM  = list(K.t = ArimaSpec[[2]]$K.t[[1]], k.t = ArimaSpec[[2]]$k.ti[[1]])
        
        
        # step 2. project kappa TS
        
        #kappaSIM = ProjectSUR.1.model.MultiPop_BIS(TS.list=kappaTS, M.list=kappaM, nSim, nAhead, constraint)
        #ProjectSUR.1.model.MultiPop_BIS=function(TS.list, M.list, nSim, nAhead, constraint){
        
        TS.list=kappaTS
        M.list=kappaM
        
        TS_n = length(TS.list)		# number of kappa time series
        nY   = length(TS.list[[1]])	# number of years in calibration period
        
        ### step 1. get input (string) systemfit
        
        formula.list=list()

        
        
        for(i in 1:TS_n){
          formula 	 = CreateFormulaString_BIS(TS=TS.list[[i]],Name=M.list[[i]], index=i)
          formula.list = append(formula.list, list(as.formula(formula)))
        }
        
        ### step 2. perform SUR regression
        
        # first fit unconstrained problem
        regression_results = systemfit(formula=formula.list, method="SUR", maxiter = 10000, methodResidCov="noDfCor")
        
        ### step 3. simulate
        
        output = TS.list	# initialize list of length TS_n
        
        for(i in 1:TS_n){ 
          if(nSim!=0){
            output[[i]]  	= array(dim=c(nY + nAhead, nSim))
          }
          if(nSim==0){
            output[[i]]  	= array(dim=c(nY + nAhead, 1))
          }
          output[[i]][1:nY,]= TS.list[[i]]
        }
        
        if(nAhead!=0){	
          if(nSim!=0){
            for(i in 1:nSim){
              errors = t( mvrnorm(n=nAhead, mu=rep(0, TS_n), Sigma=regression_results$residCov) )
              
              for(j in 1:TS_n){
                output[[j]][(nY+1):(nY+nAhead),i] = Simulate_BIS(TS.list[[j]], M.list[[j]], regression_results$eq[[j]]$coefficients, errors[j,])
              }
              
            }
          }
          if(nSim==0){
            errors = matrix(0, nrow=TS_n, ncol=nAhead)
            
            for(j in 1:TS_n){
              output[[j]][(nY+1):(nY+nAhead),1] = Simulate_BIS(TS.list[[j]], M.list[[j]], regression_results$eq[[j]]$coefficients, errors[j,])
            }
            
            
          }
        }
        
        
        kappaSIM <- list(K.t=output[[1]], k.t=output[[2]], drift = unname(regression_results$eq[[1]]$coefficients),
                         a = unname(regression_results$eq[[2]]$coefficients), C = unname(regression_results$residCov))
        
        #}
        
        
        # step 3. get results
        
        K.t.proj 	= kappaSIM$K.t
        k.t.proj 	= list(kappaSIM$k.t)
        
        names(k.t.proj) = "Males"
        
        projList_F[[1]] = list(K.t=K.t.proj, k.t=k.t.proj, drift = kappaSIM$drift, a = kappaSIM$a, C = kappaSIM$C)
        
###################################
#### plot projected parameters
###################################       
        # compute percentiles 
        vP = c(0.05, 0.5, 0.95)
        mQuantiles_M=array(NA, dim=c(3, nAhead+1) )
        
        for(k in 1:nAhead){
          mQuantiles_M[,(k+1)]=quantile(projList_M$A.x$K.t[56+k,], probs=vP) 
        }
        mQuantiles_M[,1]= quantile(projList_M$A.x$K.t[56,], probs=c(0.5))

        
        colour="black"
        # make plot
        par(fig=c(0,1,0.5,1))
        plot(yv, projList_M$A.x$K.t[1:56,1] ,ylim=c(-160,40), type="l", col=colour, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=expression(K[t]))
        # polygon 
        x	= (tail(yv,1):(tail(yv,1)+nAhead))
        xx	= c(x,rev(x))
        yy 	= c(mQuantiles_M[1,], rev(mQuantiles_M[3,]) )
        polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
        lines(x=x, y=mQuantiles_M[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
        
        
        # compute percentiles 
        vP = c(0.05, 0.5, 0.95)
        mQuantiles_M=array(NA, dim=c(3, nAhead+1) )
        
        for(k in 1:nAhead){
          mQuantiles_M[,(k+1)]=quantile(projList_M$A.x$k.t$Males[56+k,], probs=vP) 
        }
        mQuantiles_M[,1]= quantile(projList_M$A.x$k.t$Males[56,], probs=c(0.5))
        
        colour="blue"
        # make plot
        par(fig=c(0,0.5,0,0.5),new=T)
        plot(yv, projList_M$A.x$k.t$Males[1:56,1] , type="l", col=colour, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=expression(kappa[t]^M))
        # polygon 
        x	= (tail(yv,1):(tail(yv,1)+nAhead))
        xx	= c(x,rev(x))
        yy 	= c(mQuantiles_M[1,], rev(mQuantiles_M[3,]) )
        polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
        lines(x=x, y=mQuantiles_M[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
   
        
        # compute percentiles 
        vP = c(0.05, 0.5, 0.95)
        mQuantiles_F=array(NA, dim=c(3, nAhead+1) )
        
        for(k in 1:nAhead){
          mQuantiles_F[,(k+1)]=quantile(projList_F$A.x$k.t$Males[56+k,], probs=vP) 
        }
        mQuantiles_F[,1]= quantile(projList_F$A.x$k.t$Males[56,], probs=c(0.5))
        
        colour="red"
        # make plot
        par(fig=c(0.5,1,0,0.5),new=T)
        plot(yv, projList_F$A.x$k.t$Males[1:56,1] , type="l", col=colour, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=expression(kappa[t]^F))
        # polygon 
        x	= (tail(yv,1):(tail(yv,1)+nAhead))
        xx	= c(x,rev(x))
        yy 	= c(mQuantiles_F[1,], rev(mQuantiles_F[3,]) )
        polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
        lines(x=x, y=mQuantiles_F[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile

        
        
#### 	13. PROJECT MORTALITY q(t,x)
# make projections  males      
        fitList=fit_M 
        projList=projList_M
        CountrySPEC=c("Males")
        Methods=Methods       

          T 	= length(yv)		# number of years in calibration period
          nA	= length(xv)		# number of ages
          q.xt	= fitList			# initialize list of same length fitList
          

            
          fit.j	 = fitList
          proj.j = projList[[1]]
            
          # country specific mortality
          m.xti	 = array(NA,dim=c(nA,(T+nAhead),nSim))	# initialize
          m.xti[,1:T,]= t(fit.j$m.tx)				# fitted mortality
          
          for(t in (T+1):(T+nAhead)){ 	# projected death rates: use most recent death rate as starting point
            m.xti[,t,]=m.xti[,T,]*exp(fit.j$B.x %o% (proj.j$K.t[t,] - fit.j$K.t[T]) + fit.j$b.x %o% (proj.j$k.t[[1]][t,] -fit.j$k.t[T]) )
          }
          
          q.xti = list(1-exp(-m.xti))		# projected mortality rates
          
          names(q.xti) = CountrySPEC
          
          
          q.xt[[1]] = q.xti
          q.xti_M = q.xti
          
          # make projections  females      
          fitList=fit_F 
          projList=projList_F
          CountrySPEC=c("Females")
          Methods=Methods       
          
          T 	= length(yv)		# number of years in calibration period
          nA	= length(xv)		# number of ages
          q.xt_F	= fitList			# initialize list of same length fitList
          
          
          
          fit.j	 = fitList
          proj.j = projList[[1]]
          
          # country specific mortality
          m.xti	 = array(NA,dim=c(nA,(T+nAhead),nSim))	# initialize
          m.xti[,1:T,]= t(fit.j$m.tx)				# fitted mortality
          
          for(t in (T+1):(T+nAhead)){ 	# projected death rates: use most recent death rate as starting point
            m.xti[,t,]=m.xti[,T,]*exp(fit.j$B.x %o% (proj.j$K.t[t,] - fit.j$K.t[T]) + fit.j$b.x %o% (proj.j$k.t[[1]][t,] -fit.j$k.t[T]) )
          }
          
          q.xti = list(1-exp(-m.xti))		# projected mortality rates
          
          names(q.xti) = CountrySPEC
          
          
          q.xt_F[[1]] = q.xti
            
          
          
###plot projected death rates for male of 25,45,65,85 years old
          age <- c(26,46,66,86)
          par(mfrow=c(4,2))
          
          for(i in 1:4){
            colour="blue"
            obs_death_rate = matrix(ncol=1,nrow=56)
            obs_death_rate[,1] <- dtx_M[,age[i]]/etx_M[,age[i]]
            #1-exp(-mu)
            obs_qxt <- 1-exp(-obs_death_rate)
            fit_qxt <- q.xti_M$Males[age[i],1:56,2]
            qQuantiles <- q.xti_M$Males[age[i],57:96,]
            
            # compute percentiles 
            vP = c(0.05, 0.5, 0.95)
            mQuantiles_M=array(NA, dim=c(3, nAhead+1) )
            
            for(k in 1:nAhead){
              mQuantiles_M[,(k+1)]=quantile(qQuantiles[k,], probs=vP) 
            }
            mQuantiles_M[,1]= fit_qxt[56]

            
    
            bounds = c(min(mQuantiles_M),max(obs_qxt))
            plot(yv, t(obs_qxt) ,ylim=bounds, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=paste("Males: Age ",age[i]-1))
            lines(yv,t(fit_qxt),col="blue",lty=2)
            # polygon 
            x	= (tail(yv,1):(tail(yv,1)+nAhead))
            xx	= c(x,rev(x))
            yy 	= c(mQuantiles_M[1,], rev(mQuantiles_M[3,]) )
            polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
            lines(x=x, y=mQuantiles_M[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
          
            colour="red"
            obs_death_rate = matrix(ncol=1,nrow=56)
            obs_death_rate[,1] <- dtx_F[,age[i]]/etx_F[,age[i]]
            #1-exp(-mu)
            obs_qxt <- 1-exp(-obs_death_rate)
            fit_qxt <- q.xti$Females[age[i],1:56,2]
            qQuantiles <- q.xti$Females[age[i],57:96,]
            
            # compute percentiles 
            vP = c(0.05, 0.5, 0.95)
            mQuantiles_M=array(NA, dim=c(3, nAhead+1) )
            
            for(k in 1:nAhead){
              mQuantiles_M[,(k+1)]=quantile(qQuantiles[k,], probs=vP) 
            }
            mQuantiles_M[,1]= fit_qxt[56]
            
            
            
            bounds = c(min(mQuantiles_M),max(obs_qxt))
            plot(yv, t(obs_qxt) ,ylim=bounds, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=paste("Females: Age ",age[i]-1))
            lines(yv,t(fit_qxt),col="red",lty=2)
            # polygon 
            x	= (tail(yv,1):(tail(yv,1)+nAhead))
            xx	= c(x,rev(x))
            yy 	= c(mQuantiles_M[1,], rev(mQuantiles_M[3,]) )
            polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
            lines(x=x, y=mQuantiles_M[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
            
            
            }
          