################################################
### Assignment 2
################################################
################################################
### Question 1
################################################

# set working directory
setwd("/Users/tineh/Google drive/universiteit/Master of financial and actuarial engineering/ALIM/Assignment 2")

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

### Read in data
MLT_BE = read.table("mltper_1x1_BE.txt",header=TRUE,dec=".",sep="",skip=2)
FLT_BE = read.table("fltper_1x1_BE.txt",header=TRUE,dec=".",sep="",skip=2)

MMat_BE_M = matrix(data=0, ncol=nyear, nrow=nage, dimnames=list(c(cname),c(rname)) )
MMat_BE_F = matrix(data=0, ncol=nyear, nrow=nage, dimnames=list(c(cname),c(rname)) )

xMmatrix = MLT_BE[,-4:-10]
xFmatrix = FLT_BE[,-4:-10]
for(i in 1:nyear){
  for(j in 1:nage){
    MMat_BE_M[j,i] = as.numeric(as.character(subset(xMmatrix, Year==yv[i] & Age==xv[j])$mx))
    MMat_BE_F[j,i] = as.numeric(as.character(subset(xFmatrix, Year==yv[i] & Age==xv[j])$mx))
  }
}

# Avoid problem of taking logarithm of zero
MMat_BE_M[MMat_BE_M==0]=0.01
MMat_BE_F[MMat_BE_F==0]=0.01


### figuur vergelijking central death rate men and women in Belgium
yearrange <- c(1,11,21,31,41,51)
value <- c(MMat_BE_M[1:96,1])
for(i in 2:6){
  value <- c(value,MMat_BE_M[1:96,yearrange[i]])
}
for(i in 1:6){
  value <- c(value,MMat_BE_F[1:96,yearrange[i]])
}


Sex_values <- c(rep(Sex[1],6*96),rep(Sex[2],6*96))
x <- c(rep(c(0:95),2*6))
Year <- rep(1960,96)
yearvalues <- c(1970,1980,1990,2000,2010)
for(i in 1:5){
  Year <- c(Year, rep(yearvalues[i],96))
}
Year <- c(Year,Year)
value <- log(value)

df <- data.frame(x, Sex_values, value,Year)


d <- ggplot(df, aes(x=x, y=value, group=interaction(Year,Sex_values), colour=Sex_values,alpha=Year)) +
  geom_line(size=0.5) +
  background_grid(major = "xy", minor = "none") +
  labs(x="Years",y="Death rate")
grid.arrange(arrangeGrob(d))



################################################
### Question 2
################################################

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

ExpMat_BE_M[DeathMat_BE_M==0]=0.01
ExpMat_BE_F[DeathMat_BE_F==0]=0.01

dtx_M <- DeathMat_BE_M
etx_M <- ExpMat_BE_M
dtx_F <- DeathMat_BE_F
etx_F <- ExpMat_BE_F

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

#### 1. FIT LIFEMETRICS MODELS
### LEE-CARTER (Model M1)	 	
##		log m(t,x) = beta1(x) + beta2(x)*kappa2(t) + Poisson error

# fill in parameters function 
int = xv*0  
constraints = "ALL"  
exclude.coh = FALSE

  mtx_M	= dtx_M/etx_M      			# matrix of death rates
  mtx_F = dtx_F/etx_F
  qtx_M	= 1-exp(-mtx_M)   			# matrix of mortality rates
  qtx_F	= 1-exp(-mtx_F) 
  
  n = length(xv)				# number of ages
  m = length(yv)				# number of years
  
  cy= (yv[1]-xv[n]):(yv[m]-xv[1])  	# cohort approximate years of birth
  
  # Initialise parameter vectors
  beta1v_M=int
  beta2v_M=(1:n)*0
  beta3v_M=(1:n)*0			# dummy vector, this will stay at 0
  kappa2v_M=(1:m)*0
  gamma3v_M=(1:(n+m-1))*0		# dummy vector, this will stay at 0
  beta1v_F=int
  beta2v_F=(1:n)*0
  beta3v_F=(1:n)*0			# dummy vector, this will stay at 0
  kappa2v_F=(1:m)*0
  gamma3v_F=(1:(n+m-1))*0		# dummy vector, this will stay at 0
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
      beta1v_M[j]=sum(log(mtx_M[,j])*wa[,j])/sum(wa[,j])
      beta1v_F[j]=sum(log(mtx_F[,j])*wa[,j])/sum(wa[,j])      
    }
    beta2v_M[j]=1/n
    beta2v_F[j]=1/n
  }
  kappa2v_M=(m:1)-(m+1)/2
  kappa2v_F=(m:1)-(m+1)/2
  
  # Stage 1: iterate
  l0=-1000000
  l1=-999999
  iteration=0
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001

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
  
  
    while(abs(l1-l0) > 1e-20)
  {
    iteration=iteration+1
    
    l0=l1
    
    #     mhat_M=mtx*0
    #     for(i in 1:m)
    #     {
    #       mhat_M[i,]=exp(beta1v_M+beta2v_M*kappa2v_M[i]+beta3v_M*gamma3v_M[(n+i-1):i])
    #     }
    #     epsilon_M=(dtx_M-etx_M*mhat_M)/sqrt(etx_M*mhat_M)
    #     l1=sum((dtx_M*log(etx_M*mhat_M)-etx_M*mhat_M-lgamma(dtx_M+1))*wa)
    #     
    #    if(iteration==1) 
    # Stage 1B optimise over the beta2(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx_M[,j]	# actual deaths
      ev=etx_M[,j]	# exposure
      beta2v_M[j]=llmaxM2B(beta1v_M[j],beta2v_M[j],beta3v_M[j],
                         kappa2v_M,gamma3v_M[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat_M=mtx_M*0
    for(i in 1:m)
    {
      mhat_M[i,]=exp(beta1v_M+beta2v_M*kappa2v_M[i]+beta3v_M*gamma3v_M[(n+i-1):i])
    }
    epsilon_M=(dtx_M-etx_M*mhat_M)/sqrt(etx_M*mhat_M)
    l1=sum((dtx_M*log(etx_M*mhat_M)-etx_M*mhat_M-lgamma(dtx_M+1))*wa)
    #cat(l1,"->")
    
    # Stage 1D optimise over the kappa2(t)
    for(i in 1:m)
    {		 		 
      # cycle through the range of years
      dv=dtx_M[i,]	# actual deaths
      ev=etx_M[i,]	# exposure
      kappa2v_M[i]=llmaxM2D(beta1v_M,beta2v_M,beta3v_M,
                          kappa2v_M[i],gamma3v_M[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat_M=mtx_M*0
    for(i in 1:m)
    {
      mhat_M[i,]=exp(beta1v_M+beta2v_M*kappa2v_M[i]+beta3v_M*gamma3v_M[(n+i-1):i])
    }
    epsilon_M=(dtx_M-etx_M*mhat_M)/sqrt(etx_M*mhat_M)
    l1=sum((dtx_M*log(etx_M*mhat_M)-etx_M*mhat_M-lgamma(dtx_M+1))*wa)
    #cat(l1,"->")	
    
    # Now apply the constraints
    if(constraints=="ALL"){
      fac21=mean(kappa2v_M)
      fac22=sum(beta2v_M)
      kappa2v_M=fac22*(kappa2v_M-fac21)    # sum kappa2=0
      beta2v_M=beta2v_M/fac22             # sum beta2=1
      beta1v_M=beta1v_M+beta2v_M*fac22*fac21 # => adjust beta1 
    } else{
      fac22=sum(beta2v_M)
      kappa2v_M=fac22*kappa2v_M           # => adjust kappa2
      beta2v_M=beta2v_M/fac22  		  # sum beta2=1
    }
    
    mhat_M=mtx_M*0
    for(i in 1:m)
    {
      mhat_M[i,]=exp(beta1v_M+beta2v_M*kappa2v_M[i]+beta3v_M*gamma3v_M[(n+i-1):i])
    }
    epsilon_M=(dtx_M-etx_M*mhat_M)/sqrt(etx_M*mhat_M)
    l1=sum((dtx_M*log(etx_M*mhat_M)-etx_M*mhat_M-lgamma(dtx_M+1))*wa)
    #cat(l1,"->")
    
    
    # Stage 1A optimise over the beta1(x)
    if(all(int==0)){
      for(j in 1:n)
      {		 		 
        # cycle through the range of years
        wv=1	    # can be set to a vector of weights
        # to e.g. exclude duff years
        wv=wa[,j]
        s1=sum(wv*dtx_M[,j])
        s2=sum(wv*etx_M[,j]*exp(beta2v_M[j]*kappa2v_M+beta3v_M[j]*gamma3v_M[(n+1-j):(n+m-j)]))
        beta1v_M[j]=log(s1)-log(s2)
      }
    }
    
    mhat_M=mtx_M*0
    for(i in 1:m)
    {
      mhat_M[i,]=exp(beta1v_M+beta2v_M*kappa2v_M[i]+beta3v_M*gamma3v_M[(n+i-1):i])
    }
    epsilon_M=(dtx_M-etx_M*mhat_M)/sqrt(etx_M*mhat_M)
    l1=sum((dtx_M*log(etx_M*mhat_M)-etx_M*mhat_M-lgamma(dtx_M+1))*wa)
    #cat(l1,"->")
    
    
  }		 # end while loop
  
  # calculate number of parameters and deduct the number of constraints
  # also count number of parameters in 'beta1v_M'!
  if(constraints=="ALL"){
    npar_M=length(beta1v_M)+length(beta2v_M)+length(kappa2v_M)-2
  } else{
    npar_M=length(beta1v_M)+length(beta2v_M)+length(kappa2v_M)-1
  }
  
  # Calculate the BIC
  BIC_M=-2*l1+log(sum(wa))*npar_M
  
  l0=-1000000
  l1=-999999
  iteration=0
  
  while(abs(l1-l0) > 1e-20)
  {
    iteration=iteration+1
    
    l0=l1
    
    #     mhat_M=mtx*0
    #     for(i in 1:m)
    #     {
    #       mhat_M[i,]=exp(beta1v_M+beta2v_M*kappa2v_M[i]+beta3v_M*gamma3v_M[(n+i-1):i])
    #     }
    #     epsilon_M=(dtx_M-etx_M*mhat_M)/sqrt(etx_M*mhat_M)
    #     l1=sum((dtx_M*log(etx_M*mhat_M)-etx_M*mhat_M-lgamma(dtx_M+1))*wa)
    #     
    #    if(iteration==1) 
    # Stage 1B optimise over the beta2(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx_F[,j]	# actual deaths
      ev=etx_F[,j]	# exposure
      beta2v_F[j]=llmaxM2B(beta1v_F[j],beta2v_F[j],beta3v_F[j],
                           kappa2v_F,gamma3v_F[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat_F=mtx_F*0
    for(i in 1:m)
    {
      mhat_F[i,]=exp(beta1v_F+beta2v_F*kappa2v_F[i]+beta3v_F*gamma3v_F[(n+i-1):i])
    }
    epsilon_F=(dtx_F-etx_F*mhat_F)/sqrt(etx_F*mhat_F)
    l1=sum((dtx_F*log(etx_F*mhat_F)-etx_F*mhat_F-lgamma(dtx_F+1))*wa)
    #cat(l1,"->")
    
    # Stage 1D optimise over the kappa2(t)
    for(i in 1:m)
    {		 		 
      # cycle through the range of years
      dv=dtx_F[i,]	# actual deaths
      ev=etx_F[i,]	# exposure
      kappa2v_F[i]=llmaxM2D(beta1v_F,beta2v_F,beta3v_F,
                            kappa2v_F[i],gamma3v_F[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat_F=mtx_F*0
    for(i in 1:m)
    {
      mhat_F[i,]=exp(beta1v_F+beta2v_F*kappa2v_F[i]+beta3v_F*gamma3v_F[(n+i-1):i])
    }
    epsilon_F=(dtx_F-etx_F*mhat_F)/sqrt(etx_F*mhat_F)
    l1=sum((dtx_F*log(etx_F*mhat_F)-etx_F*mhat_F-lgamma(dtx_F+1))*wa)
    #cat(l1,"->")	
    
    # Now apply the constraints
    if(constraints=="ALL"){
      fac21=mean(kappa2v_F)
      fac22=sum(beta2v_F)
      kappa2v_F=fac22*(kappa2v_F-fac21)    # sum kappa2=0
      beta2v_F=beta2v_F/fac22             # sum beta2=1
      beta1v_F=beta1v_F+beta2v_F*fac22*fac21 # => adjust beta1 
    } else{
      fac22=sum(beta2v_F)
      kappa2v_F=fac22*kappa2v_F           # => adjust kappa2
      beta2v_F=beta2v_F/fac22  		  # sum beta2=1
    }
    
    mhat_F=mtx_F*0
    for(i in 1:m)
    {
      mhat_F[i,]=exp(beta1v_F+beta2v_F*kappa2v_F[i]+beta3v_F*gamma3v_F[(n+i-1):i])
    }
    epsilon_F=(dtx_F-etx_F*mhat_F)/sqrt(etx_F*mhat_F)
    l1=sum((dtx_F*log(etx_F*mhat_F)-etx_F*mhat_F-lgamma(dtx_F+1))*wa)
    #cat(l1,"->")
    
    
    # Stage 1A optimise over the beta1(x)
    if(all(int==0)){
      for(j in 1:n)
      {		 		 
        # cycle through the range of years
        wv=1	    # can be set to a vector of weights
        # to e.g. exclude duff years
        wv=wa[,j]
        s1=sum(wv*dtx_F[,j])
        s2=sum(wv*etx_F[,j]*exp(beta2v_F[j]*kappa2v_F+beta3v_F[j]*gamma3v_F[(n+1-j):(n+m-j)]))
        beta1v_F[j]=log(s1)-log(s2)
      }
    }
    
    mhat_F=mtx_F*0
    for(i in 1:m)
    {
      mhat_F[i,]=exp(beta1v_F+beta2v_F*kappa2v_F[i]+beta3v_F*gamma3v_F[(n+i-1):i])
    }
    epsilon_F=(dtx_F-etx_F*mhat_F)/sqrt(etx_F*mhat_F)
    l1=sum((dtx_F*log(etx_F*mhat_F)-etx_F*mhat_F-lgamma(dtx_F+1))*wa)
    #cat(l1,"->")
    
    
  }		 # end while loop
  
  # calculate number of parameters and deduct the number of constraints
  # also count number of parameters in 'beta1v_M'!
  if(constraints=="ALL"){
    npar_F=length(beta1v_F)+length(beta2v_F)+length(kappa2v_F)-2
  } else{
    npar_F=length(beta1v_F)+length(beta2v_F)+length(kappa2v_F)-1
  }
  
  # Calculate the BIC
  BIC_F=-2*l1+log(sum(wa))*npar_F
  
  
  
par(mfrow=c(2,3))
plot(beta1v_M,type="l",xlab="Age",ylab="BE",main=expression(A[x]^M))
plot(beta2v_M,type="l",xlab="Age",ylab="BE",main=expression(B[x]^M))
plot(yv, kappa2v_M,type="l",xlab="Year",ylab="BE",main=expression(K[t]^M))
plot(beta1v_F,type="l",xlab="Age",ylab="BE",main=expression(A[x]^F))
plot(beta2v_F,type="l",xlab="Age",ylab="BE",main=expression(B[x]^F))
plot(yv, kappa2v_F,type="l",xlab="Year",ylab="BE",main=expression(K[t]^F))


######## Vanaf hier nog man en vrouw plot maken naast elkaar
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
image(t(epsilon_M),xlab="Age",ylab="Year",col=r,main="Male",axes=F)
mtext(text=yv_fig, side=2, line=0.3, at=seq(0,1,5/length(yv)), cex=0.8,las=1)
mtext(text=xv_fig, side=1, line=0.3, at=seq(0,1,5/length(xv)), las=2, cex=0.8)
par(fig=c(0.1,0.5,0,1),new=T)
image.plot(t(epsilon_M),legend.only=T)
par(fig=c(0.52,0.95,0,1),new=T)
image(t(epsilon_F),xlab="Age",ylab="Year",col=r,main="Female",axes=F)
mtext(text=yv_fig, side=2, line=0.3, at=seq(0,1,5/length(yv)), cex=0.8,las=1)
mtext(text=xv_fig, side=1, line=0.3, at=seq(0,1,5/length(xv)), las=2, cex=0.8)
par(fig=c(0.6,1,0,1),new=T)
image.plot(t(epsilon_M),legend.only=T)


 #######################""
    # Male 
    dataplotM <- melt(epsilon_M)
    colnames(dataplotM) <- c("Year","Age","Residual")
    
    cols <- rev(rainbow(20)[-20])
    hex <- c("#FF0000", "#FFA500", "#FFFF00", "#008000", "#9999ff", "#000066")

    
    p1 <- ggplot(data = dataplotM, aes(x = dataplotM$Age, y = dataplotM$Year))  + geom_tile(aes(fill = dataplotM$Residual)) + scale_fill_gradient2() 
    p1 <- p1  + scale_fill_gradientn(colours = r, limits=c(-5,15), name="Value")
    p1 <- p1 + expand_limits(x = 0, y = 1960) + ylim(1960,2015) + xlim(0,101)

    p1 <- p1 + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    p1 <- p1 + xlab("Age") + ylab("Year") + ggtitle("Male")+theme(plot.title = element_text(lineheight=.8, face="bold"))+theme(plot.title = element_text(hjust = 0.5))
    p1
    
    
    # Female
    dataplotF <- melt(epsilon_F)
    colnames(dataplotF) <- c("Year","Age","Residual")
    
    hex <- c("#FF0000", "#FFA500", "#FFFF00", "#008000", "#9999ff", "#000066")

    
    p2 <- ggplot(data = dataplotF, aes(x = dataplotF$Age, y = dataplotF$Year))  + geom_tile(aes(fill = dataplotF$Residual)) + scale_fill_gradient2() 
    p2 <- p2  + scale_fill_gradientn(colours = r, limits=c(-5,15), name="Value")
    p2 <- p2 + expand_limits(x = 0, y = 1960) + ylim(1960,2015) + xlim(0,101)

    p2 <- p2 + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    p2 <- p2 + xlab("Age") + ylab("Year") + ggtitle("Female") +theme(plot.title = element_text(lineheight=.8, face="bold"))+theme(plot.title = element_text(hjust = 0.5))
    p2
    
    plotq2 <- grid.arrange(p1, p2, nrow = 1)
    ggsave("plotq2.png", plot = plotq2, width = 12,height = 5)


###project the parameters
x=kappa2v_M# input the time dependent parameter

order=c(0,1,0)# the Arima order (example: Random Walk)

library(forecast)
fit<-Arima(x,order=order,include.drift=TRUE)
summary(fit);
mu=coef(fit)# the drift parameter


#Project
f.years=40
project.x<-forecast(fit,h=f.years,level=c(80,95))
sim.x_M<-replicate(1000,simulate(fit,nsim=f.years)); 

t=1:length(x)
years<-1:(length(x)+f.years)
limit.y_M<-range(c(range(sim.x_M),range(x)))
limit.x<-range(years)

x=kappa2v_F# input the time dependent parameter

fit<-Arima(x,order=order,include.drift=TRUE)
summary(fit);
mu=coef(fit)# the drift parameter


#Project
f.years=40
project.x<-forecast(fit,h=f.years,level=c(80,95))
sim.x_F<-replicate(1000,simulate(fit,nsim=f.years)); 

t=1:length(x)
years<-1:(length(x)+f.years)
limit.y_F<-range(c(range(sim.x_F),range(x)))
limit.x<-range(years)

#plot the projected parameters
nSim 	 = 1000			# number of simulations
nAhead = 40	# number of years to project in future

# compute percentiles 
vP = c(0.05, 0.5, 0.95)
mQuantiles_M=array(NA, dim=c(3, nAhead+1) )
mQuantiles_F=array(NA, dim=c(3, nAhead+1) )

for(k in 1:nAhead){
  mQuantiles_M[,(k+1)]=quantile(sim.x_M[k,], probs=vP) 
  mQuantiles_F[,(k+1)]=quantile(sim.x_F[k,], probs=vP) 
}
mQuantiles_M[,1]= kappa2v_M[length(kappa2v_M)]
mQuantiles_F[,1]= kappa2v_F[length(kappa2v_F)]

bounds=c(-160,30)
colour="blue"
# make plot
par(fig=c(0,1,0,1))
plot(yv, kappa2v_M , type="l", ylim= bounds, col=colour, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=expression(K[t]))
# polygon 
x	= (tail(yv,1):(tail(yv,1)+nAhead))
xx	= c(x,rev(x))
yy 	= c(mQuantiles_M[1,], rev(mQuantiles_M[3,]) )
polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
lines(x=x, y=mQuantiles_M[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
colour="red"
lines(yv, kappa2v_F , type="l", ylim= bounds, col=colour, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="")
x	= (tail(yv,1):(tail(yv,1)+nAhead))
xx	= c(x,rev(x))
yy 	= c(mQuantiles_F[1,], rev(mQuantiles_F[3,]) )
polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
lines(x=x, y=mQuantiles_F[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile

###plot projected death rates for male of 25,45,65,85 years old
age <- c(26,46,66,86)
par(mfrow=c(4,2))

for(i in 1:4){
  colour="blue"
  obs_death_rate = matrix(ncol=1,nrow=56)
  obs_death_rate[,1] <- dtx_M[,age[i]]/etx_M[,age[i]]
  #1-exp(-mu)
  obs_qxt <- 1-exp(-obs_death_rate)
  fit_mu <- matrix(ncol = 1, nrow=56)
  fit_mu[,1] <- exp(beta1v_M[age[i]] + beta2v_M[age[i]]*kappa2v_M)
  fit_qxt <- 1-exp(-fit_mu)
  
  muQuantiles <- exp(beta1v_M[age[i]] + beta2v_M[age[i]]*mQuantiles_M)
  qQuantiles <- 1-exp(-muQuantiles)
  bounds = c(min(qQuantiles),max(obs_qxt))
  plot(yv, t(obs_qxt) ,ylim=bounds, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=paste("Males: Age ",age[i]-1))
  lines(yv,t(fit_qxt),col="blue",lty=2)
  # polygon 
  x	= (tail(yv,1):(tail(yv,1)+nAhead))
  xx	= c(x,rev(x))
  yy 	= c(qQuantiles[1,], rev(qQuantiles[3,]) )
  polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
  lines(x=x, y=qQuantiles[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
  colour="red"
  obs_death_rate = matrix(ncol=1,nrow=56)
  obs_death_rate[,1] <- dtx_F[,age[i]]/etx_F[,age[i]]
  #1-exp(-mu)
  obs_qxt <- 1-exp(-obs_death_rate)
  fit_mu <- matrix(ncol = 1, nrow=56)
  fit_mu[,1] <- exp(beta1v_F[age[i]] + beta2v_F[age[i]]*kappa2v_F)
  fit_qxt <- 1-exp(-fit_mu)
  
  muQuantiles <- exp(beta1v_F[age[i]] + beta2v_F[age[i]]*mQuantiles_F)
  qQuantiles <- 1-exp(-muQuantiles)
  bounds = c(min(qQuantiles),max(obs_qxt))
  plot(yv, t(obs_qxt) ,ylim=bounds, xlim=c(yv[1],tail(yv,1)+nAhead), xlab = "Year", ylab="",main=paste("Females: Age ",age[i]-1))
  lines(yv,t(fit_qxt),col="red",lty=2)
  # polygon 
  x	= (tail(yv,1):(tail(yv,1)+nAhead))
  xx	= c(x,rev(x))
  yy 	= c(qQuantiles[1,], rev(qQuantiles[3,]) )
  polygon(xx, yy, col=rgb(t(col2rgb(colour)/255), alpha=0.3), border=rgb(t(col2rgb(colour)/255), alpha=1), lwd=1)	# plot 5%, 95% percentile
  lines(x=x, y=qQuantiles[2,], col=rgb(t(col2rgb(colour)/255), alpha=1), lty=1, lwd=1)					# plot 50% percentile
}

