library(genlasso)
library(quadprog)
library(quantmod)
source("trendfilter3.R")

#  data 
set.seed(0)
n = 500
beta0 = sin(1:n*4*pi/n)*exp(seq(0,3,length.out = n))
y = beta0 + rnorm(n,sd=1)
plot(1:n,y)


# trend filter
n = length(y)
T = 500 # bootstrap size
Beta = matrix(NA,ncol = T,nrow = n)
Or = 3
a = trendfilter(y,ord=Or)
Lambda = 5*10^3
for(t in 1:T)
{
  cat(t,"\n")
  Beta[,t] = as.numeric(trendfilter3(y,w=rexp(n+1),lambda=Lambda, ord=Or))
}

# standard error
SD = apply(Beta,1,sd)


# plot

#pdf( file="1-filtering.pdf", height=5, width=7 )
pdf( file="1-filtering.pdf" , height=6, width=8 )
b = coef.genlasso(a,lambda = Lambda)$beta # trend estimate
#plot(1:length(y),y,xlab="position",main="Trend filtering (cubic)",pch=21,bg="gray",cex=1.5)
plot(1:length(y),y,xlab="position",main="",pch=21,bg="gray",cex=1.2, cex.axis=1.4,cex.lab=1.4, las=1)
#lines(1:length(y),b,lwd=4,col="red")
#lines(1:length(y),b+2*SD,lwd=4,col="cornflowerblue")
#lines(1:length(y),b-2*SD,lwd=4,col="cornflowerblue")
lines(1:length(y),b,lwd=2,col="red")
lines(1:length(y),b+2*SD,lwd=2,col="cornflowerblue")
lines(1:length(y),b-2*SD,lwd=2,col="cornflowerblue")

dev.off()
