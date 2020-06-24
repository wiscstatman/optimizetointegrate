# Normal Means
# Setting: 1 parameter "theta" and 1 observation "y",
# y = theta + N(0,1)
# theta has a lasso prior

# exact WBB solution
WBB = function(w1,w2,y,lambda) (y - lambda*w2/w1*sign(y))*(lambda*w2/w1<abs(y))

# calculate "posterior mean" by WBB
post.mean = function(y,N,lambda)
{
  w1 = rexp(N)
  w2 = rexp(N)
  bootstrap = NULL 
  
  for(i in 1:N)
  {
    bootstrap[i] = WBB(w1[i],w2[i],y,lambda)
  }
  
  return(mean(bootstrap))
}

# plot posterior mean (WBB) vs. y. The actual posterior mean is 
lambda = c(0.5,1,2)
y = seq(-5,5,0.1)
post = matrix(0,nrow=length(lambda),ncol = length(y))
actual.mode = matrix(0,nrow=length(lambda),ncol = length(y))
actual.mean = matrix(0,nrow=length(lambda),ncol = length(y))

N = 10000 # bootstrap size
for(i in 1:length(lambda))
{
  for(j in 1:length(y)) 
  {
    post[i,j] = post.mean(y[j],N,lambda[i])
    actual.mode[i,j] = WBB(1,1,y[j],lambda[i])
    
    n = integrate(function(x) x*exp(-(y[j]-x)^2/2-lambda[i]*abs(x)), lower = -Inf,upper = Inf)$value
    d = integrate(function(x) exp(-(y[j]-x)^2/2-lambda[i]*abs(x)), lower = -Inf,upper = Inf)$value
    actual.mean[i,j] = n/d
  }
    
}

pdf( file="1-simple.pdf" )
plot(y,y,"n",ylab="posterior mean", cex.axis=1.3 , cex.lab=1.3)
for(i in 1:length(lambda))
{
  lines(y,post[i,],col=i+1,lwd=2)
  #lines(y,actual.mode[i,],col=i+1,lwd=2,lty=3)
  lines(y,actual.mean[i,],col=i+1,lwd=2,lty=2)
}
abline(h=0)
abline(v=0)
#legend("topleft",legend=c("lambda = 0.5","lambda = 1","lambda = 2"),
#       lwd=2,col=c(1:length(lambda)+1), cex=1.3)
legend("topleft",legend= c( expression(lambda== 0.5), expression(lambda==1), expression(lambda==2) ),
       lwd=2,col=c(1:length(lambda)+1), cex=1.3)

dev.off()
