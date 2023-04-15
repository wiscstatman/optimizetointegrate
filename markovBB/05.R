

# a synthetic autoregressive example, but with goofy bimodal conditionals

set.seed(312345126)

rbimod <- function(mu, delta=2, sig=(1/4), phi=.95 )
 {
   # make a biomodal draw centered at mu
   sgn <- sample( c(-delta,delta), size=1, prob=c(1/2,1/2) ) 
   x <- phi*mu + sig*(rnorm(1)+sgn)*sqrt(abs(mu))
   x
 }

n <- 500
z <- numeric(n)
z[1] <- rnorm(1)
for( i in 2:n){ z[i] <- rbimod(z[i-1]) }





x <- z ## don't show, but match paper notation
B <- n
y <- x[2:B] - x[1:(B-1)] ## the deviations


x.from <- x[1:(B-1)]
x.to <- x[2:B]



plot(x, ylab="z", type="l",  col="blue", las=1)






library(hexbin)
h <- hexbin( cbind(x.from, x.to), xbins=50, xlab="z[from]", ylab="z[to]" )
plot(h)




# reordered by response level
oo <- order(x.from)
x.from.o <- x.from[oo]
x.to.o <- x.to[oo]

tau <- 1/1000
env <- outer(x.from.o,x.from.o,function(u,v){ exp( -(1/tau)*(u-v)^2 ) } )


# Now we realign to x.to2

oo2 <- order( x.to.o )
x.to.o2 <- x.to.o[oo2] ## sorted x.to's

env2 <- env[,oo2]


## one case
w <- env2[200,]/sum( env2[200,] )
k200 <- density( x.to.o2, weights=w )

image( sqrt(env2), xlab="from quantile", ylab="to quantile", main="final urns" )


## show how k200, for example, approximates the generative kernel at that x.from position...

## one case
w <- env2[200,]/sum( env2[200,] )
k200 <- density( x.to.o2, weights=w ) 
plot(k200,main="normalized measure at z[from] = -0.27", xlab="z[to]")


plot( x.to.o2, w, type="h", xlab="z[to]", main="discrete view of normalized measure at z[from]=-0.27" ,ylab="normalized envelope", las=1 )



zgrid <- seq( min(x), max(x), length=200 )
alpha <- zgrid
for( j in 1:length(zgrid) )
 {
   alpha[j] <-  sum( exp( -(1/tau)*(zgrid[j]-x.from)^2 )  ) 
 }
plot( zgrid, alpha, type="l", lwd=2, col="blue", las=1 , xlab="z[from]", ylab="total urn mass" )


# generate a weight vector for each 


## simplest parameter is correlation
#
# E[  (Z_i-mu) * (Z_{i-1}-mu ) ]
# = E[  (Z_{i-1}-mu) E( (Z_i-mu) | Z_{i-1} ) ] 
# = \sum_{z_{from}} pi_{z_from} ( z_from - mu ) \sum_{z_to} 
# 	w_{z.from, z.to} (z_to - mu)


## compute correlation (and other features of kernel)


nsim <- 100
rho <- numeric(nsim)
mu <- rho; vv <- rho
for( isim in 1:nsim )
 {
  W <- matrix(NA,n-1,n-1)
  # transition matrix
  for( i in 1:(n-1) ){ W[i,] <- rgamma( n-1, shape=env2[i,] ) }
  W <- W/rowSums(W)
  # power method to get stationary dist
  m <- 20
  A <- W
  for( i in 1:m ){ A <- A %*% W }
  pi <- A[1,]  ## or any row
  mu[isim] <- sum( pi*x.to.o2 )
  cmeans <- W %*% x.to.o2
  vv[isim] <- sum( pi*(x.to.o2-mu[isim] )^2  )
  rho[isim] <- sum( pi*(x.to.o2-mu[isim])*(cmeans-mu[isim]) )/vv[isim]
  print(isim)
 }
