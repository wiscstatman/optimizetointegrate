
# March 24/23

# 01 and 02 are offbase on how to the weights
# it's better if a weight associated with distance affects the next point that
# the chain reaches, not the current point

# a synthetic example; generated as AR, but with goofy bimodal conditionals

rm( list=ls() )

set.seed(312345126)


rbimod <- function(mu, delta=2, sig=(1/4), phi=.95 )
 {
   # make a biomodal draw centered at mu
   sgn <- sample( c(-delta,delta), size=1, prob=c(1/2,1/2) ) 
   #x <- phi*mu + sig*(rnorm(1)+sgn)
   x <- phi*mu + sig*(rnorm(1)+sgn)*sqrt(abs(mu))
   x
 }

B <- 1000
x <- numeric(B)
x[1] <- rnorm(1)
for( i in 2:B){ x[i] <- rbimod(x[i-1]) }
y <- x[2:B] - x[1:(B-1)] ## the deviations

tfun <- approxfun( density(y) )

##don't reorder...but keep track of where points go

x.from <- x[1:(B-1)]
x.to <- x[2:B]

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

#image( x.from.o, x.to.o2, env2 )
# rank transform
image( sqrt(env2) )

plot( x.from.o, x.to.o )

## it works
## to do
## show how k200, for example, approximates the generative kernel at that x.from position...
#
