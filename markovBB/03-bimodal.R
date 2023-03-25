
# March 24/23

# 01 and 02 are offbase on how to the weights
# it's better if a weight associated with distance affects the next point that
# the chain reaches, not the current point

# a synthetic example; generated as AR, but with goofy bimodal conditionals

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

tau <- 1/10
env <- outer(x.from,x.from,function(u,v){ exp( -(1/tau)*(u-v)^2 ) } )

oo <- order(x.from)
env2 <- env[oo,oo]

# now apply this envelope to x.to

x.to2 <- x.to[oo] # tracks the output of each ordered state in x.from

# Now we realign to x.to2

** something's not right...

oo2 <- order( x.to2 )
#env3 <- env2[,oo2]
env3 <- env2[oo2,]

#
#sx <- sort(x)
#cmean <- apply( TM, 1, function(p){ sum(p*sx) } )
#cvar <- apply( TM, 1, function(p){ sum(p*sx^2) - ( sum(p*sx))^2 } )
