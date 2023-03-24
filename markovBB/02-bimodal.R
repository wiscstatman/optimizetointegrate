
# March 23/23

# a synthetic example; generated as AR, but with goofy bimodal conditionals

set.seed(312345126)

** not working as expected

rbimod <- function(mu, delta=2, sig=(1/2) )
 {
   # make a biomodal draw centered at mu
   sgn <- sample( c(-delta,delta), size=1, prob=c(1/2,1/2) ) 
   x <- mu + sig*sgn*rnorm(1)
   x
 }

B <- 1000
x <- numeric(B)
x[1] <- rnorm(1)
phi <- 0.7
for( i in 2:B){ x[i] <- rbimod(x[i-1]) }

oo <- order(x)

env <- outer(x,x,function(u,v){ exp( -(u-v)^2 ) } )
e2 <- env[oo,oo]

umass <- rowSums(e2)
TM <- e2/umass ## 'transition matrix'

