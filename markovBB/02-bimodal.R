
# March 23/23

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

oo <- order(x)

#env <- outer(x,x,function(u,v){  tfun(u-v) } )
#env[is.na(env)] <- 0

tau <- 1/10
env <- outer(x,x,function(u,v){ exp( -(1/tau)*(u-v)^2 ) } )



e2 <- env[oo,oo]

umass <- rowSums(e2)
TM <- e2/umass ## 'transition matrix'

sx <- sort(x)
cmean <- apply( TM, 1, function(p){ sum(p*sx) } )
cvar <- apply( TM, 1, function(p){ sum(p*sx^2) - ( sum(p*sx))^2 } )

## we get the variance right for small absolute mu
