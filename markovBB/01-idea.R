
# idea from Feb 9, 2023, to adapt/extend the weighted bootstrap to Markov processes

# a synthetic example; generated as AR

B <- 1000
x <- numeric(B)
x[1] <- rnorm(1)
phi <- 0.7
for( i in 2:B){ x[i] <- phi*x[i-1] + rnorm(1) }

## Now from the analysts perspective, we get one x[i] at a time
## we start imagining there is an urn, say U_1(y) at every point y in the space (same space as x[i])
## suppose the notation is U_1(y) means at time 1, at position y, we have a measure, which we'll
## suppose for starters is alpha*N(y,1),  and let's also condition on x[1]

## Update rule: U_2(y) = U_1(y) + pointmass(x[1]) * exp( -dist( x[1], y ) )
## so all urns update, but the mass on the points diminishes the farther it is from the prior value

## Then, as a predictive matter, x[2] ~ U_2(x[1])/normalizer

## Continuing, U_n(y) = U_{n-1}(y) + pointmass( x[n-1] ) * exp( -dist( x[n-1], y) ) 

## and x[n] ~ U_n( x[n-1] )/normalizer

## After running out n, we then let alpha -> 0, to get some discrete measures, one for each data point
## hmm

## It seems like all urns, over time, will get more full; 

## I'm not sure if the urn story, which is predictive, fits with the random weight story, which is
## posterior...

## Feb 18...the urn will contain measures; limiting alpha->0 will give a discrete measure,and
## so Dirichlet distributed weights if we want posteriors...

env <- outer(x,x,function(u,v){ exp( -(u-v)^2 ) } )

# Each Urn contains point masses at all data points; but we don't give unit masses, rather we
# give masses proportional equal to exponential of negative distance...[could include a decay parameter
# like exp( -(u-v)^2 * decay ) 
# these form an envelope of sorts, 

# if we pick any point as the 'last state', the 'row' makes a measure; sort of like a TM
oo <- order(x)

e2 <- env[oo,oo]

umass <- rowSums(e2)
TM <- e2/umass ## 'transition matrix'

## cool...each row of TM now is the center of a Dirichlet dist over row-specific random weights
## it's curious that image(e2) looks nice; but image(TM) doesn't scale as well
## anyway, e2 is the envelope set of measures, so reasonable to display

## hunch; what if X is just some Markov process, with transition kernel unspecified, but imbued with
## some local regularity...I wonder if the measures will converge to the true kernel, say possibly
## if decay diminishes...
## There might also be a UQ type of theorem on posterior probability near the "kernel"...hmmm
## lots to do!


