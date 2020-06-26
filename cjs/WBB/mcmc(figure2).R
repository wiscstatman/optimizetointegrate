library(mcmc)
library(glmnet)
library(MASS)
library(reshape2)
library(ggplot2)
library(BoomSpikeSlab)

# data
set.seed(123)
x = mvrnorm(n=50,mu=c(0,0),Sigma = matrix(c(1,0.3,0.3,1),2,2))
x1 = x[,1]
x2 = x[,2]
beta = c(0,2)
y = x1*beta[1] + x2*beta[2] + rnorm(50)


# prior 
lambda = 10

# WBB
B = 5000
wbb_sample = NULL
N = length(y)
for(i in 1:B)
{
  print(i)
  w = rexp(N)
  w0 = rexp(2)
  res = glmnet(x,y,weights = w, intercept = FALSE, penalty.factor = w0)
  wbb_sample = rbind(wbb_sample, as.vector(coef(res,s=lambda/N)[-1]))
}

plot(density(wbb_sample[,1]))

# mcmc
easyMCMC <- function(target, niter, startval, proposalsd) {
  x <- matrix(0, niter,2)
  x[1,] <- startval
  for (i in 2:niter) {
    currentx <- x[i - 1,]
    proposedx = c(0,0)
    proposedx[1] <- rnorm(1, mean = currentx[1], sd = proposalsd)
    proposedx[2] <- rnorm(1, mean = currentx[2], sd = proposalsd)
    A <- target(proposedx) / target(currentx)
    print(target(currentx))
    print(target(proposedx))
    if (runif(1) < A) {
      x[i,] <- proposedx
    } else {
      x[i,] <- currentx
    }
  }
  return(x)
}

f_l1 = function(beta)
{
  exp(-0.5*sum((y-x1*beta[1]-x2*beta[2])^2) - lambda*(abs(beta[1])+abs(beta[2])))
}

f_l0 = function(beta)
{
  exp(-0.5*sum((y-x1*beta[1]-x2*beta[2])^2) - lambda*((abs(beta[1])>0)+(abs(beta[2])>0)))
}


f_l0_500 = function(beta)
{
  exp(-0.5*sum((y-x1*beta[1]-x2*beta[2])^2) - 10*lambda*((abs(beta[1])>0)+(abs(beta[2])>0)))
}



mcmc_l1sample = tail(easyMCMC(f_l1,50000,c(2,2),0.2),5000)
plot(mcmc_l1sample[,1])
plot(mcmc_l1sample[,2])


mcmc_l0sample = tail(easyMCMC(f_l0,50000,c(2,2),0.2),5000)
plot(mcmc_l0sample[,1])
plot(mcmc_l0sample[,2])

mcmc_l05sample = tail(easyMCMC(f_l0_500,50000,c(2,2),0.2),5000)
plot(mcmc_l05sample[,1])
plot(mcmc_l05sample[,2])


# mcmc_ss = tail(lm.spike(y ~ x, niter = 50000, error.distribution = "student")$beta[,-1], 5000)

# plot
sample = data.frame(rbind(wbb_sample,mcmc_l1sample,mcmc_l0sample,mcmc_l05sample))
colnames(sample) = c("beta_1",'beta_2')
rownames(sample) = NULL

method = rep(c('WBB','MCMC_L1','MCMC_L0','MCMC_L0_5'),each=5000)
sample = cbind(method, sample)
sample_df = melt(sample,id.vars = 'method')

p <- ggplot(sample_df,aes(x=value,color=method)) + 
  geom_density() + 
  facet_grid(variable~., scales="free") + 
  theme_bw() + 
  xlab("coefficient")+
  theme(text = element_text(size=21))
pdf('mcmc.pdf',height = 6,width = 8)
print(p)
dev.off()
