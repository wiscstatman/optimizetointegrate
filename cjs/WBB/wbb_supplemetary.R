# Title: Weighted Bayesian Bootstrap for Scalable Posterior Distributions
# Authors: Michael Newton, Nick Polson, Jianeng Xu
# This file contains the code for the diabetes example, in Section 3.2

# Packages --------------
library(glmnet)
library(lars)
library(ggplot2)
library(reshape2)
library(monomvn)
# The diabetes dataset is contained in "lars" package. 
# There are 442 rows and 3 columns. 
# Reference: Efron et al 2004, Least Angle Regression.
data("diabetes") 


# Preparation --------------
# number of observations
N = dim(diabetes$x)[1]
# number of variables
p = dim(diabetes$x)[2] 
# posterior sample size
B = 1000
# variable names
name = colnames(diabetes$x)
# center the data
y = diabetes$y
X = scale(diabetes$x)

# choose regularization parameter
res0 = cv.glmnet(X,y,nfolds = 10)
l = res0$lambda.min


# Posterior sampling --------------
# WBB: common prior weight
beta_common = matrix(NA,B,p)
colnames(beta_common) = name
for(t in 1:B)
{
  # generate weights and compute the corresponding solution
  res = glmnet(X,y,lambda = rexp(1)*l, weights = rexp(N))
  beta_common[t,] = as.vector(res$beta)
}

# WBB: separate prior weight
beta_separate = matrix(NA,B,p)
colnames(beta_separate) = name
for(t in 1:B)
{
  res = glmnet(X,y,lambda = rexp(1)*l, weights = rexp(N), penalty.factor = rexp(p))
  beta_separate[t,] = as.vector(res$beta)
}

# BLASSO
# generate Bayesian Lasso samples
blasso_res = blasso(X,y,T=2*B,verb = 0)
beta_blasso = tail(blasso_res$beta,B)
colnames(beta_blasso) = name


# Plot --------------
beta_common = melt(beta_common)[,-1]
beta_separate = melt(beta_separate)[,-1]
beta_blasso = melt(beta_blasso)[,-1]

beta = rbind(beta_common,beta_separate,beta_blasso)
method = c(rep("WBB_one_prior",B*p),rep("WBB",B*p),rep("BLASSO",B*p))
beta = cbind(method,beta)
colnames(beta)[2:3] = c("variable","beta")

p = ggplot(beta, aes(beta,colour = method)) +
  geom_density(alpha = 1) + 
  facet_wrap(~ variable,scales = "free",ncol = 2)+
  theme_bw()

print(p)