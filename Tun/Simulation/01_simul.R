library(monomvn) # blasso
library(MASS) # mvrnorm
library(glmnet) 
library(tidyverse)

# dimension 
p = 10

# number of true non-zero regression parameters
q = 6

# number of simulated datasets
nsimul = 500

# number of posterior samples to draw
B = 3000

# (MCMC) burn-in period
burn_in = 3000

# specify variance-covariance matrix of predictors
cov_X <- diag(p)
for(i in 1:q){
  for(j in 1:q){
    if(i != j){
      cov_X[i,j] = .3^(abs(i-j))
    }
  }
}

# specify true beta
beta_true = rep(0,p)
for(j in 1:q){
  beta_true[j] = .75 + .25*j
}

#---------------------------------------------------------------------------------
# wrap the following into a simulation function of input n, error_dist and conf_lvl

# sample size
n = 500
# error distribution
error_dist = "shifted_chisq"
# confidence level
conf_lvl = .9 # 90% confidence interval and credible interval

# user-defined function
CI_func <- function(beta_matrix, half_lvl = (1 - conf_lvl)/2){
  t(apply(beta_matrix, 2, function(x){
    unname( quantile(x, probs = c(half_lvl, 1 - half_lvl)) ) })) %>% 
    as.data.frame() %>% 
    setNames(., c("lower", "upper")) %>%
    mutate(CI_length = upper - lower,
           upper_cover = ifelse(upper >= beta_true,1,0),
           lower_cover = ifelse(lower <= beta_true,1,0)) %>%
    mutate(coverage = upper_cover * lower_cover)
}

# initialize probability of variables being selected (averaged across simulated datasets)
select_MCMC <- rep(0, p)
select_RW1 <- rep(0, p)
select_RW2 <- rep(0, p)
select_RW3 <- rep(0, p)
select_residboot <- rep(0, p)

# initialize coverage of CI (averaged across simulated datasets)
cover_MCMC <- rep(0, p)
cover_RW1 <- rep(0, p)
cover_RW2 <- rep(0, p)
cover_RW3 <- rep(0, p)
cover_residboot <- rep(0, p)

# initialize width of CI (averaged across simulated datasets)
width_MCMC <- rep(0, p)
width_RW1 <- rep(0, p)
width_RW2 <- rep(0, p)
width_RW3 <- rep(0, p)
width_residboot <- rep(0, p)

for(simul in 1:nsimul){
  # simulate dataset
  set.seed(n*B+simul)
  X = mvrnorm(n = n, mu = rep(0,p), Sigma = cov_X)
  if(error_dist=="shifted_chisq"){
    eps = rchisq(n=n,df=2)-2
  } else{
    eps = rnorm(n=n)
  }
  Y = X %*% beta_true + eps
  
  # scale Y and X
  Y = scale(Y)
  X = scale(X)
  
  # Bayesian inference
  blasso_fit <- blasso(X = X, y = Y, T = B + burn_in)
  beta_MCMC <- blasso_fit$beta[(burn_in + 1) : (B + burn_in),]
  
  # Random-weighting
  
  # first get lambda
  res0 = cv.glmnet(x = X, y = Y); plot(res0)
  lamb = res0$lambda.min
  
  beta_RW1 <- matrix(NA, nrow = B, ncol = p) # no penalty weight
  beta_RW2 <- matrix(NA, nrow = B, ncol = p) # common penalty weight
  beta_RW3 <- matrix(NA, nrow = B, ncol = p) # different penalty weights
  for(b in 1:B){
    w = rexp(n)
    wp = rexp(p)
    RW1_fit = glmnet(x = X, y = Y, weights = w)
    RW3_fit = glmnet(x = X, y = Y, weights = w, penalty.factor = wp)
    
    beta_RW1[b,] <- as.vector(coef(RW1_fit, s = lamb))[-1]
    beta_RW2[b,] <- as.vector(coef(RW1_fit, s = lamb*wp[1]))[-1]
    beta_RW3[b,] <- as.vector(coef(RW3_fit, s = lamb*mean(wp)))[-1]
   }
  
  # Residual Bootstrap
  beta_res0 <-  coef(res0, s = lamb)[-1]
  resid_res0 <- Y - X %*% beta_res0
  beta_residboot <- matrix(NA, nrow = B, ncol = p)
  for(b in 1:B){
    resid_boot <- base::sample(x = resid_res0, size = n, replace = T)
    Y_boot <- X %*% beta_res0 + resid_boot
    
    Y_boot <- scale(Y_boot)
    fit_boot <- glmnet(x = X, y = Y_boot)
    beta_residboot[b,] <- as.vector(coef(fit_boot, s = lamb))[-1]
  }
  
  # aggregate results
  
  # probability of variables being selected 
  select_MCMC <- select_MCMC + apply(beta_MCMC, 2, function(x){
    1 - sum(as.numeric(x == 0))/B})
  select_RW1 <- select_RW1 + apply(beta_RW1, 2, function(x){
    1 - sum(as.numeric(x == 0))/B})
  select_RW2 <- select_RW2 + apply(beta_RW2, 2, function(x){
    1 - sum(as.numeric(x == 0))/B})
  select_RW3 <- select_RW3 + apply(beta_RW3, 2, function(x){
    1 - sum(as.numeric(x == 0))/B})
  select_residboot <- select_residboot + apply(beta_residboot, 2, function(x){
    1 - sum(as.numeric(x == 0))/B})
  
  # quantile
  CI_MCMC <- CI_func(beta_MCMC)
  CI_RW1 <- CI_func(beta_RW1)
  CI_RW2 <- CI_func(beta_RW2)
  CI_RW3 <- CI_func(beta_RW3)
  CI_residboot <- CI_func(beta_residboot)
  
  # coverage of CI 
  cover_MCMC <- cover_MCMC + CI_MCMC$coverage
  cover_RW1 <- cover_RW1 + CI_RW1$coverage
  cover_RW2 <- cover_RW2 + CI_RW2$coverage
  cover_RW3 <- cover_RW3 + CI_RW3$coverage
  cover_residboot <- cover_residboot + CI_residboot$coverage
  
  # width of CI 
  width_MCMC <- width_MCMC + CI_MCMC$CI_length
  width_RW1 <- width_RW1 + CI_RW1$CI_length
  width_RW2 <- width_RW2 + CI_RW2$CI_length
  width_RW3 <- width_RW3 + CI_RW3$CI_length
  width_residboot <- width_residboot + CI_residboot$CI_length
}

# probability of variables being selected (averaged across simulated datasets)
select_MCMC <- select_MCMC/nsimul
select_RW1 <- select_RW1/nsimul
select_RW2 <- select_RW2/nsimul
select_RW3 <- select_RW3/nsimul
select_residboot <- select_residboot/nsimul

# coverage of CI (averaged across simulated datasets)
cover_MCMC <- cover_MCMC/nsimul
cover_RW1 <- cover_RW1/nsimul
cover_RW2 <- cover_RW2/nsimul
cover_RW3 <- cover_RW3/nsimul
cover_residboot <- cover_residboot/nsimul

# width of CI (averaged across simulated datasets)
width_MCMC <- width_MCMC/nsimul
width_RW1 <- width_RW1/nsimul
width_RW2 <- width_RW2/nsimul
width_RW3 <- width_RW3/nsimul
width_residboot <- width_residboot/nsimul
