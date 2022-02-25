library(monomvn) # blasso
library(MASS) # mvrnorm
library(glmnet) 
library(parallel)
library(doParallel)
# library(doSNOW) 
library(ggplot2)
library(gridExtra)
library(reshape2) # melt
library(xtable)

#####################################################################################
#                             some user-defined functions                           #
#####################################################################################

CV2step.func <- function(Y, X, lambda_path, cv_K = 4){
  # initiate cv MSE
  cv_MSE <- rep(NA, length(lambda_path))
  
  # training sample size
  n = length(Y)
  
  # cv partition
  part_count <- rep(floor(n/cv_K), cv_K)
  part_count[cv_K] <- n - sum(head(part_count,-1)) 
  cv_part <- rep(1:cv_K, times = part_count)
  cv_part <- sample(cv_part, size = n, replace = F)
  
  # compute CV MSE  
  for(l in 1:length(lambda_path)){
    lamb <- lambda_path[l]
    cv_SSE <- 0
    
    for(k in 1:cv_K){
      Y_cvTrain <- Y[(cv_part != k)]
      X_cvTrain <- X[(cv_part != k),]
      Y_cvTest <- Y[(cv_part == k)]
      X_cvTest <- X[(cv_part == k),]
      
      cvfit <- glmnet(x = X_cvTrain, y = Y_cvTrain, intercept = F)
      beta_cvfit <- as.vector( coef(cvfit, s = lamb)[-1] )
      if( sum(abs(beta_cvfit)) != 0 ){
        nonzero_cvfit <- as.vector( which(beta_cvfit != 0) )
        beta_cvfit[nonzero_cvfit] <- as.vector( coef(lm(Y_cvTrain ~ X_cvTrain[, nonzero_cvfit] - 1)) )
      }
      cv_SSE <- cv_SSE + sum( (Y_cvTest -  X_cvTest %*%beta_cvfit)^2 )
    }
    
    cv_MSE[l] <- cv_SSE/n
  }
  
  # get lamb.min
  return( lambda_path[which.min(cv_MSE)] )
}

MSE.func <- function(beta_mat, Y, X){
  mean( apply(beta_mat, 1, function(b){
    mean( (Y - X %*% b)^2 )
  }) )
}

MSPE.func <- function(beta_mat, Y_test, X_test){
  mean( apply(beta_mat, 1, function(b){
    mean( (Y_test - X_test %*% b)^2 )
  }) )
}

select.func <- function(beta_mat, B){
  apply(beta_mat, 2, function(b){
    1 - sum(as.numeric(b == 0))/B
  })
} 

CIbound.func <- function(beta_mat, half_lvl){
  t( apply(beta_mat, 2, function(b){
    as.vector( quantile(b, probs = c(half_lvl, 1 - half_lvl)) )
  }) )
} 

TV.func <- function(f1, f2, tv_seq){
  tv_n = length(tv_seq)
  # pmin_f =  pmin(f1,f2)
  # TV = 1 - ( diff(tv_seq) %*% (pmin_f[2:tv_n] + pmin_f[1:(tv_n-1)])/2 )
  abs_f = abs(f1 - f2)
  TV = .5*( diff(tv_seq) %*% (abs_f[2:tv_n] + abs_f[1:(tv_n-1)])/2 )
  return( TV )
}

simul.func <- function(n, simul_seed, error_dist="std_norm", covar_type="yes", p=10, q=6, store_beta="yes", conf_lvl=.9, nsimul=500, burn_in=2000, B=1000){

  # specify true beta
  beta_true = rep(0,p)
  for(j in 1:q){
    beta_true[j] = .75 + .25*j
  }
  
  # types of covariance matrix
  if(covar_type == "yes"){
    # specify variance-covariance matrix of predictors (strong irrep)
    cov_X <- diag(p)
    for(i in 1:q){
      for(j in 1:q){
        if(i != j){
          cov_X[i,j] = .3^(abs(i-j))
        }
      }
    }
  } else{
    # specify variance-covariance matrix of predictors (NOT strong irrep)
    cov_X <- diag(p)
    for(i in 1:p){
      for(j in 1:p){
        if(i != j){
          cov_X[i,j] = .5
        }
      }
    }
    for(i in 1:q){
      for(j in 1:q){
        if(i != j){
          cov_X[i,j] = .4
        }
      }
    }
  } 
  
  # range of to calculate total variation for ecdfs
  tv_range <- c(
    seq(-3.5, -0.001, by =.001),
    -.0001, .0001, 
    seq(.001, 6, by = .001)
  )
  
  # alpha/2
  half_lvl = (1 - conf_lvl)/2
  
  # initialize cache for all simulated datasets
  
  # each MSE value is averaged across B posterior samples
  MSE_MCMC <- rep(NA, nsimul)
  MSE_RW1 <- rep(NA, nsimul)
  MSE_RW2 <- rep(NA, nsimul)
  MSE_RW3 <- rep(NA, nsimul)
  MSE_ResidBoot <- rep(NA, nsimul)
  
  # each MSPE value is averaged across B posterior samples
  MSPE_MCMC <- rep(NA, nsimul)
  MSPE_RW1 <- rep(NA, nsimul)
  MSPE_RW2 <- rep(NA, nsimul)
  MSPE_RW3 <- rep(NA, nsimul)
  MSPE_ResidBoot <- rep(NA, nsimul)
  
  # probability of selecting each component of beta: 1, ..., p
  select_MCMC <- matrix(NA, nrow = nsimul, ncol = p)
  select_RW1 <- matrix(NA, nrow = nsimul, ncol = p)
  select_RW2 <- matrix(NA, nrow = nsimul, ncol = p)
  select_RW3 <- matrix(NA, nrow = nsimul, ncol = p)
  select_ResidBoot <- matrix(NA, nrow = nsimul, ncol = p)
  
  # coverage probability for each component of beta: 1, ..., p
  Cover_MCMC <- rep(0, p)
  Cover_RW1 <- rep(0, p)
  Cover_RW2 <- rep(0, p)
  Cover_RW3 <- rep(0, p)
  Cover_ResidBoot <- rep(0, p)
  
  # CI width (averaged across nsimul datasets) for each beta: 1, ..., p
  CIwidth_MCMC <- rep(0, p)
  CIwidth_RW1 <- rep(0, p)
  CIwidth_RW2 <- rep(0, p)
  CIwidth_RW3 <- rep(0, p)
  CIwidth_ResidBoot <- rep(0, p)
  
  # each value is averaged across p components of betas
  TV_RW1 <- rep(0, nsimul)
  TV_RW2 <- rep(0, nsimul)
  TV_RW3 <- rep(0, nsimul)
  TV_ResidBoot <- rep(0, nsimul)
  
  # store lambdas to compare
  lamb_1step <- rep(NA, nsimul)
  lamb_2step <- rep(NA, nsimul)
  
  # store betas
  if(store_beta == "yes"){
    beta_MCMC_cache <- array(NA, dim=c(B,p,nsimul))
    beta_others_cache <- array(NA, dim=c(B,4*p,nsimul))
  }
  
  for(dts in 1:nsimul){
    
    # simulate dataset
    set.seed(2*B*simul_seed + dts)
    
    X = mvrnorm(n = n, mu = rep(0,p), Sigma = cov_X)
    X_test = mvrnorm(n = n, mu = rep(0,p), Sigma = cov_X)
    
    if(error_dist=="shifted_chisq"){
      eps = rchisq(n=n,df=2)-2
      eps_test = rchisq(n=n,df=2)-2
    } else{
      eps = rnorm(n=n)
      eps_test = rnorm(n=n)
    }
    
    Y = X %*% beta_true + eps
    Y_test = X_test %*% beta_true + eps_test
    
    #-----------------------------------------------------------------
    # Bayesian inference
    blasso_fit <- blasso(X = X, y = Y, T = B + burn_in, icept = F, verb = 0)
    beta_MCMC <- blasso_fit$beta[(burn_in + 1) : (B + burn_in),]
    
    # performance measure for MCMC samples
    
    MSE_MCMC[dts] <- MSE.func(beta_MCMC, Y=Y, X=X) 
    MSPE_MCMC[dts] <- MSPE.func(beta_MCMC, Y_test=Y_test, X_test=X_test)
    select_MCMC[dts,] <- select.func(beta_MCMC, B=B) 
    
    CIbound_MCMC <- CIbound.func(beta_MCMC, half_lvl=half_lvl)  
    CIwidth_MCMC <- CIwidth_MCMC + as.vector(CIbound_MCMC[,2] - CIbound_MCMC[,1])
    Cover_MCMC <- Cover_MCMC + as.integer(beta_true >= CIbound_MCMC[,1]) * as.integer(beta_true <= CIbound_MCMC[,2])
    
    if(store_beta == "yes"){
      beta_MCMC_cache[,,dts] <- beta_MCMC
    }
    
    #-----------------------------------------------------------------
    # Random-weighting & Residual Bootstrap
    
    # prep for residual bootstrap
    res0 <- cv.glmnet(x = X, y = Y, nfolds = 4, intercept = F)
    beta_res0 <-  as.vector( coef(res0, s = res0$lambda.min)[-1] )
    resid_res0 <- Y - X %*% beta_res0
    resid_res0 <- resid_res0 - mean(resid_res0)
    Xbeta_res0 <- X %*% beta_res0 
    lamb_1step[dts] <- res0$lambda.min
    
    # CV for 2-step RW's
    # lamb <- CV2step.func(Y=Y, X=X, lambda_path=res0$lambda, cv_K = 4)
    lamb <- CV2step.func(Y=Y, X=X, lambda_path=exp(seq(-5,.5,by=.1)), cv_K = 4)
    lamb_2step[dts] <- lamb
    # lamb <- res0$lambda.min
    
    # Use parallel computing
    NumCore <- detectCores()-1
    registerDoParallel(cores = NumCore)
    # cl <- makeCluster(NumCore)
    # registerDoSNOW(cl) 
    beta_others <- foreach(i = 1:B, .combine = "rbind") %dopar%{
      # RW schemes 1,2,3
      w <- rexp(n)
      wp <- rexp(p)
      fit_w1 <- glmnet::glmnet(x = X, y = Y, weights = w, intercept = F)
      fit_w3 <- glmnet::glmnet(x = X, y = Y, weights = w, penalty.factor = wp, intercept = F)
      
      b_RW1 <- as.vector( coef(fit_w1, s = lamb)[-1] )
      b_RW2 <- as.vector( coef(fit_w1, s = lamb*wp[1])[-1] )
      b_RW3 <- as.vector( coef(fit_w3, s = lamb*mean(wp))[-1] )
      
      if( sum(abs(b_RW1)) != 0 ){
        nonzero_RW1 <- as.vector( which(b_RW1 != 0) )
        b_RW1[nonzero_RW1] <- as.vector( coef(lm(Y ~ X[, nonzero_RW1] - 1, weights = w)) )
      }
      
      if( sum(abs(b_RW2)) != 0 ){
        nonzero_RW2 <- as.vector( which(b_RW2 != 0) )
        b_RW2[nonzero_RW2] <- as.vector( coef(lm(Y ~ X[, nonzero_RW2] - 1, weights = w)) )
      }
      
      if( sum(abs(b_RW3)) != 0 ){
        nonzero_RW3 <- as.vector( which(b_RW3 != 0) )
        b_RW3[nonzero_RW3] <- as.vector( coef(lm(Y ~ X[, nonzero_RW3] - 1, weights = w)) )
      }
      
      # Residual bootstrap
      resid_boot <- base::sample(x = resid_res0, size = n, replace = T)
      Y_boot <- Xbeta_res0 + resid_boot
      fit_boot <- glmnet::glmnet(x = X, y = Y_boot, intercept = F)
      
      # gets all non-MCMC betas
      return(c(
        b_RW1, b_RW2, b_RW3,
        as.vector( coef(fit_boot, s = res0$lambda.min)[-1] )
      ))
    }
    # stopCluster(cl)
    registerDoSEQ()
    
    if(store_beta == "yes"){
      beta_others_cache[,,dts] <- beta_others
    }
    
    beta_RW1 <- beta_others[, 1:p]
    beta_RW2 <- beta_others[, (p+1):(2*p)]
    beta_RW3 <- beta_others[, (2*p+1):(3*p)]
    beta_ResidBoot <- beta_others[, (3*p+1):(4*p)]
    
    #-----------------------------------------------------------------
    # Performance measure for all other betas
    
    MSE_RW1[dts] <- MSE.func(beta_RW1, Y=Y, X=X) 
    MSE_RW2[dts] <- MSE.func(beta_RW2, Y=Y, X=X) 
    MSE_RW3[dts] <- MSE.func(beta_RW3, Y=Y, X=X) 
    MSE_ResidBoot[dts] <- MSE.func(beta_ResidBoot, Y=Y, X=X) 
    
    MSPE_RW1[dts] <- MSPE.func(beta_RW1, Y_test=Y_test, X_test=X_test)
    MSPE_RW2[dts] <- MSPE.func(beta_RW2, Y_test=Y_test, X_test=X_test)
    MSPE_RW3[dts] <- MSPE.func(beta_RW3, Y_test=Y_test, X_test=X_test)
    MSPE_ResidBoot[dts] <- MSPE.func(beta_ResidBoot, Y_test=Y_test, X_test=X_test)
    
    select_RW1[dts,] <- select.func(beta_RW1, B=B) 
    select_RW2[dts,] <- select.func(beta_RW2, B=B) 
    select_RW3[dts,] <- select.func(beta_RW3, B=B) 
    select_ResidBoot[dts,] <- select.func(beta_ResidBoot, B=B) 
    
    CIbound_RW1 <- CIbound.func(beta_RW1, half_lvl=half_lvl)  
    CIbound_RW2 <- CIbound.func(beta_RW2, half_lvl=half_lvl)  
    CIbound_RW3 <- CIbound.func(beta_RW3, half_lvl=half_lvl)  
    CIbound_ResidBoot <- CIbound.func(beta_ResidBoot, half_lvl=half_lvl)  
    
    CIwidth_RW1 <- CIwidth_RW1 + as.vector(CIbound_RW1[,2] - CIbound_RW1[,1])
    CIwidth_RW2 <- CIwidth_RW2 + as.vector(CIbound_RW2[,2] - CIbound_RW2[,1])
    CIwidth_RW3 <- CIwidth_RW3 + as.vector(CIbound_RW3[,2] - CIbound_RW3[,1])
    CIwidth_ResidBoot <- CIwidth_ResidBoot + as.vector(CIbound_ResidBoot[,2] - CIbound_ResidBoot[,1])
    
    Cover_RW1 <- Cover_RW1 + as.integer(beta_true >= CIbound_RW1[,1]) * as.integer(beta_true <= CIbound_RW1[,2])
    Cover_RW2 <- Cover_RW2 + as.integer(beta_true >= CIbound_RW2[,1]) * as.integer(beta_true <= CIbound_RW2[,2])
    Cover_RW3 <- Cover_RW3 + as.integer(beta_true >= CIbound_RW3[,1]) * as.integer(beta_true <= CIbound_RW3[,2])
    Cover_ResidBoot <- Cover_ResidBoot + 
      as.integer(beta_true >= CIbound_ResidBoot[,1]) * as.integer(beta_true <= CIbound_ResidBoot[,2])
    
    #-----------------------------------------------------------------
    # total variation between ecdf of beta_MCMC and ecdf of other betas
    
    for (j in 1:p){
      ecdf_MCMC <- ecdf(beta_MCMC[,j])(tv_range)
      ecdf_RW1 <- ecdf(beta_RW1[,j])(tv_range)
      ecdf_RW2 <- ecdf(beta_RW2[,j])(tv_range)
      ecdf_RW3 <- ecdf(beta_RW3[,j])(tv_range)
      ecdf_ResidBoot <- ecdf(beta_ResidBoot[,j])(tv_range)
      
      TV_RW1[dts] <- TV_RW1[dts] + TV.func(f1 = ecdf_MCMC, f2 = ecdf_RW1, tv_seq = tv_range)
      TV_RW2[dts] <- TV_RW2[dts] + TV.func(f1 = ecdf_MCMC, f2 = ecdf_RW2, tv_seq = tv_range)
      TV_RW3[dts] <- TV_RW3[dts] + TV.func(f1 = ecdf_MCMC, f2 = ecdf_RW3, tv_seq = tv_range) 
      TV_ResidBoot[dts] <- TV_ResidBoot[dts] + 
        TV.func(f1 = ecdf_MCMC, f2 = ecdf_ResidBoot, tv_seq = tv_range) 
    }
    TV_RW1[dts] <- TV_RW1[dts]/p
    TV_RW2[dts] <- TV_RW2[dts]/p
    TV_RW3[dts] <- TV_RW3[dts]/p
    TV_ResidBoot[dts] <- TV_ResidBoot[dts]/p
    
  }
  
  Cover_MCMC <- Cover_MCMC/nsimul
  Cover_RW1 <- Cover_RW1/nsimul
  Cover_RW2 <- Cover_RW2/nsimul
  Cover_RW3 <- Cover_RW3/nsimul
  Cover_ResidBoot <- Cover_ResidBoot/nsimul
  
  CIwidth_MCMC <- CIwidth_MCMC/nsimul
  CIwidth_RW1 <- CIwidth_RW1/nsimul
  CIwidth_RW2 <- CIwidth_RW2/nsimul
  CIwidth_RW3 <- CIwidth_RW3/nsimul
  CIwidth_ResidBoot <- CIwidth_ResidBoot/nsimul
  
  result <- list(
    MSE_MCMC = MSE_MCMC,
    MSE_RW1 = MSE_RW1,
    MSE_RW2 = MSE_RW2,
    MSE_RW3 = MSE_RW3,
    MSE_ResidBoot = MSE_ResidBoot,
    
    MSPE_MCMC = MSPE_MCMC,
    MSPE_RW1 = MSPE_RW1,
    MSPE_RW2 = MSPE_RW2,
    MSPE_RW3 = MSPE_RW3,
    MSPE_ResidBoot = MSPE_ResidBoot,
    
    Cover_MCMC = Cover_MCMC,
    Cover_RW1 = Cover_RW1,
    Cover_RW2 = Cover_RW2,
    Cover_RW3 = Cover_RW3,
    Cover_ResidBoot = Cover_ResidBoot,
    
    CIwidth_MCMC = CIwidth_MCMC,
    CIwidth_RW1 = CIwidth_RW1,
    CIwidth_RW2 = CIwidth_RW2,
    CIwidth_RW3 = CIwidth_RW3,
    CIwidth_ResidBoot = CIwidth_ResidBoot,
    
    TV_RW1 = TV_RW1,
    TV_RW2 = TV_RW2,
    TV_RW3 = TV_RW3,
    TV_ResidBoot = TV_ResidBoot,
    
    select_RW1 = select_RW1,
    select_RW2 = select_RW2,
    select_RW3 = select_RW3,
    select_ResidBoot = select_ResidBoot,
    select_MCMC = select_MCMC,
    
    lamb_1step = lamb_1step,
    lamb_2step = lamb_2step
  )
  
  if(store_beta == "yes"){
    result <- c(result, list(
      beta_MCMC_cache = beta_MCMC_cache,
      beta_others_cache = beta_others_cache
    ))
  }
  
  return(result)
}

#####################################################################################
#                                   run simulations                                 #
#####################################################################################

start_time <- Sys.time()
setting1 <- simul.func( n=100, simul_seed=1, error_dist="std_norm", covar_type="yes", p=10, q=6 )
Sys.time() - start_time

save(setting1, file = "Setting1.RData")
setting1_small <- setting1[!(names(setting1) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting1_small, file = "Setting1_small.RData")

#----------------------------------------------------------------------------------
rm(setting1, setting1_small); gc()
start_time <- Sys.time()
setting2 <- simul.func( n=500, simul_seed=2, error_dist="std_norm", covar_type="yes", p=10, q=6 )
Sys.time() - start_time

save(setting2, file = "Setting2.RData")
setting2_small <- setting2[!(names(setting2) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting2_small, file = "Setting2_small.RData")

#----------------------------------------------------------------------------------
rm(setting2, setting2_small); gc()
start_time <- Sys.time()
setting3 <- simul.func( n=100, simul_seed=3, error_dist="shifted_chisq", covar_type="yes", p=10, q=6 )
Sys.time() - start_time

save(setting3, file = "Setting3.RData")
setting3_small <- setting3[!(names(setting3) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting3_small, file = "Setting3_small.RData")

#----------------------------------------------------------------------------------
rm(setting3, setting3_small); gc()
start_time <- Sys.time()
setting4 <- simul.func( n=500, simul_seed=4, error_dist="shifted_chisq", covar_type="yes", p=10, q=6 )
Sys.time() - start_time

save(setting4, file = "Setting4.RData")
setting4_small <- setting4[!(names(setting4) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting4_small, file = "Setting4_small.RData")

#----------------------------------------------------------------------------------
rm(setting4, setting4_small); gc()
start_time <- Sys.time()
setting5 <- simul.func( n=100, simul_seed=5, error_dist="std_norm", covar_type="no", p=10, q=6 )
Sys.time() - start_time

save(setting5, file = "Setting5.RData")
setting5_small <- setting5[!(names(setting5) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting5_small, file = "Setting5_small.RData")

#----------------------------------------------------------------------------------
rm(setting5, setting5_small); gc()
start_time <- Sys.time()
setting6 <- simul.func( n=500, simul_seed=6, error_dist="std_norm", covar_type="no", p=10, q=6 )
Sys.time() - start_time

save(setting6, file = "Setting6.RData")
setting6_small <- setting6[!(names(setting6) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting6_small, file = "Setting6_small.RData")

#----------------------------------------------------------------------------------
rm(setting6, setting6_small); gc()
start_time <- Sys.time()
setting7 <- simul.func( n=100, simul_seed=7, error_dist="std_norm", covar_type="yes", p=50, q=6 )
Sys.time() - start_time

save(setting7, file = "Setting7.RData")
setting7_small <- setting7[!(names(setting7) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting7_small, file = "Setting7_small.RData")

#----------------------------------------------------------------------------------
rm(setting7, setting7_small); gc()
start_time <- Sys.time()
setting8 <- simul.func( n=500, simul_seed=8, error_dist="std_norm", covar_type="yes", p=50, q=6 )
Sys.time() - start_time

save(setting8, file = "Setting8.RData")
setting8_small <- setting8[!(names(setting8) %in% c("beta_MCMC_cache","beta_others_cache"))]
save(setting8_small, file = "Setting8_small.RData")

#####################################################################################
#                                  Plot and Tabulate                                #
#####################################################################################

load("Setting1_small.RData")
load("Setting2_small.RData")
load("Setting3_small.RData")
load("Setting4_small.RData")
load("Setting5_small.RData")
load("Setting6_small.RData")
load("Setting7_small.RData")
load("Setting8_small.RData")

nsimul = 500

# color palette
# correspond to MCMC, RW1, RW2, RW3, residualbootstrap
pal_simul <- c("gray50", "red", "blue","green", "orange")
pal_TV <- pal_simul[-1] # correspond to RW1, RW2, RW3, residualbootstrap

getval.func <- function(MCMC_item=1, item, list_vec = paste0("setting", 1:8,"_small")){
  vec <- NULL
  item_vec <- paste(item, c("MCMC", "RW1", "RW2", "RW3", "ResidBoot"), sep = "_")
  if(is.null(MCMC_item)){
    for(L in list_vec){
      list_L <- get(L)
      vec <- c(vec, c(
        list_L[names(list_L) == item_vec[2]][[1]],
        list_L[names(list_L) == item_vec[3]][[1]],
        list_L[names(list_L) == item_vec[4]][[1]],
        list_L[names(list_L) == item_vec[5]][[1]]
      ))
    }
  } else{
    for(L in list_vec){
      list_L <- get(L)
      vec <- c(vec, c(
        list_L[names(list_L) == item_vec[1]][[1]],
        list_L[names(list_L) == item_vec[2]][[1]],
        list_L[names(list_L) == item_vec[3]][[1]],
        list_L[names(list_L) == item_vec[4]][[1]],
        list_L[names(list_L) == item_vec[5]][[1]]
      ))
    }
  }
  return(vec)
}

#----------------------------------------------------------------------------------
# distributions of lambdas

lamb_vec <- NULL
list_vec <- paste0("setting", 1:8,"_small")
for(L in list_vec){
  list_L <- get(L)
  lamb_vec <- c(lamb_vec, c(
    list_L[names(list_L) == "lamb_1step"][[1]],
    list_L[names(list_L) == "lamb_2step"][[1]]
  ))
}
lambda_df <- data.frame(
  value = lamb_vec, 
  method = rep( rep(c("1-step","2-step"), each = nsimul), 8 ),
  Setting = rep( paste0("Setting_", c(1:8)), each = 2*nsimul )
)
ggplot(lambda_df, aes(x = value, color = method, fill = method)) +
  geom_histogram(aes(y=..density..), alpha=.5, position = "identity", bins=40) +
  facet_wrap(~Setting, ncol = 2, scales = "free") +
  labs(x = expr(paste("minimum ", lambda)), title = expr(paste("Density histograms of ",lambda," obtained from cross validation for 500 simulated datasets in 8 simulation settings"))) +
  theme(plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------------------------------
# MSE

MSE_df <- data.frame(
  value = getval.func(item = "MSE" ),
  method = rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = nsimul), 8 ),
  Setting = rep( paste0("Setting_", c(1:8)), each = 5*nsimul )
)
MSE_df$method <- factor(MSE_df$method, levels = c("MCMC", "RW1", "RW2", "RW3", "RB"))
ggplot(MSE_df, aes(x=method, y=value, fill=method)) +
  geom_boxplot(outlier.shape=".") +
  facet_wrap(~Setting, ncol = 2, scales = "free") +
  scale_fill_manual(values = pal_simul) +
  # coord_cartesian(ylim = c(0,7)) +
  labs(x="method", y= "MSE", title = "Boxplots of average MSEs for 500 simulated datasets in 8 simulation settings") + 
  theme(plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------------------------------
# MSPE

MSPE_df <- data.frame(
  value = getval.func(item = "MSPE" ),
  method = rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = nsimul), 8 ),
  Setting = rep( paste0("Setting_", c(1:8)), each = 5*nsimul )
)
MSPE_df$method <- factor(MSPE_df$method, levels = c("MCMC", "RW1", "RW2", "RW3", "RB"))
ggplot(MSPE_df, aes(x=method, y=value, fill=method)) +
  geom_boxplot(outlier.shape=".") +
  facet_wrap(~Setting, ncol = 2, scales = "free") +
  scale_fill_manual(values = pal_simul) +
  # coord_cartesian(ylim = c(0,8)) +
  labs(x="method", y= "MSPE", title = "Boxplots of average MSPEs for 500 simulated datasets in 8 simulation settings") + 
  theme(plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------------------------------
# Total Variation

TV_df <- data.frame(
  value = getval.func(MCMC_item = NULL ,item = "TV" ),
  method = rep( rep( c("RW1", "RW2", "RW3", "RB"), each = nsimul), 8 ),
  Setting = rep( paste0("Setting_", c(1:8)), each = 4*nsimul )
)
TV_df$method <- factor(TV_df$method, levels = c("RW1", "RW2", "RW3", "RB"))

# png("simul1_TV.png", width = 6, height = 6, units = "in", res = 1200)
setEPS()
postscript("simul1_TV.eps", width = 6, height = 7)
ggplot(TV_df, aes(x=method, y=value, fill=method)) +
  geom_boxplot(outlier.shape=".", width = .4) +
  facet_wrap(~Setting, ncol = 2, scales = "free") +
  scale_fill_manual(values = pal_TV) +
  # coord_cartesian(ylim = c(0,8)) +
  labs(x="method", y= expr(paste("Total variation (averaged across all ",beta,"'s)")), 
  # title = "Boxplots of total variation between ecdf(MCMC betas) and ecdf(non-MCMC betas) 
  #      for 500 simulated datasets in 8 simulation settings"
  title = "") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#---------------------------------------------------------------------------------
# Probability of selecting variable

p1 = 10
p2 = 50

select_list <- paste0("setting", 1:8,"_small")
select_vec <- paste("select", c("MCMC", "RW1", "RW2", "RW3", "ResidBoot"), sep = "_")
select_df <- data.frame()
for(L in select_list){
  list_L <- get(L)
  for(s in select_vec){
    select_df <- rbind(
      select_df, 
      melt(list_L[names(list_L) == s][[1]])[,-1]
    )
  }
}
colnames(select_df) <- c("variable", "value")
select_df$variable <- paste0("Variable_", select_df$variable)
select_df$variable <- factor(select_df$variable, levels = paste0("Variable_", 1:50))

select_df$method = c(
  rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = nsimul*p1), 6),
  rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = nsimul*p2), 2)
)
select_df$method <- factor(select_df$method, levels = c("MCMC", "RW1", "RW2", "RW3", "RB"))

select_df$Setting = c(
  rep( paste0("Setting_", c(1:6)), each = 5*nsimul*p1 ),
  rep( paste0("Setting_", c(7,8)), each = 5*nsimul*p2 )
)

# png("simul1_prob_select.png", width = 6, height = 7, units = "in", res = 1200)
setEPS()
postscript("simul1_prob_select.eps", width = 6, height = 7)
ggplot(select_df[select_df$variable %in% paste0("Variable_", c(1,7)),], aes(x = Setting, y = value)) +
  geom_boxplot(aes(fill = method), position = position_dodge2(preserve = "single"), outlier.shape=".") +
  facet_wrap(~variable, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = pal_simul) +
  labs(title = "",
    #title = "Distribution of Pr(select variable) across 500 simulated datasets",
       x = "Simulation Settings", y = "probability of selecting variable") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
dev.off()

#---------------------------------------------------------------------------------
# Coverage and average width of 90% CI

p1 = 10
p2 = 50

# get coverage data frame
# 5methods * (6settings * 10vars + 2settings * 50vars) = 800 rows
cover_df <- data.frame(
  value = getval.func(item = "Cover"),
  method = c(
    rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = p1), 6), # Settings 1-6
    rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = p2), 2)  # Settings 7&8
  ),
  variable = c(
    rep( c(1:10), 5*6),
    rep( c(1:50), 5*2)
  ),
  Setting = c(
    rep( paste0("Setting_", c(1:6)), each = 5*p1 ),
    rep( paste0("Setting_", c(7,8)), each = 5*p2 )
  )
)

# get average 90% CI data frame (each average value was taken across B=1000 posterior samples)
# 5methods * (6settings * 10vars + 2settings * 50vars) = 800 rows
aveCI_df <- data.frame(
  value = getval.func(item = "CIwidth"),
  method = c(
    rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = p1), 6), # Settings 1-6
    rep( rep( c("MCMC", "RW1", "RW2", "RW3", "RB"), each = p2), 2)  # Settings 7&8
  ),
  variable = c(
    rep( c(1:10), 5*6),
    rep( c(1:50), 5*2)
  ),
  Setting = c(
    rep( paste0("Setting_", c(1:6)), each = 5*p1 ),
    rep( paste0("Setting_", c(7,8)), each = 5*p2 )
  )
)

library(tidyverse)

# get table for coverage and CI width (each table is variable specific)
getCoverWidth_byVar.func <- function(var_j){
  cover_wider <- cover_df %>% 
    filter(variable == var_j) %>%
    mutate(val2 = round(value,3)) %>%
    pivot_wider(
      id_cols = Setting,
      names_from = method,
      values_from = val2
    )
  aveCI_wider <- aveCI_df %>% 
    filter(variable == var_j) %>%
    mutate(val2 = paste0("(",round(value,3),")")) %>%
    pivot_wider(
      id_cols = Setting,
      names_from = method,
      values_from = val2
    )
  res <- rbind(cover_wider, aveCI_wider) %>% 
    arrange(Setting) 
  return(res)
}

# get table for coverage and CI width (each table is setting specific)
getCoverWidth_bySetting.func <- function(set_j){
  cover_wider <- cover_df %>% 
    filter(Setting == paste0("Setting_",set_j)) %>%
    mutate(val2 = round(value,3)) %>%
    pivot_wider(
      id_cols = variable,
      names_from = method,
      values_from = val2
    )
  aveCI_wider <- aveCI_df %>% 
    filter(Setting == paste0("Setting_",set_j)) %>%
    mutate(val2 = paste0("(",round(value,3),")")) %>%
    pivot_wider(
      id_cols = variable,
      names_from = method,
      values_from = val2
    )
  res <- rbind(cover_wider, aveCI_wider) %>% 
    arrange(variable) 
  res <- res %>%
    mutate(variable = rep(c(seq(from=1,by=.25,length.out=6),rep(0, nrow(aveCI_wider)-6)),each = 2)) %>%
    rename(true_beta = variable)
  return(res)
}

# tabulate by variable
print(xtable(getCoverWidth_byVar.func(1)), include.rownames=FALSE)
print(xtable(getCoverWidth_byVar.func(7)), include.rownames=FALSE)

# tabulate by simulation settings 
print(xtable(getCoverWidth_bySetting.func(2)), include.rownames=FALSE)
print(xtable(getCoverWidth_bySetting.func(6)), include.rownames=FALSE)

# combine tabulation by simulation settings (similar to Das & Lahiri 2019) 
print(xtable(
  cbind( getCoverWidth_bySetting.func(1), getCoverWidth_bySetting.func(2) %>% select(-c(true_beta)) )
  ), include.rownames=FALSE)

#####################################################################################
#                                Simulation Part II: compare                        #
#####################################################################################

load("samelamb.RData")
samelamb = setting1_small
load("difflamb.RData")
difflamb = setting1_small
rm(setting1_small); gc()

p = 10

# color palette
pal2 <- pal_simul[c(2,3,4)] # for RW1,2,3

boxplot.func <- function(RW1_vec, RW2_vec, RW3_vec, nsimul, ylab){
  df <- data.frame(
      value=c(RW1_vec, RW2_vec, RW3_vec),
      method=rep( c("RW1", "RW2", "RW3"), each = nsimul)
    )
  df$method <- factor(df$method, levels = c("RW1", "RW2", "RW3"))

  ggplot(df, aes(x=method, y=value, fill=method)) +
    geom_boxplot(outlier.shape=".", width = .6) +
    labs(x="method", y=ylab) + 
    scale_fill_manual(values = pal2) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
}

#------------------------------------------------------------------------------------
# lambda

ggplot(data.frame(
  value = c(difflamb$lamb_1step, difflamb$lamb_2step),
  method = rep(c("1-step","2-step"), each = nsimul)
), aes(x = value, color = method, fill = method)) +
  geom_histogram(aes(y=..density..), alpha=.5, position = "identity", bins=70) +
  labs(x = "minimum lambda")

#------------------------------------------------------------------------------------
# Total Variation

# png(file="simul2_TV.png", width = 6.3, height = 4.5, units = "in", res = 1200)
setEPS()
postscript("simul2_TV.eps", width = 6.3, height = 4.5)
grid.arrange(
  boxplot.func(
    RW1_vec = samelamb$TV_RW1, 
    RW2_vec = samelamb$TV_RW2, 
    RW3_vec = samelamb$TV_RW3, 
    nsimul = 500, 
    ylab = expr(paste("Total variation (averaged across all ",beta,"'s)"))
  ) + labs(title = expr(paste("Using ",lambda,"'s from 1-step cv")) ) + 
    coord_cartesian(ylim = c(.002, .021)), 
  boxplot.func(
    RW1_vec = difflamb$TV_RW1, 
    RW2_vec = difflamb$TV_RW2, 
    RW3_vec = difflamb$TV_RW3, 
    nsimul = 500, 
    ylab = expr(paste("Total variation (averaged across all ",beta,"'s)"))
  ) + labs(title = expr(paste("Using ",lambda,"'s from 2-step cv")) ) + 
    coord_cartesian(ylim = c(.002, .021)), 
  ncol = 2
)
dev.off()

#------------------------------------------------------------------------------------
# Probability of selecting variable

select_df1 <- rbind(
  melt(samelamb$select_RW1)[,-1],
  melt(samelamb$select_RW2)[,-1],
  melt(samelamb$select_RW3)[,-1]
)
colnames(select_df1) <- c("variable", "value")
select_df1$method = rep( c("RW1", "RW2", "RW3"), each = nsimul*p)
select_df1$method <- factor(select_df1$method, levels = c("RW1", "RW2", "RW3"))
select_df1$variable <- paste0("Variable_", select_df1$variable)
select_df1$variable <- factor(select_df1$variable, levels = paste0("Variable_",1:10))

select_df2 <- rbind(
  melt(difflamb$select_RW1)[,-1],
  melt(difflamb$select_RW2)[,-1],
  melt(difflamb$select_RW3)[,-1]
)
colnames(select_df2) <- c("variable", "value")
select_df2$method = rep( c("RW1", "RW2", "RW3"), each = nsimul*p)
select_df2$method <- factor(select_df2$method, levels = c("RW1", "RW2", "RW3"))
select_df2$variable <- paste0("Variable_", select_df2$variable)
select_df2$variable <- factor(select_df2$variable, levels = paste0("Variable_",1:10))

# png(file="simul2_prob_select.png", width=8, height=9, units = "in", res = 1500)
setEPS()
postscript("simul2_prob_select.eps", width = 8, height = 9)
grid.arrange(
  ggplot(select_df1, aes(x=method, y=value, fill=method)) +
    geom_boxplot(outlier.shape=".", width = .5) +
    facet_wrap(~ variable,ncol = 2, scales = "free")+
    scale_fill_manual(values = pal2) +
    coord_cartesian(ylim = c(0,1)) + 
    labs(x="method", y="probability of selecting variable",
         title = expr(paste("Using ",lambda,"'s from 1-step cv"))) + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom"),
  ggplot(select_df2, aes(x=method, y=value, fill=method)) +
    geom_boxplot(outlier.shape=".", width = .5) +
    facet_wrap(~ variable,ncol = 2, scales = "free")+
    scale_fill_manual(values = pal2) +
    coord_cartesian(ylim = c(0,1)) + 
    labs(x="method", y="probability of selecting variable",
         title = expr(paste("Using ",lambda,"'s from 2-step cv"))) + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom"),
  ncol = 2
)
dev.off()

#------------------------------------------------------------------------------------

cbind(
  paste0("Variable_", 1:10),
  round((difflamb$Cover_RW1)/(samelamb$Cover_RW1),2),
  round((difflamb$Cover_RW2)/(samelamb$Cover_RW2),2),
  round((difflamb$Cover_RW3)/(samelamb$Cover_RW3),2),
  round((difflamb$CIwidth_RW1)/(samelamb$CIwidth_RW1),2),
  round((difflamb$CIwidth_RW2)/(samelamb$CIwidth_RW2),2),
  round((difflamb$CIwidth_RW3)/(samelamb$CIwidth_RW3),2)
) %>%
  data.frame() %>%
  setNames(., c("Variable", "Ratio_Cover_RW1", "Ratio_Cover_RW2", "Ratio_Cover_RW3",
                "Ratio_CIwidth_RW1", "Ratio_CIwidth_RW2", "Ratio_CIwidth_RW3")) 

#------------------------------------------------------------------------------------

# MSE
grid.arrange(
  boxplot.func(
    RW1_vec = samelamb$MSE_RW1,
    RW2_vec = samelamb$MSE_RW2,
    RW3_vec = samelamb$MSE_RW3,
    nsimul = 500,
    ylab = "MSE"
  ) + labs(title = expr(paste("Using ",lambda,"'s from 1-step cv")) ),
  boxplot.func(
    RW1_vec = difflamb$MSE_RW1,
    RW2_vec = difflamb$MSE_RW2,
    RW3_vec = difflamb$MSE_RW3,
    nsimul = 500,
    ylab = "MSE"
  ) + labs(title = expr(paste("Using ",lambda,"'s from 2-step cv")) ),
  ncol = 2
)


# MSPE
grid.arrange(
  boxplot.func(
    RW1_vec = samelamb$MSPE_RW1,
    RW2_vec = samelamb$MSPE_RW2,
    RW3_vec = samelamb$MSPE_RW3,
    nsimul = 500,
    ylab = "MSPE"
  ) + labs(title = expr(paste("Using ",lambda,"'s from 1-step cv")) ) + 
    coord_cartesian(ylim = c(.8, 1.52)),
  boxplot.func(
    RW1_vec = difflamb$MSPE_RW1,
    RW2_vec = difflamb$MSPE_RW2,
    RW3_vec = difflamb$MSPE_RW3,
    nsimul = 500,
    ylab = "MSPE"
  ) + labs(title = expr(paste("Using ",lambda,"'s from 2-step cv")) ) + 
    coord_cartesian(ylim = c(.8, 1.52)),
  ncol = 2
)
