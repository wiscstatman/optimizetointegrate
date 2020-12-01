library(glmnet)
library(monomvn) # blasso
library(basad)
library(tidyverse)
library(parallel)
library(doSNOW) # Windows

####################################################################################### 
#                                      Setup Dataset                                  #
####################################################################################### 

# read dataset from csv
train_input <- read_csv("pria_rmi_cv/file_0.csv")
for(i in 1:4){
  filename <- paste0("pria_rmi_cv/file_",i,".csv")
  train_input <- rbind( train_input, read_csv(file = filename) )
}
test_input <- read_csv("pria_prospective.csv")

# parameters
p <- 1024 # dim size
B <- 2E3 # number of MCMC and RW samples 
n_train <- nrow(train_input)
n_test <- nrow(test_input)
burn_in <- 1E3 # number of burn-in for MCMC

# function to process dataset
list_func <- function(input_df){
  # convert character string to columns of predictors and standardize column-wise
  X <- t( sapply( input_df$Fingerprints, function(x){as.numeric(str_split_fixed(x,"",p))} ) )
  colnames(X) <- paste0("x", 1:p)
  X <- scale(X)

  # center or standardize Y 
  Y <- as.numeric(scale(input_df$Keck_Pria_Continuous))
  # Y <- input_df$Keck_Pria_Continuous - mean(input_df$Keck_Pria_Continuous)

  # get list
  return(list(
    X = X,
    Y = Y
  ))
}

train_list <- list_func(train_input)
test_list <- list_func(test_input)

# check if any columns of design matrix X has 1 unique value only
# check_unique <- apply(train_list$X, 2, function(x){length(unique(x))})
# table(check_unique) # all good
# sum(as.numeric( colnames(train_list$X) == colnames(test_list$X) )) == ncol(train_list$X) # all good

####################################################################################### 
#                                     Bayesian Inference                              #
####################################################################################### 

# try bayesian Lasso (Park & Casella, 2008)
blasso_start_time <- Sys.time()
beta_blasso <- monomvn::blasso( X = train_list$X, y = train_list$Y, T = B + burn_in)
beta_blasso <- beta_blasso$beta[ (burn_in+1):(B+burn_in) , ]
blasso_end_time <- Sys.time()
blasso_proc_time <- blasso_end_time - blasso_start_time; blasso_proc_time

# try basad: shrinking and diffusing prior (Narisetty & He, 2014)
basad_start_time <- Sys.time()
fit_basad <- basad(x = train_list$X, y = train_list$Y, # alternative = T, 
                   nburn = burn_in, niter = B)
basad_end_time <- Sys.time()
basad_proc_time <- basad_end_time - basad_start_time; basad_proc_time

save(fit_basad, file = "basad.RData")
load("basad.RData")

beta_basad <- fit_basad$allB
predY_basad <- predict.basad(fit_basad, testx = test_list$X)


####################################################################################### 
#                                      Random-Weighting                               #
####################################################################################### 

# easy way out, use cv.glmnet
res0 <- cv.glmnet(x = train_list$X, y = train_list$Y)
lamb <- res0$lambda.min

#--------------------------------------------------------------------------------------
# no parallelization 

# initialize
beta_w1 <- matrix(NA, nrow = B, ncol = p) # no penalty weight
beta_w2 <- matrix(NA, nrow = B, ncol = p) # common penalty weight
beta_w3 <- matrix(NA, nrow = B, ncol = p) # different penalty weights
predY_w1 <- matrix(NA, nrow = B, ncol = n_test) # prediction based on beta_w1
predY_w2 <- matrix(NA, nrow = B, ncol = n_test) # prediction based on beta_w2
predY_w3 <- matrix(NA, nrow = B, ncol = n_test) # prediction based on beta_w3

for(b in 1:B){
  w <- rexp(n_train)
  wp <- rexp(p)
  fit_w1 <- glmnet(x = train_list$X, y = train_list$Y, weights = w)
  fit_w3 <- glmnet(x = train_list$X, y = train_list$Y, weights = w, penalty.factor = wp)
  
  beta_w1[b,] <- as.vector(coef(fit_w1, s = lamb))[-1]
  beta_w2[b,] <- as.vector(coef(fit_w1, s = lamb*wp[1]))[-1]
  beta_w3[b,] <- as.vector(coef(fit_w3, s = lamb*mean(wp)))[-1]
  
  predY_w1[b,] <- test_list$X %*% beta_w1[b,]
  predY_w2[b,] <- test_list$X %*% beta_w2[b,]
  predY_w3[b,] <- test_list$X %*% beta_w3[b,]
  
  if(b %% 100 == 0){
    print(b)
  }
}

save(lamb, beta_w1, beta_w2, beta_w3,
     predY_w1, predY_w2, predY_w3,
     file = "RW_Lasso.RData")

#--------------------------------------------------------------------------------------
# with parallelization 

RW_start_time <- Sys.time()

NumCore <- 8
cl <- makeCluster(NumCore)
registerDoSNOW(cl) 

RW_beta <- foreach(i = 1:B, .combine = "rbind") %dopar%{
  w <- rexp(n_train)
  wp <- rexp(p)
  fit_w1 <- glmnet::glmnet(x = train_list$X, y = train_list$Y, weights = w)
  fit_w3 <- glmnet::glmnet(x = train_list$X, y = train_list$Y, weights = w, penalty.factor = wp)
  return(c(
    as.vector(coef(fit_w1, s = lamb))[-1],
    as.vector(coef(fit_w1, s = lamb*wp[1]))[-1],
    as.vector(coef(fit_w3, s = lamb*mean(wp)))[-1]
  ))
}
stopCluster(cl)

RW_end_time <- Sys.time()
RW_proc_time <- RW_end_time - RW_start_time; RW_proc_time

save(lamb, RW_beta, file = "RW_Lasso.RData")

