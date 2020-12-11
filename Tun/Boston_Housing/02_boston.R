library(MASS)
library(gridExtra)
library(monomvn)
library(glmnet)
library(reshape2)
library(ggplot2)

data("Boston")

###############################################################################
#                               some user-defined functions 
###############################################################################

# plot marginal posterior density of beta
pos_dens_beta <- function(beta_MCMC, beta_RW1, beta_RW2, beta_RW3, title){
  df <- rbind(
    melt(beta_RW1)[,-1],
    melt(beta_RW2)[,-1],
    melt(beta_RW3)[,-1],
    melt(beta_MCMC)[,-1]
  )
  colnames(df) <- c("variable", "beta")
  method = c(
    rep("RW_scheme1",B*p),
    rep("RW_scheme2",B*p),
    rep("RW_scheme3",B*p),
    rep("MCMC",B*p)
  )
  df <- cbind(method, df)
  
  ggplot(df, aes(beta,colour = method)) +
    geom_density(alpha = 1) + 
    facet_wrap(~ variable,scales = "free",ncol = 2)+
    theme_bw() + 
    scale_color_manual(values = c("black", "red", "blue","green")) +
    labs(title = title) + 
    theme(plot.title = element_text(hjust = 0.5))
}

# plot marginal posterior cdf of beta
pos_cdf_beta <- function(beta_MCMC, beta_RW1, beta_RW2, beta_RW3, title){
  df <- rbind(
    melt(beta_RW1)[,-1],
    melt(beta_RW2)[,-1],
    melt(beta_RW3)[,-1],
    melt(beta_MCMC)[,-1]
  )
  colnames(df) <- c("variable", "beta")
  method = c(
    rep("RW_scheme1",B*p),
    rep("RW_scheme2",B*p),
    rep("RW_scheme3",B*p),
    rep("MCMC",B*p)
  )
  df <- cbind(method, df)
  
  ggplot(df, aes(beta,colour = method)) +
    stat_ecdf(geom = "step") + 
    facet_wrap(~ variable,scales = "free",ncol = 2)+
    theme_bw() + 
    scale_color_manual(values = c("black", "red", "blue","green")) +
    labs(title = title, y = "ecdf") + 
    theme(plot.title = element_text(hjust = 0.5))
}

# plot probability of variables being selected
plot_prob_var_select <- function(select_MCMC, select_RW1, select_RW2, select_RW3){
  prob_iii <- match( names( sort(select_MCMC, decreasing = T) ),
                     names( select_MCMC ) )
  
  plot(1:p, select_MCMC[prob_iii], type = "b", xaxt = "n", pch = 16, 
       ylab = "probability of variable being selected", main = "", xlab = "variables" )
  axis(side = 1, at = 1:p, labels = var_name[prob_iii])
  lines(1:p, select_RW1[prob_iii], type = "b", lty = 2, pch = 16, col = "red")
  lines(1:p, select_RW2[prob_iii], type = "b", lty = 3, pch = 16, col = "blue")
  lines(1:p, select_RW3[prob_iii], type = "b", lty = 4, pch = 16, col = "green")
  legend('bottomleft', bty = "n", lty = c(1:4), col = c("black", "red", "blue","green"),
         c("MCMC", "RW scheme 1", "RW scheme 2","RW scheme 3"))
}

# ggplot function for density or ecdf plot
dens_plot.func <- function(MCMC_vec, RW1_vec, RW2_vec, RW3_vec, xlab, title, plot_type){
  df <- data.frame(
    method = rep(c("MCMC", "RW_scheme1", "RW_scheme2", "RW_scheme3"), 
                 times = c(length(MCMC_vec), length(RW1_vec), length(RW2_vec), length(RW3_vec))),
    value = c(MCMC_vec, RW1_vec, RW2_vec, RW3_vec)
  )
  if(plot_type == "density"){
    ggplot(df, aes(x=value, color=method)) +
      geom_density() + 
      labs(x = xlab, title = title) +
      scale_color_manual(values = c("black", "red", "blue","green")) +
      theme(plot.title = element_text(hjust = 0.5))
  } else if(plot_type == "ecdf"){
    ggplot(df, aes(x=value, color=method)) +
      stat_ecdf(geom = "step") + 
      labs(x = xlab, title = title, y = "ecdf") +
      scale_color_manual(values = c("black", "red", "blue","green")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
}


###############################################################################
#                        split training/validation and test set 
###############################################################################

set.seed(1)
test_samp <- base::sample(1: nrow(Boston), size = ceiling(nrow(Boston)/5))
train_samp <- which( !(c(1: nrow(Boston)) %in% test_samp) )

Y = Boston$medv
X = subset(Boston, select = -c(medv))

Y_train = Y[train_samp]
X_train = X[train_samp,]
Y_test = Y[test_samp]
X_test = X[test_samp,]

# scale predictors and responses
Y_train = scale(Y_train)
X_train = scale(X_train)
Y_test = scale(Y_test)
X_test = scale(X_test)

n_train = nrow(X_train)
n_test = nrow(X_test)
p = ncol(X)
var_name = colnames(X)

rm(Y,X); gc()

B = 2000
burn_in = 2000

SSTO_train <- sum( (Y_train - mean(Y_train))^2 )
# check 
eigen(crossprod(X_train))$values # check it's invertible
eigen(crossprod(X_test))$values # check it's invertible

###############################################################################
#                                    Bayesian inference
###############################################################################

set.seed(1)
blasso_fit <- blasso(X = X_train, y = Y_train, T = B + burn_in)
beta_MCMC <- blasso_fit$beta[(burn_in + 1) : (B + burn_in),]
colnames(beta_MCMC) <- var_name

# probability of variables being selected
select_MCMC <- apply(beta_MCMC, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})

# SSE and R^2
SSE_MCMC <- apply(beta_MCMC, 1, function(b){
  sum( (Y_train - X_train %*% b)^2 )
})
R2_MCMC <- 1 - SSE_MCMC/SSTO_train

# MSPE
MSPE_MCMC <- apply(beta_MCMC, 1, function(b){
  mean( (Y_test - X_test %*% b)^2 )
})

#------------------------------------------------------------------------------
# Assess MCMC 

# check trace plots of betas
par(mfrow=c(7,2))
for(j in 1:p){
  plot(1:B, beta_MCMC[,j], type = "b", pch = ".", xlab = "iteration",  ylab = "beta",
       main = paste0("Trace plots for coeff of variable ", colnames(beta_MCMC)[j]))
}

# check traceplot of coefficient of determination
plot(1:B, R2_MCMC, type = "b", pch = ".", xlab = "iteration",  ylab = "R2",
     main = paste0("Trace plots for coefficient of determination R^2"))

###############################################################################
#                         random-weighting (one-step)
###############################################################################

# one-step

# # first get lambda
# res0 = cv.glmnet(x = X_train, y = Y_train); plot(res0)
# lamb = res0$lambda.min
# 
# beta_RW1 <- matrix(NA, nrow = B, ncol = p) # no penalty weight
# beta_RW2 <- matrix(NA, nrow = B, ncol = p) # common penalty weight
# beta_RW3 <- matrix(NA, nrow = B, ncol = p) # different penalty weights
# colnames(beta_RW1) <- var_name
# colnames(beta_RW2) <- var_name
# colnames(beta_RW3) <- var_name
# 
# set.seed(1)
# 
# # 1-step
# for(b in 1:B){
#   w = rexp(n_train)
#   wp = rexp(p)
#   RW1_fit = glmnet(x = X_train, y = Y_train, weights = w)
#   RW3_fit = glmnet(x = X_train, y = Y_train, weights = w, penalty.factor = wp)
#   
#   beta_RW1[b,] <- as.vector(coef(RW1_fit, s = lamb))[-1]
#   beta_RW2[b,] <- as.vector(coef(RW1_fit, s = lamb*wp[1]))[-1]
#   beta_RW3[b,] <- as.vector(coef(RW3_fit, s = lamb*mean(wp)))[-1]
#   
#   if(b %% 100 == 0){
#     print(b)
#   }
# }
# 
# save(lamb, beta_RW1, beta_RW2, beta_RW3, file = "Boston_RW_onestep.RData")
load("Boston_RW_onestep.RData")

###############################################################################
#                         random-weighting (two-step)
###############################################################################

# my own cv step

# grid for log-lambdas
loglamb_vec <- seq(-7, -.5, by = .05) 

# initiate cv MSE
cv_MSE <- rep(NA, length(loglamb_vec))

# K-fold cv
cv_K <- 4

# cv partition
part_count <- rep(floor(n_train/cv_K), cv_K)
part_count[cv_K] <- n_train - sum(head(part_count,-1)) 
cv_part <- rep(1:cv_K, times = part_count)
cv_part <- sample(cv_part, size = n_train, replace = F)

# compute CV MSE  
for(l in 1:length(loglamb_vec)){
  loglamb <- loglamb_vec[l]
  lamb <- exp(loglamb)
  cv_SSE <- 0
  
  for(k in 1:cv_K){
    Y_cvTrain <- scale(Y_train[(cv_part != k)])
    X_cvTrain <- scale(X_train[(cv_part != k),])
    Y_cvTest <- scale(Y_train[(cv_part == k)])
    X_cvTest <- scale(X_train[(cv_part == k),])
    
    cvfit = glmnet(x = X_cvTrain, y = Y_cvTrain)
    beta_cvfit = as.vector(coef(cvfit, s = lamb))[-1]
    if( sum(abs(beta_cvfit)) != 0 ){
      nonzero_cvfit <- as.vector( which(beta_cvfit != 0) )
      beta_cvfit[nonzero_cvfit] <- as.vector(coef(lm(Y_cvTrain ~ X_cvTrain[, nonzero_cvfit] - 1)))
    }
    cv_SSE = cv_SSE + sum( (Y_cvTest -  X_cvTest %*%beta_cvfit)^2 )
  }
  
  cv_MSE[l] <- cv_SSE/n_train
  print(l)
}

# plot cv MSE
plot( loglamb_vec, cv_MSE, type = "b", pch = 16, xlab = expr(paste(log(lambda))), ylab =  "CV MSE")

# get lamb.min
lamb.min <- exp( loglamb_vec[which.min(cv_MSE)] ); lamb.min

#------------------------------------------------------------------------------
# 2-step

# lamb <- lamb.min
# 
# beta_RW1 <- matrix(NA, nrow = B, ncol = p) # no penalty weight
# beta_RW2 <- matrix(NA, nrow = B, ncol = p) # common penalty weight
# beta_RW3 <- matrix(NA, nrow = B, ncol = p) # different penalty weights
# colnames(beta_RW1) <- var_name
# colnames(beta_RW2) <- var_name
# colnames(beta_RW3) <- var_name
# 
# set.seed(1)
# for(b in 1:B){
#   w = rexp(n_train)
#   wp = rexp(p)
#   RW1_fit = glmnet(x = X_train, y = Y_train, weights = w)
#   RW3_fit = glmnet(x = X_train, y = Y_train, weights = w, penalty.factor = wp)
#   
#   beta_RW1[b,] <- as.vector(coef(RW1_fit, s = lamb))[-1]
#   beta_RW2[b,] <- as.vector(coef(RW1_fit, s = lamb*wp[1]))[-1]
#   beta_RW3[b,] <- as.vector(coef(RW3_fit, s = lamb*mean(wp)))[-1]
#   
#   nonzero_RW1 <- as.vector( which(beta_RW1[b,] != 0) )
#   nonzero_RW2 <- as.vector( which(beta_RW2[b,] != 0) )
#   nonzero_RW3 <- as.vector( which(beta_RW3[b,] != 0) ) 
#   
#   beta_RW1[b, nonzero_RW1] <- as.vector(coef(lm(Y_train ~ X_train[, nonzero_RW1] - 1, weights = w)))
#   beta_RW2[b, nonzero_RW2] <- as.vector(coef(lm(Y_train ~ X_train[, nonzero_RW2] - 1, weights = w)))
#   beta_RW3[b, nonzero_RW3] <- as.vector(coef(lm(Y_train ~ X_train[, nonzero_RW3] - 1, weights = w)))
#   
#   if(b %% 100 == 0){
#     print(b)
#   }
# }
# 
# save(lamb, beta_RW1, beta_RW2, beta_RW3, file = "Boston_RW_twostep.RData")
load("Boston_RW_twostep.RData")

###############################################################################
#                       random-weighting (metric to compare)
###############################################################################

# probability of variables being selected
select_RW1 <- apply(beta_RW1, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})
select_RW2 <- apply(beta_RW2, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})
select_RW3 <- apply(beta_RW3, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})

# SSE
SSE_RW1 <- apply(beta_RW1 , 1, function(b){
  sum( (Y_train - X_train %*% b)^2 )
})
SSE_RW2 <- apply(beta_RW2 , 1, function(b){
  sum( (Y_train - X_train %*% b)^2 )
})
SSE_RW3 <- apply(beta_RW3 , 1, function(b){
  sum( (Y_train - X_train %*% b)^2 )
})

# R^2
R2_RW1 <- 1 - SSE_RW1/SSTO_train
R2_RW2 <- 1 - SSE_RW2/SSTO_train
R2_RW3 <- 1 - SSE_RW3/SSTO_train

# MSPE
MSPE_RW1 <- apply(beta_RW1, 1, function(b){
  mean( (Y_test - X_test %*% b)^2 )
})
MSPE_RW2 <- apply(beta_RW2, 1, function(b){
  mean( (Y_test - X_test %*% b)^2 )
})
MSPE_RW3 <- apply(beta_RW3, 1, function(b){
  mean( (Y_test - X_test %*% b)^2 )
})


###############################################################################
#                                       compare betas!
###############################################################################

# compare marginal posterior density of beta
pos_dens_beta(
  beta_MCMC = beta_MCMC, 
  beta_RW1 = beta_RW1, 
  beta_RW2 = beta_RW2, 
  beta_RW3 = beta_RW3, 
  title = "Marginal Posterior Density of Beta"
)

# compare marginal posterior cdf of beta
pos_cdf_beta(
  beta_MCMC = beta_MCMC, 
  beta_RW1 = beta_RW1, 
  beta_RW2 = beta_RW2, 
  beta_RW3 = beta_RW3, 
  title = "Marginal Posterior Distribution of Beta"
)

# probability of variables being selected
plot_prob_var_select(
  select_MCMC = select_MCMC, 
  select_RW1 = select_RW1, 
  select_RW2 = select_RW2, 
  select_RW3 = select_RW3
)

# compare density of R^2
dens_plot.func(
  MCMC_vec = R2_MCMC, 
  RW1_vec = R2_RW1, 
  RW2_vec = R2_RW2, 
  RW3_vec = R2_RW3, 
  xlab = "Coefficient of determination R^2", 
  title = "Posterior density of R^2", 
  plot_type = "density"
) + xlim(.68,.75)

# compare ecdf of R^2
dens_plot.func(
  MCMC_vec = R2_MCMC, 
  RW1_vec = R2_RW1, 
  RW2_vec = R2_RW2, 
  RW3_vec = R2_RW3, 
  xlab = "Coefficient of determination R^2", 
  title = "Posterior ecdf of R^2", 
  plot_type = "ecdf"
) + xlim(.68,.75)

# compare density of MSPE
dens_plot.func(
  MCMC_vec = MSPE_MCMC, 
  RW1_vec = MSPE_RW1, 
  RW2_vec = MSPE_RW2, 
  RW3_vec = MSPE_RW3, 
  xlab = "MSPE", 
  title = "Density of MSPE", 
  plot_type = "density"
) + xlim(.22,.38)

# compare ecdf of MSPE
dens_plot.func(
  MCMC_vec = MSPE_MCMC, 
  RW1_vec = MSPE_RW1, 
  RW2_vec = MSPE_RW2, 
  RW3_vec = MSPE_RW3, 
  xlab = "MSPE", 
  title = "ecdf of MSPE", 
  plot_type = "ecdf"
) + xlim(.22,.35)
