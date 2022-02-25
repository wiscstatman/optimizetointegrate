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

# color palette
# correspond to MCMC, RW1, RW2, RW3, RB
pal_boston <- c("black", "red", "blue","green", "orange")

# 2-step CV function
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
      
      cvfit <- glmnet(x = X_cvTrain, y = Y_cvTrain, intercept = T)
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

# plot marginal posterior density of beta
pos_dens_beta <- function(beta_MCMC, beta_RW1, beta_RW2, beta_RW3, beta_RB, title){
  df <- rbind(
    melt(beta_RW1)[,-1],
    melt(beta_RW2)[,-1],
    melt(beta_RW3)[,-1],
    melt(beta_MCMC)[,-1],
    melt(beta_RB)[,-1]
  )
  colnames(df) <- c("variable", "beta")
  method = c(
    rep("RW1",B*p),
    rep("RW2",B*p),
    rep("RW3",B*p),
    rep("MCMC",B*p),
    rep("RB",B*p)
  )
  df <- cbind(method, df)
  df$method = factor(df$method, levels=c("MCMC", "RW1", "RW2", "RW3", "RB"))
  
  ggplot(df, aes(beta,colour = method)) +
    geom_density(alpha = 1) + 
    facet_wrap(~ variable,scales = "free",ncol = 2)+
    theme_bw() + 
    scale_color_manual(values = pal_boston) +
    labs(title = title) + 
    theme(plot.title = element_text(hjust = 0.5))
}

# plot marginal posterior cdf of beta
pos_cdf_beta <- function(beta_MCMC, beta_RW1, beta_RW2, beta_RW3, beta_RB, title){
  df <- rbind(
    melt(beta_RW1)[,-1],
    melt(beta_RW2)[,-1],
    melt(beta_RW3)[,-1],
    melt(beta_MCMC)[,-1],
    melt(beta_RB)[,-1]
  )
  colnames(df) <- c("variable", "beta")
  method = c(
    rep("RW1",B*p),
    rep("RW2",B*p),
    rep("RW3",B*p),
    rep("MCMC",B*p),
    rep("RB",B*p)
  )
  df <- cbind(method, df)
  df$method = factor(df$method, levels=c("MCMC", "RW1", "RW2", "RW3", "RB"))
  
  ggplot(df, aes(beta,colour = method)) +
    stat_ecdf(geom = "step") + 
    facet_wrap(~ variable,scales = "free",ncol = 2)+
    theme_bw() + 
    scale_color_manual(values = pal_boston) +
    labs(title = title, y = "ecdf") + 
    theme(plot.title = element_text(hjust = 0.5))
}

# plot probability of variables being selected
plot_prob_var_select <- function(select_MCMC, select_RW1, select_RW2, select_RW3, select_RB){
  prob_iii <- match( names( sort(select_MCMC, decreasing = T) ),
                     names( select_MCMC ) )
  
  plot(1:p, select_MCMC[prob_iii], type = "b", xaxt = "n", pch = 16, col = pal_boston[1], 
       ylab = "probability of variable being selected", main = "", xlab = "variables" )
  axis(side = 1, at = 1:p, labels = var_name[prob_iii])
  lines(1:p, select_RW1[prob_iii], type = "b", lty = 2, pch = 16, col = pal_boston[2])
  lines(1:p, select_RW2[prob_iii], type = "b", lty = 3, pch = 16, col = pal_boston[3])
  lines(1:p, select_RW3[prob_iii], type = "b", lty = 4, pch = 16, col = pal_boston[4])
  lines(1:p, select_RB[prob_iii], type = "b", lty = 5, pch = 16, col = pal_boston[5])
  legend('bottomleft', bty = "n", lty = c(1:4), col = pal_boston,
         c("MCMC", "RW1", "RW2","RW3", "RB"))
}

# ggplot function for density or ecdf plot
dens_plot.func <- function(MCMC_vec, RW1_vec, RW2_vec, RW3_vec, RB_vec, xlab, title, plot_type){
  df <- data.frame(
    method = rep(c("MCMC", "RW1", "RW2", "RW3", "RB"), 
                 times = c(length(MCMC_vec), length(RW1_vec), length(RW2_vec), 
                           length(RW3_vec), length(RB_vec))),
    value = c(MCMC_vec, RW1_vec, RW2_vec, RW3_vec, RB_vec)
  )
  df$method = factor(df$method, levels=c("MCMC", "RW1", "RW2", "RW3", "RB"))
  if(plot_type == "density"){
    ggplot(df, aes(x=value, color=method)) +
      geom_density() + 
      labs(x = xlab, title = title) +
      scale_color_manual(values = pal_boston) +
      theme(plot.title = element_text(hjust = 0.5))
  } else if(plot_type == "ecdf"){
    ggplot(df, aes(x=value, color=method)) +
      stat_ecdf(geom = "step") + 
      labs(x = xlab, title = title, y = "ecdf") +
      scale_color_manual(values = pal_boston) +
      theme(plot.title = element_text(hjust = 0.5))
  }
}

###############################################################################
#                                      Data                                   #
###############################################################################

Y = scale(Boston$medv)
X = scale(subset(Boston, select = -c(medv)))

n = nrow(X)
p = ncol(X)
var_name = colnames(X)
B = 1000
burn_in = 2000

SSTO_train <- sum( (Y - mean(Y))^2 )
eigen(crossprod(X))$values # check it's invertible


###############################################################################
#                                    Bayesian inference
###############################################################################

set.seed(1)
blasso_fit <- blasso(X = X, y = Y, T = B + burn_in)
beta_MCMC <- blasso_fit$beta[(burn_in + 1) : (B + burn_in),]
colnames(beta_MCMC) <- var_name

# probability of variables being selected
select_MCMC <- apply(beta_MCMC, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})

# SSE and R^2
SSE_MCMC <- apply(beta_MCMC, 1, function(b){
  sum( (Y - X %*% b)^2 )
})
R2_MCMC <- 1 - SSE_MCMC/SSTO_train
MSE_MCMC <- SSE_MCMC/n

#------------------------------------------------------------------------------
# Assess MCMC 

# check trace plots of betas
par(mfrow=c(7,2))
for(j in 1:p){
  plot(1:B, beta_MCMC[,j], type = "b", pch = ".", xlab = "iteration",  ylab = "beta",
       main = paste0("Trace plots for coeff of variable ", colnames(beta_MCMC)[j]))
}
par(mfrow=c(1,1))

# check traceplot of coefficient of determination
plot(1:B, R2_MCMC, type = "b", pch = ".", xlab = "iteration",  ylab = "R2",
     main = paste0("Trace plots for coefficient of determination R^2"))

###############################################################################
#                         random-weighting (two-step)
###############################################################################

lamb <- CV2step.func(Y = Y, X = X, lambda_path = exp(seq(-7, -.5, by = .05)) , cv_K = 4)

beta_RW1 <- matrix(NA, nrow = B, ncol = p) # no penalty weight
beta_RW2 <- matrix(NA, nrow = B, ncol = p) # common penalty weight
beta_RW3 <- matrix(NA, nrow = B, ncol = p) # different penalty weights
colnames(beta_RW1) <- var_name
colnames(beta_RW2) <- var_name
colnames(beta_RW3) <- var_name

for(b in 1:B){
  w = rexp(n)
  wp = rexp(p)
  RW1_fit = glmnet(x = X, y = Y, weights = w)
  RW3_fit = glmnet(x = X, y = Y, weights = w, penalty.factor = wp)

  beta_RW1[b,] <- as.vector(coef(RW1_fit, s = lamb))[-1]
  beta_RW2[b,] <- as.vector(coef(RW1_fit, s = lamb*wp[1]))[-1]
  beta_RW3[b,] <- as.vector(coef(RW3_fit, s = lamb*mean(wp)))[-1]

  nonzero_RW1 <- as.vector( which(beta_RW1[b,] != 0) )
  nonzero_RW2 <- as.vector( which(beta_RW2[b,] != 0) )
  nonzero_RW3 <- as.vector( which(beta_RW3[b,] != 0) )

  beta_RW1[b, nonzero_RW1] <- as.vector(coef(lm(Y ~ X[, nonzero_RW1] - 1, weights = w)))
  beta_RW2[b, nonzero_RW2] <- as.vector(coef(lm(Y ~ X[, nonzero_RW2] - 1, weights = w)))
  beta_RW3[b, nonzero_RW3] <- as.vector(coef(lm(Y ~ X[, nonzero_RW3] - 1, weights = w)))

  if(b %% 100 == 0){
    print(b)
  }
}

###############################################################################
#                               residual bootstrap (RB)
###############################################################################

res0 <- cv.glmnet(x = X, y = Y, intercept = T)
beta_res0 <-  as.vector( coef(res0, s = res0$lambda.min)[-1] )
resid_res0 <- Y - X %*% beta_res0
resid_res0 <- resid_res0 - mean(resid_res0)
Xbeta_res0 <- X %*% beta_res0 
lamb_1step <- res0$lambda.min

# cache
beta_RB <- matrix(NA, nrow = B, ncol = p) 
colnames(beta_RB) <- var_name

for(b in 1:B){
  resid_boot <- base::sample(x = resid_res0, size = n, replace = T)
  Y_boot <- Xbeta_res0 + resid_boot
  fit_boot <- glmnet::glmnet(x = X, y = Y_boot, intercept = T)
  beta_RB[b,] <- as.vector( coef(fit_boot, s = res0$lambda.min)[-1] )
  
  if(b %% 100 == 0){
    print(b)
  }
}

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
select_RB <- apply(beta_RB, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})

# SSE
SSE_RW1 <- apply(beta_RW1 , 1, function(b){
  sum( (Y - X %*% b)^2 )
})
SSE_RW2 <- apply(beta_RW2 , 1, function(b){
  sum( (Y - X %*% b)^2 )
})
SSE_RW3 <- apply(beta_RW3 , 1, function(b){
  sum( (Y - X %*% b)^2 )
})
SSE_RB <- apply(beta_RB , 1, function(b){
  sum( (Y - X %*% b)^2 )
})


# R^2
R2_RW1 <- 1 - SSE_RW1/SSTO_train
R2_RW2 <- 1 - SSE_RW2/SSTO_train
R2_RW3 <- 1 - SSE_RW3/SSTO_train
R2_RB <- 1 - SSE_RB/SSTO_train

# MSE
MSE_RW1 <- SSE_RW1/n
MSE_RW2 <- SSE_RW2/n
MSE_RW3 <- SSE_RW3/n
MSE_RB <- SSE_RB/n

###############################################################################
#                                       compare betas!
###############################################################################

# rename variable
colnames(beta_MCMC)[which(colnames(beta_MCMC)=="black")] <- "Black"
colnames(beta_RW1)[which(colnames(beta_RW1)=="black")] <- "Black"
colnames(beta_RW2)[which(colnames(beta_RW2)=="black")] <- "Black"
colnames(beta_RW3)[which(colnames(beta_RW3)=="black")] <- "Black"
colnames(beta_RB)[which(colnames(beta_RB)=="black")] <- "Black"

# compare marginal posterior density of beta
pos_dens_beta(
  beta_MCMC = beta_MCMC, 
  beta_RW1 = beta_RW1, 
  beta_RW2 = beta_RW2, 
  beta_RW3 = beta_RW3, 
  beta_RB = beta_RB, 
  title = "Marginal Posterior Density of Beta"
)

# compare marginal posterior cdf of beta
# png("boston_beta_cdf.png", width = 6.5, height = 10, units = "in", res = 1200)
setEPS()
postscript("boston_beta_cdf.eps", width = 6.5, height = 10)
pos_cdf_beta(
  beta_MCMC = beta_MCMC, 
  beta_RW1 = beta_RW1, 
  beta_RW2 = beta_RW2, 
  beta_RW3 = beta_RW3, 
  beta_RB = beta_RB, 
  title = ""
) + labs(x = expr(paste(beta)), y = expr(paste("Marginal posterior/conditional cdf of ", beta)))
dev.off()

# probability of variables being selected
plot_prob_var_select(
  select_MCMC = select_MCMC, 
  select_RW1 = select_RW1, 
  select_RW2 = select_RW2, 
  select_RW3 = select_RW3,
  select_RB = select_RB
)

# compare density of R^2
dens_plot.func(
  MCMC_vec = R2_MCMC, 
  RW1_vec = R2_RW1, 
  RW2_vec = R2_RW2, 
  RW3_vec = R2_RW3, 
  RB_vec = R2_RB, 
  xlab = "Coefficient of determination R^2", 
  title = "Posterior density of R^2", 
  plot_type = "density"
) + coord_cartesian(xlim=c(.68 , .75))

# compare ecdf of R^2
dens_plot.func(
  MCMC_vec = R2_MCMC, 
  RW1_vec = R2_RW1, 
  RW2_vec = R2_RW2, 
  RW3_vec = R2_RW3, 
  RB_vec = R2_RB, 
  xlab = "Coefficient of determination R^2", 
  title = "Posterior ecdf of R^2", 
  plot_type = "ecdf"
) + coord_cartesian(xlim=c(.68,.75))

# compare density of MSE
dens_plot.func(
  MCMC_vec = MSE_MCMC, 
  RW1_vec = MSE_RW1, 
  RW2_vec = MSE_RW2, 
  RW3_vec = MSE_RW3, 
  RB_vec = MSE_RB, 
  xlab = "MSE", 
  title = "Density of MSE", 
  plot_type = "density"
) + coord_cartesian(xlim=c(.22,.38))

# compare ecdf of MSE
dens_plot.func(
  MCMC_vec = MSE_MCMC, 
  RW1_vec = MSE_RW1, 
  RW2_vec = MSE_RW2, 
  RW3_vec = MSE_RW3, 
  RB_vec = MSE_RB, 
  xlab = "MSE", 
  title = "ecdf of MSE", 
  plot_type = "ecdf"
) + coord_cartesian(xlim=c(.22,.35))


