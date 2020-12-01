library(MASS)
library(GGally)
library(gridExtra)
library(monomvn)
library(glmnet)
library(reshape2)

data("Boston")

# check data
lapply(Boston, class)
sum(is.na(Boston))
sum(duplicated(Boston))

# Exploratory data analysis
ggcorr(Boston)
summary(Boston)

png(file="ggpairs.png", width=1500, height=1500)
ggpairs(Boston)
dev.off()

#-----------------------------------------------------------------------------
# some user-defined functions 

# ggplot function for density plot
dens_plot.func <- function(MCMC_vec, RW1_vec, RW2_vec, RW3_vec, xlab){
  df <- data.frame(
    method = rep(c("MCMC", "RW_scheme1", "RW_scheme2", "RW_scheme3"), 
                 times = c(length(MCMC_vec), length(RW1_vec), length(RW2_vec), length(RW3_vec))),
    value = c(MCMC_vec, RW1_vec, RW2_vec, RW3_vec)
  )
  ggplot(df, aes(x=value, color=method)) +
    geom_density() + 
    labs(x = xlab) +
    scale_color_manual(values = c("black", "red", "blue","green"))
}

# ggplot function for density histogram
hist_plot.func <- function(MCMC_vec, RW1_vec, xlab){
  df <- data.frame(
    method = rep(c("MCMC", "RW_scheme1"), 
                 times = c(length(MCMC_vec), length(RW1_vec))),
    value = c(MCMC_vec, RW1_vec)
  )
  ggplot(df, aes(x=value, color=method, fill=method)) +
    geom_histogram(aes(y=..density..), alpha=0.5, position="identity", bins = 70) + 
    labs(x = xlab)
}

#-----------------------------------------------------------------------------
# First, let's use all data for training/validation. No test set

Y = Boston$medv
X = subset(Boston, select = -c(medv))
Y = scale(Y)
X = scale(X)

n = nrow(Boston)
p = ncol(X)
var_name = colnames(X)

B = 3000
burn_in = 2000

SSTO <- sum( (Y - mean(Y))^2 )
eigen(crossprod(X))$values # check it's invertible

#-----------------------------------------------------------------------------
# Bayesian inference

blasso_fit <- blasso(X = X, y = Y, T = B + burn_in)
beta_MCMC <- blasso_fit$beta[(burn_in + 1) : (B + burn_in),]
colnames(beta_MCMC) <- var_name

# probability of variables being selected
select_MCMC <- apply(beta_MCMC, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})

# check trace plots of betas
png(file="beta_traceplots.png", width=1000, height=2500)
par(mfrow=c(7,2))
for(j in 1:p){
  plot(1:B, beta_MCMC[,j], type = "b", pch = ".", xlab = "iteration",  ylab = "beta",
       main = paste0("Trace plots for coeff of variable ", colnames(beta_MCMC)[j]))
}
dev.off()

# check traceplot of (adjusted) coefficient of determination

SSE_MCMC <- apply(beta_MCMC, 1, function(b){
  sum( (Y - X %*% b)^2 )
})
R2_MCMC <- 1 - SSE_MCMC/SSTO
R2_adj_MCMC <- 1 - SSE_MCMC/SSTO*(n-1)/(n-1-p-1)

png(file="R2_traceplots.png", width=600, height=600)
par(mfrow=c(2,1))
plot(1:B, R2_MCMC, type = "b", pch = ".", xlab = "iteration",  ylab = "R2",
     main = paste0("Trace plots for coefficient of determination R2"))
plot(1:B, R2_adj_MCMC, type = "b", pch = ".", xlab = "iteration",  ylab = "adjusted R2",
     main = paste0("Trace plots for adjusted R2"))
dev.off()

#-----------------------------------------------------------------------------
# random-weighting

# first get lambda
res0 = cv.glmnet(x = X, y = Y); plot(res0)
lamb = res0$lambda.min
# lamb = 0.01

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
  
  if(b %% 100 == 0){
    print(b)
  }
}

# probability of variables being selected
select_RW1 <- apply(beta_RW1, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})
select_RW2 <- apply(beta_RW2, 2, function(x){
  1 - sum(as.numeric(x == 0))/B})
select_RW3 <- apply(beta_RW3, 2, function(x){
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

# (adjusted) R^2
R2_RW1 <- 1 - SSE_RW1/SSTO
R2_RW2 <- 1 - SSE_RW2/SSTO
R2_RW3 <- 1 - SSE_RW3/SSTO
R2_adj_RW1 <- 1 - SSE_RW1/SSTO*(n-1)/(n-1-p-1)
R2_adj_RW2 <- 1 - SSE_RW2/SSTO*(n-1)/(n-1-p-1)
R2_adj_RW3 <- 1 - SSE_RW3/SSTO*(n-1)/(n-1-p-1)

#-----------------------------------------------------------------------------
# plot marginal posterior distribution 

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

png(file="marginal_posterior_beta.png", width=1800, height=2500)
ggplot(df, aes(beta,colour = method)) +
  geom_density(alpha = 1) + 
  facet_wrap(~ variable,scales = "free",ncol = 2)+
  theme_bw() +
  scale_color_manual(values = c("black", "red", "blue","green"))
dev.off()

#-----------------------------------------------------------------------------
# compare probability of variables selected  

prob_iii <- match( names( sort((select_MCMC+select_RW1+select_RW2)/3) ),
                  names( (select_MCMC+select_RW1+select_RW2)/3 ) )

png(file="prob_var_select.png", width=600, height=500)
plot(1:p, select_MCMC[prob_iii], type = "b", xaxt = "n", pch = 16, 
     ylab = "probability of variable being selected", main = "", xlab = "variables" )
axis(side = 1, at = 1:p, labels = var_name[prob_iii])
lines(1:p, select_RW1[prob_iii], type = "b", lty = 2, pch = 16, col = "red")
lines(1:p, select_RW2[prob_iii], type = "b", lty = 3, pch = 16, col = "blue")
lines(1:p, select_RW3[prob_iii], type = "b", lty = 4, pch = 16, col = "green")
legend('bottomright', bty = "n", lty = c(1:4), col = c("black", "red", "blue","green"),
       c("MCMC", "RW scheme 1", "RW scheme 2","RW scheme 3"))
dev.off()

#-----------------------------------------------------------------------------
# distribution of (adjusted) R^2

png(file="R2.png", width=750, height=1500)
grid.arrange(
  dens_plot.func(
    MCMC_vec = R2_MCMC, 
    RW1_vec = R2_RW1, 
    RW2_vec = R2_RW2, 
    RW3_vec = R2_RW3,
    xlab = "Coefficient of determination R2"
  ) + xlim(.65,.75),
  dens_plot.func(
    MCMC_vec = R2_adj_MCMC, 
    RW1_vec = R2_adj_RW1, 
    RW2_vec = R2_adj_RW2, 
    RW3_vec = R2_adj_RW3,
    xlab = "Adjusted R2"
  ) + xlim(.65,.74),
  hist_plot.func(
    MCMC_vec = R2_adj_MCMC, 
    RW1_vec = R2_adj_RW1, 
    xlab = "Adjusted R2"
  )  + scale_fill_manual(values = c("black", "red")) + 
    scale_color_manual(values = c("black", "red")),
  nrow = 3
)
dev.off()

