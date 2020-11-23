library(ggplot2)
library(reshape2) # melt
library(tidyverse)

# load MCMC results
load("basad_limited.RData")

# load either one of them for RW samples
load("RW_Lasso.RData")

p <- 1024
B <- nrow(RW_beta)

# randomly select 2000 out of 10,000 MCMC samples
pick_sample <- sample(1:1E4, B)

# plot to compare lines
plot.func <- function(MCMC_vec, w1_vec, w2_vec, w3_vec, ylab, type="l"){
  iii <- match(sort(MCMC_vec),MCMC_vec)
  plot(1:length(MCMC_vec), MCMC_vec[iii], type = type, xaxt = "n", ylab = ylab, main = "", 
       xlab = paste0("variables x1 - x", p, " arranged according to MCMC samples' order of ",ylab) )
  lines(1:length(MCMC_vec), w1_vec[iii], type = type, lty = 2, col = "red")
  lines(1:length(MCMC_vec), w2_vec[iii], type = type, lty = 3, col = "blue")
  lines(1:length(MCMC_vec), w3_vec[iii], type = type, lty = 4, col = "green")
  lines(1:length(MCMC_vec), MCMC_vec[iii], type = type,lty = 1, col = "black")
  legend('topleft', lty = c(1:4), col = c("black", "red", "blue","green"),
       c("MCMC", "RW scheme 1", "RW scheme 2","RW scheme 3"))
}

plot.func2 <- function(MCMC_vec, w1_vec, w2_vec, w3_vec, ylab, type="l", ylim){
  iii <- match(sort(MCMC_vec),MCMC_vec)
  plot(1:length(MCMC_vec), MCMC_vec[iii], type = type, xaxt = "n", ylab = ylab, main = "", ylim=ylim, 
       xlab = paste0("test cases arranged according to MCMC samples' order of ",ylab) )
  lines(1:length(MCMC_vec), w1_vec[iii], type = type, lty = 2, col = "red")
  lines(1:length(MCMC_vec), w2_vec[iii], type = type, lty = 3, col = "blue")
  lines(1:length(MCMC_vec), w3_vec[iii], type = type, lty = 4, col = "green")
  lines(1:length(MCMC_vec), MCMC_vec[iii], type = type,lty = 1, col = "black")
  legend('topleft', lty = c(1:4), col = c("black", "red", "blue","green"),
         c("MCMC", "RW scheme 1", "RW scheme 2","RW scheme 3"))
}


#-------------------------------------------------------------------------------
# compare probability of variables selected 

beta_select <- apply(RW_beta,2,function(x){
  1 - sum(as.numeric(x == 0))/B
})

beta_select_w1 <- beta_select[1:p] # no penalty weight
beta_select_w2 <- beta_select[(p+1):(2*p)] # common penalty weight
beta_select_w3 <- beta_select[(2*p+1):(3*p)] # different penalty weights
basad_select <- fit_basad_limited$posteriorZ[-1]  # average of posterior probability P(Z=1)

# plot 
# arranged according to average of posterior probability P(Z=1)
plot.func(
  MCMC_vec = basad_select, 
  w1_vec = beta_select_w1,
  w2_vec = beta_select_w2, 
  w3_vec = beta_select_w3,
  ylab = "probability of variable being selected"
)

#-------------------------------------------------------------------------------
# compare model sizes

# for MCMC samples, set those betas with corresponding P(Z=1) < .5 as zero
MCMC_beta <- fit_basad_limited$allB[pick_sample,-1] * fit_basad_limited$allZ[pick_sample, -1]

size_w1 <- apply(RW_beta[,1:p],1,function(x){
  p - sum(as.numeric(x == 0))
})

size_w2 <- apply(RW_beta[,(p+1):(2*p)],1,function(x){
  p - sum(as.numeric(x == 0))
})

size_w3 <- apply(RW_beta[,(2*p+1):(3*p)],1,function(x){
  p - sum(as.numeric(x == 0))
})

size_MCMC <- apply(MCMC_beta,1,function(x){
  p - sum(as.numeric(x == 0))
})

par(mfrow=c(2,2))
hist(size_w1, breaks = 70)
hist(size_w2, breaks = 70)
hist(size_w3, breaks = 70)
hist(size_MCMC, breaks = 70)
par(mfrow=c(1,1))


#-------------------------------------------------------------------------------
# get 95% CI

CI_w1 <- t( apply(RW_beta[,1:p], 2, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(CI_length = upper - lower)

CI_w2 <- t( apply(RW_beta[,(p+1):(2*p)], 2, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(CI_length = upper - lower)

CI_w3 <- t( apply(RW_beta[,(2*p+1):(3*p)], 2, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(CI_length = upper - lower)

CI_MCMC <- t( apply(MCMC_beta, 2, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(CI_length = upper - lower)

# CI coverage

# compare length of CI
plot.func(
  MCMC_vec = CI_MCMC$CI_length, 
  w1_vec = CI_w1$CI_length,
  w2_vec = CI_w2$CI_length, 
  w3_vec = CI_w3$CI_length,
  ylab = "Length of Confidence Intervals"
)

#-------------------------------------------------------------------------------
# Pick top few variables with highest P(Z=1) and plot marginal posterior distribution 
prob_thresh <- .93

theta_ind <- which(
  (beta_select_w1 >  prob_thresh ) &
  (beta_select_w2 >  prob_thresh ) &
  (beta_select_w3 >  prob_thresh ) &
  (basad_select >  prob_thresh )
)

# top 10 variables?
colnames(RW_beta) <- rep( paste0("x", 1:p), 3)
colnames(MCMC_beta) <- paste0("x", 1:p)

df <- rbind(
    melt(RW_beta[,1:p][,theta_ind])[,-1],
    melt(RW_beta[,(p+1):(2*p)][,theta_ind])[,-1],
    melt(RW_beta[,(2*p+1):(3*p)][,theta_ind])[,-1],
    melt(MCMC_beta[,theta_ind])[,-1]
  )
colnames(df) <- c("variable", "beta")
method = c(
    rep("RW_scheme1",B*length(theta_ind)),
    rep("RW_scheme2",B*length(theta_ind)),
    rep("RW_scheme3",B*length(theta_ind)),
    rep("MCMC",B*length(theta_ind))
  )
df <- cbind(method, df)

ggplot(df, aes(beta,colour = method)) +
  geom_density(alpha = 1) + 
  facet_wrap(~ variable,scales = "free",ncol = 2)+
  theme_bw()


#-------------------------------------------------------------------------------
# predict

test_input <- read_csv("pria_prospective.csv")
test_X <- t( sapply( test_input$Fingerprints, function(x){as.numeric(str_split_fixed(x,"",p))} ) )
colnames(test_X) <- paste0("x", 1:p)
test_X <- scale(test_X)
test_Y <- as.numeric(scale(test_input$Keck_Pria_Continuous))

predY_w1 <- apply( RW_beta[,1:p], 1, function(x){
  test_X %*% x
})
predY_w2 <- apply( RW_beta[,(p+1):(2*p)], 1, function(x){
  test_X %*% x
})
predY_w3 <- apply( RW_beta[,(2*p+1):(3*p)], 1, function(x){
  test_X %*% x
})
predY_MCMC <- apply(MCMC_beta, 1, function(x){
  test_X %*% x
})

# get 95% prediction interval & width
predint_w1 <- t( apply(predY_w1, 1, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(PI_length = upper - lower, 
         cover = ifelse((test_Y >= lower) & (test_Y <= upper), 1, 0))

predint_w2 <- t( apply(predY_w2, 1, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(PI_length = upper - lower, 
         cover = ifelse((test_Y >= lower) & (test_Y <= upper), 1, 0))

predint_w3 <- t( apply(predY_w3, 1, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(PI_length = upper - lower, 
         cover = ifelse((test_Y >= lower) & (test_Y <= upper), 1, 0))

predint_MCMC <- t( apply(predY_MCMC, 1, function(x){
  quantile(x, probs = c(.025,.975))
} ) ) %>% as.data.frame() %>% 
  setNames(., c("lower", "upper")) %>%
  mutate(PI_length = upper - lower, 
         cover = ifelse((test_Y >= lower) & (test_Y <= upper), 1, 0))

# plot PI length
plot.func2(
  MCMC_vec = predint_MCMC$PI_length, 
  w1_vec = predint_w1$PI_length,
  w2_vec = predint_w2$PI_length, 
  w3_vec = predint_w3$PI_length,
  ylab = "Length of Prediction Intervals",
  ylim=c(0,2)
)


# prediction interval coverage
sum(predint_w1$cover)/length(test_Y)
sum(predint_w2$cover)/length(test_Y)
sum(predint_w3$cover)/length(test_Y)
sum(predint_MCMC$cover)/length(test_Y)

# squared prediction error 
SPE_w1 <- apply( cbind(predY_w1, test_Y), 1, function(x){
  (head(x,-1) - tail(x,1))^2
})
SPE_w2 <- apply( cbind(predY_w2, test_Y), 1, function(x){
  (head(x,-1) - tail(x,1))^2
})
SPE_w3 <- apply( cbind(predY_w3, test_Y), 1, function(x){
  (head(x,-1) - tail(x,1))^2
})
SPE_MCMC <- apply( cbind(predY_MCMC, test_Y), 1, function(x){
  (head(x,-1) - tail(x,1))^2
})

# distribution of MSPE across all test cases
MSPE_w1 <- rowMeans(SPE_w1) 
MSPE_w2 <- rowMeans(SPE_w2)
MSPE_w3 <- rowMeans(SPE_w3)
MSPE_MCMC <- rowMeans(SPE_MCMC)

par(mfrow=c(2,2))
hist(MSPE_w1, breaks = 70)
hist(MSPE_w2, breaks = 70)
hist(MSPE_w3, breaks = 70)
hist(MSPE_MCMC, breaks = 70)
par(mfrow=c(1,1))


# distribution of MSPE by each test case
test_MSPE_w1 <- colMeans(SPE_w1)
test_MSPE_w2 <- colMeans(SPE_w2)
test_MSPE_w3 <- colMeans(SPE_w3)
test_MSPE_MCMC <- colMeans(SPE_MCMC)

thresh <- 5
sum(as.numeric(test_MSPE_w1 > thresh))
sum(as.numeric(test_MSPE_w2 > thresh))
sum(as.numeric(test_MSPE_w3 > thresh))
sum(as.numeric(test_MSPE_MCMC > thresh))

which(test_MSPE_w1 > 30) == which(test_MSPE_w2 > 30)
which(test_MSPE_w2 > 30) == which(test_MSPE_w3 > 30)
which(test_MSPE_w3 > 30) == which(test_MSPE_MCMC > 30)

plot.func2(
  MCMC_vec = test_MSPE_MCMC , 
  w1_vec = test_MSPE_w1,
  w2_vec = test_MSPE_w2, 
  w3_vec = test_MSPE_w3,
  ylab = "test_case_MSPE (averaged across posterior samples)",
  ylim = c(0,10)
)
