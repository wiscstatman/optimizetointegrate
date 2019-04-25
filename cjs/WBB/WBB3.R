library(glmnet)
library(ggplot2)
library(reshape2)
library(lars)
library(monomvn)
# library(BayesBridge)

########################
##  simulation example
########################

# Data Generate
N = 50
theta = 5
x = rnorm(N,sd=10)
X = as.matrix(data.frame(rep(1,N),x))
y = rnorm(N) + x*theta + 1 # y=1+5x+noise

# Weighted Bootstrap
T = 1000

theta_w = NULL

l = 10
for(t in 1:T)
{
  # weight
  w = rexp(N+1)
  w = w/mean(w)
  
  # lasso
  res = glmnet(X,y,lambda = l*w[N+1]/N,weights = w[1:N])
  theta_w = c(theta_w,res$beta[2])
}

hist(theta_w,breaks = 20)
mean(theta_w) 


# regularization path  (takes about 2 mins)

ptm <- proc.time()
l = exp(seq(0,15,0.1))
T = 500
theta_mean = NULL
for(i in 1:length(l))
{
  
  print(i)
  theta_w2 = NULL
  
  for(t in 1:T)
  {
    # weight
    w = rexp(N+1)
    w = w/mean(w)
    
    # lasso
    res = glmnet(X,y,lambda = l[i]*w[N+1]/N,weights = w[1:N])
    theta_w2 = c(theta_w2,res$beta[2])
  }
  
  theta_mean[i] = mean(theta_w2)
}

proc.time() - ptm
plot(log(l),theta_mean,"l",
     xlab = "Log(lambda)",ylab = "Posterior mean",
     main = "Regularization path")


########################
##  diabetes example
########################
data("diabetes")
# number of observations
N = dim(diabetes$x)[1]
# number of variables
p = dim(diabetes$x)[2] 
# variable names
name = colnames(diabetes$x)

# Center the data.
y = diabetes$y - mean(diabetes$y);
X = scale(diabetes$x)



res0 = cv.glmnet(X,y)
l = res0$lambda.min

T = 2000

theta_w2 = matrix(NA,T,p) # weighted prior
colnames(theta_w2) = name
theta_w3 = matrix(NA,T,p) # weighted individual prior
colnames(theta_w3) = name

# theta_w4 = matrix(NA,T,p) # k-fold cross-validation
# colnames(theta_w4) = name

# theta_w5 = matrix(NA,T,p) # bayes brigde
# colnames(theta_w5) = name
# tt = bridge.reg.stb(y, X, nsamp=T, alpha=1,nu.shape=2.0, nu.rate=2.0)
# theta_w5 = tt$beta

for(t in 1:T)
{
  print(t)
  w = rexp(N)
  wp = rexp(p)
  res = glmnet(X,y,weights = w, penalty.factor = wp)
  res_0 = glmnet(X,y,weights = w)
  # res.cv = glmnet(X,y,weights = w)
  
  theta_w2[t,] = as.vector(coef(res_0,s=l*wp[1]))[-1]
  theta_w3[t,] = as.vector(coef(res,s=l*mean(wp)))[-1]
  # theta_w4[t,] = as.vector(coef(res.cv,s=l))[-1]
}

theta_w4 = blasso(X,y,T=T+500)
theta_w4 = theta_w4$beta[500+1:T,]
colnames(theta_w4) <- name

########################
##  posterior plot
########################

theta_wbb0 = melt(theta_w2)[,-1]
theta_wbb = melt(theta_w3)[,-1]
theta_blasso = melt(theta_w4)[,-1]
# theta_bb = melt(theta_w5)[,-1]


theta = rbind(theta_wbb0,theta_wbb,theta_blasso)
method = c(rep("WBB_one_prior",T*p),rep("WBB",T*p),rep("BLASSO",T*p))
theta = cbind(method,theta)
colnames(theta)[2:3] = c("variable","beta")


p1 = ggplot(theta, aes(beta,colour = method)) +
  geom_density(alpha = 1) + 
  facet_wrap(~ variable,scales = "free",ncol = 2)+
  theme_bw()

p1

