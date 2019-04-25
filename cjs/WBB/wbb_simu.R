library(mvtnorm)
library(monomvn)
library(glmnet)
library(microbenchmark)
library(foreach)
library(doParallel)

make_design <- function(n, p, snr, rho, case=1) {
  if(case==1)
    beta_0 = rep(1,10)
  if(case==2)
    beta_0 = c(rep(1,5),rep(10,5))
  if(case==3)
    beta_0 = rep(1,p)
  
  sig <- outer(1:p,1:p,function(x,y)rho^abs(x-y))
  x <- matrix(rnorm(n*p),n,p) %*% chol(sig*0.1)
  beta <- c(beta_0,rep(0,p-length(beta_0)))  
  mu <- as.numeric(x %*% beta)
  sigma <- sqrt(sum(mu^2) / (n * snr)) 
  y <- mu + sigma * matrix(rnorm(n), n, 1)
  
  return(list(x=x,y=y,beta=beta))
}

get_mode = function(x)
{
  d = density(x)
  return(d$x[which.max(d$y)])
}

discrepancy = function(beta,beta0)
{
  N = dim(beta)[1]
  p = dim(beta)[2]
  beta0_matrix = matrix(rep(beta0,N),N,p,byrow = TRUE)
  
  beta0_iszero = (beta0_matrix == 0)
  beta_iszero = (beta == 0)
  
  d = mean(beta0_iszero != beta_iszero)
  return(d)
}

wbb = function(x, y, lambda, b=100)
{
  N = dim(x)[1]
  p = dim(x)[2]
  theta = matrix(0,nrow=b, ncol = p)
  
  theta<- foreach(icount(b),.combine = rbind,.packages = "glmnet") %dopar%
  {
    w = rexp(N)
    res = glmnet(x,y,weights = w, penalty.factor = rexp(p), intercept=FALSE)
    as.vector(coef(res,s=lambda))[-1]
  }
  
  return(theta)
}

interval=function(x,lower,upper)
{
  sapply(1:length(x), function(i) x[i]>=lower & x[i] <= upper)
}

simu_table=function(p,n,B=500,case)
{
  WBB_beta_mse_table <- WBB_oos_mse_table <- WBB_coverage_max <- WBB_dis_table <- NULL
  Blasso_beta_mse_table <- Blasso_oos_mse_table <- Blasso_coverage_mean <- Blasso_dis_table <- NULL
  
  Time_ratio_table <- NULL
  
  for(i in 1:length(p))
  {
    cat("n=",n[i],"\t p=",p[i],"\n")
    
    Blasso_beta_mse <- Blasso_oos_mse <- Blasso_dis <- c()
    WBB_beta_mse <- WBB_oos_mse <- WBB_dis <- c()
    WBB_coverage <- NULL
    Blasso_coverage <- NULL
    time_ratio <- c()
    
    for(b in 1:B)
    {
      print(b)
      data <- make_design(n = n[i]*2, p = p[i], snr = 2, rho = 0.8,case=case)
      beta = data$beta
      
      
      bench_time <- microbenchmark({lasso <- cv.glmnet(data$x,data$y,intercept=FALSE,nfolds=3);
      lambda <- lasso$lambda.min; wbb_beta <- wbb(data$x[1:n[i],],data$y[1:n[i]],lambda,b=200)
      },
      Blasso <- blasso(data$x[1:n[i],],data$y[1:n[i]],T=500,icept=FALSE,verb = 0),times=1,unit = "s")
      
      
      Blasso$beta = Blasso$beta[301:500,]
      Blasso_mean = colMeans(Blasso$beta)
      Blasso_prediction = data$x[n[i]+1:n[i],] %*% Blasso_mean
      Blasso_beta_mse = c(Blasso_beta_mse, mean((Blasso_mean - beta)^2))
      Blasso_oos_mse = c(Blasso_oos_mse, mean((data$y[n[i]+1:n[i]]-Blasso_prediction)^2))
      Blasso_coverage = rbind(Blasso_coverage, interval(beta,apply(Blasso$beta,2,function(x)quantile(x,0.025)), apply(Blasso$beta,2,function(x)quantile(x,0.975))))
      Blasso_dis = c(Blasso_dis, discrepancy(Blasso$beta,beta))
      
      WBB_mean = colMeans(wbb_beta)
      wbb_prediction = data$x[n[i]+1:n[i],] %*% WBB_mean
      WBB_beta_mse = c(WBB_beta_mse, mean((WBB_mean - beta)^2))
      WBB_oos_mse = c(WBB_oos_mse, mean((data$y[n[i]+1:n[i]]-wbb_prediction)^2))
      WBB_coverage = rbind(WBB_coverage, interval(beta,apply(wbb_beta,2,function(x)quantile(x,0.025)),apply(wbb_beta,2,function(x)quantile(x,0.975))))
      WBB_dis = c(WBB_dis, discrepancy(wbb_beta, beta))
    }
    
    WBB_dis_table = c(WBB_dis_table,mean(WBB_dis))
    Blasso_dis_table = c(Blasso_dis_table,mean(Blasso_dis))
    Blasso_beta_mse_table = c(Blasso_beta_mse_table, mean(Blasso_beta_mse))
    WBB_beta_mse_table = c(WBB_beta_mse_table,mean(WBB_beta_mse))
    Blasso_oos_mse_table = c(Blasso_oos_mse_table,mean(Blasso_oos_mse))
    WBB_oos_mse_table = c(WBB_oos_mse_table,mean(WBB_oos_mse))
    
    WBB_coverage_mean = c(WBB_coverage_mean,mean(apply(WBB_coverage,2,mean)))
    Blasso_coverage_mean = c(Blasso_coverage_mean,mean(apply(Blasso_coverage,2,mean)))
    
    
  }
  
  t = rbind(n,
            Blasso_beta_mse_table,WBB_beta_mse_table,
            Blasso_oos_mse_table,WBB_oos_mse_table,
            Blasso_coverage_mean,WBB_coverage_mean,
            Blasso_dis_table,WBB_dis_table)
  t = as.data.frame(t)
  colnames(t) = paste0("p=",p)
  return(t)
}


numCores <- detectCores()
registerDoParallel(numCores)
# Nslots <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
# print(sprintf("%d slots were allocated", Nslots))
# cl <- makeCluster(Nslots)
# registerDoParallel(cl)

# case 5
print(5)
p = c(40, 60, 80, 100, 120)
n = p/2
res5 = simu_table(p,n,case=2)
save(file="simu_case5.Rdata",res5)

# case 6
print(6)
p = c(40, 60, 80, 100, 120)
n = p/2
res6 = simu_table(p,n,case=3)
save(file="simu_case6.Rdata",res6)

# case 1
print(1)
p = c(40, 60, 80, 100, 120)
n = rep(50,5)
res1 = simu_table(p,n,case=1)
save(file="simu_case1.Rdata",res1)

# case 2
print(2)
p = c(40, 60, 80, 100, 120)
n = rep(50,5)
res2 = simu_table(p,n,case=2)
save(file="simu_case2.Rdata",res2)

# case 3
print(3)
p = c(40, 60, 80, 100, 120)
n = rep(50,5)
res3 = simu_table(p,n,case=3)
save(file="simu_case3.Rdata",res3)

# case 4
print(4)
p = c(40, 60, 80, 100, 120)
n = p/2
res4 = simu_table(p,n,case=1)
save(file="simu_case4.Rdata",res4)

# case 5
print(5)
p = c(40, 60, 80, 100, 120)
n = p/2
res4 = simu_table(p,n,case=2)
save(file="simu_case5.Rdata",res5)

# case 6
print(6)
p = c(40, 60, 80, 100, 120)
n = p/2
res4 = simu_table(p,n,case=3)
save(file="simu_case6.Rdata",res6)
