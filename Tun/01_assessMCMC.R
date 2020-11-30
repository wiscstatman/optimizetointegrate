library(tidyverse)
library(ggplot2)

load("basad.RData")

# read dataset from csv
train_input <- read_csv("pria_rmi_cv/file_0.csv")
for(i in 1:4){
  filename <- paste0("pria_rmi_cv/file_",i,".csv")
  train_input <- rbind( train_input, read_csv(file = filename) )
}
test_input <- read_csv("pria_prospective.csv")

# parameters
p <- 1024 # dim size
n_train <- nrow(train_input)
n_test <- nrow(test_input)
burn_in <- 5E4

# function to process dataset
list_func <- function(input_df){
  # convert character string to columns of predictors and standardize column-wise
  X <- unname( t( sapply( input_df$Fingerprints, function(x){ as.numeric(str_split_fixed(x,"",p)) } ) ) )
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

#-----------------------------------------------------------------------------------------
# variable selection -- binary entries
Z <-  fit_basad$allZ[(burn_in + 1):nrow(fit_basad$allZ),-1]

# check trace plots of number of variables selected in successive iterations
plot(1:nrow(Z), rowSums(Z), type = "b", pch = ".", xlab = "iteration", 
     ylab = "Number of Variables Selected",
     main = "Trace plots for number of variables selected")

#-----------------------------------------------------------------------------------------
# check trace plots for coeffs of a few variables

# sort columns (variables) according to average Pr(variables selected)
beta <- fit_basad$allB[(burn_in + 1):nrow(fit_basad$allB),-1]
colnames(beta) <- names(fit_basad$posteriorZ[-1])
iii_beta <- match( names(sort(fit_basad$posteriorZ[-1])), colnames(beta) ) 
colsort_beta <- beta[,iii_beta]

head(fit_basad$posteriorZ[-1][iii_beta]) # check
tail(fit_basad$posteriorZ[-1][iii_beta]) # check
fit_basad$posteriorZ[-1][iii_beta][c(811,812,813)] # check

# trace plots of a few variables with low average Pr(variables selected)
plot(1:nrow(colsort_beta), colsort_beta[,1], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[1]))
plot(1:nrow(colsort_beta), colsort_beta[,2], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[2]))
plot(1:nrow(colsort_beta), colsort_beta[,3], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[3]))

# trace plots of a few variables with high average Pr(variables selected)
plot(1:nrow(colsort_beta), colsort_beta[,1024], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[1024]))
plot(1:nrow(colsort_beta), colsort_beta[,1023], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[1023]))
plot(1:nrow(colsort_beta), colsort_beta[,1022], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[1022]))

# trace plots of a few variables with medium average Pr(variables selected)
plot(1:nrow(colsort_beta), colsort_beta[,811], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[811]))
plot(1:nrow(colsort_beta), colsort_beta[,812], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[812]))
plot(1:nrow(colsort_beta), colsort_beta[,813], type = "b", pch = ".", xlab = "iteration", 
     ylab = "beta",
     main = paste0("Trace plots for coeff of variable ", colnames(colsort_beta)[813]))

rm(colsort_beta); gc()
rm(fit_basad);gc()

#-----------------------------------------------------------------------------------------
# check trace plots for coeff of determination R^2

SSTO <- sum( (train_list$Y - mean(train_list$Y))^2 )

# basad_fitted <- apply(beta , 1, function(x){
#   round( train_list$X %*% x, 4 )
# })
# Error: vector memory exhausted (limit reached?)

# basad_fitted <- matrix(NA, nrow = nrow(train_list$X), ncol = nrow(beta)) 
# Error: vector memory exhausted (limit reached?)
SSE <- rep(NA, nrow(beta))

for(i in 1:nrow(beta)){
  # basad_fitted[,i] <- round( train_list$X %*% beta[i,], 4 )
  # SSE[i] <- sum( ( train_list$Y - basad_fitted[,i] )^2 )
  basad_fitted <- train_list$X %*% beta[i,]
  SSE[i] <- sum( ( train_list$Y - basad_fitted )^2 )
  if(i %% 100 == 0){
    print(i)
  }
}

save(SSE, file = "SSE_train.RData")
# save(basad_fitted, file = "basad_fitted.RData")
# rm(basad_fitted); gc()

R2 <- 1 - SSE/SSTO

plot(1:nrow(beta), R2, type = "b", pch = ".", xlab = "iteration", 
     ylab = "Coefficient of Determination",
     main = "Trace plots for coefficient of determination R2")


#-----------------------------------------------------------------------------------------
# check if X_train and X_test are full ranked

gram_Xtrain <- crossprod(train_list$X)
gram_Xtest <- crossprod(test_list$X)
eigenval_gram_Xtrain <- eigen(gram_Xtrain)$values
eigenval_gram_Xtest <- eigen(gram_Xtest)$values
max(eigenval_gram_Xtrain) ; min(eigenval_gram_Xtrain)
max(eigenval_gram_Xtest) ; min(eigenval_gram_Xtest)
rm(gram_Xtrain,gram_Xtest); gc()
