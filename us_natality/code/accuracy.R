library("expm")
library("MCMCpack")

Accuracy <- function(res_var, full_var, nobs){
  
  
  res_var_mean <- apply(res_var, 2, mean) 
  full_var_mean <- apply(full_var, 2, mean)  
  
  res_var_cov <- cov(res_var) 
  full_var_cov <- cov(full_var) 
  
  accura_square <- (dist(rbind(res_var_mean, full_var_mean), method = "euclidean"))^2  +
    sum(diag(res_var_cov + full_var_cov - 2*(sqrtm(sqrtm(res_var_cov)%*%full_var_cov%*%sqrtm(res_var_cov)))))
  
  accuracy <- sqrt(accura_square)*sqrt(nobs)
  return (accuracy)
}  
  
  
  