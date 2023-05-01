install.packages("bootstrap")
install.packages("boot")

library(bootstrap)
library(boot)
data(scor)

lambda_hat <- eigen(cov(scor))$values
set.seed(100)
theta_hat <- lambda_hat[1] / sum(lambda_hat)
B <- 500 # number of bootstrap samples
n <- nrow(scor) # number of rows (data size)
# Bootstrap
func <- function(dat, index){
  # input: dat, data; index, a vector of the bootstrap index
  # output: theta, the estimated theta using the bootstrap sample
  x <- dat[index,]
  lambda <- eigen(cov(x))$values
  theta <- lambda[1] / sum(lambda)
  return(theta)
}
bootstrap_result <- boot(
  data = cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),
  statistic = func, R = B)

theta_b <- bootstrap_result$t
bias_boot <- mean(theta_b) - theta_hat
# the estimated bias of theta_hat, using bootstrap
se_boot <- sqrt(var(theta_b))
# the estimated standard error (se) of theta_hat, using bootstrap
# Jackknife
theta_j <- rep(0, n)
for (i in 1:n) {
  x <- scor [-i,]
  lambda <- eigen(cov(x))$values
  theta_j[i] <- lambda[1] / sum(lambda)
  # the i-th entry of theta_j is the i-th "leave-one-out" estimation of theta
}
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat)
# the estimated bias of theta_hat, using jackknife
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
# the estimated se of theta_hat, using jackknife
# 1
bias_boot
se_boot
# 2
bias_jack
se_jack
# 3
boot.ci(boot.out = bootstrap_result, conf = 0.95, type = c("perc","bca"))


