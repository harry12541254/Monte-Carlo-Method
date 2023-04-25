#1
n <- 10000
R <- 100 #模擬次數
g <- function(x){
  (x^2)*(x>1)
  } 
true <- 0.400626
integrals  <- numeric(100)
for(i in 1:R){
  x <- rnorm(n)
  integrals[i] <- mean(g(x))
}
g
theta1 <- mean(integrals)
var1 <- var(integrals)/R
print(c(theta1,var1))
plot(1:R, integrals,type= 'l',col = 'blue',xlab='Iteration', ylim=c(0.2,1),
     ylab ='theta_hat', lwd =2)
abline(a=true,b=0,col='red',lwd=2)

#2
# Importance sampling method 
g <- function(x) {(x^2/sqrt(2*pi))*exp(-x^2/2)}
f <- function(x) {exp(-(x-1))} 
# 重要性分布 shift exp(1) to 1
integrals_2  <- numeric(100)
set.seed(100)
n <- 10000
true <- 0.400626
for ( i in 1:R){
  x <- rexp(n,1) +1
  weit <- g(x) / f(x)
  integrals_2[i] <- mean(weit)
}

theta2_1 <- mean(integrals_2) 
var2_1 <- var(integrals_2)/R
print(c(theta2_1, var2_1)) 

plot(1:R, integrals_2, type="l", xlab="Number of iterations", ylab="Estimated integral value",ylim=c(0.35,0.43))

abline(a=true,b=0,col='red',lwd=2)

# self-normalized var減少但是偏誤

integrals_2_2  <- numeric(100)
for (i in 1:R){
  x <- rexp(n,1) +1
  w <- g(x)/f(x)
  p <- w /sum(w)
  x.star <- sample(x, size=n, prob=p, replace=TRUE)
  weit <- g(x.star) / f(x.star)
  integrals_2_2[i] <- mean(weit)
}
theta2_2 <- mean(integrals_2_2)
var2_2 <- var(integrals_2_2)/R
print(c(theta2_2, var2_2)) 

plot(1:R, integrals_2_2, type="l", xlab="Number of iterations", ylab="Estimated integral value",ylim=c(0.35,0.5))

abline(a=true,b=0,col='red',lwd=2)

#3
# transformation of the integrand to a interval of [0, 1]  u=1-1/x
g <- function(x) {(x^2/sqrt(2*pi))*exp(-x^2/2)}
trans_g <- function(u) {
  x <- 1/(1 - u)
  jacobian <- 1/((1 - u)^2)
  return (g(x) * jacobian)
}
integrals_3 <- numeric(100)
for(i in 1:R){
  x <- runif(n) 
  integrals_3[i] <- mean(trans_g(x))
}

theta3 <- mean(integrals_3)
var3 <- var(integrals_3)/R
print(c(theta3,var3))
plot(1:R, integrals_3, type="l", xlab="Number of iterations", ylab="Estimated integral value",ylim=c(0.35,0.43))
abline(a=true,b=0,col='red',lwd=2)

#4
n <- 10000
integrals_4 <- sigma <- numeric(100)
g <- function(x){
  (x^2)*(x>1)
}

f1 <- function(x){
  dgamma(x,shape=2,scale=2)
}
f2<- function(x){
  0.5*(x^2 + (1-x)^2)
}
theta.hat <- se <- numeric(2)

for(i in 1:R){
  x <- rnorm(n)
  L <- lm(g(x) ~ f1(x) + f2(x))
  c.star <- -c(L$coeff[2],L$coeff[3])
  mu1 <- mean(f1(x))
  mu2 <- mean(f2(x))
  integrals_4[i] <- sum(L$coefficients*c(1,mu1, mu2))
}

theta4 <- mean(integrals_4)
var4 <- var(integrals_4)/R
print(c(theta4,var4))
plot(1:R, integrals_4, type="l", xlab="Number of iterations", ylab="Estimated integral value",ylim=c(0.35,0.43))
abline(a=true,b=0,col='red',lwd=2)

#5
set.seed(1999)
g <- function(x) {(x^2/sqrt(2*pi))*exp(-x^2/2)}
stra <- function(k){
  n <- 10000
  N <- 50
  T2 <- numeric(k)
  estimates <- matrix(0,N,2)
  for (i in 1:N){
    estimates[i,1] <- 0.5-mean(g(runif(n)))
    for(j in 1:k){
      T2[j] <- 0.5-mean(g(runif(n/k, (j-1)/k, j/k)))
    }
    estimates[i,2] <- mean(T2)
  }
  return (estimates)
}

estimates_10 <- stra(k=10)
apply(estimates_10,2,mean) #0.4006631 0.4006315
apply(estimates_10,2,var) #5.398601e-07 6.902123e-09

estimates_1000 <- stra(k=1000)
apply(estimates_1000,2,mean) #0.4007636 0.4006259
apply(estimates_1000,2,var) #9.427131e-07 6.771316e-13


#6
n <- 10000
g <- function(x) {
  (x^2 * exp(-0.5*x^2) / sqrt(2*pi)) * (x>1)
}
ft <- function(t,x){
  (exp(t*x)*(1/sqrt(2*pi))*exp(-x^2/2))/exp((t^2)/2)
}
integrals6  <- numeric(100)
t <- 3

for(i in 1:R){
  u <- runif(n)
  x <- t + qnorm(u)
  fg <- g(x)/ft(t,x)
  integrals6[i] <- mean(fg)
}

theta6 <- mean(integrals6)
var6 <- var(integrals6)/R
print(c(theta6,var6))
plot(1:R, integrals6, type="l", xlab="Number of iterations", ylab="Estimated integral value",ylim=c(0.35,0.43))
abline(a=0.40062,b=0,col='red')

n <- 10000
R <- 100 #模擬次數
g <- function(x) {
  (x^2 * exp(-0.5*x^2) / sqrt(2*pi)) * (x>1)
}
ft <- function(t,x){
  (exp(t*x)*(1/sqrt(2*pi))*exp(-x^2/2))/exp((t^2)/2)
}

t_values <- seq(-10, 10, by=1) # t值從-10到10
results <- matrix(0, length(t_values), 2)

for (t_index in 1:length(t_values)) {
  t <- t_values[t_index]
  integrals <- numeric(R)
  
  for(i in 1:R){
    u <- runif(n)
    x <- t + qnorm(u)
    fg <- g(x)/ft(t,x)
    integrals[i] <- mean(fg)
  }
  
  results[t_index, 1] <- mean(integrals)
  results[t_index, 2] <- var(integrals)/n
}

plot(t_values, results[,1], type='l', col='blue', xlab='t values', ylab='Estimated value', lwd=2)
