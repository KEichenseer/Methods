N <- 50000

mu1 <- 3
sd1 <- 2
a1 <- -10


mu0 <- 5
var0 <- 1^2




skew_mu <- function(x, y, sigma, rho, mu_prior, sigma_prior) {

  rnorm(1,(n*sigma_prior^2*(mean(x)-sigma*rho*mean(y))+sigma^2*(1-rho^2)*mu_prior)/
          (n*sigma_prior^2+sigma^2*(1-rho^2)),
                 sqrt((sigma_prior^2*sigma^2*(1-rho^2))/(n*sigma_prior^2+sigma^2*(1-rho^2))) )
}



epsilon_post <- function(x, y, mu, sigma, epsilon, rho, m, S) {
  n <- length(x)
  sum(log(sigma^(-n+m+1)) + (-((1-rho^2)*S + sum(x-mu)^2 - 2*sigma*rho*sum(x-mu)*y)/
                         2*sigma^2*(1-rho^2)) + (epsilon))
}

exp(log(3^4) + 2.5 + 2)

3^4 * exp(2.5) * exp(2)

rho = (-1+exp(3))/(1+exp(3))
rho
epsilon_post(-c(-1),-3,2,log(2),3,0.9,1,6)

sigma_post <- function(x, mu, sigma, epsilon, lambda, m, S) {
  sum(2*dnorm(x,mu,sigma, log = T)+pnorm(lambda*x,lambda*mu,sigma, log = T) + ginvchi2dlog(sigma,m,S) + (epsilon))
}
(sigma_post(x, mu[i], sigma_new, epsilon_new, lambda[i-1], m, S))
sum(2*dnorm(x,mu[i],sigma_new, log = T))
sum(pnorm(lambda[i-1]*x,lambda[i-1]*mu[i],sigma[i-1], log = T))
ginvchi2dlog(sigma[i-1],m,S)
log(2*exp(0)/(1+exp(0))^2)


xi_post <- function(x, mu, sigma, lambda, xi, mu_xi, sd_xi) {
  sum(2*dnorm(x,mu,sigma, log = T)+pnorm(lambda*x,lambda*mu,sigma, log = T)+dnorm(xi,mu_xi,sd_xi, log = T) + log(2*exp(xi)/(1+exp(xi))^2))
}

lambda_post <- function(x, mu, sigma, lambda, mu_lambda, sd_lambda) {
  sum(2*dnorm(x,mu,sigma, log = T)+pnorm(lambda*x,lambda*mu,sigma, log = T)+dnorm(lambda,mu_lambda,sd_lambda, log = T))
}

x <- sn::rsn(1000,2,2,-10)
mean(sn::rsn(1000,0,2,-10))

hist(x,20)
nIter = 1000

y <- rep(NA,nIter)
mu <- rep(NA,nIter)
sigma <- rep(NA,nIter)
epsilon <- rep(NA,nIter)

lambda <- rep(NA,nIter)
rho <- rep(NA,nIter)
xi <- rep(NA,nIter)

n <- length(x)

mu[1] <- 1
sigma[1] <- 2
epsilon[1] <- log(sigma[1])
lambda[1] <- 1
rho[1] <- 0
xi[1] <- 0

mu_prior <- 0
sigma_prior <- 3

m <- 1
S <- 6

mu_lambda <- 0
sd_lambda <- 4

mu_xi <- 0
sd_xi <- 4

for(i in 2:nIter) {
  z <- (x - mu[i-1])/sigma[i-1]
  y <- truncnorm::rtruncnorm(n,0,Inf,rho[i-1]*z,sqrt(1-rho[i-1]^2))
  mu[i] <- skew_mu(x, y, sigma[i-1], rho[i-1], mu_prior, sigma_prior)
  epsilon_new <- log(sigma[i-1])+rnorm(1,0,0.1)
  sigma_new <- exp(epsilon_new)
  HRepsilon <- exp(sigma_post(x, mu[i], sigma_new, epsilon_new, lambda[i-1], m, S)-
                     sigma_post(x, mu[i], sigma[i-1], log(sigma[i-1]), lambda[i-1], m, S))
  if(HRepsilon > runif(1,0,1)) {sigma[i] <- sigma_new} else {sigma[i] <- sigma[i-1]}

 # xi_new <- xi[i-1]+rnorm(1,0,0.1)
 # rho_new <- (-1+exp(xi_new))/(1+exp(xi_new))
  lambda_new <- lambda[i-1]+rnorm(1,0,0.1)
  HRxi <- exp(lambda_post(x, mu[i], sigma[i], lambda_new, mu_lambda, sd_lambda)-
              lambda_post(x, mu[i], sigma[i], lambda[i-1], mu_lambda, sd_lambda))
  if(HRxi > runif(1,0,1)) {
    xi[i] <- xi_new
    rho[i] <- rho_new
    lambda[i] <- lambda_new
  } else {
    xi[i] <- xi[i-1]
    rho[i] <- rho[i-1]
    lambda[i] <- lambda[i-1]
  }


} # end of nIter

plot(mu)
plot(sigma)
plot(lambda)


xi_new <- 5
rho_new <- (-1+exp(xi_new))/(1+exp(xi_new))
lambda_new <- rho_new/sqrt(1-rho_new^2)
lambda_new

(2*dnorm(xseq,mu1,sigma)*pnorm(lambda*xseq,lambda*mu1,sigma))



ginvchi2d <- function(x,m,S) {
  x^(-m-1) * exp(-S/(2*x^2))
}

ginvchi2dlog <- function(x,m,S) {
  log(x^(-m-1)) + (-S/(2*x^2))
}


######


sigma_post <- function(x, mu, sigma, lambda, a, b) {
  sum(2*dnorm(x,mu,sigma, log = T)+pnorm(lambda*x,lambda*mu,sigma, log = T) + dgamma(sigma,a,b,log = T))
}

mu1 <- 3
sd1 <- 2
a1 <- 4
rho1 <-  a1/sqrt(1+a1^2)
omega1 <- sqrt((sd1^2)/(1-(2*rho1^2)/pi)) #sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * rho1


x <- sn::rsn(1000,3,2,4)
mean(sn::rsn(1000,0,2,-10))

hist(x,20)
nIter = 1000

y <- rep(NA,nIter)
mu <- rep(NA,nIter)
sigma <- rep(NA,nIter)
epsilon <- rep(NA,nIter)

lambda <- rep(NA,nIter)
rho <- rep(NA,nIter)
xi <- rep(NA,nIter)

n <- length(x)

mu[1] <- 1
sigma[1] <- 2
epsilon[1] <- log(sigma[1])
lambda[1] <- 1
rho[1] <- 0
xi[1] <- 0

mu_prior <- 0
sigma_prior <- 3

m <- 1
S <- 6

a = 1
b = 1

plot(xseq,dgamma(xseq,a,b))
mu_lambda <- 0
sd_lambda <- 4

mu_xi <- 0
sd_xi <- 4

for(i in 2:nIter) {
  z <- (x - mu[i-1])/sigma[i-1]
  y <- truncnorm::rtruncnorm(n,0,Inf,-rho[i-1]*z,sqrt(1-rho[i-1]^2))
  mu[i] <- skew_mu(x, y, sigma[i-1], -rho[i-1], mu_prior, sigma_prior)
  sigma_new <- sigma[i-1]+rnorm(1,0,0.1)
  if(sigma_new<=0) {HRsigma <- 0} else {
    HRsigma <- exp(sigma_post(x, mu[i], sigma_new, lambda[i-1], a, b)-
                     sigma_post(x, mu[i], sigma[i-1], lambda[i-1], a, b))
  }
  if(HRsigma > runif(1,0,1)) {sigma[i] <- sigma_new} else {sigma[i] <- sigma[i-1]}

  # xi_new <- xi[i-1]+rnorm(1,0,0.1)
  # rho_new <- (-1+exp(xi_new))/(1+exp(xi_new))
  lambda_new <- lambda[i-1]+rnorm(1,0,0.1)
  HRlambda <- exp(lambda_post(x, mu[i], sigma[i], lambda_new, mu_lambda, sd_lambda)-
                lambda_post(x, mu[i], sigma[i], lambda[i-1], mu_lambda, sd_lambda))
  if(HRlambda > runif(1,0,1)) {
    #xi[i] <- xi_new
    #rho[i] <- rho_new
    lambda[i] <- lambda_new
  } else {
    #xi[i] <- xi[i-1]
    #rho[i] <- rho[i-1]
    lambda[i] <- lambda[i-1]
  }
  rho[i] <- lambda[i]/sqrt(1+lambda[i]^2)


} # end of nIter

plot(mu)
plot(sigma)
plot(lambda)





xseq <- seq(0,20,0.01)

plot(xseq,ginvchi2d(xseq,1,6), type = "l")
points(xseq,7*ginvchi2d(xseq,5,10), type = "l", col = "red")
points(xseq,0.2*ginvchi2d(xseq,1,2), type = "l", col = "blue")
plot(xseq,0.2*ginvchi2d(xseq,0,10), type = "l", col = "green")



for(i in 2:N){
  z[i] <- (x - mu[i-1])/sigma

  y[i] <- truncnorm::rtruncnorm(n,0,Inf,rho*z[i],sqrt(1-rho^2))

  mu[i] <- rnorm(1,(n*var0*(mean(x)-sigma*rho*mean(y[i]))+sigma^2*(1-rho^2)*mu0)/(n*var0+sigma^2*(1-rho^2)),
                 sqrt((var0*sigma^2*(1-rho^2))/(n*var0+sigma^2*(1-rho^2))) )
}
