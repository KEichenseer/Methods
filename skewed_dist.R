library(fGarch)
library(truncnorm)


q1 <- seq(-4,4,0.01)
plot(q1,dsnorm(q1,mu1,sd1,a1), type = "l")

N <- 20000

mu1 <- 0
sd1 <- 1
a1 <- 5
sigma1 <- a1/sqrt(1+a1^2)
omega1 <- sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * sigma1


z1 <- rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
#y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
y1 <- epsilon1 + omega1*sigma1*z1 + omega1*sqrt(1-(sigma1^2))*rnorm(N,0,omega1)

hist(y1,100, xlim = c(-4,4), col = rgb(0,0,1,0.33))
hist(y1,100, xlim = c(-4,4),col = rgb(1,0,0,0.33), add = T)

hist((y1-2)/2,100, xlim = c(-4,4),col = rgb(0,1,0,0.33), add = T)

mu1+sd1*sqrt(2/pi)*sigma1

location <- mu1
mean(y1)


#### Try Stan

library("rstan") # observe startup messages
## utter failure, crashes / doesn't install

### end Stan

mu0 <- 0
sd0 <- 1

yobs <- 1
sd1 <- 1

N <- 100000
yestimate          = rnorm(N,
                             sd0^2/(sd1^2+sd0^2)*yobs + sd1^2/(sd1^2+sd0^2)*mu0,
                             sqrt(1/sd0^2 + 1/sd1^2))


hist(rnorm(N,mu0,sd0), breaks = seq(-100,100,0.1), col = rgb(0,0,0,0.25), xlim = c(-6,6), ylim = c(0,4000))
hist(rnorm(N,yobs,sd1), breaks = seq(-100,100,0.1), col = rgb(1,0,0,0.25), add = T)

hist(yestimate, breaks = seq(-100,100,0.1), col = rgb(0,0,1,0.25), add = T)


http://koreascience.or.kr/article/JAKO200504840590864.pdf

alpha = seq(-10,10,0.1)
rho = alpha/sqrt(alpha^2+1)
plot(alpha,rho)

omega1 <- sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * sigma1


N <- 50000

mu1 <- 8
sd1 <- 1
a1 <- 5


mu0 <- 1
var0 <- 1.5^2

sigma1 <- a1/sqrt(1+a1^2)
omega1 <- sd1 #sqrt((sd1^2)/(1-(2*sigma1^2)/pi)) #sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 # mu1 - omega1 * sqrt(2/pi) * sigma1

#sigma1 <- a1/sqrt(1+a1^2)
#omega1 <- 1
#epsilon1 <- 0


z1 <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
#y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
y1 <- epsilon1 + omega1*sigma1*z1 + omega1*sqrt(1-(sigma1^2))*rnorm(N,0,1)
# from https://academic.oup.com/biostatistics/article/11/2/317/268224   # eq. 2.3
#hist(y1,100)
#hist(rsnorm(5000,mu1,sd1,a1),50, add = T, col = rgb(0,1,0,0.2))

mean(y1)
sd(y1)


n =1
sigma_ori <- sd(y1)
#x <- y1[1:n]
x_ori <- mean(y1)
alpha = a1
rho = alpha/sqrt(alpha^2+1)

#x <- sigma_ori/sqrt(1-(2*rho^2)/pi)
#sigma <- x_ori - omega1 * sqrt(2/pi) * rho

sigma <- sqrt((sigma_ori^2)/(1-(2*rho^2)/pi))# sigma_ori/sqrt(1-(2*rho^2)/pi)
x <- x_ori - sigma * sqrt(2/pi) * rho

x <- 8
sigma <- 1

N <- 50000

mu1 <- 8
sd1 <- 1
a1 <- 5

sigma1 <- a1/sqrt(1+a1^2)
omega1 <- sqrt((sd1^2)/(1-(2*sigma1^2)/pi)) #sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * sigma1

x<-epsilon1
sigma <-omega1
rho <- a1/sqrt(1+a1^2)

N <- 50000
mu <- rep(NA,N)
mu[1] <- 1

z <- rep(NA,N)
y <- rep(NA,N)


for(i in 2:N){
z[i] <- (x - mu[i-1])/sigma

y[i] <- truncnorm::rtruncnorm(n,0,Inf,rho*z[i],sqrt(1-rho^2))

mu[i] <- rnorm(1,(n*var0*(mean(x)-sigma*rho*mean(y[i]))+sigma^2*(1-rho^2)*mu0)/(n*var0+sigma^2*(1-rho^2)),
               sqrt((var0*sigma^2*(1-rho^2))/(n*var0+sigma^2*(1-rho^2))) )
}

#hist(mu[40000:50000],breaks = seq(-100,100,0.1), xlim = c(-3,5), add = T, col = rgb(0,0,1,0.33))
#hist(rsnorm(N,x,sigma,-3), breaks = seq(-100,100,0.1), col = rgb(0,1,0,0.2), xlim = c(0,6), add = T)



mu1 <- mean(mu[1000:N])
mu1
sd1
a1
sigma1 <-  a1/sqrt(1+a1^2)
omega1 <- sqrt((sd1^2)/(1-(2*sigma1^2)/pi)) #sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * sigma1
z1 <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
#y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
y2 <- epsilon1 + omega1*sigma1*z1 + omega1*sqrt(1-(sigma1^2))*rnorm(N,0,1)

hist(y1,seq(-100,100,0.1), xlim = c(-2,9), ylim = c(0,2000), add = T, col = rgb(.9,.9,0,0.2))

#hist(y2,seq(-100,100,0.1), xlim = c(-3,3), col = rgb(1,0,0,0.2), add = T)

hist(rnorm(50000,mu0,sqrt(var0)),breaks = seq(-100,100,0.1),col = rgb(0,1,0,0.2),add = T)
hist(mu[1:50000],breaks = seq(-100,100,0.1), xlim = c(-3,3), add = T, col = rgb(0,0,1,0.33))



mu2 <- mu+sd1 * sqrt(2/pi) * sigma1
sd2 <- sqrt(omega1^2 - omega1^2*2*sigma1^2/pi)
hist(mu[1:50000],breaks = seq(-100,100,0.1), xlim = c(-3,3), add = T, col = rgb(1,0,0,0.33))

mean(mu2)
mean(y1)



### c.f. normal-normal
pred <- -7
sdy <- 2
yd_mu <- -3
yd_sd <- 1

y3 <- rnorm(N,
      sdy^2/(yd_sd^2+sdy^2)*yd_mu + yd_sd^2/(yd_sd^2+sdy^2)*pred,
      sqrt((sdy^2 * yd_sd^2)/(sdy^2 + yd_sd^2))) # fix this in the main sampler as well
hist(y3[40000:50000],breaks = seq(-100,100,0.05), xlim = c(-3,3), add = T, col = rgb(1,0,1,0.33))




#### JAGS implementation

model <- function() {

  ## Likelihood
  for (i in 1:n){
    y_measurement[i] ~ dnorm(y_estimate[i], 1/(y_se[i]*y_se[i]))
    y_estimate[i] ~ dnorm(mu0, sd0)
  }


  ## Priors
mu0 <- -3
sd0 <- 1

}


y_measurement <- -7
y_se <- 2
n <- 1


regression_data <- list("y_measurement","y_se","n")

fit  <- R2jags::jags(data = regression_data,
                parameters.to.save = c(
                                       "y_estimate"
                ),
                n.iter = 50000,
                n.thin = 1,
                n.chains =  2, # Other values set at default (for simplicity)
                model.file = model)

y4 <- fit$BUGSoutput$sims.list$y_estimate

hist(y4[40000:50000],breaks = seq(-100,100,0.05), xlim = c(-3,3), add = T, col = rgb(1,1,0,0.33))



model <- function() {

  ## Likelihood

    z1 ~ dnorm(0,1);T(0,)
    e1 ~ dnorm(0,1)
    y1 <- epsilon1 + omega1*rho*z1 + omega1*sqrt(1-(rho^2))*e1

    y1 ~ dnorm(y_estimate, 1/(y_se*y_se))
    y_estimate ~ dnorm(mu0, sd0)


  ## Priors
  mu0 <- -3
  sd0 <- 1

}


y_measurement <- -7
y_se <- 2
n <- 1


mu1 <- y_measurement
sd1 <- y_se
a1 <- 4
rho <- a1/sqrt(1+a1^2)
omega1 <- sd1/sqrt(1-(2*rho^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * rho



regression_data <- list("epsilon1","omega1","rho", "y_se")

fit  <- R2jags::jags(data = regression_data,
                     parameters.to.save = c(
                       "y_estimate"
                     ),
                     n.iter = 50000,
                     n.thin = 1,
                     n.chains =  2, # Other values set at default (for simplicity)
                     model.file = model)

y4 <- fit$BUGSoutput$sims.list$y_estimate

hist(y4[40000:50000],breaks = seq(-100,100,0.05), xlim = c(-3,3), add = T, col = rgb(1,1,0,0.33))



###
# Ok, this seems to work after all. Turns out that it is difficult for the prior to
# draw the estimated mu into the direction of the long tail of the skewed distribution.
# But should be ok. Let's put this into our lat grad model and triumph!

sigma2 <- sd(mu)/sqrt(1-(2*rho^2)/pi)
x2 <- mu - sigma2 * sqrt(2/pi) * rho
hist(x2, breaks = seq(-100,100,0.1), col = rgb(0,0,1,0.1), xlim = c(0,6), add = T)

0.5*(4-pi)*(mean(mu)^3)/(var(mu)^(3/2))

##################################

alpha = seq(-10,10,0.1)
rho = alpha/sqrt(alpha^2+1)
plot(alpha,rho)

mu1 <- 3
sigma1 <- 1
alpha = -2
rho = alpha/sqrt(alpha^2+1)

omega1 <- sigma1/sqrt(1-(2*rho^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * rho


mu0 <- 3
var0 <- 100



N <- 100000
mu <- rep(NA,N)
mu[1] <- 3

for(i in 2:N){
  z <- (x - mu[i-1])/omega1

  y <- truncnorm::rtruncnorm(1,0,Inf,rho*z,1-rho^2)

  mu[i] <- rnorm(1,(1*var0*(epsilon1-omega1*rho*y)+omega1^2*(1-rho^2)*mu0)/(1*var0+omega1^2*(1-rho^2)),
                 (var0*omega1^2*(1-rho^2))/(1*var0+omega1^2*(1-rho^2)) )
}

hist(mu,1000)
