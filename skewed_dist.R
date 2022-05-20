library(fGarch)
library(truncnorm)


q1 <- seq(-4,4,0.01)
plot(q1,dsnorm(q1,mu1,sd1,a1), type = "l")

N <- 10000

mu1 <- 0
sd1 <- 1
a1 <- 4
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
