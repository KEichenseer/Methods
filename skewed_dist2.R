### Skew 2

# Try this: https://www.sciencedirect.com/science/article/pii/S0888613X14000504

library()

N <- 50000

mu1 <- 8
sd1 <- 1
a1 <- 5

sigma1 <- a1/sqrt(1+a1^2)
omega1 <- sqrt((sd1^2)/(1-(2*sigma1^2)/pi)) #sd1/sqrt(1-(2*sigma1^2)/pi)
epsilon1 <- mu1 - omega1 * sqrt(2/pi) * sigma1

#sigma1 <- a1/sqrt(1+a1^2)
#omega1 <- 1
#epsilon1 <- 0


z1 <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
#y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
y1 <- epsilon1 + omega1*sigma1*z1 + omega1*sqrt(1-(sigma1^2))*rnorm(N,0,1)
# from https://academic.oup.com/biostatistics/article/11/2/317/268224   # eq. 2.3
#hist(y1,100)
#hist(rsnorm(5000,mu1,sd1,a1),50, add = T, col = rgb(0,1,0,0.2))
y2 <- fGarch::rsnorm(N,2,4,5)
y3 <- sn::rsn(N,epsilon1,omega1,a1)
mean(y1)
sd(y1)

mu <- epsilon1 # 4
sigma <- omega1 # 2
lambda <- a1

x_mu <- 1
tau <- 1


h_lambda <- lambda/sqrt(sigma^2+tau^2*(1+lambda^2))
theta0 <- (x_mu*sigma^2 + mu*tau^2)/(sigma^2+tau^2)
tau0 <- (sigma^2*tau^2)/(sigma^2+tau^2)

post_mu <- theta0 + h_lambda * (sigma*tau^2)/sqrt(sigma^2 + tau^2) * (
  dnorm(h_lambda*(sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)),0,1)/
  pnorm(h_lambda*(sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)),0,1))

post_var <- tau0 - h_lambda^2 * ((sigma^2)*tau^4)/(sigma^2+tau^2) * (
  dnorm(h_lambda*(sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)),0,1)/
    pnorm(h_lambda*(sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)),0,1)) * (
      h_lambda * (sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)) + (
          dnorm(h_lambda*(sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)),0,1)/
            pnorm(h_lambda*(sigma*(x_mu - mu)/sqrt(sigma^2 + tau^2)),0,1))
    )

  post_mu
  post_var

library(sn)

  ysn <- rsun(n=50000,theta0,matrix(tau0),matrix(lambda*tau0),-lambda*(mu-theta0),matrix(sigma^2+lambda^2*tau0))
  ysn2 <- rsun(n=50000,theta0,matrix(tau0),matrix(lambda*tau0),-lambda*(mu-theta0),matrix(sigma^2+lambda^2*tau0))

  hist(y3,seq(-100,100,0.1), xlim = c(-2,11), ylim = c(0,4000))
mean(y3)
  #hist(y2,seq(-100,100,0.1), xlim = c(-3,3), col = rgb(1,0,0,0.2), add = T)

  hist(rnorm(50000,x_mu,tau),breaks = seq(-100,100,0.1),col = rgb(0,1,0,0.2),add = T)

  #hist(rnorm(50000,post_mu,sqrt(post_var)),breaks = seq(-100,100,0.1),col = rgb(1,0,0,0.2),add = T)

  hist(ysn,breaks = seq(-100,100,0.1),col = rgb(1,0,1,0.2),add = T)

  ###
  ###
  ###

  mu2 <- mu+sd1 * sqrt(2/pi) * sigma1
  sd2 <- sqrt(omega1^2 - omega1^2*2*sigma1^2/pi)
  hist(mu2[1:50000],breaks = seq(-100,100,0.1), xlim = c(-3,3), add = T, col = rgb(0,0,1,0.33))



  #rnorm(1,(n*var0*(mean(x)-sigma*rho*mean(y[i]))+sigma^2*(1-rho^2)*mu0)/(n*var0+sigma^2*(1-rho^2)),
   #            sqrt((var0*sigma^2*(1-rho^2))/(n*var0+sigma^2*(1-rho^2))) )
