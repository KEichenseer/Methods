x = 4
sigma = 2
alpha = 0
# transform alpha
rho = alpha/sqrt(1+alpha^2)

# Prior parameters
delta = 1
tau = 2

# Gibbs sampler
N = 50000 # number of iterations
mu <- rep(NA,N) # vector to store mu
mu[1] <- 1 # initial value of mu

for(i in 2:N){
  z <- (x - mu[i-1])/sigma
  y <- truncnorm::rtruncnorm(1,0,Inf,rho*z,sqrt(1-rho^2)) # needs truncnorm package
  mu[i] <- rnorm(1,(tau^2*(x-sigma*rho*y)+sigma^2*(1-rho^2)*delta)/(tau^2+sigma^2*(1-rho^2)),
                 sqrt((tau^2*sigma^2*(1-rho^2))/(tau^2+sigma^2*(1-rho^2))))
}

N <- 50000

omega <- sqrt(sigma^2/(1-(2*rho^2)/pi))
epsilon <- x - omega * sqrt(2/pi) * rho
z <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
y <- epsilon + omega*rho*z + omega*sqrt(1-(rho^2))*rnorm(N,0,1)

par(mar = c(2.25,4,.25,.25),las = 1, mgp = c(3,0.8,0))
hist(y,seq(-100,100,0.2), xlim = c(-3.6,8.25), ylim = c(0,3500), col = rgb(0,0,1,0.2), main = NA, xlab = "", yaxt = "n")
hist(rnorm(N,delta,tau),seq(-100,100,0.2), add = T, col = rgb(0,1,0,0.2))
hist(mu,breaks = seq(-100,100,0.1), xlim = c(-3,3), add = T, col = rgb(1,0,0,0.2))
axis(2,at = seq(0,4000,500), labels = c(0,NA,1000,NA,2000,NA,3000,NA,4000))
legend("topleft", c("prior: N(1,2)", "data: SN(4,2,5)", expression("posterior of "*mu)),
       fill = c(rgb(0,1,0,0.25), rgb(0,0,1,0.25), rgb(1,0,0,0.25)), bty = "n")

hist(y2,seq(-100,100,0.05), xlim = c(-3,3), col = rgb(1,0,0,0.2), add = T)

hist(rnorm(50000,mu0,sqrt(var0)),breaks = seq(-100,100,0.1),col = rgb(0,1,0,0.2),add = T)
