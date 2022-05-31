sunpdf <- function(x,epsilon,delta,omega,gamma,lambda) {
  1/(pnorm(delta/gamma)) * dnorm(x,epsilon,omega) * pnorm((delta+lambda/omega*(x-epsilon))/sqrt(gamma^2-lambda^2/omega))
}

xseq <- seq(-1,4.5,0.01)

plot(xseq,sunpdf(xseq,1,1,2^2,5,1), type = "l")

rho <- lambda/sqrt(1+lambda^2)
omega <- sigma/sqrt(1-(2*rho^2)/pi)
epsilon <- mu - omega * sqrt(2/pi) * rho
plot(xseq,sunpdf(xseq,epsilon,0,omega,0,0), type = "l")


mu = 1
sigma = 2
lambda = 5

plot(xseq,sunpdf(xseq,epsilon = mu, delta = 1, omega = sigma^2, gamma = sigma^2*(1+lambda^2),lambda = lambda*sigma^2), type = "l")

plot(xseq,sunpdf(xseq,epsilon = mu, gamma = 0, lambda = sigma^2, delta = sigma^2*(1+lambda^2), omega = lambda*sigma^2), type = "l")


points(xseq,fGarch::dsnorm(xseq,mu,sigma,lambda), type = "l", col = "red", lty = 2)


### That one is correct:
plot(xseq,2*dnorm(xseq,mu,sigma)*pnorm(lambda*xseq,lambda*mu,sigma), type = "l", col = "red", lty = 2, yaxs = "i", add = T)

### That's the skew from http://koreascience.or.kr/article/JAKO200504840590864.pdf
points(xseq,2/sigma*dnorm((xseq-mu)/sigma)*pnorm(lambda*(xseq-mu)/sigma), type = "l", col = "blue", lty = 3)

### Good!



dn1 <- sapply(xseq, function(x) sn::dsun(x,
                                         xi = mu,
                                         Omega = matrix(sigma^2),
                                         tau = matrix(lambda*sigma^2),
                                         Delta = 0,
                                         Gamma = matrix(sigma^2 * (1+lambda^2))))

points(xseq,dn1,
                                         type = "l", col = "orange", lty = 3)

ysnd <- sapply(xseq, function(x)
  sn::dsun(x,theta0,matrix(tau0),matrix(lambda*tau0),-lambda*(mu-theta0),matrix(sigma^2+lambda^2*tau0), log = T))


h_lambda <- lambda/sqrt(sigma^2+tau^2*(1+lambda^2))
theta0 <- (x_mu*sigma^2 + mu*tau^2)/(sigma^2+tau^2)
tau0 <- (sigma^2*tau^2)/(sigma^2+tau^2)


points(xseq,sunpdf(xseq,epsilon = xi, lambda = lambda*sigma^2, omega = sigma^2, gamma = (sigma^2 * (1+lambda^2)), delta = 0), type = "l")


x = xseq
xi = mu
Omega = (sigma^2)
tau = 0
Delta = (lambda*sigma^2)
Gamma = (sigma^2 * (1+lambda^2))

tz <- (x - xi)/sqrt(Omega)
p1 <- pnorm(tau + Delta * tz, 0, Gamma - Delta^2)
p2 <- pnorm(tau, 0, Gamma)

pdfN <- dnorm(x, xi, Omega, log = FALSE)
pdfN * p1/p2

points(xseq,pdfN)


if (!(missing(Delta) & missing(Omega)) && !is.null(dp))
  stop("You cannot set both component parameters and 'dp'")
if (!is.null(dp)) {
  if (length(dp) != 5)
    stop("wrong length of non-null 'dp'")
  xi <- mu
  Omega <- matrix(sigma^2)
  Delta <- matrix(lambda * sigma^2)
  tau <- 0
  Gamma <- matrix(sigma^2 * (1+lambda^2))
}
if (!all.numeric(x, xi, Omega, Delta, tau, Gamma))
  stop("non-numeric argument(s)")
d <- dim(Omega)[1]
d
d <- 1

if (length(xi) != d | dim(Omega)[2] != d)
  stop("mismatch of dimensions")
omega <- sqrt(diag(Omega))
Omega.bar <- cov2cor(matrix(Omega))
O.inv <- solve(Omega.bar)
O.inv <- 1
m <- length(tau)
if (m == 1 & !silent)
  warning("When m=1, functions for the SN/ESN distr'n are preferable")
if (any(dim(Gamma) != c(m, m) | dim(Delta) != c(d, m)))
  stop("mismatch of dimensions")
x <- if (is.vector(x))
  matrix(x, 1, d)
else data.matrix(x)
n <- nrow(x)
if (is.vector(xi))
#  xi <- outer(rep(1, n), as.vector(matrix(xi, 1, d)))
#tz <- t(x - xi)/omega
tz <- (1.5 - xi)/omega

D.Oinv <- t(Delta) %*% O.inv
D.Oinv <- Delta
#p1 <- pmnorm(t(tau + D.Oinv %*% tz), rep(0, m), Gamma - D.Oinv %*%
#               Delta, ...)
p1 <- pnorm(tau + Delta * tz, 0, Gamma - Delta^2)
p2 <- pnorm(tau, 0, Gamma)

pdfN <- dnorm(x, xi, Omega, log = log)
if (log)
  pdfN + logb(p1) - logb(p2)
else pdfN * p1/p2
