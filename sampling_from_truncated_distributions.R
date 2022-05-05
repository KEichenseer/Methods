plot(seq(0,1,0.001),qnorm(seq(0,1,0.001),0,1))
u = 1/runif(1,0,pnorm(2,0,1))
u
1/(qnorm(u,0,1))

x <- rep(NA, 100000)
for(i in 1:100000) {
Fa = pnorm(-0.75,0,1)
Fb = pnorm(3,0,1)
v = Fa + (Fb-Fa)*runif(1,0,1)   #/* V ~ U(F(a), F(b))  */
x[i] = qnorm(v,0,1)           #/* truncated normal on [a,b] */
}
hist(x,100)

pexp()

hist(rgamma(100000,4,1)+20,100)
