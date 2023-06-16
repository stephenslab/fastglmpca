# A short script to verify my calculations for the quadratic
# approximation to exp(x) that is intended to be most accurate for
# small x.
library(pracma)
ans <- polyApprox(exp,a = -4,b = 0,n = 2)
b0 <- ans$p[3]
b1 <- ans$p[2]
b2 <- ans$p[1]
f <- ans$f
x <- seq(-4,2,length.out = 1000)
plot(x,exp(x),type = "l",lwd = 1.5)
lines(x,f(x),lwd = 1.5,col = "dodgerblue")
lines(x,b0 + b1*x + b2*x^2,lwd = 1.5,col = "red",lty = "dashed")
