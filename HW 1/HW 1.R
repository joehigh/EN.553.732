#Problem 4
#Part 2

theta<-c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
Y<-rep(1,11)
for(i in 1:11)
  Y[i]<-choose(100,57)*(theta[i]^57)*(1-theta[i])^43
print(Y)

plot(theta, Y, type = "h", main = "Problem 4, Part 2", xlab = expression(paste(theta)), 
     ylab="Pr(Y=57 | theta)")

#Part 3

x<-rep(1,11)
x1<-rep(1,11)
for(i in 1:11)
  x[i]<-(choose(100,57)*(theta[i]^57)*(1-theta[i])^43)*(1/11)
NormConstant<-1/(sum(x))
x1<-x*NormConstant
print(x1)

plot(theta, x1, type = "h", main = "Problem 4 Part 3", xlab = expression(paste(theta)),
     ylab = "Posterior")

#Part 4
f<-curve(choose(100,57)*x^57*((1-x)^43), from = 0, to = 1, main = "Problem 4, Part 4", 
         xlab = expression(paste(theta)), ylab = "posterior")

#Part 5
v<-seq(0, 1, length = 200)
z<-dbeta(v, 58, 44)
plot(v, z, type = "l", main = "Problem 4, Part 5 Beta(58,44) Distribution",
     xlab = expression(paste(theta)), ylab = "Posterior")


#Problem 5
theta_0<-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) 
n_0<-c(1,2,8,16,32)

a<-matrix(0L, nrow =length(theta_0), ncol =length(n_0)) 
b<-matrix(0L, nrow =length(theta_0), ncol =length(n_0)) 
for (i in 1:length(theta_0)) 
{for (j in 1:length(n_0)) 
{a[i,j]=theta_0[i]*n_0[j] 
b[i,j]=(1-theta_0[i])*n_0[j] 
} 
}
Pr<-matrix(0L, nrow =length(theta_0), ncol =length(n_0)) 
for (i in 1:length(theta_0))
{for (j in 1:length(n_0))
{ 
  f <- function(x) 
  {choose(100,57)*(x^57)*((1- x)^43)*(gamma(a[i,j]+b[i,j])/
                      (gamma(a[i,j])*gamma(b[i,j])))*(x^(a[ i,j]-1))*(1-x)^(b[i,j]-1)} 
  bot<-integrate(f,0, 1, rel.tol=1e-10)$value 
  top<-integrate(f,0.5, 1, rel.tol=1e-10)$value 
  Pr[i,j]<-top/bot
  
} 
}
contour(theta_0, n_0, Pr,main = "Problem 5 Countour Plot", xlab=expression(paste(theta)), 
        ylab='n_0')