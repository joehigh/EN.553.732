#Problem 4

n<-5000 

#Chosen values and Prior Parameters:
theta0 = 1
sig0 = 0.5
v0 = 1
k0 = 1
kn<-k0+n
vn<-v0+n

y<-rnorm(n, mean=0, sd = 1)
yb<-mean(y)
SS<-sum((y-yb)^2)
theta_n<-(k0*theta0+n*yb)/kn
SSn<-(v0*sig0+SS+(k0*n)*(y-theta0)^2/kn)/vn
sig <- 1/rgamma(5000,vn/2,vn*SSn/2)
theta<-rnorm(5000, theta_n, sqrt(sig/kn))
t<-rt(5000, df=vn)*sqrt(SSn/kn)+theta_n
theta_density<-density(theta)
t_dist<-density(t)
plot(theta_density)
lines(t_dist, col="red")