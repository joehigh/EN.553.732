#Problem4
#part a
glucose = read.table("glucose.dat", header = FALSE);
data=as.matrix(glucose)
data=as.numeric(data)
hist(data,breaks=seq(50,200,5),freq=FALSE,main="Problem 4 Part a")
lines(density(data))
  
#part c
y=data
set.seed(123)
n=length(y)
iter=10000
a=1
b=1
mu0=120
tau0.sq=200
sigma0.sq=1000
nu0=10
x=matrix(0,iter, n)
p=numeric(iter)
theta1=numeric(iter)
theta2=numeric(iter)
sigma1.sq=numeric(iter)
sigma2.sq=numeric(iter)
p[1]=rbeta(1, a, b)
x[1,]=rbinom(n,1,p[1])
theta1[1]=rnorm(1,mu0,sqrt(tau0.sq))
theta2[1]=rnorm(1,mu0,sqrt(tau0.sq))
sigma1.sq[1]=1/rgamma(1, nu0/2, nu0*sigma0.sq/2)
sigma2.sq[1]=1/rgamma(1, nu0/2, nu0*sigma0.sq/2)
for (i in 2:iter)
{
  for (j in 1:n) 
  {
    y1=dnorm(y[j], theta1[i-1], sqrt(sigma1.sq[i-1]))
    y2=dnorm(y[j], theta2[i-1], sqrt(sigma2.sq[i-1]))
    x[i,j]=rbinom(1,1,(p[i-1]*y1)/(p[i-1]*y1+(1-p[i-1])*y2))
  }
    
c=sum(x[i,])
p[i]=rbeta(1,a+c,b+n-c)
y_1.bar=mean(y[x[i,]==1])
mu_n=(mu0/tau0.sq+c*y_1.bar/sigma1.sq[i-1])/(1/tau0.sq+c/sigma1.sq[i-1])
tau2_n=1/(1/tau0.sq+c/sigma1.sq[i-1])
theta1[i]=rnorm(1, mu_n, sqrt(tau2_n))
    
nu_n=nu0+c
s2_n=sum((y[x[i,] == 1]-theta1[i])^2)/c
sigma2_n=(nu0*sigma0.sq+c*s2_n)/nu_n
sigma1.sq[i]=1/rgamma(1,nu_n/2,nu_n*sigma2_n/2)
    
y_2.bar=mean(y[x[i,] == 0])
mu_n=(mu0/tau0.sq+(n-c)*y_2.bar/sigma2.sq[i-1])/(1/tau0.sq+(n-c)/sigma2.sq[i-1])
tau2_n=1/(1/tau0.sq+(n-c)/sigma2.sq[i-1])
theta2[i]=rnorm(1, mu_n, sqrt(tau2_n))
    
nu_n=nu0+(n-c)
s2_n=sum((y[x[i,] == 0]-theta2[i])^2)/(n-c)
sigma2_n=(nu0*sigma0.sq+(n-c)*s2_n)/nu_n
sigma2.sq[i]=1/rgamma(1,nu_n/2,nu_n*sigma2_n/2)
}
theta_1s=rep(0,iter)
theta_2s=rep(0,iter)
for (i in 1:iter) 
{
  theta_1s[i]=min(theta1[i], theta2[i])
  theta_2s[i]=max(theta1[i], theta2[i])
}
acf(theta_1s,main="Problem 4 Part c, theta1 ")
acf(theta_2s,main="Problem 4 Part c, theta2 ")
effectiveSize(theta_1s)
effectiveSize(theta_2s)
  
#part d
x_1=rbinom(length(p), 1, p)
y_1=numeric(iter)
for (i in 1:iter)
{
  if (x_1[i] == 1)
  {
    y_1[i]=rnorm(1, theta1[i], sqrt(sigma1.sq[i]))
  } 
  else
  {
   y_1[i]=rnorm(1, theta2[i], sqrt(sigma2.sq[i]))
  }
}
hist(y_1,breaks=seq(0,300,5),freq=FALSE,main="Problem 4 Part d")
lines(density(y))