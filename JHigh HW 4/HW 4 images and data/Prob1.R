#Problem 1
#Part a
set.seed(1)
t=10000
azdiabetes = read.table("azdiabetes.dat", header = TRUE);
y=as.matrix(azdiabetes[,2])
g=length(y)
X=as.matrix(cbind(rep(1,g),azdiabetes[,-c(2,8)]))
n=nrow(X)
m=ncol(X)
nu_0=2
sigma2_0=1
Hg=(g/(g+1))*X%*%solve(t(X)%*%X)%*%t(X)
SS=t(y)%*%(diag(1,nrow=n)-Hg)%*%y
sigma2=1/rgamma(t,(nu_0+n)/2,(nu_0*sigma2_0+SS)/2)
Vb=g*solve(t(X)%*%X)/(g+1)
Eb=Vb%*%t(X)%*%y
E=matrix(rnorm(t*m,0,sqrt(sigma2)),t,m)
b=t(t(E%*%chol(Vb))+c(Eb))
for (i in 1:m) 
{
  print(quantile(b[,i],c(0.025, 0.975)))
}         
plot(density(b[,2]),col=1, xlab=names(azdiabetes)[1],main ="Problem 1 Part a")
abline(v=mean(b[,2]),col=2)
abline(v=quantile(b[,2],0.025),col=4)
abline(v=quantile(b[,2],0.975),col=4)
for (i in 2:6) 
{
  plot(density(b[,i+1]),col=1, xlab=names(azdiabetes)[i+1],main ="Problem 1 Part a")
  abline(v=mean(b[,i+1]),col=2)
  abline(v=quantile(b[,i+1],0.025),col=4)
  abline(v=quantile(b[,i+1],0.975),col=4)
}

#Part b
lpy.X=function(y, X) 
{
  n=nrow(X)
  m=ncol(X)
  g=length(y)
  nu_0=1
  sigma2_0=try(summary(lm(y ~ -1+X))$sigma^2, silent = TRUE)
  if (m == 0) 
  {
    Hg=0
    sigma2_0=mean(y^2)
  }
  else if (m > 0) 
  {
    Hg=(g/(g+1))*X%*%solve(t(X)%*%X)%*%t(X)
  }
  SS=t(y)%*%(diag(1, nrow = n)-Hg)%*%y
  -(1/2)*(n*log(pi)+m*log(1+g)+(nu_0+n)*log(nu_0*sigma2_0+SS)-nu_0*log(nu_0*sigma2_0))+
    lgamma((nu_0+n)/2)-lgamma(nu_0/2)
}
t=1000
A=matrix(NA, t, m)
B=matrix(0, t, m)
a=rep(1, m)
lpy.c=lpy.X(y, X[, a == 1, drop = FALSE])
for(i in 1:t) 
{
  for(j in sample(1:m)) 
  {
    temp=a
    temp[j]=1-temp[j]
    lpy.m=lpy.X(y, X[, temp == 1, drop = FALSE])
    r=(lpy.m-lpy.c)*(-1)^(temp[j] == 0)
    a[j]=rbinom(1, 1, 1/(1+exp(-r)))
    if (a[j] == temp[j]) 
    {
      lpy.c=lpy.m
    }
    A[i, ]=a
    Hg=(g/(g+1))*X[, A[i,] == 1, drop = FALSE]%*%solve(t(X[, A[i,] == 1, drop = FALSE])%*%X[, A[i,] == 1, drop = FALSE])%*%t(X[, A[i,] == 1, drop = FALSE])
    SS=t(y)%*%(diag(1, nrow = n)-Hg)%*%y
    sigma2=1/rgamma(t, (nu_0+n)/2, (nu_0*sigma2_0+SS)/2)
    Vb=g*solve(t(X[, A[i,] == 1, drop = FALSE])%*%X[, A[i,] == 1, drop = FALSE])/(g+1)
    Eb=Vb%*%t(X[, A[i,] == 1, drop = FALSE])%*%y
    E=matrix(rnorm(sum(A[i, ]), 0, sqrt(sigma2)), 1, sum(A[i, ]))
    B[i, A[i,] == 1]=t(t(E%*%chol(Vb))+c(Eb))
  }
}
for (i in 1:m) 
{
  print(sum(B[,i]!= 0)/t)
}
for (i in 1:m) 
{
  print(quantile(B[,i],c(0.025, 0.975)))
}         
plot(density(B[,2]),col=1, xlab=names(azdiabetes)[1],main ="Problem 1 Part b")
abline(v=mean(B[,2]),col=2)
abline(v=quantile(B[,2],0.025),col=4)
abline(v=quantile(B[,2],0.975),col=4)
for (i in 2:6) 
{
  plot(density(B[,i+1]),col=1, xlab=names(azdiabetes)[i+1],main ="Problem 1 Part b")
  abline(v=mean(B[,i+1]),col=2)
  abline(v=quantile(B[,i+1],0.025),col=4)
  abline(v=quantile(B[,i+1],0.975),col=4)
}