#Problem 2 
#part a
crime = read.table("crime.dat", header = TRUE)
y = as.matrix(crime[,1])
X = as.matrix(cbind(rep(1, dim(y)[1], 1),crime[,-1]))
g = length(y)
nu_0 = 2
sigma2_0 = 1
n = dim(X)[1]
m = dim(X)[2]
t = 10000
Hg = (g/(g+1))*X%*%solve(t(X)%*%X)%*%t(X)
SSR = t(y)%*%(diag(1, nrow = n)-Hg)%*%y
sigma2 = 1/rgamma(t, (nu_0+n)/2, (nu_0*sigma2_0+SSR)/2)
Vb = g*solve(t(X)%*%X)/(g+1)
Eb = Vb%*%t(X)%*%y
E = matrix(rnorm(t*m, 0, sqrt(sigma2)), t, m)
beta = t(t(E%*%chol(Vb))+c(Eb));
par(mfrow = c(2, 2))
beta_mean=numeric(m)
beta_l=numeric(m)
beta_u=numeric(m)
dev.off()
for (i in 1:m) {
  plot(density(beta[,i]), col = 2, lwd = 1, xlab = names(crime)[i], main = "Problem 2 part a");
  abline(v = quantile(beta[,i], c(0.025, 0.975))[1], col = 3, lwd = 2);
  abline(v = quantile(beta[,i], c(0.025, 0.975))[2], col = 3, lwd = 2);
  abline(v = mean(beta[,i]), col = 4, lwd = 2);
  beta_mean[i]=mean(beta[,i])
  beta_l[i]=quantile(beta[,i], c(0.025, 0.975))[1]
  beta_u[i]=quantile(beta[,i], c(0.025, 0.975))[2]
}
beta_mean
beta_l
beta_u
# Least squares estimates
LS = lm(y ~ M+So+Ed+Po1+Po2+LF+M.F+Pop+NW+U1+U2+GDP+Ineq+Prob+Time, data = crime)
summary(LS)
MSE = sum(LS$residuals^2) / (n - 16);
#upper 95% bound
LS$coefficients +qt(1 - 0.05 / 2, df = n - 16) * sqrt(diag(MSE * solve(t(X) %*% X)))
#lower 95% bound
LS$coefficients - qt(1 - 0.05 / 2, df = n - 16) * sqrt(diag(MSE * solve(t(X) %*% X)))


#Part b (i)

set.seed(1);
tr = sample(1:n, n/2);
te = (1:n)[-tr];
X_tr = X[tr,];
X_te = X[te,];
y_tr = y[tr,];
y_te = y[te,];
beta_ols = solve(t(X_tr)%*%X_tr)%*%t(X_tr)%*%y_tr;
y_ols = X_te%*%beta_ols;
dev.off()
plot(y_te, y_ols, xlab = "y_te", ylab = "y_ols", main='Problem 2 part b (i)', xlim = c(-2, 4), ylim = c(-2, 4), asp = 1);
abline(c(0, 0), c(1, 1), col = "red");
sum((y_ols-y_te)^2)/length(y_te)

#Part b (ii)

g = length(y_tr);
nu_0 = 2;
sigma2_0 = 1;
n = dim(X_tr)[1];
m = dim(X_tr)[2];
t = 10000;
Hg = (g/(g+1))*X_tr%*%solve(t(X_tr)%*%X_tr)%*%t(X_tr);
SSR = t(y_tr)%*%(diag(1, nrow = n)-Hg)%*%y_tr;
sigma2 = 1/rgamma(t, (nu_0+n)/2, (nu_0*sigma2_0+SSR)/2);
Vb = g*solve(t(X_tr)%*%X_tr)/(g+1);
Eb = Vb%*%t(X_tr)%*%y_tr;
E = matrix(rnorm(t*m, 0, sqrt(sigma2)), t, m);
beta_bayes = colMeans(t(t(E%*%chol(Vb))+c(Eb)));
y_bayes = X_te%*%beta_bayes;
plot(y_te, y_bayes, xlab = "y_te", ylab = "y_bayes", main='Problem 2 part b (ii)', xlim = c(-2, 4), ylim = c(-2, 4), asp = 1);
abline(c(0, 0), c(1, 1), col = "red");
sum((y_bayes-y_te)^2)/length(y_te);

#Part c

t1 = 10000;
error_ols = rep(NA, 1, t1);
error_bayes = rep(NA, 1, t1);
for (i in 1:t1) {
  set.seed(i);
  n = 47;
  tr = sample(1:n, n/2);
  te = (1:n)[-tr];
  X_tr = X[tr,];
  X_te = X[te,];
  y_tr = y[tr,];
  y_te = y[te,];
  beta_ols = solve(t(X_tr)%*%X_tr)%*%t(X_tr)%*%y_tr;
  y_ols = X_te%*%beta_ols;
  error_ols[i] = sum((y_ols-y_te)^2)/length(y_te)
  g = length(y_tr);
  nu_0 = 2;
  sigma2_0 = 1;
  n = dim(X_tr)[1];
  m = dim(X_tr)[2];
  t = 10000;
  Hg = (g/(g+1))*X_tr%*%solve(t(X_tr)%*%X_tr)%*%t(X_tr);
  SSR = t(y_tr)%*%(diag(1, nrow = n)-Hg)%*%y_tr;
  sigma2 = 1/rgamma(t, (nu_0+n)/2, (nu_0*sigma2_0+SSR)/2);
  Vb = g*solve(t(X_tr)%*%X_tr)/(g+1);
  Eb = Vb%*%t(X_tr)%*%y_tr;
  E = matrix(rnorm(t*m, 0, sqrt(sigma2)), t, m);
  beta_bayes = colMeans(t(t(E%*%chol(Vb))+c(Eb)));
  y_bayes = X_te%*%beta_bayes;
  error_bayes[i] = sum((y_bayes-y_te)^2)/length(y_te);
}
mean(error_ols)
mean(error_bayes)
