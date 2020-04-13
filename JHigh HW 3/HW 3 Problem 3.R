#Problem 3
x=c(2.3656491, 2.4952035, 1.0837817, 0.7586751, 0.8780483, 1.2765341, 1.4598699,
      0.1801679, -1.0093589, 1.4870201, -0.1193149, 0.2578262)
n = 100000
set.seed(123)
mu= numeric(n)
tau= numeric(n)
p=numeric(n)
mu[1] = rbeta(1, 2, 2)
tau[1] = rlnorm(1, 1, 10)
p[1] = 1/(sqrt(tau[1]*2*pi)^length(x))*exp(-sum((x-mu[1])^2)/(2*tau[1]))*dbeta(mu[1], 2, 2)*dlnorm(tau[1], 1, 10);
for (i in 1:n) {
  mu_star = rbeta(1, 2, 2)
  tau_star = rlnorm(1, 1, 10);
  p_star = 1/(sqrt(tau_star*2*pi)^length(x))*exp(-sum((x-mu_star)^2)/(2*tau_star))*dbeta(mu_star, 2, 2)*dlnorm(tau_star, 1, 10)
  theta = min(p_star/p[i], 1)
  U = runif(1)
  if(U < theta) {
    mu[i+1] = mu_star
    tau[i+1] = tau_star
    p[i+1] = p_star
    }
  else {
    mu[i+1] = mu[i]
    tau[i+1] = tau[i]
    p[i+1] = p[i]
    }
  }
#posterior prob of mu>=0.5
sum(mu >= 0.5)/n
hist(mu[1000:n+1])
plot(mu[1000:n+1], type='l', main="Trace plot for mu", xlab='Simulation number')
hist(tau[1000:n+1])
plot(tau[1000:n+1],type='l', main="Trace plot for tau", xlab='Simulation number')
plot(mu[1000:n+1], tau[1000:n+1], type='l', main='Trace plot for mu and tau')
acf(mu)
acf(tau)