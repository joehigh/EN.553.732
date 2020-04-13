#Problem 1
n = 10000;
x= numeric(n);
f= numeric(n);
g= numeric(n);
x= rnorm(n);
g= dnorm(x);
for (i in 1:n) {
  y= runif(1, 0, 1);
  if (y<= 0.3) 
  {
    f[i] = dbeta(x[i], 5, 2);
  } 
  else 
  {
    f[i] = dbeta(x[i], 2, 8);
  }
}
expected.value <- sum(x*f/g)/sum(f/g)

#Computing the probability that the random variable is in the interval (0.35,0.55)
for (i in 1:n) {
  temp <- runif(1, 0, 1);
  if (temp <= 0.3) {
    f[i] = rbeta(1, 5, 2);
  } else {
    f[i] = rbeta(1, 2, 8);
  }
}
prob <- sum(f >= 0.35 & f <= 0.55)/n

expected.value
prob
