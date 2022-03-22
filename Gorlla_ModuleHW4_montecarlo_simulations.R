# MODULE 4 HOMEWORK

## PROBLEM 1
#Suppose that for certain microRNA of size 20 the probability of a purine is binomially distributed with probability 0.7. Say there are 100 such microRNAs,
# each independent of the other.
#Let Y denote the average number of purine in these microRNAs. Find the probability that Y is great than 15. Please give a theoretical calculation, do NOT
# use Monte Carlo simulation to approximate. Show all the steps and formulas in your calculation.

# P(15<X<=20)= P(X<=20)-P(X<=15)
pnorm(20,mean=14,sd=sqrt(4.2)/sqrt(100))-pnorm(15,mean=14,sd=sqrt(4.2)/sqrt(100)) 
# output : [1] 5.317746e-07

## PROBLEM 2
# Two genes’ expression values follow a bivariate normal distribution. Let X and Y denote their expression values respectively. Also, assume that X has mean=9 and
# variance=3.Y has mean=10 and variance=5. The covariance between X and Y is 2.
# In a trial, 50 independent measurements of the expression values of the two genes are collected, We wish to find the probability
# , i.e., the probability that the sample mean for the second gene exceeds the sample mean of the first gene more than 0.5. 
# Conduct a Monte Carlo simulation to approximate this probability, providing a 95% confidence interval for your estimation. 
# Submit your R script for the Monte Carlo simulation, and a brief summary of the actual simulation results.

#install and load the package
install.packages("mvtnorm")
library(mvtnorm)
nsim = 1000
XmeanLess_sim <- rep(NA, nsim)
for (i in 1:nsim) {
  data.sim<-rmvnorm(50, mean=c(9, 10), sigma=matrix(c(3,2,2,5), nrow=2))
  mean.sim<-apply(data.sim, 2, mean) 
  Xmean<-mean.sim[1]
  Ymean<-mean.sim[2]
  XmeanLess_sim[i] <- (Xmean+0.5<Ymean)
}
# Mean of all the Monte Carlo simulations 
mean(XmeanLess_sim) 
# output: [1] 0.963
# 95% confidence interval for  estimation
mean(XmeanLess_sim) + c(-1,1)*1.96*sqrt(var(XmeanLess_sim)/10000)
# The mean lies between the interval 
# output: [1] 0.9592984 0.9667016 
# The mean 0.963 lies between the interval  0.9592984 0.9667016 


## PROBLEM 3
# X1 - chi-square distribution
X1 <- rchisq(10000, df = 10)
# X2 - gamma distribution
X2 <- rgamma(10000, shape = 1, rate = 2)
# X3 - t-distribution
X3 <- rt(10000, 3)
# Defining a new random variable Y using X1, X2. X3
Y <- (sqrt(X1)*X2) + 4*(X3^2)
# Monte Carlo simulation to find the mean of Y
mean(Y)
# output : [1] 14.23247
# summary of the actual simulation results
Ex1 <- 10 # mean of X1 
Ex2 <- 2 # mean of X2
Ex3 <- 0 # mean of x3
Ey <- sqrt(EX1)*EX2 + 4*(EX3^2)
Ey
# output: [1] 6.324555
# The empirical mean 14.23247 from the Monte Carlo simulation is not close to actual simulation mean 6.324555


## PROBLEM 4
# f(x)=exp(-x)exp(-exp(-x))
n <- 10000
an <- sqrt(2*log(n)) - 0.5*(log(log(n))+log(4*pi))*(2*log(n))^(-1/2)
bn <- (2*log(n))^(-1/2)
e <- double()
for (i in 1:1000) e[i] <- (max(rnorm(n))-an)/bn
# plot 1
plot(density(e),ylim=c(0,0.5))
# function 
f<-function(x){exp(-x)*exp(-exp(-x))}
# plot 2
curve(f,range(density(e)$x),add=TRUE,col = "blue")
# plot 3
curve(dnorm,add=TRUE,col = "red")
legend("topleft", legend=c("Normal Distribution","curve for the function 'f'", "Extreme Value"), lty=c(1,12,2), lwd=c(2,2,2),col=c( "red", "blue", "black"))
# Extreme value investigation. The black line (extreme value) fits to the blue line ( curve for the function "f" density of generated extreme data) much better than the red line (normal distribution)





