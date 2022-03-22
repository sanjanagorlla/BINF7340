# MODULE-5 

## problem 1
## numerical optimization of likelihood
# A random sample of size 6 from the exp(λ) distribution results in observations:  
#1.636, 0.374, 0.534, 3.015, 0.932, 0.179.
# Find the MLE on this data set in two ways: by numerical optimization of the likelihood (please include R code) and by the analytic formula.

observations<- c(1.636, 0.374, 0.534, 3.015, 0.932, 0.179);
nloglik <- function(x,theta) -sum(dexp(x = x, rate = theta, log = TRUE))
# minimizing nloglik using optimize() and extracting the values through minimum
optimize(f = nloglik,x = observations,interval = c(0,5))$minimum

#0.8995525

## Analytical formula
ana_form_exp <- 1/mean(observations)
ana_form_exp

#  0.8995502



### problem 2
# A random sample X1, X2, ………., X53 follows chi-square distribution with m degree of freedom, has sample mean  and sample standard deviation . 
# (a) Find the point estimator of m using the method of moments.
# (b) Find a one-sided 90% lower confidence interval of m. 
# Please provide the formulas and the derivations together with your numerical answer.
# (a) 100.8
# (b)
# problem 2
100.8+qt(0.1,52)*12.4/sqrt(53)
# [1] 98.58908


## problem 3
# On the Golub et al. (1999) data set, analyze the Zyxin gene expression data separately for the ALL and AML groups. 
    # (a) Find the bootstrap 95% CIs for the mean and for the variance of the gene expression in each group separately.
    # (b) Find the parametric 95% CIs for the mean and for the variance of the gene expression in each group separately. (You need to choose the appropriate approximate formula to use: z-interval, t-interval or chi-square interval.)
   #  (c)  Find the bootstrap 95% CI for the median gene expression in both groups separately.
    # (d)  Considering the CIs in parts (a)-(c), does the Zyxin gene express differently in ALL and AML patients?

######## 3(a)##########

##  Bootstrap 95% CIs for MEAN of ALL patients
data(golub, package = "multtest")
Zyxin<-golub[2124,]  #take the 2124th row, Zyxin gene
gol.fac<-factor(golub.cl, levels = 0:1, labels = c("ALL","AML")) 
x_all <- golub[2124,gol.fac=="ALL"]
#sample size n
n_all<-length(x_all)
nboot<-1000    #number of bootstrap runs
boot.xbar <- rep(NA, nboot) #empty vector to save bootstraped statistic (mean)
for (i in 1:nboot) { #do nboot iterations
  data.star <- x_all[sample(1:n_all,replace=TRUE)] #random sample n observations
  boot.xbar[i]<-mean(data.star) #i-th bootstraped mean, save in boot.xbar
}
quantile(boot.xbar,c(0.025,0.975))  #find quantiles of the bootstraped means

# -0.58522281 -0.03523272 


## Bootstrap 95% CIs for VARIANCE of ALL patients
boot.xbar <- rep(NA, nboot) #empty vector to save bootstraped statistic (mean)
for (i in 1:nboot) { #do nboot iterations
  data.star <- x_all[sample(1:n_all,replace=TRUE)] #random sample n observations
  boot.xbar[i]<-var(data.star)   #i-th bootstraped mean, save in boot.xbar
}
quantile(boot.xbar,c(0.025,0.975))  #find quantiles of the bootstraped means

# 0.3334907 0.6511065 


## Bootstrap 95% CIs for VARIANCE of AML patients
x_aml <- golub[2124,gol.fac=="AML"]
#sample size n
n_aml<-length(x_aml)
nboot<-1000    #number of bootstrap runs
boot.xbar <- rep(NA, nboot) #empty vector to save bootstraped statistic (mean)
for (i in 1:nboot) { #do nboot iterations
  data.star <- x_aml[sample(1:n_aml,replace=TRUE)] #random sample n observations
  boot.xbar[i]<-mean(data.star)  #i-th bootstraped mean, save in boot.xbar
}
quantile(boot.xbar,c(0.025,0.975))  #find quantiles of the bootstraped means

# 1.379465 1.810960 

## Bootstrap 95% CIs for VARIANCE of AML patients
boot.xbar <- rep(NA, nboot) #empty vector to save bootstraped statistic (mean)
for (i in 1:nboot) { #do nboot iterations
  data.star <- x_aml[sample(1:n_aml,replace=TRUE)] #random sample n observations
  boot.xbar[i]<-var(data.star)   #i-th bootstraped mean, save in boot.xbar
}
quantile(boot.xbar,c(0.025,0.975))  #find quantiles of the bootstraped means

# 0.05250602 0.20121013 


######## 3(b)#############################
## The parametric 95% CI for the MEAN of ALL patients:

x_all <- golub[2124, gol.fac=="ALL"] #take the 2058th row, and columns for ALL patients
n_all<-length(x_all)   #sample size n
v_all = sd(x_all)
ci.all_mean <- mean(x_all)+qt(c(0.025,0.975),df=n_all-1)*sd(x_all)/sqrt(n_all) #95% t-interval
ci.all_mean
# 0.580738750 -0.008846435
# The parametric 95% CI for the VARIANCE of ALL patients:
ci.all_var <-(n_all-1)*(v_all^2)/qchisq(c(0.025, 0.975), n_all-1)
ci.all_var
# # 0.9812951 0.3240441

# The parametric 95% CI for the MEAN of AML patients:
y_aml <- golub[2124,gol.fac=="AML"] #take the 2058th row, and columns for AML patients
n_aml<-length(y_aml)
v_aml = sd(y_aml)
ci.aml_mean <- mean(y_aml)+qt(c(0.025,0.975),df=n_aml-1)*sd(y_aml)/sqrt(n_aml) #95% t-interval
ci.aml_mean
# # 1.339698 1.833638
# The parametric 95% CI for the VARIANCE  of AML patients:
ci.aml_var <-(n_aml-1)*(v_aml^2)/qchisq(c(0.025, 0.975), n_aml-1)
ci.aml_var
# # 0.41621602 0.06597815

###################### 3(c)##############

## Bootstrap 95% CIs for MEDIAN of ALL patients
x_all <- golub[2124,gol.fac=="ALL"]
#sample size n
n_all<-length(x_all)

nboot<-1000

boot.xbar <- rep(NA, nboot)

for (i in 1:nboot) {
  
  data.star <- x_all[sample(1:n_all,replace=TRUE)]
  
  boot.xbar[i]<-median(data.star) #save bootstraped median instead
  
}

quantile(boot.xbar,c(0.025,0.975)) #quantiles of bootstraped median
# -0.73507  0.31432 

## Bootstrap 95% CIs for MEDIAN of AML patients
x_aml <- golub[2124,gol.fac=="AML"]
#sample size n
n_aml<-length(x_aml)

nboot<-1000

boot.xbar <- rep(NA, nboot)

for (i in 1:nboot) {
  
  data.star <- x_aml[sample(1:n_aml,replace=TRUE)]
  
  boot.xbar[i]<-median(data.star) #save bootstraped median instead
  
}

quantile(boot.xbar,c(0.025,0.975)) #quantiles of bootstraped median
# 1.22814 1.82829 

#################### 3(d) ##########
#  Yes，the Zyxin gene expresses differently

# The mean  Zyxin gene expresses differently in ALL and AML patients
##  Bootstrap 95% CIs for MEAN of ALL patients
# [1] -0.58522281 -0.03523272 
## Bootstrap 95% CIs for VARIANCE of ALL patients
# 0.3334907 0.6511065 
##  Bootstrap 95% CIs for MEAN of AML patients
# 1.379465 1.810960 
## Bootstrap 95% CIs for VARIANCE of AML patients
# 0.05250602 0.20121013 


# For parametric, the mean Zyxin gene expresses differently in ALL and AML patients
##  Parametric 95% CIs for MEAN of ALL patients
# -0.580738750 -0.008846435
## parametric  95% CIs for VARIANCE of ALL patients
# 0.9812951 0.3240441
## parametric 95% CIs for MEAN of AML patients
# 1.339698 1.833638
## parametric 95% CIs for VARIANCE of AML patients
# 0.41621602 0.06597815

## The median  of Zyxin gene express differently in ALL and AML patients
## Bootstrap 95% CIs for MEDIAN of ALL patients
# -0.73507  0.31432 
## Bootstrap 95% CIs for MEDIAN of AML patients
# 1.22814 1.82829 


# PROBLEM 4 
# a) Write a R-script to conduct a Monte Carlo study for the coverage probabilities of the two CIs. That is, to generate nsim=1000 such data sets from the Poisson distribution. Check the proportion of the CIs that contains the true parameter λ.
# (b) Run the Monte Carlo simulation for nsim=1000 runs, at three different parameter values: λ=0.1, λ=1 and λ=10. Report the coverage probabilities of these two CIs at each of the three parameter values.
# (c) Considering your result in part (b), which one of these two CI formulas should you use in practice? Can you explain the pattern observed in (b)? 

## 4(a)
MCsim <- function(nsim, lambda) {
  cov1<-cov2<-rep(NA,nsim) # create empty matrices to store data
  for (i in 1:nsim){
    x <- rpois(50, lambda = lambda) # taking 50 random samples from poisson dist.
    xbar <- mean(x) # mean of x
    Xsd <- sd(x) # sd of x
    # finding confidence interval for mean 
    CI1_lower_interval <- xbar + (qt(p = 0.05, df = 49) * sqrt(xbar/50))
    CI1_upper_interval <- xbar + (qt(p = 0.95, df = 49) * sqrt(xbar/50))
    CI1<- c(CI1_lower_interval, CI1_upper_interval)
    # finding confidence interval for variance
    CI2_lower_interval <- (49 * (Xsd^2))/(qchisq(p=0.95, df = 49))
    CI2_upper_interval <- (49 * (Xsd^2))/(qchisq(p=0.05, df = 49))
    CI2<- c(CI2_lower_interval, CI2_upper_interval)
    cov1[i] <- (CI1[1]<lambda)&(lambda<CI1[2])
    cov2[i] <- (CI2[1]<lambda)&(lambda<CI2[2])
  }
  
  print(paste("When lambda=", lambda, ": coverage for first CI is", mean(cov1), ", coverage for second CI is", mean(cov2), "."))}

## 4(b)
# lambda = 0.1
MCsim(1000,0.1)
# [1] "When lambda= 0.1 : coverage for first CI is 0.866 , coverage for second CI is 0.516 ."
# lambda = 1
MCsim(1000,1 )
# [1] "When lambda= 1 : coverage for first CI is 0.906 , coverage for second CI is 0.853 ."
# lambda = 10
MCsim(1000,10 )
# [1] "When lambda= 10 : coverage for first CI is 0.907 , coverage for second CI is 0.877 ."


## 4(c)
# The t-interval should be used in practice
# Considering the results from part (b), it is evident that the CI1 which used t-interval has consistent valuses when compared to cI2 coverages values  which used possion  and differs largely irrespective of the lambda values






