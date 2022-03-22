# MODULE 3 HOMEWORK

# PROBLEM 2

x = 1
f.x <- function (x) (2^x/factorial(x)) * exp(-2)
f_x <- function (x) f.x(x)
sum(f_x(x))
# 0.2706706

X_range <- c(0,1,2,3)
f.x <- function (x) (2^x/factorial(x)) * exp(-2)
f_x <- function (x) f.x(x)*(x %in% X_range)
sum(f_x(X_range))
# 0.8571235

   
#PROBLEM 4
# For Y following a binomial (n = 3, p = 0.25) distribution, compute the following: P(Y ≤ 2) =, E(Y) =, Var(Y) =
set.Y<-c(0,1,2)
sum(dbinom(set.Y,size=3,p=0.25))
# [or]
1-dbinom(3,size=3,p=0.25)
# 0.984375

# E(Y)

EX<-sum(set.Y*dbinom(set.Y,size=3,p=0.25))
EX
# 0.703125

VarX<-sum((set.Y-EX)^2*dbinom(set.Y,size=3,p=0.25))
VarX
# 0.4822655


# PROBLEM 5
pchisq(4,df=3)-pchisq(1,df=3)
# 0.5397878
EX = 3 # mean of chi-square is m 
var = 2*3 # variance of chi-square is 2m 
#The R code for Monte Carlo simulation with sample size  100000
n <- 100000
X <- rchisq(n,df=3)
X1 <-X[X>1 & X<4]
length(X1)/n
# The Monte Carlo estimate is P(1 < X < 4) = 0.5397878.
# The estimated probability is close to the theoretical value 0.53972
# I agree with the estimated value



# PROBLEM 7
# (a) What is the probability that a randomly chosen patient have the Zyxin gene expression values between 1 and 1.6?
pnorm(1.6,1.6,0.4) - pnorm(1,1.6,0.4)
# 0.4331928

# (b) Use a Monte Carlo simulation of sample size n=500,000 to estimate the probability in part (a). Give your R code, and show the value of your estimate.
x=rnorm(500000,1.6,0.4)
y=subset(x,x<=1.6&x>=1)
length(y)/500000
# 0.432752

# (c) What is the probability that exactly 2 out of 5 patients have the Zyxin gene expression values between 1 and 1.6? 
dbinom(2,size=5,p= 0.4331928)
# 0.3417185


# PROBLEM 8
## (a) Hand in a R script that calculates the mean and variance of two random variables X~F(m=2,n=5) and Y~F(m=10,n=5) from their density functions.
# mean of x
integrate(function(x) x*df(x,2,5), lower =-Inf, upper = Inf)
# 1.666667 with absolute error < 2.2e-07

# variance of x
integrate(function(y) (y-1.666667)^2*df(y,2,5),  lower = -Inf, upper = Inf)
# 13.88889 with absolute error < 5.9e-06

# mean of y
integrate(function(x) x*df(x,10,5), lower =-Inf, upper = Inf)
# 1.666667 with absolute error < 2e-08

# variance of y 
integrate(function(y) (y-1.666667)^2*df(y,10,5),  lower = -Inf, upper = Inf)
# 7.222222 with absolute error < 4e-04



## (b) Use the formula in Table 3.4.1 to calculate the means and variances directly.
# mean of x 
n = 5
x_mean = n/(n-2)
x_mean
# 1.666667

# variance of x
n = 5
m =2
variance_x = (2*(n^2)*(m+n-2))/(m*(n-2)^2*(n-4))
variance_x
# 13.88889

# mean of y
n = 5
y_mean = n/(n-2)
y_mean
# 1.666667

# variance of y
n = 5
m =10
variance_x = (2*(n^2)*(m+n-2))/(m*(n-2)^2*(n-4))
variance_x
# 7.222222


# c) X random variable function has : Mean from code = 1.666667, Variance from code = 13.88889, Y random variable function has: Mean from code =  1.666667 , Variance from code = 7.222222
# The mean and variance calculated from the r-code is same as the table calculated values 
# The answers agree with those from part(a) and part(b)



