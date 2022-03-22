# MODULE 1 HOMEWORK


#PROBLEM 1
##(a)
class (vec <-c(5,TRUE)) 
# output [1] "numeric"

##(b)
x <- 1:4
y <- 1:2
x+y
# output [1] 2 4 4 6
# An integer vector with the values 2, 4, 4, 6

##(c)
fsin<-function(x) sin(pi*x)
fsin(1)
# output [1] 1.224647e-16
# The number 1 is returned

##(d)
c(1,2) %*% t(c(1,2))
#     [,1] [,2]
#[1,]   1    2
#[2,]   2    4
# A two by two matrix

##(e)
f <- function(x) {
  g <- function(y) {
    y + z
  }
  z <- 4
  x + g(x)
}
z<-15
f(3)
# output [1] 10
# 10


#PROBLEM 2
# using R to calculate integration

# method 1 (general calculation)
x <- 1:1000 # 1 to 1000
sum(x^2) # squaring and finding the sum
# output [1] 333833500

# Alternatively we can also do 
# method 2 (using integral)
integrand <- function(x) x^2  # writing the integrand function
integrate(integrand, lower = 1, upper = 1000)
# 333333333 with absolute error < 3.7e-06



#PROBLEM 3
##(a)
X <- 2 * seq(1,20)  # a vector X of length 20, with the k th element in X = 2k, for k=1…20
X # Print out the values of X
# output [1]  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40


##(b)
Y <-rep(0,20) # a vector Y of length 20, with all elements in Y equal to 0
Y # Print out the values of Y
# output [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0


##(c)
for ( K in 1:20) {   # for loop, reassigning the value of the K-th element in Y
  if (K < 12) {      # when K < 12
    Y[K] <- cos(3*K) # reassigned as the cosine of 3K
  } else {           # when k is = 12 or k > 12
    Y[K] <- integrate( sqrt, lower=0, upper=K)$value # reassigned as the value of integral and extracting the value 
  }
}
Y # Print out the reassigned values of Y
# output [1] -0.98999250  0.96017029 -0.91113026  0.84385396 -0.75968791  0.66031671 -0.54772926  0.42417901 -0.29213881  0.15425145 -0.01327675 27.71281603 31.24811456 34.92213953 38.72983781 42.66667146 46.72853567 50.91169396 55.21272615 59.62848609








