# MODULE 3 HOMEWORK

# PROBLEM 2

X_range <- c(0,1,2,3)
f.x <- function (x) (2^x/factorial(x)) * exp(-2)
f_x <- function (x) f.x(x)*(x %in% X_range)
sum(f_x(X_range))

   