# Newton--Raphson
newtRap <- function(start,tol,f,df,details=TRUE){
  sol <- c(start)
  check <- 1
  while(check>tol){
    xn <- sol[length(sol)] - f(sol[length(sol)])/df(sol[length(sol)])
    sol <- c(sol,xn)
    check <- abs(sol[length(sol)]-sol[length(sol)-1])
if(details){print(xn)}
  }
  return(sol)
}


# The log-likelihood with the constraint imposed by a log transform
ll.expr.c <- expression(n*(xbar*theta-exp(theta)-log(1-exp(-exp(theta)))))

# Finding the derivatives symbollically
dll.expr.c <- D(ll.expr.c,"theta") # The first derivative
ddll.expr.c <- D(dll.expr.c,"theta") # The second derivative

g <- function(theta) eval(dll.expr.c) # Making a function out of the first derivative
dg <- function(theta) eval(ddll.expr.c) # Making a function out of the second derivative

# The "data"
n <- 55
xbar <- 1.56363636363636

# Finding the MLE
theta0 <- log(100)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])

theta0 <- log(1.56)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])

theta0 <- log(0.93)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])

theta0 <- log(0.4)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])

theta0 <- log(0.1)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])

theta0 <- log(0.01)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])


theta0 <- log(0.005)
theta.est <- newtRap(theta0,10^-6,g,dg,details=FALSE)
exp(theta.est[length(theta.est)])


# Finding the MLE for lambda
# A "decent" starting value. Warning a starting value of less than about 0.005 causes a very large first Newton step and can exceed R's limits and we get NaNs!
start = log(0.01)
cx = newtRap(start,10^-15,g,dg,details=FALSE) # The steps can be very small so the tolerance needs to be very small
cz = exp(cx) # the step values for Newton

# # The negative log-likelihood function and log-likelihood function.
ll=function(lambda) n*(xbar*log(lambda)-lambda-log(1-exp(-lambda)))
cy=ll(cz) 

xrep=seq(0.01,2,0.001)
yrep=ll(xrep) 
plot(xrep,yrep,type="l",xlab=expression(lambda),ylab = "Log-likelihood")

# Adding the Newton steps to the plot
points(cz,cy,pch=4)
zeroton=0:(length(cz)-1)
library(calibrate)
textxy(cz,cy,zeroton,cx=1.1)

########################################
## What if we ignore the constraint? ###
########################################

# The log likelihood without the transform
ll.expr <- expression(n*(xbar*log(lambda)-lambda-log(1-exp(-lambda))))

# Derivatives
dll.expr <- D(ll.expr,"lambda")
ddll.expr <- D(dll.expr,"lambda")
f <- function(lambda) eval(dll.expr) 
df <- function(lambda) eval(ddll.expr) 


# Finding the MLE
lambda0 <- 100
newtRap(lambda0,10^-6,f,df,details=FALSE)
# Does not work! The step is too far!

lambda0 <- 1.56
newtRap(lambda0,10^-6,f,df,details=FALSE)

lambda0 <- 0.93
newtRap(lambda0,10^-6,f,df,details=FALSE)

lambda0 <- 0.4
newtRap(lambda0,10^-6,f,df,details=FALSE)

lambda0 <- 0.1
newtRap(lambda0,10^-6,f,df,details=FALSE) 

lambda0 <- 0.01
newtRap(lambda0,10^-6,f,df,details=FALSE)

# What about an extreme starting value

lambda0 <- 10^-15
newtRap(lambda0,10^-15,f,df,details=FALSE)
# Seems to work!