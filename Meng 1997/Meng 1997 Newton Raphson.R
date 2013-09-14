# Meng 1997
# Th first part of this script does NOT impose the constraint!

# The data
data=list(x=c(0,1,2,3,4,5),nx=c(168,32,16,6,1,0))
xbar=sum(data$x*data$nx)/sum(data$nx[-1])
n=sum(data$nx[-1])
as.character(xbar)

# The log-likelihood expression. It needs to be an "expression" so that we can get R to find the derivative wrt theta symbolically.
ll.expr <- expression(n*(xbar*log(lambda)-lambda-log(1-exp(-lambda))))

# Finding the first derivative of log-likelihood symbollically
dll.expr <- D(ll.expr,"lambda")
dll.expr 

# Finding the second derivative of the log-likelihood symbollically
ddll.expr <- D(dll.expr,"lambda")
ddll.expr 

# Creating functions which can be evaluated out of our expressions for the symbolic derivatives
f=function(lambda) eval(dll.expr) # first derivative
f.=function(lambda) eval(ddll.expr) # second derivative


# Code to do Newton
# start is the starting value
# tol is a tolerance term that is used to stop Newton once the absolute difference between consecutive steps is smaller than tol
# f is the function we are looking for a root of
# f. is the first derivative of the above function

newtRap = function(start,tol,f,df){
  sol=c(start)
  check=1
  while(check>tol){ 
    xn = sol[length(sol)] - f(sol[length(sol)])/df(sol[length(sol)]) 
    sol=c(sol,xn) 
    check=abs(sol[length(sol)]-sol[length(sol)-1])
  }
  return(sol) 
}

# Finding the MLE for lambda
# A "decent" starting value. Warning a starting value of greater than about 2.157 cause the Newton step to go out of the parameter space!
start1 = 0.0000001
z=newtRap(start1,10^-20,f,f.) # The steps can be very small so the tolerance needs to be very small
z # the step values for Newton

# # The negative log-likelihood function and log-likelihood function.
ll=function(lambda) n*(xbar*log(lambda)-lambda-log(1-exp(-lambda)))
y=ll(z) 
y # the log-likelihood values


# A plot
xrep=seq(10^-20,2,0.00001)
yrep=ll(xrep) 
plot(xrep,yrep,type="l",xlab=expression(lambda),ylab = "Log-likelihood")

# Adding the Newton steps to the plot
points(z,y,pch=4)
zeroto6=0:(length(z)-1)
library(calibrate)
textxy(z,y,zeroto6,cx=0.7)



# Meng's functions
# g = function(lambda) lambda -xbar*(1-exp(-lambda))
# g. = function(lambda) 1 - xbar*exp(-lambda)
# z=newton(0.4,10^-5,g,g.)
# z

# I think Meng's algorithm fails because of some of the simplifications he has made. Not that his simplifications are incorrect, but that not simplifying helps to "protect" against a failure. 


# PART 2
# I now impose the constraint by a log transform

# The log-likelihood expression. It needs to be an "expression" so that we can get R to find the derivative wrt theta symbolically.
ll.expr.c <- expression(n*(xbar*theta-exp(theta)-log(1-exp(-exp(theta)))))


# The first derivative of log-likelihood
dll.expr.c <- D(ll.expr.c,"theta")
dll.expr.c 

# The second derivative of the log-likelihood
ddll.expr.c <- D(dll.expr.c,"theta")
ddll.expr.c #It is a real mess and looks like it would be tedious to do by hand!

# Creating functions which can be evaluated out of our expressions for the symbolic derivatives
cf=function(theta) eval(dll.expr.c) # first derivative
cf.=function(theta) eval(ddll.expr.c) # second derivative


# Finding the MLE for lambda
# A "decent" starting value. Warning a starting value of less than about 0.005 causes a very large first Newton step and can exceed R's limits and we get NaNs!
start = log(0.4)
cx = newtRap(start,10^-15,cf,cf.) # The steps can be very small so the tolerance needs to be very small
cx
cz = exp(cx) # the step values for Newton
cz

# # The negative log-likelihood function and log-likelihood function.
cy=ll(cz) 
cy # the log-likelihood values


xrep=seq(0.01,2,0.001)
yrep=ll(xrep) 
plot(xrep,yrep,type="l",xlab=expression(lambda),ylab = "Log-likelihood")

# Adding the Newton steps to the plot
points(cz,cy,pch=4)
zeroton=0:(length(cz)-1)
library(calibrate)
textxy(cz,cy,zeroton,cx=1.1)

