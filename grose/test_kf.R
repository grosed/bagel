library(functional)
library(purrr)

source("kalman_filter.R")


P <- diag(2)
x <- matrix(0,2,1)

#print(P)
#print(x)


kf <- kalman_filter(P,x)

#print(kf$P())
#print(kf$x())
#print(class(kf))


F <- diag(2)
B <- matrix(0,2,1)
Q <- matrix(0,2,2)
u <- matrix(0)

predict <- partial(kalman_filter_predict,F = F,B = B, Q = Q)

predict(kf,u)


H <- diag(2)
R <- matrix(0,2,2)
z <- matrix(0,2,1)



kalman_filter_update(kf,F,B,Q,H,R,u,z)

