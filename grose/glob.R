

library(rlist)

H <- function(t,tau)
{
  ## definition of vectors a and b
  ## taking from p.fearnhed bagel.R  ###
  a <- function(t)
  {
    return(matrix(c(1,t),nc=1))
  }

  b <- function(t,tau)
  {
    if(t <= tau)
    { 
      return(matrix(0,nc = 1,nr = 2))
    }
    else
    {
      return(matrix(c(1,t),nc = 1))
    }
  }
  if(tau == 0)
  {
     return(a(t))
  }
  else
  {
     return(rbind(a(t),b(t,tau)))
  }   
}


library(bagelR)

p <- 1.0;
p0 <- 0.9;
s <- 1.0;

blob <- new(bagelR,p0,p,s)

dat <- as.numeric(read.csv("./data/example_2.dat",header=FALSE)[[1]])

w0t <- c()
t <- 1
for(x in dat)
{
  print(t)
  taus <- blob$taus()
  taus <- c(taus,t-1)
  blob$feature_vectors(taus,Map(function(tau) return(H(t,tau)),taus))
  w0t <- c(w0t,blob$update(x))
  t <- t + 1
}


library(bagelR)

p <- 1.0;
p0 <- 0.9;
s <- 1.0;

set.seed(0)
Z <- rnorm(4000,0,1)

blob <- new(bagelR,p0,p,s)



w0t <- c()
t <- 1
for(z in Z)
{
  # print(t)
  taus <- blob$taus()
  taus <- c(taus,t-1)
  blob$feature_vectors(taus,Map(function(tau) return(H(t,tau)),taus))
  w0t <- c(w0t,blob$update(x))
  t <- t + 1
}



set.seed(0)
Z <- rnorm(4000,0,1)

blob <- new(bagelR,p0,p,s)



w0t <- c()
t <- 1
for(z in Z)
{
  # print(t)
  #taus <- blob$taus()
  #taus <- c(taus,t-1)
  #blob$feature_vectors(taus,Map(function(tau) return(H(t,tau)),taus))
  w0t <- c(w0t,blob$update(x))
  #t <- t + 1
}




