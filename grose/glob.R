

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


prior <- function(t)
{

  ## definition of prior
  delta<-c(1,1,1) ###hyper-parameter

  ##definition of matrices in THM2 at time t -- these are for tau = t-1
  sigma.beta.beta <- diag(delta[1:2]) ##does not depend on t or tau
  sigma.gamma.beta <- sigma.beta.gamma <- function(t)
  {
    return(diag(c(0,0)))
  }
  
  sigma.gamma.gamma <- function(t)
  {
    return(matrix(delta[3]*c((t-1)^2,-(t-1),-(t-1),1),nc = 2))
  }
  mu.beta <- matrix(c(0,0),nc = 1)
  ## again this is prior mean for tau = t-1
  mu.gamma <- function(t)
  {
    return(matrix(c(0,0),nc = 1))
  }
  if(t == 1)
  {
    return(list("mu" = mu.beta,"sigma" = sigma.beta.beta))
  }
  else
  {
    return(
    list("mu" = rbind(mu.beta,mu.gamma(t)),
         "sigma" = rbind(cbind(sigma.beta.beta,sigma.beta.gamma(t)),cbind(sigma.gamma.beta(t),sigma.gamma.gamma(t)))
          )
         )
   }
}



library(bagelR)

p <- 1.0;
p0 <- 0.9;
s <- 1.0;

blob <- new(bagelR,p0,p,s)

dat <- as.numeric(read.csv("./data/example_2.dat",header=FALSE)[[1]])

w0t <- c()
for(x in dat)
{
  taus <- blob$get_taus()
  t <- blob$get_time()
  taus <- c(taus,t-1)
  cat("t is : ",t,"\n")
  cat("number of particles is : ",length(taus) + 1,"\n")
  cat(taus,"\n")
  blob$set_feature_vectors(taus,Map(function(tau) return(H(t,tau)),taus))
  ts <- c(taus[-1],t)				  
  priors <- Map(prior,ts)
  prior_mus <- Map(function(x) return(x$mu), priors)
  prior_sigmas <- Map(function(x) x$sigma, priors)		
  blob$set_priors(ts,prior_mus,prior_sigmas)
  w0t <- c(w0t,blob$update(x))
}


library(bagelR)

p <- 1.0;
p0 <- 0.9;
s <- 1.0;

set.seed(0)
Z <- rnorm(3000,0,1)

blob <- new(bagelR,p0,p,s)

w0t <- c()
t <- 1
start.time <- Sys.time()
for(z in Z)
{
  taus <- blob$get_taus()
  t <- blob$get_time()
  taus <- c(taus,t-1)
  #cat("t is : ",t,"\n")
  #cat("number of particles is : ",length(taus) + 1,"\n")
  #cat(taus,"\n")
  blob$set_feature_vectors(taus,Map(function(tau) return(H(t,tau)),taus))
  ts <- c(taus[-1],t)				  
  priors <- Map(prior,ts)
  prior_mus <- Map(function(x) return(x$mu), priors)
  prior_sigmas <- Map(function(x) x$sigma, priors)		
  blob$set_priors(ts,prior_mus,prior_sigmas)
  w0t <- c(w0t,blob$update(x))
}
end.time <- Sys.time()
end.time - start.time


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




