
source("particle.R")
source("predict.R")
source("update.R")



## SET-UP FOR EXAMPLE 2





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
  if(tau == 1)
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













