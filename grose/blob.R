
library(bagelR)

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



blob <- new(bagelR,3.14)

n <- 1000
Hs <- list()	
for(i in 1:n)
{
   taus <- blob$taus()

   for(tau in taus)
   {
      list.append(Hs,matrix(c(1,2,3,4),4,1))
      # list.append(H(1000,tau))
      # blob$feature_vectors(list(matrix(c(1,2,3,4),4,1),matrix(c(1,2,3,4),4,1)))
      # Hs[[tau+1]] <- H(1000,tau)
      # glob <- H(1000,tau)
      # glob <- matrix(c(1,2,3,4),4,1)
   }
    blob$feature_vectors(Hs)	
}




