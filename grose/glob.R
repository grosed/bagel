library(bagelR)

# type
setClass("bagel_type", slots=list(H = "function",
      		       		  prior = "function",
				  p0 = "numeric",
				  p = "numeric",
				  s = "numeric",
				  n = "numeric",
				  bagel_object = "Rcpp_bagelR"))

# constructor
bagel <- function(H,prior,p0,p,s,n)
{

return(new("bagel_type",H=H,prior=prior,p0=p0,p=p,s=s,n=n,bagel_object=new(bagelR,p0,p,s,n)))	
}


# update
setMethod("update","bagel_type",
          function(object,y)
	  {
	    bagel_object <- object@bagel_object 
	    H <- object@H
	    prior <- object@prior
            taus <- bagel_object$get_taus()
  	    t <- bagel_object$get_time()
  	    taus <- c(taus,t-1)
  	    bagel_object$set_feature_vectors(taus,Map(function(tau) return(H(t,tau)),taus))
  	    ts <- c(taus[-1],t)				  
  	    priors <- Map(prior,ts)
  	    prior_mus <- Map(function(x) return(x$mu), priors)
  	    prior_sigmas <- Map(function(x) x$sigma, priors)		
  	    bagel_object$set_priors(ts,prior_mus,prior_sigmas)
	    return(bagel_object$update(x))
	  })

# weights
setGeneric("weights",function(object) standardGeneric("weights"))
setMethod("weights","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_weights())
	  })

# taus
setGeneric("taus",function(object) standardGeneric("taus"))
setMethod("taus","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_taus())
	  })

# time
setGeneric("time",function(object) standardGeneric("time"))
setMethod("time","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_time())
	  })







H <- function(t,tau)
{
  ## taking from p.fearnhed bagel.R  ###
  if(tau == 0)
  {
     M <- matrix(c(0),nc=1,nr=2)	  
     M[1,1] <- 1
     M[2,1] <- t
     return(M)
  }
  else
  {
    M <- matrix(c(0),nc=1,nr=4)	  
    M[1,1] <- 1
    M[2,1] <- t
    if(t > tau)
    {
      M[3,1] <- 1
      M[4,1] <- t
    }
  return(M)
  }
}


prior <- function(t)
{
  ## taking from p.fearnhed bagel.R  ###
  if(t == 1)
  {
     mu <- matrix(c(0),nc=1,nr=2)
     sigma <- diag(2)
     return(list("mu" = mu,"sigma" = sigma))
  }
  mu <- matrix(c(0),nc=1,nr=4)
  sigma <- diag(4)
  sigma[3,3] <- (t-1)*(t-1)
  sigma[4,3] <- sigma[3,4] <- -(t-1)
  return(list("mu" = mu,"sigma" = sigma))
}



library(bagelR)

p <- 1.0;
p0 <- 0.9;
s <- 1.0;
n <- 1000

blob <- bagel(H,prior,p0,p,s,n)
dat <- as.numeric(read.csv("./data/example_2.dat",header=FALSE)[[1]])
w0t <- c()
for(x in dat)
{
  w0t <- c(w0t,update(blob,x))
}


library(bagelR)

p <- 1.0;
p0 <- 0.9;
s <- 1.0;
n <- 1000

blob <- new(bagelR,p0,p,s,n)

dat <- as.numeric(read.csv("./data/example_2.dat",header=FALSE)[[1]])

start.time <- Sys.time()
w0t <- c()
for(x in dat)
{
  taus <- blob$get_taus()
  t <- blob$get_time()
  taus <- c(taus,t-1)
  #cat("t is : ",t,"\n")
  #cat("number of particles is : ",length(taus) + 1,"\n")
  # cat(taus,"\n")
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




