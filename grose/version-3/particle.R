
# type
setClass("particle_type", slots=list(prior = "function", H = "function", tau = "integer", post.mu = "matrix", post.sigma = "matrix"))
# constructor
particle <- function(prior,H,tau)
{
   return(new("particle_type", prior = prior, H = H, tau = tau, post.mu = matrix(), post.sigma = matrix()))
}

# type
setClass("particle_kv_type", slots=list(s = "numeric"),contains="particle_type")
# constructor
particle_kv <- function(prior,H,tau,s)
{
   return(new("particle_type", prior = prior, H = H, tua = tau, post.mu = matrix(), post.sigma = matrix(), s = s))
}


# theorem_2
setGeneric("theorem_2",function(object,...) object)
setMethod("theorem_2","particle_type",
          function(object,t)
	  {
	     if(object@tau != 0)
	     {
	       # throw an error !!
	     }
	     if(t == 1L)
	     {
	       prior.t <- object@prior(t)
	       post.mu.t <- prior.t$mu
	       post.sigma.t <- prior.t$sigma
	     }
	     else
	     {
	        prior.1 <- object@prior(1L)
   		prior.mu.1 <- prior.1$mu
   		prior.sigma.1 <- prior.1$sigma

   		prior.t <- object@prior(t)
   		prior.mu.t <- prior.t$mu
   		prior.sigma.t <- prior.t$sigma

   		d1 <- nrow(prior.mu.1)
   		d2 <- nrow(prior.sigma.1)
		n <- d1 + d2

   		sigma.beta.beta.t <- prior.sigma.t[1:d1,1:d1]
   		sigma.beta.gamma.t <- prior.sigma.t[1:d1,(d1+1):n]
   		sigma.gamma.beta.t <- prior.sigma.t[(d1+1):n,1:d1]
   		sigma.gamma.gamma.t <- prior.sigma.t[(d1+1):n,(d1+1):n]
   		mu.beta.t <- as.matrix(prior.mu.t[1:d1,1])
   		mu.gamma.t <- as.matrix(prior.mu.t[(d1+1):n,1])

   		B <- t(sigma.beta.gamma.t) %*% solve(sigma.beta.beta.t)
   
		# sigma 

   		top_left <- object@post.sigma[1:d1,1:d1]
   		top_right <- B %*% object@post.sigma[1:d1,1:d1]
   		bottom_left <- t(top_right)
   		bottom_right <- sigma.gamma.gamma.t + B %*% (object@post.sigma[1:d1,1:d1]-sigma.beta.beta.t)%*%t(B)

   		post.sigma.t <- rbind(cbind(top_left,top_right),
                                      cbind(bottom_left,bottom_right))
		# mu

   		top <- as.matrix(object@post.mu[1:d1,1])
   		bottom <- mu.gamma.t + B %*% (as.matrix(object@post.mu[1:d1,1]) - mu.beta.t)
		post.mu.t <- rbind(top,bottom)
	     }


             return(new(class(object)[1], prior = object@prior, H = object@H, tau = t - 1L, post.mu = post.mu.t, post.sigma = post.sigma.t))
	  })

# theorem_3
setGeneric("theorem_3",function(object,...) object)
setMethod("theorem_3","particle_type",
          function(object,t,y)
	  {
	    y <- matrix(c(y),1,1)
	    prior.t <- prior(t)
	    mu <- object@post.mu
	    sigma <- object@post.sigma
	    H.t <- H(t,object@tau)
	    I <- diag(1)
	    e <- y - t(H.t) %*% mu
	    Q <- (t(H.t) %*% sigma %*% H.t) + I
	    Q <- Q[1,1] # treat Q as a scalar
	    A <- sigma %*% H.t%*% solve(Q)
	    sigma <- sigma - A %*% t(A) * Q
	    mu <- mu + A %*% e
	    return(new(class(object)[1], prior = object@prior, H = object@H, tau = object@tau, post.mu = mu, post.sigma = sigma))
	  })
