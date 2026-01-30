
# type
setClass("particle_type", slots=list(prior = "function", H = "function", tau = "integer", post.mu = "matrix", post.sigma = "matrix", weight = "numeric", p0 = "numeric", p = "numeric"))
# constructor
particle <- function(prior,H,tau,p0,p)
{   
   return(new("particle_type", prior = prior, H = H, tau = tau, p0 = p0, p = p, post.mu = matrix(), post.sigma = matrix(), weight = 0.0))
}

# type
setClass("particle_kv_type", slots=list(s = "numeric"),contains="particle_type")
# constructor
particle_kv <- function(prior,H,tau,p,p0,s)
{
   return(new("particle_type", prior = prior, H = H, tua = tau, p0 = p0, p = p, post.mu = matrix(), post.sigma = matrix(), weight = 0.0, s = s))
}


# theorem_1
setGeneric("theorem_1",function(object.1,object.t,t) standardGeneric("theorem_1"))
setMethod("theorem_1",c("particle_type","particle_type","integer"),
          function(object.1,object.t,t)
	  {
	     # NOTE - IMPLEMENTED ONLY FOR UNIFORM FOR NOW
	     if(object.1@tau != 0) # has to be intial object
	     {
	       # throw an error !!
	     }
	     if(object.t@tau == object.1@tau) # initial particle
	     {
	       if(t == 1L)
	       {
		  weight <- object.1@p0  
	       }
	       else
	       {
	          weight <- object.t@weight
	       }
	       return(new(class(object.t)[1], prior = object.t@prior, H = object.t@H, tau = object.t@tau,p0 = object.t@p0, p = object.t@p,
	                                      post.mu = object.t@post.mu, post.sigma = object.t@post.sigma, weight = weight))
             }
	     if(object.t@tau == t - 1L) # latest particle
	     {
	       weight <- object.1@weight*(1.0 - object.t@p0)/(object.t@p0*(t-1L))
	     }
	     else 
	     {
	       if(t == 2)
	       {
	         weight <- (1.0 - object.t@p0) #
	       }
	       else
	       {
	         weight <- object.t@weight*(t-2L)/(t-1L)
	       }
	     }
             return(new(class(object.t)[1], prior = object.t@prior, H = object.t@H, tau = object.t@tau,p0 = object.t@p0, p = object.t@p,
	                                  post.mu = object.t@post.mu, post.sigma = object.t@post.sigma, weight = weight))
	  })



# theorem_2
setGeneric("theorem_2",function(object,t) standardGeneric("theorem_2"))
setMethod("theorem_2",c("particle_type","integer"),
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

             return(new(class(object)[1], prior = object@prior, H = object@H, tau = t - 1L,p0 = object@p0, p = object@p,
	                                  post.mu = post.mu.t, post.sigma = post.sigma.t, weight = object@weight))
	  })

# theorem_3
setGeneric("theorem_3",function(object,t,y) standardGeneric("theorem_3"))
setMethod("theorem_3",c("particle_type","integer","numeric"),
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
	    return(new(class(object)[1], prior = object@prior, H = object@H, tau = object@tau, p0 = object@p0, p = object@p,
	                                 post.mu = mu, post.sigma = sigma, weight = object@weight))
	  })



# theorem_4
setGeneric("theorem_4",function(object,t,y) standardGeneric("theorem_4"))
setMethod("theorem_4",c("particle_type","integer","numeric"),
          function(object,t,y)
	  {
             #if(object@tau == 0)
	     #{
	     #   return(object)
	     #}
             sigma.post <- object@post.sigma
	     mu.post <- object@post.mu
	     
	     H.t.tau <- object@H(t,object@tau)
	     sigma <- 1.0 # just testing - this needs to be in the devided particle - need more info on this
	     var.pred <- sigma^2 * (1.0 + t(H.t.tau) %*% sigma.post %*% H.t.tau)
	     weight <- object@weight * dnorm(y,t(H.t.tau) %*% mu.post,sqrt(var.pred))
	     
             return(new(class(object)[1], prior = object@prior, H = object@H, tau = object@tau,p0 = object@p0, p = object@p,
	                                  post.mu = object@post.mu, post.sigma = object@post.sigma, weight = weight))
	  })

