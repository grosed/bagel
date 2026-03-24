
# type
setClass("bagel_type", slots=list(particle_constructor = "function",
                                  prior = "function",
				  H = "function",
				  particles = "environment",
				  t = "numeric"))

# constructor
bagel <- function(particle_constructor,prior,H)
{
   return(new("bagel_type", particle_constructor = particle_constructor, prior = prior, H = H, particles = dict(),t = 1))
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












	




