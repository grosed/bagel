
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
	    return(bagel_object$update(y))
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


