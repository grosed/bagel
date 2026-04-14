
### NOTE - real time < streamed < online < sequential. where a < b <=> a is more restrictive than b 

### A provisional taxonomy of algorithm type (n = number of data elements processed)
### ----------------------------------------

###                latency     |   storage   |    notes
### real time        O(1)      |    O(1)     |  unbounded data, indefinite persistance, fixed latency, very limited resources (memory, processor speed, power)
### streamed         O(1)      |  O(log(n))  |  large data, long persistancy, fixed latency  
### online          O(log(n))  |     NA      |  moderate data, moderate persistancy 
### sequential       NA        |     NA      |  digital computers are sequential !! how does this term suggest any (proper) subset of all possible algorithms ?

### bagel has O(1) latency and O(n) storage so, by this criterion it is online and (trivially) sequential - but neither streamed or real time !!  


# type
setClass("bagel_type", slots=list(H = "function",
      		       		  prior = "function",
				  p0 = "numeric",
				  p = "numeric",
				  s = "numeric",
				  n = "numeric",
				  bagel_object = "Rcpp_bagelR"))

# constructor
bagel_online <- function(H,prior,p0,p,s,n)
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


# weights
setGeneric("ratios",function(object) standardGeneric("ratios"))
setMethod("ratios","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_ratios())
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
	    return(object@bagel_object$get_time()-1)
	  })



### off line interface


# type
setClass("bagel_results_type", slots=list(bagel_object = "bagel_type",
                                          feature_vector_function = "function",
      		       			  prior_function = "function",
				  	  p0 = "numeric",
				  	  p = "numeric",
				  	  noise_sd = "numeric",
				  	  max_particles = "numeric",
				          threshold = "numeric",
					  w0t = "vector",
					  y = "vector"))

# constructor
bagel_results <- function(bagel_object,
                          feature_vector_function,
                          prior_function,
			  p0,
			  p,
			  noise_sd,
			  max_particles,
			  threshold,
			  w0t,
			  y)
{
return(new("bagel_results_type",bagel_object=bagel_object,
                                feature_vector_function=feature_vector_function,
				prior_function=prior_function,
				p0=p0,
				p=p,
				noise_sd=noise_sd,
				max_particles=max_particles,
				threshold=threshold,
				w0t=w0t,
				y=y))	
}

bagel_offline <- function(feature_vector_function,
	 		  prior_function,
			  p0,
			  p,
			  noise_sd,
			  max_particles,
			  threshold,
			  y)
{
   bagel_object <- bagel_online(feature_vector_function,
                                prior_function,
			        p0,
			        p,
			        noise_sd,
			        max_particles)
   w0t <- c()
   for(yt in y)
   {
	w0 <- update(bagel_object,yt)
   	w0t <- c(w0t,w0) # log the result
	if(1.0 - w0 > threshold) { break }
   }
   
   return(bagel_results(bagel_object,feature_vector,prior_function,p0,p,noise_sd,max_particles,threshold,w0t,y))
}


# weights
setMethod("weights","bagel_results_type",
          function(object)
	  {
	    return(weights(object@bagel_object))
	  })

# taus
setMethod("taus","bagel_results_type",
          function(object)
	  {
	    return(taus(object@bagel_object))
	  })

# time
setMethod("time","bagel_results_type",
          function(object)
	  {
	    return(time(object@bagel_object))
	  })


# weights
setGeneric("weights_tau_zero",function(object) standardGeneric("weights_tau_zero"))
setMethod("weights_tau_zero","bagel_results_type",
          function(object)
	  {
	    return(object@w0t)
	  })


# weights
setMethod("ratios","bagel_results_type",
          function(object)
	  {
	    return(ratios(object@bagel_object))	
	  })

