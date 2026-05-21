

setClass("bagel_result_type", slots=list(w0ts = "numeric",weights="numeric",ratios="list",t="numeric"))

bagel_result <- function(w0ts,weights,ratios,t)
{
	return(new("bagel_result_type",w0ts=w0ts,weights=weights,ratios=ratios,t=t))
}




setClass("bagel_uv_approximate_type", slots=list(feature_vector = "function",
      		       			         prior = "function",
				                 transform = "function",
				                 p0 = "numeric",
				                 p = "numeric",
				                 nu = "numeric",
	                                         iota = "numeric",	
				                 max_particles = "numeric",
				                 bagel_object = "Rcpp_bagelR_uv_approximate"))

bagel_uv_approximate <- function(Y,feature_vector,prior,transform,p0,p,nu,iota,max_particles,threshold)
{

   bagel_object <- new("bagel_uv_approximate_type",feature_vector=feature_vector,
						   prior=prior,
			                  	   transform=transform,
			                  	   p0=p0,
			                  	   p=p,
			                  	   nu=nu,
			                  	   iota=iota,
			                  	   max_particles=max_particles,
			                  	   bagel_object=new(bagelR_uv_approximate,p0,p,nu,iota,max_particles))
  return(analyse(Y,threshold,bagel_object))

}

setClass("bagel_uv_exact_type", slots=list(feature_vector = "function",
      		       			   prior = "function",
				           transform = "function",
				           p0 = "numeric",
				           p = "numeric",
				           nu = "numeric",
	                                   iota = "numeric",	
				           max_particles = "numeric",
				           bagel_object = "Rcpp_bagelR_uv_exact"))

bagel_uv_exact <- function(Y,feature_vector,prior,transform,p0,p,nu,iota,max_particles,threshold)
{

   bagel_object <- new("bagel_uv_exact_type",feature_vector=feature_vector,
				             prior=prior,
			                     transform=transform,
			                     p0=p0,
			                     p=p,
			                     nu=nu,
			                     iota=iota,
			                     max_particles=max_particles,
			                     bagel_object=new(bagelR_uv_exact,p0,p,nu,iota,max_particles))
  return(analyse(Y,threshold,bagel_object))
}



setClass("bagel_kv_approximate_type", slots=list(feature_vector = "function",
      		       			         prior = "function",
				                 transform = "function",
				                 p0 = "numeric",
				                 p = "numeric",
				                 sigma = "numeric",
				                 max_particles = "numeric",
				                 bagel_object = "Rcpp_bagelR_kv_approximate"))

bagel_kv_approximate <- function(Y,feature_vector,prior,transform,p0,p,sigma,max_particles,threshold)
{

   bagel_object <- new("bagel_kv_approximate_type",feature_vector=feature_vector,
						   prior=prior,
			                  	   transform=transform,
			                  	   p0=p0,
			                  	   p=p,
			                  	   sigma=sigma,
			                  	   max_particles=max_particles,
			                  	   bagel_object=new(bagelR_kv_approximate,p0,p,sigma,max_particles))
  return(analyse(Y,threshold,bagel_object))

}


setClass("bagel_kv_exact_type", slots=list(feature_vector = "function",
      		       			   prior = "function",
				           transform = "function",
				           p0 = "numeric",
				           p = "numeric",
				           sigma = "numeric",
				           max_particles = "numeric",
				           bagel_object = "Rcpp_bagelR_kv_exact"))

bagel_kv_exact <- function(Y,feature_vector,prior,transform,p0,p,sigma,max_particles,threshold)
{

   bagel_object <- new("bagel_kv_exact_type",feature_vector=feature_vector,
				             prior=prior,
			                     transform=transform,
			                     p0=p0,
			                     p=p,
			                     sigma=sigma,
			                     max_particles=max_particles,
			                     bagel_object=new(bagelR_kv_exact,p0,p,sigma,max_particles))
  return(analyse(Y,threshold,bagel_object))

}



# call

analyse <- function(Y,threshold,bagel_object)
{
   w0ts <- list()
   for(y in Y)
   {
      w0t <- update(bagel_object,y)
      w0ts <- append(w0ts,w0t)
      if(1.0 - w0t > threshold)
      {
	break
      }
   }
   return(bagel_result(as.numeric(w0ts),
                       as.numeric(bagel_object@bagel_object$get_weights()),
                       # as.numeric(bagel_object@bagel_object$get_ratios()),
		       bagel_object@bagel_object$get_ratios(),
                       bagel_object@bagel_object$get_time()))  
}









update <- function(object,y)
{
        bagel_object <- object@bagel_object 
	feature_vector <- object@feature_vector
	prior <- object@prior
	transform <- object@transform
	taus <- bagel_object$get_taus()
	t <- bagel_object$get_time()
	taus <- c(taus,t-1)
	bagel_object$set_feature_vectors(taus,Map(function(tau) return(feature_vector(t,tau)),taus))
	ts <- c(taus[-1],t)
	priors <- Map(prior,ts)
	prior_mus <- Map(function(x) return(x$mu), priors)
	prior_sigmas <- Map(function(x) x$sigma, priors)
	bagel_object$set_priors(ts,prior_mus,prior_sigmas)
	bagel_object$set_transformations(taus,Map(function(tau) return(transform(tau)),taus))
	return(bagel_object$update(y))
}