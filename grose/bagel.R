library(collections)
source("particle.R")

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

	     theorem_2 <- function(t)
	     {
		if(t == 1) # this is the first particle
		{
			prior <- object@prior(t)
     			H <- object@H(t,t)
     			# create the particle
     			particle <- object@particle_constructor(prior$mu,prior$sigma,H)
			return(particle)		   			
		}
		# create a new particle from the first
		# TODO
	     }


             theorem_3 <- function(t,tau)
             {
	        # THEOREM 3
     		particle <- update(particle,y)
     		# add it to the population
     		object@particles$set(object@t,particle)
		return(object)
	     }


	     object@particles$set(t) <- theorem_2(t)

	     for(tau in 1:t)
	     {
	        object@particles$set(tau) <- the
	     }


   	     if(object@t == 1)
   	     {
	        
	        return(theorem_3())
             }
	

	  })


