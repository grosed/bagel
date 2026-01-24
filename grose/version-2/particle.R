
# type
setClass("particle_type", slots=list(mu = "matrix", sigma = "matrix"))
particle <- function(mu,sigma)
{
   return(new("particle_type", mu = mu, sigma = sigma))
}

# type
setClass("particle_kv_type", slots=list(s = "numeric"),contains="particle_type")
particle_kv <- function(mu,sigma,s)
{
   return(new("particle_type", mu = mu, sigma = sigma, s = s))
}

# update
setMethod("update","particle_type",
          function(object,H,y)
	  {
	    y <- matrix(c(y),1,1)
	    mu <- object@mu
	    sigma <- object@sigma
	    I <- diag(1)
	    e <- y - t(H) %*% mu
	    Q <- (t(H) %*% sigma %*% H) + I
	    Q <- Q[1,1] # treat Q as a scalar
	    A <- sigma %*% H%*% solve(Q)
	    sigma <- sigma - A %*% t(A) * Q
	    mu <- mu + A %*% e 
	    return(particle(mu,sigma))
	  })
