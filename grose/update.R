
setMethod("update","particle_type",
          function(object,y)
	  {
	    y <- matrix(c(y),1,1)
	    mu <- object@mu
	    sigma <- object@sigma
	    H <- object@H
	    I <- diag(1)
	    e <- y - t(H) %*% mu
	    Q <- (t(H) %*% sigma %*% H) + I
	    Q <- Q[1,1] # treat Q as a scalar
	    A <- sigma %*% H%*% solve(Q)
	    sigma <- sigma %*% A %*% t(A) * Q
	    mu <- mu + A %*% e 
	    return(particle(mu,sigma,H))
	  })