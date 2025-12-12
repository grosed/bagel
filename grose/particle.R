


setClass("particle_type", slots=list(mu = "matrix", sigma = "matrix", H = "matrix"))
particle <- function(mu,sigma,H)
{
   return(new("particle_type", mu = mu, sigma = sigma, H = H))
}


setClass("particle_kv_type", slots=list(s = "numeric"),contains="particle_type")
particle_kv <- function(mu,sigma,H,s)
{
   return(new("particle_type", mu = mu, sigma = sigma, H = H, s = s))
}







