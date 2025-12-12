source("particle.R")
source("update.R")
source("predict.R")

library(collections)


bagel_kv <- function(ptype,prior,H,sigma,y)
{
   population <- dict()
   t <- 1
   
   if(t == 1)
   {
     prior.1 <- prior(t)
     H.1 <- H(t,t)
     # create the first (t = 1) particle
     particle.1 <- ptype(prior.1$mu,prior.1$sigma,H.1)
     # initialise first particle
     particle.1 <- update(particle.1,y)
     # add it to the population
     population$set(1,particle.1)
   } 


  t <- t + 1
  particles$set(t,predict(population$get(1),prior(t),H(t,t)))

  return(particles)




}