
library(curry)

source("particle.R")
source("example-2.R")


set.seed(0)
Y <- rnorm(100)


particles <- list()

tau <- 0L
t <- 1L

#p <- particle(prior,H,tau,0.9,1.0)
p <- particle_kv(prior,H,tau,0.9,1.0,1.0)
p <- theorem_2(p,t)
p <- theorem_1(p,p,t)
p <- theorem_4(p,t,Y[t])
p <- theorem_3(p,t,Y[t])

particles <- append(particles,p)

for(t in 2:length(Y))
{
   t <- as.integer(t)
   particles <- append(particles,theorem_2(particles[[1]],t)) # add the new particle
   particles <- Map(theorem_1 %><% list(object.1 = particles[[1]],t = t),particles)
   particles <- Map(theorem_4 %><% list(t = t,y = Y[t]),particles)
   particles <- Map(theorem_3 %><% list(t = t,y = Y[t]),particles)
   sumW <- Reduce("+",Map(function(.) .@weight,particles),0)
   particles <- Map(function(p) return(set_weight(p,p@weight/sumW)),particles) 
}





print("*****")
for(p in particles)
{
  print(p@weight)
}
