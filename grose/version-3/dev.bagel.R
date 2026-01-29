

source("particle.R")
source("example-2.R")


set.seed(0)
Y <- rnorm(100)


particles <- list()

tau <- 0L
t <- 1L


p <- particle(prior,H,tau)
p <- theorem_2(p,t)
p <- theorem_3(p,t,Y[t])
particles <- append(particles,p)

for(t in 2:length(Y))
# for(t in 2:2)
{
   t <- as.integer(t)
   particles <- append(particles,theorem_2(particles[[1]],t)) # add the new particle
   for(p in 1:length(particles))
   {
      particles[[p]] <- theorem_3(particles[[p]],t,Y[t])
   }
}