source("particle.R")
source("example-2.R")


set.seed(0)
y <- rnorm(100)


tau <- 1L

p.1 <- particle(prior,H,tau)

p.1 <- theorem_2(p.1,1L)
p.1 <- theorem_3(p.1,1L,y[1])


p.2 <- theorem_2(p.1,2L)
p.2 <- theorem_3(p.2,2L,y[2])

p.3 <- theorem_2(p.1,3L)
p.3 <- theorem_3(p.3,3L,y[3])



