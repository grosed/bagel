source("particle.R")
source("example-2.R")


set.seed(0)
y <- rnorm(100)


tau <- 0L
t <- 1L

p.1 <- particle(prior,H,tau,p = 1.0,p0 = 0.9)
p.1 <- theorem_2(p.1,t)
p.1 <- theorem_1(p.1,p.1,t)
p.1 <- theorem_4(p.1,t,y[t])
p.1 <- theorem_3(p.1,t,y[t])



t <- t + 1L

p.2 <- theorem_2(p.1,t)
p.1 <- theorem_1(p.1,p.1,t)
p.2 <- theorem_1(p.1,p.2,t)
p.1 <- theorem_4(p.1,t,y[t])
p.2 <- theorem_4(p.2,t,y[t])
p.1 <- theorem_3(p.1,t,y[t])
p.2 <- theorem_3(p.2,t,y[t])



t <- t + 1L

p.3 <- theorem_2(p.1,t)
p.1 <- theorem_1(p.1,p.1,t)
p.2 <- theorem_1(p.1,p.2,t)
p.3 <- theorem_1(p.1,p.3,t)
p.1 <- theorem_4(p.1,t,y[t])
p.2 <- theorem_4(p.2,t,y[t])
p.3 <- theorem_4(p.3,t,y[t])
p.1 <- theorem_3(p.1,t,y[t])
p.2 <- theorem_3(p.2,t,y[t])
p.3 <- theorem_3(p.3,t,y[t])


t <- t + 1L

p.4 <- theorem_2(p.1,t)
p.1 <- theorem_1(p.1,p.1,t)
p.2 <- theorem_1(p.1,p.2,t)
p.3 <- theorem_1(p.1,p.3,t)
p.4 <- theorem_1(p.1,p.4,t)
p.1 <- theorem_4(p.1,t,y[t])
p.2 <- theorem_4(p.2,t,y[t])
p.3 <- theorem_4(p.3,t,y[t])
p.4 <- theorem_4(p.4,t,y[t])
p.1 <- theorem_3(p.1,t,y[t])
p.2 <- theorem_3(p.2,t,y[t])
p.4 <- theorem_3(p.4,t,y[t])









