
source("particle.R")
source("bagel.R")

source("example-2.R")

##data
set.seed(0)
y <- rnorm(100)


# initialise bagel
res <- bagel(particle,prior,H)

res <- update(res,y[1])