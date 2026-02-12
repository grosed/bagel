

library(curry)

source("particle.R")
source("example-2.R")
source("bagel.R")

set.seed(0)
# Y <- rnorm(100)

Y <- c(rnorm(100),rnorm(100) + seq(1,100,1)/20)



p0 <- 0.9
p <- 1.0
sigma <- 1.0



bagel <- bagel_kv(prior,H,p0,p,sigma)

ws <- c()

for(y in Y)
{
  bagel <- update(bagel,y)
  ws <- c(ws,weights(bagel)[1])  
}


plot(ws)



















