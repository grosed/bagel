
predict <- function(particle,prior,H)
{
  d1 <- nrow(particle@mu)
  n <- nrow(prior$mu)
  sigma.beta.beta <- prior$sigma[1:d1,1:d1]
  sigma.beta.gamma <- prior$sigma[1:d1,(d1+1):n]
  sigma.gamma.beta <- prior$sigma[(d1+1):n,1:d1]
  sigma.gamma.gamma <- prior$sigma[(d1+1):n,(d1+1):n]
  mu.beta <- prior$mu[1:d1,1]
  mu.gamma <- prior$mu[(d1+1):n,1]
  inv.sigma.beta.beta <- solve(sigma.beta.beta)

  
  return(
  particle(
  rbind(
  particle@mu,
  mu.gamma + sigma.gamma.beta %*% inv.sigma.beta.beta %*% (particle@mu - mu.beta)),
  rbind(
  cbind(
  particle@sigma,
  sigma.gamma.beta %*% inv.sigma.beta.beta %*% particle@sigma),
  cbind(
  particle@sigma %*% inv.sigma.beta.beta %*% sigma.gamma.beta,
  sigma.gamma.gamma + sigma.gamma.beta %*% inv.sigma.beta.beta %*% (particle@sigma - sigma.beta.beta ) %*% inv.sigma.beta.beta %*% sigma.beta.gamma)),
  H))
}