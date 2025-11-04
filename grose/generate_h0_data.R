

generate_h0_data <- function(n,mean,sd,tau,delta)
{
  x=rnorm(n,mean,sd)
  x[(tau+1):n]=x[(tau+1):n]+delta
  return(x)
}