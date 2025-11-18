

generate_piecewise_linear_data <- function(n, tau, beta1,beta2,gamma1,gamma2, sigma2, continuous=T) {
  x1 <- beta1+beta2*seq(1,tau,1)
  if(continuous == T){
    gamma1 = -tau*gamma2
    x2 = (beta1+gamma1)+(beta2+gamma2)*seq((tau+1),n,1)
  }else{
    x2 = (beta1+gamma1)+(beta2+gamma2)*seq((tau+1),n,1)}
  x = c(x1,x2)+rnorm(n,0,sqrt(sigma2))
  return(x)
}