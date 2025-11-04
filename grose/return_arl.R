

return_arl <- function(n,sigma2,decay, prob_change,mu0,tau0,delta0,thresh,error ='error1',K)
{
  x=rnorm(n,0,sd=sqrt(sigma2))
  return(max(OBCD_online_mean_varknown(x, decay,prob_change, mu0,tau0,delta0,sigma2,thresh=Inf,error,K)$change_trace,na.rm = TRUE))
}