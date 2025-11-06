


update_post_parameter <- function(mu_prior,Sigma_prior,sigma2,obs,w_prev)
{
	if(length(mu_prior)==1)
	{
		h=matrix(c(1), ncol = 1)
    	}
	else
	{
		h=matrix(c(1, 1), ncol = 1)
    	}
  	Q = as.numeric(sigma2+t(h)%*%Sigma_prior%*%h)
  	e = as.numeric(obs - t(h) %*% mu_prior)
  	A = (Sigma_prior %*% h)* 1/Q
  	Sigma = Sigma_prior - A%*%t(A)*Q #THE INVERSE OF THE LAMBDA
  	mu = mu_prior + A*e
  	w = w_prev * dnorm(obs, t(h)%*%mu_prior, sqrt(t(h)%*%Sigma_prior%*%h+sigma2)) #marginal likelihood
  	return(list(mu_post=mu, Sigma_post=Sigma, w=w))
}

