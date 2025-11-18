

bagel_mean <- function(x,decay,probability,mu,tau,delta,sigma,threshold,error,K)
{
    return(bagel_mean_original(x,decay,probability, mu,tau,delta,sigma,threshold,error,K))
}

bagel_mean_original <- function(x, decay,prob_change, mu0,tau0,delta0,sigma2,thresh,error='error1',K=100)
{


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



	update_prior <- function(Sig0,mu0,mu_post,Sigma_post)
	{
		Sig_beta_beta <- Sig0[1:1, 1:1]
  		Sig_beta_gamma <- Sig0[1:1, 2:2]
  		Sig_gamma_beta <- Sig0[2:2, 1:1]
  		Sig_gamma_gamma <- Sig0[2:2, 2:2]
  
		mu_beta <- mu0[1:1]
  		mu_gamma <- mu0[2:2]
  
		# Update mean
  		Sig_beta_beta_inv <- solve(Sig_beta_beta)
  		mu_gamma_updated <- mu_gamma + Sig_gamma_beta %*% Sig_beta_beta_inv %*% (mu_post - mu_beta)  # (mu_beta - mu_beta) = 0
  		mu_updated <- c(mu_post, mu_gamma_updated)  # so mu_updated = mu0_dis
  
		# Update covariance
  		Sig_gamma_gamma_updated <- Sig_gamma_gamma + Sig_gamma_beta %*% Sig_beta_beta_inv %*% 
    					(Sigma_post - Sig_beta_beta) %*% Sig_beta_beta_inv %*% Sig_beta_gamma
  		Sig_beta_gamma_updated <- Sig_gamma_beta%*%Sig_beta_beta_inv %*%Sigma_post
  		Sig_gamma_beta_updated <- Sigma_post %*% Sig_beta_beta_inv %*% Sig_beta_gamma
  		Sig_updated <- rbind(cbind(Sigma_post, Sig_beta_gamma_updated),cbind(Sig_gamma_beta_updated, Sig_gamma_gamma_updated))
  		return(list(mu=mu_updated,Sig=Sig_updated))
		}


  TT=length(x)
  be_change=NA
  prune =error_ratio= NA
  candidate=c(1)
  #for data point 1
  p=prob_change
  Sigma0=tau0
  sfstats=update_post_parameter(mu_prior = mu0[1], Sigma_prior = Sigma0, sigma2=sigma2,obs=x[1],w_prev =1-p)
  mu = list(sfstats$mu_post)
  Sigma = list(sfstats$Sigma_post)
  w = sfstats$w
  for (t in 2:TT)
  {
    new_ob = x[t]
    len=length(w)
    #current updating  
    if(t==2){g_previous=c(1,1);g=c(1,1)}else{
      g_previous=g
      q=decay^{(t-1):1}
      q=q/sum(q)
      g=c(1,q)}   
    current=sapply(1:len, function(n) update_post_parameter(mu_prior = mu[[n]], Sigma_prior = Sigma[[n]], sigma2=sigma2,obs=new_ob,w_prev=w[n]*g[n]/g_previous[n]),simplify = TRUE)
    mu = append(mu,current[1,])
    Sigma = append(Sigma,current[2,])
    w = append(w,unlist(current[3,]))
    
    # new leaf
    mu0_n=c(mu0)
    Sig0_n = tau0*matrix(data=c(1,-1,
                             -1,2),ncol=2,byrow=TRUE)
    updated_prior = update_prior(Sig0_n,mu0_n,mu[[1]],Sigma[[1]])
    mu0_n=updated_prior$mu;Sig0_n=updated_prior$Sig
    #    current = update_post_parameter(mu_prior = delta0, tau_prior = tau0, sigma2=sigma2,obs=new_ob,w_prev=w[1]*(-exp(-p*t)+exp(-p*(t-1)))/exp(-p*(t-1)))
    current = update_post_parameter(mu_prior = mu0_n, Sigma_prior = Sig0_n, sigma2=sigma2,obs=new_ob,w_prev=w[1]*p*g[length(g)]/g_previous[1]/(1-p))
    mu = append(mu,list(current$mu_post))
    Sigma = append(Sigma,list(current$Sigma_post))
    w = append(w,unlist(current$w))
    #delete old candidates, and only store current candidates
    mu=mu[-c(1:len)]
    Sigma=Sigma[-c(1:len)]
    w=w[-c(1:len)]
    #marginalize
    w=w/sum(w)
    #    factor = exp(-p)/sum(w)*factor
    be_change[t]=1-w[1]
    if(be_change[t]>thresh){return(list(uncertainty = w, stoppingtime=t, change_trace =be_change,prune=prune,candidate=c(candidate,t),error_ratio=error_ratio,mu=mu,Sigma=Sigma))}
    
    # pruning
    if(length(w)==K){
      if(error=='error1'){
        #error term is the difference between two distributions
        ind =find_min_err_index(w[2:(length(candidate))], mu[2:(length(candidate))], Sigma[2:(length(candidate))])
        error_ratio = c(error_ratio,ind$ratio)
        ind = ind$index
        print(ind)
      }else if(error =='error3'){
        ind = which.min(w[2:(length(candidate))])
        error_ratio = c(error_ratio,NA)}
      ind=ind+1
      prune = c(prune, candidate[ind])  #add to prune list
      candidate = candidate[-ind]  #remove from the searching space
      mu = mu[-ind] #remove corresponding parameters
      Sigma = Sigma[-ind]
      if(error=='error3'){
        w=w[-ind]
      }else{
        w[ind+1]=w[ind]+w[ind+1] #combine both
        w=w[-ind]}
    }else{prune=c(prune,NA)
    error_ratio=c(error_ratio,NA)}
    candidate=c(candidate, t)
  }
  return(list(uncertainty = w, stoppingtime=Inf, change_trace =be_change,prune=prune,candidate=candidate,error_ratio=error_ratio,mu=mu,Sigma=Sigma)) 
}