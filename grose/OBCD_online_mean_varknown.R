

OBCD_online_mean_varknown <- function(x, decay,prob_change, mu0,tau0,delta0,sigma2,thresh,error='error1',K=100)
{
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