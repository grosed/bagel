
bagel_slope <- function(x, p_dis,p_conti,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error,K,prior)
{

bagel_slope_call <- function(x, p_dis=0.001,p_conti=0.001,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error='error1',K=100,prior='prior1')
{

	update_post_parameter <- function(mu_prior,Sigma_prior,sigma2,obs,w_prev,t)
	{
		if(length(mu_prior)==2){h=matrix(c(1, t), ncol = 1)}else{h=matrix(c(1, t,1,t), ncol = 1)}
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
		Sig_beta_beta <- Sig0[1:2, 1:2]
  		Sig_beta_gamma <- Sig0[1:2, 3:4]
  		Sig_gamma_beta <- Sig0[3:4, 1:2]
  		Sig_gamma_gamma <- Sig0[3:4, 3:4]
  
		mu_beta <- mu0[1:2]
  		mu_gamma <- mu0[3:4]
  
		# Update mean
  		 Sig_beta_beta_inv <- solve(Sig_beta_beta)
  		 mu_gamma_updated <- mu_gamma + Sig_gamma_beta %*% Sig_beta_beta_inv %*% (mu_post - mu_beta)  # (mu_beta - mu_beta) = 0
  		 mu_updated <- c(mu_post, mu_gamma_updated)  # so mu_updated = mu0_dis
  
		# Update covariance
  		Sig_gamma_gamma_updated <- Sig_gamma_gamma + Sig_gamma_beta %*% Sig_beta_beta_inv %*% 
    		(Sigma_post - Sig_beta_beta) %*% Sig_beta_beta_inv %*% Sig_beta_gamma
  		Sig_beta_gamma_updated <- Sig_gamma_beta%*%Sig_beta_beta_inv %*%Sigma_post
  		Sig_gamma_beta_updated <- Sigma_post %*% Sig_beta_beta_inv %*% Sig_beta_gamma
  		Sig_updated <- rbind(
    		cbind(Sigma_post, Sig_beta_gamma_updated),
    		cbind(Sig_gamma_beta_updated, Sig_gamma_gamma_updated)
  		)
  		return(list(mu=mu_updated,Sig=Sig_updated))
	}





  TT=length(x)
  be_change=NA
  prune_dis = error_ratio_dis =prune_conti = error_ratio_conti= NA
  candidate_dis=candidate_conti=c(1)
  #for data point 1
  p=p_dis+p_conti
  Sigma0=diag(c(delta_1,delta_2))
  sfstats=update_post_parameter(mu_prior = mu0, Sigma_prior = Sigma0, sigma2=sigma2,obs=x[1],w_prev =1-p,t=1)
  mu_dis = mu_conti = list(sfstats$mu_post)
  Sigma_dis = Sigma_conti = list(sfstats$Sigma_post)
  w_dis =w_conti = sfstats$w
  prior_ind_dis=prior_ind_conti=c(1)
  for (t in 2:TT) {
    new_ob = x[t]
    len_dis=length(w_dis)
    len_conti=length(w_conti)
    #current updating
    # if(t==2){g_dis=c(1,decay)
    #   norm_pre_dis=decay}else{
    #     g_dis=c(1,prior_ind_dis[-length(prior_ind_dis)]*norm_pre_dis*decay,decay)
    #     norm_pre_dis=norm_pre_dis*decay+decay
    #     g_dis[-1]=g_dis[-1]/norm_pre_dis
    #   }
    if(t==2){
      g_dis_previous=c(1,1)
      g_dis=c(1,1)}else{
        g_dis_previous=g_dis
        q=decay^{(t-1):1}
        q=q/sum(q)
        g_dis=c(1,q)}
    current_dis=sapply(1:len_dis, function(n) update_post_parameter(mu_prior = mu_dis[[n]], Sigma_prior = Sigma_dis[[n]], sigma2=sigma2,obs=new_ob,w_prev=w_dis[n]*g_dis[n]/g_dis_previous[n],t=t))
    mu_dis = append(mu_dis,current_dis[1,])
    Sigma_dis = append(Sigma_dis,current_dis[2,])
    w_dis = append(w_dis,unlist(current_dis[3,]))
    
    # if(t==2){g_conti=c(1,decay)
    # norm_pre_conti=decay}else{
    #   g_conti=c(1,prior_ind_conti[-length(prior_ind_conti)]*norm_pre_conti*decay,decay)
    #   norm_pre_conti=norm_pre_conti*decay+decay
    #   g_conti[-1]=g_conti[-1]/norm_pre_conti
    # }
    # print(g_conti)
    if(t==2){
      g_conti_previous=c(1,1)
      g_conti=c(1,1)}else{
        g_conti_previous=g_conti
        q=decay^{(t-1):1}
        q=q/sum(q)
        g_conti=c(1,q)}
    current_conti=sapply(1:len_conti, function(n) update_post_parameter(mu_prior = mu_conti[[n]], Sigma_prior = Sigma_conti[[n]], sigma2=sigma2,obs=new_ob,w_prev=w_conti[n]*g_conti[n]/g_conti_previous[n],t=(t)))
    mu_conti = append(mu_conti,current_conti[1,])
    Sigma_conti = append(Sigma_conti,current_conti[2,])
    w_conti = append(w_conti,unlist(current_conti[3,]))
    
    # new leaf
    if(prior=='prior1'){
#      mu0_dis=c(mu0,-mu_dis[[1]][2]*t,0)
      mu0_dis=c(mu0,-mu0[2]*t,0)
      Sig0_dis = matrix(data=c(delta_1,0,-delta_1,0,
                               0,delta_2,0,-delta_2,
                               -delta_1,0,2*delta_1+t^2*delta_2,-t*delta_2,
                               0,-delta_2,-t*delta_2,2*delta_2),ncol=4,byrow=TRUE)
      updated_prior = update_prior(Sig0_dis,mu0_dis,mu_dis[[1]],Sigma_dis[[1]])
      mu0_dis=updated_prior$mu;Sig0_dis=updated_prior$Sig
    }else{
      mu0_dis=c(mu0,0,0)
      Sig0_dis = matrix(data=c(delta_1,0,0,0,
                               0,delta_2,0,0,
                               0,0,delta_3,0,
                               0,0,0,delta_4),ncol=4,byrow=TRUE)
    }
    current_dis = update_post_parameter(mu_prior = mu0_dis, Sigma_prior = Sig0_dis, sigma2=sigma2,obs=new_ob,w_prev=w_dis[1]*p_dis*g_dis[length(g_dis)]/g_dis_previous[1],t=t)
    #    current_dis = update_post_parameter(mu_prior = mu0_dis, Sigma_prior = Sig0_dis, sigma2=sigma2,obs=new_ob,w_prev=w_dis[1]*p_dis/(2*t-2),t=t)
    mu_dis = append(mu_dis,list(current_dis$mu_post))
    Sigma_dis = append(Sigma_dis,list(current_dis$Sigma_post))
    w_dis = append(w_dis,unlist(current_dis$w))
    
    mu0_conti=c(mu0,0,0)
    Sig0_conti = matrix(data=c(delta_1,0,0,0,
                               0,delta_2,0,0,
                               0,0,(t-1)^2*delta_4,-(t-1)*delta_4,
                               0,0,-(t-1)*delta_4,delta_4),ncol=4,byrow=TRUE)
    updated_prior = update_prior(Sig0_conti,mu0_conti,mu_conti[[1]],Sigma_conti[[1]])
    mu0_conti=updated_prior$mu;Sig0_conti=updated_prior$Sig
    current_conti = update_post_parameter(mu_prior = mu0_conti, Sigma_prior = Sig0_conti, sigma2=sigma2,obs=new_ob,w_prev=w_conti[1]*p_conti*g_conti[length(g_conti)]/g_conti_previous[1],t=t)
    #    current_conti = update_post_parameter(mu_prior = mu0_conti, Sigma_prior = Sig0_conti, sigma2=sigma2,obs=new_ob,w_prev=w_conti[1]*p_conti/(2*t-2),t=t)
    
    mu_conti = append(mu_conti,list(current_conti$mu_post))
    Sigma_conti = append(Sigma_conti,list(current_conti$Sigma_post))
    w_conti = append(w_conti,unlist(current_conti$w))
    
    #delete old candidates, and only store current candidates
    mu_dis=mu_dis[-c(1:len_dis)]
    Sigma_dis=Sigma_dis[-c(1:len_dis)]
    w_dis=w_dis[-c(1:len_dis)]
    
    mu_conti=mu_conti[-c(1:len_conti)]
    Sigma_conti=Sigma_conti[-c(1:len_conti)]
    w_conti=w_conti[-c(1:len_conti)]
    
    #marginalize
    fac = (sum(w_dis)+sum(w_conti[-1]))
    w_dis=w_dis/fac
    w_conti=w_conti/fac
    be_change[t]=1-w_dis[1]
    
    if(be_change[t]>thresh){
      return(list(uncertainty_dis = w_dis, uncertainty_conti=w_conti,stoppingtime = t, change_trace=be_change,
                  prune_dis=prune_dis, candidate_dis = c(candidate_dis,t),error_ratio_dis=error_ratio_dis,
                  prune_conti=prune_conti, candidate_conti = c(candidate_conti,t),error_ratio_conti=error_ratio_conti,
                  mu_dis = mu_dis,Sigma_dis=Sigma_dis,mu_conti=mu_conti,Sigma_conti=Sigma_conti))}
    
    # pruning
    if((length(w_dis)+length(w_conti))>=K){
      if(error=='KL'){
        #KL distance
        #        ind_conti = KL_min_index(w_conti[2:length(candidate_conti)], lapply(Sigma_conti[2:length(candidate_conti)], function(x) x+diag(c(1e-10,1e-10))), mu_conti[2:length(candidate_conti)])
        ind_conti = KL_min_index(w_conti[2:length(candidate_conti)], Sigma_conti[2:length(candidate_conti)], mu_conti[2:length(candidate_conti)],model='c')
        ind_dis = KL_min_index(w_dis[2:length(candidate_dis)], Sigma_dis[2:length(candidate_dis)], mu_dis[2:length(candidate_dis)],model='d')
      }else if(error=='HD'){
        # Hellinger bound
        ind_conti = Hellinger_min_index(w_conti[2:length(candidate_conti)], Sigma_conti[2:length(candidate_conti)], mu_conti[2:length(candidate_conti)],model='c')
        ind_dis = Hellinger_min_index(w_dis[2:length(candidate_dis)], Sigma_dis[2:length(candidate_dis)], mu_dis[2:length(candidate_dis)],model='d')
      }else if(error == 'Tylor'){
        # taylor assumption
        ind_conti = triangel_min_index(w_conti[2:length(candidate_conti)], Sigma_conti[2:length(candidate_conti)], mu_conti[2:length(candidate_conti)],model='c')
        ind_dis = triangel_min_index(w_dis[2:length(candidate_dis)], Sigma_dis[2:length(candidate_dis)], mu_dis[2:length(candidate_dis)],model='d')
      }else if(error =='test'){
        ind_conti = triangel_min_index(w_conti[2:length(candidate_conti)], Sigma_conti[2:length(candidate_conti)], mu_conti[2:length(candidate_conti)],model='c')
        ind_dis = triangel_min_index(w_dis[2:length(candidate_dis)], Sigma_dis[2:length(candidate_dis)], mu_dis[2:length(candidate_dis)],model='d')
      }else if(error =='benchmark'){
        ind_conti = which.min(w_conti[2:length(candidate_conti)])
        ind_dis = which.min(w_dis[2:length(candidate_dis)])
      }
      #discrete model
      ind_dis=ind_dis+1
      error_ratio_dis = c(error_ratio_dis,w_dis[ind_dis]/w_dis[ind_dis+1])
      prune_dis = c(prune_dis, candidate_dis[ind_dis])  #add to prune list
      candidate_dis = candidate_dis[-ind_dis]  #remove from the searching space
      mu_dis = mu_dis[-ind_dis] #remove corresponding parameters
      Sigma_dis = Sigma_dis[-ind_dis]
      if(error=='benchmark'){
        w_dis = w_dis[-ind_dis]
        prior_ind_dis = prior_ind_dis[-ind_dis]
      }else{
        w_dis[ind_dis+1]=w_dis[ind_dis]+w_dis[ind_dis+1] #combine both
        w_dis = w_dis[-ind_dis]
        #       prior_ind_dis[ind_dis+1]=prior_ind_dis[ind_dis]+prior_ind_dis[ind_dis+1] #combine both
        prior_ind_dis=prior_ind_dis[-ind_dis]}
      #continus model
      ind_conti=ind_conti+1
      error_ratio_conti = c(error_ratio_conti,w_conti[ind_conti]/w_conti[ind_conti+1])
      prune_conti = c(prune_conti, candidate_conti[ind_conti])  #add to prune list
      candidate_conti = candidate_conti[-ind_conti]  #remove from the searching space
      mu_conti = mu_conti[-ind_conti] #remove corresponding parameters
      Sigma_conti = Sigma_conti[-ind_conti]
      if(error=='benchmark'){
        w_conti = w_conti[-ind_conti]
        prior_ind_conti = prior_ind_conti[-ind_conti]}else{
          w_conti[ind_conti+1]=w_conti[ind_conti]+w_conti[ind_conti+1] #combine both
          w_conti = w_conti[-ind_conti]
          prior_ind_conti=prior_ind_conti[-ind_conti]}
    }else{
      prune_dis=c(prune_dis,NA)
      prune_conti=c(prune_conti,NA)
      error_ratio_conti=c(error_ratio_conti,NA)
      error_ratio_dis=c(error_ratio_dis,NA)}
    candidate_dis=c(candidate_dis, t)
    candidate_conti=c(candidate_conti, t)
    prior_ind_dis=c(prior_ind_dis,1)
    prior_ind_conti=c(prior_ind_conti,1)
  }
  return(list(uncertainty_dis = w_dis, uncertainty_conti=w_conti, stoppingtime=Inf,
              stoppingtime = t, change_trace=be_change,
              prune_dis=prune_dis, candidate_dis = candidate_dis,error_ratio_dis=error_ratio_dis,
              prune_conti=prune_conti, candidate_conti = candidate_conti,error_ratio_conti=error_ratio_conti,
              mu_dis = mu_dis,Sigma_dis=Sigma_dis,mu_conti=mu_conti,Sigma_conti=Sigma_conti))
}

return(bagel_slope_call(x, p_dis,p_conti,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error,K,prior))

}