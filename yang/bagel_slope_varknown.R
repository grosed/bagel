library(extraDistr)
library(Rcpp)
library(BH)
library(parallel)

#sourceCpp('~/Bayesian/slope/KL_approx.cpp')
#sourceCpp('~/Bayesian/slope/Hellinger_approx.cpp')
#sourceCpp('~/Bayesian/slope/triangler_approx.cpp')

sourceCpp('triangler_approx.cpp')
sourceCpp('KL_approx.cpp')
sourceCpp('Hellinger_approx.cpp')
update_post_parameter = function(mu_prior,Sigma_prior,sigma2,obs,w_prev,t){
  if(length(mu_prior)==2){h=matrix(c(1, t), ncol = 1)}else{h=matrix(c(1, t,1,t), ncol = 1)}
  Q = as.numeric(sigma2+t(h)%*%Sigma_prior%*%h)
  e = as.numeric(obs - t(h) %*% mu_prior)
  A = (Sigma_prior %*% h)* 1/Q
  Sigma = Sigma_prior - A%*%t(A)*Q #THE INVERSE OF THE LAMBDA
  mu = mu_prior + A*e
  w = w_prev * dnorm(obs, t(h)%*%mu_prior, sqrt(t(h)%*%Sigma_prior%*%h+sigma2)) #marginal likelihood
  return(list(mu_post=mu, Sigma_post=Sigma, w=w))
}
update_prior = function(Sig0,mu0,mu_post,Sigma_post){
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
OBCD_slope_varknown=function(x, p_dis=0.001,p_conti=0.001,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error='error1',K=100,prior='prior1'){
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
#
#################################functions----------------------------------
split_number = function(weight, ratios){
  x1 = ratios/(1+ratios)*weight
  x2 = 1/(1+ratios)*weight
  return(c(x1,x2))
}
recover_posteior =function(pruned_candidates, stored_candidate, w, ratios){
  result <- list()
  result[[1]] = lapply(seq_len(stored_candidate), function(i) i)
  for (i in 1:length(pruned_candidates)) {
    result[[i+1]] = result[[i]]
    for (j in 1:length(result[[i+1]])){
      if(pruned_candidates[i] %in% result[[i+1]][[j]]){
        pruned_index=j
      }
    }
    combined_subset = c(unlist(result[[i+1]][pruned_index]),unlist(result[[i+1]][pruned_index+1]))
    result[[i+1]][[pruned_index+1]] = combined_subset
    result[[i+1]][pruned_index] = NULL
  }
  for (i in 1:length(ratios)) {
    procedure = result[[length(result)-i+1]]
    for (j in 1:length(procedure)){
      if(pruned_candidates[length(pruned_candidates)-i+1] %in% procedure[[j]]){
        class = j
      }
    }
    w = append(w, split_number(w[class],ratios[length(ratios)-i+1]), after =class)
    w = w[-class]
  }
  return(w)}
find_cover = function(w, true_change){
  w=w/sum(w)
  lower=which(cumsum(w)>=0.025)[1]-1
  upper=which(cumsum(w)>=0.975)[1]
  return(list(cover = as.numeric(true_change>=lower & true_change <= upper), len = upper-lower))
}
run_benchmark = function(res,tau){
  ratio_int=sum(res$uncertainty_dis[-1])/sum(res$uncertainty_conti[-1])
  ratio_max=max(res$uncertainty_dis[-1])/max(res$uncertainty_conti[-1])
  if(max(res$uncertainty_dis[-1])>max(res$uncertainty_conti[-1])){
    uncertainty = res$uncertainty_dis
    prune = res$prune_dis
    candidate = res$candidate_dis
    model = 'd'
  }else{
    uncertainty = res$uncertainty_conti
    prune = res$prune_conti
    candidate = res$candidate_conti
    model = 'c'}
  if(is.infinite(res$stoppingtime)){
    post_prob_benchmark = rep(0,min(res$stoppingtime,2000))
    post_prob_benchmark[candidate] = uncertainty
  }else{
    post_prob_benchmark = rep(0,min(res$stoppingtime,2000))
    post_prob_benchmark[candidate[-1]] = uncertainty[-1]}
  post_prob_benchmark=post_prob_benchmark/sum(post_prob_benchmark)
  map_benchmark=which.max(post_prob_benchmark)
  coverage_benchmark=find_cover(post_prob_benchmark,tau)$cover
  coverage_length=find_cover(post_prob_benchmark,tau)$len
  stop_benchmark = res$stoppingtime
  return(list(map=map_benchmark,coverage=coverage_benchmark,prob=post_prob_benchmark,stop=stop_benchmark,
              length=coverage_length,model=model,ratio_int=ratio_int,ratio_max=ratio_max))
}
exact_post = function(res, tau){
  ratio_int=sum(res$uncertainty_dis[-1])/sum(res$uncertainty_conti[-1])
  ratio_max=max(res$uncertainty_dis[-1])/max(res$uncertainty_conti[-1])
  if(max(res$uncertainty_dis[-1])>max(res$uncertainty_conti[-1])){
    if(is.infinite(res$stoppingtime)){post_prob_noprune =res$uncertainty_dis}else{
      post_prob_noprune = c(0,res$uncertainty_dis[-1])}
    post_prob_noprune = post_prob_noprune/sum(post_prob_noprune)
    map_noprune = which.max(post_prob_noprune)
    coverage_noprune=find_cover(post_prob_noprune,tau)$cover
    length_noprune=find_cover(post_prob_noprune,tau)$len
    stop_noprune = min(2000,res$stoppingtime)
    return(list(model='d',map=map_noprune,coverage=coverage_noprune,prob=post_prob_noprune,stop=stop_noprune,
                ratio_int=ratio_int,ratio_max=ratio_max,length=length_noprune))
  }else{
    if(is.infinite(res$stoppingtime)){post_prob_noprune = res$uncertainty_conti}else{
      post_prob_noprune = c(0,res$uncertainty_conti[-1])}
    post_prob_noprune = post_prob_noprune/sum(post_prob_noprune)
    map_noprune = which.max(post_prob_noprune)
    coverage_noprune=find_cover(post_prob_noprune,tau)$cover
    length_noprune=find_cover(post_prob_noprune,tau)$len
    stop_noprune = min(2000,res$stoppingtime)
    return(list(model='c',map=map_noprune,coverage=coverage_noprune,prob=post_prob_noprune,stop=stop_noprune,
                ratio_int=ratio_int,ratio_max=ratio_max,length=length_noprune))
  }}
run_ours = function(res,tau){
  ratio_int=sum(res$uncertainty_dis[-1])/sum(res$uncertainty_conti[-1])
  ratio_max=max(res$uncertainty_dis[-1])/max(res$uncertainty_conti[-1])
  if(max(res$uncertainty_dis[-1])>max(res$uncertainty_conti[-1])){
    model='d'
    uncertainty = res$uncertainty_dis
    prune = res$prune_dis
    candidate = res$candidate_dis
    error_ratio = res$error_ratio_dis
  }else{
    model='c'
    uncertainty = res$uncertainty_conti
    prune = res$prune_conti
    candidate = res$candidate_conti
    error_ratio = res$error_ratio_conti}
  w=6
  w[1]=uncertainty[1]
  for (i in 2:length(uncertainty)) {
    index = which(candidate[i-1]< prune & prune<candidate[i])
    if(length(index)==0){
      w=c(w, uncertainty[i])
    }else{
      pruned_candidates=prune[index]-candidate[i-1];
      stored_candidate=candidate[i]-1-candidate[i-1];
      weight=uncertainty[i];
      ratios=error_ratio[index]
      w_1=recover_posteior(pruned_candidates,stored_candidate,weight,ratios)
      w=c(w,as.numeric(w_1))
    }
  }
  if(is.infinite(res$stoppingtime)){
    post_prob_prune = w
  }else{
    post_prob_prune = c(0,w[-1])}
  post_prob_prune=post_prob_prune/sum(post_prob_prune)
  map_prune=which.max(post_prob_prune)
  coverage_prune=find_cover(post_prob_prune,tau)$cover
  coverage_length=find_cover(post_prob_prune,tau)$len
  stop_prune = res$stoppingtime
  return(list(map=map_prune,coverage=coverage_prune,prob=post_prob_prune,stop=stop_prune,
              ratio_int=ratio_int,ratio_max=ratio_max,length=coverage_length,model=model))
}
find_cover = function(w, true_change){
  w=w/sum(w)
  cumw = cumsum(w)
  ord = order(w, decreasing=TRUE)
  cumw_sorted = cumsum(w[ord])
  cut_idx = which(cumw_sorted >= 0.95)[1]
  selected = ord[1:cut_idx]  
  return(list(cover = as.numeric(true_change %in% selected), len = length(selected)))
  
  # lower=which(cumsum(w)>=0.025)[1]-1
  # upper=which(cumsum(w)>=0.975)[1]
  #  return(list(cover = as.numeric(true_change>=lower & true_change <= upper), len = upper-lower))
}

#################################simulation----------------------------------
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


# set.seed(1)
# y=generate_piecewise_linear_data(2000,1000,0,0,-2.2,0.002,1,continuous=T)
# plot(y)
# 
# res=OBCD_slope_varknown(y,p_dis=0.1,p_conti=0.1,decay=1,mu0=c(0,0),
#                         delta_1=0.5^2,delta_2=0.5^2,delta_3=0.5,delta_4=1,
#                         sigma2=1,thresh=0.89,error ='KL',K=50000,prior='prior1')
# sum(res$uncertainty_dis[-1])
# sum(res$uncertainty_conti[-1])
# plot(res$candidate_dis[-1],res$uncertainty_dis[-1])
# plot(res$candidate_conti[-1],res$uncertainty_conti[-1])
# plot(res$change_trace)
# plot(y[1:res$stoppingtime])
# abline(v=1000)
# res$uncertainty_dis[1173]
# res$uncertainty_conti[1157]
# which.max(res$uncertainty_dis[-1])
# which.max(res$uncertainty_conti[-1])
run_simulation_add <- function(y,tau,p_dis,p_conti,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error =error,K=K,prior) {
  tem=OBCD_slope_varknown(y,p_dis,p_conti,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error =error,K=K,prior)
  if(error=='benchmark'){
    nam=run_benchmark(tem,tau)
    stop = nam$stop;map =nam$map;coverage =nam$coverage;length=nam$length;model=nam$model;ratio_int=nam$ratio_int;ratio_max=nam$ratio_max
  }else if(error =='exact'){
    nam=exact_post(tem,tau)
    stop = nam$stop;map =nam$map;coverage =nam$coverage;length=nam$length;model=nam$model;ratio_int=nam$ratio_int;ratio_max=nam$ratio_max
  }else{
    nam=run_ours(tem,tau)
    stop = nam$stop;map =nam$map;coverage =nam$coverage;length=nam$length;model=nam$model;ratio_int=nam$ratio_int;ratio_max=nam$ratio_max
  }
  return(list(stop = stop,map=map,coverage=coverage,length=length,model=model,ratio_int=ratio_int,ratio_max=ratio_max))
}


return_arl=function(n,p_dis,p_conti,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh=Inf,error =error,K=K,prior){
  x=rnorm(n,0,sd=sqrt(sigma2))
  return(max(OBCD_slope_varknown(x,p_dis,p_conti,decay,mu0,delta_1,delta_2,delta_3,delta_4,sigma2,thresh,error =error,K=K,prior)$change_trace,na.rm = TRUE))
}



simNo=as.numeric(commandArgs(trailingOnly = TRUE))

error_values <- c('KL', 'benchmark')
K_values <- c(50,100,200)
prob_values <- c(0.1)
delta_values <- c(0.125^2,0.25^2,0.5^2,1)
results <- expand.grid(K = K_values, error = error_values,prob_value=prob_values,delta=delta_values)
prepend_rows <- data.frame(K = c(130000, 130000, 130000, 130000),
                           error = c("KL", "KL", "KL", "KL"),
                           prob_value=c(0.1,0.1,0.1,0.1),
                           delta =  c(0.125^2,0.25^2,0.5^2,1),
                           stringsAsFactors = FALSE
)
results = rbind(prepend_rows,results)

pram <- results[simNo, ]
#find the threshold - continuous change
# set.seed(1)
# cat("Finding threshold", simNo, "/", nrow(results), "\n")
# at_arl <- mclapply(1:500,function(i){
#   set.seed(i)
#   return_arl(2000, p_dis=prob_value, p_conti=prob_value,decay=1,mu0=c(0,0),
#              delta_1=pram$delta,delta_2=pram$delta,delta_3=pram$delta,delta_4=pram$delta,
#              sigma2=1,thresh=Inf,error =pram$error,K=pram$K,prior='prior1')},
#   mc.set.seed = TRUE,mc.cores = detectCores() %/% 2-1)
# c=quantile(unlist(at_arl),0.95)
filename <- paste0("slope_known_thresh", simNo, ".Rdata")
#save(c, file = filename)
#cat("Testing", simNo, "/", nrow(results), "\n",c)
load(filename)
c=as.numeric(c)
set.seed(1)
res <- mclapply(1:500, function(i) {
  set.seed(i)
  y=generate_piecewise_linear_data(2000,1000,0,0,-1.75,0.002,1,continuous=F)
  run_simulation_add(y,tau=1000,p_dis=pram$prob_value, p_conti=0.1,decay=1,mu0=c(0,0),
                     delta_1=pram$delta,delta_2=pram$delta,delta_3=pram$delta,delta_4=pram$delta,
                     sigma2=1,thresh=c,error =pram$error,K=pram$K,prior='prior1')},
  mc.set.seed = TRUE,mc.cores = detectCores() %/% 2-1)
filename <- paste0("slope_known_fix_tau_dis_new_175", simNo, ".Rdata")
save(res, file = filename)
cat("Testing", simNo, "/", nrow(results), "\n",c)
set.seed(1)
res <- mclapply(1:500, function(i) {
  set.seed(i)
  y=generate_piecewise_linear_data(2000,1000,0,0,-4,0.002,1,continuous=T)
  run_simulation_add(y,tau=1000,p_dis=pram$prob_value, p_conti=0.1,decay=1,mu0=c(0,0),
                     delta_1=pram$delta,delta_2=pram$delta,delta_3=pram$delta,delta_4=pram$delta,
                     sigma2=1,thresh=c,error =pram$error,K=pram$K,prior='prior1')},
  mc.set.seed = TRUE,mc.cores = detectCores() %/% 2-1)
filename <- paste0("slope_known_fix_tau_conti_new", simNo, ".Rdata")
save(res, file = filename)

# 
# 
# 
# ###check the results------------------------------
# library(dplyr)

# outputs_table=function(res){
#   df <- bind_rows(res)
#   coverage_rate <- mean(df$coverage == 1)
#   model_rate<-mean(df$model=='c')
#   failure_rate <- sum(is.infinite(df$stop))
#   false_alarm <- sum(df$stop<999)
#   index_out <- which(df$stop<999|is.infinite(df$stop))
#   drop_out <- length(index_out)/500
#   mean_delay <- mean(df$stop[-index_out]-1000)
#   sd_delay <- sd(df$stop[-index_out]-1000)
#   mean_map <- mean(df$map[-index_out])
#   sd_map <- sd(df$map[-index_out])
#   mean_length <- mean(df$length[-index_out])
#   sd_length <- sd(df$length[-index_out])
#   model_rate_out<-mean(df$model[-index_out]=='c')
#   coverage_rate_out <- mean(df$coverage[-index_out] == 1)
#   coverage_rate_in <- mean(df$coverage[index_out] == 1)
#   ratio_int <- mean(df$ratio_int[-index_out])
#   ratio_max <- mean(df$ratio_max[-index_out])
#   result <- data.frame(coverage_rate = coverage_rate,coverage_rate_out=coverage_rate_out,ratio_int=ratio_int,ratio_max=ratio_max,model_rate=model_rate,model_rate_out=model_rate_out,drop_out=drop_out,failure_rate = failure_rate,
#                        false_alarm = false_alarm,mean_delay = mean_delay,sd_delay = sd_delay,
#                        mean_map = mean_map,sd_map = sd_map,mean_length = mean_length,
#                        sd_length = sd_length)
#   return(result)
# }
# file_list <- list.files('~/Bayesian/slope/results', pattern = "tau_dis_new_18[0-9]+\\.Rdata$", full.names = TRUE)
# tau_numbers <- as.numeric(gsub(".*dis_new_18([0-9]+)\\.Rdata$", "\\1", file_list))
# 
# file_list <- list.files('~/Bayesian/slope/results', pattern = "tau_dis_old[0-9]+\\.Rdata$", full.names = TRUE)
# tau_numbers <- as.numeric(gsub(".*dis_old([0-9]+)\\.Rdata$", "\\1", file_list))
# 
# 
# file_list <- list.files('~/Bayesian/slope/results', pattern = "tau_dis_new_15[0-9]+\\.Rdata$", full.names = TRUE)
# tau_numbers <- as.numeric(gsub(".*dis_new_15([0-9]+)\\.Rdata$", "\\1", file_list))
file_list <- list.files('~/Bayesian/slope/results', pattern = "tau_dis_new_22[0-9]+\\.Rdata$", full.names = TRUE)
tau_numbers <- as.numeric(gsub(".*dis_new_22([0-9]+)\\.Rdata$", "\\1", file_list))

file_list <- list.files('~/Bayesian/slope/results', pattern = "tau_conti_new[0-9]+\\.Rdata$", full.names = TRUE)
tau_numbers <- as.numeric(gsub(".*conti_new([0-9]+)\\.Rdata$", "\\1", file_list))
#filename <- paste0("slope_known_fix_tau_conti_new", simNo, ".Rdata")

file_list <- list.files('~/Bayesian/slope/results', pattern = "tau_dis_new_175[0-9]+\\.Rdata$", full.names = TRUE)
tau_numbers <- as.numeric(gsub(".*dis_new_175([0-9]+)\\.Rdata$", "\\1", file_list))

sorted_file_list <- file_list[order(tau_numbers)]
all_results <- list()
for (file in sorted_file_list[1:28]) {
  load(file)
  all_results <- append(all_results, list(outputs_table(res)))
}
all_results=bind_rows(all_results)
final_results=bind_cols(results,all_results)
final_results
# # #
write.csv(final_results,file = '~/Bayesian/slope/results/dis_results_new_conti.csv',row.names = FALSE)

# 
# 
# # speed ------------------------
# library(microbenchmark)
# library(ggplot2)
# 
# set.seed(1)
# # 
# k_vals <- c(20,50,100, 150, 200)
# benchmark_results <- list()
# res_exact <- microbenchmark(
#   'exact' = {
#     y=c(rnorm(1000,0,1))
#     OBCD_slope_varknown(y,p_dis=0.2,p_conti=0.4,decay=1,mu0=c(0,0),
#                                                  delta_1=0.001^2,delta_2=0.001^2,delta_3=0.001^2,delta_4=0.2^2,
#                                                  sigma2=1,thresh=Inf,error ='exact',K=10000000,prior='prior1')})
# res_exact$K <- "exact"
# benchmark_results[["exact"]] <- res_exact
# 
# for (k in k_vals) {
#   cat("Running benchmark for K =", k, "\n")
#   res_K <- microbenchmark(
#       'KL' ={    
#         y=c(rnorm(1000,0,1))
#         OBCD_slope_varknown(y,p_dis=0.2,p_conti=0.4,decay=1,mu0=c(0,0),
#                                                  delta_1=0.001^2,delta_2=0.001^2,delta_3=0.001^2,delta_4=0.2^2,
#                                                  sigma2=1,thresh=Inf,error ='KL',K=k,prior='prior1')},
#     'HD' ={
#       y=c(rnorm(1000,0,1))
#       OBCD_slope_varknown(y,p_dis=0.2,p_conti=0.4,decay=1,mu0=c(0,0),
#                           delta_1=0.001^2,delta_2=0.001^2,delta_3=0.001^2,delta_4=0.2^2,
#                           sigma2=1,thresh=Inf,error ='HD',K=k,prior='prior1')},
#     'Tylor' ={
#       y=c(rnorm(1000,0,1))
#       OBCD_slope_varknown(y,p_dis=0.2,p_conti=0.4,decay=1,mu0=c(0,0),
#                           delta_1=0.001^2,delta_2=0.001^2,delta_3=0.001^2,delta_4=0.2^2,
#                           sigma2=1,thresh=Inf,error ='Tylor',K=k,prior='prior1')},
#     'benchmark' ={
#     y=c(rnorm(1000,0,1))
#       OBCD_slope_varknown(y,p_dis=0.2,p_conti=0.4,decay=1,mu0=c(0,0),
#                           delta_1=0.001^2,delta_2=0.001^2,delta_3=0.001^2,delta_4=0.2^2,
#                           sigma2=1,thresh=Inf,error ='benchmark',K=k,prior='prior1')},
#   times=200)
#   res_K$K <- k  # tag each result with K value
#   benchmark_results[[as.character(k)]] <- res_K
# }
# 
# library(dplyr)
# library(purrr)
# summary_df <- map_dfr(names(benchmark_results), function(name) {
#   df <- as.data.frame(benchmark_results[[name]])
#   df$K <- name
#   df
# }) %>%
#   group_by(K, expr) %>%
#   summarise(mean_ms = mean(time) / 1e6/1000, .groups = "drop")  # convert ns to ms per step
# 
# 
# summary_df <- summary_df %>%
#   mutate(K_numeric = suppressWarnings(as.numeric(as.character(K))))
# 
# custom_colors <- c(
#   "exact" = "#D6CFC4",
#   "KL" = "#466CA6",
#   "HD" = "#D1AE45",
#   "Tylor" = "#872410",
#   "benchmark" = "#040203"
# )
# exact_df <- summary_df %>% filter(expr == "exact")
# non_exact_df <- summary_df %>% filter(expr != "exact")
# ggplot() +
#   geom_line(data = non_exact_df, aes(x = K_numeric, y = mean_ms, color = expr), size = 2) +
#   geom_point(data = non_exact_df, aes(x = K_numeric, y = mean_ms, color = expr), size = 2) +
#   geom_hline(data = exact_df, aes(yintercept = mean_ms), color = custom_colors["exact"], size = 2) +
#   labs(x = "M",y = "Mean Runtime (ms)",color = "Method") +
#   scale_color_manual(
#     values = custom_colors)+theme_classic(base_size = 15)+theme(legend.position = 'none')
# ggsave("runtime_slope_knownvar.pdf", width = 6, height = 4, units = "in")
