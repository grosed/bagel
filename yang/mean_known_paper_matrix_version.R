library(Rcpp)
library(parallel)
sourceCpp('find_min_err_index_list_matrix.cpp')
#sourceCpp('~/Bayesian/mean_change_var_known/find_min_err_index_list_matrix.cpp')
#sourceCpp('~/Bayesian/mean_change_var_known/approx_pnorm.cpp')
update_post_parameter = function(mu_prior,Sigma_prior,sigma2,obs,w_prev){
  if(length(mu_prior)==1){
    h=matrix(c(1), ncol = 1)}else{h=matrix(c(1, 1), ncol = 1)}
  Q = as.numeric(sigma2+t(h)%*%Sigma_prior%*%h)
  e = as.numeric(obs - t(h) %*% mu_prior)
  A = (Sigma_prior %*% h)* 1/Q
  Sigma = Sigma_prior - A%*%t(A)*Q #THE INVERSE OF THE LAMBDA
  mu = mu_prior + A*e
  w = w_prev * dnorm(obs, t(h)%*%mu_prior, sqrt(t(h)%*%Sigma_prior%*%h+sigma2)) #marginal likelihood
  return(list(mu_post=mu, Sigma_post=Sigma, w=w))
}
update_prior = function(Sig0,mu0,mu_post,Sigma_post){
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
  Sig_updated <- rbind(
    cbind(Sigma_post, Sig_beta_gamma_updated),
    cbind(Sig_gamma_beta_updated, Sig_gamma_gamma_updated)
  )
  return(list(mu=mu_updated,Sig=Sig_updated))
}

OBCD_online_mean_varknown=function(x, decay,prob_change, mu0,tau0,delta0,sigma2,thresh,error='error1',K=100){
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
  for (t in 2:TT) {
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

split_number = function(weight, ratios){
  x1 = ratios/(1+ratios)*weight
  x2 = 1/(1+ratios)*weight
  return(c(x1,x2))
}
recover_posteior =function(pruned_candidates, stored_candidate, w, ratios){
  result <- list()
  result[[1]] = lapply(seq_len(stored_candidate), function(ii) ii)
  for (ii in 1:length(pruned_candidates)) {
    result[[ii+1]] = result[[ii]]
    for (j in 1:length(result[[ii+1]])){
      if(pruned_candidates[ii] %in% result[[ii+1]][[j]]){
        pruned_index=j
      }
    }
    combined_subset = c(unlist(result[[ii+1]][pruned_index]),unlist(result[[ii+1]][pruned_index+1]))
    result[[ii+1]][[pruned_index+1]] = combined_subset
    result[[ii+1]][pruned_index] = NULL
  }
  for (iii in 1:length(ratios)) {
    procedure = result[[length(result)-iii+1]]
    for (jj in 1:length(procedure)){
      if(pruned_candidates[length(pruned_candidates)-iii+1] %in% procedure[[jj]]){
        class_index = jj
      }
    }
    w = append(w, split_number(w[class_index],ratios[length(ratios)-iii+1]), after =class_index)
    w = w[-class_index]
  }
  return(w)}
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
run_ours = function(res,tau){
  w=6
  w[1]=res$uncertainty[1]
  for (i in 2:length(res$uncertainty)) {
    index = which(res$candidate[i-1]< res$prune & res$prune<res$candidate[i])
    if(length(index)==0){
      w=c(w, res$uncertainty[i])
    }else{
      pruned_candidates=res$prune[index]-res$candidate[i-1]
      stored_candidate=res$candidate[i]-1-res$candidate[i-1]
      weight=res$uncertainty[i]
      ratios=res$error_ratio[index]
      w_1=recover_posteior(pruned_candidates,stored_candidate,weight,ratios)
      w=c(w,as.numeric(w_1))
    }
  }
  if(is.infinite(res$stoppingtime)){
    post_prob_prune = w
  }else{
    post_prob_prune = c(0,w[-1]/sum(w[-1]))}
  map_prune=which.max(post_prob_prune)
  coverage_prune=find_cover(post_prob_prune,tau)$cover
  coverage_length=find_cover(post_prob_prune,tau)$len
  stop_prune = res$stoppingtime
  mu_prune = res$mu[which.max(res$uncertainty)]
  tau_prune = res$tau[which.max(res$uncertainty)]
  return(list(map=map_prune,coverage=coverage_prune,prob=post_prob_prune,stop=stop_prune,length=coverage_length,
              mu=mu_prune,tau=tau_prune))
}
run_benchmark = function(res,tau){
  if(is.infinite(res$stoppingtime)){
    post_prob_benchmark = rep(0,min(res$stoppingtime,2000))
    post_prob_benchmark[res$candidate] = res$uncertainty
  }else{
    post_prob_benchmark = rep(0,min(res$stoppingtime,2000))
    post_prob_benchmark[res$candidate[-1]] = res$uncertainty[-1]
    post_prob_benchmark=post_prob_benchmark/sum(post_prob_benchmark)
  }
  map_benchmark = which.max(post_prob_benchmark)
  coverage_benchmark=find_cover(post_prob_benchmark,tau)$cover
  coverage_length=find_cover(post_prob_benchmark,tau)$len
  stop_benchmark = res$stoppingtime
  mu_benchmark = res$mu[which.max(res$uncertainty)]
  tau_benchmark = res$tau[which.max(res$uncertainty)]
  return(list(map=map_benchmark,coverage=coverage_benchmark,prob=post_prob_benchmark,
              stop=stop_benchmark,length=coverage_length,mu=mu_benchmark,tau=tau_benchmark))
}



generate_h0_data<- function(n,mean,sd,tau,delta){
  x=rnorm(n,mean,sd)
  x[(tau+1):n]=x[(tau+1):n]+delta
  return(x)
}



run_simulation_add <- function(x,tau,decay,prob_change,mu0,tau0,delta0,sigma2,thresh,error,K) {
  tem=OBCD_online_mean_varknown(x, decay,prob_change, mu0,tau0,delta0,sigma2,thresh,error,K)
  if(error=='error3'){
    nam=run_benchmark(tem,tau)
    stop = nam$stop;map =nam$map;coverage =nam$coverage;length=nam$length
  }else if(error =='exact'){
    if(is.infinite(tem$stoppingtime)){
      post_prob_noprune=tem$uncertainty
    }else{
      post_prob_noprune=c(0,tem$uncertainty[-1])
    }
    post_prob_noprune=post_prob_noprune/sum(post_prob_noprune)
    map_noprune = tem$candidate[which.max(post_prob_noprune)]
    coverage_noprune=find_cover(post_prob_noprune,tau)$cover
    length_noprune=find_cover(post_prob_noprune,tau)$len
    stop_noprune = tem$stoppingtime
    stop = stop_noprune;map = map_noprune;coverage = coverage_noprune;length=length_noprune
  }else{
    nam=run_ours(tem,tau)
    stop = nam$stop;map =nam$map;coverage =nam$coverage;length=nam$length
  }
  return(list(stop = stop,map=map,coverage=coverage,length=length))
} 

return_arl=function(n,sigma2,decay, prob_change,mu0,tau0,delta0,thresh,error ='error1',K){
  x=rnorm(n,0,sd=sqrt(sigma2))
  return(max(OBCD_online_mean_varknown(x, decay,prob_change, mu0,tau0,delta0,sigma2,thresh=Inf,error,K)$change_trace,na.rm = TRUE))
}



simNo=as.numeric(commandArgs(trailingOnly = TRUE))

prior_profiles <- data.frame(
  tau0 = c(0.125^2,0.25^2,0.5^2,1)
)
other_grid <- expand.grid(
  error = c("error1", "error3"),
  decay = c(1),
  K = c(50, 100, 200),
  stringsAsFactors = FALSE
)
grid_full <- merge(other_grid, prior_profiles, by = NULL)

prepend_rows <- data.frame(error = rep("exact",4),
                           decay = c(1, 1,1,1),K = rep(130000,4),
                           tau0 = c(0.125^2,0.25^2,0.5^2,1),
                           stringsAsFactors = FALSE
)
results = rbind(grid_full,prepend_rows)

prob_value=0.1
pram <- results[simNo, ]
#find the threshold
set.seed(1)
cat("Finding threshold", simNo, "/", nrow(results), "\n")

at_arl <- mclapply(1:500,function(i){
  set.seed(i)
  return_arl(n=2000,sigma2=1,decay=pram$decay, prob_change=prob_value,mu0=c(0,0),tau0=pram$tau0,delta0=0,thresh=Inf,error =pram$error,K=pram$K)
},mc.set.seed = TRUE,mc.cores = detectCores() %/% 2-1)
c=quantile(unlist(at_arl),0.95)
cat("Testing", simNo, "/", nrow(results), "\n")
#find the result
res <- mclapply(1:500, function(i) {
  set.seed(i)
  y=generate_h0_data(n=2000,mean=0,sd=1,tau=1000,delta=0.5)
  run_simulation_add(y,tau=1000,decay=pram$decay,prob_change=prob_value,mu0=c(0,0),delta0=0,
                     tau0=pram$tau0,sigma2=1,thresh=c,error =pram$error,K=pram$K)},
  mc.set.seed = TRUE,mc.cores = detectCores() %/% 2-1)
filename <- paste0("mean_known_fix_tau_np", simNo, ".Rdata")  
save(res, file = filename)  

# x=c(rnorm(500,0,1),rnorm(500,0.5,1))
# res=OBCD_online_mean_varknown(x, decay=1,prob_change=0.00001, mu0=c(0,0),tau0=1,delta0=0,sigma2=1,thresh=Inf,error='error3',K=100)
# plot(res$candidate,res$uncertainty)
# abline(v=500)
# 
# library(dplyr)
# outputs_table=function(res){
#   df <- bind_rows(res)
#   coverage_rate <- mean(df$coverage == 1)
#   failure_rate <- sum(is.infinite(df$stop))
#   false_alarm <- sum(df$stop<1000)
#   index_out <- which(df$stop<1000|is.infinite(df$stop))
#   drop_out <- length(index_out)/500
#   mean_delay <- mean(df$stop[-index_out]-1000)
#   sd_delay <- sd(df$stop[-index_out]-1000)
#   mean_map <- mean(df$map[-index_out])
#   sd_map <- sd(df$map[-index_out])
#   mean_length <- mean(df$length[-index_out])
#   sd_length <- sd(df$length[-index_out])
#   coverage_rate_out <- mean(df$coverage[-index_out] == 1)
#   coverage_rate_in <- mean(df$coverage[index_out] == 1)
#   result <- data.frame(coverage_rate = coverage_rate,coverage_rate_out=coverage_rate_out,drop_out=drop_out,failure_rate = failure_rate,
#                        false_alarm = false_alarm,mean_delay = mean_delay,sd_delay = sd_delay,
#                        mean_map = mean_map,sd_map = sd_map,mean_length = mean_length,
#                        sd_length = sd_length)
#   return(result)
# }
# file_list <- list.files('~/Bayesian/mean_change_var_known/results', pattern = "tau_np[0-9]+\\.Rdata$", full.names = TRUE)
# tau_numbers <- as.numeric(gsub(".*tau_np([0-9]+).*", "\\1", file_list))
# sorted_file_list <- file_list[order(tau_numbers)]
# all_results <- list()
# for (file in sorted_file_list) {
#   load(file)
#   all_results <- append(all_results, list(outputs_table(res)))
# }
# all_results=bind_rows(all_results)
# final_results=bind_cols(results,all_results)
# write.csv(final_results,file = '~/Bayesian/mean_change_var_known/results/known_results_np.csv',row.names = FALSE)
# 

# load('/home/yangz40/Bayesian/mean_change_var_known/results/mean_known_fix_tau_np1.Rdata')
# boxplot(bind_rows(res)$stop[bind_rows(res)$stop>999])
# mean(bind_rows(res)$stop[bind_rows(res)$stop>999])
# res1=bind_rows(res)$stop
# load('/home/yangz40/Bayesian/mean_change_var_known/results/mean_known_fix_tau_np25.Rdata')
# boxplot(bind_rows(res)$stop[bind_rows(res)$stop>999])
# mean(bind_rows(res)$stop[bind_rows(res)$stop>999])
# res2=bind_rows(res)$stop
# res1-res2
# boxplot(res1[!is.infinite(res1)&res1>1000],res2[!is.infinite(res2)&res2>1000])
# median(res1[!is.infinite(res1)&res1>1000])
# median(res2[!is.infinite(res2)&res2>1000])
# 
# library(microbenchmark)
# library(ggplot2)
# library(microbenchmark)
# library(ggplot2)
# library(ggfortify)  # for autoplot
# 
# set.seed(1)
# 
# k_vals <- c(20,50,100, 150, 200)
# benchmark_results <- list()
# res_exact <- microbenchmark(
#   'exact' = {
#     x=c(rnorm(1000,0,1))
#     OBCD_online_mean_varknown(x, decay=1,prob_change=0.00001, mu0=0,tau0=1,delta0=0,sigma2=1,thresh=Inf,error='error1',K=100000)},times=200)c
# res_exact$K <- "exact"
# benchmark_results[["exact"]] <- res_exact
# 
# for (k in k_vals) {
#   cat("Running benchmark for K =", k, "\n")
#   res_K <- microbenchmark(
#     'ours' ={
#       x=c(rnorm(1000,0,1))
#       OBCD_online_mean_varknown(x, decay=1,prob_change=0.00001, mu0=0,tau0=1,delta0=0,sigma2=1,thresh=Inf,error='error1',K=k)},
#     'benchmark' ={
#       x=c(rnorm(1000,0,1))
#       OBCD_online_mean_varknown(x, decay=1,prob_change=0.00001, mu0=0,tau0=1,delta0=0,sigma2=1,thresh=Inf,error='error3',K=k)},
#     times=200)
#   res_K$K <- k  # tag each result with K value
#   benchmark_results[[as.character(k)]] <- res_K
# }
# 
# library(dplyr)
# library(purrr)
# 
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
#   "ours_KL" = "#466CA6",
#   "ours_HD" = "#D1AE45",
#   "ours" = "#872410",
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
# ggsave("runtime_mean_knownvar.pdf", width = 6, height = 4, units = "in")
