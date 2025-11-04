

run_simulation_add <- function(x,tau,decay,prob_change,mu0,tau0,delta0,sigma2,thresh,error,K)
{
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