

run_benchmark <- function(res,tau)
{
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