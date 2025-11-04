

run_ours <- function(res,tau)
{
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