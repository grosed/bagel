

recover_posteior <- function(pruned_candidates, stored_candidate, w, ratios)
{
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

