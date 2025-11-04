

find_cover <- function(w, true_change)
{
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