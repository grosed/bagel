

split_number <- function(weight, ratios)
{
  x1 = ratios/(1+ratios)*weight
  x2 = 1/(1+ratios)*weight
  return(c(x1,x2))
}