library(Rcpp)


source("OBCD_online_mean_varknown.R")
source("find_cover.R")
source("generate_h0_data.R")
source("recover_posterior.R")
source("return_arl.R")
source("run_benchmark.R")
source("run_ours.R")
source("run_simulation_add.R")
source("split_number.R")
source("update_post_paramater.R")
source("update_prior.R")


sourceCpp('find_min_err_index_list_matrix.cpp')



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


sigma2 <- 1
n <- 500

x <- rnorm(n,0,sd=sqrt(sigma2))

prob_change <- prob_value

#result <- OBCD_online_mean_varknown(x, decay=pram$decay,prob_change=prob_value, mu0=c(0,0),tau0=pram$tau0,delta0=0,sigma2,thresh=Inf,error=pram$error,K=pram$K)


inputs <- list("basic_test"=list("x"=x,
               "decay"=pram$decay,
               "prob_change"=prob_value,
               "mu0"=c(0,0),
               "tau0"=pram$tau0,
               "delta0"=0,
               "sigma2"=sigma2,
               "thresh"=Inf,
               "error"=pram$error,
               "K"=501))
save(file="inputs",list=c("inputs")) 


result <- OBCD_online_mean_varknown(x, decay=pram$decay,prob_change=prob_value, mu0=c(0,0),tau0=pram$tau0,delta0=0,sigma2,thresh=Inf,error=pram$error,K=501)



save(file="dump",list=ls())





#at_arl <- mclapply(1:500,function(i){
#  set.seed(i)
#  return_arl(n=2000,sigma2=1,decay=pram$decay, prob_change=prob_value,mu0=c(0,0),tau0=pram$tau0,delta0=0,thresh=Inf,error =pram$error,K=pram$K)
#},mc.set.seed = TRUE,mc.cores = detectCores() %/% 2-1)



