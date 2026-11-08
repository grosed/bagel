library(bagelR)

feature <- function(t,tau) matrix(if(tau==0) 1 else c(1,1),ncol=1)
prior <- function(t) {
  if(t==1) list(mu=matrix(0,1,1),sigma=matrix(1,1,1))
  else list(mu=matrix(0,2,1),sigma=matrix(c(1,-1,-1,2),2,2))
}
transform <- function(tau) matrix(if(tau==0) 1 else c(1,1),nrow=1)
fit <- function(y, threshold = 1, ...) bagel_kv_exact(y,feature,prior,transform,.9,1,1,10,threshold,...)

extreme <- fit(c(0,1000,1001,999,1000))
stopifnot(all(is.finite(weights_tau_zero(extreme))),
          all(is.finite(weights(extreme))), abs(sum(weights(extreme))-1)<1e-12)

t <- seq_len(200)
strong <- fit(.2*sin(t)+ifelse(t>50,10,0))
stopifnot(tail(weights_tau_zero(strong),1)==0,
          bagelR::time(strong)==length(t),
          is.na(changepoint(strong)$detected))

alarm <- fit(.2*sin(t)+ifelse(t>50,10,0),threshold=.95)
stopifnot(bagelR::time(alarm)<length(t),
          tail(1-weights_tau_zero(alarm),1)>.95,
          all(head(1-weights_tau_zero(alarm),-1)<=.95),
          alarm@threshold==.95, changepoint(alarm)$detected==bagelR::time(alarm))

for(mode in c("bagel","LP-M","LP-P")) {
  zero <- bagel_kv_exact(c(0,.1,-.1,.2,0),feature,prior,transform,1,1,1,3,1,pruning=mode)
  stopifnot(all(is.finite(unlist(ratios(zero)))),all(is.finite(weights(zero))),
            weights_tau_zero(zero)==1)
}

sequential <- sequential_bagel_kv_exact(feature,prior,transform,.9,1,1,10)
bagelR::update(sequential,0)
bagelR::update(sequential,.1)
module <- sequential@bagel_object
stopifnot(abs(log(weights_tau_zero(sequential))-module$get_log_w0())<1e-14)
module$set_log_weights(c(-1000,0))
stopifnot(module$get_log_w0()== -1000,weights_tau_zero(sequential)==0)
set_weights(sequential,c(.25,.75))
stopifnot(abs(module$get_log_w0()-log(.25))<1e-14)
stopifnot(inherits(try(module$set_log_weights(c(0)),silent=TRUE),"try-error"))

failure <- try(fit(c(0,1e200)),silent=TRUE)
stopifnot(inherits(failure,"try-error"),grepl("All particle weights are zero",failure))
cat("Log-weight regression tests passed.\n")
