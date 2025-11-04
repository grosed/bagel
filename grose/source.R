library(Rcpp)
library(parallel)
sourceCpp('find_min_err_index_list_matrix.cpp')
#sourceCpp('~/Bayesian/mean_change_var_known/find_min_err_index_list_matrix.cpp')
#sourceCpp('~/Bayesian/mean_change_var_known/approx_pnorm.cpp')































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
