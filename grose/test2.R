library(bagel)



load("inputs")

result <- rlang::exec(bagel_mean,!!!inputs[[2]])

save(file="results",list=c("result"))


