library(bagel)

load("inputs")

result <- rlang::exec(bagel_mean,!!!inputs[[2]])

save(file="results",list=c("result"))


set.seed(1)
y <- generate_piecewise_linear_data(2000,1000,0,0,-2.2,0.002,1,continuous=T)



result <- bagel_slope(y,p_dis=0.1,p_conti=0.1,decay=1,mu0=c(0,0),
                           delta_1=0.5^2,delta_2=0.5^2,delta_3=0.5,delta_4=1,
                           sigma2=1,thresh=0.89,error ='KL',K=50000,prior='prior1')








