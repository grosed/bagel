####BASIC CODE FOR BAGEL WITHOUT THE RESAMPLING
####ASSUMING SIGMA IS KNOWN

###SET-UP FOR EXAMPLE 2

##definition of vectors a and b
a<-function(t){
  return(matrix(c(1,t),nc=1))
}

b<-function(t,tau){
  if(t<=tau){ 
    return(matrix(0,nc=1,nr=2))
  }else{
    return(matrix(c(1,t),nc=1))
  }
}

##definition of prior
delta<-c(1,1,1) ###hyper-parameter of prior this is defined outside the functions that use them

##definition of matrices in THM2 at time t -- these are for tau=t-1
Sigma.beta <- diag(delta[1:2]) ##does not depend on t or tau
Sigma.beta.gamma <- function(t){
  return(diag(c(0,0)))
} 
Sigma.gamma <- function(t){
  return(matrix(delta[3]*c((t-1)^2,-(t-1),-(t-1),1),nc=2))
}
mu.beta <- matrix(c(0,0),nc=1)
##again this is prior mean for tau=t-1
mu.gamma <- function(t){
  return(matrix(c(0,0),nc=1))
}

##hyper-parameters
sigma<-1 ## sd of noise
p<-1 ##parameter for prior on change location. p=1 gives uniform prior
p0 <-0.9 ##prob of  no change at any time 
##data
y<-rnorm(100)

##initialise
d1<-length(mu.beta)
d2<-length(mu.gamma(1))
n<-length(y)

###at time t we have t particles, to simplify storage we will set up an array for the parameters
##not the best approach for large n or when we have pruning
mu.post <-array(0,dim=c(d1+d1,n))
Sigma.post <-array(0,dim=c(d1+d2,d1+d2,n)) ##entry for no-change, then change at tau=1,2,...,n-1, 

mu.post[1:d1,1]<-mu.beta
Sigma.post[1:d1,1:d1,1]<-Sigma.beta

logw<-0 ##log of sum of weights -- not sure if needed to store
#t=1, our approach assumes y_1 is from the pre-change model
t <- 1
var.pred <- sigma^2*(1 + t(a(t))%*%Sigma.beta%*%a(t)) ##prediction variance
w <- dnorm(y[1],t(a(t))%*%mu.beta,sqrt(var.pred))
##normalise weights
logw<-logw+log(sum(w))
w<-w/sum(w)

##THM 3 but with t in place of t-1, some issue with THM 3 with transpose
h.i <- a(t)
e.i <- as.numeric(y[t] - t(h.i) %*% mu.post[1:d1,i])  ##I think e_i is this
Q <- as.numeric(t(h.i) %*% Sigma.post[1:d1,1:d1,i] %*% h.i +1)
A <- (1/Q)* (Sigma.post[1:d1,1:d1,i] %*% h.i)
##for no-change model only update the beta components
Sigma.post[1:d1,1:d1,i] <- Sigma.post[1:d1,1:d1,i] - Q*((A)%*%t(A)) ##I think the update for Sigma is this
mu.post[1:d1,i] <- mu.post[1:d1,i]+A*e.i

##update t
t<-t+1
while(t<=n){
  ##THM 2 to calculate prior for change at tau=t-1 -- which we store in entry t
  B <- t(Sigma.beta.gamma(t)) %*% solve(Sigma.beta)
  mu.post[,t] <- c(mu.post[1:d1,1],mu.gamma(t)+ B %*%(mu.post[1:d1,1]-mu.beta))
  Sigma.post[1:d1,1:d1,t] <- Sigma.post[1:d1,1:d1,1]
  Sigma.post[1:d1,d1+1:d2,t] <- B %*% Sigma.post[1:d1,1:d1,1]
  Sigma.post[d1+1:d2,1:d1,t] <- t(Sigma.post[1:d1,d1+1:d2,t])
  Sigma.post[d1+1:d2,d1+1:d2,t] <- Sigma.gamma(t) + B%*%(Sigma.post[1:d1,1:d1,1]-Sigma.beta)%*%t(B)
  
  ###UPDATE weights based on prior from THM1 without likelihood
  if(t==2){##at any time we have prob p0 of no change
    w<-c(p0,1-p0)
  }else{
    if(p<1){#geometric prior
      w<-c(w[1],w[2:(t-1)]*p*(1-p^(t-2))/(1-p^(t-1)),w[1]*(1-p0)*(1-p)/((1-p^(t-1))*p0))
    }else{##uniform -- we could just use this
      w<-c(w[1],w[2:(t-1)]*(t-2)/(t-1),w[1]*(1-p0)/(p0*(t-1)))
    }
  }
  ###THM 4 to update weights
  for(i in 1:t){
    if(i==1){
      h.i <- a(t)
      var.pred <- sigma^2 * (1+ t(h.i)%*%Sigma.post[1:d1,1:d1,i]%*%h.i)
      w[i] <- w[i] * dnorm(y[t],t(h.i)%*%mu.post[1:d1,i],sqrt(var.pred))
    }else{
      h.i<- c(a(t),b(t,i-1))
      var.pred <- sigma^2 * (1+ t(h.i)%*%Sigma.post[,,i]%*%h.i)
      w[i] <- w[i] * dnorm(y[t],t(h.i)%*%mu.post[,i],sqrt(var.pred))
      
    }
  }
  
  ##normalise weights
  logw<-logw+log(sum(w))
  w<-w/sum(w)
  ###THIS IS WHERE YOU WOULD CHECK IS THERE IS EVIDENCE FOR A CHANGE/STOP THE ALGORITHM
  
  ###THM 3 to update posterior parameters
  for(i in 1:t){
    if(i==1){
      h.i <- a(t)
      e.i <- as.numeric(y[t] - t(h.i) %*% mu.post[1:d1,i])  ##I think e_i is this
      Q <- as.numeric(t(h.i) %*% Sigma.post[1:d1,1:d1,i] %*% h.i +1)
      A <- (1/Q)* (Sigma.post[1:d1,1:d1,i] %*% h.i)
      ##for no-change model only update the beta components
      Sigma.post[1:d1,1:d1,i] <- Sigma.post[1:d1,1:d1,i] - Q*((A)%*%t(A)) ##I think the update for Sigma is this
      mu.post[1:d1,i] <- mu.post[1:d1,i]+A*e.i
    }else{
      h.i<- c(a(t),b(t,i-1))
      e.i <- as.numeric(y[t] - t(h.i) %*% mu.post[,i])  ##I think e_i is this
      Q <- as.numeric(t(h.i) %*% Sigma.post[,,i] %*% h.i +1)
      A <- (1/Q)* (Sigma.post[,,i] %*% h.i)
      ##for no-change model only update the beta components
      Sigma.post[,,i] <- Sigma.post[,,i] - Q*((A)%*%t(A)) ##I think the update for Sigma is this
      mu.post[,i] <- mu.post[,i]+A*e.i
    }
  }
  ##THIS IS WHERE YOU WOULD PRUNE 
  ##THE ABOVE WOULD BE MAINLY UNCHANGED EXCEPT WE HAVE USED THAT tau associate with entry i is tau=i-1 in definition of h.i
  
  ##increment t
  t<-t+1
}
