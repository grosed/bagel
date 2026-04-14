# Installing and Using the bagelR Package

## Installing bagelR

The **bagelR** package can be installed via an interactive **R** session using the following code.


```R
if(system.file(package = "bagelR") != "")
    {
      remove.packages("bagelR")
    }
library("remotes")
install_github("grosed/bagel/bagelR",force=TRUE)
```

    Removing package from ‘/home/grosed/DASS/R-packages’
    (as ‘lib’ is unspecified)
    
    Using GitHub PAT from the git credential store.
    
    Downloading GitHub repo grosed/bagel@HEAD
    


    


    Running `R CMD build`...
    


    * checking for file ‘/tmp/RtmpdoIvEV/remotes73c1c77f33c/grosed-bagel-c671288/bagelR/DESCRIPTION’ ... OK
    * preparing ‘bagelR’:
    * checking DESCRIPTION meta-information ... OK
    * cleaning src
    * installing the package to process help pages
    * saving partial Rd database
    * cleaning src
    * checking for LF line-endings in source and make files and shell scripts
    * checking for empty or unneeded directories
    * building ‘bagelR_1.1.0.tar.gz’


    Installing package into ‘/home/grosed/DASS/R-packages’
    (as ‘lib’ is unspecified)
    


## Check package version


```R
packageVersion("bagelR")
```


    [1] ‘1.1.0’


## Using bagelR - 1 : Detecting Change in Mean with Known Variance

#### Load the bagel package


```R
library(bagelR)
```

    
    Attaching package: ‘bagelR’
    
    
    The following object is masked from ‘package:stats’:
    
        time
    
    


### Example 1 - no change

#### Define the feature vector function 

These functions for the **feature vector** and the **prior** correspeond to example 1 in "**Bagel : A Fast Baysian Online Changepoint Detection Algorithm for Linear Models**" *Fearnhead et al.*


```R
feature_vector <- function(t,tau)
{
  if(tau == 0)
  {
    M <- matrix(c(0),nc=1,nr=1)       
    M[1,1] <- 1
    return(M)
  }
  else
  {
    M <- matrix(c(0),nc=1,nr=2)       
    M[1,1] <- 1
    if(t > tau)
    {
      M[2,1] <- 1
    }
    return(M)
  }
}
```

#### Define the prior function


```R

prior <- function(t)
{
  mu1 <- 1;
  delta <- c(5,1)
  if(t == 1)
  {
    mu <- matrix(c(mu1),nc=1,nr=1)
    sigma <- delta[1]*diag(1)
    return(list("mu" = mu,"sigma" = sigma))
  }
  mu <- matrix(c(mu1,0),nc=1,nr=2)
  sigma <- diag(delta)
  return(list("mu" = mu,"sigma" = sigma))
}
```

#### Generate some data - no change


```R
set.seed(0)
Y <- rnorm(1000,0,1)
plot(Y)
```


    
![png](output_16_0.png)
    


#### Set model parameters 


```R
p <- 1.0;
p0 <- 0.9;
noise_sd <- 1.0;
max_particles <- 250
threshold <- 0.8
```

#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_22_0.png)
    


### Example 2 - single change in mean

#### Generate some data - change at t = 501 


```R
set.seed(0)
Y <- c(rnorm(500,0,1),rnorm(500,1,1))
plot(Y)
```


    
![png](output_25_0.png)
    


#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_29_0.png)
    


#### Determine time at which change was detected


```R
time(results)
```


513


#### Visulaise the weights after detection of change


```R
weights <- weights(results)
taus <- taus(results)
plot(taus,weights)
```


    
![png](output_33_0.png)
    


### Determine tau for maximum weight post change detection


```R
taus[which.max(weights)]
```


501


#### Get ratios


```R
ratios <- ratios(results)
```

### Example 3 - single change in mean different prior

#### Define new prior


```R
prior <- function(t)
{
  mu1<-1;
  delta1<-1
  if(t == 1)
  {
    mu <- matrix(c(mu1),nc=1,nr=1)
    sigma <- delta1*diag(1)
    return(list("mu" = mu,"sigma" = sigma))
  }
  mu <- matrix(c(mu1,0),nc=1,nr=2)
  sigma <- delta1*matrix(c(1,-1,-1,2),nr=2,nc=2)
  return(list("mu" = mu,"sigma" = sigma))
}
```

#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_44_0.png)
    


#### Determine time at which change was detected


```R
time(results)
```


509


#### Visulaise the weights after detection of change


```R
weights <- weights(results)
taus <- taus(results)
plot(taus,weights)
```


    
![png](output_48_0.png)
    


### Determine tau for maximum weight post change detection 


```R
taus[which.max(weights)]
```


501


#### Get ratios


```R
ratios <- ratios(results)
```

## Using bagelR - 2 : Detecting Change in Slope with Known Variance

### Example 1 - no change

#### Generate some data - no change


```R
set.seed(0)
Y <- rnorm(1000,0,1)
plot(Y)
```


    
![png](output_56_0.png)
    


#### Define the feature vector function 

These functions for the **feature vector** and the **prior** correspeond to example 2 in "**Bagel : A Fast Baysian Online Changepoint Detection Algorithm for Linear Models**" *Fearnhead et al.*


```R
feature_vector <- function(t,tau)
{
  if(tau == 0)
  {
     M <- matrix(c(0),nc=1,nr=2)	  
     M[1,1] <- 1
     M[2,1] <- t
     return(M)
  }
  else
  {
    M <- matrix(c(0),nc=1,nr=4)	  
    M[1,1] <- 1
    M[2,1] <- t
    if(t > tau)
    {
      M[3,1] <- 1
      M[4,1] <- t
    }
  return(M)
  }
}
```

#### Define the prior function


```R
prior <- function(t)
{
  if(t == 1)
  {
     mu <- matrix(c(0),nc=1,nr=2)
     sigma <- diag(2)
     return(list("mu" = mu,"sigma" = sigma))
  }
  mu <- matrix(c(0),nc=1,nr=4)
  sigma <- diag(4)
  sigma[3,3] <- (t-1)*(t-1)
  sigma[4,3] <- sigma[3,4] <- -(t-1)
  return(list("mu" = mu,"sigma" = sigma))
}
```

#### Set model parameters 


```R
p <- 1.0;
p0 <- 0.9;
noise_sd <- 1.0;
max_particles <- 250
threshold <- 0.8
```

#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_67_0.png)
    


### Example 2 - single change in slope

#### Generate data with change in slope at t = 1000


```R
set.seed(0)
Y <- c(rnorm(1000,0,1),seq(1,200,1)/50 + rnorm(200,0,1))
plot(Y)
```


    
![png](output_70_0.png)
    


#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_74_0.png)
    


#### Determine time at which change was detected


```R
time(results)
```


1058


#### Visulaise the weights after detection of change


```R
weights <- weights(results)
taus <- taus(results)
plot(taus,weights)
```


    
![png](output_78_0.png)
    


#### Locate the maximum weight for tau > 0


```R
taus[which.max(weights[-1])]
```


1020


#### Get ratios


```R
ratios <- ratios(results)
```
