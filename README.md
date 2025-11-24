# Installing and using the bagel package for R

## Installing bagel

The **bagel** package can be installed via an interactive **R** session using the following code.


```R
library("remotes")
install_github("grosed/bagel/R-package",force=TRUE)
```

    Downloading GitHub repo grosed/bagel@HEAD
    
    Running `R CMD build`...
    


    * checking for file ‘/tmp/Rtmp5AmVmC/remotes19adcf658c4dc4/grosed-bagel-538a353/R-package/DESCRIPTION’ ... OK
    * preparing ‘bagel’:
    * checking DESCRIPTION meta-information ... OK
    * checking for LF line-endings in source and make files and shell scripts
    * checking for empty or unneeded directories
    * building ‘bagel_0.0.0.tar.gz’


    Installing package into ‘/home/grosedj1/bagel-env/R-packages’
    (as ‘lib’ is unspecified)
    


Note that this method will install the latest stable version of **bagel** from the **main** branch of the repository.

## Using bagel

The latest version of **bagel** can be used to detect changes in mean or slope. In both cases it is necessary to specify the (known) variance. A future release of the package will also detect changes in mean or slope when the variance is unknown.

### Example 1 - detecting a change in mean

#### load the bagel package


```R
library(bagel)
```

#### generate some test data


```R
set.seed(0)
y <- c(rnorm(250,0,1),rnorm(250,1,1))
plot(y)
```


    
![png](./images/output_11_0.png)
    


#### set model parameters


```R
decay <- 1
probability <- 0.1
mu <- c(0,0)
tau <- 0.015625
delta <- 0
sigma <- 1
threshold <- 0.2
error <- 'error1'
K <- 200
```

#### test for changes in mean


```R
result <- bagel_mean(y,decay,probability,mu,tau,delta,sigma,threshold,error,K)
```

#### location of change


```R
result$stoppingtime
```


266


### Example 2 - detecting a change in slope

#### generate test data


```R
set.seed(1)
y <- generate_piecewise_linear_data(2000,1000,0,0,-2.2,0.002,1,continuous=T)
plot(y)
```


    
![png](./images/output_20_0.png)
    


#### set model parameters


```R
p_dis <- 0.1
p_conti <- 0.1
decay <- 1
mu0 <- c(0,0)
delta_1 <- 0.5^2
delta_2 <- 0.5^2
delta_3 <- 0.5
delta_4 <- 1
sigma2 <- 1
thresh <- 0.89
error <- 'KL'
K <- 50
prior <- 'prior1'
```


```R
result <- bagel_slope(y, p_dis, p_conti, decay, mu0, delta_1, delta_2, delta_3, delta_4,
                      sigma2, thresh, error, K, prior)
```


```R
result$stoppingtime
```


1283

