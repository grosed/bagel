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
    


    * checking for file ‘/tmp/RtmpAXU8UQ/remotes7eed06a94c8fa/grosed-bagel-8d88d25/bagelR/DESCRIPTION’ ... OK
    * preparing ‘bagelR’:
    * checking DESCRIPTION meta-information ... OK
    * cleaning src
    * installing the package to process help pages
    * saving partial Rd database
    * cleaning src
    * checking for LF line-endings in source and make files and shell scripts
    * checking for empty or unneeded directories
    * building ‘bagelR_1.2.0.tar.gz’


    Installing package into ‘/home/grosed/DASS/R-packages’
    (as ‘lib’ is unspecified)
    


## Check package version


```R
packageVersion("bagelR")
```


    [1] ‘1.2.0’


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

#### Check for a changepoint


```R
changepoint(results)
```


&lt;NA&gt;


#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_24_0.png)
    


### Example 2 - single change in mean

#### Generate some data - change at t = 501 


```R
set.seed(0)
Y <- c(rnorm(500,0,1),rnorm(500,1,1))
plot(Y)
```


    
![png](output_27_0.png)
    


#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Check for a changepoint


```R
changepoint(results)
```


<dl>
	<dt>$location</dt>
		<dd>501</dd>
	<dt>$detected</dt>
		<dd>513</dd>
</dl>



#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_33_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -25
end <- 25
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_35_0.png)
    


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

#### Check for a changepoint


```R
changepoint(results)
```


<dl>
	<dt>$location</dt>
		<dd>501</dd>
	<dt>$detected</dt>
		<dd>509</dd>
</dl>



#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_44_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -25
end <- 25
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_46_0.png)
    


## Using bagelR - 2 : Detecting Change in Slope with Known Variance

### Example 1 - no change

#### Generate some data - no change


```R
set.seed(0)
Y <- rnorm(1000,0,1)
plot(Y)
```


    
![png](output_50_0.png)
    


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

#### Check for a changepoint


```R
changepoint(results)
```


&lt;NA&gt;


#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_63_0.png)
    


### Example 2 - single change in slope

#### Generate data with change in slope at t = 1000


```R
set.seed(0)
Y <- c(rnorm(1000,0,1),seq(1,200,1)/50 + rnorm(200,0,1))
plot(Y)
```


    
![png](output_66_0.png)
    


#### Offline analysis


```R
results <- bagel_offline(feature_vector,prior,p0,p,noise_sd,max_particles,threshold,Y)
```

#### Check for a changepoint


```R
changepoint(results)
```


<dl>
	<dt>$location</dt>
		<dd>1021</dd>
	<dt>$detected</dt>
		<dd>1058</dd>
</dl>



#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_72_0.png)
    



```R
cpt = changepoint(results)$location
start <- -30
end <- 30
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_73_0.png)
    


### Other methods

### Get the taus

- Returns a vector of tau values for the active (unpruned) particles


```R
taus(results)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0</li><li>1</li><li>2</li><li>3</li><li>4</li><li>5</li><li>6</li><li>7</li><li>8</li><li>9</li><li>10</li><li>11</li><li>12</li><li>13</li><li>14</li><li>15</li><li>16</li><li>17</li><li>18</li><li>19</li><li>20</li><li>21</li><li>22</li><li>23</li><li>24</li><li>25</li><li>26</li><li>27</li><li>28</li><li>29</li><li>30</li><li>31</li><li>32</li><li>33</li><li>34</li><li>35</li><li>36</li><li>37</li><li>38</li><li>39</li><li>40</li><li>41</li><li>42</li><li>99</li><li>111</li><li>128</li><li>164</li><li>180</li><li>212</li><li>270</li><li>366</li><li>378</li><li>394</li><li>418</li><li>434</li><li>450</li><li>474</li><li>522</li><li>538</li><li>546</li><li>554</li><li>562</li><li>570</li><li>578</li><li>586</li><li>594</li><li>602</li><li>610</li><li>618</li><li>626</li><li>634</li><li>642</li><li>650</li><li>658</li><li>666</li><li>674</li><li>682</li><li>690</li><li>698</li><li>706</li><li>714</li><li>722</li><li>730</li><li>738</li><li>746</li><li>754</li><li>762</li><li>770</li><li>778</li><li>786</li><li>794</li><li>802</li><li>810</li><li>818</li><li>826</li><li>834</li><li>838</li><li>842</li><li>846</li><li>850</li><li>854</li><li>858</li><li>862</li><li>866</li><li>870</li><li>874</li><li>878</li><li>882</li><li>886</li><li>890</li><li>892</li><li>894</li><li>896</li><li>898</li><li>900</li><li>902</li><li>904</li><li>906</li><li>908</li><li>910</li><li>912</li><li>914</li><li>916</li><li>918</li><li>920</li><li>922</li><li>924</li><li>926</li><li>928</li><li>930</li><li>932</li><li>934</li><li>936</li><li>938</li><li>940</li><li>942</li><li>944</li><li>945</li><li>946</li><li>947</li><li>948</li><li>949</li><li>950</li><li>951</li><li>952</li><li>953</li><li>954</li><li>955</li><li>956</li><li>957</li><li>958</li><li>959</li><li>960</li><li>961</li><li>962</li><li>963</li><li>964</li><li>965</li><li>966</li><li>967</li><li>968</li><li>969</li><li>970</li><li>971</li><li>972</li><li>973</li><li>974</li><li>975</li><li>976</li><li>977</li><li>978</li><li>979</li><li>980</li><li>981</li><li>982</li><li>983</li><li>984</li><li>985</li><li>986</li><li>987</li><li>988</li><li>989</li><li>990</li><li>991</li><li>992</li><li>993</li><li>994</li><li>995</li><li>996</li><li>997</li><li>998</li><li>999</li><li>1000</li><li>1001</li><li>1002</li><li>1003</li><li>1004</li><li>1005</li><li>1006</li><li>1007</li><li>1008</li><li>1009</li><li>1010</li><li>1011</li><li>1012</li><li>1013</li><li>1014</li><li>1015</li><li>1016</li><li>1017</li><li>1018</li><li>1019</li><li>1020</li><li>1021</li><li>1022</li><li>1023</li><li>1024</li><li>1025</li><li>1026</li><li>1027</li><li>1028</li><li>1029</li><li>1030</li><li>1031</li><li>1032</li><li>1033</li><li>1034</li><li>1035</li><li>1036</li><li>1037</li><li>1038</li><li>1039</li><li>1040</li><li>1041</li><li>1042</li><li>1043</li><li>1044</li><li>1045</li><li>1046</li><li>1047</li><li>1048</li><li>1049</li><li>1050</li><li>1051</li><li>1052</li><li>1053</li><li>1054</li><li>1055</li><li>1056</li><li>1057</li></ol>



### Get the ratios

- Returns a list of vectors containing the ratios for each active (unpruned) particle


```R
ratios <- ratios(results)
ratios[[44]]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0.0369157342762245</li><li>0.0360736863222321</li><li>0.0351171334175432</li><li>0.0340261451486438</li><li>0.032882050595927</li><li>0.0315529278212478</li><li>0.0305388558344099</li><li>0.0295831346010051</li><li>0.0225092765657438</li><li>0.021395803878344</li><li>0.0212786090318044</li><li>0.020635315826559</li><li>0.0203756645109818</li><li>0.0200234984110007</li><li>0.0197582226568702</li><li>0.0196731633454171</li><li>0.0190305428298988</li><li>0.0188476945170813</li><li>0.0185349521026749</li><li>0.0181708937053113</li><li>0.0178497480631719</li><li>0.0175585735870665</li><li>0.0172215112397935</li><li>0.0168765858823221</li><li>0.016493380085714</li><li>0.0161171573626624</li><li>0.0158138312184265</li><li>0.0151615688280308</li><li>0.0147785813612674</li><li>0.0145393785422443</li><li>0.0145127220726176</li><li>0.0141507046390281</li><li>0.0141331999648893</li><li>0.014004683834216</li><li>0.0139712479747545</li><li>0.0137231004430856</li><li>0.0134535064716414</li><li>0.01316198239605</li><li>0.0129129122425529</li><li>0.0126753580628231</li><li>0.0125541541989876</li><li>0.0124177104905518</li><li>0.0122797018343779</li><li>0.0121094147822408</li><li>0.0119520050871808</li><li>0.0117872654300719</li><li>0.011624231893999</li><li>0.0114839564380647</li><li>0.011350453500698</li><li>0.0112272845599675</li><li>0.011102185898862</li><li>0.010983511582391</li><li>0.0108601796748974</li><li>0.0107395217647639</li><li>0.0106309383577982</li><li>0.010431000559458</li><li>0.0104334442744107</li></ol>



### Get weights of active particles

- Returns a vector of weights for the currently active (unpruned) particles


```R
active_weights(results)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0.0855637856099427</li><li>5.19743464008564e-06</li><li>3.90384686711566e-06</li><li>2.71816302034523e-06</li><li>2.31712570606257e-06</li><li>2.3076920929784e-06</li><li>2.29018928131473e-06</li><li>1.79375529559154e-06</li><li>1.32377023605271e-06</li><li>1.00816104791613e-06</li><li>8.02588874135339e-07</li><li>7.38945178340599e-07</li><li>7.02243020922234e-07</li><li>6.40701366834082e-07</li><li>5.60287307749622e-07</li><li>4.90287372025317e-07</li><li>4.29985691563115e-07</li><li>3.77340586032686e-07</li><li>3.36948926343934e-07</li><li>2.98052814729899e-07</li><li>2.68684251013336e-07</li><li>2.39267787618365e-07</li><li>2.14681095512342e-07</li><li>1.95408739921179e-07</li><li>1.79410485769523e-07</li><li>1.67087773433268e-07</li><li>1.56105664484728e-07</li><li>1.47020636986539e-07</li><li>1.4012758316507e-07</li><li>1.32971968142954e-07</li><li>1.25154137011956e-07</li><li>1.1817570450591e-07</li><li>1.11695609329369e-07</li><li>1.05478385852784e-07</li><li>9.96299073013507e-08</li><li>9.40435360879632e-08</li><li>8.93648133688222e-08</li><li>8.55728969180521e-08</li><li>8.24320205513305e-08</li><li>7.93343619903732e-08</li><li>7.68399016888708e-08</li><li>7.43888166398823e-08</li><li>7.25900249677347e-08</li><li>1.92409943186823e-06</li><li>2.20225415139368e-07</li><li>2.60605291773494e-07</li><li>4.00885641989203e-07</li><li>1.49270068306689e-07</li><li>2.58855352967131e-07</li><li>4.15685024590275e-07</li><li>7.25408573312886e-07</li><li>8.69298650308779e-08</li><li>1.19813736706894e-07</li><li>1.97004706497576e-07</li><li>1.4652154836486e-07</li><li>1.60106632913877e-07</li><li>2.56228508811661e-07</li><li>6.51967289461094e-07</li><li>2.49783151682162e-07</li><li>1.33123356263719e-07</li><li>1.38126699614058e-07</li><li>1.43909540778874e-07</li><li>1.5214916585372e-07</li><li>1.61432768388641e-07</li><li>1.6705231353283e-07</li><li>1.69633901121588e-07</li><li>1.69809429181448e-07</li><li>1.74204378261086e-07</li><li>1.76995866568913e-07</li><li>1.79186720501849e-07</li><li>1.78704644799418e-07</li><li>1.79353393665684e-07</li><li>1.84826639828951e-07</li><li>1.91968115404506e-07</li><li>1.99441763829544e-07</li><li>2.07911651340076e-07</li><li>2.13576873185532e-07</li><li>2.17043744967799e-07</li><li>2.21588656583486e-07</li><li>2.29351286932304e-07</li><li>2.37773216602116e-07</li><li>2.43408612422975e-07</li><li>2.49917429951732e-07</li><li>2.50702710503117e-07</li><li>2.4804176383523e-07</li><li>2.5048392225929e-07</li><li>2.57974043017596e-07</li><li>2.76635063607054e-07</li><li>2.95001655974499e-07</li><li>3.16524994943358e-07</li><li>3.52848933433258e-07</li><li>4.10222327211216e-07</li><li>4.63603832348892e-07</li><li>5.39632111425854e-07</li><li>6.19603197583372e-07</li><li>7.02966795658939e-07</li><li>3.83447965805362e-07</li><li>4.15055883626084e-07</li><li>4.57305584505923e-07</li><li>5.18675201241884e-07</li><li>5.85455694923097e-07</li><li>6.52465283158608e-07</li><li>7.24045082408897e-07</li><li>8.25726672768395e-07</li><li>9.65555969539259e-07</li><li>1.13610195877846e-06</li><li>1.34902467408966e-06</li><li>1.58749027468961e-06</li><li>1.89734801553438e-06</li><li>2.27657089920758e-06</li><li>1.23743890280051e-06</li><li>1.32557350336651e-06</li><li>1.40413451237877e-06</li><li>1.50382243672475e-06</li><li>1.61605890154258e-06</li><li>1.72417592083199e-06</li><li>1.83731325751e-06</li><li>1.97768945956497e-06</li><li>2.1305005686952e-06</li><li>2.25717157007011e-06</li><li>2.43773390335471e-06</li><li>2.66367372938013e-06</li><li>2.91799513671797e-06</li><li>3.24113278117745e-06</li><li>3.65549458992972e-06</li><li>4.11179258898575e-06</li><li>4.71090515653807e-06</li><li>5.3779167909755e-06</li><li>6.0816269309293e-06</li><li>7.00290091889062e-06</li><li>8.33861185702633e-06</li><li>9.84736025666228e-06</li><li>1.16475070606288e-05</li><li>1.38977333833596e-05</li><li>1.68199126807267e-05</li><li>2.02993401375041e-05</li><li>2.49578108053684e-05</li><li>1.39100790009925e-05</li><li>1.53352034832266e-05</li><li>1.68607261658317e-05</li><li>1.85416415258642e-05</li><li>2.02726903735378e-05</li><li>2.2311932542649e-05</li><li>2.42571640850705e-05</li><li>2.63863525893956e-05</li><li>2.8451602873201e-05</li><li>3.09948649032785e-05</li><li>3.36934624771299e-05</li><li>3.70019087285439e-05</li><li>4.05827859184176e-05</li><li>4.45201393360191e-05</li><li>4.85324937623724e-05</li><li>5.16746338763743e-05</li><li>5.58232912069149e-05</li><li>5.96661160955412e-05</li><li>6.52576334139291e-05</li><li>7.17016260802448e-05</li><li>7.85322507172133e-05</li><li>8.67163515127849e-05</li><li>9.67416973270705e-05</li><li>0.000107957936259731</li><li>0.000120931126176391</li><li>0.00013832504163205</li><li>0.00015652561394846</li><li>0.00017756247342631</li><li>0.000203759826173803</li><li>0.000231089101894346</li><li>0.000257711641855123</li><li>0.00029250502690569</li><li>0.000330124425343699</li><li>0.000369034615656759</li><li>0.000413285035473904</li><li>0.00046073860180659</li><li>0.000511075198409551</li><li>0.000557593328917338</li><li>0.000601702228532859</li><li>0.000657095457867935</li><li>0.000723251679219596</li><li>0.000789928022237195</li><li>0.00089418768172427</li><li>0.00103731894978273</li><li>0.00119433291191719</li><li>0.00137752081087418</li><li>0.00151344049667592</li><li>0.00169035860345791</li><li>0.0019086778163354</li><li>0.00214095245251943</li><li>0.00239692665188916</li><li>0.00269134685829706</li><li>0.00308440073418274</li><li>0.00347790025905576</li><li>0.00401096958720477</li><li>0.00470983825073261</li><li>0.00546871969629779</li><li>0.00633062628368218</li><li>0.00766171495958668</li><li>0.0092996386196859</li><li>0.0109968971998205</li><li>0.0126123555473879</li><li>0.0144992565279259</li><li>0.0168523757324487</li><li>0.0198371093813007</li><li>0.0233452566293649</li><li>0.0261397358112377</li><li>0.029012015677814</li><li>0.0332943885535601</li><li>0.035725480880395</li><li>0.0392549350908904</li><li>0.0412235129514488</li><li>0.0430765779725741</li><li>0.0428704256496632</li><li>0.0440316212827238</li><li>0.0443077022339246</li><li>0.0454627834545923</li><li>0.0452731338460102</li><li>0.0428912577410694</li><li>0.0378878046277294</li><li>0.0315108279430773</li><li>0.0255184991221867</li><li>0.0219413928649518</li><li>0.019333573084094</li><li>0.0163935949860657</li><li>0.0140133388629269</li><li>0.0112304655000818</li><li>0.00857969391495955</li><li>0.00653640331801604</li><li>0.00590264514557588</li><li>0.00544448665780643</li><li>0.00510263331992912</li><li>0.00401524750285734</li><li>0.0034653706926633</li><li>0.0033657866302441</li><li>0.00285113415235325</li><li>0.00242476782003675</li><li>0.00188880856068107</li><li>0.0014839217901611</li><li>0.00108830885560241</li><li>0.000909575717240043</li><li>0.000795003267347645</li><li>0.000721623866254783</li><li>0.000732925336651902</li><li>0.000531626312952514</li><li>0.0002814028418961</li><li>0.000195651116737345</li><li>0.000138716739894694</li><li>0.000110342079516142</li><li>4.37151007557838e-05</li><li>4.24629012259596e-05</li><li>2.71380102842447e-05</li><li>1.78502254229239e-05</li></ol>


