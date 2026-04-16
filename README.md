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


## Using bagelR - 1 : Detecting Change in Mean

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

#### Offline analysis - known variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y,noise_sd)
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
    


#### Offline analysis - uknown variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y)
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


    
![png](output_30_0.png)
    


### Example 2 - single change in mean

#### Generate some data - change at t = 501 


```R
set.seed(0)
Y <- c(rnorm(500,0,1),rnorm(500,1,1))
plot(Y)
```


    
![png](output_33_0.png)
    


#### Offline analysis - known variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y,noise_sd)
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


    
![png](output_39_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -25
end <- 25
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_41_0.png)
    


#### Offline analysis - unknown variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y)
```

#### Check for a changepoint


```R
changepoint(results)
```


<dl>
	<dt>$location</dt>
		<dd>501</dd>
	<dt>$detected</dt>
		<dd>512</dd>
</dl>



#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_47_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -25
end <- 25
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_49_0.png)
    


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

#### Offline analysis - known variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y,noise_sd)
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


    
![png](output_58_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -25
end <- 25
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_60_0.png)
    


#### Offline analysis - unknown variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y)
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


    
![png](output_66_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -25
end <- 25
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_68_0.png)
    


## Using bagelR - 2 : Detecting Change in Slope

### Example 1 - no change

#### Generate some data - no change


```R
set.seed(0)
Y <- rnorm(1000,0,1)
plot(Y)
```


    
![png](output_72_0.png)
    


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

#### Offline analysis - known variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y,noise_sd)
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


    
![png](output_85_0.png)
    


#### Offline analysis - unknown variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y)
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


    
![png](output_91_0.png)
    


### Example 2 - single change in slope

#### Generate data with change in slope at t = 1000


```R
set.seed(0)
Y <- c(rnorm(1000,0,1),seq(1,200,1)/50 + rnorm(200,0,1))
plot(Y)
```


    
![png](output_94_0.png)
    


#### Offline analysis - known variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y,noise_sd)
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


    
![png](output_100_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -30
end <- 30
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_102_0.png)
    


#### Offline analysis - unknown variance


```R
results <- bagel_offline(feature_vector,prior,p0,p,max_particles,threshold,Y)
```

#### Check for a changepoint


```R
changepoint(results)
```


<dl>
	<dt>$location</dt>
		<dd>1032</dd>
	<dt>$detected</dt>
		<dd>1058</dd>
</dl>



#### Visualise weight at tau = 0


```R
plot(weights_tau_zero(results))
```


    
![png](output_108_0.png)
    


#### Plot weights in region of changepoint


```R
cpt = changepoint(results)$location
start <- -30
end <- 30
plot(seq(start,end,1),weights(results)[-1][(cpt+start):(cpt+end)])
```


    
![png](output_110_0.png)
    


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
<ol class=list-inline><li>0</li><li>1</li><li>2</li><li>3</li><li>4</li><li>5</li><li>6</li><li>7</li><li>8</li><li>9</li><li>10</li><li>11</li><li>14</li><li>16</li><li>18</li><li>21</li><li>22</li><li>23</li><li>24</li><li>25</li><li>26</li><li>27</li><li>28</li><li>30</li><li>31</li><li>32</li><li>33</li><li>34</li><li>35</li><li>36</li><li>37</li><li>38</li><li>39</li><li>40</li><li>42</li><li>43</li><li>44</li><li>45</li><li>46</li><li>47</li><li>49</li><li>50</li><li>51</li><li>52</li><li>56</li><li>58</li><li>70</li><li>74</li><li>79</li><li>86</li><li>97</li><li>101</li><li>103</li><li>104</li><li>105</li><li>106</li><li>107</li><li>109</li><li>113</li><li>118</li><li>122</li><li>123</li><li>124</li><li>125</li><li>128</li><li>135</li><li>142</li><li>146</li><li>147</li><li>148</li><li>149</li><li>150</li><li>151</li><li>152</li><li>153</li><li>158</li><li>160</li><li>161</li><li>162</li><li>185</li><li>191</li><li>204</li><li>211</li><li>213</li><li>214</li><li>218</li><li>222</li><li>224</li><li>227</li><li>229</li><li>235</li><li>240</li><li>241</li><li>242</li><li>243</li><li>244</li><li>263</li><li>265</li><li>357</li><li>365</li><li>372</li><li>380</li><li>382</li><li>384</li><li>385</li><li>386</li><li>387</li><li>389</li><li>390</li><li>403</li><li>411</li><li>424</li><li>426</li><li>427</li><li>432</li><li>434</li><li>435</li><li>436</li><li>438</li><li>453</li><li>560</li><li>608</li><li>725</li><li>759</li><li>764</li><li>767</li><li>771</li><li>781</li><li>804</li><li>816</li><li>828</li><li>835</li><li>844</li><li>852</li><li>856</li><li>858</li><li>859</li><li>861</li><li>862</li><li>870</li><li>893</li><li>901</li><li>903</li><li>906</li><li>907</li><li>913</li><li>918</li><li>920</li><li>921</li><li>922</li><li>925</li><li>928</li><li>935</li><li>943</li><li>944</li><li>945</li><li>946</li><li>947</li><li>948</li><li>949</li><li>950</li><li>951</li><li>952</li><li>953</li><li>954</li><li>955</li><li>956</li><li>957</li><li>958</li><li>959</li><li>960</li><li>961</li><li>965</li><li>967</li><li>968</li><li>971</li><li>972</li><li>973</li><li>976</li><li>978</li><li>979</li><li>980</li><li>981</li><li>982</li><li>983</li><li>984</li><li>985</li><li>986</li><li>987</li><li>988</li><li>989</li><li>993</li><li>995</li><li>997</li><li>999</li><li>1000</li><li>1001</li><li>1003</li><li>1006</li><li>1007</li><li>1008</li><li>1009</li><li>1010</li><li>1011</li><li>1012</li><li>1013</li><li>1014</li><li>1015</li><li>1016</li><li>1017</li><li>1018</li><li>1019</li><li>1020</li><li>1021</li><li>1022</li><li>1023</li><li>1024</li><li>1025</li><li>1026</li><li>1027</li><li>1028</li><li>1029</li><li>1030</li><li>1031</li><li>1032</li><li>1033</li><li>1034</li><li>1035</li><li>1036</li><li>1037</li><li>1038</li><li>1039</li><li>1040</li><li>1041</li><li>1042</li><li>1043</li><li>1044</li><li>1045</li><li>1046</li><li>1047</li><li>1048</li><li>1049</li><li>1050</li><li>1051</li><li>1052</li><li>1053</li><li>1054</li><li>1055</li><li>1056</li><li>1057</li></ol>



### Get the ratios

- Returns a list of vectors containing the ratios for each active (unpruned) particle


```R
ratios <- ratios(results)
ratios[[44]]
```


1


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
<ol class=list-inline><li>0.133126152818256</li><li>1.17098996548298e-05</li><li>7.14722350707829e-06</li><li>6.10367023166043e-06</li><li>8.12485280405636e-06</li><li>6.38782844896821e-06</li><li>1.27204834156928e-06</li><li>1.43092815226752e-06</li><li>2.05188208380923e-06</li><li>1.4080708137156e-06</li><li>1.83341750550405e-07</li><li>1.58609408727374e-07</li><li>2.58412814068089e-07</li><li>1.76717119785081e-07</li><li>1.95777028855652e-07</li><li>2.84472027388314e-07</li><li>1.28163398486342e-07</li><li>1.56778256370441e-07</li><li>1.40948785835059e-07</li><li>1.90972454419062e-07</li><li>1.92270389678651e-07</li><li>1.2782675447717e-07</li><li>1.2612619501017e-07</li><li>1.97929507424119e-07</li><li>1.47001109267379e-07</li><li>1.96483075367692e-07</li><li>2.88892444060236e-07</li><li>4.36483523232203e-07</li><li>4.83349306154401e-07</li><li>3.21873846099633e-07</li><li>2.93366049136795e-07</li><li>3.74309857098558e-07</li><li>2.68853254649256e-07</li><li>3.02624894118601e-07</li><li>2.20089064833941e-07</li><li>2.00549091061836e-07</li><li>2.61381816799501e-07</li><li>2.6149852167238e-07</li><li>2.37092363010098e-07</li><li>1.26262390136046e-07</li><li>1.97672405051065e-07</li><li>1.70057856371003e-07</li><li>2.78650895140994e-07</li><li>2.71828865062546e-07</li><li>1.40322359817444e-07</li><li>1.34773978920527e-07</li><li>1.64695059495912e-07</li><li>1.29516346170582e-07</li><li>2.48776597440567e-07</li><li>2.58847957847807e-07</li><li>5.46094778429603e-07</li><li>2.57467766670971e-07</li><li>2.28949061200577e-07</li><li>2.30544047893865e-07</li><li>2.17984531576782e-07</li><li>2.9131928361262e-07</li><li>1.71676585793455e-07</li><li>1.63654236249923e-07</li><li>1.99675519987472e-07</li><li>2.08244502025018e-07</li><li>2.93735615425763e-07</li><li>1.4377288988522e-07</li><li>1.70223640077992e-07</li><li>1.3926586973224e-07</li><li>1.42446966220921e-07</li><li>1.87205082403744e-07</li><li>1.78749579667229e-07</li><li>2.93676487381189e-07</li><li>1.28718460991831e-07</li><li>2.33417185355242e-07</li><li>3.50128443706869e-07</li><li>3.34979758018475e-07</li><li>1.64263599770506e-07</li><li>1.5563089441412e-07</li><li>2.11875665030156e-07</li><li>1.31497014649579e-07</li><li>1.72583253822788e-07</li><li>2.06723170658404e-07</li><li>2.6334572573849e-07</li><li>2.57450314187056e-07</li><li>2.13099649663135e-07</li><li>3.44740350139873e-07</li><li>1.96428075595439e-07</li><li>1.58018753034826e-07</li><li>1.29747123046555e-07</li><li>3.1905300857783e-07</li><li>2.3511329731933e-07</li><li>2.31296580436941e-07</li><li>2.44938421249196e-07</li><li>2.76233953679588e-07</li><li>1.89260811266333e-07</li><li>2.43048987457595e-07</li><li>1.65354869525355e-07</li><li>1.93117276717968e-07</li><li>1.39956750039373e-07</li><li>1.99931503639076e-07</li><li>3.37005619360818e-07</li><li>1.60823887808981e-07</li><li>3.66496210994054e-07</li><li>1.80058221231906e-07</li><li>1.8495131476063e-07</li><li>2.71888642730603e-07</li><li>1.52166662523776e-07</li><li>2.07406959392694e-07</li><li>2.10669740035206e-07</li><li>4.69760769154173e-07</li><li>6.52182862842249e-07</li><li>3.64402961696128e-07</li><li>2.54186845975598e-07</li><li>3.24943775255666e-07</li><li>1.61544967531062e-07</li><li>2.57003173844353e-07</li><li>2.70214754206207e-07</li><li>2.2637459950401e-07</li><li>4.23920587395987e-07</li><li>2.8322114483801e-07</li><li>3.30428643972373e-07</li><li>2.58750956164194e-07</li><li>2.89086708983101e-07</li><li>2.74049237532618e-07</li><li>8.44101050272715e-07</li><li>6.17598922564642e-07</li><li>1.00242897169449e-06</li><li>8.44610196785116e-07</li><li>1.12541317207676e-06</li><li>4.88268049450523e-07</li><li>1.00531893449994e-06</li><li>8.40607533559899e-07</li><li>1.63074782108577e-06</li><li>8.94299563703863e-07</li><li>8.84252414837957e-07</li><li>9.3854001441844e-07</li><li>2.19120453907265e-06</li><li>1.86573572248706e-06</li><li>1.90202187679516e-06</li><li>2.67795746810446e-06</li><li>1.32726727514234e-06</li><li>4.0301581604301e-06</li><li>3.36780000093992e-06</li><li>4.93482157282705e-06</li><li>5.79970825815224e-06</li><li>5.40837761466287e-06</li><li>4.87466443165536e-06</li><li>1.07564827967827e-05</li><li>5.40358739602619e-06</li><li>9.81077403147276e-06</li><li>1.44943254439947e-05</li><li>1.62045207200543e-05</li><li>8.16318876408154e-06</li><li>7.85447772775378e-06</li><li>1.16785556510933e-05</li><li>1.75280078480629e-05</li><li>2.43680005935932e-05</li><li>9.19739495574062e-05</li><li>3.33569742476221e-05</li><li>4.54525110031176e-05</li><li>7.35382205646925e-05</li><li>0.000131885814844327</li><li>0.000171385577330686</li><li>0.000226277196562355</li><li>0.000129453498068077</li><li>0.000216118181676151</li><li>0.000191606388098246</li><li>0.000217101943042248</li><li>0.000373585300780021</li><li>0.000494169567673899</li><li>0.00119733987989781</li><li>0.0033139493262468</li><li>0.00427212179658708</li><li>0.000373973924854536</li><li>0.000267938787271378</li><li>0.000175114013454459</li><li>0.000282078589407088</li><li>0.000390341780575317</li><li>0.000358311749941584</li><li>0.000726284750826467</li><li>0.000409703025150144</li><li>0.000384691664211269</li><li>0.000680029898370651</li><li>0.000788591976429736</li><li>0.000795715496312821</li><li>0.00127225722460704</li><li>0.00123623170684907</li><li>0.00171334246237871</li><li>0.00273269637725282</li><li>0.00520633280230676</li><li>0.00779309174837626</li><li>0.00235418514548239</li><li>0.00212126007160128</li><li>0.00409656151877303</li><li>0.00496444605817324</li><li>0.00175003648001585</li><li>0.00231285375417786</li><li>0.00329743941200314</li><li>0.00453391332242189</li><li>0.00366698562854509</li><li>0.00507909455883998</li><li>0.0069704835009334</li><li>0.0154430866223695</li><li>0.00973648274345358</li><li>0.0160195266212548</li><li>0.0268143431035235</li><li>0.0206470610759251</li><li>0.0323741228550206</li><li>0.0266874221998823</li><li>0.0122074774116735</li><li>0.0124154816423483</li><li>0.0126939905204126</li><li>0.016763803253236</li><li>0.0155709372533831</li><li>0.0132148822414063</li><li>0.0171418684696961</li><li>0.0179613780443542</li><li>0.0248004120414112</li><li>0.0297523070219905</li><li>0.0290492227458503</li><li>0.0340444824792816</li><li>0.042613679075774</li><li>0.0195723402471913</li><li>0.0184441962501854</li><li>0.0254261868960831</li><li>0.0306147535432261</li><li>0.0466201367632341</li><li>0.0862536690100928</li><li>0.0972728569775062</li><li>0.00618163440846622</li><li>0.00738039968232687</li><li>0.00822100219932533</li><li>0.00475875284603136</li><li>0.00198249221600879</li><li>0.000719376819569842</li><li>0.000704629597756609</li><li>0.000729878198116092</li><li>0.000899658486897449</li><li>0.000954793546788097</li><li>0.00116067085296595</li><li>0.000599200338686834</li><li>0.000619877217604593</li><li>0.000649050790482323</li><li>0.00056515871517072</li><li>0.000681949764943114</li><li>0.00078858115352492</li><li>0.000235048947111482</li><li>0.000228436252350799</li><li>0.000205584020046313</li><li>0.00020456548300303</li><li>1.40224142998467e-05</li><li>2.17269883052614e-05</li><li>1.99362565749796e-05</li><li>3.06153139121448e-05</li></ol>


