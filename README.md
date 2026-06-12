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
    


    * checking for file ‘/tmp/RtmpAWlt4x/remotes145e511bdb2c82/grosed-bagel-62aaf57/bagelR/DESCRIPTION’ ... OK
    * preparing ‘bagelR’:
    * checking DESCRIPTION meta-information ... OK
    * cleaning src
    * installing the package to process help pages
    * saving partial Rd database
    * cleaning src
    * checking for LF line-endings in source and make files and shell scripts
    * checking for empty or unneeded directories
    * building ‘bagelR_1.5.0.tar.gz’


    Installing package into ‘/home/grosed/DASS/R-packages’
    (as ‘lib’ is unspecified)
    


## Check package version


```R
packageVersion("bagelR")
```

[1] ‘1.5.0’


## Using bagelR

See [here](https://github.com/grosed/bagel/blob/development/notebooks/examples.ipynb) for details of and examples showing how to use **bagelR**.

