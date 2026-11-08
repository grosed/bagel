# bagelR

Fast Bayesian online changepoint detection for linear models, with an R interface and a C++ backend. The R package source is in `bagelR/`; the standalone C++ implementation is in `c++/`.

## Installation

Install the development package from GitHub:

```r
install.packages("remotes")
remotes::install_github("grosed/bagel", subdir = "bagelR")
```

Installation from source requires a C++20 toolchain. To build the introductory vignette as well:

```r
install.packages(c("knitr", "rmarkdown"))
remotes::install_github("grosed/bagel", subdir = "bagelR",
                        build_vignettes = TRUE)
vignette("bagelR", package = "bagelR")
```

Vignette building also requires Pandoc. The tutorial source is [Getting started with bagelR](bagelR/vignettes/bagelR.Rmd). It covers mean, slope and seasonal changes, autocorrelated errors, and sequential detection.

```r
library(bagelR)
packageVersion("bagelR")
```

See [Changes](#changes) for release notes.

The introduction video is distributed separately from the R package. The [video upload instructions](release-assets/bagel-introduction/README.md) describe how to publish the prepared MP4 at the address used by the vignette.

## Building and checking

From this repository directory, run:

```sh
R CMD build bagelR
R CMD check --as-cran bagelR_1.8.2.tar.gz
```

Use the version from `bagelR/DESCRIPTION` if it changes. Build and check outputs, compiled files and R session histories are excluded from version control. The package vignette is rebuilt from its R Markdown source.
