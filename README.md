# konez
K-aggregated count data regression with covariates and random effects

## How to install
First, you will need R package _devtools_. Run the following command in R, if you don't have devtools on your computer.
```R
install.packages('devtools') # You don't need this line if you already have it on you computer
library('devtools')
```
Either one of the following method can be used to install and load konez package:

* Download the zip file of konez package from GitHub, unzip it to your working directory on your computer, and run the following commands in R,
```R
devtools::install('konez')
library(konez)
```
* For this method, you don't need to download anything manually, but a live internet connection at the time of installation. Just run the following two lines of code and you are all set,
```R
devtools::install_github("HVoltBb/konez")
library(konez)
```

The functionality of this package depends on JAGS 4.X.X, and r packages: rjags, and runjags. JAGS software can be downloaded at http://mcmc-jags.sourceforge.net/, and rjags and runjags can be installed normally in R through CRAN.


If there is any questions, let me know at <eidotog@gmail.com>.

Thanks,

Can Zhou
