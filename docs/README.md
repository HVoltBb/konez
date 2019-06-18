# konez
An R package for k-aggregated count data regression with covariates and random effects

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

The functionality of this package depends on JAGS 4.X.X which can be downloaded at <http://mcmc-jags.sourceforge.net>.

## Build-in datasets
There are four datasets included in this package, and they are seabird bycatch (partial), Legionnaire's disease incidence in Singapore, Leadbeater's possum abundance and frigatebird mesting sites. These datasets were used as case studies in a manuscript currently under review in Ecological Modelling. These datasets were either extracted from published literature or downloaded directly from the publisher's data repository. Seabird bycatch data is only a partial summary, only enough to generate a histogram to view the distribution of counts. To replicate the results presented in that manuscript, however, you need to acquire the whole dataset from the original authors of the publication or from relevant authorities. The other datasets are the same as the ones used in the manuscript. Relevant references are given below. 

* Xu, H.-Y., Xie, M., Goh, T.N., 2014. Objective Bayes analysis of zero-inflated Poisson distribution with application to healthcare data. IIE Transactions 46, 843-852.
* Cunningham, R.B., Lindenmayer, D.B., 2005. Modeling count data of rare species: some statistical issues. Ecology 86, 1135-1142.
* Lindenmayer, D., Nix, H., McMahon, J., Hutchinson, M., Tanton, M., 1991. The conservation of Leadbeater's possum, Gymnobelideus leadbeateri (McCoy): a case study of the use of bioclimatic modelling. Journal of Biogeography, 371-383.
* Li, Y., Jiao, Y., Browder, J.A., 2016. Assessment of seabird bycatch in the US Atlantic pelagic longline fishery, with an extra exploration on modeling spatial variation. ICES Journal of Marine Science 73, 2687-2694.

I __do not__ either own or maintain any of the datasets mentioned above. Any data related questions should be directed toward the respective authors or relevant authorities. 




If there is any questions to me, contact at <eiDOTog@gmail.com>.

Thanks,

Can Zhou

5/21/2019
