# konez
An R package for k-aggregated count data regression and simulation with covariates and random effects

## What is k-aggregation?
K-aggregation is a modelling technique to accommodate excess ones in count data. Many count data of rare events have an excess amount of ones such that traditional distributions, e.g., Poisson, negative binomial and Conway-Maxwell-Poisson distributions, under-estimate the probability of the one count. K-aggregated transformation is a generalization of the baseline distributions, and as a result, the fit of a selected k-aggregated model can't get any worse than the baseline model. The selected k-aggregated model is at least as good as the baseline model.

If what I just said makes perfect sense to you, then you can skip reading the paper. It is not a difficult paper to read, I certainly hope not. So, __what is a k-aggregated model, anyway?__ 

The long answer is "[K-aggregated transformation of discrete distributions improves modeling count data with excess ones](https://authors.elsevier.com/a/1ZIPN15DJ~xLzr)", and the short answer is the following picture (Figure 1).

![A k-aggregated distribution](https://HVoltBb.github.io/pics/pic.png)

__Figure 1. The probability mass function of a k-aggregated distribution over the positive range.__

In that figure, I plotted the probability mass function of a k-aggregated distribution P(x) only over the positive range. The probability of zero would have a zero-inflated or hurdle structure, which would complicate the picture. So, I left P(0) out, but you can always add it back in. The k-aggregated distribution P(x) is constructed from a baseline count distribution P<sub>b</sub>(x), which has probabilities depicted as those maroon boxes. P(1) is constructed by stacking P<sub>b</sub>(1) through P<sub>b</sub>(k+1), and for i>1, P(i) is P<sub>b</sub>(k+i). Specifically, with k=0, P(x) is equivalent to P<sub>b</sub>(x). The letter "k" has no other significance than an index variable. Any other letter could have been used, but when I originally wrote the code, it was the letter "k" that was used in the loop to sum up those P<sub>b</sub>s.

## Errata 
[K-aggregated transformation of discrete distributions improves modeling count data with excess ones](https://authors.elsevier.com/a/1ZIPN15DJ~xLzr) in Ecological Modelling

1. On page 4, left column, line 3: <img src="https://latex.codecogs.com/gif.latex?\Gamma " /> is the gamma function, not the lower case <img src="https://latex.codecogs.com/gif.latex?\gamma " /> as appeared in the text. Apparently "\Gamma" was auto-corrected to "\gamma" during production. It was even rendered as <img src="https://latex.codecogs.com/gif.latex?\Delta " /> in the final proof, even through the latex is correct.  

2. The first sentence of section 2.3 should be "To model excess ones, we aggregated the initial k+1 probabilities of an original distribution g(x) over the positive range to represent the probability of a singleton outcome, and left shifted remaining probabilities to represent the probability of two or more counts, ". The original sentence wasn't clear. 

## How to install
First, you will need R package _devtools_. Run the following command in R, if you don't have devtools on your computer.
```R
install.packages('devtools') # You don't need this line if you already have it on you computer
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

The functionality of this package depends on JAGS 4.X.X, which can be downloaded at <http://mcmc-jags.sourceforge.net>.

## Build-in datasets
There are four datasets included in this package, and they are seabird bycatch (partial), Legionnaire's disease incidence in Singapore, Leadbeater's possum abundance and frigatebird nesting sites. These datasets were used as case studies in [a paper in Ecological Modelling](https://authors.elsevier.com/a/1ZIPN15DJ~xLzr).

These datasets were either extracted from published literature or downloaded directly from the publisher's data repository. Seabird bycatch data is only a partial summary, only enough to generate a histogram to view the distribution of counts. To replicate the results presented in that manuscript, however, you need to acquire the whole dataset from the original authors of the publication or from relevant authorities. The other datasets are the same as the ones used in the manuscript. Relevant references are given below. 

* Xu, H.-Y., Xie, M., Goh, T.N., 2014. Objective Bayes analysis of zero-inflated Poisson distribution with application to healthcare data. IIE Transactions 46, 843-852.
* Cunningham, R.B., Lindenmayer, D.B., 2005. Modeling count data of rare species: some statistical issues. Ecology 86, 1135-1142.
* Lindenmayer, D., Nix, H., McMahon, J., Hutchinson, M., Tanton, M., 1991. The conservation of Leadbeater's possum, Gymnobelideus leadbeateri (McCoy): a case study of the use of bioclimatic modelling. Journal of Biogeography, 371-383.
* Li, Y., Jiao, Y., Browder, J.A., 2016. Assessment of seabird bycatch in the US Atlantic pelagic longline fishery, with an extra exploration on modeling spatial variation. ICES Journal of Marine Science 73, 2687-2694.

I __do not__ either own or maintain any of the datasets mentioned above. Any data related questions should be directed toward the respective authors or relevant authorities. 

## Tutorial
I plan to add a tutorial on how to use k-aggregated models from start to finish. Sections to include are:
1. Data formating
2. Fitting a family of k-aggregated model
3. Model selection based on DIC
4. Model prediction


If there is any questions to me, contact me at <eiDOTog@gmail.com>.

Thanks,

Can Zhou

7/1/2019 (updated)
