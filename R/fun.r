#' Find k Based on DIC
#'
#' @description Fit a series of models with different ks to data, and select the k value based on DIC.
#' @param model character string, partially matched to c('poisson', 'negbinom', 'cmp', 'Tpoisson','Tnegbinom', 'Tcmp'): see 'Details'.
#' @param data numeric vector of counts.
#' @param ks numeric vector of integers.
#' @param type numeric value, specifying the type of regression model: 0: no covariates; 1: with covariates; 2: with covariates and random effects. Defaults to \code{0}.
#' @param model_par list of model parameters. \code{Xc} is the covariate dataframe for the count model, \code{Xrc} is the covariate dataframe of the random effect for the count model, \code{Xz} is the covariate dataframe for the zero/one model, \code{Xrz} is the covariate dataframe of the random effect for the zero/one model, \code{size_upper} is the upper bound of the prior of the size parameter of a negative binomial distribution (defaults to 100) and \code{maxiter} is a positive integer, specifying the number of positive term to keep in the calculation of the Conway-Maxwell-Poisson distribution (defaults to 50).
#' @param jags_par list of variables to pass to run.jags function.
#' @examples x = find_k('pois', legion) # Find k-aggregated zero-inflated Poisson model to the Legionnaires data
#' @export

find_k <- function(model=c('poisson', 'negbinom', 'cmp', 'Tpoisson','Tnegbinom', 'Tcmp'), count, ks = 0:5, type = 0, model_par = list(Xc = NA, Xz = NA, Xrc = NA, Xrz = NA, maxiter=50, size_upper=100), jags_par=list(chain = 2, sample = 1, thin = 5, method = 'rjparallel', burnin = 500, inits = inix, dic.sample = 1e3)){
  model = parse_model(model, type)
  cat('model:', model, '\n')
  listx = parse_data(count, type, model_par)
  
  runjags::runjags.options('silent.runjags'=TRUE, 'silent.jags'=TRUE)
  rjags::load.module('lecuyer')
  rjags::parallel.seeds('lecuyer::RngStream', jags_par$chain)
  dics = rep(NA, length(ks))
  names(dics) <- paste0('k=',ks)
  for (i in 1:length(ks)){
    listx$datalist$k = ks[i]
    cat("k=", ks[i], '>> ')
    fitj = do_fit(model, listx, jags_par)
    temp = runjags::extract(fitj, 'dic', n.iter = jags_par$dic.sample)
    dics[i] = sum(temp[[1]]+temp[[2]], na.rm = TRUE)

  }
  cat('Finished k =',paste0(ks,', '),'\b\b\b.\n')
  runjags::runjags.options('silent.runjags'=FALSE, 'silent.jags'=FALSE)
  return(dics)
}

#' Fit a k-aggregated model
#'
#' Fit a k-aggregated model with a fixed k.
#'
#' @param model character string, partially matched to \code{c('poisson', 'negbinom', 'cmp', 'Tpoisson','Tnegbinom', 'Tcmp')}: see 'Details'.
#' @param data numeric vector of counts.
#' @param k numeric (non-negative integer).
#' @param type numeric value, specifying the type of regression model: 0: no covariates; 1: with covariates; 2: with covariates and random effects. Defaults to \code{0}.
#' @param model_par list of model parameters. \code{Xc} is the covariate dataframe for the count model, \code{Xrc} is a dataframe of the covariate of the random effect for the count model, \code{Xz} is the covariate dataframe for the zero/one model, \code{Xrz} is the covariate dataframe of the random effect for the zero/one model, \code{size_upper} is the upper bound of the prior of the size parameter of a negative binomial distribution (defaults to 100) and \code{maxiter} is a positive integer, specifying the number of positive term to keep in the calculation of the Conway-Maxwell-Poisson distribution (defaults to 50).
#' @param jags_par list of variables to pass to run.jags function.
#' @examples fit = fit_k('n', poss$X.Lb, k=1, model_par = list(Xc = log(1+poss['X.Stags']), Xz = log(1+poss['X.Stags'])) # Fit 1-aggregated negative model to the Leadbeater's possum abundance data
#' @export
fit_k <- function(model=c('poisson', 'negbinom', 'cmp', 'Tpoisson','Tnegbinom', 'Tcmp'), count, k = 0, type = 0, model_par = list(Xc = NA, Xz = NA, Xrc = NA, Xrz = NA, maxiter=50, size_upper=100), jags_par=list(chain = 3, sample = 500, thin = 10, method = 'rjparallel', burnin = 1e3, inits = inix)){
  model = parse_model(model, type)
  cat('model:', model, '\n')
  listx = parse_data(count, type, model_par)

  options(warn = -1)
  rjags::load.module('lecuyer')
  rjags::parallel.seeds('lecuyer::RngStream', jags_par$chain)
  listx$datalist$k = k
  fitj = do_fit(model, listx, jags_par)
  options(warn = 0)
  return(fitj)
}

#' Fit a k-aggregated model
do_fit <- function(model, listx, jags_par){
  fitj = runjags::run.jags(model = eval(parse(text = model)),
                  monitor = listx$monitorlist,
                  data = listx$datalist, n.chains = jags_par$chain,
                  inits = jags_par$inits, burnin = jags_par$burnin,
                  sample = jags_par$sample, thin = jags_par$thin,
                  method = jags_par$method, summarise = FALSE)
  return(fitj)
}

#' Parse data for the jags input
parse_data <- function(count, type, model_par){
  dat = c(count)
  inx = order(dat)
  datalist = list(nz = sum(dat==0), n = length(dat), zeros=rep(0, length(dat)), k = 0, BirdNum = dat[inx], maxiter=model_par$maxiter, zero = 0, xc = matrix(0, nrow = length(dat), ncol = 1), xz = matrix(0, nrow = length(dat), ncol = 1), ncov1 = 1, ncov2 = 1, nlevel1 = 1, nlevel2 = 1, rei1=0, rei2=0, xrc = matrix(0, nrow = length(dat), ncol = 1), xrz = matrix(0, nrow = length(dat), ncol = 1), size_upper = model_par$size_upper)
  monitorlist = c('c1', 'c2', 'size', 'logmu', 'lognu')
  if(type >= 1){
    if(is.null(model_par$Xc) || is.na(model_par$Xc)) model_par[['Xc']] = matrix(0, nrow = length(dat), ncol = 1)
    if(is.null(model_par$Xz) || is.na(model_par$Xz)) model_par[['Xz']] = matrix(0, nrow = length(dat), ncol = 1)
    datalist[['xc']] = data.matrix(model_par$Xc[inx,])
    datalist[['xz']] = data.matrix(model_par$Xz[inx,])
    datalist[['ncov1']] = dim(model_par$Xz)[2]
    datalist[['ncov2']] = dim(model_par$Xc)[2]
    monitorlist = c(monitorlist, 'b1', 'b2')
  }
  if(type==2){
    if(is.null(model_par$Xrc) || is.na(model_par$Xrc)) model_par[['Xrc']] = matrix(0, nrow = length(dat), ncol = 1)
    if(is.null(model_par$Xrz) || is.na(model_par$Xrz)) model_par[['Xrz']] = matrix(0, nrow = length(dat), ncol = 1)
    datalist[['xrc']] = data.matrix(model_par$Xrc[inx,])
    datalist[['xrz']] = data.matrix(model_par$Xrz[inx,])
    datalist[['nlevel1']] = dim(model_par$Xrz)[2]
    datalist[['nlevel2']] = dim(model_par$Xrc)[2]
    datalist[['rei1']] = 1
    datalist[['rei2']] = 1
    monitorlist = c(monitorlist, 'r1', 'r2', 'var1', 'var2')
  }
  return(list(datalist= datalist, monitorlist=monitorlist))
}

#' Parse model
parse_model <- function(model, type){
  modellist = c('poissonk', 'negbinomk', 'cmpk', 'Tpoissonk','Tnegbinomk', 'Tcmpk')
  model = modellist[pmatch(model[1], modellist, nomatch = 1)]
  if(type >=1) model = paste0(model, 'x')
  return(model)
}

#' Initial value generating function
#' User should provide their own inits function.

inix = function(chain){
  return(list(c1 = runif(1, -0.6, 0.9),
              c2 = runif(1, -1, 2),
              size = runif(1, 1, 9),
              lognu = rnorm(1, 0, .1),
              var1 = runif(1, 0, 1),
              var2 = runif(1, 0, 1)
  ))
}

dkx_ = function(x, k, family, param=list(lambda=NA, size=NA, logmu=NA,lognu=NA)){
  familylist = c('pois', 'tpois', 'negbinom', 'tnegbinom', 'cmp', 'tcmp')
  family = familylist[pmatch(x = family[1], table = familylist, nomatch = 1)]

  if(family == 'pois'){
    if(!is.na(param$lambda)){
    if(x==0) return(dpois(0, param$lambda))
    else if(x==1) return(ppois(k+1, param$lambda) - dpois(0, param$lambda))
    else {
      return(dpois(x+k, param$lambda))
    }
  }else{
    stop('lambda is not specified\n')
  }}

  if(family == 'tpois'){
    if(!is.na(param$lambda)){
    if(x==0) return(0)
    else if(x==1) return((ppois(k+1, param$lambda) - dpois(0, param$lambda))/(1-dpois(0, param$lambda)))
    else {
      return(dpois(x+k, param$lambda)/(1-dpois(0, param$lambda)))
    }
  }else{
    stop('lambda is not specified\n')
  }}

  if(family == 'negbinom'){
    if(!is.na(param$size) && !is.na(param$logmu)){
    mu = exp(param$logmu)
    if(x==0) return(dnbiom(0, param$size), mu = mu)
    else if(x==1) return(pnbiom(k+1, param$size, mu = mu) - dnbinom(0, param$size, mu = mu))
    else {
      return(dnbinom(x+k, param$size, mu = mu))
    }
  }else {
    stop('Either size or logmu is not specified\n')
  }}

  if(family == 'tnegbinom'){
    if(!is.na(param$size) && !is.na(param$logmu)){
    mu = exp(param$logmu)
    if(x==0) return(0)
    else if(x==1) return((pnbinom(k+1, param$size, mu = mu) - dnbinom(0, param$size, mu = mu))/(1-dnbinom(0, param$size, mu=mu)))
    else {
      return(dnbinom(x+k, param$size, mu = mu)/(1-dnbinom(0, param$size, mu = mu)))
    }
  }else {
    stop('Either size or logmu is not specified\n')
  }}

  if(family == 'cmp'){
    if(!is.na(param$lognu) && !is.na(param$logmu)){
    if(x==0) return(dcmp_(0, param$lognu, param$logmu, maxiter))
    else if(x==1) return(pcmp_(k+1, param$lognu, param$logmu, maxiter) - dcmp_(0, param$lognu, param$logmu, maxiter))
    else {
      return(dcmp_(x+k, param$lognu, param$logmu))
    }
  }else {
    stop('Either lognu or logmu is not specified\n')
  }}

  if(family == 'tcmp'){
    if(!is.na(param$lognu) && !is.na(param$logmu)){
    if(x==0) return(0)
    else if(x==1) return((pcmp_(k+1, param$lognu, param$logmu) - dcmp_(0, param$lognu, param$logmu))/(1-dcmp_(0, param$lognu, param$logmu)))
    else {
      return(dcmp_(x+k, param$lognu, param$logmu)/(1-dcmp_(0, param$lognu, param$logmu)))
    }
  }else {
    stop('Either lognu or logmu is not specified\n')
  }}

}

dcmp_ = function(x, lognu, logmu, maxiter = 50){
  logres = exp(lognu)*(x*logmu - lfactorial(x))
  s = 1
  for(i in 1:maxiter){
    s = s + exp(exp(lognu)*(i*logmu - lfactorial(i)))
  }
  return(exp(logres - log(s)))
}

#' K-aggregated Distribution
#'
#'Density, distribution function, and random generation for the k-aggregated discrete distributions
#'@param x vector of (non-negative integers) quantiles.
#'@param q vector of quantiles.
#'@param n number of random number samples to return.
#'@param k the additional number of mass to aggregate into the single mass.
#'@param family character string partial matched from a list of base line distributions (Poisson, zero truncated Poisson, negative binomial, zero truncated negative binomial, Conway-Maxwell-Poisson, and zero truncated CMP distributions) \code{c('pois', 'tpois', 'negbinom', 'tnegbinom', 'cmp', 'tcmp')}. Defaults to 'pois'.
#'@param param list of parameters to pass to the base line function.
#'@details For different base line distributions, a different set of parameters needs to be supplied through \code{param}: for families 'pois' and 'tpois', \code{param$lambda} cannot be \code{NA}; for families 'negbinom' and 'tnegbinom', \code{param$size} and \code{param$logmu} cannot be \code{NA}; for families 'cmp' and 'tcmp', \code{param$lognu} and \code{param$logmu} cannot be \code{NA}.
#'@seealso dcmp, dpois, dNegBinom
#'@describeIn dkx Density function for a k-aggregated distribution.
#'@export
dkx = Vectorize(dkx_, vectorize.args = 'x')

pkx_ = function(q, k, family, param=list(lambda=NA, size=NA, logmu=NA,lognu=NA)){
  return(sum(dkx(0:q, k, family, param)))
}

#'@describeIn dkx Distribution function for the k-aggregated distribution.
#'@export
pkx = Vectorize(pkx_, vectorize.args = 'q')

#'@describeIn dkx Random number generator for the k-aggregated distribution.
#'@export
rkx = function(n, k, family, param=list(lambda=NA, size=NA, logmu=NA,lognu=NA)){
  familylist = c('pois', 'tpois', 'negbinom', 'tnegbinom', 'cmp', 'tcmp')
  family = familylist[pmatch(family[1], familylist, nomatch=1)]
  res = rep(0, n)
  if(family == 'pois' && !is.na(param$lambda)){
    res = rpois(n, param$lambda)
  }

  if(family == 'tpois' && !is.na(param$lambda)){
    res = rtpois(rep(log(param$lambda), n))
  }

  if(family == 'negbinom' && !is.na(param$size) && !is.na(param$logmu)){
    res = rnbinom(n, param$size, mu = exp(param$logmu))
  }

  if(family == 'tnegbinom' && !is.na(param$size) && !is.na(param$logmu)){
    res = rtnb(n, param$size, param$logmu)
  }

  if(family == 'cmp' && !is.na(param$lognu) && !is.na(param$logmu)){
    res = rcmp(logmu = param$logmu, lognu = param$lognu)
  }

  if(family == 'tcmp' && !is.na(param$lognu) && !is.na(param$logmu)){
    res = rcmpk(exp(logmu), exp(lognu), zip = FALSE)
  }

  indx = which(res > 0 & res <= k+1)
  res[indx] <- 1
  indx = which(res > 1)
  res[indx] <- res[indx] - k
  return(res)
}

rtpois <- function(logmu=1:10){
  len = length(logmu)
  T = exp(logmu)
  U = runif(len)
  t = -log(1- U*(1-exp(-T)))
  T1 = T - t

  rzl = rpois(len, T1) + 1
  return(rzl)
}

rtnb <- function(n, size, logmu){
  return(qtnb(runif(n), size, logmu))
}

qtnb <- function(p, size, logmu){
  p0 = dnbinom(0, size, mu = exp(logmu))
  logp = log(p)
  logp = logp + log(1-p0)
  pnew = exp(logp) + p0
  res = qnbinom(pnew, size, mu = exp(logmu))
  res[p<dkx_(x=1, k = 0, family = 'tn', param = list(size=size, logmu=logmu))] <- 1
  return(res)
}



#' Conway-Maxwell-Poisson Distribution
#'
#' Density, distribution function and random generaton for the Conway-Maxwell-Poisson distribution with parameters \code{logmu} and \code{lognu}.
#'@param x vector of (non-negative integer) quantiles.
#'@param q vector of quantiles.
#'@param logmu the mean using Guikema and Goffelt's (2008) parametrization.
#'@param lognu the shape parameter using Guikema and Goffelt's (2008) parametrization.
#'@param maxiter integer specifying number of term to keep in calculating the normalizing constant. Defaults to 50.
#'@references [1] Guikema, S. D., & Goffelt, J. P. (2008). A flexible count data regression model for risk analysis. Risk Analysis: An International Journal, 28(1), 213-223.
#'@describeIn dcmp Density function for the CMP distribution.
#'@export
dcmp = Vectorize(dcmp_, vectorize.args = 'x')

pcmp_ = function(q, lognu, logmu, maxiter = 50){
  return(sum(dcmp(0:q, lognu, logmu, maxiter)))
}
#'@describeIn dcmp Distribution function for the CMP distribution.
#'@export
pcmp = Vectorize(pcmp_, vectorize.args = 'q')

#'@describeIn dcmp Generates random numbers from the CMP distribution. This is a wrapper for the k-aggregated CMP random number generator.
#'@export
rcmp = function(logmu, lognu, maxiter=50){
  return(rcmpk(exp(logmu), exp(lognu), k=0, maxiter, zip = TRUE))
}

rcmpk <- function(mu=1:10, nu=1, k=0, maxiter=50, zip=F){
  len = length(mu)
  U = runif(len)
  rzl = rep(0, len)
  psi = matrix(, len, maxiter)

  for(j in 1:len){
    for(i in 1:maxiter){
      psi[j,i] <- mu[j]^((i)*nu)/factorial(i)^nu
    }
  }
  #  psi[psi <= prec]  <- 0
  if(zip){
    psi = cbind(rep(1,len),psi)
    psi = apply(psi, 1, function(x){x/sum(x)})
    cdf = apply(rbind(rep(0,len),psi), 2, cumsum)

    for(i in 1:(maxiter + 1)){
      rzl[which(U >= cdf[i,] & U <= cdf[i+1,])] <- i - 1
    }
    indx = which(rzl <= k+1 & rzl > 0)
    rzl[indx] <- 1
    rzl[-indx] <- rzl[-indx] - k
    rzl[rzl<0] <- 0
  } else {
    psi = apply(psi, 1, function(x){x/sum(x)})

    cdf = apply(rbind(rep(0,len),psi), 2, cumsum)

    for(i in 1:maxiter){
      rzl[which(U >= cdf[i,] & U <= cdf[i+1,])] <- i
    }
    indx = which(rzl <= k+1)
    rzl[indx] <- 1
    rzl[-indx] <- rzl[-indx] - k
  }
  return(rzl)
}
