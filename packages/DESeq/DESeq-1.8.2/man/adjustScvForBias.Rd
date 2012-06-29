\name{adjustScvForBias}
\Rdversion{1.1}
\alias{adjustScvForBias}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Adjust an SCV value for the bias arising when it is calculated from
   unbiased estimates of mean and variance.
}
\description{
   Assume that a small sample of i.i.d. random variables from a negative binomial
   distribution is given, and you have obtained unbiased estimates of mean and
   raw variance. Then, a new bias is introduced when the squared coefficient of 
   variation (SCV, a.k.a. dispersion) is calculated from these unbiased estimates by dividing the
   raw variance by the square of the mean. This bias can be calculated by numerical
   simulation and a pre-calculated adjustment table (or rather a fit through tabulated
   values) is supplied with the package. The present function uses this to remove the 
   bias from a raw SCV estimate.
   
   This function is used internally in \code{\link{nbinomTest}}. You will rarely need 
   to call it directly.
}
\usage{
adjustScvForBias(scv, nsamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scv}{
     An estimate for the raw squared coefficient of variation (SCV) for negative
     binomially distributed data, which has been obtained by dividing an unbiased
     estimate of the raw variance by the square of an unbiased estimate of the mean.
}
  \item{nsamples}{
     The size of the sample used in the estimation.
}
}
\value{
   an unbiased estimate of the raw SCV
}
\author{
   Simon Anders
}
\examples{
   true_mean <- 100
   true_scv <- .1
   nsamples <- 3
   res <- replicate( 1000, {
      mySample <- rnbinom( nsamples, mu=true_mean, size=1/true_scv )
      mu_est <- mean( mySample )
      raw_var_est <- var( mySample ) - mean( mySample )
      raw_scv_est <- raw_var_est / mu_est^2
      unbiased_raw_scv_est <- adjustScvForBias( raw_scv_est, 4 )
      c( raw_scv_est = raw_scv_est, unbiased_raw_scv_est = unbiased_raw_scv_est ) } )
   rowMeans( res )
}
