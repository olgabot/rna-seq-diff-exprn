\name{fitDispersionFunction}
\alias{fitDispersionFunction}
\title{Fit the mean-variance function.}
\description{
      This function fits a parametric model of the mean-dispersion
      relationship to the per-gene estimates of mean \eqn{\hat{\mu}} and
      dispersion \eqn{\hat{\alpha}}. The parametric
      model is \deqn{\alpha(\mu) = \frac{\alpha_1}{\mu} + \alpha_0,}
      where \eqn{\mu} is the mean, \eqn{\alpha} the dispersion and 
      \eqn{\alpha_1} and \eqn{\alpha_0} are two parameters.
      After this, for each exon, the maximum
      between the per-gene estimate \eqn{\hat{\alpha}} and the modelled
      value \eqn{\hat{\alpha}_1/\hat{\mu} + \hat{\alpha}_0} is
	stored in \code{fData$dispersion}.
}
\usage{
fitDispersionFunction(ecs)
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
}
\value{
	An ExonCountSet object with information of the fit included, as well as \code{fData(ecs)$dispersion} filled.
}
\examples{
	data("pasillaExons", package="pasilla")
	pasillaExons <- estimateSizeFactors( pasillaExons )
        pasillaExons <- estimateDispersions( pasillaExons )
        pasillaExons <- fitDispersionFunction( pasillaExons )
}
