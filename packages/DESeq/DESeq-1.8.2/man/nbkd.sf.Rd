\name{nbkd.sf}
\alias{nbkd.sf}
\title{
GLM family for a negative binomial with known dispersion and log link with size factors
}
\description{
A distribution family for use with \code{\link{glm}}. It describes a negative binomial
(as \code{negative.binomial} in the MASS package), but with a special link function, namely
eta[i] = log( mu[i] / sf[i] ), i.e., each count value is divided by its size factor before
the log is taken. This is used internally by \code{\link{fitNbinomGLMs}}.
}
\usage{
nbkd.sf(r, sf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
      The 'size' argument (see \code{\link{dnbinom}}), i.e., the reciprocal of
      the dispersion.
}
  \item{sf}{
      A vector of size factors.
}
}
\value{
A GLM family object.
}
\author{
Simon Anders, anders@embl.de
}

