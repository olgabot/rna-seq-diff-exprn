\name{fitNbinomGLMs}
\alias{fitNbinomGLMs}

\title{ Fit a generalized linear model (GLM) for each gene. }
\description{
  Use this function to estimate coefficients and calculate deviance from a GLM
  for each gene. The GLM uses the \code{\link{nbkd.sf}} family, with the dispersion
  estimate according to \code{getVarianceFunction(cds)}. Note that this requires
  that the variance functions were estimated with method "pooled" or "blind".
}
\usage{
fitNbinomGLMs( cds, modelFormula, glmControl=list() )
}
\arguments{
  \item{cds}{ a CountDataSet }
  \item{modelFormula}{ a formula. The left hand side must be 'count' (not 'counts'!), the right
     hand side can involve any column of \code{pData(cds)}, i.e., \code{pData(cds)}
     is used as the model frame. If you have passed just a single factor to the
     'conditions' argument of \code{\link{newCountDataSet}}, it can be referred to as
     'condition' in the formula. If you have passed a data frame to 'conditions',
     all columns of this data frame will be available. }
  \item{glmControl}{ list of additional parameters to be passed to \code{\link[stats]{glm.control}} }
}

\value{
A data frame with one row for each gene and columns as follows:
\itemize{
\item{ one column for each estimated coefficient, on a log2 scale (i.e., the natural
log reported by \code{\link{glm}} is rescaled to base 2) }
\item{ a column 'deviance', with the deviance of the fit }
\item{ a boolean column 'converged', indicating whether the fit converged }
}
Furthermore, the data frame has a scalar attribute 'df.residual' that contains
the number of residual degrees of freedom.
}
\author{
Simon Anders (sanders@fs.tum.de)
}
\seealso{
\code{\link{newCountDataSet}},\code{\link{nbinomGLMTest}}, \code{\link{nbkd.sf}}
}
\examples{
   # see nbinomGLMTest for an example
}   
