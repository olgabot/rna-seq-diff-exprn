\name{fitNbinomGLMsForMatrix}
\alias{fitNbinomGLMsForMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit negative binomial GLMs to a count matrix.
}
\description{
This is a low-level function that is wrapped by \code{\link{nbinomGLMTest}}.
}
\usage{
fitNbinomGLMsForMatrix(counts, sizeFactors, rawScv, modelFormula, 
   modelFrame, quiet = FALSE, reportLog2 = TRUE, glmControl = list() )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
     a matrix of integer counts. Rows for genes, Columns for samples.
}
  \item{sizeFactors}{
     a vector with a size factor for each column in 'counts'.
}
  \item{rawScv}{
     a vector with a raw SCV (i.e., a dispersion) for each row in 'counts'.
}
  \item{modelFormula}{
     a model formula. The left hand side should read 'count ~'.
}
  \item{modelFrame}{
     a model frame (with one row for each column in 'counts')
}
  \item{quiet}{
     whether to not print dots
}
  \item{reportLog2}{
     whether to convert reported coefficients from natural log to log2 scale
}
  \item{glmControl}{ 
     list of additional parameters to be passed to \code{\link[stats]{glm.control}} 
}
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
Simon Anders, sanders@fs.tum.de
}
\examples{
# See the code of fitNbinomGLMs for an example.
}
