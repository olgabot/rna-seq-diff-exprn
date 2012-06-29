\name{nbinomGLMTest}
\alias{nbinomGLMTest}
\title{
   Perform chi-squared tests comparing two sets of GLM fits
}
\description{
   For each gene, the function calculates a chi-square p value by simply calculating:
   \code{1 - pchisq(resReduced$deviance - resFull$deviance, attr(resReduced, 
    "df.residual") - attr(resFull, "df.residual"))}
}
\usage{
nbinomGLMTest(resFull, resReduced)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resFull, resReduced}{
     GLM fit data frames, as returned by \code{\link{fitNbinomGLMs}},
     first the full, then the reduced model. }
}
\value{
  a vector of p values
}
\author{
Simon Anders, anders@embl.de
}
\seealso{
\code{\link{fitNbinomGLMs}}
}
\examples{
cds <- makeExampleCountDataSet()[ 1:100, ]
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="pooled" )
fit1 <- fitNbinomGLMs( cds, count ~ condition )
fit0 <- fitNbinomGLMs( cds, count ~ 1 )
nbinomGLMTest( fit1, fit0 )
}

