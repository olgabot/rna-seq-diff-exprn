\name{testGeneForDEU}
\alias{testGeneForDEU}
\title{
      Test a single gene for differential exon usage.
}
\description{
    This function first fits a GLM for the null model, then a GLM for
    the full model for each exon counting bin. Then, p values are derived with
    a chi-squared test from the deviance differences between the models.
}
\usage{
testGeneForDEU( ecs, gene, formula0=NULL, formula1=NULL )
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{gene}{
     The ID of the gene to be tested for differential exon usage.
  }
  \item{formula0}{
     Formula for the null model. If NULL, the default 
     "\code{count ~ sample + exon + condition} is used.
  }
  \item{formula1}{
     Formula for the full model. If NULL, the default 
     "\code{count ~ sample + exon + condition * I(exon==exonID)} is used.
  }
}
\value{
	A data frame with columns "deviance", "df" (degrees of freedom) and pvalues from the test.
}
\details{
   The terms in the formulas must be columns of \code{design(ecs)}. In addition,
   in \code{formula1}, the variable \code{exonID} is set to the
   ID of the currently tested exon counting bin, looping through all the counting bins.
   
   The GLMs are of the negative binomial family, using the dispersions from the \code{dispersion}
   column in \code{fData(ecs)}.
   
}
\examples{
   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   pasillaExons <- estimateDispersions( pasillaExons )
   pasillaExons <- fitDispersionFunction( pasillaExons )
   testGeneForDEU(pasillaExons, "FBgn0085442")
}
\seealso{
   \code{testForDEU}
}
