\name{getBaseMeansAndVariances}
\Rdversion{1.1}
\alias{getBaseMeansAndVariances}
\title{
    Perform row-wise estimates of base-level means and variances for count data.
}
\description{
    This function is called internally by a number of other functions. You will need
    to call it directly only in very special cases.
}
\usage{
getBaseMeansAndVariances(counts, sizeFactors)
}
\arguments{
  \item{counts}{
     a matrix of data frame of count data. All the columns of this matrix will be
     considered as replicates of the same condition.
}
  \item{sizeFactors}{
     the size factors of the columns, as estimated e.g. with \code{\link{estimateSizeFactorsForMatrix}}
}
}
\value{
A data frame with one row for each row in 'counts' and two columns:
  \item{baseMean }{The base mean for each row. This is the mean of the counts after they have been divided
     by the size factors}
  \item{comp2 }{The base variance for each row. This is the variance of the counts after they have been divided
     by the size factors}
}
\details{
This function is kept for backwards compatibility. See the example below for an alternative and
more self-explanatory way to get the same data.
}

\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
head( getBaseMeansAndVariances( counts(cds), sizeFactors(cds) ) )

# You can get the same as follows
head( rowMeans( counts( cds, normalized=TRUE ) ) )
head( genefilter::rowVars( counts( cds, normalized=TRUE ) ) )
}
