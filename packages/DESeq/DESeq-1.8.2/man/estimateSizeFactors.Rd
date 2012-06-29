\name{estimateSizeFactors}
\Rdversion{1.1}
\docType{methods}
\alias{estimateSizeFactors,CountDataSet-method}
\alias{estimateSizeFactors}

\title{
   Estimate the size factors for a CountDataSet
}
\description{
   Estimate the size factors for a CountDataSet
}
\usage{
\S4method{estimateSizeFactors}{CountDataSet}( object, locfunc=median, ... )
}
\arguments{
  \item{object}{
     a CountDataSet
}
  \item{locfunc}{
      a function to compute a location for a sample. By default, the median is used. However,
      especially for low counts, the \code{\link[genefilter]{shorth}} may give better results.
}   
   \item{...}{extra arguments are ignored}
}
\details{
   You need to call this function right after \code{\link{newCountDataSet}} unless
   you have manually specified size factors.
   
   Typically, the function is called with the idiom
   
   \code{ cds <- estimateSizeFactors( cds ) }
   
   This estimates the size factors and stores the information in the object.

   Internally, the function calls \code{\link{estimateSizeFactorsForMatrix}}. See there
   for more details on the calculation.
}
\value{
   The CountDataSet passed as parameters, with the size factors filled in.
}

\author{
   Simon Anders, sanders@fs.tum.de
}
\seealso{
\code{\link{estimateSizeFactorsForMatrix}}
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
}
