\name{design}
\docType{methods}
\alias{design,ExonCountSet-method}
\alias{design<-,ExonCountSet-method}
\alias{design}
\title{
    Accessor function for the design annotation from a ExonCountSet object.
}
\description{
    The design vector is a factor or data frame that assigns to each column of the count data
    a condition (or treatment, or phenotype, or the like). This information is stored
    in the ExonCountSet's "phenoData" slot as a row.
}
\usage{
\S4method{design}{ExonCountSet}(object, drop=TRUE, asAnnotatedDataFrame=FALSE)
\S4method{design}{ExonCountSet}(object) <- value
}
\arguments{
  \item{object}{
      An ExonCountSet
   }
  \item{drop}{
      Indicates whether to return a single factor instead of a data frame in case of
      a one-way design
   }
  \item{asAnnotatedDataFrame}{
      Indicates whether the result should be presented as an AnnotatedDataFrame.
   }
  \item{value}{
      A vector or matrix with conditions for the samples, one row for each column in the count data.
   }
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   library(DEXSeq)
   data("pasillaExons", package="pasilla")
   design( pasillaExons )
}
