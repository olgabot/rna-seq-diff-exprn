\name{sizeFactors}
\Rdversion{1.1}
\docType{methods}
\alias{sizeFactors,CountDataSet-method}
\alias{sizeFactors<-,CountDataSet,numeric-method}
\alias{sizeFactors}
\title{
    Accessor functions for the 'sizeFactors' information in a
    CountDataSet object.
}
\description{
    The sizeFactors vector assigns to each column of the count data
    a value, the size factor, such that count values in the columns can be brought to
    a common scale by dividing by the corresponding size factor.
}
\usage{
\S4method{sizeFactors}{CountDataSet}(object)
\S4method{sizeFactors}{CountDataSet,numeric}(object) <- value
}
\arguments{
  \item{object}{a \code{CountDataSet} object.}
  \item{value}{a numeric vector, one size factor for each column in the count data.}
}

\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
   cds <- makeExampleCountDataSet()
   cds <- estimateSizeFactors( cds )
   sizeFactors(cds)
}
\author{
   Simon Anders, sanders@fs.tum.de
}

