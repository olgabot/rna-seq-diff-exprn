\name{sizeFactors}
\docType{methods}
\alias{sizeFactors,ExonCountSet-method}
\alias{sizeFactors<-,ExonCountSet,numeric-method}
\alias{sizeFactors}
\title{
    Accessor functions for the sizeFactors information in a ExonCountSet
}
\description{
    The sizeFactors vector assigns to each column of the count data
    a value, the size factor, such that count values in the columns can be brought to
    a common scale by dividing by the corresponding size factor.
}
\usage{
\S4method{sizeFactors}{ExonCountSet}(object)
\S4method{sizeFactors}{ExonCountSet,numeric}(object) <- value
}
\arguments{
  \item{object}{
      An ExonCountSet
}
  \item{value}{
      a vector of number, one size factor for each column in the count data
}
}

\seealso{
\code{\link{estimateSizeFactors}}
}
\examples{
   data("pasillaExons", package="pasilla")
   pasillaExons <- estimateSizeFactors( pasillaExons )
   sizeFactors(pasillaExons)
}
\author{
   Simon Anders, sanders@fs.tum.de
}

