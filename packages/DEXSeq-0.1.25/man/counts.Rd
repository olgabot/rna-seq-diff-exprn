\name{counts}
\docType{methods}
\alias{counts,ExonCountSet-method}
\alias{counts<-,ExonCountSet,matrix-method}
\alias{counts}
\title{
   Accessors for the 'counts' slot of a ExonCountSet object.
}
\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, one row for each observational unit (a counting bin, i.e., an exon or part of an exon), and
   one column for each sample.
}
\usage{
\S4method{counts}{ExonCountSet}(object, normalized=FALSE)
\S4method{counts}{ExonCountSet,matrix}(object) <- value
}
\arguments{
  \item{object}{An ExonCountSet object.}
  \item{normalized}{If TRUE, the counts will be normalized by the size factors.}
  \item{value}{An integer matrix of counts, each row corresponding to an exon and each column corresponding to a sample.}
}
\examples{
   data("pasillaExons", package="pasilla")
   head( counts( pasillaExons ) )
}   
