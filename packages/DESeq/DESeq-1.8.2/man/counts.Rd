\name{counts}
\Rdversion{1.1}
\docType{methods}
\alias{counts,CountDataSet-method}
\alias{counts<-,CountDataSet,matrix-method}
\alias{counts}
\title{
   Accessors for the 'counts' slot of a CountDataSet object.
}

\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, one row for each observational unit (gene or the like), and
   one column for each sample.
}

\usage{
\S4method{counts}{CountDataSet}(object, normalized=FALSE)
\S4method{counts}{CountDataSet,matrix}(object) <- value
}

\arguments{
  \item{object}{a \code{CountDataSet} object.}
  \item{normalized}{logical indicating whether or not to divide the
    counts by the size factors before returning.}
  \item{value}{integer matrix.}
}

\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   cds <- makeExampleCountDataSet()
   head( counts( cds ) )
}   
   
