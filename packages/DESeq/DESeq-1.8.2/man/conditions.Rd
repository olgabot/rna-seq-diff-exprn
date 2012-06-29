\name{conditions}
\Rdversion{1.1}
\docType{methods}
\alias{conditions,CountDataSet-method}
\alias{conditions<-,CountDataSet-method}
\alias{conditions}
\title{
    Accessor functions for the 'conditions' information in a
    CountDataSet object.
}
\description{
    The conditions vector is a factor that assigns to each column of the count data
    a condition (or treatment, or phenotype, or the like). This information is stored
    in the CountDataSet's "phenoData" slot as a row named "condition".
}
\usage{
\S4method{conditions}{CountDataSet}(object, ...)
\S4method{conditions}{CountDataSet}(object) <- value

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a CountDataSet}
  \item{value}{a vector of suitable length, i.e. with as many elements as
    \code{object} has samples.}
  \item{...}{should not be used for this method.}
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   cds <- makeExampleCountDataSet()
   conditions( cds )
}   
   
