\name{fitInfo}
\Rdversion{1.1}
\alias{fitInfo}
\title{
    Accessor function for the fitInfo objects in a CountDataSet
}
\description{
   After calling \code{\link{estimateDispersions}}, a CountDataSet object
   is populated with one or (in case of a ``per-condition'' estimation)
   several fitInfo objects, which can be accessed with this function.
}
\usage{
fitInfo( cds, name=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cds}{
      a CountDataSet
}
  \item{name}{
      if \code{estimateDispersion} was called with \code{method="per-condition"}
      a name hasd to specified. Try \code{ls(cds@fitInfo}.
}
}
\seealso{\code{\link{estimateDispersions}} }
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   cds <- makeExampleCountDataSet()
   cds <- estimateSizeFactors( cds )
   cds <- estimateDispersions( cds )
   str( fitInfo( cds ) )
}   

