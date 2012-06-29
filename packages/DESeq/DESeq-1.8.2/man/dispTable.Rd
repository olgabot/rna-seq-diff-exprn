\name{dispTable}
\Rdversion{1.1}
\alias{dispTable}
\alias{dispTable,CountDataSet-method}
\title{
    Accessor function for the dispTable information in a CountDataSet
}
\description{
   The dispersion table ("dispTable") is a named vector that assigns to each
   condition (as name) a dispersion column (as value). If \code{\link{nbinomTest}}
   is called to compare two conditions, say "A" and "B", DESeq looks up in the dispTable,
   which dispersion columns to use. In the standard case (see example), these
   are just the dispersions for "A" and "B", i.e., the columns \code{disp_A}
   and \code{disp_B} in \code{fData(object)}. If the "pooled" or "blind" variance
   estimation is used, all conditions are assigned the same column. 
}
\usage{
dispTable(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
      a CountDataSet
}
  \item{...}{
      further argumnts are ignored
}
}
\seealso{\code{\link{estimateDispersions}}, \code{\link{nbinomTest}} }
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
   cds <- makeExampleCountDataSet()
   cds <- estimateSizeFactors( cds )
   cds <- estimateDispersions( cds )
   dispTable( cds )
}   

