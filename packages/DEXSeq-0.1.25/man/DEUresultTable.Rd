\name{DEUresultTable}
\alias{DEUresultTable}
\title{
   Get a result table from the analysis workflow.
}
\description{
   This function returns a data frame with the summary of the results from the analysis workflow.  
   It accesses the \code{\link{fData}} slots with information of the dispersion estimates obtained
   from the function \code{\link{fitDispersionFunction}}, the p values, and adjusted p values
   obtained from the function \code{\link{testForDEU}}, and log2 fold changes obtained
   from the function \code{\link{estimatelog2FoldChanges}}.
}
\usage{
DEUresultTable(ecs)
}
\arguments{
   \item{ecs}{
      An ExonCountSet object.
   }
}
\value{
   A data frame with a summary of the analysis workflow.
}
\examples{
   \dontrun{
      data("pasillaExons", package="pasilla")
      pasillaExons <- estimateSizeFactors( pasillaExons )
      pasillaExons <- estimateDispersions( pasillaExons )
      pasillaExons <- fitDispersionFunction( pasillaExons )
      pasillaExons <- testForDEU( pasillaExons )
      res <- DEUresultTable( pasillaExons )
    }
}
