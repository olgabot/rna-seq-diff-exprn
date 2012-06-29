\name{testForDEU}
\alias{testForDEU}
\title{Test for Differential Exon Usage.}
\description{
      This function tests for differential exon usage for each of 
      the genes in the object. It stores the results in the fields 
      \code{fData(ecs)$pvalue} and \code{fData(ecs)$padjust}.
}
\usage{
testForDEU(ecs, formula0=NULL, formula1=NULL, nCores=1, quiet=FALSE, file="")
}
\arguments{
  \item{ecs}{
     An ExonCountSet object.
  }
  \item{formula0}{
     Formula for the null model to be used in the GLM fit. If no formula
     is given, the default \code{count ~ sample + exon + condition} is used.
     See below for details
  }
  \item{formula1}{
     Formula for the full model to be used in the GLM fit. If no formula
     is given, the default \code{count ~ sample + exon + condition *I (exon==exonID)} is used.
     See below for details.
  }
  \item{nCores}{
     Number of CPUcores to be used to estimate the dispersions.  
     The \code{multicore} package must be loaded to use more than 1 core.
  }
  \item{quiet}{
     If TRUE, no progress report is shown.  In case the session is 
     not an interactive session and progress report is wanted.
     Change the name of the file.
  }
  \item{file}{
      A file name to write the progress reports. If file="", output
      will be written in the standard output connection.
  }
}
\details{
   The terms in the formulas must be columns of \code{design(ecs)}. In addition,
   in \code{formula1}, the variable \code{exonID} is set to the
   ID of the currently tested exon counting bin.
   
   See \code{\link{testGeneForDEU}}, which is called for each gene, for further details.
   
}
\value{
     An ExonCountSet object with \code{fData(ecs)$pvalue} and 
     \code{fData(ecs)$padjust} data slots filled.
}
\examples{
     \dontrun{
        data("pasillaExons", package="pasilla")
        pasillaExons <- estimateSizeFactors( pasillaExons )
        pasillaExons <- estimateDispersions( pasillaExons )
        pasillaExons <- fitDispersionFunction( pasillaExons )
        pasillaExons <- testForDEU( pasillaExons )
     }
}
\seealso{
   \code{estimateExonDispersionsForModelFrame}
}
