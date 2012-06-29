\name{estimateDispersions}
\docType{methods}
\alias{estimateDispersions,ExonCountSet-method}
\alias{estimateDispersions}
\title{Estimate exon dispersions}
\description{
      This function estimates the for each counting bin
      of the ExonCountSet object a dispersion value. It stores 
      these values in \code{fData(ecs)$dispersionBeforeSharing}.  
}
\usage{
\S4method{estimateDispersions}{ExonCountSet}( object, 
         formula=count ~ sample + condition * exon, 
         initialGuess=.01, nCores=1, minCount=10, 
         maxExon=70, quiet=FALSE, file="")
}
\arguments{
  \item{object}{
     An ExonCountSet object.
   }
  \item{formula}{
     Formula used in the GLM to estimate the dispersion values.  
     The terms in the formula must be design 
     columns of the ExonCountSet object, the l.h.s. must be \code{count}.
   }
  \item{initialGuess}{
     An initial guess for the dispersion values to initiate the optimization.
   }
  \item{nCores}{
     Number of cores to be used to estimate the dispersions.  
     The \code{multicore} package must be loaded in order to spread the 
     job onto several cores.
   }
  \item{minCount}{
     Counting bins with less than 'minCount' counts (summed over all
     samples) are skipped in the tests. This reduced computation
     time, as counting bins with very few countsare unlikely to give 
     a significant signal anyway. For skipped counting bins, the
     'testable' column in fData is set to FALSE.
   }
  \item{maxExon}{
     Genes with more than 'maxExon' counting bins will be skipped in the
     test. This may be necessary because genes with very many counting bins
     may take extremely long in dispersion estimation and testing for 
     differential exon usage.
   }
   \item{quiet}{
     If TRUE, no progress report is shown.  In case the session is 
     not an interactive session and progress report is wanted, include
     a file name in the parameter "file".
   }
   \item{file}{
     A file name to write the progress reports. If file="", output
     will be written in the standart output connection.
   }
}
\value{
     An ExonCountSet object with dispersion \code{featureData(object)$dispersion_CR_est}) 
     parameters filled).
}
\details{For the dispersion estimation, we use the Cox-Reid conditional maximum
likelihood method of Gordon Smyth et al., which they devised for the \code{edgeR}
package. 
}
\examples{
     \dontrun{
        data("pasillaExons", package="pasilla")
        pasillaExons <- estimateSizeFactors( pasillaExons )
        pasillaExons <- estimateDispersions( pasillaExons )
     }
}
