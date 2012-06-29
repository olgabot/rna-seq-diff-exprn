\name{DEXSeqHTML}
\alias{DEXSeqHTML}
\title{HTML report writer}
\description{
This function generates an HTML report from the results from 
\code{\link{testForDEU}} saved in an ExonCountSet object.
It uses the information from the function 
\code{\link{DEUresultTable}} and plotting from \code{\link{plotDEXSeq}}.
This gives an easy way of exploring the results of the tests.
}
     \usage{
DEXSeqHTML(ecs, geneIDs=NULL, path="DEXSeqReport", file="testForDEU.html", 
     fitExpToVar="condition", FDR=0.1, color=NULL, color.samples=NULL,
     mart="", filter="", attributes="", extraCols=NULL)
     }
     \arguments{
       \item{ecs}{An ExonCountSet object}
       \item{geneIDs}{
           A character vector of gene identificators to be included in the report. 
           If left NULL, the genes included in the report will be the significant
           hits at the given false discovery rate. See "FDR" below.}
       \item{path}{
           A path in the system where to write the report.
       }
       \item{file}{
           The name of the html file.
       }
       \item{fitExpToVar}{
           A variable contained in the design of the ecs; the counts 
           will be fitted to this variable to get the plotting values.
           (See \code{\link{plotDEXSeq}} for details.
       }
       \item{FDR}{A false discovery rate for the result.}
       \item{color}{A vector of colors, one for each of the levels of the values of "fitExpToVar".}
       \item{color.samples}{
           A vector of colors for each of the samples. If NULL, the colors 
           of each sample will be asigned according to its corresponding 
           condition. Useful to visualize complex experimental designs.
       }
       \item{mart}{object of class Mart, created with the useMart function, with dataset specified}
       \item{filter}{ 
          Filters (ONLY ONE) that should be used in the query.  A
          possible list of filters can be retrieved using the function
          listFilters. Please note that the value of this filter will
          always be the geneIDs in the ExonCountSet object.
       }
       \item{attributes}{
          Attributes you want to retrieve.  A possible list of
          attributes can be retrieved using the biomaRt function
          listAttributes.
       }
       \item{extraCols}{
          A data frame with one or more columns to add to the report.
          For example, additional information about the genes.
       }
     }
     \value{
	   This function will write an HTML report in the directory specified by 'path'. 
           There, it will create an html file with the initial report page and a 
           directory called "files" in which SCG files with the plots and other html files are
           placed.  To see an example please visit 
           \url{http://www.embl.de/~reyes/DEXSeqReport/testForDEU.html}.
     }
     \seealso{
       \code{hwrite}
     }
     \examples{
	\dontrun{
           data("pasillaExons", package="pasilla")
           pasillaExons <- estimateSizeFactors( pasillaExons )
           pasillaExons <- estimateDispersions( pasillaExons )
           pasillaExons <- fitDispersionFunction( pasillaExons )
           pasillaExons <- testForDEU( pasillaExons )
           DEXSeqHTML( pasillaExons )
        }
     }
\keyword{ExonCountSet}
