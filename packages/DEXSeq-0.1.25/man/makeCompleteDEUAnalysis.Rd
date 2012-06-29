\name{makeCompleteDEUAnalysis}
     \alias{makeCompleteDEUAnalysis}
     \title{Complete differential exon usage analysis}
     \description{
       This function performs a complete differential exon usage analysis, calling
       all the necessary functions and giving back an \code{ExonCountSet} object 
       with p values and p adjusted values.
     }
     \usage{
makeCompleteDEUAnalysis(ecs, 
        formulaDispersion=count ~ sample + condition*exon, 
        minCount=10, maxExon=50, formula0=NULL, formula1=NULL, 
        FDR=0.1, fitExpToVar="condition", nCores=1, path=NULL, 
        color=NULL, color.samples=NULL, quiet=FALSE, file="")
     }
     \arguments{
      \item{ecs}{
          An ExonCountSet object.
      }
      \item{formulaDispersion}{
         Formula used in the glm to calculate the dispersion values.  
         The factors on the formula must be present in the design 
         columns of the ExonCountSet object.
      }
       \item{minCount}{
          Minimum number of counts on an exon for it to be considered in the tests. 
          This significantly increases the speed of the dispersion 
          estimations and testing for differential exon usage. 
          This is supported by the fact that small count exons are less likely
          of being called significant, so it should not affect the results.
       }
      \item{maxExon}{
          Genes with more exons than this value will be discarded from 
          the analysis. This is a speed issue. Currently, time of dispersion 
          estimations and testing for differential exon usage
          increases with number of exons.
      }

      \item{formula0}{
         Formula for the NULL model to fit a the glm.  The factors must be 
         present in the design columns of the ExonCountSet object. As it is 
         tested for each of the exons, a factor exonID can be added to the 
         formula, so that it will iterate over the exons of the gene fitting 
         the glm for each of them. If it is left in NULL, the default formula is 
         "count~sample+exon+condition" for the NULL model.
      }
      \item{formula1}{
         Same as formula0, but for the test model.  If it is left in NULL, the 
         default formula will be "count~sample+exon+condition*I(exon==exonID)".  
         If added a factor "exonID", it will iterate over each of the exons of 
         the geneID, e.g. If a geneID contains exons E01, E02, E03,...,EN, 
         and the function is left in the default formula, the function will fit N glms, the 
         last part of the formula will change in the iterations as follows: 
         I(exon==E01), I(exon==E02), I(exon==E03),...,I(exon==EN).
      }
       \item{FDR}{
          A false discovery rate used to indicate the significant exons.
       }
      \item{fitExpToVar}{
          A variable contained in the design annotation of the \code{ExonCountSet}, 
          the expression values will be fitted to this variable using the formula
          count~fitExpToVar*exon using a model frame obtained from the function
          \code{\link{modelFrameForGene}}.  
       }
      \item{nCores}{
          Number of cores to be used to estimate the dispersions.  
          multicore package must be loaded in order to split the job 
          in several cores.
       }
       \item{path}{
          A path in the system where to write the report from 
          \code{\link{DEXSeqHTML}}. If NULL, no report will be created.
       }
       \item{color}{
          A vector of colors for each of the levels from the factor in the design of 
          the ExonCountSet object indicated by "fitExpToVar". If path is NULL,
          this parameter will be ignored.
       }
       \item{color.samples}{
          A vector of colors for each of the samples. If NULL, the colors of each 
          sample will be asigned according to its corresponding level from "fitExpToVar". 
          This option is useful to visualize complex experimental designs. If path is NULL,
          this parameter will be ignored.
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
     An object of class ExonCountSet.
  }
  \examples{
data("pasillaExons", package="pasilla")
formuladispersion <- count ~ sample + ( exon + type ) * condition
formula0 <- count ~ sample + type * exon + condition
formula1 <- count ~ sample + type * exon + condition * I(exon == exonID)
pasillaExons <- makeCompleteDEUAnalysis(pasillaExons, 
   formulaDispersion=formuladispersion, 
   formula0=formula0, 
   formula1=formula1)
   }
\keyword{ExonCountSet}
