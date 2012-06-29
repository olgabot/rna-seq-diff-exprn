\name{estimateDispersions}
\Rdversion{1.1}
\alias{estimateDispersions,CountDataSet-method}
\alias{estimateDispersions}
\title{
    Estimate and fit dispersions for a CountDataSet.
}
\description{
    This function obtains dispersion estimates for a count data set. For each condition
    (or collectively for all conditions, see 'method' argument below) it first computes
    for each gene an empirical dispersion value (a.k.a. a raw SCV value), then fits by regression
    a dispersion-mean relationship and finally chooses for each gene a
    dispersion parameter that will be used in subsequent tests from the
    empirical and the fitted value according to the 'sharingMode'
    argument.
}
\usage{
\S4method{estimateDispersions}{CountDataSet}( object, 
   method = c( "pooled", "pooled-CR", "per-condition", "blind" ), 
   sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
   fitType = c("parametric", "local"),
   locfit_extra_args=list(), lp_extra_args=list(), 
   modelFrame = NULL, modelFormula = count ~ condition, ... )
}
\arguments{
  \item{object}{
      a \code{CountDataSet} with size factors.
  }
  \item{method}{
     There are three ways how the empirical dispersion can be computed:
          
     \itemize{
     \item \code{pooled} - Use the samples from all conditions with
     replicates to estimate a single pooled empirical dispersion value, 
     called "pooled", and assign it to all samples.

     \item \code{pooled-CR} - Fit models according to \code{modelFormula} and
     estimate the dispersion by maximizing a Cox-Reid adjusted profile
     likelihood (CR-APL). This method
     is much slower than \code{method=="pooled"} but works also with crossed
     factors (as may occur, e.g., in designs with paired samples). 
     Usually, you will need to specify the model formula, which should be the
     same as the one used later in the call to \code{nbinomFitGLMs} for fitting the full model. 

     Note: The method of using CR-APL maximization for this application has been
     developed by McCarthy, Chen and Smyth [Nucl. Acid Res., 2012 and been 
     first implemented in edgeR (in 2010). DESeq optimizes the expression
     for the CR-APL given in McCarthy et al.'s paper, but does not use the
     weigthed maximum likelihood scheme proposed there.

     \item \code{per-condition} - For each condition with replicates, compute
     a gene's empirical dispersion value by considering the data from samples for this
     condition.  For samples of unreplicated conditions, the maximum
     of empirical dispersion values from the other conditions is used.  If \code{object} has a
     multivariate design (i.e., if a data frame was passed instead of a
     factor for the \code{condition} argument in
     \code{\link{newCountDataSet}}), this method is not available. (Note: This method was called ``normal'' in
        previous versions.)

     \item \code{blind} - Ignore the sample labels and compute a 
     gene's empirical dispersion value as if all samples were replicates of a
     single condition. This can be done even if there are no biological
     replicates. This method can lead to loss of power; see
     the vignette for details. The single estimated dispersion condition
     is called "blind" and used for all samples.
    
     } }

   \item{sharingMode}{
     After the empirical dispersion values have been computed for each
     gene, a dispersion-mean relationship is fitted for sharing
     information across genes in order to reduce variability of the
     dispersion estimates. After that, for each gene, we have two values: the
     empirical value (derived only from this gene's data), and the
     fitted value (i.e., the dispersion value typical for genes with an
     average expression similar to those of this gene). The
     \code{sharingMode} argument specifies which of these two values
     will be written to the \code{featureData}'s \code{disp_} columns
     and hence will be used by the functions \code{\link{nbinomTest}}
     and \code{\link{fitNbinomGLMs}}.
     
     \itemize{     
     \item \code{fit-only} - use only the fitted value, i.e., the
     empirical value is used only as input to the fitting, and then
     ignored. Use this only with very \emph{few} replicates, and when
     you are not too concerned about false positives from dispersion outliers, i.e. genes
     with an unusually high variability.

     \item \code{maximum} - take the maximum of the two values. This is
     the conservative or prudent choice, recommended once you have at
     least three or four replicates and maybe even with only two replicates.

     \item \code{gene-est-only} - No fitting or sharing, use only the
     empirical value. This method is preferable when the number of
     replicates is large and the empirical dispersion values are
     sufficiently reliable. If the number of replicates is small, this
     option may lead to many cases where the dispersion
     of a gene is accidentally underestimated and a false positive arises in
     the subsequent testing.
     }
  }

  \item{fitType}{

     \itemize{
     \item \code{parametric} - Fit a dispersion-mean relation of the
     form \code{dispersion = asymptDisp + extraPois / mean} via a robust
     gamma-family GLM. The coefficients \code{asymptDisp} and \code{extraPois}
     are given in the attribute \code{coefficients} of the \code{dispFunc}
     in the \code{fitInfo} (see below).
     }

     \itemize{
     \item \code{local} - Use the locfit package to fit a dispersion-mean 
     relation, as described in the DESeq paper.
     }

  }

  \item{locfit_extra_args, lp_extra_args}{
      (only for \code{fitType=local})
      Options to be passed to the \code{locfit} and to the \code{lp}
      function of the locfit package. Use this to adjust the local
      fitting. For example, you may pass a value for \code{nn} different
      from the default (0.7) if the fit seems too smooth or too rough by
      setting \code{lp_extra_agrs=list(nn=0.9)}. As another example, you
      can set \code{locfit_extra_args=list(maxk=200)} if you get the
      error that locfit ran out of nodes. See the documentation of the
      \code{locfit} package for details. In most cases, you will not
      need to provide these parameters, as the defaults seem to work
      quite well.}
      
   \item{modelFrame}{
      By default, the information in \code{conditions(object)} or \code{pData(object)} is used to 
      determine which samples are replicates (see \code{\link{newCountDataSet}}).
      For \code{method="pooled"}, a data frame can be passed 
      here, and all rows that are identical in this data frame are considered
      to indicate replicate samples in \code{object}. For \code{method="pooled-CR"},
      the data frame is used in the fits. For the other methods, this argument is
      ignored.
   }

   \item{modelFormula}{
      For \code{method="pooled-CR"}, this is the formual used for the dispersion fits.
      For all other methods, this argument is ignored.
   }

   \item{...}{extra arguments are ignored}

}
\details{
   Behaviour for \code{method="per-condition"}: For each replicated condition, a list, named
   with the condition's name, is placed in the environment \code{object@fitInfo}. This list
   has five named elements: The vector \code{perGeneDispEsts} contains the
   empirical dispersions. The function \code{dispFunc}
   is the fitted function, i.e., it takes as its argument a normalized
   mean expression value and returns the corresponding
   fitted dispersion. The values fitted according to this function are
   in the third element \code{fittedDispEst}, a vector of the same
   length as \code{perGeneDispEsts}. The fourt element, \code{df},
   is an integer, indicating the number of degrees of freedom of
   the per-gene estimation. The fifth element, \code{sharingMode},
   stores the value of the \code{sharingMode} argument to
   \code{esimateDispersions}.

   Behaviour for \code{method="blind"} and \code{method="pooled"}: Only one list is produced,
   named \code{"blind"} or \code{"pooled"} and placed in \code{object@fitInfo}.
   
   For each list in the \code{fitInfo} environment, the dispersion
   values that are intended to be used in subsequent testing are computed according to
   the value of \code{sharingMode} and are placed in the
   \code{featureData} data frame, in a column 
   named with the same name, prefixed with "\code{disp_}".
   
   Then, the \code{\link{dispTable}} (see there) is filled to assign to each
   condition the appropriate dispersion column in the phenoData frame.
   
   Note: Up to DESeq version 1.4.x (Bioconductor release 2.8), this function was
   called \code{estimateVarianceFunctions}, stored its result differently and had
   did not yet have the arguments \code{sharingMode} and \code{fitType}.
   \code{estimatevarianceFunction}'s behaviour corresponded
   to the settings \code{sharingMode="fit-only"} and \code{fitType="local"}. Note that
   these are not the default, because the new functionalities sharingMode="maximum"
   and fitType="local" tend to give better results in standard cases as they are more robust
}

\value{
   The \code{CountDataSet} cds, with the slots \code{fitInfo} and
   \code{featureData} updated as described in Details.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
str( fitInfo( cds ) )
head( fData( cds ) )
}
