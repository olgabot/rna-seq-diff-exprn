setMethod("estimateSizeFactors", signature(object="CountDataSet"),
  function( object, locfunc=median, ... ) {
    if( length(list(...)) != 0 )
       warning( "in estimateSizeFactors: Ignoring extra argument(s)." )
    sizeFactors(object) <- estimateSizeFactorsForMatrix( counts(object), locfunc )
    object
  })

setMethod("estimateDispersions", signature(object="CountDataSet"),
  function( object, method = c( "pooled", "pooled-CR", "per-condition", "blind" ), 
    sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
    fitType = c( "parametric", "local" ),
    locfit_extra_args=list(), lp_extra_args=list(), 
    modelFrame = NULL, modelFormula = count ~ condition, ... )
{
   stopifnot( is( object, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(object) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   method <- match.arg( method )
   sharingMode <- match.arg( sharingMode )
   fitType <- match.arg( fitType )
   if( length(list(...)) != 0 )
      warning( "in estimateDispersions: Ignoring extra argument(s)." )
   if( object@multivariateConditions && ! method %in% c( "blind", "pooled", "pooled-CR" ) )
      stop( "You have specified multivariate conditions (i.e., passed a data frame with conditions). In this case, you cannot use method 'per-condition'." )
   if( sharingMode == "gene-est-only" )
      warning( "in estimateDispersions: sharingMode=='gene-est-only' will cause inflated numbers of false positives unless you have many replicates." )
   ## FIXME this warning should only be emitted when the number of replicates is indeed small. 
   
   # Remove results from previous fits
   fData(object) <- fData(object)[ , ! colnames(fData(object)) %in% paste( "disp", object@dispTable, sep="_" ), drop=FALSE ]
   object@dispTable <- character()
   object@fitInfo = new.env( hash=TRUE )
   
   if( method == "blind" ) {

      bmv <- getBaseMeansAndVariances( counts(object), sizeFactors(object) )
      dispsAndFunc <- estimateAndFitDispersionsFromBaseMeansAndVariances( bmv$baseMean,
         bmv$baseVar, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )   
      object@fitInfo[[ "blind" ]] <- list( 
         perGeneDispEsts = dispsAndFunc$disps,
         dispFunc = dispsAndFunc$dispFunc,
         fittedDispEsts = dispsAndFunc$dispFunc( bmv$baseMean ),
         df = ncol(counts(object)) - 1,
         sharingMode = sharingMode )
      
      if( object@multivariateConditions )
         dispTable(object) <- c( "_all" = "blind" )
      else {
         a <- rep( "blind", length( levels( conditions(object) ) ) )
         names(a) <- levels( conditions(object) )
         object@dispTable <- a }
      
   } else if( method == "per-condition" ) {
   
      replicated <- names( which( tapply( conditions(object), conditions(object), length ) > 1 ) )
      if( length( replicated ) < 1 )
         stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors." )
      nonreplicated <- names( which( tapply( conditions(object), conditions(object), length ) == 1 ) )
      overall_basemeans <- rowMeans( counts( object, normalized=TRUE ) )

      for( cond in replicated ) {
         cols <- conditions(object)==cond
         bmv <- getBaseMeansAndVariances( counts(object)[ , cols ], sizeFactors(object)[ cols ] )
         dispsAndFunc <- estimateAndFitDispersionsFromBaseMeansAndVariances( bmv$baseMean,
            bmv$baseVar, sizeFactors(object)[cols], fitType, locfit_extra_args, lp_extra_args )   
         object@fitInfo[[ cond ]] <- list( 
            perGeneDispEsts = dispsAndFunc$disps,
            dispFunc = dispsAndFunc$dispFunc,
            fittedDispEsts = dispsAndFunc$dispFunc( overall_basemeans ),     # Note that we do not use bmv$baseMean here
            df = sum(cols) - 1,
            sharingMode = sharingMode ) }
         
      object@dispTable <- sapply( levels(conditions(object)), function( cond )
            ifelse( cond %in% replicated, cond, "max" ) ) 
                        
   } else if( method == "pooled" || method == "pooled-CR" ) { 
   
      if( method == "pooled" ) {
   
         if( object@multivariateConditions ) {
            if( is.null( modelFrame ) )
               modelFrame <- pData(object)[ , colnames(pData(object)) != "sizeFactor" ]
            conds <- modelMatrixToConditionFactor( modelFrame ) }
         else
            conds <- conditions(object)
         if( !any( duplicated( conds ) ) )
            stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors." )
         bmv <- getBaseMeansAndPooledVariances( counts(object), sizeFactors(object), conds )
         baseMeans <- bmv$baseMean
         dispsAndFunc <- estimateAndFitDispersionsFromBaseMeansAndVariances( bmv$baseMean,
            bmv$baseVar, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
         df <- ncol(counts(object)) - length(unique(conds))
      
      } else {  # method == "pooled-CR"
         if( is.null( modelFrame ) )
            modelFrame <- pData(object)[ , colnames(pData(object)) != "sizeFactor", drop=FALSE ]
         baseMeans <- rowMeans( counts( object, normalized=TRUE ) )
         
         dispsAndFunc <- estimateAndFitDispersionsWithCoxReid( counts(object), modelFormula, modelFrame,
            sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )      
         df <- NA
      }

      object@fitInfo[[ "pooled" ]] <- list( 
         perGeneDispEsts = dispsAndFunc$disps,
         dispFunc = dispsAndFunc$dispFunc,
         fittedDispEsts = dispsAndFunc$dispFunc( baseMeans ),
         df = df,
         sharingMode = sharingMode )

      if( object@multivariateConditions )
         dispTable(object) <- c( "_all" = "pooled" )
      else {
         a <- rep( "pooled", length( levels( conditions(object) ) ) )
         names(a) <- levels( conditions(object) )
         object@dispTable <- a }
   
   } else
   
      stop(sprintf("Invalid method '%s'.", method))
   
   for( n in ls(object@fitInfo) )
      fData(object)[[ paste( "disp", n, sep="_" ) ]] <- 
         switch( sharingMode, 
            `fit-only`      = object@fitInfo[[ n ]]$fittedDispEsts,
            `gene-est-only` = {
               a <- object@fitInfo[[ n ]]$perGeneDispEsts
               a[ is.nan(a) ] <- 0
               pmax( a, 1e-8 ) },
            `maximum`       = pmax( object@fitInfo[[ n ]]$fittedDispEsts, object@fitInfo[[ n ]]$perGeneDispEsts, na.rm=TRUE ),
            stop(sprintf("Invalid sharingMode '%s'.", sharingMode))
         ) ## switch
        
   if( "max" %in% object@dispTable )
      fData(object)[["disp_max"]] <- do.call( pmax, 
         c( fData(object)[ , colnames(fData(object)) %in% paste( "disp", object@dispTable, sep="_" ), drop=FALSE ], na.rm=TRUE ) )
        
        
   validObject( object )
   object
})

estimateVarianceFunctions <- function( ... )
{
   stop( "The function 'estimateVarianceFunctions' has been removed. Use 'estimateDispersions' intead." )
}

varianceFitDiagnostics <- function( ... )
{
   stop( "This function has been removed. Please do not use it anymore. See the vignette for our current suggestions to check fit quality." )
}

residualsEcdfPlot <- function( ... )
{
   stop( "This function has been removed. Please do not use it anymore. See the vignette for our current suggestions to check fit quality." )
}

nbinomTest <- function( cds, condA, condB, pvals_only=FALSE, eps=NULL )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( cds@multivariateConditions )
      stop( "For CountDataSets with multivariate conditions, only the GLM-based test can be used." )
   if( all( is.na( dispTable(cds) ) ) )
      stop( "Call 'estimateDispersions' first." )

   if( dispTable(cds)[condA] == "blind" || dispTable(cds)[condB] == "blind" ) {
      if( fitInfo( cds, "blind" )$sharingMode != "fit-only" )
      warning( 'You have used \'method="blind"\' in estimateDispersion without also setting \'sharingMode="fit-only"\'. This will not yield useful results.' )
   }

   stopifnot( condA %in% levels(conditions(cds)) )  
   stopifnot( condB %in% levels(conditions(cds)) )     
   if( !is.null(eps) )
      warning( "The 'eps' argument is defunct and hence ignored." )
   
   colA <- conditions(cds)==condA
   colB <- conditions(cds)==condB

   bmv <- getBaseMeansAndVariances( counts(cds)[,colA|colB], 
      sizeFactors(cds)[colA|colB] )
   
   rawScvA <- fData(cds)[ , paste( "disp", dispTable(cds)[condA], sep="_" ) ]
   rawScvB <- fData(cds)[ , paste( "disp", dispTable(cds)[condB], sep="_" ) ]

   pval <- nbinomTestForMatrices(
      counts(cds)[,colA], 
      counts(cds)[,colB], 
      sizeFactors(cds)[colA], 
      sizeFactors(cds)[colB], 
      rawScvA, 
      rawScvB )
      
   if( pvals_only )
      pval
   else {
      bmvA <- getBaseMeansAndVariances( counts(cds)[,colA], sizeFactors(cds)[colA] )
      bmvB <- getBaseMeansAndVariances( counts(cds)[,colB], sizeFactors(cds)[colB] )
      data.frame( 
         id    = rownames( counts(cds) ),
         baseMean  = bmv$baseMean,
         baseMeanA = bmvA$baseMean,
         baseMeanB = bmvB$baseMean,
         foldChange = bmvB$baseMean / bmvA$baseMean,
         log2FoldChange = log2( bmvB$baseMean / bmvA$baseMean ), 
         pval = pval,
         padj = p.adjust( pval, method="BH" ), 
         # resVarA = cds@fitInfo[[ dispTable(cds)[condA] ]]$perGeneDispEsts / rawScvA,
         # resVarB = cds@fitInfo[[ dispTable(cds)[condB] ]]$perGeneDispEsts / rawScvB,
         stringsAsFactors = FALSE ) }
}

scvPlot <- function( ... )
{
   stop( "This function has been removed. Please do not use it anymore. See the vignette for our current suggestions to check fit quality." )
}


getVarianceStabilizedData <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( "blind" %in% ls(cds@fitInfo) )
      fitInfo <- cds@fitInfo[["blind"]]
   else if( "pooled" %in% ls(cds@fitInfo) ) 
      fitInfo <- cds@fitInfo[["pooled"]]
   else
      stop( "Use 'estimateDispersions' with 'method=\"blind\"' (or \"pooled\") before calling 'getVarianceStabilizedData'" )
   ncounts <- t( t(counts(cds)) / sizeFactors(cds) )
   if( attr( fitInfo$dispFunc, "fitType" ) == "parametric" ) {
      coefs <- attr( fitInfo$dispFunc, "coefficients" )
      vst <- function( q )
         log( (1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * q + 
           2 * sqrt( coefs["asymptDisp"] * q * ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * q ) ) )
           / ( 4 * coefs["asymptDisp"] ) ) / log(2)
      vst( ncounts )
   } else {  
      # non-parametric fit -> numerical integration
      xg <- sinh( seq( asinh(0), asinh(max(ncounts)), length.out=1000 ) )[-1]
      xim <- mean( 1/sizeFactors(cds) )
      baseVarsAtGrid <- fitInfo$dispFunc( xg ) * xg^2 + xim * xg 
      integrand <- 1 / sqrt( baseVarsAtGrid )
      splf <- splinefun( 
         asinh( ( xg[-1] + xg[-length(xg)] )/2 ), 
         cumsum( 
            ( xg[-1] - xg[-length(xg)] ) * 
            ( integrand[-1] + integrand[-length(integrand)] )/2 ) )
      h1 <- quantile( rowMeans(ncounts), .95 )
      h2 <- quantile( rowMeans(ncounts), .999 )
      eta <- ( log2(h2) - log2(h1) ) / ( splf(asinh(h2)) - splf(asinh(h1)) )
      xi <- log2(h1) - eta * splf(asinh(h1))
      tc <- sapply( colnames(counts(cds)), function(clm)
         eta * splf( asinh( ncounts[,clm] ) ) + xi )
      rownames( tc ) <- rownames( counts(cds) )
      tc
   }
}

makeExampleCountDataSet <- function( ) 
{
   ngenes <- 10000
   q0 <- rexp( ngenes, rate=1/250 )
   is_DE <- runif( ngenes ) < .3
   lfc <- rnorm( ngenes, sd=2 )
   q0A <- ifelse( is_DE, q0 * 2^(  lfc/2 ), q0 )
   q0B <- ifelse( is_DE, q0 * 2^( -lfc/2 ), q0 )
   true_sf <- c( 1., 1.3, .7, .9, 1.6 )   
   conds <- c( "A", "A", "B", "B", "B" )
   m <- t( sapply( 1:ngenes, function(i) 
      sapply( 1:5, function( j )
         rnbinom( 1, mu = true_sf[j] * ifelse( conds[j]=="A", q0A[i], q0B[i] ), 
            size = 1/.2 ) ) ) )
   colnames(m) <- c( "A1", "A2", "B1", "B2", "B3" )
   rownames(m) <- paste( "gene", 1:ngenes, 
      ifelse( is_DE, "T", "F" ), sep="_" )
   newCountDataSet( m, conds )
}

fitNbinomGLMs <- function( cds, modelFormula, glmControl=list() )
{
   stopifnot( is( cds, "CountDataSet" ) )

   if( "disp_pooled" %in% colnames( fData(cds) ) )
      disps <- fData(cds)$disp_pooled
   else if( "disp_blind" %in% colnames( fData(cds) ) ) {
      if( fitInfo( cds, "blind" )$sharingMode != "fit-only" )
         warning( 'You have used \'method="blind"\' in estimateDispersion without also setting \'sharingMode="fit-only"\'. This will not yield useful results.' )
      disps <- fData(cds)$disp_blind
   } else
      stop( "Call 'estimateDispersions' with 'method=\"pooled\"' (or 'blind') first." )

   fitNbinomGLMsForMatrix( counts(cds), sizeFactors(cds), disps, 
      modelFormula, pData(cds), glmControl=glmControl )
}

nbinomGLMTest <- function( resFull, resReduced )
   1 - pchisq( resReduced$deviance - resFull$deviance, 
   attr( resReduced, "df.residual" ) - attr( resFull, "df.residual" ) )


