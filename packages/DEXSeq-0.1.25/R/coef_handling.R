arrangeCoefs <- function( frm, mf, mm = model.matrix( frm, mf ), fit = NULL, insertValues = TRUE ) {

   if( any( attr( mm, "contrasts" ) != "contr.treatment" ) )
      stop( "Can only deal with standard 'treatment' contrasts." )   # Do I need this?
   if( is.null(fit) & insertValues )
      stop( "If fit==NULL, returnCoefValues must be FALSE" )
   if( !is.null(fit) )
      stopifnot( all( colnames(mm) == names(coefficients(fit)) ) )

   fctTbl <- attr( terms(frm), "factors" )

   coefIndicesList <- 
   lapply( seq_len(ncol(fctTbl)), function( fctTblCol ) {
      termName <- colnames(fctTbl)[ fctTblCol ]
      varsInTerm <- stringr::str_split( termName, stringr::fixed(":") )[[1]] 
      stopifnot( all( fctTbl[ varsInTerm, fctTblCol ] == 1 ) )
      stopifnot( sum( fctTbl[ , fctTblCol ] ) == length( varsInTerm ) )
      coefNames <- colnames(mm)[ attr( mm, "assign" ) == fctTblCol ]
      lvlTbl <- stringr::str_match( coefNames, 
         stringr::str_c( "^", stringr::str_c( varsInTerm, "([^:]*)", collapse=":" ), "$" ) )[ , -1, drop=FALSE ]
      stopifnot( ncol(lvlTbl) == length( varsInTerm ) )
      stopifnot( nrow(lvlTbl) == length( coefNames ) )
      if( !all( sapply( varsInTerm, function(v) is.factor(mf[[v]]) | is.character(mf[[v]]) ) ) )
         stop( "Non-factor in model frame" )

      varLevels <- lapply( varsInTerm, function(v) levels( factor( mf[[v]] ) ) ) 
      coefIndices <- array( NA_character_, dim = sapply( varLevels, length ), dimnames = varLevels )
      names( dimnames( coefIndices ) ) <- varsInTerm

      for( i in seq_len( nrow(lvlTbl) ) )
         coefIndices <- do.call( `[[<-`, c( quote(coefIndices), as.list( lvlTbl[ i, ] ), coefNames[i] ) )

      coefIndices
   } )
   names( coefIndicesList ) <- colnames( fctTbl )

   if( attr( terms(frm), "intercept" ) ) {
      a <- array( c( `(Intercept)` = "(Intercept)" ) )
      dimnames(a) <- list( `(Intercept)` = c( "(Intercept)" ) )
      coefIndicesList <- c( list( `(Intercept)` = a ), coefIndicesList )
   }

   if( !insertValues )
      ans <- coefIndicesList
   else
      ans <- lapply( coefIndicesList, function(coefIndices) {
         a <- ifelse( is.na(coefIndices), 0, coefficients(fit)[ coefIndices ] )
         attr( a, "variables" ) <- attr( coefIndices, "variables" )
         a } )
      
   lapply( ans, function(x) 
      if( is.array(x) ) 
         x 
      else { 
         y <- array( x, dim=length(x) )
         attr( y, "variables" ) <- attr( x, "variables" )
         dimnames(y) <- list( names(x) )
         y } )
}

apply2 <- function( X, MARGIN, FUN, ... ) {
   if( length(MARGIN) > 0 ) 
      apply( X, MARGIN, FUN, ... ) 
   else 
      FUN( X, ... ) }

balanceExons <- function( coefs, dispersions ) {
   stopifnot( any( sapply( coefs, function(x) 
      identical( names(dimnames(x)), "(Intercept)" ) ) ) )
   termsWithExon <- sapply( coefs, function(x) "exon" %in% names(dimnames(x)) )
   meanMainEffect <- sum( sapply( coefs[!termsWithExon], mean, na.rm=TRUE ) )
   meanExonEffects <- rowSums( sapply( coefs[termsWithExon], function(x) 
      apply2( x, "exon", mean, na.rm=TRUE ) ) )

   meanExonFittedValue <- exp( meanMainEffect + meanExonEffects )

   exonWeights <-  1 / ( dispersions + 1 / meanExonFittedValue )

   shifts <- lapply( coefs[termsWithExon], function(x) { 
      nonExonDims <- which(  names(dimnames(x)) != "exon" )
      list(
         vars = names(dimnames(x))[ nonExonDims ],
         wmeans = apply2( x, nonExonDims, weighted.mean, exonWeights) ) } )

   lapply( coefs, function(x) {
      nonExonVars <- names(dimnames(x))[ names(dimnames(x)) != "exon" ]
      if( identical( nonExonVars, "(Intercept)" ) )
         whichShift <- which( sapply( shifts, function(xx) length( xx$vars ) == 0 ) )
      else
         whichShift <- which( sapply( shifts, function(xx) identical( xx$vars, nonExonVars ) ) )
      if( length( whichShift ) == 0 )
         return( x )
      if( length( whichShift ) > 1 )
         stop( "Confused about selecting shift." )
      if( "exon" %in% names(dimnames(x)) )
         x - shifts[[ whichShift ]]$wmeans
      else
         x + shifts[[ whichShift ]]$wmeans
    } )
}         


fitAndArrangeCoefs <- function( ecs, geneID, frm = count ~ condition * exon, balanceExons = TRUE )
{
   mf <- modelFrameForGene( ecs, geneID )
   if( length(levels(mf$exon)) <= 1 )
      return( NULL )
   mm <- model.matrix( frm, mf )
   fit <- try( glmnb.fit(mm, mf$count, dispersion=mf$dispersion, offset=log(mf$sizeFactor)), silent=TRUE)
   if( is( fit, "try-error" ) )
      return( NULL )
   coefs <- arrangeCoefs( frm, mf, mm, fit )
   if( balanceExons ) 
      balanceExons( coefs, tapply( mf$dispersion, mf$exon, `[`, 1 ) )
   else
      coefs
}

getEffectsForPlotting <- function( coefs, groupingVar = "condition", averageOutExpression=FALSE )
{
   groupingExonInteraction <- which( sapply( coefs, function(x) 
      all( c( groupingVar, "exon") %in% names(dimnames(x)) ) & length(dim(x)) == 2 ) ) 
   fittedValues <- coefs[[ groupingExonInteraction ]]
   if( names(dimnames(fittedValues))[1] == "exon" )
      fittedValues <- t( fittedValues )
   stopifnot( identical( names(dimnames(fittedValues)), c( groupingVar, "exon" ) ) )
   for( x in coefs[ -groupingExonInteraction ] ) {
      if( all( c( groupingVar, "exon") %in% names(dimnames(x)) ) )
         stop( "Cannot yet deal with third-order terms." )
      if( !any( c( groupingVar, "exon") %in% names(dimnames(x)) ) ) {
         fittedValues <- fittedValues + mean( x )
      } else if( averageOutExpression & identical( names(dimnames(x)), groupingVar ) ) {
         fittedValues <- fittedValues + mean( x )
      } else if( groupingVar %in% names(dimnames(x)) ) {
         groupMeans <- apply2( x, groupingVar, mean )
         stopifnot( identical( names(groupMeans), dimnames(fittedValues)[[1]] ) )
         fittedValues <- fittedValues + groupMeans
      } else if( "exon" %in% names(dimnames(x)) ) {
         exonMeans <- apply2( x, "exon", mean )
         fittedValues <- t( t(fittedValues) + exonMeans )
      } else {
         print( x )
         stop( "Unexpected term encountered." )
      }
   }
   fittedValues
}

# Tesing code:

#load( "Graveley_unstranded_fd.rda" )
#frm <- count ~ ( condition + libType ) * exon
#i <- 0
#coefs <- new.env( hash=TRUE )
#for( geneID in unique(geneIDs(ecs)) ) {
#   coefs[[ geneID ]] <- fitAndArrangeCoefs( ecs, frm )
#   i <- i + 1
#   if( i %% 100 == 0 )
#      cat(sprintf( "%d genes processed.\n", i ))
#}

