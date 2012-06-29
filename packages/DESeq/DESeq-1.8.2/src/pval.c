#include <math.h>
#include <limits.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#define MAX( a, b ) ( a > b ? a : b )
#define MIN( a, b ) ( a < b ? a : b )

void add_from_both_sides( long int left, long int right, long int kS, double pobs, 
      double muA, double vA, double muB, double vB, double eps,
      double * res_total, double * res_obstotal)
{
   double sizeA = muA*muA / (vA-muA);
   double probA = muA/vA;
   double sizeB = muB*muB / (vB-muB);
   double probB = muB/vB;

   long int kl = left;
   long int kr = right;
   double lval = Rf_dnbinom( kl, sizeA, probA, 0 ) * Rf_dnbinom( kS-kl, sizeB, probB, 0 );
   double rval = Rf_dnbinom( kr, sizeA, probA, 0 ) * Rf_dnbinom( kS-kr, sizeB, probB, 0 );
   double prevlval = lval;
   double prevrval = rval;
   double total = lval + rval;
   double esttotalperlength = total/2;
   double obstotal = 0;   
   long int step = 1;
   long int steps = 0;
   int do_left;
   if( lval <= pobs )
      obstotal += lval;
   if( rval <= pobs )
      obstotal += rval;
   while( kl < kr ) {
      steps ++;
      if( fabs( ( prevrval - rval ) ) / prevrval > .01 )
         do_left = 1;
      else if( fabs( ( prevlval - lval ) ) / prevlval > .01 )
         do_left = 0;
      else
         do_left = lval > rval;
      if( do_left ) {
         prevlval = lval;
         if( kl + step > kr )
            step = kr - kl;
         kl += step;
         lval = Rf_dnbinom( kl, sizeA, probA, 0 ) * Rf_dnbinom( kS-kl, sizeB, probB, 0 );
         if( step == 1 )
            total += lval;   
         else
            total += MIN( lval, prevlval ) * step;   
         if( lval <= pobs ) {
            if( step == 1 )
               obstotal += lval;
            else {       
               if( prevlval <= pobs )
                  obstotal += MAX( lval, prevlval ) * step;
               else
                  obstotal += MAX( lval, prevlval ) * step * fabs( (pobs-lval) / (prevlval-lval) );
            }
         }       
         if( fabs( prevlval - lval ) / prevlval < eps )
            step = MAX( step + 1, (long int) step * 1.5 );
      } else {
         prevrval = rval;
         if( kr - step < kl )
            step = kr - kl;
         kr -= step;
         rval = Rf_dnbinom( kr, sizeA, probA, 0 ) * Rf_dnbinom( kS-kr, sizeB, probB, 0 );
         if( step == 1 )
            total += rval;   
         else
            total += MIN( rval, prevrval ) * step;   
         if( rval <= pobs ) {
            if( step == 1 )
               obstotal += rval;
            else {       
               if( prevrval <= pobs )
                  obstotal += MAX( rval, prevrval ) * step;
               else
                  obstotal += MAX( rval, prevrval ) * step * fabs( (pobs-rval) / (prevrval-rval) );
            }
         }       
         if( fabs( ( prevrval - rval ) ) / prevrval < eps )
            step = MAX( step + 1, (long int) step * 1.5 );
      }
   }
   //assert( kl == kr );
   
   *res_total = total;
   *res_obstotal = obstotal;
}   


SEXP calc_pvals( SEXP kS, SEXP pobs, SEXP muA, SEXP vA, 
   SEXP muB, SEXP vB, SEXP eps ) 
{
   double totall, totalr, obstotall, obstotalr;
   long int kexp = INTEGER(kS)[0] * REAL(muA)[0] / ( REAL(muA)[0] + REAL(muB)[0] );   
   add_from_both_sides( 0, kexp, INTEGER(kS)[0], REAL(pobs)[0], REAL(muA)[0], 
      REAL(vA)[0], REAL(muB)[0], REAL(vB)[0], REAL(eps)[0], &totall, &obstotall);
   add_from_both_sides( kexp+1, INTEGER(kS)[0], INTEGER(kS)[0], REAL(pobs)[0], REAL(muA)[0], 
      REAL(vA)[0], REAL(muB)[0], REAL(vB)[0], REAL(eps)[0], &totalr, &obstotalr);
   SEXP res = Rf_allocVector( REALSXP, 2 );
   REAL(res)[0] = totall + totalr;
   REAL(res)[1] = obstotall + obstotalr;
   return res;
}   

R_CallMethodDef callMethods[] = {
   { "calc_pvals", (DL_FUNC) &calc_pvals, 7 },
   { NULL, NULL, 0 }
};

void R_init_DESeq( DllInfo *info )
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
