#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP get_occupancy(SEXP);
int factorial(int);
int rec_fact(int);

SEXP get_occupancy(SEXP pwm) 
{
   //SEXP matrix_slot = mkChar("profileMatrix\0");
   //PrintValue(matrix_slot);
 
   // SEXP dim = GET_DIM(pwm); // --> this is NULL for my pwm. Why?
   SEXP mat = GET_SLOT(pwm, install("profileMatrix"));
   int ncol = INTEGER(GET_DIM(mat))[0];
   int nrow = INTEGER(GET_DIM(mat))[1];
   // IFDEF DEBUG
   Rprintf("ncol %d\n", ncol);
   Rprintf("nrow %d\n", nrow);
   //PrintValue(mat);
   // ENDIF
   
   SEXP res = PROTECT(allocVector(INTSXP,3));
   int *p = INTEGER(res);
   p[0] = ncol;
   p[1] = nrow;
   int i = 0;
   int j = 2;
   p[2] = (int)(REAL(mat)[i*nrow+j]); // access element at row i and col j (0 based)
   UNPROTECT(1);
   return(res);
}
