#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "total_affinity.h"

SEXP get_occupancy(SEXP);

SEXP get_occupancy(SEXP pwm) 
{
   //SEXP matrix_slot = mkChar("profileMatrix\0");
   //PrintValue(matrix_slot);
 
   // SEXP dim = GET_DIM(pwm); // --> this is NULL for my pwm. Why?
   SEXP mat = GET_SLOT(pwm, install("profileMatrix"));
   matrix_ll mat_ll = NULL;
   convert_PWMMatrix_to_matrix_ll(mat, mat_ll);
   
   //Matrix name management: does not work right now. Is it needed?
   /*
   Rprintf(GET_SLOT(pwm, install("ID")));
   const char *ID = CHAR(GET_SLOT(pwm, install("ID")));
   mat_ll->name = (char *) R_alloc(MAX_MATRIX_NAME, sizeof(char));
   strcpy(mat_ll->name, ID);
   Rprintf("ID %s\n",mat_ll->name);
   */
   
   SEXP res = PROTECT(allocVector(INTSXP,3));
   int *p = INTEGER(res);
   p[0] = 42;
   p[1] = 42;
   int i = 0;
   int j = 2;
   p[2] = (int)(REAL(mat)[i*4+j]); // access element at row i and col j (0 based)
   UNPROTECT(1);
   return(res);
}

void convert_PWMMatrix_to_matrix_ll(SEXP from, matrix_ll to)
{
   int ncol = INTEGER(GET_DIM(from))[1];
   int nrow = INTEGER(GET_DIM(from))[0];
   to = (matrix_ll) R_alloc(1, sizeof(struct matrix_ll_));
   to->ll = (double **) R_alloc(MAX_MATRIX_LENGTH + 1, sizeof(double *));
   to->llrc = (double **) R_alloc(MAX_MATRIX_LENGTH + 1, sizeof(double *));
   double **cur = NULL;
	for (cur = to->ll; cur < to->ll + ncol; cur++) {
			*cur = (double *) R_alloc(BASES+1, sizeof(double));
			
	}
	for (cur = to->llrc; cur < to->llrc + ncol; cur++) {
   		*cur = (double *) R_alloc(BASES+1, sizeof(double));
			
	}
   
   // IFDEF DEBUG
   //Rprintf("ncol %d\n", ncol);
   //Rprintf("nrow %d\n", nrow);
   //PrintValue(from);
   // ENDIF
   
   if (nrow != BASES) {
      error("Error: nrow of the matrix inside PWMMatrix object != 4");
   }
   to->length = ncol;   
   
   int i = 0;
   int j = 0;
   // We have four rows (index j) and "matrix length" columns (index i) in SEXP from
   // and we have a reversed structure in matrix_ll to with a row foreach matrix element with 4 numbers.
   // That's not true the matrix seems to be already reversed in memory? 
   // Are stored as: first column / then second column etc etc
   for (j = 0; j < ncol; j++) {
      for (i = 0; i < BASES; i++) {
         to->ll[j][i] = REAL(from)[j*BASES+i];
         //Rprintf("girato %d %d %d %g \n", i, j, j*BASES+i, REAL(from)[j*BASES+i]);
      }
   }
   
   for (i = 0; i < ncol; i++) {
      for (j = 0; j < BASES; j++) {
         Rprintf("%g [%d,%d]",to->ll[i][j],i,j);  
      }
      Rprintf("\n");
   }
   
   // are those numbers really ll already?
	// call assign cutoff! here or later? XXX TODO
}

