#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "total_affinity.h"

SEXP get_occupancy(SEXP pfm, SEXP cutoff, SEXP sequence);

SEXP get_occupancy(SEXP pfm, SEXP cutoff, SEXP sequence) 
{
   //SEXP matrix_slot = mkChar("profileMatrix\0");
   //PrintValue(matrix_slot);
   
   // Getting bg.
   SEXP bg = GET_SLOT(pfm, install("bg"));
   // XXX check on length
   double *bg_c = (double *) R_alloc(BASES, sizeof(double));
   for (int i = 0; i < BASES; i++) {
      bg_c[i] = REAL(bg)[i];
   } 
   
   // SEXP dim = GET_DIM(pwm); // --> this is NULL for my pwm. Why?
   SEXP mat = GET_SLOT(pfm, install("profileMatrix"));
   matrix_ll mat_ll = NULL;
   if (convert_PFMMatrix_to_matrix_ll(mat, &mat_ll)) {
      error("Error while converting PFMMatrix to PWM: not integer counts or wrong dimensions");
   }
   if (assign_ll(mat_ll, bg_c)) {
      error("Error while assigning (log)-likelihoods: 0 bg?");
   }
   double cutoff_c = REAL(cutoff)[0];
   if (assign_cutoff_occupancy(mat_ll, cutoff_c)) {
      error("Error while assigning cutoff to matrix");
   }
   Rprintf("cutoff %g\n", mat_ll->cutoff);

   // data@femto:~$ cp /home/data/R/x86_64-pc-linux-gnu-library/3.1/S4Vectors/include/S4Vectors_defines.h /home/data/R/x86_64-pc-linux-gnu-library/3.1/Biostrings/include/
   // WTF?
   /*
   typedef struct chars_holder {
      const char *seq;
   	int length;
   } Chars_holder;
   */

   Chars_holder seq_r = hold_XRaw(sequence);
   //PrintValue(sequence);
   int *seq_c = (int *) R_alloc(seq_r.length, sizeof(int)); // same with sample association
   int seq_length = seq_r.length;
   const char *c = seq_r.seq;
   for (int i = 0; i < seq_length; i++) {
      seq_c[i] = encode_base(DNAdecode(c[i]));
      //Rprintf("%x(%c) ", c[i], DNAdecode(c[i]));
   }
   //Rprintf("\n");
   //Rprintf("seq %s %d\n", seq_r.seq, seq_length);
   double affinity = matrix_little_window_tot(mat_ll, seq_c, seq_length);
   Rprintf("affinity %g\n", affinity);
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
   p[2] = (int)(INTEGER(mat)[i*4+j]); // access element at row i and col j (0 based)
   UNPROTECT(1);
   return res;
}

int convert_PFMMatrix_to_matrix_ll(SEXP from, matrix_ll *toptr)
{
   int ncol = INTEGER(GET_DIM(from))[1];
   int nrow = INTEGER(GET_DIM(from))[0];
   // XXX TOdo R_alloc error checking
   *toptr = (matrix_ll) R_alloc(1, sizeof(struct matrix_ll_));
   matrix_ll to = *toptr;
   to->ll = (double **) R_alloc(ncol, sizeof(double *));
   to->llrc = (double **) R_alloc(ncol, sizeof(double *));
   to->freq = (double **) R_alloc(ncol, sizeof(double *));
   double **cur = NULL;
	for (cur = to->ll; cur < to->ll + ncol; cur++) {
			*cur = (double *) R_alloc(NBASES, sizeof(double));
			
	}
	for (cur = to->llrc; cur < to->llrc + ncol; cur++) {
   		*cur = (double *) R_alloc(NBASES, sizeof(double));		
	}
   for (cur = to->freq; cur < to->freq + ncol; cur++) {
      	*cur = (double *) R_alloc(NBASES, sizeof(double));
			
	}
   // IFDEF DEBUG
   //Rprintf("ncol %d\n", ncol);
   //Rprintf("nrow %d\n", nrow);
   //PrintValue(from);
   // ENDIF
   
   if (nrow != BASES) {
      return(MATRIX_DIM_ERROR);
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
         to->freq[j][i] = (double) (INTEGER(from)[j*BASES+i]);
         //Rprintf("girato %d %d %d %g \n", i, j, j*BASES+i, REAL(from)[j*BASES+i]);
      }
   }
   
   // maybe factorize as it was in a separate function?
   double tot = 0.0;
	for (j = 0; j < to->length; j++) {
		for (i = 0; i < BASES; i++) {
			int val = (int)to->freq[j][i];
			if (val != to->freq[j][i]) {	//then it's not an integer
				return MATRIX_COUNT_ERROR;
			}
			if (to->freq[j][i] <= EEEPSILON) {
				to->freq[j][i] = 1;
			}
			tot += to->freq[j][i];
		}
		
		if (!tot) {
			return MATRIX_COUNT_ERROR;
		}
		for (i = 0; i < BASES; i++) {
			to->freq[j][i] = to->freq[j][i] / tot;
		}
		tot = 0.0;
	}
   
   /*for (i = 0; i < ncol; i++) {
      for (j = 0; j < BASES; j++) {
         Rprintf("%g [%d,%d]",to->freq[i][j],i,j);  
      }
      Rprintf("\n");
   }*/
   
   return OK;
   // are those numbers really ll already? No. Similar but not the same. Will start with them then decide:
   // p.20 (or 30) of the manual and a gnumeric with counts for a sample matrix.
   // could be ok to develop a toPWM that does what we do (pseudocounts only on 0 and log2(pwm/bg)).
   
	// call assign cutoff! here or later? XXX TODO
   
   // stupid! We need to start from the pfm, get the bg and obtain the ll ourselves!
}

int assign_ll(matrix_ll m, double *bg)
{
   int j = 0;
	int i = 0;
   int error = OK;
	for (; j < m->length; j++) {
		for (i = 0; i < BASES; i++) {
			m->ll[j][i] = ratio(m->freq[j][i], bg[i], &error);
		}
      m->ll[j][N] = NN;
	} 
   for (j= 0; j < m->length; j++) {
   	for (i = 0; i < BASES; i++) {
			m->llrc[j][i] = m->ll[(m->length) - j - 1][encoded_rc(i)];
		}
		m->llrc[j][N] = NN;
	}
   
   // DEBUG
   for (i = 0; i < m->length; i++) {
      for (j = 0; j < BASES; j++) {
         Rprintf("%g [%d,%d]",m->ll[i][j],i,j);  
      }
      Rprintf("\n");
   }
   return error;   
}


/*
    Function: assign_cutoff_occupancy

    Assign a cutoff to a matrix_ll using the fractional cutoff 
    in the given run - uses likelihoods and not log likelihoods.
    The cutoff calculations are:   
    sum (log P/Pbg) > 8/10 sum(log Pmax/Pbg)
    2^
    prod (P/Pbg) > prod (Pmax/Pbg)^8/10
    
    Parameters:
        m - matrix_ll with ll loaded and missing cutoff.
        cutoff - TODO
*/
int assign_cutoff_occupancy(matrix_ll m, double cutoff)
{
   Rprintf("len %d\n", m->length);
   int j = 0;
	int i = 0;
	double max_tot = 1;
	double max;
	for (; j < m->length; j++) {
		max = m->ll[j][0];
		for (i = 1; i < BASES; i++) {
			if (m->ll[j][i] > max) {
				max = m->ll[j][i];
			}
		}
		max_tot *= max;
	}

	if (cutoff == 0) 
		m->cutoff = 0; //a^0 = 1 but with 0 we want the same results that with total affinity
	else
		m->cutoff = pow(max_tot, cutoff); 	
      
   return OK;
}

/*  Function: ratio
    
    Returns the the ratio between numerator and denominator.
    Sets error to 1 if the ratio is non-positive or denominator is ~0.
    
    Parameters:
        n - double which will be used as numerator.
        d - double which will be used as denominator.
        error - pointer to set error flags.
        
    Return:
        ratio between n and d (will be 0 if an error occurred).
*/

double ratio(double n, double d, int *error)
{
   // FIXME XXX is this the right error code?
   if (d < EPSILON) {
      Rprintf("here\n");
      *error = BACKGROUND_FREQ_ERROR;
		return 0;
	}
	double ratio = n / d;
	if (ratio <= 0) {
      Rprintf("here2\n");
		*error = BACKGROUND_FREQ_ERROR;
		return 0;
	}	
	return ratio;
}

int encoded_rc(int n)
{
   switch (n) {
		case A:
			return T;
		case C:
			return G;
		case G:
			return C;
		case T:
			return A;
		case N:
			return N;
	}
	return -1;
}

double matrix_little_window_tot(matrix_ll m, int *seq, int seq_length)
{
   int offset = 0;
	double tot = 0;
	while (offset <= seq_length - m->length) {
		double match = get_affinity(m, seq, offset);
		if (match >= m->cutoff) {
			tot += match;
		}
		offset++;
	}
	return tot;
}

double get_affinity(matrix_ll m, int *s, int start)
{
   /* results[0] is straight total, results[1] revcomp */
   int i = 0;
	double results[2];
	results[0] = 1;
	results[1] = 1;
	while (i < m->length) {
		results[0] *= m->ll[i][s[start]];
		results[1] *= m->llrc[i][s[start]]; 
		i++;
		start++;
	}
	return (results[0] > results[1]) ? results[0] : results[1];
}

int encode_base(const char c)
{
   switch (c) {
		case 'A':
			return A;
		case 'C':
			return C;
		case 'G':
			return G;
		case 'T':
			return T;
		case 'N':
			return N;
	}
	return -1;
}
