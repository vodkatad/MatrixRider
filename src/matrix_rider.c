#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "total_affinity.h"

SEXP get_occupancy(SEXP pfm, SEXP cutoff, SEXP sequence) 
{
   // Getting bg.
   SEXP bg = GET_SLOT(pfm, install("bg"));
   // Do we need a check on length? XXX
   double *bg_c = (double *) R_alloc(BASES, sizeof(double));
   for (int i = 0; i < BASES; i++) {
      bg_c[i] = REAL(bg)[i];
   } 
   
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
   
   Chars_holder seq_r = hold_XRaw(sequence);
   int *seq_c = (int *) R_alloc(seq_r.length, sizeof(int)); 
   int seq_length = seq_r.length;
   const char *c = seq_r.ptr;
   for (int i = 0; i < seq_length; i++) {
      seq_c[i] = encode_base(DNAdecode(c[i]));
      //Rprintf("%d: %x(%c)\n", seq_c[i], c[i], DNAdecode(c[i]));
   }
   //Rprintf("\n");
   //Rprintf("seq %s %d\n", seq_r.seq, seq_length);
   double affinity = matrix_little_window_tot(mat_ll, seq_c, seq_length);
   
   SEXP res = PROTECT(allocVector(REALSXP,1));
   REAL(res)[0] = affinity;
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
   
   if (nrow != BASES) {
      return(MATRIX_DIM_ERROR);
   }
   to->length = ncol;   

   // We have four rows (index j) and "matrix length" columns (index i) in SEXP from
   // and we have a reversed structure in matrix_ll to with a row foreach matrix element with 4 numbers.
   // That's not true the matrix seems to be already reversed in memory? 
   // Are stored as: first column / then second column etc etc
   for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < BASES; i++) {
         to->freq[j][i] = (double) (INTEGER(from)[j*BASES+i]);
      }
   }
   
   return(from_counts_to_ll(to));
   // XXX Reason about using ll calculations different from ours (like TFBSTools): accept PWMMatrix?
}


int from_counts_to_ll(matrix_ll m)
{
   double tot = 0.0;
   for (int j = 0; j < m->length; j++) {
        for (int i = 0; i < BASES; i++) {
            int val = (int)m->freq[j][i];
            if (val != m->freq[j][i]) {    //then it's not an integer
                return MATRIX_COUNT_ERROR;
            }
            if (m->freq[j][i] <= EEEPSILON) {
                m->freq[j][i] = 1;
            }
            tot += m->freq[j][i];
        }
        for (int i = 0; i < BASES; i++) {
            m->freq[j][i] = m->freq[j][i] / tot;
        }
        tot = 0.0;
    }
   return OK;
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
   /*for (i = 0; i < m->length; i++) {
      for (j = 0; j < BASES; j++) {
         Rprintf("%g [%d,%d]",m->llrc[i][j],i,j);  
      }
      Rprintf("\n");
   }*/
   
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
        cutoff - the cutoff under which single position affinities will not be considered.
*/
int assign_cutoff_occupancy(matrix_ll m, double cutoff)
{
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
      *error = BACKGROUND_FREQ_ERROR;
        return 0;
    }
    double ratio = n / d;
    if (ratio <= 0) {
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
   error("Wrong argument to getSeqOccupancy, 'sequence' must be based on a restricted alphabet with only 'A','C','G','T' and 'N'");
    return -1;
}

SEXP run_tests() 
{
   int n_failedTests = runAllTests();
   SEXP res = PROTECT(allocVector(INTSXP,1));
   INTEGER(res)[0] = n_failedTests;
   UNPROTECT(1);
   return res;
}

