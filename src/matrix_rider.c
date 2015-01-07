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

void convert_PWMMatrix_to_matrix_ll(SEXP from, matrix_ll to)
{
   int ncol = INTEGER(GET_DIM(from))[0];
   int nrow = INTEGER(GET_DIM(from))[1];
   if (nrow != BASES) {
      // Error and exit badly XXX TODO
   }
   to->length = ncol;   
   
   int i = 0;
   int j = 0;
   // We have four rows (index j) and "matrix length" columns (index i) in SEXP from
   // and we have a reversed structure in matrix_ll to with a row foreach matrix element with 4 numbers.
   for (i = 0; i < ncol; i++) {
      for (j = 0; j < BASES; j++)
			to->ll[i][j] = REAL(mat)[j*nrow+i];
   }
   // are those numbers really ll already?
   // do we need to use R_alloc and not calloc/mallocs?
   // Do we need the matrix name here? probably not
	// call assign cutoff! here or later? XXX TODO
}

/*
    Function: alloc_matrixes
    
    Creates space (calloc) for all the matrixes that will be loaded.
    
    Parameters:
        m - matrix_ll ** where the matrixes will be stored.
    
    Returns:
        ERROR (1) if any calloc call failed, OK (0) otherwise.
*/
int alloc_matrixes(matrix_ll ** m)
{
   *m = (matrix_ll *) calloc(MAX_MATRIXES + 1, sizeof(matrix_ll));
	if (!(*m)) {
		free(*m);
		return ERROR;
	}
	int i = 0;
	for (i = 0; i < MAX_MATRIXES; i++) {
		(*m)[i] = (matrix_ll) calloc(1, sizeof(struct matrix_ll_));
		if (!((*m)[i])) {
			free_matrixes(*m, i);
			return ERROR;
		}
		(*m)[i]->ll =
		    (double **)calloc(MAX_MATRIX_LENGTH + 1, sizeof(double *));
		if (!((*m)[i]->ll)) {
			free_matrixes(*m, i);
			return ERROR;
		}
		(*m)[i]->llrc =
		    (double **)calloc(MAX_MATRIX_LENGTH + 1, sizeof(double *));
		if (!((*m)[i]->llrc)) {
			free_matrixes(*m, i);
			return ERROR;
		}
		(*m)[i]->freq =
		    (double **)calloc(MAX_MATRIX_LENGTH + 1, sizeof(double *));
		if (!((*m)[i]->freq)) {
			free_matrixes(*m, i);
			return ERROR;
		}
		double **cur = NULL;
		for (cur = (*m)[i]->ll; cur < (*m)[i]->ll + MAX_MATRIX_LENGTH;
		     cur++) {
			*cur = (double *)calloc(BASES+1, sizeof(double));
			if (!(*cur)) {
				free_matrixes(*m, i);
				return ERROR;
			}
		}
		for (cur = (*m)[i]->llrc; cur < (*m)[i]->llrc + MAX_MATRIX_LENGTH;
		     cur++) {
			*cur = (double *)calloc(BASES+1, sizeof(double));
			if (!(*cur)) {
				free_matrixes(*m, i);
				return ERROR;
			}
		}
		for (cur = (*m)[i]->freq;
		     cur < (*m)[i]->freq + MAX_MATRIX_LENGTH; cur++) {
			*cur = (double *)calloc(BASES, sizeof(double));
			if (!(*cur)) {
				free_matrixes(*m, i);
				return ERROR;
			}
		}
		(*m)[i]->name = (char *)calloc(MAX_MATRIX_NAME, sizeof(char));
		if (!((*m)[i]->name)) {
			free_matrixes(*m, i);
			return ERROR;
		}
	}
	return OK;
}

/*
    Function: free_matrixes
    
        Frees the space used by all the loaded matrixes.
    
    Parameters:
    
        m - matrix_ll * with the matrixes to be freed.
        loaded - number of matrixes to be freed.
*/
void free_matrixes(matrix_ll * m, int loaded)
{
	int i = 0;
	for (i = 0; i < loaded; i++) {
		free(m[i]->name);
		double **cur = NULL;
		for (cur = m[i]->ll; cur < m[i]->ll + MAX_MATRIX_LENGTH; cur++)
			free(*cur);
		free(m[i]->ll);
		for (cur = m[i]->llrc; cur < m[i]->llrc + MAX_MATRIX_LENGTH; cur++)
			free(*cur);
		free(m[i]->llrc);
		for (cur = m[i]->freq; cur < m[i]->freq + MAX_MATRIX_LENGTH; cur++)
			free(*cur);
		free(m[i]->freq);
		free(m[i]);
	}
	free(m);
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
        info - run with cutoff informations.
*/
void assign_cutoff_occupancy(matrix_ll m, run info)
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

	/* If a file with cutoffs is given we will use them, otherwise we have default
	   or given parameters equals for every matrix. */
	if (info->f_cutoffs) {
		info->cutoff = m->frac_cutoff;
		info->abs_cutoff = m->abs_cutoff;
	}
	if (info->cutoff == 0) 
		m->cutoff = 0; //a^0 = 1 but with 0 we want the same results that with total affinity
	else
		m->cutoff = pow(max_tot, info->cutoff); 	
}


/*
    Function: matrix_little_window_tot
    
    Computes affinity between a matrix and a string on all viable positions and
    returns the sum of all affinities.
    
    Parameters:
    
        m   - matrix_ll matrix that will be used.
        f   - fasta seq that will be used.
        begin - starting base to compute affinity on f.
        end - TODO XXX
    Returns:
       
       double with total affinity between the given matrix and fasta.
*/
double matrix_little_window_tot(matrix_ll m, fasta f, int begin, int end)
{
   int offset = 0;
   char *seq = f->seq;
   int l = f->length;
	if (begin > 0){
		if (begin > l){
			fprintf(stderr,	"ERROR: invalid begin in matrix_little_window_tot()\n");
			exit(1);
		}
		offset = begin;
	}
	if(end > 0){
		if(end > l){
			fprintf(stderr,	"ERROR: invalid end in matrix_little_window_tot()\n");
			exit(1);
		}
	}else{
		end = l;
	}

	double tot = 0;
	(*tot_match_noN) = 0;
	while (offset <= end - m->length) {
		int foundN = 0;
		//was get_affinities_nonLog
		double match = info->get_affinity(m, seq, offset);
		if (!info->occupancy || match >= m->cutoff) {
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