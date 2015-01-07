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

void get_matrixes(run info)
{
   if (alloc_matrixes(&(info->matrixes)))
		info->error = MEMORY_ERROR;
	matrix_ll *matrixes = info->matrixes;
	char name[MAX_MATRIX_NAME];
	char old_name[MAX_MATRIX_NAME];
	double acgt[BASES];
	strcpy(old_name, "");
	strcpy(name, "");
	int length = 0;
	int ongoing = 1;
	int i = 0;
	while (ongoing) {
		if (fscanf(info->f_matrixes, "%s\t%*d\t%lf\t%lf\t%lf\t%lf\n",
			   name, &acgt[A], &acgt[C], &acgt[G], &acgt[T]) != 5) {
#ifdef DEBUG
			printf("You fool fscanf!\n");
#endif				/* DEBUG */
			if (fgetc(info->f_matrixes) != EOF)
				info->error = MATRIX_FILE_ERROR;
			ongoing = 0;
			break;	/* XXX BAD */
		}
		if (strcmp(name, old_name)) {	/* we have a new matrix */
#ifdef DEBUG
			printf("Last was %d, now adding %s\n", length, name);
#endif				/* DEBUG */
			if (info->n_matrixes != 0) {
				(*matrixes)->length = length;	/* assign length to last matrix */
				get_fractions_from_pcounts(*matrixes, info);
				check_error(info);
				matrixes++;	/* next matrix to be loaded */
			}
			info->n_matrixes++;
			length = 0;
			strcpy((*matrixes)->name, name);
		}
		for (i = 0; i < BASES; i++)
			(*matrixes)->freq[length][i] = acgt[i];
		length++;
		if (info->n_matrixes >= MAX_MATRIXES
		    || length >= MAX_MATRIX_LENGTH) {
			ongoing = 0;
			info->error = MATRIX_LIMIT_ERROR;
		}
		strcpy(old_name, name);
	}
	(*matrixes)->length = length;	/* assign length to last matrix */
	get_fractions_from_pcounts(*matrixes, info);
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
    Function: matrix_little_window_tot
    
    Computes affinity between a matrix and a string on all viable positions and
    returns the sum of all affinities.
    
    Parameters:
    
        m   - matrix_ll matrix that will be used.
        seq - sequence versus which matrix will be compared.
   	s_length - the

    
    Returns:
       
       double with total affinity between the given matrix and fasta,
       normalized upon the number of matrix matches != 0 (i.e. without an N).
*/
double matrix_little_window_tot(matrix_ll m, fasta f, int begin, run info)
{
	int offset = 0;
    int *seq = f->seq;
	double tot = 0;
	while (offset <= f->length - m->length) {
		tot += info->get_affinities_pointer(m, seq, offset);
		offset++;
	}
	return tot;
}

double get_affinities_nonLog(matrix_ll m, int *s, int start)
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

double get_affinities_log(matrix_ll m, int *s, int start)
{
	/* results[0] is straight total, results[1] revcomp */
	int i = 0;
	double results[2];
	results[0] = 0;
	results[1] = 0;
	while (i < m->length) {
		results[0] += m->ll[i][s[start]];
		results[1] += m->llrc[i][s[start]]; 
		i++;
		start++;
	}
	return (results[0] > results[1]) ? results[0] : results[1];
}

/*  Function: log2_ratio
    
    Returns the log2ratio or only the ratio (depending on info->log2),
    sets error to 1 if the ratio is non-positive.
    Computes log2ratio if info->log2 evaluates true, ratio otherwise.
    
    Parameters:
        n - double which will be used as numerator.
        d - double which will be used as denominator.
        info - run used to set error flags and determine if log has to be applied.
        
    Return:
        double which is the log2-ratio or simple ratio between n and d.
        (will be 0 if an error occurred).
*/

double log2_ratio(double n, double d, run info)
{
   if (d < EPSILON) {
   	info->error = BACKGROUND_FREQ_ERROR;
		return 0;
	}
	double ratio = n / d;
	if (ratio <= 0) {
		info->error = BACKGROUND_FREQ_ERROR;
		return 0;
	}
	
	if (info->log2)
		return log(ratio) / log(2);
	else 
		return ratio;
}