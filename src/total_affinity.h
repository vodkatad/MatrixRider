/* total_affinity.h
//
// Copyright (C) 2008 Elena Grassi <grassi.e@gmail.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but 
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
*/

/*
    Constants: TODO finish commenting
    
    BASES - n. of bases, used to load correctly log likelood ratios.
    ERROR - constant used to return an error status, evaluates true (1).
    OK - constant used to return without error from a function, evaluates false (0).
    A - constant used to represent DNA base.
    C - constant used to represent DNA base.
    G - constant used to represent DNA base.
    T - constant used to represent DNA base.
    N - constant used to assign 0 to N when evaluating matches.
    EPSILON 0.001 - constant used for double equalities evaluation. 
    EEEPSILON 0.0000000001 - smaller constant used for double equalities evaluation. 

*/
#define BASES 4
#define NBASES 5
#define OK 0
#define A 0
#define C 1
#define G 2
#define T 3
#define N 4
#define NN 0
#define EPSILON 0.001 
#define EEEPSILON 0.0000000001 

#define MATRIX_DIM_ERROR 1
#define MATRIX_COUNT_ERROR 2
#define BACKGROUND_FREQ_ERROR 3

/*#define DEBUG
#define NDEBUG */
 
#ifndef __TOTAL_AFFINITY_H__
#define __TOTAL_AFFINITY_H__

/*
    Struct: matrix_ll is a pointer to a struct used to store info on matrixes. 
*/
typedef struct matrix_ll_ *matrix_ll;

struct matrix_ll_ {
	double **ll;
	double **llrc;
   double **freq;
   double cutoff;
	int length;
	//char *name;
};


/* R entry points: */
SEXP get_occupancy(SEXP pfm, SEXP cutoff, SEXP sequence);
SEXP run_tests();

/* Internal C functions: */
int convert_PFMMatrix_to_matrix_ll(SEXP from, matrix_ll *toptr);
int from_counts_to_ll(matrix_ll m);
int assign_ll(matrix_ll m, double *bg);
int assign_cutoff_occupancy(matrix_ll m, double cutoff);
int encoded_rc(int n);
int encode_base(const char c);
double get_affinity(matrix_ll m, int *s, int start);
double ratio(double n, double d, int *error);
double matrix_little_window_tot(matrix_ll m, int *seq, int seq_length);

#endif