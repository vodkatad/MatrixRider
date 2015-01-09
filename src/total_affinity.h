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
    Constants: Matrix and n. of bases
    
    BASES - n. of bases, used to load correctly log likelood ratios.
    MAX_MATRIX_LENGTH - maximum length of used matrixes.
    MAX_MATRIX_NAME - maximum length of matrixes' names.
    MAX_MATRIXES - maximum number of matrixes that will be loaded.
    MAX_FASTA_ID - maximum length of fasta id that will be considered.
    DEFAULT_WINDOW_SIZE - used if called without -w parameter. 
                          See run_ struct for further info.
    DEFAULT_COUNT - used if called without -c parameter.
    DEFAULT_CUTOFF - used if called without -a parameter.
    DEFAULT_ABS_CUTOFF - used if called without -p parameter.
    DEFAULT_REVCOMP - print revcomp as a default if called without -s.
    DEFAULT_LOG2 - by default it uses likelihoods and not log2-likelihoods 
                   for total affinity modes (-w != 0), if called without -l.
    DEFAULT_SLIDE - by default the sliding windows slides 1 by 1. 
            
    
    ERROR - constant used to return an error status, evaluates true (1).
    OK - constant used to return without error from a function, evaluates false (0).
    A - constant used to represent DNA base.
    C - constant used to represent DNA base.
    G - constant used to represent DNA base.
    T - constant used to represent DNA base.
    N - constant used to assign 0 to N when evaluating matches.
    EPSILON 0.001 - constant used for double equalities evaluation. 
    EEEPSILON 0.0000000001 - smaller constant used for double equalities evaluation. 
    COMMAND_LINE_ERROR 1 - internal error code.
    BACKGROUND_FREQ_ERROR - internal error code.
    MEMORY_ERROR - internal error code.
    MATRIX_COUNT_ERROR - internal error code.
    ILLEGAL_FASTA_CHAR_ERROR - internal error code.
    MATRIX_LIMIT_ERROR - internal error code.
    SHORT_FASTA_WARNING - internal error code.
    MATRIX_FILE_ERROR - internal error code.
     
    STARTING_FASTA - initial number of fasta for which storage will be done, dynamic.
    STARTING_SEQ - initial number of fasta characters for which storage will be done, dynamic.
    DEBUG - define if you want a lot of debugging infos printed, otherwise comment the defining line.
*/
#define BASES 4
#define NBASES 5
#define MAX_MATRIX_LENGTH 35
#define MAX_MATRIX_NAME 50
#define MAX_MATRIXES 2000
#define MAX_FASTA_ID 10
#define DEFAULT_COUNT 1 
#define DEFAULT_LOG2 0 
//#define ERROR 1
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

#define COMMAND_LINE_ERROR 1
#define MEMORY_ERROR 3
#define ILLEGAL_FASTA_CHAR_ERROR 5
#define MATRIX_LIMIT_ERROR 6
#define SHORT_FASTA_WARNING 7
#define MATRIX_FILE_ERROR 8
#define BG_FILE_ERROR 9
#define CUTOFF_FILE_ERROR 10
#define FREQ_ZERO_ERROR 11


#define STARTING_FASTA 50
#define STARTING_SEQ 17
/*are we wasting ram? */
/*#define DEBUG
#define NDEBUG */
 
#ifndef __TOTAL_AFFINITY_H__
#define __TOTAL_AFFINITY_H__

#define strong_assert(ASSERTION) ({\
	if ((ASSERTION)==0){\
        	fprintf (stderr,\
			"Assertion failed: " # ASSERTION\
			", function %s, file %s, line %u.\n",\
			__func__, __FILE__, __LINE__);\
			exit(1);\
	}\
})

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

/*struct error_ {
   int CODE
}*/

/*
    Struct: fasta is a pointer to a struct used to store info on fastas.  
*/
typedef struct fas *fasta;

struct fas {
    int *seq;
    char *id;
    double background[BASES];
    int length;
    int n_portions;
};

/* 
    Struct: run stores informations on this run.
*/
typedef struct run_ *run;

struct run_ {
    int counts; /* if 0 we have fractions, if 1 pseudocounts/counts */
    int log2; /* if 1 we use log2likelihood even for total affinities mode, if 0 (default) only for ss mode. */
    double (*get_affinities_pointer)(matrix_ll, int *, int); /* pointer to the function that get affinities */
    FILE *f_matrixes;
    FILE *f_fasta;
    FILE *f_background;
    int end_mfasta;
    int n_fasta;
    char *next_id;
    matrix_ll *matrixes;
    int n_matrixes;
    int error;
    int normalize_on_seq_len;
};

int convert_PFMMatrix_to_matrix_ll(SEXP from, matrix_ll *toptr);
int assign_ll(matrix_ll m, double *bg);
int assign_cutoff_occupancy(matrix_ll m, double cutoff);
int encoded_rc(int n);
double get_affinity(matrix_ll m, int *s, int start)
double ratio(double n, double d, int *error);
double matrix_little_window_tot(matrix_ll m, char *seq, int seq_length);

void free_fasta(run info);
fasta load_next_fasta(run info, double *bg);
void add_count_bg(fasta f, char c, run info); 
int encode_base(char c);
int encoded_rc(int n);
void assign_zero_bg(fasta f); 
void get_id(FILE *f, char **s, run info);
void get_seq(char *seq, int offset, int len, char *match);
void get_rc(char *seq, int offset, int len, char *match); 

void check_error(run info);
int is_atofable(char *s);
int is_atoable(char *s);
int open_file(char *argv[], int position, FILE **to_open);
void parse_command_line(char *argv[], int argc, run info);
void check_parameters(run info);
char *get_next_token(char *argv[], int argc, int next_token);
void tidy(run info);
int set_beginning_info(run *info);

void print_fasta(fasta f);
void print_matrixes(run info);
#endif