#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "CuTest.h"
#include "total_affinity.h"

/* ration, encodedrc and encode_base are really simple. Unit tests will be the last one to be added. */
/* convert_PFMMatrix_to_matrix_ll has a SEXP parameter, cannot easily test from pure C. Add a standard R test? */

matrix_ll alloc_matrix(int ncol, int nrow)
{
   matrix_ll to =(matrix_ll) R_alloc(1, sizeof(struct matrix_ll_));
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
   return(to);
}

void test_from_counts_to_ll(CuTest *tc) {
   int ncol = 5;
   matrix_ll toTest = alloc_matrix(ncol,BASES);
   toTest->length = ncol;   

   /*double a[5][BASES] = {  
      {4, 5, 2, 10} ,
      {10, 1, 1, 1} ,
      {5, 5, 5, 5},
      {1, 100, 2, 6},
      {10, 2, 2, 10}
   };
   
   matrix_ll to =(matrix_ll) R_alloc(1, sizeof(struct matrix_ll_));
   to->ll = (double **) R_alloc(ncol, sizeof(double *));
   to->llrc = (double **) R_alloc(ncol, sizeof(double *));
   double **cur = NULL;
   for (cur = to->ll; cur < to->ll + ncol; cur++) {
      	*cur = (double *) R_alloc(NBASES, sizeof(double));		
	}
	for (cur = to->llrc; cur < to->llrc + ncol; cur++) {
   		*cur = (double *) R_alloc(NBASES, sizeof(double));		
	}
   to->freq = a;*/
   toTest->freq[0][0] = 4;
   toTest->freq[0][1] = 5;
   toTest->freq[0][2] = 2;
   toTest->freq[0][3] = 10;
   toTest->freq[1][0] = 10;
   toTest->freq[1][1] = 1;
   toTest->freq[1][2] = 1;
   toTest->freq[1][3] = 1;
   toTest->freq[2][0] = 5;
   toTest->freq[2][1] = 5;
   toTest->freq[2][2] = 5;
   toTest->freq[2][3] = 5;
   int res = from_counts_to_ll(toTest);
   CuAssertDblEquals(tc, toTest->freq[0][0], 0.19047619047619047, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[0][1], 0.23809523809523808, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[0][2], 0.09523809523809523, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[0][3], 0.47619047619047616, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[1][0], 0.7692307692307693, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[1][1], 0.07692307692307693, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[1][2], 0.07692307692307693, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[1][3], 0.07692307692307693, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[2][0], 0.25, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[2][1], 0.25, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[2][2], 0.25, EPSILON);
   CuAssertDblEquals(tc, toTest->freq[2][3], 0.25, EPSILON);
}


void test_from_counts_to_ll_errorNotInt(CuTest *tc) {
   int ncol = 1;
   matrix_ll toTest = alloc_matrix(ncol,BASES);
   toTest->length = ncol;   
   toTest->freq[0][0] = 0.1;
   toTest->freq[0][1] = 0;
   toTest->freq[0][2] = 0;
   toTest->freq[0][3] = 0;
   int res = from_counts_to_ll(toTest);
   CuAssertIntEquals(tc, res, MATRIX_COUNT_ERROR);
}
    
CuSuite* MatrixRiderGetSuite();   

CuSuite* MatrixRiderGetSuite() {
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_from_counts_to_ll);
  SUITE_ADD_TEST(suite, test_from_counts_to_ll_errorNotInt);
  return suite;
}

int runAllTests();
 
int runAllTests() {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();
  CuSuiteAddSuite(suite, MatrixRiderGetSuite());
  CuSuiteRun(suite);
  return suite->failCount;
  //CuSuiteSummary(suite, output);
  //CuSuiteDetails(suite, output);
  //printf("%s\n", output->buffer);
}

