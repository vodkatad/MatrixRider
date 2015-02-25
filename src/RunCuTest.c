#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "CuTest.h"
#include "total_affinity.h"

/* ration, encodedrc and encode_base are really simple.
Unit tests will be the last one to be added. */
/* convert_PFMMatrix_to_matrix_ll has a SEXP parameter, 
cannot easily test from pure C. Add a standard R test? */

matrix_ll alloc_matrix(int ncol)
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

void test_from_counts_to_ll(CuTest *tc)
{
    int ncol = 3;
    matrix_ll toTest = alloc_matrix(ncol);
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
    from_counts_to_ll(toTest);
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


void test_from_counts_to_ll_errorNotInt(CuTest *tc)
{
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

/* Error should be tested here or in ratio or in both? */
void test_assign_ll(CuTest *tc)
{
    int ncol = 3;
    matrix_ll toTest = alloc_matrix(ncol,BASES);
    toTest->length = 3;

    toTest->freq[0][0] = 0.19047619047619047;
    toTest->freq[0][1] = 0.23809523809523808;
    toTest->freq[0][2] = 0.09523809523809523;
    toTest->freq[0][3] = 0.47619047619047616;
    toTest->freq[1][0] = 0.7692307692307693;
    toTest->freq[1][1] = 0.07692307692307693;
    toTest->freq[1][2] = 0.07692307692307693;
    toTest->freq[1][3] = 0.07692307692307693;
    toTest->freq[2][0] = 0.25;
    toTest->freq[2][1] = 0.25;
    toTest->freq[2][2] = 0.25;
    toTest->freq[2][3] = 0.25;

    double *bg = (double *) R_alloc(BASES, sizeof(double));
    bg[0] = 0.1;
    bg[1] = 0.1;
    bg[2] = 0.5;
    bg[3] = 0.3;

    int res = assign_ll(toTest, bg);

    CuAssertDblEquals(tc, toTest->ll[0][0], 1.9047619047619047, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[0][1], 2.3809523809523805, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[0][2], 0.19047619047619047, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[0][3], 1.5873015873015872, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[1][0], 7.6923076923076925, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[1][1], 0.7692307692307693, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[1][2], 0.15384615384615385, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[1][3], 0.25641025641025644, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[2][0], 2.5, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[2][1], 2.5, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[2][2], 0.5, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[2][3], 0.833333333, EPSILON);
    CuAssertDblEquals(tc, toTest->ll[2][4], NN, EPSILON);
    CuAssertDblEquals(tc, toTest->llrc[2][3], 1.9047619047619047, EPSILON); 
    // Add some more? XXX

    /* We test a ratio error here. */
    bg[3] = 0;
    res = assign_ll(toTest, bg);
    CuAssertIntEquals(tc, res, BACKGROUND_FREQ_ERROR);
}

void test_assign_cutoff_occupancy(CuTest *tc)
{
    int ncol = 3;
    matrix_ll toTest = alloc_matrix(ncol,BASES);
    toTest->length = 3;

    /* We recycle freq values from before for lazyness, 
    for the test's sake it's the same.*/
    toTest->ll[0][0] = 0.19047619047619047;
    toTest->ll[0][1] = 0.23809523809523808;
    toTest->ll[0][2] = 0.09523809523809523;
    toTest->ll[0][3] = 0.47619047619047616;
    toTest->ll[1][0] = 0.7692307692307693;
    toTest->ll[1][1] = 0.07692307692307693;
    toTest->ll[1][2] = 0.07692307692307693;
    toTest->ll[1][3] = 0.07692307692307693;
    toTest->ll[2][0] = 0.25;
    toTest->ll[2][1] = 0.25;
    toTest->ll[2][2] = 0.25;
    toTest->ll[2][3] = 0.25;

    double cutoff = 0.5;
    assign_cutoff_occupancy(toTest, cutoff);
    CuAssertDblEquals(tc, toTest->cutoff, 0.3026137663344012, EPSILON);

    cutoff = 0;
    assign_cutoff_occupancy(toTest, cutoff);
    CuAssertDblEquals(tc, toTest->cutoff, 0, EPSILON);
}

void test_get_affinity(CuTest *tc)
{
    int ncol = 2;
    matrix_ll toTest = alloc_matrix(ncol,BASES);
    toTest->length = 2;

    toTest->ll[0][0] = 0.19047619047619047;
    toTest->ll[0][1] = 0.23809523809523808;
    toTest->ll[0][2] = 0.09523809523809523;
    toTest->ll[0][3] = 0.47619047619047616;
    toTest->ll[1][0] = 0.7692307692307693;
    toTest->ll[1][1] = 0.07692307692307693;
    toTest->ll[1][2] = 0.07692307692307693;
    toTest->ll[1][3] = 0.07692307692307693;
    toTest->llrc[1][3] = 0.19047619047619047;
    toTest->llrc[1][2] = 0.23809523809523808;
    toTest->llrc[1][1] = 0.09523809523809523;
    toTest->llrc[1][0] = 0.47619047619047616;
    toTest->llrc[0][3] = 0.7692307692307693;
    toTest->llrc[0][2] = 0.07692307692307693;
    toTest->llrc[0][1] = 0.07692307692307693;
    toTest->llrc[0][0] = 0.07692307692307693;

    toTest->cutoff = 0;
    int *seq = (int *) R_alloc(ncol, sizeof(int));
    seq[0] = A; /* seq-> AA*/
    seq[1] = A;
    double res = get_affinity(toTest, seq, 0);
    CuAssertDblEquals(tc, res, 0.14652014652014653, EPSILON);
}

void test_matrix_little_window_tot(CuTest *tc)
{
    int ncol = 2;
    matrix_ll toTest = alloc_matrix(ncol,BASES);
    toTest->length = 2;

    toTest->ll[0][0] = 0.19047619047619047;
    toTest->ll[0][1] = 0.23809523809523808;
    toTest->ll[0][2] = 0.09523809523809523;
    toTest->ll[0][3] = 0.47619047619047616;
    toTest->ll[1][0] = 0.7692307692307693;
    toTest->ll[1][1] = 0.07692307692307693;
    toTest->ll[1][2] = 0.07692307692307693;
    toTest->ll[1][3] = 0.07692307692307693;
    toTest->llrc[1][3] = 0.19047619047619047;
    toTest->llrc[1][2] = 0.23809523809523808;
    toTest->llrc[1][1] = 0.09523809523809523;
    toTest->llrc[1][0] = 0.47619047619047616;
    toTest->llrc[0][3] = 0.7692307692307693;
    toTest->llrc[0][2] = 0.07692307692307693;
    toTest->llrc[0][1] = 0.07692307692307693;
    toTest->llrc[0][0] = 0.07692307692307693;

    toTest->cutoff = 0;
    int seq_l = ncol*3;
    int *seq = (int *) R_alloc(seq_l, sizeof(int));
    for (int i = 0; i < seq_l; i++) {
        seq[i]= A; /* seq-> AA*3*/
    }
    double res = matrix_little_window_tot(toTest, seq, seq_l);
    CuAssertDblEquals(tc, res, 0.7326007326007327, EPSILON);

    toTest->cutoff = 10;
    res = matrix_little_window_tot(toTest, seq, seq_l);
    CuAssertDblEquals(tc, res, 0, EPSILON);

    toTest->cutoff = 0.02;
    seq[2] = T;
    seq[3] = T;
    res = matrix_little_window_tot(toTest, seq, 4);
    CuAssertDblEquals(tc, res, 0.29304029304029305, EPSILON);
    /* Maybe add a more complex case with an affinity 
    that will be summed not only from AA.*/

}


CuSuite* MatrixRiderGetSuite();

CuSuite* MatrixRiderGetSuite() {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_from_counts_to_ll);
    SUITE_ADD_TEST(suite, test_from_counts_to_ll_errorNotInt);
    SUITE_ADD_TEST(suite, test_assign_ll);
    SUITE_ADD_TEST(suite, test_assign_cutoff_occupancy);
    SUITE_ADD_TEST(suite, test_get_affinity);
    SUITE_ADD_TEST(suite, test_matrix_little_window_tot);
    return suite;
}

int runAllTests() {
    //CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, MatrixRiderGetSuite());
    CuSuiteRun(suite);
    return suite->failCount;
    //CuSuiteSummary(suite, output);
    //CuSuiteDetails(suite, output);
    //printf("%s\n", output->buffer);
}

