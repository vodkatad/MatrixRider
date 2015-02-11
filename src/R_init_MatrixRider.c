#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "total_affinity.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
   CALLMETHOD_DEF(get_occupancy, 3),
   CALLMETHOD_DEF(run_tests, 0),
   {NULL, NULL, 0}
};

void R_init_MatrixRider(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   return;
}
