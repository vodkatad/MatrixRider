#include <R.h>
#include <Rinternals.h>

SEXP get_occupancy(SEXP);
int factorial(int);
int rec_fact(int);

SEXP get_occupancy(SEXP pwm) 
{
   int local_pwm = Rf_asInteger(pwm);
   
   //SEXP res = PROTECT(Rf_allocSExp(INTSXP));
   //*(INTEGER(res)) = factorial(local_pwm);
   
   SEXP res = PROTECT(allocVector(INTSXP,1));
   int * p = INTEGER(res);
   p[0] = factorial(local_pwm);
   
   UNPROTECT(1);
   return(res);
}

int factorial(int x)
{
   if (x < 1) {
      return 0;
   } else {
      return rec_fact(x);
   }
}

int rec_fact(int x)
{
   if (x == 1) {
      return 1;
   } else {
      return x*rec_fact(x-1);
   }
}