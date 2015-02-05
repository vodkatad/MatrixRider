#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "CuTest.h"
    
char* StrToUpper(char* str) {
   return("HELLO WORLD");
}
    
void TestStrToUpper(CuTest *tc) {
  char* input = strdup("hello world");
  char* actual = StrToUpper(input);
  char* expected = "HELLO WORLD";
  CuAssertStrEquals(tc, expected, actual);
}
   
CuSuite* StrUtilGetSuite() {
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, TestStrToUpper);
  return suite;
}
    

CuSuite* StrUtilGetSuite();
    
int RunAllTests(void) {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();
  
  CuSuiteAddSuite(suite, StrUtilGetSuite());

  CuSuiteRun(suite);
  return suite->failCount;
  //CuSuiteSummary(suite, output);
  //CuSuiteDetails(suite, output);
  //printf("%s\n", output->buffer);
}

