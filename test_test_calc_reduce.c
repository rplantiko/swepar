#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "testdata.h"
#include "test_test_calc_reduce.h"

bool test(char *filename) {
  
    FILE *testdata = fopen( filename,"rb" );
    if (!testdata) {
      char err[255];
      sprintf( err, "Can't open %s", filename);
      fail( err );
      }
    
    testHeader header;
    bool ok = true;
    while (hasMoreTestSets(testdata, &header)) {
      result r = executeTestSet(testdata);
      printf("%s: ", header.description );
      if (r.ok) printf("%d OK   ", r.ok);
      if (r.not_ok) printf("%d NOT OK", r.not_ok);
      printf("\n");
      ok = ok && (r.not_ok == 0);
      }
      
    fclose( testdata );
    
    return ok;

  }
  
static int hasMoreTestSets(FILE* testdata, testHeader* header) {
  fread( header, sizeof(testHeader), 1, testdata );
  return !feof(testdata);
  }  
  
 
static result executeTestSet(FILE* testdata) {

  result r = {0,0};  

  swissCalcHeader scHeader;  
  fread( &scHeader, sizeof(scHeader), 1, testdata );  
    
  swissCalcData scData[scHeader.numberOfRecords];  
  int numberOfRecordsRead = fread( &scData, sizeof(swissCalcData), scHeader.numberOfRecords, testdata );  
  if (numberOfRecordsRead != scHeader.numberOfRecords) {  
    printf( "Wrong number %d of records read - expected %d\n", numberOfRecordsRead, scHeader.numberOfRecords );  
    exit( EXIT_FAILURE );  
    }  
  
  double precisions[6];  
  computePrecisions( precisions, scHeader.precisions );  
      
  for (int i=0;i<scHeader.numberOfRecords;i++) {  
    bool test_ok = calcsAsExpected(  
      scData[i].jd, 
      scData[i].planet, 
      scData[i].flags, 
      scData[i].flags_ret, 
      scData[i].result, 
      precisions, 
      scData[i].msg
      );  
  
    if (!test_ok) {  
      printf( "\nDifferences from test case %d", i );  
      printf( "\nJuldat %15.10f, planet %d, iflags %d\n",scData[i].jd, scData[i].planet, scData[i].flags);  
      }  
        
    if (test_ok) r.ok++; 
    else r.not_ok++;
      
    }  
  
  return r;  
    
  }

static bool calcsAsExpected(
  const double juldate,
  const int planet,
  const int flags,
  const int flags_ret,
  const double expectedResult[],
  const double precisions[],
  const char* expectedMessage
  ) {

    double actualResult[6];
    char actualMessage[255];
     
    int return_flags = swe_calc(
      juldate,
      planet,
      flags,
      actualResult,
      actualMessage );
      
  bool ok =
    diff_check_degree( "xx[0]", actualResult[0],expectedResult[0],precisions[0]);
    
  for (int i=1;i<6;i++) {
    char var[15];
    sprintf(var,"xx[%d]",i);
    ok = ok && diff_check(var,actualResult[i],expectedResult[i],precisions[i]);
    }  
   ok = ok 
     && diff_check_char( "serr", actualMessage, expectedMessage)   
     && diff_check_int( "iflags", return_flags, flags_ret );               

   return ok;

  }

static void computePrecisions( double precisions[6], int precisionsExp[6] ) {
  for (int i=0;i<6;i++) {
    precisions[i] = defaultPrecision;
    if (precisionsExp[i] != 0) precisions[i] *= pow(10,precisionsExp[i]);
    }  
  }
  
