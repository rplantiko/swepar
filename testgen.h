#include <stdio.h>
#include <string.h>
#include "testdata.h"
#include "swephexp.h"

// --- An ephemeris test data generator, configurable with a fixture file

  // Test type  
static const int SWISS_CALC_TEST = 1;

// --- Ranges (collection of single values and intervals)
typedef struct {
  AS_BOOL between;  // Interval or single value
  int  low;
  int  high;
  } range;

#define MAXRANGE 50
typedef struct {
  int size;
  range line[MAXRANGE];
  } rangetab;

typedef struct {
  rangetab *r;
  int i;  // line count
  int j;  // interval count
  } range_iterator;
static int get_next( range_iterator *iter);

// Function how to generate a series of dates (random or fixed-increment)  
typedef void(*datefun)(double dates[], int n, double jd_from, double jd_to, double incr);  


// --- Test data, as read from a fixture
#define MAX_FIXNUM 1000  
typedef struct {
  char description[100];
  int planet;
  int iflag;
  int precisions[6];
  double jd_from, jd_to, incr;
  int n_dates;
  datefun method;
  rangetab planets;
  rangetab flags;
  } testFixture;  
  
static int generateTestsetFromFixture(char* fixtureFile,char* dataFile);
static long int doGeneralHeader( FILE* out, int type, char* description );
static void doSwissCalcHeader( FILE* out, long int pos, testFixture* test, int numberOfRecords);
static int doSwissCalcTestData( FILE* out, testFixture* test);
static int parseFixtureFile(char* fixtureFile, testFixture*, int*);
static void clear( testFixture* data );
static int parseNameValuePair( char* name, char* value, testFixture* data);
static int generateTest( FILE* out, testFixture* data);
static int parseDates( char* value, testFixture* data);
static int parsePrecisions( char* value, int* precisions );
static int parseRange( char* value, rangetab *r);
static inline int equals( char* act, char* exp) { return ! strcmp( act, exp ); }
static inline char* trim(char *s) { for (char* p = s + strlen(s)-1; *p == ' '; p-- ) *p='\0'; return s; } 
static FILE* openFileBinary( const char* file );
static void doInitializations();

static char* int_to_binary(char* b, int x);
