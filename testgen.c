#include "testgen.h"

// --- An ephemeris test data generator, configurable with a fixture file

int main( int argc, char** argv) {

  doInitializations();

  if (argc == 3)
    generateTestsetFromFixture(argv[1],argv[2]);
  else {
    printf("Wrong arg number\nSyntax: testgen [fixture file] [test data file]");
    exit(EXIT_FAILURE);
    }

}

static int generateTestsetFromFixture(char* fixtureFile, char* dataFile) {

  FILE* data = openFileBinary(dataFile);

  testFixture tests[MAX_FIXNUM];
  int size;  
  if (parseFixtureFile(fixtureFile, tests, &size) == ERR) {
    printf( "Error while reading fixture file\n" );
    return ERR;
    }
    
  for (int i=0; i<size;i++) {
    if (generateTest( data, & tests[i] )==ERR) {
      printf("Error in generation of test %d\n",i+1);
      return ERR;
      }
    }

  fclose(data);
  return OK;
  }

static int generateTest( FILE* out, testFixture* test) {
    
  long int afterHeaderPos = doGeneralHeader( out, SWISS_CALC_TEST, test->description );

  doSwissCalcHeader( out, afterHeaderPos, test, 0 );
  
  int numberOfRecords = doSwissCalcTestData( out, test );
  
// Correct SwissCalcHeader block - set number of records  
  doSwissCalcHeader( out, afterHeaderPos, test, numberOfRecords );
  
// Goto end of file afterwards  
  fseek( out, 0, SEEK_END );
    
  printf("%s : %d records generated\n", test->description, numberOfRecords );

  return OK;

  }

  
static long int doGeneralHeader( FILE* out, int type, char* description ) {
  testHeader header;
  header.type = type ;  // there may be other test types than SwissCalc
  strcpy(header.description,description);
  fwrite( &header, sizeof(header), 1, out);
  return ftell(out);
  }  
  
static void doSwissCalcHeader( FILE* out, long int pos, testFixture* test, int numberOfRecords) {
  swissCalcHeader scHeader;
  for (int i=0;i<6;i++) scHeader.precisions[i] = test->precisions[i];  
  scHeader.numberOfRecords = numberOfRecords;
  fseek( out, pos, SEEK_SET);
  fwrite( &scHeader, sizeof(scHeader), 1, out);  // will be corrected later
  }  

static int doSwissCalcTestData( FILE* out, testFixture* test) {
  int numberOfRecords = 0;
  double dates[test->n_dates];
  test->method(dates,test->n_dates,test->jd_from,test->jd_to,test->incr);
  range_iterator
    planet_range = { &test->planets, 0, -1 },
    flags_range  = { &test->flags, 0, -1 };
  swissCalcData scData;
  while( (scData.planet = get_next(&planet_range)) != -1) {
    while( (scData.flags = get_next(&flags_range)) != -1) {
      for (int i=0;i<test->n_dates;i++) {
        scData.jd = dates[i];
        numberOfRecords++;
        scData.msg[0] = '\0';
        memset( scData.result, 0, 6*sizeof( double ) ) ;
        scData.flags_ret = swe_calc( dates[i], scData.planet, scData.flags, scData.result, scData.msg );        
        if (scData.flags_ret != scData.flags) {
          char bits1[33],bits2[33];
          int_to_binary(bits1,scData.flags_ret);
          int_to_binary(bits2,scData.flags);
          printf(
            "(Warning:) Flags changed! %15.10f %d %f %s : \n%32s (actual)\n%32s (requested)\n",
            scData.jd, scData.planet,
            scData.result[0], scData.msg,
            bits1,
            bits2
            );
          }
        fwrite( & scData, sizeof( scData) , 1, out );
        }
      }
    }
  return numberOfRecords;
  }
 
  
void get_random_dates(double dates[], int n, double from, double to, double dummy) {
  for (int i=0;i<n;i++) {
    dates[i] = from + rand()*(to-from)/RAND_MAX;  // Suffices for out purpose
    }
  }
void get_date_series(double dates[], int n, double from, double to, double incr) {
  for (int i=0;i<n;i++) {
    dates[i] = from + i*incr;
    }
  }

static int parseFixtureFile(char* fixtureFile, testFixture* tests, int *size) {
  FILE* fixture = fopen( fixtureFile, "r");
  if (fixture==NULL) {
      printf("Can't open %s\n",fixtureFile);
      return ERR;
      }
  char line[255], name[50], value[205];

  *size = 0;

  testFixture* data;

  for (;fgets(line,255,fixture);) {

    if (sscanf(line," %[T]ESTCASE",value)==1) {
      data = & tests[*size];
      (*size)++;
      clear( data );
      continue;
      }

// Skip everything before the first TESTCASE tag
    if (*size == 0) continue;

// Skip comments
    if (sscanf(line," %[#]",name)) continue;

// From here on: Only name:value pairs are scanned
    if (sscanf(line,"%[^: \t] : %[^\n#] \n",name,value)==2) {
      if (parseNameValuePair( name, value, data )==ERR) {
        printf("Error when parsing '%s:%s'\n",name,value);
        return ERR;
        }
      }

    }

  return OK;

  }

static int parseNameValuePair( char* name, char* value, testFixture* data) {
    if (equals(name,"description")) {
      strcpy(data->description,trim(value));
      }     
    else if (equals(name,"flags")) {
      return parseRange(value, &data->flags);
      }
    else if (equals(name,"planets")) {
      return parseRange(value, &data->planets);
      }
    else if (equals(name,"dates")) {
      return parseDates( value, data);
      }
    else if (equals(name,"jd")) {
      if (!sscanf(value,"%lf",&data->jd_from)) return ERR;
      data->n_dates = 1;
      data->jd_to = data->jd_from + 1;
      data->incr = 1;
      data->method = get_date_series;
      }
    else if (equals(name,"precisions")) {
      return parsePrecisions(value,data->precisions);
      }    
    else {
      printf("Couldn't understand '%s':'%s'\n",name,value);
      return ERR;
      }
  return OK;
  }

static int parseRange( char* value, rangetab *r) {
  char *p0 = value,
       *p1,
       *end = p0 + strlen(value);
  range line;
  do {
    int n = sscanf(p0,"%d - %d", &line.low, &line.high);
    line.between = (n==2);
    r->line[(r->size)++]=line;
    } while(
       (p1 = strchr(p0,',')) != NULL &&
       (p0 = p1+1)           <  end  &&
       r->size               < MAXRANGE );
  return OK;
  }

static int get_next( range_iterator *iter) {
  rangetab *r = iter->r;
  int value;
  if (iter->i >= r->size) {
    iter->i = 0;
    return -1;
    }
  if (r->line[iter->i].between) {
    if (iter->j == -1) {
      iter->j = r->line[iter->i].low;
      }
    else {
      (iter->j)++;
      if (iter->j > r->line[iter->i].high) {
        (iter->i)++;
        (iter->j) = -1;
        return get_next(iter);
        }
      }
    value = iter->j;
    }
  else {
    value = r->line[iter->i].low;
    (iter->i)++;
    (iter->j) = 0;
    }

  return value;

  }

static int parseDates( char* value, testFixture* data) {
  int d_from, d_to, m_from, m_to, y_from, y_to;
  char method[50];
    int args = sscanf(value,"%d . %d . %d - %d . %d . %d , %d , %s ",
       &d_from, &m_from, &y_from, &d_to, &m_to, &y_to, &(data->n_dates), method  );

    if (args < 7) {
    printf("Not enough parameters in '%s'\n",value);
    return ERR;
    }

  if (strcmp(method,"randomly")==0) {
    data->method = get_random_dates;
    } else {
    data->method = get_date_series;
    sscanf(method, "%lf",&data->incr);
    }

  data->jd_from = swe_julday( y_from, m_from, d_from, (y_from > 1582) ? 1 : 0, 0 );
  data->jd_to = swe_julday( y_to, m_to, d_to, (y_to > 1582) ? 1 : 0, 0 );

  return OK;
  }
  
static int parsePrecisions( char* value, int* precisions ) {

    char *vp = value,
         *vend = value + strlen(value);

    int i=0;
    while ( sscanf( vp, "%d", &precisions[i]) 
       && ++i < 6
       && (vp = strchr(vp,',')) != NULL
       && ++vp < vend );
  
  return OK;
  }  

static FILE* openFileBinary( const char* file ) {
  FILE* f = fopen( file,"wb");
  if (!f) {
    printf("Can't open %s\n",file);
    exit(EXIT_FAILURE);
    }
  return f;
}

static void clear( testFixture* data) {
  * data->description = '\0';
  data->planet = data->iflag = data->n_dates = 0;
  for (int i=0;i<6;i++) data->precisions[i] = i < 3 ? 0 : 1;
  data->jd_from = data->jd_to = data->incr = 0;
  data->planets.size = 0;
  data->flags.size = 0;
  }

static char *int_to_binary(char* b, int x) {
  char *p = b;
  for (unsigned int z = 1<<31; z ; z >>= 1) {
    *p++ = x & z ? '1' : '0';
    }
  *p = '\0';
  return b;
  }

  
static void doInitializations() {
// initialize pseudo random number generator
  const int RAND_SEED = 2621411;
  srand(RAND_SEED);
 }  