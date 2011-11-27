/* test_test_calc_reduce.c */

typedef struct {
  int ok;
  int not_ok;
  } result;

bool test(char *filename);

static void fail(const char *msg);
static bool diff_check(const char *var, const double act, const double exp, const double prec);
static bool diff_check_int(const char *var, const int act, const int exp);
static bool diff_check_degree(const char *var, const double act, const double exp, const double prec);
static bool diff_check_char(const char *var, const char *act, const char *exp);

static bool calcsAsExpected(
  const double juldate, 
  const int planet, 
  const int flags, 
  const int flags_exp, 
  const double expectedResult[], 
  const double precision[], 
  const char *expectedMessage );
  
static void computePrecisions( double precisions[6], int precisionsExp[6] );                          
static result executeTestSet(FILE* testdata);
static int hasMoreTestSets(FILE* testdata, testHeader* header);
static double arcdiff(double x, double y);

 // Precision for checks on floating numbers
static double defaultPrecision = 1.0e-7;

const double DEG2ARC = 0.017453292519943295769236907684886;
const double deg360 = .36e3; 
const double deg180 = .18e3; 

static inline double arcdiff( double x, double y ) {
  double d = (x-y)/deg360+.5;  
  return (d-floor(d))*deg360-deg180;
  }

static inline void fail(const char* msg) {
  fprintf(stderr,msg);
  exit(EXIT_FAILURE);
  }

static inline bool diff_check( const char* var, const double act, const double exp, const double prec) {
  bool ok = fabs(act - exp ) < prec;  
  if (!ok) printf("\n%15s : %15.10f <> %15.10f",var,act,exp);
  return ok;
  }

static inline bool diff_check_int( const char* var, const int act, const int exp) {
  bool ok = (act == exp);  
  if (!ok) printf("\n%15s : %d <> %d",var,act,exp);
  return ok;
  }

static inline bool diff_check_degree( const char* var, const double act, const double exp, const double prec) {
  bool ok = fabs(arcdiff(act,exp) ) < prec;  
  if (!ok) printf("\n%15s : %15.10f <> %15.10f",var,act,exp);
  return ok;
  }
  
static inline bool diff_check_char( const char *var, const char* act, const char* exp) {
  bool ok = ! strcmp(act,exp);
  if (!ok) printf("\n%15s : '%s' <> '%s'",var,act,exp);
  return ok;
  }   




