// --- Data structures for testing the Swiss Ephemeris

typedef struct {
  int type; 
  char description[100];  
  } testHeader;
  
typedef struct {  
  int precisions[6];
  int numberOfRecords; 
  } swissCalcHeader;

typedef struct {
  double jd;
  int planet;
  int flags;
  int flags_ret;
  double result[6];
  char msg[255];
  } swissCalcData;  
  
  