#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "swejpl.h"
#include "swephexp.h"
#include "sweph.h"
#include "swephlib.h"

// Sometimes, it gives more readable code when using named members instead of a double[6]
typedef struct {
  double x;
  double y;
  double z;
  } rect_coord; 
  
typedef struct {
  double l;
  double b;
  double r;
  } polar_coord;  
  
typedef struct {
  rect_coord pos;
  rect_coord speed;
  } rect_data;

typedef struct {
  polar_coord pos;
  polar_coord speed;
  } polar_data;

// New API definitions ...
int se_calc( double tjd, int ipl, int iflag, double *xx, char *serr, struct swe_data *swed);

// Internal functions
static void se_calc_prepare( double tjd, int *ipl, int *iflag, double* xx, char *serr, struct swe_data *swed);
static int se_calc_map_results( int ipl, int iflag, int ephemeris_requested, double *xint, double *xout);
static void se_close( struct swe_data *swed );
static void se_set_ephe_path( char *path, struct swe_data *swed);

static int swecalc(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed);


static const int IS_PLANET = 0;
static const int IS_MOON	 = 1;
static const int IS_ANY_BODY = 2;
static const int IS_MAIN_ASTEROID = 3;

#define DO_SAVE			TRUE
#define NO_SAVE			FALSE

#define SEFLG_EPHMASK	(SEFLG_JPLEPH|SEFLG_SWIEPH|SEFLG_MOSEPH)

struct meff_ele {double r,m;};

// Function references - the poor man's objects 
#define _calcfun(f) int f(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed);
typedef _calcfun((*calcfun))
  
static _calcfun(swecalc_eclnut)  
static _calcfun(swecalc_main_planet)
static _calcfun(swecalc_sun_barycentric)
static _calcfun(swecalc_moon_moshier)
static _calcfun(swecalc_moon_swisseph)
static _calcfun(swecalc_mean_lunar_node)
static _calcfun(swecalc_true_lunar_node)
static _calcfun(swecalc_osculating_apogee)
static _calcfun(swecalc_interpolated_perigee)
static _calcfun(swecalc_interpolated_apogee)
static _calcfun(swecalc_mean_apogee)
static _calcfun(swecalc_asteroid)
static _calcfun(swecalc_fictitious)

calcfun swecalc_get_calc_type(int ipl, int iflag, char* serr);  

/****************
 * global stuff *
 ****************/
struct swe_data FAR swed = {FALSE,	/* ephe_path_is_set = FALSE */
                            FALSE,	/* jpl_file_is_open = FALSE */
                            NULL,	/* fixed stars file pointer */
			    SE_EPHE_PATH,		/* ephe path */
			    SE_FNAME_DFT,	/* JPL file name, default */
			    FALSE,	/* geopos is set, for topocentric */
			    FALSE,	/* ayanamsa is set */
			    };  // ommitting rest of struct swe_data...


/*************
 * constants *
 *************/

static char *ayanamsa_name[] = {
   "Fagan/Bradley",
   "Lahiri",
   "De Luce",
   "Raman",
   "Ushashashi",
   "Krishnamurti",
   "Djwhal Khul",
   "Yukteshwar",
   "J.N. Bhasin",
   "Babylonian/Kugler 1",
   "Babylonian/Kugler 2",
   "Babylonian/Kugler 3",
   "Babylonian/Huber",
   "Babylonian/Eta Piscium",
   "Babylonian/Aldebaran = 15 Tau",
   "Hipparchos",
   "Sassanian",
   "Galact. Center = 0 Sag",
   "J2000",
   "J1900",
   "B1950",
};
static const int FAR pnoint2jpl[]   = PNOINT2JPL;

static const int pnoext2int[] = {SEI_SUN, SEI_MOON, SEI_MERCURY, SEI_VENUS, SEI_MARS, SEI_JUPITER, SEI_SATURN, SEI_URANUS, SEI_NEPTUNE, SEI_PLUTO, 0, 0, 0, 0, SEI_EARTH, SEI_CHIRON, SEI_PHOLUS, SEI_CERES, SEI_PALLAS, SEI_JUNO, SEI_VESTA, };

static const struct meff_ele FAR eff_arr[] = {
  /*
   * r , m_eff for photon passing the sun at min distance r (fraction of Rsun)
   * the values where computed with sun_model.c, which is a classic
   * treatment of a photon passing a gravity field, multiplied by 2.
   * The sun mass distribution m(r) is from Michael Stix, The Sun, p. 47.
   */
  {1.000, 1.000000},
  {0.990, 0.999979},
  {0.980, 0.999940},
  {0.970, 0.999881},
  {0.960, 0.999811},
  {0.950, 0.999724},
  {0.940, 0.999622},
  {0.930, 0.999497},
  {0.920, 0.999354},
  {0.910, 0.999192},
  {0.900, 0.999000},
  {0.890, 0.998786},
  {0.880, 0.998535},
  {0.870, 0.998242},
  {0.860, 0.997919},
  {0.850, 0.997571},
  {0.840, 0.997198},
  {0.830, 0.996792},
  {0.820, 0.996316},
  {0.810, 0.995791},
  {0.800, 0.995226},
  {0.790, 0.994625},
  {0.780, 0.993991},
  {0.770, 0.993326},
  {0.760, 0.992598},
  {0.750, 0.991770},
  {0.740, 0.990873},
  {0.730, 0.989919},
  {0.720, 0.988912},
  {0.710, 0.987856},
  {0.700, 0.986755},
  {0.690, 0.985610},
  {0.680, 0.984398},
  {0.670, 0.982986},
  {0.660, 0.981437},
  {0.650, 0.979779},
  {0.640, 0.978024},
  {0.630, 0.976182},
  {0.620, 0.974256},
  {0.610, 0.972253},
  {0.600, 0.970174},
  {0.590, 0.968024},
  {0.580, 0.965594},
  {0.570, 0.962797},
  {0.560, 0.959758},
  {0.550, 0.956515},
  {0.540, 0.953088},
  {0.530, 0.949495},
  {0.520, 0.945741},
  {0.510, 0.941838},
  {0.500, 0.937790},
  {0.490, 0.933563},
  {0.480, 0.928668},
  {0.470, 0.923288},
  {0.460, 0.917527},
  {0.450, 0.911432},
  {0.440, 0.905035},
  {0.430, 0.898353},
  {0.420, 0.891022},
  {0.410, 0.882940},
  {0.400, 0.874312},
  {0.390, 0.865206},
  {0.380, 0.855423},
  {0.370, 0.844619},
  {0.360, 0.833074},
  {0.350, 0.820876},
  {0.340, 0.808031},
  {0.330, 0.793962},
  {0.320, 0.778931},
  {0.310, 0.763021},
  {0.300, 0.745815},
  {0.290, 0.727557},
  {0.280, 0.708234},
  {0.270, 0.687583},
  {0.260, 0.665741},
  {0.250, 0.642597},
  {0.240, 0.618252},
  {0.230, 0.592586},
  {0.220, 0.565747},
  {0.210, 0.537697},
  {0.200, 0.508554},
  {0.190, 0.478420},
  {0.180, 0.447322},
  {0.170, 0.415454},
  {0.160, 0.382892},
  {0.150, 0.349955},
  {0.140, 0.316691},
  {0.130, 0.283565},
  {0.120, 0.250431},
  {0.110, 0.218327},
  {0.100, 0.186794},
  {0.090, 0.156287},
  {0.080, 0.128421},
  {0.070, 0.102237},
  {0.060, 0.077393},
  {0.050, 0.054833},
  {0.040, 0.036361},
  {0.030, 0.020953},
  {0.020, 0.009645},
  {0.010, 0.002767},
  {0.000, 0.000000}
};


static int do_fread(void *targ, int size, int count, int corrsize, 
		    FILE *fp, int fpos, int freord, int fendian, int ifno, 
		    char *serr, struct swe_data *swed);
static int get_new_segment(double tjd, int ipli, int ifno, char *serr, struct swe_data *swed);
static int main_planet(double tjd, int ipli, int epheflag, int iflag,
		       char *serr, struct swe_data* swed);
static int main_planet_bary(double tjd, int ipli, int epheflag, int iflag, 
		AS_BOOL do_save, 
		double *xp, double *xe, double *xs, double *xm, 
		char *serr, struct swe_data *swed);
static int sweplan(double tjd, int ipli, int ifno, int iflag, AS_BOOL do_save, 
		   double *xp, double *xpe, double *xps, double *xpm,
		   char *serr, struct swe_data *swed);
static int swemoon(double tjd, int iflag, AS_BOOL do_save, double *xp, char *serr, struct swe_data *swed);
static int sweph(double tjd, int ipli, int ifno, int iflag, double *xsunb, AS_BOOL do_save, 
		double *xp, char *serr, struct swe_data *swed);
static void rot_back(int ipl, struct swe_data *swed);
static int read_const(int ifno, char *serr, struct swe_data *swed);
static void embofs(double *xemb, double *xmoon);
static int app_pos_etc_plan(int ipli, int iflag, char *serr, struct swe_data *swed);
static int app_pos_etc_plan_osc(int ipl, int ipli, int iflag, char *serr, struct swe_data *swed);
static int app_pos_etc_sun(int iflag, char *serr, struct swe_data *swed);
static int app_pos_etc_moon(int iflag, char *serr, struct swe_data *swed);
static int app_pos_etc_sbar(int iflag, char *serr, struct swe_data *swed);
extern int swi_plan_for_osc_elem(int iflag, double tjd, double *xx);
static int swid_plan_for_osc_elem(int iflag, double tjd, double *xx, struct swe_data *swed);
static int app_pos_etc_mean(int ipl, int iflag, char *serr, struct swe_data *swed);
static void nut_matrix(struct nut *nu, struct epsilon *oec); 
static void calc_epsilon(double tjd, struct epsilon *e);
static int lunar_osc_elem(double tjd, int ipl, int iflag, char *serr, struct swe_data *swed);
static int intp_apsides(double tjd, int ipl, int iflag, char *serr, struct swe_data *swed); 
static double meff(double r);
static int plaus_iflag(int iflag_in, char* serr);
static int plaus_ipl(double tjd,int ipl_in, int iflag, char* serr);
void FAR PASCAL_CONV swe_set_sid_mode(int sid_mode, double t0, double ayan_t0);
void FAR PASCAL_CONV swed_set_sid_mode(int sid_mode, double t0, double ayan_t0, struct swe_data *swed);
static int app_pos_rest(struct plan_data *pdp, int iflag, 
    double *xx, double *x2000, struct epsilon *oe, char *serr, struct swe_data *swed);

// Functions that are declared elswhere and have to be decorated with a pointer to swed
FILE *swid_fopen(int ifno, char *fname, char *ephepath, char *serr, struct swe_data *swed);
int swed_fixstar(char *star, double tjd, int iflag, double *xx, char *serr, struct swe_data *swed);
int swed_fixstar_mag(char *star, double *mag, char *serr, struct swe_data *swed);
char *swed_get_planet_name(int ipl, char *s, struct swe_data *swed);
void swed_set_topo(double geolon, double geolat, double geoalt, struct swe_data *swed);
void swed_set_sid_mode(int sid_mode, double t0, double ayan_t0, struct swe_data *swed);
double swed_get_ayanamsa(double tjd_et, struct swe_data *swed);
void swid_precess_speed(double *xx, double t, int direction, struct swe_data *swed);
void swid_deflect_light(double *xx, double dt, int iflag, struct swe_data *swed);
int swid_trop_ra2sid_lon(double *xin, double *xout, double *xoutr, int iflag, char *serr, struct swe_data *swed);
int swid_trop_ra2sid_lon_sosy(double *xin, double *xout, double *xoutr, int iflag, char *serr, struct swe_data *swed);
int swid_get_observer(double tjd, int iflag, AS_BOOL do_save, double *xobs, char *serr, struct swe_data *swed);
void swid_check_ecliptic(double tjd, struct swe_data *swed);
void swid_check_nutation(double tjd, int iflag, struct swe_data *swed);
void swid_nutate(double *xx, int iflag, AS_BOOL backward, struct swe_data *swed);
    
// Utility functions (might move to swephlib)
static inline void* vec_clear( void* x, int size) {
  memset((void*)x,0,size*sizeof(double));
  return x;
  }

static inline void* vec_init( void* x, double s, int size) {
  double* fx = x;
  if (x) for (int i=size-1;i>=0;--i) fx[i] = s;
  return x;
  }

static inline void* vec_int_init( void* x, int i0, int size) {
  int* ix = x;
  if (x) for (int i=size-1;i>=0;--i) ix[i] = i0;
  return x;
  }

static inline void* vec_copy( void* to, void* from, int n) {
  if (to) memcpy(to,from,n*sizeof(double));
  return to;
  }

static inline void* vec_add( void* to, void* from, int n) {
  double *fto = (double *) to,
         *ffrom = (double *) from;
  if (to) for (int i=n-1;i>=0;--i) fto[i] += ffrom[i];
  return fto;
  }

static inline void* vec_sub( void* from, void* y, int n) {
  double *ffrom = (double *) from,
         *fy = (double *) y;
  if (from) for (int i=n-1;i>=0;--i) ffrom[i] -= fy[i];
  return from;
  }

static inline double vec_length( void* v, int n) {
  double *fv = (double *) v,
         s = 0;
  if (v) for (int i=n-1;i>=0;--i) s += fv[i]*fv[i];
  return v ? sqrt( s ) : 0;
  }

static inline void* smul( void* to, double s, int n) {
  double *fto = (double *) to;
  if (to) for (int i=n-1;i>=0;--i) fto[i] *= s;
  return fto;
  }

// Convert the arc components (longitude, latidude) from rad to deg or vice versa  
static inline void* convert_arcs( void* x, double factor) {
  double *fx = (double *) x;
  if (x) for (int i=4;i>=0;--i) if (i!=2) fx[i] *= factor;
  return x;
  }
  
  
static inline double vec_normalize(void *v, int n) {
  double d = vec_length(v,n);
  if (d > 0) smul( v, 1./d, n);
  return d;
  }  

static inline void* sdiv( void* to, double s, int n) {
  double *fto = (double *) to;
  if (to) for (int i=n-1;i>=0;--i) fto[i] /= s;
  return fto;
  }
static inline char* msg_append(char* serr, char* text) {
  int maxLengthAppendable = AS_MAXCH - strlen(serr);
  if (maxLengthAppendable > 0) {
    if (strlen(text) >= maxLengthAppendable)
      *(text+maxLengthAppendable-1) = '\0';
    strcat(serr, text);
    }
	return serr;  
  }
  
static inline bool is_set( int iflag, int flags ) {
  return iflag & flags;
  }  

static inline bool is_not_set( int iflag, int flags ) {
  return ! (iflag & flags);
  }  

static inline void set( int* iflag, int flags ) {
  (*iflag) |= flags;
  }  
  
static inline void set_ephemeris( int* iflag, int ephemeris) {
  *iflag = ( (*iflag) & ~SEFLG_EPHMASK ) | 
           ( ephemeris & SEFLG_EPHMASK ) 
         ;
  }  

static inline bool switch_to_moshier(double tjd, int ipli, int* iflag, char* serr) {
  if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
    set_ephemeris(iflag,SEFLG_MOSEPH);
    msg_append(serr," \nusing Moshier eph.; ");
    return true;
    }
  else return false;  
  }
  
static inline int get_ipli( int ipl ) {
  /* Normalized asteroid number (the bigger asteroids have two numbers each) */
  if (ipl < SE_NPLANETS)
    return pnoext2int[ipl];
  else if (ipl <= SE_AST_OFFSET + MPC_VESTA)
    return ipl - SE_AST_OFFSET + SEI_CERES - 1;
  else
    return ipl;
  }

static inline bool is_main_asteroid( int ipli ) {
  return ipli < SE_AST_OFFSET;
  }

static inline void throw(int retc, char* msg, char *serr, struct swe_data* swed) {
  msg_append( serr, msg );
  longjmp( swed->env, retc);  
  }  



  