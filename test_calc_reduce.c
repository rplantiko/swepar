// Experimental - trying to reduce the planet calculation

#include "test_calc_reduce.h"
#include "test_test_calc_reduce.c"

int main( int argc, char ** argv ) {

  double tjd, xx[6];
  char serr[255];
  int iflag = 0, iret, ipl;
  char testfile[100] = "tcr.bin";

  if (argc == 4) {
    sscanf( argv[1], "%lf", &tjd);
    sscanf( argv[2], "%d", &ipl);
    sscanf( argv[3], "%d", &iflag);

    iret = swe_calc( tjd, ipl, iflag, xx, serr );
    printf( "Return flags: %d\n", iret);
    for (int i=0;i<6;i++) printf( "xx[%d]=%15.10f\n", i,xx[i]);
    printf( "Message: %s\n", serr);
    }
  else if (argc == 2) {
    strcpy( testfile, argv[1] );
    }

  if (argc <= 2) {
// do some tests
    bool ok = test(testfile);
    exit ( ok ? EXIT_SUCCESS : EXIT_FAILURE );
    }

  }


// --- API functions, compatibility layer
// For usage see http://astro.com/swisseph/swephprg.htm
// Dead options in this version:
//   SEFLG_JPLEPH
//   SEFLG_SPEED3

int FAR PASCAL_CONV swe_calc(double tjd, int ipl, int iflag, double *xx, char *serr) {
  return se_calc(tjd,ipl,iflag,xx,serr,&swed);
  }

int FAR PASCAL_CONV swe_calc_ut(double tjd_ut, int ipl, int iflag, double *xx, char *serr) {
  return se_calc(tjd_ut + swe_deltat(tjd_ut), ipl, iflag, xx, serr, &swed);
  }

void FAR PASCAL_CONV swe_close(void) {
  se_close(&swed);
  }

void FAR PASCAL_CONV swe_set_ephe_path(char *path) {
  se_set_ephe_path(path,&swed);
  }

/* the ayanamsa (precession in longitude)
 * according to Newcomb's definition: 360 -
 * longitude of the vernal point of t referred to the
 * ecliptic of t0.
 */
double FAR PASCAL_CONV swe_get_ayanamsa(double tjd_et)
{
  return swed_get_ayanamsa(tjd_et,&swed);
}

double FAR PASCAL_CONV swe_get_ayanamsa_ut(double tjd_ut)
{
  return swe_get_ayanamsa(tjd_ut + swe_deltat(tjd_ut));
}

void FAR PASCAL_CONV swe_set_sid_mode(int sid_mode, double t0, double ayan_t0)
{
  swed_set_sid_mode(sid_mode,t0,ayan_t0,&swed);
}

int FAR PASCAL_CONV swe_fixstar(char *star, double tjd, int iflag,
  double *xx, char *serr)
{
  return swed_fixstar(star,tjd,iflag,xx,serr,&swed);
}

int FAR PASCAL_CONV swe_fixstar_ut(char *star, double tjd_ut, int iflag,
  double *xx, char *serr)
{
  return swe_fixstar(star, tjd_ut + swe_deltat(tjd_ut), iflag, xx, serr);
}

int FAR PASCAL_CONV swe_fixstar_mag(char *star, double *mag, char *serr)
{
  return swed_fixstar_mag(star,mag,serr,&swed);
}

char *FAR PASCAL_CONV swe_get_planet_name(int ipl, char *s)
{
  return (char*)swed_get_planet_name(ipl,s,&swed); // discard 'const' by cast (only for compatibility)
}

char *FAR PASCAL_CONV swe_get_ayanamsa_name(int isidmode)
{
  return (isidmode < SE_NSIDM_PREDEF) ?
      (char*)ayanamsa_name[isidmode] : NULL;  // discard 'const' by cast (only for compatibility)
}

void FAR PASCAL_CONV swe_set_topo(double geolon, double geolat, double geoalt)
{
  swed_set_topo(geolon,geolat,geoalt,&swed);
}



  
  
// --- New API functions
int se_calc( double tjd, int ipl, int iflag, double *xx, char *serr, struct swe_data *swed) {

  int iflag_act = iflag,
      ephemeris_requested = iflag & SEFLG_EPHMASK,
      ipl_act   = ipl,
      retc;
      
  if ((retc = setjmp(swed->env))==0) {    

    se_calc_prepare( tjd, &ipl_act, &iflag_act, xx, serr, swed);
    if (iflag_act==ERR) return ERR;
    
    double xint[24];  // Default result area for ephemeris computation
    double *xp = xint;
    
    int iflag_ret = swecalc(tjd, ipl_act, iflag_act, & xp, serr, swed);
    return iflag_ret >= 0 ?
             se_calc_map_results( ipl_act, iflag_ret, ephemeris_requested, xp, xx )
             : ERR;
             
    } 
  else { 
    return retc;
    }    
}

// Prepare flags
static void se_calc_prepare( double tjd, int *ipl, int *iflag, double *xx, char *serr, struct swe_data *swed) {

// Arrays need to be allocated!
  if ((serr == NULL)||(xx==NULL)) {
    *iflag = ERR;
    return;
    }

// Reset error message
  *serr = '\0';

// Reset result
  vec_clear(xx,6);

// Check and normalize flags
  *iflag = plaus_iflag(*iflag, serr);

// Check planet
  *ipl   = plaus_ipl( tjd, *ipl, *iflag, serr);

// To be moved "nearer to the file handling" :
  /* if ephemeris flag != ephemeris flag of last call,
   * we clear the save area, to prevent swecalc() using
   * previously computed data for current calculation.
   * except with *ipl = SE_ECL_NUT which is not dependent
   * on ephemeris, and except if change is from
   * ephemeris = 0 to ephemeris = SEFLG_DEFAULTEPH
   * or vice-versa.
   */
  static int epheflag_sv = 0;
  int epheflag = *iflag & SEFLG_EPHMASK;
  if ( (epheflag_sv != epheflag) &&
       (*ipl        != SE_ECL_NUT) ) {
    if (epheflag_sv) se_close(swed);
    epheflag_sv = epheflag;
    }

  }

static int se_calc_map_results( int ipl, int iflag, int ephemeris_requested, double *xint, double *xout) {

  if (iflag == ERR) return ERR;

  double *xs; // points to the requested double[6] set of coordinates

// Equatorial (offset 12) or ecliptical coordinates (offset 0)
  xs = is_set(iflag,SEFLG_EQUATORIAL) ? xint + 12 : xint;

// Cartesian coordinates  (offset 6 (ecl.) or 18 (equatorial) )
  if (is_set(iflag,SEFLG_XYZ)) xs += 6;

// How many result values do we have to copy?
  int imax = (ipl == SE_ECL_NUT) ? 4 :
               is_set(iflag,SEFLG_SPEED) ? 6 : 3;

// Copy result; from imax on, fill with zeroes
  vec_copy( xout, xs, imax);

// If requested, convert to radians
  if (is_set(iflag,SEFLG_RADIANS))
    convert_arcs( xout, DEGTORAD);

// If no ephemeris has been specified, do not return chosen ephemeris
// (as the caller may check iflag_new == iflag_old for regular termination)

  return ephemeris_requested ? iflag : iflag & ~SEFLG_DEFAULTEPH;

  }


static int swecalc(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {

  int epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag != SEFLG_MOSEPH && !swed->ephe_path_is_set)
    swe_set_ephe_path(NULL);

  if (is_set(iflag,SEFLG_SIDEREAL) && !swed->ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);

// --- set appropriate ecliptic
  swi_check_ecliptic(tjd);
  swi_check_nutation(tjd, iflag);

// --- Select computation method
  calcfun calc = swecalc_get_calc_type( ipl, iflag, serr, swed );

// Apply computation method
  int retc = calc(tjd,ipl,iflag,x,serr,swed);
  return retc < 0 ? retc : iflag;

  }

static int swecalc_eclnut(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
//  ecliptic and nutation                  *
  (*x)[0] = swed->oec.eps + swed->nut.nutlo[1];  /* true ecliptic */
  (*x)[1] = swed->oec.eps;      /* mean ecliptic */
  (*x)[2] = swed->nut.nutlo[0];    /* nutation in longitude */
  (*x)[3] = swed->nut.nutlo[1];    /* nutation in obliquity */
  smul(*x,RADTODEG,4);
  return iflag;
  }

static int swecalc_main_planet(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
  /* internal planet number */
  int ipli = pnoext2int[ipl];
  struct plan_data* pdp = &swed->pldat[ipli];
  *x = pdp->xreturn;
  int retc = main_planet(tjd, ipli,  iflag & SEFLG_EPHMASK, iflag, serr, swed);
  return retc;
  }

static int swecalc_sun_barycentric(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
  int retc;
    /* barycentric sun must be handled separately because of
     * the following reasons:
     * ordinary planetary computations use the function
     * main_planet() and its subfunction jplplan(),
     * see further below.
     * now, these functions need the swisseph internal
     * planetary indices, where SEI_EARTH = SEI_SUN = 0.
     * therefore they don't know the difference between
     * a barycentric sun and a barycentric earth and
     * always return barycentric earth.
     * to avoid this problem, many functions would have to
     * be changed. as an alternative, we choose a more
     * separate handling. */
    struct plan_data *pedp = &swed->pldat[SEI_EARTH];
    struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];

    *x = pedp->xreturn;

    if (is_set(iflag,SEFLG_SWIEPH)) {
  /* sweplan() provides barycentric sun as a by-product in save area;
   * it is saved in swed->pldat[SEI_SUNBARY].x */
        retc = sweplan(tjd, SEI_SUN, SEI_FILE_PLANET, iflag, DO_SAVE, NULL, NULL, NULL, NULL, serr, swed);
        if (retc < 0) return retc;
        psdp->teval = tjd;
      }
    else return ERR;

    /* flags */
    if ((retc = app_pos_etc_sbar(iflag, serr,swed)) < 0) return retc;

    /* iflag has possibly changed */
    iflag = pedp->xflgs;
    /* barycentric sun is now in save area of barycentric earth.
     * (pedp->xreturn = swed->pldat[SEI_EARTH].xreturn).
     * in case a barycentric earth computation follows for the same
     * date, the planetary functions will return the barycentric
     * SUN unless we force a new computation of pedp->xreturn.
     * this can be done by initializing the save of iflag.
     */
    pedp->xflgs = -1;
    return iflag;
  }

static int swecalc_moon_moshier(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {

  struct plan_data *pdp = &swed->pldat[SEI_MOON];
  *x = pdp->xreturn;
  if ( swi_moshmoon(tjd, DO_SAVE, NULL, serr) < 0) return ERR;
  /* for hel. position, we need earth as well */
  if (swi_moshplan(tjd, SEI_EARTH, DO_SAVE, NULL, NULL, serr) < 0) return ERR;
  return app_pos_etc_moon(iflag, serr, swed);
  }

static int swecalc_moon_swisseph(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
  int retc;
  struct plan_data* pdp = &swed->pldat[SEI_MOON];
  *x = pdp->xreturn;
  if ((retc = sweplan(tjd, SEI_MOON, SEI_FILE_MOON, iflag, DO_SAVE,NULL, NULL, NULL, NULL, serr,swed)) == ERR) return ERR;

  /* if sweph file not found, switch to moshier */
  if (retc == NOT_AVAILABLE) {
    if (tjd > MOSHLUEPH_START && tjd < MOSHLUEPH_END) {
       set_ephemeris(&iflag, SEFLG_MOSEPH);
       msg_append(serr,"\nusing Moshier eph.;");
       return swecalc_moon_moshier(tjd,ipl,iflag,x,serr,swed);
       }
    else return ERR;
    }
  retc = app_pos_etc_moon(iflag, serr, swed);
  return retc < 0 ? retc : iflag;
  }

static inline void lunar_node_clear_irrelevant_coordinates(int iflag, double* x) {
/* to avoid infinitesimal deviations from latitude = 0
   that result from conversions */
  if (is_not_set(iflag, SEFLG_SIDEREAL | SEFLG_J2000 )) {
    x[1] = 0.0;  /* ecl. latitude       */
    x[4] = 0.0;  /*               speed */
    x[5] = 0.0;  /*      radial   speed */
    x[8] = 0.0;  /* z coordinate        */
    x[11] = 0.0;  /*               speed */
    }  
  }
  
static int swecalc_mean_lunar_node(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
  int retc;
  struct plan_data* ndp = &swed->nddat[SEI_MEAN_NODE];
  double *ndp_x = ndp->x;

  *x = ndp->xreturn;

  if ((retc = swi_mean_node(tjd, ndp_x, serr)) == ERR) return ERR;

  /* speed (is almost constant; variation < 0.001 arcsec) */
  if ((retc = swi_mean_node(tjd - MEAN_NODE_SPEED_INTV, ndp_x+3, serr)) == ERR) return ERR;

  ndp_x[3] = swe_difrad2n(ndp_x[0], ndp_x[3]) / MEAN_NODE_SPEED_INTV;
  ndp_x[4] = ndp_x[5] = 0;
  ndp->teval = tjd;
  ndp->xflgs = -1;

  /* lighttime etc. */
  if ((retc = app_pos_etc_mean(SEI_MEAN_NODE, iflag, serr,swed)) < 0) return retc;

  /* to avoid infinitesimal deviations from latitude = 0
     that result from conversions */
  lunar_node_clear_irrelevant_coordinates(iflag,ndp->xreturn); 

  return retc;
  }
  
static int swecalc_true_lunar_node(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
  
  int retc;
  struct plan_data* ndp = &swed->nddat[SEI_TRUE_NODE];
  *x = ndp->xreturn;
  retc = lunar_osc_elem(tjd, SEI_TRUE_NODE, iflag, serr,swed);
  iflag = ndp->xflgs;

  /* to avoid infinitesimal deviations from latitude = 0
     that result from conversions */
  lunar_node_clear_irrelevant_coordinates(iflag,ndp->xreturn);

  return retc;
  }

static int swecalc_osculating_apogee(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
    struct plan_data* ndp = &swed->nddat[SEI_OSCU_APOG];
    *x = ndp->xreturn;
    return lunar_osc_elem(tjd, SEI_OSCU_APOG, iflag, serr, swed);
  }
static int swecalc_interpolated_apogee(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
    struct plan_data* ndp = &swed->nddat[SEI_INTP_APOG];
    *x = ndp->xreturn;
    return intp_apsides(tjd, SEI_INTP_APOG, iflag, serr,swed);
  }
static int swecalc_interpolated_perigee(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
    struct plan_data* ndp = &swed->nddat[SEI_INTP_PERG];
    *x = ndp->xreturn;
    return intp_apsides(tjd, SEI_INTP_PERG, iflag, serr,swed);
  }
static int swecalc_mean_apogee(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {

    struct plan_data* ndp = &swed->nddat[SEI_MEAN_APOG];
    double *xp2 = ndp->x;

    *x = ndp->xreturn;

    if (swi_mean_apog(tjd, xp2, serr) < 0) return ERR;

    /* speed (is not constant! variation ~= several arcsec) */
    if (swi_mean_apog(tjd - MEAN_NODE_SPEED_INTV, xp2+3, serr) < 0) return ERR;

    for(int i = 0; i <= 1; i++)
      xp2[3+i] = swe_difrad2n(xp2[i], xp2[3+i]) / MEAN_NODE_SPEED_INTV;
    xp2[5] = 0;
    ndp->teval = tjd;
    ndp->xflgs = -1;
    /* lighttime etc. */
    if (app_pos_etc_mean(SEI_MEAN_APOG, iflag, serr,swed) <0) return ERR;

    /* to avoid infinitesimal deviations from r-speed = 0
     * that result from conversions */
    ndp->xreturn[5] = 0.0;  /*               speed */
    return iflag;

  }


int swecalc_asteroid(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {

    int retc;

    int ipli = get_ipli( ipl );
    int ifno = is_main_asteroid(ipli) ? SEI_FILE_MAIN_AST : SEI_FILE_ANY_AST;

    char serr_sun[AS_MAXCH];
    *serr_sun = '\0';

    /* earth and sun are also needed */
    if ((retc = main_planet(tjd, SEI_EARTH, iflag & SEFLG_EPHMASK, iflag, serr_sun, swed)) < 0) return retc;

    /* iflag (ephemeris bit) has possibly changed in main_planet() */
    iflag = swed->pldat[SEI_EARTH].xflgs;


    struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
    struct plan_data* pdp  = &swed->pldat[ is_main_asteroid(ipli) ? ipli: SEI_ANYBODY ];
    *x = pdp->xreturn;

    retc = sweph(tjd, ipli, ifno, iflag, psdp->x, DO_SAVE, NULL, serr, swed);
    if (retc == ERR || retc == NOT_AVAILABLE) return ERR;

    if ((retc = app_pos_etc_plan(ipli, iflag, serr, swed)) == ERR) return ERR;

    /* app_pos_etc_plan() might have failed, if t(light-time)
     * is beyond ephemeris range. in this case redo with Moshier
     */
    if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS) {
      if (is_not_set(iflag, SEFLG_MOSEPH)) {
        set_ephemeris(&iflag, SEFLG_MOSEPH);
        msg_append(serr, "\nusing Moshier eph.; ");
        return swecalc_asteroid(tjd,ipl,iflag,x,serr,swed);
      } else return ERR;
    }

    /* add warnings from earth/sun computation */
    if (*serr_sun != '\0') msg_append( msg_append( serr, "sun:" ), serr_sun );

    return iflag;

  }

int swecalc_fictitious(double tjd, int ipl, int iflag, double **x, char *serr, struct swe_data *swed) {
   int retc;
    /* internal planet number */
    int ipli = SEI_ANYBODY;
    struct plan_data *pdp  = &swed->pldat[ipli];
    struct plan_data *pedp = &swed->pldat[SEI_EARTH];
    struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
    *x = pdp->xreturn;

    /* the earth for geocentric position */
    retc = main_planet(tjd, SEI_EARTH, iflag & SEFLG_EPHMASK, iflag, serr, swed);
    /* iflag (ephemeris bit) has possibly changed in main_planet() */
    iflag = swed->pldat[SEI_EARTH].xflgs;
    /* planet from osculating elements */
    if (swi_osc_el_plan(tjd, pdp->x, ipl-SE_FICT_OFFSET, ipli, pedp->x, psdp->x, serr) < 0) return ERR;
    if (retc == ERR) return ERR;

    if (app_pos_etc_plan_osc(ipl, ipli, iflag, serr, swed) < 0) return ERR;

    /* app_pos_etc_plan_osc() might have failed, if t(light-time)
     * is beyond ephemeris range. in this case redo with Moshier
     */
    if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS) {
      if (is_not_set(iflag, SEFLG_MOSEPH)) {
        set_ephemeris( &iflag, SEFLG_MOSEPH );
        msg_append(serr, "\nusing Moshier eph.;" );
        return swecalc_fictitious(tjd,ipl,iflag,x,serr,swed);
        }
      else return ERR;
    }
  return iflag;
  }

// Which calculation method has to be chosen
calcfun swecalc_get_calc_type(int ipl, int iflag, char* serr, struct swe_data *swed) {

// Main planet
  if ( ( ipl == SE_SUN   && is_not_set(iflag,SEFLG_HELCTR) )
    || ( ipl == SE_EARTH && is_set(iflag, SEFLG_HELCTR | SEFLG_BARYCTR) ) 
    || ipl == SE_MERCURY
    || ipl == SE_VENUS
    || ipl == SE_MARS
    || ipl == SE_JUPITER
    || ipl == SE_SATURN
    || ipl == SE_URANUS
    || ipl == SE_NEPTUNE
    || ipl == SE_PLUTO ) return swecalc_main_planet;

// Barycentric Sun is special
  if (ipl == SE_SUN && is_set(iflag, SEFLG_BARYCTR)) return swecalc_sun_barycentric;

// Moon - distinguish requested ephemeris
  if (ipl == SE_MOON)
    return (is_set(iflag,SEFLG_MOSEPH)) ? swecalc_moon_moshier : swecalc_moon_swisseph;

// Minor planet
  if ( ipl == SE_CHIRON
    || ipl == SE_PHOLUS
    || ipl == SE_CERES
    || ipl == SE_PALLAS
    || ipl == SE_JUNO
    || ipl == SE_VESTA
    || ipl > SE_AST_OFFSET) return swecalc_asteroid;

// Fictitious planet
  if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX) return swecalc_fictitious;

// Obliquity
  if (ipl == SE_ECL_NUT) return swecalc_eclnut;

// Objects related to lunar orbit
  calcfun lunar_orbit_object = 0;
  switch (ipl) {
    case SE_MEAN_NODE:
      lunar_orbit_object = swecalc_mean_lunar_node;
      break;
    case SE_TRUE_NODE:
      lunar_orbit_object = swecalc_true_lunar_node;
      break;
    case SE_OSCU_APOG:
      lunar_orbit_object = swecalc_osculating_apogee;
      break;
    case SE_MEAN_APOG:
      lunar_orbit_object = swecalc_mean_apogee;
      break;
    case SE_INTP_PERG:
      lunar_orbit_object = swecalc_interpolated_perigee;
      break;
    case SE_INTP_APOG:
      lunar_orbit_object = swecalc_interpolated_apogee;
      break;
    }
  if (lunar_orbit_object) {
    if (!(iflag & (SEFLG_HELCTR | SEFLG_BARYCTR))) {
      return lunar_orbit_object;
      }
    }

  sprintf(serr, "illegal planet number %d.",ipl);
  throw(ERR,serr,serr,swed);
  return NULL; // will not be reached
  }



  /* closes all open files, frees space of planetary data,
 * deletes memory of all computed positions
 */

static void se_close( struct swe_data *swed ) {
  int i;
  /* close SWISSEPH files */
  for (i = 0; i < SEI_NEPHFILES; i ++) {
    if (swed->fidat[i].fptr != NULL)
      fclose(swed->fidat[i].fptr);
    memset((void *) &swed->fidat[i], 0, sizeof(struct file_data));
  }
  /* free planets data space */
  for (i = 0; i < SEI_NPLANETS; i++) {
    if (swed->pldat[i].segp != NULL) {
      free((void *) swed->pldat[i].segp);
    }
    if (swed->pldat[i].refep != NULL) {
      free((void *) swed->pldat[i].refep);
    }
    memset((void *) &swed->pldat[i], 0, sizeof(struct plan_data));
  }
  for (i = 0; i <= SE_NPLANETS; i++) /* "<=" is correct! see decl. */
    memset((void *) &swed->savedat[i], 0, sizeof(struct save_positions));
  /* clear node data space */
  for (i = 0; i < SEI_NNODE_ETC; i++) {
    memset((void *) &swed->nddat[i], 0, sizeof(struct plan_data));
  }
  memset((void *) &swed->oec, 0, sizeof(struct epsilon));
  memset((void *) &swed->oec2000, 0, sizeof(struct epsilon));
  memset((void *) &swed->nut, 0, sizeof(struct nut));
  memset((void *) &swed->nut2000, 0, sizeof(struct nut));
  memset((void *) &swed->nutv, 0, sizeof(struct nut));
  /* close JPL file */
  swed->jpl_file_is_open = FALSE;
  /* close fixed stars */
  if (swed->fixfp != NULL) {
    fclose(swed->fixfp);
    swed->fixfp = NULL;
  }
}

/* sets ephemeris file path.
 * also calls swe_close(). this makes sure that swe_calc()
 * won't return planet positions previously computed from other
 * ephemerides
 */
static void se_set_ephe_path( char *path, struct swe_data *swed) {
  int i;
  char s[AS_MAXCH];
  char *sp;
  swed->ephe_path_is_set = TRUE;
  /* close all open files and delete all planetary data */
  se_close(swed);
  /* environment variable SE_EPHE_PATH has priority */
  if ((sp = getenv("SE_EPHE_PATH")) != NULL
    && strlen(sp) != 0
    && strlen(sp) <= AS_MAXCH-1-13) {
    strcpy(s, sp);
  } else if (path == NULL) {
    strcpy(s, SE_EPHE_PATH);
  } else if (strlen(path) <= AS_MAXCH-1-13)
    strcpy(s, path);
  else
    strcpy(s, SE_EPHE_PATH);
  i = strlen(s);
  if (*(s + i - 1) != *DIR_GLUE && *s != '\0')
    strcat(s, DIR_GLUE);
  strcpy(swed->ephepath, s);
}

/* calculates obliquity of ecliptic and stores it together
 * with its date, sine, and cosine
 */
static void calc_epsilon(double tjd, struct epsilon *e)
{
    e->teps = tjd;
    e->eps = swi_epsiln(tjd);
    e->seps = sin(e->eps);
    e->ceps = cos(e->eps);
}

/* computes a main planet from any ephemeris, if it
 * has not yet been computed for this date.
 * since a geocentric position requires the earth, the
 * earth's position will be computed as well. With SWISSEPH
 * files the barycentric sun will be done as well.
 * With Moshier, the moon will be done as well.
 *
 * tjd     = julian day
 * ipli    = body number
 * epheflag  = which ephemeris? JPL, SWISSEPH, Moshier?
 * iflag  = other flags
 *
 * the geocentric apparent position of ipli (or whatever has
 * been specified in iflag) will be saved in
 * &swed->pldat[ipli].xreturn[];
 *
 * the barycentric (heliocentric with Moshier) position J2000
 * will be kept in
 * &swed->pldat[ipli].x[];
 */
static int main_planet(double tjd, int ipli, int epheflag, int iflag, char *serr, struct swe_data *swed)
{
  int retc;
  switch(epheflag) {
    case SEFLG_SWIEPH:
      /* compute barycentric planet (+ earth, sun, moon) */
      retc = sweplan(tjd, ipli, SEI_FILE_PLANET, iflag, DO_SAVE, NULL, NULL, NULL, NULL, serr, swed);
      if (retc == ERR) return ERR;
      /* if sweph file not found, switch to moshier */
      if (retc == NOT_AVAILABLE) {
        if (switch_to_moshier(tjd,ipli,&iflag,serr)) goto moshier_planet;
        else return ERR; 
        }
      
      /* geocentric, lighttime etc. */
      retc = (ipli == SEI_SUN) 
               ? app_pos_etc_sun(iflag, serr, swed) 
               : app_pos_etc_plan(ipli, iflag, serr, swed); 
      if (retc == ERR) return ERR;
      
      /* if sweph file for t(lighttime) not found, switch to moshier */
      if (retc == NOT_AVAILABLE) {
        if (switch_to_moshier(tjd,ipli,&iflag,serr)) goto moshier_planet;
        else return ERR; 
        }
      break;
    case SEFLG_MOSEPH:
      moshier_planet:
      retc = swi_moshplan(tjd, ipli, DO_SAVE, NULL, NULL, serr);/**/
      if (retc == ERR) return ERR;
      /* geocentric, lighttime etc. */
      retc = (ipli == SEI_SUN) ? 
        app_pos_etc_sun(iflag, serr, swed) : app_pos_etc_plan(ipli, iflag, serr, swed);
      if (retc == ERR) return ERR;
      break;
    default:
      break;
  }
  return iflag;
}

/* Computes a main planet from any ephemeris or returns
 * it again, if it has been computed before.
 * In barycentric equatorial position of the J2000 equinox.
 * The earth's position is computed as well. With SWISSEPH
 * and JPL ephemeris the barycentric sun is computed, too.
 * With Moshier, the moon is returned, as well.
 *
 * tjd     = julian day
 * ipli    = body number
 * epheflag  = which ephemeris? JPL, SWISSEPH, Moshier?
 * iflag  = other flags
 * xp, xe, xs, and xm are the pointers, where the program
 * either finds or stores (if not found) the barycentric
 * (heliocentric with Moshier) positions of the following
 * bodies:
 * xp    planet
 * xe    earth
 * xs    sun
 * xm    moon
 *
 * xm is used with Moshier only
 */
static int main_planet_bary(double tjd, int ipli, int epheflag, int iflag, AS_BOOL do_save,
           double *xp, double *xe, double *xs, double *xm,
           char *serr, struct swe_data *swed)
{
  int retc;
  switch(epheflag) {
    case SEFLG_SWIEPH:
      /* compute barycentric planet (+ earth, sun, moon) */
      retc = sweplan(tjd, ipli, SEI_FILE_PLANET, iflag, do_save, xp, xe, xs, xm, serr,swed);
      if (retc == ERR || retc == NOT_AVAILABLE)
  return retc;
      break;
    case SEFLG_MOSEPH:
      retc = swi_moshplan(tjd, ipli, do_save, xp, xe, serr);/**/
      if (retc == ERR) return ERR;
      vec_clear(xs,6);
      break;
    default:
      break;
  }
  return OK;
}

/* SWISSEPH
 * this routine computes heliocentric cartesian equatorial coordinates
 * of equinox 2000 of
 * geocentric moon
 *
 * tjd     julian date
 * iflag  flag
 * do_save  save J2000 position in save area pdp->x ?
 * xp    array of 6 doubles for lunar position and speed
 * serr    error string
 */
static int swemoon(double tjd, int iflag, AS_BOOL do_save, double *xpret, char *serr, struct swe_data *swed)
{
  struct plan_data *pdp = &swed->pldat[SEI_MOON];
  int speedf1, speedf2;
  double xx[6], *xp;
  if (do_save)
    xp = pdp->x;
  else
    xp = xx;
  /* if planet has already been computed for this date, return
   * if speed flag has been turned on, recompute planet */
  speedf1 = pdp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == pdp->teval
  && pdp->iephe == SEFLG_SWIEPH
  && (!speedf2 || speedf1)) {
    xp = pdp->x;
  } else {
    /* call sweph for moon */
    int retc = sweph(tjd, SEI_MOON, SEI_FILE_MOON, iflag, NULL, do_save, xp, serr,swed);
    if (retc != OK) return(retc);
    if (do_save) {
      pdp->teval = tjd;
      pdp->xflgs = -1;
      pdp->iephe = SEFLG_SWIEPH;
    }
  }
  vec_copy(xpret,xp,6);
  return(OK);
}

/* SWISSEPH
 * this function computes
 * 1. a barycentric planet
 * plus, under certain conditions,
 * 2. the barycentric sun,
 * 3. the barycentric earth, and
 * 4. the geocentric moon,
 * in barycentric cartesian equatorial coordinates J2000.
 *
 * these are the data needed for calculation of light-time etc.
 *
 * tjd     julian date
 * ipli    SEI_ planet number
 * ifno    ephemeris file number
 * do_save  write new positions in save area
 * xp    array of 6 doubles for planet's position and velocity
 * xpe                                 earth's
 * xps                                 sun's
 * xpm                                 moon's
 * serr    error string
 *
 * xp - xpm can be NULL. if do_save is TRUE, all of them can be NULL.
 * the positions will be written into the save area (swed->pldat[ipli].x)
 */
static int sweplan(double tjd, int ipli, int ifno, int iflag, AS_BOOL do_save,
       double *xpret, double *xperet, double *xpsret, double *xpmret,
       char *serr, struct swe_data *swed)
{
  int retc;
  int do_earth = FALSE, do_moon = FALSE, do_sunbary = FALSE;
  struct plan_data *pdp = &swed->pldat[ipli];
  struct plan_data *pebdp = &swed->pldat[SEI_EMB];
  struct plan_data *psbdp = &swed->pldat[SEI_SUNBARY];
  struct plan_data *pmdp = &swed->pldat[SEI_MOON];
  double xxp[6], xxm[6], xxs[6], xxe[6];
  double *xp, *xpe, *xpm, *xps;
  int speedf1, speedf2;
  /* xps (barycentric sun) may be necessary because some planets on sweph
   * file are heliocentric, other ones are barycentric. without xps,
   * the heliocentric ones cannot be returned barycentrically.
   */
  if (do_save || ipli == SEI_SUNBARY || (pdp->iflg & SEI_FLG_HELIO)
    || xpsret != NULL || is_set(iflag,SEFLG_HELCTR))
    do_sunbary = TRUE;
  if (do_save || ipli == SEI_EARTH || xperet != NULL)
    do_earth = TRUE;
  if (ipli == SEI_MOON) {
    do_earth = TRUE;
    do_sunbary = TRUE;
  }
  if (do_save || ipli == SEI_MOON || ipli == SEI_EARTH || xperet != NULL || xpmret != NULL)
    do_moon = TRUE;
  if (do_save)  {
    xp = pdp->x;
    xpe = pebdp->x;
    xps = psbdp->x;
    xpm = pmdp->x;
  } else {
    xp = xxp;
    xpe = xxe;
    xps = xxs;
    xpm = xxm;
  }
  speedf2 = iflag & SEFLG_SPEED;
  /* barycentric sun */
  if (do_sunbary) {
    speedf1 = psbdp->xflgs & SEFLG_SPEED;
    /* if planet has already been computed for this date, return
     * if speed flag has been turned on, recompute planet */
    if (tjd == psbdp->teval
      && psbdp->iephe == SEFLG_SWIEPH
      && (!speedf2 || speedf1)) {
      vec_copy(xps, psbdp->x, 6);
      } 
    else {
      retc = sweph(tjd, SEI_SUNBARY, SEI_FILE_PLANET, iflag, NULL, do_save, xps, serr, swed);/**/
      if (retc != OK) return(retc);
      }
    vec_copy(xpsret,xps,6);
    }
  /* moon */
  if (do_moon) {
    speedf1 = pmdp->xflgs & SEFLG_SPEED;
    if (tjd == pmdp->teval
    && pmdp->iephe == SEFLG_SWIEPH
    && (!speedf2 || speedf1)) {
      vec_copy(xpm,pmdp->x,6);
      } 
    else {
      retc = sweph(tjd, SEI_MOON, SEI_FILE_MOON, iflag, NULL, do_save, xpm, serr, swed);
      if (retc == ERR)  return(retc);
      /* if moon file doesn't exist, take moshier moon */
      if (swed->fidat[SEI_FILE_MOON].fptr == NULL) {
        msg_append(serr," \nusing Moshier eph. for moon; ");
        retc = swi_moshmoon(tjd, do_save, xpm, serr);
        if (retc != OK) return(retc);
      }
    }
    vec_copy(xpmret,xpm,6);
    }
    
  /* barycentric earth */
  if (do_earth) {
    speedf1 = pebdp->xflgs & SEFLG_SPEED;
    if (tjd == pebdp->teval
    && pebdp->iephe == SEFLG_SWIEPH
    && (!speedf2 || speedf1)) {
      vec_copy(xpe,pebdp->x,6);
      } 
    else {
      retc = sweph(tjd, SEI_EMB, SEI_FILE_PLANET, iflag, NULL, do_save, xpe, serr, swed);
      if (retc != OK)  return(retc);
      /* earth from emb and moon */
      embofs(xpe, xpm);
      /* speed is needed, if
       * 1. true position is being computed before applying light-time etc.
       *    this is the position saved in pdp->x.
       *    in this case, speed is needed for light-time correction.
       * 2. the speed flag has been specified.
       */
      if (xpe == pebdp->x || is_set(iflag,SEFLG_SPEED))
        embofs(xpe+3, xpm+3);
      }
    vec_copy(xperet,xpe,6);
    }
    
  if (ipli == SEI_MOON) {
    vec_copy(xp,xpm,6);
    } 
  else if (ipli == SEI_EARTH) {
    vec_copy(xp,xpe,6);
    } 
  else if (ipli == SEI_SUN) {
    vec_copy(xp,xps,6);
    } 
  else {
    /* planet */
    speedf1 = pdp->xflgs & SEFLG_SPEED;
    if (tjd == pdp->teval
    && pdp->iephe == SEFLG_SWIEPH
    && (!speedf2 || speedf1)) {
      vec_copy(xp,pdp->x,6);
      return(OK);
      } 
    else {
      retc = sweph(tjd, ipli, ifno, iflag, NULL, do_save, xp, serr, swed);
      if (retc != OK)  return(retc);
      /* if planet is heliocentric, it must be transformed to barycentric */
      if (pdp->iflg & SEI_FLG_HELIO) {
      /* now barycentric planet */
        vec_add(xp,xps, do_save || is_set(iflag,SEFLG_SPEED) ? 6 : 3);
      }
    }
  }
  vec_copy(xpret,xp,6);
  return(OK);
}

/*
 * this function looks for an ephemeris file,
 * opens it, if not yet open,
 * reads constants, if not yet read,
 * computes a planet, if not yet computed
 * attention: asteroids are heliocentric
 *            other planets barycentric
 *
 * tjd     julian date
 * ipli    SEI_ planet number
 * ifno    ephemeris file number
 * xsunb  INPUT (!) array of 6 doubles containing barycentric sun
 *              (must be given with asteroids)
 * do_save  boolean: save result in save area
 * xp    return array of 6 doubles for planet's position
 * serr    error string
 */
static int sweph(double tjd, int ipli, int ifno, int iflag, double *xsunb, AS_BOOL do_save, double *xpret, char *serr, struct swe_data *swed)
{
  int i, ipl, retc, subdirlen;
  char s[AS_MAXCH], subdirnam[AS_MAXCH], fname[AS_MAXCH], *sp;
  double t, tsv;
  double xemb[6], xx[6], *xp;
  struct plan_data *pdp;
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
  struct file_data *fdp = &swed->fidat[ifno];
  int speedf1, speedf2;
  AS_BOOL need_speed;
  ipl = ipli;
  if (ipli > SE_AST_OFFSET)
    ipl = SEI_ANYBODY;
  pdp = &swed->pldat[ipl];
  if (do_save)
    xp = pdp->x;
  else
    xp = xx;
  /* if planet has already been computed for this date, return.
   * if speed flag has been turned on, recompute planet */
  speedf1 = pdp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == pdp->teval
  && pdp->iephe == SEFLG_SWIEPH
  && (!speedf2 || speedf1)
        && ipl < SEI_ANYBODY) {
    vec_copy(xpret,pdp->x,6);
    return(OK);
  }
  /******************************
   * get correct ephemeris file *
   ******************************/
  if (fdp->fptr != NULL) {
    /* if tjd is beyond file range, close old file.
     * if new asteroid, close old file. */
    if (tjd < fdp->tfstart || tjd > fdp->tfend
      || (ipl == SEI_ANYBODY && ipli != pdp->ibdy)) {
      fclose(fdp->fptr);
      fdp->fptr = NULL;
      if (pdp->refep != NULL)
  free((void *) pdp->refep);
      pdp->refep = NULL;
      if (pdp->segp != NULL)
  free((void *) pdp->segp);
      pdp->segp = NULL;
    }
  }
  /* if sweph file not open, find and open it */
  if (fdp->fptr == NULL) {
    swi_gen_filename(tjd, ipli, fname);
    strcpy(subdirnam, fname);
    sp = strrchr(subdirnam, (int) *DIR_GLUE);
    if (sp != NULL) {
      *sp = '\0';
      subdirlen = strlen(subdirnam);
    } else {
      subdirlen = 0;
    }
    strcpy(s, fname);
again:
    fdp->fptr = swid_fopen(ifno, s, swed->ephepath, serr, swed);
    if (fdp->fptr == NULL) {
      /*
       * if it is a numbered asteroid file, try also for short files (..s.se1)
       * On the second try, the inserted 's' will be seen and not tried again.
       */
      if (ipli > SE_AST_OFFSET) {
  char *spp;
  spp = strchr(s, '.');
  if (spp > s && *(spp-1) != 's') {  /* no 's' before '.' ? */
    sprintf(spp, "s.%s", SE_FILE_SUFFIX);  /* insert an 's' */
    goto again;
  }
  /*
   * if we still have 'ast0' etc. in front of the filename,
   * we remove it now, remove the 's' also,
   * and try in the main ephemeris directory instead of the
   * asteroid subdirectory.
   */
        spp--;  /* point to the character before '.' which must be a 's' */
  swi_strcpy(spp, spp + 1);  /* remove the s */
  if (subdirlen > 0 && strncmp(s, subdirnam, (size_t) subdirlen) == 0) {
    swi_strcpy(s, s + subdirlen + 1);  /* remove "ast0/" etc. */
    goto again;
  }
      }
      return(NOT_AVAILABLE);
    }
    /* during the search error messages may have been built, delete them */
    if (serr != NULL) *serr = '\0';
    read_const(ifno, serr, swed);

  }

  /* if first ephemeris file (J-3000), it might start a mars period
   * after -3000. if last ephemeris file (J3000), it might end a
   * 4000-day-period before 3000. */
  if (tjd < fdp->tfstart || tjd > fdp->tfend) {
    if (serr != NULL) {
      if (tjd < fdp->tfstart)
  sprintf(s, "jd %f < Swiss Eph. lower limit %f;",
      tjd, fdp->tfstart);
      else
  sprintf(s, "jd %f > Swiss Eph. upper limit %f;",
      tjd, fdp->tfend);
      if (strlen(serr) + strlen(s) < AS_MAXCH)
  strcat(serr, s);
    }
    return(NOT_AVAILABLE);
  }
  /******************************
   * get planet's position
   ******************************/
  /* get new segment, if necessary */
  if (pdp->segp == NULL || tjd < pdp->tseg0 || tjd > pdp->tseg1) {
    get_new_segment(tjd, ipl, ifno, serr, swed);
    /* rotate cheby coeffs back to equatorial system.
     * if necessary, add reference orbit. */
    if (pdp->iflg & SEI_FLG_ROTATE)
      rot_back(ipl,swed); /**/
    else
      pdp->neval = pdp->ncoe;
  }
  /* evaluate chebyshew polynomial for tjd */
  t = (tjd - pdp->tseg0) / pdp->dseg;
  t = t * 2 - 1;
  /* speed is needed, if
   * 1. true position is being computed before applying light-time etc.
   *    this is the position saved in pdp->x.
   *    in this case, speed is needed for light-time correction.
   * 2. the speed flag has been specified.
   */
  need_speed = (do_save || is_set(iflag,SEFLG_SPEED));
  for (i = 0; i <= 2; i++) {
    xp[i]  = swi_echeb (t, pdp->segp+(i*pdp->ncoe), pdp->neval);
    if (need_speed)
      xp[i+3] = swi_edcheb(t, pdp->segp+(i*pdp->ncoe), pdp->neval) / pdp->dseg * 2;
    else
      xp[i+3] = 0;  /* von Alois als billiger fix, evtl. illegal */
  }
  /* if planet wanted is barycentric sun and must be computed
   * from heliocentric earth and barycentric earth: the
   * computation above gives heliocentric earth, therefore we
   * have to compute barycentric earth and subtract heliocentric
   * earth from it. this may be necessary with calls from
   * sweplan() and from app_pos_etc_sun() (light-time). */
  if (ipl == SEI_SUNBARY && (pdp->iflg & SEI_FLG_EMBHEL)) {
    /* sweph() calls sweph() !!! for EMB.
     * Attention: a new calculation must be forced in any case.
     * Otherwise EARTH (instead of EMB) will possibly taken from
     * save area.
     * to force new computation, set pedp->teval = 0 and restore it
     * after call of sweph(EMB).
     */
    tsv = pedp->teval;
    pedp->teval = 0;
    retc = sweph(tjd, SEI_EMB, ifno, iflag | SEFLG_SPEED, NULL, NO_SAVE, xemb, serr, swed);
    if (retc != OK)
      return(retc);
    pedp->teval = tsv;
    // xp[] = xemb[] - xp[]
    int n = need_speed ? 6 : 3;    
    vec_add( smul(xp,-1,n), xemb, n);
    }
    
  /* asteroids are heliocentric.
   * if JPL or SWISSEPH, convert to barycentric */
  if (is_set(iflag,SEFLG_JPLEPH | SEFLG_SWIEPH)) {
    if (ipl >= SEI_ANYBODY) {
      int n = need_speed ? 6 : 3;    
      vec_add(xp, xsunb, n);
      }
    }
    
  if (do_save) {
    pdp->teval = tjd;
    pdp->xflgs = -1;  /* do new computation of light-time etc. */
    if (ifno == SEI_FILE_PLANET || ifno == SEI_FILE_MOON)
      pdp->iephe = SEFLG_SWIEPH;/**/
    else
      pdp->iephe = psdp->iephe;
  }
  vec_copy(xpret,xp,6);
  return(OK);

  }

/*
 * Alois 2.12.98: inserted error message generation for file not found
 */
FILE *swi_fopen(int ifno, char *fname, char *ephepath, char *serr) {
  return swid_fopen(ifno,fname,ephepath,serr,&swed);
  }
FILE *swid_fopen(int ifno, char *fname, char *ephepath, char *serr, struct swe_data *swed)
{
  int np, i, j;
  FILE *fp = NULL;
  char *fnamp, fn[AS_MAXCH];
  char *cpos[20];
  char s[2 * AS_MAXCH], *s1 = s + AS_MAXCH;  /* a little trick */
  if (ifno >= 0) {
    fnamp = swed->fidat[ifno].fnam;
  } else {
    fnamp = fn;
  }
  strcpy(s1, ephepath);
  np = swi_cutstr(s1, PATH_SEPARATOR, cpos, 20);
  for (i = 0; i < np; i++) {
    strcpy(s, cpos[i]);
    if (strcmp(s, ".") == 0)  /* current directory */
      *s = '\0';
    else {
      j = strlen(s);
      if (*(s + j - 1) != *DIR_GLUE && *s != '\0')
  strcat(s, DIR_GLUE);
    }
    strcpy(fnamp, s);
    if (strlen(fnamp) + strlen(fname) < AS_MAXCH)
      strcat(fnamp, fname);
    else {
      if (serr != NULL)
  sprintf(serr, "error: file path and name must be shorter than %d.", AS_MAXCH);
      return NULL;
    }
    fp = fopen(fnamp, BFILE_R_ACCESS);
    if (fp != NULL)
      return fp;
  }
  sprintf(s, "SwissEph file '%s' not found in PATH '%s'", fname, ephepath);
  s[AS_MAXCH-1] = '\0';    /* s may be longer then AS_MAXCH */
  throw( ERR, s, serr, swed );
  return NULL;  // will never be reached
}

/* converts planets from barycentric to geocentric,
 * apparent positions
 * precession and nutation
 * according to flags
 * ipli    planet number
 * iflag  flags
 * serr         error string
 */
static int app_pos_etc_plan(int ipli, int iflag, char *serr, struct swe_data *swed)
{
  int i, j, niter, retc = OK;
  int ifno, ibody;
  int flg1, flg2;
  double xx[6], dx[3], dt, t, dtsave_for_defl;
  double xobs[6], xobs2[6];
  double xearth[6], xsun[6];
  double xxsp[6], xxsv[6];
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *pdp;
  struct epsilon *oe = &swed->oec2000;
  int epheflag = iflag & SEFLG_EPHMASK;
  t = dtsave_for_defl = 0;  /* dummy assignment to silence gcc */
  /* ephemeris file */
  if (ipli > SE_AST_OFFSET) {
    ifno = SEI_FILE_ANY_AST;
    ibody = IS_ANY_BODY;
    pdp = &swed->pldat[SEI_ANYBODY];
  } else if (ipli == SEI_CHIRON
      || ipli == SEI_PHOLUS
      || ipli == SEI_CERES
      || ipli == SEI_PALLAS
      || ipli == SEI_JUNO
      || ipli == SEI_VESTA) {
    ifno = SEI_FILE_MAIN_AST;
    ibody = IS_MAIN_ASTEROID;
    pdp = &swed->pldat[ipli];
  } else {
    ifno = SEI_FILE_PLANET;
    ibody = IS_PLANET;
    pdp = &swed->pldat[ipli];
  }
  /* if the same conversions have already been done for the same
   * date, then return */
  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pdp->xflgs = iflag;
    pdp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }
  /* the conversions will be done with xx[]. */
  vec_copy(xx,pdp->x,6);
  
  /* if heliocentric position is wanted */
  if (is_set(iflag, SEFLG_HELCTR)) {
    if (pdp->iephe == SEFLG_JPLEPH || pdp->iephe == SEFLG_SWIEPH)
      vec_sub(xx,swed->pldat[SEI_SUNBARY].x,6);
    }
    
  /************************************
   * observer: geocenter or topocenter
   ************************************/
  /* if topocentric position is wanted  */
  if (is_set(iflag, SEFLG_TOPOCTR)) {    
    if (swed->topd.teval != pedp->teval
      || swed->topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag, DO_SAVE, xobs, serr) != OK)
        return ERR;
      } else {
      vec_copy(xobs,swed->topd.xobs,6);
      }  
    /* barycentric position of observer */
    vec_add(xobs,pedp->x,6);
    
    } 
  else {
    /* barycentric position of geocenter */
    vec_copy(xobs,pedp->x,6);
    }
  /*******************************
   * light-time geocentric       *
   *******************************/
  if (is_not_set(iflag,SEFLG_TRUEPOS)) {
    /* number of iterations - 1 */
    if (pdp->iephe == SEFLG_JPLEPH || pdp->iephe == SEFLG_SWIEPH)
      niter = 1;
    else   /* SEFLG_MOSEPH or planet from osculating elements */
      niter = 0;
    if (is_set(iflag,SEFLG_SPEED)) {
      /*
       * Apparent speed is influenced by the fact that dt changes with
       * motion. This makes a difference of several hundredths of an
       * arc second. To take this into account, we compute
       * 1. true position - apparent position at time t - 1.
       * 2. true position - apparent position at time t.
       * 3. the difference between the two is the part of the daily motion
       * that results from the change of dt.
       */
      for (i = 0; i <= 2; i++)
        xxsv[i] = xxsp[i] = xx[i] - xx[i+3];

      for (j = 0; j <= niter; j++) {
        for (i = 0; i <= 2; i++) {
          dx[i] = xxsp[i];
          if (is_not_set(iflag,SEFLG_HELCTR | SEFLG_BARYCTR))
            dx[i] -= (xobs[i] - xobs[i+3]);
          }
  /* new dt */
        dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;        
        for (i = 0; i <= 2; i++)   /* rough apparent position at t-1 */
          xxsp[i] = xxsv[i] - dt * pdp->x[i+3];
        }
  /* true position - apparent position at time t-1 */
      for (i = 0; i <= 2; i++)
        xxsp[i] = xxsv[i] - xxsp[i];
      }
    /* dt and t(apparent) */
    for (j = 0; j <= niter; j++) {
      for (i = 0; i <= 2; i++) {
        dx[i] = xx[i];
        if (is_not_set(iflag,SEFLG_HELCTR | SEFLG_BARYCTR ))
          dx[i] -= xobs[i];
        }
      dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
      /* new t */
      t = pdp->teval - dt;
      dtsave_for_defl = dt;
      for (i = 0; i <= 2; i++)     /* rough apparent position at t*/
        xx[i] = pdp->x[i] - dt * pdp->x[i+3];
      }
    /* part of daily motion resulting from change of dt */
    if (is_set(iflag,SEFLG_SPEED))
      for (i = 0; i <= 2; i++)
        xxsp[i] = pdp->x[i] - xx[i] - xxsp[i];
        
    /* new position, accounting for light-time (accurate) */
    switch(epheflag) {
      case SEFLG_SWIEPH:
        if (ibody == IS_PLANET)
          retc = sweplan(t, ipli, ifno, iflag, NO_SAVE, xx, xearth, xsun, NULL, serr, swed);
        else {     /*asteroid*/
          retc = sweplan(t, SEI_EARTH, SEI_FILE_PLANET, iflag, NO_SAVE, xearth, NULL, xsun, NULL, serr, swed);        
          if (retc == OK) retc = sweph(t, ipli, ifno, iflag, xsun, NO_SAVE, xx, serr, swed);
          }
        if (retc != OK) return(retc);
        break;
      case SEFLG_MOSEPH:
      default:
  /*
   * with moshier or other ephemerides, subtraction of dt * speed
   * is sufficient (has been done in light-time iteration above)
   */
 /* if speed flag is true, we call swi_moshplan() for new t.
   * this does not increase position precision,
   * but speed precision, which becomes better than 0.01"/day.
   * for precise speed, we need earth as well.
   */
        if (iflag & SEFLG_SPEED
          && !(iflag & (SEFLG_HELCTR | SEFLG_BARYCTR))) {
          if (ibody == IS_PLANET)
            retc = swi_moshplan(t, ipli, NO_SAVE, xxsv, xearth, serr);
          else {    /* if asteroid */
            retc = sweph(t, ipli, ifno, iflag, NULL, NO_SAVE, xxsv, serr, swed);
            if (retc == OK)
              retc = swi_moshplan(t, SEI_EARTH, NO_SAVE, xearth, xearth, serr);
            }
          if (retc != OK) return(retc);
        /* only speed is taken from this computation, otherwise position
         * calculations with and without speed would not agree. The difference
         * would be about 0.01", which is far below the intrinsic error of the
         * moshier ephemeris.
         */
          vec_copy(xx+3,xxsv+3,3); 
          }
        break;
      }
      
    if (is_set(iflag, SEFLG_HELCTR)) {
      if (pdp->iephe == SEFLG_JPLEPH || pdp->iephe == SEFLG_SWIEPH)
        vec_sub( xx, swed->pldat[SEI_SUNBARY].x, 6);
      }
      
    if (is_set(iflag,SEFLG_SPEED)) {
      /* observer position for t(light-time) */
      if (is_set(iflag,SEFLG_TOPOCTR)) {
        if (swi_get_observer(t, iflag, NO_SAVE, xobs2, serr) != OK) return ERR;
        vec_add(xobs2,xearth,6);
      } else {
        vec_copy(xobs2,xearth,6);
      }
    }
  }
  /*******************************
   * conversion to geocenter     *
   *******************************/
  if (is_not_set(iflag, SEFLG_HELCTR | SEFLG_BARYCTR)) {
    /* subtract earth */
    vec_sub(xx,xobs,6);
    if (is_not_set(iflag, SEFLG_TRUEPOS)) {
      /*
       * Apparent speed is also influenced by
       * the change of dt during motion.
       * Neglect of this would result in an error of several 0.01"
       */
      if (is_set(iflag,SEFLG_SPEED))
        vec_sub(xx+3,xxsp,3);
      }
    }
  if (is_not_set(iflag, SEFLG_SPEED)) vec_clear(xx+3,3);
    
  /************************************
   * relativistic deflection of light *
   ************************************/
  if (is_not_set(iflag, SEFLG_TRUEPOS | SEFLG_NOGDEFL))
    /* SEFLG_NOGDEFL is on, if SEFLG_HELCTR or SEFLG_BARYCTR */
    swi_deflect_light(xx, dtsave_for_defl, iflag);
    
  /**********************************
   * 'annual' aberration of light   *
   **********************************/
  if (is_not_set(iflag, SEFLG_TRUEPOS | SEFLG_NOABERR)) {
    /* SEFLG_NOABERR is on, if SEFLG_HELCTR or SEFLG_BARYCTR */
    swi_aberr_light(xx, xobs, iflag);
    /*
     * Apparent speed is also influenced by
     * the difference of speed of the earth between t and t-dt.
     * Neglecting this would involve an error of several 0.1"
     */
    if (is_set(iflag,SEFLG_SPEED))
      vec_sub( 
        vec_add(xx+3, xobs+3, 3), 
        xobs2+3, 
        3);
    }
    
  if (is_not_set(iflag, SEFLG_SPEED)) vec_clear(xx+3,3);

  /* ICRS to J2000 */
  if (is_not_set(iflag, SEFLG_ICRS) && swed->jpldenum >= 403) {
    swi_bias(xx, iflag, FALSE);
  }/**/
  
  /* save J2000 coordinates; required for sidereal positions */
  vec_copy(xxsv,xx,6);

  /************************************************
   * precession, equator 2000 -> equator of date *
   ************************************************/
  if (is_not_set(iflag, SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, J2000_TO_J);
    if (is_set(iflag,SEFLG_SPEED)) swi_precess_speed(xx, pdp->teval, J2000_TO_J);
    oe = &swed->oec;
  } else
    oe = &swed->oec2000;
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr, swed);
}

static int app_pos_rest(struct plan_data *pdp, int iflag,
                        double *xx, double *x2000,
                        struct epsilon *oe, char *serr, struct swe_data *swed)
{

  /************************************************
   * nutation                                     *
   ************************************************/
  if (is_not_set(iflag, SEFLG_NONUT))
    swi_nutate(xx, iflag, FALSE);
    
  /* now we have equatorial cartesian coordinates; save them */
    vec_copy(pdp->xreturn+18,xx,6);

  /************************************************
   * transformation to ecliptic.                  *
   * with sidereal calc. this will be overwritten *
   * afterwards.                                  *
   ************************************************/
  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  if (is_set(iflag, SEFLG_SPEED))
    swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
  if (is_not_set(iflag, SEFLG_NONUT)) {
    swi_coortrf2(xx, xx, swed->nut.snut, swed->nut.cnut);
    if (is_set(iflag, SEFLG_SPEED))
      swi_coortrf2(xx+3, xx+3, swed->nut.snut, swed->nut.cnut);
  }
  
  /* now we have ecliptic cartesian coordinates */
  vec_copy(pdp->xreturn+6,xx,6);

  /************************************
   * sidereal positions               *
   ************************************/
  if (is_set(iflag,SEFLG_SIDEREAL)) {
    /* project onto ecliptic t0 */
    if (swed->sidd.sid_mode & SE_SIDBIT_ECL_T0) {
      if (swi_trop_ra2sid_lon(x2000, pdp->xreturn+6, pdp->xreturn+18, iflag, serr) != OK)
  return ERR;
    /* project onto solar system equator */
    } else if (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
      if (swi_trop_ra2sid_lon_sosy(x2000, pdp->xreturn+6, pdp->xreturn+18, iflag, serr) != OK)
  return ERR;
    } else {
    /* traditional algorithm */
      swi_cartpol_sp(pdp->xreturn+6, pdp->xreturn);
      pdp->xreturn[0] -= swe_get_ayanamsa(pdp->teval) * DEGTORAD;
      swi_polcart_sp(pdp->xreturn, pdp->xreturn+6);
    }
  }
  /************************************************
   * transformation to polar coordinates          *
   ************************************************/
  swi_cartpol_sp(pdp->xreturn+18, pdp->xreturn+12);
  swi_cartpol_sp(pdp->xreturn+6, pdp->xreturn);
  
  /**********************
   * radians to degrees *
   **********************/
  smul( pdp->xreturn,    RADTODEG, 2);
  smul( pdp->xreturn+3,  RADTODEG, 2);
  smul( pdp->xreturn+12, RADTODEG, 2);
  smul( pdp->xreturn+15, RADTODEG, 2);

  /* save, what has been done */
  pdp->xflgs = iflag;
  pdp->iephe = iflag & SEFLG_EPHMASK;
  return OK;
}

void FAR PASCAL_CONV swed_set_sid_mode(int sid_mode, double t0, double ayan_t0, struct swe_data *swed)
{
  struct sid_data *sip = &swed->sidd;
  sip->sid_mode = sid_mode;
  if (sid_mode >= SE_SIDBITS)
    sid_mode %= SE_SIDBITS;
  /* standard equinoxes: positions always referred to ecliptic of t0 */
  if (sid_mode == SE_SIDM_J2000
    || sid_mode == SE_SIDM_J1900
    || sid_mode == SE_SIDM_B1950) {
    sip->sid_mode |= SE_SIDBIT_ECL_T0;
  }
  if (sid_mode >= SE_NSIDM_PREDEF && sid_mode != SE_SIDM_USER) {
    sip->sid_mode = sid_mode = SE_SIDM_FAGAN_BRADLEY;
  }
  swed->ayana_is_set = TRUE;
  if (sid_mode == SE_SIDM_USER) {
    sip->t0 = t0;
    sip->ayan_t0 = ayan_t0;
  } else {
    sip->t0 = ayanamsa[sid_mode].t0;
    sip->ayan_t0 = ayanamsa[sid_mode].ayan_t0;
  }
  swi_force_app_pos_etc(swed);
}

/* the ayanamsa (precession in longitude)
 * according to Newcomb's definition: 360 -
 * longitude of the vernal point of t referred to the
 * ecliptic of t0.
 */
double FAR PASCAL_CONV swed_get_ayanamsa(double tjd_et, struct swe_data *swed)
{
  double x[6], eps;
  struct sid_data *sip = &swed->sidd;
  if (!swed->ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);
  /* vernal point (tjd), cartesian */
  x[0] = 1;
  x[1] = x[2] = 0;
  /* to J2000 */
  if (tjd_et != J2000)
    swi_precess(x, tjd_et, J_TO_J2000);
  /* to t0 */
  swi_precess(x, sip->t0, J2000_TO_J);
  /* to ecliptic */
  eps = swi_epsiln(sip->t0);
  swi_coortrf(x, x, eps);
  /* to polar */
  swi_cartpol(x, x);
  /* subtract initial value of ayanamsa */
  x[0] = x[0] * RADTODEG - sip->ayan_t0;
  /* get ayanamsa */
  return swe_degnorm(-x[0]);
}

/*
 * input coordinates are J2000, cartesian.
 * xout   ecliptical sidereal position
 * xoutr   equatorial sidereal position
 */
int swi_trop_ra2sid_lon(double *xin, double *xout, double *xoutr, int iflag, char *serr)
{
  return swid_trop_ra2sid_lon(xin,xout,xoutr,iflag,serr,&swed);
}

int swid_trop_ra2sid_lon(double *xin, double *xout, double *xoutr, int iflag, char *serr, struct swe_data *swed)
{
  double x[6];

  struct sid_data *sip = &swed->sidd;
  struct epsilon oectmp;
  
  vec_copy(x,xin,6);

  if (sip->t0 != J2000) {
    swi_precess(x, sip->t0, J2000_TO_J);
    swi_precess(x+3, sip->t0, J2000_TO_J);  /* speed */
    }
  
  vec_copy(xoutr,x,6);

  calc_epsilon(swed->sidd.t0, &oectmp);
  swi_coortrf2(x, x, oectmp.seps, oectmp.ceps);
  if (is_set(iflag, SEFLG_SPEED))
    swi_coortrf2(x+3, x+3, oectmp.seps, oectmp.ceps);
    
  /* to polar coordinates */
  swi_cartpol_sp(x, x);
  
  /* subtract ayan_t0 */  
  x[0] -= sip->ayan_t0 * DEGTORAD;
  
  /* back to cartesian */
  swi_polcart_sp(x, xout);
  return OK;

  }

/*
 * input coordinates are J2000, cartesian.
 * xout   ecliptical sidereal position
 * xoutr   equatorial sidereal position
 */
int swi_trop_ra2sid_lon_sosy(double *xin, double *xout, double *xoutr, int iflag, char *serr)
{
  return swid_trop_ra2sid_lon_sosy(xin,xout,xoutr,iflag,serr,&swed); 
}

int swid_trop_ra2sid_lon_sosy(double *xin, double *xout, double *xoutr, int iflag, char *serr, struct swe_data *swed)
{
  double x[6], x0[6];

  struct sid_data *sip = &swed->sidd;
  struct epsilon *oe = &swed->oec2000;
  double plane_node = SSY_PLANE_NODE_E2000;
  double plane_incl = SSY_PLANE_INCL;
  
  vec_copy(x,xin,6);

  /* planet to ecliptic 2000 */
  swi_coortrf2(x, x, oe->seps, oe->ceps);
  if (is_set(iflag,SEFLG_SPEED))
    swi_coortrf2(x+3, x+3, oe->seps, oe->ceps);
  /* to polar coordinates */
  swi_cartpol_sp(x, x);
  /* to solar system equator */
  x[0] -= plane_node;
  swi_polcart_sp(x, x);
  swi_coortrf(x, x, plane_incl);
  swi_coortrf(x+3, x+3, plane_incl);
  swi_cartpol_sp(x, x);
  /* zero point of t0 in J2000 system */
  x0[0] = 1;
  x0[1] = x0[2] = 0;
  if (sip->t0 != J2000)
    swi_precess(x0, sip->t0, J_TO_J2000);
  /* zero point to ecliptic 2000 */
  swi_coortrf2(x0, x0, oe->seps, oe->ceps);
  /* to polar coordinates */
  swi_cartpol(x0, x0);
  /* to solar system equator */
  x0[0] -= plane_node;
  swi_polcart(x0, x0);
  swi_coortrf(x0, x0, plane_incl);
  swi_cartpol(x0, x0);
  /* measure planet from zero point */
  x[0] -= x0[0];
  x[0] *= RADTODEG;
  /* subtract ayan_t0 */
  x[0] -= sip->ayan_t0;
  x[0] = swe_degnorm(x[0]) * DEGTORAD;
  /* back to cartesian */
  swi_polcart_sp(x, xout);
  return OK;
}

/* converts planets from barycentric to geocentric,
 * apparent positions
 * precession and nutation
 * according to flags
 * ipli    planet number
 * iflag  flags
 */
static int app_pos_etc_plan_osc(int ipl, int ipli, int iflag, char *serr, struct swe_data *swed)
{
  int i, j, niter, retc;
  double xx[6], dx[3], dt, dtsave_for_defl;
  double xearth[6], xsun[6], xmoon[6];
  double xxsv[6], xxsp[6], xobs[6], xobs2[6];
  double t;
  struct plan_data *pdp = &swed->pldat[ipli];
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed->oec2000;
  int epheflag = SEFLG_DEFAULTEPH;
  dt = dtsave_for_defl = 0;  /* dummy assign to silence gcc */
  if (is_set(iflag,SEFLG_MOSEPH))
    epheflag = SEFLG_MOSEPH;
  else if (is_set(iflag,SEFLG_SWIEPH))
    epheflag = SEFLG_SWIEPH;
  else if (is_set(iflag,SEFLG_JPLEPH))
    epheflag = SEFLG_JPLEPH;
  /* the conversions will be done with xx[]. */
  vec_copy(xx,pdp->x,6);

  /************************************
   * barycentric position is required *
   ************************************/
  /* = heliocentric position with Moshier ephemeris */
  /************************************
   * observer: geocenter or topocenter
   ************************************/
  /* if topocentric position is wanted  */
  if (is_set(iflag, SEFLG_TOPOCTR)) {
    if (swed->topd.teval != pedp->teval
      || swed->topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag, DO_SAVE, xobs, serr) != OK)
        return ERR;
        } 
      else {
        vec_copy(xobs,swed->topd.xobs,6);
        }
      /* barycentric position of observer */
    vec_add(xobs,pedp->x,6);
    } 
  else 
    if (is_set(iflag, SEFLG_BARYCTR)) {
      vec_clear(xobs,6); 
      } 
    else if (is_set(iflag, SEFLG_HELCTR)) {
      if (is_set(iflag, SEFLG_MOSEPH)) {
        vec_clear(xobs,6);
        } 
      else {
        vec_copy(xobs,psdp->x,6);
        } 
      }  
    else {
      vec_copy( xobs, pedp->x, 6);
      }
    
  
  /*******************************
   * light-time                  *
   *******************************/
  if (is_not_set(iflag, SEFLG_TRUEPOS)) {
    niter = 1;
    if (is_set(iflag, SEFLG_SPEED)) {
      /*
       * Apparent speed is influenced by the fact that dt changes with
       * motion. This makes a difference of several hundredths of an
       * arc second. To take this into account, we compute
       * 1. true position - apparent position at time t - 1.
       * 2. true position - apparent position at time t.
       * 3. the difference between the two is the daily motion resulting from
       * the change of dt.
       */
       
      for (i = 0; i <= 2; i++)
        xxsv[i] = xxsp[i] = xx[i] - xx[i+3];
      for (j = 0; j <= niter; j++) {
        vec_copy(dx,xxsp,3);
        if (is_not_set(iflag, SEFLG_HELCTR | SEFLG_BARYCTR))
          vec_sub(vec_add(dx,xobs+3,3),xobs,3);
          
  /* new dt */
        dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
        for (i = 0; i <= 2; i++)
          xxsp[i] = xxsv[i] - dt * pdp->x[i+3];/* rough apparent position */
        }
        
      /* true position - apparent position at time t-1 */
      for (i = 0; i <= 2; i++) xxsp[i] = xxsv[i] - xxsp[i];
    }
    /* dt and t(apparent) */
    for (j = 0; j <= niter; j++) {
      vec_copy(dx,xx,3);
      if (is_not_set( iflag, SEFLG_HELCTR | SEFLG_BARYCTR )) 
        vec_sub(dx,xobs,3);
        
      /* new dt */
        dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
        dtsave_for_defl = dt;
      /* new position: subtract t * speed
       */
        for (i = 0; i <= 2; i++) {
          xx[i] = pdp->x[i] - dt * pdp->x[i+3];/**/
          xx[i+3] = pdp->x[i+3];
        }
    }
    if (is_set(iflag, SEFLG_SPEED)) {
      /* part of daily motion resulting from change of dt */
      for (i = 0; i <= 2; i++)
        xxsp[i] = pdp->x[i] - xx[i] - xxsp[i];
      
      t = pdp->teval - dt;
      /* for accuracy in speed, we will need earth as well */
      retc = main_planet_bary(t, SEI_EARTH, epheflag, iflag, NO_SAVE, xearth, xearth, xsun, xmoon, serr, swed);
      if (swi_osc_el_plan(t, xx, ipl-SE_FICT_OFFSET, ipli, xearth, xsun, serr) != OK)  return ERR;
      if (retc != OK)  return(retc);
      
      if (is_set(iflag, SEFLG_TOPOCTR)) {
        if (swi_get_observer(t, iflag, NO_SAVE, xobs2, serr) != OK) return ERR;
        vec_add(xobs2,xearth,6);
      } else {
        vec_copy(xobs2,xearth,6);
      }      
    }    
  }
  
  /*******************************
   * conversion to geocenter     *
   *******************************/
  vec_sub(xx,xobs,6);

  if (is_not_set(iflag, SEFLG_TRUEPOS)) {
    /*
     * Apparent speed is also influenced by
     * the change of dt during motion.
     * Neglect of this would result in an error of several 0.01"
     */
    if (is_set(iflag, SEFLG_SPEED))
      vec_sub(xx+3,xxsp,3);
    }
  
  if (is_not_set(iflag, SEFLG_SPEED))
    vec_clear(xx+3,3);

  /************************************
   * relativistic deflection of light *
   ************************************/
  if (is_not_set(iflag, SEFLG_TRUEPOS | SEFLG_NOGDEFL))
  /* SEFLG_NOGDEFL is on, if SEFLG_HELCTR or SEFLG_BARYCTR */
    swi_deflect_light(xx, dtsave_for_defl, iflag);
  
  /**********************************
   * 'annual' aberration of light   *
   **********************************/
  if (is_not_set(iflag, SEFLG_TRUEPOS | SEFLG_NOABERR)) {
    /* SEFLG_NOABERR is on, if SEFLG_HELCTR or SEFLG_BARYCTR */
    swi_aberr_light(xx, xobs, iflag);    
    /*
     * Apparent speed is also influenced by
     * the difference of speed of the earth between t and t-dt.
     * Neglecting this would involve an error of several 0.1"
     */
    if (is_set(iflag, SEFLG_SPEED))
      for (i = 3; i <= 5; i++)
        xx[i] += xobs[i] - xobs2[i];
    }
    
  
  /* save J2000 coordinates; required for sidereal positions */
  vec_copy(xxsv,xx,6);

  /************************************************
   * precession, equator 2000 -> equator of date *
   ************************************************/
  if (is_not_set(iflag, SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, J2000_TO_J);
    if (is_set(iflag, SEFLG_SPEED))
      swi_precess_speed(xx, pdp->teval, J2000_TO_J);
    oe = &swed->oec;
    } 
  else
    oe = &swed->oec2000;
  
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr, swed);

  }

/* influence of precession on speed
 * xx    position and speed of planet in equatorial cartesian
 *    coordinates */
void swi_precess_speed(double *xx, double t, int direction)
{
  swid_precess_speed(xx,t,direction,&swed);
}

void swid_precess_speed(double *xx, double t, int direction, struct swe_data *swed)
{
  struct epsilon *oe;
  double fac;
  double tprec = (t - J2000) / 36525.0;
  if (direction == J2000_TO_J) {
    fac = 1;
    oe = &swed->oec;
  } else {
    fac = -1;
    oe = &swed->oec2000;
  }
  /* first correct rotation.
   * this costs some sines and cosines, but neglect might
   * involve an error > 1"/day */
  swi_precess(xx+3, t, direction);
  /* then add 0.137"/day */
  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
  swi_cartpol_sp(xx, xx);
  xx[3] += (50.290966 + 0.0222226 * tprec) / 3600 / 365.25 * DEGTORAD * fac;
      /* formula from Montenbruck, German 1994, p. 18 */
  swi_polcart_sp(xx, xx);
  swi_coortrf2(xx, xx, -oe->seps, oe->ceps);
  swi_coortrf2(xx+3, xx+3, -oe->seps, oe->ceps);
}

/* multiplies cartesian equatorial coordinates with previously
 * calculated nutation matrix. also corrects speed.
 */
void swi_nutate(double *xx, int iflag, AS_BOOL backward)
{
  swid_nutate(xx,iflag,backward,&swed);
}

void swid_nutate(double *xx, int iflag, AS_BOOL backward, struct swe_data *swed)
{
  int i;
  double x[6], xv[6];
  for (i = 0; i <= 2; i++) {
    if (backward)
      x[i] = xx[0] * swed->nut.matrix[i][0] +
       xx[1] * swed->nut.matrix[i][1] +
       xx[2] * swed->nut.matrix[i][2];
    else
      x[i] = xx[0] * swed->nut.matrix[0][i] +
       xx[1] * swed->nut.matrix[1][i] +
       xx[2] * swed->nut.matrix[2][i];
  }
  if (is_set(iflag,SEFLG_SPEED)) {
    /* correct speed:
     * first correct rotation */
    for (i = 0; i <= 2; i++) {
      if (backward)
  x[i+3] = xx[3] * swed->nut.matrix[i][0] +
     xx[4] * swed->nut.matrix[i][1] +
     xx[5] * swed->nut.matrix[i][2];
      else
  x[i+3] = xx[3] * swed->nut.matrix[0][i] +
     xx[4] * swed->nut.matrix[1][i] +
     xx[5] * swed->nut.matrix[2][i];
    }
    /* then apparent motion due to change of nutation during day.
     * this makes a difference of 0.01" */
    for (i = 0; i <= 2; i++) {
      if (backward)
  xv[i] = xx[0] * swed->nutv.matrix[i][0] +
         xx[1] * swed->nutv.matrix[i][1] +
         xx[2] * swed->nutv.matrix[i][2];
      else
  xv[i] = xx[0] * swed->nutv.matrix[0][i] +
         xx[1] * swed->nutv.matrix[1][i] +
         xx[2] * swed->nutv.matrix[2][i];
      /* new speed */
      xx[3+i] = x[3+i] + (x[i] - xv[i]) / NUT_SPEED_INTV;
    }
  }
  /* new position */
  vec_copy(xx,x,3);
}

/* computes 'annual' aberration
 * xx    planet's position accounted for light-time
 *              and gravitational light deflection
 * xe      earth's position and speed
 */
void swi_aberr_light(double *xx, double *xe, int iflag) {
  int i;
  double xxs[6], v[6], u[6], ru;
  double xx2[6], dx1, dx2;
  double b_1, f1, f2;
  double v2;
  double intv = PLAN_SPEED_INTV;
  
  vec_copy(xxs,xx,6);
  vec_copy(u,xx,6);

  ru = vec_length(u,3);
  for (i = 0; i <= 2; i++)
    v[i] = xe[i+3] / 24.0 / 3600.0 / CLIGHT * AUNIT;
  v2 = square_sum(v);
  b_1 = sqrt(1 - v2);
  f1 = dot_prod(u, v) / ru;
  f2 = 1.0 + f1 / (1.0 + b_1);
  for (i = 0; i <= 2; i++)
    xx[i] = (b_1*xx[i] + f2*ru*v[i]) / (1.0 + f1);
  if (is_set(iflag,SEFLG_SPEED)) {
    /* correction of speed
     * the influence of aberration on apparent velocity can
     * reach 0.4"/day
     */
    for (i = 0; i <= 2; i++)
      u[i] = xxs[i] - intv * xxs[i+3];
    ru = vec_length(u,3);
    f1 = dot_prod(u, v) / ru;
    f2 = 1.0 + f1 / (1.0 + b_1);
    for (i = 0; i <= 2; i++)
      xx2[i] = (b_1*u[i] + f2*ru*v[i]) / (1.0 + f1);
    for (i = 0; i <= 2; i++) {
      dx1 = xx[i] - xxs[i];
      dx2 = xx2[i] - u[i];
      dx1 -= dx2;
      xx[i+3] += dx1 / intv;
    }
  }
}

/* computes relativistic light deflection by the sun
 * ipli   sweph internal planet number
 * xx    planet's position accounted for light-time
 * dt    dt of light-time
 */
void swi_deflect_light(double *xx, double dt, int iflag)
{
  swid_deflect_light(xx,dt,iflag,&swed);
}

void swid_deflect_light(double *xx, double dt, int iflag, struct swe_data *swed)
{
  int i;
  double xx2[6];
  double u[6], e[6], q[6], ru, re, uq, ue, qe, g1, g2;
  double xx3[6], dx1, dx2, dtsp;
  double xsun[6], xearth[6];
  double sina, sin_sunr, meff_fact;
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
  int iephe = pedp->iephe;
  vec_copy(xearth,pedp->x,6);
  if (is_set(iflag, SEFLG_TOPOCTR))
    vec_add(xearth,swed->topd.xobs,6);

  /* U = planetbary(t-tau) - earthbary(t) = planetgeo */
  vec_copy(u,xx,3);

  /* Eh = earthbary(t) - sunbary(t) = earthhel */
  if (iephe == SEFLG_JPLEPH || iephe == SEFLG_SWIEPH)
    for (i = 0; i <= 2; i++)
      e[i] = xearth[i] - psdp->x[i];
  else
    vec_copy(e,xearth,3);

  /* Q = planetbary(t-tau) - sunbary(t-tau) = 'planethel' */
  /* first compute sunbary(t-tau) for */
  if (iephe == SEFLG_JPLEPH || iephe == SEFLG_SWIEPH) {
    for (i = 0; i <= 2; i++)
      /* this is sufficient precision */
      xsun[i] = psdp->x[i] - dt * psdp->x[i+3];
    vec_copy(xsun+3,psdp->x+3,3);
  } else {
    vec_copy(xsun,psdp->x,6);
  }
  for (i = 0; i <= 2; i++)
    q[i] = xx[i] + xearth[i] - xsun[i];

    
  ru = vec_normalize(u,3);
       vec_normalize(q,3);
  re = vec_normalize(e,3);

  uq = dot_prod(u,q);
  ue = dot_prod(u,e);
  qe = dot_prod(q,e);

  /* When a planet approaches the center of the sun in superior
   * conjunction, the formula for the deflection angle as given
   * in Expl. Suppl. p. 136 cannot be used. The deflection seems
   * to increase rapidly towards infinity. The reason is that the
   * formula considers the sun as a point mass. AA recommends to
   * set deflection = 0 in such a case.
   * However, to get a continous motion, we modify the formula
   * for a non-point-mass, taking into account the mass distribution
   * within the sun. For more info, s. meff().
   */
  sina = sqrt(1 - ue * ue);  /* sin(angle) between sun and planet */
  sin_sunr = SUN_RADIUS / re;   /* sine of sun radius (= sun radius) */
  if (sina < sin_sunr)
    meff_fact = meff(sina / sin_sunr);
  else
    meff_fact = 1;
  g1 = 2.0 * HELGRAVCONST * meff_fact / CLIGHT / CLIGHT / AUNIT / re;
  g2 = 1.0 + qe;
  /* compute deflected position */
  for (i = 0; i <= 2; i++)
    xx2[i] = ru * (u[i] + g1/g2 * (uq * e[i] - ue * q[i]));
  if (is_set(iflag,SEFLG_SPEED)) {
    /* correction of speed
     * influence of light deflection on a planet's apparent speed:
     * for an outer planet at the solar limb with
     * |v(planet) - v(sun)| = 1 degree, this makes a difference of 7"/day.
     * if the planet is within the solar disc, the difference may increase
     * to 30" or more.
     * e.g. mercury at j2434871.45:
     *  distance from sun     45"
     *  1. speed without deflection     2d10'10".4034
     *    2. speed with deflection        2d10'42".8460 (-speed flag)
     *    3. speed with deflection        2d10'43".4824 (< 3 positions/
     *                 -speed3 flag)
     * 3. is not very precise. Smaller dt would give result closer to 2.,
     * but will probably never be as good as 2, unless int doubles are
     * used. (try also j2434871.46!!)
     * however, in such a case speed changes rapidly. before being
     * passed by the sun, the planet accelerates, and after the sun
     * has passed it slows down. some time later it regains 'normal'
     * speed.
     * to compute speed, we do the same calculation as above with
     * slightly different u, e, q, and find out the difference in
     * deflection.
     */
    dtsp = -DEFL_SPEED_INTV;
    /* U = planetbary(t-tau) - earthbary(t) = planetgeo */
    for (i = 0; i <= 2; i++)
      u[i] = xx[i] - dtsp * xx[i+3];
    /* Eh = earthbary(t) - sunbary(t) = earthhel */
    if (iephe == SEFLG_JPLEPH || iephe == SEFLG_SWIEPH) {
      for (i = 0; i <= 2; i++)
        e[i] = xearth[i] - psdp->x[i] -
               dtsp * (xearth[i+3] - psdp->x[i+3]);
      } 
    else
      for (i = 0; i <= 2; i++)
        e[i] = xearth[i] - dtsp * xearth[i+3];
        
    /* Q = planetbary(t-tau) - sunbary(t-tau) = 'planethel' */
    for (i = 0; i <= 2; i++)
      q[i] = u[i] + xearth[i] - xsun[i] -
       dtsp * (xearth[i+3] - xsun[i+3]);

    ru = vec_normalize(u,3);
         vec_normalize(q,3);
    re = vec_normalize(e,3);

    uq = dot_prod(u,q);
    ue = dot_prod(u,e);
    qe = dot_prod(q,e);

    sina = sqrt(1 - ue * ue);  /* sin(angle) between sun and planet */
    sin_sunr = SUN_RADIUS / re;   /* sine of sun radius (= sun radius) */
    
    meff_fact = (sina < sin_sunr) ?  meff(sina / sin_sunr) : 1;
      
    g1 = 2.0 * HELGRAVCONST * meff_fact / CLIGHT / CLIGHT / AUNIT / re;
    g2 = 1.0 + qe;
    for (i = 0; i <= 2; i++)
      xx3[i] = ru * (u[i] + g1/g2 * (uq * e[i] - ue * q[i]));

    for (i = 0; i <= 2; i++) {
      dx1 = xx2[i] - xx[i];
      dx2 = xx3[i] - u[i] * ru;
      dx1 -= dx2;
      xx[i+3] += dx1 / dtsp;
    }

  } /* endif speed */
  
  /* deflected position */
  vec_copy(xx,xx2,3);

  }

/* converts the sun from barycentric to geocentric,
 *          the earth from barycentric to heliocentric
 * computes
 * apparent position,
 * precession, and nutation
 * according to flags
 * iflag  flags
 * serr         error string
 */
static int app_pos_etc_sun(int iflag, char *serr, struct swe_data* swed)
{
  int j, niter, retc = OK;
  int flg1, flg2;
  double xx[6], xxsv[6], dx[3], dt, t;
  double xearth[6], xsun[6], xobs[6];
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed->oec2000;
  /* if the same conversions have already been done for the same
   * date, then return */
  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pedp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pedp->xflgs = iflag;
    pedp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }
  /************************************
   * observer: geocenter or topocenter
   ************************************/
  /* if topocentric position is wanted  */
  if (is_set(iflag,SEFLG_TOPOCTR)) {
    if (swed->topd.teval != pedp->teval
      || swed->topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag, DO_SAVE, xobs, serr) != OK) 
        return ERR;
        } 
      else {
        vec_copy(xobs,swed->topd.xobs,6);
        }
    /* barycentric position of observer */
    vec_add(xobs, pedp->x, 6);
    } 
  else {
    /* barycentric position of geocenter */
    vec_copy( xobs, pedp->x, 6);
    }
  /***************************************
   * true heliocentric position of earth *
   ***************************************/
  if (pedp->iephe == SEFLG_MOSEPH || is_set(iflag,SEFLG_BARYCTR))
    vec_copy(xx,xobs,6);
  else
    vec_sub( vec_copy(xx,xobs,6), psdp->x, 6 );

  /*******************************
   * light-time                  *
   *******************************/
  if (is_not_set(iflag,SEFLG_TRUEPOS)) {
    /* number of iterations - 1
     * the following if() does the following:
     * with jpl and swiss ephemeris:
     *   with geocentric computation of sun:
     *     light-time correction of barycentric sun position.
     *   with heliocentric or barycentric computation of earth:
     *     light-time correction of barycentric earth position.
     * with moshier ephemeris (heliocentric!!!):
     *   with geocentric computation of sun:
     *     nothing! (aberration will be done later)
     *   with heliocentric or barycentric computation of earth:
     *     light-time correction of heliocentric earth position.
     */
    if (pedp->iephe == SEFLG_JPLEPH || pedp->iephe == SEFLG_SWIEPH
      || is_set(iflag,SEFLG_HELCTR | SEFLG_BARYCTR)) {
      
      vec_copy(xearth,xobs,6);
      if (pedp->iephe == SEFLG_MOSEPH) 
          vec_clear(xsun,6);
        else
          vec_copy(xsun,psdp->x,6);
      niter = 1;  /* # of iterations */
      for (j = 0; j <= niter; j++) {
  /* distance earth-sun */
        vec_copy(dx,xearth,3);
        if (is_not_set(iflag,SEFLG_BARYCTR))
          vec_sub( dx, xsun, 3);
        
  /* new t */
        dt = sqrt(square_sum(dx)) * AUNIT / CLIGHT / 86400.0;
        t = pedp->teval - dt;
  /* new position */
        switch(pedp->iephe) {
          /* if geocentric sun, new sun at t'
           * if heliocentric or barycentric earth, new earth at t' */
          case SEFLG_SWIEPH:
            /*
              retc = sweph(t, SEI_SUN, SEI_FILE_PLANET, iflag, NULL, NO_SAVE, xearth, serr);
            */
            if (is_set(iflag,SEFLG_HELCTR | SEFLG_BARYCTR))
              retc = sweplan(t, SEI_EARTH, SEI_FILE_PLANET, iflag, NO_SAVE, xearth, NULL, xsun, NULL, serr, swed);
                  else
              retc = sweph(t, SEI_SUNBARY, SEI_FILE_PLANET, iflag, NULL, NO_SAVE, xsun, serr, swed);
            break;
          case SEFLG_MOSEPH:
            if (is_set(iflag,SEFLG_HELCTR | SEFLG_BARYCTR))
              retc = swi_moshplan(t, SEI_EARTH, NO_SAVE, xearth, xearth, serr);
            /* with moshier there is no barycentric sun */
            break;
                default:
            retc = ERR;
            break;
          }
        if (retc != OK) return(retc);
        }  
      /* apparent heliocentric earth */
      vec_copy(xx,xearth,6);
      if (is_not_set(iflag,SEFLG_BARYCTR))
        vec_sub(xx,xsun,6);
    }
  }
  
  if (is_not_set(iflag,SEFLG_SPEED))
    vec_clear(xx+3,3);

  /*******************************
   * conversion to geocenter     *
   *******************************/
  if (is_not_set(iflag,SEFLG_HELCTR | SEFLG_BARYCTR))
    smul(xx,-1,6);

  /**********************************
   * 'annual' aberration of light   *
   **********************************/
  if (is_not_set(iflag,SEFLG_TRUEPOS | SEFLG_NOABERR)) {
    /* SEFLG_NOABERR is on, if SEFLG_HELCTR or SEFLG_BARYCTR */
    swi_aberr_light(xx, xobs, iflag);
  }
  if (is_not_set(iflag,SEFLG_SPEED)) vec_clear( xx+3,3 );
    
  /* ICRS to J2000 */
  if (is_not_set(iflag,SEFLG_ICRS) && swed->jpldenum >= 403) {
    swi_bias(xx, iflag, FALSE);
  }/**/
  /* save J2000 coordinates; required for sidereal positions */
  vec_copy(xxsv,xx,6);

  /************************************************
   * precession, equator 2000 -> equator of date *
   ************************************************/
  if (is_not_set(iflag,SEFLG_J2000)) {
    swi_precess(xx, pedp->teval, J2000_TO_J);/**/
    if (is_set(iflag,SEFLG_SPEED))
      swi_precess_speed(xx, pedp->teval, J2000_TO_J);/**/
    oe = &swed->oec;
  } else
    oe = &swed->oec2000;
  return app_pos_rest(pedp, iflag, xx, xxsv, oe, serr, swed);
}


/* transforms the position of the moon:
 * heliocentric position
 * barycentric position
 * astrometric position
 * apparent position
 * precession and nutation
 *
 * note:
 * for apparent positions, we consider the earth-moon
 * system as independant.
 * for astrometric positions (SEFLG_NOABERR), we
 * consider the motions of the earth and the moon
 * related to the solar system barycenter.
 */
static int app_pos_etc_moon(int iflag, char *serr, struct swe_data *swed)
{
  int i;
  int flg1, flg2;
  double xx[6], xxsv[6], xobs[6], xxm[6], xs[6], xe[6], xobs2[6], dt;
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
  struct plan_data *pdp = &swed->pldat[SEI_MOON];
  struct epsilon *oe = &swed->oec;
  double t;
  int retc;
  /* if the same conversions have already been done for the same
   * date, then return */
  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pdp->xflgs = iflag;
    pdp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }
  /* the conversions will be done with xx[]. */
  vec_copy(xx,pdp->x,6);
  vec_copy(xxm,xx,6);

  /***********************************
   * to solar system barycentric
   ***********************************/
  vec_add(xx,pedp->x,6);

  /*******************************
   * observer
   *******************************/
  if (is_set(iflag,SEFLG_TOPOCTR)) {
    if (swed->topd.teval != pdp->teval
      || swed->topd.teval == 0) {
      if (swi_get_observer(pdp->teval, iflag, DO_SAVE, xobs, NULL) != OK) return ERR;
      } 
    else {
      vec_copy(xobs,swed->topd.xobs,6);
      }
    vec_sub(xxm,xobs,6);
    vec_add(xobs,pedp->x,6);
    } 
  else if (is_set(iflag,SEFLG_BARYCTR)) {
    vec_clear(xobs,6);
    vec_add(xxm,pedp->x,6);
  } else if (is_set(iflag,SEFLG_HELCTR)) {
    vec_copy(xobs,psdp->x,6);
    for (i = 0; i <= 5; i++)
      xxm[i] += pedp->x[i] - psdp->x[i];
  } else {
    vec_copy(xobs,pedp->x,6);
  }
  /*******************************
   * light-time                  *
   *******************************/
  if (is_set(iflag,SEFLG_TRUEPOS) == 0) {
    dt = sqrt(square_sum(xxm)) * AUNIT / CLIGHT / 86400.0;
    t = pdp->teval - dt;
    switch(pdp->iephe) {
      case SEFLG_SWIEPH:
        retc = sweplan(t, SEI_MOON, SEI_FILE_MOON, iflag, NO_SAVE, xx, xe, xs, NULL, serr, swed);
        if (retc != OK) return(retc);
        vec_add(xx,xe,6);
        break;
      case SEFLG_MOSEPH:
        /* this method results in an error of a milliarcsec in speed */
        for (i = 0; i <= 2; i++) {
          xx[i] -= dt * xx[i+3];
          xe[i] = pedp->x[i] - dt * pedp->x[i+3];
          xe[i+3] = pedp->x[i+3];
          xs[i] = 0;
          xs[i+3] = 0;
          }
        break;
      }
    if (is_set(iflag,SEFLG_TOPOCTR)) {
      if (swi_get_observer(t, iflag, NO_SAVE, xobs2, NULL) != OK) return ERR;
      vec_add(xobs2,xe,6);
    } else if (is_set(iflag,SEFLG_BARYCTR)) {
      vec_clear(xobs2,6);
    } else if (is_set(iflag,SEFLG_HELCTR)) {
      vec_copy(xobs2,xs,6);
    } else {
      vec_copy(xobs2,xe,6);
    }
  }
  /*************************
   * to correct center
   *************************/
   vec_sub(xx,xobs,6);

   /**********************************
   * 'annual' aberration of light   *
   **********************************/
  if (is_not_set(iflag,SEFLG_TRUEPOS | SEFLG_NOABERR)) {
    /* SEFLG_NOABERR is on, if SEFLG_HELCTR or SEFLG_BARYCTR */
    swi_aberr_light(xx, xobs, iflag);
    /*
     * Apparent speed is also influenced by
     * the difference of speed of the earth between t and t-dt.
     * Neglecting this would lead to an error of several 0.1"
     */
    if (is_set(iflag,SEFLG_SPEED))
      for (i = 3; i <= 5; i++)
        xx[i] += xobs[i] - xobs2[i];
  }

  /* if !speedflag, speed = 0 */
  if (is_not_set(iflag,SEFLG_SPEED))
    vec_clear(xx+3,3);

  /* ICRS to J2000 */
  if (is_not_set(iflag,SEFLG_ICRS) && swed->jpldenum >= 403) {
    swi_bias(xx, iflag, FALSE);
  }/**/
  /* save J2000 coordinates; required for sidereal positions */
  for (i = 0; i <= 5; i++)
    xxsv[i] = xx[i];
  /************************************************
   * precession, equator 2000 -> equator of date *
   ************************************************/
  if (is_not_set(iflag,SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, J2000_TO_J);
    if (is_set(iflag,SEFLG_SPEED))
      swi_precess_speed(xx, pdp->teval, J2000_TO_J);
    oe = &swed->oec;
  } else
    oe = &swed->oec2000;
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr, swed);
}

/* transforms the position of the barycentric sun:
 * precession and nutation
 * according to flags
 * iflag  flags
 * serr         error string
 */
static int app_pos_etc_sbar(int iflag, char *serr, struct swe_data *swed)
{
  int i;
  double xx[6], xxsv[6], dt;
  struct plan_data *psdp = &swed->pldat[SEI_EARTH];
  struct plan_data *psbdp = &swed->pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed->oec;
  /* the conversions will be done with xx[]. */
  vec_copy(xx,psbdp->x,6);
  /**************
   * light-time *
   **************/
  if (is_not_set(iflag,SEFLG_TRUEPOS)) {
    dt = sqrt(square_sum(xx)) * AUNIT / CLIGHT / 86400.0;
    for (i = 0; i <= 2; i++)
      xx[i] -= dt * xx[i+3];  /* apparent position */
  }
  if (is_not_set(iflag,SEFLG_SPEED)) vec_clear(xx+3,3);

  /* ICRS to J2000 */
  if (is_not_set(iflag,SEFLG_ICRS) && swed->jpldenum >= 403) {
    swi_bias(xx, iflag, FALSE);
  }
  /* save J2000 coordinates; required for sidereal positions */
  vec_copy(xxsv,xx,6);

  /************************************************
   * precession, equator 2000 -> equator of date *
   ************************************************/
  if (is_not_set(iflag,SEFLG_J2000)) {
    swi_precess(xx, psbdp->teval, J2000_TO_J);
    if (is_set(iflag,SEFLG_SPEED))
      swi_precess_speed(xx, psbdp->teval, J2000_TO_J);
    oe = &swed->oec;
  } else
    oe = &swed->oec2000;
  return app_pos_rest(psdp, iflag, xx, xxsv, oe, serr, swed);
}

/* transforms position of mean lunar node or apogee:
 * input is polar coordinates in mean ecliptic of date.
 * output is, according to iflag:
 * position accounted for light-time
 * position referred to J2000 (i.e. precession subtracted)
 * position with nutation
 * equatorial coordinates
 * cartesian coordinates
 * heliocentric position is not allowed ??????????????
 *         DAS WAERE ZIEMLICH AUFWENDIG. SONNE UND ERDE MUESSTEN
 *         SCHON VORHANDEN SEIN!
 * ipl    bodynumber (SE_MEAN_NODE or SE_MEAN_APOG)
 * iflag  flags
 * serr         error string
 */
static int app_pos_etc_mean(int ipl, int iflag, char *serr, struct swe_data *swed) {
  int flg1, flg2;
  double xx[6], xxsv[6];
  struct plan_data *pdp = &swed->nddat[ipl];
  struct epsilon *oe;
  /* if the same conversions have already been done for the same
   * date, then return */
  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = pdp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  if (flg1 == flg2) {
    pdp->xflgs = iflag;
    pdp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }
  vec_copy(xx,pdp->x,6);

  /* cartesian equatorial coordinates */
  swi_polcart_sp(xx, xx);
  swi_coortrf2(xx, xx, -swed->oec.seps, swed->oec.ceps);
  swi_coortrf2(xx+3, xx+3, -swed->oec.seps, swed->oec.ceps);
  if (is_not_set(iflag,SEFLG_SPEED)) vec_clear(xx+3,3);

  /* J2000 coordinates; required for sidereal positions */
  if ((is_set(iflag,SEFLG_SIDEREAL)
    && (swed->sidd.sid_mode & SE_SIDBIT_ECL_T0))
      || (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE)) {
    vec_copy(xxsv,xx,6);
    /* xxsv is not J2000 yet! */
    if (pdp->teval != J2000) {
      swi_precess(xxsv, pdp->teval, J_TO_J2000);
      if (is_set(iflag,SEFLG_SPEED))
        swi_precess_speed(xxsv, pdp->teval, J_TO_J2000);
    }
  }
  /*****************************************************
   * if no precession, equator of date -> equator 2000 *
   *****************************************************/
  if (is_set(iflag,SEFLG_J2000)) {
    swi_precess(xx, pdp->teval, J_TO_J2000);
    if (is_set(iflag,SEFLG_SPEED))
      swi_precess_speed(xx, pdp->teval, J_TO_J2000);
    oe = &swed->oec2000;
  } else
    oe = &swed->oec;
  return app_pos_rest(pdp, iflag, xx, xxsv, oe, serr, swed);
}

/* fetch chebyshew coefficients from sweph file for
 * tjd     time
 * ipli    planet number
 * ifno    file number
 * serr    error string
 */
static void get_new_segment(double tjd, int ipli, int ifno, char *serr, struct swe_data *swed) {

  int i, j, k, m, n, o, icoord;
  int iseg;
  int fpos;
  int nsizes, nsize[6];
  int nco;
  int idbl;
  unsigned char c[4];
  struct plan_data *pdp = &swed->pldat[ipli];
  struct file_data *fdp = &swed->fidat[ifno];
  FILE *fp = fdp->fptr;
  uint32 longs[MAXORD+1];
  
  struct file_context fc = set_file_context(ifno,serr,swed);
      
  /* compute segment number */
  iseg = (int) ((tjd - pdp->tfstart) / pdp->dseg);
  pdp->tseg0 = pdp->tfstart + iseg * pdp->dseg;
  pdp->tseg1 = pdp->tseg0 + pdp->dseg;

  /* get file position of coefficients from file */
  fpos = pdp->lndx0 + iseg * 3;
  fseek(fp, fpos, SEEK_SET);
  do_fread((void *) &fpos, 3, 1, 4, &fc);
  fseek(fp, fpos, SEEK_SET);
  
  /* clear space of chebyshew coefficients */
  if (pdp->segp == NULL)
    pdp->segp = (double *) malloc((size_t) pdp->ncoe * 3 * 8);
    
  memset((void *) pdp->segp, 0, (size_t) pdp->ncoe * 3 * 8);
  
  /* read coefficients for 3 coordinates */
  for (icoord = 0; icoord < 3; icoord++) {
    idbl = icoord * pdp->ncoe;
    /* first read header */
    /* first bit indicates number of sizes of packed coefficients */
    do_fread((void *) &c[0], 1, 2, 1, &fc);
    if (c[0] & 128) {
      nsizes = 6;
      do_fread((void *) (c+2), 1, 2, 1, &fc);
      nsize[0] = (int) c[1] / 16;
      nsize[1] = (int) c[1] % 16;
      nsize[2] = (int) c[2] / 16;
      nsize[3] = (int) c[2] % 16;
      nsize[4] = (int) c[3] / 16;
      nsize[5] = (int) c[3] % 16;
      nco = nsize[0] + nsize[1] + nsize[2] + nsize[3] + nsize[4] + nsize[5];
    } else {
      nsizes = 4;
      nsize[0] = (int) c[0] / 16;
      nsize[1] = (int) c[0] % 16;
      nsize[2] = (int) c[1] / 16;
      nsize[3] = (int) c[1] % 16;
      nco = nsize[0] + nsize[1] + nsize[2] + nsize[3];
    }
    /* there may not be more coefficients than interpolation
     * order + 1 */
    if (nco > pdp->ncoe) throw_coeff_error(nco,pdp,&fc);

    /* now unpack */
    for (i = 0; i < nsizes; i++) {
      if (nsize[i] == 0) continue;
      if (i < 4) {
        j = (4 - i);
        k = nsize[i];
        do_fread((void *) &longs[0], j, k, 4, &fc);
        for (m = 0; m < k; m++, idbl++) {
          if (longs[m] & 1)   /* will be negative */
            pdp->segp[idbl] = -(((longs[m]+1) / 2) / 1e+9 * pdp->rmax / 2);
          else
            pdp->segp[idbl] = (longs[m] / 2) / 1e+9 * pdp->rmax / 2;
          }
        } 
      else if (i == 4) {    /* half byte packing */
        j = 1;
        k = (nsize[i] + 1) / 2;
        do_fread((void *) longs, j, k, 4, &fc);
        for (m = 0, j = 0;
             m < k && j < nsize[i];
             m++) {
          for (n = 0, o = 16;
               n < 2 && j < nsize[i];
               n++, j++, idbl++, longs[m] %= o, o /= 16) {
            if (longs[m] & o)
              pdp->segp[idbl] = -(((longs[m]+o) / o / 2) * pdp->rmax / 2 / 1e+9);
            else
              pdp->segp[idbl] = (longs[m] / o / 2) * pdp->rmax / 2 / 1e+9;
            } 
          }
        } 
      else if (i == 5) {    /* quarter byte packing */
        j = 1;
        k = (nsize[i] + 3) / 4;
        do_fread((void *) longs, j, k, 4, &fc);
        for (m = 0, j = 0;
             m < k && j < nsize[i];
             m++) {
          for (n = 0, o = 64;
               n < 4 && j < nsize[i];
               n++, j++, idbl++, longs[m] %= o, o /= 4) {
            if (longs[m] & o)
              pdp->segp[idbl] = -(((longs[m]+o) / o / 2) * pdp->rmax / 2 / 1e+9);
            else
              pdp->segp[idbl] = (longs[m] / o / 2) * pdp->rmax / 2 / 1e+9;
            }
          }
        }
      }
    }
  }

static struct file_context set_file_context( int ifno, char* serr, struct swe_data *swed) {  
  struct file_context fc;
  struct file_data *fdp = &swed->fidat[ifno];
  fc.ifno = ifno;
  fc.serr = serr;
  fc.swed = swed;
  fc.fp   = fdp->fptr;
  fc.fdp  = fdp;    
// reord and endian are undefined before read_const()  
  fc.freord = (int) fdp->iflg & SEI_FILE_REORD;
  fc.fendian = (int) fdp->iflg & SEI_FILE_LITENDIAN;
  fc.fpos = SEI_CURR_FPOS;    
  return fc;
  }
  

static void throw_coeff_error(int nco,struct plan_data *pdp,struct file_context *fc) {
  char msg[AS_MAXCH];
  sprintf(msg,"error in ephemeris file: %d coefficients instead of %d. ", nco, pdp->ncoe);
  if (strlen(msg) + strlen(fc->fdp->fnam) < AS_MAXCH - 1)
    sprintf(msg, "error in ephemeris file %s: %d coefficients instead of %d. ", fc->fdp->fnam, nco, pdp->ncoe);
  free(pdp->segp);
  pdp->segp = NULL;
  throw_file_error(fc->fp,msg,fc->serr,fc->swed);
  }    

/* SWISSEPH
 * reads constants on ephemeris file
 * ifno         file #
 * serr         error string
 */
static void read_const(int ifno, char *serr, struct swe_data *swed)
{
  char s[AS_MAXCH*2];
  char sastnam[41];
  int lastnam = 19;
  int ipli, kpl;
  int lng;
  uint32 ulng;
  int flen, fpos;
  short nplan;
  double doubles[20];
  struct plan_data *pdp;
  struct file_data *fdp = &swed->fidat[ifno];
  int nbytes_ipl = 2;
  FILE *fp = fdp->fptr;
  
  struct file_context fc = set_file_context(ifno,serr,swed);
    
// version number (unused so far) */
  read_line(s,fdp,serr,swed);
  fdp->fversion = int_from_string(s,fdp,serr,swed);
  
  /*************************************
   * correct file name?                *
   *************************************/
  read_line(s,fdp,serr,swed);  
  check_filename(s,fdp,serr,swed);
  
   /*************************************
   * copyright                         *
   *************************************/
  read_line(s,fdp,serr,swed);  

  /****************************************
   * orbital elements, if single asteroid *
   ****************************************/
  if (ifno == SEI_FILE_ANY_AST) {
    read_line(s,fdp,serr,swed);
    orbital_elements_single_asteroid(s,sastnam,lastnam,fdp,serr,swed);
    }
  
  /*************************************
   * one int for test of byte order   *
   *************************************/
  determine_byte_order( &fc ); 

  /*************************************
   * length of file correct?           *
   *************************************/
  do_fread((void *) &lng, 4, 1, 4, &fc); // fp, SEI_CURR_FPOS, freord, fendian, ifno, serr, swed);
  fpos = ftell(fp);
  if (fseek(fp, 0L, SEEK_END) != 0) throw_file_damaged( fdp, serr, swed);
  flen = ftell(fp);
  if (lng != flen) throw_file_damaged( fdp, serr, swed);
  fseek(fp,fpos,SEEK_SET); // Restore current position
  
  /**********************************************************
   * DE number of JPL ephemeris which this file is based on *
   **********************************************************/
  do_fread((void *) &fdp->sweph_denum, 4, 1, 4, &fc); 
  swed->jpldenum = fdp->sweph_denum;
  /*************************************
   * start and end epoch of file       *
   *************************************/
  do_fread((void *) &fdp->tfstart, 8, 1, 8, &fc); 
  do_fread((void *) &fdp->tfend, 8, 1, 8, &fc);
    
  /*************************************
   * how many planets are in file?     *
   *************************************/
  do_fread((void *) &nplan, 2, 1, 2, &fc);
  if (nplan > 256) {
    nbytes_ipl = 4;
    nplan %= 256;
  }
  if (nplan < 1 || nplan > 20)
    throw_file_damaged( fdp, serr, swed);
  fdp->npl = nplan;
  /* which ones?                       */
  do_fread((void *) fdp->ipl, nbytes_ipl, (int) nplan, sizeof(int), &fc);

  /*************************************
   * asteroid name                     *
   *************************************/
  if (ifno == SEI_FILE_ANY_AST) {
    read_astnam(fdp->astnam, sastnam, lastnam, &fc);
    }
    
  /*************************************
   * check CRC                         *
   *************************************/
  fpos = ftell(fp);
  /* read CRC from file */
  do_fread((void *) &ulng, 4, 1, 4, &fc);

  /* read check area from file */
  fseek(fp, 0L, SEEK_SET);
  /* must check that defined length of s is less than fpos */
  if (fpos - 1 > 2 * AS_MAXCH)
    throw_file_damaged( fdp, serr, swed);
  if (fread((void *) s, (size_t) fpos, 1, fp) != 1)
    throw_file_damaged( fdp, serr, swed);
  if (swi_crc32((unsigned char *) s, (int) fpos) != ulng)
    throw_file_damaged( fdp, serr, swed);
    /*printf("crc %d %d\n", ulng2, ulng);*/
  fseek(fp, fpos+4, SEEK_SET);
  /*************************************
   * read general constants            *
   *************************************/
  /* clight, aunit, helgravconst, ratme, sunradius
   * these constants are currently not in use */
  do_fread((void *) doubles, 8, 5, 8, &fc);

  swed->gcdat.clight       = doubles[0];
  swed->gcdat.aunit        = doubles[1];
  swed->gcdat.helgravconst = doubles[2];
  swed->gcdat.ratme        = doubles[3];
  swed->gcdat.sunradius    = doubles[4];

  /*************************************
   * read constants of planets         *
   *************************************/
  for (kpl = 0; kpl < fdp->npl; kpl++) {
    /* get SEI_ planet number */
    ipli = fdp->ipl[kpl];
    pdp = &swed->pldat[ipli >= SE_AST_OFFSET ? SEI_ANYBODY:ipli];
    pdp->ibdy = ipli;
    /* file position of planet's index */
    do_fread((void *) &pdp->lndx0, 4, 1, 4, &fc);

    /* flags: helio/geocentric, rotation, reference ellipse */
    do_fread((void *) &pdp->iflg, 1, 1, sizeof(int), &fc);

    /* number of chebyshew coefficients / segment  */
    /* = interpolation order +1                    */
    do_fread((void *) &pdp->ncoe, 1, 1, sizeof(int), &fc);

    /* rmax = normalisation factor */
    do_fread((void *) &lng, 4, 1, 4, &fc);
    pdp->rmax = lng / 1000.0;

    /* start and end epoch of planetary ephemeris,   */
    /* segment length, and orbital elements          */
    do_fread((void *) doubles, 8, 10, 8, &fc);

    pdp->tfstart  = doubles[0];
    pdp->tfend    = doubles[1];
    pdp->dseg     = doubles[2];
    pdp->nndx     = (int) ((doubles[1] - doubles[0] + 0.1) /doubles[2]);
    pdp->telem    = doubles[3];
    pdp->prot     = doubles[4];
    pdp->dprot    = doubles[5];
    pdp->qrot     = doubles[6];
    pdp->dqrot    = doubles[7];
    pdp->peri     = doubles[8];
    pdp->dperi    = doubles[9];
    /* alloc space for chebyshew coefficients */
    /* if reference ellipse is used, read its coefficients */
    if (pdp->iflg & SEI_FLG_ELLIPSE) {
      if (pdp->refep != NULL) { /* if switch to other eph. file */
        free((void *) pdp->refep);
        if (pdp->segp != NULL) {
          free((void *) pdp->segp);     /* array of coefficients of */
          pdp->segp = NULL;     /* ephemeris segment        */
        }
      }
      pdp->refep = (double *) malloc((size_t) pdp->ncoe * 2 * 8);
      do_fread((void *) pdp->refep, 8, 2*pdp->ncoe, 8, &fc);
      
    }  
  }
}

static void read_astnam( char* astnam, char* sastnam, int lastnam, struct file_context *fc) {
    char sastno[12];
    char s[AS_MAXCH*2];
    /* name of asteroid is taken from orbital elements record
     * read above */
    int j = 4;  /* old astorb.dat had only 4 characters for MPC# */
    while (sastnam[j] != ' ' && j < 10)  /* new astorb.dat has 5 */
      j++;
    strncpy(sastno, sastnam, j);
    sastno[j] = '\0';
    int i = (int) atol(sastno);
    if (i == fc->fdp->ipl[0] - SE_AST_OFFSET) {
      /* element record is from bowell database */
      strncpy(astnam, sastnam+j+1, lastnam);
      /* overread old ast. name field */
      if (fread((void *) s, 30, 1, fc->fp) != 1)
        throw_file_damaged( fc->fdp, fc->serr, fc->swed);
    } else {
      /* older elements record structure: the name
       * is taken from old name field */
      if (fread((void *) astnam, 30, 1, fc->fp) != 1)
        throw_file_damaged( fc->fdp, fc->serr, fc->swed);
    }
    /* in worst case strlen of not null terminated area! */
    i = strlen(astnam) - 1;
    if (i < 0) i = 0;
    char *sp = astnam + i;
    while(*sp == ' ') {
      sp--;
    }
    sp[1] = '\0';
  }

static void read_line(char* s, struct file_data* fdp, char* serr, struct swe_data* swed) {
  char *sp;
  sp = fgets(s, AS_MAXCH*2, fdp->fptr);
  if (sp == NULL || strstr(sp, "\r\n") == NULL)
    throw_file_damaged(fdp,serr,swed);
  for(sp = strchr(s, '\r');sp>s && isspace(*sp); --sp);  // trim trailing whitespace
  *(sp+1) = '\0';    
  }

 static void orbital_elements_single_asteroid(char* s,char *sastnam,int lastnam,struct file_data* fdp, char* serr, struct swe_data* swed) {
   char *sp;
   char s2[AS_MAXCH];
   /* MPC number and name; will be analyzed below:
     * search "asteroid name" */
    while(*sp == ' ') sp++;
    while(isdigit(*sp)) sp++;
    sp++;
    int i = sp - s;
    strncpy(sastnam, sp, lastnam+i);
    *(sastnam+lastnam+i) = '\0';
    /* save elements, they are required for swe_plan_pheno() */
    strcpy(swed->astelem, s);
    /* required for magnitude */
    swed->ast_H = atof(s + 35 + i);
    swed->ast_G = atof(s + 42 + i);
    if (swed->ast_G == 0) swed->ast_G = 0.15;
    /* diameter in kilometers, not always given: */
    strncpy(s2, s+51+i, 7);
    *(s2 + 7) = '\0';
    swed->ast_diam = atof(s2);
    if (swed->ast_diam == 0) {
      /* estimate the diameter from magnitude; assume albedo = 0.15 */
      swed->ast_diam = 1329/sqrt(0.15) * pow(10, -0.2 * swed->ast_H);
    }
  }  

static void determine_byte_order( struct file_context* fc) {  
  
  int testendian;
  int lng;
  char *c, c2, *sp;

  if (fread((void *) &testendian, 4, 1, fc->fp) != 1)
    throw_file_damaged( fc->fdp, fc->serr, fc->swed);
  /* is byte order correct?            */
  if (testendian == SEI_FILE_TEST_ENDIAN)
    fc->freord = SEI_FILE_NOREORD;
  else {
    fc->freord = SEI_FILE_REORD;
    sp = (char *) &lng;
    c = (char *) &testendian;
    for (int i = 0; i < 4; i++)
      *(sp+i) = *(c+3-i);
    if (lng != SEI_FILE_TEST_ENDIAN)
      throw_file_damaged( fc->fdp, fc->serr, fc->swed);
      /* printf("%d  %x\n", lng, lng);*/
  }
  /* is file bigendian or littlendian?
   * test first byte of test integer, which is highest if bigendian */
  c = (char *) &testendian;
  c2 = SEI_FILE_TEST_ENDIAN / 16777216L;
  if (*c == c2)
    fc->fendian = SEI_FILE_BIGENDIAN;
  else
    fc->fendian = SEI_FILE_LITENDIAN;
  fc->fdp->iflg = (int) fc->freord | fc->fendian;
  }
  
  static int int_from_string(char* s, struct file_data* fdp, char* serr, struct swe_data* swed) {
  char *sp = s;
  while (isdigit((int) *sp) == 0 && *sp != '\0') sp++;
  if (*sp == '\0') throw_file_damaged( fdp, serr, swed);
  return atoi(sp);
  }  
  
static void check_filename(char* name_exp, struct file_data* fdp, char* serr, struct swe_data* swed) {
  char name_act[AS_MAXCH], msg[AS_MAXCH];
  filename_without_path( name_act, fdp->fnam);
  if (!equals_ignore_case(name_exp,name_act)) {
    sprintf(msg, "Ephemeris file name '%s' wrong; rename '%s' ", name_act, name_exp);
    throw_file_error(fdp->fptr,msg,serr,swed);
    }  
  }

static bool equals_ignore_case(char* s1, char*s2) {
  int n = strlen(s1);
  if (n != strlen(s2)) 
    return false;
  else 
    for (int i=n-1;i>=0;--i) {
      if (tolower(s1[i])!=tolower(s2[i])) return false;
      }  
  return true;
  }
  
static void filename_without_path( char* name, char* fullname) {
  char *sp = strrchr(fullname, (int) *DIR_GLUE);
  if (sp == NULL)
    sp = fullname;
  else
    sp++;
  strcpy(name, sp);  
  }  
  
static void throw_file_damaged(struct file_data *fdp, char* serr, struct swe_data *swed) {
  char *serr_file_damage = "Ephemeris file %s is damaged.";
  char msg[AS_MAXCH];
  int errmsglen = strlen(serr_file_damage) + strlen(fdp->fnam);
  sprintf(msg, serr_file_damage, (errmsglen < AS_MAXCH) ? fdp->fnam : "" );
  throw_file_error(fdp->fptr,msg,serr,swed);
  }

static void throw_file_error(FILE *fp, char* msg, char* serr, struct swe_data *swed) {
  fclose(fp);
  throw(ERR,msg,serr,swed);
  }

/* SWISSEPH
 * reads from a file and, if necessary, reorders bytes
 * targ   target pointer
 * size    size of item to be read
 * count  number of items
 * corrsize  in what size should it be returned
 *    (e.g. 3 byte int -> 4 byte int)
 * fp    file pointer
 * fpos    file position: if (fpos >= 0) then fseek
 * freord  reorder bytes or no
 * fendian  little/bigendian
 * ifno    file number
 * serr    error string
 */
static void do_fread(void *trg, int size, int count, int corrsize, struct file_context *fc) {
  unsigned char space[1000];
  unsigned char *targ = (unsigned char *) trg;

  int totsize = size * count;
  if (fc->fpos >= 0)
    fseek(fc->fp, fc->fpos, SEEK_SET);
  /* if no byte reorder has to be done, and read size == return size */
  if (!fc->freord && size == corrsize) {
    if (fread((void *) targ, (size_t) totsize, 1, fc->fp) == 0) {
      throw_file_damaged(fc->fdp,fc->serr,fc->swed);
      } 
    } 
  else {
    if (fread((void *) &space[0], (size_t) totsize, 1, fc->fp) == 0)
      throw_file_damaged(fc->fdp,fc->serr,fc->swed);
    if (size != corrsize) {
      memset((void *) targ, 0, (size_t) count * corrsize);
      }
    for(int i = 0; i < count; i++) {
      for (int j = size-1; j >= 0; j--) {
        int k = (fc->freord) ? size-j-1 : j;
        if (size != corrsize)
          if ((fc->fendian == SEI_FILE_BIGENDIAN && !fc->freord) ||
              (fc->fendian == SEI_FILE_LITENDIAN &&  fc->freord))
            k += corrsize - size;
        targ[i*corrsize+k] = space[i*size+j];
      }
    }
  }
}

/* SWISSEPH
 * adds reference orbit to chebyshew series (if SEI_FLG_ELLIPSE),
 * rotates series to mean equinox of J2000
 *
 * ipli    planet number
 */
static void rot_back(int ipli, struct swe_data *swed)
{
  int i;
  double t, tdiff;
  double qav, pav, dn;
  double omtild, com, som, cosih2;
  double x[MAXORD+1][3];
  double uix[3], uiy[3], uiz[3];
  double xrot, yrot, zrot;
  double *chcfx, *chcfy, *chcfz;
  double *refepx, *refepy;
  double seps2000 = swed->oec2000.seps;
  double ceps2000 = swed->oec2000.ceps;
  struct plan_data *pdp = &swed->pldat[ipli];
  int nco = pdp->ncoe;
  t = pdp->tseg0 + pdp->dseg / 2;
  chcfx = pdp->segp;
  chcfy = chcfx + nco;
  chcfz = chcfx + 2 * nco;
  refepx = pdp->refep;
  refepy = refepx + nco;
  tdiff= (t - pdp->telem) / 365250.0;
  if (ipli == SEI_MOON) {
    dn = pdp->prot + tdiff * pdp->dprot;
    i = (int) (dn / TWOPI);
    dn -= i * TWOPI;
    qav = (pdp->qrot + tdiff * pdp->dqrot) * cos(dn);
    pav = (pdp->qrot + tdiff * pdp->dqrot) * sin(dn);
  } else {
    qav = pdp->qrot + tdiff * pdp->dqrot;
    pav = pdp->prot + tdiff * pdp->dprot;
  }
  /*calculate cosine and sine of average perihelion longitude. */
  for (i = 0; i < nco; i++) {
    x[i][0] = chcfx[i];
    x[i][1] = chcfy[i];
    x[i][2] = chcfz[i];
  }
  if (pdp->iflg & SEI_FLG_ELLIPSE) {
    omtild = pdp->peri + tdiff * pdp->dperi;
    i = (int) (omtild / TWOPI);
    omtild -= i * TWOPI;
    com = cos(omtild);
    som = sin(omtild);
    /*add reference orbit.  */
    for (i = 0; i < nco; i++) {
      x[i][0] = chcfx[i] + com * refepx[i] - som * refepy[i];
      x[i][1] = chcfy[i] + com * refepy[i] + som * refepx[i];
    }
  }
  /* construct right handed orthonormal system with first axis along
     origin of longitudes and third axis along angular momentum
     this uses the standard formulas for equinoctal variables
     (see papers by broucke and by cefola).      */
  cosih2 = 1.0 / (1.0 + qav * qav + pav * pav);
  /*     calculate orbit pole. */
  uiz[0] = 2.0 * pav * cosih2;
  uiz[1] = -2.0 * qav * cosih2;
  uiz[2] = (1.0 - qav * qav - pav * pav) * cosih2;
  /*     calculate origin of longitudes vector. */
  uix[0] = (1.0 + qav * qav - pav * pav) * cosih2;
  uix[1] = 2.0 * qav * pav * cosih2;
  uix[2] = -2.0 * pav * cosih2;
  /*     calculate vector in orbital plane orthogonal to origin of
        longitudes.                                               */
  uiy[0] =2.0 * qav * pav * cosih2;
  uiy[1] =(1.0 - qav * qav + pav * pav) * cosih2;
  uiy[2] =2.0 * qav * cosih2;
  /*     rotate to actual orientation in space.         */
  for (i = 0; i < nco; i++) {
    xrot = x[i][0] * uix[0] + x[i][1] * uiy[0] + x[i][2] * uiz[0];
    yrot = x[i][0] * uix[1] + x[i][1] * uiy[1] + x[i][2] * uiz[1];
    zrot = x[i][0] * uix[2] + x[i][1] * uiy[2] + x[i][2] * uiz[2];
    if (fabs(xrot) + fabs(yrot) + fabs(zrot) >= 1e-14)
      pdp->neval = i;
    x[i][0] = xrot;
    x[i][1] = yrot;
    x[i][2] = zrot;
    if (ipli == SEI_MOON) {
      /* rotate to j2000 equator */
      x[i][1] = ceps2000 * yrot - seps2000 * zrot;
      x[i][2] = seps2000 * yrot + ceps2000 * zrot;
    }
  }
  for (i = 0; i < nco; i++) {
    chcfx[i] = x[i][0];
    chcfy[i] = x[i][1];
    chcfz[i] = x[i][2];
  }
}

/* Adjust position from Earth-Moon barycenter to Earth
 *
 * xemb = hel./bar. position or velocity vectors of emb (input)
 *                                                  earth (output)
 * xmoon= geocentric position or velocity vector of moon
 */
static void embofs(double *xemb, double *xmoon)
{
  int i;
  for (i = 0; i <= 2; i++)
    xemb[i] -= xmoon[i] / (EARTH_MOON_MRAT + 1.0);
}

/* calculates the nutation matrix
 * nu    pointer to nutation data structure
 * oe    pointer to epsilon data structure
 */
static void nut_matrix(struct nut *nu, struct epsilon *oe)
{
  double psi, eps;
  double sinpsi, cospsi, sineps, coseps, sineps0, coseps0;
  psi = nu->nutlo[0];
  eps = oe->eps + nu->nutlo[1];
  sinpsi = sin(psi);
  cospsi = cos(psi);
  sineps0 = oe->seps;
  coseps0 = oe->ceps;
  sineps = sin(eps);
  coseps = cos(eps);
  nu->matrix[0][0] = cospsi;
  nu->matrix[0][1] = sinpsi * coseps;
  nu->matrix[0][2] = sinpsi * sineps;
  nu->matrix[1][0] = -sinpsi * coseps0;
  nu->matrix[1][1] = cospsi * coseps * coseps0 + sineps * sineps0;
  nu->matrix[1][2] = cospsi * sineps * coseps0 - coseps * sineps0;
  nu->matrix[2][0] = -sinpsi * sineps0;
  nu->matrix[2][1] = cospsi * coseps * sineps0 - sineps * coseps0;
  nu->matrix[2][2] = cospsi * sineps * sineps0 + coseps * coseps0;
}

/* lunar osculating elements, i.e.
 * osculating node ('true' node) and
 * osculating apogee ('black moon', 'lilith').
 * tjd    julian day
 * ipl    body number, i.e. SEI_TRUE_NODE or SEI_OSCU_APOG
 * iflag  flags (which ephemeris, nutation, etc.)
 * serr    error string
 *
 * definitions and remarks:
 * the osculating node and the osculating apogee are defined
 * as the orbital elements of the momentary lunar orbit.
 * their advantage is that when the moon crosses the ecliptic,
 * it is really at the osculating node, and when it passes
 * its greatest distance from earth it is really at the
 * osculating apogee. with the mean elements this is not
 * the case. (some define the apogee as the second focus of
 * the lunar ellipse. but, as seen from the geocenter, both
 * points are in the same direction.)
 * problems:
 * the osculating apogee is given in the 'New International
 * Ephemerides' (Editions St. Michel) as the 'True Lilith'.
 * however, this name is misleading. this point is based on
 * the idea that the lunar orbit can be approximated by an
 * ellipse.
 * arguments against this:
 * 1. this procedure considers celestial motions as two body
 *    problems. this is quite good for planets, but not for
 *    the moon. the strong gravitational attraction of the sun
 *    destroys the idea of an ellipse.
 * 2. the NIE 'True Lilith' has strong oscillations around the
 *    mean one with an amplitude of about 30 degrees. however,
 *    when the moon is in apogee, its distance from the mean
 *    apogee never exceeds 5 degrees.
 * besides, the computation of NIE is INACCURATE. the mistake
 * reaches 20 arc minutes.
 * According to Santoni, the point was calculated using 'les 58
 * premiers termes correctifs au Perigee moyen' published by
 * Chapront and Chapront-Touze. And he adds: "Nous constatons
 * que meme en utilisant ces 58 termes CORRECTIFS, l'erreur peut
 * atteindre 0,5d!" (p. 13) We avoid this error, computing the
 * orbital elements directly from the position and the speed vector.
 *
 * how about the node? it is less problematic, because we
 * we needn't derive it from an orbital ellipse. we can say:
 * the axis of the osculating nodes is the intersection line of
 * the actual orbital plane of the moon and the plane of the
 * ecliptic. or: the osculating nodes are the intersections of
 * the two great circles representing the momentary apparent
 * orbit of the moon and the ecliptic. in this way they make
 * some sense. then, the nodes are really an axis, and they
 * have no geocentric distance. however, in this routine
 * we give a distance derived from the osculating ellipse.
 * the node could also be defined as the intersection axis
 * of the lunar orbital plane and the solar orbital plane,
 * which is not precisely identical to the ecliptic. this
 * would make a difference of several arcseconds.
 *
 * is it possible to keep the idea of a continuously moving
 * apogee that is exact at the moment when the moon passes
 * its greatest distance from earth?
 * to achieve this, we would probably have to interpolate between
 * the actual apogees.
 * the nodes could also be computed by interpolation. the resulting
 * nodes would deviate from the so-called 'true node' by less than
 * 30 arc minutes.
 *
 * sidereal and j2000 true node are first computed for the ecliptic
 * of epoch and then precessed to ecliptic of t0(ayanamsa) or J2000.
 * there is another procedure that computes the node for the ecliptic
 * of t0(ayanamsa) or J2000. it is excluded by
 * #ifdef SID_TNODE_FROM_ECL_T0
 */
static int lunar_osc_elem(double tjd, int ipl, int iflag, char *serr, struct swe_data *swed)
{
  int i, j, istart;
  int epheflag = SEFLG_DEFAULTEPH;
  int retc = ERR;
  int flg1, flg2;
  struct plan_data *ndp, *ndnp, *ndap;
  struct epsilon *oe;
  double speed_intv = NODE_CALC_INTV;  /* to silence gcc warning */
  double a, b;
  double xpos[3][6], xx[3][6], xxa[3][6], xnorm[6], r[6];
  double rxy, rxyz, t, dt, fac, sgn;
  double sinnode, cosnode, sinincl, cosincl, sinu, cosu, sinE, cosE;
  double uu, ny, sema, ecce, Gmsm, c2, v2, pp;
  int speedf1, speedf2;
#ifdef SID_TNODE_FROM_ECL_T0
  struct sid_data *sip = &swed->sidd;
  struct epsilon oectmp;
  if (is_set(iflag,SEFLG_SIDEREAL)) {
    calc_epsilon(sip->t0, &oectmp);
    oe = &oectmp;
  } else if (is_set(iflag,SEFLG_J2000))
    oe = &swed->oec2000;
  else
#endif
    oe = &swed->oec;
  ndp = &swed->nddat[ipl];
  /* if elements have already been computed for this date, return
   * if speed flag has been turned on, recompute */
  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = ndp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  speedf1 = ndp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == ndp->teval
  && tjd != 0
  && flg1 == flg2
  && (!speedf2 || speedf1)) {
    ndp->xflgs = iflag;
    ndp->iephe = iflag & SEFLG_EPHMASK;
    return OK;
  }
  /* the geocentric position vector and the speed vector of the
   * moon make up the lunar orbital plane. the position vector
   * of the node is along the intersection line of the orbital
   * plane and the plane of the ecliptic.
   * to calculate the osculating node, we need one lunar position
   * with speed.
   * to calculate the speed of the osculating node, we need
   * three lunar positions and the speed of each of them.
   * this is relatively cheap, if the jpl-moon or the swisseph
   * moon is used. with the moshier moon this is much more
   * expensive, because then we need 9 lunar positions for
   * three speeds. but one position and speed can normally
   * be taken from swed->pldat[moon], which corresponds to
   * three moshier moon calculations.
   * the same is also true for the osculating apogee: we need
   * three lunar positions and speeds.
   */
  /*********************************************
   * now three lunar positions with speeds     *
   *********************************************/
  if (is_set(iflag,SEFLG_MOSEPH))
    epheflag = SEFLG_MOSEPH;
  else if (is_set(iflag,SEFLG_SWIEPH))
    epheflag = SEFLG_SWIEPH;
  else if (is_set(iflag,SEFLG_JPLEPH))
    epheflag = SEFLG_JPLEPH;
  /* there may be a moon of wrong ephemeris in save area
   * force new computation: */
  swed->pldat[SEI_MOON].teval = 0;
  if (is_set(iflag,SEFLG_SPEED))
    istart = 0;
  else
    istart = 2;
  if (serr != NULL)
    *serr = '\0';
  three_positions:
  switch(epheflag) {
    case SEFLG_SWIEPH:
      speed_intv = NODE_CALC_INTV;
      for (i = istart; i <= 2; i++) {
  if (i == 0)
    t = tjd - speed_intv;
        else if (i == 1)
    t = tjd + speed_intv;
        else
    t = tjd;
  retc = swemoon(t, iflag | SEFLG_SPEED, NO_SAVE, xpos[i], serr, swed);/**/
  if (retc == ERR)
    return(ERR);
  /* light-time-corrected moon for apparent node (~ 0.006") */
  if (is_set(iflag,SEFLG_TRUEPOS) == 0 && retc >= OK) {
    dt = sqrt(square_sum(xpos[i])) * AUNIT / CLIGHT / 86400.0;
    retc = swemoon(t-dt, iflag | SEFLG_SPEED, NO_SAVE, xpos[i], serr, swed);/**/
    if (retc == ERR)
      return(ERR);
        }
  if (retc == NOT_AVAILABLE) {
    if (tjd > MOSHPLEPH_START && tjd < MOSHPLEPH_END) {
      iflag = (iflag & ~SEFLG_SWIEPH) | SEFLG_MOSEPH;
      epheflag = SEFLG_MOSEPH;
      if (serr != NULL && strlen(serr) + 30 < AS_MAXCH)
        strcat(serr, " \nusing Moshier eph.; ");
      break;
    } else
      return ERR;
  }
  /* precession and nutation etc. */
  retc = swi_plan_for_osc_elem(iflag|SEFLG_SPEED, t, xpos[i]); /* retc is always ok */
      }
      break;
    case SEFLG_MOSEPH:
      /* with moshier moon, we need a greater speed_intv, because here the
       * node and apogee oscillate wildly within small intervals */
      speed_intv = NODE_CALC_INTV_MOSH;
      for (i = istart; i <= 2; i++) {
  if (i == 0)
    t = tjd - speed_intv;
        else if (i == 1)
    t = tjd + speed_intv;
        else
    t = tjd;
  retc = swi_moshmoon(t, NO_SAVE, xpos[i], serr);/**/
  if (retc == ERR)
    return(retc);
  /* precession and nutation etc. */
  retc = swi_plan_for_osc_elem(iflag|SEFLG_SPEED, t, xpos[i]); /* retc is always ok */
      }
      break;
    default:
      break;
  }
  if (retc == NOT_AVAILABLE || retc == BEYOND_EPH_LIMITS)
    goto three_positions;
  /*********************************************
   * node with speed                           *
   *********************************************/
  /* node is always needed, even if apogee is wanted */
  ndnp = &swed->nddat[SEI_TRUE_NODE];
  /* three nodes */
  for (i = istart; i <= 2; i++) {
    if (fabs(xpos[i][5]) < 1e-15)
      xpos[i][5] = 1e-15;
    fac = xpos[i][2] / xpos[i][5];
    sgn = xpos[i][5] / fabs(xpos[i][5]);
    for (j = 0; j <= 2; j++)
      xx[i][j] = (xpos[i][j] - fac * xpos[i][j+3]) * sgn;
  }
  /* now we have the correct direction of the node, the
   * intersection of the lunar plane and the ecliptic plane.
   * the distance is the distance of the point where the tangent
   * of the lunar motion penetrates the ecliptic plane.
   * this can be very large, e.g. j2415080.37372.
   * below, a new distance will be derived from the osculating
   * ellipse.
   */
  /* save position and speed */
  for (i = 0; i <= 2; i++) {
    ndnp->x[i] = xx[2][i];
    if (is_set(iflag,SEFLG_SPEED)) {
      b = (xx[1][i] - xx[0][i]) / 2;
      a = (xx[1][i] + xx[0][i]) / 2 - xx[2][i];
      ndnp->x[i+3] = (2 * a + b) / speed_intv;
    } else
      ndnp->x[i+3] = 0;
    ndnp->teval = tjd;
    ndnp->iephe = epheflag;
  }
  /************************************************************
   * apogee with speed                                        *
   * must be computed anyway to get the node's distance       *
   ************************************************************/
  ndap = &swed->nddat[SEI_OSCU_APOG];
  Gmsm = GEOGCONST * (1 + 1 / EARTH_MOON_MRAT) /AUNIT/AUNIT/AUNIT*86400.0*86400.0;
  /* three apogees */
  for (i = istart; i <= 2; i++) {
    /* node */
    rxy =  sqrt(xx[i][0] * xx[i][0] + xx[i][1] * xx[i][1]);
    cosnode = xx[i][0] / rxy;
    sinnode = xx[i][1] / rxy;
    /* inclination */
    swi_cross_prod(xpos[i], xpos[i]+3, xnorm);
    rxy =  xnorm[0] * xnorm[0] + xnorm[1] * xnorm[1];
    c2 = (rxy + xnorm[2] * xnorm[2]);
    rxyz = sqrt(c2);
    rxy = sqrt(rxy);
    sinincl = rxy / rxyz;
    cosincl = sqrt(1 - sinincl * sinincl);
    /* argument of latitude */
    cosu = xpos[i][0] * cosnode + xpos[i][1] * sinnode;
    sinu = xpos[i][2] / sinincl;
    uu = atan2(sinu, cosu);
    /* semi-axis */
    rxyz = sqrt(square_sum(xpos[i]));
    v2 = square_sum((xpos[i]+3));
    sema = 1 / (2 / rxyz - v2 / Gmsm);
    /* eccentricity */
    pp = c2 / Gmsm;
    ecce = sqrt(1 - pp / sema);
    /* eccentric anomaly */
    cosE = 1 / ecce * (1 - rxyz / sema);
    sinE = 1 / ecce / sqrt(sema * Gmsm) * dot_prod(xpos[i], (xpos[i]+3));
    /* true anomaly */
    ny = 2 * atan(sqrt((1+ecce)/(1-ecce)) * sinE / (1 + cosE));
    /* distance of apogee from ascending node */
    xxa[i][0] = swi_mod2PI(uu - ny + PI);
    xxa[i][1] = 0;      /* latitude */
    xxa[i][2] = sema * (1 + ecce);  /* distance */
    /* transformation to ecliptic coordinates */
    swi_polcart(xxa[i], xxa[i]);
    swi_coortrf2(xxa[i], xxa[i], -sinincl, cosincl);
    swi_cartpol(xxa[i], xxa[i]);
    /* adding node, we get apogee in ecl. coord. */
    xxa[i][0] += atan2(sinnode, cosnode);
    swi_polcart(xxa[i], xxa[i]);
    /* new distance of node from orbital ellipse:
     * true anomaly of node: */
    ny = swi_mod2PI(ny - uu);
    /* eccentric anomaly */
    cosE = cos(2 * atan(tan(ny / 2) / sqrt((1+ecce) / (1-ecce))));
    /* new distance */
    r[0] = sema * (1 - ecce * cosE);
    /* old node distance */
    r[1] = sqrt(square_sum(xx[i]));
    /* correct length of position vector */
    for (j = 0; j <= 2; j++)
      xx[i][j] *= r[0] / r[1];
  }
  /* save position and speed */
  for (i = 0; i <= 2; i++) {
    /* apogee */
    ndap->x[i] = xxa[2][i];
    if (is_set(iflag,SEFLG_SPEED))
      ndap->x[i+3] = (xxa[1][i] - xxa[0][i]) / speed_intv / 2;
    else
      ndap->x[i+3] = 0;
    ndap->teval = tjd;
    ndap->iephe = epheflag;
    /* node */
    ndnp->x[i] = xx[2][i];
    if (is_set(iflag,SEFLG_SPEED))
      ndnp->x[i+3] = (xx[1][i] - xx[0][i]) / speed_intv / 2;/**/
    else
      ndnp->x[i+3] = 0;
  }
  /**********************************************************************
   * precession and nutation have already been taken into account
   * because the computation is on the basis of lunar positions
   * that have gone through swi_plan_for_osc_elem.
   * light-time is already contained in lunar positions.
   * now compute polar and equatorial coordinates:
   **********************************************************************/
  for (j = 0; j <= 1; j++) {
    double x[6];
    if (j == 0)
      ndp = &swed->nddat[SEI_TRUE_NODE];
    else
      ndp = &swed->nddat[SEI_OSCU_APOG];
    memset((void *) ndp->xreturn, 0, 24 * sizeof(double));
    /* cartesian ecliptic */
    vec_copy(ndp->xreturn+6,ndp->x,6);
    /* polar ecliptic */
    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
    /* cartesian equatorial */
    swi_coortrf2(ndp->xreturn+6, ndp->xreturn+18, -oe->seps, oe->ceps);
    if (is_set(iflag,SEFLG_SPEED))
      swi_coortrf2(ndp->xreturn+9, ndp->xreturn+21, -oe->seps, oe->ceps);
#ifdef SID_TNODE_FROM_ECL_T0
    /* sideral: we return NORMAL equatorial coordinates, there are no
     * sidereal ones */
    if (is_set(iflag,SEFLG_SIDEREAL)) {
      /* to J2000 */
      swi_precess(ndp->xreturn+18, sip->t0, J_TO_J2000);
      if (is_set(iflag,SEFLG_SPEED))
  swi_precess_speed(ndp->xreturn+21, sip->t0, J_TO_J2000);
      if (is_not_set(iflag,SEFLG_J2000)) {
  /* to tjd */
  swi_precess(ndp->xreturn+18, tjd, J2000_TO_J);
  if (is_set(iflag,SEFLG_SPEED))
    swi_precess_speed(ndp->xreturn+21, tjd, J2000_TO_J);
      }
    }
#endif
    if (is_not_set(iflag,SEFLG_NONUT)) {
      swi_coortrf2(ndp->xreturn+18, ndp->xreturn+18, -swed->nut.snut, swed->nut.cnut);
      if (is_set(iflag,SEFLG_SPEED))
  swi_coortrf2(ndp->xreturn+21, ndp->xreturn+21, -swed->nut.snut, swed->nut.cnut);
    }
    /* polar equatorial */
    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    ndp->xflgs = iflag;
    ndp->iephe = iflag & SEFLG_EPHMASK;
#ifdef SID_TNODE_FROM_ECL_T0
    /* node and apogee are already referred to t0;
     * nothing has to be done */
#else
    if (is_set(iflag,SEFLG_SIDEREAL)) {
      /* node and apogee are referred to t;
       * the ecliptic position must be transformed to t0 */
      /* rigorous algorithm */
      if ((swed->sidd.sid_mode & SE_SIDBIT_ECL_T0)
        || (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE)) {
  vec_copy(x,ndp->xreturn+18,6);

  /* remove nutation */
  if (is_not_set(iflag,SEFLG_NONUT))
    swi_nutate(x, iflag, TRUE);
  /* precess to J2000 */
  swi_precess(x, tjd, J_TO_J2000);
  if (is_set(iflag,SEFLG_SPEED))
    swi_precess_speed(x, tjd, J_TO_J2000);
        if (swed->sidd.sid_mode & SE_SIDBIT_ECL_T0)
    swi_trop_ra2sid_lon(x, ndp->xreturn+6, ndp->xreturn+18, iflag, NULL);
        /* project onto solar system equator */
        else if (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE)
          swi_trop_ra2sid_lon_sosy(x, ndp->xreturn+6, ndp->xreturn+18, iflag, NULL);
  /* to polar */
  swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
        swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
      /* traditional algorithm;
       * this is a bit clumsy, but allows us to keep the
       * sidereal code together */
      } else {
  swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
  ndp->xreturn[0] -= swe_get_ayanamsa(ndp->teval) * DEGTORAD;
  swi_polcart_sp(ndp->xreturn, ndp->xreturn+6);
      }
    } else if (is_set(iflag,SEFLG_J2000)) {
      /* node and apogee are referred to t;
       * the ecliptic position must be transformed to J2000 */
      vec_copy(x,ndp->xreturn+18,6);
      /* precess to J2000 */
      swi_precess(x, tjd, J_TO_J2000);
      if (is_set(iflag,SEFLG_SPEED))
        swi_precess_speed(x, tjd, J_TO_J2000);
      vec_copy(ndp->xreturn+18,x,6);  
      swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
      swi_coortrf2(ndp->xreturn+18, ndp->xreturn+6, swed->oec2000.seps, swed->oec2000.ceps);
      if (is_set(iflag,SEFLG_SPEED))
        swi_coortrf2(ndp->xreturn+21, ndp->xreturn+9, swed->oec2000.seps, swed->oec2000.ceps);
      swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
    }
#endif
    /**********************
     * radians to degrees *
     **********************/
    convert_arcs( ndp->xreturn,    RADTODEG);
    convert_arcs( ndp->xreturn+12, RADTODEG);
 
    ndp->xreturn[0]  = swe_degnorm(ndp->xreturn[0]);
    ndp->xreturn[12] = swe_degnorm(ndp->xreturn[12]);

  }
  return OK;
}

/* lunar osculating elements, i.e.
 */
static int intp_apsides(double tjd, int ipl, int iflag, char *serr, struct swe_data *swed)
{
  int i;
  int flg1, flg2;
  struct plan_data *ndp;
  struct epsilon *oe;
  struct nut *nut;
  double speed_intv = 0.1;
  double t, dt;
  double xpos[3][6], xx[6], x[6];
  int speedf1, speedf2;
  oe = &swed->oec;
  nut = &swed->nut;
  ndp = &swed->nddat[ipl];
  /* if same calculation was done before, return
   * if speed flag has been turned on, recompute */
  flg1 = iflag & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  flg2 = ndp->xflgs & ~SEFLG_EQUATORIAL & ~SEFLG_XYZ;
  speedf1 = ndp->xflgs & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (tjd == ndp->teval
  && tjd != 0
  && flg1 == flg2
  && (!speedf2 || speedf1)) {
    ndp->xflgs = iflag;
    ndp->iephe = iflag & SEFLG_MOSEPH;
    return OK;
  }
  /*********************************************
   * now three apsides *
   *********************************************/
  for (t = tjd - speed_intv, i = 0; i < 3; t += speed_intv, i++) {
    if (is_not_set(iflag, SEFLG_SPEED) && i != 1) continue;
    swi_intp_apsides(t, xpos[i], ipl);
  }
  /************************************************************
   * apsis with speed                                         *
   ************************************************************/
  for (i = 0; i < 3; i++) {
    xx[i] = xpos[1][i];
    xx[i+3] = 0;
  }
  if (is_set(iflag,SEFLG_SPEED)) {
    xx[3] = swe_difrad2n(xpos[2][0], xpos[0][0]) / speed_intv / 2.0;
    xx[4] = (xpos[2][1] - xpos[0][1]) / speed_intv / 2.0;
    xx[5] = (xpos[2][2] - xpos[0][2]) / speed_intv / 2.0;
  }
  memset((void *) ndp->xreturn, 0, 24 * sizeof(double));
  /* ecliptic polar to cartesian */
  swi_polcart_sp(xx, xx);
  /* light-time */
  if (is_not_set(iflag,SEFLG_TRUEPOS)) {
    dt = sqrt(square_sum(xx)) * AUNIT / CLIGHT / 86400.0;
    for (i = 1; i < 3; i++)
      xx[i] -= dt * xx[i+3];
  }
  vec_copy(ndp->xreturn+6,xx,6);

  /* equatorial cartesian */
  swi_coortrf2(ndp->xreturn+6, ndp->xreturn+18, -oe->seps, oe->ceps);
  if (is_set(iflag,SEFLG_SPEED))
    swi_coortrf2(ndp->xreturn+9, ndp->xreturn+21, -oe->seps, oe->ceps);
  ndp->teval = tjd;
  ndp->xflgs = iflag;
  ndp->iephe = iflag & SEFLG_EPHMASK;
  if (is_set(iflag,SEFLG_SIDEREAL)) {
    /* apogee is referred to t;
     * the ecliptic position must be transformed to t0 */
    /* rigorous algorithm */
    if ((swed->sidd.sid_mode & SE_SIDBIT_ECL_T0)
  || (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE)) {
      vec_copy(x,ndp->xreturn+18,6);
      /* precess to J2000 */
      swi_precess(x, tjd, J_TO_J2000);
      if (is_set(iflag,SEFLG_SPEED))
  swi_precess_speed(x, tjd, J_TO_J2000);
      if (swed->sidd.sid_mode & SE_SIDBIT_ECL_T0)
  swi_trop_ra2sid_lon(x, ndp->xreturn+6, ndp->xreturn+18, iflag, NULL);
      /* project onto solar system equator */
      else if (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE)
  swi_trop_ra2sid_lon_sosy(x, ndp->xreturn+6, ndp->xreturn+18, iflag, NULL);
      /* to polar */
      swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
      swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    } else {
    /* traditional algorithm */
      swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
      ndp->xreturn[0] -= swe_get_ayanamsa(ndp->teval) * DEGTORAD;
      swi_polcart_sp(ndp->xreturn, ndp->xreturn+6);
      swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    }
  } else if (is_set(iflag,SEFLG_J2000)) {
    /* node and apogee are referred to t;
     * the ecliptic position must be transformed to J2000 */
    vec_copy(x,ndp->xreturn+18,6);
    /* precess to J2000 */
    swi_precess(x, tjd, J_TO_J2000);
    if (is_set(iflag,SEFLG_SPEED))
      swi_precess_speed(x, tjd, J_TO_J2000);
    vec_copy(ndp->xreturn+18,x,6);      
    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    swi_coortrf2(ndp->xreturn+18, ndp->xreturn+6, swed->oec2000.seps, swed->oec2000.ceps);
    if (is_set(iflag,SEFLG_SPEED))
      swi_coortrf2(ndp->xreturn+21, ndp->xreturn+9, swed->oec2000.seps, swed->oec2000.ceps);
    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
  } else {
    /* tropical ecliptic positions */
    /* precession has already been taken into account, but not nutation */
    if (is_not_set(iflag,SEFLG_NONUT)) {
      swi_nutate(ndp->xreturn+18, iflag, FALSE);
    }
    /* equatorial polar */
    swi_cartpol_sp(ndp->xreturn+18, ndp->xreturn+12);
    /* ecliptic cartesian */
    swi_coortrf2(ndp->xreturn+18, ndp->xreturn+6, oe->seps, oe->ceps);
    if (is_set(iflag,SEFLG_SPEED))
      swi_coortrf2(ndp->xreturn+21, ndp->xreturn+9, oe->seps, oe->ceps);
    if (is_not_set(iflag,SEFLG_NONUT)) {
      swi_coortrf2(ndp->xreturn+6, ndp->xreturn+6, nut->snut, nut->cnut);
      if (is_set(iflag,SEFLG_SPEED))
        swi_coortrf2(ndp->xreturn+9, ndp->xreturn+9, nut->snut, nut->cnut);
      }
    /* ecliptic polar */
    swi_cartpol_sp(ndp->xreturn+6, ndp->xreturn);
  }
  /**********************
   * radians to degrees *
   **********************/
   convert_arcs( ndp->xreturn,    RADTODEG);
   convert_arcs( ndp->xreturn+12, RADTODEG);
  
   ndp->xreturn[0] = swe_degnorm(ndp->xreturn[0]);
   ndp->xreturn[12] = swe_degnorm(ndp->xreturn[12]);
  
  return OK;
}

/* transforms the position of the moon in a way we can use it
 * for calculation of osculating node and apogee:
 * precession and nutation (attention to speed vector!)
 * according to flags
 * iflag  flags
 * tjd          time for which the element is computed
 *              i.e. date of ecliptic
 * xx           array equatorial cartesian position and speed
 * serr         error string
 */
int swi_plan_for_osc_elem(int iflag, double tjd, double *xx)
{
  return swid_plan_for_osc_elem(iflag,tjd,xx,&swed);
}

int swid_plan_for_osc_elem(int iflag, double tjd, double *xx, struct swe_data *swed)
{
  int i;
  double x[6];
  struct nut nuttmp;
  struct nut *nutp = &nuttmp;  /* dummy assign, to silence gcc warning */
  struct epsilon *oe = &swed->oec;
  struct epsilon oectmp;
  /* ICRS to J2000 */
  if (is_not_set(iflag,SEFLG_ICRS) && swed->jpldenum >= 403) {
    swi_bias(xx, iflag, FALSE);
  }/**/
  /************************************************
   * precession, equator 2000 -> equator of date  *
   * attention: speed vector has to be rotated,   *
   * but daily precession 0.137" may not be added!*/
#ifdef SID_TNODE_FROM_ECL_T0
  struct sid_data *sip = &swed->sidd;
  /* For sidereal calculation we need node refered*
   * to ecliptic of t0 of ayanamsa                *
   ************************************************/
  if (is_set(iflag,SEFLG_SIDEREAL)) {
    tjd = sip->t0;
    swi_precess(xx, tjd, J2000_TO_J);
    swi_precess(xx+3, tjd, J2000_TO_J);
    calc_epsilon(tjd, &oectmp);
    oe = &oectmp;
  } else if (is_not_set(iflag,SEFLG_J2000)) {
#endif
    swi_precess(xx, tjd, J2000_TO_J);
    swi_precess(xx+3, tjd, J2000_TO_J);
    /* epsilon */
    if (tjd == swed->oec.teps)
      oe = &swed->oec;
    else if (tjd == J2000)
      oe = &swed->oec2000;
    else  {
      calc_epsilon(tjd, &oectmp);
      oe = &oectmp;
    }
#ifdef SID_TNODE_FROM_ECL_T0
  } else  /* if SEFLG_J2000 */
    oe = &swed->oec2000;
#endif
  /************************************************
   * nutation                                     *
   * again: speed vector must be rotated, but not *
   * added 'speed' of nutation                    *
   ************************************************/
  if (is_not_set(iflag,SEFLG_NONUT)) {
    if (tjd == swed->nut.tnut)
      nutp = &swed->nut;
    else if (tjd == J2000)
      nutp = &swed->nut2000;
    else if (tjd == swed->nutv.tnut)
      nutp = &swed->nutv;
    else {
      nutp = &nuttmp;
      swi_nutation(tjd, nutp->nutlo);
      nutp->tnut = tjd;
      nutp->snut = sin(nutp->nutlo[1]);
      nutp->cnut = cos(nutp->nutlo[1]);
      nut_matrix(nutp, oe);
    }
    for (i = 0; i <= 2; i++)
      x[i] = xx[0] * nutp->matrix[0][i] +
       xx[1] * nutp->matrix[1][i] +
       xx[2] * nutp->matrix[2][i];
    /* speed:
     * rotation only */
    for (i = 0; i <= 2; i++)
      x[i+3] = xx[3] * nutp->matrix[0][i] +
         xx[4] * nutp->matrix[1][i] +
         xx[5] * nutp->matrix[2][i];
    vec_copy(xx,x,6);
  }
  /************************************************
   * transformation to ecliptic                   *
   ************************************************/
  swi_coortrf2(xx, xx, oe->seps, oe->ceps);
  swi_coortrf2(xx+3, xx+3, oe->seps, oe->ceps);
#ifdef SID_TNODE_FROM_ECL_T0
  if (is_set(iflag,SEFLG_SIDEREAL)) {
    /* subtract ayan_t0 */
    swi_cartpol_sp(xx, xx);
    xx[0] -= sip->ayan_t0;
    swi_polcart_sp(xx, xx);
  } else
#endif
  if (is_not_set(iflag,SEFLG_NONUT)) {
    swi_coortrf2(xx, xx, nutp->snut, nutp->cnut);
    swi_coortrf2(xx+3, xx+3, nutp->snut, nutp->cnut);
  }
  return(OK);
}

static double meff(double r)
{
  double f, m;
  int i;
  if (r <= 0)
    return 0.0;
  else if (r >= 1)
    return 1.0;
  for (i = 0; eff_arr[i].r > r; i++)
    ;  /* empty body */
  f = (r - eff_arr[i-1].r) / (eff_arr[i].r - eff_arr[i-1].r);
  m = eff_arr[i-1].m + f * (eff_arr[i].m - eff_arr[i-1].m);
  return m;
}

void swi_check_ecliptic(double tjd)
{
  swid_check_ecliptic(tjd, &swed);
}

void swid_check_ecliptic(double tjd, struct swe_data *swed)
{
  if (swed->oec2000.teps != J2000) {
    calc_epsilon(J2000, &swed->oec2000);
  }
  if (tjd == J2000) {
    swed->oec.teps = swed->oec2000.teps;
    swed->oec.eps = swed->oec2000.eps;
    swed->oec.seps = swed->oec2000.seps;
    swed->oec.ceps = swed->oec2000.ceps;
    return;
  }
  if (swed->oec.teps != tjd || tjd == 0) {
    calc_epsilon(tjd, &swed->oec);
  }
}

/* computes nutation, if it is wanted and has not yet been computed.
 * if speed flag has been turned on since last computation,
 * nutation is recomputed */
void swi_check_nutation(double tjd, int iflag)
{
 swid_check_nutation(tjd,iflag,&swed);
}
void swid_check_nutation(double tjd, int iflag, struct swe_data *swed)
{
  int speedf1, speedf2;
  static int nutflag = 0;
  double t;
  speedf1 = nutflag & SEFLG_SPEED;
  speedf2 = iflag & SEFLG_SPEED;
  if (is_not_set(iflag,SEFLG_NONUT)
  && (tjd != swed->nut.tnut || tjd == 0
  || (!speedf1 && speedf2))) {
    swi_nutation(tjd, swed->nut.nutlo);
    swed->nut.tnut = tjd;
    swed->nut.snut = sin(swed->nut.nutlo[1]);
    swed->nut.cnut = cos(swed->nut.nutlo[1]);
    nutflag = iflag;
    nut_matrix(&swed->nut, &swed->oec);
    if (is_set(iflag,SEFLG_SPEED)) {
      /* once more for 'speed' of nutation, which is needed for
       * planetary speeds */
      t = tjd - NUT_SPEED_INTV;
      swi_nutation(t, swed->nutv.nutlo);
      swed->nutv.tnut = t;
      swed->nutv.snut = sin(swed->nutv.nutlo[1]);
      swed->nutv.cnut = cos(swed->nutv.nutlo[1]);
      nut_matrix(&swed->nutv, &swed->oec);
    }
  }
}

static int plaus_ipl(double tjd, int ipl_in, int iflag, char* serr) {
  int ipl = ipl_in;

  /* function calls for Pluto with asteroid number 134340
   * are treated as calls for Pluto as main body SE_PLUTO.
   * Reason: Our numerical integrator takes into account Pluto
   * perturbation and therefore crashes with body 134340 Pluto. */
  if (ipl == SE_AST_OFFSET + 134340) ipl = SE_PLUTO;

  if (ipl == SE_CHIRON && (tjd < CHIRON_START || tjd > CHIRON_END)) {
    if (serr != NULL)
      sprintf(serr, "Chiron's ephemeris is restricted to JD %8.1f - JD %8.1f",
       CHIRON_START, CHIRON_END);
      return ERR;
      }
  if (ipl == SE_PHOLUS && tjd < PHOLUS_START) {
    if (serr != NULL)
      sprintf(serr,
       "Pholus's ephemeris is restricted to the time after JD %8.1f",
       PHOLUS_START);
    return ERR;
    }

  return ipl;
  }

static int plaus_iflag(int iflag_in, char* serr)
{

  int iflag = iflag_in;

  if (is_set(iflag,SEFLG_SPEED3)) {
    strcpy(serr,"Warning: SEFLG_SPEED3 not supported in this version - using SEFLG_SPEED");
    iflag = ( iflag & ~SEFLG_SPEED3 ) | SEFLG_SPEED;
  }

  /* high precision speed prevails fast speed */
  if (is_set(iflag,SEFLG_SPEED3 & SEFLG_SPEED))
    iflag = iflag & ~SEFLG_SPEED3;

  /* cartesian flag excludes radians flag */
  if (is_set(iflag,SEFLG_XYZ & SEFLG_RADIANS))
    iflag = iflag & ~SEFLG_RADIANS;

  /* if topocentric bit, turn helio- and barycentric bits off */
  if (is_set(iflag,SEFLG_TOPOCTR))
    iflag = iflag & ~(SEFLG_HELCTR | SEFLG_BARYCTR);
    
  /* if heliocentric bit, turn aberration and deflection off */
  if (is_set(iflag,SEFLG_HELCTR))
    iflag |= SEFLG_NOABERR | SEFLG_NOGDEFL; 
    
  /* same, if barycentric bit */
  if (is_set(iflag,SEFLG_BARYCTR))
    iflag |= SEFLG_NOABERR | SEFLG_NOGDEFL; 
    
  /* if no_precession bit is set, set also no_nutation bit */
  if (is_set(iflag,SEFLG_J2000))
    iflag |= SEFLG_NONUT;
    
  /* if sidereal bit is set, set also no_nutation bit */
  if (is_set(iflag,SEFLG_SIDEREAL))
    iflag |= SEFLG_NONUT;
    
  /* if truepos is set, turn off grav. defl. and aberration */
  if (is_set(iflag,SEFLG_TRUEPOS))
    iflag |= (SEFLG_NOGDEFL | SEFLG_NOABERR);

  /* Determine a unique ephemeris version to use */
  int epheflag = iflag & SEFLG_EPHMASK;
  if (epheflag & SEFLG_JPLEPH) {
    strcpy(serr,"Warning: SEFLG_JPLEPH not supported in this version - using SEFLG_SWIEPH");
    epheflag = SEFLG_SWIEPH;
    }
  else if (epheflag & SEFLG_MOSEPH)
    epheflag = SEFLG_MOSEPH;
  else
    epheflag = SEFLG_SWIEPH; // includes default case (epheflag == 0)
  iflag = (iflag & ~SEFLG_EPHMASK) | epheflag;

  /* no barycentric calculations with Moshier ephemeris */
  if (is_set(iflag,SEFLG_BARYCTR & SEFLG_MOSEPH)) {
    msg_append(serr,"Barycentric Moshier not supported / resetting barycentric flag");
    iflag &= ~SEFLG_BARYCTR;
  }

  return iflag;
}

/**********************************************************
 * get fixstar positions
 * parameters:
 * star   name of star or line number in star file
 *    (start from 1, don't count comment).
 *        If no error occurs, the name of the star is returned
 *          in the format trad_name, nomeclat_name
 *
 * tjd     absolute julian day
 * iflag  s. swecalc(); speed bit does not function
 * x    pointer for returning the ecliptic coordinates
 * serr    error return string
**********************************************************/
int FAR PASCAL_CONV swed_fixstar(char *star, double tjd, int iflag,
  double *xx, char *serr, struct swe_data *swed)
{
  int i;
  int star_nr = 0;
  AS_BOOL  isnomclat = FALSE;
  size_t cmplen;
  double x[6], xxsv[6], xobs[6], *xpo = NULL;
  char *cpos[20];
  char sstar[SE_MAX_STNAME + 1];
  char fstar[SE_MAX_STNAME + 1];
  static char slast_stardata[AS_MAXCH];
  static char slast_starname[AS_MAXCH];
  char s[AS_MAXCH + 20], *sp, *sp2;  /* 20 byte for SE_STARFILE */
  double ra_s, ra_pm, de_pm, ra, de, t, cosra, cosde, sinra, sinde;
  double ra_h, ra_m, de_d, de_m, de_s;
  char *sde_d;
  double epoch, radv, parall, u;
  int line = 0;
  int fline = 0;
  struct plan_data *pedp = &swed->pldat[SEI_EARTH];
  struct plan_data *psdp = &swed->pldat[SEI_SUNBARY];
  struct epsilon *oe = &swed->oec2000;
  int retc;
  int epheflag, iflgsave;
  iflag |= SEFLG_SPEED; /* we need this in order to work correctly */
  iflgsave = iflag;
  if (serr != NULL)
    *serr = '\0';
  iflag = plaus_iflag(iflag, serr);
  if (iflag & SEFLG_SIDEREAL && !swed->ayana_is_set)
    swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY, 0, 0);
  epheflag = iflag & SEFLG_EPHMASK;
  /******************************************
   * obliquity of ecliptic 2000 and of date *
   ******************************************/
  swi_check_ecliptic(tjd);
  /******************************************
   * nutation                               *
   ******************************************/
  swi_check_nutation(tjd, iflag);
  strncpy(sstar, star, SE_MAX_STNAME);
  sstar[SE_MAX_STNAME] = '\0';
  if (*sstar == ',') {
    isnomclat = TRUE;
  } else if (isdigit((int) *sstar)) {
    star_nr = atoi(sstar);
  } else {
    /* traditional name of star to lower case */
    for (sp = sstar; *sp != '\0'; sp++)
      *sp = tolower((int) *sp);
    if ((sp = strchr(sstar, ',')) != NULL)
      *sp = '\0';
  }
  /*swi_right_trim(sstar);*/
  while ((sp = strchr(sstar, ' ')) != NULL)
    swi_strcpy(sp, sp+1);
  cmplen = strlen(sstar);
  if (cmplen == 0) {
    if (serr != NULL)
      sprintf(serr, "swe_fixstar(): star name empty");
    retc = ERR;
    goto return_err;
  }
  /* star elements from last call: */
  if (*slast_stardata != '\0' && strcmp(slast_starname, sstar) == 0) {
    strcpy(s, slast_stardata);
    goto found;
  }
  /******************************************************
   * Star file
   * close to the beginning, a few stars selected by Astrodienst.
   * These can be accessed by giving their number instead of a name.
   * All other stars can be accessed by name.
   * Comment lines start with # and are ignored.
   ******************************************************/
  if (swed->fixfp == NULL) {
    if ((swed->fixfp = swid_fopen(SEI_FILE_FIXSTAR, SE_STARFILE, swed->ephepath, serr, swed)) == NULL) {
      swed->is_old_starfile = TRUE;
      if ((swed->fixfp = swid_fopen(SEI_FILE_FIXSTAR, SE_STARFILE_OLD, swed->ephepath, NULL, swed)) == NULL) {
  retc = ERR;
  goto return_err;
      }
    }
  }
  rewind(swed->fixfp);
  while (fgets(s, AS_MAXCH, swed->fixfp) != NULL) {
    fline++;
    if (*s == '#') continue;
    line++;
    if (star_nr == line)
      goto found;
    else if (star_nr > 0)
      continue;
    if ((sp = strchr(s, ',')) == NULL) {
      if (serr != NULL)
  sprintf(serr, "star file %s damaged at line %d", SE_STARFILE, fline);
      retc = ERR;
      goto return_err;
    }
    if (isnomclat) {
      if (strncmp(sp, sstar, cmplen) == 0)
        goto found;
      else
        continue;
    }
    *sp = '\0';  /* cut off first field */
    strncpy(fstar, s, SE_MAX_STNAME);
    *sp = ',';
    fstar[SE_MAX_STNAME] = '\0';  /* force termination */
    /*swi_right_trim(fstar);*/
    while ((sp = strchr(fstar, ' ')) != NULL)
      swi_strcpy(sp, sp+1);
    i = strlen(fstar);
    if (i < (int) cmplen)
      continue;
    for (sp2 = fstar; *sp2 != '\0'; sp2++) {
      *sp2 = tolower((int) *sp2);
    }
    if (strncmp(fstar, sstar, cmplen) == 0)
      goto found;
  }
  if (serr != NULL && strlen(star) < AS_MAXCH - 20)
    sprintf(serr, "star %s not found", star);
  retc = ERR;
  goto return_err;
  found:
  strcpy(slast_stardata, s);
  strcpy(slast_starname, sstar);
  i = swi_cutstr(s, ",", cpos, 20);
  swi_right_trim(cpos[0]);
  swi_right_trim(cpos[1]);
  if (i < 13) {
    if (serr != NULL)
      sprintf(serr, "data of star '%s,%s' incomplete", cpos[0], cpos[1]);
    retc = ERR;
    goto return_err;
  }
  epoch = atof(cpos[2]);
  ra_h = atof(cpos[3]);
  ra_m = atof(cpos[4]);
  ra_s = atof(cpos[5]);
  de_d = atof(cpos[6]);
  sde_d = cpos[6];
  de_m = atof(cpos[7]);
  de_s = atof(cpos[8]);
  ra_pm = atof(cpos[9]);
  de_pm = atof(cpos[10]);
  radv = atof(cpos[11]);
  parall = atof(cpos[12]);
  /* return trad. name, nomeclature name */
  if (strlen(cpos[0]) > SE_MAX_STNAME)
    cpos[0][SE_MAX_STNAME] = '\0';
  if (strlen(cpos[1]) > SE_MAX_STNAME-1)
    cpos[1][SE_MAX_STNAME-1] = '\0';
  sprintf(star, "%s,%s", cpos[0], cpos[1]);
  /****************************************
   * position and speed (equinox)
   ****************************************/
  /* ra and de in degrees */
  ra = (ra_s / 3600.0 + ra_m / 60.0 + ra_h) * 15.0;
  if (strchr(sde_d, '-') == NULL)
    de = de_s / 3600.0 + de_m / 60.0 + de_d;
  else
    de = -de_s / 3600.0 - de_m / 60.0 + de_d;
  /* speed in ra and de, degrees per century */
  if (swed->is_old_starfile == TRUE) {
    ra_pm = ra_pm * 15 / 3600.0;
    de_pm = de_pm / 3600.0;
  } else {
    ra_pm = ra_pm / 10.0 / 3600.0;
    de_pm = de_pm / 10.0 / 3600.0;
    parall /= 1000.0;
  }
  /* parallax, degrees */
  if (parall > 1)
    parall = (1 / parall / 3600.0);
  else
    parall /= 3600;
  /* radial velocity in AU per century */
  radv *= KM_S_TO_AU_CTY;

  /* radians */
  ra *= DEGTORAD;
  de *= DEGTORAD;
  ra_pm *= DEGTORAD;
  de_pm *= DEGTORAD;
  ra_pm /= cos(de); /* catalogues give proper motion in RA as great circle */
  parall *= DEGTORAD;
  x[0] = ra;
  x[1] = de;
  x[2] = 1;  /* -> unit vector */
  /* cartesian */
  swi_polcart(x, x);
  /*space motion vector */
  cosra = cos(ra);
  cosde = cos(de);
  sinra = sin(ra);
  sinde = sin(de);
  x[3] = -ra_pm * cosde * sinra - de_pm * sinde * cosra
      + radv * parall * cosde * cosra;
  x[4] = ra_pm * cosde * cosra - de_pm * sinde * sinra
      + radv * parall * cosde * sinra;
  x[5] = de_pm * cosde + radv * parall * sinde;
  x[3] /= 36525;
  x[4] /= 36525;
  x[5] /= 36525;
  /******************************************
   * FK5
   ******************************************/
  if (epoch == 1950) {
    swi_FK4_FK5(x, B1950);
    swi_precess(x, B1950, J_TO_J2000);
    swi_precess(x+3, B1950, J_TO_J2000);
  }
  /* FK5 to ICRS, if jpl ephemeris is referred to ICRS
   * With data that are already ICRS, epoch = 0 */
  if (epoch != 0) {
    swi_icrs2fk5(x, iflag, TRUE);
    /* with ephemerides < DE403, we now convert to J2000 */
    if (swed->jpldenum < 403)
      swi_bias(x, iflag, FALSE);
  }
  /****************************************************
   * earth/sun
   * for parallax, light deflection, and aberration,
   ****************************************************/
   if (is_not_set(iflag,SEFLG_BARYCTR) && 
       ( is_not_set(iflag,SEFLG_HELCTR) || is_not_set(iflag, SEFLG_MOSEPH) ) ) {
       
    if ((retc = main_planet(tjd, SEI_EARTH, epheflag, iflag, serr, swed)) != OK) {
      /*retc = ERR;
      goto return_err;*/
      iflag &= ~(SEFLG_TOPOCTR|SEFLG_HELCTR);
      /* on error, we provide barycentric position: */
      iflag |= SEFLG_BARYCTR | SEFLG_TRUEPOS | SEFLG_NOGDEFL;
      retc = iflag;
    } else {
      /* iflag (ephemeris bit) may have changed in main_planet() */
      iflag = swed->pldat[SEI_EARTH].xflgs;
    }
  }
  /************************************
   * observer: geocenter or topocenter
   ************************************/
  /* if topocentric position is wanted  */
  if (is_set(iflag,SEFLG_TOPOCTR)) {
    if (swed->topd.teval != pedp->teval
      || swed->topd.teval == 0) {
      if (swi_get_observer(pedp->teval, iflag, DO_SAVE, xobs, serr) != OK)
        return ERR;
    } else {
      vec_copy(xobs,swed->topd.xobs,6);
    }
    /* barycentric position of observer */
    vec_add(xobs,pedp->x,6);
  } else if (is_not_set(iflag,SEFLG_BARYCTR) && (is_not_set(iflag,SEFLG_HELCTR) ||is_not_set(iflag,SEFLG_MOSEPH))) {
    /* barycentric position of geocenter */
    vec_copy(xobs,pedp->x,6);
  }
  /************************************
   * position and speed at tjd        *
   ************************************/
  if (epoch == 1950)
    t= (tjd - B1950);  /* days since 1950.0 */
  else /* epoch == 2000 */
    t= (tjd - J2000);  /* days since 2000.0 */
  /* for parallax */
  if (is_set(iflag,SEFLG_HELCTR & SEFLG_MOSEPH))
    xpo = NULL;    /* no parallax, if moshier and heliocentric */
  else if (is_set(iflag,SEFLG_HELCTR))
    xpo = psdp->x;
  else if (is_set(iflag,SEFLG_BARYCTR))
    xpo = NULL;    /* no parallax, if barycentric */
  else
    xpo = xobs;
  if (xpo == NULL) {
    for (i = 0; i <= 2; i++)
      x[i] += t * x[i+3];
  } else {
    for (i = 0; i <= 2; i++) {
      x[i] += t * x[i+3] - parall * xpo[i];
      x[i+3] -= parall * xpo[i+3];
    }
  }
  
  /************************************
   * relativistic deflection of light *
   ************************************/
  smul(x,10000,6); /* great distance, to allow
       * algorithm used with planets */
       
  if (is_not_set(iflag, SEFLG_TRUEPOS | SEFLG_NOGDEFL)) {
    swi_deflect_light(x, 0, iflag & SEFLG_SPEED);
  }
  
  /**********************************
   * 'annual' aberration of light   *
   * speed is incorrect !!!         *
   **********************************/
  if (is_not_set(iflag,SEFLG_TRUEPOS | SEFLG_NOABERR))
    swi_aberr_light(x, xpo, iflag & SEFLG_SPEED);
    
  /* ICRS to J2000 */
  if (is_not_set(iflag,SEFLG_ICRS) && (swed->jpldenum >= 403 || is_set(iflag, SEFLG_BARYCTR) )) 
    swi_bias(x, iflag, FALSE);

  /* save J2000 coordinates; required for sidereal positions */
  vec_copy(xxsv,x,6);

  /************************************************
   * precession, equator 2000 -> equator of date *
   ************************************************/

  if (is_not_set(iflag, SEFLG_J2000)) {
    swi_precess(x, tjd, J2000_TO_J);
    if (is_set(iflag,SEFLG_SPEED))
      swi_precess_speed(x, tjd, J2000_TO_J);
    oe = &swed->oec;
  } else
    oe = &swed->oec2000;
  /************************************************
   * nutation                                     *
   ************************************************/
  if (is_not_set(iflag,SEFLG_NONUT))
    swi_nutate(x, 0, FALSE);

   /************************************************
   * unit vector (distance = 1)                   *
   ************************************************/
  u = sqrt(square_sum(x));
  sdiv(x,u,6);
  u = sqrt(square_sum(xxsv));
  sdiv(xxsv,u,6);

  /************************************************
   * set speed = 0, because not correct (aberration)
   ************************************************/
  vec_clear(x+3,3);
  vec_clear(xxsv+3,3);

  /************************************************
   * transformation to ecliptic.                  *
   * with sidereal calc. this will be overwritten *
   * afterwards.                                  *
   ************************************************/
  if (is_not_set(iflag, SEFLG_EQUATORIAL)) {
    swi_coortrf2(x, x, oe->seps, oe->ceps);
    if (is_set(iflag,SEFLG_SPEED))
      swi_coortrf2(x+3, x+3, oe->seps, oe->ceps);
    if (is_not_set(iflag,SEFLG_NONUT)) {
      swi_coortrf2(x, x, swed->nut.snut, swed->nut.cnut);
      if (is_set(iflag,SEFLG_SPEED))
        swi_coortrf2(x+3, x+3, swed->nut.snut, swed->nut.cnut);
      }
    }
    
  /************************************
   * sidereal positions               *
   ************************************/
  if (is_set(iflag,SEFLG_SIDEREAL)) {
    /* rigorous algorithm */
    if (swed->sidd.sid_mode & SE_SIDBIT_ECL_T0) {
      if (swi_trop_ra2sid_lon(xxsv, x, xxsv, iflag, serr) != OK) return ERR;
      if (is_set(iflag,SEFLG_EQUATORIAL))
        vec_copy(x,xxsv,6);
      } 
    /* project onto solar system equator */
    else if (swed->sidd.sid_mode & SE_SIDBIT_SSY_PLANE) {
      if (swi_trop_ra2sid_lon_sosy(xxsv, x, xxsv, iflag, serr) != OK) return ERR;
      if (is_set(iflag,SEFLG_EQUATORIAL))
        vec_copy(x,xxsv,6);
    /* traditional algorithm */
    } else {
      swi_cartpol_sp(x, x);
      x[0] -= swe_get_ayanamsa(tjd) * DEGTORAD;
      swi_polcart_sp(x, x);
    }
  }
  /************************************************
   * transformation to polar coordinates          *
   ************************************************/
  if (is_not_set(iflag, SEFLG_XYZ))
    swi_cartpol_sp(x, x);
    
  /**********************
   * radians to degrees *
   **********************/
  if (is_not_set(iflag, SEFLG_RADIANS | SEFLG_XYZ)) {
    convert_arcs(x,RADTODEG);
  }

  vec_copy(xx,x,6);
  
  /* if no ephemeris has been specified, do not return chosen ephemeris */
  if (is_not_set(iflgsave, SEFLG_EPHMASK))
    iflag = iflag & ~SEFLG_DEFAULTEPH;
  iflag = iflag & ~SEFLG_SPEED;
  return iflag;
  return_err:
  vec_clear(xx,6);
  return retc;
}

/**********************************************************
 * get fixstar magnitude
 * parameters:
 * star   name of star or line number in star file
 *    (start from 1, don't count comment).
 *        If no error occurs, the name of the star is returned
 *          in the format trad_name, nomeclat_name
 *
 * mag     pointer to a double, for star magnitude
 * serr    error return string
**********************************************************/
int FAR PASCAL_CONV swed_fixstar_mag(char *star, double *mag, char *serr, struct swe_data *swed)
{
  int i;
  int star_nr = 0;
  AS_BOOL  isnomclat = FALSE;
  size_t cmplen;
  char *cpos[20];
  char sstar[SE_MAX_STNAME + 1];
  char fstar[SE_MAX_STNAME + 1];
  char s[AS_MAXCH + 20], *sp, *sp2;  /* 20 byte for SE_STARFILE */
  int line = 0;
  int fline = 0;
  int retc;
  if (serr != NULL)
    *serr = '\0';
  /******************************************************
   * Star file
   * close to the beginning, a few stars selected by Astrodienst.
   * These can be accessed by giving their number instead of a name.
   * All other stars can be accessed by name.
   * Comment lines start with # and are ignored.
   ******************************************************/
  if (swed->fixfp == NULL) {
    if ((swed->fixfp = swid_fopen(SEI_FILE_FIXSTAR, SE_STARFILE, swed->ephepath, serr,swed)) == NULL) {
      retc = ERR;
      goto return_err;
    }
  }
  rewind(swed->fixfp);
  strncpy(sstar, star, SE_MAX_STNAME);
  sstar[SE_MAX_STNAME] = '\0';
  if (*sstar == ',') {
    isnomclat = TRUE;
  } else if (isdigit((int) *sstar)) {
    star_nr = atoi(sstar);
  } else {
    /* traditional name of star to lower case */
    for (sp = sstar; *sp != '\0'; sp++)
      *sp = tolower((int) *sp);
    if ((sp = strchr(sstar, ',')) != NULL)
      *sp = '\0';
  }
  swi_right_trim(sstar);
  cmplen = strlen(sstar);
  if (cmplen == 0) {
    if (serr != NULL)
      sprintf(serr, "swe_fixstar_mag(): star name empty");
    retc = ERR;
    goto return_err;
  }
  while (fgets(s, AS_MAXCH, swed->fixfp) != NULL) {
    fline++;
    if (*s == '#') continue;
    line++;
    if (star_nr == line)
      goto found;
    else if (star_nr > 0)
      continue;
    if ((sp = strchr(s, ',')) == NULL) {
      if (serr != NULL)
  sprintf(serr, "star file %s damaged at line %d", SE_STARFILE, fline);
      retc = ERR;
      goto return_err;
    }
    if (isnomclat) {
      if (strncmp(sp, sstar, cmplen) == 0)
        goto found;
      else
        continue;
    }
    *sp = '\0';  /* cut off first field */
    strncpy(fstar, s, SE_MAX_STNAME);
    *sp = ',';
    fstar[SE_MAX_STNAME] = '\0';  /* force termination */
    swi_right_trim(fstar);
    i = strlen(fstar);
    if (i < (int) cmplen)
      continue;
    for (sp2 = fstar; *sp2 != '\0'; sp2++) {
      *sp2 = tolower((int) *sp2);
    }
    if (strncmp(fstar, sstar, cmplen) == 0)
      goto found;
  }
  if (serr != NULL && strlen(star) < AS_MAXCH - 20)
    sprintf(serr, "star %s not found", star);
  retc = ERR;
  goto return_err;
  found:
  i = swi_cutstr(s, ",", cpos, 20);
  swi_right_trim(cpos[0]);
  swi_right_trim(cpos[1]);
  if (i < 13) {
    if (serr != NULL)
      sprintf(serr, "data of star '%s,%s' incomplete", cpos[0], cpos[1]);
    retc = ERR;
    goto return_err;
  }
  *mag = atof(cpos[13]);
  /* return trad. name, nomeclature name */
  if (strlen(cpos[0]) > SE_MAX_STNAME)
    cpos[0][SE_MAX_STNAME] = '\0';
  if (strlen(cpos[1]) > SE_MAX_STNAME-1)
    cpos[1][SE_MAX_STNAME-1] = '\0';
  sprintf(star, "%s,%s", cpos[0], cpos[1]);
  return OK;
  return_err:
  *mag = 0;
  return retc;
}

const char *FAR PASCAL_CONV swed_get_planet_name(int ipl, char *s, struct swe_data *swed)
{
  int i;
  int retc;
  double xp[6];
  /* function calls for Pluto with asteroid number 134340
   * are treated as calls for Pluto as main body SE_PLUTO */
  if (ipl == SE_AST_OFFSET + 134340)
    ipl = SE_PLUTO;
  if (ipl != 0 && ipl == swed->i_saved_planet_name) {
    strcpy(s, swed->saved_planet_name);
    return s;
  }
  switch(ipl) {
    case SE_SUN:
      strcpy(s, SE_NAME_SUN);
      break;
    case SE_MOON:
      strcpy(s, SE_NAME_MOON);
      break;
    case SE_MERCURY:
      strcpy(s, SE_NAME_MERCURY);
      break;
    case SE_VENUS:
      strcpy(s, SE_NAME_VENUS);
      break;
    case SE_MARS:
      strcpy(s, SE_NAME_MARS);
      break;
    case SE_JUPITER:
      strcpy(s, SE_NAME_JUPITER);
      break;
    case SE_SATURN:
      strcpy(s, SE_NAME_SATURN);
      break;
    case SE_URANUS:
      strcpy(s, SE_NAME_URANUS);
      break;
    case SE_NEPTUNE:
      strcpy(s, SE_NAME_NEPTUNE);
      break;
    case SE_PLUTO:
      strcpy(s, SE_NAME_PLUTO);
      break;
    case SE_MEAN_NODE:
      strcpy(s, SE_NAME_MEAN_NODE);
      break;
    case SE_TRUE_NODE:
      strcpy(s, SE_NAME_TRUE_NODE);
      break;
    case SE_MEAN_APOG:
      strcpy(s, SE_NAME_MEAN_APOG);
      break;
    case SE_OSCU_APOG:
      strcpy(s, SE_NAME_OSCU_APOG);
      break;
    case SE_INTP_APOG:
      strcpy(s, SE_NAME_INTP_APOG);
      break;
    case SE_INTP_PERG:
      strcpy(s, SE_NAME_INTP_PERG);
      break;
    case SE_EARTH:
      strcpy(s, SE_NAME_EARTH);
      break;
    case SE_CHIRON:
    case SE_AST_OFFSET + MPC_CHIRON:
      strcpy(s, SE_NAME_CHIRON);
      break;
    case SE_PHOLUS:
    case SE_AST_OFFSET + MPC_PHOLUS:
      strcpy(s, SE_NAME_PHOLUS);
      break;
    case SE_CERES:
    case SE_AST_OFFSET + MPC_CERES:
      strcpy(s, SE_NAME_CERES);
      break;
    case SE_PALLAS:
    case SE_AST_OFFSET + MPC_PALLAS:
      strcpy(s, SE_NAME_PALLAS);
      break;
    case SE_JUNO:
    case SE_AST_OFFSET + MPC_JUNO:
      strcpy(s, SE_NAME_JUNO);
      break;
    case SE_VESTA:
    case SE_AST_OFFSET + MPC_VESTA:
      strcpy(s, SE_NAME_VESTA);
      break;
    default:
      /* fictitious planets */
      if (ipl >= SE_FICT_OFFSET && ipl <= SE_FICT_MAX) {
        swi_get_fict_name(ipl - SE_FICT_OFFSET, s);
        break;
      }
      /* asteroids */
      if (ipl > SE_AST_OFFSET) {
  /* if name is already available */
  if (ipl == swed->fidat[SEI_FILE_ANY_AST].ipl[0])
    strcpy(s, swed->fidat[SEI_FILE_ANY_AST].astnam);
        /* else try to get it from ephemeris file */
  else {
    retc = sweph(J2000, ipl, SEI_FILE_ANY_AST, 0, NULL, NO_SAVE, xp, NULL,swed);
    if (retc != ERR && retc != NOT_AVAILABLE)
      strcpy(s, swed->fidat[SEI_FILE_ANY_AST].astnam);
    else
      sprintf(s, "%d: not found", ipl - SE_AST_OFFSET);
  }
        /* If there is a provisional designation only in ephemeris file,
         * we look for a name in seasnam.txt, which can be updated by
         * the user.
         * Some old ephemeris files return a '?' in the first position.
         * There are still a couple of unnamed bodies that got their
         * provisional designation before 1925, when the current method
         * of provisional designations was introduced. They have an 'A'
         * as the first character, e.g. A924 RC.
         * The file seasnam.txt may contain comments starting with '#'.
         * There must be at least two columns:
         * 1. asteroid catalog number
         * 2. asteroid name
         * The asteroid number may or may not be in brackets
         */
        if (s[0] == '?' || isdigit((int) s[1])) {
          int ipli = (int) (ipl - SE_AST_OFFSET), iplf = 0;
          FILE *fp;
          char si[AS_MAXCH], *sp, *sp2;
          if ((fp = swid_fopen(-1, SE_ASTNAMFILE, swed->ephepath, NULL, swed)) != NULL) {
            while(ipli != iplf && (sp = fgets(si, AS_MAXCH, fp)) != NULL) {
              while (*sp == ' ' || *sp == '\t'
                     || *sp == '(' || *sp == '[' || *sp == '{')
                sp++;
              if (*sp == '#' || *sp == '\r' || *sp == '\n' || *sp == '\0')
                continue;
              /* catalog number of body of current line */
              iplf = atoi(sp);
              if (ipli != iplf)
                continue;
              /* set pointer after catalog number */
              sp = strpbrk(sp, " \t");
              if (sp == NULL)
                continue; /* there is no name */
              while (*sp == ' ' || *sp == '\t')
                sp++;
              sp2 = strpbrk(sp, "#\r\n");
              if (sp2 != NULL)
                *sp2 = '\0';
              if (*sp == '\0')
                continue;
              swi_right_trim(sp);
              strcpy(s, sp);
            }
            fclose(fp);
          }
        }
      } else  {
  i = ipl;
  sprintf(s, "%d", i);
      }
      break;
  }
  if (strlen(s) < 80) {
    swed->i_saved_planet_name = ipl;
    strcpy(swed->saved_planet_name, s);
  }
  return s;
}

/* set geographic position and altitude of observer */
void FAR PASCAL_CONV swed_set_topo(double geolon, double geolat, double geoalt, struct swe_data *swed)
{
  swed->topd.geolon = geolon;
  swed->topd.geolat = geolat;
  swed->topd.geoalt = geoalt;
  swed->geopos_is_set = TRUE;
  /* to force new calculation of observer position vector */
  swed->topd.teval = 0;
  /* to force new calculation of light-time etc.
   */
  swi_force_app_pos_etc(swed);
}

void swi_force_app_pos_etc(struct swe_data *swed)
{
  for (int i = 0; i < SEI_NPLANETS; i++) swed->pldat[i].xflgs = -1;
  for (int i = 0; i < SEI_NNODE_ETC; i++) swed->nddat[i].xflgs = -1;
  for (int i = 0; i < SE_NPLANETS; i++) {
    swed->savedat[i].tsave = 0;
    swed->savedat[i].iflgsave = -1;
  }
}

int swi_get_observer(double tjd, int iflag,
  AS_BOOL do_save, double *xobs, char *serr)
{
  return swid_get_observer(tjd,iflag,do_save,xobs,serr,&swed);
}

int swid_get_observer(double tjd, int iflag,
  AS_BOOL do_save, double *xobs, char *serr, struct swe_data *swed)
{
  double sidt, delt, tjd_ut, eps, nut, nutlo[2];
  double f = EARTH_OBLATENESS;
  double re = EARTH_RADIUS;
  double cosfi, sinfi, cc, ss, cosl, sinl, h;
  if (!swed->geopos_is_set) {
    throw(ERR,"geographic position has not been set",serr,swed);
  }
  /* geocentric position of observer depends on sidereal time,
   * which depends on UT.
   * compute UT from ET. this UT will be slightly different
   * from the user's UT, but this difference is extremely small.
   */
  delt = swe_deltat(tjd);
  tjd_ut = tjd - delt;
  if (swed->oec.teps == tjd && swed->nut.tnut == tjd) {
    eps = swed->oec.eps;
    nutlo[1] = swed->nut.nutlo[1];
    nutlo[0] = swed->nut.nutlo[0];
  } else {
    eps = swi_epsiln(tjd);
    if (is_not_set(iflag, SEFLG_NONUT))
      swi_nutation(tjd, nutlo);
  }
  if (is_set(iflag, SEFLG_NONUT)) {
    nut = 0;
  } else {
    eps += nutlo[1];
    nut = nutlo[0];
  }
  /* mean or apparent sidereal time, depending on whether or
   * not SEFLG_NONUT is set */
  sidt = swe_sidtime0(tjd_ut, eps, nut);
  sidt *= 15;  /* in degrees */
  /* length of position and speed vectors;
   * the height above sea level must be taken into account.
   * with the moon, an altitude of 3000 m makes a difference
   * of about 2 arc seconds.
   * height is referred to the average sea level. however,
   * the spheroid (geoid), which is defined by the average
   * sea level (or rather by all points of same gravitational
   * potential), is of irregular shape and cannot easily
   * be taken into account. therefore, we refer height to
   * the surface of the ellipsoid. the resulting error
   * is below 500 m, i.e. 0.2 - 0.3 arc seconds with the moon.
   */
  cosfi = cos(swed->topd.geolat * DEGTORAD);
  sinfi = sin(swed->topd.geolat * DEGTORAD);
  cc= 1 / sqrt(cosfi * cosfi + (1-f) * (1-f) * sinfi * sinfi);
  ss= (1-f) * (1-f) * cc;
  /* neglect polar motion (displacement of a few meters), as long as
   * we use the earth ellipsoid */
  /* ... */
  /* add sidereal time */
  cosl = cos((swed->topd.geolon + sidt) * DEGTORAD);
  sinl = sin((swed->topd.geolon + sidt) * DEGTORAD);
  h = swed->topd.geoalt;
  xobs[0] = (re * cc + h) * cosfi * cosl;
  xobs[1] = (re * cc + h) * cosfi * sinl;
  xobs[2] = (re * ss + h) * sinfi;
  /* polar coordinates */
  swi_cartpol(xobs, xobs);
  /* speed */
  xobs[3] = EARTH_ROT_SPEED;
  xobs[4] = xobs[5] = 0;
  swi_polcart_sp(xobs, xobs);
  /* to AUNIT */
  sdiv(xobs,AUNIT,6);
  /* subtract nutation, set backward flag */
  if (is_not_set(iflag,SEFLG_NONUT)) {
    swi_coortrf2(xobs, xobs, -swed->nut.snut, swed->nut.cnut);
    if (is_set(iflag,SEFLG_SPEED))
      swi_coortrf2(xobs+3, xobs+3, -swed->nut.snut, swed->nut.cnut);
    swi_nutate(xobs, iflag, TRUE);
  }
  /* precess to J2000 */
  swi_precess(xobs, tjd, J_TO_J2000);
  if (is_set(iflag, SEFLG_SPEED))
     swi_precess_speed(xobs, tjd, J_TO_J2000);
  /* neglect frame bias (displacement of 45cm) */
  /* ... */
  /* save */
  if (do_save) {
    vec_copy(swed->topd.xobs,xobs,6);
    swed->topd.teval = tjd;
    swed->topd.tjd_ut = tjd_ut;  /* -> save area */
  }
  return OK;
}

/* Equation of Time
 *
 * The function returns the difference between
 * local apparent and local mean time in days.
 * E = LAT - LMT
 * Input variable tjd is ET.
 * Algorithm according to Meeus, German, p. 190ff.
 */
int FAR PASCAL_CONV swe_time_equ(double tjd, double *E, char *serr)
{
  double L0, dpsi, eps, x[6], nutlo[2];
  double tau = (tjd - J2000) / 365250;
  L0 = 280.4664567
       + swe_degnorm(tau * 360007.6982779)
       + tau *
           tau  * (  0.03032028 +
           tau  * (  1 / 49931 -
           tau  * (  1 / 15299 -
           tau  *    1 / 1988000 ) ) );
  swi_nutation(tjd, nutlo);
  eps = (swi_epsiln(tjd) + nutlo[1]) * RADTODEG;
  dpsi = nutlo[0] * RADTODEG;
  if (swe_calc(tjd, SE_SUN, SEFLG_EQUATORIAL, x, serr) == ERR)
    return ERR;
  *E = swe_degnorm(L0 - 0.0057183 - x[0] + dpsi * cos(eps * DEGTORAD));
  if (*E > 180)
    *E -= 360;
  *E *= 4 / 1440.0;
  return OK;
}
