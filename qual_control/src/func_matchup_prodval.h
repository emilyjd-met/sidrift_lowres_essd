#ifndef FUNC_MATCHUP_PRODVAL_H
#define FUNC_MATCHUP_PRODVAL_H

/* some definitions and macros */
#define USE3DCOL   0   /* special value for triggering the use of 3D collocation */
#define USE2DCOL   333 /* special value for triggering the use of 2D collocation */

#define NNRADIUS   40.0 /* radius for searching for nearest neighbour [km] */
#define BLWEIGHT   0.99 /* minimum interpolation weight for accepting an interpolation */
#define MAXTDIFF   3    /* temporal collocation at T0 are acccepted between ]-MAXTDIFF:+MAXTDIFF[ (hours) */
#define MAXDDIFF   1    /* value of Delta = (T1-D0) are acccepted between ]-MAXDDIFF:+MAXDDIFF[ (hours) */

/* object to hold a drift matchup */
typedef struct {
   /* validation data */
   char id[STATION_NAME_LENGTH+1];
   char network[STATION_NETW_LENGTH+1];
   char source[STATION_SRC_LENGTH+1];
   buoyRecord val[2];
   /* prod data */
   float latB_prod;
   float lonB_prod;
   float latE_prod;
   float lonE_prod;
   fmsec1970 tB_prod;
   fmsec1970 tE_prod;
   size_t i_prod;
   size_t j_prod;
   float sX_prod;
   float sY_prod;
   float cXY_prod;
   short flag_prod;
   /* collocation information */
   float distance;
   float dtB;
   float dtE;
} matchProdVal;

void printMatchProdVal(matchProdVal *);
int selectMatchups(dbl_list *allMatchups);

int matchProdWithVal(
      /* inputs */
      int temporal_collocation_scheme,
      int use_nearest_neighbour,
      fmsec1970 refT0, fmsec1970 refT1, /* 'reference' time are only used with the -M flag on and are seconds since 1970-01-01 */
      size_t prodDims[],float A[2], float B[2],PJ *proj,
      float *lat_start_fprod, float *lon_start_fprod, 
      fmsec1970 *t0_fprod, fmsec1970 *t1_fprod, 
      float *driftX_fprod, float *driftY_fprod, short *flags_fprod, float fillvalf_fprod,
      float *sigdX_fprod, float *sigdY_fprod, float *corrdXdY_fprod,
      size_t nbStations, buoyTrajectory trajs[],
      /* outputs */
      dbl_list *allMatchups);

#endif /* FUNC_MATCHUP_PRODVAL_H */
