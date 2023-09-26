#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <projects.h>
#include <fmutil.h>
#include <errorcodes.h>
#include <ice_common.h>
#include <useproj.h>
#include "icedrift_common.h"
#include "icedrift_flags.h"
#include "icedrift_model.h"
#include "icedrift_prepost.h"
#include "icedrift_solve_core.h"

// Global variables
char *reportFile;
char *out_area;
double pattern_radius[NBPATTERNS];  /* [km]    */
double max_icevelocity; /* [m/sec] */
char *out_projstr;
double out_Ax, out_Bx, out_Ay, out_By;
char  *img_projstr;
short twghtStart,twghtEnd;

/* Output fields (globals from Python). */
float *driftX;
float *driftY;
short *pflag;
float *sigdX;
float *sigdY;
float *corrdXdY;
short *uflag;
float *length;
float *dir;
float *outlonB;
float *outlatB;
float *outlonE;
float *outlatE;
float *outfc;
short *snavg;
float *avgX;
float *avgY;
float *length_avg;
float *length_diff;
float *stdX;
float *stdY;
short *patternIndex;
 
/* local variables */
static PJ *out_pj;

  
static void setaPrioriUnknowns(double *sxpr, double **cixpr);
static void setStartPoint(double x[]);

int setPattern(double pattern_radius, short nbWaveBands, PJ *pj, double Ax,
	       double Bx, double Ay, double By,size_t img_dims[3],
	       double *img_lon, double *img_lat,
	       size_t *lpattern_size, size_t *lpattern_center,
	       short **lpattern_mask, long *lpattern_windexes[3],
	       double ***lpattern_img, short ***lpattern_isvalid,
	       double ***lpattern_img2, short ***lpattern_isvalid2, 
	       short **lpattern_icemask,short **lpattern_icemask2);

void set_constraint_parameters(double mult, double pow) {
   constraint_multiplier = mult;
   constraint_power      = pow;
}

void set_sigmoidlength(double len) {
   if (len < 0)
      sigmoid_length = maxdriftdistance;
   else
      sigmoid_length = len;
}

void setBestKnowledge(double lat, double lon) {
   bestKnowledge_lat = lat;
   bestKnowledge_lon = lon;
}


void initmod(int *n, double x[]) {

   int np,nu;
   unsigned int tote = NDRIFTPIXELS * NUNKNOWNS; /* should be '*n' but can we trust n?*/
   unsigned int e;
   int ret;
   
   double timediff;

   /* allocate memory */
   /* *************** */
   /* a-priori pdf */
   xpr   = fmMalloc(tote * sizeof(double));
   xdif  = fmMalloc(tote * sizeof(double));
   sxpr  = fmMalloc(tote * sizeof(*sxpr));

   /* read and set values */
   /* ******************* */
   /* set constraint parameters for evaluating the cost function */
   set_constraint_parameters(1.e6,4);

   /* intialize the out projection object (PROJ4)  */
   out_pj = NULL;
   out_pj = pj_init_plus(out_projstr);
   if (!out_pj) {
      printf("ERROR (initmod) could not initialize the PROJ4 output remapping object.\n");
      return 1;
   }

   /* intialize the img projection object (PROJ4) */
   img_pj = NULL;
   img_pj = pj_init_plus(img_projstr);
   if (!img_pj) {
      printf("ERROR (initmod) could not initialize the PROJ4 image file remapping object.\n");
      return 1;
   }
   
   /* precompute the pattern shape */
   for (size_t pat = 0 ; pat < NBPATTERNS ; pat++) {
      ret = setPattern(pattern_radius[pat], nbWaveBands,img_pj, img_Ax,
		       img_Bx, img_Ay, img_By, img_dims, img_lon, img_lat,
		       &(pattern_size[pat]), &(pattern_center[pat]),
		       &(pattern_mask[pat]), pattern_windexes[pat],
		       &(pattern_img[pat]), &(pattern_isvalid[pat]),
		       &(pattern_img2[pat]), &(pattern_isvalid2[pat]),
		       &(pattern_icemask[pat]), &(pattern_icemask2[pat]));
      if (ret) {
	fprintf(stderr, "ERROR (initmod) while setting the %u pattern shape.\n", pat);
	exit(EXIT_FAILURE);
      }
   }

   /* a-priori PDF */
   setaPrioriUnknowns(sxpr,cixpr);

   /* 'eval' or 'starting point' */
   setStartPoint(x);

   /* dump information to stdout */
   /* ************************** */
   for (short pattern = 0 ; pattern < NBPATTERNS ; pattern++) {
      printf("\tRadius for pattern #%1d is %4.1fkm\n", pattern,
	     pattern_radius[pattern]);
   }
   printf("\tMaximum drift distance is %5.2fkm\n",maxdriftdistance);

   FFLUSH;

}

int prepare_icedriftProduct(int *n, double x[], double fc[],
			    short processingflag[], short pattern[],
			    float leng[], float dire[], float latB[],
			    float lonB[], float latE[], float lonE[],
			    size_t navg[], float xavg[], float yavg[],
			    float lenavg[], float lendiff[], float xstd[],
			    float ystd[], double x_stddev[],
			    double y_stddev[], double xy_correl[],
			    short uncertaintyflag[]) {

   /* Local vars */ 
   int ret;
   fmsec1970 *tstart;
   fmsec1970 *tend;

   float fillvalf = (float) UNDEFNC_FLOAT;
   int   fillvali = (int)   UNDEFNC_INT;
   short fillvals = (short) UNDEFNC_SHORT;

   /* Allocating temporary arrays */
   tstart = fmMalloc(out_dims[TDIM]*sizeof(*tstart));
   tend   = fmMalloc(out_dims[TDIM]*sizeof(*tend));

   /* Currently filling arrays with fill values in C as this is probably more efficient than Python.*/
   for (size_t e = 0 ; e < out_dims[TDIM] ; e++) {
      driftX[e] = driftY[e] = outfc[e] = fillvalf;
      length[e] = dir[e] = fillvalf;
      outlonB[e] = outlatB[e] = fillvalf;
      outlonE[e] = outlatE[e] = fillvalf;
      avgX[e]    = avgY[e]    = fillvalf;
      stdX[e]    = stdY[e]    = fillvalf;
      length_avg[e] = length_diff[e] = fillvalf;
      tstart[e] = tend[e] = -1;
      patternIndex[e] = fillvals;
      pflag[e] = fillvals;
      snavg[e] = fillvals;
   }

   if (uflag) {
      for (size_t e = 0 ; e < out_dims[TDIM] ; e++) {
         sigdX[e]  = sigdY[e] = corrdXdY[e] = fillvalf;
         uflag[e]  = fillvals;
      }
   }

   /* fill the fields with the estimated drift vectors */
   for (size_t e = 0 ; e < NDRIFTPIXELS ; e++) {

      unsigned int out_write = owcs[e];

      pflag[out_write]        = processingflag[e];
      if (uncertaintyflag) {
         uflag[out_write]        = uncertaintyflag[e];
      }
      if ( (processingflag[e] == ICEDRIFT_OK) || 
	    (processingflag[e] == ICEDRIFT_CORRECT_BY_NEIGHBOURS) ) {

	 driftX[out_write] = (float)x[e*NUNKNOWNS+XDRIFT_IDX];
	 driftY[out_write] = (float)x[e*NUNKNOWNS+YDRIFT_IDX];
	 length[out_write] = (float)leng[e];
	 dir[out_write]    = (float)dire[e];
	 outlatE[out_write]= latE[e];
	 outlonE[out_write]= lonE[e];
	 outlatB[out_write]= latB[e];
	 outlonB[out_write]= lonB[e];
	 outfc[out_write]  = (float)fc[e];
	 patternIndex[out_write] = pattern[e];
     if (pattern[e]>0) {
        pflag[out_write]    = ICEDRIFT_SMALLER_PATTERN;
     }
	 if ( ! ((xavg[e] == 1000) && (yavg[e] == 1000)) ) {
	    avgX[out_write]        = xavg[e];
	    avgY[out_write]        = yavg[e];
	    length_avg[out_write]  = lenavg[e];
	    length_diff[out_write] = lendiff[e];
	 }
	 if ( ! ((xstd[e] == 1000) && (ystd[e] == 1000)) ) {
	    stdX[out_write]        = xstd[e];
	    stdY[out_write]        = ystd[e];
	 }
     snavg[out_write] = navg[e];
     
     if (uncertaintyflag && (uncertaintyflag[e] == ICEDRIFTPOST_OK)) {
        sigdX[out_write] = x_stddev[e];
        sigdY[out_write] = y_stddev[e];
        corrdXdY[out_write] = xy_correl[e];
     }

      }
   }
   
   /* Allow to return, the arrays needed for writing out in Python are written. */

   return ret;
}


void setaPrioriUnknowns(double *sxpr, double **cixpr) {

   /* for now we use hard-coded values */
   int np;
   unsigned int e;

   for ( np = 0 ; np < NDRIFTPIXELS ; np++ ) {
      e = np*NUNKNOWNS;
      xpr[e+XDRIFT_IDX] = 0.; sxpr[e+XDRIFT_IDX] = 3*maxdriftdistance; 
      xpr[e+YDRIFT_IDX] = 0.; sxpr[e+YDRIFT_IDX] = 3*maxdriftdistance; 
   }
}

void setStartPoint(double x[]) {

   /* copy the xpr[] values into x[] */
   unsigned int e;
   for ( e = 0 ; e < NDRIFTPIXELS*NUNKNOWNS ; e++ ) {
      x[e] = xpr[e];
   }

}

int setPattern(double pattern_radius, short nbWaveBands, PJ *pj, double Ax,
	       double Bx, double Ay, double By,size_t img_dims[3],
	       double *img_lon, double *img_lat, size_t *lpattern_size,
	       size_t *lpattern_center, short **lpattern_mask,
	       long *lpattern_windexes[3], double ***lpattern_img,
	       short ***lpattern_isvalid, double ***lpattern_img2,
	       short ***lpattern_isvalid2, short **lpattern_icemask,
	       short **lpattern_icemask2) {

   /* the 'correlation' pattern is pre-computed as a disc
    * of radius pattern_radius. */
   int ret;
   int pass;
   
   double latcentral,loncentral;
   if (strstr(out_area,"nh")) {
      latcentral = 90.;
      loncentral = 0.;
   } else if (strstr(out_area,"sh")) {
      latcentral = -90;
      loncentral = 0;
   } else {
      fprintf(stderr,"ERROR (%s) Unknown grid when setting the block-pattern shape.\n",__func__);
      return 1;
   }

   double xPOL,yPOL;
   ret = remap_ll2xy(latcentral,loncentral,pj,Ax,Bx,Ay,By,&xPOL,&yPOL,1);

   int xcentralPixel = (int)xPOL;
   int ycentralPixel = (int)yPOL;

   long centralElem   = fmivec(xcentralPixel,ycentralPixel,img_dims[XDIM]);
   loncentral = img_lon[centralElem];
   latcentral = img_lat[centralElem];
  
   /* with the image resolution we can design a rectangle which encompasses the wished disc */
   float minA = (Ax<Ay?Ax:Ay);
   int largeRegionRadius = (int)ceil(pattern_radius/minA);
   if (largeRegionRadius < 1) largeRegionRadius = 1;
   
   unsigned int elem;
   for ( pass = 0 ; pass <= 1 ; pass ++ ) {
      if (pass == 1) {
         *lpattern_size = elem;
         lpattern_windexes[XDIM] = fmMalloc(*lpattern_size * sizeof(long));
         lpattern_windexes[YDIM] = fmMalloc(*lpattern_size * sizeof(long));
         lpattern_windexes[TDIM] = fmMalloc(*lpattern_size * sizeof(long));
         *lpattern_mask     = fmMalloc(*lpattern_size * sizeof(short));
      }
 
      elem = 0;
      for ( int x = -largeRegionRadius ; x <= largeRegionRadius ; x ++ ) {
         //if (pass == 1) printf("\t\t");
         for ( int y = -largeRegionRadius ; y <= largeRegionRadius ; y ++ ) {
            long index = fmivec(xcentralPixel+x,ycentralPixel+y,img_dims[XDIM]);
            double d;
            
            compute_distance(img_lat[index],img_lon[index],latcentral,loncentral,&d);
	      
	   
	        if ( (x == 0) && (y == 0) ) {
               *lpattern_center = elem;
	           //printf("Center pixel is %u (%fkm)\n",*lpattern_center,d);
	        }

	        if ( d <= pattern_radius ) {
               if (pass == 1) {
                  (*lpattern_mask)[elem] = 1;
                  lpattern_windexes[XDIM][elem] = x;
                  lpattern_windexes[YDIM][elem] = y;
                  lpattern_windexes[TDIM][elem] = index - centralElem;
		 
                  /*
                     if ( x==0 && y==0 ) {
		                printf("o ");
		             } else {
		                printf("x ");
		             }
		          */
		          //printf("{%+01ld/%+01ld} ",lpattern_windexes[XDIM][elem],lpattern_windexes[YDIM][elem]);
               }
	           elem++;
            } else {
               ;
               /*
               if (pass == 1) {
		          printf(". ");
		          //printf("._____. ");
	           }
	           */
	        }
         }
 	     //if ( pass == 1 ) printf("\n");
 	     //if ( pass == 1 ) printf("\n\n\n\n");
      }
   }


   /* allocate the other pattern_ memory */
   *lpattern_icemask  = fmMalloc(*lpattern_size * sizeof(short));
   *lpattern_icemask2 = fmMalloc(*lpattern_size * sizeof(short));
   *lpattern_img      = fmMalloc(nbWaveBands * sizeof(double *));
   *lpattern_isvalid  = fmMalloc(nbWaveBands * sizeof(short *));
   *lpattern_img2     = fmMalloc(nbWaveBands * sizeof(double *));
   *lpattern_isvalid2 = fmMalloc(nbWaveBands * sizeof(short *));
   for (short b = 0 ; b < nbWaveBands ; b++) {
      (*lpattern_img)[b]      = fmMalloc(*lpattern_size * sizeof(double));
      (*lpattern_isvalid)[b]  = fmMalloc(*lpattern_size * sizeof(short));
      (*lpattern_img2)[b]     = fmMalloc(*lpattern_size * sizeof(double));
      (*lpattern_isvalid2)[b] = fmMalloc(*lpattern_size * sizeof(short));
   }


   return 0;
}


int setNeighbourhoodPattern(double pattern_radius, size_t *lpattern_size,
			    size_t *lpattern_center, short **lpattern_mask,
			    long *lpattern_windexes[3]) {

   /* the 'neighborhood' pattern is pre-computed as a disc of radius pattern_radius. */
   int ret;
   int pass;

   /* set local variables from global ones */
   PJ    *pj = out_pj;
   double Ax = out_Ax;
   double Ay = out_Ay;
   double Bx = out_Bx;
   double By = out_By;
   size_t *dims = out_dims;
   double *lon  = olon;
   double *lat  = olat;
   
     double latcentral,loncentral;
   if (strstr(out_area,"nh")) {
      latcentral = 90.;
      loncentral = 0.;
   } else if (strstr(out_area,"sh")) {
      latcentral = -90;
      loncentral = 0;
   } else {
      fprintf(stderr,"ERROR (%s) Unknown grid when setting the neighbouring-area shape.\n",__func__);
      return 1;
   }
   double xPOL,yPOL;
   ret = remap_ll2xy(latcentral,loncentral,pj,Ax,Bx,Ay,By,&xPOL,&yPOL,1);
   
   int xcentralPixel = (int)xPOL;
   int ycentralPixel = (int)yPOL;
   
   long centralElem   = fmivec(xcentralPixel,ycentralPixel,dims[XDIM]);
   
   loncentral = lon[centralElem];
   latcentral = lat[centralElem];
  
   /* with the grid resolution we can design a rectangle which encompases the wished disc */
   float minA = (Ax<Ay?Ax:Ay);
   int largeRegionRadius = (int)ceil(pattern_radius/minA);
   if (largeRegionRadius < 1) largeRegionRadius = 1;
   

   unsigned int elem;
   for ( pass = 0 ; pass <= 1 ; pass ++ ) {
      if (pass == 1) {
	*lpattern_size = elem;
 	lpattern_windexes[XDIM] = fmMalloc(*lpattern_size * sizeof(long));
 	lpattern_windexes[YDIM] = fmMalloc(*lpattern_size * sizeof(long));
 	lpattern_windexes[TDIM] = fmMalloc(*lpattern_size * sizeof(long));
 	*lpattern_mask     = fmMalloc(*lpattern_size * sizeof(short));
      }
 
      elem = 0;
      for ( int x = -largeRegionRadius ; x <= largeRegionRadius ; x ++ ) {
 	if (pass == 1) printf(" ");
 	for ( int y = -largeRegionRadius ; y <= largeRegionRadius ; y ++ ) {
 	   long index = fmivec(xcentralPixel+x,ycentralPixel+y,dims[XDIM]);
 	   double d;

 	   compute_distance(lat[index],lon[index],latcentral,loncentral,&d);
	   
	   if ( (x == 0) && (y == 0) ) {
	      *lpattern_center = elem;
	      //printf("Center pixel is %u (%fkm)\n",*lpattern_center,d);
	   }

	   if ( d <= pattern_radius ) {
	      if (pass == 1) {
		 (*lpattern_mask)[elem] = ((x==y)&&(x==0)?0:1); /* exclude the center to only keep the neighbours */
		 lpattern_windexes[XDIM][elem] = x;
		 lpattern_windexes[YDIM][elem] = y;
		 lpattern_windexes[TDIM][elem] = index - centralElem;
		 if ((*lpattern_mask)[elem]) {
		    printf("{%+01ld/%+01ld} ",lpattern_windexes[XDIM][elem],lpattern_windexes[YDIM][elem]);
		 } else {
		    printf(".xxxxx. ");
		 }
	      }
	      elem++;
	   } else {
	      ;
	      if (pass == 1) {
		 printf("._____. ");
	      }
	   }

 	}
 	if ( pass == 1 ) printf("\n\n\n\n");
      }
   }

   printf("\tRadius of the neighbourhood filtering pattern is %5.1fkm ([%ld,%ld])\n",
	 pattern_radius,lpattern_windexes[XDIM][0],lpattern_windexes[YDIM][0]);


   return 0;
}


