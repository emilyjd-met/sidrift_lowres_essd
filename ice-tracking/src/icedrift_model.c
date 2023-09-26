

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <projects.h>
#include "icedrift_common.h"
#include "icedrift_model.h"
#include "icedrift_flags.h"

static size_t gnuplot_image_counter = 0;

/* Number of location where we want to estimate the ice drift.*/
int   NDRIFTPIXELS;
short nbWaveBands;

/* Global variables to communicate between the various routines and the cost function */
double *xpr;
double *xdif;
double *sxpr;
double **cixpr;
unsigned int *owcs;
unsigned int *iwcs;
double   *olon;
double   *olat;
size_t   out_dims[3];
double   bestKnowledge_lat;
double   bestKnowledge_lon;

long     center_coord[3];


PJ         *img_pj;
size_t     img_dims[3]; /* [XDIM,YDIM,TDIM] */
double     *img_lat,*img_lon;
double     img_Ax;
double     img_Bx;
double     img_Ay;
double     img_By;
double     **obs[2];
short      *icelandmask[2];
short      **TCflag[2];
long       timeref[2];

size_t     pattern_center[NBPATTERNS];
size_t     pattern_size[NBPATTERNS];
short      *pattern_mask[NBPATTERNS];
short      *pattern_icemask[NBPATTERNS];
short      *pattern_icemask2[NBPATTERNS];
long       *pattern_windexes[NBPATTERNS][3];
double     **pattern_img[NBPATTERNS];
short      **pattern_isvalid[NBPATTERNS];
double     **pattern_img2[NBPATTERNS];
short      **pattern_isvalid2[NBPATTERNS];

double *lengths;
double *directions;

/* constraint_multiplier: global variable to hold a (huge) multiplier 
 * that "soft constraints" will use. 
 */
double constraint_multiplier; 
/* constraint_power: global variable to hold a polynomial power
 * that "soft constraints" will use. 
 */
double constraint_power;

/* maximum distance [km] for the drift distance */
double maxdriftdistance;

/* length of cut-off for sigmoid step function */
double sigmoid_length;

/* power for the sigmoid step function */
double sigmoid_power = 10;

double pattern_Alpha[MAX_NBWAVEBANDS];
double pattern_Beta[MAX_NBWAVEBANDS];

double distance_sigmoid(double X, double Xlim) {
   double ret = 1. / (1. + exp((X - Xlim)*sigmoid_power));
   //printf("SIGMOID: x: %f dist: %f, pow: %f ====> ret: %f\n",X,Xlim,sigmoid_power,ret);
   return ret;
}

void mytrunc(double x, long *xtrunc) {
   if (x >= 0) {
      *xtrunc = (long)floor(x);
   } else {
      *xtrunc = (long)ceil(x);
   }
}

void locfmivec(long *e, long x, long y, unsigned long nx) {
   *e = y*nx+x;
}

void locfmijmap(long elem, unsigned long nx, long *x, long *y) {
   /* integers ratio with casting : floor(). */
   *y = elem/nx; 
   *x = elem - (*y*nx);
}

/* n is the number of drift vector 4-uple we should test. */
void unknownsAreValid(int n, double *x, int *isvalid, double *penalty, double dlen) {

   long e1;
   double xdrift,ydrift,dist;

   double startPenalty = *penalty;

   /* compute the cost due to "soft constraints" on unknowns validity range */
   /* ********************************************************************* */
   /* implement a maximum distance drift. drifts with longer distance will be penalized */
   for (e1 = 0 ; e1 < n ; e1++) {

      if (dlen > maxdriftdistance) {
         *penalty += constraint_multiplier * pow(dlen - maxdriftdistance,constraint_power);
      }
   }

   if (*penalty > startPenalty)
      *isvalid = 0;
   else
      *isvalid = 1;


}

void load_subimage(int *ret, double **image, short **TCflag, short *icelandmask, size_t img_dims[3],
      long img_coords[3], size_t pattern_size, 
      short *pattern_mask, long *pattern_windexes[3], double **pattern_img, short **pattern_isvalid, short *pattern_icemask) {

   long e,wcorner,xpi,ypi;

   /* x,y coordinates of the central pixel */
   for ( e = 0 ; e < pattern_size ; e++) {

      /*
         for (short b = 0 ; b < nbWaveBands ; b++)
         pattern_img[b][e]     = 0;
         */

      if (pattern_mask[e]) {
         /* by default, we use pixels which are validated by the mask. */
         /* the 'mask' gives the pattern shape (a-priori of any local condition)
          * and the 'isvalid' gives the pattern shape, once all limiting effects
          * of e.g. image border, missing pixels, etc... are taken into account
          */

         /* x,y coordinates of the neighbouring pixels in the pattern */
         xpi = img_coords[XDIM]+pattern_windexes[XDIM][e];
         ypi = img_coords[YDIM]+pattern_windexes[YDIM][e];
         if ( pointIsInGrid((double)(xpi),(double)(ypi),img_dims) ) {
            locfmivec(&wcorner,xpi,ypi,img_dims[XDIM]);
            pattern_icemask[e] = icelandmask[wcorner];
            for (short b = 0 ; b < nbWaveBands ; b++) {
               pattern_img[b][e] = image[b][wcorner];
               pattern_isvalid[b][e] = TCflag[b][wcorner];
            }
         } else {
            //printf("%s -> Found %ld (%ld,%ld) which is outside the image grid\n",__func__,e,xpi,ypi);
            for (short b = 0 ; b < nbWaveBands ; b++) {
               pattern_isvalid[b][e] = TCIMAGE_OUTSIDE_GRID;
            }
         }
      }

   }

   *ret = 0;
}

/* compute_shifted_image: newimage depends on xdrift and ydrift */
void compute_shifted_image(int *ret,double xdrift, double ydrift,
      double **image,short **TCflag, short *icelandmask, size_t img_dims[3],
      long coord[3],
      size_t pattern_size, short *pattern_mask, long *pattern_windexes[3], 
      double **pattern_img, short **pattern_isvalid, short *pattern_icemask,
      double *olat, double *olon) {

   /* compute the value of each pixel in the drifted window */
   long e;
   long xpi,ypi;
   long xr,yr;
   short corner;
   long centralIndexes[4];
   long xcorners[4],ycorners[4],wcorner;
   double weights[4];
   double lweight[MAX_NBWAVEBANDS];  
   short  corner_flags[4], corner_iceflags[4];

   /* if the proposed drift is (dx,dy), we want to load subimage with (-dx,-dy) */
   xdrift *= -1; ydrift *= -1;

   /* x a real number. 
    *    
    *    x = xtrunc + xsign * xeps 
    *
    * with xtrunc = trunc(x)
    *      xsign  = sign(x)
    *      xeps   = abs(x - trunc(x))
    *
    * same for y
    */
   long   xtrunc, ytrunc;
   double xeps, yeps;
   short  xsign, ysign;

   mytrunc(xdrift,&xtrunc); 
   xeps  = fabs(xdrift-xtrunc);
   xsign  = 1;
   if (xdrift < 0) xsign = -1;

   mytrunc(ydrift,&ytrunc); 
   yeps  = fabs(ydrift-ytrunc);
   ysign  = 1;
   if (ydrift < 0) ysign = -1;

   weights[0] = (1. - xeps) * (1. - yeps);
   weights[1] = xeps        * (1. - yeps);
   weights[2] = (1. - xeps) * yeps;
   weights[3] = xeps        * yeps;

   /*
      printf("In load_shifted_image: pattern_size is %u. wc=%u x=%f, y=%f\n",pattern_size,wc,xdrift,ydrift);
      printf("In load_shifted_image: xtrunc is %ld. xeps is %f\n",xtrunc,xeps);
      */

   for ( e = 0 ; e < pattern_size ; e++) {

      if (pattern_mask[e]) {

         /* initialisation */
         for (short b = 0 ; b < nbWaveBands ; b++) {
            lweight[b]            = 0;
            pattern_img[b][e]     = 0; 
         }

         /* compute from which position the current pixel is coming from (integer drift) */
         xr = coord[XDIM] + pattern_windexes[XDIM][e] + xtrunc;
         yr = coord[YDIM] + pattern_windexes[YDIM][e] + ytrunc;

         /* compute the x and y coordinates of the 4 'corner' pixels */
         xcorners[0] = xr;       ycorners[0] = yr;
         xcorners[1] = xr+xsign; ycorners[1] = yr;
         xcorners[2] = xr;       ycorners[2] = yr+ysign;
         xcorners[3] = xr+xsign; ycorners[3] = yr+ysign;

         /* compute the weighted interpolation at current pixel */
         for (corner = 0 ; corner < 4 ; corner ++) {
            if (weights[corner]) {
               //printf("\tpixel %d (%d/4): (%ld,%ld) -> %f\n",e,corner,xcorners[corner],ycorners[corner],weights[corner]);
               /* make sure the corner pixel comes from inside the original image */
               if ( pointIsInGrid((double)(xcorners[corner]),(double)(ycorners[corner]),img_dims) ) {
                  locfmivec(&wcorner,xcorners[corner],ycorners[corner],img_dims[XDIM]);
                  /* only use those pixels that have a valid value */
                  corner_iceflags[corner] = icelandmask[wcorner];
                  for (short b = 0 ; b < nbWaveBands ; b++) {
                     corner_flags[corner] = TCflag[b][wcorner];
                     //if ((corner_iceflags[corner] == 3) && (corner_flags[corner] == TCIMAGE_OK)) {}
                     if ((corner_iceflags[corner] == 3 || corner_iceflags[corner] == 2) \
                           && (corner_flags[corner] == TCIMAGE_OK)) {
                        /* this pixel is on ice AND with a valid observation: we use it */
                        pattern_img[b][e] += weights[corner] * image[b][wcorner];
                        lweight[b] += weights[corner];
                     }
                  }
               } else {
                  corner_flags[corner] = TCIMAGE_OUTSIDE_GRID;
               }
            }
         }
         /* normalize the weighted average. lweight is 1 if all the 4 corners were used */
         /* from the corners flags, we also decide if this drift is valid or not */
         for (short b = 0 ; b < nbWaveBands ; b++) {
            if (lweight[b] > 0.75) {
               pattern_img[b][e] /= lweight[b];
               pattern_isvalid[b][e] = TCIMAGE_OK;
            } else {
               pattern_isvalid[b][e] = TCIMAGE_FAILED;
            }
         }
      }
   }


   /* compute the (lat,lon) of the center of the central point of the pattern, once drifted */
   xdrift += coord[XDIM];
   ydrift += coord[YDIM];
   remap_xy2ll(xdrift,ydrift,img_pj,img_Ax,img_Bx,img_Ay,img_By,olat,olon);

   *ret = 0;
}


void applyForwardModel_onevector(double *x,
      int          NDRIFTPIXELS,short nbWaveBands,
      size_t       pix,
      long         coord[3],
      double **img, short **TCflag, short *icelandmask, size_t img_dims[3], 
      size_t pattern_size, short *pattern_mask, long *pattern_windexes[3], 
      double **pattern_img, short **pattern_isvalid,short *pattern_icemask,
      double *olat, double *olon) {

   int ret;
   long e,eaux,eobs,p,wc;
   double xdrift;
   double ydrift;
   double xdrift_p;
   double ydrift_p;
   double alpha;
   double beta;


   //printf("ApplyForwardModel_onevector::begin.\n");

   e = pix;
   eaux = e*NUNKNOWNS;
   //printf("\te is %ld. eaux is %ld\n",e,eaux);

   /* get the unknowns from the x[] vector */
   xdrift = x[eaux + XDRIFT_IDX];
   ydrift = x[eaux + YDRIFT_IDX];
   /*
      alpha  = x[eaux + ALPHA_IDX];
      beta   = x[eaux + BETA_IDX];
      */

   /* load the drifted (by xdrift,ydrift) pattern */
   xdrift_p = xdrift/img_Ax;
   ydrift_p = ydrift/img_Ay;

   //printf("\tunknowns are {%f,%f,%f,%f} (passed to compute_shifted_image: %f,%f)\n",xdrift,ydrift,alpha,beta,xdrift_p,ydrift_p);
   compute_shifted_image(&ret,xdrift_p,ydrift_p,
         img,TCflag,icelandmask, img_dims,
         coord,
         pattern_size,pattern_mask,pattern_windexes,pattern_img,pattern_isvalid,pattern_icemask,
         olat,olon);


   /* transform the intensity (by alpha,beta) of the pattern */
   /*
      for ( p = 0 ; p < pattern_size ; p ++) {
      for (short b = 0 ; b < nbWaveBands ; b++) {
      if (pattern_isvalid[b][p]) {
      pattern_img[b][p] = alpha * pattern_img[b][p] + beta;
      } else {
      ;//printf("element %ld in the pattern is not valid.\n",p);
      }
      }
      }
      */


   /* put the pattern in place in ocalc[] */
   /*
      eobs = 0;
      for ( p = 0 ; p < pattern_size ; p ++) {
      if (pattern_isvalid[p]) {
      ocalc_vals[p] = pattern_img[p];
      }
      }
      */

   //printf("ApplyForwardModel_onevector::end.\n");

}

void model_onevector_sqdiff_linearScaled(int *n, double *x, double *fc, long e, short p) {

   int ret;

   double rfc    = 0;
   double badfc  = 10000.;

   long pix;
   double *xloc;
   int isvalid;

   long e1,o,wind;
   long xi,yi,xpi,ypi;
   unsigned int nobs,nobs2;

   double olat,olon,dlen;

   /* xloc: access point to the 4 parameters for one drift vector: */
   pix  = e*NUNKNOWNS;
   xloc = &(x[pix]);

   double xdrift_try = xloc[XDRIFT_IDX];
   double ydrift_try = xloc[YDRIFT_IDX];

   /* compute the cost due to the observation/model mismatch */
   /* ****************************************************** */
   /* simulate the start_image when transformed by the current values of the unknown */

   //printf("center_coord is (%ld,%ld,%ld)\n",center_coord[XDIM],center_coord[YDIM],center_coord[TDIM]);

   applyForwardModel_onevector(x,
         NDRIFTPIXELS,nbWaveBands,e,center_coord,
         obs[BEG], TCflag[BEG], icelandmask[BEG],img_dims, 
         pattern_size[p], pattern_mask[p], pattern_windexes[p], pattern_img[p], pattern_isvalid[p],pattern_icemask[p],
         &olat, &olon);


   size_t cptElems[MAX_NBWAVEBANDS];
   size_t cptValid[MAX_NBWAVEBANDS];
   for (short b = 0 ; b < nbWaveBands ; b++) {
      cptElems[b] = 0;
      cptValid[b] = 0;
      //printf("%s: [%d] alpha = %f and beta = %f\n",__func__,b,pattern_Alpha[b],pattern_Beta[b]);
   }

   for (o = 0 ; o < pattern_size[p] ; o++) {
      if ( pattern_mask[p][o] ) {
         for (short b = 0 ; b < nbWaveBands ; b++) {
            cptElems[b]++;
            if ( (pattern_isvalid2[p][b][o] == TCIMAGE_OK) && (pattern_isvalid[p][b][o] == TCIMAGE_OK) ) {
               cptValid[b] ++;
            }
         }
      }
   }
   size_t totValid = 0;
   for (short b = 0 ; b < nbWaveBands ; b++) {
      //printf("\tCPTVALID[%d] is %u (%u)\n",b,cptValid[b],cptElems[b]);
      if (cptValid[b] >= 0.75*cptElems[b]) {
         totValid ++;
      }
   }
   //printf("TOTVALID is %u\n",totValid);

   double obs_fcs[MAX_NBWAVEBANDS];
   if (totValid==0) {
      /* not enough valid pixels in the end image (all bands). just stop now */
      //printf("Not enough valid elements, exit\n");
      *fc = badfc;
   } else {


      /* compare this modified image to the end-of-drift image */
      /* SQUARED DIFFERENCE WITH IMPOSED ALPHA AND BETA */

      double obs_fc = 0;
      for (short b = 0 ; b < nbWaveBands ; b++) {

         double sqdiff = 0;
         size_t nbvals = 0;
         for (o = 0 ; o < pattern_size[p] ; o++) {
            if ( pattern_mask[p][o] && (pattern_isvalid2[p][b][o] == TCIMAGE_OK) && (pattern_isvalid[p][b][o] == TCIMAGE_OK) ) {
               double v1 = pattern_img[p][b][o]; double v2 = pattern_img2[p][b][o];
               if (isnan(v1) || isnan(v2)) continue;
               sqdiff += pow(v2 - (pattern_Alpha[b]*v1 + pattern_Beta[b]),2.);
               nbvals++;
            }
         }
         if (nbvals <= 1) {
            /* no valid pixel in the images. just go to next band */
            continue;
         }

         sqdiff /= nbvals;
         obs_fc += sqdiff;
         obs_fcs[b] = sqdiff;
      }
      obs_fc /= nbWaveBands;
      obs_fc *= 0.5;

      wind = center_coord[TDIM];
      /* compute the length of the difference vector between the 
       * best knowledge so far and the current estimate */
      /* from (olat,olon): compute the drift length */
      compute_distance(bestKnowledge_lat,bestKnowledge_lon,olat,olon,&dlen);

      /* use dlen to modify the correlation factor if dlen is too long. We use a sigmoid
       * function in order to : 
       * a) keep the correlation unchanged if dlen < len_limit, but 
       * b) force it towards -1 otherwise
       */
      double sigmoid_weight = distance_sigmoid(dlen,sigmoid_length);
      //printf("FC BEFORE:  %f\n",obs_fc);
      double dcost=0.;
      if (dlen > sigmoid_length)
         dcost = badfc*pow(dlen-sigmoid_length,4.);

      /*
         if (sigmoid_weight<1)
         printf("FC  DCOST=%f and SIGMOID_WEIGHT=%f\n",dcost,sigmoid_weight);
         */

      if (sigmoid_weight == 0) {
         obs_fc = DBL_MAX;
      } else {
         obs_fc /= sigmoid_weight;
      }

      *fc = obs_fc;

   }

}

void model_onevector_corr(int *n, double *x, double *fc, long e, short p) {

   int ret;

   double corr    = 0;
   double badcorr = -1;

   long pix;
   double *xloc;
   int isvalid;

   long e1,o,wind;
   long xi,yi,xpi,ypi;
   unsigned int nobs,nobs2;

   double olat,olon,dlen;

   /* xloc: access point to the 4 parameters for one drift vector: */
   pix  = e*NUNKNOWNS;
   xloc = &(x[pix]);

   double xdrift_try = xloc[XDRIFT_IDX];
   double ydrift_try = xloc[YDRIFT_IDX];

   /* compute the cost due to the observation/model mismatch */
   /* ****************************************************** */
   /* simulate the start_image when transformed by the current values of the unknown */

   //printf("center_coord is (%ld,%ld,%ld)\n",center_coord[XDIM],center_coord[YDIM],center_coord[TDIM]);

   applyForwardModel_onevector(x,
         NDRIFTPIXELS,nbWaveBands,e,center_coord,
         obs[BEG], TCflag[BEG], icelandmask[BEG],img_dims, 
         pattern_size[p], pattern_mask[p], pattern_windexes[p], pattern_img[p], pattern_isvalid[p],pattern_icemask[p],
         &olat, &olon);


   size_t cptElems[MAX_NBWAVEBANDS];
   size_t cptValid[MAX_NBWAVEBANDS];
   for (short b = 0 ; b < nbWaveBands ; b++) {
      cptElems[b] = 0;
      cptValid[b] = 0;
      pattern_Alpha[b] = -99.;
      pattern_Beta[b]  = -99.;
   }

   for (o = 0 ; o < pattern_size[p] ; o++) {
      if ( pattern_mask[p][o] ) {
         for (short b = 0 ; b < nbWaveBands ; b++) {
            cptElems[b]++;
            if ( (pattern_isvalid2[p][b][o] == TCIMAGE_OK) && (pattern_isvalid[p][b][o] == TCIMAGE_OK) ) {
               cptValid[b] ++;
            }
         }
      }
   }
   size_t totValid = 0;
   for (short b = 0 ; b < nbWaveBands ; b++) {
      //printf("\tCPTVALID[%d] is %u (%u)\n",b,cptValid[b],cptElems[b]);
      if (cptValid[b] >= 0.75*cptElems[b]) {
         totValid ++;
      }
   }
   //printf("TOTVALID is %u\n",totValid);

   double obs_corrs[MAX_NBWAVEBANDS];
   if (totValid==0) {
      /* not enough valid pixels in the end image (all bands). just stop now */
      //printf("Not enough valid elements, exit\n");
      *fc = badcorr;
   } else {


      /* compare this modified image to the end-of-drift image */
      /* CORRELATION */

      /* June 2008 : the one-pass correlation routine seems quite unstable. Implement a
       * two-pass one: */

      double obs_corr = 0;
      for (short b = 0 ; b < nbWaveBands ; b++) {

         double meanv,meanv2;
         double stdv,stdv2;

         meanv = meanv2 = stdv = stdv2 = 0;

         size_t nbvals = 0;
         for (o = 0 ; o < pattern_size[p] ; o++) {
            if ( pattern_mask[p][o] && (pattern_isvalid2[p][b][o] == TCIMAGE_OK) && (pattern_isvalid[p][b][o] == TCIMAGE_OK) ) {
               double v1 = pattern_img[p][b][o]; double v2 = pattern_img2[p][b][o];
               if (isnan(v1) || isnan(v2)) continue;
               meanv  += v1;
               meanv2 += v2;
               stdv   += v1*v1;
               stdv2  += v2*v2;

               if (isnan(meanv)) {
                  printf("NAN (meanv): b is %d, v1 is %f\n",b,v1);
               }
               nbvals++;

            }
         }
         //   printf("NBVALS[%d] is %d\n",b,nbvals);
         if (nbvals <= 1) {
            /* no valid pixel in the images. just go to next band */
            continue;
         }
         meanv  /= nbvals;
         meanv2 /= nbvals;
         stdv   /= nbvals;
         stdv2  /= nbvals;
         stdv   -= meanv*meanv; 
         stdv2  -= meanv2*meanv2;
         stdv    = sqrt(stdv);
         stdv2   = sqrt(stdv2);


         double correlation = 0;
         for (o = 0 ; o < pattern_size[p] ; o++) {
            if ( pattern_mask[p][o] && (pattern_isvalid2[p][b][o] == TCIMAGE_OK) && (pattern_isvalid[p][b][o] == TCIMAGE_OK) ) {
               double v1 = pattern_img[p][b][o]; double v2 = pattern_img2[p][b][o];
               if (isnan(v1) || isnan(v2)) continue;
               correlation += (v1-meanv)*(v2-meanv2);
            }
         }
         correlation /= nbvals;
         pattern_Alpha[b] = correlation / pow(stdv,2.);
         pattern_Beta[b]  = meanv2 - pattern_Alpha[b]*meanv;

         correlation /= stdv * stdv2;
         obs_corr += correlation;
         obs_corrs[b] = correlation;
         // printf("CORR[%d] is %f\n",b,correlation);
         if (isnan(correlation)) {
            printf("NAN! nbvals is %u. %f,%f,%f,%f\n",nbvals,meanv,meanv2,stdv,stdv2);
         }
      }
      obs_corr /= nbWaveBands;

      wind = center_coord[TDIM];
      /* compute the length of the difference vector between the 
       * best knowledge so far and the current estimate */
      /* from (olat,olon): compute the drift length */
      compute_distance(bestKnowledge_lat,bestKnowledge_lon,olat,olon,&dlen);
      //printf("DISTANCE IS %f (%f,%f)->(%f,%f) (%f)\n",dlen,img_lat[wind],img_lon[wind],olat,olon,maxdriftdistance);

      /* use dlen to modify the correlation factor if dlen is too long. We use a sigmoid
       * function in order to : 
       * a) keep the correlation unchanged if dlen < len_limit, but 
       * b) force it towards -1 otherwise
       */
      double sigmoid_weight = distance_sigmoid(dlen,sigmoid_length);
      //printf("Correlation BEFORE: %f (w: %f)\n",obs_corr,sigmoid_weight);
      obs_corr = (obs_corr + 1) * sigmoid_weight - 1;
      //printf("Correlation AFTER:  %f\n",obs_corr);

      *fc = obs_corr;

   }

#undef GNUPLOT_MODEL
#ifdef GNUPLOT_MODEL

   short b = 0;

   /* write the images to ascii file */
   char DataFileName[100];
   sprintf(DataFileName,"/tmp/windows/window_%05d.dat",gnuplot_image_counter);
   FILE *dataFile = fopen(DataFileName,"w");
   long firstX = pattern_windexes[p][XDIM][0];
   long firstY = firstX;
   long xx,yy,xpi2,ypi2,wc;
   long el = 0;
   for ( xx = firstX ; xx <= -firstX ; xx++ ) {
      for ( yy = firstY ; yy <= -firstY ; yy++ ) {
         //printf("xx:%ld yy:%ld windexes[XDIM]=%ld windexes[YDIM]=%ld\n",xx,yy,pattern_windexes[p][XDIM][el],pattern_windexes[p][YDIM][el]);
         if ((el < pattern_size[p]) && (pattern_windexes[p][XDIM][el] == xx) && (pattern_windexes[p][YDIM][el] == yy)) {
            /* x,y coordinates of the neighbouring pixels in the pattern */
            xpi2 = center_coord[XDIM]+pattern_windexes[p][XDIM][el];
            ypi2 = center_coord[YDIM]+pattern_windexes[p][YDIM][el];
            if ( pointIsInGrid((double)(xpi2),(double)(ypi2),img_dims) ) {
               locfmivec(&wc,xpi2,ypi2,img_dims[XDIM]);
               fprintf(dataFile,"%ld %f %f %f %f %f %f %d %d %d %d %d\n",el,xx+0.5,yy+0.5,
                     pattern_img2[p][b][el],pattern_img2[p][b+1][el],
                     pattern_img[p][b][el], pattern_img[p][b+1][el],
                     pattern_mask[p][el],
                     (pattern_isvalid2[p][b][el] == TCIMAGE_OK), (pattern_isvalid2[p][b+1][el]==TCIMAGE_OK),
                     (pattern_isvalid[p][b][el] == TCIMAGE_OK),  (pattern_isvalid[p][b+1][el]==TCIMAGE_OK)
                     );
               //fprintf(dataFile,"%ld %f %f %d\n",el,xx+0.5,yy+0.5,pointIsOnIce(wc,icelandmask[END]));
               el++;
               //if ((el < pattern_size))
               //printf("Now search for %d (%d %d)\n",el,pattern_windexes[XDIM][el],pattern_windexes[YDIM][el]);
            } else {
               //printf("Point %d %d is not in grid!\n",xpi2,ypi2);
            }
         } else {
            //printf("not egal. write 0.\n");
            fprintf(dataFile,"%ld %f %f 0. 0. 0. 0. 0 0 0 0 0\n",el,xx+0.5,yy+0.5);
         }
      }
      fprintf(dataFile,"\n");
   }
   fclose(dataFile);

   /* append the tested point */
   char trajecFileName[32];
   sprintf(trajecFileName,"/tmp/windows/%s","trajec.dat");
   FILE *trajecFile = fopen(trajecFileName,"a");
   fprintf(trajecFile,"%05u %+06.3f %+06.3f 0.5 %f\n",gnuplot_image_counter,xdrift_try,ydrift_try,*fc);
   fclose(trajecFile);

   /* write commands for plotting by gnuplot */
   FILE *cmdFile = fopen("/tmp/commandPlotWind1.gnuplot","w");
   fprintf(cmdFile,"set terminal postscript eps colour\n");
   fprintf(cmdFile,"set output '/tmp/windows/window_%05u.eps'\n",gnuplot_image_counter);
   fprintf(cmdFile,"set size 2.0,1.0 \n");
   /* images */
   fprintf(cmdFile,"set origin 0.0,0.0\n");
   fprintf(cmdFile,"set multiplot\n");
   fprintf(cmdFile,"set xlabel 'X'\n");
   fprintf(cmdFile,"set ylabel 'Y'\n");
   fprintf(cmdFile,"set nokey\n");
   fprintf(cmdFile,"set palette gray\n");
   fprintf(cmdFile,"set xrange [%ld:%ld]\n",firstX,-firstX);
   fprintf(cmdFile,"set yrange [%ld:%ld]\n",firstY,-firstY);
   fprintf(cmdFile,"set cbrange [-3:3]\n");
   fprintf(cmdFile,"set pm3d map corners2color c4\n");
   fprintf(cmdFile,"unset colorbox\n");

   fprintf(cmdFile,"set size   0.5,0.5 \n");
   fprintf(cmdFile,"set origin 0.0,0.0\n");
   if (totValid != 0) {
      fprintf(cmdFile,"set title 'DAYTWO : %+06.4f, B0'\n",obs_corrs[0]);
   } else {
      fprintf(cmdFile,"set title 'DAYTWO : INVALID, B0'\n");
   }
   fprintf(cmdFile,"splot '%s' using 2:3:4\n",DataFileName);
   fprintf(cmdFile,"set origin 0.5,0.0\n");
   if (totValid != 0) {
      fprintf(cmdFile,"set title 'DAYTWO : %+06.4f, B1'\n",obs_corrs[1]);
   } else {
      fprintf(cmdFile,"set title 'DAYTWO : INVALID, B1'\n");
   }
   fprintf(cmdFile,"splot '%s' using 2:3:5\n",DataFileName);
   fprintf(cmdFile,"set origin 0.0,0.5\n");
   fprintf(cmdFile,"set title 'DAYONE(%+05.2f,%+05.2f), B0'\n",xdrift_try,ydrift_try);
   fprintf(cmdFile,"splot '%s' using 2:3:6\n",DataFileName);
   fprintf(cmdFile,"set origin 0.5,0.5\n");
   fprintf(cmdFile,"set title 'DAYONE(%+05.2f,%+05.2f), B1'\n",xdrift_try,ydrift_try);
   fprintf(cmdFile,"splot '%s' using 2:3:7\n",DataFileName);
   /* trajectory */
   fprintf(cmdFile,"set origin 1.0,0.0\n");
   fprintf(cmdFile,"set palette color\n");
   fprintf(cmdFile,"set size   1.0,1.0 \n");
   fprintf(cmdFile,"set xlabel 'dX'\n");
   fprintf(cmdFile,"set ylabel 'dY'\n");
   fprintf(cmdFile,"set nokey\n");
   fprintf(cmdFile,"set xrange [%f:%f]\n",-1.2*maxdriftdistance,1.2*maxdriftdistance);
   fprintf(cmdFile,"set yrange [%f:%f]\n",-1.2*maxdriftdistance,1.2*maxdriftdistance);
   fprintf(cmdFile,"set parametric\n");
   fprintf(cmdFile,"plot [-pi:pi] %f*cos(t),%f*sin(t) with lines lw 3\n",maxdriftdistance,maxdriftdistance);
   fprintf(cmdFile,"unset parametric\n");
   fprintf(cmdFile,"plot '%s' using 2:3:4 with circles fill solid palette frac ($5*0.5)\n",trajecFileName);
   fprintf(cmdFile,"unset multiplot\n");


   fclose(cmdFile);

   /* launch gunplot */
   ret = system("gnuplot /tmp/commandPlotWind1.gnuplot");
   gnuplot_image_counter++;

#endif /* GNUPLOT */

}

#define TMPREMOVE
void model_onevector(int *n, double *x, double *fc, long e) {
   ;
#ifndef TMPREMOVE
   int ret;

   double xcost = 0;
   double ocost = 0;
   double ccost = 0;

   long pix;
   double *xloc;
   int isvalid;
   long cptValid;

   long e1,o,wind;
   long xi,yi,xpi,ypi;
   unsigned int nobs,nobs2;

   double olat,olon,dlen;

   /* xloc: access point to the 4 parameters for one drift vector: */
   pix  = e*NUNKNOWNS;
   xloc = &(x[pix]);


   /* compute the cost due to the a-priori on unknowns */
   /* ************************************************ */
   for ( e1 = 0 ; e1 < NUNKNOWNS ; e1++ ) {
      xdif[e1] = x[pix+e1] - xpr[pix+e1];
   }
   for ( e1 = 0 ; e1 < NUNKNOWNS ; e1 ++) {
      //for ( e2 = 0 ; e2 < tote ; e2 ++) {
      //	 xcost += xdif[e1]*cixpr[e1][e2]*xdif[e2]; 
      //}
      xcost += xdif[e1]*(1./pow(sxpr[pix+e1],2.))*xdif[e1]; 
   }

   /* compute the cost due to the observation/model mismatch */
   /* ****************************************************** */
   /* simulate the start_image when transformed by the current values of the unknown */

#undef GNUPLOT
#ifdef GNUPLOT

   /* write the correlation map to ascii file */
   FILE *dataFile = fopen("/tmp/window1.dat","w");
   long firstX = pattern_windexes[XDIM][0];
   long firstY = firstX;
   long xx,yy,xpi2,ypi2,wc;
   long el = 0;
   for ( xx = firstX ; xx <= -firstX ; xx++ ) {
      for ( yy = firstY ; yy <= -firstY ; yy++ ) {
         if ((el < pattern_size) && (pattern_windexes[XDIM][el] == xx) && (pattern_windexes[YDIM][el] == yy)) {
            /* x,y coordinates of the neighbouring pixels in the pattern */
            xpi2 = center_coord[XDIM]+pattern_windexes[XDIM][el];
            ypi2 = center_coord[YDIM]+pattern_windexes[YDIM][el];
            if ( pointIsInGrid((double)(xpi2),(double)(ypi2),img_dims) ) {
               locfmivec(&wc,xpi2,ypi2,img_dims[XDIM]);
               fprintf(dataFile,"%d %f %f %f %f\n",el,xx+0.5,yy+0.5,obs[END][wc],pattern_img2[el]);
               //fprintf(dataFile,"%d %f %f %d\n",el,xx+0.5,yy+0.5,pointIsOnIce(wc,icelandmask[END]));
               el++;
               //if ((el < pattern_size))
               //printf("Now search for %d (%d %d)\n",el,pattern_windexes[XDIM][el],pattern_windexes[YDIM][el]);
            } else {
               //printf("Point %d %d is not in grid!\n",xpi2,ypi2);
            }
         } else {
            //printf("not egal. write 0.\n");
            fprintf(dataFile,"%d %f %f %f\n",el,xx+0.5,yy+0.5,0.);
         }
      }
      fprintf(dataFile,"\n");
   }
   fclose(dataFile);
   /* write commands for plotting by gnuplot */
   FILE *cmdFile = fopen("/tmp/commandPlotWind1.gnuplot","w");
   fprintf(cmdFile,"set terminal postscript eps 25 colour\n");
   fprintf(cmdFile,"set output '/tmp/window1.eps'\n");
   fprintf(cmdFile,"set xlabel 'X'\n");
   fprintf(cmdFile,"set ylabel 'Y'\n");
   fprintf(cmdFile,"set nokey\n");
   fprintf(cmdFile,"set palette gray\n");
   fprintf(cmdFile,"set xrange [%d:%d]\n",firstX,-firstX);
   fprintf(cmdFile,"set yrange [%d:%d]\n",firstY,-firstY);
   fprintf(cmdFile,"set pm3d map corners2color c4\n");
   fprintf(cmdFile,"set size square\n");
   fprintf(cmdFile,"splot '/tmp/window1.dat' using 2:3:4\n");
   fclose(cmdFile);

   /* launch gunplot */
   ret = system("gnuplot /tmp/commandPlotWind1.gnuplot");
   printf("Image /tmp/window1.eps is ready.\n"); 



#endif /* GNUPLOT */


   applyForwardModel_onevector(x,
         NDRIFTPIXELS,nbWaveBands,e,center_coord,
         obs[BEG], TCflag[BEG], icelandmask[BEG],img_dims, 
         pattern_size, pattern_mask, pattern_windexes, pattern_img, pattern_isvalid,
         &olat, &olon);

   cptValid=0;
   for (o = 0 ; o < pattern_size ; o++) {
      if ( pattern_mask[o] && pattern_isvalid2[o] && pattern_isvalid[o] ) {
         cptValid ++;
      }
   }
   if (!cptValid) {
      /* no valid pixel in the end image. just stop now */
      *fc = 10000;
      return;
   }

   /* from (olat,olon): compute the drift length */
   wind = center_coord[TDIM];
   compute_distance(img_lat[wind],img_lon[wind],olat,olon,&dlen);

   /* decide if the unknown values are valid (ccost = soft constraints) */
   /* ***************************************************************** */
   unknownsAreValid(1,xloc,&isvalid,&ccost,dlen);

   /* compare this modified image to the end-of-drift image */
   //printf("##########\n");
   for (o = 0 ; o < pattern_size ; o++) {
      //printf("Pattern pixel %ld/%u: (%d,%d,%d): ",o,pattern_size,pattern_mask[o],pattern_isvalid[o],pattern_isvalid2[o]);
      if ( pattern_mask[o] && pattern_isvalid[o] && pattern_isvalid2[o] ) {
         odif[o] = pattern_img[o] - pattern_img2[o];
         //printf("Comparing %f with %f.\n",pattern_img[o],pattern_img2[o]);
      } else {
         //printf("Is invalid.\n");
      }
   }
   //printf("##########\n");

   for (o = 0 ; o < pattern_size ; o++) {
      if ( pattern_mask[o] && pattern_isvalid[o] && pattern_isvalid2[o] ) {
         ocost += odif[o]* (1./pow(0.5,2.)) *odif[o];
      }
   }

   *fc = 0.5 * (xcost + ocost + ccost);

#endif /* TMPREMOVE */
}


void model(int *n, double *x,double *fc) {


   long e,e1,e2,o,o1,o2;
   long pix;
   double vfc, ifc, tfc;
   long wind;
   double vciobs;
   int isvalid;
   double *xloc;

   double xcost,xcost_tot;
   double ocost,ocost_tot;
   double ccost,ccost_tot;


   unsigned int tote = NDRIFTPIXELS*NUNKNOWNS;
   unsigned int nobs;

   vfc = 0;
   for ( e = 0 ; e < NDRIFTPIXELS ; e++ ) {
      model_onevector(n,x,&ifc,e);
      vfc += ifc;
   }

   /* constraints on the geographical regularity of the fields */
   /************************************************************/
   /*                      (none)                              */
   tfc = 0;

   /* total cost function */
   /* ******************* */
   *fc = vfc + tfc;

}

