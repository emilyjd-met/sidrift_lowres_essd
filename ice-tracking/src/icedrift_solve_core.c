
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <projects.h>
#include <fmutil.h>
#include <errorcodes.h>
#include "icedrift_flags.h"
#include "icedrift_common.h"
#include "icedrift_prepost.h"
#include "icedrift_model.h"
#include "icedrift_solve_filter.h"
#include "optimization_simplex.h"
#include "icedrift_solve_common.h"
#include "icedrift_solve_core.h"
#include <signal.h>

#undef USE_NWP_FILTER

/* global variables for communications with the I/O in icedrift_prepost.c */
char *OptimMetric;
double radiusNeighbours;

/* global variables for communications with the 'wrappers' in icedrift_solve_common.c */
int     xdim;
double *x;
size_t  currentVector;
short   currentPattern;

static PJ     *out_pj;

char progname[] = "icedrift_solve_core";

int findBestX_cont (size_t p, short pat, double *startPoints[NUNKNOWNS],
		    double xcent, double ycent, double maxdist, double x[],
		    double xbest[], double bestval[], short processingflag[] );

int compute_Posterior_Uncertainties(size_t p, short pat, double xbest[],
				    double bestval[], double *dx_stddev,
				    double *dy_stddev, double *xy_correl,
				    short uncertaintyflag[]);

int core(void) {

   int ret;
   int n;
   double *xbest, *xloc;
   double *x_stddev, *y_stddev, *xy_corr;
   short  *processingflag, *uncertaintyflag;
   double *bestval;

   fmlogmsg(progname,"Start processing.");

   int numberUnprocReported=0;
  
   size_t ndims = NUNKNOWNS;

   xdim = NUNKNOWNS*NDRIFTPIXELS;
      
   /* ======= START CORE ====== */

   /* open the log file (whose name we got from the parameter file) */
   printf("\tLog file is <%s>\n",reportFile);
   FILE *logf = fopen(reportFile,"w");
   if (!logf) {
      fmerrmsg(progname,"Cannot open the log file for writing.");
      exit(OSISAF_ERROR_OTHER);
   }
   
   /* allocate some 'starting points' to the problem */
   double **startpoints;
   startpoints = fmMalloc((ndims+1)*sizeof(double *));
   for (size_t nv = 0 ; nv <= ndims ; nv++) {
      startpoints[nv] = fmMalloc(ndims*sizeof(double));
   }

   x        = fmMalloc(sizeof(double)*xdim);
   xbest    = fmMalloc(sizeof(double)*xdim);
   x_stddev = fmMalloc(sizeof(double)*xdim);
   y_stddev = fmMalloc(sizeof(double)*xdim);
   xy_corr  = fmMalloc(sizeof(double)*xdim);
   for (int i = 0 ; i < xdim ; i++) {
      xbest[i]    = -2500.;
      x_stddev[i] = -2500.;
      y_stddev[i] = -2500.;
      xy_corr[i]  = -2500.; 
   }
   processingflag  = fmMalloc(sizeof(*processingflag)*xdim);
   uncertaintyflag = fmMalloc(sizeof(*uncertaintyflag)*xdim);

   /* Calling initmod to set up some arrays and precompute the pattern shape */
   initmod(&xdim,x);

   char optiC;
   int (*model_wrapper)(size_t,double *,double *);
   ret = choose_model(OptimMetric,&optiC,&model_wrapper);
   if (ret) {
      fprintf(stderr,"Unable to choose optimisation problem from MSD and MCC.\n");
      exit(OSISAF_ERROR_CMDVALUE);
   }
   bestval = fmMalloc(sizeof(double) * NDRIFTPIXELS);
   for (size_t p = 0 ; p < NDRIFTPIXELS ; p++ ) {
      bestval[p] = DBL_MAX;
   }
   
   /* startpoints: */
   /* we have 2 dimensions (xdrift and ydrift) so we need 3 vertex for our simplex */
   double startPointsRadius = maxdriftdistance*0.25;
   double stepang = 360/3.;
   startpoints[0][XDRIFT_IDX] = startPointsRadius*cos(0*stepang*DEG_TO_RAD);
   startpoints[0][YDRIFT_IDX] = startPointsRadius*sin(0*stepang*DEG_TO_RAD); 
   startpoints[1][XDRIFT_IDX] = startPointsRadius*cos(1*stepang*DEG_TO_RAD);
   startpoints[1][YDRIFT_IDX] = startPointsRadius*sin(1*stepang*DEG_TO_RAD);
   startpoints[2][XDRIFT_IDX] = startPointsRadius*cos(2*stepang*DEG_TO_RAD);
   startpoints[2][YDRIFT_IDX] = startPointsRadius*sin(2*stepang*DEG_TO_RAD);
   
   /* prepare and initialize the variables for the filtering steps */
   short nb_neighbours_limit = 5;
   double correlation_limit = 0.3;
   double correlation_limit_for_avg = 0.3;
   double local_stddev_limit = +99999999.; // km : in practice such a large value deactivate the test
   size_t nbNeighbours,neighboursCenter;
   short *neighboursMask;
   long  *neighboursIndexes[3];
   ret = setNeighbourhoodPattern(radiusNeighbours,&nbNeighbours,&neighboursCenter,&neighboursMask,neighboursIndexes);
   if (ret) {
      fmerrmsg(progname,"Cannot set the neighborhood pattern of radius %f km. \n",radiusNeighbours);
      exit(OSISAF_ERROR_OTHER);
   }
   double *neighboursX  = fmMalloc(nbNeighbours*sizeof(double));
   double *neighboursY  = fmMalloc(nbNeighbours*sizeof(double));

   /* initialize all the drift position to 'unprocessed' (to be processed) */
   size_t *tmppointer;
   size_t *tobeprocessed = fmMalloc(NDRIFTPIXELS*sizeof(size_t));
   short  *pattern       = fmMalloc(NDRIFTPIXELS*sizeof(short));
   size_t nb_tobeprocessed = 0;
   for (size_t p = 0 ; p < NDRIFTPIXELS ; p++) {
      processingflag[p]  = ICEDRIFT_UNPROCESSED;
      uncertaintyflag[p] = ICEDRIFTPOST_UNPROCESSED;
      tobeprocessed[p]  = p;
      pattern[p]        = -1;
      nb_tobeprocessed++;
   }
   short atleast_one_unprocessed = (nb_tobeprocessed!=0);
   if (!atleast_one_unprocessed) goto finish_processing;

   size_t *new_tobeprocessed = fmMalloc(NDRIFTPIXELS*sizeof(size_t));
   size_t new_nbtobeprocessed;
   /* go through all the locations and flag those which are too close to the borders of the image grid (or outside) */
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      short pat = 0;
      size_t p = tobeprocessed[q];
      long xi,yi;

      if (iwcs[p] == img_dims[TDIM]) {
         /* the center of the pattern is outside the image grid. flag and to not use */
         processingflag[p] = ICEDRIFT_OUTSIDE_IMGBORDER;
         continue;
      }
      locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
      center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];
      pattern[p] = pat;

      load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
            center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
            pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

      long cptOut=0;
      for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
         if ( pattern_mask[pat][o] && (pattern_isvalid[pat][0][o] == TCIMAGE_OUTSIDE_GRID) ) {
            /* found one point in the pattern that is outside the image grid */
            cptOut++;
            break;
         }
      }
      if (cptOut != 0) {
         processingflag[p] = ICEDRIFT_CLOSETO_IMGBORDER;
         continue;
      }

      /* store the points who passed all the tests */
      new_tobeprocessed[new_nbtobeprocessed] = p;
      new_nbtobeprocessed++;
   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;

   /* go through all the locations and flag those whose central point is over land */
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];

      if (pointIsOnLand(iwcs[p],icelandmask[BEG])) {
         processingflag[p] = ICEDRIFT_CENTER_OVER_LAND;
      } else {
         new_tobeprocessed[new_nbtobeprocessed] = p;
         new_nbtobeprocessed++;
      }
   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;

   /* go through all the remaining locations and flag those for which
    * all pixels in the subimage are 1) nodata or 2) over the ocean or land */
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {

      size_t p = tobeprocessed[q];
      long xi,yi;
      locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
      center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];

      /* First work with the largest pattern (index 0). */
      short pat = 0;
      pattern[p] = pat;

      load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
            center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
            pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

      load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
            center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
            pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);        
      long cptUnproc=0;
      long cptMissing=0;
      long cptIce=0;
      long cptPoints=0;
      for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
         if ( pattern_mask[pat][o] ) {
            cptPoints++;
            if ( (pattern_isvalid[pat][0][o] == TCIMAGE_NODATA) || (pattern_isvalid2[pat][0][o] == TCIMAGE_NODATA) ) {
               cptMissing++;
            }
            if (((pattern_icemask[pat][o] == 3) || (pattern_icemask[pat][o] == 2)) && 
                  ((pattern_icemask2[pat][o] == 3) || (pattern_icemask2[pat][o] == 2))) {
               cptIce++;
            }
            /*
            if ((pattern_icemask[pat][o] == 3) && (pattern_icemask2[pat][o] == 3)) {
               cptIce++;
            }
            */
            if ( (pattern_isvalid[pat][0][o] == TCIMAGE_UNPROCESSED) || (pattern_isvalid2[pat][0][o] == TCIMAGE_UNPROCESSED) ) {
               cptUnproc++;
            }
         }
      }
      if ( cptIce == 0 ) {
         processingflag[p] = ICEDRIFT_NOICE;
         continue;
      } 

      /* no everything is ice. */
      if ( cptIce < cptPoints ) { 

         if (pat != (NBPATTERNS-1)) {
            /* try with the second pattern */
            pat = 1;
            pattern[p] = pat;
            load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
                  center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
                  pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

            load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
                  center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
                  pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);     
            cptMissing=0;
            cptUnproc=0;
            cptIce=0;
            cptPoints=0;
            for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
               if ( pattern_mask[pat][o] ) {
                  cptPoints++;
                  if ( (pattern_isvalid[pat][0][o] == TCIMAGE_NODATA) || (pattern_isvalid2[pat][0][o] == TCIMAGE_NODATA) ) {
                     cptMissing++;
                  }
                  if (((pattern_icemask[pat][o] == 3) || (pattern_icemask[pat][o] == 2)) && 
                        ((pattern_icemask2[pat][o] == 3) || (pattern_icemask2[pat][o] == 2))) {
                     cptIce++;
                  }
                  /*
                  if ((pattern_icemask[pat][o] == 3) && (pattern_icemask2[pat][o] == 3)) {
                     cptIce++;
                  }
                  */
                  if ( (pattern_isvalid[pat][0][o] == TCIMAGE_UNPROCESSED) 
                        || (pattern_isvalid2[pat][0][o] == TCIMAGE_UNPROCESSED) ) {
                     cptUnproc++;
                  }
               }
            }
         } 
      } 
      /* check for ice */

      /* if not enough ice is found, flag and skip this drift location */
      if (cptIce != cptPoints) {
         processingflag[p] = ICEDRIFT_CLOSETO_COAST_OR_EDGE;
         continue;
      } 
      /* => if we make it here, all the pixels in the pattern (index 'pat') are supposedly ice */
      /*    we can now check if there are missing data */
      if ( cptUnproc != 0 ) {
         /* this is bad: it shows that we avoided to process (tc or laplacian) an area which would 
          * * be usefull. Use a special flag and report to log file
          * */
         if ( numberUnprocReported < 10 ) {
            fprintf(stderr,"WARNING (%s) Found an interesting UNPROCESSED area. Check the flags to TC daily images\n",progname);
         } else if (numberUnprocReported == 10) {
            fprintf(stderr,"WARNING (%s) Stop reporting 'Unprocessed area' (%d).\n",progname,numberUnprocReported);
         } else {
            ;
         }
         numberUnprocReported++;
         processingflag[p] = ICEDRIFT_CLOSETO_UNPROCESSED;
         continue;
      }

      if ( cptMissing != 0 ) {

         /* there are some missing data. Can we try reducing the pattern? */
         if (pat != (NBPATTERNS-1)) {
            /* load the smaller pattern and check again */
            pat = 1;
            pattern[p] = pat;

            load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
                  center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
                  pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

            load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
                  center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
                  pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);      
            cptMissing=0;
            cptPoints=0;
            for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
               if ( pattern_mask[pat][o] ) {
                  cptPoints++;
                  if ( (pattern_isvalid[pat][0][o] == TCIMAGE_NODATA) || (pattern_isvalid2[pat][0][o] == TCIMAGE_NODATA) ) {
                     cptMissing++;
                  }
               }
            }
         }
         if (cptMissing != 0) {
            /* still some missing points. Cannot process. */
            processingflag[p] = ICEDRIFT_CLOSETO_MISSING;
            continue;
         } 
         /* else, we try to keep this location */
      } 


      /* keep those who passed the previous tests (and record their pattern number) */
      new_tobeprocessed[new_nbtobeprocessed] = p;
      new_nbtobeprocessed++;

   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;

   /* process until no drift location is unprocessed */ 
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {

      size_t p = tobeprocessed[q];

      /*
         printf ("P: %u, Q: %u\n",p,q);
         if (!(p%500))
         printf("p: %u (q = %u), FLAG = %d\n",p,q,processingflag[p]);
         */

      short pat = pattern[p];

      if (!( (pat >= 0) && (pat <NBPATTERNS) )) {
         printf("ERROR with pattern. (p: %u q:%u)\n",p,q);
         processingflag[p] = ICEDRIFT_FAILS;
         continue;
      }

      /* load the subimage in end_drift image */
      long xi,yi;
      locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
      center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];

      load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
            center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
            pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);

      if ((img_lat == NULL) || (img_lon == NULL)) {
         printf("ERROR with img_lat (%p) or img_lon (%p)\n",img_lat,img_lon);
      }

      /* set the correlation function parameters (sigmoid) */
      setBestKnowledge(img_lat[iwcs[p]],img_lon[iwcs[p]]);

      /* perform the optimization */
      int oret = findBestX_cont (p,pat,startpoints,0.,0.,maxdriftdistance,
            x,xbest,bestval,processingflag);

      /* compute the posterior uncertainties */
      if (processingflag[p] == ICEDRIFT_OK) {
         double sdevx,sdevy,xycorr;
         int pret = compute_Posterior_Uncertainties(p,pat,xbest,bestval,&sdevx,&sdevy,&xycorr,uncertaintyflag);
         if (uncertaintyflag[p] == ICEDRIFTPOST_OK) {
            x_stddev[p] = sdevx;
            y_stddev[p] = sdevy;
            xy_corr[p]  = xycorr;
         }
      } else {
         uncertaintyflag[p] = ICEDRIFTPOST_NOVECTOR;
      }


      new_tobeprocessed[new_nbtobeprocessed] = p;
      new_nbtobeprocessed++;

   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;

   /* post-process and transform the drift result */
   float *latB = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *lonB = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *latE = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *lonE = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *leng = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *dire = fmMalloc(NDRIFTPIXELS * sizeof(float));

   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];
      long xcoordp,ycoordp;
      locfmijmap(p,out_dims[XDIM],&xcoordp,&ycoordp);

      double xdrift = xbest[p*NUNKNOWNS + XDRIFT_IDX];
      double ydrift = xbest[p*NUNKNOWNS + YDRIFT_IDX];   

      /* compute the latitude and longitude of start point */
      latB[p] = img_lat[iwcs[p]]; lonB[p] = img_lon[iwcs[p]];
      /* compute the latitude and longitude of final point */ 
      double lat,lon;
      long xcoordi,ycoordi;
      locfmijmap(iwcs[p],img_dims[XDIM],&xcoordi,&ycoordi);
      remap_xy2ll(-(xdrift/img_Ax)+xcoordi,-(ydrift/img_Ay)+ycoordi,img_pj,img_Ax,img_Bx,img_Ay,img_By,&lat,&lon);
      latE[p] = lat; lonE[p] = lon;

      /*
         printf("ICEDRIFT is x:%f, y:%f. (%f,%f)(%f,%f) -> (%f,%f)(%f %f)\n",
         xdrift,ydrift,
         (float)xcoordi,(float)ycoordi,
         latB[p],lonB[p],
         -(xdrift/img_Ax)+xcoordi,-(ydrift/img_Ay)+ycoordi,
         latE[p],lonE[p]);
         */

      /* compute the length and direction of the arrows */ 
      double dlen,ddir;
      compute_distance(latB[p],lonB[p],latE[p],lonE[p],&dlen);
      compute_directionToNorth(latB[p],lonB[p],latE[p],lonE[p],&ddir);
      leng[p] = dlen; dire[p] = ddir;
   }

   /* compute the local average drift vector */ 
   float  *xDrift = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *yDrift = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *xavg = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *yavg = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *xstd = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *ystd = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *avgLat = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *avgLon = fmMalloc(NDRIFTPIXELS * sizeof(float));
   size_t *navg = fmMalloc(NDRIFTPIXELS * sizeof(size_t));

   float *len_avg  = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *len_diff = fmMalloc(NDRIFTPIXELS * sizeof(float));

   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {

      size_t p = tobeprocessed[q];
      long xcoordp,ycoordp;
      locfmijmap(p,out_dims[XDIM],&xcoordp,&ycoordp);

      size_t nbValidNeighbours   = 0;
      size_t nbInvalidNeighbours = 0;
      for (size_t n = 0 ; n < nbNeighbours ; n++) {
         if (neighboursMask[n]) {

            long xneighbour = xcoordp + neighboursIndexes[XDIM][n];
            long yneighbour = ycoordp + neighboursIndexes[YDIM][n];

            if ( !pointIsInGrid((double)(xneighbour),(double)(yneighbour),out_dims) ) continue;

            long windexNeighbour;
            locfmivec(&windexNeighbour,xneighbour,yneighbour,out_dims[XDIM]);

            if ( (processingflag[windexNeighbour] == ICEDRIFT_OK) && 
                  (bestval[windexNeighbour] >= correlation_limit_for_avg) ) {
               neighboursX[nbValidNeighbours] = xbest[windexNeighbour*NUNKNOWNS + XDRIFT_IDX];
               neighboursY[nbValidNeighbours] = xbest[windexNeighbour*NUNKNOWNS + YDRIFT_IDX];
               nbValidNeighbours++;
            } else if (processingflag[windexNeighbour] == ICEDRIFT_CLOSETO_MISSING) {
               ;
            } else {
               nbInvalidNeighbours++;
            }
         }
      }

      navg[p] = nbValidNeighbours;
      xavg[p] = 1000; yavg[p] = 1000;
      xstd[p] = 1000; ystd[p] = 1000;
      double xavgP,yavgP,xstdP,ystdP;
      ret = compute_mean_vector(&xavgP,&yavgP,&xstdP,&ystdP,nbValidNeighbours,neighboursX,neighboursY);
      xavg[p] = xavgP;yavg[p] = yavgP;
      xstd[p] = xstdP;ystd[p] = ystdP;
      if ((nbValidNeighbours >= nb_neighbours_limit) && (xstdP < local_stddev_limit) && (ystdP < local_stddev_limit)) {
         /* we have a local average drift vector we can trust */

         /* get the lat and lon associated with the average drift */
         long xcoordi,ycoordi;
         double latAvg,lonAvg;
         locfmijmap(iwcs[p],img_dims[XDIM],&xcoordi,&ycoordi);
         remap_xy2ll(-(xavgP/img_Ax)+xcoordi,-(yavgP/img_Ay)+ycoordi,img_pj,img_Ax,img_Bx,img_Ay,img_By,&latAvg,&lonAvg);
         avgLat[p] = latAvg;
         avgLon[p] = lonAvg;

         /* compute the length of the average drift */
         double lenavg;
         compute_distance(latB[p],lonB[p],latAvg,lonAvg,&lenavg);

         /* compute the length of the difference vector: avg - drift */
         double lendiff;
         compute_distance(latAvg,lonAvg,latE[p],lonE[p],&lendiff);

         /* store those values for writing in the file */
         len_avg[p]  = lenavg;
         len_diff[p] = lendiff;
      } else {
         len_diff[p] = -1.; 
      }
   }

   /* filtering with neighbours. */
   double difflenLimit = 10;
   int invalid_exists = 1;
   size_t nbcorr = 0;
   while (invalid_exists) {
      /* sort by order of diff_length */
      ret = sort_icedrift_bylength(nb_tobeprocessed,tobeprocessed,len_diff);
      /* test the first (and thus worst) candidate */
      size_t p = tobeprocessed[0];
      if (len_diff[p] > difflenLimit) { 

         if (processingflag[p] == ICEDRIFT_REFUSED_BY_NEIGHBOURS) {
            /* correction already attempted here */
            len_diff[p] = -3;
            continue;
         }

         nbcorr++;
         //printf("Start a new inversion for location p:%u (len_diff[p]=%f)\n",p,len_diff[p]);
         /* correct (modify or delete this vector) */
         processingflag[p]  = ICEDRIFT_REFUSED_BY_NEIGHBOURS; 
         uncertaintyflag[p] = ICEDRIFTPOST_NOVECTOR;
         double currentBestX = xbest[p*NUNKNOWNS+XDRIFT_IDX];
         double currentBestY = xbest[p*NUNKNOWNS+YDRIFT_IDX];
         double currentBestC = bestval[p];

         /* load the subimage in end_drift image */
         short pat = pattern[p];
         long xi,yi;
         locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
         center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];
         load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
               center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
               pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);

         /* modify the sigmoid  */
         set_sigmoidlength(difflenLimit);
         /* set the correlation function parameters (sigmoid center) */
         setBestKnowledge(avgLat[p],avgLon[p]);

         /* perform the optimization */
         int oret = findBestX_cont (p,pat,startpoints,xavg[p],yavg[p],difflenLimit,
               x,xbest,bestval,processingflag);

         if ( (processingflag[p] != ICEDRIFT_OK) || (bestval[p] < correlation_limit) ) {
            /* Optimization failed, or does not result in a good enough cross-correlation. Remove the vector. */
            processingflag[p]  = ICEDRIFT_REFUSED_BY_NEIGHBOURS;
            uncertaintyflag[p] = ICEDRIFTPOST_NOVECTOR;
         } else {
            processingflag[p]               = ICEDRIFT_CORRECT_BY_NEIGHBOURS;

            /* compute the posterior uncertainties */
            double sdevx,sdevy,xycorr;
            int pret = compute_Posterior_Uncertainties(p,pat,xbest,bestval,&sdevx,&sdevy,&xycorr,uncertaintyflag);
            if (uncertaintyflag[p] == ICEDRIFTPOST_OK) {
               x_stddev[p] = sdevx;
               y_stddev[p] = sdevy;
               xy_corr[p]  = xycorr;
            }

         }

         /* recompute the average for all the neighbouring pixels (which do not include 'p')*/
         long xcoordp,ycoordp;
         locfmijmap(p,out_dims[XDIM],&xcoordp,&ycoordp);
         for (size_t n = 0 ; n < nbNeighbours ; n++) {

            if (neighboursMask[n]) {

               long xneighbour = xcoordp + neighboursIndexes[XDIM][n];
               long yneighbour = ycoordp + neighboursIndexes[YDIM][n];

               if ( !pointIsInGrid((double)(xneighbour),(double)(yneighbour),out_dims) ) continue;

               long windexNeighbour;
               locfmivec(&windexNeighbour,xneighbour,yneighbour,out_dims[XDIM]);
               /* windexNeighbour is the absolute 1D index for the n-th neighbour of 'p' */
               size_t nbValidNeighbours   = 0;
               size_t nbInvalidNeighbours = 0;
               /* recompute the average around windexNeighbour */
               for (size_t np = 0 ; np < nbNeighbours ; np++) {
                  if (neighboursMask[np]) {

                     long xneighbour_np_n = xneighbour + neighboursIndexes[XDIM][np];
                     long yneighbour_np_n = yneighbour + neighboursIndexes[YDIM][np];

                     if ( !pointIsInGrid((double)(xneighbour_np_n),(double)(yneighbour_np_n),out_dims) ) continue;

                     long windexNeighbour_np;
                     locfmivec(&windexNeighbour_np,xneighbour_np_n,yneighbour_np_n,out_dims[XDIM]);

                     if (  ((processingflag[windexNeighbour_np] == ICEDRIFT_OK) ||
                              (processingflag[windexNeighbour_np] == ICEDRIFT_CORRECT_BY_NEIGHBOURS))
                           && (bestval[windexNeighbour_np] >= correlation_limit_for_avg) ) {
                        neighboursX[nbValidNeighbours] = xbest[windexNeighbour_np*NUNKNOWNS + XDRIFT_IDX];
                        neighboursY[nbValidNeighbours] = xbest[windexNeighbour_np*NUNKNOWNS + YDRIFT_IDX];
                        nbValidNeighbours++;
                     } else if (processingflag[windexNeighbour_np] == ICEDRIFT_CLOSETO_MISSING) {
                        ;
                     } else {
                        nbInvalidNeighbours++;
                     }
                  }
               }
               navg[windexNeighbour] = nbValidNeighbours;
               xavg[windexNeighbour] = 1000; yavg[windexNeighbour] = 1000;
               xstd[windexNeighbour] = 1000; ystd[windexNeighbour] = 1000;
               double xavgP,yavgP,xstdP,ystdP;
               ret = compute_mean_vector(&xavgP,&yavgP,&xstdP,&ystdP,nbValidNeighbours,neighboursX,neighboursY);
               xavg[windexNeighbour] = xavgP;yavg[windexNeighbour] = yavgP;
               xstd[windexNeighbour] = xstdP;ystd[windexNeighbour] = ystdP;
               if ((nbValidNeighbours >= nb_neighbours_limit) && (xstd[windexNeighbour] < local_stddev_limit) && (ystd[windexNeighbour] < local_stddev_limit)) {

                  /* get the lat and lon associated with the average drift */
                  long xcoordi,ycoordi;
                  double latAvg,lonAvg;
                  locfmijmap(iwcs[windexNeighbour],img_dims[XDIM],&xcoordi,&ycoordi);
                  remap_xy2ll(-(xavgP/img_Ax)+xcoordi,-(yavgP/img_Ay)+ycoordi,img_pj,img_Ax,img_Bx,img_Ay,img_By,&latAvg,&lonAvg);
                  avgLat[windexNeighbour] = latAvg;
                  avgLon[windexNeighbour] = lonAvg;

                  /* compute the length of the average drift */
                  double lenavg;
                  compute_distance(latB[windexNeighbour],lonB[windexNeighbour],latAvg,lonAvg,&lenavg);

                  /* compute the length of the difference vector: avg - drift */
                  double lendiff;
                  compute_distance(latAvg,lonAvg,latE[windexNeighbour],lonE[windexNeighbour],&lendiff);

                  /* store those values for writing in the file */
                  len_avg[windexNeighbour]  = lenavg;
                  len_diff[windexNeighbour] = lendiff;

               } else {
                  len_diff[windexNeighbour] = -2.; 
               }
            }
         }

         /* transform the new vector (in 'p') to len,dir, lat and lon and recompute its difference to the average in 'p' */
         if (processingflag[p] == ICEDRIFT_CORRECT_BY_NEIGHBOURS) {
            double xdrift = xbest[p*NUNKNOWNS+XDRIFT_IDX];
            double ydrift = xbest[p*NUNKNOWNS+YDRIFT_IDX];;
            double lat,lon;
            long xcoordi,ycoordi;
            locfmijmap(iwcs[p],img_dims[XDIM],&xcoordi,&ycoordi);
            remap_xy2ll(-(xdrift/img_Ax)+xcoordi,-(ydrift/img_Ay)+ycoordi,img_pj,img_Ax,img_Bx,img_Ay,img_By,&lat,&lon);
            latE[p] = lat; lonE[p] = lon;

            /* compute the length and direction of the arrow */ 
            double dlen,ddir;
            compute_distance(latB[p],lonB[p],latE[p],lonE[p],&dlen);
            compute_directionToNorth(latB[p],lonB[p],latE[p],lonE[p],&ddir);
            leng[p] = dlen; dire[p] = ddir;

            /* compute length of the difference vector to average (the average vector in 'p' is not modified) */
            double lendiff;
            compute_distance(avgLat[p],avgLon[p],latE[p],lonE[p],&lendiff);

            /* store those values for writing in the file */
            len_diff[p] = lendiff;
            //printf("p:%u. New distance to average is %f\n",p,len_diff[p]);
         } else {
            len_diff[p] = -3;
            //printf("p:%u. Set len_diff to %f\n",p,len_diff[p]);
         }

      } else {
         /* the worst vector passes the test. So we can stop filtering. */
         invalid_exists=0;
      }
   }

   for (size_t p = 0 ; p < NDRIFTPIXELS ; p++) {
      if ( (processingflag[p] != ICEDRIFT_OK) && (processingflag[p] != ICEDRIFT_CORRECT_BY_NEIGHBOURS) ) continue;
      //if ( (navg[p] < nb_neighbours_limit) || (xavg[p] == 1000) || (len_diff[p] < 0) ) {
      if ( (navg[p] < nb_neighbours_limit) || (xavg[p] == 1000) ) {
         processingflag[p]  = ICEDRIFT_NOAVERAGE;
         uncertaintyflag[p] = ICEDRIFTPOST_NOVECTOR;
      }
   }

finish_processing:
  /* last level of filtering: remove the vectors for which we still have no average vector as
    * well as those with too bad a correlation. */
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];
      if ((processingflag[p] != ICEDRIFT_OK) && (processingflag[p] != ICEDRIFT_CORRECT_BY_NEIGHBOURS)) continue;
      if (bestval[p] < correlation_limit) {
         processingflag[p]  = ICEDRIFT_LOWCORRELATION;
         uncertaintyflag[p] = ICEDRIFTPOST_NOVECTOR;
      }
   }

   /* Calling the prep code (writing will be done in Python) */
   ret = prepare_icedriftProduct(&xdim, xbest, bestval, processingflag,
				 pattern, leng, dire, latB, lonB, latE,
				 lonE, navg, xavg, yavg, len_avg, len_diff,
				 xstd, ystd, x_stddev, y_stddev, xy_corr,
				 uncertaintyflag);
   if (ret) {
      printf("error in preparing output\n");
   }
   
   fclose(logf);

   free(x);
   free(xbest);
   for (size_t nv = 0 ; nv <= ndims ; nv++) {
      free(startpoints[nv]);
   }

   fmlogmsg(progname,"End processing.");
   
   /* ======= END CORE ====== */

   return EXIT_SUCCESS;
}

int findBestX_cont (size_t p, short pat, double *startPoints[NUNKNOWNS], double xcent, double ycent, double maxdist,
      double x[], double xbest[], double bestval[], short processingflag[] ) {

   /* 
    * sample the correlation function around the (xcent,ycent) point to
    * find good starting points 
    */
   double bestCorrs[3],bestXs[3],bestYs[3];
   bestCorrs[0] = bestCorrs[1] = bestCorrs[2] = -1;

   size_t sta = p*NUNKNOWNS;
   double stepang = 45.; 
   double steplen = 10.;
   double sampLen = 5.;
   short  nbSampLen = 0;
   while (sampLen < maxdist) {

      double sampAng = (stepang*0.5)*(nbSampLen%2);
      while (sampAng < 360) {

         /* try this point */
         double val;
         x[sta+XDRIFT_IDX] = xcent+sampLen*cos(sampAng*DEG_TO_RAD);
         x[sta+YDRIFT_IDX] = ycent+sampLen*sin(sampAng*DEG_TO_RAD);
         model_onevector_corr(&xdim,x,&val,p,pat);
         short t;
         for (t=0 ; t < 3 ; t++) {
            if (val > bestCorrs[t]) break;
         }
         if (t <= 2) {
            for (short t2=2 ; t2 > t ; t2--) {
               bestCorrs[t2] = bestCorrs[t2-1];
               bestXs[t2]    = bestXs[t2-1];
               bestYs[t2]    = bestYs[t2-1];
            }
            bestCorrs[t] = val;
            bestXs[t]    = x[sta+XDRIFT_IDX];
            bestYs[t]    = x[sta+YDRIFT_IDX];
         }
         sampAng+=stepang;
      }
      sampLen+=steplen;
      nbSampLen++;
   }

   /* perform the SIMPLEX optimisation, starting from the best 3 points */
   //printf("Start points: \n");
   for (short t = 0 ; t < 3 ; t++) {
      startPoints[t][XDRIFT_IDX] = bestXs[t];
      startPoints[t][YDRIFT_IDX] = bestYs[t];
      // printf("\t%f,%f => %f\n",startPoints[t][XDRIFT_IDX],startPoints[t][YDRIFT_IDX],bestCorrs[t]);
   }
   currentVector  = p;
   currentPattern = pat;
   double xloc[NUNKNOWNS];
   double bestscore;
   int oret; 
   oret = findBestX(NUNKNOWNS,xloc,&bestscore,startPoints,'>',model_wrapper_CC);
   if (oret) {
      processingflag[p] = ICEDRIFT_FAILS;
   } else {
      xbest[sta+XDRIFT_IDX] = xloc[XDRIFT_IDX];
      xbest[sta+YDRIFT_IDX] = xloc[YDRIFT_IDX];
      processingflag[p] = ICEDRIFT_OK;
   }
   bestval[p] = bestscore;

   return 0;
}

/* this routine computes the a-posteriori uncertainty estimates from the Hessian of the correlation
   function around the best estimate. */
int compute_Posterior_Uncertainties(size_t p, short pat, double xbest[], double bestval[], 
      double *dx_stddev, double *dy_stddev, double *xy_correl, short uncertaintyflag[]) {

   int oret,err;


   /* length of the local perturbations */
   double epsilon = 0.0001;

   /* initialize the global variables for communications with the wrappers in icedrift_solve_common */
   currentVector  = p;
   currentPattern = pat;

   /* save the values that are in the x vector */ 
   size_t currLoc = currentVector*NUNKNOWNS; 
   double xd_sav  = xbest[currLoc+XDRIFT_IDX];
   double yd_sav  = xbest[currLoc+YDRIFT_IDX];

   /* have a NUNKNOWNS-element vector to hold the local perturbations */
   double xLoc[NUNKNOWNS];

   /* we use a 1 node approximation to the second derivative terms */
   double d2fdx2, d2fdy2; /* (F(t+h) - 2*F(t) + F(t-h))/(h*h) */
   double d2fdxdy; /* (F(x+h,y+k) - F(x+h,y-k) - F(x-h,y+k) + F(x-h,y-k)) / (4*h*k) */

   oret = 0;

   //fprintf(stdout,"\n%s\n",__func__);

   /* compute F(x,y) */
   double Fxy;
   xLoc[XDRIFT_IDX] = xd_sav;
   xLoc[YDRIFT_IDX] = yd_sav;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxy);

   //printf("%s: %u: %f %f\n",__func__,currentVector,Fxy,bestval[p]);

   /* compute F(x+h,y) and F(x-h,y) */
   double Fxphy,Fxmhy;
   xLoc[XDRIFT_IDX] = xd_sav + epsilon;
   xLoc[YDRIFT_IDX] = yd_sav;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxphy);
   xLoc[XDRIFT_IDX] = xd_sav - epsilon;
   xLoc[YDRIFT_IDX] = yd_sav;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxmhy);

   /* compute F(x,y+h) and F(x,y-h) */
   double Fxyph,Fxymh;
   xLoc[XDRIFT_IDX] = xd_sav;
   xLoc[YDRIFT_IDX] = yd_sav + epsilon;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxyph);
   xLoc[XDRIFT_IDX] = xd_sav;
   xLoc[YDRIFT_IDX] = yd_sav - epsilon;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxymh);

   /* compute F(x+h,y+h), F(x+h,y-h), F(x-h,y+h) and F(x-h,y-h) */
   double Fxphyph,Fxphymh,Fxmhyph,Fxmhymh;
   xLoc[XDRIFT_IDX] = xd_sav + epsilon;
   xLoc[YDRIFT_IDX] = yd_sav + epsilon;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxphyph);
   xLoc[XDRIFT_IDX] = xd_sav + epsilon;
   xLoc[YDRIFT_IDX] = yd_sav - epsilon;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxphymh);
   xLoc[XDRIFT_IDX] = xd_sav - epsilon;
   xLoc[YDRIFT_IDX] = yd_sav + epsilon;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxmhyph);
   xLoc[XDRIFT_IDX] = xd_sav - epsilon;
   xLoc[YDRIFT_IDX] = yd_sav - epsilon;
   oret += model_wrapper_SLSD(NUNKNOWNS,xLoc,&Fxmhymh);

   if ( oret > 0 ) {
      fprintf(stderr,"ERROR (%s) [p:%u] In evaluating the local perturbations\n",__func__,p);
      uncertaintyflag[p] = ICEDRIFTPOST_DIF;
      return 1;
   }


   /* final computations for hessian */
   double epsilon_square = epsilon*epsilon;
   d2fdx2  = (Fxphy - 2.*Fxy + Fxmhy)/epsilon_square;
   d2fdy2  = (Fxyph - 2.*Fxy + Fxymh)/epsilon_square;
   d2fdxdy = (Fxphy + Fxmhy + Fxyph + Fxymh - 2*Fxy - Fxphyph - Fxmhymh)/(2.*epsilon_square) ;
   //d2fdxdy = (Fxphyph - Fxphymh - Fxmhyph + Fxmhymh)/(4.*epsilon_square);



   /* inverse the 2-by-2 hessian matrix to get the 2-by-2 covariance matrix */
   double dx_var, dy_var, dxdy_covar;
   double delta = d2fdx2 * d2fdy2 - d2fdxdy * d2fdxdy;
   if (delta == 0) {
      fprintf(stderr,"ERROR (%s) [p:%u] Determinant is zero\n",__func__,p);
      uncertaintyflag[p] = ICEDRIFTPOST_INV;
      return 1;
   } 

   dx_var = d2fdy2 / delta;
   dy_var = d2fdx2 / delta;
   dxdy_covar = - d2fdxdy / delta;
   if ( (dx_var <= 0) || (dy_var <= 0) ) {
      fprintf(stdout,"ERROR (%s) [p:%u] Negative variances (%f,%f)\n",__func__,p,dx_var,dy_var);
      /*
         uncertaintyflag[p] = ICEDRIFTPOST_NEG;
         printf("X-H:%f X+H:%f Y-H:%f Y+H:%f\n",(xd_sav-epsilon),(xd_sav+epsilon),(yd_sav-epsilon),(yd_sav+epsilon));
         printf("%g %g %g\n",Fxmhyph/Fxy,Fxyph/Fxy,Fxmhyph/Fxy);
         printf("%g %g %g\n",Fxmhy/Fxy,Fxy,Fxphy/Fxy);
         printf("%g %g %g\n",Fxmhymh/Fxy,Fxymh/Fxy,Fxphymh/Fxy);
         printf("D2F/DX2=%e, D2F/DY2=%e, D2F/DXDY=%e\n",d2fdx2,d2fdy2,d2fdxdy);
         printf("Determinant=%e\n",delta);
         printf("Var(X)=%e, Var(Y)=%e\n",dx_var,dy_var);
         */
      return 1;
   }
   /* return to calling routine */
   *dx_stddev = sqrt(dx_var);
   *dy_stddev = sqrt(dy_var);
   *xy_correl = dxdy_covar / ((*dx_stddev) * (*dy_stddev));


   uncertaintyflag[p] = ICEDRIFTPOST_OK;
   return 0;
}
