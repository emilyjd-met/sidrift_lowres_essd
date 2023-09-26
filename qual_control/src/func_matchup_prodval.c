
#include <stdlib.h>
#include <stdio.h>
#include <projects.h>
#include <dbl_list.h>
#include <rsprod.h>
#include <fmutil.h>
#include <read_obsfile.h>
#include "icedrift_common.h"
#include "func_drifters_trajectories.h"
#include "func_matchup_prodval.h"

//#define MATCHUP_VERBOSE
#undef MATCHUP_VERBOSE

int matchups_byIJ(matchProdVal *a, matchProdVal *b);
int matchups_byValSrc(matchProdVal *a, matchProdVal *b);
int matchups_byTimeError(matchProdVal *a, matchProdVal *b);
int matchups_byDistance(matchProdVal *a, matchProdVal *b);

void printMatchProdVal(matchProdVal *this) {
   Obsdata *valB = &((this->val[0]).data);
   Obsdata *valE = &((this->val[1]).data);
   
   char      dB_val[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec1970 tB_val = (this->val[0]).time;
   fmsec19702CFepoch(tB_val,dB_val);
   char      dE_val[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec1970 tE_val = (this->val[1]).time;
   fmsec19702CFepoch(tE_val,dE_val);

   char dB_prod[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec19702CFepoch(this->tB_prod,dB_prod);
   char dE_prod[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec19702CFepoch(this->tE_prod,dE_prod);
   
   printf("\tVAL: %s [%s %s] (%s %.3f %.3f) -> (%s %.3f %.3f)\n",
         this->id,this->network,this->source,
         dB_val,valB->lat,valB->lon,
         dE_val,valE->lat,valE->lon);
   printf("\tPROD: %03d-%03d (%s %.3f %.3f) -> (%s %.3f %.3f)\n",
         this->i_prod, this->j_prod,
         dB_prod,this->latB_prod,this->lonB_prod,
         dE_prod,this->latE_prod,this->lonE_prod);
   printf("\tCOLL: %.2f km | %+04.2f h | %+04.2f h\n",
         this->distance,this->dtB,this->dtE);
}

int selectMatchups(dbl_list *allMatchups) {


   /* sort the matchups by: 
    * 1) I,J coords in the OSISAF product, 
    * 2) src of position data record (sar,gps,argos,etc...)
    * 3) temporal collocation skills
    * 4) spatial collocation distances */
   void sort_matchups(dbl_list *allMatchups);
   sort_matchups(allMatchups);

   /* now 'uniq' the list and remove the matchups that have the same I and J */
   int matchups_haveSameIJ(matchProdVal *a, matchProdVal *b);
   long nb_deleted = dbl_list_uniq(allMatchups,&matchups_haveSameIJ,NULL);

   /* TODO: CHECK THIS! */
   //int matchups_haveCloseIJ(matchProdVal *a, matchProdVal *b);  
   //nb_deleted += dbl_list_uniq(allMatchups,&matchups_haveCloseIJ,NULL);

   /* final rank of the matchups before we write them to file */
   int ret = dbl_list_sort(allMatchups,&matchups_byValSrc,NULL,NULL);

   return 0;
}

void sort_matchups(dbl_list *allMatchups) {
   int matchups_total_ranking(matchProdVal *a, matchProdVal *b);
   int ret = dbl_list_sort(allMatchups,&matchups_total_ranking,NULL,NULL);
}

/* RANKING ROUTINES */
int matchups_total_ranking(matchProdVal *a, matchProdVal *b) {
   int retcode;
   if ((retcode = matchups_byIJ(a,b)))        return retcode;
   if ((retcode = matchups_byValSrc(a,b)))    return retcode;
   if ((retcode = matchups_byTimeError(a,b))) return retcode;
   if ((retcode = matchups_byDistance(a,b)))  return retcode;
   return 0;
}

int matchups_byIJ(matchProdVal *a, matchProdVal *b) {
   /* rank by I then by J */
   if (a->i_prod < b->i_prod) {
      return -1;
   } else if (a->i_prod > b->i_prod) {
      return +1;
   } else {
      /* a and b have same I */
      if (a->j_prod < b->j_prod) {
         return -1;
      } else if (a->j_prod > b->j_prod) {
         return +1;
      } else {
         return 0;
      }
   }
}

int matchups_byValSrc(matchProdVal *a, matchProdVal *b) {
   Obsdata *aVal = &((a->val[0]).data);
   Obsdata *bVal = &((b->val[0]).data);
   /* gps > sar > argos */
   if (!strcmp(aVal->stDesc,bVal->stDesc)) return 0;
   if (!strcmp(aVal->stDesc,"gps")) return -1;
   else if (!strcmp(bVal->stDesc,"gps")) return +1;
   else if (!strcmp(aVal->stDesc,"wsm")) return -1;
   else if (!strcmp(bVal->stDesc,"wsm")) return +1;
   else if (!strcmp(aVal->stDesc,"argos")) return -1;
   else if (!strcmp(bVal->stDesc,"argos")) return +1;
   else {
      fprintf(stderr,"WARNING (%s) Know nothing about source '%s' and '%s'\n",__func__,aVal->stDesc,bVal->stDesc);
      return 0;
   }
}

int matchups_byTimeError(matchProdVal *a, matchProdVal *b) {
   /* take the max of the time collocation error (in abs value) */
   float a_dT = (fabs(a->dtB)>fabs(a->dtE)?fabs(a->dtB):fabs(a->dtE));
   float b_dT = (fabs(b->dtB)>fabs(b->dtE)?fabs(b->dtB):fabs(b->dtE));
   /* those with lower time collocation errors are ranked first */
   if (a_dT == b_dT) return 0;
   if (a_dT  < b_dT) return -1;
   return 1;
}

int matchups_byDistance(matchProdVal *a, matchProdVal *b) {
   if (a->distance == b->distance) return 0;
   if (a->distance  < b->distance) return -1;
   return 1;
}

int matchups_haveSameIJ(matchProdVal *a, matchProdVal *b) { 
   return ((a->i_prod == b->i_prod) && (a->j_prod == b->j_prod));
}

int matchups_haveCloseIJ(matchProdVal *a, matchProdVal *b) {
   int yes = ((fabs((int)a->i_prod - (int)b->i_prod) <= 1) && (fabs((int)a->j_prod - (int)b->j_prod) <= 1));
   /* printf("%s: test %03d-%03d with %03d-%03d and answer %d\n",
         __func__,a->i_prod,a->j_prod,b->i_prod,b->j_prod,yes); */
   return yes;
}

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
      dbl_list *allMatchups) {

   int ret;

   /* configuration for the collocation */
   double nearest_neighbour_radius    = NNRADIUS;
   double minimum_bilinear_weight     = BLWEIGHT;
   double maximum_abstime_difference  = MAXTDIFF;
   double maximum_duration_difference = MAXDDIFF;

   fmlogmsg(progname,"Configuration for collocation : T0 ]-%.1f hours:+%.1f hours[, Delta ]-%.1f hours:+%.1f hours[, %s (%.2f)",
         maximum_abstime_difference,maximum_abstime_difference,
         maximum_duration_difference,maximum_duration_difference,
         (use_nearest_neighbour?"NN":"INT"),
         (use_nearest_neighbour?nearest_neighbour_radius:minimum_bilinear_weight));

   double TimeLimitSeconds     = maximum_abstime_difference  * 60 * 60; /* [sec] */
   double DurationLimitSeconds = maximum_duration_difference * 60 * 60; /* [sec] */

   /* Loop through the product grid and record the smallest/largest start times and duration */
   fmsec1970 min_prod_T0 =  refT1;
   fmsec1970 max_prod_T0 =  refT0;
   for (size_t e = 0 ; e < prodDims[2] ; e++) {
      if ( driftX_fprod[e] == fillvalf_fprod ) continue;
      if ( t0_fprod[e] < min_prod_T0 ) {
         min_prod_T0 = t0_fprod[e];
      } else if ( t0_fprod[e] > max_prod_T0 ) {
         max_prod_T0 = t0_fprod[e];
      }
   }

#ifdef MATCHUP_VERBOSE
   char datestrMINT0[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec19702CFepoch(min_prod_T0,datestrMINT0);
   char datestrMAXT0[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec19702CFepoch(max_prod_T0,datestrMAXT0);
   printf("MIN and MAX of T0 in product grid are <%s> and <%s>\n",
         datestrMINT0,datestrMAXT0);
#endif 

   /* Loop through the validation data and perform collocation with product */
   size_t nb_valid_matchups = 0u;
   for (size_t s = 0 ; s < nbStations ; s++) {

      /*
      if (strcmp(trajs[s].id,"itp24")) {
         continue;
      }
      */

#ifdef MATCHUP_VERBOSE
      printf("Start collocation for station %s:\n",trajs[s].id);
#endif /* MATCHUP_VERBOSE */

#ifdef MATCHUP_VERBOSE
      if ( strstr("gcp",trajs[s].source) ) {
         printf("... it is a GCP\n");
      }
#endif

      char tmpdatestrIS[FMUTIL_CFEPOCH_LENGTH+1];
      char tmpdatestr0[FMUTIL_CFEPOCH_LENGTH+1];
      char tmpdatestr1[FMUTIL_CFEPOCH_LENGTH+1];

      /* Implement a 3-dimensional search (lat,lon,time) between the trajectory and the product's (lat,lon,t0) grid */
      dbl_list *records = trajs[s].records;
      dbl_node *startRecord = records->head;
      int recc = 0;
      float bestDistInit = 999999999.;
      float bestDistT0 =  bestDistInit;
      float Duration_Prd;
      buoyRecord *bestValRecT0;
      long  bestCornersT0[4];
      float bestWeightsT0[4];
      short bestClosestCornerT0;
      /* search for the best T0 record in the trajectory */
#ifdef MATCHUP_VERBOSE
      printf("SEARCH FOR T0 match\n");
#endif
      do {
         recc++;
         buoyRecord *rec = startRecord->c;
         Obsdata *o = &(rec->data);

         /* get date&time for record */
         double insitu_time = (double)rec->time;

         if ( strstr("gcp",trajs[s].source) ) {
            insitu_time = min_prod_T0;
         }

#ifdef MATCHUP_VERBOSE
         fmsec19702CFepoch((fmsec1970)round(insitu_time),tmpdatestrIS);
#endif

         /* if record is not corresponding to the min/max values for the product grid, skip it */
         if ( (min_prod_T0 - insitu_time) >  TimeLimitSeconds ) {
#ifdef MATCHUP_VERBOSE
            printf("Time stamp <%s> too far from MIN of product T0 in grid <%s>. Skip.\n",tmpdatestrIS,datestrMINT0);
#endif
            goto next_val_record;
         }

         if ( (insitu_time - max_prod_T0) >  TimeLimitSeconds ) {
#ifdef MATCHUP_VERBOSE
            printf("Time stamp <%s> too far from MAX of product T0 in grid <%s>. Too late for this trajectory.\n",tmpdatestrIS,datestrMAXT0);
#endif
            break;
         }

         /* remap each record of the trajectory in the product grid */
         double xo,yo;
         ret = remap_ll2xy(o->lat,o->lon,proj,A[XDIM],B[XDIM],A[YDIM],B[YDIM],&xo,&yo,0);

         /* get the product start time (t0) around this position (Nearest Neighbour or BiLinear)*/
         int xlow  = floor(xo); 
         int xhig  = ceil(xo);
         int ylow  = floor(yo); 
         int yhig  = ceil(yo);

#ifdef MATCHUP_VERBOSE
         printf ("%d In situ Time: %s Pos: (%f %f) -> (%f,%f) (%d,%d,%u),(%d,%d,%u)\n",
               recc,tmpdatestrIS,
               o->lat,o->lon,xo,yo,xlow,xhig,prodDims[XDIM],ylow,yhig,prodDims[YDIM]);
#endif


         double  t0Interp  = 0;
         double  t1Interp  = 0;

         long  interpCorners[4];
         float interpWeights[4];
         float interpDist[4];
         short interpValid[4];
         if ( !((xlow < 0) || (xlow >= prodDims[XDIM]) ||
                  (xhig < 0) || (xhig >= prodDims[XDIM]) ||
                  (ylow < 0) || (ylow >= prodDims[YDIM]) ||
                  (yhig < 0) || (yhig >= prodDims[YDIM])) )  {

            double xeps = xo - xlow; 
            double yeps = yo - ylow; 

            interpCorners[0] = fmivec(xlow,ylow,prodDims[XDIM]);
            interpCorners[1] = fmivec(xlow,yhig,prodDims[XDIM]);
            interpCorners[2] = fmivec(xhig,ylow,prodDims[XDIM]);
            interpCorners[3] = fmivec(xhig,yhig,prodDims[XDIM]);
            interpWeights[0] = (1. - xeps) * (1. - yeps);
            interpWeights[1] = (1. - xeps) * yeps;
            interpWeights[2] = xeps        * (1. - yeps);
            interpWeights[3] = xeps        * yeps;
            double  totWeight = 0;
            short   nbValidCorners = 0;
            for (short c = 0 ; c < 4 ; c++) {
               size_t ind    = interpCorners[c];
               double weight = interpWeights[c];
               if (driftX_fprod[ind] != fillvalf_fprod) {
                  interpValid[c] = 1;
                  totWeight += weight;
                  nbValidCorners++;
               } else {
                  interpValid[c] = 0;
               }
            }
            if (totWeight == 0) {
#ifdef MATCHUP_VERBOSE
               printf("All neighbours are missing (totWeight)\n");
#endif
               goto next_val_record;
            }
            if (use_nearest_neighbour && (nbValidCorners != 4)) {
#ifdef MATCHUP_VERBOSE
               printf("Found only %d valid corners in NN collocation: no valid match.\n",nbValidCorners);
#endif
               goto next_val_record;
            }

            short closestValidCorner = -1;
            float dist2closestValidCorner = bestDistInit;
            for (short c = 0 ; c < 4 ; c++) {
               double dist;
               size_t ind    = interpCorners[c];
               compute_distance(lat_start_fprod[ind],lon_start_fprod[ind],o->lat,o->lon,&dist);
               interpDist[c] = dist;

               if (!interpValid[c]) continue;
               if (dist < dist2closestValidCorner) {
                  dist2closestValidCorner = dist;
                  closestValidCorner      = c;
               }
            }
            if (closestValidCorner == -1) {
#ifdef MATCHUP_VERBOSE
               printf("All neighbours are not valid (minimum distance)\n");
#endif
               goto next_val_record;
            }


            if (use_nearest_neighbour) {
               /* nearest neighbour collocation */
               if (dist2closestValidCorner < nearest_neighbour_radius) {
                  size_t ind    = interpCorners[closestValidCorner];
                  t0Interp      = t0_fprod[ind];
                  t1Interp      = t1_fprod[ind];
#ifdef MATCHUP_VERBOSE
                  printf("Nearest Neighbour : (%.3f %.3f)<%.2f km>[%d]\n",lat_start_fprod[ind],lon_start_fprod[ind],dist2closestValidCorner,
                        closestValidCorner);
#endif
               } else {
#ifdef MATCHUP_VERBOSE
                  printf("Using a nearest neighbour approach but closest product grid is at %.3f km\n",dist2closestValidCorner);
#endif
                  goto next_val_record;
               }
            } else {
               /* bi-linear interpolation */
               if (totWeight > minimum_bilinear_weight) {
                  for (short c = 0 ; c < 4 ; c++) {
                     size_t ind    = interpCorners[c];
                     double weight = interpWeights[c];
                     if (interpValid[c]) {
                        t0Interp  += weight * t0_fprod[ind];
                        t1Interp  += weight * t1_fprod[ind];
                     }
                  }
                  t0Interp /= totWeight;
                  t1Interp /= totWeight;
               } else {
#ifdef MATCHUP_VERBOSE
                  printf("Using a bi-linear approach but totWeights is %.3f\n",totWeight);
#endif
                  goto next_val_record;
               }
            }


#ifdef MATCHUP_VERBOSE
            fmsec19702CFepoch((fmsec1970)round(t0Interp),tmpdatestr0);
            fmsec19702CFepoch((fmsec1970)round(t1Interp),tmpdatestr1);
            printf("Associated product times are T0: <%s>(%.1f) and T1: <%s>(%.1f)\n",tmpdatestr0,t0Interp,tmpdatestr1,t1Interp);
#endif

            /* if the -M flag was set, we force the in-situ records to be found close to 12UTC (central time of the product) */
            if (temporal_collocation_scheme == USE2DCOL) {
               t0Interp = refT0;
               t1Interp = refT1;
#ifdef MATCHUP_VERBOSE
               fmsec19702CFepoch((fmsec1970)round(t0Interp),tmpdatestr0);
               fmsec19702CFepoch((fmsec1970)round(t1Interp),tmpdatestr1);
               printf("2D collocation. Product times are forced to T0: <%s>(%.1f) and T1: <%s>(%.1f)\n",tmpdatestr0,t0Interp,tmpdatestr1,t1Interp);
#endif
            } else if (temporal_collocation_scheme == USE3DCOL) {
#ifdef MATCHUP_VERBOSE
               printf("3D collocation. Use product time ASIS.\n");
#endif
            } else {
               double time_delay = temporal_collocation_scheme * 60. * 60.;
               t0Interp += time_delay;
               t1Interp += time_delay;
#ifdef MATCHUP_VERBOSE
               fmsec19702CFepoch((fmsec1970)round(t0Interp),tmpdatestr0);
               fmsec19702CFepoch((fmsec1970)round(t1Interp),tmpdatestr1);
               printf("3D collocation with time delay (%02d hours). Product times are forced to T0: <%s>(%.1f) and T1: <%s>(%.1f)\n",
                     temporal_collocation_scheme,
                     tmpdatestr0,t0Interp,tmpdatestr1,t1Interp);
#endif
            } 

            /* save two best records in the trajectory : the closest to t0 and t1 */
            float DistT0 = insitu_time-t0Interp;
            if (fabs(DistT0) < fabs(bestDistT0)) {
#ifdef MATCHUP_VERBOSE
               printf("Best T0 so far (Time difference is %.2f hours)\n",DistT0/60./60.);
#endif
               bestDistT0 = DistT0;
               bestValRecT0  = rec;
               bestClosestCornerT0 = closestValidCorner;
               Duration_Prd = t1Interp - t0Interp; /* duration (time-span) of the sat drift vector */
               for (short c = 0 ; c < 4 ; c ++) {
                  bestCornersT0[c] = interpCorners[c];
                  bestWeightsT0[c] = interpWeights[c];
               }
            } else if (DistT0 > 0) {
#ifdef MATCHUP_VERBOSE
               printf("Insitu Time is later than product time (d: %.2f hours). Can stop searching now.\n",DistT0/60./60.);
#endif
               recc--;
               break;
            }
         } 
next_val_record:
         startRecord = startRecord->n;
      } while (startRecord != records->head);
      /* check if there was good enough spatial collocation */
      if ( (bestDistT0 == bestDistInit) ) {
#ifdef MATCHUP_VERBOSE
         fprintf(stdout,"WARNING (%s) No spatial match between <%s> and the product.\n",__func__,trajs[s].id);
#endif
         continue;
      }
      /* check that we have a good enough temporal match for T0 */
      if ( (fabs(bestDistT0) > TimeLimitSeconds) ) {
#ifdef MATCHUP_VERBOSE
         fprintf(stdout,"WARNING (%s) No temporal match at T0 (diff to t0: %.2f hours) for object <%s>\n",
               __func__,bestDistT0/60./60.,trajs[s].id);
#endif
         continue;
      }

      buoyRecord *bestValRecT1;
      /* if the validation source is a GCP no need to test for duration, we have a match */
      if ( strstr("gcp",trajs[s].source) ) {
#ifdef MATCHUP_VERBOSE
          printf("This is a GCP point: we are done.\n");
#endif
          bestValRecT1 = bestValRecT0;
          goto good_matchup;
      }

      /* Continue browsing the list for finding the best T1 (no spatial matching this time) */
#ifdef MATCHUP_VERBOSE
      printf("Now collocate on drift duration (Prod = %.2f)\n",Duration_Prd/60./60.);
#endif

      float bestDuration =  bestDistInit;
      do {
         recc++;
         buoyRecord *rec = startRecord->c;
         Obsdata *o = &(rec->data);

         double insitu_time = (double)rec->time;

#ifdef MATCHUP_VERBOSE
         fmsec19702CFepoch((fmsec1970)round(insitu_time),tmpdatestrIS);
         printf ("%d In situ Time: %s(%.1f) Pos: (%.3f %.3f)\n",recc,tmpdatestrIS,insitu_time,o->lat,o->lon);
#endif

         float Duration_Ref = insitu_time - bestValRecT0->time;
         if (fabs(Duration_Ref - Duration_Prd) < fabs(bestDuration - Duration_Prd)) {
#ifdef MATCHUP_VERBOSE
            printf("Best Duration so far (Duration difference is %.2f hours)\n",(Duration_Ref-Duration_Prd)/60./60.);
#endif
            bestDuration = Duration_Ref;
            bestValRecT1 = rec;
         } else if (Duration_Ref > Duration_Prd) {
#ifdef MATCHUP_VERBOSE
            printf("Insitu Duration (%.2f hours) is later than product duration (%.2f hours). Can stop searching now.\n",Duration_Ref/60./60.,Duration_Prd/60./60.);
#endif
            break;
         }
         startRecord = startRecord->n;
      } while (startRecord != records->head);
      /* check if there was good enough temporal collocation for the duration of the drift */
      if ( fabs(bestDuration - Duration_Prd) > DurationLimitSeconds ) {
#ifdef MATCHUP_VERBOSE
         fprintf(stdout,"WARNING (%s) No temporal match for Duration = T1-T0 (diff on Delta = %.2f hours) for object <%s>\n",
               __func__,(bestDuration-Duration_Prd)/60./60.,trajs[s].id);
#endif
         continue;
      }

good_matchup:
      ;

      /* =================================================================== */
      /* If we arrive here, then we have a valid matchup and we can store it */
      /* =================================================================== */

      matchProdVal oneMatchup;
      /* VAL 
       * this is easy: we just have to get to store the two position record objects. */
      sprintf(oneMatchup.id,"%s",trajs[s].id);
      sprintf(oneMatchup.network,"%s",trajs[s].network);
      sprintf(oneMatchup.source,"%s",trajs[s].source);
      memcpy(&(oneMatchup.val[0]),bestValRecT0,sizeof(*bestValRecT0));
      memcpy(&(oneMatchup.val[1]),bestValRecT1,sizeof(*bestValRecT1));
      /* PROD 
       * more difficult: we must do the collocation and then store two position record objects */
      double totWeight     = 0.;
      double totUWeight    = 0.;
      float  xpos_prod     = 0.;
      float  ypos_prod     = 0.;
      float  dX_prod       = 0.;
      float  dY_prod       = 0.;
      float  sigdX_prod    = 0.;
      float  sigdY_prod    = 0.;
      float  corrdXdY_prod = 0.;

      /* implement the nearest neighbour collocation by modifying the weights of the bilinear formula */
      if (use_nearest_neighbour) {
#ifdef MATCHUP_VERBOSE
         printf("Nearest neighbour collocation\n");
#endif
         for (size_t c = 0 ; c < 4 ; c++) {
            if (c == bestClosestCornerT0) {
               bestWeightsT0[c] = 1;
            } else {
               bestWeightsT0[c] = 0;
            }
         }
      }

      double t0_double = 0.;
      double t1_double = 0.;
#ifdef MATCHUP_VERBOSE
      printf("Corners in product grid are:\n");
#endif
      for (size_t c = 0 ; c < 4 ; c++) {
         size_t ind    = bestCornersT0[c];
         double weight = bestWeightsT0[c];

#ifdef MATCHUP_VERBOSE
         printf("\t%u %u %f %s\n",c,ind,weight,(driftX_fprod[ind] != fillvalf_fprod?"valid":"missing"));
#endif
         long i,j;
         fmijmap(ind,prodDims[XDIM],&j,&i);
         xpos_prod += weight * i ;
         ypos_prod += weight * j ;

         if (driftX_fprod[ind] != fillvalf_fprod) {
            t0_double  += weight * t0_fprod[ind];
            t1_double  += weight * t1_fprod[ind];
            dX_prod    += weight * driftX_fprod[ind];
            dY_prod    += weight * driftY_fprod[ind];
            if ((sigdX_fprod) && (sigdX_fprod[ind] != fillvalf_fprod)) {
               sigdX_prod    += weight * sigdX_fprod[ind];
               sigdY_prod    += weight * sigdY_fprod[ind];
               corrdXdY_prod += weight * corrdXdY_fprod[ind];
               totUWeight += weight;
            }
            totWeight += weight;
         }
      }
      xpos_prod            /= totWeight;
      ypos_prod            /= totWeight;
      t0_double            /= totWeight;
      t1_double            /= totWeight;
      dX_prod              /= totWeight;
      dY_prod              /= totWeight;
      fmsec1970 t0_prod     = (fmsec1970)(round(t0_double));
      fmsec1970 t1_prod     = (fmsec1970)(round(t1_double));
      if (totUWeight >= minimum_bilinear_weight) {
         sigdX_prod    /= totUWeight;
         sigdY_prod    /= totUWeight;
         corrdXdY_prod /= totUWeight;
      } else {
         sigdX_prod     = -1.;
         sigdY_prod     = -1.;
         corrdXdY_prod  = -1.;
      }
      /* get lat_b,lon_b and lat_e,lon_e out of the x,y positions */
      double latB_prod,lonB_prod;
      remap_xy2ll(xpos_prod,ypos_prod,proj,A[XDIM],B[XDIM],A[YDIM],B[YDIM],&latB_prod,&lonB_prod);
      double latE_prod,lonE_prod;
      remap_xy2ll(xpos_prod+dX_prod/A[XDIM],ypos_prod+dY_prod/A[YDIM],proj,A[XDIM],B[XDIM],A[YDIM],B[YDIM],&latE_prod,&lonE_prod);

      /* store all in the matchup object */
      oneMatchup.latB_prod = latB_prod;
      oneMatchup.lonB_prod = lonB_prod;
      oneMatchup.latE_prod = latE_prod;
      oneMatchup.lonE_prod = lonE_prod;
      oneMatchup.tB_prod   = t0_prod;
      oneMatchup.tE_prod   = t1_prod;
      long e = bestCornersT0[bestClosestCornerT0];
      long i,j;
      fmijmap(e,prodDims[XDIM],&j,&i);
      oneMatchup.i_prod    = i;
      oneMatchup.j_prod    = j;
      oneMatchup.sX_prod   = sigdX_prod;
      oneMatchup.sY_prod   = sigdY_prod;
      oneMatchup.cXY_prod  = corrdXdY_prod;
      oneMatchup.flag_prod = flags_fprod[e];
      /* compute the collocation distance at T0 and time differences */
      double coll_dist;
      compute_distance(oneMatchup.val[0].data.lat,oneMatchup.val[0].data.lon,oneMatchup.latB_prod,oneMatchup.lonB_prod,&coll_dist);
      oneMatchup.distance  = coll_dist;
      oneMatchup.dtB  = (oneMatchup.val[0].time - oneMatchup.tB_prod)/60./60.;
      oneMatchup.dtE  = (oneMatchup.val[1].time - oneMatchup.tE_prod)/60./60.;
     
#ifdef MATCHUP_VERBOSE
      fmsec19702CFepoch(t0_prod,tmpdatestr0);
      fmsec19702CFepoch(t1_prod,tmpdatestr1);
      printf("Product drift %s (%.3f %.3f) and %s (%.3f %.3f)\n",
            tmpdatestr0,latB_prod,lonB_prod,tmpdatestr1,latE_prod,lonE_prod);
      fmsec19702CFepoch(oneMatchup.val[0].time,tmpdatestr0);
      fmsec19702CFepoch(oneMatchup.val[1].time,tmpdatestr1);
      printf("Buoy drift    %s (%.3f %.3f) and %s (%.3f %.3f)\n",
            tmpdatestr0,oneMatchup.val[0].data.lat,oneMatchup.val[0].data.lon,tmpdatestr1,oneMatchup.val[1].data.lat,oneMatchup.val[1].data.lon);
#endif

      /* store the matcup object in the list of all matchups */
      ret = dbl_list_addNode(allMatchups,dbl_node_createCopyContent(allMatchups,&oneMatchup,NULL),LIST_POSITION_LAST);

      nb_valid_matchups++;
   }

   return 0;
}
