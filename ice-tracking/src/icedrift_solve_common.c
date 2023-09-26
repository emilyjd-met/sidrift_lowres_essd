
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fmutil.h>
#include <projects.h>
#include "icedrift_common.h"
#include "icedrift_model.h"
#include "icedrift_prepost.h"

/* global variables from icedrift_solve_{meth}.c */
extern int    xdim;
extern double *x;
extern size_t currentVector;
extern short  currentPattern;

int model_wrapper_CC(size_t nx,double *y,double *res) {
   
   size_t loc = currentVector*NUNKNOWNS;
   x[loc+XDRIFT_IDX] = y[0];
   x[loc+YDRIFT_IDX] = y[1];
   /*
   x[loc+ALPHA_IDX]  = 1;
   x[loc+BETA_IDX]   = 0;
   */

   *res = 0;
      
   model_onevector_corr(&xdim,x,res,(long)currentVector,currentPattern);

   return 0;
}

int model_wrapper_SD(size_t nx,double *y,double *res) {
   
   size_t loc = currentVector*NUNKNOWNS;
   x[loc+XDRIFT_IDX] = y[0];
   x[loc+YDRIFT_IDX] = y[1];
   /*
   x[loc+ALPHA_IDX]  = 1;
   x[loc+BETA_IDX]   = y[3];
   */

   *res = 0;
   model_onevector(&xdim,x,res,(long)currentVector);

   return 0;
}

/* squared linearly scaled difference */
int model_wrapper_SLSD(size_t nx,double *y,double *res) {
   
   size_t loc = currentVector*NUNKNOWNS;
   x[loc+XDRIFT_IDX] = y[0];
   x[loc+YDRIFT_IDX] = y[1];

   *res = 0;
   model_onevector_sqdiff_linearScaled(&xdim,x,res,(long)currentVector,currentPattern);

   return 0;
}

int choose_model(char *optimMetric,char *opti,int (**func)(size_t,double *,double *)) {

   if (!strcmp(optimMetric,"SD")) {
      /* minimixe the square difference */
      *opti = '<';
      *func = &model_wrapper_SD;
   } else if (!strcmp(optimMetric,"CC")) {
      /* maximize the cross correlation */
      *opti = '>';
      *func = &model_wrapper_CC;
   } else {
      fprintf(stderr,"ERROR (%s) Unknown optimization metric (%s). Can be %s or %s\n",__func__,optimMetric,"SD","CC");
      return 1;
   }

   return 0;

}

struct pair { size_t ind; float len; };

int byLength(struct pair *p1, struct pair *p2) {
   return (p1->len < p2->len);
}

int sort_icedrift_bylength(int nb, size_t indexes[], float leng[]) {

   struct pair *pairs = fmMalloc(nb * sizeof(struct pair));
   for (size_t p = 0 ; p < nb ; p++) {
      pairs[p].ind = indexes[p];
      pairs[p].len = leng[pairs[p].ind];
   }
   qsort(pairs,nb,sizeof(struct pair),&byLength);
   for (size_t p = 0 ; p < nb ; p++) {
      indexes[p] = pairs[p].ind;
   }

   return 0;

}

