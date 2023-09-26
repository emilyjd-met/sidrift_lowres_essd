

#include <stdio.h>
#include <stdlib.h>
#include <fmutil.h>
#include <math.h>

#include "icedrift_solve_filter.h"

int compute_mean_vector(double *x,double *y,double *xstd, double *ystd, size_t nbvecs, double *xs, double *ys) {

   *x = *xstd = 0;
   *y = *ystd = 0;
   for (size_t n = 0 ; n < nbvecs ; n++) {
      *x += xs[n];
      *y += ys[n];
      *xstd += pow(xs[n],2.);
      *ystd += pow(ys[n],2.);
   }
   *x /= nbvecs;
   *y /= nbvecs;
   *xstd /= nbvecs;
   *ystd /= nbvecs;

   *xstd = sqrt(*xstd - pow(*x,2.));
   *ystd = sqrt(*ystd - pow(*y,2.));

   return 0;

}
