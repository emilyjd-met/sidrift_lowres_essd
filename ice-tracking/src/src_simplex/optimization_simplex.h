#ifndef MINIMISATION_SIMPLEX_H
#define MINIMISATION_SIMPLEX_H

#define FTOL 1e-9

#define FOUND_MINIMUM     0
#define ERROR             1
#define REACHED_ITERA_MAX 2

int findBestX(size_t ndim, double *x, double *score, double *startingPoints[],
      char opti,int (*func)(size_t nx,double *x,double *f));


#endif
