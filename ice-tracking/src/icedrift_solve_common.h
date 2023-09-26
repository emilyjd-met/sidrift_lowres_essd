
#ifndef ICEDRIFT_SOLVE_COMMON_H
#define ICEDRIFT_SOLVE_COMMON_H


int model_wrapper_CC(size_t nx,double *y,double *res);
int model_wrapper_SD(size_t nx,double *y,double *res);
int model_wrapper_SLSD(size_t nx,double *y,double *res);
int choose_model(char *optimMetric,char *opti,int (**func)(size_t,double *,double *));
int sort_icedrift_bylength(int nb, size_t indexes[], float leng[]);



#endif /* ICEDRIFT_SOLVE_COMMON_H */
