#ifndef ICEDRIFT_SOLVE_CORE_H
#define ICEDRIFT_SOLVE_CORE_H

/* global variables passed from CFFI */
extern double **obs[2];
extern short **TCflag[2];
extern short *icelandmask[2];
extern double *img_lat;
extern double *img_lon;
extern double Ax;
extern double Bx;
extern double Ay;
extern double By;
extern size_t img_dims[3];
extern char *img_projstr;
extern short nbWaveBands;
extern short twghtStart;
extern short twghtEnd;
extern double maxdriftdistance;
extern double sigmoid_length;
extern char *OptimMetric;
extern double pattern_radius[NBPATTERNS];
extern double radiusNeighbours;
extern int NDRIFTPIXELS;
extern char *out_area;
extern char *out_projstr;
extern double out_Ax;
extern double out_Bx;
extern double out_Ay;
extern double out_By;
extern size_t out_dims[3];
extern double *olat;
extern double *olon;
extern unsigned int *owcs;
extern unsigned int *iwcs;
extern float *driftX;
extern float *driftY;
extern long *t0; 
extern long *t1;
extern short *pflag;
extern float *sigdX;
extern float *sigdY;
extern float *corrdXdY;
extern short *uflag;
extern float *length;
extern float *dir;
extern float *outlonB;
extern float *outlatB;
extern float *outlonE;
extern float *outlatE;
extern float *outfc;
extern short *snavg;
extern float *avgX;
extern float *avgY;
extern float *length_avg;
extern float *length_diff;
extern float *stdX;
extern float *stdY;
extern short *patternIndex;
extern char *reportFile;

int core (void);

#endif 
