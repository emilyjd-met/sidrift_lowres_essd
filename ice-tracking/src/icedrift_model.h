
#ifndef ICEDRIFT_MODEL_H
#define ICEDRIFT_MODEL_H


/* UNKNOWNS FOR EACH DRIFT LOCATIONS ON THE GRID */
/* The drift estimate at each drift-grid location is made of 4 parameters: 
 * - (x,y) is the drift vector. Both x and y have unit km and are relative to 
 *         the x and y dimensions of the satellite image 
 * - (a,b) are for adjusting the intensity of the images.
 *
 * These 2 pairs might have regularity constraints from one drift-grid location
 * to the next but for different reasons.
 * - (x,y) because of the physical properties of the drift in the compact ice-pack.
 * - (a,b) because the main source of intensity variation in time are the atmospheric
 *         conditions (Temperature, Water Vapour, etc.) which should be spatially smooth.
 */
#define XDRIFT_IDX   0
#define YDRIFT_IDX   1
//#define ALPHA_IDX    2 
//#define BETA_IDX     3
#define NUNKNOWNS    2 /* x,y,alpha,beta */


/* NUMBER OF GRID LOCATIONS WHERE WE WANT AN ESTIMATE OF THE DRIFT */
extern int   NDRIFTPIXELS;
extern short nbWaveBands;
#define NBPATTERNS      2
#define MAX_NBWAVEBANDS 6

/* Global variables to communicate between the various routines and the cost function */
extern double *xpr;
extern double *xdif;
extern double *sxpr;
extern double **cixpr;
extern unsigned int *owcs;
extern unsigned int *iwcs;
extern double   *olon;
extern double   *olat;
extern size_t   out_dims[3];

extern long     center_coord[3];

extern PJ        *img_pj;
extern size_t    img_dims[3];
extern double    *img_lat,*img_lon;
extern double    img_Ax;
extern double    img_Bx;
extern double    img_Ay;
extern double    img_By;
extern double    **obs[2];
extern short     *icelandmask[2];
extern double    *ocalc;
extern long      *ocalc_wind;
extern double    *odif;
extern double    *ciobs_diag;
extern short     **TCflag[2];
extern short     *ocalcIsValid;
extern long      timeref[2];

extern size_t     pattern_center[NBPATTERNS];
extern size_t     pattern_size[NBPATTERNS];
extern short      *pattern_mask[NBPATTERNS];
extern short      *pattern_icemask[NBPATTERNS];
extern short      *pattern_icemask2[NBPATTERNS];
extern long       *pattern_windexes[NBPATTERNS][3];
extern double     **pattern_img[NBPATTERNS];
extern short      **pattern_isvalid[NBPATTERNS];
extern double     **pattern_img2[NBPATTERNS];
extern short      **pattern_isvalid2[NBPATTERNS];
double *directions;

extern double constraint_multiplier; 
extern double constraint_power;
extern double sigmoid_length;
extern double maxdriftdistance;
extern double bestKnowledge_lat;
extern double bestKnowledge_lon;

extern double pattern_Alpha[MAX_NBWAVEBANDS];
extern double pattern_Beta[MAX_NBWAVEBANDS];

void compute_shifted_image(int *ret,double xdrift, double ydrift,
		  double **image,short **TCflag, short *icelandmask, size_t img_dims[3],
		  long coord[3],
		  size_t pattern_size, short *pattern_mask, long *pattern_windexes[3], 
		  double **pattern_img, short **pattern_isvalid, short *pattern_icemask,
		  double *olat, double *olon);

void load_subimage(int *ret, double **image, short **TCflag, short *icelandmask, size_t img_dims[3],
      long img_coords[3], size_t pattern_size, 
      short *pattern_mask, long *pattern_windexes[3], double **pattern_img, short **pattern_isvalid, short *pattern_icemask);


void locfmivec(long *e, long x, long y, unsigned long nx);
void locfmijmap(long elem, unsigned long nx, long *x, long *y);


/* prototype for cost function */
void model_onevector(int *n, double *unknown,double *fc,long e);
void model_onevector_corr(int *n, double *unknown,double *fc, long e, short pat);
void model_onevector_sqdiff_linearScaled(int *n, double *unknown,double *fc, long e, short pat);

#endif /* ICEDRIFT_MODEL_H */
