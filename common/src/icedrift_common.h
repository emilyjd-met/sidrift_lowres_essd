
#ifndef ICEDRIFT_COMMON_H
#define ICEDRIFT_COMMON_H


#ifndef PI
#define PI 3.14159265358979323846
#endif

/* Time 'named' dimensions */
#define BEG  0
#define END  1

/* Image and product 'named' dimensions */
#define XDIM 0
#define YDIM 1
#define TDIM 2

#define REARTH 6371

#define myfmijmap(e,nx,x,y) fmijmap(e,nx,y,x)
#define mmax(a,b) (a>b?a:b)

extern char progname[];

void convert_XYdrift_into_LengthDir (const double xpos, const double ypos, double xdrift, double ydrift, 
      const double NorthPole_x, const double NorthPole_y, double *Length, double *Dir);

int remap_ll2xy(double lat, double lon, PJ *pj, double Ax, double Bx, double Ay, double By, double *x, double *y, short round);
int remap_xy2ll(double x,double y, PJ *pj, double Ax, double Bx, double Ay, double By, double *lat, double *lon);
int remap_xy2xy(double *x,double *y,short round, 
      PJ *pj1, double Ax1, double Bx1, double Ay1, double By1,
      PJ *pj2, double Ax2, double Bx2, double Ay2, double By2);
void compute_distance(double lat0, double lon0, double lat1, double lon1, double *dist);
void compute_directionToNorth(double lat0, double lon0, double lat1, double lon1, double *dtn);
int pointIsInGrid(double xf, double yf, size_t grid_dims[]);
int pointIsInGridExtent(double xf, double yf, double grid_extent[/*4*/], double margin);
int pointIsOnIce(long pix, short *icelandmask);
int pointIsOnLand(long pix, short *icelandmask);
void transform_IceDrift_fromXY(double lat0, double lon0, double xdrift, double ydrift,
      float Ax,float Bx, float Ay, float By, PJ *proj,
      double *plat1, double *plon1, double *plength, double *pdir, double *pdriftNS, double *pdriftEW);
void transform_IceDrift_fromZA(double lat0, double lon0, double Zdrift, double Adrift,
      float Ax,float Bx, float Ay, float By, PJ *proj,
      double *plat1, double *plon1, double *plength, double *pdir, double *pxdrift, double *pydrift);
#endif /* ICEDRIFT_COMMON_H */
