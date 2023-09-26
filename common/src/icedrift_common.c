
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <fmutil.h>
#include <projects.h>
#include <fmutil.h>
#include "icedrift_common.h"

void compute_distance(double lat0, double lon0, double lat1, double lon1, double *dist) {

   lat0 *= DEG_TO_RAD;
   lon0 *= DEG_TO_RAD;
   lat1 *= DEG_TO_RAD;
   lon1 *= DEG_TO_RAD;

   double dlat = 0.5*(lat1-lat0);
   double dlon = 0.5*(lon1-lon0);
	   
   double a = sin(dlat)*sin(dlat) + cos(lat0)*cos(lat1)*sin(dlon)*sin(dlon);
   double c = 2.*atan2(sqrt(a),sqrt(1-a));
   *dist    = REARTH * c;

}

void compute_directionToNorth(double lat0, double lon0, double lat1, double lon1, double *dtn) {

   lat0 *= DEG_TO_RAD;
   lon0 *= DEG_TO_RAD;
   lat1 *= DEG_TO_RAD;
   lon1 *= DEG_TO_RAD;

   double dlat = 0.5*(lat1-lat0);
   double dlon = (lon1-lon0);
	    
   /* retrieves the initial bearing when travelling from lon0,lat0 to lon1,lat1 */
   /* atan2 returns angles between -pi and pi. */
   double direction = atan2( sin(dlon)*cos(lat1) ,
	       cos(lat0)*sin(lat1) - sin(lat0)*cos(lat1)*cos(dlon) );
   /* 0 is north. pi (-pi) is south. */
   direction  = fmPI + direction;
   *dtn      = direction * RAD_TO_DEG;

}
 
/* remap with pj (PROJ4 object) from lat,lon to x,y in a grid defined by Ax,Ay,Bx,By */
/* WARNING: no checks on pj == NULL */
int remap_ll2xy(double lat, double lon, PJ *pj, double Ax, double Bx, double Ay, double By, double *x, double *y, short round) {
   
   projUV G;

   G.u = lon * DEG_TO_RAD;
   G.v = lat * DEG_TO_RAD;
   G   = pj_fwd(G, pj);
   G.u /= 1000.;
   G.v /= 1000.;
	
   *x  =  (G.u - Bx)/Ax;
   *y  = -(G.v - By)/Ay;
   if (round) {
      *x = rint(*x);
      *y = rint(*y);
   }

   return 0;
}

/* remap with pj (PROJ4 object) from x,y to lat,lon [deg] in a grid defined by Ax,Ay,Bx,By */
/* WARNING: no checks on pj == NULL */
int remap_xy2ll(double x,double y, PJ *pj, double Ax, double Bx, double Ay, double By, double *lat, double *lon) {
   
   projUV G;

   G.u  = ( +x * Ax + Bx ) * 1000.;
   G.v  = ( -y * Ay + By ) * 1000.;
   G    = pj_inv(G, pj);
   *lon = G.u * RAD_TO_DEG;
   *lat = G.v * RAD_TO_DEG;

   return 0;
}

/* remap from x,y in pj1 (Ax1,Bx1,Ay1,By1) to x,y in pj2 (Ax2,Bx2,Ay2,By2) */
/* WARNING: no checks on pj == NULL */
int remap_xy2xy(double *x,double *y,short round, 
      PJ *pj1, double Ax1, double Bx1, double Ay1, double By1,
      PJ *pj2, double Ax2, double Bx2, double Ay2, double By2) {
   
   int ret;
   double lat,lon;

   //printf("remap_xy2xy: PJ1 is (%f,%f,%f,%f)\n",Ax1,Bx1,Ay1,By1);
   //printf("remap_xy2xy: PJ2 is (%f,%f,%f,%f)\n",Ax2,Bx2,Ay2,By2);

   //printf("remap_xy2xy: in PJ1, x=%f, y=%f | ",*x,*y);
   
   ret = remap_xy2ll(*x,*y,pj1,Ax1,Bx1,Ay1,By1,&lat,&lon);
   //printf("lat =%f, lon=%f | ",lat,lon);
   ret = remap_ll2xy(lat,lon,pj2,Ax2,Bx2,Ay2,By2,x,y,round);
   
   //printf("in PJ2, x=%f, y=%f\n",*x,*y);

   return 0;
}

/* pointIsInGrid: xf,yf in _pixel_ fractions, grid_dims in pixels [0..maxX][0..maxY] */
int pointIsInGrid(double xf, double yf, size_t grid_dims[]) {

   int ret = 1;
   if ( xf < 0 )
      ret = 0;
   else if ( yf < 0 )
      ret = 0;
   else if ( xf >= grid_dims[XDIM] )
      ret = 0;
   else if ( yf >= grid_dims[YDIM] )
      ret = 0;

   return (ret);
}

/* pointIsInGridExtent: xf, yf in grid's axis unit (e.g. km) and the same for grid_extent.
 *    margin must be >= 0, in same unit as the others.
 *    grid_extent[4] is [x_low,x_high,y_low,y_high] */
int pointIsInGridExtent(double xf, double yf, double grid_extent[/*4*/], double margin) {

   assert(margin >= 0);

   int ret = 1;
   if (      xf < (grid_extent[0]-margin) )
      ret = 0;
   else if ( yf < (grid_extent[2]-margin) )
      ret = 0;
   else if ( xf > (grid_extent[1]+margin) )
      ret = 0;
   else if ( yf > (grid_extent[3]+margin) )
      ret = 0;

   return (ret);
}

int pointIsOnIce(long pix, short *icelandmask) {

   return (icelandmask[pix] == 3);

}

int pointIsOnLand(long pix, short *icelandmask) {

   return (icelandmask[pix] >= 9); // 9: land, 10: coast

}

void transform_IceDrift_fromXY(double lat0, double lon0, double xdrift, double ydrift,
      float Ax,float Bx, float Ay, float By, PJ *proj,
      double *plat1, double *plon1, double *plength, double *pdir, double *pdriftNS, double *pdriftEW) {

   /* compute the grid coordinate of start point */ 
   double x0,y0;
   remap_ll2xy(lat0,lon0,proj,Ax,Bx,Ay,By,&x0,&y0,0);

   /* deduce grid coordinate of final point */
   double x1,y1;
   x1 = x0 + xdrift/Ax;
   y1 = y0 + ydrift/Ay;

   /* compute the latitude and longitude of final point */ 
   double lat1,lon1;
   remap_xy2ll(x1,y1,proj,Ax,Bx,Ay,By,&lat1,&lon1);

   /*
   printf("ICEDRIFT is x:%f, y:%f. (%f,%f)(%f,%f) -> (%f,%f)(%f %f)\n",
	 xdrift,ydrift,
	 x0,y0,lat0,lon0,
	 x1,y1,lat1,lon1);
	 */

   /* compute the length and direction of the arrows */ 
   double dlen,ddir;
   compute_distance(lat0,lon0,lat1,lon1,&dlen);
   compute_directionToNorth(lat0,lon0,lat1,lon1,&ddir);

   /* compute the East-West and North-South components */
   double driftNS = dlen*cos(ddir*DEG_TO_RAD);
   double driftEW = dlen*sin(ddir*DEG_TO_RAD);

   // printf("LENGTH: %fkm, DIR: %fdeg, driftNS: %fkm, driftEW: %fkm\n",dlen,ddir,driftNS,driftEW);

   *plat1    = lat1;
   *plon1    = lon1;
   *plength  = dlen;
   *pdir     = ddir;
   *pdriftNS = driftNS;
   *pdriftEW = driftEW;

}

void transform_IceDrift_fromZA(double lat0, double lon0, double Zdrift, double Adrift,
      float Ax,float Bx, float Ay, float By, PJ *proj,
      double *plat1, double *plon1, double *plength, double *pdir, double *pxdrift, double *pydrift) {

   lat0 *= DEG_TO_RAD;
   lon0 *= DEG_TO_RAD;

   /* compute length and direction of the drift vector */
   double dlen,ddir;
   dlen = sqrt ( Zdrift*Zdrift + Adrift*Adrift );
   ddir = atan2 ( Adrift, Zdrift );
   ddir+= fmPI;

   /* compute lat1 and lon1 */
   double dlat1,dlon1;
   dlat1 = asin(sin(lat0)*cos(dlen/REARTH) + cos(lat0)*sin(dlen/REARTH)*cos(ddir));
   dlon1 = lon0 + atan2(sin(ddir)*sin(dlen/REARTH)*cos(lat0), cos(dlen/REARTH) - sin(lat0)*sin(dlat1));

   /* remap the start and end positions */
   double x0,y0;
   remap_ll2xy(lat0,lon0,proj,Ax,Bx,Ay,By,&x0,&y0,0);
   double x1,y1;
   remap_ll2xy(dlat1,dlon1,proj,Ax,Bx,Ay,By,&x1,&y1,0);

   /* compute xdrift and ydrift */
   double dxdrift,dydrift;
   dxdrift = (x1-x0) * Ax;
   dydrift = (y1-y0) * Ay;

   /* tranform all (output) angular values back to degrees */ 
   ddir  *= RAD_TO_DEG;
   dlat1 *= RAD_TO_DEG;
   dlon1 *= RAD_TO_DEG;

   printf("LENGTH: %fkm, DIR: %fdeg, driftX: %fkm, driftY: %fkm\n",dlen,ddir,dxdrift,dydrift);

   *plat1   = dlat1;
   *plon1   = dlon1;
   *plength = dlen;
   *pdir    = ddir;
   *pxdrift = dxdrift;
   *pydrift = dydrift;

}
