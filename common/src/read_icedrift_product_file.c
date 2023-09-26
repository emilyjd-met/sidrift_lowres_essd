
/*
 * PURPOSE:
 * To read from ice-drift satellite products as processed by the OSISAF chain.
 * 
 * BUGS:
 * None known.
 * 
 * AUTHOR: 
 * Thomas Lavergne, 10.10.2008 (extracted from matchup_validation_with_product.c)
 *
 * MODIFIED:
 * Thomas Lavergne, 10.10.2008 Add reading the processing flags
 * Thomas Lavergne, 10.03.2009 Add reading the uncertainty estimates (optional)
 * Thomas Lavergne, 31.03.2010 Adapted to also read the products in 'final' format
 *
 * CVS:
 * $Id$
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <fmutil.h>
#include <netcdf.h>
#include <projects.h>
#include <rsprod.h>
#define  USE_CF_ROUTINES
#include <osisaf_ice_netcdf.h>
#include <icedrift_common.h>
#include <icedrift_finalNCformat_conventions.h>
#include <icedrift_filenames.h>
#include "read_icedrift_product_file.h"

extern char sigdXname[10];
extern char sigdYname[10];

int readProductFile(char infname[],     /* full name and path to icedrift product file (IN)*/ 
      fmsec1970 refT0, fmsec1970 refT1, /* time references for the product. 
					   The time information in the netcdf low res product file are 
					   only relative to the DAY0 and DAY1 so that we need this reference time.
					   Might be un-needed for AVHRR files. (IN) */
      size_t out_dims[3],               /* (OUT) grid info: float out_dims[3] with: 
					          out_dims[0] = nx = number_of_pixels_x, 
					          out_dims[1] = ny and
						  out_dims[2] = nx*ny */
      float A[2], float B[2],           /* (OUT) more grid info: 
					          A[0] = reso_x_pixel (km), A[1] = reso_y_pixel (km) 
                                                  B[0] = x_coord_corner_pixel (km), B[1] = y_coord_corner_pixel (km) */
      char *projstr, PJ **proj,         /* (OUT) projstr: PROJ4 proj string 
                                           (OUT) proj: initialized PJ object. */
      fmsec1970 **t0, fmsec1970 **t1,            /* (OUT) array of out_dims[2] elements, containing t0 and t1 */
      float **lat,float **lon,          /* (OUT) array of out_dims[2] elements, containing lat0 and lon0 */
      float **driftX, float **driftY,   /* (OUT) array of out_dims[2] elements, containing driftX and driftY 
					   (drift vector on the product grid, in km) */
      short **flags,                    /* (OUT) array of out_dims[2] elements, containing a flag information (type use/not-use)*/
      float *fillvalue,                  /* (OUT) fillvalue for the float datasets (especially driftX and driftY) */
      float **sigdX, float **sigdY, float **corrdXdY, short **uflags) {

  int ret,ncret,projret;
   int ncid;


   /* open the netcdf file */
   ret = nc_open(infname,NC_NOWRITE,&ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot nc_open() file for reading:\n\t%s\n",__func__,infname);
      return 1;
   }

   /* check if the filename corresponds to a 'intermediate' or 'final' product */
   short fmt_type = 0;
   char *endDirPart,*basename;
   fmbasename(infname,&endDirPart,&basename);
   if (strstr(basename,LRSID_FILEN_SUFFIX)) {
      /* final format */
      fmt_type = LRSID_FILEFORMAT_FINL;
   } else if (strstr(basename,"icedrift_")) {
      fmt_type = LRSID_FILEFORMAT_PROC;
   } else {
      fprintf(stderr,"ERROR (%s) Cannot guess the formatting ('proc' or 'final') for netcdf file:\n\t%s\n",__func__,infname);
      ret = nc_close(ncid);
      return 1;
   }
   char *t_name;
   char *drift_name;
   char *flag_name;
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      t_name = fmMalloc(strlen("dtx")+1);
      drift_name = fmMalloc(strlen("dX")+1);
      flag_name = fmMalloc(strlen("status_flag")+1);
   } else {
      t_name = fmMalloc(strlen("tx")+1);
      drift_name = fmMalloc(strlen("driftX")+1);
      flag_name = fmMalloc(strlen("flag")+1);
   }


   
   /* load the time bounds: */
   rsprod_field *FTbounds;
   if (rsprod_field_loadFromNetCDF(&FTbounds,"time_bnds",ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the 'time_bnds' dataset.\n",__func__);
      return 1;
   }
   double *FTbounds_data;
   if (rsprod_data_getType(FTbounds->data) != RSPROD_DOUBLE) {
      fprintf(stderr,"ERROR (%s) The 'time_bnds' dataset is not of type DOUBLE.\n",__func__);
      return 1;
   }
   if ( (*FTbounds->data->methods->accessValues)(FTbounds->data,&FTbounds_data) )  {
      fprintf(stderr,"ERROR (%s) Error in accessing the 'time_bounds' values.\n",__func__);
      return 1;
   }
   //printf ("Time bounds in the file are : %f (%ld) %f (%ld)\n",FTbounds_data[0],refT0,FTbounds_data[1],refT1);
   /* tranform those time bounds to reference 1-1-1970 */
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      fmsec1970 refTime = CFepoch2fmsec1970(LRSID_TIMEREF);
      FTbounds_data[0] += refTime;
      FTbounds_data[1] += refTime;
   }
   /* check that the value corresponds to the one give as input */
   if ((((fmsec1970)round(FTbounds_data[0])) != refT0) ||
	 (((fmsec1970)round(FTbounds_data[1])) != refT1)) {
      fprintf(stderr,"ERROR (%s) The time_bnds (time_bounds) values do not correspond those of refT0 and refT1\n",__func__);
      return 1;
   }

   /* lat */
   rsprod_field *Flat;
   if (rsprod_field_loadFromNetCDF(&Flat,"lat",ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the 'lat' dataset.\n",__func__);
      return 1;
   }
   float *Flat_data;
   if (rsprod_data_getType(Flat->data) != RSPROD_FLOAT) {
      fprintf(stderr,"ERROR (%s) The 'lat' dataset is not of type FLOAT.\n",__func__);
      return 1;
   }
   if ( (*Flat->data->methods->accessValues)(Flat->data,&Flat_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the 'lat' values.\n",__func__);
      return 1;
   }
   /* lon */
   rsprod_field *Flon;
   if (rsprod_field_loadFromNetCDF(&Flon,"lon",ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the 'lon' dataset.\n",__func__);
      return 1;
   }
   float *Flon_data;
   if (rsprod_data_getType(Flon->data) != RSPROD_FLOAT) {
      fprintf(stderr,"ERROR (%s) The 'lon' dataset is not of type FLOAT.\n",__func__);
      return 1;
   }
   if ( (*Flon->data->methods->accessValues)(Flon->data,&Flon_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the 'lon' values.\n",__func__);
      return 1;
   }
   /* t0 */
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      sprintf(t_name,"dt0");
   } else {
      sprintf(t_name,"t0");
   }
   rsprod_field *Ft0;
   if (rsprod_field_loadFromNetCDF(&Ft0,t_name,ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,t_name);
      return 1;
   }
   int *Ft0_data;
   if (rsprod_data_getType(Ft0->data) != RSPROD_INT) {
      fprintf(stderr,"ERROR (%s) The '%s' dataset is not of type INT.\n",__func__,t_name);
      return 1;
   }
   if ( (*Ft0->data->methods->accessValues)(Ft0->data,&Ft0_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,t_name);
      return 1;
   }
   /* t1 */
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      sprintf(t_name,"dt1");
   } else {
      sprintf(t_name,"t1");
   }
   rsprod_field *Ft1;
   if (rsprod_field_loadFromNetCDF(&Ft1,t_name,ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,t_name);
      return 1;
   }
   int *Ft1_data;
   if (rsprod_data_getType(Ft1->data) != RSPROD_INT) {
      fprintf(stderr,"ERROR (%s) The '%s' dataset is not of type INT.\n",__func__,t_name);
      return 1;
   }
   if ( (*Ft1->data->methods->accessValues)(Ft1->data,&Ft1_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,t_name);
      return 1;
   }

   /* driftX */
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      sprintf(drift_name,"dX");
   } else {
      sprintf(drift_name,"driftX");
   }
   rsprod_field *FdriftX;
   if (rsprod_field_loadFromNetCDF(&FdriftX,drift_name,ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,drift_name);
      return 1;
   }
   float *FdriftX_data;
   if (rsprod_data_getType(FdriftX->data) != RSPROD_FLOAT) {
      fprintf(stderr,"ERROR (%s) The '%s' dataset is not of type FLOAT.\n",__func__,drift_name);
      return 1;
   }
   if ( (*FdriftX->data->methods->accessValues)(FdriftX->data,&FdriftX_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,drift_name);
      return 1;
   }
   /* driftY */
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      sprintf(drift_name,"dY");
   } else {
      sprintf(drift_name,"driftY");
   }
   rsprod_field *FdriftY;
   if (rsprod_field_loadFromNetCDF(&FdriftY,drift_name,ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,drift_name);
      return 1;
   }
   float *FdriftY_data;
   if (rsprod_data_getType(FdriftY->data) != RSPROD_FLOAT) {
      fprintf(stderr,"ERROR (%s) The '%s' dataset is not of type FLOAT.\n",__func__,drift_name);
      return 1;
   }
   if ( (*FdriftY->data->methods->accessValues)(FdriftY->data,&FdriftY_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,drift_name);
      return 1;
   }

   /* flags */
   int flag_type;
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      sprintf(flag_name,"status_flag");
      flag_type = RSPROD_BYTE;
   } else {
      sprintf(flag_name,"flag");
      flag_type = RSPROD_SHORT;
   }
   rsprod_field *Fflags;
   if (rsprod_field_loadFromNetCDF(&Fflags,flag_name,ncid)) {
      fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,flag_name);
      return 1;
   }
   short *Fflags_data;
   if (rsprod_data_getType(Fflags->data) != flag_type) {
      fprintf(stderr,"ERROR (%s) The '%s' dataset is not of expected type (is %s).\n",
            __func__,flag_name, TypeName[rsprod_data_getType(Fflags->data)]);
      return 1;
   }
   if ( (*Fflags->data->methods->accessValues)(Fflags->data,&Fflags_data) )  {
      fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,flag_name);
      return 1;
   }
   
   /* test if we can find the uncertainties datasets */
   rsprod_field *FsigdX,*FsigdY,*FcorrdXdY,*Fuflags;

   if (strlen(sigdXname) == 0) 
      strcpy(sigdXname,"sigdX");
   if (strlen(sigdYname) == 0) 
      strcpy(sigdYname,"sigdY");
      
   fprintf(stdout,"INFO (%s) Look for uncertainties as <%s> <%s>\n",__func__,sigdXname,sigdYname);

   short haveUncertainties = 0;
   int varid;
   int ncretx = nc_inq_varid (ncid,sigdXname,&varid);
   int ncrety = nc_inq_varid (ncid,sigdYname,&varid);
   //printf("Search for dataset <%s> returned %d (%d). Varid is %d.\n",sigdXname,ncret,NC_NOERR,varid);
   if ( (ncretx == NC_NOERR) && (ncrety == NC_NOERR)) 
      haveUncertainties=1;

   if (haveUncertainties) {
      /* sigdX */
      if (rsprod_field_loadFromNetCDF(&FsigdX,sigdXname,ncid)) {
         fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,sigdXname);
         return 1;
      }
      float *FsigdX_data;
      if (rsprod_data_getType(FsigdX->data) != RSPROD_FLOAT) {
         fprintf(stderr,"ERROR (%s) The '%s' dataset is not of type FLOAT.\n",__func__,sigdXname);
         return 1;
      }
      if ( (*FsigdX->data->methods->accessValues)(FsigdX->data,&FsigdX_data) )  {
         fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,sigdXname);
         return 1;
      }
      /* sigdY */
      if (rsprod_field_loadFromNetCDF(&FsigdY,sigdYname,ncid)) {
         fprintf(stderr,"ERROR (%s) While loading the '%s' dataset.\n",__func__,sigdYname);
         return 1;
      }
      float *FsigdY_data;
      if (rsprod_data_getType(FsigdY->data) != RSPROD_FLOAT) {
         fprintf(stderr,"ERROR (%s) The '%s' dataset is not of type FLOAT.\n",__func__,sigdYname);
         return 1;
      }
      if ( (*FsigdY->data->methods->accessValues)(FsigdY->data,&FsigdY_data) )  {
         fprintf(stderr,"ERROR (%s) error in accessing the '%s' values.\n",__func__,sigdYname);
         return 1;
      }
      /* corrdXdY */
      if (rsprod_field_loadFromNetCDF(&FcorrdXdY,"corrdXdY",ncid)) {
         fprintf(stderr,"ERROR (%s) While loading the 'corrdXdY' dataset.\n",__func__);
         return 1;
      }
      float *FcorrdXdY_data;
      if (rsprod_data_getType(FcorrdXdY->data) != RSPROD_FLOAT) {
         fprintf(stderr,"ERROR (%s) The 'corrdXdY' dataset is not of type FLOAT.\n",__func__);
         return 1;
      }
      if ( (*FcorrdXdY->data->methods->accessValues)(FcorrdXdY->data,&FcorrdXdY_data) )  {
         fprintf(stderr,"ERROR (%s) error in accessing the 'corrdXdY' values.\n",__func__);
         return 1;
      }
      /* uflags */
      if (rsprod_field_loadFromNetCDF(&Fuflags,"uflag",ncid)) {
         fprintf(stderr,"ERROR (%s) While loading the 'uflag' dataset.\n",__func__);
         return 1;
      }
      short *Fuflags_data;
      if (rsprod_data_getType(Fuflags->data) != RSPROD_SHORT) {
         fprintf(stderr,"ERROR (%s) The 'uflags' dataset is not of type SHORT.\n",__func__);
         return 1;
      }
      if ( (*Fuflags->data->methods->accessValues)(Fuflags->data,&Fuflags_data) )  {
         fprintf(stderr,"ERROR (%s) error in accessing the 'uflags' values.\n",__func__);
         return 1;
      }
   } else {
      fprintf(stderr,"WARNING (%s) Skip reading uncertainties <%s> and <%s>.\n",__func__,sigdXname,sigdYname);
   }

   /* allocate the output arrays and transfer the data */
   rsprod_dimensions *dims = Flat->dims; 
   if (dims->nbdims != 2) {
      fprintf(stderr,"ERROR (%s) Fields have %u dimensions. Was expecting 2.\n",__func__,dims->nbdims);
      return 1;
   }
   for (short d = 0 ; d < 2 ; d ++) {
      if (!strcmp(dims->name[d],"xc")) out_dims[0] = dims->length[d];
      else if (!strcmp(dims->name[d],"yc")) out_dims[1] = dims->length[d];
   }
   out_dims[2] = out_dims[0] * out_dims[1];

   if (rsprod_data_getVoidStar(FdriftX->data,(void **)driftX)) {
      fprintf(stderr,"ERROR (%s) While copying 'driftX' data.\n",__func__);
      return 1;
   }
   if (rsprod_data_getVoidStar(FdriftY->data,(void **)driftY)) {
      fprintf(stderr,"ERROR (%s) While copying 'driftX' data.\n",__func__);
      return 1;
   }
   if (rsprod_data_getVoidStar(Flat->data,(void **)lat)) {
      fprintf(stderr,"ERROR (%s) While copying 'lat' data.\n",__func__);
      return 1;
   }
   if (rsprod_data_getVoidStar(Flon->data,(void **)lon)) {
      fprintf(stderr,"ERROR (%s) While copying 'lon' data.\n",__func__);
      return 1;
   }
   if (rsprod_data_getVoidStar(Fflags->data,(void **)flags)) {
      fprintf(stderr,"ERROR (%s) While copying 'flags' data.\n",__func__);
      return 1;
   }
   if (haveUncertainties) {
      if (rsprod_data_getVoidStar(FsigdX->data,(void **)sigdX)) {
         fprintf(stderr,"ERROR (%s) While copying 'sigdX' data.\n",__func__);
         return 1;
      }
      if (rsprod_data_getVoidStar(FsigdY->data,(void **)sigdY)) {
         fprintf(stderr,"ERROR (%s) While copying 'sigdY' data.\n",__func__);
         return 1;
      }
      if (rsprod_data_getVoidStar(FcorrdXdY->data,(void **)corrdXdY)) {
         fprintf(stderr,"ERROR (%s) While copying 'corrdXdY' data.\n",__func__);
         return 1;
      }
      if (rsprod_data_getVoidStar(Fuflags->data,(void **)uflags)) {
         fprintf(stderr,"ERROR (%s) While copying 'uflags' data.\n",__func__);
         return 1;
      } 

   } else {
      *sigdX = *sigdY = *corrdXdY = NULL;
      *uflags = NULL;
   }

   /* get the fillvalue out of the drift datasets */
   if ( rsprod_attributes_accessValue_Float(FdriftX->attr,"_FillValue",fillvalue) ) {
      fprintf(stderr,"ERROR (%s) cannot access float fill value for field 'driftX'.\n",__func__);
      return 1;
   }
   /* get the fillvalue out of the time datasets */
   int fillvalt;
   if ( rsprod_attributes_accessValue_Int(Ft0->attr,"_FillValue",&fillvalt) ) {
      fprintf(stderr,"ERROR (%s) cannot access int fill value for field 't0'.\n",__func__);
      return 1;
   }

   /* correct the T0 and T1 time so that we get 'seconds since 1970-01-01' */
   *t0 = fmMalloc(sizeof(**t0) * out_dims[TDIM]);
   *t1 = fmMalloc(sizeof(**t1) * out_dims[TDIM]);
   for (size_t e = 0 ; e < out_dims[TDIM] ; e ++) {
      if (FdriftX_data[e] != *fillvalue) {
	 (*t0)[e] = Ft0_data[e] + refT0;
	 (*t1)[e] = Ft1_data[e] + refT1;
      } else {
	 (*t0)[e] = (*t1)[e] = fillvalt;
      }
   }

   /* get the grid and projection information */
   size_t lengths[2];
   short  wTime = 0;
   fmsec1970 timeVal[1];
   char timeUnit[128];memset(timeUnit,0,128);
   projret = strstr(infname, "ease");
   if (projret) {
     ret = decode_CF1_struct(ncid,lengths,A,B,projstr,wTime,timeVal,timeUnit,PROJOBJECT_NAME_EASE);
     if (ret) {
       fprintf(stderr,"ERROR (%s) cannot access grid/proj information from <%s>.\n",__func__,PROJOBJECT_NAME_EASE);
       return 1;
     }
   } else {
     ret = decode_CF1_struct(ncid,lengths,A,B,projstr,wTime,timeVal,timeUnit,PROJOBJECT_NAME_POLSTERE);
     if (ret) {
       fprintf(stderr,"ERROR (%s) cannot access grid/proj information from <%s>.\n",__func__,PROJOBJECT_NAME_POLSTERE);
      return 1;
     }
   }
   /* initialize the proj object */
   *proj = NULL;
   *proj = pj_init_plus(projstr); 


   /* Close NetCDF file */
   ret = nc_close(ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot nc_open() file for reading:\n\t%s\n",__func__,infname);
      return 1;
   }
   /* delete the rsprod_fields */
   rsprod_field_delete(Flat); rsprod_field_delete(Flon);
   rsprod_field_delete(Ft0); rsprod_field_delete(Ft1);
   rsprod_field_delete(FdriftX); rsprod_field_delete(FdriftY);
   rsprod_field_delete(Fflags);
   if (haveUncertainties) {
      rsprod_field_delete(FsigdX); rsprod_field_delete(FsigdY);rsprod_field_delete(FcorrdXdY); 
      rsprod_field_delete(Fuflags);
   }

   return 0;
}
