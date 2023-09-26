/* *****************************************************************
 * COPYRIGHT:
 * EUMETSAT
 *
 * PRODUCED BY:
 * Norwegian Meteorological Institute (met.no)
 * Research and Development Department
 * P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
 *
 * This SW was developed by met.no and the Danish Meteorological
 * Institute (DMI) within the context of the Co-operation Agreement
 * for the development of a pilot SAF on Ocean and Sea Ice.
 * *****************************************************************/ 

/*
 * NAME:
 *   apply_summerhold_filter
 *
 * PURPOSE:
 *   To prepare a 'lev3' filtered file with an empty grid if the 
 *   product date is in the Summer hold period (May 1st -> September 30th).
 *   For the rest of the year, the lev3 is merely a copy of the lev2.
 *
 * REQUIREMENTS:
 *
 * INPUT:
 *   Required parameters are:
 *   -i <inputFile>    :   lev2 product file as output from icedrift_solve_simplex.
 *   -o <outputDir>    :   directory where the output file will be placed
 *   
 * OUTPUT:
 *
 * NOTES:
 * NA
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 *   Thomas Lavergne (met.no)
 *
 * MODIFIED:
 *    Thomas Lavergne, met.no/FoU, 07.05.2009   :   First version
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <rsprod.h>
#include <fmutil.h>
#include <errorcodes.h>
#include "directoryname.h"
#include "icedrift_filenames.h"
#include "icedrift_flags.h"
#include "func_summerhold.h"

/* Global variables */
char progname[] = "apply_summerhold_filter";

/* Local declarations */
static void usage(void); 

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++  MAIN DRIVER PROGRAM                                       ++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int main(int argc, char *argv[]) {
   
   int ret;

   int   iflg,oflg;
   char *infile,*outdir;

   iflg = oflg = 0;
   extern char *optarg;
   while ((ret = getopt(argc, argv, "i:o:")) != EOF) {
	  switch (ret) {
		 case 'i':
			infile = fmMalloc(strlen(optarg)+1);
            if (!strcpy(infile, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
			iflg++;
			break;
		 case 'o':
			outdir = fmMalloc(strlen(optarg)+1);
			if (!strcpy(outdir, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
			oflg++;
			break;
		 default:
			usage();
			break;
	  }
   }
   if (!iflg || !oflg) usage();

   fmlogmsg(progname,"Start run");

   /* analyse the input file name to read the product date and the satellite */
   char *Instrument;
   char *Platform;
   char *WaveBands;
   char *method;
   char *levf;
   char *area;
   char *startdate;
   char *enddate;
   ret = driftproduct_split_filename_proc(infile,&Instrument,&Platform,&WaveBands,&method,&levf,&area,&startdate,&enddate);
   if (ret) {
	  fprintf(stderr,"ERROR (%s) Input file is not a valid ice drift netcdf product file\n\t%s\n",__func__,infile);
	  exit(OSISAF_ERROR_CMDVALUE);
   }

   /* open the input file for reading */
   int incid;
   ret = nc_open(infile,NC_NOWRITE,&incid);
   if (ret != NC_NOERR) {
	  fprintf(stderr,"ERROR (%s) Cannot open netCDF file for reading.\n\t%s\n",__func__,infile);
	  exit(OSISAF_ERROR_SYSIO);
   }
   
   /* form the name of the output file : change 'lev2' in 'lev3' */
   char *endInDir, *basename;
   fmbasename(infile,&endInDir,&basename);
   char *outfile = fmMalloc(strlen(basename)+1);
   if (!strcpy(outfile,basename)) {
	  fprintf(stderr,"ERROR (%s) Cannot copy string <%s>",__func__,basename);
	  exit(OSISAF_ERROR_SYSMEM);
   }
   char *lev2string = strstr(outfile,"lev2");
   if (!lev2string) {
	  fprintf(stderr,"ERROR (%s) Input file is not a 'lev2' product file.\n\t<%s>\n",__func__,basename);
	  exit(OSISAF_ERROR_CMDVALUE);
   }
   memcpy(lev2string,"lev3",strlen("lev3"));
   char *path2_outfile = fmMalloc(strlen(outdir)+1+strlen(outfile)+1);
   sprintf(path2_outfile,"%s/%s",outdir,outfile);

   /* open the output file for writing */
   int oncid;
   ret = nc_create(path2_outfile,NC_CLOBBER,&oncid);
   if (ret != NC_NOERR) {
	  fprintf(stderr,"ERROR (%s) Cannot open netCDF file for writing.\n\t<%s>\n",__func__,path2_outfile);
	  exit(OSISAF_ERROR_SYSIO);
   }

   /* decide if we are in the summer hold period */
   int productdate_in_summerhold;
   ret = test_summer_hold(ymdh2fmsec1970_alt(startdate,0),area,Instrument,WaveBands,&productdate_in_summerhold);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Invalid product date %s\n",__func__,startdate);
	  exit(OSISAF_ERROR_OTHER);
   }
   if (productdate_in_summerhold) {
	  fmlogmsg(progname,"Date %s is in %s summer hold period. Empty the product grid.",startdate, area);
   } else {
	  fmlogmsg(progname,"Date %s is *not* in %s summer hold period. Copy the product.",startdate, area);
   }

   /* loop trough all the datasets and copy them to the output netcdf file */
   int nvars;
   ret = nc_inq_nvars(incid, &nvars);
   if (ret != NC_NOERR) {
	  fprintf(stderr,"ERROR (%s) Cannot inquire the number of variables in the input netCDF file.\n",__func__);
	  exit(OSISAF_ERROR_OTHER);
   }
   for (int var = 0 ; var < nvars ; var++) {
	  char varname[NC_MAX_NAME+1];
	  ret = nc_inq_varname(incid,var,varname);
	  if (ret != NC_NOERR) {
		 fprintf(stderr,"ERROR (%s) No variable with varid %d. Skip it.\n",__func__,var);
		 continue;
	  }
	  rsprod_field *field;
	  ret = rsprod_field_loadFromNetCDF(&field,varname,incid);
	  if (ret) {
		 fprintf(stderr,"ERROR (%s) Unable to load variable %s from input file. Skip it.\n",__func__,varname);
		 continue;
	  }
	  /* if we are in summer hold, modify the content of some of the datasets */
	  if (productdate_in_summerhold) {
		 if (strcmp(varname,"lat") && strcmp(varname,"lon") && strcmp(varname,"xc") && strcmp(varname,"yc") 
			   && strcmp(varname,"Polar_Stereographic_Grid") && strcmp(varname,"Lambert_Azimuthal_Grid") && strcmp(varname,"time") && strcmp(varname,"time_bnds")) {
			/* get the number of elements */
			size_t nbVals  = rsprod_data_getNbvalues(field->data);
			if ( !strcmp(varname,"flag") ) {
			   printf("VERBOSE (%s) Change content of %s (line %d)\n",progname,varname,__LINE__);
			   /* we change the flags values everywhere the ice drift could have been valid 
				* (but not the land, open water or missing areas) 
				*/
			   short *flags;
			   if ( (*field->data->methods->accessValues)(field->data,&flags)) {
				  fprintf(stderr,"ERROR (%s) Unable to access the values in %s.\n",__func__,varname);
				  exit(OSISAF_ERROR_OTHER);
			   }
			   for (size_t e = 0 ; e < nbVals ; e++) {
				  if ( (flags[e] == ICEDRIFT_OUTSIDE_IMGBORDER)   ||
				   (flags[e] == ICEDRIFT_CLOSETO_IMGBORDER)     ||
				   (flags[e] == ICEDRIFT_CLOSETO_MISSING)       ||
				   (flags[e] == ICEDRIFT_CLOSETO_UNPROCESSED)   ||
				   (flags[e] == ICEDRIFT_CENTER_OVER_LAND)      ||
				   (flags[e] == ICEDRIFT_NOICE) ) {
					 ; /* do nothing */
				  } else {
					 /* put the special flag for summer hold */
					 flags[e] = ICEDRIFT_SUMMERHOLD;
				  }
			   }
			} else if ( !strcmp(varname,"uflag") ) {
			   printf("VERBOSE (%s) Change content of %s (line %d)\n",progname,varname,__LINE__);
			   /* special case for the posterior uncertainty flag */
			   short fv = ICEDRIFTPOST_NOVECTOR;
			   if ( (*field->data->methods->initValues)(field->data,&fv) ) {
				  fprintf(stderr,"ERROR (%s) Unable to init the values in %s. Skip it.\n",__func__,varname);
				  continue;
			   }
			} else {
			   printf("VERBOSE (%s) Change content of %s (line %d)\n",progname,varname,__LINE__);
			   /* for the other datasets : read the fillvalue and erase all valid data */
			   if (rsprod_data_getType(field->data) == RSPROD_FLOAT) {
				  float fv;
				  if ( rsprod_attributes_accessValue_Float(field->attr,"_FillValue",&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to access fill value for %s. Skip it.\n",__func__,varname);
					 continue;
				  }
				  if ( (*field->data->methods->initValues)(field->data,&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to init the values in %s. Skip it.\n",__func__,varname);
					 continue;
				  }
			   } else if (rsprod_data_getType(field->data) == RSPROD_DOUBLE) {
				  double fv;
				  if ( rsprod_attributes_accessValue_Double(field->attr,"_FillValue",&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to access fill value for %s. Skip it.\n",__func__,varname);
					 continue;
				  }
				  if ( (*field->data->methods->initValues)(field->data,&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to init the values in %s. Skip it.\n",__func__,varname);
					 continue;
				  }
			   } else if (rsprod_data_getType(field->data) == RSPROD_SHORT) {
				  short fv;
				  if ( rsprod_attributes_accessValue_Short(field->attr,"_FillValue",&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to access fill value for %s. Skip it.\n",__func__,varname);
					 continue;
				  }
				  if ( (*field->data->methods->initValues)(field->data,&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to init the values in %s. Skip it.\n",__func__,varname);
					 continue;
				  }
			   } else if (rsprod_data_getType(field->data) == RSPROD_INT) {
				  int fv;
				  if ( rsprod_attributes_accessValue_Int(field->attr,"_FillValue",&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to access fill value for %s. Skip it.\n",__func__,varname);
					 continue;
				  }
				  if ( (*field->data->methods->initValues)(field->data,&fv) ) {
					 fprintf(stderr,"ERROR (%s) Unable to init the values in %s. Skip it.\n",__func__,varname);
					 continue;
				  }
			   } else {
				  fprintf(stderr,"ERROR (%s) This data type is not foreseen (%d). Skip %s.\n",__func__,
						rsprod_data_getType(field->data),varname);
				  continue;
			   }
			}
		 }
	  }
	  ret = rsprod_field_writeToNetCDF(field,oncid);
	  if (ret) {
		 fprintf(stderr,"ERROR (%s) Unable to copy variable %s to output file. Skip it.\n",__func__,varname);
		 continue;
	  }
	  rsprod_field_delete(field);
   }
   /* load the global attributes */
   rsprod_attributes *globAttr;
   ret = rsprod_attributes_loadFromNetCDF(&globAttr,incid,NC_GLOBAL);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot load the global attributes of input file.\n",__func__);
	  exit(OSISAF_ERROR_OTHER);
   }
   /* modify the 'filename' attribute (TODO) */
   /* copy global attributes into output file */
   ret = rsprod_attributes_writeToNetCDF(globAttr,oncid,NC_GLOBAL);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot copy the global attributes to output file.\n",__func__);
	  exit(OSISAF_ERROR_OTHER);
   }

   nc_close(incid);
   nc_close(oncid);
   
   fmlogmsg(progname,"Ouput file is <%s>",path2_outfile);

   fmlogmsg(progname,"Stop run");
   exit(OSISAF_EXIT_CORRECT);
}

void usage(void) {
   fprintf(stdout,"USAGE:\n\t%s -i <input_lev2_File> -o <outputDir>\n",progname);
   exit(OSISAF_ERROR_CMDSYNTAX);
}
