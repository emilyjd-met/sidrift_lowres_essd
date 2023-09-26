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
 *   laplacian_preprocessing.c 
 *
 * PURPOSE:
 *   To implement the Laplacian filtering of daily maps
 *
 * REQUIREMENTS:
 *
 * INPUT:
 *   
 * OUTPUT:
 *   The Laplacian datasets are appended to the input netCDF file (with -A).
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
 *   Thomas Lavergne, met.no/FoU, 07.03.2009:  Change '-' to '_' for names of netCDF datasets as '-' is invalid 
 *                                             (does not pass in NC3.6.3).
 *   Thomas Lavergne, met.no/FoU, 29.04.2010:  BUG: now handle missing data in ice_edge product
 *
 */ 


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <fmutil.h>
#include <rsprod.h>
#include <errorcodes.h>
#include <icedrift_filenames.h>
#include <icedrift_flags.h>

char progname[] = "laplacian_preprocessing";

int isIceObs(short ice) {
   
   int isice = 1;

   if ( ice >= 9 ) /* land and coasts */
      isice = 0;
   //else if ( (ice == 1) || (ice == 2) ) /* water and low concentration ice */
   else if ( (ice == 1) ) /* water */ 
      isice = 0;
   else if ( (ice == 0) ) /* missing value in edge product */
      isice = 0;
   
   return isice;
}

int isValidIceObs(float obs,float fillv,short ice) {

   int isvalid = 1;
   if ( obs == fillv ) 
      isvalid = 0;
   else if (!isIceObs(ice)) {
      isvalid = 0;
   }

   return isvalid;

}

void usage(void) {
   fprintf(stderr,"%s -i <infile> (-o <outfile> | -A) [-m <sizeList>] [-s <icefile>]\n",progname);
   fprintf(stderr,"\t-A append the computed field(s) to the input file (incompatible with -o).\n");
   fprintf(stderr,"\t-m <sizeList>    : comma separated list of sizes for the median filter.\n");
   fprintf(stderr,"\t-a <sizeList>    : comma separated list of sizes for the mean filter. (incompatible with -m)\n");
   fprintf(stderr,"\t-s <icefile>     : optional independent ice mask.\n");
   exit(OSISAF_ERROR_CMDSYNTAX);
}


int main (int argc, char *argv[]) {

   int ret;

   int  iflg,oflg,mflg,aflg,Aflg,sflg;
   char *infile,*outfile,*icefile;
   char **filterSizeStr;
   unsigned int nbFilters;
   extern char *optarg;
   iflg = oflg = aflg = mflg = Aflg = sflg = 0;
   while ((ret = getopt(argc, argv, "i:o:m:a:As:")) != EOF) {
      switch (ret) {
	 case 'i':
	    infile = fmMalloc(strlen(optarg)+1);
        if (!strcpy(infile, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
        iflg++;
        break;
	 case 'o':
	    outfile = fmMalloc(strlen(optarg)+1);
        if (!strcpy(outfile, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
        oflg++;
        break;
	 case 'A':
	    Aflg++;
	    break;
	 case 'm':
	    mflg++;
	    fmstrsplit(optarg,",",&nbFilters,&filterSizeStr);
	    break;
	 case 'a':
	    aflg++;
	    fmstrsplit(optarg,",",&nbFilters,&filterSizeStr);
	    break;
	 case 's':
	    sflg++;
        icefile = fmMalloc(strlen(optarg)+1);
        if (!strcpy(icefile, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
	    break;
	 default:
	    usage();
	    break;
      }
   }
   if (!iflg) usage();
   if (aflg && mflg) usage();
   if (oflg && Aflg) usage();
   if (!oflg && !Aflg) usage();

   int *filterSizes, *halfFilterSizes;
   if (mflg || aflg) {
      filterSizes = fmMalloc(nbFilters*sizeof(int));
      halfFilterSizes = fmMalloc(nbFilters*sizeof(int));
      for (size_t i = 0 ; i < nbFilters ; i++) {
	 char *endPtr;
	 filterSizes[i] = strtol(filterSizeStr[i],&endPtr,10);
	 if ((*endPtr != '\0')) 
	    usage();
	 float halfSize = (filterSizes[i] - 1)/2.0;
	 halfFilterSizes[i] = (int)rintf(halfSize);
	 if ( (filterSizes[i] <= 0) || (halfSize < 0) || (halfFilterSizes[i] != (int)halfSize)) {
	    fprintf(stderr,"ERROR (%s) the filter sizes (-m or -a) must be positive, odd numbers and >= 1\n",progname);
	    return OSISAF_ERROR_CMDVALUE;
	 }
      }
   } else {
      nbFilters = 1;
   }
  

   rsprod_field   **obs,**oflag,*ice;

   short *icedata = NULL;

   if (sflg) {
      printf("\tLoad the ice and ocean masks from separate file.\n");
   
      int ncid;
      /* Open the input netCDF file */
      if (nc_open(icefile, NC_NOWRITE, &ncid)) {
         fprintf(stderr,"ERROR (%s) cannot open nc file <%s>.\n",progname,icefile);
         exit(OSISAF_ERROR_CMDVALUE);
      }
      if (rsprod_field_loadFromNetCDF(&ice,"ice_edge",ncid)) {
	 fmerrmsg(progname,"icemask file does not contain <%s>.","ice_edge");
	 return (OSISAF_ERROR_OTHER);
      }
      nc_close(ncid);

      if ( (*ice->data->methods->accessValues)(ice->data,&icedata) ) {
	 fmerrmsg(progname,"the ice_edge mask is not of type short.");
	 return(OSISAF_ERROR_OTHER);
      }

   }
   
   printf("\tLoad the observations from the input file.\n");

   /* split the tcimage filename to access the information */
   char  *Instrument, *Platform, *WaveBands, *area, *centralTime;
   short twght;
   ret = tcimage_split_filename(infile,&twght,&Instrument,&Platform,&WaveBands,&area,&centralTime);

   char   **waveBands;
   unsigned int   nbWaveBands;
   fmstrsplit(WaveBands,Ssep,&nbWaveBands,&waveBands);

   int ncid;
   /* Open the input netCDF file */
   if (nc_open(infile, NC_NOWRITE, &ncid)) {
      fprintf(stderr,"ERROR (%s) cannot open nc file <%s>.\n",progname,infile);
      exit(OSISAF_ERROR_CMDVALUE);
   }
   obs   = fmMalloc(sizeof(rsprod_field *) * nbWaveBands);
   for (short b = 0 ; b < nbWaveBands ; b++) {
      char shortName[32];
      sprintf(shortName,"%s",waveBands[b]);
      if (rsprod_field_loadFromNetCDF(&(obs[b]),shortName,ncid)) {
	 fprintf(stderr,"ERROR (%s) inputfile does not contain <%s>.\n",progname, shortName);
	 return (OSISAF_ERROR_OTHER);
      }
   }

   nc_close(ncid);

   /* get the 'grid mapping' parameter for the observations */
   char gridMapping[] = "crs";
   char coords[] = "lat lon";

   /* Get a direct (typed) pointer to the various data sets */
   if ( (obs[0])->dims->nbdims != 3 ) {
      fprintf(stderr,"ERROR (%s) the observation does not have 3 dimensions (but %u).\n",progname,(obs[0])->dims->nbdims);
      return(OSISAF_ERROR_OTHER);
   }
   float **obsdata;
   obsdata = fmMalloc(sizeof(float *) * nbWaveBands);
   for (short b = 0 ; b < nbWaveBands ; b++) {
      if ( (*(obs[b])->data->methods->accessValues)((obs[b])->data,&(obsdata[b])) ) {
	 fprintf(stderr,"ERROR (%s) the observation is not of type float.\n",progname);
	 return(OSISAF_ERROR_OTHER);
      }
   }
   float obsFillValuedata;
   if ( rsprod_attributes_accessValue_Float((obs[0])->attr,"_FillValue",&obsFillValuedata) ) {
      fprintf(stderr,"ERROR (%s) problem accessing the _FillValue or missing_value for %s\n",progname,(obs[0])->name);
      return(OSISAF_ERROR_OTHER);
   }

   long matIndex[3];

   float **laplacian     = fmMalloc(sizeof(float *)*nbWaveBands);
   short **laplacian_flag = fmMalloc(sizeof(short *)*nbWaveBands);
   size_t *nbelems   = fmMalloc(sizeof(*nbelems)*nbWaveBands);
   for (short b = 0 ; b < nbWaveBands ; b++) {
      nbelems[b] = (obs[b])->dims->totelems;
      laplacian[b]     = fmMalloc(sizeof(float) * nbelems[b]);
      laplacian_flag[b] = fmMalloc(sizeof(short) * nbelems[b]);
   }

   short **oflagdata;
   oflagdata = fmMalloc(sizeof(short *) * nbWaveBands);
   for (short b = 0 ; b < nbWaveBands ; b++) {
      oflagdata[b] = fmMalloc(sizeof( short ) * nbelems[b]);
      for (unsigned long el = 0 ; el < nbelems[b] ; el++) {
         if ( obsdata[b][el] == obsFillValuedata ) {
            oflagdata[b][el] = TCIMAGE_NODATA;
         } else {
            oflagdata[b][el] = TCIMAGE_OK;
         }
      }
   }

   short defaultIce = 3;
   /* compute the Laplacian field approximation */
   matIndex[0] = 0; /* the 'time' dimension */
   for (short b = 0 ; b < nbWaveBands ; b++) {
      
      printf("\tCompute Laplacian field for %s.\n",waveBands[b]);

      for (unsigned long el = 0 ; el < nbelems[b] ; el++) {

	 laplacian[b][el] = obsFillValuedata;

	 short icemaskvalue = defaultIce;
	 if (icedata) {
            icemaskvalue = icedata[el];
         }

	 if ( (oflagdata[b][el] != TCIMAGE_OK) || ( ! isIceObs(icemaskvalue) ) ) {
	    laplacian_flag[b][el] = TCIMAGE_UNPROCESSED;
            laplacian[b][el] = obsFillValuedata;
	    continue;
         }

	 int cptP,cptM,cptU;
	 float termP,termM;
	 cptP=0;cptM=0;
	 termP=0;termM=0;
	 cptU=0;
	 for (int j = -2 ; j <= 2 ; j++) {
	    for (int i = -2 ; i <= 2 ; i++) {
	       if ((i == 0) && (j == 0)) continue;

	       matIndex[1] = j;
	       matIndex[2] = i;
	       unsigned int celem = vecIndexRelative(obs[b]->dims,el,3,matIndex);
	       if (celem == nbelems[b]) break;

	       short icemaskvaluen = defaultIce;
	       if (icedata) icemaskvaluen = icedata[celem];
	       if (oflagdata[b][celem] == TCIMAGE_NODATA) { cptU++; continue; }
	       else if (oflagdata[b][celem] != TCIMAGE_OK) continue;
	       if (icemaskvalue && !(isIceObs(icemaskvaluen))) continue;  /* if the center is on ice, then only use ice observations */

	       int absi = (i>=0?i:-i);int absj = (j>=0?j:-j);
	       if ( (absi <= 1) && (absj <= 1) ) {
		  cptP++;
		  termP += obsdata[b][celem];
	       } else {
		  cptM++;
		  termM += obsdata[b][celem];
	       }
	    }
	 }
	 if (cptU > 1) {
	    laplacian_flag[b][el] = TCIMAGE_NODATA;
	 } else if ((cptP >= 5) && (cptM >= 9)) {
	    laplacian[b][el]     = termP/cptP - termM/cptM;
	    laplacian_flag[b][el] = TCIMAGE_OK;
	 } else {
	    laplacian_flag[b][el] = TCIMAGE_UNPROCESSED;
	 }

      }
   }

   /* if needed, apply the filters */
   float ***flaplacians      = fmMalloc(sizeof(float **) * nbWaveBands);
   short ***flaplacians_flag = fmMalloc(sizeof(short **) * nbWaveBands);
   for (short b = 0 ; b < nbWaveBands ; b++) {
      flaplacians[b]      = fmMalloc(sizeof(float *)*nbFilters);
      flaplacians_flag[b] = fmMalloc(sizeof(short *)*nbFilters);
   }

   if ( mflg || aflg ) {

      for (short b = 0 ; b < nbWaveBands ; b++) {
   
	 printf("\tCompute the spatially filtered Laplacian field(s) for %s [%s filter]:\n",
	       waveBands[b],(mflg?"median":"mean"));

	 for (size_t i = 0 ; i < nbFilters ; i++) {
	    
	    int filtersize = filterSizes[i];
	    int ifilter    = halfFilterSizes[i];

	    if (filtersize == 1) {
	       /* no real averaging. Used to force the program to output the unfiltered field */
	       flaplacians[b][i]      = laplacian[b];
	       flaplacians_flag[b][i] = laplacian_flag[b];
	       continue;
	    }

	    printf("\t\tFilter fields with %s x %s window.\n",filterSizeStr[i],filterSizeStr[i]);

	    float  *cflaplacian      = fmMalloc(sizeof(float)*nbelems[b]);
	    short  *cflaplacian_flag = fmMalloc(sizeof(short)*nbelems[b]);
	    size_t *nvalues     = fmMalloc(sizeof(size_t)*nbelems[b]);
	    flaplacians[b][i]      = cflaplacian;
	    flaplacians_flag[b][i] = cflaplacian_flag;
	   
	    long nbsearch = filtersize*filtersize;
	    float **values = fmMalloc(nbsearch*sizeof(float *));
	    short **flags  = fmMalloc(nbsearch*sizeof(short *));

	    matIndex[0] = 0; /* 'time' is the first dimension */
	    for (unsigned long el = 0 ; el < nbelems[b] ; el++) {
	      
	       /* first we fill the field with fillvalues */
	       cflaplacian[el] = obsFillValuedata;

	       /* compute an average of all valid surrounding laplacian values */
	       /* the search area has size filtersize x filtersize */
	       unsigned int nvals = 0;
	       double avg = 0;
	       for (int j = -ifilter ; j <= ifilter ; j++) {
		  for (int i = -ifilter ; i <= ifilter ; i++) {
		     matIndex[1] = j;
		     matIndex[2] = i;
		     unsigned int celem = vecIndexRelative(obs[b]->dims,el,3,matIndex);
		     if (celem == nbelems[b]) break;

		     short *flgp = &(laplacian_flag[b][celem]);
		     flags[nvals] = flgp;
		     if (*flgp != TCIMAGE_OK) continue;

		     float *valp = &(laplacian[b][celem]);
	       
		     values[nvals] = valp;
		     avg += *valp;
		     nvals++; 
		  }
	       }
	       /* store the value ONLY if more than half of the search window was valid */
	       if (nvals >= 0.5*nbsearch) {
		  avg /= nvals;
		  cflaplacian[el] = avg;
		  nvalues[el]     = nvals;
		  cflaplacian_flag[el] = TCIMAGE_OK;
	       } else {
		  cflaplacian_flag[el] = TCIMAGE_UNPROCESSED;
	       }
	    }
	 
	    if (mflg) {
	       /* if we have enough valid laplacian values around, select the 
		* closest to the calculated average and use it for
		* Laplacian[el] */
		for (unsigned long el = 0 ; el < nbelems[b] ; el++) {

		   if ( cflaplacian_flag[el] == TCIMAGE_OK ) {
		      size_t nvals = nvalues[el];
		      double mindist = 1000;
		      float  selectedValue;
		      double avg = cflaplacian[el];
		      for (int i = 0 ; i < nvals ; i++) {
			 double dist = (avg - *(values[i]));
			 dist = (dist > 0 ? dist : -dist);
			 if (dist < mindist) {
			    mindist = dist;
			    selectedValue = *(values[i]);
			 }
		      }
		      cflaplacian[el] = selectedValue;
		   }
		}
	    }

	    free(values);
	    free(nvalues);

	 }
      }
   } else {
      for (short b = 0 ; b < nbWaveBands ; b++) {
	 flaplacians[b][0] = laplacian[b];
	 flaplacians_flag[b][0] = laplacian_flag[b];
      }
   }

      
   /* create the laplacian FIELDs (RSPROD) */
   rsprod_field ***lapFields = fmMalloc(nbWaveBands * sizeof(rsprod_field **));
   rsprod_field ***flgFields = fmMalloc(nbWaveBands * sizeof(rsprod_field **));
   for (short b = 0 ; b < nbWaveBands ; b++) {
      lapFields[b] = fmMalloc(nbFilters*sizeof(rsprod_field *));
      flgFields[b] = fmMalloc(nbFilters*sizeof(rsprod_field *));

      for (size_t nfilter = 0 ; nfilter < nbFilters ; nfilter ++) {

	 /* names (short and long) */
	 char shortName[32],shortNameFlag[32];
	 char longName[1024],longNameFlag[1024];
	 if ((!mflg && !aflg) 
	       || 
	    ((mflg || aflg) && (filterSizes[nfilter]==1))) {
	    sprintf(shortName,"%s_Lap",waveBands[b]);
	    sprintf(longName,"Laplacian of TCobs(%s) - Unfiltered",waveBands[b]);
	    sprintf(shortNameFlag,"%s_Lap_flag",waveBands[b]);
	    sprintf(longNameFlag,"Flags for Laplacian of TCobs(%s) - Unfiltered",waveBands[b]);
	 } else {
	    sprintf(shortName,"%s_Lap%s",waveBands[b],filterSizeStr[nfilter]);
	    sprintf(longName,"Laplacian of TCobs(%s) - %sx%s %s filter",waveBands[b],
		  filterSizeStr[nfilter],filterSizeStr[nfilter],(mflg?"median":"mean"));
	    sprintf(shortNameFlag,"%s_Lap%s_flag",waveBands[b],filterSizeStr[nfilter]);
	    sprintf(longNameFlag,"Flags for Laplacian of TCobs(%s) - %sx%s %s filter",waveBands[b],
		  filterSizeStr[nfilter],filterSizeStr[nfilter],(mflg?"median":"mean"));}

	 /* create */
	 if (rsprod_field_createStandard(&(lapFields[b][nfilter]),shortName,RSPROD_FLOAT,
	       obs[b]->dims->totelems,obs[b]->dims,longName,NULL,"1",&obsFillValuedata,NULL,NULL,coords,
	       gridMapping,flaplacians[b][nfilter])) {
	    fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",progname,shortName);
	    return OSISAF_ERROR_OTHER;
	 }
	 if (rsprod_field_createStandard(&(flgFields[b][nfilter]),shortNameFlag,RSPROD_SHORT,
	       obs[b]->dims->totelems,obs[b]->dims,longNameFlag,NULL,"1",NULL,NULL,NULL,coords,
	       gridMapping,flaplacians_flag[b][nfilter])) {
	    fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",progname,shortName);
	    return OSISAF_ERROR_OTHER;
	 }

      }
   }


   printf("\tWrite Laplacian fields to file.\n");

   if (Aflg) {
      if ((ret = nc_open(infile,NC_WRITE,&ncid)) != NC_NOERR) {
	 fprintf(stderr,"ERROR (%s) Cannot open nc file %s. <%s>\n",progname,infile,nc_strerror(ret));
	 return OSISAF_ERROR_OTHER;
      }
   } else {
      if ((ret = nc_create(outfile,NC_CLOBBER,&ncid)) != NC_NOERR) {
	 fprintf(stderr,"ERROR (%s) Cannot open nc file %s. <%s>\n",progname,outfile,nc_strerror(ret));
	 return OSISAF_ERROR_OTHER;
      }
      for (short b = 0 ; b < nbWaveBands ; b++) {
	 ret = rsprod_field_writeToNetCDF(obs[b],ncid);
	 if (ret) {
	    fprintf(stderr,"ERROR (%s) Cannot write Field %s to netcdf file.\n",progname,obs[b]->name);
	    return OSISAF_ERROR_OTHER;
	 } 
     fprintf(stderr,"VERBOSE (%s) Wrote field %s to netcdf file.\n",progname,obs[b]->name);
      }
   }

   for (short b = 0 ; b < nbWaveBands ; b++) {
      for (size_t nfilter = 0 ; nfilter < nbFilters ; nfilter ++) {

	 ret = rsprod_field_writeToNetCDF(lapFields[b][nfilter],ncid);
	 if (ret) {
	    fprintf(stderr,"ERROR (%s) Cannot write Field %s to netcdf file.\n",progname,lapFields[b][nfilter]->name);
	    return OSISAF_ERROR_OTHER;
	 }
	 rsprod_field_delete(lapFields[b][nfilter]);
	 
	 ret = rsprod_field_writeToNetCDF(flgFields[b][nfilter],ncid);
	 if (ret) {
	    fprintf(stderr,"ERROR (%s) Cannot write Field %s to netcdf file.\n",progname,flgFields[b][nfilter]->name);
	    return OSISAF_ERROR_OTHER;
	 }
	 rsprod_field_delete(flgFields[b][nfilter]);
      }
   }

   if (sflg) {
      /* put the ice and ocean mask for later drift processing */
      ret = rsprod_field_writeToNetCDF(ice,ncid);
      if (ret) {
         fprintf(stderr,"ERROR (%s) Cannot write Field %s to netcdf file.\n",progname,ice->name);
         return OSISAF_ERROR_OTHER;
      }
   }

   nc_close(ncid);
   
   for (short b = 0 ; b < nbWaveBands ; b++) {
      free(lapFields[b]);
      free(flgFields[b]);
      rsprod_field_delete(obs[b]);
   }
   if (sflg)
      rsprod_field_delete(ice);


   return OSISAF_EXIT_CORRECT;     
}

