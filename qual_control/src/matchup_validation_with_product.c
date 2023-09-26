
/*
 * PURPOSE:
 * To read from ice-drift satellite products as processed by the OSISAF
 * chain and collocate vectors with available in-situ estimates, extracted
 * via the "extract_drifters_trajectories" software.
 * 
 * INPUT:
 *    o name of the netcdf file containing validation trajectories
 *    o name of a   netcdf file containing the OSISAF ice drift product
 *
 * OUTPUT:
 *    o an ascii file containing all the matchup data, one line per validation trajectories.
 * 
 * BUGS:
 * None known.
 * 
 * AUTHOR: 
 * Thomas Lavergne, 19.02.2008
 *
 * MODIFIED:
 * Thomas Lavergne, met.no/FoU, 05.05.2008      : add the reading and reporting of time period (t0,t1)
 * Thomas Lavergne, met.no/FoU, 01.10.2008      : adapt to read full trajectories from the validation file
 *                                                and select the start and end records closest to the product's
 *                                                points, in space _and_ time.
 * Thomas Lavergne, met.no/FoU, 02.10.2008      : add the -M flag for triggering the 2D collocation
 * Thomas Lavergne, met.no/FoU, 11.03.2009      : also read the uncertainty estimates from the product file
 * Thomas Lavergne, met.no/FoU, 28.10.2009      : implement a nearest-neighbour collocation (to be in line with Phil Hwang, working
 *                                                   on validating Roberto's SAR ice drift at SAMS).
 * Thomas Lavergne, met.no/FoU, 29.10.2009      : add the -D flag for testing the influence of the accurate time collocation.
 * Thomas Lavergne, met.no/FoU, 29.10.2009      : Change the output format for storing the matchup
 * Thomas Lavergne, met.no/FoU, 30.10.2009      : add the -S (NN or BL) flag for choosing between NearestNeighbour and BiLinear interpolation
 *                                                   in the spatial domain.
 * Thomas Lavergne, met.no/FoU, 08.04.2010      : Make NN collocation the default choice.
 * Thomas Lavergne, met.no/FoU, 18.08.2010      : Start adapting to read validation data extracted from DTU's WSM SAR product
 * Thomas Lavergne, met.no/FoU, 19.10.2010      : Migrate to netCDF output
 * Thomas Lavergne, met.no/FoU, 11.08.2011      : Exit with success if there is no matchup
 * Thomas Lavergne, met.no/FoU, 02.11.2011      : Impose to find 4 valid corners in NN collocation
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <projects.h>
#include <netcdf.h>
#include <osisaf_ice_netcdf.h>
#include <ice_common.h>
#include <fmutil.h>
#include <read_obsfile.h>
#include <rsprod.h>
#include <icedrift_common.h>
#include <icedrift_filenames.h>
#include <dbl_list.h>
#include <errorcodes.h>
#include "read_icedrift_product_file.h"
#include "func_drifters_trajectories.h"
#include "func_matchup_prodval.h"


int readValFile(char *filename, 
      size_t *nbValidStations,
      buoyTrajectory **Trajs);

int transformMatchupList(
      /* input */
      dbl_list *allMatchups,fmsec1970 refTime,PJ *pj,
      /* output VAL */
      float **valB_lat,float **valB_lon,int **valB_time,float **valE_lat,float **valE_lon,int **valE_time,
      float **val_dX,float **val_dY,float **val_len,float **val_dir,char ***val_id,char ***val_netw,char ***val_src,
      /* output PRD */
      float **prdB_lat,float **prdB_lon,int **prdB_time,float **prdE_lat,float **prdE_lon,int **prdE_time,
      float **prd_dX,float **prd_dY,float **prd_len,float **prd_dir,short **prd_I,short **prd_J,short **prd_flag,
      float **prd_sX,float **prd_sY,float **prd_cXY,
      /* output COL */
      float **col_dist, float **col_dtB, float **col_dtE, float **col_dDur,
      size_t *nbM);

int prepare_Matchups_nc(
      rsprod_file *outFile,
      size_t nbMatch,
      fmsec1970 refTime_matchup, fmsec1970 refTime0, fmsec1970 refTime1,
      float *valB_lat,float *valB_lon,int *valB_time,float *valE_lat,float *valE_lon,int *valE_time,
      float *val_dX,float *val_dY,float *val_len,float *val_dir,char **val_id,char **val_netw,char **val_src,
      float *prdB_lat,float *prdB_lon,int *prdB_time,float *prdE_lat,float *prdE_lon,int *prdE_time,
      float *prd_dX,float *prd_dY,float *prd_len,float *prd_dir,short *prd_I,short *prd_J,short *prd_flag,
      float *prd_sX,float *prd_sY,float *prd_cXY,
      float *col_dist, float *col_dtB, float *col_dtE, float *col_dDur);

char  progname[]         = "matchup_insitu_with_product";

void usage(void) 
{
   fprintf(stdout,"\tSYNTAX: %s -v <valFiles> -p <productFile> -o <outputFile> [-D <dT>] [-M] [-S <sm>]\n", progname);
   fprintf(stdout,"\t-v <valFiles>    : comma-separated list of files with validation data (netCDF, in-situ and SAR)\n");
   fprintf(stdout,"\t-p <productFile> : file with the OSISAF ice drift product.\n");
   fprintf(stdout,"\t-o <outputFile>  : name of netCDF file for storing the matchup info.\n");
   fprintf(stdout,"\t-M               : if set, the collocation will be performed ignoring the t0 and t1 datasets, taking 12utc as product time.\n");
   fprintf(stdout,"\t-D <dT>          : if set, the collocation will be performed at t0+dT and t1+dT. dT is in (whole) hours.\n");
   fprintf(stdout,"\t-S <sm>          : specify the space collocation method: NN => NearestNeighbour and BL => BiLinear. Default is NN.\n");
   exit(OSISAF_ERROR_OTHER);
}

int main(int argc, char *argv[]) 
{
   extern char *optarg;
   int ret, i, vflg, pflg, oflg, Mflg, Dflg,Sflg;
   char *listOfValFiles,*prodFile,*outFile,*spatialMethod;
   float dTcol;

   /* Interprete commandline arguments */
   vflg = pflg = oflg = Mflg = Dflg = Sflg = 0;
   while ((ret = getopt(argc, argv, "v:p:o:S:D:Mh")) != EOF) {
      switch (ret) {
         case 'v':
            listOfValFiles = fmMalloc(strlen(optarg)+1);
            if (!strcpy(listOfValFiles, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
            vflg++;
            break;
         case 'p':
            prodFile = fmMalloc(strlen(optarg)+1);
            if (!strcpy(prodFile, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
            pflg++;
            break;
         case 'S':
            spatialMethod = fmMalloc(strlen(optarg)+1);
            if (!strcpy(spatialMethod, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
            Sflg++;
            break;
         case 'o':
            outFile = fmMalloc(strlen(optarg)+1);
            if (!strcpy(outFile, optarg)) exit(OSISAF_ERROR_CMDSYNTAX);
            oflg++;
            break;
         case 'M':
            Mflg++;
            break;
         case 'D':
            dTcol = atof(optarg);
            Dflg++;
            break;
         case 'h':
         default:
            usage();
      }
   }
   if ( !vflg || !pflg || !oflg ) {
      fmerrmsg(progname,"One of the input parameters is missing.");
      usage();
   }
   if ( Mflg && Dflg) {
      fprintf(stderr,"ERROR (%s) The -M and -D cannot be used at the same time.\n",progname);
      usage();
   }

   char *o_endDirPart,*o_basename;
   fmbasename(outFile,&o_endDirPart,&o_basename);
   if (!strstr(o_basename,".nc")) {
      fprintf(stderr,"ERROR (%s) The output file name (-o) must be netCDF\n",progname);
      usage();
   }

   int temporal_collocation_scheme = USE3DCOL;
   if (Mflg) temporal_collocation_scheme = USE2DCOL;
   else if (Dflg) temporal_collocation_scheme = round(dTcol);

   int use_nearest_neighbour   = 1;
   if (Sflg) {
      if (!strcmp(spatialMethod,"NN")) {
         use_nearest_neighbour = 1;
      } else if (!strcmp(spatialMethod,"BL")) {
         use_nearest_neighbour = 0;
      } else {
         fprintf(stderr,"ERROR (%s) Unknown parameter to the -S flag. Should be NN or BL.\n",progname);
         exit(OSISAF_ERROR_CMDSYNTAX);
      }
   }
   if (use_nearest_neighbour) {
      printf("VERBOSE (%s) Use Nearest Neighbour spatial collocation.\n",progname);
   } else {
      printf("VERBOSE (%s) Use Bi-Linear spatial collocation.\n",progname);
   }

   /* Analyse the product file name to extract the reference dates for start and end images */
   char *endDirPart,*basename;
   fmbasename(prodFile,&endDirPart,&basename);
   char *Instrument,*Platform,*WaveBands,*drift_method,*Area,*Pstartdate,*Penddate,*levFilter;
   int fmt_type = -9;
   if (is_mdriftproduct_filename(basename,&fmt_type) && (fmt_type == LRSID_FILEFORMAT_PROC)) {
      ret = mdriftproduct_split_filename(basename,&Instrument,&Platform,&drift_method,&levFilter,
            &Area,&Pstartdate,&Penddate);
   } else if (is_driftproduct_filename(basename,&fmt_type) && (fmt_type == LRSID_FILEFORMAT_PROC)) {
      ret = driftproduct_split_filename_proc(basename,&Instrument,&Platform,&WaveBands,&drift_method,
            &levFilter,&Area,&Pstartdate,&Penddate);
   } else if (is_driftproduct_filename(basename,&fmt_type) && (fmt_type == LRSID_FILEFORMAT_FINL)) {
      char *GridInfo;
      char *InstrumentAndPlatform;
      ret = driftproduct_split_filename_final(basename,&Area,&GridInfo, 
                  &InstrumentAndPlatform, &Pstartdate, &Penddate);
      /* cut the dates to YYYYMMDDHH: */
      Pstartdate[10] = '\0';
      Penddate[10]   = '\0';
      unsigned int nbToks;
      char **Toks;
      fmstrsplit(InstrumentAndPlatform,"-",&nbToks,&Toks);
      sprintf(Instrument,"%s",Toks[0]);
      sprintf(Platform,"%s",Toks[1]);
   }
   if ( (fmt_type != LRSID_FILEFORMAT_PROC && fmt_type != LRSID_FILEFORMAT_FINL) || ret) {
      fmerrmsg(progname,"The product file (-p) is not a valid ice drift file.");
      exit(OSISAF_ERROR_OTHER);
   }
   fmsec1970 refT0, refT1;
   if ( strlen(Pstartdate) == 10 ) {
      refT0 = ymdh2fmsec1970_alt(Pstartdate,0);
      refT1 = ymdh2fmsec1970_alt(Penddate,0);
   } else if ( strlen(Pstartdate) == 14+1 ) {
      Pstartdate[14] = '\0';
      Penddate[14]   = '\0';
      refT0 = ymdhms2fmsec1970_alt(Pstartdate,0);
      refT1 = ymdhms2fmsec1970_alt(Penddate,0);
   } else {
      fmerrmsg(progname,"The start and end date is not in an expected format <%s><%s>\n",Pstartdate,Penddate);
      exit(OSISAF_ERROR_OTHER);
   }

   char tstr0[FMUTIL_CFEPOCH_LENGTH+1]; 
   char tstr1[FMUTIL_CFEPOCH_LENGTH+1]; 
   fmsec19702CFepoch(refT0,tstr0);
   fmsec19702CFepoch(refT1,tstr1);
   printf("From filename: refT0=<%s> refT1=<%s>\n",tstr0,tstr1);

   fprintf(stdout,"INFO (%s) The product file is from <%s>-<%s>\n",progname,Instrument,Platform);

   /* read-in the in-situ netcdf file and extract the datasets of interest */
   /* go through the list of validation files (-v) and load the trajectories */
   size_t nbStations          =  0u;
   buoyTrajectory *validTrajs = NULL;

   unsigned int nbValFiles;
   char **ValFiles;
   fmstrsplit(listOfValFiles,",",&nbValFiles,&ValFiles);
   fmlogmsg(progname,"Read validation files:");
   for (size_t vf = 0 ; vf < nbValFiles ; vf++) {
      fprintf(stdout,"\t%s\n",ValFiles[vf]);
      ret = readValFile(ValFiles[vf],&nbStations,&validTrajs);
      if (ret) {
         fmerrmsg(progname,"While reading the validation file <%s>.",ValFiles[vf]);
         exit(OSISAF_ERROR_OTHER);
      }
   }

   printf("Found %u validation trajectories in input files\n",nbStations);
   /*
   for (size_t st = 0 ; st < nbStations ; st++) {
      printBuoyTrajectory_short(&(validTrajs[st]));
   }
   */


   /* read-in the osisaf product netcdf file and extract the datasets of interest */
   /* at this point, the datasets are the full maps. We do the geographical collocation 
    * at a later stage. */
   size_t prodDims[3];
   float  A[2],B[2];
   PJ    *pj;
   float *lat_start_fprod, *lon_start_fprod;
   fmsec1970 *t0_fprod, *t1_fprod;
   float *driftX_fprod, *driftY_fprod;
   short *flags_fprod;
   float *sigdX_fprod, *sigdY_fprod, *corrdXdY_fprod;
   short *uflags_fprod;
   float fillvalf_fprod;
   char  projstr[1028];
   memset(projstr,'\0',1028);
   sigdX_fprod = sigdY_fprod = NULL;
   if (!strcmp(Instrument,"multi")) {
      /* for a multi-sensor product, use the default uncertainty variables */
      sigdXname[0] = sigdYname[0] = '\0';
   } else {
      /* for a single-sensor product, use different uncertainty variables, 
       * depending on the temporal collocation scheme */
      if (temporal_collocation_scheme == USE2DCOL) {
         sprintf(sigdXname,"fsX_12utc");
         sprintf(sigdYname,"fsY_12utc");
      } else {
         sprintf(sigdXname,"fsX");
         sprintf(sigdYname,"fsY");
      }
   }
   ret = readProductFile(prodFile,refT0,refT1,prodDims,A,B,projstr,&pj,
         &t0_fprod,&t1_fprod,
         &lat_start_fprod, &lon_start_fprod,
         &driftX_fprod, &driftY_fprod,&flags_fprod,
         &fillvalf_fprod,
         &sigdX_fprod, &sigdY_fprod, &corrdXdY_fprod, &uflags_fprod);
   if (ret) {
      fmerrmsg(progname,"While reading the product file <%s>.",prodFile);
      exit(OSISAF_ERROR_OTHER);
   }
   if (fmt_type == LRSID_FILEFORMAT_FINL) {
      fmerrmsg(progname,"There is an issue with reading status_flag from final formatted file. Abort <%s>.",prodFile);
      exit(OSISAF_ERROR_OTHER);
   }
   free(uflags_fprod); /* we do not use it for now */
   /*
   if (!strcmp(levFilter,"lev3") && ((sigdX_fprod == NULL) || (sigdY_fprod == NULL))) {
      fmerrmsg(progname,"Something went wrong: we could not read uncertainties as <%s> and <%s>.",sigdXname,sigdYname);
      exit(OSISAF_ERROR_OTHER);
   }
   */


   /* Collocation. The aim here is to select/interpolate the OSISAF icedrift
    * at the same location (*and* time) than the buoy data.
    */
   fmlogmsg(progname,"Perform collocation against file:\n\t%s",prodFile);
   dbl_list *allMatchups = dbl_list_init(sizeof(matchProdVal));
   ret = matchProdWithVal(temporal_collocation_scheme,use_nearest_neighbour,refT0,refT1,
         prodDims,A,B,pj,lat_start_fprod,lon_start_fprod,t0_fprod,t1_fprod,driftX_fprod,driftY_fprod,flags_fprod,fillvalf_fprod,
         sigdX_fprod, sigdY_fprod, corrdXdY_fprod, 
         nbStations,validTrajs,
         allMatchups);
   if (ret) {
      fmerrmsg(progname,"While collocating the OSISAF product at validation locations.");
      exit(OSISAF_ERROR_OTHER);
   }

   /* Select only some (the best) collocation data pairs */
   ret = selectMatchups(allMatchups);
   if (ret) {
      fmerrmsg(progname,"While selecting best available matchups.");
      exit(OSISAF_ERROR_OTHER);
   }
   fprintf(stdout,"VERBOSE (%s) Found %u matchups!\n",progname, allMatchups->nbnodes);
   if (allMatchups->nbnodes == 0) {
      /* Found no matchups. This is not an error. */
      fmlogmsg(progname,"Exit now, before creating an empty matchup file.");
      exit(OSISAF_EXIT_CORRECT);
   }

   /* transform the list of matchups into several arrays */
   float *valB_lat;float *valB_lon;int *valB_time;float *valE_lat;float *valE_lon;int *valE_time;
   float *val_dX;float *val_dY;float *val_len;float *val_dir;char **val_id;char **val_netw;char **val_src;
   float *prdB_lat;float *prdB_lon;int *prdB_time;float *prdE_lat;float *prdE_lon;int *prdE_time;
   float *prd_dX;float *prd_dY;float *prd_len;float *prd_dir;short *prd_I;short *prd_J;short *prd_flag;
   float *prd_sX,*prd_sY,*prd_cXY;
   float *col_dist; float *col_dtB; float *col_dtE; float *col_dDur;
   size_t nbM;
   /* the 'time' variables are all expressed as seconds since refTM */
   fmsec1970 refTM = ymdh2fmsec1970("1978010100",0);
   ret = transformMatchupList(
         allMatchups,refTM,pj,
         &valB_lat,&valB_lon,&valB_time,&valE_lat,&valE_lon,&valE_time,
         &val_dX,&val_dY,&val_len,&val_dir,&val_id,&val_netw,&val_src,
         &prdB_lat,&prdB_lon,&prdB_time,&prdE_lat,&prdE_lon,&prdE_time,
         &prd_dX,&prd_dY,&prd_len,&prd_dir,&prd_I,&prd_J,&prd_flag,
         &prd_sX,&prd_sY,&prd_cXY,
         &col_dist,&col_dtB,&col_dtE,&col_dDur,
         &nbM);
   if (ret) {
      fmerrmsg(progname,"While transforming the list of matchups into arrays.");
      exit(OSISAF_ERROR_OTHER);
   }
   
   //dbl_list_print(allMatchups,&printMatchProdVal);

   /* write those matchups to a netCDF file */
   rsprod_file *matchupFile ;
   if (rsprod_file_create(&matchupFile,NULL,NULL)) {
      fprintf(stderr,"ERROR (%s) Did not manage to create the file object.\n",progname);
      exit(OSISAF_ERROR_OTHER);
   }
   ret = prepare_Matchups_nc(matchupFile,
         nbM,refTM,refT0,refT1,
         valB_lat,valB_lon,valB_time,valE_lat,valE_lon,valE_time,
         val_dX,val_dY,val_len,val_dir,val_id,val_netw,val_src,
         prdB_lat,prdB_lon,prdB_time,prdE_lat,prdE_lon,prdE_time,
         prd_dX,prd_dY,prd_len,prd_dir,prd_I,prd_J,prd_flag,
         prd_sX,prd_sY,prd_cXY,
         col_dist, col_dtB, col_dtE, col_dDur);
   if (ret) {
      fmerrmsg(progname,"While preparing matchups in a netCDF file object.");
      exit(OSISAF_ERROR_OTHER);
   }

   /* Add some extra global attributes */
   rsprod_attributes *globalAttributes = matchupFile->glob_attr;
   rsprod_attr *attr;
   ret = 0;
#define N_STR "history"
   /* get current date+time */
   fmtime fmtoday;
   time_t now;
   fmsec1970 today = time(&now);
   tofmtime(today,&fmtoday);
   char nowDate[64];
   sprintf(nowDate,"created %04d-%02d-%02d",fmtoday.fm_year,fmtoday.fm_mon,fmtoday.fm_mday);
   ret += rsprod_attr_createWithCopyValues(&attr,N_STR,RSPROD_CHAR,strlen(nowDate),nowDate) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
#undef N_STR
#define N_STR "satellite_product"
#define V_STR "OSI-405-"
   ret += rsprod_attr_createWithCopyValues(&attr,N_STR,RSPROD_CHAR,strlen(V_STR),V_STR) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr); 

   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot define further global attributes for the output file.\n",progname);
      return 1;
   }

   int ncid;
   /* open/create new netCDF file */
   ret = nc_create(outFile,NC_CLOBBER,&ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot open file for NetCDF writing\n\t%s\n",progname,outFile);
      exit(OSISAF_ERROR_OTHER);
   }

   /* write the file object to netCDF */
   librsprod_write_nullchar   = 0;
   librsprod_pad_strings_with = '~';
   ret = rsprod_file_writeToNetCDF(matchupFile,ncid);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot write to netCDF file\n\t%s\n",progname,outFile);
      exit(OSISAF_ERROR_OTHER);
   }

   rsprod_file_delete(matchupFile);

   /* close the netCDF file */
   ret = nc_close(ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot close file after NetCDF writing\n\t%s\n",progname,outFile);
      exit(OSISAF_ERROR_OTHER);
   }

   fmlogmsg(progname,"Done. Results are in file %s",outFile);


   exit(OSISAF_EXIT_CORRECT);
}

int readValFile(char *filename, 
      size_t *nbValidStations,
      buoyTrajectory **Trajs) {

   int ret;
   char errmsg[] = "readValFile: ";

   /* Open netCDF file */
   int ncid;
   ret = nc_open(filename,NC_NOCLOBBER,&ncid);check_ncstatus(ncid, ret, errmsg); 

   /* 
    *     LOAD THE FIELDS OF INTEREST
    */
   /* name: */
   rsprod_field *Fname;
   if (rsprod_field_loadFromNetCDF(&Fname,"id",ncid)) {
      fprintf(stderr,"%sWhile loading the 'id' dataset.\n",errmsg);
      return 1;
   }
   char *Fname_data;
   if (rsprod_data_getType(Fname->data) != RSPROD_CHAR) {
      fprintf(stderr,"%sThe 'id' dataset is not of type CHAR.\n",errmsg);
      return 1;
   }
   if ( (*Fname->data->methods->accessValues)(Fname->data,&Fname_data) )  {
      fprintf(stderr,"%serror in accessing the 'id' values.\n",errmsg);
      return 1;
   }
   size_t nbstations   = Fname->dims->length[0];
   size_t nbchars_name = Fname->dims->length[1];

   /* network: */
   rsprod_field *Fnetwork;
   if (rsprod_field_loadFromNetCDF(&Fnetwork,"network",ncid)) {
      fprintf(stderr,"%sWhile loading the 'network' dataset.\n",errmsg);
      return 1;
   }
   char *Fnetwork_data;
   if (rsprod_data_getType(Fnetwork->data) != RSPROD_CHAR) {
      fprintf(stderr,"%sThe 'network' dataset is not of type CHAR.\n",errmsg);
      return 1;
   }
   if ( (*Fnetwork->data->methods->accessValues)(Fnetwork->data,&Fnetwork_data) )  {
      fprintf(stderr,"%serror in accessing the 'network' values.\n",errmsg);
      return 1;
   }
   size_t nbchars_network = Fnetwork->dims->length[Fnetwork->dims->nbdims-1];

   /* source: */
   rsprod_field *Fsource;
   if (rsprod_field_loadFromNetCDF(&Fsource,"source",ncid)) {
      fprintf(stderr,"%sWhile loading the 'source' dataset.\n",errmsg);
      return 1;
   }
   char *Fsource_data;
   if (rsprod_data_getType(Fsource->data) != RSPROD_CHAR) {
      fprintf(stderr,"%sThe 'source' dataset is not of type CHAR.\n",errmsg);
      return 1;
   }
   if ( (*Fsource->data->methods->accessValues)(Fsource->data,&Fsource_data) )  {
      fprintf(stderr,"%serror in accessing the 'source' values.\n",errmsg);
      return 1;
   }
   size_t nbchars_source = Fsource->dims->length[Fsource->dims->nbdims-1];

   /* reformat the name and network string dataset information */
   char **Fname_data2,**Fnetwork_data2,**Fsource_data2;
   if (datastring_format(Fname_data,nbstations,nbchars_name,&Fname_data2)) {
      fprintf(stderr,"%serror in reformatting the 'name' values.\n",errmsg);
      return 1;
   }
   if (datastring_format(Fnetwork_data,nbstations,nbchars_network,&Fnetwork_data2) ) {
      fprintf(stderr,"%serror in reformatting the 'network' values.\n",errmsg);
      return 1;
   }
   if (datastring_format(Fsource_data,nbstations,nbchars_source,&Fsource_data2) ) {
      fprintf(stderr,"%serror in reformatting the 'source' values.\n",errmsg);
      return 1;
   }

   /* time: */
   rsprod_field *Ftime;
   if (rsprod_field_loadFromNetCDF(&Ftime,"time",ncid)) {
      fprintf(stderr,"%sWhile loading the 'time' dataset.\n",errmsg);
      return 1;
   }
   int *Ftime_data;
   if (rsprod_data_getType(Ftime->data) != RSPROD_INT) {
      fprintf(stderr,"%sThe 'time' dataset is not of type INT.\n",errmsg);
      return 1;
   }
   if ( (*Ftime->data->methods->accessValues)(Ftime->data,&Ftime_data) )  {
      fprintf(stderr,"%serror in accessing the 'time' values.\n",errmsg);
      return 1;
   }
   size_t maxnbtimes   = Ftime->dims->length[1];
   int timefillval;
   if ( rsprod_attributes_accessValue_Int(Ftime->attr,"_FillValue",&timefillval) ) { 
      fprintf(stderr,"%serror in accessing the 'time' _FillValue attribute.\n",errmsg);
      return 1;
   }

   /* lat: */
   rsprod_field *Flat;
   if (rsprod_field_loadFromNetCDF(&Flat,"lat",ncid)) {
      fprintf(stderr,"%sWhile loading the 'lat' dataset.\n",errmsg);
      return 1;
   }
   float *Flat_data;
   if (rsprod_data_getType(Flat->data) != RSPROD_FLOAT) {
      fprintf(stderr,"%sThe 'lon' dataset is not of type FLOAT.\n",errmsg);
      return 1;
   }
   if ( (*Flat->data->methods->accessValues)(Flat->data,&Flat_data) )  {
      fprintf(stderr,"%serror in accessing the 'lat' values.\n",errmsg);
      return 1;
   }

   /* lon: */
   rsprod_field *Flon;
   if (rsprod_field_loadFromNetCDF(&Flon,"lon",ncid)) {
      fprintf(stderr,"%sWhile loading the 'lon' dataset.\n",errmsg);
      return 1;
   }
   float *Flon_data;
   if (rsprod_data_getType(Flon->data) != RSPROD_FLOAT) {
      fprintf(stderr,"%sThe 'lon' dataset is not of type FLOAT.\n",errmsg);
      return 1;
   }
   if ( (*Flon->data->methods->accessValues)(Flon->data,&Flon_data) )  {
      fprintf(stderr,"%serror in accessing the 'lon' values.\n",errmsg);
      return 1;
   }

   /* get the reference_date_and_time attribute (truly start_date_and_time) */
   rsprod_attributes *globattr;
   ret = rsprod_attributes_loadFromNetCDF(&globattr,ncid,NC_GLOBAL);
   if (ret) {
      fprintf(stderr,"%serror cannot find and load grobal attributes.\n",errmsg);
      return 1;
   }
   char *isoRefDateAndTime;
   if ( rsprod_attributes_accessValues_Char(globattr,"start_date_and_time",&isoRefDateAndTime) ) { 
      fprintf(stderr,"%serror in accessing the grobal <start_date_and_time> attribute.\n",errmsg);
      return 1;
   }
   fmsec1970 refTime = CFepoch2fmsec1970(isoRefDateAndTime);

   /* Go through the netcdf arrays and fill in the array of trajectories */
   /* the array of trajectories we get from input is appended with the new ones */
   buoyTrajectory *trajectories = *Trajs;
   trajectories = realloc(trajectories,(*nbValidStations+nbstations) * sizeof(buoyTrajectory));
   if (!trajectories) {
      fprintf(stderr,"ERROR (%s) Problem with memory allocation\n",__func__);
      return 1;
   }
   for (size_t s = 0 ; s < nbstations ; s++) {
      buoyTrajectory *t = &(trajectories[*nbValidStations+s]);

      strcpy(t->id,Fname_data2[s]);
      strcpy(t->network,Fnetwork_data2[s]);
      strcpy(t->source,Fsource_data2[s]);
      t->records = dbl_list_init(sizeof(buoyRecord));
      long st = s*maxnbtimes;
      long rn = 0;
      while ( (Ftime_data[st+rn] != timefillval) && (rn < maxnbtimes) ) {
         //printf("%ld %d %d\n",st,Ftime_data[st+rn],timefillval);
         Obsdata *obs = fmMalloc(sizeof(*obs));

         strcpy(obs->stID,t->id);
         strcpy(obs->stType,t->network);
         strcpy(obs->stDesc,t->source);
         obs->lat   = Flat_data[st+rn];
         obs->lon   = Flon_data[st+rn];
         fmsec1970 obstime = Ftime_data[st+rn] + refTime; // get seconds since 1970-01-01
         fmtime fmt;
         tofmtime(obstime,&fmt);
         obs->year  = fmt.fm_year;
         obs->month = fmt.fm_mon;
         obs->day   = fmt.fm_mday;
         obs->hour  = fmt.fm_hour;
         obs->min   = fmt.fm_min;

         /*
            printf("%04d-%02d-%02d %02d:%02d:%02d %d %d (%ld)\n",
            obs->year,obs->month,obs->day,obs->hour,obs->min,0,
            Ftime_data[st+rn],obstime,refTime);
          */

         buoyRecord buoyRec;
         setBuoyRecord(&buoyRec,obs,obstime);
         if (dbl_list_addNode(t->records,dbl_node_createCopyContent(t->records,&buoyRec,NULL),LIST_POSITION_LAST)) {
            fprintf(stderr,"WARNING (%s) Did not manage to add the new node (%s).\n",__func__,"struct buoyRecord");
            continue;
         }

         rn++;
      }

   }


   /* Close netCDF file */
   ret = nc_close(ncid); check_ncstatus(ncid,ret,errmsg);

   /* Free the memory */
   rsprod_field_delete(Fname);
   rsprod_field_delete(Fnetwork);
   rsprod_field_delete(Ftime);
   rsprod_field_delete(Flat);rsprod_field_delete(Flon);

   /* Set the output parameters */ 
   *nbValidStations += nbstations;
   *Trajs            = trajectories;

   return 0;
}

int transformMatchupList(
      /* input */
      dbl_list *allMatchups,fmsec1970 refTime,PJ *pj,
      /* output VAL */
      float **valB_lat,float **valB_lon,int **valB_time,float **valE_lat,float **valE_lon,int **valE_time,
      float **val_dX,float **val_dY,float **val_len,float **val_dir,char ***val_id,char ***val_netw,char ***val_src,
      /* output PRD */
      float **prdB_lat,float **prdB_lon,int **prdB_time,float **prdE_lat,float **prdE_lon,int **prdE_time,
      float **prd_dX,float **prd_dY,float **prd_len,float **prd_dir,short **prd_I,short **prd_J,short **prd_flag,
      float **prd_sX,float **prd_sY,float **prd_cXY,
      /* output COL */
      float **col_dist, float **col_dtB, float **col_dtE, float **col_dDur,
      size_t *nbM) {

   /* number of matchups to be processed */
   size_t nbMatch = allMatchups->nbnodes;
   if (!nbMatch) {
      fprintf(stderr,"WARNING (%s) List of matchups is empty.\n",__func__);
      return 0; /* not an error */
   }

   /* go through the list of matchups and get the 
    * maximum length of some string arrays */
   size_t maxnamelength  = 0;
   size_t maxnetwlength  = 0;
   size_t maxsrclength   = 0;
   dbl_node *node = allMatchups->head;
   do {
      matchProdVal *M = node->c;
      size_t namelength = strlen(M->id); 
      size_t netwlength = strlen(M->network); 
      size_t srclength  = strlen(M->source); 
      node = node->n;
      if (namelength > maxnamelength)
         maxnamelength = namelength;
      if (netwlength > maxnetwlength)
         maxnetwlength = netwlength;
      if (srclength  > maxsrclength)
         maxsrclength = srclength;
   } while (node != allMatchups->head);

   /* allocate all output arrays */
   *valB_lat  = fmMalloc(sizeof(**valB_lat)*nbMatch);
   *valB_lon  = fmMalloc(sizeof(**valB_lon)*nbMatch);
   *valB_time = fmMalloc(sizeof(**valB_time)*nbMatch);
   *valE_lat  = fmMalloc(sizeof(**valE_lat)*nbMatch);
   *valE_lon  = fmMalloc(sizeof(**valE_lon)*nbMatch);
   *valE_time = fmMalloc(sizeof(**valE_time)*nbMatch);
   *val_dX    = fmMalloc(sizeof(**val_dX)*nbMatch);
   *val_dY    = fmMalloc(sizeof(**val_dY)*nbMatch);
   *val_len   = fmMalloc(sizeof(**val_len)*nbMatch);
   *val_dir   = fmMalloc(sizeof(**val_dir)*nbMatch);
   *val_id    = fmMalloc(sizeof(**val_id)*nbMatch);
   *val_netw  = fmMalloc(sizeof(**val_netw)*nbMatch);
   *val_src   = fmMalloc(sizeof(**val_src)*nbMatch);
   *prdB_lat  = fmMalloc(sizeof(**prdB_lat)*nbMatch);
   *prdB_lon  = fmMalloc(sizeof(**prdB_lon)*nbMatch);
   *prdB_time = fmMalloc(sizeof(**prdB_time)*nbMatch);
   *prdE_lat  = fmMalloc(sizeof(**prdE_lat)*nbMatch);
   *prdE_lon  = fmMalloc(sizeof(**prdE_lon)*nbMatch);
   *prdE_time = fmMalloc(sizeof(**prdE_time)*nbMatch);
   *prd_dX    = fmMalloc(sizeof(**prd_dX)*nbMatch);
   *prd_dY    = fmMalloc(sizeof(**prd_dY)*nbMatch);
   *prd_sX    = fmMalloc(sizeof(**prd_sX)*nbMatch);
   *prd_sY    = fmMalloc(sizeof(**prd_sY)*nbMatch);
   *prd_cXY   = fmMalloc(sizeof(**prd_cXY)*nbMatch);
   *prd_len   = fmMalloc(sizeof(**prd_len)*nbMatch);
   *prd_dir   = fmMalloc(sizeof(**prd_dir)*nbMatch);
   *prd_flag  = fmMalloc(sizeof(**prd_flag)*nbMatch);
   *col_dist  = fmMalloc(sizeof(**col_dist)*nbMatch);
   *col_dtB   = fmMalloc(sizeof(**col_dtB)*nbMatch);
   *col_dtE   = fmMalloc(sizeof(**col_dtE)*nbMatch);
   *col_dDur  = fmMalloc(sizeof(**col_dDur)*nbMatch);

   /* allocate and initialize charater arrays */
   for (size_t m = 0 ; m < nbMatch ; m++) {
      (*val_id)[m]  = fmMalloc(maxnamelength+1);
      memset((*val_id)[m],(int)'\0',maxnamelength+1);
      (*val_netw)[m] = fmMalloc(maxnetwlength+1);
      memset((*val_netw)[m],(int)'\0',maxnetwlength+1);
      (*val_src)[m] = fmMalloc(maxsrclength+1);
      memset((*val_src)[m],(int)'\0',maxsrclength+1);
   }

   /* fill output arrays with the content of the list */
   size_t m = 0;
   /*dbl_node * */node = allMatchups->head;
   do {
      matchProdVal *M = node->c;

      /* some intermediate variables */
      double xB,yB,xE,yE,dlen,ddir;

      /* VAL */
      strcpy((*val_id)[m],M->id);
      strcpy((*val_netw)[m],M->network);
      strcpy((*val_src)[m],M->source);
      (*valB_lat)[m]  = (M->val[0]).data.lat;
      (*valB_lon)[m]  = (M->val[0]).data.lon;
      (*valB_time)[m] = (M->val[0]).time - refTime;
      (*valE_lat)[m]  = (M->val[1]).data.lat;
      (*valE_lon)[m]  = (M->val[1]).data.lon;
      (*valE_time)[m] = (M->val[1]).time - refTime;
      float val_Dur   = (*valE_time)[m] - (*valB_time)[m]; 
      remap_ll2xy((*valB_lat)[m],(*valB_lon)[m],pj,1.,0,1.,0.,&xB,&yB,0);
      remap_ll2xy((*valE_lat)[m],(*valE_lon)[m],pj,1.,0,1.,0.,&xE,&yE,0);
      (*val_dX)[m]    = xE - xB;
      (*val_dY)[m]    = yE - yB;
      compute_distance((*valB_lat)[m],(*valB_lon)[m],(*valE_lat)[m],(*valE_lon)[m],&dlen);
      compute_directionToNorth((*valB_lat)[m],(*valB_lon)[m],(*valE_lat)[m],(*valE_lon)[m],&ddir);
      (*val_len)[m]   = dlen;
      (*val_dir)[m]   = ddir;

      /* PRD */ 
      (*prdB_lat)[m]  = M->latB_prod;
      (*prdB_lon)[m]  = M->lonB_prod;
      (*prdB_time)[m] = M->tB_prod - refTime;
      (*prdE_lat)[m]  = M->latE_prod;
      (*prdE_lon)[m]  = M->lonE_prod;
      (*prdE_time)[m] = M->tE_prod - refTime;
      float prd_Dur   = (*prdE_time)[m] - (*prdB_time)[m];
      remap_ll2xy((*prdB_lat)[m],(*prdB_lon)[m],pj,1.,0,1.,0.,&xB,&yB,0);
      remap_ll2xy((*prdE_lat)[m],(*prdE_lon)[m],pj,1.,0,1.,0.,&xE,&yE,0);
      (*prd_dX)[m]    = xE - xB;
      (*prd_dY)[m]    = yE - yB;
      compute_distance((*prdB_lat)[m],(*prdB_lon)[m],(*prdE_lat)[m],(*prdE_lon)[m],&dlen);
      compute_directionToNorth((*prdB_lat)[m],(*prdB_lon)[m],(*prdE_lat)[m],(*prdE_lon)[m],&ddir);
      (*prd_len)[m]   = dlen;
      (*prd_dir)[m]   = ddir;
      (*prd_sX)[m]    = M->sX_prod;
      (*prd_sY)[m]    = M->sY_prod;
      (*prd_cXY)[m]   = M->cXY_prod;
      (*prd_flag)[m]  = M->flag_prod;

      /* COLLOCATION */
      double cdist;
      compute_distance((*valB_lat)[m],(*valB_lon)[m],(*prdB_lat)[m],(*prdB_lon)[m],&cdist);
      (*col_dist)[m]  = cdist;
      /* all collocation temporal stats are kept as hours */
      (*col_dtB)[m]   = ((*valB_time)[m] - (*prdB_time)[m]) / 60. / 60.;
      (*col_dtE)[m]   = ((*valE_time)[m] - (*prdE_time)[m]) / 60. / 60.;
      (*col_dDur)[m]  = (val_Dur - prd_Dur) / 60. / 60.;

      node = node->n;m++;
   } while (node != allMatchups->head);

   *nbM = m;
   return 0;
}


int prepare_Matchups_nc(rsprod_file *matchupFile,
      size_t nbMatch,
      fmsec1970 refTime_matchup, fmsec1970 refTime0, fmsec1970 refTime1,
      float *valB_lat,float *valB_lon,int *valB_time,float *valE_lat,float *valE_lon,int *valE_time,
      float *val_dX,float *val_dY,float *val_len,float *val_dir,char **val_id,char **val_netw,char **val_src,
      float *prdB_lat,float *prdB_lon,int *prdB_time,float *prdE_lat,float *prdE_lon,int *prdE_time,
      float *prd_dX,float *prd_dY,float *prd_len,float *prd_dir,short *prd_I,short *prd_J,
      short *prd_flag,
      float *prd_sX,float *prd_sY,float *prd_cXY,
      float *col_dist, float *col_dtB, float *col_dtE, float *col_dDur) {
   
   int ret;

   /* dimension object: matchup */
   rsprod_dimensions *dims;
   {
      char *dimnames[1];
      dimnames[0] = fmMalloc(strlen("matchup")+1); strcpy(dimnames[0],"matchup");
      size_t dimlengths[1];
      dimlengths[0]  = nbMatch;
      short  dimorders[1];
      dimorders[0] = 1;
      ret = rsprod_dims_create(&dims,1,dimnames,dimlengths,dimorders);
      if (ret) {
         fprintf(stderr,"ERROR (%s) could not create the Dimension object for output file.\n",__func__);
         return 1;
      }
   }

   float fillval_f = -999999.;
   short fillval_s = -99;

   char time_unit[100];
   char refTime_string[FMUTIL_CFEPOCH_LENGTH+1]; 
   fmsec19702CFepoch(refTime_matchup,refTime_string);
   sprintf(time_unit,"seconds since %s",refTime_string);
   /* fields */

   /* ===================================== */
   /* FIELDS DESCRIBING THE VALIDATION DATA */
   /* ===================================== */
   rsprod_field *Fval_id;
   if (rsprod_field_createStandard(&Fval_id,"val_id",RSPROD_CHAR,nbMatch,dims,"id for validation data",NULL,NULL,
            NULL,NULL,NULL,0,NULL,val_id)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_id");
      return 1;
   }
   rsprod_field *Fval_netw;
   if (rsprod_field_createStandard(&Fval_netw,"val_network",RSPROD_CHAR,nbMatch,dims,"network for validation data",NULL,NULL,
            NULL,NULL,NULL,0,NULL,val_netw)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_netw");
      return 1;
   }
   rsprod_field *Fval_src;
   if (rsprod_field_createStandard(&Fval_src,"val_source",RSPROD_CHAR,nbMatch,dims,"geo-position technique for validation data",NULL,NULL,
            NULL,NULL,NULL,0,NULL,val_src)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_src");
      return 1;
   }
   rsprod_field *FvalB_lat;
   if (rsprod_field_createStandard(&FvalB_lat,"valB_lat",RSPROD_FLOAT,nbMatch,dims,"latitude at start of drift for validation data",
            NULL,"degrees north",
            &fillval_f,NULL,NULL,0,NULL,valB_lat)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"valB_lat");
      return 1;
   }
   rsprod_field *FvalB_lon;
   if (rsprod_field_createStandard(&FvalB_lon,"valB_lon",RSPROD_FLOAT,nbMatch,dims,"longitude at start of drift for validation data",
            NULL,"degrees east",
            &fillval_f,NULL,NULL,0,NULL,valB_lon)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"valB_lon");
      return 1;
   }
   rsprod_field *FvalB_time;
   if (rsprod_field_createStandard(&FvalB_time,"valB_time",RSPROD_INT,nbMatch,dims,"time start of drift for validation data",
            NULL,time_unit,
            &fillval_f,NULL,NULL,0,NULL,valB_time)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"valB_time");
      return 1;
   }
   rsprod_field *FvalE_lat;
   if (rsprod_field_createStandard(&FvalE_lat,"valE_lat",RSPROD_FLOAT,nbMatch,dims,"latitude at start of drift for validation data",
            NULL,"degrees north",
            &fillval_f,NULL,NULL,0,NULL,valE_lat)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"valE_lat");
      return 1;
   }
   rsprod_field *FvalE_lon;
   if (rsprod_field_createStandard(&FvalE_lon,"valE_lon",RSPROD_FLOAT,nbMatch,dims,"longitude at start of drift for validation data",
            NULL,"degrees east",
            &fillval_f,NULL,NULL,0,NULL,valE_lon)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"valE_lon");
      return 1;
   }
   rsprod_field *FvalE_time;
   if (rsprod_field_createStandard(&FvalE_time,"valE_time",RSPROD_INT,nbMatch,dims,"time stop of drift for validation data",
            NULL,time_unit,
            &fillval_f,NULL,NULL,0,NULL,valE_time)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"valE_time");
      return 1;
   }
   rsprod_field *Fval_dX;
   if (rsprod_field_createStandard(&Fval_dX,"val_dX",RSPROD_FLOAT,nbMatch,dims,"x-component of validation drift vector",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,val_dX)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_dX");
      return 1;
   }
   rsprod_field *Fval_dY;
   if (rsprod_field_createStandard(&Fval_dY,"val_dY",RSPROD_FLOAT,nbMatch,dims,"y-component of validation drift vector",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,val_dY)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_dY");
      return 1;
   }
   rsprod_field *Fval_len;
   if (rsprod_field_createStandard(&Fval_len,"val_len",RSPROD_FLOAT,nbMatch,dims,"magnitude of validation drift vector",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,val_len)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_len");
      return 1;
   }
   rsprod_field *Fval_dir;
   if (rsprod_field_createStandard(&Fval_dir,"val_dir",RSPROD_FLOAT,nbMatch,dims,"direction of validation drift vector",
            NULL,"degrees to north",
            &fillval_f,NULL,NULL,0,NULL,val_dir)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"val_dir");
      return 1;
   }

   /* ============================= */
   /* FIELDS DESCRIBING THE PRODUCT */
   /* ============================= */
   rsprod_field *FprdB_lat;
   if (rsprod_field_createStandard(&FprdB_lat,"prdB_lat",RSPROD_FLOAT,nbMatch,dims,"latitude at start of drift for osisaf data",
            NULL,"degrees north",
            &fillval_f,NULL,NULL,0,NULL,prdB_lat)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prdB_lat");
      return 1;
   }
   rsprod_field *FprdB_lon;
   if (rsprod_field_createStandard(&FprdB_lon,"prdB_lon",RSPROD_FLOAT,nbMatch,dims,"longitude at start of drift for osisaf data",
            NULL,"degrees east",
            &fillval_f,NULL,NULL,0,NULL,prdB_lon)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prdB_lon");
      return 1;
   }
   rsprod_field *FprdB_time;
   if (rsprod_field_createStandard(&FprdB_time,"prdB_time",RSPROD_INT,nbMatch,dims,"time start of drift for osisaf data",
            NULL,time_unit,
            &fillval_f,NULL,NULL,0,NULL,prdB_time)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prdB_time");
      return 1;
   }
   rsprod_field *FprdE_lat;
   if (rsprod_field_createStandard(&FprdE_lat,"prdE_lat",RSPROD_FLOAT,nbMatch,dims,"latitude at start of drift for osisaf data",
            NULL,"degrees north",
            &fillval_f,NULL,NULL,0,NULL,prdE_lat)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prdE_lat");
      return 1;
   }
   rsprod_field *FprdE_lon;
   if (rsprod_field_createStandard(&FprdE_lon,"prdE_lon",RSPROD_FLOAT,nbMatch,dims,"longitude at start of drift for osisaf data",
            NULL,"degrees east",
            &fillval_f,NULL,NULL,0,NULL,prdE_lon)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prdE_lon");
      return 1;
   }
   rsprod_field *FprdE_time;
   if (rsprod_field_createStandard(&FprdE_time,"prdE_time",RSPROD_INT,nbMatch,dims,"time stop of drift for osisaf data",
            NULL,time_unit,
            &fillval_f,NULL,NULL,0,NULL,prdE_time)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prdE_time");
      return 1;
   }
   rsprod_field *Fprd_dX;
   if (rsprod_field_createStandard(&Fprd_dX,"prd_dX",RSPROD_FLOAT,nbMatch,dims,"x-component of osisaf drift vector",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,prd_dX)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_dX");
      return 1;
   }
   rsprod_field *Fprd_dY;
   if (rsprod_field_createStandard(&Fprd_dY,"prd_dY",RSPROD_FLOAT,nbMatch,dims,"y-component of osisaf drift vector",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,prd_dY)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_dY");
      return 1;
   }
   rsprod_field *Fprd_len;
   if (rsprod_field_createStandard(&Fprd_len,"prd_len",RSPROD_FLOAT,nbMatch,dims,"magnitude of osisaf drift vector",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,prd_len)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_len");
      return 1;
   }
   rsprod_field *Fprd_dir;
   if (rsprod_field_createStandard(&Fprd_dir,"prd_dir",RSPROD_FLOAT,nbMatch,dims,"direction of osisaf drift vector",
            NULL,"degrees to north",
            &fillval_f,NULL,NULL,0,NULL,prd_dir)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_dir");
      return 1;
   }
   rsprod_field *Fprd_flag;
   if (rsprod_field_createStandard(&Fprd_flag,"prd_flag",RSPROD_SHORT,nbMatch,dims,"status_flag of osisaf drift vector",
            NULL,"1",
            &fillval_s,NULL,NULL,0,NULL,prd_flag)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_flag");
      return 1;
   }
   rsprod_field *Fprd_sX;
   if (rsprod_field_createStandard(&Fprd_sX,"prd_sX",RSPROD_FLOAT,nbMatch,dims,"x-component of osisaf drift uncertainty",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,prd_sX)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_sX");
      return 1;
   }
   rsprod_field *Fprd_sY;
   if (rsprod_field_createStandard(&Fprd_sY,"prd_sY",RSPROD_FLOAT,nbMatch,dims,"y-component of osisaf drift uncertainty",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,prd_sY)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_sY");
      return 1;
   }
   rsprod_field *Fprd_cXY;
   if (rsprod_field_createStandard(&Fprd_cXY,"prd_cXY",RSPROD_FLOAT,nbMatch,dims,
            "X-Y correlation measure of drift uncertainty",
            NULL,"1",
            &fillval_f,NULL,NULL,0,NULL,prd_cXY)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"prd_cXY");
      return 1;
   }
   /* ===================================== */
   /* FIELDS DESCRIBING THE COLLOCATION     */
   /* ===================================== */
   //float *col_dist, float *col_dtB, float *col_dtE
   rsprod_field *Fcol_dist;
   if (rsprod_field_createStandard(&Fcol_dist,"col_dist",RSPROD_FLOAT,nbMatch,dims,"distance from REF to PROD at start point",
            NULL,"km",
            &fillval_f,NULL,NULL,0,NULL,col_dist)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"col_dist");
      return 1;
   }
   rsprod_field *Fcol_dtB;
   if (rsprod_field_createStandard(&Fcol_dtB,"col_dtB",RSPROD_FLOAT,nbMatch,dims,"diffence in start time T0 (REF minus PROD)",
            NULL,"hours",
            &fillval_f,NULL,NULL,0,NULL,col_dtB)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"col_dtB");
      return 1;
   }
   rsprod_field *Fcol_dtE;
   if (rsprod_field_createStandard(&Fcol_dtE,"col_dtE",RSPROD_FLOAT,nbMatch,dims,"diffence in stop time T1 (REF minus PROD)",
            NULL,"hours",
            &fillval_f,NULL,NULL,0,NULL,col_dtE)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"col_dtE");
      return 1;
   }
   rsprod_field *Fcol_dDur;
   if (rsprod_field_createStandard(&Fcol_dDur,"col_dDur",RSPROD_FLOAT,nbMatch,dims,"diffence in duration (T1-T0) (REF minus PROD)",
            NULL,"hours",
            &fillval_f,NULL,NULL,0,NULL,col_dDur)) {
      fprintf(stderr,"ERROR (%s) Cannot create the %s FIELD object.\n",__func__,"col_dDur");
      return 1;
   }


   /* define some global attributes */
   rsprod_attributes *globalAttributes;
   rsprod_attributes_create(&globalAttributes);
   ret = 0;
   rsprod_attr *attr;
#define N_STR "title"
#define V_STR "Matchup data base for validation of OSI SAF sea ice drift products"
   ret += rsprod_attr_createWithCopyValues(&attr,N_STR,RSPROD_CHAR,strlen(V_STR),V_STR) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr); 
#undef N_STR
#undef V_STR
#define N_STR "start_date_and_time"
   char start_date_cf[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec19702CFepoch(refTime0,start_date_cf);
   ret += rsprod_attr_createWithCopyValues(&attr,N_STR,RSPROD_CHAR,strlen(start_date_cf),start_date_cf) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr); 
#undef N_STR
#define N_STR "end_date_and_time"
   char end_date_cf[FMUTIL_CFEPOCH_LENGTH+1];
   fmsec19702CFepoch(refTime1,end_date_cf);
   ret += rsprod_attr_createWithCopyValues(&attr,N_STR,RSPROD_CHAR,strlen(end_date_cf),end_date_cf) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr); 
#undef N_STR
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot define the global attribute objects\n",__func__);
      return 1;
   }

   /* assign the global attributes to the file object */
   matchupFile->glob_attr = globalAttributes;

   /* assign all the fields to the file object */ 
   ret  = 0;
   ret += rsprod_file_addDataset(matchupFile,Fval_id);
   ret += rsprod_file_addDataset(matchupFile,Fval_netw);
   ret += rsprod_file_addDataset(matchupFile,Fval_src);
   ret += rsprod_file_addDataset(matchupFile,FvalB_lat);
   ret += rsprod_file_addDataset(matchupFile,FvalB_lon);
   ret += rsprod_file_addDataset(matchupFile,FvalB_time);
   ret += rsprod_file_addDataset(matchupFile,FvalE_lat);
   ret += rsprod_file_addDataset(matchupFile,FvalE_lon);
   ret += rsprod_file_addDataset(matchupFile,FvalE_time);
   ret += rsprod_file_addDataset(matchupFile,Fval_dX);
   ret += rsprod_file_addDataset(matchupFile,Fval_dY);
   ret += rsprod_file_addDataset(matchupFile,Fval_len);
   ret += rsprod_file_addDataset(matchupFile,Fval_dir);
   ret += rsprod_file_addDataset(matchupFile,FprdB_lat);
   ret += rsprod_file_addDataset(matchupFile,FprdB_lon);
   ret += rsprod_file_addDataset(matchupFile,FprdB_time);
   ret += rsprod_file_addDataset(matchupFile,FprdE_lat);
   ret += rsprod_file_addDataset(matchupFile,FprdE_lon);
   ret += rsprod_file_addDataset(matchupFile,FprdE_time);
   ret += rsprod_file_addDataset(matchupFile,Fprd_dX);
   ret += rsprod_file_addDataset(matchupFile,Fprd_dY);
   ret += rsprod_file_addDataset(matchupFile,Fprd_sX);
   ret += rsprod_file_addDataset(matchupFile,Fprd_sY);
   ret += rsprod_file_addDataset(matchupFile,Fprd_cXY);
   ret += rsprod_file_addDataset(matchupFile,Fprd_len);
   ret += rsprod_file_addDataset(matchupFile,Fprd_dir);
   ret += rsprod_file_addDataset(matchupFile,Fcol_dist);
   ret += rsprod_file_addDataset(matchupFile,Fprd_flag);
   ret += rsprod_file_addDataset(matchupFile,Fcol_dtB);
   ret += rsprod_file_addDataset(matchupFile,Fcol_dtE);
   ret += rsprod_file_addDataset(matchupFile,Fcol_dDur);
   if (ret) {
      fprintf(stderr,"ERROR (%s) could not add the datasets to the file object.\n",__func__);
      return 1;
   } 

   return 0;
}


