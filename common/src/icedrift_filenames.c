

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fmutil.h>
#include <icedrift_instruments.h>
#include <icedrift_finalNCformat_conventions.h>
#include "icedrift_filenames.h"

/* two levels of separators (global variables) */
char Fsep[] = "_";
char Ssep[] = "-";   
char Tsep[] = "+";   

char ncExt[] = "nc";
char tcimage_prefix[] = "tc";
static char tcimage_tweight[] = "wght";

static char driftproduct_prefix[] = "icedrift";
static char mdriftproduct_id[]    = "multi";
static char tavgdrift_prefix[]    = "tavg";
static char dformproduct_prefix[] = "icedform";

int isValid_InstrumentAndPlatform(char InstrumentAndPlatform[],char **Instrument, char **Platform) {

   int ret;

   /* split on the separator */
   unsigned int nbTokens;
   char **tokens;
   fmstrsplit(InstrumentAndPlatform,Ssep,&nbTokens,&tokens);
   if ( (nbTokens != 2) ) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted Instrument+Platform string <%s>.\n",__func__,InstrumentAndPlatform);
      return 0;
   } 
   *Instrument = tokens[0];
   *Platform   = tokens[1];
   /* check that we know the instrument and platform */
   if (!compatibleInstrumentAndPlatform(*Instrument,*Platform)) {
      fprintf(stderr,"ERROR (%s) Incompatible Instrument (%s) and Platform (%s).\n",__func__,*Instrument,*Platform);
      return 0;
   } 
   free(tokens);

   return 1;
}

int isValid_InstrumentAndPlatformAndChannel(char InstrumentAndPlatformAndChannel[],
      char **Instrument, char **Platform,char **Channels) {

   int ret;

   /* split on the separator */
   unsigned int nbTokens;
   char **tokens;
   fmstrsplit(InstrumentAndPlatformAndChannel,Ssep,&nbTokens,&tokens);
   if ( (nbTokens != 3) && (nbTokens != 2)) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted Instrument+Platform+Channel string <%s>.\n",__func__,InstrumentAndPlatformAndChannel);
      return 0;
   } 
   *Instrument = tokens[0];
   *Platform   = tokens[1];
   if (nbTokens == 3) {
      *Channels   = tokens[2];
   } else {
      *Channels   = NULL;
   }
   /* check that we know the instrument and platform */
   if (!compatibleInstrumentAndPlatform(*Instrument,*Platform)) {
      fprintf(stderr,"ERROR (%s) Incompatible Instrument (%s) and Platform (%s).\n",__func__,*Instrument,*Platform);
      return 0;
   } 
   free(tokens);

   return 1;
}


void tcimage_build_filename(short Tweight,char *Instrument, char *Platform, char *WaveBands, char *area, char *centralTime,char **fname) {
   char weightStr[10];
   sprintf(weightStr,"%s%s",(Tweight?"":"no"),tcimage_tweight);
   size_t nbchars = strlen(tcimage_prefix) + 1 + strlen(weightStr) + 1 + strlen(Instrument) + 1 + strlen(Platform) + 1 + 
      strlen(WaveBands) + 1 +
      strlen(area) + 1 + strlen(centralTime) + 1 + strlen(ncExt);

   *fname = fmMalloc(nbchars+1);
   sprintf(*fname,"%s%s""%s%s""%s%s%s%s""%s%s""%s%s""%s.%s",
	 tcimage_prefix,Fsep,
	 weightStr,Fsep,
	 Instrument,Ssep,Platform,Fsep,
	 WaveBands,Fsep,
	 area,Fsep,
	 centralTime,ncExt);
}

int tcimage_split_filename(char *fname, short *Tweight, char **Instrument, char **Platform, char **WaveBands, char **area, char **centralTime) {

   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   /* split on the separator */
   unsigned int nbTokens;
   char **tokens;
   fmstrsplit(basename,Fsep,&nbTokens,&tokens);
   if ( (nbTokens != 6) || (strcmp(tokens[0],tcimage_prefix)) ) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted file name <%s>.\n",__func__,basename);
      return 1;
   }

   /* The last token contains the file extension (.nc) which should be removed */
   char *lastdot = strrchr(tokens[nbTokens-1],'.');
   *lastdot = '\0';

   /* Further break the Intrument and Platform */
   char **tokens2;
   fmstrsplit(tokens[2],Ssep,&nbTokens,&tokens2);
   if (nbTokens != 2) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted Instrument+Platform name <%s>.\n",__func__,tokens[2]);
      return 1;
   }

   /*  Assign the parameters to the tokens */
   *Instrument  = tokens2[0];
   *Platform    = tokens2[1];
   *WaveBands   = tokens[3];
   *area        = tokens[4];
   *centralTime = tokens[5];

   *Tweight = 0;
   if (!strcmp(tokens[1],"wght")) *Tweight = 1;

   /* release the memory holding the array of array of chars (the strings are still allocated) */
   free(tokens);
   free(tokens2);

   return 0;
}
int is_tcimage_filename(char *fname) {
   
   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   return (!strncmp(basename,tcimage_prefix,strlen(tcimage_prefix)) && 
	 strstr(basename,".nc"));
}

char *driftproduct_build_description(char *Instrument, char *Platform, char *area, char *WaveBands, char *method, char *flevel, 
      char *startdate, short starttwght, char *enddate, short endtwght) {

   char *wavelength_info;
   if (strcmp(Instrument,"multi")) {
      if (!WaveBands) {
         fprintf(stderr,"ERROR (%s) The WaveBand information is NULL for non multi product.\n",__func__);
         return NULL;
      }
      wavelength_info = fmMalloc(strlen(WaveBands)+1+1);
      sprintf(wavelength_info,"%s%s",WaveBands,Fsep);
   } else {
      wavelength_info = NULL;
   } 

   char *flevel_info;
   if (strcmp(Instrument,"multi")) {
      if (!flevel) {
         fprintf(stderr,"ERROR (%s) The Flevel information is NULL for non multi product.\n",__func__);
         return NULL;
      }
      flevel_info = fmMalloc(strlen(flevel)+1+1);
      sprintf(flevel_info,"%s%s",flevel,Fsep);
   } else {
      flevel_info = NULL;
   } 

   char *descr = fmMalloc(sizeof(char)*200);
   sprintf(descr,"%s%s%s%s%s%s%s%s%s%s%s%c%s%s%c",
	 Instrument,Ssep,Platform,Fsep,
	 (wavelength_info?wavelength_info:""),
	 method,Fsep,
	 (flevel_info?flevel_info:""),
	 area,Fsep,
	 startdate,(starttwght?'w':'n'),Fsep,
	 enddate,(endtwght?'w':'n'));

   return descr;
}

void driftproduct_build_filename_final(char *Prefix, char *Area, char *GridInfo, char *InstrumentAndPlatform, 
      char *startdate, char *enddate, char FFsep, char FSsep, char **fname) {

   *fname = fmMalloc(1024);

   sprintf(*fname,"%s%c%s%c%s%c%s%c%s00%c%s00.nc",
     Prefix,FFsep,
     Area,FFsep,
     GridInfo,FFsep,
     InstrumentAndPlatform,FFsep,
     startdate,FSsep,
     enddate);

}

void driftproduct_build_filename_proc(char *Instrument, char *Platform, char *area, char *WaveBands, char *method, char *flevel, 
      char *startdate, short starttwght, char *enddate, short endtwght, short is_temporal_average, char **fname) {

   char *descr = driftproduct_build_description(Instrument,Platform,area,WaveBands,method,flevel,
         startdate,starttwght,enddate,endtwght);
   if (!descr) {
      fprintf(stderr,"ERROR (%s) Cannot build the product description string.\n",__func__);
      *fname = NULL;
      return;
   } 

   char *timeprefix_info = NULL;
   if (is_temporal_average) {
      timeprefix_info = fmMalloc(strlen(tavgdrift_prefix)+1+1);
      sprintf(timeprefix_info,"%s%s",tavgdrift_prefix,Fsep);
   }

   *fname = fmMalloc((timeprefix_info?strlen(timeprefix_info):0)+
         strlen(driftproduct_prefix)+strlen(Fsep)+strlen(descr)+1+strlen(ncExt)+1);
   sprintf(*fname,"%s%s%s%s.%s",(timeprefix_info?timeprefix_info:""),
         driftproduct_prefix,Fsep,descr,ncExt);
   free(descr);

}

int is_driftproduct_filename(char *fname, int *fmt_type) {
   
   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   /* check if the file is netcdf */
   if (!strstr(basename,".nc"))
      return 0;

   /* check if the file is an icedrift product with 'proc' format */
   if (!strncmp(basename,driftproduct_prefix,strlen(driftproduct_prefix))) {
      *fmt_type = LRSID_FILEFORMAT_PROC;
      return 1;
   }

   /* check if the file is an icedrift product with 'final' format */
   if (!strncmp(basename,LRSID_FILEN_SUFFIX,strlen(LRSID_FILEN_SUFFIX))) {
      *fmt_type = LRSID_FILEFORMAT_FINL;
      return 1;
   }

   return 0;

}

int driftproduct_split_filename_proc(char *fname, char **Instrument, char **Platform, 
      char **WaveBands,char **method,char **levFilter,char **area, char **startdate, char **enddate) {
   
   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   /* split on the separator */
   unsigned int nbTokens;
   char **tokens;
   fmstrsplit(basename,Fsep,&nbTokens,&tokens);

   if ( (nbTokens != 8) || (strcmp(tokens[0],driftproduct_prefix)) ) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted file name <%s>.\n",__func__,basename);
      return 1;
   }

   /* The last token contains the file extension (.nc) which should be removed */
   char *lastdot = strrchr(tokens[nbTokens-1],'.');
   *lastdot = '\0';

   /* Further break the Intrument and Platform */
   char **tokens2;
   fmstrsplit(tokens[1],Ssep,&nbTokens,&tokens2);
   if (nbTokens != 2) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted Instrument+Platform name <%s>.\n",__func__,tokens[1]);
      return 1;
   }
   /*  Assign the parameters to the tokens */
   *Instrument  = tokens2[0];
   *Platform    = tokens2[1];
   *WaveBands   = tokens[2];
   *method      = tokens[3];
   /* room for the filtering level indication */
   *levFilter   = tokens[4];
   *area        = tokens[5];
   *startdate   = tokens[6];
   *enddate     = tokens[7];

   /* the startdate and enddate are to be processed further has they have an extra character to
    * indicate if the time composite image was weighted ('w') or not ('n') */
   short startTweight = 0;
   if (strlen(*startdate) == 11) {
      char *weightChar = &((*startdate)[10]);
      if (*weightChar == 'w') {
	 startTweight = 1;
      }
      (*startdate)[10] = '\0';
   }

   short endTweight = 0;
   if (strlen(*enddate) == 11) {
      char *weightChar = &((*enddate)[10]);
      if (*weightChar == 'w') {
	 endTweight = 1;
      }
      (*enddate)[10] = '\0';
   }
   /* release the memory holding the array of array of chars (the strings are still allocated) */
   free(tokens);

   return 0;

}

int driftproduct_split_filename_final(char *fname, char **Area, char **GridInfo, 
      char **InstrumentAndPlatform, char **startdate, char **enddate) {
      
   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   /* split on the separator */
   unsigned int nbTokens;
   char **tokens;
   fmstrsplit(basename,Fsep,&nbTokens,&tokens);

   if ( (nbTokens != 6) ) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted file name <%s>.\n",__func__,basename);
      return 1;
   }

   /* The last token contains the file extension (.nc) which should be removed */
   char *lastdot = strrchr(tokens[nbTokens-1],'.');
   *lastdot = '\0';

   //ice_drift_nh_polstere-625_ssmi-f15_201004031200-201004051200.nc

   /*  Assign the parameters to the tokens */
   *Area                   = tokens[2];
   *GridInfo               = tokens[3];
   *InstrumentAndPlatform  = tokens[4];

   /* Further break the Date information */
   char **tokens3;
   fmstrsplit(tokens[5],Ssep,&nbTokens,&tokens3);
   if (nbTokens != 2) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted Start+End date <%s>.\n",__func__,tokens[5]);
      return 1;
   }
   *startdate   = tokens3[0];
   *enddate     = tokens3[1];

   /* release the memory holding the array of array of chars (the strings are still allocated) */
   free(tokens);

   return 0;
}


void micedrift_build_filename(char *mergemethod, char *driftmethod, char *levfilter, char *area, char *startdate, 
      short starttwght, char *enddate, short endtwght, char **fname) {

   size_t nbchars = 0;
   nbchars += strlen(driftproduct_prefix) + 1;
   nbchars += strlen(mdriftproduct_id) + 1;
   nbchars += strlen(mergemethod) + 1;
   nbchars += strlen(levfilter) + 1;
   nbchars += strlen(driftmethod) + 1;
   nbchars += strlen(area) + 1;
   nbchars += strlen(startdate) + 1 + 1;
   nbchars += strlen(enddate) + 1;
   nbchars += 1 + strlen(ncExt);

   char cstarttwght, cendtwght;
   if (starttwght) cstarttwght = 'w';
   else cstarttwght = 'n';
   if (endtwght) cendtwght = 'w';
   else cendtwght = 'n';

   *fname = fmMalloc(nbchars+1);
   sprintf(*fname,"%s%s%s%s%s%s%s%s%s%s%s%s%s%c%s%s%c.%s",driftproduct_prefix,Fsep,mdriftproduct_id,Ssep,mergemethod,Fsep,
         driftmethod,Fsep,levfilter,Fsep,area,Fsep,
	 startdate,cstarttwght,Fsep,enddate,cendtwght,ncExt);

}

int mdriftproduct_split_filename(char *fname, char **Instrument, char **mergemethod, 
      char **driftmethod, char **levFilter, char **area, char **startdate, char **enddate) {
   
   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   /* split on the separator */
   unsigned int nbTokens;
   char **tokens;
   fmstrsplit(basename,Fsep,&nbTokens,&tokens);

   if ( (nbTokens != 7) || (strcmp(tokens[0],driftproduct_prefix)) ) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted file name <%s>.\n",__func__,basename);
      return 1;
   }

   /* The last token contains the file extension (.nc) which should be removed */
   char *lastdot = strrchr(tokens[nbTokens-1],'.');
   *lastdot = '\0';

   /* Have the Instrument field to reflect the 'multi' string */
   /* fmstrsplit */
   unsigned int nbStrings;
   char **Strings;
   fmstrsplit(tokens[1],Ssep,&nbStrings,&Strings);
   if (strcmp(Strings[0],mdriftproduct_id)) {
      fprintf(stderr,"ERROR (%s) Not a MULTI file <%s>.\n",__func__,basename);
      return 1;
   }
   if (nbStrings != 2) {
      fprintf(stderr,"ERROR (%s) Wrongly formatted file name <%s>. The 'multi' part misses the merging method.\n",__func__,basename);
      return 1;
   }
   *Instrument       = Strings[0];
   *mergemethod      = Strings[1];

   /*  Assign the parameters to the tokens */
   *driftmethod      = tokens[2];
   *levFilter        = tokens[3];
   *area             = tokens[4];
   *startdate        = tokens[5];
   *enddate          = tokens[6];

   /* the startdate and enddate are to be processed further has they have an extra character to
    * indicate if the time composite image was weighted ('w') or not ('n') */
   short startTweight = 0;
   if (strlen(*startdate) == 11) {
      char *weightChar = &((*startdate)[10]);
      if (*weightChar == 'w') {
	 startTweight = 1;
      }
      (*startdate)[10] = '\0';
   }

   short endTweight = 0;
   if (strlen(*enddate) == 11) {
      char *weightChar = &((*enddate)[10]);
      if (*weightChar == 'w') {
	 endTweight = 1;
      }
      (*enddate)[10] = '\0';
   }
   /* release the memory holding the array of array of chars (the strings are still allocated) */
   free(tokens);
   free(Strings);

   return 0;

}
int is_mdriftproduct_filename(char *fname, int *fmt_type) {
   
   /* fmbasename will remove (if present) the directory part in filename */
   char *endDirPart, *basename;
   fmbasename(fname,&endDirPart,&basename);

   return (is_driftproduct_filename(fname,fmt_type) &&
	 strstr(basename,mdriftproduct_id));

}

void dformproduct_build_filename_proc(char *Instrument, char *Platform, char *area, char *WaveBands, char *method, char *flevel, 
      char *startdate, short starttwght, char *enddate, short endtwght, char **fname) {

   char *descr = driftproduct_build_description(Instrument,Platform,area,WaveBands,method,flevel,
         startdate,starttwght,enddate,endtwght);
   if (!descr) {
      fprintf(stderr,"ERROR (%s) Cannot build the product description string.\n",__func__);
      *fname = NULL;
      return;
   } 

   *fname = fmMalloc(strlen(dformproduct_prefix)+strlen(Fsep)+strlen(descr)+1+strlen(ncExt)+1);
   sprintf(*fname,"%s%s%s.%s",dformproduct_prefix,Fsep,descr,ncExt);
   free(descr);

   return;

}

