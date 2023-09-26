
#include <stdio.h>
#include <string.h>
#include <icedrift_instruments.h>

/* Global variables for each instrument */
/* http://en.wikipedia.org/wiki/Special_sensor_microwave/imager */
/* SSM/I */
char  *PossibleWaveBandsSSMI[]   = {"bt19v","bt19h","bt22v","bt37v","bt37h","bt85v","bt85h"};
double PossibleFOVMAppertSSMI[]  = {    70.,    70.,    60.,    38.,    38.,    16.,    16.}; /* 38 -> 65 */
double PossibleRadNoiseSSMI[]    = {   0.45,   0.42,   0.74,   0.37,   0.38,   0.69,   0.73};
char  *DefaultWaveBandsSSMI[]    = {"bt85v","bt85h"};
/* SSMIS */
char  *PossibleWaveBandsSSMIS[]  = {"bt19v","bt19h","bt37v","bt37h","bt91v","bt91h"};
double PossibleFOVMAppertSSMIS[] = {    70.,    70.,    38.,    38.,    16.,    16.};
double PossibleRadNoiseSSMIS[]   = {     1.,     1.,     1.,     1.,     1.,     1.};
char  *DefaultWaveBandsSSMIS[]   = {"bt91v","bt91h"};
/* AMSR */
char  *PossibleWaveBandsAMSR[]   = {"bt19v","bt19h","bt37v","bt37h","bt89v","bt89h"};
double PossibleFOVMAppertAMSR[]  = {    27.,    27.,    14.,    14.,     7.,     7.};
double PossibleRadNoiseAMSR[]    = {     1.,     1.,     1.,     1.,     1.,     1.};
char  *DefaultWaveBandsAMSR[]    = {"bt37v","bt37h"};
/* ASCAT */
char  *PossibleWaveBandsASCAT[]  = {"sigma0"};
double PossibleFOVMAppertASCAT[] = {     25.};
double PossibleRadNoiseASCAT[]   = {     1.};
char  *DefaultWaveBandsASCAT[]   = {"sigma0"};
/* SEAWINDS */
char  *PossibleWaveBandsSEAWINDS[]  = {"sigma0H","sigma0V"};
double PossibleFOVMAppertSEAWINDS[] = {     25.};
double PossibleRadNoiseSEAWINDS[]   = {     1.};
char  *DefaultWaveBandsSEAWINDS[]   = {"sigma0H","sigma0V"};
/* MIRAS/SMOS */
char  *PossibleWaveBandsMIRAS[]  = {"btL"};
double PossibleFOVMAppertMIRAS[] = {  15.};
double PossibleRadNoiseMIRAS[]   = {   1.};
char  *DefaultWaveBandsMIRAS[]   = {"btL"};
/* CIMR */
char  *PossibleWaveBandsCIMR[]   = {"bt19v","bt19h","bt37v","bt37h"};
double PossibleFOVMAppertCIMR[]  = { 14.,    14.,     7.,     7.};
double PossibleRadNoiseCIMR[]    = {  1.,     1.,     1.,     1.};
char  *DefaultWaveBandsCIMR[]    = {"bt37v","bt37h"};
/* MULTI */
char  *PossibleWaveBandsMULTI[]  = {NULL};
double PossibleFOVMAppertMULTI[] = {0.};
double PossibleRadNoiseMULTI[]   = {0.};
char  *DefaultWaveBandsMULTI[]   = {NULL};

/* all in an array */
char **PossibleBands[] = {PossibleWaveBandsSSMI,PossibleWaveBandsSSMIS,PossibleWaveBandsAMSR,PossibleWaveBandsASCAT,PossibleWaveBandsSEAWINDS,PossibleWaveBandsMIRAS,PossibleWaveBandsAMSR,PossibleWaveBandsCIMR,PossibleWaveBandsMULTI};
double *PossibleMApp[] = {PossibleFOVMAppertSSMI,PossibleFOVMAppertSSMIS,PossibleFOVMAppertAMSR,PossibleFOVMAppertASCAT,PossibleFOVMAppertSEAWINDS,PossibleFOVMAppertMIRAS,PossibleFOVMAppertAMSR,PossibleFOVMAppertCIMR,PossibleFOVMAppertMULTI};
double *PossibleRNoise[] = {PossibleRadNoiseSSMI,PossibleRadNoiseSSMIS,PossibleRadNoiseAMSR,PossibleRadNoiseASCAT,PossibleRadNoiseSEAWINDS,PossibleRadNoiseMIRAS,PossibleRadNoiseAMSR,PossibleRadNoiseCIMR,PossibleRadNoiseMULTI};
char **DefaultBands[] =  {DefaultWaveBandsSSMI,DefaultWaveBandsSSMIS,DefaultWaveBandsAMSR,DefaultWaveBandsASCAT,DefaultWaveBandsSEAWINDS,DefaultWaveBandsMIRAS,DefaultWaveBandsAMSR,DefaultWaveBandsCIMR,DefaultWaveBandsMULTI};

short  nbPossibleBands[] = {sizeof(PossibleWaveBandsSSMI)/sizeof(*PossibleWaveBandsSSMI),\
                            sizeof(PossibleWaveBandsSSMIS)/sizeof(*PossibleWaveBandsSSMIS),\
                            sizeof(PossibleWaveBandsAMSR)/sizeof(*PossibleWaveBandsAMSR),\
                            sizeof(PossibleWaveBandsASCAT)/sizeof(*PossibleWaveBandsASCAT),\
                            sizeof(PossibleWaveBandsSEAWINDS)/sizeof(*PossibleWaveBandsSEAWINDS),\
                            sizeof(PossibleWaveBandsMIRAS)/sizeof(*PossibleWaveBandsMIRAS),\
                            sizeof(PossibleWaveBandsAMSR)/sizeof(*PossibleWaveBandsAMSR),\
                            sizeof(PossibleWaveBandsCIMR)/sizeof(*PossibleWaveBandsCIMR),\
                            0};

short  nbDefaultBands[] = {sizeof(DefaultWaveBandsSSMI)/sizeof(*DefaultWaveBandsSSMI),\
                            sizeof(DefaultWaveBandsSSMIS)/sizeof(*DefaultWaveBandsSSMIS),\
                            sizeof(DefaultWaveBandsAMSR)/sizeof(*DefaultWaveBandsAMSR),\
                            sizeof(DefaultWaveBandsASCAT)/sizeof(*DefaultWaveBandsASCAT),\
                            sizeof(DefaultWaveBandsSEAWINDS)/sizeof(*DefaultWaveBandsSEAWINDS),\
                            sizeof(DefaultWaveBandsMIRAS)/sizeof(*DefaultWaveBandsMIRAS),\
                            sizeof(DefaultWaveBandsAMSR)/sizeof(*DefaultWaveBandsAMSR),\
                            sizeof(DefaultWaveBandsCIMR)/sizeof(*DefaultWaveBandsCIMR),\
                            0};

/* routines to test which instrument a file corresponds to */
int is_SSMIS(char *name) {
   return !strcmp(name,"ssmis");
}
int is_SSMI(char *name) {
   return !strcmp(name,"ssmi");
}
int is_AMSR(char *name) {
   return !strcmp(name,"amsr");
}
int is_AMSR2(char *name) {
   return !strcmp(name,"amsr2");
}
int is_ASCAT(char *name) {
   return !strcmp(name,"ascat");
}
int is_SEAWINDS(char *name) {
   return !strcmp(name,"seawinds");
}
int is_MIRAS(char *name) {
   return !strcmp(name,"miras");
}
int is_CIMR(char *name) {
   return !strcmp(name,"cimr");
}
int is_MULTI(char *name) {
   return !strcmp(name,"multi");
}

int instrumentType(char *name) {

   int ret = INSTRUMENT_NDEF;
   
   if (is_SSMI(name))
      ret = INSTRUMENT_SSMI;
   else if (is_SSMIS(name))
      ret = INSTRUMENT_SSMIS;
   else if (is_AMSR(name))
      ret = INSTRUMENT_AMSR;
   else if (is_AMSR2(name))
      ret = INSTRUMENT_AMSR2;
   else if (is_ASCAT(name))
      ret = INSTRUMENT_ASCAT;
   else if (is_SEAWINDS(name))
      ret = INSTRUMENT_SEAWINDS;
   else if (is_MIRAS(name))
      ret = INSTRUMENT_MIRAS;
   else if (is_CIMR(name))
      ret = INSTRUMENT_CIMR;
   else if (is_MULTI(name))
      ret = INSTRUMENT_MULTI;

   return ret;

}

int get_default_WaveBands(char *Instrument, short *nbdefBands, char ***defBands) {
   int inst = instrumentType(Instrument);
   if ( inst == INSTRUMENT_NDEF) {
      printf("ERROR (%s) Unknown instrument (%s)\n",__func__,Instrument);
      return 0;
   }
   *defBands       = DefaultBands[inst];
   *nbdefBands     = nbDefaultBands[inst];
   return 1;
}

int compatibleInstrumentAndPlatform(char *Instrument, char *Platform) {
   
   int inst = instrumentType(Instrument);
   if ( inst == INSTRUMENT_NDEF) {
      printf("ERROR (%s) Unknown instrument (%s)\n",__func__,Instrument);
      return 0;
   }
   switch (inst) {
      case INSTRUMENT_SSMI:
         {
            int dmspN,ok;
            ok = sscanf(Platform,"f%d",&dmspN);
            if (ok != 1) {
               fprintf(stderr,"ERROR (%s) Not a DMSP platform (%s) for an SSM/I instrument\n",__func__,Platform);
               return 0;
            }
            if (dmspN >= 16) {
               fprintf(stderr,"ERROR (%s) Platform is f%02d (max for SSM/I is f15).\n",__func__,dmspN);
               return 0;
            }
         }
         break;
      case INSTRUMENT_SSMIS:
         {
            int dmspN,ok;
            ok = sscanf(Platform,"f%d",&dmspN);
            if (ok != 1) {
               fprintf(stderr,"ERROR (%s) Not a DMSP platform (%s) for an SSMIS instrument\n",__func__,Platform);
               return 0;
            }
            if (dmspN < 16) {
               fprintf(stderr,"ERROR (%s) Platform is f%02d (min for SSMIS is f16).\n",__func__,dmspN);
               return 0;
            }
         }
         break;
      case INSTRUMENT_AMSR:
         if (strcmp(Platform,"aqua")) {
            fprintf(stderr,"ERROR (%s) Not Aqua platform (%s) for AMSR instrument\n",__func__,Platform);
            return 0;
         }
         break;
      case INSTRUMENT_AMSR2:
         if (!strstr(Platform,"gw")) {
            fprintf(stderr,"ERROR (%s) Not gw? platform (%s) for AMSR2 instrument\n",__func__,Platform);
            return 0;
         }
         break;
      case INSTRUMENT_ASCAT:
         if (!strstr(Platform,"metop")) {
            fprintf(stderr,"ERROR (%s) Not metop? platform (%s) for ASCAT instrument\n",__func__,Platform);
            return 0;
         }
         break;
      case INSTRUMENT_SEAWINDS:
         if (strcmp(Platform,"qscat")) {
            fprintf(stderr,"ERROR (%s) Not qscat platform (%s) for SEAWINDS instrument\n",__func__,Platform);
            return 0;
         }
         break;
      case INSTRUMENT_MIRAS:
         if (strcmp(Platform,"smos")) {
            fprintf(stderr,"ERROR (%s) Not smos platform (%s) for MIRAS instrument\n",__func__,Platform);
            return 0;
         }
         break;
      case INSTRUMENT_CIMR:
         if (strcmp(Platform,"s9a")) {
            fprintf(stderr,"ERROR (%s) Not S9a platform (%s) for CIMR instrument\n",__func__,Platform);
            return 0;
         }
         break;
      case INSTRUMENT_MULTI:
         if (strcmp(Platform,"sel") && strcmp(Platform,"oi")) {
            fprintf(stderr,"ERROR (%s) Not 'oi' or 'sel' method (%s) for MULTI instrument\n",__func__,Platform);
            return 0;
         }
   }

   return 1;
}

int compatibleInstrumentAndPlatformAndWaveBand(char *Instrument, char *Platform, char **WaveBands, short nbWaveBands) {

   if (!compatibleInstrumentAndPlatform(Instrument,Platform)) {
      printf("ERROR (%s) Instrument (%s) and Platform (%s) are not compatible.\n",__func__,Instrument,Platform);
      return 0;
   }
   
   int instType = instrumentType(Instrument);
   if ( instType == INSTRUMENT_NDEF) {
      printf("ERROR (%s) Unknown instrument (%s)\n",__func__,Instrument);
      return 0;
   }

   short found_one_issue = 0;
   for ( short b = 0 ; b < nbWaveBands ; b++ ) {
      for ( size_t pb = 0 ; pb < nbPossibleBands[instType] ; pb++ ) {
         if (strcmp(WaveBands[b],PossibleBands[instType][pb])) {
            found_one_issue = 1;
            printf("ERROR (%s) Channel (%s) is not compatible with %s instrument.\n",__func__,WaveBands[b],Instrument);
            break;
         }
      }
   }
   if (found_one_issue) {
      return 0;
   }

   return 1;
}

/* 
 * DEPRECATED
int instrumentLongName(char *name,char *longName,short *format) {

   int ret = INSTRUMENT_NDEF;

   if ( (*format < TAGFORMAT_CODE) || (*format >= TAGFORMAT_NUMBER) ) {
      *format = TAGFORMAT_CODE;
   }
   
   if (is_SSMI(name)) {
      switch (*format) {
     case TAGFORMAT_CODE:
        sprintf(longName,"ssmi");
        break;
     case TAGFORMAT_LONG:
        sprintf(longName,"SSM/I");
        break;
         case TAGFORMAT_GCMD: 
        sprintf(longName,"SSM/I > Special Sensor Microwave/Imager");
        break;
      }
      ret = INSTRUMENT_SSMI;
   } else if (is_AMSR(name)) {
      switch (*format) {
     case TAGFORMAT_CODE:
        sprintf(longName,"amsr");
        break;
     case TAGFORMAT_LONG:
        sprintf(longName,"AMSR-E");
        break;
         case TAGFORMAT_GCMD: 
        sprintf(longName,"AMSR-E > Advanced Microwave Scanning Radiometer-EOS");
        break;
      }
      ret = INSTRUMENT_AMSR;
   } else if (is_ASCAT(name)) {
      switch (*format) {
     case TAGFORMAT_CODE:
        sprintf(longName,"ascat");
        break;
     case TAGFORMAT_LONG:
        sprintf(longName,"ASCAT");
        break;
         case TAGFORMAT_GCMD: 
        sprintf(longName,"ASCAT > Advanced SCATterometer");
        break;
      }
      ret = INSTRUMENT_ASCAT;
   }

   return ret;
}

int instrumentPlatform(char *name, char *platform, short *format) {

   int ret = INSTRUMENT_NDEF;
   if ( (*format < TAGFORMAT_CODE) || (*format >= TAGFORMAT_NUMBER) ) {
      *format = TAGFORMAT_CODE;
   }

   if (is_SSMI(name)) {
      int ssmi_platform;
      int nb = sscanf(name,"ssmi-f%d",&ssmi_platform);
      if (nb != 1) {
     return (ret);
      }
      switch(*format) {
     case TAGFORMAT_CODE:
        sprintf(platform,"f%d",ssmi_platform);
        break;
     case TAGFORMAT_LONG:
        sprintf(platform,"DMSP-F%02d",ssmi_platform);
        break;
         case TAGFORMAT_GCMD: 
        {
        short gcmd_dmsp_type = 2;
        if (ssmi_platform <= 4) 
           gcmd_dmsp_type = 1;
        sprintf(platform,"DMSP 5D-%d/F%d >  Defense Meteorological Satellite Program-F%d",gcmd_dmsp_type,ssmi_platform,ssmi_platform);
        break;
        }
      }
      ret = INSTRUMENT_SSMI;
   } else if (is_AMSR(name)) {
      switch(*format) {
     case TAGFORMAT_CODE:
        sprintf(platform,"aqua");
        break;
     case TAGFORMAT_LONG:
        sprintf(platform,"EOS-Aqua");
        break;
         case TAGFORMAT_GCMD: 
        sprintf(platform,"AQUA > Earth Observing System, AQUA");
        break;
      }
      ret = INSTRUMENT_AMSR;
   } else if (is_ASCAT(name)) {
      switch(*format) {
     case TAGFORMAT_CODE:
        sprintf(platform,"metopA");
        break;
     case TAGFORMAT_LONG:
        sprintf(platform,"METOP-A");
        break;
         case TAGFORMAT_GCMD: 
        sprintf(platform,"METOP-A > Meteorological Operational Satellite - A");
        break;
      }
      ret = INSTRUMENT_ASCAT;
   }

   return ret;
}
*/

int get_channel_apperture_m(char *Instrument, char *Platform, char **WaveBands, short nbWaveBands, 
      double *apps) {

   if (!compatibleInstrumentAndPlatformAndWaveBand(Instrument,Platform,WaveBands,nbWaveBands)) {
      fprintf(stderr,"ERRPR (%s) Something is wrong with this Instrument/Platform/WaveBands combinations",
            __func__);
      return 0;
   }
   int instType = instrumentType(Instrument);
   if ( instType == INSTRUMENT_MULTI ) {
      printf("ERROR (%s) Multi sensor %s does not have an apperture\n",__func__,Instrument);
      return 0;
   }

   for ( short b = 0 ; b < nbWaveBands ; b++ ) {
      apps[b] = PossibleMApp[instType][b] * 1000.; /* km -> m */
   }
   return 1;

}
