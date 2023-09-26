
#ifndef ICEDRIFT_FILENAMES
#define ICEDRIFT_FILENAMES

#define LRSID_FILEFORMAT_UNKN -1
#define LRSID_FILEFORMAT_PROC  0
#define LRSID_FILEFORMAT_FINL  1

extern char Fsep[];
extern char Ssep[];

extern char ncExt[];
extern char tcimage_prefix[];

int isValid_InstrumentAndPlatform(char InstrumentAndPlatform[],char **Instrument, char **Platform);
int isValid_InstrumentAndPlatformAndChannel(char InstrumentAndPlatform[],char **Instrument, char **Platform, char **Channels);

void tcimage_build_filename(short Tweight,char *Instrument, char *Platform, char *WaveBands, char *area, char *centralTime,char **fname);
int  tcimage_split_filename(char *fname, short *Tweight, char **Instrument, char **Platform, char **WaveBands, char **area, char **centralTime);
int  is_tcimage_filename(char *fname);

char *driftproduct_build_description(char *Instrument, char *Platform, char *area, char *WaveBands, char *method, char *flevel, 
      char *startdate, short startwght, char *enddate, short endwght);
void driftproduct_build_filename_proc(char *Instrument, char *Platform, char *area, char *WaveBands, char *method, char *flevel, 
      char *startdate, short starttwght, char *enddate, short endtwght, short is_temporal_average, char **fname);
void driftproduct_build_filename_final(char *Prefix, char *Area, char *GridInfo, char *InstrumentAndPlatform, 
      char *startdate, char *enddate, char FFsep, char FSsep, char **fname);
int is_driftproduct_filename(char *fname,int *fmt_type);
int driftproduct_split_filename_proc(char *fname, char **Instrument, char **Platform, char **WaveBands,char **method,
      char **levFilter,char **area, char **startdate, char **enddate);
int driftproduct_split_filename_final(char *fname, char **Area, char **GridInfo, char **InstrumentAndPlatform, 
      char **startdate, char **enddate);
int mdriftproduct_split_filename(char *fname, char **Instrument, char **mergemethod, 
      char **driftmethod, char **levFilter, char **area, char **startdate, char **enddate);
void micedrift_build_filename(char *mergemethod, char *driftmethod, char *levfilter, char *area, char *startdate, 
      short starttwght, char *enddate, short endtwght, char **fname);
int is_mdriftproduct_filename(char *fname,int *fmt_type);

void dformproduct_build_filename_proc(char *Instrument, char *Platform, char *area, char *WaveBands, char *method, char *flevel, 
      char *startdate, short starttwght, char *enddate, short endtwght, char **fname);

#endif /* ICEDRIFT_FILENAMES */
