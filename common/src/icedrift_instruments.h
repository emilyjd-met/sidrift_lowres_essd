
#ifndef ICEDRIFT_INSTRUMENTS_H
#define ICEDRIFT_INSTRUMENTS_H

#define INSTRUMENT_NDEF  -1
#define INSTRUMENT_SSMI   0
#define INSTRUMENT_SSMIS  1
#define INSTRUMENT_AMSR   2
#define INSTRUMENT_ASCAT  3
#define INSTRUMENT_SEAWINDS  4
#define INSTRUMENT_MIRAS  5
#define INSTRUMENT_AMSR2  6
#define INSTRUMENT_CIMR   7
#define INSTRUMENT_MULTI  8
#define INSTRUMENT_NUMBER 9

#define TAGFORMAT_NDEF   -1
#define TAGFORMAT_CODE    0
#define TAGFORMAT_LONG    1
#define TAGFORMAT_GCMD    2
#define TAGFORMAT_NUMBER  3


extern char *PossibleWaveBandsSSMI[];
extern double PossibleFOVMAppertSSMI[];
extern char *PossibleWaveBandsSSMIS[];
extern double PossibleFOVMAppertSSMIS[];
extern char *PossibleWaveBandsAMSR[];
extern double PossibleFOVMAppertAMSR[];
extern char *PossibleWaveBandsASCAT[];
extern double PossibleFOVMAppertASCAT[];
extern char *PossibleWaveBandsQSCAT[];
extern double PossibleFOVMAppertQSCAT[];
extern char *PossibleWaveBandsMIRAS[];
extern double PossibleFOVMAppertMIRAS[];
extern double PossibleFOVMAppertCIMR[];
extern char *PossibleWaveBandsMULTI[];
extern double PossibleFOVMAppertMULTI[];
extern char **PossibleBands[];
extern double *PossibleMApp[];
extern double *PossibleRNoise[];
extern short nbPossibleBands[];

extern char *DefaultWaveBandsSSMI[];
extern char *DefaultWaveBandsSSMIS[];
extern char *DefaultWaveBandsAMSR[];
extern char *DefaultWaveBandsASCAT[];
extern char *DefaultWaveBandsQSCAT[];
extern char *DefaultWaveBandsMIRAS[];
extern char *DefaultWaveBandsCIMR[];
extern char *DefaultWaveBandsMULTI[];
extern char **DefaultBands[];
extern short nbDefaultBands[];

/* routines to test which instrument a file corresponds to */
int is_SSMI(char *name);
int is_AMSR(char *name);
int is_ASCAT(char *name);
int is_SEAWINDS(char *name);
int is_MIRAS(char *name);
int is_CIMR(char *name);
int is_MULTI(char *name);
int instrumentType(char *name);

int compatibleInstrumentAndPlatform(char *Instrument, char *Platform);
int compatibleInstrumentAndPlatformAndWaveBand(char *Instrument, char *Platform, char **WaveBands, short nbWaveBands);
int get_default_WaveBands(char *Intrument, short *nbdefBands, char ***defBands);

/* routines to format the instrument/platform string */
/* DEPRECATED 
int instrumentLongName(char *name, char *longName, short *fmt);
int instrumentPlatform(char *name, char *platform, short *fmt);
*/

/* routines to obtain FOV apperture information */
int get_channel_apperture_m(char *Instrument, char *Platform, char **WaveBands, short nbWaveBands, double *apps);

#endif /* ICEDRIFT_INSTRUMENTS_H */
