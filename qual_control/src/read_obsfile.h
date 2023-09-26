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
 * read_obsfile.h
 *
 * PURPOSE:
 * Header file related to reading of ASCII observation data.
 *
 * REQUIREMENTS:
 * NA
 *
 * INPUT:
 * NA
 *
 * OUTPUT:
 * NA
 *
 * NOTES:
 * NA
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * Steinar Eastwood, met.no, 03.12.2007
 *
 * MODIFIED:
 * NA
 *
 * CVS_ID:
 * $Id: read_obsfile.h 1425 2007-12-03 10:51:37Z steinare $
 */ 

/* Values to indicate file type */
#define FTYPEUDEF "NA"
#define FTYPEDRAU "DRIBU"        /* met.no buoy observation file */
#define FTYPESYNO "SYNOP"        /* met.no synop observation file (ships ++) */
#define FTYPESHIP "FNMOCSHIP"    /* FNMOC ship observation file */
#define FTYPEPROF "FNMOCPROF"    /* FNMOC profile observation file */

#define STRLENSTID 16        /* Lenght of string stID in Obsdata */
#define STRLENSTTYPE 16      /* Lenght of string stType and stDesc in Obsdata */
#define STRLENSATF 32        /* Lenght of string satFile in MDBdata */


typedef struct {
  char stID[STRLENSTID];
  char stType[STRLENSTTYPE];
  char stDesc[STRLENSTTYPE];
  float lat;
  float lon;
  unsigned int year;
  unsigned int month;
  unsigned int day;
  unsigned int hour;
  unsigned int min;
  float Tw;
  float TTT;
  float Td;
  float PPPP;
  int ff;
  int dd;
  float RRR;
  int E;
  int sss;
  int N;
  int Nh;
  int Cl;
  int Cm;
  int Ch;
  int VV;
  int ww;
  int W1;
  int W2;
  int ci;
  int Si;
  int bi;
  int Di;
  int zi;
  int qualCode;
  int measLevel;
  int measDepth;
} Obsdata;


/* Function prototyping */

int read_obsdata(char *obsfile, Obsdata **obs, int *numObs, char *filetype);

int read_obsdata_fnmoc(char *obsfile, Obsdata **obs, int *numObs);

int read_obsdata_old(char *obsfile, Obsdata **obs, int *numObs);

int read_QCobsdata(char *obsfile1, char *obsfile2, Obsdata **obs, int *numObs);
