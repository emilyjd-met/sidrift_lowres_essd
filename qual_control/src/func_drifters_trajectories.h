#ifndef FUNC_DRIFTERS_TRAJECTORIES_H
#define FUNC_DRIFTERS_TRAJECTORIES_H

/* maximum number of characters that will
 * be written in the netcdf file for the name
 * of the station. */
#define STATION_NAME_LENGTH 10
#define STATION_NETW_LENGTH 10
#define STATION_SRC_LENGTH  10


typedef struct {
   /* drift */
   char      id[STATION_NAME_LENGTH+1];
   char      network[STATION_NETW_LENGTH+1];
   float     start_lat; float start_lon; fmsec1970 start_time;
   float     end_lat;   float end_lon;   fmsec1970 end_time;
   float     total_distance;
   /* remapped drift */
   float     xdrift;    float ydrift;
   /* interpolation information */
   short     inGrid;
   size_t    interpCorners[4];
   float     interpWeights[4];
   /* remapping information  */
   PJ       *proj;
   float    *Ax,*Ay,*Bx,*By;
   int      *numX,*numY;
} buoyDrift;

typedef struct {
   char      id[STATION_NAME_LENGTH+1];
   char      network[STATION_NETW_LENGTH+1];
   char      source[STATION_SRC_LENGTH+1];
   dbl_list *records;
   double    total_distance ;
   short     valid;
} buoyTrajectory;

typedef struct {
   Obsdata   data;
   fmsec1970 time;
   short     valid;
} buoyRecord;


void computeRemappedDrift(buoyDrift *bd);

/* support functions for the buoyRecord datasets */
void printBuoyRecord(buoyRecord *b);
void setBuoyRecord(buoyRecord *b, Obsdata *obs, fmsec1970 time);
int BuoyRecord_byTimeStamp(buoyRecord *a,buoyRecord *b);
int buoyRecords_haveSameTimeStamp(buoyRecord *a,buoyRecord *b);

/* support functions for the buoyTrajectory datasets */
void setBuoyTrajectory(buoyTrajectory *t,char *name,char *network,char *source,dbl_list *oneTrajectory);
void printBuoyTrajectory(buoyTrajectory *t);
void printBuoyTrajectory_short(buoyTrajectory *t);
int  copyBuoyTrajectory(buoyTrajectory *to,buoyTrajectory *from);
void deleteBuoyTrajectory(buoyTrajectory *t);
int  haveSameIDs(buoyTrajectory *t1,buoyTrajectory *t2);
int  trajectoryIsNotAlwaysAboveLimit(buoyTrajectory *t);
int  trajectoryHasOnlyFewRecords(buoyTrajectory *t);
void computeTotalDriftDistance(buoyTrajectory *t);
void sortRecords_byTimeStamp(buoyTrajectory *t);
void uniqEachTrajectory(buoyTrajectory *t);
void filter_trajectory_onVelocities(buoyTrajectory *t);
int trajectoryHasInvalidatedRecords (buoyTrajectory *t);
int trajectoryIsNotRepresentativeOfTimeSpan(buoyTrajectory *t);
int buoyDidNotMove(buoyTrajectory *t);
void keepOnlyFirstAndLastRecord(buoyTrajectory *t);
struct {
   char id[32];
} buoyHasNotWantedID_params;
int buoyHasNotWantedID(buoyTrajectory *t);


/* support function for the buoyDrift datasets */
void transformDisplacement2Drift(buoyDrift *bd,buoyTrajectory *bt);
void setMappingParameters(buoyDrift *bd, PJ *proj, float *Ax, float *Ay, float *Bx, float *By, int *numX,int *numY);
void printBuoyDrift(buoyDrift *bd);

/* global variables used for communication between caller and routines in this file */
float usedLatLimit;
fmsec1970 start_date_and_time;
fmsec1970 end_date_and_time;

#endif /* FUNC_DRIFTERS_TRAJECTORIES_H */
