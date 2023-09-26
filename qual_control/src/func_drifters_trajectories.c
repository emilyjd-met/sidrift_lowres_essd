
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <read_obsfile.h>
#include <errorcodes.h>
#include <projects.h>
#include <fmutil.h>
#include <dbl_list.h>
#include "icedrift_common.h"
#include "func_drifters_trajectories.h"


/* support functions for the buoyRecord datasets */
void printBuoyRecord(buoyRecord *b) {

   double lat = (b->data).lat;
   double lon = (b->data).lon;
   
   char latc='N';if (lat < 0) latc = 'S';
   char lonc='E';if (lon < 0) lonc = 'W';

   char timeStr[FMUTIL_ISODATETIME_LENGTH+1];
   fmsec19702isodatetime(b->time,timeStr);
   printf("(%c) <Pos: {%05.2f%c,%06.2f%c}; Time: %s (%ld)>",(b->valid?'G':'B'),lat,latc,lon,lonc,timeStr,b->time);

}
void setBuoyRecord(buoyRecord *b, Obsdata *data, fmsec1970 time) {

   /* store the observation data (copy) */
   memcpy(&(b->data),data,sizeof(Obsdata));
   b->time = time;
   b->valid = 1;

}
int latIsBelowLimit(buoyRecord *b) {
   int ret;
   extern float usedLatLimit;
   if (usedLatLimit > 0) {
      ret = (b->data.lat <= usedLatLimit);
   } else {
      ret = (b->data.lat >= usedLatLimit);
   }
   //printf("In %s: %f %f %f %d\n",__func__,b->data.lon, b->data.lat,usedLatLimit,ret);
   return ret;
}
int buoyRecord_isInvalid(buoyRecord *b) {
   return !(b->valid);
}

int buoyRecords_haveSameTimeStamp(buoyRecord *a,buoyRecord *b) {
   if (a->time == b->time) 
      return 1;
   return 0;
}

int BuoyRecord_byTimeStamp(buoyRecord *a,buoyRecord *b) {
   if (a->time == b->time) return 0;
   if (a->time < b->time) return -1;
   return 1;
}
double compute_distance_betweenTwoRecords(buoyRecord *firstRecord,buoyRecord *secondRecord,short only_with_valid) {
   Obsdata *P1 = &(firstRecord->data);
   Obsdata *P2 = &(secondRecord->data);

   double distance;
   if ( only_with_valid ) {
      if ( (!firstRecord->valid) || (!secondRecord->valid) ) 
	 return 0;
   }

   compute_distance(P1->lat,P1->lon,P2->lat,P2->lon,&distance);

   return distance;
}
double compute_timediff_betweenTwoRecords(buoyRecord *firstRecord,buoyRecord *secondRecord) {
   return (secondRecord->time - firstRecord->time);
}
     

/* support functions for the buoyTrajectory datasets */
void setBuoyTrajectory(buoyTrajectory *t,char *name,char *network,char *source,dbl_list *oneTraj) {

   //printf("in %s, name is <%s><%s><%s>\n",__func__,name,network,source);

   sprintf(t->id,"%s",name);
   sprintf(t->network,"%s",network);
   sprintf(t->source,"%s",source);
   t->records        = oneTraj;
   t->total_distance = 0;
   t->valid          = 1;

}
void printBuoyTrajectory(buoyTrajectory *t) {
   printf("<Buoy %s (%s/%s) Dist: %f km>\n",t->id,t->network,t->source,t->total_distance);
   dbl_list_print(t->records,&printBuoyRecord);
}
void printBuoyTrajectory_short(buoyTrajectory *t) {
   printf("<Buoy %s (%s/%s) nb_records: %d>\n",t->id,t->network,t->source,t->records->nbnodes);
}
int copyBuoyTrajectory(buoyTrajectory *to,buoyTrajectory *from) {

   strcpy(to->id,  from->id);
   strcpy(to->network,from->network);
   strcpy(to->source,from->source);
   to->records = dbl_list_copy(from->records,NULL);
   to->total_distance = from->total_distance;
   to->valid          = from->valid;
   if ( !(to->records) )
      return 1;
   return 0;

}

void deleteBuoyTrajectory(buoyTrajectory *t) {
   dbl_list_delete(t->records,NULL);
}

void sortRecords_byTimeStamp(buoyTrajectory *t) {
   int ret = dbl_list_sort(t->records,&BuoyRecord_byTimeStamp,NULL,NULL);
}

void uniqEachTrajectory(buoyTrajectory *t) {
  dbl_list_uniq(t->records,&buoyRecords_haveSameTimeStamp,NULL);
}

void computeTotalDriftDistance(buoyTrajectory *t) {

   double totalDriftDistance  = 0.;
   double totalDriftDistance2 = 0.;
   dbl_node *firstRecord = t->records->head;
   dbl_node *secondRecord;
   do {
      secondRecord = firstRecord->n;

      double dist = compute_distance_betweenTwoRecords(firstRecord->c,secondRecord->c,1);
      totalDriftDistance += dist;
      double dist2 = compute_distance_betweenTwoRecords(firstRecord->c,secondRecord->c,0);
      totalDriftDistance2 += dist;

      firstRecord = secondRecord;
      //printf("Distance : %f %f %f\n",dist,totalDriftDistance,totalDriftDistance2);
   } while (secondRecord->n != t->records->head);

   t->total_distance = totalDriftDistance;
}

int haveSameIDs(buoyTrajectory *t1,buoyTrajectory *t2) {
   return (!strcmp(t1->id,t2->id));
}

int trajectoryIsNotAlwaysAboveLimit(buoyTrajectory *t) {
   dbl_node *invalid_record = dbl_list_findIf(t->records,NULL,&latIsBelowLimit);
   return (invalid_record!=NULL);
}


int buoyHasNotWantedID(buoyTrajectory *t) {
   return (strcmp(t->id,buoyHasNotWantedID_params.id));
}

int trajectoryHasOnlyFewRecords(buoyTrajectory *t) {
   return (t->records->nbnodes <= 1);
}

int computeVelocities(dbl_list *records, double *vels, double *avg, double *std) {


   double avg_vel = 0;
   double std_vel = 0;

   /* compute the velocities*/
   size_t segN = 1;
   dbl_node *first = records->head;
   dbl_node *end   = first->p;
   dbl_node *sec;
   do {
      sec = first->n;

      double dist  = compute_distance_betweenTwoRecords(first->c,sec->c,0);
      dist        *= 1000; /* transform to meters */
      double tdiff = compute_timediff_betweenTwoRecords(first->c,sec->c);

      if (tdiff <= 0) {
	 fprintf(stderr,"ERROR (%s) the list is not ordered chronologically (or multiples exist). Tdiff = %f\n",__func__,tdiff);
	 return 1;
      } 
      double vel = dist/tdiff;
      vels[segN-1] = vel;
      //printf("\t Segment: %u, Velocity is %f\n",segN,vels[segN-1]);

      avg_vel += vel;
      std_vel += pow(vel,2.);

      first = sec;segN++;
   } while (sec != end);
  
   /* normalize the average value and compute the standard deviation */
   avg_vel /= segN;
   std_vel /= segN;
   std_vel  = - pow(avg_vel,2.) + std_vel; 
   std_vel  = sqrt(std_vel);

   *avg = avg_vel;
   *std = std_vel;

   return 0;
}

int analyseVelocities(dbl_list *records, double velLimit, short forward, short *valid, size_t *nbValid) {

   //printf("ANALYSE %s\n",(forward>0?"FORWARD":"BACKWARD"));


   /* initialize: first point is valid:*/
   if (forward>0)
      valid[0] = 1;
   else
      valid[records->nbnodes-1] = 1;
   *nbValid = 1;

   /* compute the velocities*/
   int    nbS      = (forward>0?1:records->nbnodes-2);
   int    incS     = (forward>0?+1:-1);
   dbl_node *first = (forward>0?records->head:records->head->p);
   dbl_node *end   = (forward>0?first->p:first->n);
   dbl_node *sec;
   do {
      sec = first;

      short LocValid = 0;
      while (!LocValid) { 
	 sec = (forward>0?sec->n:sec->p);
	 //printf("search for valid second.\n");
	 double dist  = compute_distance_betweenTwoRecords(first->c,sec->c,0);
	 dist        *= 1000; /* transform to meters */
	 double tdiff = compute_timediff_betweenTwoRecords(first->c,sec->c);
	 tdiff       *= (forward>0?+1:-1);
	 double vel = dist/tdiff;

	 if (vel > velLimit) {
	    LocValid = 0;
	 } else {
	    LocValid = 1;
	    (*nbValid)++;
	 }
	 valid[nbS] = LocValid;

	 /*
	 printf("P0 (%ld) -> P1 (%ld) V=%f. [%d] is %s (%d)\n",
	       ((buoyRecord *)(first->c))->time,
	       ((buoyRecord *)(sec->c))->time,
	       vel,nbS,(valid[nbS]?"GOOD":"BAD"),*nbValid);
	       */

	 if (sec == end)
	    break;

	 nbS+=incS;
      }


      //printf("Move first.\n");
      first = sec;
   } while (sec != end);

   return 0;

}

void filter_trajectory_onVelocities(buoyTrajectory *t) {
   int ret;

   if (t->records->nbnodes <= 1) {
      fprintf(stderr,"ERROR (%s) Trajectory is empty or has only 1 record.\n",__func__);
      t->valid = 0;
      return;
   }

   size_t nbSegments = t->records->nbnodes - 1 ;

   double *velocities = fmMalloc(nbSegments * sizeof(double));
   double avg,std;
//   printf("Compute velocities for %s\n",t->id);
   //dbl_list_print(t->records,&printBuoyRecord);
   if (computeVelocities(t->records,velocities,&avg,&std)) {
      fprintf(stderr,"ERROR (%s) could not compute velocities for trajectory.\n",__func__);
      t->valid = 0;
      return;
   }
   double Vlimit = avg+3*std;
  
   size_t nbvalid_forward;
   short *valid_forward  = fmMalloc(t->records->nbnodes * sizeof(short));
   size_t nbvalid_backward;
   short *valid_backward = fmMalloc(t->records->nbnodes * sizeof(short));

   ret  = 0;
   ret += analyseVelocities(t->records,Vlimit,+1,valid_forward, &nbvalid_forward);
   ret += analyseVelocities(t->records,Vlimit,-1,valid_backward,&nbvalid_backward);
   if (ret) {
      fprintf(stderr,"ERROR (%s) could not analyse and filter the velocities.\n",__func__);
      t->valid = 0;
      return;
   }

   size_t ind=0;
   dbl_node *start = t->records->head;
   do {
      buoyRecord *b = start->c;
      if (valid_forward[ind] && valid_backward[ind]) {
	 b->valid = 1; 
      } else {
	 b->valid = 0;
      }
      if (valid_forward[ind] && (nbvalid_forward > nbvalid_backward)) {
	 b->valid = 1;
      } 
      if (valid_backward[ind] && (nbvalid_backward > nbvalid_forward)) {
	 b->valid = 1;
      }
      ind++;
      start = start->n;
   } while (start != t->records->head);

   //dbl_list_print(t->records,&printBuoyRecord);

   free(velocities);
   free(valid_forward);free(valid_backward);

}

int trajectoryHasInvalidatedRecords (buoyTrajectory *t) {
   dbl_node *invalid_record = dbl_list_findIf(t->records,NULL,&buoyRecord_isInvalid);
   return (invalid_record!=NULL);
}

int trajectoryIsNotRepresentativeOfTimeSpan(buoyTrajectory *t) {
   dbl_list *records = t->records;
   buoyRecord *firstRecord = records->head->c;
   buoyRecord *lastRecord = records->head->p->c;

   fmsec1970 timeFirstRecord = firstRecord->time;
   fmsec1970 timeLastRecord  = lastRecord->time;

   extern fmsec1970 start_date_and_time, end_date_and_time;

   float timeRange = 1 * 60 * 60; /* 1 hour */

   if ( (abs(timeFirstRecord - start_date_and_time) <= timeRange) &&
	 (abs(timeLastRecord - end_date_and_time) <= timeRange) ) {
      return 0;
   } else {
      return 1;
   }

}
int buoyDidNotMove(buoyTrajectory *t) {
   return (t->total_distance == 0);
}

void keepOnlyFirstAndLastRecord (buoyTrajectory *t) {

   dbl_list *list  = t->records;
   dbl_node *start = list->head->n;
   dbl_node *end   = list->head->p;
   dbl_node *curr  = start;
   do {
      dbl_node *next = curr->n;
      dbl_list_removeNode(list,curr,NULL);
      curr = next;
   } while (curr != end);

}

/* support function for the buoyDrift datasets */
void transformDisplacement2Drift(buoyDrift *bd,buoyTrajectory *bt) {
   dbl_list *records = bt->records;
   buoyRecord *firstRecord = records->head->c;
   buoyRecord *lastRecord  = records->head->p->c;

   bd->start_lat = (firstRecord->data).lat;
   bd->start_lon = (firstRecord->data).lon;
   bd->start_time= firstRecord->time;
   bd->end_lat   = (lastRecord->data).lat;
   bd->end_lon   = (lastRecord->data).lon;
   bd->end_time  = lastRecord->time;
   sprintf(bd->id,"%s",bt->id);
   sprintf(bd->network,"%s",bt->network);
   bd->total_distance = bt->total_distance;
}

void setMappingParameters(buoyDrift *bd, PJ *proj, float *Ax, float *Ay, float *Bx, float *By,int *numX, int *numY) {
   bd->proj = proj;
   bd->Ax   = Ax;
   bd->Ay   = Ay;
   bd->Bx   = Bx;
   bd->By   = By;
   bd->numX = numX;
   bd->numY = numY;
}

void printBuoyDrift(buoyDrift *bd) {

   char start_timeStr[FMUTIL_ISODATETIME_LENGTH+1];
   fmsec19702isodatetime(bd->start_time,start_timeStr);
   char end_timeStr[FMUTIL_ISODATETIME_LENGTH+1];
   fmsec19702isodatetime(bd->end_time,end_timeStr);

   printf("<(%s,%f,%f) => (%s,%f,%f) [%d,%f,%f,%f]>",
	 start_timeStr,bd->start_lat,bd->start_lon,
	 end_timeStr,  bd->end_lat,  bd->end_lon,bd->inGrid,bd->xdrift,bd->ydrift,bd->total_distance);
}


void computeRemappedDrift(buoyDrift *bd) {

   double xstart,ystart,xend,yend;
   /* remap the start point */
   remap_ll2xy(bd->start_lat,bd->start_lon,bd->proj,*(bd->Ax),*(bd->Bx),*(bd->Ay),*(bd->By),&xstart,&ystart,0);
   /* remap the end point */
   remap_ll2xy(bd->end_lat,bd->end_lon,bd->proj,*(bd->Ax),*(bd->Bx),*(bd->Ay),*(bd->By),&xend,&yend,0);

   /* store the 'map-dependent' drift */
   bd->xdrift = *(bd->Ax) * (xend - xstart);
   bd->ydrift = *(bd->Ay) * (yend - ystart);

   printf("\t%s\tSTART: (%f,%f) => (%f %f) | END: (%f,%f) => (%f %f)\n",
	 bd->id,
	 bd->start_lat,bd->start_lon,xstart,ystart,
	 bd->end_lat,bd->end_lon,xend,yend);

   printf("\t%s\tAx:%f Bx:%f Ay:%f By:%f NX:%u NY:%u\n",bd->id,*(bd->Ax),*(bd->Bx),*(bd->Ay),*(bd->By),*(bd->numX),*(bd->numY));

   /* compute the interpolation corners and weights (only start position) */

   //ystart  = *(bd->numY) - 1 - ystart;
  
   int xlow  = floor(xstart); 
   int xhig  = ceil(xstart);
   int ylow  = floor(ystart); 
   int yhig  = ceil(ystart);
   int numX  = *(bd->numX);
   int numY  = *(bd->numY);
   if ( !((xlow < 0) || (xlow >= numX) ||
	  (xhig < 0) || (xhig >= numX) ||
	  (ylow < 0) || (ylow >= numY) ||
	  (yhig < 0) || (yhig >= numY)) )  {

      printf("Compute the interpolation parameters\n");

      bd->inGrid = 1;
      double xeps = xstart - xlow; 
      double yeps = ystart - ylow; 
      
      bd->interpCorners[0] = fmivec(xlow,ylow,numX);
      bd->interpCorners[1] = fmivec(xlow,yhig,numX);
      bd->interpCorners[2] = fmivec(xhig,ylow,numX);
      bd->interpCorners[3] = fmivec(xhig,yhig,numX);
      bd->interpWeights[0] = (1. - xeps) * (1. - yeps);
      bd->interpWeights[1] = (1. - xeps) * yeps;
      bd->interpWeights[2] = xeps        * (1. - yeps);
      bd->interpWeights[3] = xeps        * yeps;

      /*
	   startX[station] =  xstart;
	   startY[station] =  ystart;
	   driftX[station] =  driftx;
	   driftY[station] =  drifty;
	   */
   } else {
      printf("Point (%f %f) is not in the grid\n",xstart,ystart);
      bd->inGrid = 0;
   }

}
