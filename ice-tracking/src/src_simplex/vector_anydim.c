

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "memory.h"
#include "vector_anydim.h"

Vect *NewVect(size_t nx,double x[],double score)
{

   if (!nx) return NULL;

   Vect *new;
   new     = xMalloc(sizeof(Vect));
   new->nx = nx;
   size_t size = nx*sizeof(double);
   new->x  = xMalloc(size);
   if (x)
      memcpy(new->x,x,size);
   new->score = score;

   return new;
}


Vect *CopyVect(Vect *ori)
{
   Vect *copy;
   copy = NewVect(ori->nx,ori->x,ori->score);
   
   return copy;
}

void AssignVect(Vect *v1,const Vect *v2)
{
   v1->nx = v2->nx;
   memcpy(v1->x,v2->x,sizeof(double)*v2->nx);
   v1->score = v2->score;
}

void ExtractX(const Vect *v,double *x) {
   memcpy(x,v->x,sizeof(double)*(v->nx));
}


void FreeVect(Vect *v)
{
   free(v->x);
   free(v);
}

void PrintVect(const Vect *v)
{
   /* first part (the params) should be delegated to the physical X */
   /* param values: */
   printf("X: {");
   for (size_t n = 0 ; n < v->nx ; n ++) {
      printf("%f,",v->x[n]);
   }
   /* finish with score */
   printf("\b}. Score=%.10f\n",v->score);
}

void PrintCoordVect(const Vect *v)
{
   /* we should also delegate to the phisical X */
   printf("X: {");
   for (size_t n = 0 ; n < v->nx ; n ++) {
      printf("%f,",v->x[n]);
   }
   printf("\b}.\n");

}

int ValidVect(const Vect *v) 
{
   /* should delegate to the physical X */
   return 1;
}

int EvaluateVect(Vect *v, int (*func)(size_t n,double *x,double *f)) {
   /* we should delegate to the phisical X or directely to the model */

   if (!ValidVect((const Vect *)v)) {
      printf("+++++++++++++++\n");
      printf("***** BIG ERROR : evaluate a non valid point : ");PrintVect((const Vect *)v);
      printf("+++++++++++++++\n");
   }
   
   int  everything_ok;
 
   everything_ok = (*func)(v->nx,v->x,&(v->score));
   
   return everything_ok;
}

/* compute the cartesian centroid of nvecs vectors, all with mid->nx dimensions */
void GetCentroid(size_t nvecs,const Vect **vs,Vect *mid) {

   /* loop over the dimensions */
   for (size_t n = 0 ; n < mid->nx ; n ++) {
      mid->x[n] = 0.;
      /* loop over the vectors */
      for (size_t nv = 0 ; nv < nvecs ; nv ++) {
	 mid->x[n] += vs[nv]->x[n];
      }
      mid->x[n] /= mid->nx;
   }

}

/* compute the difference of two vectors */
void SubVect(const Vect *v1,const Vect *v2,Vect *sub) {

   for (size_t n = 0 ; n < sub->nx ; n++) {
      sub->x[n] = v1->x[n] - v2->x[n];
   }
}

/* shift an origin point via a scaled direction. The output overwrites the 'direction' */
void AddScaleVect(const Vect *origin, Vect *direction, const double lambda)
{
   for (size_t n = 0 ; n < origin->nx ; n++) {
      direction->x[n] = origin->x[n] + lambda*direction->x[n];
   }
}

int better_if_lower (double s1,double s2) { return(s1 < s2); }
int better_if_higher(double s1,double s2) { return(s1 > s2); }

/* rank the nv vectors by descending order of score value,
 * so that v[0].score < v[1].score < ... < v[n-1].score */
void RankByScores(size_t nv, Vect **vx, Vect *swap,int (*ranker)(double s1,double s2)) {

   /* tbs: element to be processed */
   for (size_t n2 = 0 ; n2 < (nv-1) ; n2++) {
   
      /* search the greatest score in the rest of the array */
      Vect *sv = vx[n2];
      for (size_t n1 = n2+1 ; n1 < nv ; n1++) {
         //if (vx[n1]->score < sv->score)
         if ((*ranker)(vx[n1]->score,sv->score))
	    sv = vx[n1];
      }
      /* swap the best elem in sub array with the first element */
      AssignVect(swap,(const Vect *)vx[n2]);
      AssignVect(vx[n2],(const Vect *)sv);
      AssignVect(sv,(const Vect *)swap);

   }

}
