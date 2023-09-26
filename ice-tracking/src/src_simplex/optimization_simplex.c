
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "memory.h"
#include "vector_anydim.h"
#include "optimization_simplex.h"

#define ITERA_MAX 1000u
#define FEPS   1e-6

#define INSIDE_CONTRACTION  1
#define OUTSIDE_CONTRACTION 2
#define SHRINK_CONTRACTION  3


static int check_tol(const double,const double);

int findBestX(size_t ndim, double *x, double *score, double *startingPoints[],
      char opti, int (*func)(size_t nx,double *x,double *f))
{

   int success=ERROR;
   /* prepare for a minimization or a maximization */
   int (*isBetter)(double,double);
   double initScore;
   switch(opti) {
      case '<': /* minimization */
#ifdef MINIM_DEBUG
   printf("\nSetting up for a MINIMIZATION problem.\n\n");
#endif
	 initScore = VERY_HIGH_VALUE;
         isBetter = &better_if_lower;
	 break;
      case '>': /* maximization */
#ifdef MINIM_DEBUG
   printf("\nSetting up for a MAXIMIZATION problem.\n\n");
#endif	 
         initScore = VERY_LOW_VALUE;
         isBetter = &better_if_higher;
         break;
      default: /* unknown */
	 fprintf(stderr,"ERROR (findBestX): optiChar is '%c': minimization ('<'), maximization ('>')\n",opti);
	 return success;
   }

   short contraction;
   int ok;


   size_t nv = ndim+1;  /* number of vertices: if x has 2 dimensions, we need 3 vertices */
   Vect **vertices;

   /* indexes to find the worst, second worst and best element in the ranked array of nv vertices */
   size_t worst  = nv-1;
   size_t sworst = nv-2;
   size_t best   = 0;

   Vect *swap        = NewVect(ndim,NULL,initScore);
   Vect *centroid    = NewVect(ndim,NULL,initScore);
   Vect *trans_vec   = NewVect(ndim,NULL,initScore);
   Vect *refl_point  = NewVect(ndim,NULL,initScore);
   Vect *expd_point  = NewVect(ndim,NULL,initScore);
   Vect *cont_point  = NewVect(ndim,NULL,initScore);

   /* allocate the (ndim+1) vertices of the simplex */
   vertices = xMalloc(sizeof(Vect *)*nv);

   /* initialize the (ndim+1) vertices of the simplex */
   for (size_t v = 0 ; v < nv ; v++) {
      vertices[v] = NewVect(ndim,startingPoints[v],initScore);
   }

#ifdef MINIM_DEBUG
   printf("Initialize the search with simplex:\n");
   for (int i=0 ; i<nv ; i++) {
      PrintVect(vertices[i]);
   }
#endif

   /* Evaluate each vertex */
   for (int i=0 ; i<nv ; i++) {
      ok=EvaluateVect(vertices[i],func);
      if (ok != 0) {
         printf("ERROR: could not evaluate vertex #%d.\n",i);
         goto ErrorFinish;
      }
   }



   /* Check the convergence condition at end of the loop */
   unsigned int count_iterations=0;
   while (1) {
      
#ifdef MINIM_DEBUG
      printf("**************** Iteration #%u\n",count_iterations);
#endif

#ifdef MINIM_DEBUG
      printf("Before ranking the vertices:\n");
      for (int i=0 ; i<nv ; i++) {
	 PrintVect(vertices[i]);
      }
#endif
      
      /* Rank the vertices according to there score */
      RankByScores(nv,vertices,swap,isBetter);
      
#ifdef MINIM_DEBUG
      printf("\tRanked vertices:\n");
      for (int i=0 ; i<nv ; i++) {
	 PrintVect(vertices[i]);
      }
      printf("\tBest   is  : ");PrintVect(vertices[best]);
      printf("\tSWorst is  : ");PrintVect(vertices[sworst]);
      printf("\tWorst  is  : ");PrintVect(vertices[worst]);
#endif
     
      /* check the convergence as the score distances between worst and best */
      if (check_tol((const double)(vertices[worst]->score),(const double)(vertices[best]->score))) {
	 success=FOUND_MINIMUM;
	 break;
      }
     
      /* get the centroid point of all-except-the-worst. We use the fact that the worst is known
       * to be the last element in the array */
      GetCentroid(nv-1,(const Vect **)vertices,centroid);
#ifdef MINIM_DEBUG
      printf("\tCentroid : ");PrintCoordVect((const Vect *)centroid);
#endif
      /* reflec the worst vertex through the centroid point */
      /* a) compute the 'shift vector' from worst towards centroid */
      SubVect((const Vect *)centroid,(const Vect *)(vertices[worst]),trans_vec);
      
      /* b) compute the reflected point*/
      AssignVect(refl_point,trans_vec);
      AddScaleVect((const Vect *)(vertices[worst]),refl_point,2.0); 
      
      /* if the reflected point is not valid, we will try an inside contraction instead */
      if (ValidVect((const Vect *)refl_point)) {
	 /* evaluate this new point */
         ok=EvaluateVect(refl_point,func);
         if (ok != 0) {
            printf("ERROR: could not evaluate a vertex.\n");
            goto ErrorFinish;
         }
#ifdef MINIM_DEBUG
         printf("\tReflected point : ");PrintVect((const Vect *)refl_point);
#endif
	 /* if the reflected point improves on second worst (but *not* on best),
	  * we store it and go to next iteration. */
	 //if ( (refl_point->score >= vertices[best]->score) && (refl_point->score < vertices[sworst]->score) ) {
	 if ( !(*isBetter)(refl_point->score,vertices[best]->score) 
	       && (*isBetter)(refl_point->score,vertices[sworst]->score) ) {
	    AssignVect((vertices[worst]),(const Vect *)refl_point);
#ifdef MINIM_DEBUG
	    printf("\tKeep 'Reflected'.\n");
#endif
	    goto NextIteration;
	 } /* endif  */
	 else
	 /* if the reflected point improves on the current best point, try a step further : expansion */
	 if ((*isBetter)(refl_point->score,vertices[best]->score)) {
	    AssignVect(expd_point,trans_vec);
            AddScaleVect((const Vect *)(vertices[worst]),expd_point,3.0);
#ifdef MINIM_DEBUG
            printf("\t\tReflected point is better than best : EXPANSION to ");PrintCoordVect((const Vect *)expd_point);
#endif          

	    if (ValidVect((const Vect *)expd_point)) {
	       /* evaluate this new expansion point */
               ok=EvaluateVect(expd_point,func);
               if (ok != 0) {
                  printf("ERROR: could not evaluate a vertex.\n");
                  goto ErrorFinish;
               }

#ifdef MINIM_DEBUG
               printf("\t\tValid EXPANSION : ");PrintVect((const Vect *)expd_point);
#endif           
	       if ( (*isBetter)(expd_point->score,refl_point->score) ) {
		  /* expansion point is better than the reflected point */
#ifdef MINIM_DEBUG
                  printf("\t\tKeep 'Expansion' point.\n");
#endif           
		  AssignVect(vertices[worst],(const Vect *)expd_point);
		  goto NextIteration;
	       }
	    } 
	    
	    /* the new expansion point is not valid OR was found worse than the reflected point */
#ifdef MINIM_DEBUG
            printf("\tWe keep the valid reflection point.\n");
#endif           
	    AssignVect(vertices[worst],(const Vect *)refl_point);
	    goto NextIteration;
	 
	 } /* endif. */ 
	 else
	 /* the reflected point is worse that the second worst: */
	 if (!(*isBetter)(refl_point->score,vertices[sworst]->score)) {

	    /* go for a contraction */

	    double contraction_factor;
	    /* choose between INSIDE or OUTSIDE contraction */
	    if ( (*isBetter)(refl_point->score,vertices[worst]->score) ) {
	       /* go for an outisde contraction */
	       contraction = OUTSIDE_CONTRACTION;
	       contraction_factor = 1.5;
	    } else {
	       /* go for an inside contraction */
	       contraction = INSIDE_CONTRACTION;
	       contraction_factor = 0.5;
	    }

	    /* do and test the contraction */
	    AssignVect(cont_point,trans_vec);
	    AddScaleVect((const Vect *)vertices[worst],cont_point,contraction_factor); 
#ifdef MINIM_DEBUG
	    printf("\t\tReflected point is worse than sworst but better than worst: %s CONTRACTION to ",
		  (contraction==OUTSIDE_CONTRACTION?"OUTSIDE":"INSIDE"));
	    PrintCoordVect((const Vect *)cont_point);
#endif          

	    if (ValidVect((const Vect *)cont_point)) {
	       /* evaluate this new expansion point */
	       ok=EvaluateVect(cont_point,func);
	       if (ok != 0) {
		  printf("ERROR: could not evaluate a vertex.\n");
		  goto ErrorFinish;
	       }

#ifdef MINIM_DEBUG
	       printf("\t\tValid CONTRACTION : ");PrintVect((const Vect *)cont_point);
#endif          
	       
	       if (contraction == OUTSIDE_CONTRACTION) {
		  if ( (*isBetter)(cont_point->score,refl_point->score)
			|| (cont_point->score == refl_point->score) ) {
#ifdef MINIM_DEBUG
		     printf("\t\tKeep OUTSIDE CONTRACTION point.\n");
#endif           
		     AssignVect(vertices[worst],(const Vect *)cont_point);
		     goto NextIteration;
		  }
	       } else { 
		  if ( (*isBetter)(cont_point->score,vertices[worst]->score) ) {
#ifdef MINIM_DEBUG
		     printf("\t\tKeep INSIDE CONTRACTION point.\n");
#endif           
		     AssignVect(vertices[worst],(const Vect *)cont_point);
		     goto NextIteration;
		  }
	       }
	    }
	 } /* endif contractions */
      } /* endif validVect(refl_point) */
      
      /* We arrive here because (exclusive) :
       * 1) the refl_point was not valid
       * 2) a contraction was unsuccessfull (non valid or not improving)
       *
       * In both case, we try a multidimensional SHRINK
       */
      contraction=SHRINK_CONTRACTION;

#ifdef MINIM_DEBUG
      printf("\tPerform a SHRINK contraction.\n");
#endif           
      /* Try a multidimensional SHRINK: modify all points but the best */
      for (size_t p = 1 ; p < nv ; p++) {
	 SubVect((const Vect *)vertices[p],(const Vect *)vertices[best],vertices[p]);
	 AddScaleVect((const Vect *)vertices[best],vertices[p],0.5);
	 /* here, we should check the validity of each newly created point */
      }
      
NextIteration:
      count_iterations++;
      
      if (count_iterations == ITERA_MAX) {
         printf("The maximum number of iterations (%u) is reached. The best point up to now is ",count_iterations);
	 PrintVect((const Vect *)vertices[best]);
         success=REACHED_ITERA_MAX;
	 break; 
      }
   }


   if (success == FOUND_MINIMUM) {
      ;//printf("The minimum was found (after %u iterations): ",count_iterations);PrintVect((const Vect *)vertices[best]);
   }

   
Finish:


   ExtractX((const Vect *)vertices[best],x);
   *score = vertices[best]->score;

   /* free the nv vertices */
   for (int i=0; i<nv ; i++) FreeVect(vertices[i]);
   
   /* free the temporary points */
   FreeVect(swap);
   FreeVect(centroid);
   FreeVect(trans_vec);
   FreeVect(refl_point);
   FreeVect(expd_point);
   FreeVect(cont_point);

   return success;
ErrorFinish:
   success=1;
   goto Finish;
}

int check_tol(const double fworst,const double fbest) {
   double delta = fabs(fworst - fbest);
   double accuracy = (fworst + fbest)*FTOL;

   return (delta < (accuracy + FEPS));
}

