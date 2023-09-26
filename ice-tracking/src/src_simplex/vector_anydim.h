
#ifndef VECTOR_ANYDIM_H
#define VECTOR_ANYDIM_H

struct vector_ {
   size_t nx;
   double *x;
   double score;
};

typedef struct vector_ Vect;

extern Vect *NewVect(size_t,double *,double);
extern void  AllocVect(size_t,Vect *);
extern Vect *CopyVect(Vect *);
extern int   EvaluateVect(Vect *,int (*func)(size_t,double *,double *));
extern void  AssignVect(Vect *,const Vect *);
extern void  PrintVect(const Vect *);
extern void  PrintCoordVect(const Vect *);
extern int   ValidVect(const Vect *);
extern void  FreeVect(Vect *);
extern void  ExtractX(const Vect *,double *);
extern void  GetCentroid(size_t,const Vect **,Vect *);
extern void  SubVect(const Vect *,const Vect *,Vect *);
extern void  AddScaleVect(const Vect *,Vect *,const double);
extern void  RankByScores(size_t, Vect **,Vect *,int (*)(double,double));
/* ranking routines */
extern int better_if_lower (double s1,double s2);
extern int better_if_higher(double s1,double s2);

#endif /* VECTOR_ANYDIM_H */
