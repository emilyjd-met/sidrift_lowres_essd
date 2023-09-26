

#define VERY_HIGH_VALUE  1.e50
#define VERY_LOW_VALUE  -1.e50

#define deg2rad(x) (double)(x*PI/180.)
#define rad2deg(x) (double)(x*180./PI)

#define EPSILON 0.001
#define are_equals(x,y) (fabs(x-y)<=EPSILON*fabs(x))

#ifndef MINIM_DEBUG
#undef MINIM_DEBUG
#endif
