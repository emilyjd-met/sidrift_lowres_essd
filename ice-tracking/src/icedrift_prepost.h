#ifndef ICEDRIFT_PREPOST_H
#define ICEDRIFT_PREPOST_H

#define FFLUSH fflush(stdout)

/* Fill values copied from osibinio.h (Note: these differs from 
   osisaf_ice_netcdf.h which uses reference to INT32_MAX) */
#define UNDEFNC_FLOAT  -1.0e10  
#define UNDEFNC_DOUBLE -1.0e10
#define UNDEFNC_INT    -1e10 
#define UNDEFNC_SHORT  -32767

extern void initmod(int *n, double x[]);

void set_sigmoidlength(double len);

void setBestKnowledge(double lat, double lon);

int setNeighbourhoodPattern(double pattern_radius,size_t *lpattern_size,size_t *lpattern_center,short **lpattern_mask,long *lpattern_windexes[3]);

#endif 
