#ifndef READ_ICEDRIFT_PRODUCT_FILE_H
#define READ_ICEDRIFT_PRODUCT_FILE_H

char sigdXname[10];
char sigdYname[10];

int readProductFile(char infname[],fmsec1970 refT0, fmsec1970 refT1, size_t out_dims[], 
      float A[2], float B[2], char *projstr, PJ **proj,
      fmsec1970 **t0,fmsec1970 **t1,
      float **lat,float **lon,
      float **driftX, float **driftY, short **flags, float *fillvalf,
      float **sigdX, float **sigdY, float **corrdXdY, short **uflags);


#endif /* READ_ICEDRIFT_PRODUCT_FILE_H */
