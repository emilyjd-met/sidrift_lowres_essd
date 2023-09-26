
#include <stdio.h>
#include <stdlib.h>

#include "memory.h"


void *xMalloc(size_t size) {
   void *ret;
   if (!(ret=malloc(size))) {
      fprintf(stderr,"Memory error: could not allocate %lu bytes.\n",(unsigned long)size);
      exit(EXIT_FAILURE);
   }
   return ret;   
}

void *xCalloc(size_t nelem,size_t size) {
   void *ret;
   size_t totsize;

   if ((totsize=(nelem*size)) < size) {
      fprintf(stderr,"Memory error: The size asked to Calloc is too large for size_t type.\n");
      exit(EXIT_FAILURE);
   }
   if (!(ret=calloc(nelem,size))) {
      fprintf(stderr,"Memory error: could not allocate %lu bytes.\n",(unsigned long)totsize);
      exit(EXIT_FAILURE);
   }
   return ret;
}
