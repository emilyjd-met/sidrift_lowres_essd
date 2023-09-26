#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "fmutil.h"

int hasTrailingSlash (const char *string) {

   char *lastSlash=strrchr(string,'/');
   return (lastSlash == &(string[strlen(string)-1]));

}

void deleteTrailingSlash (char **string) {

   while ( hasTrailingSlash(*string) ) {
      int  len = strlen(*string);
      char *tmp = fmMalloc((len-1+1)*sizeof(char));
      strncpy(tmp,*string,len-1);
      tmp[len-1]='\0';
      free(*string);
      *string = tmp;
   }

}

void appendSlash (char **string) {

   int  len = strlen(*string);
   char *tmp = fmMalloc((len+1+1)*sizeof(char));
   sprintf(tmp,"%s/",*string);
   free(*string);
   *string = tmp;

}
