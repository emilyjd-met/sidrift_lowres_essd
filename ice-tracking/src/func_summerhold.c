
#include <stdio.h>
#include <fmutil.h>
#include "icedrift_instruments.h"
#include "func_summerhold.h"

int test_summer_hold(fmsec1970 startdate, char *area, char *instr, char *channels, int *productdate_in_summer_hold) {

   int ret;

   if (is_MIRAS(instr)) {
      *productdate_in_summer_hold = 0;
      return 0;
   } else if (is_AMSR2(instr) || is_AMSR(instr)) {
      if (strstr(channels,"bt19")) {
         *productdate_in_summer_hold = 0;
         return 0;
      }
   }

   /* transform the given time */
   fmtime fmt;
   tofmtime(startdate,&fmt);

   int year = fmt.fm_year;

   fmtime last_before_summer;
   fmtime first_after_summer;

   if (strstr(area,"nh")) {
      last_before_summer.fm_year = year;
      first_after_summer.fm_year = year;

      last_before_summer.fm_mon   = 5;
      last_before_summer.fm_mday  = 1; 
      last_before_summer.fm_hour  = 0;
      last_before_summer.fm_min   = 0;
      last_before_summer.fm_sec   = 0;

      first_after_summer.fm_mon   = 10;
      first_after_summer.fm_mday  = 1; 
      first_after_summer.fm_hour  = 0;
      first_after_summer.fm_min   = 0;
      first_after_summer.fm_sec   = 0;
   } else if (strstr(area,"sh")) {
      last_before_summer.fm_year = year;
      first_after_summer.fm_year = year;

      first_after_summer.fm_mon   = 4;
      first_after_summer.fm_mday  = 1; 
      first_after_summer.fm_hour  = 0;
      first_after_summer.fm_min   = 0;
      first_after_summer.fm_sec   = 0;

      last_before_summer.fm_mon   = 11;
      last_before_summer.fm_mday  = 1; 
      last_before_summer.fm_hour  = 0;
      last_before_summer.fm_min   = 0;
      last_before_summer.fm_sec   = 0;
   } else {
      fprintf(stderr,"ERROR (%s) Unknown area <%s>\n",__func__,area);
      return 1;
   }

   fmsec1970 last_before_summer_time = tofmsec1970(last_before_summer);
   if (last_before_summer_time == -1) {
      fprintf(stderr,"ERROR (%s) Cannot convert 'last_before_summer' date. Year is %d\n",__func__,year);
      return 1;
   }
   fmsec1970 first_after_summer_time = tofmsec1970(first_after_summer);

   if (first_after_summer_time == -1) {
      fprintf(stderr,"ERROR (%s) Cannot convert 'first_after_summer' date. Year is %d\n",__func__,year);
      return 1;
   }

   if (strstr(area,"nh")) {
      if ((startdate >= last_before_summer_time) && (startdate <= first_after_summer_time)) {
         *productdate_in_summer_hold = 1;
      } else {
         *productdate_in_summer_hold = 0;
      }
   } else {
      if ((startdate >= last_before_summer_time) || (startdate <= first_after_summer_time)) {
         *productdate_in_summer_hold = 1;
      } else {
         *productdate_in_summer_hold = 0;
      }
   }



   return 0;

}

#undef MAINTEST
#ifdef MAINTEST
int main (void) {

   char firstdate[] = "2007010112";
   char lastdate[]  = "2009123112";

   fmsec1970 firsttime = ymdh2fmsec1970_alt(firstdate,0);
   fmsec1970 lasttime  = ymdh2fmsec1970_alt(lastdate,0);

   fmsec1970 currtime = firsttime;
   while (currtime <= lasttime) {

	  char currdate[FMUTIL_CFEPOCH_LENGTH+1];
	  fmsec19702CFepoch(currtime,currdate);

      int ret;
	  int productdate_in_summer_hold_nh;
	  ret = test_summer_hold(currtime,"nh",0,&productdate_in_summer_hold_nh);
	  if (ret) {
		 printf("%s. Skip this date\n",currdate);continue;
	  }
	  int productdate_in_summer_hold_sh;
	  ret = test_summer_hold(currtime,"sh",0,&productdate_in_summer_hold_sh);
	  if (ret) {
		 printf("%s. Skip this date\n",currdate);continue;
	  }

	  printf("%s %d %d\n",currdate,!productdate_in_summer_hold_nh,!productdate_in_summer_hold_sh);

	  currtime += 24*60*60; /* +1 day */
   }

}
#endif /* MAINTEST */
