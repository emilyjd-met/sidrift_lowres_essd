
#ifndef ICEDRIFT_FLAGS_H
#define ICEDRIFT_FLAGS_H

/* FLAG VALUES FOR DAILY IMAGES (INTERNAL) */
#define TCIMAGE_OUTSIDE_GRID              -2
#define TCIMAGE_NODATA                    -1
#define TCIMAGE_OK                         0
#define TCIMAGE_UNPROCESSED                1
#define TCIMAGE_FAILED                     2

/* FLAG VALUES FOR DAILY PRODUCTS (INTERNAL) */
#define ICEDRIFT_UNPROCESSED              -1
#define ICEDRIFT_OK                        0
#define ICEDRIFT_OUTSIDE_IMGBORDER         1
#define ICEDRIFT_CLOSETO_IMGBORDER         2
#define ICEDRIFT_CENTER_OVER_LAND          3
#define ICEDRIFT_NOICE                     4
#define ICEDRIFT_CLOSETO_COAST_OR_EDGE     5
#define ICEDRIFT_CLOSETO_MISSING           6
#define ICEDRIFT_CLOSETO_UNPROCESSED       7
#define ICEDRIFT_OPTIMISATION_FAILS        8
#define ICEDRIFT_FAILS                     9
#define ICEDRIFT_LOWCORRELATION           10
#define ICEDRIFT_TOOLONG                  11
#define ICEDRIFT_REFUSED_BY_NEIGHBOURS    12
#define ICEDRIFT_CORRECT_BY_NEIGHBOURS    13
#define ICEDRIFT_NOAVERAGE                14
#define ICEDRIFT_SUMMERHOLD               15
#define ICEDRIFT_OIINTERP                 16
#define ICEDRIFT_SMALLER_PATTERN          17
#define ICEDRIFT_MASKED_NWP               18
#define ICEDRIFT_NBFLAGS                  19 

/* FLAG VALUES FOR DAILY UNCERTAINTIES (INTERNAL) */
#define ICEDRIFTPOST_UNPROCESSED    -1
#define ICEDRIFTPOST_OK              0
#define ICEDRIFTPOST_DIF             1
#define ICEDRIFTPOST_INV             2
#define ICEDRIFTPOST_NEG             3
#define ICEDRIFTPOST_NOVECTOR        ICEDRIFTPOST_UNPROCESSED

#endif /* ICEDRIFT_FLAGS_H */

