"""CFFI building code to wrap the ice drift core C code into a library
accessible from a Python wrapper function."""

import os
from cffi import FFI
ffibuilder = FFI()

# Code paths
here = os.path.dirname(os.path.abspath(__file__))
osisaf_path = os.path.join(here, '../../../osisaf-hl-sw')
icedrift_path = os.path.join(here, '../../output/icedrift')
proj_path = os.path.join(here, '../../data/proj-4.5.0_red/src')

ffibuilder.set_source("_idcore", # name of the output C extension
"""
    # include <stdio.h>
    # include <stdlib.h>
    # include <string.h>
    # include <float.h>
    # include <math.h>
    # include <projects.h>
    # include <fmutil.h>
    # include <assert.h>
    # include <errorcodes.h>
    # include <useproj.h>
    # include "icedrift_flags.h"
    # include "icedrift_common.h"
    # include "icedrift_prepost.h"
    # include "icedrift_model.h"
    # include "icedrift_solve_filter.h"
    # include "optimization_simplex.h"
    # include "vector_anydim.h"
    # include "memory.h"
    # include "icedrift_solve_common.h"
    # include "ice_common.h"
    # include "proj_api.h"
    # include "icedrift_solve_core.h"

""",
    include_dirs=[os.path.join(icedrift_path, 'ice-tracking/src'),
                  os.path.join(icedrift_path,
                      'ice-tracking/src/src_simplex'),
                  os.path.join(icedrift_path,
                      'common/src'),
                  os.path.join(osisaf_path,
                      'OSI_HL_Ice/common/libicecommon/include'),
                  os.path.join(osisaf_path,
                      'OSI_HL_common/commondefs'),
                  os.path.join(osisaf_path,
                      'OSI_HL_AUX/libs/libfmutil/include'),
                  os.path.join(osisaf_path,
                               'OSI_HL_AUX/libs/libfmutil/src'),
                  os.path.join(osisaf_path,
                      'OSI_HL_AUX/OSI_HL_FUNC/useproj/include'),
                  # Reduced version of proj4 code
                  proj_path],

    sources=['icedrift_solve_core.c', 'icedrift_prepost.c',
             'icedrift_solve_common.c', 'icedrift_model.c',
             'icedrift_solve_filter.c',
             os.path.join(icedrift_path,
                 'ice-tracking/src/src_simplex/optimization_simplex.c'),
             os.path.join(icedrift_path,
                 'ice-tracking/src/src_simplex/vector_anydim.c'),
             os.path.join(icedrift_path,
                 'ice-tracking/src/src_simplex/memory.c'),
             os.path.join(icedrift_path,
                 'common/src/icedrift_common.c'),
             os.path.join(osisaf_path,
                 'OSI_HL_AUX/libs/libfmutil/src/fmtime.c'),
             os.path.join(osisaf_path,
                 'OSI_HL_AUX/libs/libfmutil/src/fmsolar.c'),
             os.path.join(osisaf_path,
                 'OSI_HL_AUX/libs/libfmutil/src/fmangleconversion.c'),
             os.path.join(osisaf_path,
                 'OSI_HL_AUX/libs/libfmutil/src/fmstorage.c'),
             os.path.join(osisaf_path,
                 'OSI_HL_AUX/libs/libfmutil/src/fmstrings.c'),
             os.path.join(osisaf_path,
                 'OSI_HL_AUX/libs/libfmutil/src/fmerrmsg.c'),
             # Only the proj4 code which it complains without
             os.path.join(proj_path, 'adjlon.c'),
             os.path.join(proj_path, 'dmstor.c'),
             os.path.join(proj_path, 'pj_auth.c'),
             os.path.join(proj_path, 'pj_datums.c'),
             os.path.join(proj_path, 'pj_datum_set.c'),
             os.path.join(proj_path, 'pj_ellps.c'),
             os.path.join(proj_path, 'pj_ell_set.c'),
             os.path.join(proj_path, 'pj_errno.c'),
             os.path.join(proj_path, 'pj_fwd.c'),
             os.path.join(proj_path, 'pj_init.c'),
             os.path.join(proj_path, 'pj_inv.c'),
             os.path.join(proj_path, 'PJ_laea.c'),
             os.path.join(proj_path, 'pj_latlong.c'),
             os.path.join(proj_path, 'pj_list.c'),
             os.path.join(proj_path, 'pj_malloc.c'),
             os.path.join(proj_path, 'pj_open_lib.c'),
             os.path.join(proj_path, 'pj_param.c'),
             os.path.join(proj_path, 'pj_qsfn.c'),
             os.path.join(proj_path, 'PJ_stere.c'),
             os.path.join(proj_path, 'pj_tsfn.c'),
             os.path.join(proj_path, 'pj_units.c')],

    libraries=['c'])


# cdef() expects a single string declaring the C types, functions and
# globals, in valid C syntax.
ffibuilder.cdef("""
                extern double **obs[2];
                extern short **TCflag[2];
                extern short *icelandmask[2];
                extern double *img_lat;
                extern double *img_lon;
                extern double img_Ax;
                extern double img_Bx;
                extern double img_Ay;
                extern double img_By;
                extern size_t img_dims[3];
                extern char *img_projstr;
                extern short nbWaveBands;
                extern short twghtStart;
                extern short twghtEnd;
                extern double maxdriftdistance;
                extern double sigmoid_length;
                extern char *OptimMetric;
                extern double pattern_radius[2];
                extern double radiusNeighbours;
                extern int NDRIFTPIXELS;
                extern char *out_area;
                extern char *out_projstr;
                extern double out_Ax;
                extern double out_Bx;
                extern double out_Ay;
                extern double out_By;
                extern size_t out_dims[3];
                extern double *olat;
                extern double *olon;
                extern unsigned int *owcs;
                extern unsigned int *iwcs;
                extern float *driftX;
                extern float *driftY;
                extern short *pflag;
                extern float *sigdX;
                extern float *sigdY;
                extern float *corrdXdY;
                extern short *uflag;
                extern float *length;
                extern float *dir;
                extern float *outlonB;
                extern float *outlatB;
                extern float *outlonE;
                extern float *outlatE;
                extern float *outfc;
                extern short *snavg;
                extern float *avgX;
                extern float *avgY;
                extern float *length_avg;
                extern float *length_diff;
                extern float *stdX;
                extern float *stdY;
                extern short *patternIndex;
                extern char *reportFile;

                int core(void);
                void *malloc(size_t size);
                void free(void *ptr);
                """)


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
