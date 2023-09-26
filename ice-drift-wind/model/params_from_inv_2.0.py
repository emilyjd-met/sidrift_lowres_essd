# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 01:58:17 2017

@author: javimozo
"""

import os
import sys
import re
from datetime import timedelta
import argparse
from argparse import RawDescriptionHelpFormatter
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from numpy.linalg import inv
from pyresample import geometry, utils, AreaDefinition, kd_tree
try:
    from pyresample.area_config import parse_area_file
except ImportError as e:
    from pyresample.utils import parse_area_file
from time_periods import period_def
from write_params_athc import write_params
sys.path.append(os.path.join(os.path.dirname(__file__), '../../merge/src'))
from icedrift_flags import id_flags

area_long = {'osi405_nh': 'nh-polstere-625',
             'osi405_sh': 'sh-polstere-625',
             'osi455_nh': 'nh-ease2-750',
             'osi455_sh': 'sh-ease2-750'}

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')


def parse_args():

    valid_period = ["summer", "winter", "jan", "feb", "mar", "apr", "may",
                    "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    valid_icefmt = ['proc', 'final']

    p = argparse.ArgumentParser('params_from_inv',
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-s', '--start_yr',
                   help='starting date (integer of format YYYY) for computing '
                   'parameters, can be a comma-separated list')
    p.add_argument('-e', '--stop_yr',
                   help='ending date (integer of format YYYY) for computing '
                   'parameters, if start_yr is a comma-separated list, this '
                   'should be too')
    p.add_argument('-p', '--period', choices=valid_period,
                   help='Month or season, valid choices are {}'.format(valid_period))
    p.add_argument('-a', '--area',
                   help='area for which the parameters are calculated '
                   'e.g. "osi405_nh"')
    p.add_argument('-d', '--duration', default=2,
                   help='Duration of ice drift in days, default=2')
    p.add_argument('-z', '--inst', default='multi-oi',
                   help='Platform and instrument, e.g. ssmi-f10, amsr2-gw1 for single sensor, can be a comma-separated list. Default is multi-oi')
    p.add_argument('-l', '--wavelen', default='19',
                   help='Wavelength if single sensor, default is 19GHz. If instrument is a comma-separated list, this should be too.')
    p.add_argument('-w', '--wind_dir', help='where to look for wind files')
    p.add_argument('-i', '--ice_dir', help='where to look for ice files, must be a comma-separated list if instrument is a comma-separated list')
    p.add_argument('-n', '--nodriftmask', required=False, default=None,
                   help='Location of file to mask outside maximum ice drift region')
    p.add_argument('-f', '--icefmt', required=False, default='proc',
                   choices=valid_icefmt,
                   help='Format of ice files, default is proc, valid choices are {}'.format(valid_icefmt))
    p.add_argument('-o', '--out_dir', help='where to write output file')

    args = p.parse_args()

    return args


def load_files(start_yr, stop_yr, period, area, duration, wind_dir, ice_dir,
               inst='multi-oi', wavelen='19', icefmt='proc'):

    wind_files = []
    ice_files = []

    startyrlist = start_yr.split(',')
    stopyrlist = stop_yr.split(',')
    yearlist = []
    for starty, stopy in zip(startyrlist, stopyrlist):
        yearlist.extend(list(range(int(starty), int(stopy)+1)))

    instlist = inst.split(',')
    wavelenlist = wavelen.split(',')
    icedirlist = ice_dir.split(',')
    for year in yearlist:

        start_date, end_date = period_def(period, year)

        # The dates are taken as the final dates of the drift period
        ddate = start_date
        while ddate <= end_date:

            d1 = ddate.strftime("%Y%m%d")
            mdate = ddate + timedelta(days=-duration) # Start of period
            d0 = mdate.strftime("%Y%m%d")
            ddate += timedelta(days=1) # Iterator date

            windf = 'NWP_{}_aggr_{}12-{}12.nc'.format(area, d0, d1)
            windf = os.path.join(wind_dir, str(d1[0:4]), str(d1[4:6]), windf)

            # Finding the ice file(s)
            icefls = []
            for i, ins in enumerate(instlist):
                if ins == 'multi-oi':
                    if icefmt == 'final':
                        icef = 'icedrift_{}_multi-oi_{}1200-{}1200.nc'.format(
                            area_long[area], d0, d1)
                    else:
                        icef = 'icedrift_multi-oi_simplex_lev3_{}_{}12w_{}12w.nc'.format(area_long[area], d0, d1)
                else:
                    icef = 'icedrift_{}_tb{}v-Lap+tb{}h-Lap_simplex_lev2_{}_{}12w_{}12w.nc'.format(ins, wavelenlist[i], wavelenlist[i], area_long[area], d0, d1)

                icef = os.path.join(icedirlist[i], str(d1[0:4]), str(d1[4:6]),
                                    str(d1[6:8]), icef)
                if os.path.isfile(icef):
                    icefls.append(icef)

            # Both wind and ice files should be present for each date
            if os.path.isfile(windf) and icefls:
                for icef in icefls:
                    wind_files.append(windf)
                    ice_files.append(icef)
            elif not os.path.isfile(windf) and icefls:
                print("WARNING: Wind file {} not present but ice file(s) {} present: Skipping".format(windf, icefls))
            elif os.path.isfile(windf) and not icefls:
                print("WARNING: Ice file(s) {} not present but wind file {} present: Skipping".format(icefls, windf))
            else:
                print("WARNING: Neither Wind file {} nor ice file(s) {} present: Skipping".format(windf, icefls))

    # print("wind files = ", wind_files)
    # print("ice files = ", ice_files)

    return (wind_files, ice_files, yearlist[0], yearlist[-1])


def load_data(wind_files, ice_files, duration, icefmt='proc'):
    wind_files_size = np.size(wind_files)
    ice_files_size = np.size(ice_files)

    with Dataset(wind_files[0], mode='r') as datum1:
        xdim = datum1.dimensions['xc'].size
        ydim = datum1.dimensions['yc'].size
    # 2D array to save LOCAL (X) WIND VELOCITY values
    wu = np.zeros([xdim * ydim, wind_files_size], dtype=np.float)
    # 2D array to save MERIDIONAL (Y) WIND VELOCITY values
    wv = np.zeros([xdim * ydim, wind_files_size], dtype=np.float)
    # 2D array to save LOCAL (X) ICE VELOCITY values
    iu = ma.zeros([xdim * ydim, ice_files_size], dtype=np.float)
    # 2D array to save MERIDIONAL (X) ICE VELOCITY values
    iv = ma.zeros([xdim * ydim, ice_files_size], dtype=np.float)

    # Counters
    n = -1
    nn = -1

    for windf in wind_files:
        print("Processing wind file: ", windf)
        with Dataset(windf, mode='r') as datum1:
            # ERAI
#            lat1 = datum1.variables['lat1'][:]
#            lon1 = datum1.variables['lon1'][:]
#            wu_data = datum1.variables['u_wind_avg'][:]
#            wv_data = datum1.variables['v_wind_avg'][:]
            # ERA5
            lat1 = datum1.variables['lat'][:]
            lon1 = datum1.variables['lon'][:]
            wu_data = datum1.variables['uwind_avg'][:]
            wv_data = datum1.variables['vwind_avg'][:]

        gsiz = lat1.size
        gshp = lat1.shape

        n += 1
        wu[:,n] = np.ravel(wu_data)
        wv[:,n] = np.ravel(wv_data)

    for icef in ice_files:
        print("Processing ice file: ", icef)
        with Dataset(icef, mode='r') as datum2:
#            lat1 = datum2.variables['lat1'][:]
#            lon1 = datum2.variables['lon1'][:]
            lat1 = datum2.variables['lat'][:]
            lon1 = datum2.variables['lon'][:]
            if 'dX' in datum2.variables:
                iu_data = datum2.variables['dX'][:]
                try:
                    iv_data = datum2.variables['dY_v1p4'][:]
                except KeyError:
                    iv_data = datum2.variables['dY'][:]
            elif 'driftX' in datum2.variables:
                iu_data = datum2.variables['driftX'][:]
                iv_data = datum2.variables['driftY'][:]

        # Proc files have the y-drift inverted
        if icefmt == 'proc':
            iv_data = iv_data * -1.0

# UNITS
#        # Convert to m/s from km/duration in days before transforming to columns
#        iu_data = iu_data[0, :, :] * (1000.0 / (2 * 86400))
#        iv_data = iv_data[0, :, :] * (1000.0 / (2 * 86400))
        iu_data = iu_data[0, :, :] * (1000.0 / (duration * 86400))
        iv_data = iv_data[0, :, :] * (1000.0 / (duration * 86400))
#        iu_data = iu_data[0, :, :]
#        iv_data = iv_data[0, :, :]
        nn += 1
        # Saves LOCAL (X) ICE VELOCITY values in columns
        iu[:,nn] = ma.ravel(iu_data)
        # Saves MERIDIONAL (Y) ICE VELOCITY values in columns
        iv[:,nn] = ma.ravel(iv_data)

    return (wu, wv, iu, iv, xdim, ydim)


def inversion(wu, wv, iu, iv, period, xdim, ydim, hemi):

    # 2D array to save A and C values as they're computed for each grid point
    # (rows) for every year (columns)
    c = np.zeros([xdim * ydim], dtype=np.complex)
    a = np.zeros([xdim * ydim], dtype=np.complex)
    # 2D array to save LOCAL (X) residual values
    rms_res_r = np.zeros([xdim * ydim],dtype=np.float)
    # 2D array to save MERIDIONAL (Y) residual values
    rms_res_i = np.zeros([xdim * ydim],dtype=np.float)

    # Initialise the arrays with nans
    c[:].real = np.nan
    c[:].imag = np.nan
    a[:].real = np.nan
    a[:].imag = np.nan
    rms_res_r[:] = np.nan
    rms_res_i[:] = np.nan

    # Loop over every gridpoint
    for element in range(ydim * xdim):

        # Extracts one by one every row (which represents wind values
        # for each grid point)
        wu_pt = wu[element,:]
        wv_pt = wv[element,:]
        # Extracts one by one every row (which represents ice drift
        # values for each grid point)
        iu_pt = iu[element,:]
        iv_pt = iv[element,:]

        # Form matrices
        wc = np.array(wu_pt, dtype=np.complex)
        wc.imag = (wv_pt)
        ic = np.array(iu_pt, dtype=np.complex)
        ic.imag = (iv_pt)

        # Ignores the masked values
        k = ma.getmask(iu_pt)
        wc = wc[~k]
        ic = ic[~k]

        # Minimum length of time series to compute for NH:
        # 5 in summer (counted as May-Sept in NH, )
        # 10 in winter (counted as Oct-Apr)
        summernh = ['summer', 'may', 'jun', 'jul', 'aug', 'sep']
        winternh = ['winter', 'oct', 'nov', 'dec', 'jan', 'feb', 'mar', 'apr']
        summersh = ['summer', 'oct', 'nov', 'dec', 'jan', 'feb']
        wintersh = ['winter', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep']
        if hemi == 'nh':
            winterlist = winternh
            summerlist = summernh
        elif hemi == 'sh':
            winterlist = wintersh
            summerlist = summersh
        if period in summerlist:
            min_len_ts = 5
        elif period in winterlist:
            min_len_ts = 10
        else:
            raise(ValueError, "Unrecognised time period: {}".format(period))

        if len(ic) > min_len_ts:
            idmat = np.eye(len(wc))
            # Transforms the horizontal array into a vertical one
            ic = ic[:, np.newaxis]
            # Attaches an array of ones to the list of complex
            # velocity values
            wc = [wc, np.ones_like(wc, dtype=np.complex)]
            # Converts a list of a 1D list and 1D array into a
            # 2 row matrix
            wc = np.mat(wc)
            # Transposes the matrix into a 2 column one (1 x 2)
            wc = wc.T
            # Transforms it into a Hermitian (conjugate
            # transpose) matrix
            wch = wc.getH()
            # Computes the inversion 1x2 matrix with values
            # of A and C (in complex form)
            ac = (inv(wch * wc) * wch * ic)

            # Saves A values by columns
            a[element] = ac[0,0]
            # Appends C values by columns
            c[element] = ac[1,0]

            # Calculates the residuals
            res = (idmat - (wc * (inv(wch * wc)) * wch)) * ic
            rms_res_r[element] = np.sqrt(np.mean(np.square(res.real)))
            rms_res_i[element] = np.sqrt(np.mean(np.square(res.imag)))

    return (a, c, rms_res_r, rms_res_i)


def postprocess(a, c, rms_res_r, rms_res_i, xdim, ydim):

    # |A| informs of coupling between wind and ice (and internal ice stresses)
    abs_a = np.absolute(a)
    # Reshapes the A parameter 1D array into 2D
    abs_a = np.reshape(abs_a, (ydim, xdim))
    # Angle between wind and ice motion vectors
    theta_a = np.angle(a, deg=True)
    theta_a = np.reshape(theta_a, (ydim, xdim))

    # Normalized vectors
    modulus = abs(c)
    c_mod = c/modulus
    c_mod = np.reshape(c_mod, (ydim, xdim))

    # Reshaping to 2-D
    a_2d = np.reshape(a, (ydim, xdim))
    c_2d = np.reshape(c, (ydim, xdim))
    rms_res_r_2d = np.reshape(rms_res_r, (ydim, xdim))
    rms_res_i_2d = np.reshape(rms_res_i, (ydim, xdim))

    # Decomposing the real and imaginary parts into their separate arrays
    real_c = c_2d.real
    imag_c = c_2d.imag
    real_a = a_2d.real
    imag_a = a_2d.imag
    real_c_mod = c_mod.real
    imag_c_mod = c_mod.imag

    return (real_a, imag_a, real_c, imag_c, real_c_mod, imag_c_mod, abs_a,
            theta_a, rms_res_r_2d, rms_res_i_2d)


def grid_landmask(monthnum, area):

    # Area definition of the parameter file
    grid_def_file = os.path.join(os.path.dirname(__file__), 'plots/grids_py.def')
    area_def_list = parse_area_file(grid_def_file)
    area_defs = {}
    for x in area_def_list:
        area_defs[x.area_id] = x
    if area == 'osi405_nh':
        param_area_def = area_defs['nh625']
    elif area == 'osi405_sh':
        param_area_def = area_defs['sh625']
    elif area == 'osi455_nh':
        param_area_def = area_defs['nhease750']
    elif area == 'osi455_sh':
        param_area_def = area_defs['shease750']

    # Read in and regrid a land mask
    landmaskbyarea = {'osi405_nh': (os.path.join(datadir, 'landmask/oceanmask_nh_stere_100.nc'), 'oceanmask', 0, 3),
                      'osi405_sh': (os.path.join(datadir, 'landmask/oceanmask_sh_stere_100.nc'), 'oceanmask', 0, 3),
                      'osi455_nh': (os.path.join(datadir, 'landmask/LandOceanLakeMask_nh-ease2-125.nc'), 'climatology', 0, 'mask'),
                      'osi455_sh': (os.path.join(datadir, 'landmask/LandOceanLakeMask_sh-ease2-125.nc'), 'climatology', 0, 'mask')}
    lmask = landmaskbyarea[area][0]
    datavar = landmaskbyarea[area][1]
    oceanval = landmaskbyarea[area][2]
    landval = landmaskbyarea[area][3]

    with Dataset(lmask, 'r') as dataset:

        # Area definition of the land mask
        grid_mapping = dataset.variables[datavar].__dict__['grid_mapping']
        proj4_string = dataset.variables[grid_mapping].__dict__['proj4_string']
        try:
            proj_dict = utils._proj4.proj4_str_to_dict(proj4_string)
        except:
            proj_dict = utils.proj4.proj4_str_to_dict(proj4_string)
        area_extent = (float(dataset['xc'][0]*1000.),
                       float(dataset['yc'][-1]*1000.),
                       float(dataset['xc'][-1]*1000.),
                       float(dataset['yc'][0]*1000.))
        lmask_area_def = AreaDefinition('data', 'data', 'data', proj_dict,
                                        dataset['xc'].shape[0],
                                        dataset['yc'].shape[0],
                                        area_extent)
        try:
            lmask_data = dataset[datavar][monthnum-1, :, :]
        except:
            lmask_data = dataset[datavar][:, :]

    # Regridding
    default_nbs = 1
    default_rad = 25000
    input_index, output_index, index_array, distance_array = \
        kd_tree.get_neighbour_info(lmask_area_def, param_area_def,
                                      default_rad, neighbours=default_nbs,
                                      reduce_data=False)

    lmask = kd_tree.get_sample_from_neighbour_info('nn',
            param_area_def.shape, lmask_data, input_index, output_index,
            index_array, distance_array=distance_array)

    ocean = lmask == oceanval
    if landval == 'mask':
        land = lmask.mask
    else:
        land = lmask == landval
#    poss_ice = np.logical_and(lmask != oceanval, lmask != landval)
    poss_ice = np.logical_and(~ocean, ~land)
    lmask[ocean] = 0
    lmask[land] = 1
    lmask[poss_ice] = 2

    return lmask, param_area_def


def extrapolate(real_a, imag_a, real_c, imag_c, real_c_mod, imag_c_mod, abs_a,
                theta_a, period, area, nodriftmask=None):

    fillval = 9.95e36

    fields = {0: real_a,
              1: imag_a,
              2: real_c,
              3: imag_c,
              4: real_c_mod,
              5: imag_c_mod,
              6: abs_a,
              7: theta_a
          }

    flagfields = {}
    gapfillfields = {}

    # Reading in the landmask
    startd, endd = period_def(period, 2020) # The year doesn't matter
    monthnum = startd.strftime('%m')
    monthnum = int(monthnum)
    lmask, area_def = grid_landmask(monthnum, area)
    lons, lats = area_def.get_lonlats()

    # Reading in the "no drift" mask
    if nodriftmask is not None:
        with Dataset(nodriftmask, 'r') as dataset:
            nodriftpix = dataset['maxdriftmask'][:, :] == 0

    # General flag fields
#    landpix = lmask == 3
#    noicepix = lmask == 0
    landpix = lmask == 1
    noicepix = lmask == 0
    fillpix = np.logical_or(landpix, noicepix)
    if nodriftmask is not None:
        fillpix = np.logical_or(fillpix, nodriftpix)
    basic_flag = np.empty_like(fields[1], dtype=np.int8)
    # Flag values are defined in accordance with those from cross-correlation
    # icedrift
    if nodriftmask is not None:
        basic_flag[nodriftpix] = 20
    basic_flag[landpix] = 3
    basic_flag[noicepix] = 4

    # Extrapolation parameters for gapfill region and data smoothing
    roi_gf = 1000000
    sig_gf = 0.5 * roi_gf
    roi_data = 125000
    sig_data = 0.5 * roi_data

    # Looping over each parameter
    for i in range(8):

        # Specific flag fields
        datapix = fields[i] <= fillval
        # Mask the pole hole
        datapix = np.logical_and(datapix, lats < 89.5)
        emptypix = np.logical_and(~fillpix, ~datapix)
        flagfields[i] = basic_flag
        flagfields[i][datapix] = id_flags['icedrift_ok']
        flagfields[i][emptypix] = id_flags['icedrift_windgapfill']

        # Swath coverage information
        datalons = lons[datapix]
        datalats = lats[datapix]
        coverage = geometry.SwathDefinition(lons=datalons, lats=datalats)

        # Extrapolation
        datacover = fields[i][datapix]
        gapfillfields[i] = kd_tree.resample_gauss(coverage, datacover,
                                                  area_def,
                                                  radius_of_influence=roi_gf,
                                                  sigmas=sig_gf)
        smooth_data = kd_tree.resample_gauss(coverage, datacover, area_def,
                                             radius_of_influence=roi_data,
                                             sigmas=sig_data)
        gapfillfields[i][datapix] = smooth_data[datapix]
        gapfillfields[i] = ma.array(gapfillfields[i])
        gapfillfields[i].mask = fillpix


    # Checking the flag fields are all the same
    flagsame = True
    for i in [x + 1 for x in range(7)]:
        if not np.all(flagfields[i] == flagfields[0]):
            print("flagfields[{}] differs from flagfields[0]".format(i))
            flagsame = False
    if flagsame:
        print("Flags match.")

    return gapfillfields[0], gapfillfields[1], gapfillfields[2], gapfillfields[3], gapfillfields[4], gapfillfields[5], gapfillfields[6], gapfillfields[7], flagfields[0]


def compute_params(start_yr, stop_yr, period, area, duration, wind_dir,
                   ice_dir, inst='multi-oi', wavelen='19', icefmt='proc',
                   nodriftmask=None, out_dir='.'):

    if re.search('nh', area):
        hemi = 'nh'
    elif re.search('sh', area):
        hemi = 'sh'
    else:
        print("WARNING: Unrecognised hemisphere from area {}. Defaulting "
              "to 'nh'".format(area))
        hemi = 'nh'

    wind_files, ice_files, starty, stopy = load_files(start_yr, stop_yr,
                                                      period, area, duration,
                                                      wind_dir, ice_dir,
                                                      inst=inst,
                                                      wavelen=wavelen,
                                                      icefmt=icefmt)
    wu, wv, iu, iv, xdim, ydim = load_data(wind_files, ice_files, duration)
    a, c, rms_res_r, rms_res_i = inversion(wu, wv, iu, iv, period, xdim, ydim,
                                           hemi)
    real_a, imag_a, real_c, imag_c, real_c_mod, imag_c_mod, abs_a, theta_a, rms_res_r_2d, rms_res_i_2d = postprocess(a, c, rms_res_r, rms_res_i, xdim, ydim)

    real_a_gapfill, imag_a_gapfill, real_c_gapfill, imag_c_gapfill, real_c_mod_gapfill, imag_c_mod_gapfill, abs_a_gapfill, theta_a_gapfill, flags = extrapolate(real_a, imag_a, real_c, imag_c, real_c_mod, imag_c_mod, abs_a, theta_a, period, area, nodriftmask=nodriftmask)

    dfile = write_params(starty, stopy, period, duration, area, real_a,
                         imag_a, real_c, imag_c, real_c_mod, imag_c_mod,
                         abs_a, theta_a, rms_res_r_2d, rms_res_i_2d,
                         real_a_gapfill, imag_a_gapfill, real_c_gapfill,
                         imag_c_gapfill, real_c_mod_gapfill,
                         imag_c_mod_gapfill, abs_a_gapfill, theta_a_gapfill,
                         flags, out_dir)
    print('Output parameter file: {}'.format(dfile))


def main():

    args = parse_args()

    start_yr = args.start_yr
    stop_yr = args.stop_yr
    period = args.period
    area = args.area
    duration = int(args.duration)
    inst = args.inst
    wavelen = args.wavelen
    wind_dir = args.wind_dir
    ice_dir = args.ice_dir
    nodriftmask = args.nodriftmask
    icefmt = args.icefmt
    out_dir = args.out_dir

    compute_params(start_yr, stop_yr, period, area, duration, wind_dir,
                   ice_dir, inst=inst, wavelen=wavelen, icefmt=icefmt,
                   nodriftmask=nodriftmask, out_dir=out_dir)


if __name__ == '__main__':

    main()
