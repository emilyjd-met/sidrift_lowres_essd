"""
Created on Thu Jun 14 16:48:36 2017

@author: javier
"""

import os
import sys
import re
from datetime import datetime, timedelta, time
import argparse
from argparse import RawDescriptionHelpFormatter
from netCDF4 import Dataset
import numpy as np
from numpy import ma
from pyresample import geometry, utils, AreaDefinition
try:
    from pyresample.area_config import parse_area_file
except ImportError as e:
    from pyresample.utils import parse_area_file
import cartopy
import cartopy.crs as ccrs

here = os.path.dirname(os.path.abspath(__file__))
codedir = os.path.join(here, '../..')
sys.path.append(os.path.join(codedir, 'ice-tracking/src'))
from write_icedrift_product_file import write_icedrift

def read_args():
    '''Read and pre-process user input'''


    valid_area = ['osi405_nh', 'osi405_sh', 'osi455_nh', 'osi455_sh']

    p = argparse.ArgumentParser("icedrift_from_winds",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('enddate',
                   help='End date (YYYYMMDD) for running free-drift model')
    p.add_argument('-d', '--duration', default=2,
                   help='Duration of the ice drift.')
    p.add_argument('-a', 'area', choices=valid_area,
                   help='Area for which the parameters are calculated, '
                   'valid choices are {}'.format(valid_area))
    p.add_argument('-w', '--wind_dir',
                   help='Input directory containing wind files.', default=None)
    p.add_argument('-p', '--param_file', default=None,
                   help='Parameter file from inverse calc.')
    p.add_argument('-o', '--out_dir', default='./case',
                   help='Where to write ice drift files.')
    p.add_argument('-s', '--suffix', default=None,
                   help='Optional suffix for the output file')
    args = p.parse_args()

    return args


def find_ccrs(area, area_def):

    if re.search('nh', area):
        hemi = 'nh'
    elif re.search('sh', area):
        hemi = 'sh'
    else:
        raise ValueError("Unrecognised hemisphere in area {}".format(area))
    if re.search('ease', area) or re.search('osi455', area):
        coordsys = 'ease'
    else:
        coordsys = 'polstere'

    try:
        proj_dict = utils._proj4.proj4_str_to_dict(str(area_def.proj4_string))
    except:
        proj_dict = utils.proj4.proj4_str_to_dict(str(area_def.proj4_string))

    if coordsys == 'polstere':
        try:
            data_globe = ccrs.Globe(semimajor_axis=proj_dict['a'],
                                    semiminor_axis=proj_dict['b'])
        except:
            data_globe = ccrs.Globe(semimajor_axis=proj_dict['a'],
                                    inverse_flattening=proj_dict['rf'])
        if hemi == 'nh':
            area_ccrs = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                              globe=data_globe)
        else:
            area_ccrs = ccrs.SouthPolarStereo(central_longitude=0.0,
                                              globe=data_globe)

    else:
        if hemi == 'nh':
            area_ccrs = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=90,
                false_easting=0, false_northing=0)
        else:
            area_ccrs = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=-90,
                false_easting=0, false_northing=0)

    return area_ccrs


def compute_drift(enddate, duration, area, wind_dir, param_file, out_dir='.',
                  suffix=None):

    # Create filename
    edate = datetime.combine(enddate, time(12))
    sdate = edate - timedelta(days=duration)
    windf = 'NWP_{}_aggr_{:%Y%m%d}12-{:%Y%m%d}12.nc'.format(area, sdate, edate)
    windf = os.path.join(wind_dir, datetime.strftime(edate, '%Y'),
                         datetime.strftime(edate, '%m'), windf)

    # Load information from wind file
    with Dataset(windf, mode='r') as datum1:
        # ERA5
        lat1 = datum1.variables['lat'][:]
        lon1 = datum1.variables['lon'][:]
        wu = datum1.variables['uwind_avg'][:]
        wv = datum1.variables['vwind_avg'][:]

    gsiz = lat1.size
    gshp = lat1.shape

    # Applying the model requires complex numbers
    wc_2d = np.array(wu, dtype=np.complex)
    wc_2d.imag = wv

    wc = np.zeros([gsiz,], dtype=np.complex)
    # Initialise with nans
    wc[:,].real = np.nan
    wc[:,].imag = np.nan

    wc[:,] = wc_2d.flat

    # Load the free-drift parameters
    with Dataset(param_file, mode='r') as datum2:
        a_real = datum2.variables['A_real_gapfill'][:]
        a_imag = datum2.variables['A_imag_gapfill'][:]
        c_real = datum2.variables['C_real_gapfill'][:]
        c_imag = datum2.variables['C_imag_gapfill'][:]
        flags = datum2.variables['flags'][:]

    # Recombine real and imaginary parts of both model parameters
    a = np.array(a_real, dtype=np.complex)
    a.imag = a_imag
    c = np.array(c_real, dtype=np.complex)
    c.imag = c_imag

    # 2D array to save A and C values
    ac = np.zeros((gsiz, 2), dtype=np.complex)
    ac[:,:].real = np.nan                           # changes zeros to nans
    ac[:,:].imag = np.nan
    ac[:,0] = a.flat
    ac[:,1] = c.flat

    # 1D array to save iceDrift values
    id_arr = np.zeros([gsiz,], dtype=np.complex)
    id_arr[:,].real = np.nan                        # changes zeros to nans
    id_arr[:,].imag = np.nan

    # Loop through the grid cells and apply the model
    for i in range(gsiz,):         # loop over every gridpoint
         # Extracts one by one every row (which represents wind values
        # for each grid point)
        wc_pt = wc[i,]
        # Attach an array of ones to the list of complex velocity values
        wc_pt = [wc_pt, np.ones_like(wc_pt, dtype=np.complex)]
        wc_pt = np.mat(wc_pt)

        ac_pt = ac[i,:]
        ac_pt = np.mat(ac_pt)
        ac_pt = ac_pt.T

        id_pt = wc_pt * ac_pt
        id_arr[i,] = id_pt[0,0]

    ic = id_arr
    ic = np.reshape(ic, gshp)

    # UNITS
    # Get to the dX dY components
    # tspan_s is seconds in duration days
    tspan_s = duration * 86400.
    dx = np.ma.masked_invalid(ic.real * (tspan_s / 1000.0))
    dy = np.ma.masked_invalid(ic.imag * (tspan_s / 1000.0))

    # Mask also where the inversion parameters flags show "no ice" (4) or
    # land (3)
    flagnoice = np.logical_or(flags == 3, flags == 4)
    dx.mask = np.logical_or(dx.mask, flagnoice)
    dy.mask = np.logical_or(dy.mask, flagnoice)

    # Masking very very low values, and a couple of pixels near the pole
    # with very high values
    dx.mask[dx > 1000.] = True
    dx.mask[dx < -1000.] = True
    dy.mask[dy > 1000.] = True
    dy.mask[dy < -1000.] = True

    # Finding the area definition
    grid_def_file = os.path.join(os.path.dirname(__file__), 'grids_py.def')
    area_def_list = parse_area_file(grid_def_file)
    area_defs = {}
    for x in area_def_list:
        area_defs[x.area_id] = x
    if area == 'osi405_nh':
        area_def = area_defs['nh625']
    elif area == 'osi405_sh':
        area_def = area_defs['sh625']
    elif area == 'osi455_nh':
        area_def = area_defs['nhease750']
    elif area == 'osi455_sh':
        area_def = area_defs['shease750']
    projstr = str(area_def.proj4_string)
    area_ccrs = find_ccrs(area, area_def)

    if area == 'osi405_nh':
        gridname = 'nh-polstere-625'
    elif area == 'osi405_sh':
        gridname = 'sh-polstere-625'
    elif area == 'osi455_nh':
        gridname = 'nh-ease2-750'
    elif area == 'osi455_sh':
        gridname = 'sh-ease2-750'

    uncert_dx = 0.1 * dx
    uncert_dy = 0.1 * dy

    # Time variables
    fillvalf = -1.0e-10
    sizearr = dx.shape
    datamask = dx.mask
    timestart = np.full(sizearr, fillvalf)
    timestart = ma.array(timestart, mask=datamask)
    timestart[~datamask] = 0
    timeend = np.full(sizearr, fillvalf)
    timeend = ma.array(timeend, mask=datamask)
    timeend[~datamask] = 0

    # Beginning and ending lats and lons
    pc = ccrs.PlateCarree()
    lons, lats = area_def.get_lonlats()
    olatb = np.full(sizearr, fillvalf)
    olatb = ma.array(olatb, mask=datamask)
    olatb[~datamask] = lats[~datamask]
    olonb = np.full(sizearr, fillvalf)
    olonb = ma.array(olonb, mask=datamask)
    olonb[~datamask] = lons[~datamask]
    olate = np.full(sizearr, fillvalf)
    olate = ma.array(olate, mask=datamask)
    olone = np.full(sizearr, fillvalf)
    olone = ma.array(olone, mask=datamask)
    dist = np.full(sizearr, fillvalf)
    dist = ma.array(dist, mask=datamask)
    dirn = np.full(sizearr, fillvalf)
    dirn = ma.array(dirn, mask=datamask)
    points = ma.where(datamask==False)
    for ix, iy in zip(points[0], points[1]):
        lon0 = lons[ix, iy]
        lat0 = lats[ix, iy]
        x0, y0 = area_ccrs.transform_point(lon0, lat0, pc)
        x1 = x0 + (dx[ix, iy] * 1000.)
        y1 = y0 - (dy[ix, iy] * 1000.)
        lon1, lat1 = pc.transform_point(x1, y1, area_ccrs)
        olate[ix, iy] = lat1
        olone[ix, iy] = lon1

    # Writing the NetCDF file out
    blank0 = np.zeros(sizearr)
    blank1 = np.full(sizearr, 1)
    blankfill = np.full(sizearr, -32767)
    dfile = write_icedrift(area_def,
                                  None,
                                  None,
                                  None,
                                  sdate,
                                  edate,
                                  gridname,
                                  projstr,
                                  out_dir,
                                  dx,
                                  dy,
                                  timestart,
                                  timeend,
                                  flags,
                                  uncert_dx,
                                  uncert_dy,
                                  blank0,
                                  blankfill,
                                  dist,
                                  dirn,
                                  olonb,
                                  olatb,
                                  olone,
                                  olate,
                                  blank1,
                                  None,
                                  None,
                                  None,
                                  None,
                                  None,
                                  None,
                                  None,
                                  blank0,
                                  method='wind',
                                  suffix=suffix)

    return dfile


def main():

    args = read_args()
    enddate = args.enddate
    duration = int(args.duration)
    area = args.area
    wind_dir = args.wind_dir
    param_file = args.param_file
    out_dir = args.out_dir
    suffix = args.suffix

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    enddate = datetime.strptime(enddate, '%Y%m%d')
    compute_drift(enddate, duration, area, wind_dir, param_file, out_dir,
                  suffix)


if __name__ == '__main__':

    main()
