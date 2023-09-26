# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 23:08:07 2017

@author: admin
"""

import os
import re
from netCDF4 import Dataset
import numpy as np
#from pyresample import geometry
try:
    from pyresample import parse_area_file
except ImportError as e:
    from pyresample.utils import parse_area_file

shapes = {'nh':(177,119), 'sh':(131,125)}


def write_params(start_yr, end_yr, period, duration, area, real_a, imag_a,
                 real_c, imag_c, real_c_mod, imag_c_mod, abs_a, theta_a,
                 rms_res_r_2d, rms_res_i_2d, real_a_gapfill, imag_a_gapfill,
                 real_c_gapfill, imag_c_gapfill, real_c_mod_gapfill,
                 imag_c_mod_gapfill, abs_a_gapfill, theta_a_gapfill, flags,
                 out_dir):
    """
    start_ yr = initial year of period
    end_yr = final year of period
    season = 'summer', 'winter', 'jan', 'feb', 'mar', 'apr', 'may', 'jun'
             'jul, 'aug', 'sep', 'oct', 'nov' or 'dec'
    area = e.g. osi455_nh
    real_a = real component of the complex A parameter
    imag_a = imaginary component of the complex A parameter
    real_c = real component of the complex C parameter
    imag_c = imaginary component of the complex C parameter
    real_c_mod = normalized real component of the complex C parameter
    imag_c_mod = normalized imaginary component of the complex C parameter
    abs_a = absolute value of complex A parameter
    theta_a = angle value of complex A parameter
    rms_res_r_2d = yearly-averaged Root Mean Squared of real component of
    residuals
    rms_res_i_2d = yearly-averaged Root Mean Squared of imaginary component
    of residuals
    odir = output directory
    """

    # Area definition to write out
    grid_def_file = os.path.join(os.path.dirname(__file__), 'plots/grids_py.def')
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

    fname = 'inv_params_{}_{}-{}_{}_{}day.nc'.format(area, start_yr, end_yr, period, duration)
    fname = os.path.join(out_dir,fname)

    #adef = get_areadef(hemis)
    nx = area_def.width
    ny = area_def.height
    lons, lats = area_def.get_lonlats()

    with Dataset(fname, 'w', format='NETCDF3_CLASSIC') as dataset:

        # Dimensions
        dimx = dataset.createDimension('xc', nx)
        dimy = dataset.createDimension('yc', ny)
#        dimt = dataset.createDimension('time', 1)
#        dimnv = dataset.createDimension('nv', 2)

        # Auxiliary (coordinate) variables
        # The "if" clause deals with "+no_defs" etc
        pdict = dict([el.split('=') for el in area_def.proj4_string.split()
                      if (len(el.split('=')) == 2)])
        if re.search('ease', area_def.area_id):
            crsname = 'Lambert_Azimuthal_Grid'
            gridmapname = 'lambert_azimuthal_grid'
        else:
            crsname = 'Polar_Stereographic_Grid'
            gridmapname = 'polar_stereographic'
        crs = dataset.createVariable(crsname, np.int32)
        crs.grid_mapping_name = gridmapname
        crs.straight_vertical_longitude_from_pole = np.float32(pdict['+lon_0'])
        crs.latitude_of_projection_origin = np.float32(pdict['+lat_0'])
        if '+lat_ts' in pdict.keys():
            crs.standard_parallel = np.float32(pdict['+lat_ts'])
        crs.false_easting = 0.0
        crs.false_northing = 0.0
        if '+a' in pdict.keys():
            crs.semi_major_axis = np.float32(pdict['+a'])
        if '+b' in pdict.keys():
            crs.semi_minor_axis = np.float32(pdict['+b'])
            crs_type = 'b'
        elif '+rf' in pdict.keys():
            crs.inverse_flattening = np.float32(pdict['+rf'])
            crs_type = 'rf'
        if '+datum' in pdict.keys():
            crs.reference_ellipsoid_name = pdict['+datum']
        crs.proj4_string = str(area_def.proj4_string)

        xc = dataset.createVariable('xc', np.float64, ('xc',))
        xc.axis = "X"
        xc.units = "km"
        xc.long_name = "x coordinate of projection (eastings)"
        xc.standard_name = "projection_x_coordinate"
        xc[:] = area_def.projection_x_coords * 0.001

        yc = dataset.createVariable('yc', np.float64, ('yc',))
        yc.axis = "Y"
        yc.units = "km"
        yc.long_name = "y coordinate of projection (northings)"
        yc.standard_name = "projection_y_coordinate"
        yc[:] = area_def.projection_y_coords * 0.001

        lat = dataset.createVariable('lat', np.float32, ('yc', 'xc'))
        lat.long_name = "latitude coordinate"
        lat.standard_name = "latitude"
        lat.units = "degrees_north"
        lat[:] = lats

        lon = dataset.createVariable('lon', np.float32, ('yc', 'xc'))
        lon.long_name = "longitude coordinate"
        lon.standard_name = "longitude"
        lon.units = "degrees_east"
        lon[:] = lons

        # data
        real_a_var = dataset.createVariable('A_real', np.float32,
                                            ('yc', 'xc'),)
        real_a_var.long_name = "real component of the parameter A"
        real_a_var.standard_name = "A real"
        real_a_var.units = " "
        real_a_var.grid_mapping = crsname
        real_a_var.coordinates = "lat lon"
        real_a_var[:] = real_a

        imag_a_var = dataset.createVariable('A_imag', np.float32,
                                            ('yc', 'xc'),)
        imag_a_var.long_name = "imaginary component of the parameter A"
        imag_a_var.standard_name = "A imag"
        imag_a_var.units = " "
        imag_a_var.grid_mapping = crsname
        imag_a_var.coordinates = "lat lon"
        imag_a_var[:] = imag_a

        real_c_var = dataset.createVariable('C_real', np.float32,
                                            ('yc', 'xc'),)
        real_c_var.long_name = "real component of the parameter C"
        real_c_var.standard_name = "C real"
        real_c_var.units = " "
        real_c_var.grid_mapping = crsname
        real_c_var.coordinates = "lat lon"
        real_c_var[:] = real_c

        imag_c_var = dataset.createVariable('C_imag', np.float32,
                                            ('yc', 'xc'),)
        imag_c_var.long_name = "imaginary component of the parameter C"
        imag_c_var.standard_name = "C imag"
        imag_c_var.units = " "
        imag_c_var.grid_mapping = crsname
        imag_c_var.coordinates = "lat lon"
        imag_c_var[:] = imag_c

        abs_a_var = dataset.createVariable('absA', np.float32,
                                           ('yc', 'xc'),)
        abs_a_var.long_name = "absolute value of the parameter A"
        abs_a_var.standard_name = "abs A"
        abs_a_var.units = " "
        abs_a_var.grid_mapping = crsname
        abs_a_var.coordinates = "lat lon"
        abs_a_var[:] = abs_a

        theta_a_var = dataset.createVariable('ThetaA', np.float32,
                                             ('yc', 'xc'),)
        theta_a_var.long_name = "angle value of the parameter A"
        theta_a_var.standard_name = "Theta A"
        theta_a_var.units = "degrees"
        theta_a_var.grid_mapping = crsname
        theta_a_var.coordinates = "lat lon"
        theta_a_var[:] = theta_a

        real_c_mod_var = dataset.createVariable('c_real', np.float32,
                                                ('yc', 'xc'),)
        real_c_mod_var.long_name = "normalized real component of the parameter C"
        real_c_mod_var.standard_name = "c real"
        real_c_mod_var.units = " "
        real_c_mod_var.grid_mapping = crsname
        real_c_mod_var.coordinates = "lat lon"
        real_c_mod_var[:] = real_c_mod

        imag_c_mod_var = dataset.createVariable('c_imag', np.float32,
                                                ('yc', 'xc'),)
        imag_c_mod_var.long_name = "normalized imaginary component of the parameter C"
        imag_c_mod_var.standard_name = "c imag"
        imag_c_mod_var.units = " "
        imag_c_mod_var.grid_mapping = crsname
        imag_c_mod_var.coordinates = "lat lon"
        imag_c_mod_var[:] = imag_c_mod

        rms_real_var = dataset.createVariable('RMS_Res_real', np.float32,
                                              ('yc','xc'),)
        rms_real_var.long_name = "Root mean square of the real component of residuals"
        rms_real_var.standard_name = "RMS Res real"
        rms_real_var.units = " "
        rms_real_var.grid_mapping = crsname
        rms_real_var.coordinates = "lat lon"
        rms_real_var[:] = rms_res_r_2d

        rms_imag_var = dataset.createVariable('RMS_Res_imag', np.float32,
                                              ('yc','xc'),)
        rms_imag_var.long_name = "Root mean square of the imaginary component of residuals"
        rms_imag_var.standard_name = "RMS Res imag"
        rms_imag_var.units = " "
        rms_imag_var.grid_mapping = crsname
        rms_imag_var.coordinates = "lat lon"
        rms_imag_var[:] = rms_res_i_2d

        real_a_gf_var = dataset.createVariable('A_real_gapfill', np.float32,
                                               ('yc', 'xc'),)
        real_a_gf_var.long_name = "real component of the parameter A, gapfilled"
        real_a_gf_var.standard_name = "A real gapfill"
        real_a_gf_var.units = " "
        real_a_gf_var.grid_mapping = crsname
        real_a_gf_var.coordinates = "lat lon"
        real_a_gf_var[:] = real_a_gapfill

        imag_a_gf_var = dataset.createVariable('A_imag_gapfill', np.float32,
                                               ('yc', 'xc'),)
        imag_a_gf_var.long_name = "imaginary component of the parameter A, gapfilled"
        imag_a_gf_var.standard_name = "A imag gapfill"
        imag_a_gf_var.units = " "
        imag_a_gf_var.grid_mapping = crsname
        imag_a_gf_var.coordinates = "lat lon"
        imag_a_gf_var[:] = imag_a_gapfill

        real_c_gf_var = dataset.createVariable('C_real_gapfill', np.float32, ('yc', 'xc'),)
        real_c_gf_var.long_name = "real component of the parameter C, gapfilled"
        real_c_gf_var.standard_name = "C real gapfill"
        real_c_gf_var.units = " "
        real_c_gf_var.grid_mapping = crsname
        real_c_gf_var.coordinates = "lat lon"
        real_c_gf_var[:] = real_c_gapfill

        imag_c_gf_var = dataset.createVariable('C_imag_gapfill', np.float32,
                                               ('yc', 'xc'),)
        imag_c_gf_var.long_name = "imaginary component of the parameter C, gapfilled"
        imag_c_gf_var.standard_name = "C imag gapfill"
        imag_c_gf_var.units = " "
        imag_c_gf_var.grid_mapping = crsname
        imag_c_gf_var.coordinates = "lat lon"
        imag_c_gf_var[:] = imag_c_gapfill

        abs_a_gf_var = dataset.createVariable('absA_gapfill', np.float32, ('yc', 'xc'),)
        abs_a_gf_var.long_name = "absolute value of the parameter A, gapfilled"
        abs_a_gf_var.standard_name = "abs A gapfill"
        abs_a_gf_var.units = " "
        abs_a_gf_var.grid_mapping = crsname
        abs_a_gf_var.coordinates = "lat lon"
        abs_a_gf_var[:] = abs_a_gapfill

        theta_a_gf_var = dataset.createVariable('ThetaA_gapfill', np.float32,
                                                ('yc', 'xc'),)
        theta_a_gf_var.long_name = "angle value of the parameter A, gapfilled"
        theta_a_gf_var.standard_name = "Theta A gapfill"
        theta_a_gf_var.units = "degrees"
        theta_a_gf_var.grid_mapping = crsname
        theta_a_gf_var.coordinates = "lat lon"
        theta_a_gf_var[:] = theta_a_gapfill

        real_c_mod_gf_var = dataset.createVariable('c_real_gapfill',
                                                   np.float32, ('yc', 'xc'),)
        real_c_mod_gf_var.long_name = "normalized real component of the parameter C, gapfilled"
        real_c_mod_gf_var.standard_name = "c real gapfill"
        real_c_mod_gf_var.units = " "
        real_c_mod_gf_var.grid_mapping = crsname
        real_c_mod_gf_var.coordinates = "lat lon"
        real_c_mod_gf_var[:] = real_c_mod_gapfill

        imag_c_mod_gf_var = dataset.createVariable('c_imag_gapfill',
                                                   np.float32, ('yc', 'xc'),)
        imag_c_mod_gf_var.long_name = "normalized imaginary component of the parameter C, gapfilled"
        imag_c_mod_gf_var.standard_name = "c imag gapfill"
        imag_c_mod_gf_var.units = " "
        imag_c_mod_gf_var.grid_mapping = crsname
        imag_c_mod_gf_var.coordinates = "lat lon"
        imag_c_mod_gf_var[:] = imag_c_mod_gapfill

        flag_var = dataset.createVariable('flags', np.int8, ('yc', 'xc'),)
        flag_var.long_name = "status flags "
        flag_var.standard_name = "status flags"
        flag_var.units = " "
        flag_var.grid_mapping = crsname
        flag_var.coordinates = "lat lon"
        flag_var[:] = flags


    return fname
