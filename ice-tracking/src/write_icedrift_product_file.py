"""Function to write a NetCDF file from the ice drift output, modified
from WriteLRSID.py by T. Lavergne."""

import os
import re
import numpy as np
from netCDF4 import Dataset, date2num

here = os.path.dirname(os.path.abspath(__file__))
codedir = os.path.join(here, '../..')

sys.path.append(os.path.join(codedir, 'ice-tracking/src'))
from helper_fcns import fmtdate, round_sig

fill_values = {np.float32:-1e10, np.int16:-32767, np.int32:2147483647}

def write_icedrift(area_def,
                   inst, plat, channels, start_date, end_date,
                   out_area, out_projstr, out_dir,
                   driftx, drifty, t0, t1, flag,
                   sigdx, sigdy, corrdxdy,
                   uflag, length, direction,
                   outlonb, outlatb, outlone, outlate,
                   outfc, snavg, avgx, avgy,
                   length_avg, length_diff,
                   stdx, stdy, patternIndex, method=None, suffix=None):

    """
    area_def      Output area definition
    inst          Instrument name
    plat          Platform
    channels      Full list of channels, including "_Lap" if present
    start_date    Datetime object for start date
    end_date      Datetime object for end date
    out_area      Code for area, e.g. 'nh625'
    out_projstr   Projection of output grid
    out_dir       Directory for output
    driftx        x-component of drift
    drifty        y-component of drift
    t0            Start time
    t1            End time
    flag          Quality/Processing Flag
    sigdx         Standard deviation on driftx
    sigdy         Standard deviation on drifty
    corrdxdy      Correlation of errors in driftx and drifty
    uflag         Uncertainty flag
    length        Drift length
    direction     Direction of the drift
    outlonb       Longitude of drift point at start
    outlatb       Latitude of drift point at start
    outlone       Longitude of drift point at end
    outlate       Latitude of drift point at end
    outfc         Cross-correlation at best estimate
    snavg         Number of neighbour vectors to compute average
    avgx          Local average drift in x-direction
    avgy          Local average drift in y-direction
    length_avg    Length of local average drift
    length_diff   Length of the difference vector between the drift and average
    stdx          Local variability drift in x-direction
    stdy          Local variability drift in y-direction
    patternIndex  Size of the correlation pattern
    """

    start_date_fmt = fmtdate(start_date)
    end_date_fmt = fmtdate(end_date)

    # Concatenating channel list
    if channels is not None:
        channelsr = sorted(channels, reverse=True)
        channel_list = channelsr[0].replace('_', '-')
        for i in range(1, len(channels)):
            channel_list += '+' + channelsr[i].replace('_', '-')

    # Create output name
    if suffix is None:
        suffix = ''
    else:
        suffix = '_{}'.format(suffix)
    if method == 'merge':
        fname = ('icedrift_multi-oi_simplex_lev3_{}_{}w_{}w{}.nc'.format(
            out_area, start_date_fmt, end_date_fmt, suffix))
    elif method == 'wind':
        fname = ('icedrift_ecmwf-era5_none_wind_lev2_{}_{}w_{}w{}.nc'.format(
            out_area, start_date_fmt, end_date_fmt, suffix))
    else:
        fname = ('icedrift_{}-{}_{}_simplex_lev2_{}_{}w_{}w{}.nc'.format(
            inst, plat, channel_list, out_area, start_date_fmt, end_date_fmt,
            suffix))
    fnameandpath = os.path.join(out_dir, fname)

    # Checking size of the area definition
#    if area_def.pixel_size_x != 62500 or area_def.pixel_size_y != 62500:
#        raise ValueError("The grid parameters do not match a grid spacing "
#                         "of 62.5km! (x:{} y:{})".format(
#                         area_def.pixel_size_x, area_def.pixel_size_y))

    # Fetching parameters from the area definition
    ny, nx = area_def.shape
    lons, lats = area_def.get_lonlats()
    (out_llx, out_lly, out_urx, out_ury) = area_def.area_extent
    unitsstr = [x for x in area_def.proj4_string.split()
                if '+units=' in x][0][7:]
    if unitsstr == 'm':
        unitconv = 0.001
    else:
        unitconv = 1.0
    out_Ax = round_sig(area_def.resolution[0]) * unitconv
    out_Ay = round_sig(area_def.resolution[1]) * unitconv
    out_Bx = out_llx * unitconv
    out_By = out_ury * unitconv

    with Dataset(fnameandpath, 'w', format='NETCDF3_CLASSIC') as dataset:

        # dimension and auxiliary datasets
        dimx = dataset.createDimension('xc', nx)
        dimy = dataset.createDimension('yc', ny)
        dimt = dataset.createDimension('time', 1)
        dimnv = dataset.createDimension('nv', 2)

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

        tunits = "seconds since 1970-01-01 00:00:00"
        time = dataset.createVariable('time', np.float64, ('time',))
        time.axis = "T"
        time.long_name = "reference time of product"
        time.standard_name = "time"
        time.units = tunits
        time.calendar = "standard"
        time.bounds = "time_bnds"
        time[:] = date2num(end_date, tunits)

        tbounds = dataset.createVariable('time_bnds',
                                         np.float64, ('time', 'nv'))
        tbounds.units = tunits
        tbounds[:] = [date2num(start_date, tunits), date2num(end_date, tunits)]

        xc = dataset.createVariable('xc', np.float64, ('xc'))
        xc.axis = "X"
        xc.units = "km"
        xc.long_name = "x coordinate of projection (eastings)"
        xc.standard_name = "projection_x_coordinate"
        xc[:] = area_def.projection_x_coords * unitconv

        yc = dataset.createVariable('yc', np.float64, ('yc'))
        yc.axis = "Y"
        yc.units = "km"
        yc.long_name = "y coordinate of projection (northings)"
        yc.standard_name = "projection_y_coordinate"
        yc[:] = area_def.projection_y_coords * unitconv

        # Data

        var_lat = dataset.createVariable('lat', np.float32, ('yc', 'xc'))
        var_lat.long_name = "latitude coordinate"
        var_lat.standard_name = "latitude"
        var_lat.units = "degrees_north"
        var_lat[:] = lats

        var_lon = dataset.createVariable('lon', np.float32, ('yc', 'xc'))
        var_lon.long_name = "longitude coordinate"
        var_lon.standard_name = "longitude"
        var_lon.units = "degrees_east"
        var_lon[:] = lons

        var_driftx = dataset.createVariable('driftX', np.float32,
                                            ('time', 'yc', 'xc'),
                                            fill_value=fill_values[np.float32])
        var_driftx.long_name = "drift along X axis of grid"
        var_driftx.standard_name = "drift_x"
        var_driftx.units = "km"
        var_driftx.grid_mapping = crsname
        var_driftx.coordinates = "lat lon"
        var_driftx[:] = driftx

        var_drifty = dataset.createVariable('driftY', np.float32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float32])
        var_drifty.long_name = "drift along Y axis of grid"
        var_drifty.standard_name = "drift_x"
        var_drifty.units = "km"
        var_drifty.grid_mapping = crsname
        var_drifty.coordinates = "lat lon"
        var_drifty[:] = drifty

        var_t0 = dataset.createVariable('t0', np.int32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.int32])
        var_t0.long_name = "time of start pixel"
        var_t0.standard_name = "time_start"
        var_t0.units = "seconds since {:%Y-%m-%dT%H:%M:%SZ}".format(start_date)
        var_t0.grid_mapping = crsname
        var_t0.coordinates = "lat lon"
        var_t0[:] = t0

        var_t1 = dataset.createVariable('t1', np.int32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.int32])
        var_t1.long_name = "time of end pixel"
        var_t1.standard_name = "time_end"
        var_t1.units = "seconds since {:%Y-%m-%dT%H:%M:%SZ}".format(end_date)
        var_t1.grid_mapping = crsname
        var_t1.coordinates = "lat lon"
        var_t1[:] = t1

        var_flag = dataset.createVariable('flag', np.int16,
                                          ('time', 'yc', 'xc'),
                                          fill_value=fill_values[np.int16])
        var_flag.long_name = "Quality/Processing Flag"
        var_flag.standard_name = "qual_flag"
        var_flag.units = "1"
        var_flag.grid_mapping = crsname
        var_flag.coordinates = "lat lon"
        var_flag[:] = flag

        var_sigdx = dataset.createVariable('sigdX', np.float32,
                                           ('time', 'yc', 'xc'),
                                           fill_value=fill_values[np.float32])
        var_sigdx.long_name = "standard deviation on driftX"
        var_sigdx.standard_name = "stddevX"
        var_sigdx.units = "km"
        var_sigdx.grid_mapping = crsname
        var_sigdx.coordinates = "lat lon"
        var_sigdx[:] = sigdx

        var_sigdy = dataset.createVariable('sigdY', np.float32,
                                           ('time', 'yc', 'xc'),
                                           fill_value=fill_values[np.float32])
        var_sigdy.long_name = "standard deviation on driftY"
        var_sigdy.standard_name = "stddevY"
        var_sigdy.units = "km"
        var_sigdy.grid_mapping = crsname
        var_sigdy.coordinates = "lat lon"
        var_sigdy[:] = sigdy

        var_corrdxdy = dataset.createVariable('corrdXdY', np.float32,
                                              ('time', 'yc', 'xc'),
                                              fill_value=fill_values[np.float32])
        var_corrdxdy.long_name = "correlation of errors in driftX and driftY"
        var_corrdxdy.standard_name = "corrdXdY"
        var_corrdxdy.units = "km"
        var_corrdxdy.grid_mapping = crsname
        var_corrdxdy.coordinates = "lat lon"
        var_corrdxdy[:] = corrdxdy

        var_uflag = dataset.createVariable('uflag', np.int16,
                                           ('time', 'yc', 'xc'),
                                           fill_value=fill_values[np.int16])
        var_uflag.long_name = "Uncertainty Flag"
        var_uflag.standard_name = "uncert_flag"
        var_uflag.units = "1"
        var_uflag.grid_mapping = crsname
        var_uflag.coordinates = "lat lon"
        var_uflag[:] = uflag

        var_length = dataset.createVariable('length', np.float32,
                                            ('time', 'yc', 'xc'),
                                            fill_value=fill_values[np.float32])
        var_length.long_name = "length of the drift"
        var_length.standard_name = "length"
        var_length.units = "km"
        var_length.grid_mapping = crsname
        var_length.coordinates = "lat lon"
        var_length[:] = length

        var_direction = dataset.createVariable('dir', np.float32,
                                               ('time', 'yc', 'xc'),
                                               fill_value=fill_values[np.float32])
        var_direction.long_name = "direction of the drift"
        var_direction.standard_name = "direction"
        var_direction.units = "degrees from North"
        var_direction.grid_mapping = crsname
        var_direction.coordinates = "lat lon"
        var_direction[:] = direction

        var_outlatb = dataset.createVariable('latStart', np.float32,
                                             ('time', 'yc', 'xc'),
                                             fill_value=fill_values[np.float32])
        var_outlatb.long_name = "latitude of start of drift point"
        var_outlatb.standard_name = "lat_start"
        var_outlatb.units = "degrees"
        var_outlatb.grid_mapping = crsname
        var_outlatb.coordinates = "lat lon"
        var_outlatb[:] = outlatb

        var_outlonb = dataset.createVariable('lonStart', np.float32,
                                             ('time', 'yc', 'xc'),
                                             fill_value=fill_values[np.float32])
        var_outlonb.long_name = "longitude of start of drift point"
        var_outlonb.standard_name = "lon_start"
        var_outlonb.units = "degrees"
        var_outlonb.grid_mapping = crsname
        var_outlonb.coordinates = "lat lon"
        var_outlonb[:] = outlonb

        var_outlate = dataset.createVariable('latEnd', np.float32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float32])
        var_outlate.long_name = "latitude of end of drift point"
        var_outlate.standard_name = "lat_end"
        var_outlate.units = "degrees"
        var_outlate.grid_mapping = crsname
        var_outlate.coordinates = "lat lon"
        var_outlate[:] = outlate

        var_outlone = dataset.createVariable('lonEnd', np.float32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float32])
        var_outlone.long_name = "longitude of end of drift point"
        var_outlone.standard_name = "lon_end"
        var_outlone.units = "degrees"
        var_outlone.grid_mapping = crsname
        var_outlone.coordinates = "lat lon"
        var_outlone[:] = outlone

        var_corr = dataset.createVariable('corr', np.float32,
                                          ('time', 'yc', 'xc'),
                                          fill_value=fill_values[np.float32])
        var_corr.long_name = "cross-correlation at best estimate"
        var_corr.standard_name = "cross-correlation"
        var_corr.units = "1"
        var_corr.grid_mapping = crsname
        var_corr.coordinates = "lat lon"
        var_corr[:] = outfc

        if snavg is not None:
            var_navg = dataset.createVariable('navg', np.int16,
                                              ('time', 'yc', 'xc'),
                                              fill_value=fill_values[np.int16])
            var_navg.long_name = "number of neighbour vectors to compute average"
            var_navg.standard_name = "n_avg_drift"
            var_navg.units = "1"
            var_navg.grid_mapping = crsname
            var_navg.coordinates = "lat lon"
            var_navg[:] = snavg

        if avgx is not None:
            var_avgx = dataset.createVariable('avgX', np.float32,
                                              ('time', 'yc', 'xc'),
                                              fill_value=fill_values[np.float32])
            var_avgx.long_name = "local average drift along X axis of grid"
            var_avgx.standard_name = "avg_drift_x"
            var_avgx.units = "km"
            var_avgx.grid_mapping = crsname
            var_avgx.coordinates = "lat lon"
            var_avgx[:] = avgx

        if avgy is not None:
            var_avgy = dataset.createVariable('avgY', np.float32,
                                              ('time', 'yc', 'xc'),
                                              fill_value=fill_values[np.float32])
            var_avgy.long_name = "local average drift along Y axis of grid"
            var_avgy.standard_name = "avg_drift_y"
            var_avgy.units = "km"
            var_avgy.grid_mapping = crsname
            var_avgy.coordinates = "lat lon"
            var_avgy[:] = avgy

        if length_avg is not None:
            var_avglen = dataset.createVariable('avglen', np.float32,
                                                ('time', 'yc', 'xc'),
                                                fill_value=fill_values[np.float32])
            var_avglen.long_name = "length of local average drift"
            var_avglen.standard_name = "avg_len"
            var_avglen.units = "km"
            var_avglen.grid_mapping = crsname
            var_avglen.coordinates = "lat lon"
            var_avglen[:] = length_avg

        if length_diff is not None:
            var_difflen = dataset.createVariable('difflen', np.float32,
                                                 ('time', 'yc', 'xc'),
                                                 fill_value=fill_values[np.float32])
            var_difflen.long_name = "length of difference vector AVG-DRIFT"
            var_difflen.standard_name = "diff_len"
            var_difflen.units = "km"
            var_difflen.grid_mapping = crsname
            var_difflen.coordinates = "lat lon"
            var_difflen[:] = length_diff

        if stdx is not None:
            var_stdx = dataset.createVariable('stdX', np.float32,
                                              ('time', 'yc', 'xc'),
                                              fill_value=fill_values[np.float32])
            var_stdx.long_name = "local variability drift along X axis of grid"
            var_stdx.standard_name = "std_drift_x"
            var_stdx.units = "km"
            var_stdx.grid_mapping = crsname
            var_stdx.coordinates = "lat lon"
            var_stdx[:] = stdx

        if stdy is not None:
            var_stdy = dataset.createVariable('stdY', np.float32,
                                              ('time', 'yc', 'xc'),
                                              fill_value=fill_values[np.float32])
            var_stdy.long_name = "local variability drift along Y axis of grid"
            var_stdy.standard_name = "std_drift_y"
            var_stdy.units = "km"
            var_stdy.grid_mapping = crsname
            var_stdy.coordinates = "lat lon"
            var_stdy[:] = stdy

        var_pattern = dataset.createVariable('pattern', np.int16,
                                             ('time', 'yc', 'xc'),
                                             fill_value=fill_values[np.int16])
        var_pattern.long_name = "Size of the correlation pattern"
        var_pattern.standard_name = "pattern_index"
        var_pattern.units = "1"
        var_pattern.grid_mapping = crsname
        var_pattern.coordinates = "lat lon"
        var_pattern[:] = patternIndex

        # Global Attributes
        dataset.filename = fname
        dataset.Conventions = "CF-1.0"
        if method == 'merge':
            dataset.instrument = 'multi-oi'
            dataset.platform = 'multi-oi'
        elif method == 'wind':
            dataset.instrument = 'era5'
            dataset.platform = 'wind'
        else:
            dataset.instrument = inst
            dataset.platform = plat
        dataset.processing_software = "icedrift_solve_core"
        dataset.start_date_and_time = "{:%Y-%m-%dT%H:%M:%SZ}".format(start_date)
        dataset.end_date_and_time = "{:%Y-%m-%dT%H:%M:%SZ}".format(end_date)

        # It is assumed for now that the grid file is in the local directory
        dataset.remapping_gridfile = os.path.join(os.getcwd(), 'grids_py.def')
        dataset.remapping_gridname = out_area

        # The projstr is in km not m, so this is filtered to use only
        # required projstr params
        projstr_list = out_projstr.split()
        c_out_projstr = ''
        if re.search('ease', area_def.area_id):
            req_proj_items = ['proj', 'datum', 'lat_0', 'lon_0', 'x_0', 'y_0']
        else:
            if crs_type == 'b':
                req_proj_items = ['proj', 'a', 'b', 'lat_0', 'lat_ts', 'lon_0']
            elif crs_type == 'rf':
                req_proj_items = ['proj', 'a', 'rf', 'lat_0', 'lat_ts', 'lon_0']
        for proj_item in req_proj_items:
            c_out_projstr += ([x for x in projstr_list if
                               ('+' + proj_item + '=') in x][0]).strip() + ' '
        dataset.remapping_projstr = c_out_projstr.strip()

        # The gridparams written out here are C-standard ones (so pixel
        # centres, not pixel corners)
        c_rmp_gdpms = np.asarray([out_Ax, out_By - (0.5 * out_Ay), out_Ay,
                                  out_Bx + (0.5 * out_Ax)], dtype=np.float32)
        c_remapping_gridparams = (c_rmp_gdpms[0], c_rmp_gdpms[1],
                                  c_rmp_gdpms[2], c_rmp_gdpms[3])
        dataset.remapping_gridparams = c_remapping_gridparams

    print('Writing out {}'.format(fnameandpath))

    return fnameandpath
