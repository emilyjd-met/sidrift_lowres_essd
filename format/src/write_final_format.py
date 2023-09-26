'''Final formatting to create CF compliant netCDF4 file for distribution.'''

import os
import sys
import re
import subprocess
import argparse
from argparse import RawDescriptionHelpFormatter
from datetime import datetime
import numpy as np
from netCDF4 import Dataset, date2num
try:
    import simplejson as json
except ImportError:
    import json
from collections import OrderedDict

sys.path.append(os.path.join(os.path.dirname(__file__), '../../merge/src'))
from read_icedrift_product_file import read_drift
from icedrift_flags import id_flags
from icedrift_final_flags import id_final_flags, flag_matrix
sys.path.append(os.path.join(os.path.dirname(__file__), '../../ice-tracking/src'))

def get_tracking_id():
    cmd  = 'uuidgen -r'
    uuid = subprocess.check_output(cmd.split()).decode('utf-8').replace('\n','')
    return uuid


fill_values = {np.float32:-1e10, np.int16:-32767, np.int32:2147483647,
               np.int8:-127, np.float64:-1e10}

sensor_n = {'ssmi':'SSM/I',
            'ssmis':'SSMIS',
            'amsr':'AMSR-E',
            'amsr2':'AMSR2',
            'multi-oi': 'SSM/I,SSMIS,AMSR-E,AMSR2'}
platform_n = {'ssmi-f10':'DMSP-F10',
              'ssmi-f11':'DMSP-F11',
              'ssmi-f13':'DMSP-F13',
              'ssmi-f14':'DMSP-F14',
              'ssmi-f15':'DMSP-F15',
              'ssmis-f16':'DMSP-F16',
              'ssmis-f17':'DMSP-F17',
              'ssmis-f18':'DMSP-F18',
              'amsr-aq':'Aqua',
              'amsr2-gw1':'GCOM-W1',
              'multi-oi': 'DMSP-F<08,10,11,13,14,15>,DMSP-F<16,17,18>,Aqua,GCOM-W1'}

source_n = {'ssmi': "FCDR of SMMR / SSMI / SSMIS Brightness Temperatures R4 (beta) from EUMETSAT Climate Monitoring SAF",
            'ssmis': "FCDR of SMMR / SSMI / SSMIS Brightness Temperatures R4 (beta) from EUMETSAT Climate Monitoring SAF",
            'amsr': "FCDR of AMSR-E Brightness Temperatures V3 (doi: 10.5067/AMSR-E/AE_L2A.003)",
            'amsr2': "AMSR2 L1R data from JAXA",
            'wind': "C3S/ECMWF ERA5 data"}
#source_n['multi-oi'] = ',\n '.join([source_n['ssmi'], source_n['amsr'],
source_n['multi-oi'] = ', '.join([source_n['ssmi'], source_n['amsr'],
                                    source_n['amsr2'], source_n['wind']])

hs_dict = {'nh':'Northern Hemisphere', 'sh':'Southern Hemisphere'}


def parse_args():

    p = argparse.ArgumentParser("write_final_format",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--infile', required=True,
                   help="Input file to convert from proc to final format")
    p.add_argument('-o', '--outdir', required=True,
                   help="Output directory")
    p.add_argument('-w', '--verinfo', required=False, default='cdr-v1p0',
                   help="Dataset and version information")
    p.add_argument('-r', '--reproc', action='store_true', default=False,
                   help="Write output files to YYYY/MM directories")
    p.add_argument('-f', '--force', action='store_true', default=False,
                   help="Force overwrite if the formatted file already exists")

    args = p.parse_args()

    return args


def convert_flags(proc_flags):
    '''Convert proc format flags to final format flags'''

    fflags = np.empty_like(proc_flags)
    for flg in flag_matrix:
        fflags[proc_flags == id_flags[flg]] = id_final_flags[flag_matrix[flg]][0]

    return fflags


def write_final_format(infile, outdir, verinfo='cdr-v1p0', reproc=False,
                       force=False):
    '''Code to read an ice drift file in proc format and output to final
    format'''

    # Read the proc format file
    multi = False
    if 'multi-oi' in infile:
        multi = True
    driftdata = read_drift(infile, short=True, multi=multi)

    # Create the output name
    inst = driftdata['globals']['instrument']
    if inst == 'multi-oi':
        platinst = inst
    elif inst == 'wind':
        platinst = inst
    else:
        platinst = '{}-{}'.format(inst, driftdata['globals']['platform'])

    gridname = driftdata['globals']['remapping_gridname'].replace('h-', 'h_')
    hemi = gridname[:2]
    sdatestr = driftdata['globals']['start_date_and_time']
    sdate = datetime.strptime(sdatestr, '%Y-%m-%dT%H:%M:%SZ')
    sdatefmtlng = datetime.strftime(sdate, '%Y%m%d%H%M')
    edatestr = driftdata['globals']['end_date_and_time']
    edate = datetime.strptime(edatestr, '%Y-%m-%dT%H:%M:%SZ')
    edatefmtlng = datetime.strftime(edate, '%Y%m%d%H%M')
    dur = (edate - sdate).days
    if dur == 1:
        durstr = '24h'
    elif dur == 2:
        durstr = '48h'
    else:
        raise ValueError("Unrecognised duration in days: {}".format(dur))

    if platinst == 'multi-oi':
        outbase = 'ice_drift_{}_{}_{}-{}.nc'.format(gridname, verinfo,
                                                    durstr, edatefmtlng)
    else:
        outbase = 'ice_drift_{}_{}-{}_{}-{}.nc'.format(gridname, verinfo,
                                                       platinst, durstr,
                                                       edatefmtlng)
    if reproc:
        outdir = os.path.join(outdir, datetime.strftime(edate, '%Y'),
                              datetime.strftime(edate, '%m'))
    outname = os.path.join(outdir, outbase)
    if os.path.isfile(outname) and not force:
        print("{} already exists and force is not set. Exiting.".format(outname))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Convert the proc format flags to final format
    fflags = convert_flags(driftdata['flag'])

    # The dY must have the sign reversed to match with the y-axis
    driftdata['driftY'] = -1.0 * driftdata['driftY']

    # Read in the json file
    dirname = os.path.dirname(__file__)
    if verinfo == 'cdr-v1p0':
        json_file = os.path.join(dirname, 'ice_drift_cdr_v1_metadata.json')
    else:
        raise ValueError("Unrecognised version without metadata available: {}".format(verinfo))
    if not os.path.isfile(json_file):
        raise FileNotFoundError('Unable to find {}\n'.format(json_file))
    with open(json_file, 'r') as jsf:
        json_dict = json.load(jsf, object_pairs_hook=OrderedDict)


    # Write out file
    with Dataset(outname, 'w', format='NETCDF4') as dataset:

        # Dimension and auxiliary datasets
        dimt = dataset.createDimension('time', None)
        dimx = dataset.createDimension('xc', len(driftdata['xc']))
        dimy = dataset.createDimension('yc', len(driftdata['yc']))
        dimnv = dataset.createDimension('nv', 2)

        # Retrieving the projection info and creating this variable
        # The "if" clause deals with "+no_defs" etc
        pdict = dict([el.split('=') for el in driftdata['proj4_string'].split()
                      if len(el.split('=')) == 2])
        if re.search('Lambert', driftdata['grid_mapping']):
            crsname = 'Lambert_Azimuthal_Equal_Area'
            gridmapname = 'lambert_azimuthal_equal_area'
        else:
            crsname = 'Polar_Stereographic_Grid'
            gridmapname = 'polar_stereographic'
        crs = dataset.createVariable(crsname, np.int32, zlib=True)
        crs.grid_mapping_name = gridmapname
        if crs.grid_mapping_name == 'lambert_azimuthal_equal_area':
            crs.longitude_of_projection_origin = np.float32(pdict['+lon_0'])
        else:
            crs.straight_vertical_longitude_from_pole = np.float32(pdict['+lon_0'])
        crs.latitude_of_projection_origin = np.float32(pdict['+lat_0'])
        if '+lat_ts' in pdict.keys():
            crs.standard_parallel = np.float32(pdict['+lat_ts'])
        crs.false_easting = 0.0
        crs.false_northing = 0.0
        if crsname == 'Lambert_Azimuthal_Equal_Area' and '+datum=WGS84' in driftdata['proj4_string'].split():
            crs.semi_major_axis = 6378137.0,
            crs.inverse_flattening = 298.257223563,
        if '+a' in pdict.keys():
            crs.semi_major_axis = np.float32(pdict['+a'])
        if '+b' in pdict.keys():
            crs.semi_minor_axis = np.float32(pdict['+b'])
            crs_type = 'b'
        elif '+rf' in pdict.keys():
            crs.inverse_flattening = np.float32(pdict['+rf'])
            crs_type = 'rf'
        if '+datum' in pdict.keys() and not crs.grid_mapping_name == 'lambert_azimuthal_equal_area':
            crs.reference_ellipsoid_name = pdict['+datum']
        crs.proj4_string = driftdata['proj4_string']

        # Dimension information

        tunits = "seconds since 1978-01-01 00:00:00"
        vtime = dataset.createVariable('time', np.float64, ('time'),
                                       zlib=True)
        vtime.axis = "T"
        vtime.long_name = "reference time of product"
        vtime.standard_name = "time"
        vtime.units = tunits
        vtime.calendar = "standard"
        vtime.bounds = "time_bnds"
        vtime.comment = "The time scalar dimension holds the date of _end_ of motion."
        vtime[:] = date2num(edate, tunits)
        tbounds = dataset.createVariable('time_bnds',
                                         np.float64, ('time', 'nv'),
                                         zlib=True)
        tbounds.units = tunits
        tbounds[0, :] = [date2num(sdate, tunits), date2num(edate, tunits)]

        xc = dataset.createVariable('xc', np.float64, ('xc'),
                                    zlib=True)
        xc.axis = "X"
        xc.units = "km"
        xc.long_name = "x coordinate of projection (eastings)"
        xc.standard_name = "projection_x_coordinate"
        xc[:] = driftdata['xc']

        yc = dataset.createVariable('yc', np.float64, ('yc'),
                                    zlib=True)
        yc.axis = "Y"
        yc.units = "km"
        yc.long_name = "y coordinate of projection (northings)"
        yc.standard_name = "projection_y_coordinate"
        yc[:] = driftdata['yc']

        # Variables

        var_lat = dataset.createVariable('lat', np.float32, ('yc', 'xc'),
                                         zlib=True)
        var_lat.long_name = "latitude coordinate"
        var_lat.standard_name = "latitude"
        var_lat.units = "degrees_north"
        var_lat[:] = driftdata['lat']

        var_lon = dataset.createVariable('lon', np.float32, ('yc', 'xc'),
                                         zlib=True)
        var_lon.long_name = "longitude coordinate"
        var_lon.standard_name = "longitude"
        var_lon.units = "degrees_east"
        var_lon[:] = driftdata['lon']

        ddt0 = np.float64(driftdata['t0'])
        # Converting to seconds since 1978
        ddt0 = ddt0 + date2num(sdate, tunits)
        ddt0[ddt0 == fill_values[np.int32]] = fill_values[np.float64]
        ddt0[ddt0.mask] = fill_values[np.float32]
        var_t0 = dataset.createVariable('t0', np.float64,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float64],
                                        zlib=True)
        var_t0.long_name = "time for start of displacement"
        var_t0.units = tunits
        var_t0.grid_mapping = crsname
        var_t0.coordinates = "lat lon"
        var_t0[0, :] = ddt0

        ddt1 = np.float64(driftdata['t1'])
        # Converting to seconds since 1978
        ddt1 = ddt1 + date2num(edate, tunits)
        ddt1[ddt1 == fill_values[np.int32]] = fill_values[np.float64]
        var_t1 = dataset.createVariable('t1', np.float64,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float64],
                                        zlib=True)
        var_t1.long_name = "time for end of displacement"
        var_t1.units = tunits
        var_t1.grid_mapping = crsname
        var_t1.coordinates = "lat lon"
        var_t1[0, :] = ddt1

        var_lat1 = dataset.createVariable('lat1', np.float32,
                                          ('time', 'yc', 'xc'),
                                          fill_value=fill_values[np.float32],
                                          zlib=True)
        var_lat1.long_name = "latitude at end of displacement"
        var_lat1.units = "degrees_north"
        var_lat1.grid_mapping = crsname
        var_lat1.coordinates = "lat lon"
        var_lat1[0, :] = driftdata['latEnd']

        var_lon1 = dataset.createVariable('lon1', np.float32,
                                          ('time', 'yc', 'xc'),
                                          fill_value=fill_values[np.float32],
                                          zlib=True)
        var_lon1.long_name = "longitude at end of displacement"
        var_lon1.units = "degrees_east"
        var_lon1.grid_mapping = crsname
        var_lon1.coordinates = "lat lon"
        var_lon1[0, :] = driftdata['lonEnd']

        var_dx = dataset.createVariable('dX', np.float32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float32],
                                        zlib=True)
        var_dx.long_name = "component of the displacement along the x axis of the grid"
        var_dx.standard_name = "sea_ice_x_displacement"
        var_dx.units = "km"
        var_dx.grid_mapping = crsname
        var_dx.coordinates = "lat lon"
        var_dx.ancillary_variables = "uncert_dX_and_dY status_flag"
        var_dx[0, :] = driftdata['driftX']

        var_dy = dataset.createVariable('dY', np.float32,
                                        ('time', 'yc', 'xc'),
                                        fill_value=fill_values[np.float32],
                                        zlib=True)
        var_dy.long_name = "component of the displacement along the y axis of the grid"
        var_dy.standard_name = "sea_ice_y_displacement"
        var_dy.units = "km"
        var_dy.grid_mapping = crsname
        var_dy.coordinates = "lat lon"
        var_dy.ancillary_variables = "uncert_dX_and_dY status_flag"
        var_dy[0, :] = driftdata['driftY']

        flag_vals = [np.int8(v[0]) for v in id_final_flags.values()]
        flag_meanings = ' '.join([v[1] for v in id_final_flags.values()])
        flag_desc = ['{} -> {}'.format(v[0], v[2])
                     for v in id_final_flags.values()]
        var_flag = dataset.createVariable('status_flag', np.int8,
                                          ('time', 'yc', 'xc'),
                                          zlib=True)
        var_flag.long_name = "rejection and quality level flag"
        var_flag.standard_name = "status_flag"
        var_flag.valid_min = min(flag_vals)
        var_flag.valid_max = max(flag_vals)
        var_flag.grid_mapping = crsname
        var_flag.coordinates = "lat lon"
        var_flag.flag_values = flag_vals
        var_flag.flag_meanings = flag_meanings
        var_flag.flag_descriptions = ', '.join(flag_desc)
        var_flag[0, :] = fflags

        var_uncert = dataset.createVariable('uncert_dX_and_dY', np.float32,
                                            ('time', 'yc', 'xc'),
                                            fill_value=fill_values[np.float32],
                                            zlib=True)
        var_uncert.long_name = "uncertainty (1 standard deviation) of dX and dY components of the drift vector"
        var_uncert.standard_name = "sea_ice_x_displacement standard_error"
        var_uncert.units = "km"
        var_uncert.coverage_content_type = "qualityInformation"
        var_uncert.grid_mapping = crsname
        var_uncert.coordinates = "lat lon"
        if multi:
            var_uncert[0, :] = driftdata['sigdX']
        else:
            var_uncert[0, :] = driftdata['fsX']

        # Global Attributes
        for attr in list(json_dict['attributes'].keys()):
            if attr == "keywords":
                setattr(dataset, attr,
                        json_dict['attributes'][attr].format(hs_dict[hemi]))
            elif attr == "geospatial_lat_max":
                norlat = np.max(driftdata['lat'])
                if norlat > 89.5:
                    norlat = 90.
                setattr(dataset, attr, norlat)
            elif attr == "geospatial_lat_min":
                soulat = np.min(driftdata['lat'])
                if soulat < -89.5:
                    soulat = -90.
                setattr(dataset, attr, soulat)
            elif attr == "geospatial_lon_max":
                easlon = np.max(driftdata['lon'])
                if easlon > 179.5:
                    easlon = 180.
                setattr(dataset, attr, easlon)
            elif attr == "geospatial_lon_min":
                weslon = np.min(driftdata['lon'])
                if weslon < -179.5:
                    weslon = -180.
                setattr(dataset, attr, weslon)
            elif attr == "sensor":
                setattr(dataset, attr, sensor_n[inst])
            elif attr == "platform":
                setattr(dataset, attr, platform_n[platinst])
            elif attr == "source":
                setattr(dataset, attr, source_n[inst])
            elif attr == 'time_coverage_start':
                setattr(dataset, attr,
                        json_dict['attributes'][attr].format(d=sdate))
            elif attr == 'time_coverage_end':
                setattr(dataset, attr,
                        json_dict['attributes'][attr].format(d=edate))
            elif attr in ['history', 'date_created']:
                setattr(dataset, attr,
                        json_dict['attributes'][attr].format(
                            d=datetime.utcnow()))
            elif attr == 'geospatial_bounds_crs':
                setattr(dataset, attr,
                        json_dict['attributes'][attr][gridname])
            elif ( attr == 'tracking_id' ):
                setattr(dataset, attr, get_tracking_id())
            else:
                setattr(dataset, attr, json_dict['attributes'][attr])


    print('Writing out {}'.format(outname))

    return outname


def main():

    args = parse_args()
    infile = args.infile
    outdir = args.outdir
    verinfo = args.verinfo
    reproc = args.reproc
    force = args.force

    outfile = write_final_format(infile, outdir, verinfo=verinfo,
                                 reproc=reproc, force=force)


if __name__ == '__main__':

    main()
