"""Function to concatenate icedrift parameter files to a per-hemisphere output"""

import os
import re
import subprocess
import argparse
from argparse import RawDescriptionHelpFormatter
from glob import glob
try:
    import simplejson as json
except ImportError:
    import json
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset, date2num
from collections import OrderedDict

fill_values = {np.float32:-1e10, np.int16:-32767, np.int32:2147483647}

patt = "inv_params_osi455_{}_{}-{}_{}_1day.nc"

monlist = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep',
           'oct', 'nov', 'dec']

gvardict = {
    'grid': 'Lambert_Azimuthal_Grid',
    'xc': 'xc',
    'yc': 'yc',
    'lat': 'lat',
    'lon': 'lon',
}

pvardict = {
#    'a_real': 'A_real',
#    'a_imag': 'A_imag',
#    'c_real': 'C_real',
#    'c_imag': 'C_imag',
#    'thetaa': 'ThetaA',
#    'absa': 'absA',
#    'c_real_2': 'c_real',
#    'c_imag_2': 'c_imag',
    'rms_real': 'RMS_Res_real',
    'rms_imag': 'RMS_Res_imag',
    'a_real_gf': 'A_real_gapfill',
    'a_imag_gf': 'A_imag_gapfill',
    'c_real_gf': 'C_real_gapfill',
    'c_imag_gf': 'C_imag_gapfill',
    'absa_gf': 'absA_gapfill',
    'thetaa_gf': 'ThetaA_gapfill',
#    'c_real_2_gf': 'c_real_gapfill',
#    'c_imag_2_gf': 'c_imag_gapfill',
    'flags': 'flags',
}

hs_dict = {'nh':'Northern Hemisphere', 'sh':'Southern Hemisphere'}
gridname = {'nh': 'nh_ease2-750', 'sh': 'sh_ease2-750'}

def get_tracking_id():
    cmd  = 'uuidgen -r'
    uuid = subprocess.check_output(cmd.split()).decode('utf-8').replace('\n','')
    return uuid


def parse_args():

    valid_period = ["summer", "winter", "jan", "feb", "mar", "apr", "may",
                    "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    valid_icefmt = ['proc', 'final']

    p = argparse.ArgumentParser('concatenate_paramfiles',
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--indir',
                   help='Input directory for parameter files')
    p.add_argument('-o', '--outdir',
                   help='Output directory for concatenated parameter files')

    args = p.parse_args()

    return args


def read_paramfile(paramfile):

    pfile = {}
    pfile['data'] = {}
    pfile['metadata'] = {}

    # Load the parameters
    with Dataset(paramfile, mode='r') as dataset:
        for key, value in gvardict.items():
            pfile['data'][key] = dataset.variables[value][:]
            pfile['metadata'][key] = dataset.variables[value].__dict__
        for key, value in pvardict.items():
            pfile['data'][key] = dataset.variables[value][:]
            pfile['metadata'][key] = dataset.variables[value].__dict__

    return pfile


def concatenate_pfiles(indir, outdir):
    """Concatenate the monthly parameter files with metadata added """

    # Read in the json file for global metadata
    dirname = os.path.dirname(__file__)
    json_file = os.path.join(dirname, '../../format/src',
                             'ice_drift_cdr_v1_metadata.json')
    if not os.path.isfile(json_file):
        raise FileNotFoundError('Unable to find {}\n'.format(json_file))
    with open(json_file, 'r') as jsf:
        json_dict = json.load(jsf, object_pairs_hook=OrderedDict)

    # Glob all the files from the input directory that match the pattern
    gpatt = patt.replace("{}", "*")
    globlist = glob(os.path.join(indir, gpatt))

    # One hemisphere at a time
    for hemi in ['nh', 'sh']:

        # Find the sublist of matching files for the hemisphere
        hgloblist = [x for x in globlist if hemi in x]

        # Cycle over the months in order
        for i, month in enumerate(monlist):

            # Assume there is only one per-month files for each hemisphere
            monthf = [x for x in hgloblist if month in x][0]
            monthdata = read_paramfile(monthf)

            # Initialise an array to hold the info
            if i == 0:
                data = {}
                metadata = {}

                data['startdates'] = []
                data['enddates'] = []

                # Copying the general variables
                for key in gvardict:
                    data[key] = monthdata['data'][key]
                    metadata[key] = monthdata['metadata'][key]

                # Grabbing a data array and defining the shape for the
                # overall data arrays
                firstkey = list(pvardict.keys())[0]
                pshape = monthdata['data'][firstkey].shape
                dshape = (12, pshape[0], pshape[1])

                # Initialising the data arrays and saving the metadata
                for key in pvardict:
                    data[key] = np.zeros(dshape,
                                         dtype=monthdata['data'][key].dtype)
                    metadata[key] = monthdata['metadata'][key]

            # Add the data to the arrays
            for key in pvardict:
                data[key][i, :] = monthdata['data'][key]

            # Check the dates in the filename to find the earliest and
            # latest dates - the last date is set to 23:59:59 on the last
            # day of the month
            rmatch = re.match("^(.*)(\d{4})-(\d{4})(.*)$",
                              os.path.basename(monthf))
            sdate = datetime(int(rmatch[2]), i + 1, 1, 0, 0, 0)
            if i == 11:
                edateplus = datetime(int(rmatch[3]) + 1, 1, 1, 0, 0, 0)
            else:
                edateplus = datetime(int(rmatch[3]), i + 2, 1, 0, 0, 0)
            edate = edateplus - timedelta(seconds=1)
            data['startdates'].append(sdate)
            data['enddates'].append(edate)

        # Create output name
        mindate = min(data['startdates'])
        maxdate = max(data['enddates'])
        mindatestr = datetime.strftime(mindate, '%Y%m')
        maxdatestr = datetime.strftime(maxdate, '%Y%m')
        fname = 'inv_params_osi455_{}_{}-{}_1day.nc'.format(hemi, mindatestr,
                                                            maxdatestr)
        fnameandpath = os.path.join(outdir, fname)

        # Write the file
        with Dataset(fnameandpath, 'w') as dataset:

            # dimension and auxiliary datasets
            dimx = dataset.createDimension('xc', dshape[1])
            dimy = dataset.createDimension('yc', dshape[2])
            dimm = dataset.createDimension('month', dshape[0])
            dimnv = dataset.createDimension('nv', 2)
            dimshape = ('month', 'yc', 'xc')

            newvar = {}

            # Retrieving the projection info and creating this variable
            # The "if" clause deals with "+no_defs" etc
            pdict = dict([el.split('=')
                          for el in metadata['grid']['proj4_string'].split()
                          if len(el.split('=')) == 2])
            if re.search('lambert', metadata['grid']['grid_mapping_name']):
                crsname = 'Lambert_Azimuthal_Equal_Area'
                gridmapname = 'lambert_azimuthal_equal_area'
            else:
                crsname = 'Polar_Stereographic_Grid'
                gridmapname = 'polar_stereographic'
            newvar['crs'] = dataset.createVariable(crsname, np.int32)
            newvar['crs'].grid_mapping_name = gridmapname
            if newvar['crs'].grid_mapping_name == 'lambert_azimuthal_equal_area':
                newvar['crs'].longitude_of_projection_origin = np.float32(pdict['+lon_0'])
            else:
                newvar['crs'].straight_vertical_longitude_from_pole = np.float32(pdict['+lon_0'])
            newvar['crs'].latitude_of_projection_origin = np.float32(pdict['+lat_0'])
            if '+lat_ts' in pdict.keys():
                newvar['crs'].standard_parallel = np.float32(pdict['+lat_ts'])
            newvar['crs'].false_easting = 0.0
            newvar['crs'].false_northing = 0.0
            if crsname == 'Lambert_Azimuthal_Equal_Area' and '+datum=WGS84' in metadata['grid']['proj4_string'].split():
                newvar['crs'].semi_major_axis = 6378137.0,
                newvar['crs'].inverse_flattening = 298.257223563,
            if '+a' in pdict.keys():
                newvar['crs'].semi_major_axis = np.float32(pdict['+a'])
            if '+b' in pdict.keys():
                newvar['crs'].semi_minor_axis = np.float32(pdict['+b'])
                crs_type = 'b'
            elif '+rf' in pdict.keys():
                newvar['crs'].inverse_flattening = np.float32(pdict['+rf'])
                crs_type = 'rf'
            if '+datum' in pdict.keys() and not newvar['crs'].grid_mapping_name == 'lambert_azimuthal_equal_area':
                newvar['crs'].reference_ellipsoid_name = pdict['+datum']

            # The projstr is in km not m, so this is filtered to use only
            # required projstr params
            out_projstr_list = []
            if crsname == 'Lambert_Azimuthal_Equal_Area':
                req_proj_items = ['+proj', '+datum', '+lat_0', '+lon_0',
                                  '+x_0', '+y_0']
            else:
                if crs_type == 'b':
                    req_proj_items = ['+proj', '+a', '+b', '+lat_0', '+lat_ts',
                                      '+lon_0']
                elif crs_type == 'rf':
                    req_proj_items = ['+proj', '+a', '+rf', '+lat_0',
                                      '+lat_ts', '+lon_0']
            for proj_item in req_proj_items:
                if proj_item in pdict.keys():
                    out_projstr_list.append('{}={}'.format(proj_item,
                                                           pdict[proj_item]))
            out_projstr = ' '.join(out_projstr_list)

            newvar['crs'].proj4_string = out_projstr

            newvar['month'] = dataset.createVariable('month', np.int16,
                                                     ('month'))
            newvar['month'].long_name = "Index of month"
            newvar['month'][:] = [x + x for x in range(dshape[0])]

            newvar['monthname'] = dataset.createVariable('monthname', str,
                                                         ('month'))
            newvar['monthname'].long_name = "Name of month"
            for i, mon in enumerate(monlist):
                newvar['monthname'][i] = mon

            # Creating the time bounds array
            tunits = "seconds since 1970-01-01 00:00:00"
            sd = [date2num(x, tunits) for x in data['startdates']]
            ed = [date2num(x, tunits) for x in data['enddates']]
            timebounds = np.array([sd, ed]).T
            newvar['tbounds'] = dataset.createVariable('time_bnds',
                                             np.float64, ('month', 'nv'))
            newvar['tbounds'].units = tunits
            newvar['tbounds'][:] = timebounds

            newvar['xc'] = dataset.createVariable('xc', np.float64, ('xc'))
            newvar['xc'].axis = "X"
            newvar['xc'].units = "km"
            newvar['xc'].long_name = "x coordinate of projection (eastings)"
            newvar['xc'].standard_name = "projection_x_coordinate"
            newvar['xc'][:] = data['xc']

            newvar['yc'] = dataset.createVariable('yc', np.float64, ('yc'))
            newvar['yc'].axis = "Y"
            newvar['yc'].units = "km"
            newvar['yc'].long_name = "y coordinate of projection (northings)"
            newvar['yc'].standard_name = "projection_y_coordinate"
            newvar['yc'][:] = data['yc']

            newvar['lat'] = dataset.createVariable('lat', np.float32,
                                                   ('yc', 'xc'))
            newvar['lat'].long_name = "latitude coordinate"
            newvar['lat'].standard_name = "latitude"
            newvar['lat'].units = "degrees_north"
            newvar['lat'][:] = data['lat']

            newvar['lon'] = dataset.createVariable('lon', np.float32,
                                                   ('yc', 'xc'))
            newvar['lon'].long_name = "longitude coordinate"
            newvar['lon'].standard_name = "longitude"
            newvar['lon'].units = "degrees_east"
            newvar['lon'][:] = data['lon']

            for key, value in pvardict.items():
                newvar[key] = dataset.createVariable(value, data[key].dtype,
                    dimshape)
                for mkey, mval in metadata[key].items():
                    if mkey == 'grid_mapping':
                        newvar[key].setncattr(mkey, crsname)
                    elif mkey != 'standard_name':
                        newvar[key].setncattr(mkey, mval)
                newvar[key][:] = data[key]

            # Global Attributes

            for attr in list(json_dict['attributes'].keys()):
                if attr == 'title':
                    setattr(dataset, attr,
                            "Per-month parameters for wind-driven free-drift model of sea-ice drift, accompanies the Global Sea Ice Drift Climate Data Record Version 1 from the EUMETSAT OSI SAF")
                elif attr == 'summary':
                    setattr(dataset, attr,
                            "Per-month parameters for wind-driven free-drift model of sea-ice drift obtained by solving the inverse free-drift problem with the ECMWF/C3S ERA5 wind reanalysis and coarse resolution passive microwave satellite data from AMSR-E and AMSR2. This dataset was generated by the EUMETSAT Ocean and Sea Ice Satellite Application Facility (OSI SAF) and accompanies the Global Sea Ice Drift Climate Data Record Version 1 with doi 0.15770/EUM SAF OSI 0012")
                elif attr == 'algorithm':
                    setattr(dataset, attr, "solving for inverse model with least squares")
                elif attr == 'time_coverage_duration':
                    setattr(dataset, attr, "P1M")
                elif attr == 'time_coverage_resolution':
                    setattr(dataset, attr, "P1M")
                elif attr == "keywords":
                    setattr(dataset, attr,
                            json_dict['attributes'][attr].format(hs_dict[hemi]))
                elif attr == "northernmost_latitude":
                    norlat = np.max(data['lat'])
                    if norlat > 89.5:
                        norlat = 90.
                    setattr(dataset, attr, norlat)
                elif attr == "southernmost_latitude":
                    soulat = np.min(data['lat'])
                    if soulat < -89.5:
                        soulat = -90.
                    setattr(dataset, attr, soulat)
                elif attr == "easternmost_longitude":
                    easlon = np.max(data['lon'])
                    if easlon > 179.5:
                        easlon = 180.
                    setattr(dataset, attr, easlon)
                elif attr == "westernmost_longitude":
                    weslon = np.min(data['lon'])
                    if weslon < -179.5:
                        weslon = -180.
                    setattr(dataset, attr, weslon)
                elif attr == "sensor":
                    setattr(dataset, attr, 'AMSR-E,AMSR2')
                elif attr == "platform":
                    setattr(dataset, attr, 'Aqua,GCOM-W1')
                elif attr == "source":
                    setattr(dataset, attr, 'FCDR of AMSR-E Brightness Temperatures V3 (doi: 10.5067/AMSR-E/AE_L2A.003), AMSR2 L1R data from JAXA, C3S/ECMWF ERA5 data')
                elif attr == 'time_coverage_start':
                    setattr(dataset, attr,
                            json_dict['attributes'][attr].format(d=mindate))
                elif attr == 'time_coverage_end':
                    setattr(dataset, attr,
                            json_dict['attributes'][attr].format(d=maxdate))
                elif attr in ['history', 'date_created']:
                    setattr(dataset, attr,
                            json_dict['attributes'][attr].format(
                            d=datetime.utcnow()))
                elif attr == 'geospatial_bounds_crs':
                    setattr(dataset, attr,
                            json_dict['attributes'][attr][gridname[hemi]])
                elif ( attr == 'tracking_id' ):
                    setattr(dataset, attr, get_tracking_id())
                else:
                    setattr(dataset, attr, json_dict['attributes'][attr])

        print('Writing out {}'.format(fnameandpath))

    return


def main():

    args = parse_args()
    indir = args.indir
    outdir = args.outdir

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    concatenate_pfiles(indir, outdir)


if __name__ == '__main__':

    main()
