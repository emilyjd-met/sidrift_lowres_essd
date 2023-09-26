'''Make a mask of the maximum ice drift from all files in a directory glob'''

import os
import sys
import argparse
import glob
import re
from datetime import datetime
import numpy as np
from netCDF4 import Dataset, date2num

sys.path.append(os.path.join(os.path.dirname(__file__), '../../merge/src'))
from read_icedrift_product_file import read_drift


def read_args():
    '''Read and pre-process user input'''

    parser = argparse.ArgumentParser(description='Ice drift using winds')

    parser.add_argument('-i', '--indirs', required=True,
                        help='Comma-separated list of input directories')

    parser.add_argument('-o', '--outfile', required=True,
                        help='Output filepath')

    parser.add_argument('-z', '--hemi', required=False, default='nh',
                        help="Hemisphere, 'nh' or 'sh'")

    args = parser.parse_args()

    return args


def make_max_drift_mask(indirs, outfile, hemi='nh'):

    indirlist = indirs.split(',')
    infiles = []
    for indir in indirlist:
        infiles.extend(glob.glob(os.path.join(indir, '*/*/*/icedrift*_{}*nc'.format(hemi))))

    ## Make a starting grid
    #with Dataset(infiles[0], 'r') as refdataset:
    #    dx = refdataset['driftX'][0, :, :]
    #    refshape = dx.shape
    #mask = np.zeros(refshape, dtype=np.int8)

    # Find information from the first file
    driftdata = read_drift(infiles[0])
    area_def = driftdata['area_def']
    ny, nx = area_def.shape
    lons, lats = area_def.get_lonlats()

    # Make a starting grid
    mask = np.zeros((ny, nx), dtype=np.int8)

    start_date = datetime.strptime('20100101', '%Y%m%d')
    end_date = datetime.strptime('20100102', '%Y%m%d')
    # Find the maximum drift mask
    for infile in infiles:
        with Dataset(infile, 'r') as dataset:
            tbnds = dataset['time_bnds'][:]
            lbnd = datetime.fromtimestamp(tbnds[0][0])
            ubnd = datetime.fromtimestamp(tbnds[0][1])
            if lbnd < start_date:
                start_date = lbnd
            if ubnd > end_date:
                end_date = ubnd
            dx = dataset['driftX'][0, :, :]
            mask[~dx.mask] = 1

    # Writing out the file
    with Dataset(outfile, 'w') as dataset:
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
        xc[:] = area_def.projection_x_coords

        yc = dataset.createVariable('yc', np.float64, ('yc'))
        yc.axis = "Y"
        yc.units = "km"
        yc.long_name = "y coordinate of projection (northings)"
        yc.standard_name = "projection_y_coordinate"
        yc[:] = area_def.projection_y_coords

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

        msk = dataset.createVariable('maxdriftmask', np.int8, ('yc', 'xc'))
        msk.long_name = "Maximum ice drift mask"
        msk.grid_mapping = crsname
        msk.coordinates = "lat lon"
        msk[:] = mask

    print("Wrote {}".format(outfile))


def main():

    args = read_args()
    indirs = args.indirs
    outfile = args.outfile
    hemi = args.hemi

    make_max_drift_mask(indirs, outfile, hemi=hemi)


if __name__ == '__main__':

    main()
