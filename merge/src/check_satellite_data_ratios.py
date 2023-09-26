'''Check number of valid ice data per day and plot graph of ratio with ice area, for purpose of deciding threshold to say a satellite image has enough valid data'''

import os
import sys
import argparse
import glob
import re
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset, date2num, num2date
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')

sys.path.append(os.path.join(here, '../../merge/src'))
from read_icedrift_product_file import read_drift


def read_args():
    '''Read and pre-process user input'''

    parser = argparse.ArgumentParser(description='Ice drift using winds')

    parser.add_argument('-s', '--startdate', required=False, default='20201231',
                        help='Start date for the record')

    parser.add_argument('-e', '--enddate', required=False, default='19910101',
                        help='End date for the record')

    parser.add_argument('-d', '--duration', required=False, default=1,
                        help='Duration of each ice drift')

    parser.add_argument('-z', '--hemi', required=False, default='nh',
                        help="Hemisphere, 'nh' or 'sh'")

    parser.add_argument('-i', '--indir', required=True,
                        help='Top-level input directory')

    parser.add_argument('-o', '--outfile', required=True,
                        help='Output textfile filepath')

    parser.add_argument('-p', '--plotfile', required=True,
                        help='Output plotfile filepath')

    args = parser.parse_args()

    return args


def check_satellite_data(startdate, enddate, duration, indir, outfile,
                         plotfile, hemi='nh'):

    instdict = {'ssmi-f10': '90',
                'ssmi-f11': '90',
                'ssmi-f13': '90',
                'ssmi-f14': '90',
                'ssmi-f15': '90',
                'ssmis-f16': '90',
                'ssmis-f17': '90',
                'ssmis-f18': '90',
                'amsr-aq': '37',
                'amsr2-gw1': '37'}

    sdate = datetime.strptime(startdate, '%Y%m%d')
    edate = datetime.strptime(enddate, '%Y%m%d')

    patt = 'icedrift_{}_tb{}v-Lap+tb{}h-Lap_simplex_lev2_{}-ease2-750_{}12w_{}12w.nc'

    # Find information from the first file
    fdate = edate - timedelta(days=int(duration))
    inst = list(instdict.keys())[-1]
    fname = patt.format(inst, instdict[inst], instdict[inst], hemi,
                        datetime.strftime(fdate, '%Y%m%d'),
                        datetime.strftime(edate, '%Y%m%d'))
    fpath = os.path.join(indir, inst, datetime.strftime(edate, '%Y'),
                         datetime.strftime(edate, '%m'),
                         datetime.strftime(edate, '%d'), fname)
    driftdata = read_drift(fpath)
    area_def = driftdata['area_def']
    ny, nx = area_def.shape

    # Make a pole hole mask
    if hemi == 'nh':
        with Dataset(fpath) as dataset:
            startdate = dataset['time_bnds'][:][0][0]
            enddate = dataset['time_bnds'][:][0][1]
            tunits = dataset['time_bnds'].units
            phsdate = num2date(startdate, units=tunits)
            phedate = num2date(enddate, units=tunits)
            polehole = np.zeros_like(dataset['flag'][0, :, :])
            # Lat >= 86 seems to give a pole hole 1 pixel bigger than
            # SSMIs for an EASE2 grid
            polehole[dataset['lat'][:] >= 86.] = 1
        write_testfile(os.path.join(datadir, 'info/polehole.nc'), area_def, phsdate, phedate, polehole, testfield2=driftdata['driftX'])

    # Initialising the list of all counts
    ddatesl = []
    driftsl = []
    ratiosl = []
    instsl = []
    fpathsl = []

    with open(outfile, 'a+') as outf:

        # Cyling through the dates
        ddate = sdate
        while ddate < edate:

            # Incrementing the date
            ddate =  ddate + timedelta(days=1)
            mdate = ddate - timedelta(days=int(duration))
            # Check for season and throw out summer
            # Summer is Jun - Sep in NH
            # Summer is Nov - Feb in SH

            month = int(datetime.strftime(ddate, '%m'))
            if hemi == 'nh':
                if month >= 6 and month <= 9:
                    continue
            else:
                if month >= 11 or month <= 2:
                    continue

            # Find the filelist
            infiles = []
            ilist = []
            for inst in instdict.keys():
                fname = patt.format(inst, instdict[inst], instdict[inst],
                                    hemi, datetime.strftime(mdate, '%Y%m%d'),
                                    datetime.strftime(ddate, '%Y%m%d'))
                fpath = os.path.join(indir, inst,
                                     datetime.strftime(ddate, '%Y'),
                                     datetime.strftime(ddate, '%m'),
                                     datetime.strftime(ddate, '%d'), fname)
                if os.path.isfile(fpath):
                    driftgrid = np.zeros((ny, nx), dtype=np.int8)
                    icegrid = np.zeros((ny, nx), dtype=np.int8)
                    with Dataset(fpath, 'r') as dataset:
                        # Find where there is ice drift
                        dx = dataset['driftX'][0, :, :]
                        driftgrid[~dx.mask] = 1
                        # Not considering the pole hole in the NH
                        if hemi == 'nh':
                            driftgrid[polehole] = 0
                        driftnum = np.sum(driftgrid)
                        # Find where there is ice
                        flag = dataset['flag'][0, :, :]
                        icegrid[np.logical_or(flag <= 0, flag > 5)] = 1
                        # Not considering the pole hole in the NH
                        if hemi == 'nh':
                            icegrid[polehole] = 0
                        icenum = np.sum(icegrid)
                        # Calculate the ratio of drift vectors to ice pixels
                        ratio = driftnum / icenum
                    outf.write("{}\t{}\t{}\t{}\n".format(ddate, driftnum,
                                                         ratio, inst))

                    ddatesl.append(ddate)
                    driftsl.append(driftnum)
                    ratiosl.append(ratio)
                    instsl.append(inst)
                    fpathsl.append(fpath)

    ddates  = [x for y, x in sorted(zip(ratiosl, ddatesl))]
    drifts  = [x for y, x in sorted(zip(ratiosl, driftsl))]
    insts  = [x for y, x in sorted(zip(ratiosl, instsl))]
    fpaths  = [x for y, x in sorted(zip(ratiosl, fpathsl))]
    ratios = sorted(ratiosl)

    dirn = os.path.dirname(outfile)
    basen = os.path.basename(outfile)
    sortfile = os.path.join(dirn, 'sorted_{}'.format(basen))
    with open(sortfile, 'a+') as sortf:
        for i in range(len(ratios)):
            sortf.write("{}\t{}\t{}\t{}\n".format(ddates[i], drifts[i],
                                                  ratios[i], insts[i]))

    nccmd1 = os.path.join(dirn, '{}_rp0-p3.sh'.format(basen[:-4]))
    nccmd2 = os.path.join(dirn, '{}_rp3-p4.sh'.format(basen[:-4]))
    nccmd3 = os.path.join(dirn, '{}_rp4-p5.sh'.format(basen[:-4]))
    ncdict = {nccmd1: (0.0, 0.3),
              nccmd2: (0.3, 0.4),
              nccmd3: (0.4, 0.5)}
    for outnc in ncdict.keys():
        with open(outnc, 'a+') as outn:
            outn.write("ncview ")
            for i, rat in enumerate(ratios):
                if rat > ncdict[outnc][0] and rat <= ncdict[outnc][1]:
                    outn.write("{} ".format(fpaths[i]))


    # Plotting a histogram
    fig = plt.figure(figsize=(8,6))
    plt.hist(ratios, density=True, bins=30)
    plt.ylabel('Number of days')
    plt.xlabel('Ratios')
    plt.savefig(plotfile)


def write_testfile(testncname, area_def, sdate, edate, testfield1,
                   testfield2=None):
    '''NOTE: Be careful of data types.'''

    ny, nx = area_def.shape
    lons, lats = area_def.get_lonlats()

    # Writing out the file
    with Dataset(testncname, 'w') as dataset:
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
        time[:] = date2num(edate, tunits)

        tbounds = dataset.createVariable('time_bnds',
                                         np.float64, ('time', 'nv'))
        tbounds.units = tunits
        tbounds[:] = [date2num(sdate, tunits), date2num(edate, tunits)]

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

        field1 = dataset.createVariable('field1', np.int8, ('yc', 'xc'))
        field1.long_name = "Test field 1"
        field1.grid_mapping = crsname
        field1.coordinates = "lat lon"
        field1[:] = testfield1

        if testfield2 is not None:
            field2 = dataset.createVariable('field2', np.float32, ('yc', 'xc'))
            field2.long_name = "Test field 2"
            field2.grid_mapping = crsname
            field2.coordinates = "lat lon"
            field2[:] = testfield2

    print("Wrote {}".format(testncname))


def main():

    args = read_args()
    startdate = args.startdate
    enddate = args.enddate
    duration = args.duration
    hemi = args.hemi
    indir = args.indir
    outfile = args.outfile
    plotfile = args.plotfile

    check_satellite_data(startdate, enddate, duration, indir, outfile,
                         plotfile, hemi=hemi)


if __name__ == '__main__':

    main()
