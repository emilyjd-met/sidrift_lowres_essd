'''Determine which single buoy trajectory files in a directory have data
   covering a particular area and time and write these trajectories to
   an output NetCDF file'''

import os
import sys
import re
import subprocess
from datetime import datetime, timedelta
import argparse
from argparse import RawDescriptionHelpFormatter
from glob import glob
import numpy as np
import pandas as pd
import xarray
from netCDF4 import Dataset, date2num


def parse_args():

    valid_areas = ['nh', 'sh']

    p = argparse.ArgumentParser("select_traj_2netcdf",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('indir',
                   help="Directory containing single buoy trajectory files")
    p.add_argument('-s', '--start_date', required=False, default=None,
                   help="Start date for which to select the trajectories. Can be of form YYYYMMDD in which case trajectories will be selected from YYYYMMDD 00:00:00, or YYYYMMDDHH, YYYYMMDDHHMM, YYYYMMDDHHMMSS")
    p.add_argument('-e', '--end_date', required=False, default=None,
                   help="End date for which to select the trajectories. Can be of form YYYYMMDD in which case trajectories will be selected up to YYYYMMDD+1day 00:00:00, or YYYYMMDDHH, YYYYMMDDHHMM, YYYYMMDDHHMMSS, in which case the exact time given will be used with 00 minutes and seconds if not specified. If not given, the start date plus three days will be used.")
    p.add_argument('-d', '--duration', required=False, default=2,
                   help="Duration in days to select the trajectories over, default is 2 days. If both start_date and end_date are set, the duration is calculated from these.")
    p.add_argument('-a', '--area', required=False, default='nh',
                   choices=valid_areas,
                   help="Area for which to select trajectories, default='nh', choices={}".format(valid_areas))
    p.add_argument('-o', '--output', required=False, default='.',
                   help="Either an output directory in which to store the output NetCDF which will be given a standard name, or the full output filepath.")
    p.add_argument('-r', '--reproc', default=False, action='store_true',
                   help="If set, store files in YYYY/mm subdirectories.")
    p.add_argument('-v', '--verbose', required=False, action='store_true',
                   help="If set, write more to the console.")
    args = p.parse_args()

    return args


def id_name(longid, netw):
    '''The buoy names must be unique and 9 characters.

    NOTE: This adds a 'z' to the end of the ID name because the matchups code deletes the last character.
    '''

    # Take suffixes off antsid buoys
    if longid[-4:] in ['_0pT', '_2pT', '_3pT', '_6pT']:
        longid = longid[:-4]

    fixlen = 8
    longid = longid.replace('-', '').replace('_', '')

    # Changing name for the IABP 3000* and 3002* buoys
    if netw == 'iabp' and longid[0:3] == '300' and len(longid) >= 13:
        longid = '{}{}'.format('j',longid[7:-1])

    # Adding the network name or part of it to the short IDs
    lenid = len(longid)
    if lenid <= fixlen:
        longid = '{}{}'.format(netw[0:fixlen-lenid], longid)
        longid = '{}z'.format(longid)
        if len(longid) < fixlen + 1:
            return longid.ljust(fixlen + 1)
        return longid

    longidend = longid[-fixlen:]
    # Special cases
    if longidend == '00000000':
        return '{}z'.format(longid[0:fixlen])

    return '{}z'.format(longidend)


def select_trajectories(indir, start_date, end_date, duration, area, output,
                        reproc=False, verbose=True):

    # Checking variables
    if start_date is None and end_date is None:
        print("Must specify a start date with -s or an end date with -e")
        sys.exit()

    # Setting the start and end dates
    valid_datefmt = {8: '%Y%m%d',
                     10: '%Y%m%d%H',
                     12: '%Y%m%d%H%M',
                     14: '%Y%m%d%H%M%S'}
    if start_date is not None:
        try:
            datefmt = valid_datefmt[len(start_date)]
        except:
            raise ValueError("start_date should be a date string with YYYYMMDD, YYYYMMDDHH, YYYYMMDDHHMM or YYYYMMDDHHMMSS, found {}".format(start_date))
        sdate = datetime.strptime(start_date, datefmt)

    if end_date is not None:
        try:
            datefmt = valid_datefmt[len(end_date)]
        except:
            raise ValueError("end_date should be a date string with YYYYMMDD, YYYYMMDDHH, YYYYMMDDHHMM or YYYYMMDDHHMMSS, found {}".format(end_date))
        edate = datetime.strptime(end_date, datefmt)
        # If the end_date is specified as %Y%m%d, one day is added to include
        # all of the specified day in the data
        if datefmt == '%Y%m%d':
            edate = edate + timedelta(days=1)

    # If both start date and end date are set, the duration is overwritten
    if start_date is not None and end_date is not None:
        duration = (edate - sdate).days
    # If either start_date or end_date are not set, then these are calculated
    elif start_date is None:
        sdate = edate - timedelta(days=duration)
    elif end_date is None:
        edate = sdate + timedelta(days=duration)

    # Filelist of possible single buoy trajectories
    flist = glob(os.path.join(indir, '*.nc'))

    if area == 'nh':
        long_area = "Northern Hemisphere"
    elif area == 'sh':
        long_area = "Southern Hemisphere"

    # Setting up storage for the matching trajectories
    trajdict = {}
    # Time, lat, lon arrays of dimension (trajectory_id, record) with
    # oversized dimensions to start
    dtarr = np.full((1000, 10000), np.nan, dtype='datetime64[s]')
    latarr = np.full((1000, 10000), np.nan, dtype=np.float64)
    lonarr = np.full((1000, 10000), np.nan, dtype=np.float64)
    # Cycle through all *nc files in the directory
    trajid = 0
    maxreclen = 0

    for fl in flist:
        if verbose:
            print("Checking buoy file {}...".format(fl))
        with Dataset(fl, 'r') as dataset:

            # Filtering out files which do not match the area or which
            # have no data overlapping the required datetime period
            if dataset.area != long_area:
                #print("Area doesn't match.")
                continue
            f_start = datetime.strptime(dataset.time_coverage_start,
                                        "%Y-%m-%dT%H:%M:%SZ")
            if f_start > edate:
                #print("Data too late.")
                continue
            f_end = datetime.strptime(dataset.time_coverage_end,
                                        "%Y-%m-%dT%H:%M:%SZ")

            if f_end < sdate:
                #print("Data too early.")
                continue

            # With the remaining files, read the time/lat/lon data and
            # filter for datetime. Note that the datetime was ordered in
            # previous code
            dt = [np.datetime64(datetime.fromtimestamp(t), 's')
                  for t in dataset['time'][:]]
            dt = np.array(dt)
            lat = dataset['lat'][:]
            lon = dataset['lon'][:]
            filt = np.logical_and(dt >= sdate, dt <= edate)
            dtfilt = dt[filt]
            latfilt = lat[filt]
            lonfilt = lon[filt]

            # Check that the record length is not zero or only 1 (which
            # can happen if there are no records in the time period, even
            # if the start and end dates are consistent)
            reclen = len(dtfilt)
            if reclen <= 1:
                continue

            if verbose:
                print("Matching data found with file {}".format(fl))

            dtarr[trajid, 0:reclen] = dtfilt
            latarr[trajid, 0:reclen] = latfilt
            lonarr[trajid, 0:reclen] = lonfilt
            # Add the buoy ID and network to a dictionary
            if dataset['source'][:][0] == 0:
                src = 'argo'
            elif dataset['source'][:][0] == 1:
                src = 'gps'
            elif dataset['source'][:][0] == 2:
                src = 'irid'
            else:
                raise ValueError("Unrecognised source value {}".format(dataset['source'][:][0]))
            trajdict[trajid] = (dataset.id, dataset.network, src)

            if reclen > maxreclen:
                maxreclen = reclen
            trajid += 1

    # Trim the arrays to the number of trajectories and the maximum
    # record length
    dtarr = dtarr[0:trajid, 0:maxreclen]
    latarr = latarr[0:trajid, 0:maxreclen]
    lonarr = lonarr[0:trajid, 0:maxreclen]

    # Check that there is valid data
    if len(dtarr) == 0:
        raise ValueError("No valid trajectory data found for this time period and area")

    # lon from -180 to +180 to match original files
    lonarr[lonarr > 180.] = lonarr[lonarr > 180.] - 360.

    # Determining the min and max latitudes for the NetCDF output and
    # then setting the required fill value
    min_lat = np.nanmin(latarr)
    max_lat = np.nanmax(latarr)
    min_lon = np.nanmin(lonarr)
    max_lon = np.nanmax(lonarr)
    latarr[np.isnan(latarr)] = -1.e10
    lonarr[np.isnan(lonarr)] = -1.e10

    # Converting the buoy IDs and network names to character strings
    idstrings = [id_name(v[0], v[1]) for k, v in trajdict.items()]
    id = np.array([[char for char in idstr] for idstr in idstrings], dtype='S1')
    netstrings = [v[1].ljust(4) for k, v in trajdict.items()]
    network = np.array([[char for char in netstr] for netstr in netstrings],
                       dtype='S1')
    srcstrings = [v[2].ljust(4) for k, v in trajdict.items()]
    source = np.array([[char for char in srcstr] for srcstr in srcstrings],
                       dtype='S1')

    xdata = xarray.Dataset(
        data_vars=dict(
            id=(["station", "nch_0"], id,
                dict(long_name="id of each station")),
            network=(["station", "nch_1"], network,
                     dict(long_name="network for each station")),
            source=(["station", "nch_2"], source,
                     dict(long_name="technique for geo-positioning for each station")),
            time=(["station", "record"], dtarr,
                  dict(_FillValue=-2147483648)),
            lat=(["station", "record"], latarr,
                 dict(long_name="latitude of position record",
                      units="degrees north", _FillValue=-1.e10)),
            lon=(["station", "record"], lonarr,
                 dict(long_name="longitude of position record",
                      units="degrees east", _FillValue=-1.e10)),
        ),
        attrs=dict(start_date_and_time=datetime.strftime(sdate,
                                                         "%Y-%m-%d %H:%M:%S"),
                   end_date_and_time=datetime.strftime(edate,
                                                       "%Y-%m-%d %H:%M:%S"),
                   area=long_area,
                   northernmost_latitude=max_lat,
                   southernmost_latitude=min_lat,
                   easternmost_longitude=max_lon,
                   westernmost_longitude=min_lon
               ),
    )

    if verbose:
        print("trajdict = ", trajdict)
        print("===========================")
        print(xdata)

    # Assembling the output name for the NetCDF file
    if output.endswith('nc'):
        outname = output
    else:
        if reproc:
            output = os.path.join(output, datetime.strftime(edate, '%Y'),
                                  datetime.strftime(edate, '%m'))
        if not os.path.isdir(output):
            os.makedirs(output)
        fbase = "insitu-validation-{}-{}-{}.nc".format(area,
            datetime.strftime(sdate, '%Y%m%d%H'),
            datetime.strftime(edate, '%Y%m%d%H'))
        outname = os.path.join(output, fbase)
    tmpname = "{}_tmp.nc".format(outname[:-3])

    # Writing out the NetCDF file. TODO: For now this needs to be NetCDF3
    # classic to interface with the C code. If the C code is updated/
    # replaced, update this to NetCDF4.
    print("Writing file to {}".format(outname))
    # WARNING! The C matchups code expects the time in units of seconds-
    # since-start-time, so be sure to use this. The C code takes the
    # start time from the globa variable start_date_and_time
    sunits = "seconds since {}".format(datetime.strftime(sdate,
                                                        '%Y-%m-%d %H:%M:%S'))
    xdata.to_netcdf(tmpname, format="NETCDF3_CLASSIC",
                    encoding={"time": {"dtype": "int32", "units": sunits},
                              "lat": {"dtype": "float32"},
                              "lon": {"dtype": "float32"}})

    # Call NCO to get the char strings formatted exactly like the original.
    # TODO: Find better method.
    subprocess.check_call(['ncwa', '-O', '-h', '-a', 'string1', tmpname,
                           outname])
    subprocess.check_call(['rm', '-fr', tmpname])


def main():

    args = parse_args()
    indir = args.indir
    start_date = args.start_date
    end_date = args.end_date
    duration = int(args.duration)
    area = args.area
    output = args.output
    reproc = args.reproc
    verbose = args.verbose

    select_trajectories(indir, start_date, end_date, duration, area, output,
                        reproc, verbose=verbose)


if __name__ == '__main__':

    main()
