'''Convert all types of buoy txt files to a common format NetCDF file'''

# TODO - Need a function which examines the lat/lon positions and determines
# the source (argos/gps) based on how erratic the data is

import os
import pathlib
import re
import sys
import traceback
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray
import argparse
from argparse import RawDescriptionHelpFormatter
import csv
from glob import glob
from netCDF4 import Dataset, stringtochar, num2date
from pyresample import utils, AreaDefinition, SwathDefinition, kd_tree, spherical
import cartopy
import cartopy.crs as ccrs

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')


def parse_args():

    default_buoydir = os.path.join(datadir, 'SID_BUOY_DATA/')
    valid_buoycats = ['aari', 'antsid', 'argos', 'awi', 'bbb', 'crrel',
                      'hudson', 'iabp', 'itp', 'pipers', 'sams', 'sedna',
                      'simba', 'sip', 'tara']
    default_outdir = '.'

    p = argparse.ArgumentParser("convert_buoy_txt2netcdf",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('buoyname',
                   help="Name of the buoy to fetch data for")
    p.add_argument('buoycat', choices=valid_buoycats,
                   help="Name of the buoy catalogue to fetch data from, valid options are {}".format(valid_buoycats))    
    p.add_argument('--buoydir', required=False,
                   default = default_buoydir,
                   help="Top-level directory for buoy data, default {}".format(default_buoydir))
    p.add_argument('--outdir', required=False, default=default_outdir,
                   help="Output directory for the NetCDF file if not specified, default {}".format(default_outdir))
    p.add_argument('--outname', required=False, default=None,
                   help="Output filename for the NetCDF file if not specified")
    p.add_argument('-f', '--force', action='store_true', default=False,
                   help="Force creation of new NetCDF files if they exist. Default behaviour is that they will not be created if they exist and have a timestamp newer than the text file")

    args = p.parse_args()

    return args


def netcdf_attr_dict(buoycat):
    '''Return NetCDF attributes specific to the buoy catalogue'''

    ncatt = {'aari': {'catname': 'Arctic and Antarctic Research Institute'},
             'antsid': {'catname': 'Atlas of Antarctic Sea Ice Drift 2004'},
             'awi': {'catname': 'Alfred Wegener Institute buoys'},
             'bbb': {'catname': 'Baffin Bay buoys'},
             'argos': {'catname': 'Argos Damocles buoys'},
             'crrel': {'catname': 'Cold Regions Research and Engineering Laboratory'},
             'hudson': {'catname': 'Hudson Bay GPS collar'},
             'iabp': {'catname': 'International Arctic Buoy Programme'},
             'itp': {'catname': 'Ice-Tethered Profiler buoys from the Woods Oceanographic Institute'},
             'pipers': {'catname': 'Polynyas, Ice Production, and seasonal Evolution in the Ross Sea expedition'},
             'sams': {'catname': 'Scottish Association for Marine Science'},
             'sedna': {'catname': 'Sea Ice Experiment - Dynamic Nature of the Arctic'},
             'simba': {'catname': 'Sea Ice Mass Balance in the Antarctic'},
             'tara': {'catname': 'Tara'},
             'sip': {'catname': 'SIP'},
    }
    return ncatt[buoycat]


def read_aari(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the Arctic and Antarctic Research Institute'''

    try:
        bname = '{}.vector.txt'.format(buoyname)
        bfile = os.path.join(buoydir, 'aari', bname)
        assert os.path.isfile(bfile)
    except AssertionError:
        bname = '{}.drift.1h.txt'.format(buoyname)
        bfile = os.path.join(buoydir, 'aari', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('lat', 'f8'), ('lon', 'f8'), ('index', 'i4'), ('x', 'f8'),
                 ('y', 'f8'), ('date', '<U10'), ('time', '<U8'), ('sog', 'f8'),
                 ('cog', 'i4'), ('track', 'f8'), ('days', 'f8'), ('dt', 'f8')]

    # Use invalid_raise as there is a file which has one line with the
    # incorrect number of columns
    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, comments='#',
                             invalid_raise=False)
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    fixed = datetime.strptime('00:00:00', '%H:%M:%S')
    try:
        datecol = [datetime.strptime(d, '%Y.%m.%d')
                   + (datetime.strptime(t, '%H:%M:%S') - fixed)
                   for d, t in zip(data['date'], data['time'])]
    except:
        datecol = [datetime.strptime(d, '%Y-%m-%d')
                   + (datetime.strptime(t, '%H:%M:%S') - fixed)
                   for d, t in zip(data['date'], data['time'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def read_antsid(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the Atlas of Antarctic Sea Ice Drift 2004'''

    bname = '{}.txt'.format(buoyname)
    bfile = os.path.join(buoydir, 'AtlasofAntarcticSeaIceDrift-2004', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('id','f8'), ('year','f8'), ('doy','f8'), ('dom','f8'),
                 ('mon','f8'), ('lon','f8'), ('lat','f8'), ('p','f8'),
                 ('tair','f8')]

    # Note: Use unicode-escape as some files have non-ascii characters
    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, comments='!',
                             encoding= 'unicode_escape')
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    datecol = [datetime(int(y), 1, 1) + timedelta(d)
               for y, d in zip(data['year'], data['doy'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # Replace bad data
    pddata[pddata == 999.999] = np.nan
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def skip_no_loc_lines(fl):
    '''Code snippet to skip lines in AWI files where no GPS location'''

    for line in fl:
        if not b"GPS position better than 100m" in line:
            # Don't yield this line, we don't want it
            continue
        yield line


def read_awi(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the AWI buoys'''

    bnames = ['ANT-XIX_2_{}_buoy_data.tab'.format(buoyname),
              'ANT-XXI_4_{}_buoy_data.tab'.format(buoyname)]
    bfiles = [os.path.join(buoydir, 'pangaea_ipab', bname) for bname in bnames]
    bfiles = [bf for bf in bfiles if os.path.isfile(bf)]

    # Checking if the file creation should be skipped
    if out_ts is not None:
        for bfile in bfiles:
            pname = pathlib.Path(fl)
            txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
            if out_ts > txt_ts:
                return None, bfiles

    dtypefile = [('datetime','<U16'), ('lat','f8'), ('lon','f8')]

    try:
        data = np.array([], dtype=dtypefile)
        for bfile in bfiles:
            # NOTE: Open file in binary mode for use with generator
            with open(bfile, 'rb') as fl:
                skipped = skip_no_loc_lines(fl)
                # Just read the first 3 columns
                fldata = np.genfromtxt(skipped, usecols=[0,1,2], dtype=dtypefile)
                data = np.append(data, fldata, axis=0)
    except:
        raise IOError('Files {} are not readable.'.format(bfiles))

    # Finding the date column
    datecol = [datetime.strptime(d, '%Y-%m-%dT%H:%M') for d in data['datetime']]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (GPS here)
    source = [1] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfiles


def read_bbb(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the Baffin Bay buoy'''

    bname = 'bbb_baffin_bay_{}.txt'.format(buoyname)
    bfile = os.path.join(buoydir, 'bbBuoy', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('lat','f8'), ('lon','f8'), ('index', 'i4'),
                 ('date','<U10'), ('time','<U8'), ('p','f8'),
                 ('q','f8')]

    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, comments='!')
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    fixed = datetime.strptime('00:00:00', '%H:%M:%S')
    datecol = [datetime.strptime(d, '%Y.%m.%d')
               + (datetime.strptime(t, '%H:%M') - fixed)
               for d, t in zip(data['date'], data['time'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def read_argos(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of Argos Damocles buoys'''

    bname = 'argos_damocles{}.dat'.format(buoyname)
    bfile = os.path.join(buoydir, 'buoyarray', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('no','i4'), ('buoy','i8'), ('deziday','f8'), ('year','i4'),
                 ('mon','i4'), ('day','i4'), ('hour','i4'), ('min','i4'),
                 ('sec','i4'), ('lat','f8'), ('lon','f8')]

    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, skip_header=13)
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    datecol = [datetime(year=y, month=m, day=d, hour=h, minute=mi, second=s)
               for y, m, d, h, mi, s in
               zip(data['year'], data['mon'], data['day'], data['hour'],
                   data['min'], data['sec'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # Replace bad data
    pddata[pddata == -999.99902] = np.nan
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def read_crrel(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the crrel buoys'''

    bname = '{}_cleanPos.csv'.format(buoyname)
    bfile = os.path.join(buoydir, 'crrel', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

#    dtypefile = [('date','<U18'), ('lonstr','<U10'), ('latstr','<U10'),
#                 ('gps','<U3'),
#                 ('c1','<U10'), ('c2','<U10'), ('c3','<U10'), ('c4','<U10'),
#                 ('c5','<U10'), ('c6','<U10'), ('c7','<U10'), ('c8','<U10'),
#                 ('c9','<U10'), ('c10','<U10'), ('c11','<U10'), ('c12','<U10'), 
#                 ('c13','<U10'), ('c14','<U10'), ('c15','<U10'),
#                 ('c16','<U10'), ('c17','<U10'), ('c18','<U10'),
#                 ('c19','<U10'), ('c20','<U10'), ('c21','<U10'),
#                 ('c22','<U10'), ('c23','<U10'), ('c24','<U10'),
#                 ('c25','<U10'), ('c26','<U10'), ('c27','<U10'),
#                 ('c28','<U10'), ('c29','<U10'), ('c30','<U10'),
#                 ('c31','<U10'), ('c32','<U10'), ('c33','<U10'),
#                 ('c34','<U10'), ('c35','<U10'), ('c36','<U10'),
#                 ('c37','<U10')]
#
#    try:
#        data = np.genfromtxt(bfile, dtype=dtypefile, delimiter=',',
#                             converters={0:lambda x: x.replace('\"','')})
#    except:
#        raise IOError('File {} is not readable.'.format(bfile))

#    data['latstr'] = [re.sub('"', '', x).strip() for x in data['latstr']]
#    data['lonstr'] = [re.sub('"', '', x).strip() for x in data['lonstr']]
#    data['lat'] = [float(y) for y in data['latstr']]
#    data['lon'] = [float(y) for y in data['lonstr']]

    # Messy way of reading a CSV file, but could not get genfromtxt to work -
    # issue with bytes/string type to convert "s away
    try:
        dates = []
        lats = []
        lons = []
        sources = []
        with open(bfile) as csvfile:
            csvdata = csv.reader(csvfile)
            for row in csvdata:
                dates.append(row[0])
                lats.append(row[1])
                lons.append(row[2])
                sources.append(row[3])
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    datecol = [datetime.strptime(dt, '%m/%d/%Y %H:%M') for dt in dates]

    # lons > 0 and filtering for valid lonlats
    # lons = [float(l.strip()) for l in lons]
    # lons = [(l + 360. if l < 0. else l) for l in lons]
    # lats = [float(l.strip()) for l in lats]
    lons_filt = []
    lats_filt = []
    datecol_filt = []
    sources_filt = []
    for i in range(len(datecol)):
        # Use try/except loops to skip over empty entries like ''
        try:
            lat = float(lats[i].strip())
        except:
            break
        try:
            lon = float(lons[i].strip())
        except:
            break
        if lon < 0:
            lon = lon + 360.
        if not (lon >= 0.) and (lon < 360.):
            break
        if not (lat >= -90.) and (lat <= 90.):
            break
        if re.search("gps", sources[i], re.IGNORECASE):
            src = 1
        else:
            src = 0
        lons_filt.append(lon)
        lats_filt.append(lat)
        datecol_filt.append(datecol[i])
        sources_filt.append(src)

    # Creating the data dictionary
    data = {}
    data['datecol'] = datecol_filt
    data['lat'] = lats_filt
    data['lon'] = lons_filt
    data['source'] = sources_filt

    # Create the dataframe%
    pddata = pd.DataFrame(data, index=data['datecol'], columns=('lat', 'lon',
                                                                'source'))
    # Replace bad data
    pddata[pddata == 999.999] = np.nan
    # The datetime column is also needed
    pddata.insert(0, 'datetime', data['datecol'])

    return pddata, bfile


def read_hudson(buoyname, buoycat, buoydir, out_ts):
    '''The Hudson Bay data is in eastings/northings in the NAD83(CSRS) / Teranet Ontario Lambert projection (EPSG:5321) coordinate system'''

    bfile = os.path.join(buoydir, 'hudson', 'Drifting_collar_location_data.tab')

    dtypefile = [('buoyname', '<U8'), ('year', 'f4'), ('date', '<U10'),
                 ('time', '<U5'), ('posx', 'f8'), ('posy', 'f8')]

    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, skip_header=1)
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    fixed = datetime.strptime('00:00:00', '%H:%M:%S')
    datecol = [datetime.strptime(d, '%Y-%m-%d')
               + (datetime.strptime(t, '%H:%M') - fixed)
               for d, t in zip(data['date'], data['time'])]

    # Define EPSG:5321 coordinate system
    proj4_str = "+proj=lcc +lat_1=44.5 +lat_2=54.5 +lat_0=0 +lon_0=-84 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +type=crs" # +units=m +no_defs"
    try:
        proj_dict = utils._proj4.proj4_str_to_dict(proj4_str)
    except:
        proj_dict = utils.proj4.proj4_str_to_dict(proj4_str)
    area_extent = [-3454415.91, 6848343.70, 1381376.16, 10828392.64]
    epsg5321_areadef = AreaDefinition('epsg5321', 'epdg5321', 'epsg5321',
                                      proj_dict, 1000, 1000, area_extent)
    epsg5321_ccrs = epsg5321_areadef.to_cartopy_crs()

    # Transforming the coordinates to lat/lon
    lat = []
    lon = []
    for i in range(len(data['posx'])):
        lonp, latp = ccrs.PlateCarree().transform_point(data['posx'][i], data['posy'][i], src_crs=epsg5321_ccrs)
        lat.append(latp)
        lon.append(lonp)

    # Creating the data dictionary
    data = {}
    data['datecol'] = datecol
    data['lat'] = lat
    data['lon'] = lon
    data['source'] = [0] * len(datecol)

    # Create the dataframe%
    pddata = pd.DataFrame(data, index=data['datecol'], columns=('lat', 'lon',
                                                                'source'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', data['datecol'])

    return pddata, bfile


def read_iabp(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the IABP buoys from the IABP_archive dir'''

    bname = '{}.ll.pos'.format(buoyname)
    flist = glob(os.path.join(buoydir, 'IABP_archive', '*', '3HOURLY', bname))

    # Checking if the file creation should be skipped
    if out_ts is not None:
        for fl in flist:
            pname = pathlib.Path(fl)
            txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
            if out_ts > txt_ts:
                return None, fl

    lenstart = len(os.path.join(buoydir, 'IABP_archive'))
    lenend = len(os.path.join('3HOURLY', bname))
#    years = [x[lenstart+1:-(lenend+1)] for x in flist]

    dtypeshfile = [('id','f8'), ('fracdoy','f8'), ('lat','f8'), ('lon','f8')]
    dtypefile = [('id','f8'), ('fracdoy','f8'), ('lat','f8'), ('lon','f8'),
                 ('year', 'i4')]

    try:
#        data = np.array([], dtype=dtypefile)
#        for fl in flist:
#            year = fl[lenstart+1:-(lenend+1)]
#            data = np.append(data, np.genfromtxt(fl, dtype=dtypefile), axis=0)
        data = np.array([], dtype=dtypefile)
        for fl in flist:
            year = int(fl[lenstart+1:-(lenend+1)])
            fldata = np.genfromtxt(fl, dtype=dtypeshfile)
            flyrdata = [(x[0], x[1], x[2], x[3], year) for x in fldata]
            flyrdatanp = np.array(flyrdata, dtype=dtypefile)
            data = np.append(data, flyrdatanp, axis=0)
    except:
        raise IOError('Files {} are not readable.'.format(flist))

    # Finding the date column
    datecol = [datetime(y, 1, 1) + timedelta(d - 1)
               for y, d in zip(data['year'], data['fracdoy'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, str(flist)


def read_iabp_new(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the IABP buoys from the IABP_3HOURLY_DATA dir'''

    bnamefr = '{}.dat'.format(buoyname)
    bname3h = '{}.csv'.format(buoyname)
    flist = []
    flist.extend(glob(os.path.join(buoydir, 'IABP_3HOURLY_DATA/iabp.apl.uw.edu/Data_Products/BUOY_DATA/FULL_RESOLUTION_DATA/', '*', '*', bnamefr)))
    flist.extend(glob(os.path.join(buoydir, 'IABP_3HOURLY_DATA/iabp.apl.uw.edu/Data_Products/BUOY_DATA/3HOURLY_DATA/', '*', bname3h)))

    print("flist = ", flist)
    # Checking if the file creation should be skipped
    if out_ts is not None:
        for fl in flist:
            pname = pathlib.Path(fl)
            txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
            if out_ts > txt_ts:
                return None, fl

    # Different formats of these files
    # (3-hourly dat files exist, but these are the same)
    # 3-hourly csv
    # BuoyID,Year,Hour,DOY,Lat,Lon + 0-3 sci output
    # Full-res dat
    # BuoyID,Year,Hour,Min,DOY,POS_DOY,Lat,Lon + 0-3 sci output

    dtypefile3h = [('id','f8'), ('year', 'i4'), ('hour', 'i4'),
                   ('fracdoy','f8'), ('lat','f8'), ('lon','f8')]
    dtypefilefr = [('id','f8'), ('year', 'i4'), ('hour', 'i4'),
                   ('min', 'i4'), ('fracdoy','f8'), ('posdoy', 'f8'),
                   ('lat','f8'), ('lon','f8')]

    dtypefile3h1s = dtypefile3h.copy()
    dtypefile3h1s.append(('sci1', 'f8'))
    dtypefile3h2s = dtypefile3h1s.copy()
    dtypefile3h2s.append(('sci2', 'f8'))
    dtypefile3h3s = dtypefile3h2s.copy()
    dtypefile3h3s.append(('sci3', 'f8'))
    dtypefilefr1s = dtypefilefr.copy()
    dtypefilefr1s.append(('sci1', 'f8'))
    dtypefilefr2s = dtypefilefr1s.copy()
    dtypefilefr2s.append(('sci2', 'f8'))
    dtypefilefr3s = dtypefilefr2s.copy()
    dtypefilefr3s.append(('sci3', 'f8'))
    dtypedict = {}
    dtypedict['3hour'] = {6: dtypefile3h,
                          7: dtypefile3h1s,
                          8: dtypefile3h2s,
                          9: dtypefile3h3s}
    dtypedict['fullres'] = {8: dtypefilefr,
                            9: dtypefilefr1s,
                            10: dtypefilefr2s,
                            11: dtypefilefr3s}

    try:
        data = np.array([], dtype=dtypefilefr)
        for fl in flist:
            with open(fl, 'r') as flf:
                fieldnum = len(flf.readline().split(','))
            if 'FULL_RESOLUTION' in fl:
                dtypefile = dtypedict['fullres'][fieldnum]
                flrddata = np.genfromtxt(fl, dtype=dtypefile,
                                         delimiter=',', skip_header=1)
                # Just take the first 8 fields
                fldata = [(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                           x[7]) for x in flrddata]
                fldatanp = np.array(fldata, dtype=dtypefilefr)
                data = np.append(data, fldatanp, axis=0)
            else:
                dtypefile = dtypedict['3hour'][fieldnum]
                flrddata = np.genfromtxt(fl, dtype=dtypefile,
                                         delimiter=',', skip_header=1)
                # Add minutes and posdoy (copy of fracdoy here) to 3h data
                fldata = [(x[0], x[1], x[2], '00', x[3], x[3], x[4],
                           x[5]) for x in flrddata]
                fldatanp = np.array(fldata, dtype=dtypefilefr)
                data = np.append(data, fldatanp, axis=0)
#            print("-- data.shape = ", data.shape)

    except:
        raise IOError('Files {} are not readable.'.format(flist))

#    print("== data.shape = ", data.shape)

    # Fix the years where they are only 2-digit
    for i in range(len(data)):
        if data[i][1] >= 0 and data[i][1] <= 25:
            data[i][1] += 2000
        elif data[i][1] >= 60 and data[i][1] <= 99:
            data[i][1] += 1900

    # Some of the lines have -999 in the date column, skip these lines
    # Also some have -90 or -180 in the lat and lon columns, also skip
    poplines = []
    for i in range(len(data)):
        if (data[i][1] < 1960
            or data[i][1] > 2025
            or int(data[i][6]) == -90
            or int(data[i][7]) == -180
            or int(data[i][7]) == 180):
            poplines.append(i)
    poplines.reverse()
    if poplines:
        for i in poplines:
            data = np.delete(data, (i))

    # Finding the date column
    datecol = [datetime(y, 1, 1) + timedelta(days = d - 1)
               for y, d in zip(data['year'], data['fracdoy'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, str(flist)


def read_itp(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the ITP buoys'''

    bname = '{}rawlocs.dat'.format(buoyname)
    bfile = os.path.join(buoydir, 'itp', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('year','i4'), ('fracdoy','f8'), ('lon','f8'), ('lat','f8')]

    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, comments='%')
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    datecol = [datetime(y, 1, 1) + timedelta(d - 1)
               for y, d in zip(data['year'], data['fracdoy'])]

    # lons > 0
    data['lon'] = [(l + 360. if l < 0. else l) for l in data['lon']]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def read_pipers(buoyname, buoycat, buoydir, out_ts):
    '''Read the nc files of the PIPERS buoys'''

    bfile = os.path.join(buoydir, 'pipers', 'PIPERS_NIWA_waves_ice_2017.nc')

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('time', 'i8'), ('lon','f8'), ('lat','f8')]
    data = np.array([], dtype=dtypefile)

    try:
        with Dataset(bfile, 'r') as dataset:
            buoy = dataset['buoy'][:]
            time = dataset['time'][:]
            lat = dataset['lat'][:]
            lon = dataset['lon'][:]
            tunits = dataset['time'].units
        buoynames = np.array([str(b) for b in buoy.data])
        idx = np.where(buoynames == buoyname)[0][0]
        blat = lat[idx, :]
        blon = lon[idx, :]
        msk = np.logical_or(blat.mask, blon.mask)
        ftime = time[~msk]
        flat = blat[~msk]
        flon = blon[~msk]

        for i in range(len(ftime)):
            single = np.array([(ftime[i], flon[i], flat[i])], dtype=dtypefile)
            data = np.append(data, single, axis=0)
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    sunits = str(tunits)
    unitpatt = "seconds since (\d{2}) (\d{2}) (\d{4}) (\d{2}):(\d{2}) UTC"
    m = re.match(unitpatt, sunits)
    unitstr = "seconds since {}-{}-{} {}:{} UTC"
    units = unitstr.format(m[3], m[2], m[1], m[4], m[5])
    datecol = [num2date(int(d), units) for d in data['time']]

    # lons > 0
    data['lon'] = [(l + 360. if l < 0. else l) for l in data['lon']]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def read_sams(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the SAMS buoys'''

    bname = '{}_qcgps.txt'.format(buoyname)
    if buoyname.startswith('Access'):
        bfile = os.path.join(buoydir, 'sams/ACCESS', bname)
    elif buoyname.startswith('Icebell'):
        bfile = os.path.join(buoydir, 'sams/ICEBELL', bname)
    elif buoyname.startswith('Kopri'):
        bfile = os.path.join(buoydir, 'sams/KOPRI', bname)
    elif buoyname.startswith('Naacos'):
        bfile = os.path.join(buoydir, 'sams/NAACOS', bname)
    elif buoyname.startswith('SIP'):
        bfile = os.path.join(buoydir, 'sams/SIP', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('qc','i4'), ('date','<U10'), ('time','<U8'), ('speed','f8'),
                 ('lat','f8'), ('lon','f8')]
    dtypeout = [('qc','i4'), ('datetime','<U10'), ('lat','f8'), ('lon','f8')]

    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, comments='#')
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    fixed = datetime.strptime('00:00:00', '%H:%M:%S')
    datecol = [datetime.strptime(d, '%Y-%m-%d')
               + (datetime.strptime(t, '%H:%M:%S') - fixed)
               for d, t in zip(data['date'], data['time'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon', 'qc'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)
    # Remove the entries with qc=0, and then remove the qc column
    pddata = pddata[pddata.qc != 0]
    del pddata['qc']

    return pddata, bfile


def read_sedna(buoyname, buoycat, buoydir, out_ts):
    '''Read the NetCDF files of the Sedna buoys'''

    bnamevel = 'Velocity/SEDNA_Buoy_{}_Velocity_2Months.nc'.format(buoyname)
    bnamestrn = 'Strain/SEDNA_Buoy_{}_GPS_2Months.nc'.format(buoyname)
    try:
        bfile = os.path.join(buoydir, 'sedna', bnamevel)
        assert os.path.isfile(bfile)
    except AssertionError:
        bfile = os.path.join(buoydir, 'sedna', bnamestrn)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    with Dataset(bfile, 'r') as dataset:
        lats = dataset['latitude'][:]
        lons = dataset['longitude'][:]
        timebytes = dataset['time_string'][:]
        sourcestr = dataset.source
        commentstr = dataset.comment

    # Finding the date column
    timestr = [''.join(x.astype(str)) for x in timebytes]
    datecol = [datetime.strptime(t, '%d-%b-%Y %H:%M:%S') for t in timestr]

    # Creating a source column
    if (re.search("gps", sourcestr, re.IGNORECASE) or
        re.search("gps", commentstr, re.IGNORECASE)):
        source = [1] * len(datecol)
    else:
        source = [0] * len(datecol)

    # Creating the data dictionary
    data = {}
    data['datecol'] = datecol
    data['lat'] = lats
    data['lon'] = lons
    data['source'] = source

    # Create the dataframe%
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon', 'source'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)

    return pddata, bfile


def read_simba(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the SIMBA buoys'''

    bname = '{}.Position.dat'.format(buoyname)
    bfile = os.path.join(buoydir, 'SIMBA', bname)

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    dtypefile = [('year','i4'), ('doy','i4'), ('mon','i4'), ('day','i4'),
                 ('hour','i4'), ('min','i4'), ('dectime','f8'), ('lon','f8'),
                 ('lat','f8')]

    try:
        data = np.genfromtxt(bfile, dtype=dtypefile, skip_header=1)
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    datecol = [datetime(year=y, month=m, day=d, hour=h, minute=mi)
               for y, m, d, h, mi in
               zip(data['year'], data['mon'], data['day'], data['hour'],
                   data['min'])]

    # lons > 0
    data['lon'] = [(l + 360. if l < 0. else l) for l in data['lon']]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def read_tara(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the tara buoys'''

    bname = 'tara_argos_drift.ori.*'
    flist = glob(os.path.join(buoydir, 'tara/argos', bname))

    dtypeshfile = [('id','i4'), ('latstr','<U8'), ('lonstr','<U8'),
                   ('qi','i4'), ('dateind', '<U18')]
    dtypefile = [('id','<U5'), ('latstr','<U8'), ('lonstr','<U8'),
                 ('lat', 'f8'), ('lon', 'f8'), ('qi','i4'),
                 ('datestr', '<U18'), ('year', '<U4')]

    # Checking if the file creation should be skipped
    if out_ts is not None:
        for fl in flist:
            try:
                pname = pathlib.Path(fl)
            except:
                pass
            txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
            if out_ts > txt_ts:
                return None, fl

    try:
        data = np.array([], dtype=dtypefile)
        for fl in flist:
            year = fl[-4:]
            fldata = np.genfromtxt(fl, dtype=dtypeshfile)
            flyrdata = [(x[0], x[1], x[2], 0.0, 0.0, x[3], x[4], year)
                        for x in fldata]
            flyrdatanp = np.array(flyrdata, dtype=dtypefile)
            data = np.append(data, flyrdatanp, axis=0)
    except:
        raise IOError('Files {} are not readable.'.format(flist))

    # Delete the rows from the pandas dataframe which do not match the named
    # buoy (since buoy name not in the filename)
    data = [d for d in data if d[0] == buoyname]
    data = np.array(data, dtype=dtypefile)

    # Edit the lats and lons to be integers
    for i, lon in enumerate(data['lonstr']):
        if lon.endswith('E'):
            data['lon'][i] = float(lon.replace('E', ''))
        elif lon.endswith('W'):
            data['lon'][i] = 360. - float(lon.replace('W', ''))
        else:
            data['lon'][i] = np.nan
    for i, lat in enumerate(data['latstr']):
        if lat.endswith('N'):
            data['lat'][i] = float(lat.replace('N', ''))
        elif lat.endswith('S'):
            data['lat'][i] = 0. - float(lat.replace('S', ''))
        else:
            data['lat'][i] = np.nan

    # Finding the date column
    fixed = datetime.strptime('00:00:00', '%H:%M:%S')
    datecol = [datetime.strptime(y, '%Y')
               + (datetime.strptime(u[-8:], '%j/%H%M') - fixed)
               for u, y in zip(data['datestr'], data['year'])]

    # Create the dataframe
    pddata = pd.DataFrame(data, index=datecol, columns=('lat', 'lon', 'qi'))
    # The datetime column is also needed
    pddata.insert(0, 'datetime', datecol)
    # Add the source column (assume argos, for now this looks at the
    # tara/argos subdir)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)
    # Remove the entries with qi!=3, and then remove the qi column
    pddata = pddata[pddata.qi == 3]
    del pddata['qi']

    return pddata, str(flist)


def read_sip(buoyname, buoycat, buoydir, out_ts):
    '''Read the text files of the SIP buoys (www.seaiceportal.de)'''

    bname = os.path.join(buoydir, 'www.seaiceportal.de',
                         '{}_{}_(proc|gps).csv'.format(buoyname, '\d{15}'))
    flist = glob(os.path.join(buoydir, 'www.seaiceportal.de', '*.csv'))

    for fl in flist:
        try:
            bfile = re.match(bname, fl).group(0)
            break
        except AttributeError:
            pass

    # Checking if the file creation should be skipped
    if out_ts is not None:
        try:
            pname = pathlib.Path(bfile)
        except:
            raise ValueError('File {} is not found'.format(bfile))
        txt_ts = datetime.fromtimestamp(pname.stat().st_mtime)
        if out_ts > txt_ts:
            return None, bfile

    # Messy way of reading a CSV file, but could not get genfromtxt to work -
    # issue with bytes/string type to convert "s away
    try:
        dates = []
        lats = []
        lons = []
        with open(bfile) as csvfile:
            csvdata = csv.reader(csvfile)
            # Skip the one line header
            next(csvdata)
            for row in csvdata:
                dates.append(row[0])
                lats.append(row[1])
                lons.append(row[2])
    except:
        raise IOError('File {} is not readable.'.format(bfile))

    # Finding the date column
    datecol = [datetime.strptime(dt, '%Y-%m-%dT%H:%M:%S') for dt in dates]
    # lons > 0 and filtering for valid lonlats
    # lons = [float(l.strip()) for l in lons]
    # lons = [(l + 360. if l < 0. else l) for l in lons]
    # lats = [float(l.strip()) for l in lats]
    lons_filt = []
    lats_filt = []
    datecol_filt = []
    for i in range(len(datecol)):
        # Use try/except loops to skip over empty entries like ''
        try:
            lat = float(lats[i].strip())
        except:
            break
        try:
            lon = float(lons[i].strip())
        except:
            break
        if lon < 0:
            lon = lon + 360.
        if not (lon >= 0.) and (lon < 360.):
            break
        if not (lat >= -90.) and (lat <= 90.):
            break
        lons_filt.append(lon)
        lats_filt.append(lat)
        datecol_filt.append(datecol[i])

    # Creating the data dictionary
    data = {}
    data['datecol'] = datecol_filt
    data['lat'] = lats_filt
    data['lon'] = lons_filt

    # Create the dataframe%
    pddata = pd.DataFrame(data, index=data['datecol'], columns=('lat', 'lon'))
    # Replace bad data
    pddata[pddata == 999.999] = np.nan
    # The datetime column is also needed
    pddata.insert(0, 'datetime', data['datecol'])
    # Add the source column (assume argos)
    source = [0] * len(datecol)
    pddata.insert(3, 'source', source)

    return pddata, bfile


def remove_invalid_lonlat(pddata):
    '''Remove any lines where lat or lon are Inf or Nan, or lon < 0 lon < 360,
    or lat < -90 lat > 90. If lon is in range -180 - 0, then this should be
    readjusted'''

    pddata = pddata[pddata.lat != np.inf]
    pddata = pddata[pddata.lon != np.inf]
    pddata = pddata[pddata.lat != np.nan] # This numpy nan does not seem
    pddata = pddata[pddata.lon != np.nan] # to work for pandas
    pddata = pddata[pddata.lon.notna()]
    pddata = pddata[pddata.lat.notna()]

    # Adjusting lons to range 0 - 360.
    ind = (pddata.lon < 0.) & (pddata.lon >= -180.)
    pddata.lon[ind] = pddata.lon[ind] + 360.

    pddata = pddata[pddata.lon >= 0.]
    pddata = pddata[pddata.lon < 360.]
    pddata = pddata[pddata.lat >= -90.]
    pddata = pddata[pddata.lat <= 90.]

    return pddata


def compute_distances_and_vels(pddata):
    '''Calculation of the distances and velocities between consecutive points
    in the pandas dataframe'''

#    if lat1 > 0.:
#        proj4_params = {'proj': 'stere',
#                        'lat_0': 90.,
#                        'lat_ts' : 70.,
#                        'lon_0': -45.0,
#                        'a': 6378273,
#                        'b': 6356889.44891}
#        area_def = AreaDefinition('plot', 'plot', 'plot', proj4_params,
#                                  1120, 760,
#                                  [-3850000, 3740000, 5840000, -5350000])
#    else:
#        area_def = {'proj': 'stere',
#                        'lat_0': -90.,
#                        'lat_ts' : -70.,
#                        'lon_0': 0.0,
#                        'a': 6378273,
#                        'b': 6356889.44891}
#        area_def = AreaDefinition('plot', 'plot', 'plot', proj4_params,
#                                  830, 790,
#                                  [-3950000, 3940000, 4340000, -3950000])
#
#    lons = (lon1, lon2)
#    lats = (lat1, lat2)
#    swath_def = SwathDefinition(lons=lons, lats=lats)

    # Initialise the distances and velocity arrays
    distances = [np.nan]
    velocities = [np.nan]
    sc1 = spherical.SCoordinate(pddata.lon[0], pddata.lat[0])
    t1 = pddata.datetime[0]
    for i in range(len(pddata.lat) - 1):
        sc2 = spherical.SCoordinate(pddata.lon[i+1], pddata.lat[i+1])
        t2 = pddata.datetime[i+1]
        dist = sc2.distance(sc1) * 1000.
        interval = (t2 - t1).total_seconds()
        vel = dist / interval
        distances.append(dist)
        velocities.append(vel)
        sc1 = sc2
        t1 = t2

    pddata.insert(4, 'distance', distances)
    pddata.insert(5, 'velocity', velocities)

    # Change the first item in the distance and velocity columns to the
    # means so they don't get thrown out
    pddata.distance[0] = pddata.distance.mean()
    pddata.velocity[0] = pddata.velocity.mean()

    return pddata


def sort_and_filter(pddata):
    '''Filter and sort the pandas dataframe of timestamp, lat, lon'''

    # NOTE: The original C code seems to remove all trajectories with
    # even 1 datapoint in the wrong hemisphere - such cases are raised
    # as a ValueError while setting metadata

    # Sort chronologically
    pddata = pddata.sort_index()

    # Remove any duplicate records (this seems to be done in the original
    # by removing any records which have the same timestamp)
    pddata = pddata.drop_duplicates(subset='datetime', keep=False)

    # Compute the distance moved and velocities between each consecutive
    # pair of points in the pandas dataframe
    pddata = compute_distances_and_vels(pddata)

    # Find data points at which buoys did not move, or buoys had velocities
    # which differ from the mean velocity by more than +/-3 std
    meanvel = pddata.velocity.mean()
    stdvel = pddata.velocity.std()
    # Buoys which didn't move
    errdist = pddata.distance == 0
    # Buoys with erroneous velocities.
    #errvel = np.logical_or(pddata.velocity > meanvel + (3 * stdvel),
    #                       pddata.velocity < meanvel - (3 * stdvel))
    errvel = pddata.velocity > meanvel + (3 * stdvel)

    # Finding the combined arrays of both these errors
    err = np.logical_or(errdist, errvel)
    # Calculation to remove both the array points which have the non-moving
    # buoy or the erroneous velocity
    errshift = err.copy()
    for i in range(len(pddata.lat) - 1):
        errshift[i] = errshift[i + 1]
    errshift[-1] = False
    totalerr = np.logical_or(err, errshift)

    # Use this error column to remove the bad data
    pddata.insert(6, 'totalerr', totalerr)
    pddata = pddata[pddata.totalerr != True]
    del pddata['totalerr']
    del pddata['distance']
    del pddata['velocity']

    # The original C code removes all records which have only one point in
    # the required time frame, but the time range is not selected yet -
    # moved to select_traj_2netcdf.py

    return pddata


def buoy_file_to_netcdf(buoyname, buoycat, buoydir, outdir, outname, force):

    # Finding the output name
    if outname is None:
        outname = '{}_{}.nc'.format(buoycat, buoyname)
    outfname = os.path.join(outdir, outname)
    # Checking if this file exists and what the datestamp is. If the force
    # option is given, this datestamp checking will not be used but the
    # file created anyway
    out_ts = None
    if not force and os.path.isfile(outfname):
        pname = pathlib.Path(outfname)
        out_ts = datetime.fromtimestamp(pname.stat().st_mtime)

    # Read the trajectory from the text file
    if buoycat == 'aari':
        pddata, orig_file = read_aari(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'antsid':
        pddata, orig_file = read_antsid(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'awi':
        pddata, orig_file = read_awi(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'bbb':
        pddata, orig_file = read_bbb(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'argos':
        pddata, orig_file = read_argos(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'crrel':
        pddata, orig_file = read_crrel(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'hudson':
        pddata, orig_file = read_hudson(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'iabp':
        # Use the new IABP reading
        pddata, orig_file = read_iabp_new(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'itp':
        pddata, orig_file = read_itp(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'pipers':
        pddata, orig_file = read_pipers(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'sams':
        pddata, orig_file = read_sams(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'sedna':
        pddata, orig_file = read_sedna(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'simba':
        pddata, orig_file = read_simba(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'tara':
        pddata, orig_file = read_tara(buoyname, buoycat, buoydir, out_ts)
    elif buoycat == 'sip':
        pddata, orig_file = read_sip(buoyname, buoycat, buoydir, out_ts)
    else:
        raise ValueError("Unrecognised buoy catalogue {}".format(buoycat))

    # Dictionary for network names
    networks = {'aari': 'AARI',
                'antsid': 'ANT',
                'awi': 'AWI',
                'bbb': 'BBB',
                'argos': 'ARG',
                'crrel': 'CRRE',
                'hudson': 'HUDS',
                'iabp': 'IABP',
                'itp': 'ITP',
                'pipers': 'PIP',
                'sams': 'SAMS',
                'sedna': 'SEDN',
                'simba': 'SIMB',
                'tara': 'TARA',
                'sip': 'SIP'}

    # If None is returned from the textfile reading, this is because it
    # was skipped
    if isinstance(pddata, pd.DataFrame):
        pass
    else:
        if pddata == None:
            print("{} already exists with a newer timestamp than the textfile, skipping.".format(outname))
            sys.exit(4)

    # Check for an empty dataframe and exit
    if pddata.empty:
        try:
            raise ValueError("WARNING: Dataframe for buoy {} catalogue {} "
                             "is empty, continuing".format(buoyname, buoycat))
        except:
            traceback.print_exc()
            sys.exit(3)

    # Filter the pandas dataframe for bad lats and lons
    pddata = remove_invalid_lonlat(pddata)

    # Sort and further filtering (call multiple times as sometimes when the
    # buoys get stuck they wibble a bit, so need to filter for distance 0
    # mutiple times)
    # TODO - this was trialled to try to get rid of stuck buoys by looping
    # but still issue with buoys "wiggling" when stuck. So need to check for
    # "wiggling" in small locus
    lpddata = len(pddata)
    stillsort = True
    while stillsort:
        pddata = sort_and_filter(pddata)
        nlpddata = len(pddata)
        if nlpddata == lpddata:
            stillsort = False
        lpddata = nlpddata

    if pddata.empty:
        try:
            raise ValueError("WARNING: Dataframe for buoy {} catalogue {} "
                             "is empty, continuing".format(buoyname, buoycat))
        except:
            traceback.print_exc()
            sys.exit(3)

    # Calculate seconds since 1970-01-01 and add this column
    reftime = datetime.strptime('19700101', '%Y%m%d')
    time = [(d - reftime).total_seconds() for d in pddata.datetime]
    pddata.insert(0, 'time', time)

    # Setting some metadata parameters
    cov_start_pt = pddata['datetime'].min()
    coverage_start = cov_start_pt.strftime('%Y-%m-%dT%H:%M:%SZ')
    cov_end_pt = pddata['datetime'].max()
    coverage_end = cov_end_pt.strftime('%Y-%m-%dT%H:%M:%SZ')
    min_lat = pddata['lat'].min()
    max_lat = pddata['lat'].max()
    min_lon = pddata['lon'].min()
    max_lon = pddata['lon'].max()
    if (max_lat > 0) and (min_lat > 0):
        area = 'Northern Hemisphere'
    elif (max_lat < 0) and (min_lat < 0):
        area = 'Southern Hemisphere'
    else:
        raise ValueError("Trajectory appears to cross the equator: min latitude {}, max latitude {}".format(min_lat, max_lat))

    # Removing the datetime column (we use time/seconds since 1970 instead)
    del pddata['datetime']

    # Modify the index
    pddata.set_index('time', inplace=True)

    # Create xarray Dataset from pandas dataframe
    xr = xarray.Dataset.from_dataframe(pddata)

    # Add variable attribute metadata
    xr['time'].attrs={'standard_name': 'time',
                      'long_name':'time',
                      'units': 'seconds since 1970-01-01 00:00:00'}
#    xr['datetime'].attrs={'standard_name': 'time',
#                          'long_name':'time'}
    xr['lat'].attrs={'standard_name': 'latitude',
                     'long_name':'latitude',
                     'units': 'degrees_north'}
    xr['lon'].attrs={'standard_name': 'longitude',
                     'long_name':'longitude',
                     'units': 'degrees_east'}
    xr['source'].attrs={'long_name':'source of data point',
                        'description': '0: argos, 1: gps, 2:iridium'}

    # Add global attribute metadata
    bmeta = netcdf_attr_dict(buoycat)
    xr.attrs={'Conventions':'CF-1.7',
              'title':'Buoy data from {}'.format(bmeta['catname']),
              'summary':'NetCDF conversion of buoy data from an existing text file, meant as an intermediate processing step', 
              'featureType': 'trajectory',
              'geospatial_lat_min' : min_lat,
              'geospatial_lat_max' : max_lat,
              'geospatial_lon_min' : min_lon,
              'geospatial_lon_max' : max_lon,
              'time_coverage_start' : coverage_start,
              'time_coverage_end' : coverage_end,
              'area' : area,
              'id' : buoyname,
              'network' : networks[buoycat],
              'history': '',
              'date_created': datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
              'project': ''.format(bmeta['catname']),
              'source': 'convert_buoy_txt2netcdf.py',
              'upstream_file': orig_file,
    }
#    print(xr)

    # Save to netCDF
    print("Writing file to {}".format(outfname))
    xr.to_netcdf(outfname, format="NETCDF4")

    # Add the trajectory variable
    name_strlen = len(buoyname)
    with Dataset(outfname, 'r+', format='NETCDF4_CLASSIC') as dataset:
        dataset.createDimension('name_strlen', name_strlen)
        traj = dataset.createVariable('trajectory', 'S1', 'name_strlen')
        traj[:] = stringtochar(np.array([buoyname], 'S{}'.format(name_strlen)))
        setattr(dataset['trajectory'], 'long_name', buoyname)
        setattr(dataset['trajectory'], 'cf_role', 'trajectory_id')


if __name__ == '__main__':

    args = parse_args()
    buoyname = args.buoyname
    buoycat = args.buoycat
    buoydir = args.buoydir
    outdir = args.outdir
    outname = args.outname
    force = args.force

    buoy_file_to_netcdf(buoyname, buoycat, buoydir, outdir, outname, force)
