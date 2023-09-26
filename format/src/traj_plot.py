'''Simple code to plot buoy trajectories'''

import sys
import os.path
from glob import glob
from datetime import datetime, timedelta
import argparse
from argparse import RawDescriptionHelpFormatter

import numpy as np
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.lines as mlines
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import cmocean

from grid_info import region_params, valid_regions


def parse_args():

    p = argparse.ArgumentParser("traj_plot",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('trajinp',
                   help="Either a NetCDF file with selected trajectories or a directory with the individual NetCDF trajectory files")
    p.add_argument('-r', '--region', required=False, choices=valid_regions,
                   default='nh',
                   help="Region to plot")
    p.add_argument('-s', '--startdate', required=False, default='19910101',
                   help="Start date of the time range to plot")
    p.add_argument('-e', '--enddate', required=False, default='20201231',
                   help="End date of the timerange to plot")
    p.add_argument('-o', '--output', required=False, default='.',
                   help="Either an output directory or a full filepath")
    p.add_argument('-l', '--label_traj', action='store_true', default=False,
                   help="Label trajectories with buoy ID")
    p.add_argument('-c', '--col_time', action='store_true', default=False,
                   help="Color trajectory with time (default is by program)")


    args = p.parse_args()

    return args


def val_get_netw_color(netw):
    netw = netw.lower()
    # NH
    if netw in ['aari', 'AAR']:
        color = 'red'
    elif netw in ['argos', 'argo', 'ARG']:
        color = 'orange'
    elif netw in ['bbb', 'BBB']:
        color = 'yellow'
    elif netw in ['crrel', 'crre', 'CRR']:
        color = 'lime'
    elif netw in ['dtu', 'DTU']:
        color = 'mint'
    elif netw in ['huds', 'HUD']:
        color = 'olive'
    elif netw in ['iabp', 'IAB']:
        color = 'green'
    elif netw in ['itp', 'ITP']:
        color = 'blue'
    elif netw in ['sams', 'SAM']:
        color = 'cyan'
    elif netw in ['sedna', 'sedn', 'SED']:
        color = 'purple'
    elif netw in ['simba', 'simb', 'SIM']:
        color = 'lavendar'
    elif netw in ['sip', 'SIP']:
        color = 'magenta'
    elif netw in ['tara', 'TAR']:
        color = 'teal'
    # SH
    elif netw in ['antsid', 'ants', 'ANT']:
        color = 'maroon'
    # Undefined
    else:
        color = 'black'
        print("WARNING: Unsupported network value <{}>. Default to {}."
              "".format(netw, color))

    return color


def traj_read(ncfile, sdate, edate, hemi):

    ncdata = {}

    # Reading from the NetCDF file
    with Dataset(ncfile, 'r') as dataset:

        try:
            ncdata['tc_start'] = datetime.strptime(dataset.start_date_and_time,
                                                   '%Y-%m-%d %H:%M:%S')
        except:
            ncdata['tc_start'] = datetime.strptime(dataset.time_coverage_start,
                                                   '%Y-%m-%dT%H:%M:%SZ')
        try:
            ncdata['tc_end'] = datetime.strptime(dataset.end_date_and_time,
                                                 '%Y-%m-%d %H:%M:%S')
        except:
            ncdata['tc_end'] = datetime.strptime(dataset.time_coverage_end,
                                                 '%Y-%m-%dT%H:%M:%SZ')

        # Check that the file covers the correct hemisphere and return None
        # if not
        if dataset.area == 'Northern Hemisphere':
            area = 'nh'
        elif dataset.area == 'Southern Hemisphere':
            area = 'sh'
        else:
            raise ValueError('Unrecognised region in ncfile')
        if not area == hemi:
            return None

        # Check that the file covers the required time period and return None
        # if not
        if ncdata['tc_start'] >= edate or ncdata['tc_end'] <= sdate:
            return None

        # Read the variable data
        ncdata['name'] = dataset.variables['trajectory'].long_name
        timevarsec = dataset.variables['time'][:]
        latvar = np.array(dataset.variables['lat'][:])
        lonvar = np.array(dataset.variables['lon'][:])

        # Other useful data
        ncdata['id'] = dataset.id
        ncdata['network'] = dataset.network

    # Convert seconds since 1970-01-01 to a python datetime
    timevar = np.array([datetime.fromtimestamp(t) for t in timevarsec])

    # Select only the relevent trajectory data
    select = np.logical_and(timevar >= sdate, timevar <= edate)
    ncdata['time'] = timevar[select]
    ncdata['lat'] = latvar[select]
    ncdata['lon'] = lonvar[select]

    # Create a timedelta variable
    ncdata['timedelta'] = np.array([(t - sdate).total_seconds()
                                    for t in ncdata['time']])
    ncdata['period'] = (edate - sdate).total_seconds() / 86400.
#    ncdata['lapsed'] = np.array([t / period for t in ncdata['timedelta']])
    # Use lapsed in terms of number of days
    # Seconds in a day = 86400
    ncdata['lapsed'] = np.array([t / 86400 for t in ncdata['timedelta']])

    # Other useful data
    ncdata['filename'] = ncfile

    return ncdata


def comb_trajfile_read(trajinp):

    # Reading data from the NetCDF file
    with Dataset(trajinp, 'r') as dataset:

        sdate = datetime.strptime(dataset.start_date_and_time,
                                  '%Y-%m-%d %H:%M:%S')
        edate = datetime.strptime(dataset.end_date_and_time,
                                  '%Y-%m-%d %H:%M:%S')
        period = (edate - sdate).total_seconds() / 86400.
        secs_to_1970 = (sdate
                        - datetime.strptime('19700101000000',
                                            '%Y%m%d%H%M%S')).total_seconds()

        slen = dataset.dimensions['station'].size

        ids = dataset['id'][:]
        networks = dataset['network'][:]
        times_sec = dataset['time'][:]
        lats = dataset['lat'][:]
        lons = dataset['lon'][:]

    # Compiling the data into a dictionary
    tjs = {}
    for i in range(slen):
        tj = {}
        tj['id'] = ''.join([b.decode() for b in ids[i, :]]).strip()
        tj['network'] = ''.join([n.decode() for n in networks[i, :]]).strip()

        # Create a timedelta variable
        tj['timedelta'] = np.array(times_sec[i, :])
        tj['lapsed'] = np.array([t / 86400 for t in tj['timedelta']])
        tj['period'] = period

        # Convert seconds since start date to seconds since 1970-01-01
        # and then to a python datetime
        tj['time'] = np.array([datetime.fromtimestamp(t + secs_to_1970)
                               for t in tj['timedelta']])

        tj['lat'] = lats[i, :]
        tj['lon'] = lons[i, :]
        tjs[i] = tj

    return tjs, sdate, edate


def traj_plot(trajinp, region, startdate, enddate, output, label_traj,
              col_time):


    # Finding the region parameters
    rp = region_params(region)
    if rp['lllat'] > 0.:
        hemi = 'nh'
    else:
        hemi = 'sh'
    if region.startswith('ease'):
        gridtype = 'ease'
    elif region.startswith('pol'):
        gridtype = 'polstere'
    else:
        raise ValueError("Unrecognised grid type for region {} - "
                         "expected to start with 'ease' or 'pol'".format(region))
    # Define grid based on region
    if gridtype == 'polstere':
        if hemi == 'nh':
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': 90.,
                                 'lat_ts' : 70.,
                                 'lon_0': -45.0,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.NorthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)
        else:
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': -90.,
                                 'lat_ts' : -70.,
                                 'lon_0': 0.,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.SouthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)
    elif gridtype == 'ease':
        if hemi == 'nh':
            plot_crs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                      central_latitude=90,
                                                      false_easting=0,
                                                      false_northing=0)
        else:
            plot_crs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                      central_latitude=-90,
                                                      false_easting=0,
                                                      false_northing=0)
    else:
        raise ValueError("Unrecognised region {}".format(region))

    pc = ccrs.PlateCarree()

    # Defining the plot corners
    llx, lly = plot_crs.transform_point(rp['lllon'], rp['lllat'], src_crs=pc)
    urx, ury = plot_crs.transform_point(rp['urlon'], rp['urlat'], src_crs=pc)

    # The projection keyword determines how the plot will look
    plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=plot_crs)

    # Setting the region
    ax.set_extent([llx, urx, lly, ury], crs=plot_crs)

    # Converting the start and end dates to Python datetimes (using the
    # entire day of the second date)
    sdate = datetime.strptime(startdate, '%Y%m%d')
    edate = datetime.strptime(enddate + '235959', '%Y%m%d%H%M%S')

    # Reading the trajectories
    if os.path.isdir(trajinp):

        index = 0
        tjs = {}
        flist = glob(os.path.join(trajinp, '*.nc'))
        #print("flist = ", flist)
        for fl in flist:
            print("Checking file: {}".format(fl))

            #try:
            trajdata = traj_read(fl, sdate, edate, hemi)
            #except:
            #    print("Problem reading trajectory from {}".format(fl))
            #    continue
            if trajdata is not None:
                tjs[index] = trajdata
                index += 1
                print("Data found in file...")

    elif os.path.isfile(trajinp):
        # Check if this is a single trajectory or a buoy collection
        single = False
        with Dataset(trajinp, 'r') as dataset:
            if 'source' in dataset.__dict__.keys():
                if dataset.source == "convert_buoy_txt2netcdf.py":
                    single = True
                    sdate = datetime.strptime(dataset.time_coverage_start,
                                              '%Y-%m-%dT%H:%M:%SZ')
                    edate = datetime.strptime(dataset.time_coverage_end,
                                              '%Y-%m-%dT%H:%M:%SZ')
        if single:
            tjs = {}
            tjs[0] = traj_read(trajinp, sdate, edate, hemi)
            #sdate = tjs[0]['tc_start']
            #edate = tjs[0]['tc_end']
        else:
            tjs, sdate, edate = comb_trajfile_read(trajinp)

    else:
        raise ValueError("{} is not a valid file or directory".format(trajinp))

    # Checking if there is valid data
    if len(tjs) == 0:
        print("No valid trajectory data to plot, exiting.")
        sys.exit(0)

    # Setting the colourmap
    cmap_time = cm.get_cmap('Spectral_r')
    period = tjs[0]['period']
    cb = ax.imshow(np.array([[0.0, 0.0], [0.0, 0.0]]), cmap=cmap_time,
                   vmin=0, vmax=period)

    # Plotting the trajectories
    netw_keys = []
    for i, tj in enumerate(tjs):

        # Default is to use network color
        netw = tjs[i]['network']
        col = val_get_netw_color(netw)

# PLOTTING VECTORS
#        for j in range(tjs[i]['lon'].size):
#
#            # If colouring by datetime, set this here per point
#            if col_time:
#                col = cmap_time(tjs[i]['lapsed'][j] / period)
#
#            x1, y1 = plot_crs.transform_point(tjs[i]['lon'][j],
#                                              tjs[i]['lat'][j], src_crs=pc)
#
#            # Label
#            if label_traj and j == 0:
#                label = '{}-{}'.format(tjs[i]['id'], tjs[i]['network'])
#                if ((x1 > llx) and (x1 < urx) and (y1 > lly) and (y1 < ury)):
#                    ax.text(x1, y1, label, fontsize='small')
#
#            # Plot the drift
#            if j > 0:
#                ax.plot([x0, x1], [y0, y1], color=col)
#                if not netw in netw_keys:
#                    netw_keys.append(netw)
#            x0 = x1
#            y0 = y1

        for j in range(tjs[i]['lon'].size):

            # If colouring by datetime, set this here per point
            if col_time:
                col = cmap_time(tjs[i]['lapsed'][j] / period)

            x0, y0 = plot_crs.transform_point(tjs[i]['lon'][j],
                                              tjs[i]['lat'][j], src_crs=pc)

            # Label
            if label_traj and j == 0:
                label = '{}-{}'.format(tjs[i]['id'], tjs[i]['network'])
                #if ((x1 > llx) and (x1 < urx) and (y1 > lly) and (y1 < ury)):
                ax.text(x0, y0, label, fontsize='small')

            # Plot the drift
            ax.plot(x0, y0, '-o', markersize=2, color=col)
            if not netw in netw_keys:
                netw_keys.append(netw)


    if col_time:
        # For the lines coloured by time, add the colourbar
        if period <= 10:
            cmap_lvl = np.arange(0, period)
        elif period > 10 and period <= 100:
            cmap_lvl = np.arange(0, period, 10)
        elif period > 100 and period <= 1000:
            cmap_lvl = np.arange(0, period, 100)
        else:
            cmap_lvl = np.arange(0, period, 1000)
        plt.colorbar(cb, ticks=cmap_lvl,
                     orientation='horizontal', pad=0.05, shrink=0.4
        ).set_label('Elapsed time (days)')
    else:
        # For the lines coloured by network, add the key
        handles = []
        for n in netw_keys:
            col = val_get_netw_color(n)
            handles.append(mlines.Line2D([], [], color=col, label=n))
        ax.legend(handles=handles)

    # Coastlines
#    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land',
#                                                '50m', edgecolor='black',
#                                                facecolor='white'))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land',
                                                '50m', edgecolor='grey',
                                                facecolor='grey'))
    ice_shelves = cfeature.NaturalEarthFeature(
        category='physical',
        name='antarctic_ice_shelves_polys',
        scale='10m',
        facecolor='grey',
        edgecolor='grey')
    ax.add_feature(ice_shelves

)


    # Title
    title = '{}-{} {}'.format(datetime.strftime(sdate, '%Y%m%d'),
                              datetime.strftime(edate, '%Y%m%d'),
                              rp['long_name'])
    plt.title(title, fontsize='small')

    # If the path given is only a directory, add a filename. Otherwise
    # we assume it is a valid filepath
    if os.path.isdir(output):
        name = 'trajectories.png'
        output = os.path.join(output,name)

    plt.savefig(output,bbox_inches='tight')
#    plt.show()
    plt.close()
    print("Figure is in {}".format(output))


def main():

    args = parse_args()
    trajinp = args.trajinp
    region = args.region
    startdate = args.startdate
    enddate = args.enddate
    output = args.output
    label_traj = args.label_traj
    col_time = args.col_time

    traj_plot(trajinp, region, startdate, enddate, output, label_traj, col_time)


if __name__ == '__main__':

    main()
