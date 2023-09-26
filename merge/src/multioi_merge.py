import os
import sys
import math
import argparse
from argparse import RawDescriptionHelpFormatter
from glob import glob
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
from numpy import ma
import cartopy
import cartopy.crs as ccrs
from pyresample import geometry, kd_tree, spherical
from merge_periods import merge_season, merge_frac

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')

from read_icedrift_product_file import read_drift
sys.path.append(os.path.join(here, '../../ice-tracking/src'))
from write_icedrift_product_file import write_icedrift
from add_uncert_dict import ud
from icedrift_flags import id_flags

max_drift_mask = {'nh': os.path.join(datadir, 'paramfiles/max_ice_drift_mask_nh.nc'),
                  'sh': os.path.join(datadir, 'paramfiles/max_ice_drift_mask_sh.nc')}

max_allowed_uncert = 50.
lat_max_discard = 87.5

min_sat_pix_ratio = {'nh': 0.4, 'sh': 0.4}
pole_hole_mask_lat = 86.


def parse_args():

    valid_grids = ['nh-polstere-625', 'sh-polstere-625', 'nh-ease2-750',
                   'sh-ease2-750']
    valid_levels = ['lev2', 'lev3']
    valid_tbname = ['bt', 'tb']

    p = argparse.ArgumentParser("multioi_merge",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--indirs', required=True,
                   help="Comma-separated list of input directories")
    p.add_argument('-o', '--output', required=True,
                   help="Output directory")
    p.add_argument('-s', '--startdate', required=True,
                   help="Start date of the files to combine, in either YYmmdd or YYmmddHH format. If no hour is specified, defaults to 12H.")
    p.add_argument('-e', '--enddate', required=True,
                   help="End date of the files to combine, in either YYmmdd or YYmmddHH format. If no hour is specified, defaults to 12H.")
    p.add_argument('-I', '--sensors', required=True,
                   help="List of sensors to combine")
    p.add_argument('-g', '--grid', required=False, default='nh_polstere-625',
                   choices=valid_grids,
                   help="Grids to combine, valid grids are {}, default grid 'nh_polstere-625'".format(valid_grids))
#    p.add_argument('--reproc', action='store_true',
#                   help='expect input files in YYYY/MM/DD dirs, and place results in YYYY/MM/DD ones.',)
    p.add_argument('--reproc', action='store_true',
                   help='Expect input files in inst/YYYY/MM/DD dirs, and place results in YYYY/MM ones.',)
    p.add_argument('-x', '--windname', required=False, default=None,
                   help='The name of the wind drift sub-directory',)
    p.add_argument('-l', '--level', required=False, default='lev2',
                   choices=valid_levels,
                   help="Data level, valid choices are {}, default level 'lev3'")
    p.add_argument('-tb', '--tbname', required=False, default='bt',
                   choices=valid_tbname,
                   help="Name of brightness temperature, default it 'bt', valid choices are {}".format(valid_tbname))
    p.add_argument('-c', '--chan', required=False, default=None,
                   help="Name of passive microwave channel, default is '37' for AMSR2 and '90' for SSMI")
    p.add_argument('-f', '--fixedgsf', action='store_true', default=False,
                   help="Fix the uncertainties of the merge product to a static value where they are determined by gapfilling (rather than gapfill the uncertainties)")
    args = p.parse_args()

    return args


def icedrift_fname(sensor, level, sdate, edate, grid, tbname='bt', chan=None):

    sdatef = datetime.strftime(sdate, '%Y%m%d%H')
    edatef = datetime.strftime(edate, '%Y%m%d%H')

    # Default channels for AMSR2 and SSMIS
    if chan is None:
        if sensor.startswith('amsr'):
            chan = '37'
        elif sensor.startswith('ssmi'):
            chan='90'

    if sensor.startswith('amsr') or sensor.startswith('ssmi'):
        fname = 'icedrift_{}_{}{}v-Lap+{}{}h-Lap_simplex_{}_{}_{}w_{}w.nc'.format(sensor, tbname, chan, tbname, chan, level, grid, sdatef,
                  edatef)
    elif sensor.startswith('ascat'):
        fname = 'icedrift_{}_sigma0-Lap_simplex_{}_{}_{}w_{}w.nc'.format(sensor, level, grid, sdatef, edatef)
    elif sensor.startswith('wind'):
        fname = 'icedrift_ecmwf-era5_none_{}_{}_{}_{}w_{}w.nc'.format(sensor, level,
                                                          grid, sdatef, edatef)
    else:
        raise ValueError("Unrecognised sensor {}".format(sensor))

    return fname


def compute_distance(lat0, lon0, lat1, lon1):

    rearth = 6371.
    deg_to_rad = math.pi / 180.
    lat0 *= deg_to_rad
    lon0 *= deg_to_rad
    lat1 *= deg_to_rad
    lon1 *= deg_to_rad
    dlat = 0.5 * (lat1 - lat0)
    dlon = 0.5 * (lon1 - lon0)
    a = (math.sin(dlat) * math.sin(dlat)) + (
        math.cos(lat0) * math.cos(lat1) * math.sin(dlon) * math.sin(dlon))
    dist = rearth * 2. * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return dist


def compute_directiontonorth(lat0, lon0, lat1, lon1):

    deg_to_rad = math.pi / 180.
    lat0 *= deg_to_rad
    lon0 *= deg_to_rad
    lat1 *= deg_to_rad
    lon1 *= deg_to_rad
    dlat = 0.5 * (lat1 - lat0)
    dlon = lon1 - lon0
    dirn = math.atan2(math.sin(dlon) * math.cos(lat1),
                      (math.cos(lat0) * math.sin(lat1))
                      - (math.sin(lat0) * math.cos(lat1) * math.cos(dlon)))
    dirn += math.pi
    dirn /= deg_to_rad

    return dirn


def check_sat_data(sensors, driftdata, hemi):
    '''Check that there are enough pixels in the satellite data fields and
    throw these out if not'''

    ny, nx = driftdata[sensors[0]]['area_def'].shape

    # Make a mask for the pole hole
    polehole_mask = np.zeros((ny, nx), dtype=np.int8)
    polehole_mask = driftdata[sensors[0]]['lat'][:] >= pole_hole_mask_lat

    remlist = []
    for sens in sensors:
        if sens != 'wind':
            driftgrid = np.zeros((ny, nx), dtype=np.int8)
            icegrid = np.zeros((ny, nx), dtype=np.int8)
            # Find where there is ice drift
            dx = driftdata[sens]['driftX'][:]
            driftgrid[~dx.mask] = 1
            # Not considering the pole hole in the NH
            if hemi == 'nh':
                driftgrid[polehole_mask] = 0
            driftnum = np.sum(driftgrid)
            # Find where there is ice
            flag = driftdata[sens]['flag'][:]
            icegrid[np.logical_or(flag <= 0, flag > 5)] = 1
            # Not considering the pole hole in the NH
            if hemi == 'nh':
                icegrid[polehole_mask] = 0
            icenum = np.sum(icegrid)
            # Calculate the ratio of drift vectors to ice pixels
            ratio = driftnum / icenum
            # Remove this sensor if there are not enough drift vectors
            # compared to the number of ice pixels
            if ratio < min_sat_pix_ratio[hemi]:
                remlist.append(sens)
                print("WARNING: Removing data from sensor {} as not enough drift vectors".format(sens))

    for rs in remlist:
        sensors.pop(sensors.index(rs))

    return sensors, driftdata


def mask_wind_to_sat(sensors, driftdata):
    '''Mask the wind drift data to the shape of the valid satellite data'''

    valid_data = [id_flags['icedrift_ok'],
                  id_flags['icedrift_correct_by_neighbours'],
                  id_flags['icedrift_smaller_pattern']]

    # Find a list of the satellite sensors only
    satsensors = [x for x in sensors if not x.startswith('wind')]

    # Initialise a datacube
    refshape = driftdata[satsensors[0]]['flag'][:, :].shape
    refflags = np.zeros((len(satsensors), refshape[0], refshape[1]),
                        dtype=np.int8)

    # Read the flags from all the reference single sensor icedrift files
    # into the datacube
    for i, sens in enumerate(satsensors):
        refflags[i, :, :] = driftdata[sens]['flag'][:, :]

    # The land should be the same in every image
    land = np.zeros((refshape[0], refshape[1]), dtype=np.int8)
    land[refflags[0, :, :] == id_flags['icedrift_center_over_land']] = 1

    # For noice, we want the area where none of the single sensor files
    # have valid data
    driftveccube = np.zeros((len(satsensors), refshape[0], refshape[1]),
                            dtype=np.int8)
    for val in valid_data:
        driftveccube[:, :, :][refflags == val] = 1
    driftvec = np.any(driftveccube, axis=0)

    # Apply the masking to the wind data
    windmask = driftdata['wind']['driftX'].mask
    newmask = np.logical_or(windmask, ~driftvec)
    driftdata['wind']['driftX'][newmask] = np.nan
    driftdata['wind']['driftX'].mask = newmask
    driftdata['wind']['driftY'][newmask] = np.nan
    driftdata['wind']['driftY'].mask = newmask
    driftdata['wind']['fsX_12utc'][newmask] = np.nan
    driftdata['wind']['fsX_12utc'].mask = newmask
    driftdata['wind']['fsY_12utc'][newmask] = np.nan
    driftdata['wind']['fsY_12utc'].mask = newmask

    return driftdata, driftvec, windmask, newmask


def find_masking(sensors, driftdata, sizearr, nofilllist):
    '''Find the masking for the datacube'''

    masking = {}
    nofillmask = {}
    for i, sensor in enumerate(sensors):
        masking[sensor] = np.ma.zeros(sizearr)
        nofillmask[sensor] = np.ma.zeros(sizearr)
        # Where drift or uncertainty is masked
        masking[sensor][:][driftdata[sensor]['driftX'].mask] = 1
        masking[sensor][:][driftdata[sensor]['driftY'].mask] = 1
        masking[sensor][:][driftdata[sensor]['fsX_12utc'].mask] = 1
        # Where drift or uncertainty is NaN
        masking[sensor][:][driftdata[sensor]['driftX'] == np.nan] = 1
        masking[sensor][:][driftdata[sensor]['driftY'] == np.nan] = 1
        masking[sensor][:][driftdata[sensor]['fsX_12utc'] == np.nan] = 1
        # Where uncertainty is very high
        masking[sensor][:][driftdata[sensor]['fsX_12utc'] >= max_allowed_uncert] = 1
        # Where latitude is high and it is not a nominal vector - UNLESS
        # we only have a single wind file, then we don't want to do this
        if len(sensors) == 1 and 'wind' in sensors:
            print("SKIPPING HIGH LAT MASKING FOR SINGLE WIND FIELD")
        else:
            masking[sensor][:][np.logical_and(
                driftdata[sensor]['lat'] >= lat_max_discard,
                driftdata[sensor]['flag'] != id_flags['icedrift_ok'])] = 1
        # Where there is land/coast/border
        for nof in nofilllist:
            nofillmask[sensor][:][driftdata[sensor]['flag'] == nof] = 1
        masking[sensor][:][nofillmask[sensor][:] == 1] = 1
            #masking[sensor][:][driftdata[sensor]['flag'] == nof] = 1

    return masking, nofillmask


def multioi_merge(indirs, output, startdate, enddate, sensors, grid, level,
                  reproc, windname=None, tbname='bt', chan=None, wghts=None,
                  osi405yr=False, fixedgsf=False):

    fillvalf = -1.0e-10

    hemi = grid[0:2]

    # If no hours are specified, this defaults to 12h
    if len(startdate) == 8:
        startdate = '{}12'.format(startdate)
    elif len(startdate) != 10:
        raise ValueError("startdate must be in format YYmmdd or YYmmddHH, "
                         "not {}".format(startdate))
    if len(enddate) == 8:
        enddate = '{}12'.format(enddate)
    elif len(enddate) != 10:
        raise ValueError("enddate must be in format YYmmdd or YYmmddHH, "
                         "not {}".format(enddate))
    sdate = datetime.strptime(startdate, '%Y%m%d%H')
    edate = datetime.strptime(enddate, '%Y%m%d%H')

    # Find list of files
    sensors = sensors.split(',')
    fnames = {}
    for sensor in sensors:
        fnames[sensor] = icedrift_fname(sensor, level, sdate, edate, grid,
                                        tbname=tbname, chan=chan)

    fpaths = {}
    if reproc:
        inpdirs = []
        for sensor in sensors:
            sens_dirs = [os.path.join(idir, sensor, enddate[0:4],
                                     enddate[4:6],
                                     enddate[6:8]) for idir in indirs.split(',')]
            inpdirs.extend(sens_dirs)
        if windname:
            inpdirs.extend([os.path.join(idir, windname, enddate[0:4],
                                         enddate[4:6],
                                         enddate[6:8]) for idir in indirs.split(',')])
    else:
        inpdirs = indirs.split(',')

    for inpdir in inpdirs:
        for sensor in sensors:
            fpath = os.path.join(inpdir, fnames[sensor])
            if os.path.isfile(fpath):
                fpaths[sensor] = fpath

    # Check that some files to merge were found
    remlist = []
    for sensor in sensors:
        if not sensor in fpaths:
            print("WARNING: No data found for sensor {}".format(sensor))
            remlist.append(sensor)
    for rs in remlist:
        sensors.pop(sensors.index(rs))
    if len(fpaths) < 1:
        raise ValueError("No data found to merge.")
    else:
        print("Processing files...", fpaths)

    # Read in the files
    driftdata = {}
    remlist = []
    for sensor in sensors:
        try:
            driftdata[sensor] = read_drift(fpaths[sensor], short=True)
        except:
            print("WARNING: File {} could not be read in".format(fpaths[sensor]))
            remlist.append(sensor)
    for rs in remlist:
        sensors.pop(sensors.index(rs))

    # Check for no data and exit if so
    if len(driftdata) == 0:
        raise ValueError("No files to merge could be successfully read")

    # Check the number of pixels in the satellite data fields. If there
    # is too little data, files are removed. There should still be a wind
    # file which will then be used.
    sensors, driftdata = check_sat_data(sensors, driftdata, hemi)

    # Initialise
    numfiles = len(sensors)
    sizearr = driftdata[sensors[0]]['driftX'].shape
    id_land = id_flags['icedrift_center_over_land']
    id_ccoe = id_flags['icedrift_closeto_coast_or_edge']
    id_noice = id_flags['icedrift_noice']
    id_outbord = id_flags['icedrift_outside_imgborder']
    id_closebord = id_flags['icedrift_closeto_imgborder']
    id_maxdriftmask = id_flags['icedrift_maxdriftmask']
    nofilllist = [id_land, id_ccoe, id_noice, id_outbord, id_closebord,
                  id_maxdriftmask]

    # Mask where data is not going to contribute to the merged product
    masking, nofillmask = find_masking(sensors, driftdata, sizearr, nofilllist)

    # If there are sensors which shall not contribute due to being entirely
    # masked, then they are removed from the list
    remlist = []
    for i, sensor in enumerate(sensors):
        if np.all(masking[sensor]):
           remlist.append(sensor)
    for rs in remlist:
        sensors.pop(sensors.index(rs))
    # Rerun the masking code
    masking, nofillmask = find_masking(sensors, driftdata, sizearr, nofilllist)

    # Now check that we have not removed ALL the data (which can happen in
    # winter), and if so, it is necessary to reinstate the wind drift
    # with it's base uncertainty value
    if len(sensors) == 0:
        print("WARNING: No good data found for any sensor, reinstating wind")
        sensors = ['wind']

        # We use a high rather than very high uncertainty value
        windvh = ud['cdr']['veryhigh']['wind'][0]
        windh = ud['cdr']['high']['wind'][0]
        driftdata['wind']['fsX_12utc'][driftdata['wind']['fsX_12utc'] >= windvh] = windh

        ## Finding the uncertainty value for the wind in this case
        #season, _ = merge_season(edate, hemi)
        ## Flag 0 is used, this is anyway taken as the uncertainty from
        ## all the flags
        #w_uncert = ud['cdr'][hemi]['wind'][season][0]
        #driftdata['wind']['fsX_12utc'][driftdata['wind']['fsX_12utc'] > w_uncert] = w_uncert

        # Rerun the code to determine the masking
        masking, nofillmask = find_masking(sensors, driftdata, sizearr,
                                           nofilllist)
        if np.all(masking['wind']):
            raise ValueError("No good data found for any sensor, failing.")
    print("Sensors used in the merging: {}".format(sensors))

    # If there is both a wind drift file and at least one satellite
    # drift file, then it is necessary to trim the wind data to the
    # shape of the satellite data.
    if 'wind' in sensors and len(sensors) > 1:
        driftdata, driftvec, windmask, newmask = mask_wind_to_sat(sensors, driftdata)
        # Rerun the masking
        masking, nofillmask = find_masking(sensors, driftdata, sizearr,
                                           nofilllist)

    # Initialise
    stacksize = (len(sensors), sizearr[0], sizearr[1])
    funcwght = np.ma.zeros(stacksize)
    funcXs = np.ma.zeros(stacksize)
    funcYs = np.ma.zeros(stacksize)
    flags = np.ma.zeros(stacksize)
    nofills = np.ma.zeros(stacksize)

    # Apply the masking
    for i, sensor in enumerate(sensors):

        # Apply the masking to the drift and uncertainty and flags
        driftdata[sensor]['driftX'].mask = masking[sensor][:]
        driftdata[sensor]['driftY'].mask = masking[sensor][:]
        driftdata[sensor]['fsX_12utc'].mask = masking[sensor][:]

        # Only need to read in the lats and lons once
        if i == 0:
            lats = driftdata[sensor]['lat']
            lons = driftdata[sensor]['lon']
            area_def = driftdata[sensor]['area_def']
            fillval = driftdata[sensor]['driftX_fv']

        # For the uncertainties, the X uncertainties corrected to 12UTC are
        # used (the X and Y uncertainties are considered to be the same)
        funcwght[i, :, :] = 1 / (driftdata[sensor]['fsX_12utc']
                                 * driftdata[sensor]['fsX_12utc'])

        funcXs[i, :, :] = driftdata[sensor]['driftX'] * funcwght[i, :, :]
        funcYs[i, :, :] = driftdata[sensor]['driftY'] * funcwght[i, :, :]

        # And also "masking" any flags where the sensor has been taken out
        # of the final product and this does not fall into a "no fill"
        # region (land, coast, etc). This is dine by setting the flag to a
        # very high value since the minimum value is later taken for the flag
        flagfld = driftdata[sensor]['flag'].copy()
        flgmask = np.logical_and(masking[sensor][:] == 1,
                                 nofillmask[sensor][:] == 0)
        flagfld[flgmask] = 100
        flags[i, :, :] = flagfld
        nofills[i, :, :] = flags[i, :, :] == nofilllist[0]
        for j in range(len(nofilllist)-1):
            nofills[i, :, :] = np.logical_or(nofills[i, :, :],
                                             flags[i, :, :] == nofilllist[j+1])

    sumfuncw = np.nansum(funcwght, axis=0)
    sumfuncX = np.nansum(funcXs, axis=0)
    sumfuncY = np.nansum(funcYs, axis=0)

    avgx = sumfuncX / sumfuncw
    avgy = sumfuncY / sumfuncw

#    # Count the number of sensors (and take case of a mask if there
#    # is one). It is assumed that masks/NaNs match between the driftX,
#    # driftY and uncertainty arrays
#    nnmask = np.invert(np.isfinite(funcXs)) | funcXs.mask
#    sarr = np.ma.MaskedArray(funcXs, mask=nnmask)
#    sumsens = sarr.count(axis=0)

    # sig = alpha * sqrt( 1 / (SUM_i (1 / sig_i^2)))
    alpha = 1.5
    uncert = alpha * np.sqrt(1 / (sumfuncw))

    nofillpix = np.any(nofills, axis=0)

    # Gap interpolation
    # NOTE: The pyresample "sigma" sig_gf is actually related to the gaussian
    # rms c by sig_gf^2 = 2 * c^2. In the ATBD the standard deviation of the
    # gaussian is 200km.
    sig_gf = 200000. * math.sqrt(2)
    # In the ATBD, the region of interest is +/-4 pix
    roi_gf = 300000. # ease2, each pixel=75km

    gapfillfields={}
    # Looping over each parameter
    fields = [avgx, avgy, uncert]
    for i in range(len(fields)):

        # Specific flag fields
        datapix = ~fields[i].mask
        # Mask the pole hole
        datapix = np.logical_and(datapix, lats < 89.5)

        # Swath coverage information
        datalons = lons[datapix]
        datalats = lats[datapix]
        coverage = geometry.SwathDefinition(lons=datalons, lats=datalats)

        # Extrapolation
        datacover = fields[i][datapix]
        gapfillfields[i] = kd_tree.resample_gauss(coverage, datacover,
                                                  area_def,
                                                  radius_of_influence=roi_gf,
                                                  sigmas=sig_gf)
        # If we want to fix the status flag for the gapfilled uncertainty
        # then wipe out the gapfilled field here
        if fixedgsf and i == 2:
            gapfillfields[i].fill(ud['cdr'][hemi]['multi-oi'][16])
        gapfillfields[i][datapix] = fields[i][datapix]
        gapfillfields[i] = ma.array(gapfillfields[i])
        gapfillfields[i].mask = nofillpix

    avgx_gf = gapfillfields[0]
    avgy_gf = gapfillfields[1]
    uncert_gf = gapfillfields[2]

    # Flag combination
    # Where the product is successfully merged, use flag icedrift_ok
    # Where the product is interpolated, use flag icedrift_oiinterp
    # Elsewhere, use the minimum flag of the combined maps
    flag = np.nanmin(flags, axis=0)
    flag[np.logical_and(avgx.mask, ~avgx_gf.mask)] = id_flags['icedrift_oiinterp']

    # Changing the icedrift_ok flag to reflect the seasonal filling with wind
    # drift
    if 'wind' in sensors:
        if len(sensors) == 1:
            flag[flag == id_flags['icedrift_ok']] = id_flags['icedrift_fillbywind']
            # The uncertainty should then match the measured uncertainty of
            # the wind drift for that season
            season, _ = merge_season(edate, hemi)
            # Flag 0 is used, this is anyway taken as the uncertainty from
            # all the flags
            wfill_uncert = ud['cdr'][hemi]['wind'][season][0] * alpha
            uncert_gf[flag == id_flags['icedrift_fillbywind']] = wfill_uncert
            uncert_gf[flag == id_flags['icedrift_windgapfill']] = wfill_uncert

        elif len(sensors) > 1:
            good_icedrift = np.zeros_like(flag)
            for gid in [id_flags['icedrift_ok'],
                        id_flags['icedrift_correct_by_neighbours'],
                        id_flags['icedrift_smaller_pattern']]:
                good_icedrift = np.logical_or(good_icedrift, flag == gid)
            flag[good_icedrift] = id_flags['icedrift_blendwithwind']

    # Masking any drift outside the maximum drift mask, as we do not want
    # coastal and fjord drift
    maxdds = Dataset(max_drift_mask[hemi], 'r')
    maxdriftmask = maxdds['maxdriftmask'][:]
    extramsk = np.logical_and(avgx_gf.mask == 0, maxdriftmask == 0)
    print("EXTRA {} PIXELS MASKED BY MAX ICE DRIFT MASK".format(np.sum(extramsk)))
    avgx_gf.mask[extramsk] = 1
    avgy_gf.mask[extramsk] = 1
    uncert_gf.mask[extramsk] = 1
    flag[extramsk] = id_flags['icedrift_maxdriftmask']

    # And finally, check where the interpolation has not worked since the
    # pixels are too far away. If there are pixels where dX = dY = 0 and
    # the pixel is interpolated, fill with wind and flag this.
    # NOTE: The wind file will have been trimmed to match the satellite
    # data so needs to be read in again
    interpfail = np.logical_and(avgx_gf == 0, avgy_gf == 0)
    interpfail = np.logical_and(interpfail, flag == id_flags['icedrift_oiinterp'])
    print("{} INTERP FAILED PIXELS SET TO WIND DRIFT".format(np.sum(interpfail)))

    if np.sum(interpfail) > 0:
        try:
            driftdata['windorig'] = read_drift(fpaths['wind'], short=True)
        except:
            print("WARNING: File {} could not be read in".format(fpaths[sensor]))        
        avgx_gf[interpfail] = driftdata['windorig']['driftX'][interpfail]
        avgy_gf[interpfail] = driftdata['windorig']['driftY'][interpfail]
        # Finding the uncertainty value for the wind in this case
        season, _ = merge_season(edate, hemi)
        # Flag 0 is used, this is anyway taken as the uncertainty from
        # all the flags
        wfill_uncert = ud['cdr'][hemi]['wind'][season][0] * alpha
        uncert_gf[interpfail] = wfill_uncert
        # The flag is either fillbywind or the gapfilled wind parameter file
        # flag
        gapfl = np.zeros_like(driftdata['wind']['flag'][:])
        gapfl = driftdata['wind']['flag'][:] == id_flags['icedrift_windgapfill']
        interpf_gfp = np.logical_and(interpfail, gapfl)
        interpf_nogfp = np.logical_and(interpfail, gapfl == 0)
        flag[interpf_gfp] = id_flags['icedrift_windgapfill']
        flag[interpf_nogfp] = id_flags['icedrift_fillbywind']

    # Time variables
    datamask = avgx_gf.mask
    t0 = np.full(sizearr, fillvalf)
    t0 = ma.array(t0, mask=datamask)
    t0[~datamask] = 0
    t1 = np.full(sizearr, fillvalf)
    t1 = ma.array(t1, mask=datamask)
    t1[~datamask] = 0

    # Beginning and ending lats and lons
    pc = ccrs.PlateCarree()
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
        x0, y0 = driftdata[sensors[0]]['data_ccrs'].transform_point(lon0, lat0, pc)
        x1 = x0 + (avgx_gf[ix, iy] * 1000.)
        y1 = y0 - (avgy_gf[ix, iy] * 1000.)
        lon1, lat1 = pc.transform_point(x1, y1, driftdata[sensors[0]]['data_ccrs'])
        olate[ix, iy] = lat1
        olone[ix, iy] = lon1

        # Distance and direction calculation
       # TODO - Can these be done with python functions from e.g.
        # pyresample or cartopy?
        dist[ix, iy] = compute_distance(lat0, lon0, lat1, lon1)
        dirn[ix, iy] = compute_directiontonorth(lat0, lon0, lat1, lon1)

    # If reproc, add the inst/date subdirs
    if reproc:
        output = os.path.join(output, enddate[0:4], enddate[4:6], enddate[6:8])
        if not os.path.exists(output):
            os.makedirs(output)

    # Writing the final file out
    blank0 = np.zeros(sizearr)
    blank1 = np.full(sizearr, 1)
    blankfill = np.full(sizearr, -32767)
    projstr = str(driftdata[sensors[0]]['data_crs'])
    projstr.replace('_PROJ4Projection(', '')
    projstr.replace(')', '')
    write_icedrift(driftdata[sensors[0]]['area_def'],
                          None,
                          None,
                          None,
                          sdate,
                          edate,
                          grid,
                          projstr,
                          output,
                          avgx_gf,
                          avgy_gf,
                          t0,
                          t1,
                          flag,
                          uncert_gf,
                          uncert_gf,
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
                          method='merge')


def main():

    args = parse_args()
    indirs = args.indirs
    output = args.output
    startdate = args.startdate
    enddate = args.enddate
    sensors = args.sensors
    grid = args.grid
    level = args.level
    tbname = args.tbname
    chan = args.chan
    reproc = args.reproc
    windname = args.windname
    fixedgsf = args.fixedgsf

    multioi_merge(indirs, output, startdate, enddate, sensors, grid, level,
                  reproc, windname=windname, tbname=tbname, chan=chan,
                  fixedgsf=fixedgsf)


if __name__ == '__main__':

    main()
