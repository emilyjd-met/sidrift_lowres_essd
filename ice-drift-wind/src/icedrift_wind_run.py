import os
import sys
import shutil
import argparse
from glob import glob
from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'model'))
from time_periods import param_find
from remap_nwp import get_aggr_fields
from icedrift_from_winds import compute_drift

def read_args():
    '''Read and pre-process user input'''

    valid_area = ['osi405_nh', 'osi405_sh', 'osi455_nh', 'osi455_sh']
    valid_source = ['ecmwf', 'era5']
    valid_icefmt = ['proc', 'final']
    valid_icemasktype = ['winter', 'summer']

    parser = argparse.ArgumentParser(description='Ice drift using winds')

    parser.add_argument('-e', '--enddate', required=False,
                        help='End date of drift as YYYYmmdd (default is today)')

    parser.add_argument('-d', '--duration', required=False, default=2,
                        help='Time in days over which to measure the drift')

    parser.add_argument('-m', '--mode', required=False, default='oper',
                        help='Mode of processing, "oper" or "reproc"')

    parser.add_argument('-s', '--source', required=False, default='emcwf',
                        choices=valid_source,
                        help='Mode of processing, valid options: {}'.format(valid_source))

    parser.add_argument('-a', '--area', required=True, choices=valid_area,
                        help='Area grid, valid options: {}'.format(valid_area))

    parser.add_argument('-f', '--fimexdir', required=True,
                        help='Parameter directory for fimex regridding')

    parser.add_argument('-p', '--pardir', required=True,
                        help='Parameter directory for files to retrieve drift '
                        'from wind')

    parser.add_argument('--pstart', required=False, default='2013',
                        help='Start year for parameter file')

    parser.add_argument('--pend', required=False, default='2020',
                        help='End year for parameter file')

    parser.add_argument('-i', '--inwinddir', required=True,
                        help='Input directory for input wind files')

    parser.add_argument('-z', '--tmpdir', required=True,
                        help='Temporary working directory')

    parser.add_argument('-w', '--aggwinddir', required=True,
                        help='Directory to store aggregated wind files')

    parser.add_argument('-l', '--icemaskloc', required=True,
                        help='Directories where single sensor files containing ice masks are located, in a comma-separated list')

    parser.add_argument('-t', '--icemasktype', required=False,
                        default='winter', choices=valid_icemasktype,
                        help="Type of masking to use, the winter masking masks more than the summer masking, default is 'winter', valid choices {}".format(valid_icemasktype))

    parser.add_argument('-o', '--outdir', required=True,
                        help='Output directory')

    parser.add_argument('-fmt', '--icefmt', required=False, default='proc',
                   choices=valid_icefmt,
                   help='Format of ice files, default is proc, valid choices are {}'.format(valid_icefmt))

    parser.add_argument('-v', '--verbose', help='Increase output verbosity',
                        action='store_true')

    args = parser.parse_args()

    # Determining which date to run for
    if args.enddate is None:
        # If no date input, use "today" as default
        args.enddate = datetime.utcnow
    else:
        try:
            # Make sure input date is valid
            args.enddate = datetime.strptime(args.enddate, '%Y%m%d')
        except ValueError:
            print("\nERROR - Invalid date string, need YYYYmmdd.\n")
            sys.exit(1)
    print("\n - Run {} ({:%Y%m%d})\n".format(os.path.basename(__file__),
                                             args.enddate))

    return args


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def paramfilefind(area, pstart, pend, mon, duration, pardir):

    # If the parameter file is found with the passed years, return directly
    paramfile = os.path.join(pardir,
                             'inv_params_{}_{}-{}_{}_{}day.nc'
                             ''.format(area, pstart, pend, mon, duration))
    if os.path.isfile(paramfile):
        return paramfile

    # Otherwise try years shifted by +1
    pstartp = str(int(pstart) + 1)
    pendp = str(int(pend) + 1)
    paramfile = os.path.join(pardir,
                             'inv_params_{}_{}-{}_{}_{}day.nc'
                             ''.format(area, pstartp, pendp, mon,
                                       duration))
    if os.path.isfile(paramfile):
        return paramfile

    # Then years shifted by -1
    pstartm = str(int(pstart) - 1)
    pendm = str(int(pend) - 1)
    paramfile = os.path.join(pardir,
                             'inv_params_{}_{}-{}_{}_{}day.nc'
                             ''.format(area, pstartm, pendm, mon,
                                       duration))
    if os.path.isfile(paramfile):
        return paramfile

    # If none of these parameter files are found, return none
    return None


def ice_drift_wind(enddate, duration, mode, source, area, fimexdir, pardir,
                   pstart, pend, inwinddir, tmpdir, aggwinddir, icemaskloc,
                   outdir, verbose, icemasktype='winter', icefmt='proc'):
    '''Run the ice drift from winds.'''

    # Setting paths
    try:
        os.environ['OSI_ROOT']
    except KeyError:
        dirname = os.path.abspath(os.path.dirname(__file__))
        os.environ['OSI_ROOT'] = os.path.abspath(os.path.join(dirname,
                                                              '../../../..'))
        print("OSI_ROOT is set locally as {}".format(os.environ['OSI_ROOT']))

    # Check if the wind file exists
    startdate = enddate - timedelta(days=duration)
    windf = 'NWP_{}_aggr_{:%Y%m%d}12-{:%Y%m%d}12.nc'.format(area, startdate,
                                                            enddate)
    windf = os.path.join(aggwinddir,
                         '{}'.format(datetime.strftime(enddate, '%Y')),
                         '{}'.format(datetime.strftime(enddate, '%m')), windf)
    print("WINDF = ", windf)
    if not os.path.isfile(windf):
        try:
            aggr_file = get_aggr_fields(area, enddate-timedelta(days=duration),
                                        enddate, source=source, mode=mode,
                                        fimexdir=fimexdir, winddir=inwinddir,
                                        tmpdir=tmpdir, outdir=aggwinddir,
                                        verbose=verbose)
            print("Aggregated wind file {} created.".format(aggr_file))
        except:
            print("WARNING: Aggregated wind file {} does not exist and could not be created.".format(windf))
    else:
        print("Aggregated wind file {} already exists: not creating.".format(windf))

    # Selection of parameter files (summer/winter)
#    summer_st = datetime.strptime(pdate.strftime('%Y') + '0601', '%Y%m%d')
#    summer_end = datetime.strptime(pdate.strftime('%Y') + '0930', '%Y%m%d')
#    if (pdate >= summer_st) and (pdate <= summer_end):
#        paramfile = os.path.join(pardir,
#                                 'inversion_parameters_2013-2020_summer.nc')
#    else:
#        paramfile = os.path.join(pardir,
#                                 'inversion_parameters_2014-2020_winter.nc')

    # Selection of parameter files for a monthly average
    mon1, weight1, mon2, weight2 = param_find(enddate)
    paramfile1 = paramfilefind(area, pstart, pend, mon1, duration, pardir)
    paramfile2 = paramfilefind(area, pstart, pend, mon2, duration, pardir)
    print("For date {}, using parameter files:".format(datetime.strftime(enddate, '%Y%m%d')))
    print("\t {} : weighting {}".format(paramfile1, weight1))
    print("\t {} : weighting {}".format(paramfile2, weight2))

    # The ice drift from winds is run for each parameter file
    dfile1 = compute_drift(enddate, duration, area, aggwinddir, paramfile1,
                           out_dir=tmpdir, suffix=mon1)
    print("dfile1 = ", dfile1)
    dfile2 = compute_drift(enddate, duration, area, aggwinddir, paramfile2,
                           out_dir=tmpdir, suffix=mon2)
    print("dfile2 = ", dfile2)

    # Merging of drift files
    dfile = merge_drift(dfile1, dfile2, weight1, weight2, enddate, outdir,
                        icefmt=icefmt)

    # Applying the icemask from a reference single sensor file
    dfile = apply_icemask(dfile, startdate, enddate, area, icemaskloc,
                          icemasktype=icemasktype)


def merge_drift(dfile1, dfile2, weight1, weight2, enddate, outdir,
                icefmt='proc'):

    # Create an output name
    outdir = os.path.join(outdir, datetime.strftime(enddate, '%Y'),
                         datetime.strftime(enddate, '%m'),
                         datetime.strftime(enddate, '%d'))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    dbasename = '{}.nc'.format(os.path.basename(dfile1)[:-7])
    dfile = os.path.join(outdir, dbasename)

    # Copy the first file to the output name
    shutil.copy(dfile1, dfile)

    # Read in the drifts from each file
    with Dataset(dfile1) as dataset1:
        driftx1 = dataset1['driftX'][0, :, :]
        drifty1 = dataset1['driftY'][0, :, :]
    with Dataset(dfile2) as dataset2:
        driftx2 = dataset2['driftX'][0, :, :]
        drifty2 = dataset2['driftY'][0, :, :]

    # Masking by the fill value
    fillval = 9.96920996e+36
    driftx1.mask = np.logical_or(driftx1.mask, driftx1 >= fillval)
    drifty1.mask = np.logical_or(drifty1.mask, drifty1 >= fillval)
    driftx2.mask = np.logical_or(driftx2.mask, driftx2 >= fillval)
    drifty2.mask = np.logical_or(drifty2.mask, drifty2 >= fillval)

    # Create new results arrays
    driftx = np.full_like(driftx1, 0)
    driftx = ma.array(driftx)
    drifty = np.full_like(drifty1, 0)
    driftx = ma.array(drifty)
    # Merge
    driftx = ((weight1 * driftx1) + (weight2 * driftx2)) / (weight1 + weight2)
    drifty = ((weight1 * drifty1) + (weight2 * drifty2)) / (weight1 + weight2)

    # TODO How to merge other parameters apart from the drift

    # If the file should be a proc format file, invert y
    if icefmt == 'proc':
        drifty = drifty * -1.0

    # Write the merged drifts to the output file
    with Dataset(dfile, 'r+') as dataset:
        dataset['driftX'][0, :, :] = driftx
        dataset['driftY'][0, :, :] = drifty

    # Remove temporary files
#    if os.path.isfile(dfile1):
#        os.remove(dfile1)
#    if os.path.isfile(dfile2):
#        os.remove(dfile2)

    print("Output drift file written to {}".format(dfile))

    return dfile


def apply_icemask(dfile, startdate, enddate, area, icemaskloc,
                  icemasktype='winter'):
    '''Locate single sensor ice drift files, compile an ice mask and
    apply it.
    The mask returned has a landmask and a noice mask.
    The noice mask is either:
    Full -
    All except flags 0 (icedrift OK), 13 (corrected by neighbours) and
    17 (smaller pattern) are flagged.
    Summer -
    Only the noice (4) and close_to_coast_or_edge (5) are flagged.
    '''

    # Setting the data which is valid
    if icemasktype == 'winter':
        valid_data = [0, 13, 17]
    elif icemasktype == 'summer':
        valid_data = [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

    # Cycle over the filepaths and find the single sensor satellite files
    icemaskfiles = []
    locs = icemaskloc.split(',')
    for loc in locs:
        icemaskfiles.extend(find_icemask_files(startdate, enddate, area, loc))
    # If no icemask is found, find the one closest in time
    if not icemaskfiles:
        sd = startdate
        ed = enddate
        tdelt = 1
        while not icemaskfiles:
            sd = startdate - timedelta(days=tdelt)
            ed = enddate - timedelta(days=tdelt)
            for loc in locs:
                icemaskfiles.extend(find_icemask_files(sd, ed, area, loc))
            if icemaskfiles:
                break
            else:
                sd = startdate + timedelta(days=tdelt)
                ed = enddate + timedelta(days=tdelt)
                for loc in locs:
                    icemaskfiles.extend(find_icemask_files(sd, ed, area, loc))
                if icemaskfiles:
                    break
                else:
                    tdelt += 1

    print("Masking non-ice and land with a mask from single sensor references:")
    for imf in icemaskfiles:
          print("\t{}".format(imf))

    # Initialise a datacube
    with Dataset(icemaskfiles[0], 'r') as refdataset:
        flagf = refdataset['flag'][0, :, :]
        refshape = flagf.shape
        refflags = np.zeros((len(icemaskfiles), refshape[0],
                             refshape[1]), dtype=np.int8)

    # Read the flags from all the reference single sensor icedrift files
    # into the datacube
    for i, imf in enumerate(icemaskfiles):
        with Dataset(imf, 'r') as refdataset:
            refflags[i, :, :] = refdataset['flag'][0, :, :]

    # The land should be the same in every image
    land = np.zeros((refshape[0], refshape[1]), dtype=bool)
    land[refflags[0, :, :] == 3] = 1

    # For noice, we want the area where none of the single sensor files
    # have valid data
    icecube = np.zeros((len(icemaskfiles), refshape[0], refshape[1]),
                         dtype=np.int8)
    for val in valid_data:
        icecube[:, :, :][refflags == val] = 1
    ice = np.any(icecube, axis=0)

    # Now applying the flagging in the wind drift file
    with Dataset(dfile, 'r+') as dataset:

        # Reading the wind ice drift flags
        flag = dataset['flag'][0, :, :]

        # Setting the land (3) and noice (4)
        flag[ice == 0] = 4
        flag[land] = 3
        refmask = np.logical_or(flag == 3, flag == 4)

        # And writing the flags back to the icedrift file
        dataset['flag'][0, :, :] = flag

        # Applying the masking to the other variables
        varlist = ['driftX', 'driftY', 't0', 't1', 'sigdX', 'sigdY',
                   'length', 'dir', 'latStart', 'lonStart', 'latEnd',
                   'lonEnd']
        for var in varlist:
            fv = dataset[var]._FillValue
            datalayer = dataset[var][0, :, :]
            datalayer[refmask] = fv
            dataset[var][0, :, :] = datalayer



def find_icemask_files(startdate, enddate, area, icemaskloc):
    '''Search a directory structure for a single sensor ice drift files
    for the matching dates and return the first one found'''

    icemaskfiles = []

    # Making a list of possible filenames
    griddict = {'osi405_nh': 'nh-polstere-625',
                'osi405_sh': 'sh-polstere-625',
                'osi455_nh': 'nh-ease2-750',
                'osi455_sh': 'sh-ease2-750'}
    grid = griddict[area]
    inst_wb = [('amsr2-gw1', '37'),
               ('amsr-aq', '37'),
               ('ssmis-f18', '90'),
               ('ssmis-f17', '90'),
               ('ssmis-f16', '90'),
               ('ssmi-f15', '90'),
               ('ssmi-f14', '90'),
               ('ssmi-f13', '90'),
               ('ssmi-f11', '90'),
               ('ssmi-f10', '90')]
    patt = 'icedrift_{}_tb{}v-Lap+tb{}h-Lap_simplex_lev2_{}_{}12w_{}12w.nc'
    fnamelist = []
    for iwb in inst_wb:
        fnamelist.append(patt.format(iwb[0], iwb[1], iwb[1], grid,
                                     datetime.strftime(startdate, '%Y%m%d'),
                                     datetime.strftime(enddate, '%Y%m%d')))

    # Use a path including the date
    globlist = glob(os.path.join(icemaskloc, datetime.strftime(enddate, '%Y'),
                                 datetime.strftime(enddate, '%m'),
                                 datetime.strftime(enddate, '%d'), '*.nc'))
    for gl in globlist:
        if os.path.basename(gl) in fnamelist:
            icemaskfiles.append(gl)

    # We need a check on the files - if they have no good data in the
    # flag field, they should be thrown out
    badlist = []
    for imf in icemaskfiles:
        with Dataset(imf, 'r') as refdataset:
            flags = refdataset['flag'][0, :, :]
            if np.max(flags) < 4 or np.min(flags) > 0:
                badlist.append(imf)
    for badfile in badlist:
        icemaskfiles.remove(badfile)

    return icemaskfiles


def main():

    # Read in and check user input
    args = read_args()
    enddate = args.enddate
    duration = int(args.duration)
    mode = args.mode
    source = args.source
    area = args.area
    fimexdir = args.fimexdir
    pardir = args.pardir
    pstart = args.pstart
    pend = args.pend
    inwinddir = args.inwinddir
    tmpdir = args.tmpdir
    aggwinddir = args.aggwinddir
    icemaskloc = args.icemaskloc
    icemasktype = args.icemasktype
    outdir = args.outdir
    icefmt = args.icefmt
    verbose = args.verbose

    ice_drift_wind(enddate, duration, mode, source, area, fimexdir, pardir,
                   pstart, pend, inwinddir, tmpdir, aggwinddir, icemaskloc,
                   outdir, verbose, icemasktype=icemasktype, icefmt=icefmt)


if __name__ == '__main__':

    main()
