import sys
import os
from subprocess import call
import argparse
from datetime import datetime, timedelta
import glob
from netCDF4 import Dataset
from fix_line import fix_line

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')

nwproot = os.path.join(datadir, 'era5/raw')
ecmwf_raw = os.path.join(datadir, 'ecmwf/grib/')
ecmwf_reproc = os.path.join(datadir, 'ECMWF_oper_archive')
fimexdir = os.path.join(here, '../par')
tmpdir = os.path.join(datadir, 'temp')
area = 'osi405_nh'
time_step = 1.0
tsFormat = "%Y%m%d"


def read_args():
    '''Read and pre-process user input'''

    valid_area = ['osi405_nh', 'osi405_sh', 'osi455_nh', 'osi455_sh', 'polstere-125_nh', 'polstere-125_sh']
    valid_source = ['ecmwf', 'era5']
    valid_mode = ['oper', 'reproc']

    parser = argparse.ArgumentParser(description='Read in user-defined vars')
    parser.add_argument('-e', '--enddate', required=True,
                        help='Date as YYYYmmdd or YYYYmmddHH. If no hours'
                        'are specified, this defaults to 12h.')

    parser.add_argument('-d', '--duration', required=False, default=2,
                        help='Duration over which to create the remapped '
                        'wind files (default is 2)')

    parser.add_argument('-a', '--area', required=True, choices=valid_area,
                        help='Area grid, valid options: {}'.format(valid_area))

    parser.add_argument('-s', '--source', required=False,
                        choices=valid_source,
                        help='Wind files to process, valid options: {}'
                        ''.format(valid_source))

    parser.add_argument('-m', '--mode', required=False, choices=valid_mode,
                        help='Mode of processing, valid options: {}'
                        ''.format(valid_mode))

    parser.add_argument('-x', '--fimexdir', required=False,
                        help='Parameter directory for fimex information')

    parser.add_argument('-w', '--winddir', required=False,
                        help='Input directory for wind files')

    parser.add_argument('-t', '--tmpdir', required=False, default='./tmp',
                        help='Temporary working directory')

    parser.add_argument('-o', '--outdir', required=True,
                        help='Output directory')

    parser.add_argument('-r', '--reproc', action='store_true', default=False,
                        help='Put the output wind files into YYYY/mm folders')

    parser.add_argument('-f', '--force', action='store_true', default=False,
                        help='Force creation if file already exists')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase output verbosity')

    args = parser.parse_args()


    # Make sure input date is valid
    try:
        edate = datetime.strptime(args.enddate, '%Y%m%d')
    except ValueError:
        try:
            edate = datetime.strptime(args.enddate, '%Y%m%d%H')
        except ValueError:
            print("ERROR: Invalid date string, need YYYYmmdd or YYYYmmddHH.")
            sys.exit(1)

    print("\n - Run {} ({:%Y%m%d})\n".format(os.path.basename(__file__),
                                             edate))

    return args


def find_ecmwf_files(wind_date, ecmwfdir, mode='oper'):

    """
    Copied and modified from nwp.py by Atle Sørenson

    wind_date:      Date of required ECMWF files

    We pick the 'freshest' files valid for the datetime, ranking by
    analysis time.
    """

    ana_dates = []

    # ECMWF filenames are of the form
    # N1S{xxxxxxxx}{yyyyyyyy}1
    # xxxxxxxx is the analysis time in mmddHHMM format
    # yyyyyyyy is the file time in mmddHHMM format
    if mode == 'oper':
        patt = '{}/N1S????????{}????1'.format(ecmwf_raw,
                                              wind_date.strftime('%m%d'))
    elif mode == 'reproc':
        patt = '{}/{}/osisaf_seaice_oper_{}????'.format(ecmwf_reproc,
                wind_date.strftime('%Y'), wind_date.strftime('%Y%m%d'))
    else:
        raise ValueError("Unrecognised mode: {}".format(mode))

    date_files = sorted(glob.glob(patt))
    for fname in date_files:

        # Extract the analysis time associated with these files
        yearnow = datetime.utcnow().year
        if mode == 'oper':
            ana_date = datetime.strptime('{}{}'.format(yearnow,
                                         os.path.basename(fname)[3:9]),
                                         '%Y%m%d%H')
        elif mode == 'reproc':
            ana_date = datetime.strptime(os.path.basename(fname)[19:29],
                                         '%Y%m%d%H')
        # Handle the year change
        ana_dates.append(ana_date)

    # Find the latest analysis time to pick the 'freshest' files. Return
    # a warning if no files are found
    if not ana_dates:
        print('No ECMWF files covering time {} found.'.format(wind_date))
        sys.exit(1)
    else:
        latest_ana_date = sorted(ana_dates)[-1]

        # Find the set of files from that analysis time
        if mode == 'oper':
            patt = '{}/N1S{}{}????1'.format(ecmwf_raw,
                    latest_ana_date.strftime('%m%d%H%M'),
                    wind_date.strftime('%m%d'))
        elif mode == 'reproc':
            patt = '{}/{}/osisaf_seaice_oper_{}??'.format(ecmwf_reproc,
                    latest_ana_date.strftime('%Y'),
                    latest_ana_date.strftime('%Y%m%d%H'))
        date_files = sorted(glob.glob(patt))

        # Check that some files were found
        if not date_files:
            print('No ECMWF files covering time {} found.'.format(wind_date))
            sys.exit(1)

    return date_files, patt


def get_aggr_fields(area, date0, date1, source='ecmwf', mode='oper',
                    fimexdir=fimexdir, winddir=ecmwf_raw, tmpdir=tmpdir,
                    outdir=tmpdir, reproc=False, force=False, verbose=False):

    print("winddir = ", winddir)
    print("mode = ", mode)
    print("source = ", source)
    print("outdir = ", outdir)

    # TODO - This code assumes that the wind variables are unpacked in the
    # input NetCDF, i.e. no scale_factor and no add_offset

    d = date0
    tmpfiles = []
    tmpfiles2 = []

    # Setting up the output directory
    if reproc:
        outdir = os.path.join(outdir, date1.strftime("%Y"),
                              date1.strftime("%m"))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Defining the wind variables according to how the file is created
    if source == 'era5':
        wind_vars = ('u10', 'v10')
    elif source == 'ecmwf':
        wind_vars = ('x_wind_10m', 'y_wind_10m')

    while d <= date1:
        ts = d.strftime(tsFormat)
        outbase = 'NWP_{}_{}.nc'.format(area, ts)
        outfile = os.path.join(tmpdir, outbase)
        outfile2 = outfile[:-3] + '_tmp.nc'
        ts1 = date0.strftime("%Y%m%d%H")
        ts2 = date1.strftime("%Y%m%d%H")
        outbase3 = 'NWP_{}_aggr_{}-{}.nc'.format(area,ts1,ts2)
        outfile3 = os.path.join(outdir, outbase3)

        if os.path.exists(outfile3):
            if not force:
                print('File {} exists and force flag is not set, exiting.'.format(outfile3))
                sys.exit()
            else:
                os.remove(outfile3)
                if os.path.exists(outfile):
                    os.remove(outfile)

        cmd = "fimex"
        cmd += " -c {}/fimex_remap_{}.cfg".format(fimexdir, area)

        if source == 'era5':
            inpbase = 'era5_{}.nc'.format(ts)
            inpfile = os.path.join(winddir, str(d.year), inpbase)
            print("Reading from {}".format(inpfile))
            cmd += " --input.config {}/fimex_input_nc.ncml".format(fimexdir)
            cmd += " --input.file {}".format(inpfile)
        elif source == 'ecmwf':
            grib_files, patt = find_ecmwf_files(d, winddir, mode)
            cmd += " --input.type grib"
            cmd += " --input.config {}/fimex_input_grib.ncml".format(fimexdir)
            cmd += " --input.file glob:{}".format(patt)
        else:
            raise ValueError("Unrecognised source {}\n".format(source))

        cmd += " --process.rotateVector.direction grid"
        cmd += " --process.rotateVector.x {}".format(wind_vars[0])
        cmd += " --process.rotateVector.y {}".format(wind_vars[1])

        cmd += " --output.file {}".format(outfile)
        cmd += " --extract.reduceTime.start {}".format(date0.strftime("%Y-%m-%dT%H:%M:%S"))
        cmd += " --extract.reduceTime.end {}".format(date1.strftime("%Y-%m-%dT%H:%M:%S"))
        if verbose:
            print(cmd)
        call(cmd.split())

        # There is a weird artifact caused by a fimex bug for the grib files.
        # I am currently using a quick script to fix this
        if source == 'ecmwf':
            fix_line(outfile)

        # The time dimension is not currently the record dimension, so this
        # needs to be fixed with ncks
        cmd = 'ncks -O --mk_rec_dmn time'
        cmd += ' {}'.format(outfile)
        cmd += ' -o {}'.format(outfile2)
        if verbose:
            print(cmd)
        call(cmd.split())

        # Then need to unpack the variables as the scale factors and offsets
        # are not the same between the files. Overwrite the first file
        cmd = 'ncpdq -U -O'
        cmd += ' {}'.format(outfile2)
        cmd += ' -o {}'.format(outfile)
        if verbose:
            print(cmd)
        call(cmd.split())

        tmpfiles.append(outfile)
        tmpfiles2.append(outfile2)
        d += timedelta(hours=24)

    # Aggregate the individual files into 1 file
    cmd = "ncrcat -v {},{} {} {}".format(wind_vars[0], wind_vars[1],
                                         ' '.join(tmpfiles), outfile3)
    if verbose:
        print(cmd)
    call(cmd.split())

    grp = Dataset(outfile3,"r+")

    # compute time-averaged u_wind and v_wind
    for v in wind_vars:
        nc_wind = grp.variables[v]

        # mean
        nc_wind_av_var = grp.createVariable('{}_avg'.format(v), 'f4',
                                        nc_wind.dimensions[1:])
        nc_wind_av_var[:, :] = nc_wind[:, :, :].mean(axis=0)
        nc_wind_av_var.long_name = "wind component {} averaged over 48 hours".format(v)
        nc_wind_av_var.grid_mapping = nc_wind.grid_mapping
        nc_wind_av_var.coordinates  = nc_wind.coordinates
        nc_wind_av_var.units = "m/s"

        # stddev
        nc_wind_std_var = grp.createVariable('{}_std'.format(v), 'f4',
                                             nc_wind.dimensions[1:])
        nc_wind_std_var[:, :] = nc_wind[:, :, :].std(axis=0)
        nc_wind_std_var.long_name = "wind component {} standard-deviation over 48 hours".format(v)
        nc_wind_std_var.grid_mapping = nc_wind.grid_mapping
        nc_wind_std_var.coordinates  = nc_wind.coordinates
        nc_wind_std_var.units = "m/s"

    grp.close()

    # remove some variables (mainly the time-dependent fields: we only
    # keep the time-averaged ones)
    cmd = "ncks -O -x -v {},{} {} {}".format(wind_vars[0], wind_vars[1],
                                             outfile3, outfile3)
    if verbose:
        print(cmd)
    call(cmd.split())

    # Rename some variables (especially that fimex for some reason
    # names xc and yc as latitude and longitude)
    cmd = 'ncrename -O'
    cmd += ' -d latitude,yc -v latitude,yc'
    cmd += ' -d longitude,xc -v longitude,xc'
    cmd += ' -v {}_avg,uwind_avg -v {}_std,uwind_std'.format(wind_vars[0],
                                                             wind_vars[0])
    cmd += ' -v {}_avg,vwind_avg -v {}_std,vwind_std'.format(wind_vars[1],
                                                             wind_vars[1])
    cmd += ' {}'.format(outfile3)
    if verbose:
        print(cmd)
    call(cmd.split())

    # Remove temporary files
    for f in tmpfiles:
        os.remove(f)
    for f in tmpfiles2:
        os.remove(f)

    return outfile3


def main():

    # Read in and check user input
    args = read_args()
    enddate = args.enddate
    duration = int(args.duration)
    area = args.area
    source = args.source
    mode = args.mode
    fimexdir = args.fimexdir
    winddir = args.winddir
    tmpdir = args.tmpdir
    outdir = args.outdir
    reproc = args.reproc
    force = args.force
    verbose = args.verbose

    # If no hours are specified, this defaults to 12h
    if len(enddate) == 8:
        enddate = '{}12'.format(enddate)
    elif len(enddate) != 10:
        raise ValueError("enddate must be in format YYmmdd or YYmmddHH, "
                         "not {}".format(enddate))
    edate = datetime.strptime(enddate, '%Y%m%d%H')
    sdate = edate - timedelta(days=duration)

    resFile = get_aggr_fields(area, sdate, edate, source=source, mode=mode,
                              fimexdir=fimexdir, winddir=winddir,
                              tmpdir=tmpdir, outdir=outdir, reproc=reproc,
                              force=force, verbose=verbose)

    print("Aggregated NWP file ready as {}".format(resFile))


if __name__ == '__main__':

    main()
