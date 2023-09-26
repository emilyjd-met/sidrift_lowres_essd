import sys
import os
from subprocess import call
import argparse
from datetime import datetime, timedelta


def read_args():
    '''Read and pre-process user input'''

    parser = argparse.ArgumentParser(description='Read in user-defined vars')
    parser.add_argument('-s', '--startdate', required=True,
                        help='Date as YYYYmmdd, i.e. 20141202')
    parser.add_argument('-e', '--enddate', required=True,
                        help='Date as YYYYmmdd, i.e. 20141203')
    parser.add_argument('-w', '--winddir', required=True,
                        help='Input directory for wind files')
    parser.add_argument('-o', '--outdir', required=True,
                        help='Output directory')
    parser.add_argument('-t', '--tmpdir', required=False,
                        help='Temporary working directory')
    args = parser.parse_args()

    return args


def preproc_era5(sdate, edate, winddir, outdir, tmpdir):

    windvars = ['u10', 'v10']

    d = sdate
    while d <= edate:
        # Find input and output paths
        outdiry = os.path.join(outdir, str(d.year))
        if not os.path.exists(outdiry):
            os.makedirs(outdiry)
        inpbase = 'era5_{}.nc'.format(d.strftime('%Y%m%d'))
        inpfile = os.path.join(winddir, str(d.year), inpbase)
        tmpbase = 'era5_{}_tmp.nc'.format(d.strftime('%Y%m%d'))
        tmpfile = os.path.join(tmpdir, tmpbase)
        outfile = os.path.join(outdiry, inpbase)

        # Copies the wind files across and removes unnecessary variables
        cmd = "ncks -O -v {},{} {} {}".format(windvars[0], windvars[1],
                                              inpfile, tmpfile)
        call(cmd.split())

        # Unpacks the wind variables
        cmd = "ncpdq -U {} {}".format(tmpfile, outfile)
        call(cmd.split())

        # Remove the temporary file
        os.remove(tmpfile)

        d += timedelta(days=1)


def main():

    # Read in and check user input
    args = read_args()
    startdate = args.startdate
    enddate = args.enddate
    winddir = args.winddir
    outdir = args.outdir
    if args.tmpdir is not None:
        tmpdir = args.tmpdir
    else:
        tmpdir = outdir

    sdate = datetime.strptime(startdate, '%Y%m%d')
    edate = datetime.strptime(enddate, '%Y%m%d')

    preproc_era5(sdate, edate, winddir, outdir, tmpdir)


if __name__ == '__main__':

    main()
