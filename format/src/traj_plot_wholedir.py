'''Code snippet to plot trajectory maps of individual buoy NetCDF files'''

import os
import glob
from subprocess import check_call, CalledProcessError
import argparse
from argparse import RawDescriptionHelpFormatter
from netCDF4 import Dataset

here = os.path.dirname(os.path.abspath(__file__))


def parse_args():

    p = argparse.ArgumentParser("traj_plot_wholedir",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--indir', required=True,
                   help="Filepath to directory of netcdf trajectory files")
    p.add_argument('-p', '--plotdir', required=True,
                   help="Filepath to directory to store maps")
    p.add_argument('-w', '--Wide', action='store_true', default=False,
                   help="Plot the wide angle version.")
    args = p.parse_args()

    return args


def traj_plot_wholedir(indir, plotdir, wide=False):

    flist = glob.glob(os.path.join(indir, '*.nc'))
    for fl in flist:

        # Find the region
        with Dataset(fl, 'r') as dataset:
            hemi = dataset.area
        if hemi == 'Northern Hemisphere':
            if wide:
                reg = 'nh-ease-wide'
            else:
                reg = 'nh-ease'
        elif hemi == 'Southern Hemisphere':
            if wide:
                reg = 'sh-ease-wide'
            else:
                reg = 'sh-ease'
        else:
            raise ValueError('No region found from file {}'.format(fl))

        # Create the output name
        basen = os.path.basename(fl)
        plotfile = os.path.join(plotdir, '{}_{}.png'.format(basen[:-3], reg))

        # Run the plotting
        cmd = "python {} {} -o {} -r {} -c".format(
            os.path.join(here, 'traj_plot.py'), fl, plotfile, reg)
        print(cmd)
        try:
            ret = check_call(cmd, shell=True)
        except CalledProcessError as e:
            print('Error: {}'.format(e))


if __name__ == '__main__':

    args = parse_args()
    indir = args.indir
    plotdir = args.plotdir
    wide = args.wide

    traj_plot_wholedir(indir, plotdir, wide=wide)
