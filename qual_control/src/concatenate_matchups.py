'''Small script to concatenate daily matchup files to a larger file'''

import os.path
import shutil
import argparse
from argparse import RawDescriptionHelpFormatter
from glob import glob
import subprocess
from datetime import datetime
import xarray as xr

namedict = {'multi-oi': 'multi-oi_simplex_lev3',
            'amsr-aq': 'amsr-aq_tb37v-Lap+tb37h-Lap_simplex_lev2',
            'amsr2-gw1': 'amsr2-gw1_tb37v-Lap+tb37h-Lap_simplex_lev2',
            'ssmi-f10': 'ssmi-f10_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmi-f11': 'ssmi-f11_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmi-f13': 'ssmi-f13_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmi-f14': 'ssmi-f14_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmi-f15': 'ssmi-f15_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmis-f16': 'ssmis-f16_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmis-f17': 'ssmis-f17_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'ssmis-f18': 'ssmis-f18_tb90v-Lap+tb90h-Lap_simplex_lev2',
            'wind': 'ecmwf-era5_none_wind_lev2'}



def parse_args():

    valid_sensors = ['multi-oi', 'amsr-aq', 'amsr2-gw1', 'ssmi-f10',
                     'ssmi-f11', 'ssmi-f13', 'ssmi-f14', 'ssmi-f15',
                     'ssmis-f16', 'ssmis-f17', 'ssmis-f18',
                     'wind']

    p = argparse.ArgumentParser("concatenate_matchups",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-y', '--year', required=True,
                   help="Year to concatenate")
    p.add_argument('-m', '--month', required=False, default=None,
                   help="Month to concatenate, if not set then a whole "
                   "year is concatenated")
    p.add_argument('-s', '--sensor', required=False, default='multi-oi',
                   choices=valid_sensors,
                   help="Sensor whose matchup files to concatenate")
    p.add_argument('-i', '--indir', required=True,
                   help="Directory to concatenate the files from")
    p.add_argument('-o', '--outdir', required=True,
                   help="Directory to write the output concatenated file to")

    args = p.parse_args()

    return args


def concatenate_matchups(year, month, sensor, indir, outdir):

    # Loop over the hemispheres
    for hemi in ['nh', 'sh']:

        # Determine pattern for yearly files, if month is not specified
        if month is None:
            pdate = datetime.strptime('{}0101'.format(year), '%Y%m%d')
            patt = "{}/{:%Y}/*/matchup_{}_{}-ease2-750_*12w_*12w_NN3Dcol.nc".format(indir, pdate, namedict[sensor], hemi)
        # Otherwise determine pattern for monthly files
        else:
            pdate = datetime.strptime('{}{}01'.format(year, month.zfill(2)),
                                      '%Y%m%d')
            patt = "{}/{:%Y}/{:%m}/matchup_{}_{}-ease2-750_*12w_*12w_NN3Dcol.nc".format(indir, pdate, pdate, namedict[sensor], hemi)

        # Glob the files, check uniqueness and sort
        flist = glob(patt)
        uniqueflist = list(set(flist))
        flist = sorted(uniqueflist)

        # Use xarray to concatenate these files
        dataset = xr.open_mfdataset(flist, combine='nested',
                                    concat_dim='matchup')

        ## Take the first file and add an unlimited fimension
        #firstfile = flist.pop(0)
        #unlimfile = "{}/temp/first_matchup_file_{}_{}_{:%Y%m}.nc".format(outdir, sensor, hemi, pdate)
        #ncks_cmd = "ncks --no_abc -O --mk_rec_dmn matchup {} {}".format(firstfile, unlimfile)
        #subprocess.check_call(ncks_cmd.split())

        # Determining filename
        if month is None:
            outfile = "{}/matchup_{}_{}-ease2-750_{:%Y}_NN3Dcol.nc".format(outdir, sensor, hemi, pdate)
        else:
            outfile = "{}/matchup_{}_{}-ease2-750_{:%Y%m}_NN3Dcol.nc".format(outdir, sensor, hemi, pdate)            

        ## Append all the files to this first one
        #strlist = ' '.join(flist)
        #ncrcat_cmd = "ncrcat {} {} {}".format(unlimfile, strlist, outfile)
        #subprocess.check_call(ncrcat_cmd.split())

        # Write out dataset
        dataset.to_netcdf(path=outfile, mode='w')

        print("Written {}".format(outfile))

        # Close dataset
        dataset.close()


def main():

    args = parse_args()
    year = args.year
    month = args.month
    sensor = args.sensor
    indir = args.indir
    outdir = args.outdir

    concatenate_matchups(year, month, sensor, indir, outdir)


if __name__ == '__main__':

    main()
