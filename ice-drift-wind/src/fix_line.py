import shutil
import argparse
import numpy as np
import math
import datetime
from netCDF4 import Dataset


def p_parseargs():
    """Read in the user-defined arguments"""

    parser = argparse.ArgumentParser(description='Read in user-defined vars')

    parser.add_argument('-i', required=True, dest='infile', type=str,
                        help='''Filename of the file to fix''')
    parser.add_argument('-o', required=False, dest='outfile', type=str,
                        default=None,
                        help='''Output filename for the fixed file''')

    args = parser.parse_args()

    return vars(args)


def fix_line(infile, outfile=None):
    """Fix the inverted line from i=60, j>=92 in x_wind_10m and y_wind_10m"""

    print("Fixing the fimex bug at lat=-45")
    if outfile is not None:
        shutil.copy(infile, outfile)
    else:
        outfile = infile

    with Dataset(outfile, 'r+') as dataset:

        lon = dataset['lon'][:]
        # Finding where lon = -45
        x45, y45 = np.where(np.logical_and(lon > -45.001, lon < -44.999))
#        print("lat/lon: ", min(x45), max(x45), min(y45), max(y45))
        # Cycle over the wind components
        for wind in ['x_wind_10m', 'y_wind_10m']:
            print("Wind component:", wind)
            # Cycle over the timesteps
            for z in range(len(dataset[wind][:, 0, 0, 0])):
                print("Time: ", datetime.datetime.fromtimestamp(dataset['time'][z]).strftime("%Y%m%d%H"))
                # Check the differences to the neighbouring columns.
                # The column with lon=-45 is checked against the neighbouring
                # columns with + and - one pixel, while the two neighbouring
                # columns are checked against each other as a baseline "good"
                # variation
                diff_sum_sq_testp = 0
                diff_sum_sq_testm = 0
                diff_sum_sq_good = 0
                samples = 0
                for i in range(len(x45)):
                    samples += 1
                    diff_sum_sq_testp += ((dataset[wind][z,0,x45[i],y45[i]] -
                                           dataset[wind][z,0,x45[i],y45[i]+1]) *
                                          (dataset[wind][z,0,x45[i],y45[i]] -
                                           dataset[wind][z,0,x45[i],y45[i]+1]))
                    diff_sum_sq_testm += ((dataset[wind][z,0,x45[i],y45[i]] -
                                           dataset[wind][z,0,x45[i],y45[i]-1]) *
                                          (dataset[wind][z,0,x45[i],y45[i]] -
                                           dataset[wind][z,0,x45[i],y45[i]-1]))
                    diff_sum_sq_good += ((dataset[wind][z,0,x45[i],y45[i]+1] -
                                          dataset[wind][z,0,x45[i],y45[i]-1]) *
                                         (dataset[wind][z,0,x45[i],y45[i]+1] -
                                          dataset[wind][z,0,x45[i],y45[i]-1]))

                std_testp = math.sqrt(diff_sum_sq_testp / samples)
                std_testm = math.sqrt(diff_sum_sq_testm / samples)
                std_good = math.sqrt(diff_sum_sq_good / samples)
#                print("std_testp = ", std_testp)
#                print("std_testm = ", std_testm)
#                print("std_good = ", std_good)

                if ((std_testp > 2 * std_good) and (std_testm > 2 * std_good)):
                    print("Fixing...")
                    for i in range(len(x45)):
                        dataset[wind][z, 0, x45[i], y45[i]] *= -1.
#                    print("...done")

#        dataset[wind][:,0,92:,60] *= -1.



if __name__ == '__main__':

    uservars = p_parseargs()
    infile = uservars.get('infile')
    outfile = uservars.get('outfile')

    fix_line(infile, outfile)
