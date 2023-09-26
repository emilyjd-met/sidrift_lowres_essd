from subprocess import check_call, CalledProcessError
import os
import sys
import re
from glob import glob
import argparse
from argparse import RawDescriptionHelpFormatter
sys.path.append(os.path.abspath('.'))

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')
outdir = os.path.join(here, '../../output')


def parse_args():

    default_buoydir = os.path.join(datadir, 'SID_BUOY_DATA/')
    default_outdir = os.path.join(outdir, 'individual_buoys')

    p = argparse.ArgumentParser("convert_buoy_txt2netcdf",
                                formatter_class=RawDescriptionHelpFormatter)

    p.add_argument('--buoydir', required=False, default = default_buoydir,
                   help="Top-level directory for buoy data, default {}".format(default_buoydir))
    p.add_argument('-o', '--outdir', required=False, default=default_outdir,
                   help="Output directory for the NetCDF files, default {}".format(default_outdir))
    p.add_argument('-f', '--force', action='store_true', default=False,
                   help="Force creation of new NetCDF files if they exist. Default behaviour is that they will not be created if they exist and have a timestamp newer than the text file")
    args = p.parse_args()

    return args


def main(buoydir, outdir, force):

    namedict = {1: ('aari', os.path.join(buoydir, 'aari', '*.drift.1h.txt'), 1),
                2: ('aari', os.path.join(buoydir, 'aari', '*.vector.txt'), 1),
                3: ('antsid', os.path.join(buoydir, 'AtlasofAntarcticSeaIceDrift-2004', '*.txt'), 1),
                4: ('awi', os.path.join(buoydir, 'pangaea_ipab', 'ANT-*_*_*_buoy_data.tab'), 3),
                5: ('bbb', os.path.join(buoydir, 'bbBuoy', 'bbb_baffin_bay_*.txt'), 1), 
                6: ('argos', os.path.join(buoydir, 'buoyarray', 'argos_damocles*.dat'), 1), 
                7: ('crrel', os.path.join(buoydir, 'crrel', '*_cleanPos.csv'), 1),
####                8: ('iabp', os.path.join(buoydir, 'IABP_archive', '*', '3HOURLY', '*.ll.pos'), 2),
                8: ('iabp', os.path.join(buoydir, 'IABP_3HOURLY_DATA/iabp.apl.uw.edu/Data_Products/BUOY_DATA/FULL_RESOLUTION_DATA', '*', '*', '*.dat'), 3),
                9: ('iabp', os.path.join(buoydir, 'IABP_3HOURLY_DATA/iabp.apl.uw.edu/Data_Products/BUOY_DATA/3HOURLY_DATA', '*', '*.csv'), 2),

                10: ('itp', os.path.join(buoydir, 'itp', '*rawlocs.dat'), 1),
                11: ('sams', os.path.join(buoydir, 'sams/ACCESS', '*_qcgps.txt'), 1),
                12: ('sams', os.path.join(buoydir, 'sams/ICEBELL', '*_qcgps.txt'), 1),
                13: ('sams', os.path.join(buoydir, 'sams/KOPRI', '*_qcgps.txt'), 1),
                14: ('sams', os.path.join(buoydir, 'sams/NAACOS', '*_qcgps.txt'), 1),
                15: ('sams', os.path.join(buoydir, 'sams/SIP', '*_qcgps.txt'), 1),
                16: ('sedna', os.path.join(buoydir, 'sedna/Strain', 'SEDNA_Buoy_*_GPS_2Months.nc'), 1),
                17: ('sedna', os.path.join(buoydir, 'sedna/Velocity', 'SEDNA_Buoy_*_Velocity_2Months.nc'), 1),
                18: ('simba', os.path.join(buoydir, 'SIMBA', '*.Position.dat'), 1),
                19: ('sip', os.path.join(buoydir, 'www.seaiceportal.de', '*_*_proc.csv'), 1),
                20: ('sip', os.path.join(buoydir, 'www.seaiceportal.de', '*_*_gps.csv'), 1),
            }

    failures = []
    emptydata = []
    for key in namedict:
        # Find all the files in the directory and loop over them
        flist = glob(namedict[key][1])
        for fl in flist:
            # Create a regex search string and extract the buoy name
            ss = namedict[key][1].replace('*', '(.+)')
            buoyname = re.search(ss, fl).group(namedict[key][2])

            # Call the code which writes the nc file
            cmd = "python {} {} {} --outdir {}".format(
                os.path.join(here, 'convert_buoy_txt2netcdf.py'),
                buoyname, namedict[key][0], outdir)
            if force:
                cmd = '{} -f'.format(cmd)
#            print("Running command: {}".format(cmd))
            try:
                ret = check_call(cmd, shell=True)
            except CalledProcessError as e:
                if e.returncode == 4:
                    pass
                else:
                    print("*****************************************************")
                    print("*****************************************************")
                    print("******************** WARNING ************************")
                    if e.returncode == 3:
                        print("** No good data found in file {}".format(fl))
                        emptydata.append(fl)
                    else:
                        print("** File {} failed".format(fl))
                        failures.append(fl)
                    print("*****************************************************")
                    print("*****************************************************")

    print("=============================================")
    print("Data files containing no good data:")
    for fe in emptydata:
        print(fe)
    print("List of failures:")
    for fl in failures:
        print(fl)
    print("=============================================")

    # Tara is fixed for now (buoyname in file)
    cmd = "python {} 53975 tara --outdir {}".format(
        os.path.join(here, 'convert_buoy_txt2netcdf.py'), outdir)
    try:
        check_call(cmd, shell=True)
    except:
        print("Tara buoy file failed")

    # Hudson is also fixed for now (buoyname in file)
    buoylist = ['X19909', 'X10561', 'X12273', 'X12702', 'X17042', 'X19351',
                'X12753', 'X19135', 'X19250', 'X33246', 'X33267', 'X10354',
                'X11456', 'X12502', 'X17288', 'X33519', 'X33445', 'X33556',
                'X12777', 'X33149', 'X33591']
    for buoy in buoylist:
        cmd = "python {} {} hudson --outdir {}".format(
            os.path.join(here, 'convert_buoy_txt2netcdf.py'), buoy, outdir)
        try:
            check_call(cmd, shell=True)
        except:
            print("Hudson Bay buoy {} failed".format(buoy))

    # PIPERS is also fixed for now (buoyname in file)
    buoylist = ['10421', '10422', '10423', '10424', '10425', '10426',
                '10427', '10428', '10429', '10430', '10431', '10432',
                '10433', '10434']
    for buoy in buoylist:
        cmd = "python {} {} pipers --outdir {}".format(
            os.path.join(here, 'convert_buoy_txt2netcdf.py'), buoy, outdir)
        try:
            check_call(cmd, shell=True)
        except:
            print("PIPERS buoy {} failed".format(buoy))



if __name__ == '__main__':

    args = parse_args()
    buoydir = args.buoydir
    outdir = args.outdir
    force = args.force

    main(buoydir, outdir, force)
