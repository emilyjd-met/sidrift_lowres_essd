#!/usr/bin/env python3

import os
import sys
from datetime import datetime, date, time, timedelta
from netCDF4 import Dataset
import glob
from subprocess import check_call, check_output, CalledProcessError
import itertools
import pyresample as pr
import numpy as np
from itertools import combinations
from icedrift_wrapper import icedrift_wrapper

import uuid

here = os.path.dirname(os.path.abspath(__file__))
par = os.path.join(here, '../par')
datadir = os.path.join(here, '../../data')

# we use some modules from the pmr_gridding package, which must be checked out at the same level as
# ice-tracking/ with the command
# git clone https://gitlab.met.no/thomasl/pmr_gridding.git
sid_rrdp_sftw = os.path.join(here, '/../../pmr_gridding/')
sys.path.append(sid_rrdp_sftw)
import gridding.utils as gr_utils
import gridding.io_handler as gr_io

def find_concv2_file(area,dt,ext='nc',strict_date=True):

    dates = [dt,]
    if not strict_date:
        for marg in range(1,3):
            dates.append(dt + timedelta(days=marg),)
            dates.append(dt - timedelta(days=marg),)

    for ds in dates:
        found = False
        if ds >= datetime.today() - timedelta(16):
            # dm1
            cfile = 'ice_conc_{}_ease2-250_dm1-amsr2_{:%Y%m%d}1200.{}'.format(area,ds,ext)
            cdir = os.path.join(datadir, 'ice/conc/dm1/L4')
        elif ds.date() >= date(2016,1,1):
            # icdr
            cfile = 'ice_conc_{}_ease2-250_icdr-v2p0_{:%Y%m%d}1200.{}'.format(area,ds,ext)
            cdir = os.path.join(datadir, 'ice/conc-cont-reproc/v2p0/{}/{:%Y/%m}'.format({'nc':'','png':'quicklooks'}[ext],ds))
        else:
            # cdr
            cfile = 'ice_conc_{}_ease2-250_cdr-v2p0_{:%Y%m%d}1200.{}'.format(area,ds,ext)
            cdir = os.path.join(datadir, 'ice/conc/osi-450/{}/{:%Y/%m}'.format({'nc':'','png':'quicklooks'}[ext],ds))

        for cd in cdir:
            fn = cd + '/' + cfile
            if os.path.exists(fn):
                found = True
                break
            else:
                pass
                ##print("{}: file not found {}".format(ds,fn))

        if found:
            break

    if not found:
        raise ValueError("Did not find {} file for {} {}".format(ext,area,ds))

    return fn

def guess_hemisphere(aid):

    # which hemisphere corresponds to the grid ?
    nh_grids = ('nh','ao',)
    sh_grids = ('sh','wed','wd')

    hemis = None
    for gr in nh_grids:
        if aid.startswith(gr):
            hemis = 'nh'
            break

    if hemis is None:
        for gr in sh_grids:
            if aid.startswith(gr):
                hemis = 'sh'
                break

    if hemis is None:
        raise ValueError("Could not determine which hemisphere corresponds to {}".format(aid))

    return hemis


def remap_iceconc_file(outdir, dt, adef_o, force=False):

    # create outdir if needed
    if not os.path.exists(outdir):
        print("Create directory {} to hold ice masks".format(outdir,))
        os.mkdir(outdir,)

    # guess hemisphere
    hemis = guess_hemisphere(adef_o.area_id)

    # name of the output concentration file
    iceconc_o = os.path.join(outdir,"icemask-multi-{}-{:%Y%m%d}12.nc".format(adef_o.area_id, dt))
    if os.path.exists(iceconc_o) and not force:
        print("File {} exists already".format(iceconc_o,))
        return iceconc_o

    # find input 10km file (we use OSI-450 + OSI-430-b + DM1-SSMIS SIC files)
    iceconc_i = find_concv2_file( hemis, dt)

    print("Create ice mask from SIC file {}".format(iceconc_i))

    # read from netCDF file (and apply masking)
    with Dataset(iceconc_i) as nc:
        sic = nc.variables['ice_conc'][0,:]
        flg = nc.variables['status_flag'][0,:]

    sic_limit = 40
    edge = np.zeros(sic.shape,dtype='i2')
    edge[sic >= sic_limit] = 3
    edge[sic < sic_limit] = 1
    edge[flg == 1] = 9 # land

    adef_i = gr_utils.load_grid_defs_from_OSISAF_ncCF_file(iceconc_i)

    # resample
    edge_o = pr.kd_tree.resample_nearest(adef_i, edge, adef_o,
                                        radius_of_influence=30000,
                                        fill_value=None,reduce_data=False)
    edge_o = edge_o[np.newaxis, ...]

    # write to netCDF
    maps = dict()
    maps[adef_o.area_id] = dict()
    maps[adef_o.area_id]['_grid'] = adef_o
    maps['_what_to_map'] = ['ice_edge',]
    maps[adef_o.area_id]['ice_edge'] = edge_o
    # write the file structure to netCDF/CF
    roi_files = gr_io.out_nc(outdir, dt, 'iceconc', maps, ['ice_conc',],
                             skip_lonlat=True, nc_format='NETCDF3_CLASSIC')

    # rename
    os.rename(roi_files[0],iceconc_o)

    return iceconc_o

def preproc_img(img_file, icemask, btvar):

    # a list of all the bt channel variables in the file
    print("Preproc image {}".format(img_file))
    with Dataset(img_file) as _:
        chns = sorted ( set ([ v.split('_')[0] for v in list(_.variables.keys()) if v.startswith(btvar) ]) )

    # rename the image file if needed
    sname    = os.path.basename(img_file).split('_')
    if sname[0] == 'tc' and sname[1] != 'wght':

        all_chns = '-'.join(chns,)
        if len(sname) == 5:
            new_name = '_'.join((sname[0],'wght',sname[1],all_chns,sname[2],sname[3],sname[4]))
        else:
            new_name = '_'.join((sname[0],'wght',sname[1],all_chns,sname[2],sname[3]))
        new_name = os.path.join(os.path.dirname(img_file),new_name)
        print("Rename {} in {}".format(img_file, new_name))
        if os.path.exists(new_name):
            os.remove(new_name)
        os.rename(img_file,new_name)
        img_file = new_name

    # we must have a NetCDF3 file (bc of laplacian preprocessing)

    # check what type has the file now
    cmd = 'ncdump -k {}'.format(img_file)
    try:
        nc_type = check_output(cmd, shell=True, universal_newlines=True)
        nc_type = nc_type.rstrip("\n")
    except CalledProcessError as e:
        print("Failed with {}".format(cmd,))
        raise(e)

    if nc_type != 'classic':
        print("Convert {} file to NetCDF3.".format(nc_type))
        # convert to netcdf3 with nccopy
        img_file_nc3 = img_file.replace('.nc', '.nc3')
        cmd = 'nccopy -k classic {} {}'.format(img_file, img_file_nc3)
        try:
            check_call(cmd, shell=True,)
        except CalledProcessError as e:
            print("Failed with {}".format(cmd,))
            if os.path.exists(img_file_nc3):
                os.remove(img_file_nc3)
            raise(e)
        # now overwrite image file with the netCDF3 version
        if os.path.exists(img_file_nc3):
            os.remove(img_file)
            os.rename(img_file_nc3, img_file)

    # check if laplacian variables are already there
    if args.force:
        run_laplacian = True
    else:
        run_laplacian = False
        with Dataset(img_file) as _:
            for chn in chns:
                try:
                    _.variables[chn+'_Lap']
                except KeyError as k:
                    # missing Laplacian, must re-run
                    run_laplacian = True
                    break

    # check if the icemask file is newer than the image file
    #if ( os.path.getmtime(icemask) > os.path.getmtime(img_file) ):
    #    run_laplacian = True

    # running the laplacian pre-processing requires cleaning the
    #   netCDF file of all the Lap-related variables
    if run_laplacian:
        # Use ncks to remove the Lap-related variables
        # First, find out if some Lap-related variables are in the file
        with Dataset(img_file) as _:
            variables = _.variables.keys()
        lap_vars = []
        for _ in variables:
            if '_Lap' in _:
                lap_vars.append(_)
            if 'ice_edge' in _:
                lap_vars.append(_)
        if len(lap_vars) > 0:
            # If yes, remove them with ncks. We must use a temporary _cleaned.nc file
            img_file_cleaned = img_file.replace('.nc', '_cleaned.nc')
            cmd = 'ncks -h -x -v {} {} {}'.format(','.join(lap_vars,), img_file, img_file_cleaned)
            try:
                check_call(cmd, shell=True)
            except CalledProcessError as e:
                print("Failed with {}".format(cmd,))
                if os.path.exists(img_file_cleaned):
                    os.remove(img_file_cleaned)
                    raise(e)

            # now overwrite image file with the cleaned version
            if os.path.exists(img_file_cleaned):
                os.remove(img_file)
                os.rename(img_file_cleaned, img_file)

        # Use C software to compute the Laplacian fields
        cmd = "{}/laplacian_preprocessing -i {} -A -s {}".format(here,img_file, icemask)
        try:
            check_call(cmd, shell=True)
        except CalledProcessError:
            print("Failed with {}".format(cmd,))
    else:
        print("all laplacian variables are in place, do not re-run preprocessing")

def run_cmcc(tt0, tt1, area, swt_dir_t0, swt_dir_t1, out_dir,
             instr, channels, plot=False, tmp_dir='/tmp', update=False,
             radius=75., rad_neigh=100.):

    # get a uuid for this run. This will allow to process dates in parallel.
    tag = uuid.uuid4()

    # create the text parameter file that icedrift_solve_simplex awaits
    paramf = os.path.join(tmp_dir,'param_cmcc.{}.txt'.format(tag,))

    # Defining a log file for the C code
    logf = os.path.join(tmp_dir,'cmcc.{}.log'.format(tag,))

    # Guess name of output drift product (lev2), filename is of the form
    # icedrift_ssmis-f18_bt90v-Lap+bt90h-Lap_simplex_lev2_nh625_20130128213708w_20130130004758w.nc
    chnls = [ch.replace('_','-') for ch in channels]
    chns  = '+'.join(chnls)
    tfmt = '%Y%m%d%H%M%S'
    if tt0.time() == time(12):
        tfmt = "%Y%m%d%H"
    lev2g = 'icedrift_{}_{}_simplex_lev2_{}625_{}w_{}w.nc'.format(instr,chns,area[:2],tt0.strftime(tfmt),tt1.strftime(tfmt))
    lev2g = os.path.join(out_dir,lev2g)
    if update and os.path.exists(lev2g):
        print("Product file exists already... do not re-run ({})".format(lev2g,))
        return

# Call the python wrapper which calls the C core ice drift code
    print("\nCalling C core ice drift code with the parameters:")
    print("Start date:", tt0)
    print("End date:", tt1)
    print("Directory of start file:", swt_dir_t0)
    print("Directory of end file:", swt_dir_t1)
    print("Output directory:", out_dir)
    print("C log file:", logf)
    print("Instrument and platform:", instr)
    print("Area:", area)
    print("Channels:", channels)
    print("Radius of pattern:", radius)
    print("Radius to look for neighbours:", rad_neigh)
    fname = icedrift_wrapper(tt0, tt1, swt_dir_t0, swt_dir_t1, out_dir,
                             logf, instr, area, channels, radius, rad_neigh)

    # plot quicklooks (unless de-activated)
    if not args.no_QLs:
        zlist = {'nh':('arc',), 'sh':('ant',)}[area[0:2]]
        for z in zlist:
            qln = fname.replace('.nc','_{}.png'.format(z,))
            cmd = "python ../../format/src/ql_figure.py -d {} -o {} -r {}".format(fname, qln, z)
            try:
                check_call(cmd, shell=True,)
            except CalledProcessError as e:
                print("Failed with {}".format(cmd))
                raise(e)

    # cleanup temp files
    for tmp_f in (paramf, logf):
        if os.path.exists(tmp_f):
            #print "Clean temporary file {}".format(tmp_f)
            os.remove(tmp_f,)
            pass


if __name__ == '__main__':

    import argparse
    from argparse import RawDescriptionHelpFormatter

    p = argparse.ArgumentParser("run_cmcc",formatter_class=RawDescriptionHelpFormatter)
    apT0 = p.add_argument('START_DATE', help='datetime (YYYMMDD or YYYYMMDDHHMSS) for start image')
    apT1 = p.add_argument('STOP_DATE', help='datetime (YYYMMDD or YYYYMMDDHHMSS) for stop image. Can also be D+<offset> with <offset> a number of days after START_DATE')
    p.add_argument('INSTR', help='instrument (e.g. ssmis-f17)')
    p.add_argument('AREA', help='area code name (e.g. ao125)')
    p.add_argument('--inp_dir',help='where to look for image files.',default='.')
    p.add_argument('--reproc_instr',action='store_true',help='expect input files in INSTR/ dirs, and place results in INSTR/ ones.',)
    p.add_argument('--reproc_date',action='store_true',help='expect input files in YYYY/MM dirs, and place results in YYYY/MM/DD ones.',)
#    p.add_argument('--reproc_date',action='store_true',help='expect input files in YYYY/MM dirs, and place results in YYYY/MM ones.',)
    p.add_argument('--reproc','-R',action='store_true',help='same as --reproc_instr --reproc_date',)
    p.add_argument('--laplacian',help='Use Laplacian-preprocessed images',action='store_true')
    p.add_argument('--out_dir',help='where to write ice drift files.',default='/tmp')
    p.add_argument('--icemask_dir',help='where to write icemasks.',default=None)
    p.add_argument('--no-QLs',help='do not create quicklook images.',action='store_true')
    p.add_argument('--update',help='only run ice drift products that are missing from directory.',action='store_true')
    p.add_argument('--wave',help='wavelength to compute ice drift from',default=None)
    p.add_argument('--btvar',help='use "bt" as the brightness temp variable name rather than "tb"',action='store_true')
    p.add_argument('--radius',help='radius of the sub-window (in km, defaults to 75.0km)'.format(),default=75.,type=float)
    p.add_argument('--rad_neigh',help='radius to search for neighbours, defaults to 100km which finds 8 neighbours on the regular polstere grid', default=100.,type=float)
    p.add_argument('--tspan',help='limit the allowed tspans (in hours, defaults to no limit: all tspans are processed)',default=None)
    p.add_argument('--force',help='do not rely on existing intermediate files: redo everything',action='store_true')
    p.add_argument('--only_preproc',help='stop after the pre-processing step (ice mask, lapacian,...): no ice drift',action='store_true')
    apAS = p.add_argument('--all-swaths',help='Do swath-to-swath drift using all swaths from START_DATE and STOP_DATE.',action='store_true')
    args = p.parse_args()

    # ========================================================
    # decode / process command line parameters
    # ========================================================

    #handle the --reproc flag
    if args.reproc:
        args.reproc_instr= True
        args.reproc_date = True

    #icemask_dir
    if args.icemask_dir is None:
        imask_inp_dir   = os.path.join(args.inp_dir,'icemask')
    else:
        imask_inp_dir   = args.icemask_dir

    # inp_dir
    swath_inp_dir   = args.inp_dir
    if args.reproc_instr:
        swath_inp_dir = os.path.join(swath_inp_dir, args.INSTR)

    # decode datetime information
    t0 = None
    daily_t0 = False
    dt_fmts = ("%Y%m%d","%Y%m%d%H%M%S")
    for i_fmt, dt_fmt in enumerate(dt_fmts,):
        try:
            t0 = datetime.strptime(args.START_DATE,dt_fmt)
            if i_fmt == 0:
                daily_t0 = True
            break
        except ValueError:
            pass

    if t0 is None:
        raise argparse.ArgumentError(apT0, "Wrong format for START_DATE")

    try:
        t1 = datetime.strptime(args.STOP_DATE,dt_fmts[i_fmt])
    except ValueError:
        if args.STOP_DATE.startswith('D+'):
            try:
                offset_days = int(args.STOP_DATE[2:])
                t1 = t0 + timedelta(days=offset_days)
            except Exception:
                raise argparse.ArgumentError(apT1, "Wrong format for STOP_DATE (must be same as START_DATE or D+<days>)")
        else:
            raise argparse.ArgumentError(apT1, "Wrong format for STOP_DATE (must be same as START_DATE or D+<days>)")

    if args.reproc_date:
#        swath_inp_dir_t0 = os.path.join(swath_inp_dir,'{:%Y/%m/%d}'.format(t0))
#        swath_inp_dir_t1 = os.path.join(swath_inp_dir,'{:%Y/%m/%d}'.format(t1))
        swath_inp_dir_t0 = os.path.join(swath_inp_dir,'{:%Y/%m}'.format(t0))
        swath_inp_dir_t1 = os.path.join(swath_inp_dir,'{:%Y/%m}'.format(t1))
        imask_inp_dir_t0 = os.path.join(imask_inp_dir,'{:%Y/%m}'.format(t0))
        imask_inp_dir_t1 = os.path.join(imask_inp_dir,'{:%Y/%m}'.format(t1))
    else:
        swath_inp_dir_t0 = swath_inp_dir
        swath_inp_dir_t1 = swath_inp_dir
        imask_inp_dir_t0 = imask_inp_dir
        imask_inp_dir_t1 = imask_inp_dir


    if args.all_swaths and not daily_t0:
        raise argparse.ArgumentError(apAS, "--all_swaths works only with YYYYMMDD dates")

    if daily_t0:
        t0 = datetime.combine( t0.date(), time(12) )
        t1 = datetime.combine( t1.date(), time(12) )

    # we need also a datetime object for the day (e.g. to find the ice_conc_ file)
    t0_12 = datetime.combine( t0.date(), time(12) )
    t1_12 = datetime.combine( t1.date(), time(12) )

    # sanity check on datetimes
    if daily_t0:
        if args.all_swaths:
            if t1 < t0:
                raise argparse.ArgumentError(apT0,"STOP_DATE must come after START_DATE")
        else:
            if t1 <= t0:
                raise argparse.ArgumentError(apT0,"STOP_DATE must come after START_DATE")
    else:
        if t1 <= t0:
            raise argparse.ArgumentError(apT0,"STOP_DATE must come after START_DATE")

    # wavelengths parameters
    if args.btvar:
        btvar = 'bt'
    else:
        btvar = 'tb'
    if not args.wave:
        if args.INSTR.startswith('amsr'):
            f = (btvar + '37',)
        elif args.INSTR.startswith('cimr'):
            f = (btvar + '37',)
        elif args.INSTR.startswith('ssmi'):
            f = (btvar + '90',)
        else:
            sys.exit("no default wavelelength for instrument {}".format(args.INSTR))
    else:
        f = args.wave.split(',')


    combis = []
    for wf in f:
        if not wf.startswith(btvar):
            sys.exit("The --wave parameter should only hold {}XX (got {})".format(btvar, wf))
        vpol = wf+'v'
        hpol = wf+'h'
        if args.laplacian:
            vpol += '_Lap'
            hpol += '_Lap'
        combis.append(hpol)
        combis.append(vpol)
    combis = [combis,]
    print("Wavelengths for the sea-ice drift computations: {}".format(combis))


    # ========================================================
    # find all the image files we want to use as start/stop
    # ========================================================
    if daily_t0:
        tpatt_fmt = '{:%Y%m%d}12'
        if args.all_swaths:
            tpatt_fmt = '{:%Y%m%d}*'
    else:
        tpatt_fmt = '{:%Y%m%d%H%M%S}'

    # start day
    fpatt1 = ('tc_{}_{}_'+tpatt_fmt+'.nc').format(args.INSTR, args.AREA, t0)
    fpatt2 = ('tc_wght_{}_*_{}_'+tpatt_fmt+'.nc').format(args.INSTR, args.AREA, t0)

    for fpatt in (fpatt1, fpatt2):
        patt = os.path.join(swath_inp_dir_t0,fpatt)
        t0_files = glob.glob(patt)
        # remove the daily_map from the list if --all-swaths was used
        if args.all_swaths:
            t0_12_fext = '{:%Y%m%d}12.nc'.format(t0.date())
            for i_,_ in enumerate(t0_files,):
                if t0_12_fext in _:
                    t0_files.pop(i_)
                    break
        if len(t0_files) != 0:
            break

    if len(t0_files) == 0:
        sys.exit("Found no {} in {}".format(fpatt,swath_inp_dir_t0))

    # end day
    fpatt1 = ('tc_{}_{}_'+tpatt_fmt+'.nc').format(args.INSTR, args.AREA, t1)
    fpatt2 = ('tc_wght_{}_*_{}_'+tpatt_fmt+'.nc').format(args.INSTR, args.AREA, t1)
    for fpatt in (fpatt1, fpatt2):
        patt = os.path.join(swath_inp_dir_t1,fpatt)
        t1_files = glob.glob(patt)
        # remove the daily_map from the list if --all-swaths was used
        if args.all_swaths:
            t1_12_fext = '{:%Y%m%d}12.nc'.format(t1.date())
            for i_,_ in enumerate(t1_files,):
                if t1_12_fext in _:
                    t1_files.pop(i_)
                    break
        if len(t1_files) != 0:
            break
    if len(t1_files) == 0:
        sys.exit("found no {} in {}".format(fpatt,swath_inp_dir_t1))

    # ========================================================
    # pre-processing of the images (Laplacian, icemask)
    # ========================================================

    # Check we have an ice_mask file, otherwise create it
    #    we need an area_def for the target grid
    adef = gr_utils.load_grid_defs_from_OSISAF_ncCF_file(t0_files[0])
    adef.area_id = args.AREA
    icemasks = dict()
    icemasks[t0_12] = remap_iceconc_file(imask_inp_dir, t0_12, adef, force=args.force)
    if t1_12 != t0_12:
        icemasks[t1_12] = remap_iceconc_file(imask_inp_dir, t1_12, adef, force=args.force)

    # Run through the t0 files and apply the pre-processing filters (e.g. Laplacian)
    for f in t0_files :
        preproc_img(f, icemasks[t0_12], btvar)

    # pre-process only the t1 files that were not already processed as part of t0_files
    #     (happens if start and stop dates are the same)
    for f in list(set(t1_files).difference(t0_files)):
        preproc_img(f, icemasks[t1_12], btvar)

    #if --only_preproc was used, we are done
    if args.only_preproc:
        print("Stop here, since --only_preproc was used")
        sys.exit(0)

    # ========================================================
    # sea-ice drift computations
    # ========================================================

    # Extract timestamps from the filenames
    if daily_t0 and not args.all_swaths:
        t0_timestamps = sorted([datetime.strptime((os.path.basename(f).split('_')[-1])[:-3],'%Y%m%d%H') for f in t0_files])
        t1_timestamps = sorted([datetime.strptime((os.path.basename(f).split('_')[-1])[:-3],'%Y%m%d%H') for f in t1_files])
    else:
        t0_timestamps = sorted([datetime.strptime((os.path.basename(f).split('_')[-1])[:-3],'%Y%m%d%H%M%S') for f in t0_files])
        t1_timestamps = sorted([datetime.strptime((os.path.basename(f).split('_')[-1])[:-3],'%Y%m%d%H%M%S') for f in t1_files])

    # Cartesian product of the two lists is what we want to run as ice drift
    t0t1s = list(itertools.product(t0_timestamps,t1_timestamps))

    # screen/filter on a minimum separation time (also that t0 must be before t1)
    list_t0t1s = []

    if args.tspan is None:
        # all tspans are allowed
        for tt0,tt1 in t0t1s:
            if tt1 > tt0:
                list_t0t1s.append((tt0,tt1))
    else:
        tsp = float(args.tspan)
        tspans = (tsp-1.,tsp+1.)
        for tt0,tt1 in t0t1s:
            if tt1 > tt0:
                tdiff_h = (tt1-tt0).total_seconds() / 60. / 60.
                if tdiff_h>=tspans[0] and tdiff_h<=tspans[1]:
                    list_t0t1s.append((tt0,tt1))

    if len(list_t0t1s) == 0:
        sys.exit("Found no suitable start/stop image pairs.")

    # ready to run: prepare output directories
    d = args.out_dir
    if args.reproc_instr:
        d = os.path.join(d,args.INSTR)
    if not os.path.exists(d):
        print("Make output directory {}".format(d,))
        os.makedirs(d)

    # start looping over the start/stop pairs, and run ice drift processing
    for tt0,tt1 in list_t0t1s:

        outd = d
        if args.reproc_date:
            outd = os.path.join(d,'{:%Y/%m/%d}'.format(tt1))
        if not os.path.exists(outd):
            os.makedirs(outd)

        for ch in combis:
            print("Compute sea-ice drift for {}".format(ch))
            run_cmcc(tt0, tt1, args.AREA, swath_inp_dir_t0, swath_inp_dir_t1,
                     outd, args.INSTR, ch, tmp_dir=args.out_dir,
                     update=args.update, radius=args.radius,
                     rad_neigh=args.rad_neigh)


    sys.exit('Done with run_cmcc.py')
