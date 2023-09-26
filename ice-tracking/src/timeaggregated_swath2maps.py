#!/usr/bin/env python3

import sys
import os
from glob import glob
from datetime import datetime, date, time
import numpy as np
import numpy.ma as ma

import pyresample as pr
import netCDF4
from netCDF4 import Dataset, date2num
from math import sqrt

def read_swaths(swath,varn,):

    # check types
    if isinstance(varn, str):
        varn = [varn,]
    if isinstance(swath, str):
        swath = [swath,]

    # read all swath files, appending
    ret = dict()
    for sw in swath:
        # get read access to swath file
        try:
            grp = Dataset(sw,"r")
        except Exception as ex:
            print("Cannot open netcdf {} for reading. Skip it. {}".format(sw,ex,))
            continue

        for var in varn:
            # get read access to variable
            try:
                varo = grp.variables[var]
                varo.set_auto_maskandscale(True)
            except KeyError as e:
                raise ValueError("Variable {} not found in {}. Cannot recover.".format(var,sw,))

            # read (force masked array)
            dat = varo[:]
            try:
                dat.mask
            except AttributeError:
                dat = ma.array(dat,mask=np.zeros(dat.shape).astype('bool'))

            # store (append or new array)
            try:
                ret[var] = ma.array(np.append(ret[var],dat),mask=np.append(ret[var].mask,dat.mask))
            except KeyError:
                ret[var] = dat

        # close swath file
        grp.close()

    # sanith check for what we read
    shape = ret[varn[0]].shape
    for var in varn[1:]:
        if ret[var].shape != shape:
            raise ValueError("Shape mismatch for {} (got {} but expected {})",format(var,ret[var].shape,shape))

    # internal, book-keeping data begin with '_'
    ret['_files'] = swath
    ret['_shape'] = shape
    ret['_varn']  = varn
    return ret

def plot_sig_loc(testP,tg,all_inc,all_sig,all_azim=None,all_kp=None,ref_ang=40.):
    from matplotlib import pylab as plt

    def fit_corr(ang, sig, ref_ang=45.):
        lin_fit = np.polyfit(ang,sig,1)
        sig_corr = np.polyval(lin_fit,ref_ang)
        return sig_corr, lin_fit[0],np.polyval(lin_fit,ang)

    def fixed_corr(ang, sig, ref_ang=45.):
        dist = ref_ang - ang
        slope = -0.19
        sig_corr = sig + (-0.19 * dist)
        closest = np.argmin(abs(dist))
        return sig_corr[closest], slope

    all_inc = all_inc.compressed()
    all_sig = all_sig.compressed()
    if all_kp is not None:
        all_kp = all_kp.compressed()
    if all_azim is not None:
        all_azim = all_azim.compressed()

    if all_kp is None:
        plt.scatter(all_inc,all_sig,marker='o',facecolor='none',)
    else:
        plt.errorbar(all_inc,all_sig,fmt='o',yerr=all_kp,color='k')

    ax = plt.gca()
    xmin = 25.
    xmax = 65.

    sig_corr_fixed, slope_fixed = fixed_corr(all_inc,all_sig,ref_ang=ref_ang)
    sig_corr_fit, slope_fit, all_sigs_fit = fit_corr(all_inc,all_sig,ref_ang=ref_ang)

    plt.scatter(ref_ang,sig_corr_fixed,marker='x',color='red')
    plt.scatter(ref_ang,sig_corr_fit,marker='+',color='blue')
    plt.scatter(all_inc,all_sigs_fit,marker='+',color='k')
    if all_azim is not None:
        for iaz,az in enumerate(all_azim):
            plt.text(all_inc[iaz]+0.25,all_sig[iaz],'{:.1f}'.format(az,),fontsize=8)
    ax.set_xlabel('Incidence Angle')
    ax.set_xlim(xmin,xmax)
    ax.set_ylabel(r'$\sigma^0$ [dB]')
    ax.set_title('{} {}'.format(tg.area_id,testP,))
    ax.text(0.05,0.95,'slope = {:.3f}'.format(slope_fit,),transform=ax.transAxes)
    ax.text(0.05,0.90,r'$\sigma_c$ = {:.3f} dB'.format(sig_corr_fit,),transform=ax.transAxes)
    oname = '/tmp/test.png'
    plt.savefig(oname,bbox_inches='tight')
    print('{} is ready'.format(oname,))

def map_swaths(data,grid_defs):

    # find and apply a common mask, remove all missing data, and flatten
    common_mask = np.zeros(data['_shape']).astype('bool')
    for w in data['_varn']:
        common_mask = np.logical_or(common_mask,data[w].mask)

    for w in data['_varn']:
        data[w].mask = common_mask
        data[w] = data[w].compressed()
        #print w, data[w].min(), data[w].max()


    # transform data into a staggered array
    data_arr = np.empty((data[w].size,len(data['_what_to_map'])))
    for iw,w in enumerate(data['_what_to_map']):
       data_arr[:,iw] = data[w]

    # prepare swath for pyresample remapping
    swath_def = pr.geometry.SwathDefinition(lons=data['lon'], lats=data['lat'])

    sig_wf = 0.75 * 12500
    wf = lambda r: np.exp(-0.5*(r/sig_wf)**2)
    wf = [wf,]*len(data['_what_to_map'])

    # remap to target grids
    maps = dict()
    for tg in grid_defs:
        maps[tg.area_id] = dict()
        valid_input_index, valid_output_index, index_array, distance_array = \
                pr.kd_tree.get_neighbour_info(swath_def,tg,radius_of_influence=12500.*sqrt(2.),neighbours=50)

        # from inside pr.kdtree.get_sample_from_neighbour_info()
        # ==================================================================================
        output_shape = tg.shape
        for iw,w in enumerate(data['_what_to_map']):
            if data[w].ndim > 2 and data[w].shape[0] * data[w].shape[1] == valid_input_index.size:
                data[w] = data[w].reshape(data[w].shape[0] * data[w].shape[1], data[w].shape[2])
                print("reshaped {} to {}".format(w,data[w].shape))
        valid_input_size = valid_input_index.sum()
        valid_output_size = valid_output_index.sum()
        input_size = valid_input_size
        if len(output_shape) > 1:
            output_size = output_shape[0] * output_shape[1]
        else:
            output_size = output_shape[0]
        neighbours = index_array.shape[1]
        new_data = dict()
        for iw,w in enumerate(data['_what_to_map']):
            new_data[w] = data[w][valid_input_index]

        #Get neighbours and masks of valid indices
        ch_neighbour_list = []
        index_mask_list = []
        all_sig = np.empty((output_size,3*neighbours))
        all_inc = np.empty((output_size,3*neighbours))
        all_kp  = np.empty((output_size,3*neighbours))
        all_msk = np.empty((output_size,3*neighbours)).astype('bool')
        for i in range(neighbours): # Iterate over number of neighbours
            # Make working copy neighbour index and 
            # set out of bounds indices to False
            index_ni = index_array[:, i].copy()
            index_mask_ni = (index_ni == input_size)
            index_ni[index_mask_ni] = False

            # Get channel data for the corresponing indices
            for iv,view in enumerate(('fore','mid','back'),):
                all_sig[:,3*i+iv] = new_data['sigma0_'+view][index_ni]
                all_kp[:,3*i+iv]  = new_data['kp_'+view][index_ni]
                all_inc[:,3*i+iv] = new_data['inc_ang_'+view][index_ni]
                all_msk[:,3*i+iv] = index_mask_ni

        # ==================================================================================

        all_sig = ma.array(all_sig,mask=all_msk)
        all_kp  = ma.array(all_kp,mask=all_msk)
        all_inc = ma.array(all_inc,mask=all_msk)

        # http://stackoverflow.com/questions/28237428/fully-vectorise-numpy-polyfit 
        n = (~all_sig.mask).sum(axis=1)
        x = all_inc
        y = all_sig
        slp_corr = (n*(x*y).sum(axis=1) - x.sum(axis=1)*y.sum(axis=1)) / (n*(x*x).sum(axis=1) - x.sum(axis=1)*x.sum(axis=1))
        int_corr = y.mean(axis=1) - slp_corr * x.mean(axis=1)
        sig_corr = slp_corr * 40.0 + int_corr
        sig_sdev = y.std(axis=1)
        nb_corr  = n / 3.

        nb_corr  = nb_corr.reshape(tg.shape)
        sig_corr = sig_corr.reshape(tg.shape)
        slp_corr = slp_corr.reshape(tg.shape)
        sig_sdev = sig_sdev.reshape(tg.shape)
        #kp_corr  = kp_corr.reshape(tg.shape)

        sig_corr = ma.array(sig_corr,mask=(nb_corr<2))
        slp_corr = ma.array(slp_corr,mask=(nb_corr<2))
        sig_sdev = ma.array(sig_sdev,mask=(nb_corr<2))
        #kp_corr  = ma.array(kp_corr,mask=(nb_corr==0))

        # transfer back in the ouput dict
        maps[tg.area_id]['sigma_40'] = sig_corr
        maps[tg.area_id]['sigma_slope'] = slp_corr
        #maps[tg.area_id]['kappa_40'] = kp_corr
        maps[tg.area_id]['sigma_nb'] = nb_corr
        maps[tg.area_id]['sigma_40_flg'] = np.where(maps[tg.area_id]['sigma_nb']==0,-1,0)
        maps[tg.area_id]['sigma_40_sdev'] = sig_sdev
        maps[tg.area_id]['sigma_40_dtime'] = sig_corr.copy()
        maps[tg.area_id]['sigma_40_dtime'].data[:] = 12.
        maps[tg.area_id]['_grid'] = tg

    # return
    return maps
    
def out_nc_one(outn,day,instr,mapo):
    grp = Dataset(outn,"w",format='NETCDF3_CLASSIC')

    areadef = mapo['_grid']
    lon, lat = areadef.get_lonlats()
    midday = datetime.combine(day,time(12,))
    
    # create Dimensions
    grp.createDimension('time',1)
    grp.createDimension('xc',areadef.x_size)
    grp.createDimension('yc',areadef.y_size)

    # create variables associated with dimensions

    # - time
    varn = grp.createVariable('time','f8',('time',),)
    time_unit = 'seconds since 1970-01-01 00:00:00'
    setattr(varn, 'axis', 'T')
    setattr(varn, 'units', time_unit)
    setattr(varn, 'calendar', 'standard')
    setattr(varn, 'standard_name', 'time')
    setattr(varn, 'long_name', 'reference time of map')
    varn[:] = date2num(midday,units=time_unit)

    # - xc
    varn = grp.createVariable('xc','f8',('xc',))
    setattr(varn,'axis','X')
    setattr(varn,'units','km')
    setattr(varn,'long_name','x coordinate of projection (eastings)')
    setattr(varn,'standard_name','projection_x_coordinate')
    varn[:] = areadef.projection_x_coords / 1000.

    # - yc
    varn = grp.createVariable('yc','f8',('yc',))
    setattr(varn,'axis','Y')
    setattr(varn,'units','km')
    setattr(varn,'long_name','y coordinate of projection (eastings)')
    setattr(varn,'standard_name','projection_y_coordinate')
    varn[:] = areadef.projection_y_coords / 1000.

    # - crs
    varn = grp.createVariable('Polar_Stereographic_Grid','i4')
    proj_str = ' '.join(["+%s=%s"%(str(k), str(areadef.proj_dict[k]))
                               for k in sorted(areadef.proj_dict.keys())])
    setattr(varn,'proj4_string',proj_str)

    # create all other variables

    # - lat
    varn =  grp.createVariable('lat','f4',('yc','xc',))
    setattr(varn,'units','degrees_north')
    setattr(varn,'long_name','latitude coordinate')
    setattr(varn,'standard_name','latitude')
    varn[:] = lat

    # - lon
    varn =  grp.createVariable('lon','f4',('yc','xc',))
    setattr(varn,'units','degrees_east')
    setattr(varn,'long_name','longitude coordinate')
    setattr(varn,'standard_name','longitude')
    varn[:] = lon

    # variables
    for k in list(mapo.keys()):
        if k[0] == '_': continue
        if k == 'sigma_40':
            varname = 'sigma0_TCobs'
            dtype = 'f4'
            units = 'dB'
        elif k == 'sigma_40_flg':
            varname = 'sigma0_TCobs_flag'
            dtype = 'i2'
            units = '1'
        elif k == 'sigma_40_sdev':
            varname = 'sigma0_TCdev'
            dtype = 'f4'
            units = 'dB'
        elif k == 'sigma_40_dtime':
            varname = 'sigma0_avgTime'
            dtype = 'f4'
            units = 'hours since start_time'
        else:
            continue

        varn = grp.createVariable(varname,dtype,('time','yc','xc',),
                                  fill_value=netCDF4.default_fillvals[dtype])
        setattr(varn,'grid_mapping','crs')
        setattr(varn,'coordinates','lat lon')
        setattr(varn,'units',units)
        varn[0,:] = mapo[k].astype(dtype)

    grp.close()

def out_nc(outdir,day,instr,maps,tweight=False):
        
    twstr = 'nowght'
    if tweight:
        twstr = 'wght'

    outns = []
    for aid in list(maps.keys()):
        outn = 'tc_{}_{}_sigma0_{}_{}12.nc'.format(twstr,instr,aid,day.strftime('%Y%m%d'),)
        outn = os.path.join(outdir,outn,)
        #try:
        out_nc_one(outn,day,instr,maps[aid])
        outns.append(outn,)
        #except Exception as ex:
        #    print "Something wrong while writing {} ({})".format(outn,ex,)
        #    outns.append('error',)
        #    continue

    return outns

if __name__ == '__main__':

    import argparse

    grid_def_file = './grids.def'

    # Parse command line parameters 
    p = argparse.ArgumentParser()
    p.add_argument('DATE',help='Gather all swath for that day (YYYYMMDD)')
    p.add_argument('INSTR',help='Intrument and platform (e.g. ascat-metopA)')
    p.add_argument('GRIDS',help='Name of target grids (see --grid-def-file)')
    p.add_argument('INDIR',help='Directory to search for swath files (see -R)')
    p.add_argument('OUTDIR',help='Directory to write the gridded map into')
    p.add_argument('-R',dest='REPROC',action='store_true',
                   help='If set, will search swath files in <indir>/YYYY/MM/DD sub-directories')
    p.add_argument('--grid_def_file',dest='gdef',default=grid_def_file,
                   help='Alternative grid definition file (defaults to {})'.format(grid_def_file,))
    p.add_argument('-L',dest='logf',default=None,
                   help='Log file for the run')
    p.add_argument('-T',dest='TWEIGHT',action='store_true',
                   help='Currently does nothing but modify name of output file')
    args = p.parse_args()

    if args.logf is not None:
        sys.stdout = open(args.logf, 'w')

    # sanity check
    try:
        datestr = args.DATE
        if len(datestr) > len('YYYYMMDD'):
            datestr = datestr[:len('YYYYMMDD')]
        DATE = datetime.strptime(datestr,"%Y%m%d").date()
    except ValueError as ve:
        sys.exit("Error in decoding datestring {} ({})".format(args.DATE,ve,))

    try:
        instr, platf = args.INSTR.split('-')
    except ValueError:
        sys.exit('Invalid INSTR value (expect <instr>-<platform>): {}'.format(args.INSTR,))

    # read grid definition file
    grid_defs = [pr.utils.parse_area_file(args.gdef,tg)[0] for tg in args.GRIDS.split(',')]

    # per instrument/sensor configuration
    if instr == 'ascat':
        platf = platf.lower()
        allowed_ascat = ('metopa','metopb','metopc',)
        if platf not in allowed_ascat:
            sys.exit('Invalid INSTR value: ascat is only on-board {}'.format(allowed_ascat,))
        patt = 'ascat_{}_*_{}_*_125.nc'.format(DATE.strftime('%Y%m%d'),platf,)
        what_to_map = ['sigma0_fore','sigma0_mid','sigma0_back','inc_ang_fore','inc_ang_mid','inc_ang_back',
                      'azim_ang_fore','azim_ang_mid','azim_ang_back','kp_fore','kp_mid','kp_back']
        aux_to_map = ['lat','lon']
        fit_in_dB = True
    else:
        sys.exit('Know nothing about instrument {}'.format(instr,))

    # Find all swath files
    if args.REPROC:
        patt = '{}/{}'.format(DATE.strftime('%Y/%m/%d'),patt,)

    swath = glob(os.path.join(args.INDIR,patt))
    if len(swath) == 0:
        sys.exit("Found no swath for {} {}".format(args.DATE,args.INSTR,))

    # Read
    data = read_swaths(swath,what_to_map+aux_to_map)
    data['_what_to_map'] = what_to_map
    data['_aux_to_map']   = aux_to_map

    # Pre-process
    for view in ('fore','mid','back'):
        if not fit_in_dB:
            data['sigma0_'+view] = 10. ** (data['sigma0_'+view] * 0.1)
        data['kp_'+view] = 0.01 * data['kp_'+view] * data['sigma0_'+view]

    # Remap and time-composite
    maps = map_swaths(data,grid_defs,)

    # Post-process
    for tg in grid_defs:
        if not fit_in_dB:
            maps[tg.area_id]['sigma_40_lin'] = maps[tg.area_id]['sigma_40']
            maps[tg.area_id]['sigma_40'] = 10. * np.log10(maps[tg.area_id]['sigma_40_lin'])
        try:
            maps[tg.area_id]['kappa_40'] = maps[tg.area_id]['kappa_40'] / maps[tg.area_id]['sigma_40']
            if not fit_in_dB:
                maps[tg.area_id]['kappa_40_lin'] = maps[tg.area_id]['kappa_40']
                maps[tg.area_id]['kappa_40'] = 100. * 10. * np.log10(maps[tg.area_id]['kappa_40_lin'] / maps[tg.area_id]['sigma_40_lin'])
        except KeyError:
            pass

    # Write output files (1 per target grid)
    outs = out_nc(args.OUTDIR,DATE,args.INSTR,maps,tweight=args.TWEIGHT)

    # Conclude
    print("OUTF<{}>".format(outs[0],))
    print("Wrote {} remapped file(s): {}".format(len(outs,),outs))


