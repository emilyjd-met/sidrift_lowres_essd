#!/usr/bin/env python3

import os
import sys
from netCDF4 import Dataset, num2date, date2num
import argparse
import numpy
from datetime import date, datetime, timedelta
from add_uncert_dict import ud, fitDT0dict

sys.path.append(os.path.join(os.path.dirname(__file__), '../../merge/src'))
from merge_periods import merge_season, merge_frac


def get_season(dt, area):

    if area == 'nh':
        ret = 'winter'
        # summer = may, jun, jul, aug, sep
        if dt.month >= 6 and dt.month <= 8: ret = 'summer'
        elif dt.month == 5: ret = 'spring'
        elif dt.month == 9: ret = 'autumn'
        return ret
    elif area == 'sh':
        ret = 'winter'
        # summer = nov, dec, jan, feb, mar
        if dt.month >= 12 or dt.month <= 2: ret = 'summer'
        elif dt.month == 11: ret = 'spring'
        elif dt.month == 3: ret = 'autumn'
        return ret
    else:
        raise NotImplementedError("Only done for NH and SH")

def timedelta_total_seconds(td):
    """ Implements timedelta.total_seconds() as available from Py2.7 """
    try:
        ret = td.total_seconds()
    except NameError:
        ret = (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
    return ret

def get_uncerts_fitDT0(prod, mode='nrt'):

    if mode == 'nrt':
        if prod.startswith('amsr'):
            fitDT0coeffs = fitDT0dict['nrt']['amsr']
        elif prod.startswith('ssmi'):
            fitDT0coeffs = fitDT0dict['nrt']['ssmi']
        elif prod.startswith('ascat'):
            fitDT0coeffs = fitDT0dict['nrt']['ascat']
        else:
            raise ValueError('Know nothing about the uncertainties of {}'.format(prod))
    elif mode == 'cdr':
        fitDT0coeffs = fitDT0dict['cdr']['all']

    return numpy.poly1d(fitDT0coeffs)


def get_uncerts_flag(prod, area='nh',season='winter', mode='nrt'):

    # Finds and returns a dictionary as follows
    # {0:uncert_nominal,    # nominal
    #  13:uncert_correct,   # correct by neighbours
    #  17:uncert_smaller,   # smaller pattern
    #  30:uncert_nominal}   # wind icedrift from gapfilled param file

    prodind = prod
    if prod.startswith('ssmi'):
        prodind = 'ssmi'
    elif prod.startswith('amsr'):
        prodind = prodind.replace('tb', 'bt')
    elif prod.startswith('ascat'):
        prodind = 'ascat'
    elif prod.startswith('wind') or prod.startswith('ecmwf'):
        prodind = 'wind'

    try:
        uncert_dict = ud[mode][area][prodind][season]
    except:
        raise ValueError("No entry in uncertainty dictionary with mode {} - area {} - prod index {} - season {}: ".format(mode, area, prodind, season))

    return uncert_dict


def raise_uncertainty_with_dt0(t0,ref_t0,product, mode='nrt'):

    # Prepare the polynomial fit
    pol = get_uncerts_fitDT0(product, mode=mode)

    # Compute dt0 in hours (the sign does not matter)
    dt0s = numpy.ma.empty(t0.shape)
    dt0s.mask = t0.mask
    if dt0s.count() == 0:
        corr = numpy.ma.zeros(t0.shape)
        corr.mask = t0.mask
    else:
        dt0s[~t0.mask] = date2num(t0[~t0.mask],'hours since {}'.format(ref_t0.isoformat()))

        # restrict to validity domain of the polynomial fit
        dt0s = numpy.clip(dt0s,-10.,+10.)

        # Estimate the polynomial fit correction
        corr = pol(dt0s)

    return corr

def add_uncertainty_fields(fname, mode='nrt', method='ramp'):

    # open the file for read/write
    try:
        grp = Dataset(fname,"a")
    except Exception as e:
        raise ValueError("Cannot open file {} with 'writing' mode ({})".format(fname,e))

    # check it is a 'proc' file
    try:
        d = grp.variables['driftX']
    except KeyError as k:
        raise ValueError("File does not seem to be a 'proc' LRSID file (does not contain {})".format(k,))

    # guess which instrument sensor it was processed from (using the filename)
    fn = os.path.basename(fname,)
    fn_toks = fn.split('_')
    prod = fn_toks[1]
    if prod == 'multi-oi':
        raise ValueError("Can only process single-sensor LRSID files (got a {} one)".format(prod,))
    # Create a product name with wavebands (except for wind drift where
    # product is just 'wind')
    if prod != 'wind' and prod != 'ecmwf-era5':
        wbands = fn_toks[2] # (bt19v-Lap+bt19h-Lap)
        wband = wbands.split('-')[0] # (bt19v)
        if wband.startswith('bt') or wband.startswith('tb'):
            wband = wband[:-1] # remove the polarization
        prod = "{}-{}-{}".format(grp.instrument,grp.platform,wband)

#    area = fn_toks[5][:2]
    if 'nh' in [x[:2] for x in fn_toks]:
        area = 'nh'
    elif 'sh' in [x[:2] for x in fn_toks]:
        area = 'sh'
    else:
        raise ValueError("No recognised area in filename {}".format(fname))

    # load the status_flag, and driftX, and time
    sflag = grp.variables['flag'][:]
    dX = grp.variables['driftX'][:]
    time_t0 = num2date(grp.variables['time_bnds'][0,0],
                       grp.variables['time_bnds'].units)
    t0 = numpy.ma.array(num2date(grp.variables['t0'][0,:],
        'seconds since {}'.format(time_t0.isoformat())))
    t0.mask = dX.mask[0,:]

    if mode == 'nrt':
        pass
        # TODO: Implement this
        #season=get_season(time_season, area)
    elif mode == 'cdr':
        date_season = num2date(grp.variables['time_bnds'][0,1],
                               grp.variables['time_bnds'].units).strftime('%Y%m%d%H')
        date_season = datetime.strptime(date_season, '%Y%m%d%H')
        season, _ = merge_season(date_season, area)

    uncert_values = dict()
    if mode == 'nrt':
        if season == 'winter' or season == 'summer':
            # This will get the uncertainties right for both summer and winter
            uncert_values = get_uncerts_flag(prod, area=area, season=season,
                                             mode=mode)
        else:
            # Transition between the two levels
            winter_uncert = get_uncerts_flag(prod, area=area, season='winter',
                                             mode=mode)
            summer_uncert = get_uncerts_flag(prod, area=area, season='summer',
                                             mode=mode)
            if season == 'early_summer':
                u1 = winter_uncert
                u2 = summer_uncert
            elif season == 'late_summer':
                u1 = summer_uncert
                u2 = winter_uncert
            else:
                raise ValueError("Unknown season {}".format(season,))

            w = time_t0.day / 30.
            for k in list(winter_uncert.keys()):
                uncert_values[k] = (1. - w) * u1[k] + w * u2[k]
            #print season,w,k,u1[k],u2[k],uncert_values[k]

    elif mode == 'cdr':

        # These are fixed spring and autumn season uncertainties
        # (NOT USED)
        if method == 'fixseason':
            uncert_values = get_uncerts_flag(prod, area=area, season=season,
                                             mode=mode)
        # In winter, the wind values are set to 'veryhigh' so they do not
        # contribute, and in summer the same for the satellite data.
        elif method == 'ramp':

            windprod = False
            nomseas = 'winter'
            if prod.startswith('wind') or prod.startswith('ecmwf'):
                windprod = True
                nomseas = 'summer'

            nomflags = get_uncerts_flag(prod, area=area, season=nomseas,
                                        mode=mode)

            # In winter the wind does not contribute (unless no satellite data)
            if season == 'winter':
                if windprod:
                    uncert_values = ud['cdr']['veryhigh']['wind']
                else:
                    uncert_values = nomflags

            # In summer, the satellite data is not used
            elif season == 'summer':
                if windprod:
                    uncert_values = nomflags
                else:
                    uncert_values = ud['cdr']['veryhigh']['sat']

            # In spring and autumn there is a transition between two levels,
            # 'high' and the nominal value
            elif season == 'spring' or season == 'autumn':
                if season == 'spring':
                    if windprod:
                        u1 = ud['cdr']['high']['wind']
                        u2 = nomflags
                    else:
                        u1 = nomflags
                        u2 = ud['cdr']['high']['sat']
                elif season == 'autumn':
                    if windprod:
                        u1 = nomflags
                        u2 = ud['cdr']['high']['wind']
                    else:
                        u1 = ud['cdr']['high']['sat']
                        u2 = nomflags

                # Merge fraction uses the exact number of days in the month
                # Also different for NRT, the date used is the end date
                # not the start date
                wfrac, sfrac = merge_frac(date_season, area)
                if season == 'spring':
                    frac = sfrac
                elif season == 'autumn':
                    frac = wfrac
                for k in list(nomflags.keys()):
                    uncert_values[k] = (frac * (u2[k] - u1[k])) + u1[k]

    # prepare an array for the uncertainties based on the flag values (same mask as dX)
    uncert = numpy.ma.empty_like(dX)
    uncert.data[~uncert.mask] = -1
    for flgs in list(uncert_values.keys()):
        uncert[sflag==flgs] = uncert_values[flgs]

    # check that all the flags were assigned
    if (uncert==-1).sum() > 0:
        sys.stderr.write("WARNING: could not set all the uncertainty values: some flags were un-expected")

    # prepare a new array for the uncertainties to be used for 12utc drift (e.g. NN2Dcol, multi-oi,...)
    uncert12 = uncert.copy()
    # Note: wind drift is not changed by this, as it is 12UTC
    uncert12 += raise_uncertainty_with_dt0(t0,time_t0,prod, mode=mode)

    # write the new variables to file
    nameX = 'fsX'
    nameY = 'fsY'
    if not nameX in list(grp.variables.keys()):
        grp.createVariable(nameX,dX.dtype,('time','yc','xc'),fill_value=dX.fill_value)
        setattr(grp.variables[nameX],'long_name','Uncertainties (1 stddev) of driftX, based on flag values')
        setattr(grp.variables[nameX],'grid_mapping',grp.variables['driftX'].grid_mapping)
        setattr(grp.variables[nameX],'coordinates',grp.variables['driftX'].coordinates)

    if not nameY in list(grp.variables.keys()):
        grp.createVariable(nameY,dX.dtype,('time','yc','xc'),fill_value=dX.fill_value)
        setattr(grp.variables[nameY],'long_name','Uncertainties (1 stddev) of driftY, based on flag values')
        setattr(grp.variables[nameY],'grid_mapping',grp.variables['driftX'].grid_mapping)
        setattr(grp.variables[nameY],'coordinates',grp.variables['driftX'].coordinates)

    grp.variables[nameX].set_auto_maskandscale(True)
    grp.variables[nameX][0,:] = uncert.copy()
    grp.variables[nameY].set_auto_maskandscale(True)
    grp.variables[nameY][0,:] = uncert.copy()

    nameX_12 = 'fsX_12utc'
    nameY_12 = 'fsY_12utc'
    if not nameX_12 in list(grp.variables.keys()):
        grp.createVariable(nameX_12,dX.dtype,('time','yc','xc'),fill_value=dX.fill_value)
        setattr(grp.variables[nameX_12],'long_name','Uncertainties (1 stddev) of driftX, based on flag values')
        setattr(grp.variables[nameX_12],'grid_mapping',grp.variables['driftX'].grid_mapping)
        setattr(grp.variables[nameX_12],'coordinates',grp.variables['driftX'].coordinates)

    if not nameY_12 in list(grp.variables.keys()):
        grp.createVariable(nameY_12,dX.dtype,('time','yc','xc'),fill_value=dX.fill_value)
        setattr(grp.variables[nameY_12],'long_name','Uncertainties (1 stddev) of driftY, based on flag values')
        setattr(grp.variables[nameY_12],'grid_mapping',grp.variables['driftX'].grid_mapping)
        setattr(grp.variables[nameY_12],'coordinates',grp.variables['driftX'].coordinates)

    grp.variables[nameX_12].set_auto_maskandscale(True)
    grp.variables[nameX_12][0,:] = uncert12.copy()
    grp.variables[nameY_12].set_auto_maskandscale(True)
    grp.variables[nameY_12][0,:] = uncert12.copy()


    grp.close()

if __name__ == '__main__':

    import sys, traceback

    valid_modes = ['nrt', 'cdr']
    valid_methods= ['ramp', 'fixseason']
    p = argparse.ArgumentParser('add_uncertainties.py')
    p.add_argument('INPF',help="Path to a 'proc' ice drift file")
    p.add_argument('-m', '--mode', default='nrt', choices=valid_modes,
                   help="Mode of processing, default is 'nrt', valid modes are {}".format(valid_modes))
    p.add_argument('-q', '--method', default='ramp', choices=valid_methods,
                   help="Method to calculate uncertainty for CDR, default is 'ramp', valid methods are {}".format(valid_methods))
    args = p.parse_args()

    try:
        add_uncertainty_fields(args.INPF, mode=args.mode, method=args.method)
    except Exception as ex:
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        print('-'*60)
        sys.exit("Failed with {}".format(args.INPF,))
