import sys
import matchups_io as mio
import matchups_selector as msel
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
import find_matchups_files as mfil
import scipy as sp
from scipy import ma
from scipy.stats import t as statst
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import cm
from matplotlib.patches import Ellipse
from datetime import datetime, date, timedelta, time
from math import sqrt, ceil
from collections import Counter
from itertools import chain
import cmocean
from sklearn import linear_model, datasets

from matplotlib import dates as mdates
from matplotlib import gridspec

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# global variables to set labels
tstep  = 0.05
tstart = 0.95
tx     = 0.03
fs     = 14

#cX = 'red'
#cY = 'blue'
cX = '#4363d8'
cY = '#ffe119'

dpi = 110

def titled(mydatelist, pos=0):

    try:
        titledate = mydatelist[pos].isoformat()
    except:
        titledate = mydatelist[0][pos].isoformat()

    return titledate


def timedelta_total_seconds(td):
    """ Implements timedelta.total_seconds() as available from Py2.7 """
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

color_scheme_validation_source  = {'gps': 'magenta', 'argos': 'red', 'wsm': 'green'}
#color_scheme_validation_network = {'ITP': 'magenta',
#                                   'AARI': 'blue',
#                                   'IABP': 'lightgrey',
#                                   'Tara-GPS': 'orange',
#                                   'BBArray': 'lightgrey',
#                                   'DTU': 'lightblue',
#                                   'SAMS': 'green',
#                                   'CRREL': 'teal',
#                                   'SIP': 'slateblue'}
color_scheme_validation_network = {'AARI': 'Red',
                                   'AAR': 'Red',
                                   'ANT': 'Maroon',
                                   'ARG': 'Orange',
                                   'BBArray': 'Yellow',
                                   'BBB': 'Yellow',
                                   'CRREL': 'Lime',
                                   'CRR': 'Lime',
                                   'DTU': 'Mint',
                                   'HUD': 'Olive',
                                   'IABP': 'Green',
                                   'IAB': 'Green',
                                   'ITP': 'Blue',
                                   'SAMS': 'Cyan',
                                   'SAM': 'Cyan',
                                   'SED': 'Purple',
                                   'SIM': 'Lavender',
                                   'SIP': 'Magenta',
                                   'Tara-GPS': 'Teal',
                                   'TAR': 'Teal'}



flag_dict = {20:(17,'small_patt'),
             21:(13,'corr_neighb'),
             22:(16,'interp.'),
             30:(0,'nominal'),}


def get_color(m, scheme='source', start_period=None, end_period=None):
    """
    Get colours for plotting matchups, depending on some colour scheme
    """

    legend = None
    cblab   = None
    cbticks = None
    if scheme == 'source' or scheme == 'network':
        cmap = None
        cmin = None
        cmax = None
        legend = dict()
        if scheme == 'source':
            color_scheme = color_scheme_validation_source
            apply_to = 'val.src'
        elif scheme == 'network':
            color_scheme = color_scheme_validation_network
            apply_to = 'val.net'

        colors = ma.array(['black']*m.get_nb_matchups(),mask=m.mask,dtype='S10').astype(str)
        for k,c in color_scheme.items():
            ls = msel.nocase_label_selector(k,apply_to)
            colorize_with_c = ls.select(m) * ~m.mask
            colors[colorize_with_c] = c
            if sum(colorize_with_c) > 0:
                legend[k] = c

    elif scheme == 'time':
        t0s = np.array(m.val.t0)
        if m.val.t0.mask is ma.nomask:
            t_index = np.array([True]*m.get_nb_matchups())
        else:
            t_index = ~m.val.t0.mask

        if isinstance(start_period, list):
            sp0 = start_period[0]
        else:
            sp0 = start_period
        if isinstance(end_period, list):
            ep0 = end_period[-1]
        else:
            ep0 = end_period
        #min_t = t0s[t_index].min()
        #max_t = t0s[t_index].max()
        min_t = datetime.combine(sp0, time(0))
        max_t = datetime.combine(ep0, time(0)) + timedelta(days=1)
        tdelta_tot = float(timedelta_total_seconds(max_t - min_t))
        colors = sp.array([timedelta_total_seconds(t - min_t)/tdelta_tot for t in t0s])
        colors = ma.array(colors, mask = m.val.t0.mask)
        cmin   = 0
        cmax   = 1
        cmap = cm.get_cmap('plasma')
    elif scheme == 'mismatch_len':
        cmap   = cm.RdBu_r
        colors = m.prd.len - m.val.len
        cmin   = -15.0
        cmax   = -cmin
        cbticks = [-15.,-10.,-5.]
        cbticks = cbticks + [0] + [abs(t) for t in cbticks[::-1]]
        cblab   = 'prd_len - val_len [km]'
    elif scheme == 'mismatch_dt':
        cmap   = cm.RdBu_r
        colors = ma.array([dt.total_seconds()/60./60. for dt in m.coll.diff_dt.data],mask=m.coll.diff_dt.mask)
        cmin   = -5
        cmax   = -cmin
        cbticks = [-5.,-2.5,]
        cbticks = cbticks + [0] + [abs(t) for t in cbticks[::-1]]
        cblab   = 'prd_dT - val_dT [h]'
    elif scheme == 'uncert':
        cmap   = cm.jet
        colors = (m.prd.sX + m.prd.sY)*0.5
        cmin   = 0.
        cmax   = colors.max()
        cblab = 'Product Uncertainty [km]'
    else:
        raise ValueError("Coloring scheme %s is not supported" % (scheme,))


    symbols = 'o'
    siz = 8
    return colors,siz, symbols, cmap, cmin, cmax, legend, cbticks, cblab

def compute_correls_statistics(m1,m2):
    """
       compute main matchup statistics (bias, RMSE, correlation, etc...) for
       a pair of matchup objects.
    """
    from scipy import stats
    mstats = dict()

    # Use the compute_matchup_statistics() on each matchup, and transfer to the mstats
    mstats1 = compute_matchup_statistics(m1)
    mstats2 = compute_matchup_statistics(m2)
    for k in mstats1.keys():
        mstats[k+'_1'] = mstats1[k]
    for k in mstats2.keys():
        mstats[k+'_2'] = mstats2[k]

    # now compute the cross-matchups statistics
    if mstats['N_1'] != mstats['N_2']:
        raise ValueError("Something went wrong with masking/selecting in the matchups! m1 has {} while m2 has {}".format(mstats['N_1'],mstats['N_2']))

    if mstats['N_1'] == 0:
        return mstats

    # bias and RMSE
    mismatch_x_1 = (m1.prd.dX - m1.val.dX)
    mismatch_y_1 = (m1.prd.dY - m1.val.dY)
    mismatch_x_2 = (m2.prd.dX - m2.val.dX)
    mismatch_y_2 = (m2.prd.dY - m2.val.dY)
    covarM_x_12 = ma.cov(mismatch_x_1,mismatch_x_2)
    covarM_y_12 = ma.cov(mismatch_y_1,mismatch_y_2)

    mstats['corr_dx_12'] = covarM_x_12[0,1] / (sqrt(covarM_x_12[0,0]) * sqrt(covarM_x_12[1,1]))
    mstats['corr_dy_12'] = covarM_y_12[0,1] / (sqrt(covarM_y_12[0,0]) * sqrt(covarM_y_12[1,1]))

    return mstats

def compute_matchup_statistics(m):
    """
       compute main matchup statistics (bias, RMSE, correlation, etc...) for
       a matchup object, taking into account the selection masks, if any.
    """
    from scipy import stats
    mstats = dict()

    # number of matchups
    mstats['N'] = m.prd.dX.count()
    if mstats['N'] != m.get_nb_selected_matchups():
        raise ValueError("Something went wrong with masking/selecting in the matchups!")

    if m.get_nb_selected_matchups() == 0:
        return mstats

    # bias and RMSE
    mismatch_x = m.prd.dX - m.val.dX
    mismatch_y = m.prd.dY - m.val.dY

    covarM = ma.cov(mismatch_x,mismatch_y)

    mstats['bias_dx'] = mismatch_x.mean()
    mstats['bias_dy'] = mismatch_y.mean()
    mstats['sdev_dx'] = sqrt(covarM[0,0])
    mstats['sdev_dy'] = sqrt(covarM[1,1])
    mstats['corr_dxdy'] = covarM[0,1] / (mstats['sdev_dx'] * mstats['sdev_dy'])

    # parameters for draing the error ellipse
    eigvals, eigvect = np.linalg.eig(covarM)
    mstats['ell_length'] = sqrt(eigvals[0])
    mstats['ell_height'] = sqrt(eigvals[1])
    mstats['ell_alpha'] = 0

    # linear regression for concatenated dX and dY (a,b,corr)
    xy_prd = ma.concatenate((m.prd.dX,m.prd.dY)).compressed()
    xy_val = ma.concatenate((m.val.dX,m.val.dY)).compressed()
    val_n = len(xy_val)
    prd_n = len(xy_prd)
    if val_n != prd_n:
        raise ValueError("Something is wrong with masking and ma.compressed()")
    slope, intercept, r_value, p_value, std_err = stats.linregress(xy_val,xy_prd)

    mstats['a'] = slope
    mstats['b'] = intercept
    mstats['corr'] = r_value

    # linear regression for dX (a,b,corr)
    x_prd = m.prd.dX.compressed()
    x_val = m.val.dX.compressed()
    val_n = len(x_val)
    prd_n = len(x_prd)
    if val_n != prd_n:
        raise ValueError("Something is wrong with masking and ma.compressed()")
    res = stats.linregress(x_val,x_prd)
    # Two-sided inverse Students t-distribution
    # p - probability, df - degrees of freedom
    tinv = lambda p, df: abs(statst.ppf(p/2, df))
    ts = tinv(0.05, len(x_val)-2)
    slope_err = ts * res.stderr
    intercept_err = ts * res.intercept_stderr
    mstats['a_dx'] = res.slope
    mstats['b_dx'] = res.intercept
    mstats['corr_dx'] = res.rvalue
    mstats['a_dx_e95'] = slope_err
    mstats['b_dx_e95'] = intercept_err

    # linear regression for dY (a,b,corr)
    y_prd = m.prd.dY.compressed()
    y_val = m.val.dY.compressed()
    val_n = len(y_val)
    prd_n = len(y_prd)
    if val_n != prd_n:
        raise ValueError("Something is wrong with masking and ma.compressed()")
    res = stats.linregress(y_val,y_prd)
    # Two-sided inverse Students t-distribution
    # p - probability, df - degrees of freedom
    tinv = lambda p, df: abs(statst.ppf(p/2, df))
    ts = tinv(0.05, len(y_val)-2)
    slope_err = ts * res.stderr
    intercept_err = ts * res.intercept_stderr
    mstats['a_dy'] = res.slope
    mstats['b_dy'] = res.intercept
    mstats['corr_dy'] = res.rvalue
    mstats['a_dy_e95'] = slope_err
    mstats['b_dy_e95'] = intercept_err

    return mstats

def dump_matchup_to_ascii(m,outfile,outliers=False):

    mismatch_lim = 0.
    if outliers:
        mismatch_lim = 15.

    with open(outfile,"w") as fp:

        header = " ".join(('val_id','val_netw','val_t0','val_lon0','val_lat0','val_dX', 'val_dY', 'val_len',\
                           'prd_t0','prd_lon0','prd_lat0','prd_dX', 'prd_dY','diff_len','diff_t0','prd_flg',))
        fp.write("# "+header+'\n')

        for lm in range(len(m)):
            if (m.mask is not ma.nomask) and (m.mask[lm]):
                continue
            misX = m.prd.dX[lm] - m.val.dX[lm]
            misY = m.prd.dY[lm] - m.val.dY[lm]
            if (abs(misX) >= mismatch_lim or abs(misY) >= mismatch_lim):
                line = " ".join((
                    m.val.id[lm],
                    m.val.net[lm],
                    m.val.t0[lm].strftime('%Y-%m-%dT%H:%M:%SZ'),
                    str(m.val.lon0[lm]),str(m.val.lat0[lm]),
                    str(m.val.dX[lm]),str(m.val.dY[lm]),
                    str(m.val.len[lm]),
                    m.prd.t0[lm].strftime('%Y-%m-%dT%H:%M:%SZ'),
                    str(m.prd.lon0[lm]),str(m.prd.lat0[lm]),
                    str(m.prd.dX[lm]),str(m.prd.dY[lm]),
                    str(m.prd.len[lm] - m.val.len[lm]),
                    str(m.coll.diff_t0[lm].total_seconds()/60./60.),
                   "%d" % (m.prd.flg[lm]),))
                fp.write(line+'\n')


def plot_geo_coverage(m,date1,date2,outfile,product='na',channels=None,zone='all',color_scheme='mismatch_len',buoylabs=None,):

    titlefont = 24
    tickfont = 20

    vals_lon_1 = m.val.lon0
    vals_lat_1 = m.val.lat0

    # prepare the map (cartopy)
    if m.area == 'nh':
        ax = plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_latitude=+90.0,))
        plotcrs = ccrs.LambertAzimuthalEqualArea(central_latitude=+90.0,)
        xlims = (-2500000,2000000)
        ylims = (-1800000,2700000)
    elif m.area == 'sh':
        ax = plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_latitude=-90.0,))
#        xlims = (-3000000,+0)
#        ylims = (+500000,+3000000)
        xlims = (-3000000,3500000)
        ylims = (-3000000,3500000)
        plotcrs = ccrs.LambertAzimuthalEqualArea(central_latitude=-90.0,)
        ice_shelves = cfeature.NaturalEarthFeature(
            category='physical',
            name='antarctic_ice_shelves_polys',
            scale='10m',
            facecolor='grey',
            edgecolor='grey')
        ax.add_feature(ice_shelves, )
    else:
        raise ValueError("Only know about nh and sh matchups")

    plt.gcf().set_size_inches((10,10))
    ax.add_feature(cartopy.feature.LAND, color='gray')

    ax.set_global()
    ax.gridlines()

    cl, siz, mk, cmap, cmin, cmax, legend, cbticks, cblab = get_color(
        m, color_scheme, start_period=date1[0], end_period=date2[-1])

    cl = cl.compressed()

    points = ax.scatter(vals_lon_1.compressed(), vals_lat_1.compressed(),
                        s=siz, c=cl, marker=mk, cmap=cmap, vmin=cmin,
                        vmax=cmax, transform=ccrs.PlateCarree())

    if buoylabs:
        for i in range(len(buoylabs[0])):
            bpx, bpy =  plotcrs.transform_point(buoylabs[0][i],
                                                buoylabs[1][i],
                                                src_crs=ccrs.PlateCarree())
            ax.text(bpx, bpy, buoylabs[2][i], fontsize='xx-small')

    figratio = ax.get_ylim()[1] / ax.get_xlim()[1]
    cbar_orient = 'horizontal'
    tit_func = plt.title
    if figratio > 1.:
        cbar_orient = 'vertical'
        tit_func = ax.set_ylabel

    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)

    points.set_linewidth(0.1)

    tit_func("Validation drifters for %s\n%s (%s -> %s)" % (product,m.area.upper(),titled(date1),titled(date2, pos=-1)), fontsize=titlefont)

    cb = None
    if color_scheme.startswith('mismatch_'):
        cb = plt.colorbar(points,orientation=cbar_orient,shrink=0.5,pad=0.05,ticks=cbticks)
    elif color_scheme == 'time':
        cb = plt.colorbar(points,shrink=0.5,pad=0.05,orientation=cbar_orient,ticks=[])
        cblab = '{} to {}'.format(titled(date1),titled(date2,pos=-1))
    elif color_scheme == 'network':
        markers = []
        labels = []
        for src in legend.keys():
            markers.append(mlines.Line2D(range(1), range(1), color="white", marker=mk, markerfacecolor=legend[src]))
            labels.append(src)
            ax.legend(markers,labels,scatterpoints=1,scatteryoffsets=[0.5],markerscale=1,prop={'size':fs},\
                      loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=len(markers), borderaxespad=0.)

    if cb is not None:
        if (figratio > 1):
            cb.ax.set_ylabel(cblab, fontsize=tickfont)
            cb.ax.yaxis.set_label_position('left')
        else:
            cb.ax.set_xlabel(cblab, fontsize=tickfont)
            cb.ax.xaxis.set_label_position('top')

    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()
    plt.close()

    return outfile

def plot_magnitude_histograms(m,date1,date2,outfile,as_velocity=True):
    prds = np.sort(m.prd.len)
    vals = np.sort(m.val.len)

    unit = 'km'
    if as_velocity:
        convf = 1000. / (2. * 24 * 60. * 60.)
        prds *=  convf
        vals *=  convf
        unit = 'm/s'

    vmin = 0
    vmax = vals.max()
    vbin = 30.

    print("Max of validation: {} {}".format(vmax,unit))

    y0 = 0
    y1 = 5
    marg = (y1-y0)*0.25


    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    Vn, Vbins = np.histogram(vals,bins=vbin,range=(vmin,vmax))
    Pn, Pbins = np.histogram(prds,bins=vbin,range=(vmin,vmax))
    Vc = (Vbins[:-1] + Vbins[1:]) / 2
    Pc = (Pbins[:-1] + Pbins[1:]) / 2
    ax.bar(Vc,Vn,zs=y0,zdir='y',align='center',color='g',width=0.8*(Vbins[1]-Vbins[0]))
    ax.bar(Pc,Pn,zs=y1,zdir='y',align='center',color='r',width=0.8*(Pbins[1]-Pbins[0]))
    ax.set_xlim(vmin,vmax)
    ax.set_ylim(y0-marg,y1+marg)
    xlab = 'Magnitude of 48h drift [{}]'.format(unit,)
    if as_velocity:
        xlab = 'Mean velocity in 48h drift [{}]'.format(unit,)
    ax.set_xlabel(xlab)
    ax.set_zlabel('Occurrence')
    #ax.set_ylabel('Y - Class')
    ax.set_yticks((y0,y1))
    ax.set_yticks(())
    #ax.set_yticklabels(('Validation','Product'))
#    ax.text((vmin+vmax)*0.5,y0-0.5*marg,0.,'Validation',zdir='x',ha='center',color='g')
    ax.text((vmin+vmax)*0.5,y0-0.5*marg,0.,'Buoy',zdir='x',ha='center',color='g')
    ax.text((vmin+vmax)*0.5,y1-0.5*marg,0.,'Product',zdir='x',ha='center',color='r')
    plt.title("Histogram of 48h drift magnitude\n%s (%s -> %s)" % (m.area.upper(),titled(date1),titled(date2,pos=-1)))
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

    return outfile


def plot_qq_plot(m,date1,date2,outfile,product='na',channels=None,status_flag=-1):
    prds = np.sort(m.prd.len)
    vals = np.sort(m.val.len)

    vmin = 0
    vmax = 100

    mstats = compute_matchup_statistics(m)

    plt.plot(vals,prds,'r+')
    plt.plot([vmin,vmax],[vmin,vmax],'k-')
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim(vmin,vmax)
    ax.set_ylim(vmin,vmax)
    #ax.set_xlabel('Validation data [km]')
    ax.set_xlabel('Buoy vector [km]')
    ax.set_ylabel('Product vector [km]')
    ax.set_title("QQ-plot of 48h drift magnitude\n%s (%s -> %s)" % (m.area.upper(),titled(date1),titled(date2,pos=-1)))
    ty = tstart
    product_str = product.upper()
    if channels:
        product_str += ' (' + channels.upper() + ')'
    ax.text(tx,ty,product_str,transform=ax.transAxes,fontsize=fs); ty -= tstep
    ax.text(tx,ty,"N = %d" % (mstats['N']),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    if (status_flag >= 0):
        txt = 'Only Flag {}'.format(flag_dict[str(status_flag)])
    else:
        txt = ''
        #txt = 'All Flags'
    ax.text(tx,ty,txt,
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

    return outfile

def plot_monthly_timeseries(m,date1,date2,outfile,product='na',status_flag=-1):

    area  = m.area
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator(interval=2)  # every 2 months
    DateFmt = mdates.DateFormatter('%Y-%m')

    # collect the monthly validation statistics
    T = []
    monthly_stats = []
    insitu_net = []
    locdate1 = [date1[0],]
    locdate2 = [date2[-1],]
    sub_periods = mfil.split_period_in_months(locdate1,locdate2)
    for sp in sub_periods:
        sav_mask = m.mask
        try:
            s0 = sp[0]
            s1 = sp[1]+timedelta(days=1)-timedelta(seconds=0.1)
            psel = msel.period_selector(s0,s1,'prd.t0')
            m.apply_selectors(psel)
            c = Counter(m.val.net[~m.mask])
            mstats = compute_matchup_statistics(m)
            monthly_stats.append(mstats)
            insitu_net.append(c)
            T.append(s0)
        except Exception as ex:
            print("Warning: caught ex while computing statistics for {}->{}: {}".format(s0.date(),s1.date(),ex,))

        m.mask = sav_mask
        m._applymask()

    uniq_insitu_netw = set(list(chain(*[itw.keys() for itw in insitu_net])))

    # transform these in arrays
    msk = np.zeros(len(monthly_stats)).astype('bool')
    Time= ma.empty(len(monthly_stats),dtype=type(T[0]))
    TimeNum = ma.empty(len(monthly_stats),dtype='float')
    N   = ma.empty(len(monthly_stats),dtype='int')
    Sx  = ma.empty(len(monthly_stats))
    Sy  = ma.empty_like(Sx)
    Bx  = ma.empty_like(Sx)
    By  = ma.empty_like(Sx)
    for mon in range(len(N)):
        Time[mon] = T[mon] + timedelta(days=15)
        TimeNum[mon] = mdates.date2num(Time[mon])
        N[mon]  = monthly_stats[mon]['N']
        if N[mon] == 0:
            msk[mon] = True
        else:
            Sx[mon] = monthly_stats[mon]['sdev_dx']
            Sy[mon] = monthly_stats[mon]['sdev_dy']
            Bx[mon] = monthly_stats[mon]['bias_dx']
            By[mon] = monthly_stats[mon]['bias_dy']

    #apply mask
    Time.mask = msk
    TimeNum.mask = msk
    N.mask    = msk
    Sx.mask   = msk
    Sy.mask   = msk
    Bx.mask   = msk
    By.mask   = msk

    #additional mask on the statistics if N too low
    min_N = 25
    Sx.mask[N<min_N]   = True
    Sy.mask[N<min_N]   = True
    Bx.mask[N<min_N]   = True
    By.mask[N<min_N]   = True

    max_S = 10.
    min_B = -2.

    # needed to avoid strange behaviour of pylab when some
    #   points are outside the plot area
    Sx.data[Sx>max_S] = max_S
    Sy.data[Sy>max_S] = max_S
    Bx.data[Bx<min_B] = min_B
    By.data[By<min_B] = min_B

    min_T = mdates.num2date(TimeNum.min(),).date()
    max_T = mdates.num2date(TimeNum.max(),).date()
    # plot
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax1 = plt.subplot(gs[0])
    xx = ax1.plot_date(TimeNum,Sx,ls='-',color=cX,label='sdev dX',lw=2)
    yy = ax1.plot_date(TimeNum,Sy,ls='-',color=cY,label='sdev dY',lw=2)
    ax1.plot_date(TimeNum,Bx,ls='-.',color=cX,label='bias dX')
    ax1.plot_date(TimeNum,By,ls='-.',color=cY,label='bias dY')
    ax1.plot([date1[0],date2[-1]],[0.,0.],ls=':',color='k',label=None)


    ty = 0.90
    ty = 0.90
    ax1.text(tx,ty,product.upper(),transform=ax1.transAxes,fontsize=fs); ty -= 2*tstep
    if (status_flag >= 0):
        txt = 'Flag {}:{}'.format(status_flag,flag_dict[status_flag][1])
    else:
        txt = ''
        #txt = 'All Flags'
    ax1.text(tx,ty,txt,
             fontsize=fs,transform=ax1.transAxes); ty -= tstep

    ax1.set_ylabel('Statistics of mismatch [km]')
    ax1.set_title("Monthly validation statistics\n{} ({} -> {})".format(m.area.upper(),titled(date1),titled(date2,pos=-1)))
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(DateFmt)
    ax1.legend(loc='upper center',bbox_to_anchor=[0.5,0.],ncol=4,prop={'size':fs})

    datelim1 = date1[0]
    datelim2 = date2[-1]
    if isinstance(datelim1, list):
        datelim1 = datelim1[0]
        datelim2 = datelim2[-1]
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_xlim(datelim1, datelim2+timedelta(days=1))
    ax1.set_ylim(min_B,max_S)

    ax2 = plt.subplot(gs[1],sharex=ax1)
    ax2.plot_date(TimeNum,N,ls='-',color='k',lw=2,label='all')
    for netw in uniq_insitu_netw:
        N_netw = [inw[netw] for inw in insitu_net]
        ax2.plot(TimeNum,N_netw,ls='-',
                 color=color_scheme_validation_network[netw.upper()],
                 lw=2,label=netw)

    ax2.set_ylabel('N')
    maxN = 1000
    if m.area.upper() == 'SH':
        maxN = 200
    ax2.legend(loc='upper center',bbox_to_anchor=[0.5,1.],ncol=len(uniq_insitu_netw)+1,prop={'size':int(ceil(fs/1.5))})
    ax2.xaxis.set_major_formatter(DateFmt)
    ax2.set_xlim(datelim1, datelim2+timedelta(days=1))
    ax2.set_ylim(0,maxN)
    fig.autofmt_xdate()

    # add vertical shading for summer
#    if area == 'nh':
#        summer = [date(1970,5,1),date(1970,6,1),date(1970,9,1),date(1970,10,1)]
#    else:
#        summer = [date(1970,11,1),date(1970,12,1),date(1970,3,1)-timedelta(days=1),date(1970,3,31)]

#    import calendar
#    lmin_T = date(min_T.year,min_T.month,1)
#    lmax_T = date(max_T.year,max_T.month,calendar.monthrange(max_T.year,max_T.month)[1])
#    lmin_T = date1[0]
#    lmax_T = date2[-1]
#    for ax in [ax1, ax2,]:
#        if area == 'nh':
#            for yyyy in range(min_T.year,max_T.year+1):
#                ax.axvspan(max([lmin_T,date(yyyy,summer[0].month,summer[0].day)]),
#                            min([date(yyyy,summer[3].month,summer[3].day),lmax_T]),
#                            facecolor='0.85',edgecolor='none',lw=0)
#                ax.axvspan(max([lmin_T,date(yyyy,summer[1].month,summer[1].day)]),
#                            min([date(yyyy,summer[2].month,summer[2].day),lmax_T]),
#                            facecolor='0.75',edgecolor='none',lw=0)
#        else:
#            for yyyy in range(min_T.year-1,max_T.year+1):
#                ax.axvspan(max([lmin_T,date(yyyy,summer[0].month,summer[0].day)]),
#                           min([date(yyyy+1,summer[3].month,summer[3].day),lmax_T]),
#                           facecolor='0.85',edgecolor='none',lw=0)
#                ax.axvspan(max([lmin_T,date(yyyy,summer[1].month,summer[1].day)]),
#                            min([date(yyyy+1,summer[2].month,summer[2].day),lmax_T]),
#                            facecolor='0.75',edgecolor='none',lw=0)

#    # final force of the xlim (in case the summer shading changed things)
#    ax1.set_xlim(datelim1, datelim2+timedelta(days=1))
#    ax2.set_xlim(datelim1, datelim2+timedelta(days=1))

    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

    return outfile

def plot_correl_scatterplot(m1,m2,date1,date2,outfile,product1,product2,status_flag=-1):

    # now get the matchups that are common to both m1 and m2, which means they have the same uuid
    m1,m2 = mio.get_common_matchups(m1,m2)
    print("Now {} {} matchups".format(m1.count(),m2.count()))

    prds_dX_1 = m1.prd.dX.compressed()
    vals_dX_1 = m1.val.dX.compressed()
    prds_dY_1 = m1.prd.dY.compressed()
    vals_dY_1 = m1.val.dY.compressed()

    prds_dX_2 = m2.prd.dX.compressed()
    vals_dX_2 = m2.val.dX.compressed()
    prds_dY_2 = m2.prd.dY.compressed()
    vals_dY_2 = m2.val.dY.compressed()

    mstats = compute_correls_statistics(m1,m2)

    eX_1 = prds_dX_1 - vals_dX_1
    eY_1 = prds_dY_1 - vals_dY_1
    eX_2 = prds_dX_2 - vals_dX_2
    eY_2 = prds_dY_2 - vals_dY_2

    vmin = -20.
    vmax = +20.

    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    nx  = ('X','Y')

    for ix,ax in enumerate((ax1, ax2)):
        ax.set_aspect('equal')
        ax.plot([vmin,vmax],[0.,0.],'k-.')
        ax.plot([0.,0.],[vmin,vmax],'k-.')
        var_x = eval('e{}_1'.format(nx[ix]))
        var_y = eval('e{}_2'.format(nx[ix]))

        xx = ax.scatter(var_x,var_y,marker='+',label=None,s=20,lw=0.5)
        ax.set_xlim(vmin,vmax)
        ax.set_ylim(vmin,vmax)
        ax.set_xlabel('Diff in d{} for {} [km]'.format(nx[ix],product1))
        ax.set_ylabel('Diff in d{} for {} [km]'.format(nx[ix],product2))
        ty = tstart
        ax.text(tx,ty,"N = %d" % (mstats['N_1']),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
        ax.text(tx,ty,"sdev_{} {} = {:.2f}".format(nx[ix],product1,mstats['sdev_d{}_1'.format(nx[ix].lower())]),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
        ax.text(tx,ty,"sdev_{} {} = {:.2f}".format(nx[ix],product2,mstats['sdev_d{}_2'.format(nx[ix].lower())]),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
        ax.text(tx,ty,"corr_{} ({},{}) = {:.2f}".format(nx[ix],product1,product2,mstats['corr_d{}_12'.format(nx[ix].lower())]),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
        if (status_flag >= 0):
            txt = 'Flag {}:{}'.format(status_flag,flag_dict[status_flag][1])
        else:
            txt = ''
            #txt = 'All Flags'
        ax.text(tx,ty,txt,
             fontsize=fs,transform=ax.transAxes); ty -= tstep

    fig.suptitle("Correlation of mismatch {} vs {}\n{} ({} -> {})".format(
        product1, product2, m1.area.upper(), titled(date1), titled(date2,pos=-1)))
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

    return outfile

def plot_eXeY_scatterplot(m,date1,date2,outfile,product='na',channels=None,status_flag=-1, datestr=None):

    prds_dX = m.prd.dX
    vals_dX = m.val.dX
    prds_dY = m.prd.dY
    vals_dY = m.val.dY

    eX = prds_dX - vals_dX
    eY = prds_dY - vals_dY

    vmin = -20
    vmax = 20

    try:
        mstats = compute_matchup_statistics(m)
    except:
        print("No matchup stats could be calculated")
        sys.exit(0)
    if mstats['N'] == 0:
        print("No points available")
        sys.exit(0)

    fig = plt.figure(figsize=(4.15,4.15))
    ax = fig.add_subplot(1,1,1)
# This method doesn't work
#    ax = plt.gca()
    ax.set_aspect('equal')
    ax.plot([vmin,vmax],[0.,0.],'k-.')
    ax.plot([0.,0.],[vmin,vmax],'k-.')

    xx = ax.scatter(eX,eY,marker='+',label=None,s=20,lw=0.5)

#    sig_muls = (1.5,2.5,)
#    sig_cols = ('red','green')
#    for m in range(len(sig_muls)):
#        e = Ellipse(xy=(mstats['bias_dx'], mstats['bias_dy']),
#                    width=sig_muls[m]*mstats['sdev_dx'],
#                    height=sig_muls[m]*mstats['sdev_dy'],)
#        ax.add_artist(e)
#        e.set_clip_box(ax.bbox)
#        e.set_fill(False)
#        e.set_edgecolor(sig_cols[m])

##    ax.scatter([mstats['bias_dx']],[mstats['bias_dy']],marker='+',label=None,s=20,lw=3,c='black')
    ty = tstart
    product_str = product.upper()
    if channels:
        product_str += ' (' + channels.upper() + ')'
    ax.text(tx,ty,product_str,transform=ax.transAxes,fontsize=fs); ty -= tstep
    ax.text(tx,ty,"bias [km] = (%+05.3f,%+05.3f)" % (mstats['bias_dx'],mstats['bias_dy']),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    ax.text(tx,ty,"sdev [km] = (%05.3f,%05.3f)" % (mstats['sdev_dx'],mstats['sdev_dy']),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    ax.text(tx,ty,"corr(eX,eY) = %05.3f" % (mstats['corr_dxdy']),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    ax.text(tx,ty,"N = %d" % (mstats['N']),
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    if (status_flag >= 0):
        txt = 'Only Flag {}'.format(str(status_flag))
    else:
        txt = ''
        #txt = 'All Flags'
    ax.text(tx,ty,txt,
             fontsize=fs,transform=ax.transAxes); ty -= tstep
    ax.set_xlim(vmin,vmax)
    ax.set_ylim(vmin,vmax)
    ax.set_xlabel('Deviation in dX (Prd - Val) [km]')
    ax.set_ylabel('Deviation in dY (Prd - Val) [km]')
    if datestr is None:
        ax.set_title("Bi-dimensional PDF of errors in dX and dY\n%s (%s -> %s)" % (m.area.upper(),titled(date1),titled(date2,pos=-1)))
    else:
        ax.set_title("Bi-dimensional PDF of errors in dX and dY\n{} ({})".format(m.area.upper(), datestr))
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()
    plt.close()

    return outfile

def plot_dXdY_scatterplot(m,date1,date2,outfile,product='na',channels=None,status_flag=-1, hexplot=False, datestr=None):

    fsdxdy = 20
    titlefont = 24
    tickfont = 20

#    prds_dX = m.prd.dX
#    vals_dX = m.val.dX
#    prds_dY = m.prd.dY
#    vals_dY = m.val.dY

    prds_dX = m.prd.dX.compressed()
    vals_dX = m.val.dX.compressed()
    prds_dY = m.prd.dY.compressed()
    vals_dY = m.val.dY.compressed()

    vmin = -100.
    vmax = +100.

    vmin = -50.
    vmax = +50.

    prd_relangle = m.prd.relangle

    try:
        mstats = compute_matchup_statistics(m)
    except:
        print("No matchup stats could be calculated")
        sys.exit(1)
    if mstats['N'] == 0:
        print("No points available")
        sys.exit(1)
    #print("mstats['N'] = ", mstats['N'])
    ppercent = mstats['N'] * 0.0001

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    ax.plot([vmin,vmax],[vmin,vmax],'k-.')
    ax.set_aspect('equal')

#    # Robustly fit linear model with RANSAC algorithm
#    horiz = np.arange(-50., 50.)[:, np.newaxis]
#    ransac_x = linear_model.RANSACRegressor()
#    ransac_x.fit(np.array(vals_dX)[:, np.newaxis], np.array(prds_dX))
#    inlier_mask_x = ransac_x.inlier_mask_
#    outlier_mask_x = np.logical_not(inlier_mask_x)
#    x_vert = ransac_x.predict(horiz)
#    ransac_y = linear_model.RANSACRegressor()
#    ransac_y.fit(np.array(vals_dY)[:, np.newaxis], np.array(prds_dY))
#    inlier_mask_y = ransac_y.inlier_mask_
#    outlier_mask_y = np.logical_not(inlier_mask_y)
#    y_vert = ransac_y.predict(horiz)
#    ax.plot(horiz, x_vert, c=cX, ls='-', lw=3)
#    ax.plot(horiz, y_vert, c=cY, ls='-', lw=3)

#    print("Estimated coefficients RANSAC dX: {}".format(ransac_x.estimator_.coef_))
#    print("Estimated coefficients RANSAC dY: {}".format(ransac_y.estimator_.coef_))

    # Finding the confidence limits
#    nrange = 100
#    vrange = np.array([vmin + (x * (vmax - vmin) / nrange) for x in range(nrange)])
#    ulimarr = np.zeros((nrange, 2))
#    ulimarr[:, 0] = (vrange * (mstats['a_dx'] + mstats['a_dx_e95'])) + mstats['b_dx'] + mstats['b_dx_e95']
#    ulimarr[:, 1] = (vrange * (mstats['a_dx'] - mstats['a_dx_e95'])) + mstats['b_dx'] + mstats['b_dx_e95']
#    ulimx = np.max(ulimarr, axis=1)
#    llimarr = np.zeros((nrange, 2))
#    llimarr[:, 0] = (vrange * (mstats['a_dx'] + mstats['a_dx_e95'])) + mstats['b_dx'] - mstats['b_dx_e95']
#    llimarr[:, 1] = (vrange * (mstats['a_dx'] - mstats['a_dx_e95'])) + mstats['b_dx'] - mstats['b_dx_e95']
#    llimx = np.min(llimarr, axis=1)

#    ulimarr = np.zeros((nrange, 2))
#    ulimarr[:, 0] = (vrange * (mstats['a_dy'] + mstats['a_dy_e95'])) + mstats['b_dy'] + mstats['b_dy_e95']
#    ulimarr[:, 1] = (vrange * (mstats['a_dy'] - mstats['a_dy_e95'])) + mstats['b_dy'] + mstats['b_dy_e95']
#    ulimy = np.max(ulimarr, axis=1)
#    llimarr = np.zeros((nrange, 2))
#    llimarr[:, 0] = (vrange * (mstats['a_dy'] + mstats['a_dy_e95'])) + mstats['b_dy'] - mstats['b_dy_e95']
#    llimarr[:, 1] = (vrange * (mstats['a_dy'] - mstats['a_dy_e95'])) + mstats['b_dy'] - mstats['b_dy_e95']
#    llimy = np.min(llimarr, axis=1)

#    ax.plot([vmin,vmax],np.array([vmin,vmax])*mstats['a_dx'] + mstats['b_dx'],c=cX,ls='-',lw=2)
#    ax.plot([vmin,vmax],np.array([vmin,vmax])*mstats['a_dy'] + mstats['b_dy'],c=cY,ls='-',lw=2)

#    print("vrange = ", vrange)
#    print("llimx = ", llimx)
#    print("ulimx = ", ulimx)
#    print("llimy = ", llimy)
#    print("ulimy = ", ulimy)
#    ax.fill_between(vrange, llimx, ulimx, color=cX, alpha=0.5, edgecolor="")
#    ax.fill_between(vrange, llimy, ulimy, color=cY, alpha=0.5, edgecolor="")
#    ax.fill_between(vrange, llimx, y2=ulimx, color=cX)#, alpha=0.5, edgecolor="")
#    ax.fill_between(vrange, llimy, y2=ulimy, color=cY)#, alpha=0.5, edgecolor="")

    if not hexplot:
        xx = ax.scatter(vals_dX,prds_dX,marker='+',c=cX,label=None,s=20,lw=0.5)
        yy = ax.scatter(vals_dY,prds_dY,marker='+',c=cY,label=None,s=20,lw=0.5)
        ax.legend([xx,yy],['dX','dY'],scatterpoints=1,loc=4,scatteryoffsets=[0.5],markerscale=2,prop={'size':fs})
    else:
#        vals_dXYs = np.ma.concatenate([vals_dX, vals_dY])
#        prds_dXYs = np.ma.concatenate([prds_dX, prds_dY])
        vals_dXYs = np.concatenate([vals_dX, vals_dY])
        prds_dXYs = np.concatenate([prds_dX, prds_dY])
#        hbcmap = cmocean.cm.dense
        hbcmap = cmocean.cm.haline
        hb = ax.hexbin(vals_dXYs, prds_dXYs, gridsize=50, bins='log',
                       mincnt=ppercent, cmap=hbcmap)
        cbaxes = ax.inset_axes([0.2, 0.05, 0.6, 0.03])
        cbar = plt.colorbar(hb, cax=cbaxes, orientation='horizontal')
        cbar.ax.tick_params(labelsize=tickfont)

    #xx = ax.scatter(vals_dX,prds_dX,marker='o',c=m.prd.relangle,label=None,s=30,lw=0.5,cmap=cm.inferno,vmin=0,vmax=180)
    #yy = ax.scatter(vals_dY,prds_dY,marker='^',c=m.prd.relangle,label=None,s=30,lw=0.5,cmap=cm.inferno,vmin=0,vmax=180)

    #xx = ax.scatter(vals_dX,prds_dX,marker='o',c=m.coll.dist,label=None,s=30,lw=0.5,cmap=cm.inferno,vmin=0,vmax=30)
    #yy = ax.scatter(vals_dY,prds_dY,marker='^',c=m.coll.dist,label=None,s=30,lw=0.5,cmap=cm.inferno,vmin=0,vmax=30)

    #xx = ax.scatter(vals_dX,prds_dX,marker='o',c=m.prd.flg,label=None,s=30,lw=0.5,cmap=cm.tab20c,vmin=0,vmax=20)
    #yy = ax.scatter(vals_dY,prds_dY,marker='^',c=m.prd.flg,label=None,s=30,lw=0.5,cmap=cm.tab20c,vmin=0,vmax=20)

    ty = tstart
    product_str = product.upper()
    if product_str.startswith('WIND'):
        product_str = 'FREE-DRIFT MODEL'
    if product_str.startswith('MRG'):
        product_str = 'MULTI-OI'
    if product_str.startswith('AMSRE'):
        product_str = 'AMSR-E'
    if product_str.endswith('_NH') or product_str.endswith('_SH'):
        product_str = product_str[:-3]
    if product_str.endswith('_TB37') or product_str.endswith('_TB19'):
        product_str = product_str[:-5]
    if channels and channels != [None] and channels != 'none' and channels != 'None':
        product_str += ' (' + channels.upper() + ')'
    ax.text(tx,ty,product_str,transform=ax.transAxes,fontsize=fsdxdy); ty -= tstep
    if mstats['N'] > 0:
        ax.text(tx,ty,"bias [km] = (%+05.3f, %+05.3f)" % (mstats['bias_dx'],mstats['bias_dy']),
             fontsize=fsdxdy,transform=ax.transAxes); ty -= tstep
        ax.text(tx,ty,"sdev [km] = (%05.3f, %05.3f)" % (mstats['sdev_dx'],mstats['sdev_dy']),
             fontsize=fsdxdy,transform=ax.transAxes); ty -= tstep
#    ax.text(tx,ty,"linfit dX = (a:%04.2f, b:%+05.3f km, r:%04.2f)" % (mstats['a_dx'],mstats['b_dx'], mstats['corr_dx']),
#             fontsize=fsdxdy,transform=ax.transAxes); ty -= tstep
#    ax.text(tx,ty,"linfit dY = (a:%04.2f, b:%+05.3f km, r:%04.2f)" % (mstats['a_dy'],mstats['b_dy'], mstats['corr_dy']),
#             fontsize=fsdxdy,transform=ax.transAxes); ty -= tstep
    ax.text(tx,ty,"N = %d" % (mstats['N']),
             fontsize=fsdxdy,transform=ax.transAxes); ty -= tstep
    if (status_flag >= 0):
#        txt = 'Only Flag {}'.format(flag_dict[str(status_flag)])
        txt = 'Only Flag {}'.format(str(status_flag))
    else:
        txt = ''
        #txt = 'All Flags'
    ax.text(tx,ty,txt,
             fontsize=fsdxdy,transform=ax.transAxes); ty -= tstep

    ax.xaxis.set_major_locator(MultipleLocator(12.5))
    ax.yaxis.set_major_locator(MultipleLocator(12.5))

    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.tick_params(labelsize=tickfont)

    ax.set_xlim(vmin,vmax)
    ax.set_ylim(vmin,vmax)
#    ax.set_xlabel('Validation data [km]', fontsize=titlefont)
    ax.set_xlabel('Buoy vector [km]', fontsize=titlefont)
    ax.set_ylabel('Product vector [km]', fontsize=titlefont)
    if datestr is None:
        ax.set_title("Scatterplot of drift components\n%s (%s -> %s)" % (m.area.upper(),titled(date1),titled(date2,pos=-1)), fontsize=titlefont)
    else:
        ax.set_title("Scatterplot of drift components\n{} ({})".format(m.area.upper(), datestr), fontsize=titlefont)
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.close()

    return outfile, mstats

def decimal_hours(T):
    return T.hour+T.minute/60.+T.second/3600.

def plot_timecoll_plot(m,date1,date2,outfile,product='na',status_flag=-1):
    prds_dX = m.prd.dX.compressed()
    vals_dX = m.val.dX.compressed()
    prds_dY = m.prd.dY.compressed()
    vals_dY = m.val.dY.compressed()
    diff_t0 = m.coll.diff_t0.compressed()

    hmin  = -9
    hmax  = -hmin

    hours = np.linspace(hmin, hmax, num=2*hmax+1)
    mismatch_msk = np.zeros(hours.size-1).astype('bool')
    mismatch_dX_med = ma.zeros(hours.size-1)
    mismatch_dX_std = ma.zeros(hours.size-1)
    mismatch_dY_med = ma.zeros(hours.size-1)
    mismatch_dY_std = ma.zeros(hours.size-1)

    for i,hrs in enumerate(hours[:-1]):
        lb = timedelta(hours=hours[i])
        lu = timedelta(hours=hours[i+1])
        indx = (diff_t0>=lb) * (diff_t0<lu)
        #print i, lb, lu, indx.sum()
        if sum(indx) == 0:
            mismatch_msk[i] = True
        else:
            mismatch_dX = prds_dX[indx] - vals_dX[indx]
            mismatch_dX_med[i] = mismatch_dX.mean()
            mismatch_dX_std[i] = mismatch_dX.std()
            mismatch_dY = prds_dY[indx] - vals_dY[indx]
            mismatch_dY_med[i] = mismatch_dY.mean()
            mismatch_dY_std[i] = mismatch_dY.std()

    mismatch_dX_med.mask = mismatch_msk
    mismatch_dX_std.mask = mismatch_msk
    mismatch_dY_med.mask = mismatch_msk
    mismatch_dY_std.mask = mismatch_msk

    hours = hours+0.5*(hours[1]-hours[0])
    hours = hours[:-1]

    #polynomial fit
    x_fit = hours[~mismatch_msk]
    y_fit = (mismatch_dX_std[~mismatch_msk] + mismatch_dY_std[~mismatch_msk])*0.5
    c_fit = np.polyfit(x_fit,y_fit,2)
    print(product.upper(), c_fit)
    r_pol = np.poly1d(c_fit)

    mmin = -1
    mmax = 6

    plt.plot([hmin,hmax],[0.,0.],'k-.')
    plt.plot([0,0],[mmin,mmax],'k-.')
    sxx = plt.plot(hours,mismatch_dX_std,ls='-',color=cX,lw=2)
    syy = plt.plot(hours,mismatch_dY_std,ls='-',color=cY,lw=2)
    bxx = plt.plot(hours,mismatch_dX_med,ls='-.',color=cX,lw=2)
    byy = plt.plot(hours,mismatch_dY_med,ls='-.',color=cY,lw=2)
    plt.plot(hours,r_pol(hours),ls='-',color='k',lw=1)
    ax = plt.gca()
    ax.set_xlabel('Time mis-collocation [hours]')
    ax.set_ylabel('Mismatch in drift component [km]')
    ax.set_xlim(hmin,hmax)
    ax.set_ylim(mmin,mmax)

    tx = 0.5 + 0.03
    ty = tstart
    ax.text(tx,ty,product.upper(),transform=ax.transAxes,fontsize=fs); ty -= tstep
    if (status_flag >= 0):
        txt = 'Flag {}:{}'.format(status_flag,flag_dict[status_flag][1])
    else:
        txt = ''
        #txt = 'All Flags'
    ax.text(tx,ty,txt,
             fontsize=fs,transform=ax.transAxes); ty -= tstep

    ax.legend([sxx,syy,bxx,byy],['sdev dX','sdev dY','bias dX','bias dY'],
              scatterpoints=1,loc=4,scatteryoffsets=[0.5],markerscale=2,prop={'size':fs})
    ax.set_title("Time Collocation plot\n%s (%s -> %s)" % (m.area.upper(),titled(date1),titled(date2,pos=-1)))
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

    return outfile

def plot_timecollcorrel_plot(m1,m2,date1,date2,outfile,product1,product2,status_flag=-1):
    diff_t0_1 = m1.coll.diff_t0
    diff_t0_2 = m2.coll.diff_t0

    hmin  = -9
    hmax  = -hmin

    hours = np.linspace(hmin, hmax, num=2*hmax+1)
    nh = hours.size-1
    mismatch_msk = np.zeros((nh,nh)).astype('bool')
    mismatch_dX_corr = ma.zeros((nh,nh))
    mismatch_dY_corr = ma.zeros((nh,nh))

    for i1,hrs1 in enumerate(hours[:-1]):
        lb1 = timedelta(hours=hours[i1])
        lu1 = timedelta(hours=hours[i1+1])
        indx1 = (diff_t0_1>=lb1) * (diff_t0_1<lu1)
        wher1 = np.ma.where(indx1)[0]
        #print diff_t0_1.mask[wher1]
        for i2,hrs2 in enumerate(hours[:-1]):
            lb2 = timedelta(hours=hours[i2])
            lu2 = timedelta(hours=hours[i2+1])
            print('[{}:{}]x[{}:{}]'.format(lb1,lu1,lb2,lu2))
            indx2 = (diff_t0_2>=lb2) * (diff_t0_2<lu2)
            wher2 = np.ma.where(indx2)[0]
            #print diff_t0_2.mask[wher2]
            if sum(indx1) == 0 or sum(indx2) == 0:
                mismatch_msk[i1,i2] = True
            else:
                sm1  = m1.index(wher1)
                sm2  = m2.index(wher2)
                sdt1 = [dt.total_seconds()/60./60. for dt in sm1.coll.diff_t0]
                sdt2 = [dt.total_seconds()/60./60. for dt in sm2.coll.diff_t0]

                #dump_matchup_to_ascii(sm1,'/tmp/ascii_dump_1_b.txt')
                #dump_matchup_to_ascii(sm2,'/tmp/ascii_dump_2_b.txt')
                #print "xxdiff /tmp/ascii_dump_1_b.txt /tmp/ascii_dump_2_b.txt"

                #print "Before get_common: {} {} matchups".format(sm1.count(),sm2.count())
                sm1,sm2 = mio.get_common_matchups(sm1,sm2)
                #print "Now {} {} matchups".format(sm1.count(),sm2.count())

                #dump_matchup_to_ascii(sm1,'/tmp/ascii_dump_1_a.txt')
                #dump_matchup_to_ascii(sm2,'/tmp/ascii_dump_2_a.txt')
                #print "xxdiff /tmp/ascii_dump_1_a.txt /tmp/ascii_dump_2_a.txt"

                mismatch_dX_1 = sm1.prd.dX - sm1.val.dX
                mismatch_dY_1 = sm1.prd.dY - sm1.val.dY
                mismatch_dX_2 = sm2.prd.dX - sm2.val.dX
                mismatch_dY_2 = sm2.prd.dY - sm2.val.dY

                covarM_x_12 = np.cov(mismatch_dX_1,mismatch_dX_2)
                mismatch_dX_corr[i1,i2] = covarM_x_12[0,1] / (sqrt(covarM_x_12[0,0]) * sqrt(covarM_x_12[1,1]))
                covarM_y_12 = np.cov(mismatch_dY_1,mismatch_dY_2)
                mismatch_dY_corr[i1,i2] = covarM_y_12[0,1] / (sqrt(covarM_y_12[0,0]) * sqrt(covarM_y_12[1,1]))

                del sm1
                del sm2

    mismatch_dX_corr.mask = mismatch_msk
    mismatch_dY_corr.mask = mismatch_msk

    hours = hours+0.5*(hours[1]-hours[0])
    hours = hours[:-1]

    mmin = -1
    mmax = 1
    plt.plot([hmin,hmax],[0.,0.],'k-.')
    plt.plot([0,0],[mmin,mmax],'k-.')

    for i1,hrs1 in enumerate(hours[:-1]):
        plt.plot(hours-hours[i1],mismatch_dX_corr[i1,:])

    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def plot_uncertainty_diagram(m,date1,date2,outfile,product='na',status_flag=-1, hemi='nh'):

    legendfont = 16
    labelfont = 20
    titlefont = 22
    tickfont = 20
    percfont = 12

    prds_dX = m.prd.dX.compressed()
    vals_dX = m.val.dX.compressed()
    prds_dY = m.prd.dY.compressed()
    vals_dY = m.val.dY.compressed()
    prds_sX = m.prd.sX.compressed()
    prds_sY = m.prd.sY.compressed()

    if m.area != 'nh' and m.area != 'sh':
        raise ValueError("Know nothing of area {}".format(m.area))

    #mult = 1.0
    #prds_sX *= mult
    #prds_sY *= mult

    mstats = compute_matchup_statistics(m)

    vals_sX = 0.0
    vals_sY = 0.0

    lim_chis = 12.
    #lim_chis = 4.
    chis_dX = ((prds_dX - vals_dX)/(prds_sX+vals_sX))**2
    chis_dX = chis_dX[chis_dX<=lim_chis]
    chi_square_dX = chis_dX.sum() / chis_dX.size
    chis_dY = ((prds_dY - vals_dY)/(prds_sY+vals_sY))**2
    chis_dY = chis_dY[chis_dY<=lim_chis]
    chi_square_dY = chis_dY.sum() / chis_dY.size

    smin = 0.
#    smax = 12.
    smax = 12.
    mmin = 0.
#    mmax = 12.
    mmax = 12.

    if hemi == 'nh':
        smin = 0.
        smax = 4.
        mmin = 0.
        mmax = 4.
    else:
        smin = 0.
        smax = 6.
        mmin = 0.
        mmax = 6.

    edges = np.linspace(smin, smax, num=int((smax * 2)+1))
#    print("edges = ", edges)
#    print("prds_sX = ", prds_sX)

    mismatch_dX_msk = np.zeros(edges.size-1).astype('bool')
    edges_dX = ma.zeros(edges.size-1)
    mismatch_dX_med = ma.zeros(edges.size-1)
    mismatch_dX_std = ma.zeros(edges.size-1)
    mismatch_dX_num = ma.zeros(edges.size-1)
    mismatch_dY_msk = np.zeros(edges.size-1).astype('bool')
    edges_dY = ma.zeros(edges.size-1)
    mismatch_dY_med = ma.zeros(edges.size-1)
    mismatch_dY_std = ma.zeros(edges.size-1)
    mismatch_dY_num = ma.zeros(edges.size-1)
    for i,ed in enumerate(edges[:-1]):
        lb = edges[i]
        ub = edges[i+1]
        indx_sX = (prds_sX>=lb) * (prds_sX<ub)
        mismatch_dX_num[i] = sum(indx_sX)
        if sum(indx_sX) <= 50:
            mismatch_dX_msk[i] = True
        else:
            mismatch_dX = prds_dX[indx_sX] - vals_dX[indx_sX]
            mismatch_dX_med[i] = mismatch_dX.mean()
            mismatch_dX_std[i] = mismatch_dX.std()
            edges_dX[i] = (prds_sX[indx_sX]).mean()

        indx_sY = (prds_sY>=lb) * (prds_sY<ub)
        mismatch_dY_num[i] = sum(indx_sY)
        if sum(indx_sY) <= 50:
            mismatch_dY_msk[i] = True
        else:
            mismatch_dY = prds_dY[indx_sY] - vals_dY[indx_sY]
            mismatch_dY_med[i] = mismatch_dY.mean()
            mismatch_dY_std[i] = mismatch_dY.std()
            edges_dY[i] = (prds_sY[indx_sY]).mean()

    edges_dX = edges_dX[~mismatch_dX_msk]
    mismatch_dX_med = mismatch_dX_med[~mismatch_dX_msk]
    mismatch_dX_std = mismatch_dX_std[~mismatch_dX_msk]
    mismatch_dX_num = mismatch_dX_num[~mismatch_dX_msk]
    edges_dY = edges_dY[~mismatch_dY_msk]
    mismatch_dY_med = mismatch_dY_med[~mismatch_dY_msk]
    mismatch_dY_std = mismatch_dY_std[~mismatch_dY_msk]
    mismatch_dY_num = mismatch_dY_num[~mismatch_dY_msk]

    plt.figure(figsize=(6,6))
    plt.plot([smin,smax],[mmin+vals_sX,+mmax+vals_sX],'k-.')
    plt.plot([smin,smax],[mmin-vals_sX,-mmax-vals_sX],'k-.')
    ax = plt.gca()
    ax.plot([smin,smax],[0.,0.],'k-.')
    if len(edges_dX) > 0 and len(edges_dY) > 0:
        xx = ax.errorbar(edges_dX,mismatch_dX_med,yerr=mismatch_dX_std,fmt='o',color=cX)
        yy = ax.errorbar(edges_dY,mismatch_dY_med,yerr=mismatch_dY_std,fmt='o',color=cY)
        for n in range(len(edges_dX)):
            frac_dX = 100.*mismatch_dX_num[n]/mstats['N']
            frac_dY = 100.*mismatch_dY_num[n]/mstats['N']
            pos = (edges_dX[n],max(mismatch_dX_med[n]+mismatch_dX_std[n],mismatch_dY_med[n]+mismatch_dY_std[n])+0.3)#+0.5)
            ax.text(pos[0],pos[1],'{:.0f}%'.format(frac_dX,),
                    color='k',rotation=90,ha='center',va='bottom',fontsize=percfont)
            #ax.text(edges_dY[n]+0.5,mismatch_dY_med[n]+mismatch_dY_std[n]+0.1,'{:.0f}%'.format(frac_dY,),
            #        color=cY,rotation=+90,ha='left',va='top')
    ty = tstart
    ax.text(tx,ty,product.upper(),transform=ax.transAxes,fontsize=legendfont); ty -= tstep
    ax.text(tx,ty,"N = %d" % (mstats['N']),
             fontsize=legendfont,transform=ax.transAxes); ty -= tstep
    ax.text(tx,ty,r"$\chi^2$ = ({:.2f},{:.2f})".format(chi_square_dX,chi_square_dY),
             fontsize=legendfont,transform=ax.transAxes); ty -= tstep
    if (status_flag >= 0):
        #txt = 'Flag {}:{}'.format(status_flag,flag_dict[status_flag][1])
        txt = 'Flag {}'.format(str(status_flag))
    else:
        txt = ''
        #txt = 'All Flags'
    ax.text(tx,ty,txt,
             fontsize=legendfont,transform=ax.transAxes); ty -= tstep

    ax.tick_params(labelsize=tickfont)
#    ax.set_xlim(smin,smax)
#    ax.set_ylim(-mmax,mmax)
    if hemi == 'nh':
        ax.set_xlim(0,4)
        ax.set_ylim(-4,4)
    else:
        ax.set_xlim(0,6)
        ax.set_ylim(-7,7)
    ax.set_xlabel(r'Product uncertainty (1$\sigma$) [km]', fontsize=labelfont)
#    ax.set_ylabel(r'Validation mismatch (1$\sigma$) [km]', fontsize=labelfont)
    ax.set_ylabel(r'Buoy mismatch (1$\sigma$) [km]', fontsize=labelfont)
    if len(edges_dX) > 0 and len(edges_dY) > 0:
        ax.legend([xx,yy],['dX','dY'],scatterpoints=0,loc=3,scatteryoffsets=[0.5],markerscale=2,prop={'size':legendfont})
    ax.set_title("Uncertainty diagram\n%s (%s -> %s)" % (m.area.upper(),titled(date1),titled(date2,pos=-1)), fontsize=titlefont)

    #np.random.shuffle(chis_dX)
    #chis_dX = chis_dX[:10*(chis_dX.size/10)]
    #for e in range(chis_dX.size/10):
    #    print chis_dX[e*10:(e+1)*10].size
    #rect = [0.0,0.0,0.4,0.4]
    #chis_m = 5
    #ax2 = add_subplot_axes(ax,rect)
    #_, _, patches = ax2.hist(chis_dX, 20, normed=1, range=(0.,chis_m), histtype='stepfilled',facecolor=cX,alpha=0.5)
    #_, _, patches = ax2.hist(chis_dY, 20, normed=1, range=(0.,chis_m), histtype='stepfilled',facecolor=cY,alpha=0.5)
    #ax2.set_ylim(0,chis_m)

    #rect = [0.0,0.0,0.4,0.4]
    plt.savefig(outfile, bbox_inches='tight',dpi=dpi)
    plt.clf()

    return outfile


def _write_one_statistics_line(fh,mstats,pstr,latex=False,last_column=''):

    if mstats['N'] > 0:
        b_x  = "{:+.3f}".format(mstats['bias_dx'],)
        b_y  = "{:+.3f}".format(mstats['bias_dy'],)
        s_x  = "{:.3f}".format(mstats['sdev_dx'],)
        s_y  = "{:.3f}".format(mstats['sdev_dy'],)
        a    = "{:.2f}".format(mstats['a'],)
        b    = "{:.3f}".format(mstats['b'],)
        corr = "{:.2f}".format(mstats['corr'],)
    else:
        b_x  = b_y = s_x = s_y = a = b = corr = '--'

    if latex:
        fh.write("{} & {} & {} & {} & {} & {} & {} & {}\n".format(b_x,b_y,s_x,s_y,a,b,corr,mstats['N']))
    else:
        fh.write(("<tr align=center><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>"+\
                  "<td>{}</td><td>{}</td><td>{}</td><td>{}</td>{}</tr>\n").format(pstr,b_x,b_y,s_x,
                            s_y,a,b,corr,mstats['N'],last_column))

def write_statistics_table(m,date1,date2,outfile,product='na',channels=None, split_in_months=True,reverse_months=False,
                           link_to_figs=None,full_html=False,embed_pics=False,period_string='na',latex=False):
    """
        Compute the matchup statistics (bias, RMSE, correlation, etc...) and
        format them to HTML (or LaTeX) in outfile.
    """
    if latex:
        # de-activate all the linking to pics
        link_to_figs=None
        full_html=False
        embed_pics=False

    if split_in_months:
        sub_periods = mfil.split_period_in_months(date1,date2,reverse=reverse_months)
    else:
        sub_periods = zip(date1, date2);

    # possibly add a last_column with links to figs
    last_column=""
    if link_to_figs is not None:
        print("Handle last column")
        urls = []
        for f in link_to_figs.keys():
            urls.append("<a href=%s>%s</a>" % (link_to_figs[f],f,))
        last_column="<td rowspan=%d>%s</td>" % (len(sub_periods),', '.join(urls))

    nbcolumns = 9
    if len(last_column) > 0:
        nbcolumns += 1

    product_str = product.upper()
    if channels:
        product_str += ' (' + channels.upper() + ')'

    with open(outfile,"w") as fh:
        if full_html:
            fh.write("<!DOCTYPE html>\n<html>\n<body>\n")
            fh.write("<p>\n<table border=1>\n")
            product_line = "<th colspan={}>{}</th>".format(nbcolumns,product_str)
            fh.write("<tr>%s</tr>\n" % (product_line,))
            header_line = "<th>Period</th><th>b(X) [km]</th><th>b(Y) [km]</th>"
            header_line += "<th>s(X) [km]</th><th>s(Y) [km]</th><th>a</th><th>b [km]</th><th>corr</th><th>N</th>"
            if len(last_column) > 0:
                header_line += "<th width=70>figs</th>"
            fh.write("<tr>%s</tr>\n" % (header_line,))

        if split_in_months:
            for sp in sub_periods:
                sav_mask = m.mask
                try:
                    s0 = sp[0]
                    s1 = sp[1]+timedelta(days=1)-timedelta(seconds=0.1)
                    psel = msel.period_selector(s0,s1,'prd.t0')
                    m.apply_selectors(psel)
                    mstats = compute_matchup_statistics(m)
                    pstr = sp[0].strftime("%b%Y")
                    _write_one_statistics_line(fh,mstats,pstr,latex,last_column)
                except Exception as ex:
                    raise Exception("caught ex while computing statistics: {}".format(ex,))

                # last column is only written if on the 1st line of the table (rowspan)
                last_column=''
                # restore mask before datetime selection
                m.mask = sav_mask
                m._applymask()
                sav_mask = m.mask

        else: # NOT split_in_months:
            # we create a combined_selector with OR operator with all the time periods.
            # This is then evaluated only once.
            psels = []
            for sp in sub_periods:
                s0 = sp[0]
                s1 = sp[1]+timedelta(days=1)-timedelta(seconds=0.1)
                psels.append(msel.period_selector(s0,s1,'prd.t0'))

            psel = psels[0]
            if len(psels) >= 2:
                psel = msel.combined_selector(psels[0],psels[1],'or')
                for p in psels[2:]:
                    psel = msel.combined_selector(psel,p,'or')

            m.apply_selectors(psel)
            mstats = compute_matchup_statistics(m)
            _write_one_statistics_line(fh,mstats,period_string,latex,last_column)

        if full_html:
            fh.write("</table>\n</p>\n")

        if embed_pics and link_to_figs is not None:
            fh.write("<p>\n<table border=1>\n")
            for f in link_to_figs.keys():
                fh.write("<tr><td><a href={0}><img src={0} width={1} hspace=0 vspace=0></a></td></tr>\n".format(link_to_figs[f],'75%'))
            fh.write("</table>\n")

        if full_html:
            fh.write("</body>\n</html>\n")

if __name__ == '__main__':
    import doctest
    doctest.testmod()
