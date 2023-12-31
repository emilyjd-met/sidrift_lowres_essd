import os
import sys
import math
import argparse
from argparse import RawDescriptionHelpFormatter
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib
import matplotlib.style
matplotlib.use('Agg')
matplotlib.style.use('classic')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyresample import geometry, utils, AreaDefinition
import cmocean

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../format/src'))
from grid_info import region_params

tickfont = 20


def parse_args():

    p = argparse.ArgumentParser("plot_ncparam",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--infile', required=True,
                   help="Input file with field to plot")
    p.add_argument('-p', '--param', required=True,
                   help="Parameter to plot")
    p.add_argument('-p2', '--param2', required=False, default=None,
                   help="Second parameter to average with first parameter")
    p.add_argument('-o', '--outdir', required=True,
                   help="Output plot directory")
    p.add_argument('-c', '--cline', action='store_true', default=False,
                   help="Plot a contour line around the un-gapfilled data")
    args = p.parse_args()

    return args


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def cmap_norm_to_lims(cmap, mapmin, mapmax, vmin, vmax):
    '''Return a colourmap that will be normalised from mapmin to mapmax (so
    as to get zero in the right place) but which is cropped to vmin and vmax'''

    if vmin < mapmin or vmax > mapmax:
        raise ValueError("vmin and vmax must be within the bounds of mapmin and mapmax")
    minval = (vmin - mapmin) / (mapmax - mapmin)
    maxval = (vmax - mapmin) / (mapmax - mapmin)
    new_cmap = truncate_colormap(cmap, minval=minval, maxval=maxval)

    return new_cmap


def colbar_flags_discrete():

    # The flags we were interested in are
    # 0 - Good data
    # 1 - Gapfilled
    # 2 - Land
    # 3 - No ice
    # 4 - No drift vectors

    colorlist = ["firebrick", "goldenrod", "grey", "midnightblue", "teal"]
    cmap_flags = colors.ListedColormap(colorlist)

    return cmap_flags


def plot_ncparam(infile, param, param2=None, outdir='.', cline=False):

    with Dataset(infile, 'r') as dataset:
        try:
            grid_mapping = dataset.variables[param].__dict__['grid_mapping']
        except:
            grid_mapping = dataset.variables['A_real'].__dict__['grid_mapping']

        proj4_string = dataset.variables[grid_mapping].__dict__['proj4_string']
        data_proj_dict = utils.proj4.proj4_str_to_dict(proj4_string)
#        data_proj_dict = utils._proj4.proj4_str_to_dict(proj4_string)
        xc = dataset['xc']
        yc = dataset['yc']
        flip = False
        # Flip if necessary
        if yc[-1] < yc[0]:
            flip = True
            yc = np.flip(yc)
        data_area_extent = (float(xc[0]*1000.),
                            float(yc[-1]*1000.),
                            float(xc[-1]*1000.),
                            float(yc[0]*1000.))
        mpl_data_extent = (data_area_extent[0],
                           data_area_extent[2],
                           data_area_extent[3],
                           data_area_extent[1])
        #print("mpl_data_extent = ", mpl_data_extent)
        data_area_def = AreaDefinition('data', 'data', 'data',
                                       data_proj_dict,
                                       xc.shape[0],
                                       yc.shape[0],
                                       data_area_extent)

        if grid_mapping in ['Polar_Stereographic_Grid',
                            'projection_stere']:
            data_globe = ccrs.Globe(semimajor_axis=dataset['proj_dict']['a'],
                                    semiminor_axis=dataset['proj_dict']['b'])
            if dataset['lat'][0, 0] > 0:
                data_ccrs = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                                  globe=data_globe)
                region = 'nh'
                hemi = 'nh'
            else:
                data_ccrs = ccrs.SouthPolarStereo(central_longitude=0.0,
                                                  globe=data_globe)
                region = 'sh'
                hemi = 'sh'
        elif grid_mapping in ['LambertAzimuthalEqualArea',
                              'Lambert_Azimuthal_Equal_Area',
                              'Lambert_Azimuthal_Grid',
                              'projection_laea']:
            if dataset['lat'][0, 0] > 0:
                data_ccrs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                           central_latitude=90,
                                                           false_easting=0, 
                                                           false_northing=0)
                region = 'nh-ease-wide'
                hemi = 'nh'
            else:
                data_ccrs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                           central_latitude=-90,
                                                           false_easting=0,
                                                           false_northing=0)
                region = 'sh-ease-wide'
                hemi = 'sh'
        else:
            raise ValueError("Unrecognised grid mapping {}".format(
                grid_mapping))

        data = dataset.variables[param][:, :]
        dmask = dataset.variables[param][:, :].mask
        if param2 is not None:
            data2 = dataset.variables[param2][:, :]
            dmask2 = dataset.variables[param2][:, :].mask

        # Oh, no, some of the data is FLOAT nans, not numpy nans...
        if dmask.sum() == 0:
            dmask_new = np.zeros_like(data, dtype=np.int8)
            for iy, ix in np.ndindex(data.shape):
                if math.isnan(data[iy, ix]):
                    dmask_new[iy, ix] = 1
            dmask = dmask_new
        if param2 is not None:
            if dmask2.sum() == 0:
                dmask2_new = np.zeros_like(data2, dtype=np.int8)
                for iy, ix in np.ndindex(data2.shape):
                    if math.isnan(data2[iy, ix]):
                        dmask2_new[iy, ix] = 1
                dmask2 = dmask2_new
            dmask = np.logical_and(dmask, dmask2)

        if param2 is not None:
            data = np.mean(np.array([data, data2]), axis=0)
            data = np.ma.array(data)

        if param != 'flags':
            data[dmask] = np.nan
            data.mask = dmask
            count = data.count()
            thresh = math.ceil(count / 100.)
            data1d = data.flatten()
            data1dclean = [x for x in data1d if x != np.nan]
            datasort = sorted(data1dclean)
            print("-- min(data) = ", np.nanmin(data))
            print("-- max(data) = ", np.nanmax(data))
            print("-- lower threshold data (99%) = ", datasort[thresh - 1])
            print("-- upper threshold data (99%) = ", datasort[-thresh])

        if param == 'flags':
            # Make a plotting array with
            # 0 - OK (0)
            # 1 - Gapfill (19)
            # 2 - Land (3)
            # 3 - No ice (4)
            # 4 - No drift (20)
            data[data == 19] = 1
            data[data == 3] = 2
            data[data == 4] = 3
            data[data == 20] = 4

        # Information for the contour line
        ff = dataset.variables['flags'][:, :]
        ffmask = dataset.variables['flags'][:, :].mask
        ff[ff == 0] = 100
        ff[ff != 100] = 0
        if ffmask:
            ff.mask = False
            ff[ffmask] = 0
        if flip:
            ff = np.flip(ff, axis=0)
            
    # Plotting setup
    colbar_type = 'cont'
    inv_yax = False
    if param in  ['A_real', 'A_imag', 'A_real_gapfill', 'A_imag_gapfill']:
        datacmap = cmocean.cm.haline
        vmin = -0.02
        vmax = 0.05
        data_cmap_lvl  = [-0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05]
        norm = colors.TwoSlopeNorm(vmin=-0.05, vcenter=0, vmax=vmax)
    if param in  ['C_real', 'C_imag', 'C_real_gapfill', 'C_imag_gapfill']:
        datacmap = cmocean.cm.balance
        vmin = -0.06
        vmax = 0.06
        data_cmap_lvl  = [-0.06, -0.03, 0.0, 0.03, 0.06]
        norm = colors.TwoSlopeNorm(vmin=-0.25, vcenter=0, vmax=vmax)
    if param in  ['RMS_Res_real', 'RMS_Res_imag']:
        datacmap = cmocean.cm.haline
        vmin = 0.0
        vmax = 0.1
        data_cmap_lvl  = [0, 0.02, 0.04, 0.06, 0.08, 0.1]
        norm = None
    if param in  ['absA', 'absA_gapfill']:
        datacmap = cmocean.cm.haline
        vmin = 0.00
        vmax = 0.03
        data_cmap_lvl  = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
        norm = None
    if param in  ['ThetaA', 'ThetaA_gapfill']:
        if hemi == 'nh':
#            datacmap = cmocean.cm.balance
#            vmin = -60.
#            vmax = 20.
#            data_cmap_lvl  = [-60, -50, -40, -30, -20, -10, 0, 10, 20]
#            datacmap = cmocean.cm.balance
            vmin = -60.
            vmax = 20.
            data_cmap_lvl  = [-60, -50, -40, -30, -20, -10, 0, 10, 20]
            norm = None
            datacmap = cmap_norm_to_lims(cmocean.cm.balance, -60, 60, 
                                         vmin, vmax)
        else:
#            datacmap = cmocean.cm.balance
#            vmin = 10.
#            vmax = 50.
#            data_cmap_lvl  = [10, 20, 30, 40, 50]
#            datacmap = cmocean.cm.balance_r
#            vmin = -20.
#            vmax = 60.
#            data_cmap_lvl  = [-20, -10, 0, 10, 20, 30, 40, 50, 60]
            vmin = -20.
            vmax = 60.
            data_cmap_lvl  = [-20, -10, 0, 10, 20, 30, 40, 50, 60]
            norm = None
            inv_yax = True
            datacmap = cmap_norm_to_lims(cmocean.cm.balance_r, -60, 60, vmin, vmax)
    if param in  ['lon']:
        datacmap = cmocean.cm.balance
        vmin = -180
        vmax = 180
        data_cmap_lvl  = [-150, -100, -50, 0, 50, 100, 150]
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    if param in  ['flags']:
        datacmap = colbar_flags_discrete()
        vmin = 0
        vmax = 5
        data_cmap_lvl  = [0.5, 1.5, 2.5, 3.5, 4.5]
        cmap_labs = ['Param', 'Gapfilled', 'Land', 'No ice', 'No drift']
        norm = None
        colbar_type = 'discrete'

    colbar_label = ""
    
    # Fetching the region parameters from the file and converting 
    # from lat/lon
    pc = ccrs.PlateCarree()
    rp = region_params(region)
    llx, lly = data_ccrs.transform_point(rp['lllon'], rp['lllat'], 
                                         src_crs=pc)
    urx, ury = data_ccrs.transform_point(rp['urlon'], rp['urlat'], 
                                         src_crs=pc)
    plot_data_extent = [llx, urx, lly, ury]

    # Setting up the plot
    colbar = True
    if colbar:
        plt.figure(figsize=(8, 7))
    else:
        plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=data_ccrs)
    ax.set_extent(plot_data_extent, crs=data_ccrs)

    # Plotting the data
    datacbar = ax.imshow(data, cmap=datacmap, extent=mpl_data_extent,
                         transform=data_ccrs, vmin=vmin, vmax=vmax,
                         norm=norm, interpolation='none')
    # Adding a contour line
    if cline:
        cs = ax.contour(ff, [95], colors='red', linewidths=2, 
                        transform=data_ccrs, extent=mpl_data_extent)

    if colbar:
        if colbar_type == 'discrete':
            print("Am I here?")
            cb = plt.colorbar(datacbar, ticks=data_cmap_lvl,
                         orientation='vertical', pad=0.05, shrink=0.8)
            cb.ax.set_yticklabels(cmap_labs)
            cb.set_label(colbar_label)
            cb.ax.tick_params(labelsize=tickfont)
        else:
            print("HERE!!!!!!")
            plt.yticks(fontsize=tickfont)
            plt.colorbar(datacbar, ticks=data_cmap_lvl,
                         orientation='vertical', pad=0.05, shrink=0.8
                     ).set_label(colbar_label)
            #ticklabs = cb.ax.get_yticklabels()
            #cbar.ax.set_yticklabels(ticklabs, fontsize=tickfont)
            #plt.xticks(fontsize=tickfont)
            #plt.xticks(fontsize=tickfont)

    if param != 'flags':
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land',
                                                    '50m', edgecolor='grey',
                                                    facecolor='grey'))
        if hemi == 'sh':
            ax.add_feature(cfeature.NaturalEarthFeature('physical',
                    'antarctic_ice_shelves_polys', '10m', facecolor='grey',
                                                        edgecolor='grey'))

    if param2 is not None:
        param = '{}_{}'.format(param, param2)
    outbname = "{}_{}.png".format(os.path.basename(infile)[:-3], param)
    outname = os.path.join(outdir, outbname)
    plt.savefig(outname, bbox_inches='tight')
    plt.close()
    print("Figure is in {}".format(outname))


if __name__ == '__main__':

    args = parse_args()
    infile = args.infile
    param = args.param
    param2 = args.param2
    outdir = args.outdir
    cline = args.cline

    plot_ncparam(infile, param, param2=param2, outdir=outdir, cline=cline)
