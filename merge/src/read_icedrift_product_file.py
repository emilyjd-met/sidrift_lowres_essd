from netCDF4 import Dataset
from pyresample import geometry, utils, AreaDefinition
import numpy as np
from numpy import ma
import cartopy
import cartopy.crs as ccrs

def read_drift(driftfile, short=False, multi=False):
    '''Reading the drift data from a level 3 product file'''

    var_list = ['driftX', 'driftY', 't0', 't1', 'flag', 'sigdX', 'sigdY',
                'corrdXdY', 'uflag', 'length', 'dir', 'latStart', 'lonStart',
                'latEnd', 'lonEnd', 'corr', 'pattern']
    if not multi:
        var_list.extend(['fsX', 'fsY', 'fsX_12utc', 'fsY_12utc'])
    if not short:
        var_list.extend(['navg', 'avgX', 'avgY', 'avglen', 'difflen', 'stdX',
                         'stdY'])
    driftdata = {}

    with Dataset(driftfile, 'r') as dataset:
        driftdata['grid_mapping'] = dataset.variables['driftX'].__dict__['grid_mapping']
        driftdata['proj4_string'] = dataset.variables[driftdata['grid_mapping']].__dict__['proj4_string']
        try:
            driftdata['proj_dict'] = utils._proj4.proj4_str_to_dict(driftdata['proj4_string'])
        except:
            driftdata['proj_dict'] = utils.proj4.proj4_str_to_dict(driftdata['proj4_string'])
        driftdata['xc'] = dataset['xc'][:]
        driftdata['yc'] = dataset['yc'][:]
        driftdata['gridname'] = dataset.remapping_gridname

        driftdata['lon'] = dataset.variables['lon'][:]
        driftdata['lat'] = dataset.variables['lat'][:]

        for var in var_list:
            fv = dataset.variables[var]._FillValue
            vardata = dataset[var][:]
            vardata = vardata[0, :, :]
            driftdata[var] = ma.array(vardata, fill_value=fv)
            driftdata[var].mask = driftdata[var].data == fv
            fvvar = '{}_fv'.format(var)
            driftdata[fvvar] = fv

        # Read the global variables
        driftdata['globals'] = {}
        for item in dataset.__dict__:
            driftdata['globals'][item] = dataset.__dict__[item]


    # Grid spacing
    sorted_xc = np.sort(driftdata['xc'])
    sorted_yc = np.sort(driftdata['yc'])
    smallest_xc = sorted_xc[0]
    second_smallest_xc = sorted_xc[1]
    smallest_yc = sorted_yc[0]
    second_smallest_yc = sorted_yc[1]
    driftdata['ax'] = second_smallest_xc - smallest_xc
    driftdata['ay'] = second_smallest_yc - smallest_yc

    scale = 1000.
    llx = (driftdata['xc'][0] - (0.5 * driftdata['ax'])) * scale
    lly = (driftdata['yc'][-1] - (0.5 * driftdata['ay'])) * scale
    urx = (driftdata['xc'][-1] + (0.5 * driftdata['ax'])) * scale
    ury = (driftdata['yc'][0] + (0.5 * driftdata['ay'])) * scale
    driftdata['area_extent'] = (float(llx),
                                float(lly),
                                float(urx),
                                float(ury))
    driftdata['mpl_extent'] = (driftdata['area_extent'][0],
                            driftdata['area_extent'][2],
                            driftdata['area_extent'][3],
                            driftdata['area_extent'][1])
    # Shifting by half a pixel for divergences and convergences, which
    # are calculated between pixels
    driftdata['mpl_extent_divs'] = (driftdata['area_extent'][0]
                                 + (0.5 * 1000. * driftdata['ax']),
                                 driftdata['area_extent'][2]
                                 + (0.5 * 1000. * driftdata['ax']),
                                 driftdata['area_extent'][3]
                                 + (0.5 * 1000. * driftdata['ay']),
                                 driftdata['area_extent'][1]
                                 + (0.5 * 1000. * driftdata['ay']))
    driftdata['area_def'] = AreaDefinition(driftdata['gridname'],
                                           driftdata['gridname'],
                                           driftdata['gridname'],
                                           driftdata['proj_dict'],
                                           driftdata['xc'].shape[0],
                                           driftdata['yc'].shape[0],
                                           driftdata['area_extent'])
    driftdata['data_crs'] = driftdata['area_def'].to_cartopy_crs()


    if driftdata['grid_mapping'] == 'Polar_Stereographic_Grid':
        data_globe = ccrs.Globe(semimajor_axis=driftdata['proj_dict']['a'],
                                semiminor_axis=driftdata['proj_dict']['b'])
        if driftdata['lat'][0, 0] > 0:
            driftdata['data_ccrs'] = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                                        globe=data_globe)
        else:
            driftdata['data_ccrs'] = ccrs.SouthPolarStereo(central_longitude=0.0,
                                                        globe=data_globe)

    elif driftdata['grid_mapping'] in ['LambertAzimuthalEqualArea',
                                       'Lambert_Azimuthal_Equal_Area',
                                       'Lambert_Azimuthal_Grid']:
        if driftdata['lat'][0, 0] > 0:
            driftdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=90,
                false_easting=0, false_northing=0)
        else:
            driftdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=-90,
                false_easting=0, false_northing=0)
    else:
        raise ValueError("Unrecognised grid mapping {}".format(driftdata['grid_mapping']))


    return driftdata
