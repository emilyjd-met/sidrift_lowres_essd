import os
import re
import math
import fnmatch

import numpy as np
from pyresample import parse_area_file, geometry, image
from netCDF4 import Dataset

from _idcore import ffi, lib

from helper_fcns import fmtdate, corr_time_from_unit, round_sig
from checkinst import checkinst
from write_icedrift_product_file import write_icedrift


def icedrift_wrapper(start_date, end_date, start_dir, end_dir, out_dir,
                     c_logfile, instr, area, channels, radius, rad_neigh):
    """Wrapper to do I/O for C core ice drift code"""

    # Split instrument into instrument and platform and check compatibility
    inst, plat = instr.split('-')
    checkinst(inst, plat)

    # Checking for too many wavebands
    MAX_NBWAVEBANDS = 6
    if len(channels) > MAX_NBWAVEBANDS:
        raise ValueError("Too many wavebands ({}) for the processing. "
                         "Maximum is set to MAX_WAVEBANDS = {}".format(
                         len(channels), MAX_NBWAVEBANDS))

    # Checking the area
    if re.match('nh', area):
        if re.search('ease', area):
            out_area = 'nhease750'
            out_area_name = 'nh-ease2-750'
        else:
            out_area = 'nh625'
            out_area_name = 'nh-polstere-625'
    elif re.match('sh', area):
        if re.search('ease', area):
            out_area = 'shease750'
            out_area_name = 'sh-ease2-750'
        else:
            out_area = 'sh625'
            out_area_name = 'sh-polstere-625'
    else:
        raise ValueError("Unknown input grid {} (should be 'nhXX' or "
                         "'shXX')".format(area))

    # Prepare for wavelength filtering, note that channels are in sorted
    # order in run_cmcc
    chan_sort = sorted(channels)
    wbs = ''
    for i, wb in enumerate(chan_sort):
        if i > 0:
            wbs = wbs + '-'
        m = re.match('(.*)_(.*)', wb)
        if m:
            wbs = wbs + m.group(1)
        else:
            wbs = wbs + (wb)

    # Find the required files
    lncfile_start_pat = ('tc_wght_' + instr + '_' + '*' + wbs + '*' + '_'
                         + area + '_' + fmtdate(start_date) + '.nc')
    lncfile_end_pat = ('tc_wght_' + instr + '_' + '*' + wbs + '*' + '_'
                       + area + '_' + fmtdate(end_date) + '.nc')

    startfiles = fnmatch.filter(os.listdir(start_dir), lncfile_start_pat)
    if not startfiles:
        raise FileNotFoundError('No file matching the pattern {} in directory '
                                '{}'.format(lncfile_start_pat, start_dir))
    endfiles = fnmatch.filter(os.listdir(end_dir), lncfile_end_pat)
    if not endfiles:
        raise FileNotFoundError('No file matching the pattern {} in directory '
                                '{}'.format(lncfile_end_pat, end_dir))

    lncfile_start = os.path.join(start_dir, startfiles[0])
    lncfile_end = os.path.join(end_dir, endfiles[0])
    if not (os.path.exists(lncfile_start) and os.path.isfile(lncfile_start)):
        raise FileNotFoundError('Start file {} does not exist.'.format(
                                lncfile_start))
    elif not (os.path.exists(lncfile_end) and os.path.isfile(lncfile_end)):
        raise FileNotFoundError('End file {} does not exist.'.format(
                                lncfile_end))
    else:
        print("\nDaily maps found ({}):\n".format(instr))
        print("\tDay 1 : {} \n".format(lncfile_start))
        print("\tDay 2 : {} \n".format(lncfile_end))
        print("\n")
    filelist = [lncfile_start, lncfile_end]

    # Reading in the images

    # Array lengths
    filenum = len(filelist)
    wbnum = len(channels)

    # XDIM, YDIM and TDIM are defined this way on the C-side in
    # OSI_HL_Ice/drift/lowres/common/src/icedrift_common.h
    XDIM = 0
    YDIM = 1
    TDIM = 2

    # Shape and coordinate information of input files - these can be
    # taken from one file, but the other should be checked to make sure
    # they are the same
    # A set of parameters from the end file for grid checking

    with Dataset(filelist[1], 'r') as enddataset:

        # Check what the CRS is called
        if 'Polar_Stereographic_Grid' in enddataset.variables:
            crs = 'Polar_Stereographic_Grid'
        elif 'crs' in enddataset.variables:
            crs = 'crs'

        if 'area' in enddataset.ncattrs():
            end_img_area_id = enddataset.area
        elif 'grid' in enddataset.ncattrs():
            end_img_area_id = enddataset.grid
        else:
            end_img_area_id = enddataset['crs'].area_id
        end_img_nx = enddataset.dimensions['xc'].size
        end_img_ny = enddataset.dimensions['yc'].size
        end_img_projstr = enddataset[crs].proj4_string
        end_img_x0 = np.array(enddataset['xc'])[0]
        end_img_y0 = np.array(enddataset['yc'])[0]

    # Most of the grid info is taken from the start file
    with Dataset(filelist[0], 'r') as dataset:
        if 'area' in dataset.ncattrs():
            img_area_id = dataset.area
        elif 'grid' in dataset.ncattrs():
            img_area_id = dataset.grid
        else:
            img_area_id = dataset['crs'].area_id
        img_nx = dataset.dimensions['xc'].size
        img_ny = dataset.dimensions['yc'].size
        img_projstr = dataset[crs].proj4_string
        img_x0 = np.array(dataset['xc'])[0]
        img_y0 = np.array(dataset['yc'])[0]

        # Checking the grids match
        gridcheck = {"grid": [img_area_id, end_img_area_id],
                     "dimensions['xc'].size": [img_nx, end_img_nx],
                     "dimensions['yc'].size": [img_ny, end_img_ny],
                     "[{}].proj4_string".format(crs): [img_projstr,
                                                       end_img_projstr],
                     "['xc'][0]": [img_x0, end_img_x0],
                     "['yc'][0]": [img_y0, end_img_y0]}
        for key in gridcheck:
            if (gridcheck[key][0] != gridcheck[key][1]):
                raise ValueError("Mismatch in grids of start and end files: "
                                 "Parameter {} has value {} in the start "
                                 "file and value {} in the end file".format(
                                 key, gridcheck[key][0], gridcheck[key][1]))


        # Removing problematic parts of the projstr which prevent the C
        # code from parsing it
        ips = img_projstr.split()
        ips = [x for x in ips if not re.match('\+no_defs', x)]
        ips = [x for x in ips if not re.match('\+type=crs', x)]
        img_projstr = ' '.join(ips)

        # Dimensions of grid
        img_size = img_nx * img_ny
        img_dims = np.zeros(3, dtype=np.uint)
        img_dims[XDIM] = img_nx
        img_dims[YDIM] = img_ny
        img_dims[TDIM] = img_size

        # The ice drift code requires the grid values in km, but the output
        # file uses them in m, so stick to m for now
        if dataset['xc'].units == 'km':
            xunits = 1000.
        else:
            xunits = 1.
        if dataset['yc'].units == 'km':
            yunits = 1000.
        else:
            yunits = 1.
        # It is necessary to adjust for the fact that pyresample uses
        # the grid corners and C uses the centres of the pixels at the
        # grid corners
        img_llx = np.array(dataset['xc']).astype(np.float64)[0] * xunits
        img_ury = np.array(dataset['yc']).astype(np.float64)[0] * yunits
        img_urx = np.array(dataset['xc']).astype(np.float64)[-1] * xunits
        img_lly = np.array(dataset['yc']).astype(np.float64)[-1] * yunits
        img_xres = (img_urx - img_llx) / (img_nx - 1)
        img_yres = (img_ury - img_lly) / (img_ny - 1)
        img_llx -= 0.5 * img_xres
        img_ury += 0.5 * img_yres
        img_urx += 0.5 * img_xres
        img_lly -= 0.5 * img_yres

        # Setting up an area definition for the input file
        img_area_def = geometry.AreaDefinition(img_area_id, img_area_id,
                                               img_area_id, img_projstr,
                                               img_nx, img_ny,
                                               (img_llx, img_lly,
                                                img_urx, img_ury))
        # Setting up clean, contiguous arrays in numpy
        lats = np.zeros((img_size), dtype=np.float64)
        lons = np.zeros((img_size), dtype=np.float64)
        images = np.zeros((filenum, wbnum, img_size), dtype=np.float64)
        flags = np.zeros((filenum, wbnum, img_size), dtype=np.short)
        avtimes = np.zeros((filenum, wbnum, img_size), dtype=np.float64)
        iceedges = np.zeros((filenum, img_size), dtype=np.uint8)

        # Input file lats/lons
        lats = np.reshape(np.array(dataset['lat']).astype(np.float64), img_size)
        lons = np.reshape(np.array(dataset['lon']).astype(np.float64), img_size)


    # Cycle over the images and wavebands, reading in the images
    timeref = [None] * 2
    for k, fname in enumerate(filelist):

        # The ice edge mask is only needed per file, not per waveband
        with Dataset(fname, 'r') as dataset:
            iceedges[k] = np.reshape(np.array(dataset['ice_edge']).astype(np.uint8), img_size)

            for j, channel in enumerate(chan_sort):
                # Check waveband name
                m = re.match('(.*)_(.*)', channel)
                if m:
                    wb = m.group(1)
                else:
                    wb = channel

                # NetCDF variables to fetch
                imagename = channel
                images[k][j] = np.reshape(np.array(dataset[imagename]
                                          ).astype(np.float64), img_size)
                flagname = channel + '_flag'
                flagname = channel + '_flag'
                flags[k][j] = np.reshape(np.array(dataset[flagname]
                                          ).astype(np.short), img_size)
                # The average time variable can be constant across wavebands
                if (channel[:-4] + '_avgTime' in dataset.variables):
                    avtimename = channel[:-4] + '_avgTime'
                    avtimes[k][j] = np.reshape(np.array(dataset[avtimename]
                                                    ).astype(np.float64),
                                               img_size)
                elif ('dtime' in dataset.variables):
                    avtimename = 'dtime'
                    avtimes[k][j] = np.reshape(np.array(dataset[avtimename]
                                                    ).astype(np.float64),
                                               img_size)
                elif ('time' in dataset.variables):
                    avtime = dataset['time'][0].astype(np.float64)
                    avtimes[k][j] = avtime
                else:
                    print("WARNING: Issue reading average times")

            # Finding the time since 1970 for each file
            time_unit = dataset['time'].units
            time_val = dataset['time'][0].data.item()
            timeref[k] = corr_time_from_unit(time_val, time_unit)

    # Filling in the rest of the lpattern_radius array
    NBPATTERNS = 2
    pattern_radius = np.zeros((NBPATTERNS), dtype=np.float64)
    pattern_radius[0] = radius
    for i in range(1, NBPATTERNS):
        pradius = pattern_radius[i-1] * 0.5
        if pradius < 20:
            pradius = 20.
        pattern_radius[i] = pradius

    # Compute and set the maxdriftdistance variable
    max_icevelocity = 0.45 * 0.001 * 60. *60. *24. # in km/day
    # Time difference between the two files
    timediff = (end_date - start_date).days
    if timediff < 1:
        timediff = 1
    maxdriftdistance = max_icevelocity * timediff # in km
    sigmoidlength = maxdriftdistance

    # Read grid definition file - it is assumed for now that this is in
    # the local directory
    workdir = os.path.dirname(os.path.abspath(__file__))
    griddeffilename = 'grids_py.def'
    griddeffile = os.path.join(workdir, '../par', griddeffilename)
    if not (os.path.exists(griddeffile) and os.path.isfile(griddeffile)):
        raise FileNotFoundError("{} is not found.".format(griddeffile))
    out_area_def = parse_area_file(griddeffile, out_area_name)[0]

    out_dims = np.zeros(3, dtype=np.uint)
    out_dims[XDIM] = out_area_def.width
    out_dims[YDIM] = out_area_def.height
    out_dims[TDIM] = out_area_def.size
    out_size = out_area_def.size
    (out_llx, out_lly, out_urx, out_ury) = out_area_def.area_extent
    out_projstr = out_area_def.proj_str
    # Removing problematic parts of the projstr which prevent the C
    # code from parsing it
    ops = out_projstr.split()
    ops = [x for x in ops if not re.match('\+no_defs', x)]
    ops = [x for x in ops if not re.match('\+type=crs', x)]
    out_projstr = ' '.join(ops)

    # Calculating the Ax, Ay, Bx, By required for the ice drift code.
    # For the C code, these should be in km
    out_unitsstr = [x for x in out_area_def.proj4_string.split()
                    if '+units=' in x][0][7:]
    if out_unitsstr == 'm':
        out_unitconv = 0.001
    else:
        out_unitconv = 1.0
    out_Ax = round_sig(out_area_def.resolution[0]) * out_unitconv
    out_Ay = round_sig(out_area_def.resolution[1]) * out_unitconv
    # The C grid is half a pixel smaller as it uses pixel centres
    out_Bx = (out_llx + (0.5 * out_area_def.resolution[0])) * out_unitconv
    out_By = (out_ury - (0.5 * out_area_def.resolution[1])) * out_unitconv

    img_unitsstr = [x for x in img_area_def.proj4_string.split()
                    if '+units=' in x][0][7:]
    if img_unitsstr == 'm':
        img_unitconv = 0.001
    else:
        img_unitconv = 1.0

    img_Ax = round_sig(img_area_def.resolution[0]) * img_unitconv
    img_Ay = round_sig(img_area_def.resolution[1]) * img_unitconv
    # The C grid is half a pixel smaller as it uses pixel centres
    img_Bx = (img_llx + (0.5 * img_area_def.resolution[0])) * img_unitconv
    img_By = (img_ury - (0.5 * img_area_def.resolution[1])) * img_unitconv

    # Create the output latitude and longitude arrays
    olons, olats = out_area_def.get_lonlats()
    out_lats = np.reshape(np.array(olats).astype(np.float64), out_size)
    out_lons = np.reshape(np.array(olons).astype(np.float64), out_size)

    # Creating an array of output image indices
    owcs = np.zeros((out_size), dtype=np.uint32)
    owcs = np.array(range(out_size)).astype(np.uint32)

    # Calculating iwcs per individual pixel lat/lon in the same style
    # as the original C code. This results in a regularly spaced array,
    # but does not precisely match the original C code output (which
    # has some irregularities)
    # First find the x and y indices corresponding to the 1-D indices

#    print("======================================")
#    print("Input grid: ")
#    print("x,y,lat,lon = ", 0, 0, img_area_def.get_lonlat(0, 0))
#    print("x,y,lat,lon = ", 0, 607, img_area_def.get_lonlat(0, 607))
#    print("x,y,lat,lon = ", 895, 0, img_area_def.get_lonlat(895, 0))
#    print("x,y,lat,lon = ", 895, 607, img_area_def.get_lonlat(895, 607))
#    print("======================================")
#    print("Output grid: ")
#    print("x,y,lat,lon = ", 0, 0, out_area_def.get_lonlat(0, 0))
#    print("x,y,lat,lon = ", 0, 118, out_area_def.get_lonlat(0, 118))
#    print("x,y,lat,lon = ", 176, 0, out_area_def.get_lonlat(176, 0))
#    print("x,y,lat,lon = ", 176, 118, out_area_def.get_lonlat(176, 118))
#    print("======================================")

    iwcs = np.zeros((out_size), dtype=np.uint32)
    for i in range(len(owcs)):
        row = math.floor(owcs[i] / out_dims[XDIM])
        col = int(owcs[i] - (row * out_dims[XDIM]))
        lon = out_lons[i]
        lat = out_lats[i]
        (newrow, newcol) = img_area_def.get_xy_from_lonlat(lon, lat)
        try:
            iwcs[i] = int((newcol * img_nx) + newrow)
        except:
            iwcs[i] = int(out_size)

#    # Alternative method of calculating iwcs. This results in an
#    # array with some irregularities.
#    # Relating the input image indices to the output image indices
#    # according to the area definitions
#    iwcs = np.zeros((out_size), dtype=np.uint32)
#    iwcs_inds = np.array([x for x in range(img_size)], dtype=np.uint32)
#    iwcs_data = np.reshape(iwcs_inds, (img_ny, img_nx))
#    # Radius of influence chosen at 10km, 1km is too small but 5km looks
#    # OK, 10km chosen for safety margin
#    iwcs_con = image.ImageContainerNearest(iwcs_data, img_area_def,
#                                           radius_of_influence=10000)
#    iwcs_2d = iwcs_con.resample(out_area_def).image_data
#    iwcs = np.array(iwcs_2d).flatten()

    # Create output arrays for the drift
    driftx = np.zeros((out_size), dtype=np.single)
    drifty = np.zeros((out_size), dtype=np.single)
    t0 = np.zeros((out_size), dtype=np.int_)
    t1 = np.zeros((out_size), dtype=np.int_)
    flag = np.zeros((out_size), dtype=np.short)
    sigdx = np.zeros((out_size), dtype=np.single)
    sigdy = np.zeros((out_size), dtype=np.single)
    corrdxdy = np.zeros((out_size), dtype=np.single)
    uflag = np.zeros((out_size), dtype=np.short)
    length = np.zeros((out_size), dtype=np.single)
    direction = np.zeros((out_size), dtype=np.single)
    outlonb = np.zeros((out_size), dtype=np.single)
    outlatb = np.zeros((out_size), dtype=np.single)
    outlone = np.zeros((out_size), dtype=np.single)
    outlate = np.zeros((out_size), dtype=np.single)
    outfc = np.zeros((out_size), dtype=np.single)
    snavg = np.zeros((out_size), dtype=np.short)
    avgx = np.zeros((out_size), dtype=np.single)
    avgy = np.zeros((out_size), dtype=np.single)
    length_avg = np.zeros((out_size), dtype=np.single)
    length_diff = np.zeros((out_size), dtype=np.single)
    stdx = np.zeros((out_size), dtype=np.single)
    stdy = np.zeros((out_size), dtype=np.single)
    patternindex = np.zeros((out_size), dtype=np.short)

    # Prepare and set the globals

    # Define an allocator
    allocator = ffi.new_allocator(lib.malloc, None)

    # Define a list to avoid freeing memory when the cffi python
    # container goes out of scope
    owned_memory = []

    # STRINGS ====================================================

    iproj_bytestr = img_projstr.encode()
    iproj = allocator("char[]", iproj_bytestr)
    owned_memory.append(iproj)
    lib.img_projstr = iproj

    om_bytestr = 'CC'.encode()
    om = allocator("char[]", om_bytestr)
    owned_memory.append(om)
    lib.OptimMetric = om

    oarea_bytestr = out_area.encode()
    oarea = allocator("char[]", oarea_bytestr)
    owned_memory.append(oarea)
    lib.out_area = oarea

    oproj_bytestr = out_projstr.encode()
    oproj = allocator("char[]", oproj_bytestr)
    owned_memory.append(oproj)
    lib.out_projstr = oproj

    clogf_bytestr = c_logfile.encode()
    clogf = allocator("char[]", clogf_bytestr)
    owned_memory.append(clogf)
    lib.reportFile = clogf

    # SINGLE NUMERICAL VALUES =====================================

    # Input image dimensions. Ax and Ay are pixel size, and Bx and By
    # are the coords of the upper left pixel
    lib.img_Ax = img_Ax
    lib.img_Bx = img_Bx
    lib.img_Ay = img_Ay
    lib.img_By = img_By

    # Number of wavebands
    lib.nbWaveBands = wbnum

    # Weighting, pass 1 for now (placeholder in ice-drift code)
    lib.twghtStart = 1
    lib.twghtEnd = 1

    lib.maxdriftdistance = maxdriftdistance
    lib.sigmoid_length = sigmoidlength

    # Number of points in output grid
    lib.NDRIFTPIXELS = out_size

    # Output image dimensions. Ax and Ay are pixel size, and Bx and By
    # are the coords of the upper left pixel
    lib.out_Ax = out_Ax
    lib.out_Bx = out_Bx
    lib.out_Ay = out_Ay
    lib.out_By = out_By

    # Radius to look for neighbours
    lib.radiusNeighbours = rad_neigh

    # DATA ARRAYS =================================================

    # Start and end images
    for k in range(filenum):
        obs_k = allocator("double*[]", wbnum)
        owned_memory.append(obs_k)
        lib.obs[k] = obs_k
        for j in range(wbnum):
            obs_kj = allocator("double[]", img_size)
            owned_memory.append(obs_kj)
            lib.obs[k][j] = obs_kj
            for i in range(img_size):
                lib.obs[k][j][i] = images[k][j][i]

    # Flag for images
    for k in range(filenum):
        tcflag_k = allocator("short*[]", wbnum)
        owned_memory.append(tcflag_k)
        lib.TCflag[k] = tcflag_k
        for j in range(wbnum):
            tcflag_kj = allocator("short[]", img_size)
            owned_memory.append(tcflag_kj)
            lib.TCflag[k][j] = tcflag_kj
            for i in range(img_size):
                lib.TCflag[k][j][i] = flags[k][j][i]

    # Ice edge mask
    for k in range(filenum):
        iceedges_k = allocator("short[]", img_size)
        owned_memory.append(iceedges_k)
        lib.icelandmask[k] = iceedges_k
        for i in range(img_size):
            lib.icelandmask[k][i] = iceedges[k][i]

    # Input and output image dimensions
    # The size of these have already been defined in the build, so no
    # need to save the memory
    for k in range(3):
        lib.img_dims[k] = img_dims[k]
    for k in range(3):
        lib.out_dims[k] = out_dims[k]

    # Latitude and longitude of input image
    img_lat_save = allocator("double[]", img_size)
    owned_memory.append(img_lat_save)
    lib.img_lat = img_lat_save
    for i in range(img_size):
        lib.img_lat[i] = lats[i]
    img_lon_save = allocator("double[]", img_size)
    owned_memory.append(img_lon_save)
    lib.img_lon = img_lon_save
    for i in range(img_size):
        lib.img_lon[i] = lons[i]

    # Latitude and longitude of output image
    olat_save = allocator("double[]", out_size)
    owned_memory.append(olat_save)
    lib.olat = olat_save
    for i in range(out_size):
        lib.olat[i] = out_lats[i]
    olon_save = allocator("double[]", out_size)
    owned_memory.append(olon_save)
    lib.olon = olon_save
    for i in range(out_size):
        lib.olon[i] = out_lons[i]

    # Starting pattern radius
    for k in range(NBPATTERNS):
        lib.pattern_radius[k] = pattern_radius[k]

    # Input and output world coordinate system extraction indices
    iwcs_save = allocator("unsigned int[]", out_size)
    owned_memory.append(iwcs_save)
    lib.iwcs = iwcs_save
    for i in range(out_size):
        lib.iwcs[i] = iwcs[i]
    owcs_save = allocator("unsigned int[]", out_size)
    owned_memory.append(owcs_save)
    lib.owcs = owcs_save
    for i in range(out_size):
        lib.owcs[i] = owcs[i]

    # Output arrays
    driftX_save = allocator("float[]", out_size)
    owned_memory.append(driftX_save)
    lib.driftX = driftX_save
    driftY_save = allocator("float[]", out_size)
    owned_memory.append(driftY_save)
    lib.driftY = driftY_save
    pflag_save = allocator("short[]", out_size)
    owned_memory.append(pflag_save)
    lib.pflag = pflag_save
    sigdX_save = allocator("float[]", out_size)
    owned_memory.append(sigdX_save)
    lib.sigdX = sigdX_save
    sigdY_save = allocator("float[]", out_size)
    owned_memory.append(sigdY_save)
    lib.sigdY = sigdY_save
    corrdXdY_save = allocator("float[]", out_size)
    owned_memory.append(corrdXdY_save)
    lib.corrdXdY = corrdXdY_save
    uflag_save = allocator("short[]", out_size)
    owned_memory.append(uflag_save)
    lib.uflag = uflag_save
    length_save = allocator("float[]", out_size)
    owned_memory.append(length_save)
    lib.length = length_save
    dir_save = allocator("float[]", out_size)
    owned_memory.append(dir_save)
    lib.dir = dir_save
    outlonB_save = allocator("float[]", out_size)
    owned_memory.append(outlonB_save)
    lib.outlonB = outlonB_save
    outlatB_save = allocator("float[]", out_size)
    owned_memory.append(outlatB_save)
    lib.outlatB = outlatB_save
    outlonE_save = allocator("float[]", out_size)
    owned_memory.append(outlonE_save)
    lib.outlonE = outlonE_save
    outlatE_save = allocator("float[]", out_size)
    owned_memory.append(outlatE_save)
    lib.outlatE = outlatE_save
    outfc_save = allocator("float[]", out_size)
    owned_memory.append(outfc_save)
    lib.outfc = outfc_save
    snavg_save = allocator("short[]", out_size)
    owned_memory.append(snavg_save)
    lib.snavg = snavg_save
    avgX_save = allocator("float[]", out_size)
    owned_memory.append(avgX_save)
    lib.avgX = avgX_save
    avgY_save = allocator("float[]", out_size)
    owned_memory.append(avgY_save)
    lib.avgY = avgY_save
    length_avg_save = allocator("float[]", out_size)
    owned_memory.append(length_avg_save)
    lib.length_avg = length_avg_save
    length_diff_save = allocator("float[]", out_size)
    owned_memory.append(length_diff_save)
    lib.length_diff = length_diff_save
    stdX_save = allocator("float[]", out_size)
    owned_memory.append(stdX_save)
    lib.stdX = stdX_save
    stdY_save = allocator("float[]", out_size)
    owned_memory.append(stdY_save)
    lib.stdY = stdY_save
    patternIndex_save = allocator("short[]", out_size)
    owned_memory.append(patternIndex_save)
    lib.patternIndex = patternIndex_save

    # Call the model
    print("Calling the C core code...")
    _ = lib.core()

    # Putting the output back into numpy arrays
    for i in range(out_size):
        driftx[i] = lib.driftX[i]
        drifty[i] = lib.driftY[i]
        flag[i] = lib.pflag[i]
        sigdx[i] = lib.sigdX[i]
        sigdy[i] = lib.sigdY[i]
        corrdxdy[i] = lib.corrdXdY[i]
        uflag[i] = lib.uflag[i]
        length[i] = lib.length[i]
        direction[i] = lib.dir[i]
        outlonb[i] = lib.outlonB[i]
        outlatb[i] = lib.outlatB[i]
        outlone[i] = lib.outlonE[i]
        outlate[i] = lib.outlatE[i]
        outfc[i] = lib.outfc[i]
        snavg[i] = lib.snavg[i]
        avgx[i] = lib.avgX[i]
        avgy[i] = lib.avgY[i]
        length_avg[i] = lib.length_avg[i]
        length_diff[i] = lib.length_diff[i]
        stdx[i] = lib.stdX[i]
        stdy[i] = lib.stdY[i]
        patternindex[i] = lib.patternIndex[i]

    # Creating the output times, masking according to the ice drift
    UNDEFNC_FLOAT = -1.0e10
    UNDEFNC_INT = -1e10
    drift_mask = np.zeros(out_size)
    t0_vals = np.full((out_size), UNDEFNC_INT, dtype=np.intc)
    t1_vals = np.full((out_size), UNDEFNC_INT, dtype=np.intc)
    for i in range(out_size):
        if driftx[i] == UNDEFNC_FLOAT:
            drift_mask[i] = 1
        else:
            t0_vals[owcs[i]] = avtimes[0][0][iwcs[i]]
            t1_vals[owcs[i]] = avtimes[1][0][iwcs[i]]
    t0 = np.ma.masked_array(t0_vals, mask=drift_mask)
    t1 = np.ma.masked_array(t1_vals, mask=drift_mask)

    # =================================================================
    # Writing out
    w_driftx = driftx.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_drifty = drifty.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_t0 = t0.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_t1 = t1.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_flag = flag.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_sigdx = sigdx.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_sigdy = sigdy.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_corrdxdy = corrdxdy.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_uflag = uflag.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_length = length.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_direction = direction.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_outlonb = outlonb.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_outlatb = outlatb.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_outlone = outlone.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_outlate = outlate.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_outfc = outfc.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_snavg = snavg.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_avgx = avgx.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_avgy = avgy.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_length_avg = length_avg.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_length_diff = length_diff.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_stdx = stdx.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_stdy = stdy.reshape((out_dims[YDIM], out_dims[XDIM]))
    w_patternindex = patternindex.reshape((out_dims[YDIM], out_dims[XDIM]))

    fname = write_icedrift(out_area_def,
                           inst, plat, channels, start_date, end_date,
                           out_area_name, out_projstr, out_dir,
                           w_driftx, w_drifty, w_t0, w_t1, w_flag,
                           w_sigdx, w_sigdy, w_corrdxdy,
                           w_uflag, w_length, w_direction,
                           w_outlonb, w_outlatb, w_outlone, w_outlate,
                           w_outfc, w_snavg, w_avgx, w_avgy,
                           w_length_avg, w_length_diff,
                           w_stdx, w_stdy, w_patternindex)

    return fname
