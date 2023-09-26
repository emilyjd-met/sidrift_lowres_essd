import os
import re
import datetime as dt
from netCDF4 import getlibversion

platinstdict = {
    'ssmi-f10': ("SSMI (Special Sensor Microwave Imager)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmi-f11': ("SSMI (Special Sensor Microwave Imager)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmi-f13': ("SSMI (Special Sensor Microwave Imager)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmi-f14': ("SSMI (Special Sensor Microwave Imager)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmi-f15': ("SSMI (Special Sensor Microwave Imager)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmis-f16': ("SSMI/S (Special Sensor Microwave Imager/Sounder)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmis-f17': ("SSMI/S (Special Sensor Microwave Imager/Sounder)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'ssmis-f18': ("SSMI/S (Special Sensor Microwave Imager/Sounder)",
                  "DMSP (Defense Meteorological Satellite Program)"),
    'amsr-aq': ("Advanced Microwave Scanning Radiometer for EOS (AMSR-E)",
                "JAXA Advanced Earth Observation Satellite-II 'Midori II' (ADEOS-II)"),
    'amsr2-gw1': ("AMSR2 (Advanced Microwave Scanning Radiometer 2)",
                  "JAXA Global Change Observation Mission 1st Water"),
    'multi-oi': ("Multi-sensor analysis",
                 "Multi-sensor analysis"),
    'wind': ("ERA5 Wind Data",
             "ERA5 Wind Data")
}


def global_data(hemi, platinst):
    '''Compile the global variables for the NetCDF file'''

    gcmdk_area = 'ANA_GBLATTRVAL_GCMDKEYWORDS_REGION_' + hemi.upper()
    areaname = 'ANA_GBLATTRVAL_REGION_' + hemi.upper()

    try:
        osiroot = os.environ['OSI_ROOT']
    except KeyError:
        print("Please export OSI_ROOT path to OSISAF software OSI_HL_Ice directory")
    absdir = os.path.dirname(os.path.abspath(__file__))

    # Reading the C headers with the global information into a lookup dictionary
    gdict = {}
    head_list = [os.path.join(osiroot,
                 'analysis_nc/src/analysis_product_conventions.h'),
                 os.path.join(osiroot, '../osisaf_hl_versions.h')]
    read_defs_from_c_headers(gdict, head_list)

    # Forming a dictionary of the required global information
    glob_data = {
        'title': gdict['ICEANA_LRSID_GBLATTRVAL_TITLE'],
        'product_name': gdict['ICEANA_LRSID_GBLATTRVAL_NAME'],
        'product_status': gdict['ICEANA_LRSID_GBLATTRVAL_STATUS'],
        'abstract': open(os.path.join(absdir, 'lrsid_cdr_abstract.txt'),
                         'r').read(),
        'topiccategory': gdict['ICEANA_GBLATTRVAL_TOPICS'],
        'keywords': (gdict['ICEANA_LRSID_GBLATTRVAL_KEYWORDS'] + ',' +
                     gdict['ICEANA_GBLATTRVAL_KEYWORDS']),
        'gcmd_keywords': (gdict['ICEANA_LRSID_GBLATTRVAL_GCMDKEYWORDS_PROD']
                          + '\n' +
                          gdict['ICEANA_GBLATTRVAL_GCMDKEYWORDS_VLOC']
                          + '\n' +
                          gdict['ICEANA_GBLATTRVAL_GCMDKEYWORDS_INST']),

        'activity_type': gdict['ANA_GBLATTRVAL_ACTIVITY'].strip(),
        'area': gdict[areaname],
        'instrument_type': platinstdict[platinst][0],
        'platform_name': platinstdict[platinst][1],
        'project_name': gdict['ICEANA_GBLATTRVAL_PROJNAME'],
        'institution': gdict['ICEANA_GBLATTRVAL_INST'],
        'PI_name': gdict['ICEANA_LRSID_GBLATTRVAL_PI'],
        'contact': gdict['ICEANA_GBLATTRVAL_CONTACT'],
        'distribution_statement': gdict['ICEANA_GBLATTRVAL_CONTACT'],
        'copyright_statement': 'Copyright {} EUMETSAT'.format(dt.date.today().year),
        'references': ("Product User Manual for OSI-405-c, Lavergne, v1.8, 2016\n" + "Validation Report for OSI-405-c, Lavergne, v5, 2016\n" + "Algorithm Theoretical Basis Document for OSI-405-c, Lavergne, v1.3, 2016\n" +
                       gdict['ICEANA_GBLATTRVAL_WEB_CENTRAL']),
        'history': "{:04d}-{:02d}-{:02d} creation".format(dt.date.today().year,
                                                          dt.date.today().month,
                                                          dt.date.today().day),
        'product_version': gdict['ICEANA_LRSID_GBLATTRVAL_VER_PROD'],
        'software_version': gdict['ICEANA_GBLATTRVAL_VER_SOFT'],
        'netcdf_version': getlibversion().split()[0],
        'Conventions': gdict['ANA_GBLATTRVAL_CONV']
    }

    return glob_data


def read_defs_from_c_headers(gdict, head_list):
    '''Read definitions from a C header file and append to a Python
    dictionary'''

    # Reading all header files into a list
    head_data = []
    for head_file in head_list:
        with open(head_file, 'r') as headf:
            head_data.extend(headf.readlines())

    # Reading all definition lines into a dictionary
    def_expr = re.compile("^#define\s+([A-Z0-9_]+)\s+(.+)$")
    comment_val = re.compile("^(.+)\s+/\*.+\*/\s*$")
    for line in head_data:
        def_match = re.match(def_expr, line)
        if def_match:
            key = def_match.group(1)
            value = def_match.group(2)
            # Stripping off any comments
            com_match = re.match(comment_val, value)
            if com_match:
                value = com_match.group(1)
            # Stripping whitespace
            value = value.strip()
            # Stripping off quotes
            value = value.rstrip('"').lstrip('"')
            # Setting in the dictionary
            gdict[key] = value

    # Consolidating the global dictionary
    cons_flg = True
    while cons_flg:
        for key, value in gdict.items():
            cons_flg = False
            if value in gdict.keys():
                gdict[key] = gdict[value]
                cons_flg = True
