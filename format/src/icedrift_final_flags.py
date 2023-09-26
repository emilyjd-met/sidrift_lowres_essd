'''Flags for the final format of the LR SID product and conversion between final and proc flags'''

#Flags for the operational LR sea ice drift product
id_final_flags = {

# REJECTION FLAG : UNPROCESSED
    'flag_value_missing': (0, 'missing_input_data',
        "missing satellite image data"),
    'flag_value_land': (1, 'over_land',
        "grid point is over land"),
    'flag_value_noice': (2, 'no_ice',
        "grid point is not over sufficient ice"),
    'flag_value_coast': (3, 'close_to_coast_or_edge',
        "grid point is too close to coast or edge"),
    'flag_value_summer': (4, 'summer_period',
        "unreliable vector was removed because in summer period"),
    'flag_value_masked_outside_icedrift': (5, 'masked_outside_icedrift',
        "ice drift from winds masked since value is outside region of icedrift from CMCC"),

# REJECTION FLAG : ERRONEOUS
    'flag_value_fail': (10, 'processing_failed',
        "optimization of the correlation (CMCC) failed"),
    'flag_value_lowmcc': (11, 'too_low_correlation',
        "vector removed because too low cross-correlation"),
    'flag_value_no_neighbour': (12, 'not_enough_neighbours',
        "vector removed because not enough neighbouring vectors"),
    'flag_value_filt_neighbour': (13, 'filtered_by_neighbours',
        "vector removed because motion is not consistent with neighbours"),

#QUALITY LEVEL FLAG : USE WITH CARE
    'flag_value_smallpatt': (20, 'smaller_pattern',
        "vector processed using a smaller matching window"),
    'flag_value_corr_neighbour': (21, 'corrected_by_neighbours',
        "vector corrected using the neighbouring vectors"),
    'flag_value_interp': (22, 'interpolated',
        "vector interpolated from the neighbours (multi-sensor product)"),
    'flag_value_windgapfill': (23, 'gapfilled_wind_parameter',
        "wind drift parameter file gap filled"),
    'flag_value_fillbywind': (24, 'filled_by_wind',
        "vector replaced by wind drift due to lack of satellite data"),
    'flag_value_blendwithwind': (25, 'blended_with_wind',
        "vector comprised of blended satellite and wind drift"),


#QUALITY LEVEL FLAG : NOMINAL
    'flag_value_nominal': (30, 'nominal_quality',
        "vector was retrieved with nominal algorithm"),
}

# Matrix to match the proc flags up with the final format flags
flag_matrix = {
    'icedrift_ok': 'flag_value_nominal',
    'icedrift_outside_imgborder': 'flag_value_missing',
    'icedrift_closeto_imgborder': 'flag_value_missing',
    'icedrift_center_over_land': 'flag_value_land',
    'icedrift_noice': 'flag_value_noice',
    'icedrift_closeto_coast_or_edge': 'flag_value_coast',
    'icedrift_closeto_missing': 'flag_value_missing',
    'icedrift_closeto_unprocessed': 'flag_value_missing',
    'icedrift_optimisation_fails': 'flag_value_fail',
    'icedrift_fails': 'flag_value_fail',
    'icedrift_lowcorrelation': 'flag_value_lowmcc',
    'icedrift_toolong': 'flag_value_missing', # This flag was missed in the original matrix
    'icedrift_refused_by_neighbours': 'flag_value_filt_neighbour',
    'icedrift_correct_by_neighbours': 'flag_value_corr_neighbour',
    'icedrift_noaverage': 'flag_value_no_neighbour',
    'icedrift_summerhold': 'flag_value_summer',
    'icedrift_oiinterp': 'flag_value_interp',
    'icedrift_smaller_pattern': 'flag_value_smallpatt',
    'icedrift_masked_nwp': 'flag_value_missing', # This flag was missed in the original matrix
    'icedrift_windgapfill': 'flag_value_windgapfill',
    'icedrift_maxdriftmask': 'flag_value_masked_outside_icedrift',
    'icedrift_fillbywind': 'flag_value_fillbywind',
    'icedrift_blendwithwind': 'flag_value_blendwithwind'
}
