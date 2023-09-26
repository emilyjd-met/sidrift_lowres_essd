'''Dictionary of uncertainties for NRT, CDR, etc
The uncertainties are specified differently for sensor, season and
status flag.
Later the uncertainties are mofified by an addition according to the
difference to the central time.

The polynomial fit coefficients for correcting to 12UTC are supplied
in a separate dictionary.
'''

ud = {}

# NRT NEAR-REAL TIME PRODUCT

ud['nrt'] = {}
ud['nrt']['nh'] = {}
ud['nrt']['nh']['amsr2-gw1-bt37'] = {}
ud['nrt']['nh']['amsr2-gw1-bt37']['winter'] = {0: 2.1,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['nh']['amsr2-gw1-bt37']['summer'] = {0: 7.5,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['nh']['amsr2-gw1-bt19'] = {}
ud['nrt']['nh']['amsr2-gw1-bt19']['winter'] = {0: 2.5,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['nh']['amsr2-gw1-bt19']['summer'] = {0: 5.0,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['nh']['amsr-aqua-bt37'] = {}
ud['nrt']['nh']['amsr-aqua-bt37']['winter'] = {0: 2.5,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['nh']['amsr-aqua-bt37']['summer'] = {0: 8.5,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['nh']['amsr-aqua-bt19'] = {}
ud['nrt']['nh']['amsr-aqua-bt19']['winter'] = {0: 3.0,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['nh']['amsr-aqua-bt19']['summer'] = {0: 6.0,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['nh']['ssmi'] = {}
ud['nrt']['nh']['ssmi']['winter'] = {0: 3.5,
                                     13: 9.0,
                                     17: 5.0}
ud['nrt']['nh']['ssmi']['summer'] = {0: 100.0,
                                     13: 100.0,
                                     17: 100.0}
ud['nrt']['nh']['ascat'] = {}
ud['nrt']['nh']['ascat']['winter'] = {0: 3.9,
                                      13: 9.0,
                                      17: 5.0}
ud['nrt']['nh']['ascat']['summer'] = {0: 100.0,
                                      13: 100.0,
                                      17: 100.0}
ud['nrt']['nh']['wind'] = {}
ud['nrt']['nh']['wind']['winter'] = {0: 3.0,
                                     19: 3.0}
ud['nrt']['nh']['wind']['summer'] = {0: 5.0,
                                     19: 5.0}

ud['nrt']['sh'] = {}
ud['nrt']['sh']['amsr2-gw1-bt37'] = {}
ud['nrt']['sh']['amsr2-gw1-bt37']['winter'] = {0: 2.7,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['sh']['amsr2-gw1-bt37']['summer'] = {0: 9.8,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['sh']['amsr2-gw1-bt19'] = {}
ud['nrt']['sh']['amsr2-gw1-bt19']['winter'] = {0: 3.0,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['sh']['amsr2-gw1-bt19']['summer'] = {0: 6.0,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['sh']['amsr-aqua-bt37'] = {}
ud['nrt']['sh']['amsr-aqua-bt37']['winter'] = {0: 3.25,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['sh']['amsr-aqua-bt37']['summer'] = {0: 11.05,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['sh']['amsr-aqua-bt19'] = {}
ud['nrt']['sh']['amsr-aqua-bt19']['winter'] = {0: 3.6,
                                               13: 9.0,
                                               17: 5.0}
ud['nrt']['sh']['amsr-aqua-bt19']['summer'] = {0: 7.2,
                                               13: 100.0,
                                               17: 100.0}
ud['nrt']['sh']['ssmi'] = {}
ud['nrt']['sh']['ssmi']['winter'] = {0: 5.3,
                                     13: 9.0,
                                     17: 5.0}
ud['nrt']['sh']['ssmi']['summer'] = {0: 100.0,
                                     13: 100.0,
                                     17: 100.0}
ud['nrt']['sh']['ascat'] = {}
ud['nrt']['sh']['ascat']['winter'] = {0: 3.1,
                                      13: 9.0,
                                      17: 5.0}
ud['nrt']['sh']['ascat']['summer'] = {0: 100.0,
                                      13: 100.0,
                                      17: 100.0}
ud['nrt']['sh']['wind'] = {}
ud['nrt']['sh']['wind']['winter'] = {0: 3.0,
                                     19: 3.0}
ud['nrt']['sh']['wind']['summer'] = {0: 5.0,
                                     19: 5.0}

# CDR CLIMATE DATA RECORD

#highval = 15.
highval = 10.
veryhighval = 100.

ud['cdr'] = {}

ud['cdr']['high'] = {}
ud['cdr']['high']['wind'] = {0: highval,
                             19: highval}
ud['cdr']['high']['sat'] = {0: highval,
                            13: highval,
                            17: highval}
ud['cdr']['veryhigh'] = {}
ud['cdr']['veryhigh']['wind'] = {0: veryhighval,
                                 19: veryhighval}
ud['cdr']['veryhigh']['sat'] = {0: veryhighval,
                                13: veryhighval,
                                17: veryhighval}

ud['cdr']['nh'] = {}
ud['cdr']['nh']['multi-oi'] = {0: 2.4,
                               16: 3.2}

ud['cdr']['nh']['amsr2-gw1-bt37'] = {}
ud['cdr']['nh']['amsr2-gw1-bt37']['winter'] = {0: 1.7, # 1.4
                                               13: 8.1, # 8.0
                                               17: 3.3} # 3.7
ud['cdr']['nh']['amsr-aq-bt37'] = {}
ud['cdr']['nh']['amsr-aq-bt37']['winter'] = {0: 1.7, # 1.4
                                               13: 8.1, # 8.0
                                               17: 3.3} # 3.7
ud['cdr']['nh']['ssmi'] = {}
ud['cdr']['nh']['ssmi']['winter'] = {0: 2.3,
                                     13: 8.0,
                                     17: 3.7}
ud['cdr']['nh']['wind'] = {}
ud['cdr']['nh']['wind']['summer'] = {0: 2.6,
                                     19: 1.8}
ud['cdr']['nh']['wind']['autumn'] = {0: 2.9,    # Only used for gapfilling,
                                     19: 2.9}   # use overall, not per-flag
ud['cdr']['nh']['wind']['winter'] = {0: 2.7,    # Only used for gapfilling,
                                     19: 2.7}   # use overall, not per-flag
ud['cdr']['nh']['wind']['spring'] = {0: 2.3,    # Only used for gapfilling,
                                     19: 2.3}   # use overall, not per-flag

ud['cdr']['sh'] = {}
ud['cdr']['sh']['multi-oi'] = {0: 3.8,
                               16: 6.9}

ud['cdr']['sh']['amsr2-gw1-bt37'] = {}
ud['cdr']['sh']['amsr2-gw1-bt37']['winter'] = {0: 2.8, # 2.9
                                               13: 8.3,
                                               17: 5.3}
ud['cdr']['sh']['amsr-aq-bt37'] = {}
ud['cdr']['sh']['amsr-aq-bt37']['winter'] = {0: 2.8, # 2.9
                                             13: 8.3,
                                             17: 5.3}
ud['cdr']['sh']['ssmi'] = {}
ud['cdr']['sh']['ssmi']['winter'] = {0: 3.6,
                                     13: 8.7,
                                     17: 6.2}
ud['cdr']['sh']['wind'] = {}
ud['cdr']['sh']['wind']['summer'] = {0: 3.0,
                                     19: 5.2}
ud['cdr']['sh']['wind']['autumn'] = {0: 3.0,    # Only used for gapfilling,
                                     19: 3.0}   # use overall, not per-flag
ud['cdr']['sh']['wind']['winter'] = {0: 3.0,    # Only used for gapfilling,
                                     19: 3.8}   # use overall, not per-flag
ud['cdr']['sh']['wind']['spring'] = {0: 3.0,    # Only used for gapfilling,
                                     19: 3.7}   # use overall, not per-flag



#ud['cdr'] = {}
#
#ud['cdr']['high'] = {}
#ud['cdr']['high']['wind'] = {0: 15.,
#                             19: 15.}
#ud['cdr']['high']['sat'] = {0: 15.,
#                            13: 15.,
#                            17: 15.}
#ud['cdr']['veryhigh'] = {}
#ud['cdr']['veryhigh']['wind'] = {0: 100.,
#                                 19: 100.}
#ud['cdr']['veryhigh']['sat'] = {0: 100.,
#                                13: 100.,
#                                17: 100.}
#
#ud['cdr']['nh'] = {}
#ud['cdr']['nh']['multi-oi'] = {16: 3.3}
#
#ud['cdr']['nh']['amsr2-gw1-bt37'] = {}
#ud['cdr']['nh']['amsr2-gw1-bt37']['winter'] = {0: 1.5,
#                                               13: 8.0,
#                                               17: 3.8}
#ud['cdr']['nh']['amsr2-gw1-bt37']['spring'] = {0: 1.9,
#                                               13: 7.0,
#                                               17: 2.7}
#ud['cdr']['nh']['amsr2-gw1-bt37']['summer'] = {0: 3.5,
#                                               13: 8.1,
#                                               17: 4.3}
#ud['cdr']['nh']['amsr2-gw1-bt37']['autumn'] = {0: 2.1,
#                                               13: 7.4,
#                                               17: 4.3}
#
#ud['cdr']['nh']['amsr2-gw1-bt19'] = {}
#ud['cdr']['nh']['amsr2-gw1-bt19']['winter'] = {0: 2.0,
#                                               13: 7.6,
#                                               17: 3.2}
#ud['cdr']['nh']['amsr2-gw1-bt19']['spring'] = {0: 2.3,
#                                               13: 8.5,
#                                               17: 3.4}
#ud['cdr']['nh']['amsr2-gw1-bt19']['summer'] = {0: 3.2,
#                                               13: 8.0,
#                                               17: 4.3}
#ud['cdr']['nh']['amsr2-gw1-bt19']['autumn'] = {0: 2.4,
#                                               13: 7.0,
#                                               17: 4.5}
#
#ud['cdr']['nh']['amsr-aq-bt37'] = {}
#ud['cdr']['nh']['amsr-aq-bt37']['winter'] = {0: 1.1,
#                                               13: 8.0, # Fill value
#                                               17: 1.7}
#ud['cdr']['nh']['amsr-aq-bt37']['spring'] = {0: 2.0,
#                                               13: 7.4,
#                                               17: 2.3}
#ud['cdr']['nh']['amsr-aq-bt37']['summer'] = {0: 3.7,
#                                               13: 8.5,
#                                               17: 4.5}
#ud['cdr']['nh']['amsr-aq-bt37']['autumn'] = {0: 3.3,
#                                               13: 8.2,
#                                               17: 5.0}
#
#ud['cdr']['nh']['amsr-aq-bt19'] = {}
#ud['cdr']['nh']['amsr-aq-bt19']['winter'] = {0: 2.6,
#                                               13: 8.6,
#                                               17: 2.6}
#ud['cdr']['nh']['amsr-aq-bt19']['spring'] = {0: 2.3,
#                                               13: 8.3,
#                                               17: 3.1}
#ud['cdr']['nh']['amsr-aq-bt19']['summer'] = {0: 3.6,
#                                               13: 8.5,
#                                               17: 4.6}
#ud['cdr']['nh']['amsr-aq-bt19']['autumn'] = {0: 3.4,
#                                               13: 9.3,
#                                               17: 5.5}
#
#ud['cdr']['nh']['ssmi'] = {}
#ud['cdr']['nh']['ssmi']['winter'] = {0: 2.3,
#                                     13: 8.0,
#                                     17: 3.7}
#ud['cdr']['nh']['ssmi']['spring'] = {0: 3.2,
#                                     13: 8.3,
#                                     17: 4.1}
#ud['cdr']['nh']['ssmi']['summer'] = {0: 5.4,
#                                     13: 10.7,
#                                     17: 6.4}
#ud['cdr']['nh']['ssmi']['autumn'] = {0: 3.5,
#                                     13: 8.0,
#                                     17: 4.9}
#
#ud['cdr']['nh']['wind'] = {}
#ud['cdr']['nh']['wind']['winter'] = {0: 3.1,
#                                     19: 3.0}
#ud['cdr']['nh']['wind']['spring'] = {0: 2.4,
#                                     19: 2.0}
#ud['cdr']['nh']['wind']['summer'] = {0: 2.6,
#                                     19: 2}
#ud['cdr']['nh']['wind']['autumn'] = {0: 3.3,
#                                     19: 3.5}
#
#ud['cdr']['sh'] = {}
#ud['cdr']['sh']['multi-oi'] = {16: 6.9}
#
#ud['cdr']['sh']['amsr2-gw1-bt37'] = {}
#ud['cdr']['sh']['amsr2-gw1-bt37']['winter'] = {0: 2.9,
#                                               13: 8.3,
#                                               17: 5.3}
#ud['cdr']['sh']['amsr2-gw1-bt37']['spring'] = {0: 2.7,
#                                               13: 15.0, # Fill value
#                                               17: 2.0}
#ud['cdr']['sh']['amsr2-gw1-bt37']['summer'] = {0: 3.3,
#                                               13: 7.0,
#                                               17: 5.4}
#ud['cdr']['sh']['amsr2-gw1-bt37']['autumn'] = {0: 3.3,
#                                               13: 7.7,
#                                               17: 5.9}
#
#ud['cdr']['sh']['amsr2-gw1-bt19'] = {}
#ud['cdr']['sh']['amsr2-gw1-bt19']['winter'] = {0: 3.1,
#                                               13: 8.8,
#                                               17: 6.3}
#ud['cdr']['sh']['amsr2-gw1-bt19']['spring'] = {0: 3.1,
#                                               13: 15.0, # Fill value
#                                               17: 2.7}
#ud['cdr']['sh']['amsr2-gw1-bt19']['summer'] = {0: 3.5,
#                                               13: 7.6,
#                                               17: 5.3}
#ud['cdr']['sh']['amsr2-gw1-bt19']['autumn'] = {0: 4.3,
#                                               13: 8.3,
#                                               17: 4.1}
#
#ud['cdr']['sh']['amsr-aq-bt37'] = {}
#ud['cdr']['sh']['amsr-aq-bt37']['winter'] = {0: 2.9, # Fill value
#                                             13: 8.3, # Fill value
#                                             17: 5.3} # Fill value
#ud['cdr']['sh']['amsr-aq-bt37']['spring'] = {0: 3.2,
#                                             13: 15.0, # Fill value
#                                             17: 15.0} # Fill value
#ud['cdr']['sh']['amsr-aq-bt37']['summer'] = {0: 3.2,
#                                             13: 3.8,
#                                             17: 15.0} # Fill value
#ud['cdr']['sh']['amsr-aq-bt37']['autumn'] = {0: 15.0, # Fill value
#                                             13: 15.0, # Fill value
#                                             17: 15.0} # Fill value
#
#ud['cdr']['sh']['amsr-aq-bt19'] = {}
#ud['cdr']['sh']['amsr-aq-bt19']['winter'] = {0: 3.1,
#                                             13: 15.0, # Fill value
#                                             17: 15.0} # Fill value
#ud['cdr']['sh']['amsr-aq-bt19']['spring'] = {0: 3.9,
#                                             13: 15.0, # Fill value
#                                             17: 15.0} # Fill value
#ud['cdr']['sh']['amsr-aq-bt19']['summer'] = {0: 3.7,
#                                             13: 15.0, # Fill value
#                                             17: 15.0} # Fill value
#ud['cdr']['sh']['amsr-aq-bt19']['autumn'] = {0: 15.0, # Fill value
#                                             13: 15.0, # Fill value
#                                             17: 15.0} # Fill value
#
#ud['cdr']['sh']['ssmi'] = {}
#ud['cdr']['sh']['ssmi']['winter'] = {0: 3.6,
#                                     13: 6.0,
#                                     17: 8.7}
#ud['cdr']['sh']['ssmi']['spring'] = {0: 4.1,
#                                     13: 7.8,
#                                     17: 4.4}
#ud['cdr']['sh']['ssmi']['summer'] = {0: 4.5,
#                                     13: 9.0,
#                                     17: 6.2}
#ud['cdr']['sh']['ssmi']['autumn'] = {0: 3.4,
#                                     13: 9.1,
#                                     17: 4.4}
#
#ud['cdr']['sh']['wind'] = {}
#ud['cdr']['sh']['wind']['winter'] = {0: 3.8,
#                                     19: 9.0}
#ud['cdr']['sh']['wind']['spring'] = {0: 3.6,
#                                     19: 6.0}
#ud['cdr']['sh']['wind']['summer'] = {0: 3.0,
#                                     19: 5.2}
#ud['cdr']['sh']['wind']['autumn'] = {0: 3.4,
#                                     19: 6.0}


# fitDT0 dictionary
fitDT0dict = {}
fitDT0dict['nrt'] = {}
fitDT0dict['nrt']['amsr'] = [0.01948381, -0.01426112, 0.]
fitDT0dict['nrt']['ssmi'] = [0.01422614, -0.00535288, 0.]
fitDT0dict['nrt']['ascat'] = [0.01029225, 0.0052347, 0.]
fitDT0dict['cdr'] = {}
fitDT0dict['cdr']['all'] = [0.015, -0.005, 0.]
