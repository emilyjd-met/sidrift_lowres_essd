import os
from netCDF4 import Dataset
import numpy as np


here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, '../../data')
outdir = os.path.join(here, '../../output')

paramf = {'nh': os.path.join(datadir, 'paramfiles/inv_params_osi455_nh_200301-202012_1day.nc'),
          'sh': os.path.join(datadir, 'paramfiles/inv_params_osi455_sh_200208-202007_1day.nc')}

outfile = os.path.join(outdir, 'latex_stats_paramfiles.tex')

mondict = {0: 'Jan',
           1: 'Feb',
           2: 'Mar',
           3: 'Apr',
           4: 'May',
           5: 'Jun',
           6: 'Jul',
           7: 'Aug',
           8: 'Sep',
           9: 'Oct',
           10: 'Nov',
           11: 'Dec'}


def main():

    varlist = ['flags', 'absA_gapfill', 'ThetaA_gapfill', 'RMS_Res_real',
               'RMS_Res_imag', 'C_real_gapfill', 'C_imag_gapfill', 'lat', 'lon']
    stats = {}
    data = {}

    for hemi in ['nh', 'sh']:

        stats[hemi] = {}
        data[hemi] = {}

        with Dataset(paramf[hemi], 'r') as dataset:

            for var in varlist:
                data[hemi][var] = dataset.variables[var][:]

        # Finding data shape
        shx, shy = data['nh']['RMS_Res_real'][0, :, :].shape
        rmsdtype = data['nh']['RMS_Res_real'].dtype

        # Finding the per-month flags
        flagok = np.zeros((12, shx, shy), dtype=data[hemi]['flags'][:].dtype)
        for mon in range(12):
            flagok[mon, :, :] = data[hemi]['flags'][mon, :, :] == 0

        # Calculating the geostrophic speeds
        cdtype = np.double
        data[hemi]['C_gapfill'] = np.zeros((12, shx, shy), dtype=cdtype)

        # In addition to the other flagging, we want to flag any NaNs
        # in the geostropic currents
        flagokgeo = flagok.copy()
        for mon in range(12):
            flagokgeo[mon, :, :] = np.logical_and(flagok[mon, :, :], ~data[hemi]['C_real_gapfill'][mon, :, :].mask)
            flagokgeo[mon, :, :] = np.logical_and(flagok[mon, :, :], ~data[hemi]['C_imag_gapfill'][mon, :, :].mask)
        data[hemi]['C_real_gapfill'].mask = 0
        data[hemi]['C_imag_gapfill'].mask = 0

        for mon in range(12):
            # Set the gapfilled and sea/land areas to 0 to avoid nasty numbers
            data[hemi]['C_real_gapfill'][mon, :, :][flagokgeo[mon, :, :] == 0] = 0
            data[hemi]['C_imag_gapfill'][mon, :, :][flagokgeo[mon, :, :] == 0] = 0

            ccomb = np.zeros((2, shx, shy), dtype=cdtype)
            cc1 = data[hemi]['C_real_gapfill'][mon, :, :]
            ccomb[0, :, :] = cc1 * cc1
            cc2 =  data[hemi]['C_imag_gapfill'][mon, :, :]
            ccomb[1, :, :] = cc2 * cc2
            ccombsum = np.nansum(ccomb, axis=0)
            data[hemi]['C_gapfill'][mon, :, :] = np.sqrt(ccombsum)

        # Averaging the RMS var
        data[hemi]['RMS_Res'] = np.zeros((12, shx, shy), dtype=rmsdtype)
        for mon in range(12):
            rmsav = np.zeros((2, shx, shy), dtype=rmsdtype)
            rmsav[0, :, :] = data[hemi]['RMS_Res_real'][mon, :, :]
            rmsav[1, :, :] = data[hemi]['RMS_Res_real'][mon, :, :]
            data[hemi]['RMS_Res'][mon, :, :] = np.nanmean(rmsav, axis=0)

        for var in ['absA_gapfill', 'ThetaA_gapfill', 'RMS_Res', 'C_gapfill']:
            stats[hemi][var] = {}

            for mon in range(12):
                # Calculate the per-month stats where the flag shows
                # no gapfilling
                if var == 'C_gapfill':
                    field = data[hemi][var][mon, :, :][flagokgeo[mon, :, :] == 1]
                else:
                    field = data[hemi][var][mon, :, :][flagok[mon, :, :] == 1]
                stats[hemi][var][mon] = np.nanmean(field)
                # Changing the absA variables to percentage
                if var == 'absA_gapfill':
                    stats[hemi][var][mon] = 100. * stats[hemi][var][mon]

    # Write out the latex file
    with open(outfile, 'a+') as outf:

        outf.write("\\begin{table}[htb]\n")
        outf.write("\\begin{center}\n")
        outf.write("\\begin{tabular}{|l||c|c|c||c|c|c|}\n")
        outf.write("\hline\n")
        outf.write("{} & \multicolumn{3}{|c||}{Northern hemisphere} & \multicolumn{3}{|c|}{Southern hemisphere} \\\\ \n")
        outf.write("\hline\n")
        outf.write("Month & $\\langle|A|\\rangle$ & $\\langle\\theta \\rangle$ & $\\langle U_{wg} \\rangle$ & $\\langle|A|\\rangle$ & $\\langle\\theta\\rangle$ & $\\langle U_{wg} \\rangle$ \\\\ \n")
        outf.write("{} & / \\% & / degree & / ms$^{-1}$ & / \\% & / degree & / ms$^{-1}$ \\\\ \n")
        outf.write("\hline\n")

        for mon in range(12):
            outf.write("{} & {} & {} & {} & {} & {} & {} \\\\ \n".format(
                mondict[mon],
                "{:+.1f}".format(stats['nh']['absA_gapfill'][mon]),
                "{:+.1f}".format(stats['nh']['ThetaA_gapfill'][mon]),
                "{:+.3f}".format(stats['nh']['C_gapfill'][mon]),
                "{:+.1f}".format(stats['sh']['absA_gapfill'][mon]),
                "{:+.1f}".format(stats['sh']['ThetaA_gapfill'][mon]),
                "{:+.3f}".format(stats['sh']['C_gapfill'][mon])))
        outf.write("\hline\n")
        outf.write("\end{tabular}\n")
        outf.write("\end{center}\n")
        outf.write("\caption[Parameter statistics]{}\n")
        outf.write("\label{tab:stats:param}\n")
        outf.write("\end{table}\n")
        outf.write("\n\n")

#    # OLD VERSION WITH RMS ERRORS
#    # Write out the latex file
#    with open(outfile, 'a+') as outf:
#
#        outf.write("\\begin{table}[htb]\n")
#        outf.write("\\begin{center}\n")
#        outf.write("\\begin{tabular}{|l||c|c|c||c|c|c|}\n")
#        outf.write("\hline\n")
#        outf.write("{} & \multicolumn{3}{|c||}{Northern hemisphere} & \multicolumn{3}{|c|}{Southern hemisphere} \\\\ \n")
#        outf.write("\hline\n")
#        outf.write("Month & $\\langle|A|\\rangle$ & $\\langle\\theta \\rangle$ & $\\langle \\varepsilon(dX) + \\varepsilon(dY) \\rangle$ & $\\langle|A|\\rangle$ & $\\langle\\theta\\rangle$ & $\\langle \\varepsilon(dX) + \\varepsilon(dY) \\rangle$ \\\\ \n")
#        outf.write("{} & / \\% & / degree & / ms$^{-1}$ & / \\% & / degree & / ms$^{-1}$ \\\\ \n")
#        outf.write("\hline\n")
#
#        for mon in range(12):
#            outf.write("{} & {} & {} & {} & {} & {} & {} \\\\ \n".format(
#                mondict[mon],
#                "{:+.1f}".format(stats['nh']['absA_gapfill'][mon]),
#                "{:+.1f}".format(stats['nh']['ThetaA_gapfill'][mon]),
#                "{:+.3f}".format(stats['nh']['RMS_Res'][mon]),
#                "{:+.1f}".format(stats['sh']['absA_gapfill'][mon]),
#                "{:+.1f}".format(stats['sh']['ThetaA_gapfill'][mon]),
#                "{:+.3f}".format(stats['sh']['RMS_Res'][mon])))
#        outf.write("\hline\n")
#        outf.write("\end{tabular}\n")
#        outf.write("\end{center}\n")
#        outf.write("\caption[Parameter statistics]{}\n")
#        outf.write("\label{tab:stats:param}\n")
#        outf.write("\end{table}\n")
#        outf.write("\n\n")


if __name__ == '__main__':

    main()
