# Processing software for the EUMETSAT OSI SAF Sea Ice Drift Climate Data Record v1 (OSI-455)

This repository holds the processing software used at the Norwegian Meteorological Institute in 2020-2022 to prepare the Global Sea Ice Drift Climate Data Record v1 of the EUMETSAT OSI SAF (OSI-455), based on Passive Microwave Radiometer data (SSM/I, SSMIS, AMSR-E, and AMSR2) and ERA5 winds.

The sea-ice drift CDR and its documentation are accessible as:
https://osi-saf.eumetsat.int/products/osi-455 

The manuscript introducing the sea-ice drift CDR is at Earth System Science Data (ESSD):

Lavergne, T. and Down, E.: A Climate Data Record of Year-Round Global Sea Ice Drift from the EUMETSAT OSI SAF, Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2023-40, in review, 2023.

## Other related resources:
* The software to reproduce the figures and tables in the ESSD paper (https://doi.org/10.5281/zenodo.7595043).
* Example software (notebook) to transform and rotate the CDR vectors into other projections (https://doi.org/10.5281/zenodo.8315155).
* The DOI landing page at EUMETSAT (http://dx.doi.org/10.15770/EUM_SAF_OSI_0012).

## Caveats and limitations:

The processing software in this repository is to help users of the data record to gain some insights about how the algorithms were structured and implemented. This supports _transparent_ and _open_ science.

We evertheless acknowledge that the software will be cumbersome for anyone but us to install and run. The software in this repository is not enough for _reproducible_ science.
