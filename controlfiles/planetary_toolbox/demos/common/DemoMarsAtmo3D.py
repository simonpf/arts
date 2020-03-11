################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for deriving Martian (atmospheric) data from the     #
# arts-xml-data package and convert it to the common spatial grids             #
# (p/lat/lon_grid), such that they can be applied in radiative transfer        #
# calculations.                                                                #
# It is for a 3D atmosphere (for 1D use DemoMarsAtmo1D.arts instead).          #
#                                                                              #
# It provides following output:                                                #
#   atmosphere_dim    as the WSV                                               #
#   p_grid            as the WSV                                               #
#   lat_grid          as the WSV                                               #
#   lon_grid          as the WSV                                               #
#   z_field           as the WSV                                               #
#   t_field           as the WSV                                               #
#   vmr_field         as the WSV                                               #
#   wind_u/v/w_field  as the WSV                                               #
#   abs_species       as the WSV                                               #
#                                                                              #
# The user is supposed to select (atmospheric case, species to include) from   #
# lists. Details of setting rules are given at the place of the settings.      #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/common/createvars.arts                                            #
#                                                                              #
# This template makes use of the following include files                       #
#   includes/mars/atmo_mars.arts                                               #
#   includes/mars/getatmo_mars.arts                                            #
#   includes/common/getgrids_3Dsimple.arts                                     #
#   includes/common/makeatmo3D.arts                                            #
#   includes/mars/getwind_mars.arts                                            #
#   includes/common/makefield3D.arts                                           #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
# set up name arrays and the like for selections
ws.execute_controlfile("planetary_toolbox/includes/mars/atmo_mars.arts")
# NOT modify
# prepare the variables for the atmosphere case & species selections
ws.IndexCreate("Ls")
ws.IndexCreate("daytime")
ws.IndexCreate("dust")
ws.IndexCreate("solar")
ws.ArrayOfIndexCreate("basespecies")
ws.ArrayOfIndexCreate("h2ospecies")
ws.ArrayOfIndexCreate("ch4species")
ws.ArrayOfIndexCreate("Necase")
ws.ArrayOfIndexCreate("vertwind")
ws.ArrayOfIndexCreate("NSwind")
ws.ArrayOfIndexCreate("EWwind")
ws.NumericCreate("latmin")
ws.NumericCreate("latmax")
ws.NumericCreate("lonmin")
ws.NumericCreate("lonmax")
ws.IndexCreate("nlat")
ws.IndexCreate("nlon")
ws.IndexCreate("auxfield_zeropad")
ws.IndexCreate("vmr_zeropad")
ws.IndexCreate("interp_order")
# do NOT modify
# set atmospheric dimensionality to 1D (for 3D use another template!)
ws.AtmosphereSet3D()
# only MODIFY if you know, what you are doing (else the default setting should
#  be fine).
#
# interpolation order for atmospheric data
ws.IndexSet(ws.interp_order, 1)
# assume species-vmr to be zero at pressures/altitudes that are not covered by
#  the species' profile data?
ws.IndexSet(ws.vmr_zeropad, 1)
# as above, but for other data (namely: wind)
ws.IndexSet(ws.auxfield_zeropad, 1)
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# Define limits for vertical grid (in terms of pressure)
# ---
# The grid itself is taken from the data (z_field).
# Setting limits to very low and high values will preserve the data grid limit
#  at the respective end.
# For Mars, surface pressure is about 8e2, pressure @ 100km around 1e-2, and
#  pressure @ 300km (highest altitude on Mars, where we have data) around 1e-8.
ws.NumericSet(ws.pmin, 0.01)
ws.NumericSet(ws.pmax, 1e99)
# Define limits for lat/lon grids
# ---
# allowed for latitudes: [-90,90]
# allowed for longitudes:
#    values in [-360,360]
#    lon_max - lon_min <= 360
# Note: In case of global coverage (lon_max - lon_min == 360), cyclicity rules
# are applied (e.g., the measurement is allowed to switch over the 360-0 border
# meridian; equality of x and 360+x is recognised).
ws.NumericSet(ws.latmin, -90.0)
ws.NumericSet(ws.latmax, 90.0)
ws.NumericSet(ws.lonmin, -180.0)
ws.NumericSet(ws.lonmax, 180.0)
# Define number of points in lat/lon grids
# ---
# Minimum: nlat = nlon = 2
# Minimum is perfectly sufficient for homogeneous fields, but problematic when
#  fields have some structure.
# NOTE: The lat/lon grids with this spacing will also be applied for the
#  surface. That is, use more than the minimum spacing if the surface fields
#  are inhomogeneous (which is true for surface altitude, temperature, and
#  refractive index for Mars; finest spacing there is 5deg, which corresponds to
#  nlat = 37, nlon = 73 for a global grid).
ws.IndexSet(ws.nlat, 37)
ws.IndexSet(ws.nlon, 73)
# ===========================================
# Select the atmospheric scenario to be used
# ---
# Mars atmospheric scenarios are available for 4 different seasons, 2 daytimes,
# 3 dust loadings, and 3 solar activities. You chose all of them separately, and
# the resulting full scenario (name, location) is derived from these components.
# Season: (northern) spring (Ls0), summer (Ls90), fall (Ls180), winter (Ls270)
#                       0        ,    1         ,   2         ,    3
ws.IndexSet(ws.Ls, 2)
# Daytime: day, night
#           0 ,   1
ws.IndexSet(ws.daytime, 1)
# Dust Loading: low, medium, high
#                0 ,   1   ,  2
ws.IndexSet(ws.dust, 0)
# Solar activity: low, average, high
#                  0 ,    1   ,  2
ws.IndexSet(ws.solar, 1)
# ===========================================
# Select the trace gases (and possible sub-scenarios) to be used
# ---
# Basic species
# ---
# refers to species with only one version here. no sub-options/further
#  specifications required
# NOTE: if you select CO2-CO2 CIA here, you have to use an appropriate set for
#  the abs_xsec_agenda.
#
# Select ALL species you like to take into account.
### CO, CO2, CO2-CO2 CIA, CO2-CO2 PWR, H2, H2O2, H2S, HCl, N2, N2O, O, O2, ...
#    0,  1,   2,           3,           4,  5,    6,   7,   8,  9, 10, 11, ...
# ... O3, OCS, OH, SO2
# ... 12, 13,  14,  15
ws.ArrayOfIndexSet(ws.basespecies, [0, 1, 2, 7, 9, 12, 15])
# Species with sub-scenarios
# ---
# only set UP TO ONE for each species (else the species will be included
#  several times, which usually does not make sense).
# EMPTY selection de-selects the whole species.
### CH4: low, standard
#         0 ,    1
# select UP TO ONE
ws.ArrayOfIndexSet(ws.ch4species, [1])
# Electron density: depending on solar zenith angle.
# NOTE that SZA<90 are only available for day data, SZA>90 for night data
#  (non-matching selection will lead to a runtime error!). That is, COORDINATE
#  selection here with your selection of DAYTIME.
#       day       day       day       day       night
### SZA0-30, SZA30-50, SZA50-70, SZA70-90, SZA120-180
#         0,        1,        2,        3,          4
# select UP TO ONE
ws.ArrayOfIndexSet(ws.Necase, [4])
# Select species with separate isotopologue profiles available
# ---
# Here it is ok to select more than one entry. General species (aka 'all')
# selects all (remaining, i.e., not yet selected) abs lines of the species. That
# is, general species (i.e., highest index) shall be LAST in selection.
### H2O: HDO (162), all (remaining)
#          0      ,   1
# select as many AS YOU WANT, but 'all' index has to be in LAST position
ws.ArrayOfIndexSet(ws.h2ospecies, [1])
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# do NOT modify
# now, let the prepared include files do the actual work:
# (a) read in the raw atmosphere including all selected species
ws.execute_controlfile("planetary_toolbox/includes/mars/getatmo_mars.arts")
# (b) get the common grids for the atmosphere
ws.execute_controlfile("planetary_toolbox/includes/common/getgrids_3Dsimple.arts")
# (c) do the conversion from raw data with individual grids to the common grids
ws.execute_controlfile("planetary_toolbox/includes/common/makeatmo3D.arts")
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# Get non-abs_species data: wind
# ---
# Only select ONE element per wind component (else the latter will overwrite the
#  earlier).
### vertical wind: standard
#                     0
# select EXACTLY ONE.
# also, you need to UNCOMMENT lines further below if you want vertical wind.
ws.ArrayOfIndexSet(ws.vertwind, [0])
### N-S wind: standard
#                0
# select EXACTLY ONE.
# also, you need to UNCOMMENT lines further below if you want N-S wind.
ws.ArrayOfIndexSet(ws.NSwind, [0])
### E-W wind: standard
#                0
# select EXACTLY ONE.
# also, you need to UNCOMMENT lines further below if you want E-W wind.
ws.ArrayOfIndexSet(ws.EWwind, [0])
# if you WANT any wind (but only then), UNCOMMENT the next two lines
ws.execute_controlfile("planetary_toolbox/includes/mars/getwind_mars.arts")
ws.IndexSet(ws.abs_f_interp_order, 1)
# if you want VERTICAL winds, UNCOMMENT the following 3 lines
ws.Copy(ws.rawfield, ws.wind_w_raw)
ws.execute_controlfile("planetary_toolbox/includes/common/makefield3D.arts")
ws.Copy(ws.wind_w_field, ws.finalfield)
# if you want N-S winds, UNCOMMENT the following 3 lines
ws.Copy(ws.rawfield, ws.wind_v_raw)
ws.execute_controlfile("planetary_toolbox/includes/common/makefield3D.arts")
ws.Copy(ws.wind_v_field, ws.finalfield)
# if you want E-W winds, UNCOMMENT the following 3 lines
ws.Copy(ws.rawfield, ws.wind_u_raw)
ws.execute_controlfile("planetary_toolbox/includes/common/makefield3D.arts")
ws.Copy(ws.wind_u_field, ws.finalfield)
################################################################################
# END USER SETTINGS                                                            #
################################################################################
