################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for deriving Mars cloud data from the arts-xml-data  #
# package and convert it to the common spatial grids (p_grid), such that they  #
# can be applied in radiative transfer calculations. It is for a 1D atmosphere.#
#                                                                              #
# It provides following output:                                                #
#   pnd_field         as the WSV                                               #
#   scat_data         as the WSV                                               #
#   cloudbox_on       as the WSV                                               #
#   cloudbox_limits   as the WSV                                               #
#                                                                              #
# The user is supposed to select cloud types and scenarios from lists.         #
# Details of setting rules are given at the place of the settings.             #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
# Files to be included before this file:                                       #
#   demos/common/DemoMarsAtmo1D.arts      for atmospheric scenario selection   #
#                                                                              #
# It requires the following input:                                             #
#   f_grid            as the WSV                                               #
#                                                                              #
# This template makes use of the following include files                       #
#   includes/mars/clouds_mars.arts                                             #
#   includes/mars/getclouds_mars.arts                                          #
#   includes/common/makeclouds1D.arts                                          #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
# set up name arrays and the like for selections
ws.execute_controlfile("planetary_toolbox/includes/mars/clouds_mars.arts")
# NOT modify
# prepare the variables for the cloud case (pnd & ssd) selections
ws.IndexCreate("dustRI")
ws.IndexCreate("co2RI")
ws.IndexCreate("h2oRI")
ws.ArrayOfIndexCreate("dustcase")
ws.ArrayOfIndexCreate("co2case")
ws.ArrayOfIndexCreate("h2ocase")
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# ---
# Select the different cloud types to be included
# ---
# Select (up to) one case per cloud type - you are able to choose more (and case
#  will run), but this likely results in nonsense. However, it is ok to combine
#  different "cloud" types (dust, CO2 ice clouds, H2Oice clouds).
###
# Martian dust
# ---
# dust scenario
###
# particle size:     small      medium large very-large
# dust profile:  verywell-mixed   =      =       =
#                      0       ,  1   ,  2  ,    3  , ...
#
# particle size:    medium         =        =      =            =            =
# dust profile: verywell-mixed well-mixed mixed confined very-confined higly-confined
#                ... (1)      ,    4     ,  5  ,   6    ,       7     ,      8
#
# select (up to) one
ws.ArrayOfIndexSet(ws.dustcase, [0])
# dust refractive index
###
# minimum-absorption standard maximum-absorption
#         0         ,    1   ,        2
ws.IndexSet(ws.dustRI, 1)
# CO2 ice clouds
# ---
# cloud scenario
###
# cloud type: mesospheric-day polarnight-ch1      =            =
# profile:     gaussian       layer-at-20km  layer-at-8km tower-to-8km
#                 0          ,      1       ,     2      ,     3      , ...
#
# cloud type: polarnight-ch4       =
# profile:     layer-at-10km layer-at-4km
#             ...    4      ,      5
#
# select (up to) one
ws.ArrayOfIndexSet(ws.co2case, [3])
# CO2 ice refractive index
###
# minimum-absorption  best-estimate  maximum-absorption
#         0         ,      1       ,        2
ws.IndexSet(ws.co2RI, 1)
# H2O ice clouds
# ---
# cloud scenario
###
# cloud type:      type1             =              type2
# profile:    high-wide-layer low-narrow-layer high-wide-layer
#                    0       ,       1        ,       2
#
# select (up to) one
ws.ArrayOfIndexSet(ws.h2ocase, [2])
# H2O ice refractive index is fairly well known, hence no choice here
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# do NOT modify
ws.IndexSet(ws.h2oRI, 3)
# do NOT modify
# now, let the prepared include files do the actual work:
# (a) prepare the pnd field and single scattering data of all particle types/
#      cloud layers to include for digestion into respective ARTS WSVs
ws.execute_controlfile("planetary_toolbox/includes/mars/getclouds_mars.arts")
# (b) get the data into the actual WSVs and do conversion of pnd field from raw
#      data with individual grids to the common p/lat/lon_grids
ws.execute_controlfile("planetary_toolbox/includes/common/makeclouds1D.arts")
