# DEFINITIONS:  -*-sh-*-
#
############
# Venus specific settings
#
############
#
# Authors: Jana Mendrok
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
#
# Isotopologue ratios
#
ws.ReadXML(ws.isotopologue_ratios, "planets/Venus/isotopratio_Venus.xml")
#
# Reference ellipsoid (a spherical ellipsoid must be used for 1D)
#
ws.refellipsoidVenus(ws.refellipsoid, "Sphere")
#
# Weight of dry air [g/mol]
# (needed for hydrostatic equilibrium calculations)
# source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
#
ws.NumericSet(ws.molarmass_dry_air, 43.45)
#
# Gravity
# (needed for hydrostatic equilibrium calculations)
#
@arts_agenda
def g0_agenda(ws):
    ws.Ignore(ws.lon)
    ws.Ignore(ws.lat)
    ws.g0Venus()


ws.g0_agenda = g0_agenda

#
# Sidereal rotation period (−243.0185 Earth days)
#
ws.NumericSet(ws.planet_rotation_period, -20997000.0)
