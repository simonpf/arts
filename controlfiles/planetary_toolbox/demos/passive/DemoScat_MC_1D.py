################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for doing scattering calculations of passive         #
# measurements using the Monte Carlo scattering solver. It is for a 1D         #
# atmosphere.                                                                  #
# The file is supposed to be used as an include file in a full radiative       #
# transfer case. Atmospheric scenario, surface and sensor settings, etc. have  #
# to be done in the calling file before this file is included (hence,          #
# executed).                                                                   #
# This file requires the following input parameters:                           #
#   atmosphere_dim    as the WSV                                               #
#   f_grid            as the WSV                                               #
#   iy_unit           as the WSV                                               #
#   stokes_dim        as the WSV                                               #
#   p_grid            as the WSV                                               #
#   z_field           as the WSV                                               #
#   t_field           as the WSV                                               #
#   vmr_field         as the WSV                                               #
#   pnd_field         as the WSV                                               #
#   scat_data         as the WSV                                               #
#   cloudbox_on       as the WSV                                               #
#   cloudbox_limits   as the WSV                                               #
#   abs_species       as the WSV                                               #
#   z_surface         as the WSV                                               #
#   t_surface         as the WSV                                               #
#   rte_pos           as the WSV                                               #
#   allzang                     (Vector)  Sensor viewing angles                #
#   bad_partition_functions_ok  (Index)   Partition function handling flag     #
#                                                                              #
# It provides following OUTPUT (written to file):                              #
#   iy         as the WSV                                                      #
#               radiance; units selectable                                     #
#   iy_aux     as the WSV                                                      #
#               auxiliary output parameters (particularly of along-the-path    #
#               type), selectable variety                                      #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
# This template does not makes use of further include files.                   #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
ws.NumericCreate("ntmp")
ws.VectorCreate("vtmp")
ws.IndexCreate("nlat")
ws.IndexCreate("nlon")
################################################################################
# START USER SETTINGS - Modify settings according to you wishes                #
################################################################################
# ---
# Convergence limits for scattering solution
# ---
# Whichever is reached first stops the iteration (i.e., the calculation for more
#  "photons"). Setting a parameter to -1 means this paramter will be ignored in
#  the convergence decision.
# ---
# Statistical error of solution [W / (m2 sr Hz)]
# Note: dTb_planck=0.1K roughly correspond to 5e-22, 5e-20, 5e-18 W/(m2 sr Hz)
#  at 5, 50, and 500GHz, respectively.
ws.NumericSet(ws.mc_std_err, 5e-21)
# Maximum run time of solver (per monochromatic frequency)
ws.IndexSet(ws.mc_max_time, -1)
# Maximum number of iterations (i.e., "photons")
ws.IndexSet(ws.mc_max_iter, 250000)
# ---
# Define (auxiliary) data output
# ---
# Uncomment all parameters you want as auxiliary output (i.e., in addition to
#  total radiance/brigthness temperature). For meaning of each paramters see
#  online-doc of the WSM selected for iy_main_agenda (here: iyMC).
# NOTE: Last element NOT to be followed by comma.
# NOTE: Only use "Absorption, species X" up to the number of entries in
#  abs_species (clearsky calculations in Venus have at maximum 19 abs_species
#  entries, i.e. highest valid index is 18).
# NOTE: Only use "PND, type Y" up to the number of entries in scat_data.
#  Planet toolbox particle calculations have at maximum 10 particle entries,
#  i.e., highest valid index is 9).
# ---
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Error (uncorrelated)"])
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# only MODIFY if you know, what you are doing (else the default setting should
#  be fine).
#####
# setting agendas needed for RT calc (there are alternative settings, though)
#####
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# do NOT modify
#####
# use MC scattering module
#####
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__ScattMC)
# MC specific settings
#####
ws.MCSetSeedFromTime()
ws.mc_antennaSetPencilBeam()
# adaption due to MC requires 3D
###
# 1) set lat/lon components of platform position to dummies
ws.VectorSetConstant(ws.vtmp, 2, 0.0)
ws.Append(ws.rte_pos, ws.vtmp)
# 2a) blow up the clear-sky atmosphere
ws.AtmosphereSet3D()
ws.VectorLinSpace(ws.lat_grid, -30.0, 50.0, 5.0)
ws.VectorLinSpace(ws.lon_grid, -30.0, 30.0, 5.0)
ws.AtmFieldsExpand1D()
# 2b-1) adapt surface altitude to 3D
ws.nelemGet(ws.nlat, ws.lat_grid)
ws.nelemGet(ws.nlon, ws.lon_grid)
ws.VectorExtractFromMatrix(ws.vtmp, ws.z_surface, 0, "row")
ws.Extract(ws.ntmp, ws.vtmp, 0)
ws.MatrixSetConstant(ws.z_surface, ws.nlat, ws.nlon, ws.ntmp)
# 2b-2) adapt surface temperature to 3D
# note: It's necessary only for surface cases B-3x, no effects on others.
ws.VectorExtractFromMatrix(ws.vtmp, ws.t_surface, 0, "row")
ws.Extract(ws.ntmp, ws.vtmp, 0)
ws.MatrixSetConstant(ws.t_surface, ws.nlat, ws.nlon, ws.ntmp)
# 2c) blow up the pnd_field to 3D. but only in a reduced region of the
#  full 3D atmo we consider (cloudbox needs sufficient clear-sky around!)
ws.Extract(ws.itmp, ws.cloudbox_limits, 0)
ws.Extract(ws.pmax_cb, ws.p_grid, ws.itmp)
ws.Extract(ws.itmp, ws.cloudbox_limits, 1)
ws.Extract(ws.pmin_cb, ws.p_grid, ws.itmp)
ws.cloudboxSetManually(
    p1=ws.pmax_cb, p2=ws.pmin_cb, lat1=-10.0, lat2=30.0, lon1=-10.0, lon2=10.0
)
ws.pnd_fieldExpand1D()
# no Jacobians
#####
ws.jacobianOff()
# the checks necessary for full RT calc
#####
ws.atmfields_checkedCalc(bad_partition_functions_ok=ws.bad_partition_functions_ok)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc(scat_data=ws.scat_data_raw, scat_data_type="raw")
ws.propmat_clearsky_agenda_checkedCalc()
# and the RT calc
#####
ws.Copy(ws.strtmp, ws.iy_unit)
ws.StringSet(ws.iy_unit, "1")
ws.NumericCreate("za")
ws.AgendaCreate("forloop_agenda_angles")


@arts_agenda
def forloop_agenda_angles(ws):
    ws.Extract(ws.za, ws.allzang, ws.forloop_index)
    ws.rte_losSet(za=ws.za, aa=0.0)
    ws.Print(ws.rte_los, 0)
    ws.iyCalc()
    ws.iyApplyUnit(iy_unit=ws.strtmp)
    ws.WriteXMLIndexed(ws.output_file_format, ws.forloop_index, ws.iy, "", 0)
    ws.WriteXMLIndexed(ws.output_file_format, ws.forloop_index, ws.iy_aux, "", 0)


ws.forloop_agenda_angles = forloop_agenda_angles

ws.IndexCreate("nangles")
ws.nelemGet(ws.nangles, ws.allzang)
ws.IndexStepDown(ws.nangles, ws.nangles)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_angles)
ws.ForLoop(ws.forloop_agenda, 0, ws.nangles, 1)
