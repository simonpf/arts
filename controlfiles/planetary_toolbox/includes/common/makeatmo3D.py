################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file interpolates raw 1D basic atmospheric data (z_field_raw,           #
# t_field_raw, vmr_field_raw) to the calculation grids (p_grid) for 3D         #
# atmosphere output.                                                           #
#                                                                              #
# This file expects the following input parameters:                            #
#   p_grid             as the WSV                                              #
#   lat_grid           as the WSV                                              #
#   lon_grid           as the WSV                                              #
#   z_field_raw        as the WSV                                              #
#   t_field_raw        as the WSV                                              #
#   vmr_field_raw      as the WSV                                              #
#   interp_order       (Index)         Grid interpolation order                #
#   vmr_zeropad        (Index)         Flag, whether to fill VMR at            #
#                                       non-covered profile regions with zeros #
# Output:                                                                      #
#   z_field            as the WSV                                              #
#   t_field            as the WSV                                              #
#   vmr_field          as the WSV                                              #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# here, we do LTE only so far, but need to initialize NLTE t-field accordingly
ws.Touch(ws.nlte_field_raw)
# this brings z_, t_, and vmr_field_raw to p/lat/lon_grid for 2D or 3D atmosphere
ws.AtmFieldsCalcExpand1D(interp_order=ws.interp_order, vmr_zeropadding=ws.vmr_zeropad)
