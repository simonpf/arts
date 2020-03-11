################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file interpolates raw 1D GriddedField3 data, e.g., wind component       #
# fields to the calculation grids (p_grid) for 3D atmosphere output.           #
#                                                                              #
# This file expects the following input parameters:                            #
#   p_grid             as the WSV                                              #
#   lat_grid           as the WSV                                              #
#   lon_grid           as the WSV                                              #
#   rawfield           (GriddedField3) a raw atmospheric field (note: the      #
#                                       respective field, e.g. wind_w_raw, has #
#                                       to be copied to rawfield BEFORE        #
#                                       including the file at hand.            #
#   interp_order       (Index)         Grid interpolation order                #
#   vmr_zeropad        (Index)         Flag, whether to fill VMR at            #
#                                       non-covered profile regions with zeros #
# Output:                                                                      #
#   finalfield         (Tensor3)       the atmospheric field interpolated to   #
#                                       p/lat/lon_grid, ready for doing RT     #
#                                       calculations with it.                  #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.GriddedFieldLatLonExpand(ws.rawfield, ws.rawfield)
ws.GriddedFieldLatLonRegrid(ws.rawfield, ws.lat_grid, ws.lon_grid, ws.rawfield, 1)
ws.GriddedFieldPRegrid(ws.rawfield, ws.p_grid, ws.rawfield, 1, ws.auxfield_zeropad)
ws.FieldFromGriddedField(
    ws.finalfield, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.rawfield
)
