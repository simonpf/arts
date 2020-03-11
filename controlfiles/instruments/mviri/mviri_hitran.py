# Included by mviri_fast.arts and mviri_reference.arts

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Read HITRAN catalog:
# abs_linesReadFromHitran2004( abs_lines,
#                             hitran_file,
#                             2.0022234e+12,
#                             2.01e+12 )
ws.abs_linesReadFromHitran(
    ws.abs_lines, ws.hitran_file, 100000000000.0, 100000000000000.0
)
ws.abs_lines_per_speciesCreateFromLines()
# write out for HITRAN-independent use in test case runs
# WriteXML( "zascii", abs_lines, "../../testdata/abs_lines_IR.xml.gz" )
# WARNING: If you redefine abs_species, and want to do a line-by-line
# calculation, you also have to call
# abs_lines_per_speciesCreateFromLines again.
