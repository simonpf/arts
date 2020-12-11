from pathlib import Path
import pyarts
from pyarts.workspace import Workspace
SSDB_PATH = Path.home() / ("Dendrite/SSDB/SSDB_EUMETSAT/SSD_ARO/"
                           "AzimuthallyRandom/Ice/SingleCrystals/"
                           "Pristine/PlateType1_Id9/"
                           "AzimuthallyRandom_beta060.0deg")
PARTICLE_PATH = (SSDB_PATH / "Dmax00013um_Dveq00010um_Mass4.79983e-13kg.nc")

#
# Read single scattering particle.
#

ws = Workspace()
ws.ScatteringParticleCreate("scattering_particle")
ws.scattering_particleReadFromARTSSSDB(ws.scattering_particle, PARTICLE_PATH)
ws.scattering_particle.print()
ws.WriteXML("binary", ws.scattering_particle, "particle.xml", 0)

#
# Read scattering habit.
#

ws.ArrayOfScatteringParticleCreate("scattering_particles")
ws.particle_habitReadFromARTSSSDB(ws.scattering_particles, SSDB_PATH)
ws.scattering_particles.print()
