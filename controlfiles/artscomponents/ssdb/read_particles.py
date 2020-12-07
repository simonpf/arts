from pathlib import Path
import pyarts
from pyarts.workspace import Workspace
SSDB_PATH = Path.home() / ("Dendrite/SSDB/SSDB_EUMETSAT/SSD_ARO/"
                           "AzimuthallyRandom/Ice/SingleCrystals/"
                           "Pristine/PlateType1_Id9")
PARTICLE_PATH = SSDB_PATH / "AzimuthallyRandom_beta010.0deg"

#
# Read single scattering particle.
#

ws = Workspace()
ws.ScatteringParticleCreate("scattering_particle")
ws.scattering_particleReadFromARTSSDB(ws.scattering_particle, PARTICLE_PATH)
ws.scattering_particle.print()

#
# Read scattering habit.
#

ws.ArrayOfScatteringParticleCreate("scattering_particles")
ws.scattering_particleReadFromARTSSDB(ws.scattering_particles, SSDB_PATH)
ws.scattering_particles.print()
