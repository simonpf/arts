import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.xml import save

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/planet_earth.arts")

# Dimensionality of the atmosphere
ws.AtmosphereSet1D()

# Read data created by setup_test.m
ws.ReadXML(ws.p_grid, "test_data/p_grid.xml")
ws.ReadXML(ws.t_field, "test_data/t_field.xml")
ws.ReadXML(ws.z_field, "test_data/z_field.xml")
ws.ReadXML(ws.vmr_field, "test_data/vmr_field.xml")
ws.ReadXML(ws.particle_bulkprop_field, "test_data/particle_bulkprop_field")
ws.ReadXML(ws.particle_bulkprop_names, "test_data/particle_bulkprop_names")

ws.VectorSet(ws.f_grid, np.array([3.15e+10, 1.65e+11, 6.66e+11]))
ws.IndexSet(ws.stokes_dim, 2)
ws.VectorSet(ws.aa_grid, np.array([0.0]))
ws.VectorSet(ws.za_grid, np.array([0.0]))

#
# Define scattering species.
#

ws.ReadXML(ws.scat_data_raw, "test_data/scat_data.xml")
ws.ReadXML(ws.scat_meta, "test_data/scat_meta.xml")

## Rain
#@arts_agenda
#def psd_agenda_rain(ws):
#    ws.scattering_habitGetParticleSizes(x_unit="dveq")
#    ws.Copy(ws.psd_size_grid, ws.scat_species_x)
#    ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
#    ws.psdWangEtAl16(t_min=273, t_max=999.0)
#    ws.pndFromPsdBasic()
#
#ws.scattering_speciesAddScatteringHabit(name="rain",
#                                        scattering_data=ws.scat_data_raw.value[0],
#                                        scattering_meta_data=ws.scat_meta.value[0],
#                                        psd_agenda=psd_agenda_rain,
#                                        pnd_agenda_input=["RWC"])

# ICE
@arts_agenda
def psd_agenda_ice(ws):
    ws.scattering_habitGetParticleSizes(x_unit="dveq")
    ws.Copy(ws.psd_size_grid, ws.scat_species_x)
    ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
    ws.psdMcFarquaharHeymsfield97(t_min=10.0, t_max=273.0, t_min_psd=210.0)
    ws.pndFromPsdBasic()

ws.scattering_speciesAddScatteringHabit(name="ice",
                                       scattering_data=ws.scat_data_raw.value[1],
                                       scattering_meta_data=ws.scat_meta.value[1],
                                       psd_agenda=psd_agenda_ice,
                                       pnd_agenda_input=["IWC"])

#
# Calculate bulk properties.
#

ws.Tensor6Create("absorption_coeff")
ws.scattering_speciesCalcBulkAbsorptionCoeff(bulk_absorption_coeff=ws.absorption_coeff)
#ws.Tensor6Create("absorption_coeff_ref")
#ws.ReadXML(out=ws.absorption_coeff_)
ws.MatrixCreate("absorption_coeff_ref")
ws.ReadXML(out=ws.absorption_coeff_ref, filename="bulk_absorption_stokes_1.xml")


ac = ws.absorption_coeff.value[:, 0, 0, :, 0, 0]

ac_ref = ws.absorption_coeff_ref.value
