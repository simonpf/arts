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

#Rain
@arts_agenda
def psd_agenda_rain(ws):
   ws.scattering_habitGetParticleSizes(x_unit="dveq")
   ws.Copy(ws.psd_size_grid, ws.scat_species_x)
   ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
   ws.psdWangEtAl16(t_min=273, t_max=999.0)
   ws.pndFromPsdBasic()

ws.scattering_speciesAddScatteringHabit(name="rain",
                                       scattering_data=ws.scat_data_raw.value[0],
                                       scattering_meta_data=ws.scat_meta.value[0],
                                       psd_agenda=psd_agenda_rain,
                                       pnd_agenda_input=["RWC"])

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

# Absorption coefficient
ws.Tensor6Create("absorption_coeff")
ws.scattering_speciesCalcBulkAbsorptionCoeff(bulk_absorption_coeff=ws.absorption_coeff)
ws.MatrixCreate("absorption_coeff_ref")
ws.ReadXML(out=ws.absorption_coeff_ref, filename="bulk_absorption_stokes_1.xml")

absorption_coeff = ws.absorption_coeff.value[:, 0, 0, :, 0, 0]
absorption_coeff_ref = ws.absorption_coeff_ref.value
assert np.all(np.isclose(absorption_coeff, absorption_coeff_ref.T))

# Extinction coefficient
ws.Tensor6Create("extinction_coeff")
ws.scattering_speciesCalcBulkExtinctionCoeff(bulk_extinction_coeff=ws.extinction_coeff)
ws.MatrixCreate("extinction_coeff_ref")
ws.ReadXML(out=ws.absorption_coeff_ref, filename="bulk_extinction_stokes_1.xml")

extinction_coeff = ws.extinction_coeff.value[:, 0, 0, :, 0, 0]
extinction_coeff_ref = ws.absorption_coeff_ref.value
assert np.all(np.isclose(extinction_coeff, extinction_coeff_ref.T))

# Legendre coefficients

ws.Tensor7Create("legendre_coeffs")
ws.scattering_speciesCalcPhaseFunctionLegendreCoeffs(legendre_coeffs=ws.legendre_coeffs,
                                                     n_coeffs=9)
ws.Tensor3Create("legendre_coeffs_ref")
ws.ReadXML(out=ws.legendre_coeffs_ref, filename="legendre_coeffs.xml")

legendre_coeffs = ws.legendre_coeffs.value[:, 0, 0, :, 0, 0]
legendre_coeffs = 0.5  * (legendre_coeffs[1:] + legendre_coeffs[:-1])
legendre_coeffs = legendre_coeffs[::-1]
legendre_coeffs /= legendre_coeffs[:, :, [0]]
legendre_coeffs = np.nan_to_num(legendre_coeffs, nan=0.0)

legendre_coeffs_ref = np.transpose(ws.legendre_coeffs_ref.value, [1, 0, 2])
assert np.all(np.isclose(legendre_coeffs, legendre_coeffs_ref, atol=0.2))

# Absorption vector

ws.stokes_dim = 1
ws.Tensor7Create("absorption_vector")
ws.scattering_speciesCalcAbsorptionVector(bulk_absorption_vector=ws.absorption_vector)
absorption_vector = ws.absorption_vector.value[:, 0, 0, -1, 0, 0, 0]
absorption_vector = 0.5 * (absorption_vector[1:] + absorption_vector[:-1])
ws.Tensor5Create("absorption_vector_ref")
ws.ReadXML(out=ws.absorption_vector_ref, filename="absorption_vector_stokes_1.xml")
absorption_vector_ref = np.copy(ws.absorption_vector_ref.value[0, :, 0, 0, 0])
assert np.all(np.isclose(absorption_vector, absorption_vector_ref))

ws.stokes_dim = 2
ws.scattering_speciesCalcAbsorptionVector(bulk_absorption_vector=ws.absorption_vector)
absorption_vector = ws.absorption_vector.value[:, 0, 0, -1, 0, 0, :]
absorption_vector = 0.5 * (absorption_vector[1:] + absorption_vector[:-1])
ws.Tensor5Create("absorption_vector_ref")
ws.ReadXML(out=ws.absorption_vector_ref, filename="absorption_vector_stokes_2.xml")
absorption_vector_ref = ws.absorption_vector_ref.value[0, :, 0, 0, :]
assert np.all(np.isclose(absorption_vector, absorption_vector_ref))
# Extinction matrix

ws.stokes_dim = 1
ws.Tensor7Create("extinction_matrix")
ws.scattering_speciesCalcExtinctionMatrix(bulk_extinction_matrix=ws.extinction_matrix)
extinction_matrix = ws.extinction_matrix.value[:, 0, 0, -1, 0, 0, 0]
extinction_matrix = 0.5 * (extinction_matrix[1:] + extinction_matrix[:-1])
ws.Tensor6Create("extinction_matrix_ref")
ws.ReadXML(out=ws.extinction_matrix_ref, filename="extinction_matrices_stokes_1.xml")
extinction_matrix_ref = ws.extinction_matrix_ref.value[0, :, 0, 0, 0, 0]
assert np.all(np.isclose(extinction_matrix, extinction_matrix_ref))

ws.stokes_dim = 2
ws.scattering_speciesCalcExtinctionMatrix(bulk_extinction_matrix=ws.extinction_matrix)
extinction_matrix = ws.extinction_matrix.value[:, 0, 0, -1, 0, 0, :]
extinction_matrix = 0.5 * (extinction_matrix[1:] + extinction_matrix[:-1])
extinction_matrix = extinction_matrix.reshape(-1, 2, 2)
ws.ReadXML(out=ws.extinction_matrix_ref, filename="extinction_matrices_stokes_2.xml")
extinction_matrix_ref = ws.extinction_matrix_ref.value[0, :, 0, 0, :]
assert np.all(np.isclose(extinction_matrix, extinction_matrix_ref))

ws.Tensor6Create("extinction_matrix_ref")
ws.ReadXML(out=ws.extinction_matrix_ref, filename="extinction_matrices_stokes_1.xml")


# Scattering matrix

#
# For the tests of the scattering data we must use scattering data that
# has no temperature dependence, since this was how scattering matrices
# were computed in ARTS previously.
#

ws.ReadXML(ws.scat_data_raw, "test_data/scat_data_single_t.xml")
ws.ReadXML(ws.za_grid, "za_grid.xml")
ws.scattering_species = []
ws.scattering_speciesAddScatteringHabit(name="rain",
                                       scattering_data=ws.scat_data_raw.value[0],
                                       scattering_meta_data=ws.scat_meta.value[0],
                                       psd_agenda=psd_agenda_rain,
                                       pnd_agenda_input=["RWC"])
ws.scattering_speciesAddScatteringHabit(name="ice",
                                       scattering_data=ws.scat_data_raw.value[1],
                                       scattering_meta_data=ws.scat_meta.value[1],
                                       psd_agenda=psd_agenda_ice,
                                       pnd_agenda_input=["IWC"])

ws.scattering_speciesCalcExtinctionMatrix(bulk_extinction_matrix=ws.extinction_matrix)
n_mu = len(ws.za_grid.value) // 2
extinction_matrix = ws.extinction_matrix.value[:, 0, 0, -1, 0, n_mu-1::-1, 0]
ws.scattering_speciesCalcAbsorptionVector(bulk_absorption_vector=ws.absorption_vector)
absorption_vector = ws.absorption_vector.value[:, 0, 0, -1, 0, n_mu-1::-1, 0]
scattering_coeff = extinction_matrix - absorption_vector
scattering_coeff = np.expand_dims(scattering_coeff, [1, 3])


ws.stokes_dim = 1
ws.Tensor7Create("scattering_matrix")
ws.scattering_speciesCalcScatteringMatrix(bulk_scattering_matrix=ws.scattering_matrix,
                                          quadrature_type="lobatto",
                                          quadrature_degree=ws.za_grid.value.size)
#ws.scattering_speciesCalcScatteringMatrix(bulk_scattering_matrix=ws.scattering_matrix)
scattering_matrix = ws.scattering_matrix.value[:, -1, :, 0, :, :, 0]
n_mu = scattering_matrix.shape[1] // 2
scattering_matrix = scattering_matrix[:, n_mu-1::-1, n_mu-1::-1]
scattering_matrix = scattering_matrix * scattering_coeff
scattering_matrix = 0.5 * (scattering_matrix[1:] + scattering_matrix[:-1])

ws.Tensor6Create("scattering_matrix_ref")
ws.ReadXML(out=ws.scattering_matrix_ref, filename="scattering_matrices_stokes_1.xml")
scattering_matrix_ref = ws.scattering_matrix_ref.value[:, 0, :, 0, :,]
assert np.all(np.isclose(scattering_matrix, scattering_matrix_ref, rtol=1e-3))


ws.stokes_dim = 2
ws.Tensor7Create("scattering_matrix")
ws.scattering_speciesCalcScatteringMatrix(bulk_scattering_matrix=ws.scattering_matrix,
                                          quadrature_type="lobatto",
                                          quadrature_degree=ws.za_grid.value.size)
#ws.scattering_speciesCalcScatteringMatrix(bulk_scattering_matrix=ws.scattering_matrix)
scattering_matrix = ws.scattering_matrix.value[:, -1, :, 0, :, :, 0]
n_mu = scattering_matrix.shape[1] // 2
scattering_matrix = scattering_matrix[:, n_mu-1::-1, n_mu-1::-1]
scattering_matrix = 0.5 * (scattering_matrix[1:] + scattering_matrix[:-1])
scattering_matrix *= scattering_coeff.reshape(-1, 1, 1, 1)

ws.ReadXML(out=ws.scattering_matrix_ref, filename="scattering_matrices_stokes_2.xml")
scattering_matrix_ref = ws.scattering_matrix_ref.value[:, 0, :, 0, :,]
assert np.all(np.isclose(scattering_matrix, scattering_matrix_ref))
