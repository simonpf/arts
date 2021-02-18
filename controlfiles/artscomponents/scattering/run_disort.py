import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.xml import save

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")

# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)

# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)

# Blackbody surface
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.VectorSet(ws.surface_scalar_reflectivity, np.array([0.]))

# Standard ppath calculations
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)

# Radiative transfer agendas
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.iy_cloudbox_agenda, ws.iy_cloudbox_agenda__QuarticInterpField)

# Absorption species
ws.abs_speciesSet(species=['N2-SelfContStandardType', 'O2-PWR93', 'H2O-PWR98'])

# No line data needed here
ws.abs_lines_per_speciesSetEmpty()

# Dimensionality of the atmosphere
ws.AtmosphereSet1D()

# Brigtness temperatures used
ws.StringSet(ws.iy_unit, "PlanckBT")

# Various things not used
ws.ArrayOfStringSet(ws.iy_aux_vars, [])

ws.jacobianOff()
# Read data created by setup_test.m
ws.ReadXML(ws.p_grid, "test_data/p_grid.xml")
ws.ReadXML(ws.t_field, "test_data/t_field.xml")
ws.ReadXML(ws.z_field, "test_data/z_field.xml")
ws.ReadXML(ws.vmr_field, "test_data/vmr_field.xml")
ws.ReadXML(ws.particle_bulkprop_field, "test_data/particle_bulkprop_field")
ws.ReadXML(ws.particle_bulkprop_names, "test_data/particle_bulkprop_names")

#ws.ReadXML(ws.scat_data_raw, "test_data/scat_data.xml")
#ws.ReadXML(ws.scat_meta, "test_data/scat_meta.xml")

ws.ReadXML(ws.scat_data_raw, "test_data/scat_data.xml")
#ws.scat_data_raw = ws.scat_data_raw.value[1:]
ws.ReadXML(ws.scat_meta, "test_data/scat_meta.xml")
#ws.scat_meta = ws.scat_meta.value[1:]

ws.StringCreate("species_id_string")


# Scat species 0
ws.StringSet(ws.species_id_string, "RWC")
ws.ArrayOfStringSet(ws.pnd_agenda_input_names, ['RWC'])
@arts_agenda
def pnd_agenda_array(ws):
   ws.ScatSpeciesSizeMassInfo(species_index=ws.agenda_array_index, x_unit="dveq")
   ws.Copy(ws.psd_size_grid, ws.scat_species_x)
   ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
   ws.psdWangEtAl16(t_min=273.0, t_max=999.0)
   ws.pndFromPsdBasic()
ws.Append(ws.pnd_agenda_array, pnd_agenda_array)
ws.Append(ws.scat_species, ws.species_id_string)
ws.Append(ws.pnd_agenda_array_input_names, ws.pnd_agenda_input_names)

@arts_agenda
def psd_agenda(ws):
   ws.scattering_habitGetParticleSizes(x_unit="dveq")
   ws.Copy(ws.psd_size_grid, ws.scat_species_x)
   ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
   ws.psdWangEtAl16(t_min=273, t_max=999.0)
   ws.pndFromPsdBasic()

ws.scattering_speciesAddScatteringHabit(name="rain",
                                       scattering_data=ws.scat_data_raw.value[0],
                                       scattering_meta_data=ws.scat_meta.value[0],
                                       psd_agenda=psd_agenda,
                                       pnd_agenda_input=["RWC"])


# Scat species 1
ws.StringSet(ws.species_id_string, "IWC")
ws.ArrayOfStringSet(ws.pnd_agenda_input_names, ['IWC'])
@arts_agenda
def pnd_agenda_array(ws):
    ws.ScatSpeciesSizeMassInfo(species_index=ws.agenda_array_index, x_unit="dveq", x_fit_start=0.0001)
    ws.Copy(ws.psd_size_grid, ws.scat_species_x)
    ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
    ws.psdMcFarquaharHeymsfield97(t_min=10.0, t_max=273.0, t_min_psd=210.0)
    ws.pndFromPsdBasic()

ws.Append(ws.pnd_agenda_array, pnd_agenda_array)
ws.Append(ws.scat_species, ws.species_id_string)
ws.Append(ws.pnd_agenda_array_input_names, ws.pnd_agenda_input_names)

@arts_agenda
def psd_agenda(ws):
    ws.scattering_habitGetParticleSizes(x_unit="dveq")
    ws.Copy(ws.psd_size_grid, ws.scat_species_x)
    ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
    ws.psdMcFarquaharHeymsfield97(t_min=10.0, t_max=273.0, t_min_psd=210.0)
    ws.pndFromPsdBasic()

ws.scattering_speciesAddScatteringHabit(name="ice",
                                       scattering_data=ws.scat_data_raw.value[-1],
                                       scattering_meta_data=ws.scat_meta.value[-1],
                                       psd_agenda=psd_agenda,
                                       pnd_agenda_input=["IWC"])

# RT4 creates own grirs, so we need copies of the ones created above
ws.VectorCreate("za_grid_copy")
ws.VectorCreate("aa_grid_copy")

#
# Hybrid requires that ppath_lmax is not too high
ws.NumericSet(ws.ppath_lmax, 100.0)
ws.AgendaCreate("iy_hybrid_agenda")

@arts_agenda
def iy_hybrid_agenda(ws):
    ws.Ignore(ws.iy_id)
    ws.ppathCalc(cloudbox_on=0)
    ws.iyHybrid()
    # The line below is just temporary
    ws.Touch(ws.iy_aux)
ws.iy_hybrid_agenda = iy_hybrid_agenda

# Versions of y for various calculations
ws.VectorCreate("y_doit")
ws.VectorCreate("y_disort")
ws.VectorCreate("y_hybrid")
ws.VectorCreate("y_rt4")
#
ws.StringCreate("message")
# Perform some basic checks
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
# Intitial settings for tests
ws.IndexSet(ws.stokes_dim, 1)
# Scattering data tailored to these frequencies, so don't change!
ws.VectorSet(ws.f_grid, np.array([3.15e+10, 1.65e+11, 6.66e+11]))
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.MatrixSet(ws.sensor_pos, np.array([[20000.],
       [20000.],
       [10000.],
       [ 5000.]]))
ws.MatrixSet(ws.sensor_los, np.array([[180.],
       [130.],
       [160.],
       [ 20.]]))
# Some stuff that depends on the settings above
ws.sensorOff()
ws.atmgeom_checkedCalc()
ws.sensor_checkedCalc()
ws.scat_dataCalc()
# We need here a higher sca_mat_threshold. Otherwise there is an error for RWC
# and 668 GHz
ws.scat_data_checkedCalc(sca_mat_threshold=0.25)
#
ws.VectorExtractFromMatrix(ws.rtp_pos, ws.z_surface, 0, "row")
ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
# Test 1: Cloubox with no scattering
# ---------------------------------------------------------------------
ws.cloudboxSetFullAtm()
ws.pnd_fieldCalcFromParticleBulkProps()
ws.cloudbox_checkedCalc()
#
#
# Angular grids for DOIT and DISORT
ws.DOAngularGridsSet(N_za_grid=38, N_aa_grid=37)

ws.DisortCalc( pfct_method = "interpolate", nstreams=8)
ws.yCalc( y=ws.y_disort )
#save(ws.y_disort.value, "y_disort.xml", format="binary")


#
# Extract simplified scattering data
#

#ws.Copy( ws.za_grid_copy, ws.za_grid )
#ws.Copy( ws.aa_grid_copy, ws.aa_grid )
#ws.RT4Calc( nstreams = 16,
#            auto_inc_nstreams=32,
#            robust=1,
#            quad_type = "l",
#         pfct_aa_grid_size=37,
#            pfct_method = "interpolate" )
#ws.yCalc( y = y_rt4 )
#ws.Copy( ws.za_grid, ws.za_grid_copy )
#ws.Copy( ws.aa_grid, ws.aa_grid_copy )


#scat_data = [ws.scat_data_raw.value[0]]
#n = len(scat_data[0]) // 2
#scat_data[0] = scat_data[0][n-1 : n+2]
#
#scat_data[0][1].pha_mat_data = scat_data[0][1].pha_mat_data[:, :1, ...]
#scat_data[0][1].ext_mat_data = scat_data[0][1].ext_mat_data[:, :1, ...]
#scat_data[0][1].abs_vec_data = scat_data[0][1].abs_vec_data[:, :1, ...]
#scat_data[0][1].T_grid = scat_data[0][1].T_grid[:1]
#scat_data[0][0] = scat_data[0][1]
#scat_data[0][2] = scat_data[0][1]
#
#meta_data = [ws.scat_meta.value[0]]
#meta_data[0] = meta_data[0][n-1:n+2]
#
#from pyarts.xml import save
#save(scat_data, "test_data/scat_data_simple.xml")
#save(meta_data, "test_data/scat_meta_simple.xml")
#
#import matplotlib.pyplot as plt
#res = np.loadtxt("res.txt")
