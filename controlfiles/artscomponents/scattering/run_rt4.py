import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

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

ws.ReadXML(ws.scat_data_raw, "test_data/scat_data_single_t.xml")

ws.scat_data_raw = ws.scat_data_raw.value
ws.ReadXML(ws.scat_meta, "test_data/scat_meta.xml")
ws.scat_meta = ws.scat_meta.value

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

 #Scat species 1
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
                                      scattering_data=ws.scat_data_raw.value[1],
                                      scattering_meta_data=ws.scat_meta.value[1],
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
ws.IndexSet(ws.stokes_dim, 2)
# Scattering data tailored to these frequencies, so don't change!
ws.VectorSet(ws.f_grid, np.array([3.15e+10, 1.65e+11, 6.66e+11]))
ws.VectorSet(ws.f_grid, np.array([6.66e+11]))
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
#
# Extract simplified scattering data
#
ws.DOAngularGridsSet( N_za_grid=38, N_aa_grid=37 )

ws.Copy( ws.za_grid_copy, ws.za_grid )
ws.Copy( ws.aa_grid_copy, ws.aa_grid )
ws.RT4Calc( nstreams = 32,
            auto_inc_nstreams=0,
            robust=1,
            quad_type = "l",
            pfct_aa_grid_size=37,
            pfct_method = "interpolate" )
ws.yCalc( y = ws.y_rt4 )
ws.Copy( ws.za_grid, ws.za_grid_copy )
ws.Copy( ws.aa_grid, ws.aa_grid_copy )




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

sd = ws.scat_data.value
pf = sd[0][0].pha_mat_data[0, 0, :, 0, 0, 0, 0]
y_0 = pf
x_0 = sd[0][0].za_grid
x = ws.za_grid.value[16-1::-1]
x_i = np.linspace(0, np.pi, 32)
mu_i = -np.cos(x_i)
y = np.interp(x, x_0, y_0)

s = np.array([
4.72144e-06
    ,5.27914e-06
    ,5.62248e-06
    ,5.94807e-06
    ,6.53917e-06
    ,7.76194e-06
    ,9.35437e-06
    ,8.98869e-06
    ,1.75357e-05
    ,1.71107e-05
    ,2.31398e-05
    ,7.95722e-05
    ,2.62074e-05
    ,0.000436247
    ,0.0025088
    ,0.00461849
])

s_i = np.concatenate([s[::-1], s])

s_0 = np.array([
0.00025996
    ,0.000172796
    ,4.75674e-05
    ,1.55194e-05
    ,1.69143e-05
    ,1.06839e-05
    ,9.44047e-06
    ,8.88711e-06
    ,7.52712e-06
    ,7.3433e-06
    ,6.96316e-06
    ,6.63064e-06
    ,6.43655e-06
    ,6.40781e-06
    ,6.28984e-06
    ,5.98151e-06
    ])

s_0_i = np.concatenate([s_0[::-1], s_0])

s_ref = np.array([
5.903e-06
    ,6.60026e-06
    ,7.02951e-06
    ,7.43659e-06
    ,8.17561e-06
    ,9.70438e-06
    ,1.16953e-05
    ,1.12381e-05
    ,2.19241e-05
    ,2.13927e-05
    ,2.89305e-05
    ,9.94853e-05
    ,3.27659e-05
    ,0.000545418
    ,0.00313664
    ,0.00577427
])
s_ref_i = np.concatenate([s_ref[::-1], s_ref])

s_0_ref = np.array([
    0.000256527
    ,0.000170514
    ,4.69392e-05
    ,1.53145e-05
    ,1.66909e-05
    ,1.05428e-05
    ,9.31581e-06
    ,8.76976e-06
    ,7.42773e-06
    ,7.24633e-06
    ,6.87121e-06
    ,6.54308e-06
    ,6.35155e-06
    ,6.3232e-06
    ,6.20678e-06
    ,5.90252e-06
    ])
s_0_ref_i = np.concatenate([s_0_ref[::-1], s_0_ref])

mu = np.cos(np.deg2rad(x))
#
from pyarts.xml import load,save
import numpy as np
sd = load("test_data/scat_data_single_t.xml")


sca_coeff_0 = sd[0][-1].ext_mat_data - sd[0][-1].abs_vec_data
def reduce_t(data):
    t_grid = data.T_grid
    i = t_grid.size // 2
    t_grid = t_grid[[i]]

    print(data)
    data.T_grid = t_grid
    data.pha_mat_data = data.pha_mat_data[:, i:i+1]
    data.ext_mat_data = data.ext_mat_data[:, i:i+1]
    data.abs_vec_data = data.abs_vec_data[:, i:i+1]

    data.za_grid = sd[0][0].za_grid
    sca_coeff = data.ext_mat_data - data.abs_vec_data
    data.pha_mat_data = np.copy(sd[0][-1].pha_mat_data)
    data.pha_mat_data *= (sca_coeff / sca_coeff_0)[..., np.newaxis, np.newaxis]


for h in sd:
    for p in h:
        reduce_t(p)
        print(p.T_grid.size)
