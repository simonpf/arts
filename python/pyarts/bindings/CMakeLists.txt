include_directories (
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/3rdparty/
  ${CMAKE_SOURCE_DIR}/3rdparty/Eigen_unsupported
  ${CMAKE_SOURCE_DIR}/3rdparty/scattering/include
  ${CMAKE_SOURCE_DIR}/3rdparty/netcdfhpp/include
  ${CMAKE_SOURCE_DIR}/3rdparty/scattering/ext/shtns
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(SOURCE_FILES
  bindings.cc
  scattering_particle.cc
  single_scattering_data.cc
  xml_io.cc
  messages.cc)
pybind11_add_module(cxx ${SOURCE_FILES})
target_link_libraries(cxx PUBLIC arts_api)
