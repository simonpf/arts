# - Find the SHTns library
#
# SHTNS_FOUND       - system has shtns
# SHTNS_INCLUDE_DIR - the shtns include directory
# SHTNS_LIBRARIES   - Link these to use shtns
#
# Copyright (c) 2018, Oliver Lemke, <olemke@core-dump.info>

if (FFTW_FOUND AND NOT NO_SHTNS)
  find_path (SHTNS_INCLUDE_DIR shtns.h)

  find_library (SHTNS_LIBRARY NAMES shtns)

  set (SHTNS_LIBRARIES ${SHTNS_LIBRARY})

  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (
    SHTNS DEFAULT_MSG
    SHTNS_LIBRARIES
    SHTNS_INCLUDE_DIR
    )

  mark_as_advanced (SHTNS_INCLUDE_DIR SHTNS_LIBRARY)

  if (SHTNS_FOUND)
    set (ENABLE_SHTNS 1)
  endif()
endif()

