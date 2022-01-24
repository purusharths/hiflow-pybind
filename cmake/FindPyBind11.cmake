# - Try to find PyBind11

find_path (PyBind11 pybind11.h DOC "Pybind11 Include Directory")

set(VCL_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/contrib/VCL)

IF(EXISTS ${PB11_DIR}/pybind11.h)
  SET(PB11_FOUND YES)
  find_path(PB11_INCLUDE_DIR pybind11.h PATH ${PB11_DIR})
ELSE(EXISTS ${PB11_DIR}/pybind11.h)
  SET(PB11_FOUND NO)
ENDIF(EXISTS ${PB11_DIR}/pybind11.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Vcl DEFAULT_MSG VCL_INCLUDE_DIR)
