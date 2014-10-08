message(STATUS "Searching for WRF source tree...")
find_path(
  WRF_SRC_PATH
  phys/module_mp_kessler.F
  HINTS /Users/dorota/WRFV3 
  ENV WRF_SRC_PATH
)
if (NOT WRF_SRC_PATH)
  message(FATAL_ERROR "WRF source tree not found! (you can specify it with WRF_SRC_PATH env. var.)")
else()
  message(STATUS "WRF located in: ${WRF_SRC_PATH}")
endif()

add_library(kessler SHARED ${WRF_SRC_PATH}/phys/module_mp_kessler.F kessler_wrap.f90)
set_target_properties(kessler PROPERTIES COMPILE_FLAGS "-ffree-form -fdefault-real-8")
set_target_properties(kessler PROPERTIES SUFFIX ".so")
