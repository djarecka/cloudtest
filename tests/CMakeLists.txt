# TODO: search for pytest
# TODO: search for cffi
# TODO: search for numpy
# TODO: search for libcloudphxx

foreach(test 
  "blk_1m_libcloudphxx;-s;${CMAKE_SOURCE_DIR}/tests/blk_1m_pytest.py;--libname=libcloudphxx_blk_1m_pytest" 
  "blk_1m_wrfkessler;-s;${CMAKE_SOURCE_DIR}/tests/blk_1m_pytest.py;--libname=wrfkessler_blk_1m_pytest"
  "blk_1m_eqs;-s;${CMAKE_SOURCE_DIR}/tests/blk_1m_eqs_pytest.py"
)
  list(GET test 0 id)
  list(REMOVE_AT test 0)
  add_test(
    NAME ${id}
    COMMAND ${PYTHON_EXECUTABLE} "-m" "pytest" ${test}
  )
  set(PYTHONPATH "${CMAKE_SOURCE_DIR}")
  set(PYTHONPATH "${PYTHONPATH}:${CMAKE_SOURCE_DIR}/libcloudphxx")
  set(PYTHONPATH "${PYTHONPATH}:${CMAKE_SOURCE_DIR}/wrf_microphys/kessler")
  set_tests_properties(${id} PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHONPATH};LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/wrf_microphys/kessler")
endforeach()
