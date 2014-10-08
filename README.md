Work in progress!

Codes for testing microphysical schemes used in atmospheric numerical models.

Simple examples of usage:

> cd tests

1-MOMENT SCHEMES

** libcloudph 1-moment
  - >python2.7 -m pytest -s blk_1m_pytest.py --libname=libcloudphxx_blk_1m_pytest

** kessler from WRF
 - >cd ../wrf_microphys/kessler/ (TODO!)
 - >python2.7 -m pytest -s ../../tests/blk_1m_pytest.py --libname=wrfkessler_blk_1m_pytest

** thermodynamic equation used in libcloudph 1-moment
 - python2.7 -m pytest -s blk_1m_eqs_pytest.py


=============================================


2-MOMENT SCHEMES

** libcloudph 2-moment
 - >python2.7 -m pytest -s blk_2m_pytest.py --libname=libcloudphxx_blk_2m_pytest

** 2-moment Morrison (NOT from WRF!)
 - >cd ../wrf_microphys/morrison_2momNc/ (TODO!)
 - >python2.7 -m pytest -s ../../tests/blk_2m_pytest.py  --libname=wrfmorrison_blk_2m_pytest

