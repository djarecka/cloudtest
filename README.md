Work in progress!

Codes for testing microphysical scheme used in atmospheric numerical models.

Simple examples of usage:

1-MOMENT SCHEMES

** libcloudph 1-moment
 - be sure that in blk_1m_pytest adj_cellwise is imported from libcloudphxx_blk_1m_pytes (TODO!)
 - >python2.7 -m pytest -s blk_1m_pytest.py

** kessler from WRF
 -  be sure that in blk_1m_pytest adj_cellwise is imported  from wrfkessler_blk_1m_pytest (TODO!)
 - >cd wrf_microphys/kessler/ (TODO!)
 - >python2.7 -m pytest -s ../../blk_1m_pytest.py

** thermodynamic equation used in libcloudph 1-moment
 - python2.7 -m pytest -s blk_1m_eqs_pytest.py


=============================================


2-MOMENT SCHEMES

** libcloudph 2-moment
 - be sure that in blk_2m_pytest adj_cellwise is imported from libcloudphxx_blk_2m_pytest (TODO!)
 - >python2.7 -m pytest -s blk_2m_pytest.py

** 2-moment Morrison (NOT from WRF!)
 -  be sure that in blk_2m_pytest adj_cellwise is imported  from wrfmorrison_blk_2m_pytest (TODO!)
 - >cd wrf_microphys/morrison_2momNc/ (TODO!)
 - >python2.7 -m pytest -s ../../blk_2m_pytest.py

