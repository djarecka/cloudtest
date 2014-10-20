The Morrison 2moment microphysical scheme written in Fortran95 and used in WRF is called from Python using the ccfi library.

Preparing the fortran file:
  * getting the WRF model from http://www2.mmm.ucar.edu/wrf/users/download/get_sources.html
  * the Morrison schme is in the `phys` directory: `module_mp_morr_two_moment.F`
  * running cpp to obtain .f90 file (TODO):
   
    `$ cpp -P -traditional -DEM_CORE=1 -DNMM_CORE=0 -DCOAMPS_CORE=0 -DDA_CORE=0 -DIWORDSIZE=4 -DDWORDSIZE=8 -DRWORDSIZE=4 -DLWORDSIZE=4 -DNONSTANDARD_SYSTEM_FUNC -DWRF_USE_CLM -DNO_IEEE_MODULE  -DDM_PARALLEL -DSTUBMPI module_mp_morr_two_moment.F > module_mp_morr_two_moment.f90`
   
   * adjusting the fortran file to cffi/python codes:

    `$ patch < module_mp_morr_two_moment_adjustment.txt`

Making a share library - libmorrison_2mom.so (using gfortran and double precision):

    $ gfortran -fdefault-real-8 -c -fPIC module_mp_morr_two_moment.f90 -o module_mp_morr_two_moment.o
    $ gfortran -fdefault-real-8 -shared -fPIC module_mp_morr_two_moment.o morrison_2mom_wrap.f90 -o libmorrison_2mom.so

Simple example of usage:
 
    $ python morrison_2momNc_call.py