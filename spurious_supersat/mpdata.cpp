#include <boost/python.hpp>
namespace bp = boost::python;

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
using namespace libmpdataxx;

using arr_t = blitz::Array<double, 1>;
using py_ptr_t = long; // TODO: acquire it using some decltype()

void mpdata_cpp(
  const arr_t &arr, 
  const double &C, 
  const int &nt
) {
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
    enum { opts = opts::fct | opts::abs };
  };
  
  using slv_t = solvers::mpdata<ct_params_t>;
  typename slv_t::rt_params_t p;
  p.grid_size = { arr.extent()[0] };
  concurr::serial<slv_t, bcond::open, bcond::open> run(p);
  run.advector() = C;
  run.advectee() = arr;
  run.advance(nt);
  arr(blitz::Range::all()) = run.advectee();
}

void mpdata_py(
  const bp::numeric::array &arg, 
  const double &C, 
  const int &nt
) {
  mpdata_cpp(
    arr_t(                                                             
      // pointer to the data                                                  
      reinterpret_cast<typename arr_t::T_numtype*>(                           
	(py_ptr_t)bp::extract<py_ptr_t>(arg.attr("ctypes").attr("data"))         
      ),                                                                      
      // length of the array (regardless of the original dimensionality, we do 1D)
      blitz::shape(bp::extract<long>(arg.attr("size"))),                      
      // ensure Blitz++ does not try to free the memory when done             
      blitz::neverDeleteData                                                  
    ),
    C,
    nt   
  );
}

BOOST_PYTHON_MODULE(libmpdata)
{
  bp::numeric::array::set_module_and_type("numpy", "ndarray");
  bp::def("mpdata", &mpdata_py);
}
