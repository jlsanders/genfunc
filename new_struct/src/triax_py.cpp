#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface.h"
#include <gsl/gsl_poly.h>
#include "utils.h"
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"
#include "orbit.h"
#include "aa.h"
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
using namespace boost::python;

namespace {
	boost::python::numeric::array Stack_Triax_Forces(boost::python::numeric::array f){
		StackelTriaxial T(3.61/500.,-30.,-20.);
		VecDoub x(3,0); for(int i=0;i<3;i++)x[i]=extract<double>(f[i]);
		VecDoub v = T.Forces(x);
    	return boost::python::numeric::array(make_tuple(v[0],v[1],v[2]));
	}
	double Stack_Triax_H(boost::python::numeric::array f){
		StackelTriaxial T(3.61/500.,-30.,-20.);
		VecDoub x(6,0);for(int i=0;i<6;i++)x[i]=extract<double>(f[i]);
		return T.H(x);
	}
	boost::python::numeric::array Stack_Triax_Actions(boost::python::numeric::array f){
		StackelTriaxial T(3.61/500.,-30.,-20.);
		VecDoub x(6,0);for(int i=0;i<6;i++)x[i]=extract<double>(f[i]);
		Actions_TriaxialStackel AA(&T);
		VecDoub v = AA.actions(x);
    	return boost::python::numeric::array(make_tuple(v[0],v[1],v[2]));
	}

	boost::python::numeric::array Stack_Triax_Freqs(boost::python::numeric::array f){
		StackelTriaxial T(3.61/500.,-30.,-20.);
		VecDoub x(6,0);for(int i=0;i<6;i++)x[i]=extract<double>(f[i]);
		Actions_TriaxialStackel AA(&T);
		VecDoub v = AA.actions(x,1);
    	return boost::python::numeric::array(make_tuple(v[3],v[4],v[5]));
	}
}

BOOST_PYTHON_MODULE(triax_py){
	boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    def("Stack_Triax_Forces", Stack_Triax_Forces);
    def("Stack_Triax_H", Stack_Triax_H);
    def("Stack_Triax_Actions", Stack_Triax_Actions);
    def("Stack_Triax_Freqs", Stack_Triax_Freqs);
}	