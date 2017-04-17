/*
    ANETO C++ Library - Arbitrary precisioN solvEr with pseudo-specTral MethOds
    Project homepage: https://github.com/DSantosO/anetolib
    Copyright (c) 2017 Daniel Santos-Olivan and Carlos F. Sopuerta
    **************************************************************************
    ANETO is under GNU General Public License ("GPL").

    GNU General Public License ("GPL") copyright permissions statement:
    This file is part of ANETO

    ANETO is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ANETO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <cstddef>
#include <chrono>
#include <boost/format.hpp>

#if QUAD == 1
	#include <boost/multiprecision/float128.hpp>
#endif
#if QD == 1
	#include <qd/fpu.h>
	#include <qd/dd_real.h>
	#include <qd/qd_real.h>
#endif

#include "tests/include/mpreal.h"

#include "aneto/spectral_domain.hpp"
#include "tests/include/auxiliars.hpp"


template<typename T> static T check_error_derivative (bool verbose, size_t spec_points, std::string start_str)
{
		aneto::spectral_domain<T> spec(spec_points);

		T *function, *num_der, *theor_der;
		function    = new T[spec_points + 1];
		theor_der = new T[spec_points + 1];
		num_der   = new T[spec_points + 1];


		/* Cosinus */
		for (size_t i = 0; i <= spec_points; i++) {
			function[i] = log(((T)1.000) + spec.get_X(i) * spec.get_X(i))/((T)(2.000));
			theor_der[i] = (spec.get_X(i))/(((T)1.000) + spec.get_X(i) * spec.get_X(i));
		}

		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < 1; i++)
			spec.compute_1st_der(function, num_der);
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


		T max_error = compute_max_error(num_der, theor_der, spec_points);

		if (verbose) std::clog << start_str << " \t" << toString_5e(max_error) << " \t"  << duration /1. << std::endl;

		delete[] function;
		delete[] theor_der;
		delete[] num_der;

		return max_error;
}





int singledom_derivatives(bool verbose, size_t spec_points)
{

	int ERROR;

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "########################## Test Singledomain:  Derivative Precision Comparison ##########################" << std::endl;
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl<<  "# Using a single domain of " << spec_points << " points." << std::endl ;
	if (verbose) std::clog << "#bits  decimal     error        time(microsec)" << std::endl;

	if (verbose) {
		check_error_derivative<float>(verbose, spec_points,  "#float   008");
		check_error_derivative<double>(verbose, spec_points, "#double  016");
		check_error_derivative<long double>(verbose, spec_points, "#longdo  016");
		#if QUAD == 1
			check_error_derivative<boost::multiprecision::float128>(verbose, spec_points, "#quad    034");
		#endif
		#if QD == 1
			check_error_derivative<dd_real>(verbose, spec_points, "#dd_real 032");
			check_error_derivative<qd_real>(verbose, spec_points, "#qd_real 064");
		#endif


// 		for (size_t dec_prec = 5; dec_prec <= 120; dec_prec+=5) {
		for (size_t dec_prec = 30; dec_prec <= 120; dec_prec+=30) {
			mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));
			check_error_derivative<mpfr::mpreal>(verbose, spec_points, (boost::format("%04i    %04i") % mpfr_get_default_prec() % dec_prec).str() );
		}
	}


	if (verbose) std::clog  << std::endl << "# Using a single domain of " << 400 << " points." << std::endl;
	if (verbose) std::clog << "#bits  decimal     error        time(microsec)" << std::endl;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(150));
	mpfr::mpreal error150 = check_error_derivative<mpfr::mpreal>(verbose, 400, (boost::format("%04i    %04i") % mpfr_get_default_prec() % 150).str() );

	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);
	if (error150 < 1e-144) {
		std::cout << "----- Test Singledomain Derivative Comparison   ('singledom_derivatives'):           OK!\n";
		ERROR = 0;
	}
	else {
		std::cout << "----- Test Singledomain Derivative Comparison   ('singledom_derivatives'):           Fail! Error in 499 bit integral is greater than expected.\n";
		ERROR = 1;
	}
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << std::endl;


	return ERROR;









}
