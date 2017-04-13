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

#include <boost/format.hpp>
#include <chrono>

#if QUAD == 1
	#include <boost/multiprecision/float128.hpp>
#endif
#include "tests/include/mpreal.h"

#include "aneto/spectral_domain.hpp"
#include <tests/include/auxiliars.hpp>


template<typename T> static T check_error_integral (bool verbose, size_t spec_points, std::string start_str)
{
		aneto::spectral_domain<T> spec(spec_points);

		T *function, *num_integ, *theor_integ;
		function    = new T[spec_points + 1];
		theor_integ = new T[spec_points + 1];
		num_integ   = new T[spec_points + 1];


		/* Cosinus */
		for (size_t i = 0; i <= spec_points; i++) {
			function[i] = cos(spec.get_X(i));
			theor_integ[i] = sin(spec.get_X(i)) - sin(spec.get_X(0)) ;
		}

		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		spec.compute_integral(function, num_integ);
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


		T max_error = compute_max_error(num_integ, theor_integ, spec_points);

		if (verbose) std::clog << start_str << " \t" << toString_5e(max_error) << " \t"  << duration << std::endl;

		delete[] function;
		delete[] theor_integ;
		delete[] num_integ;

		return max_error;
}




int singledom_integral_comparison(bool verbose, size_t spec_points)
{
	int ERROR;

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "########################## Test Singledomain:  Integrals Precision Comparison ##########################" << std::endl;
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << "# Using a single domain of " << spec_points << " points." << std::endl;
	if (verbose) std::clog << "#bits  decimal     error        time(microsec)" << std::endl;

	if (verbose) {
		check_error_integral<float>(verbose, spec_points,  "#float   008");
		check_error_integral<double>(verbose, spec_points, "#double  016");
		check_error_integral<long double>(verbose, spec_points, "#longdo  016");
		#if QUAD == 1
			check_error_integral<boost::multiprecision::float128>(verbose, spec_points, "#quad    034");
		#endif

		for (size_t dec_prec = 5; dec_prec <= 120; dec_prec+=5) {
			mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));
			check_error_integral<mpfr::mpreal>(verbose, spec_points, (boost::format("%04i    %04i") % mpfr_get_default_prec() % dec_prec).str() );
		}
	}


	if (verbose) std::clog  << std::endl << "# Using a single domain of " << 340 << " points." << std::endl;
	if (verbose) std::clog << "#bits  decimal     error        time(microsec)" << std::endl;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(150));
	mpfr::mpreal error150 = check_error_integral<mpfr::mpreal>(verbose, 340, (boost::format("%04i    %04i") % mpfr_get_default_prec() % 150).str() );

	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);
	if (error150 < 1e-148) {
		std::cout << "----- Test Singledomain Integral Comparison     ('singledom_integral_comparison'):   OK!\n";
		ERROR = 0;
	}
	else {
		std::cout << "----- Test Singledomain Integral Comparison     ('singledom_integral_comparison'):   Fail! Error in 499 bit integral is greater than expected.\n";
		ERROR = 1;
	}
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << std::endl;


	return ERROR;

}
