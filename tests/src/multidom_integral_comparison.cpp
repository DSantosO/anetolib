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
#include <fstream>
#include <cstddef>
#include <chrono>
#include <boost/format.hpp>

#if QUAD == 1
	#include <boost/multiprecision/float128.hpp>
#endif
#include "tests/include/mpreal.h"

#include "aneto/multidomain.hpp"
#include "tests/include/auxiliars.hpp"


template<typename T> static T check_error_integral(bool verbose, size_t dom_number, size_t dom_point, std::string start_str)
{
		aneto::multidomain<T> mesh(dom_number, dom_point, 0.000, +1.00);

		T *function, *num_integ, *theor_integ;
		function    = new T[mesh.number_points()];
		theor_integ = new T[mesh.number_points()];
		num_integ   = new T[mesh.number_points()];

		for (size_t i = 0; i < mesh.number_points(); i++) {
			function[i]    = cos(mesh.xp(i));
			theor_integ[i] = sin(mesh.xp(i));
		}

		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < 1; i++)
			mesh.comp_integral(function, num_integ, aneto::INTEG_OPTION::CUSTOM_BC, 0.00, 0.000);
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1.;

		T max_error = compute_max_error(num_integ, theor_integ, mesh.number_points()-1);
		if (verbose) std::clog << start_str << " \t" << toString_5e(max_error) << " \t"  << duration << std::endl;

		delete [] function;
		delete [] num_integ;
		delete [] theor_integ;

		return max_error;
}




int multidom_integral_comparison(bool verbose, size_t spec_domains = 14, size_t spec_points = 63)
{
	int ERROR;

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "########################## Test Multidomain:  Integrals Precision Comparison ##########################" << std::endl;
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << "# Using a multigrid of [" << spec_domains << ", " << spec_points << "]" << std::endl;
	if (verbose) std::clog << "#bits  decimal     error        time(microsec)" << std::endl;


	if (verbose) {
		check_error_integral<float>(verbose, spec_domains, spec_points,  "#float   008");
		check_error_integral<double>(verbose, spec_domains, spec_points, "#double  016");
		check_error_integral<long double>(verbose, spec_domains, spec_points, "#longdo  016");
		#if QUAD == 1
			check_error_integral<boost::multiprecision::float128>(verbose, spec_domains, spec_points, "#quad    034");
		#endif

		for (size_t bit_prec = 5; bit_prec <= 500; bit_prec+=15) {
			mpfr::mpreal::set_default_prec(bit_prec);
			check_error_integral<mpfr::mpreal>(verbose, spec_domains, spec_points, (boost::format("%04i    %04i") % mpfr_get_default_prec() % mpfr::bits2digits(bit_prec)).str() );
		}
	}

	if (verbose) std::clog << std::endl << "# Using a multigrid of [" << 14 << ", " << 63 << "]" << std::endl;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(150));
	mpfr::mpreal error150 = check_error_integral<mpfr::mpreal>(verbose, 14, 63, (boost::format("%04i    %04i") % mpfr_get_default_prec() % 150).str() );

	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);
	if (error150 < 1e-149) {
		std::cout << "----- Test Multidomain Integral Comparison      ('multidom_integral_comparison'):    OK!\n";
		ERROR = 0;
	}
	else {
		std::cout << "----- Test Multidomain Integral Comparison      ('multidom_integral_comparison'):    Fail! Error in 499 bit integral is greater than expected.\n";
		ERROR = 1;
	}
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << std::endl;

	return ERROR;

}


