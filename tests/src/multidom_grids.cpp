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

	if (dom_number != 14)
		abort();

	T boundaries [13] = {(T)0.05, (T)0.10, (T)0.15, (T)0.20, (T)0.25, (T) 0.30,
		                 (T)0.35, (T)0.45, (T)0.55, (T)0.65, (T)0.75, (T) 0.85, (T) 0.95};

	aneto::multidomain<T> mesh_unif(dom_number, dom_point, 0.000, +1.00);
	T *f_unif, *ni_unif, *ti_unif;
	f_unif    = new T[mesh_unif.number_points()];
	ti_unif   = new T[mesh_unif.number_points()];
	ni_unif   = new T[mesh_unif.number_points()];

	aneto::multidomain<T> mesh_defB(dom_number, dom_point, 0.000, +1.00, +1,
					aneto::DERGRID_OPTION::SINGLE, aneto::GRID_OPTION::DEF_BNDS, boundaries);
	T *f_defB, *ni_defB, *ti_defB;
	f_defB    = new T[mesh_defB.number_points()];
	ti_defB   = new T[mesh_defB.number_points()];
	ni_defB   = new T[mesh_defB.number_points()];

	for (size_t i = 0; i < mesh_unif.number_points(); i++) {
		f_unif[i]  = cos(mesh_unif.xp(i));
		ti_unif[i] = sin(mesh_unif.xp(i));

		f_defB[i]  = cos(mesh_defB.xp(i));
		ti_defB[i] = sin(mesh_defB.xp(i));
	}

	mesh_unif.comp_integral(f_unif, ni_unif, aneto::INTEG_OPTION::CUSTOM_BC, 0.00, 0.000);
	mesh_defB.comp_integral(f_defB, ni_defB, aneto::INTEG_OPTION::CUSTOM_BC, 0.00, 0.000);

	T max_error_unif = compute_max_error(ni_unif, ti_unif, mesh_unif.number_points()-1);
	T max_error_defB = compute_max_error(ni_defB, ti_defB, mesh_defB.number_points()-1);

	if (verbose) std::clog << start_str << " \t" << toString_5e(max_error_unif) << " \t";
	if (verbose) std::clog << toString_5e(max_error_defB) << " \t";
	if (verbose) std::clog << toString_5e(max_error_defB / max_error_unif) << " \t";
	if (verbose) std::clog << toString_5e(max_error_defB - max_error_unif) << std::endl;


	if (verbose) {
		std::ofstream file(start_str + "_multidom_grids.dat");
		for (size_t i = 0; i < mesh_unif.number_points(); i++){
			file << toString(mesh_unif.xp(i)) << "  "  << toString(f_unif[i]) << "  "  << toString(ni_unif[i]) << "  ";
		    file << toString(mesh_defB.xp(i)) << "  "  << toString(f_defB[i]) << "  "  << toString(ni_defB[i]) << std::endl;
		}
		file.close();
	}

	delete [] f_unif;
	delete [] ni_unif;
	delete [] ti_unif;
	delete [] f_defB;
	delete [] ni_defB;
	delete [] ti_defB;

	return cabs(max_error_unif - max_error_defB);
}




int multidom_grids(bool verbose, size_t spec_domains = 14, size_t spec_points = 63)
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

		for (size_t dec_prec = 5; dec_prec <= 120; dec_prec+=5) {
			mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));
			check_error_integral<mpfr::mpreal>(verbose, spec_domains, spec_points, (boost::format("%04i    %04i") % mpfr_get_default_prec() % dec_prec).str() );
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


