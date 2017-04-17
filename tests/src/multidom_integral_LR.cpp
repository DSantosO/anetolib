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

#include "tests/include/mpreal.h"

#include "aneto/spectral_domain.hpp"
#include "aneto/multidomain.hpp"
#include <tests/include/auxiliars.hpp>

template<typename T> static T integral_comparison (bool verbose, size_t dom_number, size_t domain_points, std::string filename)
{
	T *function, *num_integ_left, *num_integ_right, *num_integ_R_thr_L, *num_integ_custom;

	aneto::multidomain<T> mesh(dom_number, domain_points, -1, +1);
	aneto::spectral_domain<T> spec_grid = mesh.get_spectral_domain(0);

	function  = new T[mesh.number_points()];
	num_integ_left     = new T[mesh.number_points()];
	num_integ_right  = new T[mesh.number_points()];
	num_integ_R_thr_L   = new T[mesh.number_points()];
	num_integ_custom   = new T[mesh.number_points()];


	/* Exponential */
	for (size_t i = 0; i < mesh.number_points(); i++) {
// 		function[i]  = exp(- tan(mesh.xp(i) * spec_grid.pi()/((T)4) + spec_grid.pi()/((T)4) ));
			function[i]  = exp(- mesh.xp(i));
	}

	mesh.comp_integral(function, num_integ_left);
	mesh.comp_integral(function, num_integ_right, aneto::INTEG_OPTION::RIGHT_BC);
	mesh.comp_integral(function, num_integ_custom, aneto::INTEG_OPTION::CUSTOM_BC);

	for (size_t i = 0; i < mesh.number_points(); i++) {
		num_integ_R_thr_L[i] = num_integ_left[mesh.number_points() -1] - num_integ_left[i];
	}


	if (verbose) {
		std::ofstream ofile (filename);
		ofile << "# Col 1: coordinate x = [-1:+1]"  << std::endl;
		ofile << "# Col 2: function   f(x) = exp(- tan( - x * PI/4  + PI/4))"  << std::endl;
		ofile << "# Col 3: Integral Left:  IL"  << std::endl;
		ofile << "# Col 4: IL(pi/2) - IL(x)"  << std::endl;
		ofile << "# Col 5: Integral Right."  << std::endl;

		for (size_t i = 0; i < mesh.number_points(); i++) {
			ofile << toString(mesh.xp(i)) << "  ";
			ofile << toString(function[i]) << "  ";
			ofile << toString(num_integ_left[i]) << "  ";
			ofile << toString(num_integ_R_thr_L[i]) << "  ";
			ofile << toString(num_integ_right[i])  << "  " ;
			ofile << toString(num_integ_custom[i]) << std::endl;
		}
		ofile.close();

		std::clog << "# File generated at: " << filename << std::endl;
		std::clog << "# Col 1: coordinate x = [-1:+1]"  << std::endl;
		std::clog << "# Col 2: function   f(x) = exp(- tan( - x * PI/4  + PI/4))"  << std::endl;
		std::clog << "# Col 3: Integral Left:  IL"  << std::endl;
		std::clog << "# Col 4: IL(pi/2) - IL(x)"  << std::endl;
		std::clog << "# Col 5: Integral Right."  << std::endl;
	}

	/* Error to return */
	T max_error = compute_max_error<T>(num_integ_R_thr_L, num_integ_right, mesh.number_points() - 1);

	delete[] function;
	delete[] num_integ_left;
	delete[] num_integ_right;
	delete[] num_integ_R_thr_L;

	return max_error;
}




int multidom_integral_LR(bool verbose, size_t spec_domains = 14, size_t spec_points = 63)
{
	int ERROR;

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "############################## Test Multidomain:  Left / Right Integrals ##############################" << std::endl;
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "# Using a multigrid of [" << spec_domains << ", " << spec_points << "]" << std::endl << std::endl;

	double error = integral_comparison<double>(verbose, spec_domains, spec_points, "multidom_integral_LR.dat");

	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);
	if (error < 1e-15) {
		std::cout << "----- Test Multidomain Integral Left / Right    ('multidom_integral_LR'):            OK!\n";
		ERROR = 0;
	}
	else{
		std::cout << "----- Test Multidomain Integral Left / Right    ('multidom_integral_LR'):            Fail! Difference between both computations is greater than expected.\n";
		ERROR = 1;
	}
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << std::endl;

	return ERROR;

}
