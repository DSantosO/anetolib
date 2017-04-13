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
#include <vector>
#include <boost/format.hpp>

#include "tests/include/mpreal.h"

#include "aneto/spectral_domain.hpp"
#include "aneto/multidomain.hpp"
#include "tests/include/auxiliars.hpp"


template<typename T> std::vector<T> check_error_dual(bool verbose, size_t dom_number, size_t dom_points, std::string start_str, bool return_check = false)
{

	std::ofstream file;

	T *dom_func, *dom_num_der, *dom_theor_der, *dom_num_d2, *dom_theor_d2;
	T err_norm_d, err_dual_d, err_norm_d2, err_dual_d2;

	aneto::multidomain<T> multdom(dom_number, dom_points, 0.00, +1.00);
	aneto::multidomain<T> multdom_dual(dom_number, dom_points, 0.00, +1.00, 1, aneto::DERGRID_OPTION::DUAL);

	dom_func      = new T[multdom.number_points()];
	dom_theor_der = new T[multdom.number_points()];
	dom_num_der   = new T[multdom.number_points()];

	dom_theor_d2 = new T[multdom.number_points()];
	dom_num_d2   = new T[multdom.number_points()];




	for (size_t i = 0; i < multdom.number_points(); i++) {
		dom_func[i] = exp(-multdom.xp(i) * multdom.xp(i));
		dom_theor_der[i] = - ((T)2.00) * multdom.xp(i) * exp(-multdom.xp(i) * multdom.xp(i));
		dom_theor_d2[i] = - (T)2 * exp(-multdom.xp(i) * multdom.xp(i)) + (T)4 * exp(-multdom.xp(i) * multdom.xp(i)) * multdom.xp(i) * multdom.xp(i);
	}

	multdom.comp_derivative(dom_func, dom_num_der);
	multdom.comp_derivative(dom_num_der, dom_num_d2);

	err_norm_d  = compute_max_error<T>(dom_num_der, dom_theor_der, multdom.number_points()-10, 10);
	err_norm_d2 = compute_max_error<T>(dom_num_d2, dom_theor_d2, multdom.number_points()-10, 10);


	if (verbose) {
		file.open(start_str + "_multidom.dat");
		for (size_t i = 0; i < multdom.number_points(); i++) {
			file << toString(multdom.xp(i)) << "  ";
			file << toString(dom_theor_der[i]) << "  ";
			file << toString(dom_num_der[i]) << "  ";
			file << toString(dom_theor_d2[i]) << "  ";
			file << toString(dom_num_d2[i]) << "  ";
			file << toString(cabs<T>(dom_num_der[i] - dom_theor_der[i])) <<  "  ";
			file << toString(cabs<T>(dom_num_d2[i] - dom_theor_d2[i]));
			file << std::endl;
		}
		file.close();
	}


	multdom_dual.comp_derivative(dom_func, dom_num_der);
	multdom_dual.comp_derivative(dom_num_der, dom_num_d2);

	err_dual_d  = compute_max_error<T>(dom_num_der, dom_theor_der, multdom.number_points()-10, 10);
	err_dual_d2 = compute_max_error<T>(dom_num_d2, dom_theor_d2, multdom.number_points()-10, 10);

	if (verbose) std::clog << start_str << "  \t";
	if (verbose) std::clog << toString_5e(err_norm_d) << "  ";
	if (verbose) std::clog << toString_5e(err_dual_d)  << "  ";
	if (verbose) std::clog << toString_5e(err_norm_d / err_dual_d) << "  ";
	if (verbose) std::clog << toString_5e(err_norm_d2) << "  ";
	if (verbose) std::clog << toString_5e(err_dual_d2) << "  ";
	if (verbose) std::clog << toString_5e(err_norm_d2 / err_dual_d2) << "  ";
	if (verbose) std::clog <<  std::endl;



	if (verbose) {
		file.open(start_str + "_dual.dat");
		for (size_t i = 0; i < multdom_dual.number_points(); i++) {
			file << toString(multdom.xp(i)) << "  ";
			file << toString(dom_theor_der[i]) << "  ";
			file << toString(dom_num_der[i]) << "  ";
			file << toString(dom_theor_d2[i]) << "  ";
			file << toString(dom_num_d2[i]) << "  ";
			file << toString(cabs<T>(dom_num_der[i] - dom_theor_der[i])) <<  "  ";
			file << toString(cabs<T>(dom_num_d2[i] - dom_theor_d2[i]));
			file << std::endl;
		}
		file.close();
	}

	delete[] dom_func;
	delete[] dom_theor_der;
	delete[] dom_num_der;
	delete[] dom_theor_d2;
	delete[] dom_num_d2;


	/* Checked values */
	if (return_check == true){
		std::vector<T> vec(2);
		vec[0] = err_norm_d / err_dual_d;
		vec[1] = err_norm_d2 / err_dual_d2;
		return vec;
	}

	return std::vector<T>();

}






int multidom_derivative_dual(bool verbose, size_t dom_number = 14, size_t dom_points = 63)
{
	int ERROR;

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "################################# Test Multidomain:  Derivative Dual  #################################" << std::endl;
	print_sep_line(verbose, std::clog);

	if (verbose) std::clog << std::endl << "##type [domain, points]    errD(Norm)   errD(Dual)   Improv      errD2(Norm)   errD(Dual)    Improv" << std::endl;
	if (verbose) std::clog << std::endl << "###### Default Examples" << std::endl;
	check_error_dual<double>(verbose, 10, 12, "double_D10_N12");

	check_error_dual<double>(verbose, 10, 24, "double_D10_N24");
	check_error_dual<double>(verbose, 10, 36, "double_D10_N32");
	check_error_dual<double>(verbose, 10, 48, "double_D10_N48");

	std::vector<double> check_1060 = check_error_dual<double>(verbose, 10, 60, "double_D10_N60", true);


	if (verbose) std::clog <<  std::endl << "###### User Example" << std::endl;
	std::string str_start = (boost::format("double_D%d_N%d") % dom_number % dom_points ).str();
	check_error_dual<double>(verbose, dom_number, dom_points, str_start);

	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);

	if (check_1060[0] > 40 && check_1060[1] > 900)
	{
		std::cout << "----- Test Multidomain Derivative Dual          ('multidom_derivative_dual'):        OK!\n";
		ERROR = 0;
	}
	else {
		std::cout << "----- Test Multidomain Derivative Dual          ('multidom_derivative_dual'):        Fail! Improvements in double [10,60] case are less than expected.\n";
		ERROR = 1;
	}

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << std::endl;


	return ERROR;

}
