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


#include <fstream>
#include <vector>
#include <boost/format.hpp>

#include "tests/include/mpreal.h"

#include "aneto/spectral_domain.hpp"
#include "aneto/multidomain.hpp"
#include "tests/include/auxiliars.hpp"

using mpfr::mpreal;

template<typename T> std::vector<T> check_error_precision(bool verbose, size_t dom_number, size_t dom_points, std::string start_str, bool return_check = false)
{

	std::ofstream file;

	T *dom_func, *dom_theor_der, *dom_theor_d2;
	T *dom_norm_der, *dom_norm_d2, *dom_HP_der, *dom_HP_d2;
	T err_norm_d, err_HP_d, err_norm_d2, err_HP_d2;

	aneto::multidomain<T> multdom(dom_number, dom_points, 0.00, +1.00);


	dom_func      = new T[multdom.number_points()];
	dom_theor_der = new T[multdom.number_points()];
	dom_theor_d2 = new T[multdom.number_points()];

	dom_norm_der  = new T[multdom.number_points()];
	dom_norm_d2   = new T[multdom.number_points()];
	dom_HP_der    = new T[multdom.number_points()];
	dom_HP_d2     = new T[multdom.number_points()];


	mpfr::mpreal::set_default_prec(mpfr::digits2bits(300));
	aneto::multidomain<mpreal> multdom_HP(dom_number, dom_points, 0.00, +1.00);
	mpreal *intermid_fun, *intermid_der, *intermid_d2;
	intermid_fun  = new mpreal[multdom_HP.number_points()];
	intermid_der  = new mpreal[multdom_HP.number_points()];
	intermid_d2   = new mpreal[multdom_HP.number_points()];


	for (size_t i = 0; i < multdom.number_points(); i++) {
// 		dom_func[i] = exp(-multdom.xp(i) * multdom.xp(i));
		dom_func[i] = multdom.xp(i);
		dom_theor_der[i] = - ((T)2.00) * multdom.xp(i) * exp(-multdom.xp(i) * multdom.xp(i));
		dom_theor_d2[i] = - (T)2 * exp(-multdom.xp(i) * multdom.xp(i)) + (T)4 * exp(-multdom.xp(i) * multdom.xp(i)) * multdom.xp(i) * multdom.xp(i);
	}

	multdom.comp_derivative(dom_func, dom_norm_der);
	multdom.comp_derivative(dom_norm_der, dom_norm_d2);

	err_norm_d  = compute_max_error<T>(dom_norm_der, dom_theor_der, multdom.number_points()-10, 10);
	err_norm_d2 = compute_max_error<T>(dom_norm_d2, dom_theor_d2, multdom.number_points()-10, 10);


	if (verbose) {
		file.open(start_str + "_normal.dat");
		for (size_t i = 0; i < multdom.number_points(); i++)
			file << toString(multdom.xp(i)) << "  " << toString(cabs<T>(dom_norm_der[i] - dom_theor_der[i])) <<  "  " << toString(cabs<T>(dom_norm_d2[i] - dom_theor_d2[i])) << std::endl;
		file.close();
	}


	for (size_t i = 0; i < multdom.number_points(); i++)
		intermid_fun[i] = multdom_HP.xp(i);
// 		intermid_fun[i] = (mpreal) dom_func[i];

	multdom_HP.comp_derivative(intermid_fun, intermid_der);
	multdom_HP.comp_derivative(intermid_der, intermid_d2);

	for (size_t i = 0; i < multdom.number_points(); i++) {
		dom_HP_der[i] = intermid_der[i].toDouble();
		dom_HP_d2[i] = intermid_d2[i].toDouble();
	}


	err_HP_d  = compute_max_error<T>(dom_HP_der, dom_theor_der, multdom.number_points()-10, 10);
	err_HP_d2 = compute_max_error<T>(dom_HP_d2,  dom_theor_d2,  multdom.number_points()-10, 10);

	if (verbose) std::clog << start_str << "  \t";
	if (verbose) std::clog << toString_5e(err_norm_d) << "  ";
	if (verbose) std::clog << toString_5e(err_HP_d)  << "  ";
	if (verbose) std::clog << toString_5e(err_norm_d / err_HP_d) << "  ";
	if (verbose) std::clog << toString_5e(err_norm_d2) << "  ";
	if (verbose) std::clog << toString_5e(err_HP_d2) << "  ";
	if (verbose) std::clog << toString_5e(err_norm_d2 / err_HP_d2) << "  ";
	if (verbose) std::clog <<  std::endl;



	if (verbose) {
		file.open(start_str + "_HP.dat");
		for (size_t i = 0; i < multdom_HP.number_points(); i++)
			file << toString(multdom.xp(i)) << "  " << toString(cabs<T>(dom_HP_der[i] - dom_theor_der[i])) << "  " << toString(cabs<T>(dom_HP_d2[i] - dom_theor_d2[i])) << std::endl;
		file.close();
	}

	if (verbose) {
		file.open(start_str + "_TEST.dat");
		for (size_t i = 0; i < multdom_HP.number_points(); i++)
			file << toString(multdom.xp(i)) << "  " << toString(dom_norm_der[i]) << "  " << toString(dom_HP_der[i]) << "  " << toString(intermid_der[i]) << std::endl;
		file.close();
	}

	delete[] dom_func;
	delete[] dom_theor_der;
	delete[] dom_theor_d2;
	delete[] dom_norm_der;
	delete[] dom_norm_d2;
	delete[] dom_HP_der;
	delete[] dom_HP_d2;
	delete[] intermid_fun;
	delete[] intermid_der;
	delete[] intermid_d2;

	/* Checked values */
	if (return_check == true){
		std::vector<T> vec(2);
		vec[0] = err_norm_d / err_HP_d;
		vec[1] = err_norm_d2 / err_HP_d2;
		return vec;
	}

	return std::vector<T>();

}






int multidom_derivative_precision(bool verbose, size_t dom_number = 14, size_t dom_points = 63)
{
	int ERROR;

	std::clog << "Warning!!! UNSUPPORTED!!!!!" << std::endl;


	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "################################# Test Multidomain:  Derivative Precision  #################################" << std::endl;
	print_sep_line(verbose, std::clog);

	if (verbose) std::clog << std::endl << "##type [domain, points]    errD(Norm)   errD(Dual)   Improv      errD2(Norm)   errD(Dual)    Improv" << std::endl;
	if (verbose) std::clog << std::endl << "###### Default Examples" << std::endl;
	check_error_precision<double>(verbose, 10, 12, "double_D10_N12");

	check_error_precision<double>(verbose, 10, 24, "double_D10_N24");
	check_error_precision<double>(verbose, 10, 36, "double_D10_N32");
	check_error_precision<double>(verbose, 10, 48, "double_D10_N48");

	std::vector<double> check_1060 = check_error_precision<double>(verbose, 10, 60, "double_D10_N60", true);


	if (verbose) std::clog <<  std::endl << "###### User Example" << std::endl;
	std::string str_start = (boost::format("double_D%d_N%d") % dom_number % dom_points ).str();
	check_error_precision<double>(verbose, dom_number, dom_points, str_start);

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
