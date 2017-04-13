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

#include "tests/include/mpreal.h"
#include "tests/include/auxiliars.hpp"

#include "aneto/spectral_domain.hpp"
#include "aneto/multidomain.hpp"


template <typename T> int check_interpolation(bool verbose, size_t spec_points, size_t dom_num, size_t dom_points, size_t dec_prec, std::string str_start)
{
	int exp_error_f =  -(dec_prec-2);
	int exp_error_d =  -(dec_prec-4);
	int exp_error_d2 =  -(dec_prec-7);

	aneto::spectral_domain<T> spec(spec_points);
	aneto::multidomain<T> multidom(dom_num, dom_points, -1.00, +1.00);


	T *r_spec = new T [spec_points + 1];
	T *r_mult = new T [multidom.number_points()];
	T *spec_coeffs = new T [spec_points + 1];

	T x_interp, fun_x, der_x, dv2_x;
	T fun_mul, fun_mul2;
	T fun_spec, der_spec, dv2_spec;
	T fun_coeff_spec;

	for (size_t i = 0; i <= spec_points; i++)
		r_spec[i] = - log( ((T)0.5) * (spec.get_X(i) + (T)1.5));

	for (size_t i = 0; i < multidom.number_points(); i++)
		r_mult[i] = - log( ((T)0.5) * (multidom.xp(i) + (T)1.5));


	x_interp = 0.55;

	fun_x = - log( ((T)0.5) * (x_interp + (T)1.5));
	der_x = - (T)0.5 / ( (T)0.5 * x_interp + (T)0.75);
	dv2_x = + (T)0.25 / ( (T)0.5 * x_interp + (T)0.75) / ( (T)0.5 * x_interp + (T)0.75);


	spec.compute_spectral_func_der(r_spec, aneto::JAC_OPTION::SPECTRAL, NULL, aneto::COEFF_OPTION::D2);
	fun_spec = spec.get_fun(x_interp);
	der_spec = spec.get_der(x_interp);
	dv2_spec = spec.get_d2(x_interp);

	spec.get_spectral_coeffs(r_spec, spec_coeffs);
	fun_coeff_spec = spec.interpolate_with_spec_coeffs(x_interp, spec_coeffs);

	fun_mul =  multidom.interpolate(r_mult, x_interp, false);
	fun_mul2 = multidom.interpolate(r_mult, x_interp, true);



	if (verbose) std::clog << str_start           << "   ";
	if (verbose) std::clog << toString_5e(cabs(fun_x-fun_spec))  << "   ";
	if (verbose) std::clog << toString_5e(cabs(fun_x-fun_coeff_spec))   << "   ";
	if (verbose) std::clog << toString_5e(cabs(fun_x-fun_mul))   << "   ";
	if (verbose) std::clog << toString_5e(cabs(fun_x-fun_mul2))   << "   ";

	if (verbose) std::clog << toString_5e(cabs(der_spec - der_x))  << "   ";
	if (verbose) std::clog << toString_5e(cabs(dv2_spec - dv2_x))  << "   ";

	T error_fun = (cabs(fun_x-fun_spec) + cabs(fun_x-fun_coeff_spec) +  \
				 cabs(fun_x-fun_mul) + cabs(fun_x-fun_mul2) ) / 4.;


	delete[] r_spec;
	delete[] r_mult;

	if (error_fun < pow((T)10.00, exp_error_f) && cabs(der_spec - der_x) < pow((T)10.00, exp_error_d) && cabs(dv2_spec - dv2_x) < pow((T)10.00, exp_error_d2)) {
		if (verbose) std::clog << "OK!" << std::endl;
		return 0;
	}
	else {
		if (verbose) std::clog << "Fail!" << std::endl;
		return 1;
	}

}




int interpolation (bool verbose, size_t user_spec_points, size_t user_domains, size_t user_dom_points, size_t user_dec_prec)
{

	size_t spec_points = 340;
	size_t domains = 100;
	size_t dom_points = 64;
	size_t dec_prec = 100;
	std::string str_start;
	size_t TOT_ERROR = 0;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "############################# Test Multidomain(m)/Single(s):  Interpolation #############################" << std::endl;
	print_sep_line(verbose, std::clog);

	if (verbose) std::clog << "#type   bits    decim    N(s)    D(m)    N(m)     f_spec_nor     f_spec_coeff     f_mul_dir       f_mul_stor     f_spec_der     f_spec_der2\n";


	str_start = (boost::format("mpfr:   %04ib    %04i    %04i    %04i    %04i") \
	            % mpfr_get_default_prec() % dec_prec % spec_points % domains % dom_points).str();
	TOT_ERROR += check_interpolation<mpfr::mpreal>(verbose, spec_points, domains, dom_points, dec_prec, str_start);

	str_start = (boost::format("mpfr:   %04ib    %04i    %04i    %04i    %04i") \
	            % mpfr_get_default_prec() % dec_prec % spec_points % domains % dom_points).str();
	TOT_ERROR += check_interpolation<mpfr::mpreal>(verbose, spec_points, domains, dom_points, dec_prec, str_start);

	str_start = (boost::format("doub:   %04ib    %04i    %04i    %04i    %04i") \
	            % 53 % 16 % spec_points % domains % dom_points).str();
	TOT_ERROR += check_interpolation<double>(verbose, spec_points, domains, dom_points, 16, str_start);


	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);

	if (TOT_ERROR == 0) {
		std::cout << "----- Test Interpolation                        ('interpolation'):                   OK!\n";
	}
	else {
		std::cout << "----- Test Interpolation                        ('interpolation'):                   Failed!  Some of the interpolation finded have greater error than expected.\n";
	}

	print_sep_line(verbose, std::clog);

	if(verbose == false)
		return TOT_ERROR;

	/* User defined tests */
// 	spec_points = user_spec_points;
// 	domains = user_domains;
// 	dom_points = user_dom_points;
// 	dec_prec = user_dec_prec;

// 	mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));
//
// 	print_sep_line(verbose, std::clog);
// 	if (verbose) std::clog << "#User defined tests\n";
// 	if (verbose) std::clog << "#type   bits    decim    N(s)    D(m)    N(m)     sol(theo)       sol(s)         sol (m)       err(s)         err(m)  \n";
// 	str_start = (boost::format("mpfr:   %04ib    %04i    %04i    %04i    %04i") \
// 	            % mpfr_get_default_prec() % dec_prec % spec_points % domains % dom_points).str();
// 	check_root_finder<mpfr::mpreal>(verbose, spec_points, domains, dom_points, dec_prec, str_start);
//
// 	str_start = (boost::format("doub:   %04ib    %04i    %04i    %04i    %04i") \
// 	            % 53 % 16 % spec_points % domains % dom_points).str();
// 	check_root_finder<double>(verbose, spec_points, domains, dom_points, 16, str_start);
// 	print_sep_line(verbose, std::clog);
// 	if (verbose) std::clog << std::endl << std::endl;



	return TOT_ERROR;







}
