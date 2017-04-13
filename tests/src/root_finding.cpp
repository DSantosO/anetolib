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


template <typename T> int check_root_finder(bool verbose, size_t spec_points, size_t dom_num, size_t dom_points, size_t dec_prec, std::string str_start)
{
	int exp_error =  -(dec_prec-2);

	aneto::spectral_domain<T> spec(spec_points);
	aneto::multidomain<T> multidom(dom_num, dom_points, -1.00, +1.00);


	T *r_spec = new T [spec_points + 1];
	T *r_mult = new T [multidom.number_points()];
	T value, sol, sol_mul, sol_spec;

	for (size_t i = 0; i <= spec_points; i++)
		r_spec[i] = - log(.5*(spec.get_X(i) + 1.5));

	for (size_t i = 0; i < multidom.number_points(); i++)
		r_mult[i] = - log(.5*(multidom.xp(i) + 1.5));


	value = 0.80;
	sol = ((T)2.00) * exp(-value) - (T)(1.5000);
	sol_mul = multidom.root_finder_decreasing(value, r_mult);
	sol_spec = spec.root_finder(value, r_spec);

	if (verbose) std::clog << str_start           << "   ";
	if (verbose) std::clog << toString_5e(sol)       << "   ";
	if (verbose) std::clog << toString_5e(sol_spec)  << "   ";
	if (verbose) std::clog << toString_5e(sol_mul)   << "   ";
	if (verbose) std::clog << toString_5e(cabs(sol-sol_spec)) << "   ";
	if (verbose) std::clog << toString_5e(cabs(sol-sol_mul))  << "   ";

	delete[] r_spec;
	delete[] r_mult;

	if (cabs(sol-sol_spec) < pow((T)10.00, exp_error) && cabs(sol-sol_mul) < pow((T)10.00, exp_error)) {
		if (verbose) std::clog << "OK!" << std::endl;
		return 0;
	}
	else {
		if (verbose) std::clog << "Fail!" << std::endl;
		return 1;
	}

}




int root_finding (bool verbose, size_t user_spec_points, size_t user_domains, size_t user_dom_points, size_t user_dec_prec)
{
	size_t spec_points = 340;
	size_t domains = 100;
	size_t dom_points = 64;
	size_t dec_prec = 100;
	std::string str_start;
	size_t TOT_ERROR = 0;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "############################# Test Multidomain(m)/Single(s):  Root Finder #############################" << std::endl;
	print_sep_line(verbose, std::clog);
// 	if (verbose) std::clog
	if (verbose) std::clog << "#type   bits    decim    N(s)    D(m)    N(m)     sol(theo)       sol(s)         sol (m)       err(s)         err(m)  \n";


	str_start = (boost::format("mpfr:   %04ib    %04i    %04i    %04i    %04i") \
	            % mpfr_get_default_prec() % dec_prec % spec_points % domains % dom_points).str();
	TOT_ERROR += check_root_finder<mpfr::mpreal>(verbose, spec_points, domains, dom_points, dec_prec, str_start);

	str_start = (boost::format("mpfr:   %04ib    %04i    %04i    %04i    %04i") \
	            % mpfr_get_default_prec() % dec_prec % spec_points % domains % dom_points).str();
	TOT_ERROR += check_root_finder<mpfr::mpreal>(verbose, spec_points, domains, dom_points, dec_prec, str_start);

	str_start = (boost::format("doub:   %04ib    %04i    %04i    %04i    %04i") \
	            % 53 % 16 % spec_points % domains % dom_points).str();
	TOT_ERROR += check_root_finder<double>(verbose, spec_points, domains, dom_points, 16, str_start);


	print_sep_line(verbose, std::clog);
	print_sep_line(verbose, std::clog);
	if (TOT_ERROR == 0) {
		std::cout << "----- Test Root Finding                         ('root_finding'):                    OK!\n";
	}
	else {
		std::cout << "----- Test Root Finding                         ('root_finding'):                    Failed!  Some of the roots finded have greater error than expected.\n";
	}
	print_sep_line(verbose, std::clog);

	if(verbose == false)
		return TOT_ERROR;


	/* User defined tests */
	spec_points = user_spec_points;
	domains = user_domains;
	dom_points = user_dom_points;
	dec_prec = user_dec_prec;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(dec_prec));

	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << "#User defined tests\n";
	if (verbose) std::clog << "#type   bits    decim    N(s)    D(m)    N(m)     sol(theo)       sol(s)         sol (m)       err(s)         err(m)  \n";
	str_start = (boost::format("mpfr:   %04ib    %04i    %04i    %04i    %04i") \
	            % mpfr_get_default_prec() % dec_prec % spec_points % domains % dom_points).str();
	check_root_finder<mpfr::mpreal>(verbose, spec_points, domains, dom_points, dec_prec, str_start);

	str_start = (boost::format("doub:   %04ib    %04i    %04i    %04i    %04i") \
	            % 53 % 16 % spec_points % domains % dom_points).str();
	check_root_finder<double>(verbose, spec_points, domains, dom_points, 16, str_start);
	print_sep_line(verbose, std::clog);
	if (verbose) std::clog << std::endl << std::endl;




	return TOT_ERROR;







}
