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
#include <chrono>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "tests/include/mpreal.h"

#include "aneto/spectral_domain.hpp"
#include "aneto/multidomain.hpp"
#include "tests/include/tests.hpp"


namespace po = boost::program_options;

int main(int argc, char **argv)
{
	size_t dom_number;
	size_t dom_points;
	size_t spec_points;
	size_t dec_prec;
	std::string test;
	bool verbose;
	bool files;
	size_t TOT_ERROR = 0;


	// Declare the supported options.
	po::options_description desc("ANETO lib Tests");
	desc.add_options()
		("help,h", "Show this help")
		("verbose,v", "Print all the information about the tests.")
 		("test,t", po::value<std::string>(&test)->default_value("all"), "Type of test to execute")
		("spec_point,S", po::value<size_t>(&spec_points)->default_value(240), "Number of Points for the single domain tests")
		("dom_number,D", po::value<size_t>(&dom_number)->default_value(14), "Number of Domains")
		("dom_points,N", po::value<size_t>(&dom_points)->default_value(63), "Number of Points per domain")
		("dec_prec,p", po::value<size_t>(&dec_prec)->default_value(100), "Significant digits for the tests that admits it.")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		exit(1);
	}

	if (vm.count("verbose"))
		verbose = true;
	else
		verbose = false;

	if (test == "all") {
		TOT_ERROR += multidom_derivative_dual(verbose, dom_number, dom_points);
		TOT_ERROR += multidom_integral_comparison(verbose, dom_number, dom_points);
		TOT_ERROR += multidom_integral_LR(verbose, dom_number, dom_points);
		TOT_ERROR += singledom_derivatives(verbose, spec_points);
		TOT_ERROR += singledom_integral_comparison(verbose, spec_points);
		TOT_ERROR += singledom_integral_LR(verbose, spec_points);
		TOT_ERROR += singledom_spectral_fft(verbose, spec_points);
		TOT_ERROR += root_finding(verbose, spec_points, dom_number, dom_points, dec_prec);
		TOT_ERROR += interpolation(verbose, spec_points, dom_number, dom_points, dec_prec);
	}
	else if (test == "multidom") {
		TOT_ERROR += multidom_derivative_dual(verbose, dom_number, dom_points);
		TOT_ERROR += multidom_integral_comparison(verbose, dom_number, dom_points);
		TOT_ERROR += multidom_integral_LR(verbose, dom_number, dom_points);
		TOT_ERROR += root_finding(verbose, spec_points, dom_number, dom_points, dec_prec);
		TOT_ERROR += interpolation(verbose, spec_points, dom_number, dom_points, dec_prec);
	}
	else if (test == "singledom") {
		TOT_ERROR += singledom_derivatives(verbose, spec_points);
		TOT_ERROR += singledom_integral_comparison(verbose, spec_points);
		TOT_ERROR += singledom_integral_LR(verbose, spec_points);
		TOT_ERROR += singledom_spectral_fft(verbose, spec_points);
		TOT_ERROR += root_finding(verbose, spec_points, dom_number, dom_points, dec_prec);
		TOT_ERROR += interpolation(verbose, spec_points, dom_number, dom_points, dec_prec);
	}
	else if (test == "multidom_integral_LR") {
		TOT_ERROR += multidom_integral_LR(verbose, dom_number, dom_points);
	}
	else if (test == "multidom_integral_comparison") {
		TOT_ERROR += multidom_integral_comparison(verbose, dom_number, dom_points);
	}
	else if (test == "multidom_derivative_dual") {
		TOT_ERROR += multidom_derivative_dual(verbose, dom_number, dom_points);
	}
	else if (test == "multidom_derivative_precision") {
		TOT_ERROR += multidom_derivative_precision(verbose, dom_number, dom_points);
	}
	else if (test == "multidom_grids") {
		TOT_ERROR += multidom_grids(verbose, dom_number, dom_points);
	}
	else if (test == "root_finding") {
		TOT_ERROR += root_finding(verbose, spec_points, dom_number, dom_points, dec_prec);
	}
	else if (test == "singledom_derivatives"){
		TOT_ERROR += singledom_derivatives(verbose, spec_points);
	}
	else if (test == "singledom_integral_LR"){
		TOT_ERROR += singledom_integral_LR(verbose, spec_points);
	}
	else if (test == "singledom_integral_comparison"){
		TOT_ERROR += singledom_integral_comparison(verbose, spec_points);
	}
	else if (test == "singledom_spectral_fft"){
		TOT_ERROR += singledom_spectral_fft(verbose, spec_points);
	}
	else if (test == "interpolation"){
		TOT_ERROR += interpolation(verbose, spec_points, dom_number, dom_points, dec_prec);
	}
	else {
		std::clog << "Error, option " << test << " not valid\n";
		std::exit(1);
	}

	if (TOT_ERROR == 0) {
		std::clog << std::endl <<  "All tests are satisfactory! Congratulations!" << std::endl;
	}
	else {
		std::clog << std::endl;
		std::clog << "Some of the test have failed." << std::endl;
		std::clog << "You can rerun the specific test with '-t name_test' and using '-v' to verbose all the information." << std::endl;
	}


}
