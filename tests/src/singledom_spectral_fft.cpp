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


#include "tests/include/mpreal.h"
#include "aneto/spectral_domain.hpp"
#include "boost/format.hpp"
#include <chrono>
#include <iostream>
#include <fstream>

#include "eigen3/unsupported/Eigen/FFT"

#include <tests/include/auxiliars.hpp>



using mpfr::mpreal;
using namespace std;

template <typename T> int check_FFT(bool verbose, size_t spec_points, size_t dec_prec, string start_str, size_t iterations = 1, bool only_fft = false)
{
	T *function, *coeffs_mat;
	T *coeffs_fft;


	function   = new T [spec_points + 1];
	coeffs_mat = new T [spec_points + 1];
	coeffs_fft = new T [spec_points + 1];

	aneto::spectral_domain<T> spec_mat(spec_points, aneto::TRANSF_OPTION::MATRIX);
	aneto::spectral_domain<T> spec_fft(spec_points, aneto::TRANSF_OPTION::FFT);


	for (size_t i = 0; i <= spec_points; i++)
		function[i] = sin(spec_mat.get_X(i)) + cos(spec_mat.get_X(i));



	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	if (only_fft == false) {
		for (size_t i = 0; i < iterations; i++)
		spec_mat.get_spectral_coeffs(function, coeffs_mat);
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto dur_matrix = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / iterations;

	t1 = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < iterations; i++) {
		spec_fft.get_spectral_coeffs(function, coeffs_fft);
	}
	t2 = std::chrono::high_resolution_clock::now();
	auto dur_fft = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / iterations;

	if (verbose) {

		std::ofstream  outfile(start_str + "_fft.dat");

		for (size_t i = 0; i <= spec_points; i++) {

			outfile << i << "  ";
			outfile << toString(spec_mat.get_X(i)) << "  ";
			outfile << toString(function[i]) << "  ";
			outfile << toString(coeffs_mat[i]) << "  ";
			outfile << toString(coeffs_fft[i]) << "  ";
			outfile << std::endl;
		}
	}

	T max_error = compute_max_error(coeffs_mat, coeffs_fft, spec_points);


	if (only_fft == false) {
		if (verbose) std::clog << start_str << "\t\t" << toString_5e(max_error) << " \t\t"  << dur_matrix << "\t\t" << dur_fft << "\t\t";
	}
	else {
		if (verbose) std::clog << start_str << "\t\t" << "No_Comparison" << "\t\t"  << "\" \"  " << "\t\t" << dur_fft << "\n" ;
		return 0;
	}

	if (spec_points < 200) {
		if (max_error < pow((T)10.00, -(dec_prec-2))){
			if (verbose) std::clog << "Ok!" << std::endl;
			return 0;
		}
		else {
			std::clog << "Fail error!" << std::endl;
			return 2;
		}
	}else {
		if (max_error < pow((T)10.00, -(dec_prec-2)) && dur_matrix > dur_fft  ){
			if (verbose) std::clog << "Ok!" << std::endl;
			return 0;
		}
		else if (max_error > pow((T)10.00, -(dec_prec-2))){
			if (verbose) std::clog << "Fail error!" << std::endl;
			return 2;
		}else {
			if (verbose) std::clog << "Ok!   ";
			if (verbose) std::clog << "Time of the fft transf is bigger than the matrix transf. ";
			if (verbose) std::clog << "If N is a power of 2, check that because there is probably something wrong!" << std::endl;
			return 1;
		}

	}

}



int singledom_spectral_fft(bool verbose, size_t spec_points)
{
	string str;


	if (verbose)
	{
		print_sep_line(verbose, std::clog);
		if (verbose) std::clog << "#################### Test Singledomain:  Spectral Coefficients Matrix / FFT ########################" << std::endl;
		print_sep_line(verbose, std::clog);
		if (verbose) std::clog << std::endl << "# Double precision." << std::endl;
		if (verbose) std::clog << "points\t\terror    \t\tt_mat(ms)\tt_fft(ms)\tTest "<< std::endl;

		for (size_t i = 2; i <=6; i++) {
			size_t poi = pow(2, i);

			str = (boost::format("%04d") % poi ).str();

			check_FFT<double>(verbose, poi, 16, str, 100);
		}
		for (size_t i = 7; i <=12; i++) {
			size_t poi = pow(2, i);

			str = (boost::format("%04d") % poi ).str();

			check_FFT<double>(verbose, poi, 16, str, 1);
		}
		for (size_t i = 13; i <=14; i++) {
			size_t poi = pow(2, i);

			str = (boost::format("%04d") % poi ).str();

			check_FFT<double>(verbose, poi, 16, str, 1, true);
		}

		std::clog << std::endl;
		print_sep_line(verbose, std::clog);
		if (verbose) std::clog << "# MPFR - 300b precision." << std::endl;
		if (verbose) std::clog << "points\t\terror\t\tt_mat(ms)\t\tt_fft(ms)"<< std::endl;

		mpfr::mpreal::set_default_prec(300);

		for (size_t i = 2; i <=6; i++) {
			size_t poi = pow(2, i);

			str = (boost::format("%04d") % poi ).str();

			check_FFT<mpreal>(verbose, poi, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), str, 100);
		}
		for (size_t i = 7; i <=9; i++) {
			size_t poi = pow(2, i);

			str = (boost::format("%04d") % poi ).str();

			check_FFT<mpreal>(verbose, poi, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), str, 1);
		}
		for (size_t i = 10; i <=11; i++) {
			size_t poi = pow(2, i);

			str = (boost::format("%04d") % poi ).str();

			check_FFT<mpreal>(verbose, poi, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), str, 1, true);
		}


		if (verbose) std::clog << std::endl << "# Using a single domain of " << spec_points << " points." << std::endl;

	}

	int ERROR = 0;
	size_t poi = 1024;
	mpreal::set_default_prec(300);

	ERROR += check_FFT<double>(verbose, poi, 16, (boost::format("%04d") % poi ).str());
	ERROR += check_FFT<mpreal>(verbose, poi, mpfr::bits2digits(mpfr::mpreal::get_default_prec()), (boost::format("%04d") % poi ).str(), 1);

	if (ERROR == 0) {
		std::cout << "----- Test Singledomain Spect Coeffs FFT        ('singledom_spectral_fft'):          OK!\n";
		return 0;
	}
	else{
		std::cout << "----- Test Singledomain Spect Coeffs FFT        ('singledom_spectral_fft'):          Fail! Error / time of computation have failed.\n";
		return 1;
	}



}
