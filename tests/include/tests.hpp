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


#ifndef _TESTS_HPP
#define _TESTS_HPP


int multidom_derivative_dual(bool verbose, size_t dom_number = 10, size_t dom_points = 47);
int multidom_derivative_precision(bool verbose, size_t dom_number = 14, size_t dom_points = 63);
int multidom_integral_comparison(bool verbose, size_t spec_domains = 14, size_t spec_points = 63);
int multidom_integral_LR(bool verbose, size_t spec_domains = 14, size_t spec_points = 63);
int multidom_grids(bool verbose, size_t spec_domains = 14, size_t spec_points = 63);
int root_finding (bool verbose, size_t user_spec_points, size_t user_domains, size_t user_dom_points, size_t user_dec_prec);
int singledom_derivatives(bool verbose, size_t spec_points);
int singledom_integral_LR(bool verbose, size_t spec_points);
int singledom_integral_comparison(bool verbose, size_t spec_points);
int singledom_spectral_fft(bool verbose, size_t spec_points);
int interpolation (bool verbose, size_t user_spec_points, size_t user_domains, size_t user_dom_points, size_t user_dec_prec);

#endif /* _TESTS_HPP */
