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


#ifndef _AUXILIARS_HPP
#define _AUXILIARS_HPP

#include <iostream>
#include <boost/multiprecision/mpfr.hpp>

#include "boost/format.hpp"
#include "mpreal.h"


static inline void print_sep_line(bool verbose, std::ostream &out)
{
	if (verbose) out << "#######################################################################################################" << std::endl;
}


template<typename T> T cabs(T x)
{
	if (x > 0)
		return x;
	else
		return -x;
}

template<typename T> T compute_max_error(T *num, T *theo, size_t points, size_t first = 0)
{
	T temp;
	T max_err = 0;

	for(size_t i = first; i < points; i++) {
		temp = cabs(num[i] - theo[i]);

		if (temp > max_err)
			max_err = temp;
	}
	return cabs(max_err);
}


static inline std::string toString(float number)
{
	return (boost::format("%.8e") % number).str();
}

static inline std::string toString(double number)
{
	return (boost::format("%.16e") % number).str();
}

static inline std::string toString(long double number)
{
	return (boost::format("%.32e") % number).str();
}

static inline std::string toString(mpfr::mpreal number)
{
	return number.toString();
}

static inline std::string toString_5e(float number)
{
	return (boost::format("%+.5e") % number).str();
}

static inline std::string toString_5e(double number)
{
	return (boost::format("%+.5e") % number).str();
}

static inline std::string toString_5e(long double number)
{
	return (boost::format("%+.5e") % number).str();
}

static inline std::string toString_5e(mpfr::mpreal number)
{
	return number.toString("%+.5Re");
}
static inline std::string toString_5e(boost::multiprecision::mpf_float number)
{
	return number.str();
}


static inline std::string toString_f(float number)
{
	return (boost::format("%+.3f") % number).str();
}
static inline std::string toString_f(double number)
{
	return (boost::format("%+.3f") % number).str();
}
static inline std::string toString_f(long double number)
{
	return (boost::format("%+.3f") % number).str();
}
static inline std::string toString_f(mpfr::mpreal number)
{
	return number.toString("%+.3Rf");
}

static inline std::string toString_f(boost::multiprecision::mpf_float number)
{
	return number.str();
}





#if QUAD == 1
	#include <boost/multiprecision/float128.hpp>

	static inline std::string toString(boost::multiprecision::float128 number)
	{
		return (boost::format("%.32e") % number).str();
	}

	static inline std::string toString_5e(boost::multiprecision::float128 number)
	{
		return (boost::format("%.5e") % number).str();
	}

#endif



#endif
