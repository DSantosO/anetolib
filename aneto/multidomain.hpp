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


#ifndef _MULTIDOMAIN_HPP
#define _MULTIDOMAIN_HPP

#include <iostream>

#include "aneto/spectral_domain.hpp"


#ifndef OMP
	#define OMP 0
#endif

#if OMP == 1
	#include <omp.h>
// #include "extrae_user_events.h"
#endif

namespace aneto {


	/**
	 * Mode that select the the position of the internal domains of the grid.
	 */
	enum class GRID_OPTION
	{
		UNIFORM,  	 ///< All domains have the same size.
		DEF_BNDS  	 ///< User defined boundaries for the domains.
	};

	/**
	 * Derivatives can be computed using one set of pseudospectral domains
	 * or two of them. In the latest, the domains are centred in the nodes
	 * of the original domain. The final derivatives are computed using a mixed
	 * of both grids that reduces notably the error near the boundaries of the
	 * domains.
	 */
	enum class DERGRID_OPTION
	{
		SINGLE,		///< Derivatives are computed in one domain.
		DUAL		///< Derivatives are computed in two overlapping domains.
		            ///<     This reduces the error near the domain boundaries
	};





	/** \brief Template of the multidomain structure defined a set of spectral_domain.
	 *
	 * Class that create a multidomain structure with all the individual
	 * domains being defined with the spectral_domain class. \n
	 *
	 * The multidomain can be created directly with all the options or just created
	 * empty and be initialised afterwards. Once it is initialise the structure
	 * can not be modify but can be cleaned up and initialised again.\n
	 *
	 * The size of the grid needs to be specify at the moment of initialisation. \n
	 * The position of the internal domains can be specified as an option that
	 * can take the value GRID_OPTION::UNIFORM for constant size domains and also
	 * GRID_OPTION::DEF_BNDS if we want to specify the position of the boundaries.
	 *
	 * In the initialisation also can be choosen the mode of performing derivatives.
	 * DERGRID_OPTION::SINGLE will compute normal derivatives and DERGRID_OPTION::DUAL
	 * will create a dual grid structure that can increase the accuracy of the derivatives
	 * near the domain boundaries in two orders of magnitude. \n
	 *
	 * The number_threads option for select to do the operations in parallel is not yet
	 * supported and only admit the default value of one.
	 *
	 * \n\n
	 * The template can be used with any general type where
	 * the usual operations + - * / are defined.\n
	 * Also a cos / acos functions are required for the specific type.\n
	 * \tparam T type to be used for the spectral grid
	 */
	template <typename T> class multidomain
	{
	private:
		bool m_allocated;
		bool m_der_dual;


		/* OpenMP Variables */
		size_t m_num_threads;


		T *m_offset;

		/* Multidomain variables */
		T *m_main_xp;							///< Array that stores the physical coordinate.
		T *m_main_dxp_dX;						///< Array that stores the Jacobian
		T *m_main_dX_dxp;						///< Array that stores the Jacobian (inverse)

		/* Dual Grid variables */
		T *m_dual_xp;						///< Array that stores the physical coordinate of the dual multidomain.
		T *m_dual_dxp_dX;					///< Array that stores the Jacobian of the dual multidomain.
		T *m_dual_dX_dxp;					///< Array that stores the Jacobian (inverse) of the dual multidomain.
		T *m_partit;						///< Array that stores the weights of main/dual grid for computing the derivative
										///< Dt(xp) = partit(x) * Dmain (xp) + (1 + partit(x)) * Ddual (xp)

		spectral_domain<T> *m_dualDomain;		///< Grid for the single dual domain.

		T m_zero;
		T m_one;
		T m_two;
		T m_half;

		size_t m_total_points;
				// total points
		size_t m_domain_points;						// points per domain
		size_t m_domain_number;						// total domains

		T *m_dom_nodes;
		T *m_aux1;
		T *m_aux2;

		spectral_domain<T> *m_lobGrid;
	///< Grid for a single domain.
		size_t *m_dom_beg;
		size_t *m_dom_end;

	public:


		/** Void constructor.
		 *
		 * A multidomain is created empty.
		 * It needs to be initialised in order to be used.
		 */
		multidomain() : m_allocated(false), m_total_points(0), m_domain_points(0), m_domain_number(0) {}

		/** Constructor.
		 *
		 * It creates a multidomain with D domains and N+1 grid points each.
		 *
		 * \param p_domain_number Number of domains
		 * \param p_domain_points Number of points in each domain (p_domain_points + 1)
		 * \param left_bndry      left global boundary. Default 0.
		 * \param right_bndry     right global boundary. Default +1.
		 * \param number_threads  Number of OpenMP threads. Warning! now, only supports one.
		 * \param dergrid_option  Option for derivatives.
		 *                        DERGRID_OPTION::SINGLE for normal multidomain.
		 *                        DERGRID_OPTION::DUAL   for a dual grid.
		 * \param grid_option     Option for the grid structure.
		 *                        GRID_OPTION::UNIFORM:  Equally distributed domains
		 *                        GRID_OPTION::DEF_BNDS: Internal boundaries can be defined by the user.
		 * \param boundaries      Array with the boundaries of the domains.
		 *                        It suppose to have D-1 components starting with the
		 *                        boundary between the first and second domains, not including the global boundaries.
		 *                        If the boundaries are not ordered correctly the run will be aborted.
		 *                        Can be left empty or passing a NULL points if it is not needed.
		 */
		multidomain(size_t p_domain_number, size_t p_domain_points, T left_bndry = 0,
					T right_bndry = +1, size_t number_threads = +1,
					DERGRID_OPTION dergrid_option = DERGRID_OPTION::SINGLE,
					GRID_OPTION grid_option = GRID_OPTION::UNIFORM, T *boundaries = NULL) {

			initialise(p_domain_number, p_domain_points, left_bndry, right_bndry, number_threads, \
			           dergrid_option, grid_option, boundaries);
		}

		/**
		 * It initialises the multidomain with D domains and N grid points each.
		 * The exact range of the coordinate can be chosen.
		 *
		 * \param p_domain_number Number of domains
		 * \param p_domain_points Number of points in each domain (p_domain_points + 1)
		 * \param left_bndry      left global boundary. Default 0.
		 * \param right_bndry     right global boundary. Default +1.
		 * \param number_threads  Number of OpenMP threads. Warning! now, only supports one.
		 * \param dergrid_option  Option for derivatives.
		 *                        DERGRID_OPTION::SINGLE for normal multidomain.
		 *                        DERGRID_OPTION::DUAL   for a dual grid.
		 * \param grid_option     Option for the grid structure.
		 *                        GRID_OPTION::UNIFORM:  Equally distributed domains
		 *                        GRID_OPTION::DEF_BNDS: Internal boundaries can be defined by the user.
		 * \param boundaries      Array with the boundaries of the domains.
		 *                        It suppose to have D-1 components starting with the
		 *                        boundary between the first and second domains, not including the global boundaries.
		 *                        If the boundaries are not ordered correctly the run will be aborted.
		 *                        Can be left empty or passing a NULL points if it is not needed.
		 */
		void initialise(size_t p_domain_number, size_t p_domain_points, T left_bndry = 0,
						T right_bndry = +1, size_t number_threads = +1,
						DERGRID_OPTION dergrid_option = DERGRID_OPTION::SINGLE,
						GRID_OPTION grid_option = GRID_OPTION::UNIFORM, T *boundaries = NULL) {

			if (grid_option == GRID_OPTION::DEF_BNDS && boundaries == NULL) {
				std::clog << "The bndries have to be given with the GRID_OPTION::DEF_BNDS. Check that, please." << std::endl;
				abort();
			}

			/* Checking entered values */
			if (p_domain_number > 0) {
				m_domain_number = p_domain_number;
			}
			else {
				std::clog << "Please introduced a number of domains greater than zero." << std::endl;
				abort();
			}

			if (p_domain_points >= 3) {
				m_domain_points = p_domain_points;
			}
			else {
				std::clog << "Please introduced a number of points per domain greater than three." << std::endl;
				abort();
			}
			if (number_threads != 1) {
// 				std::clog << "Warning!! OpenMP implementation is not fully supported yet. Please use carefully." << std::endl;
// 				abort();
			}

			/* Setting other internal variables */
			m_total_points = m_domain_number * (m_domain_points + 1);
			m_der_dual = (dergrid_option==DERGRID_OPTION::DUAL) ? true : false;




			_configure_omp(number_threads);

			_allocate_memory();
			_set_constants();

			m_allocated = true;

			if (grid_option == GRID_OPTION::UNIFORM)
				_set_grid_uniform(left_bndry, right_bndry);
			else if (grid_option == GRID_OPTION::DEF_BNDS)
				_set_grid_def_bnds(left_bndry, right_bndry, boundaries);
			else {
				std::clog << "The grid option is not valid. Check that, please." << std::endl;
				abort();
			}



			if (dergrid_option == DERGRID_OPTION::DUAL)
				_set_dual_grid();

			if (grid_option != GRID_OPTION::UNIFORM && dergrid_option == DERGRID_OPTION::DUAL){
				std::clog << "Check comp_derivative with dual grid. It is going to fail. " << std::endl;
				abort();
			}

		}

		/**
		 * Returns the total number of points.\n
		 * total = domain_number * (domain_points + 1);
		 *
		 */
		inline size_t number_points() const
		{
			return m_total_points;
		}

		/**
		 * Returns the number of points per domain.
		 *
		 */
		inline size_t dom_points() const
		{
			return m_domain_points;
		}

		/**
		 * Returns the number of domains.
		 *
		 */
		inline size_t dom_number() const
		{
			return m_domain_number;
		}

		/**
		 * Returns the value of the multidomain coordinate in the collocation point i.
		 * \param i the global index
		 */
		inline T xp(size_t i) const
		{
			return m_main_xp[i];
		}

		/**
		 * Returns the value of the multidomain coordinate in the collocation point i.
		 * \param dom the domain
		 * \param k the index in the domain
		 */
		inline T xp(size_t dom, size_t k) const
		{
			return m_main_xp[get_id(dom,k)];
		}

		/**
		 * Returns the value of the multidomain jacobian in the point k fo the domain dom.
		 * \param i the global index
		 */
		inline T dxp_dX(size_t i) const
		{
			return m_main_dxp_dX[i];
		}

		/**
		 * Returns the value of the multidomain jacobian in the point k fo the domain dom.
		 * \param dom the domain
		 * \param k the index in the domain
		 */
		inline T dxp_dX(size_t dom, size_t k) const
		{
			return m_main_dxp_dX[get_id(dom,k)];
		}

		/**
		 * Returns the value of the multidomain jacobian (inverse) in the point k fo the domain dom.
		 * \param i the global index
		 */
		inline T dX_dxp(size_t i) const
		{
			return m_main_dX_dxp[i];
		}

		/**
		 * Returns the value of the multidomain jacobian (inverse) in the point k fo the domain dom.
		 * \param dom the domain
		 * \param k the index in the domain
		 */
		inline T dX_dxp(size_t dom, size_t k) const
		{
			return m_main_dX_dxp[get_id(dom,k)];
		}


		/**
		 * Returns the left boundary of a domain.
		 * \param dom the domain
		 */
		inline T get_l_bndry(size_t dom) const
		{
			return m_dom_nodes[dom];
		}

		/**
		 * Returns the right boundary of a domain.
		 * \param dom the domain
		 *
		 */
		inline T get_r_bndry(size_t dom)
		{
			return m_dom_nodes[dom+1];
		}

		/**
		 * Returns the left global boundary.
		 */
		inline T get_l_glob_bndry()
		{
			return m_dom_nodes[0];
		}

		/**
		 * Returns the right global boundary.
		 */
		inline T get_r_glob_bndry()
		{
			return m_dom_nodes[m_domain_number];
		}

		/**
		 * Returns the index array of the point i of the domain a.\n
		 * The input is not checked.
		 * \param dom  the domain
		 * \param poin the point
		 */
		inline size_t get_id (const size_t dom, const size_t poin) const
		{
			return dom * (m_domain_points + 1) + poin;
		}

		/**
		 * Returns the domain that correspond with a certain global index.\n
		 * The input is not checked.
		 * \param id index of the array.
		 */
		inline size_t get_dom(const size_t id) const
		{
			return id / (m_domain_points + 1);
		}

		/**
		 * Returns number of point in the spectral grid that correspond with a certain global index.\n
		 * The input is not checked.
		 * \param id index of the array.
		 */
		inline size_t get_point(const size_t id) const
		{
			return id - get_dom(id) * (m_domain_points + 1);
		}

		/**
		 * Returns the index array of the point i of the domain a.\n
		 * The input is not checked.
		 * \param dom  the domain
		 * \param poin the point
		 */
		inline T* get_pointer_dom (T *pointer_a0, size_t dom)
		{
			return &(pointer_a0[get_id(dom, 0)]);
		}

		/**
		 * Returns the domain to which correspond a given coordinate.\n
		 * The method will fail if the given coordinate is outside the grid.
		 * \param x_i coordinate of the grid.
		 */
		inline size_t get_dom_coord(T x_i) const
		{
			size_t i_minu, i_plus, i_midd;

			i_minu = 0;
			i_plus = m_domain_number;

			// check v is in range
			if (x_i < m_dom_nodes[i_minu]) {
				std::clog << "Error!!  v = " <<  x_i << " is below the computing range (x_L: " << m_dom_nodes[i_minu] << ")" << std::endl;
				abort();
			}
			if (x_i > m_dom_nodes[i_plus]) {
				std::clog << "Error!!  v = " <<  x_i << " is above the computing range (x_R: " << m_dom_nodes[i_plus] << ")" << std::endl;
				abort();
			}

			while(true) {
				i_midd = (i_plus + i_minu)/2;


				if ( (i_plus - i_minu) == 1)
					return i_minu;
				else if (x_i < m_dom_nodes[i_midd])
					i_plus = i_midd;
				else
					i_minu = i_midd;
			}
		}

		/**
		 * Return the spectral coordinate in the the domain that the given coordinate is in.\n
		 * The method will fail if the given coordinate is outside the grid.
		 * \param x_i coordinate of the grid.
		 */
		inline T get_X(T x_i) {
			size_t dom = get_dom_coord(x_i);

			return (m_two * x_i - (m_dom_nodes[dom] + m_dom_nodes[dom+1]) )/(m_dom_nodes[dom+1] - m_dom_nodes[dom]);
		}

		/**
		 * Grid coordinate that corresponds with the spectral coordinate X of the domain dom.\n
		 * \param dom Domain we are interested in.
		 * \param X   Spectral coordinate.
		 */
		inline T get_xp(size_t dom, T X) {
			return m_half * ((m_dom_nodes[dom] + m_dom_nodes[dom+1]) + (m_dom_nodes[dom+1] - m_dom_nodes[dom]) * X);
		}

		/**
		 * Compute the value of a function in an arbitrary function thought spectral interpolation.\n
		 * Careful with the use_spec_stored option because the class don't control if the
		 * stored values correspond with the same function or domain of the current call.
		 *
		 * \param function         Function we want to interpolate.
		 * \param x_value          Coordinate where we want the function value.
		 * \param use_spec_stored  If activated, the interpolation will use the stored values
		 *                         from a previous interpolation.
		 *
		 */
		T interpolate(T *function, T x_value, bool use_spec_stored = false)
		{
			size_t dom = get_dom_coord(x_value);
			T spec_X = get_X(x_value);

			if (use_spec_stored)
				return m_lobGrid[0].get_fun(spec_X);
			else
				return m_lobGrid[0].interpolate(spec_X, get_pointer_dom(function, dom));
		}


		/**
		 * It find the coordinate where the function f(x) have a given value
		 * in a given domain.
		 * \f{eqnarray*}{
		 *     f(f_xr) - \mathrm{f_xr} = 0;
		 * \f}
		 *
		 * It requires that the function f(x) is monotonically decreasing.
		 * In this case the root, if exists is unique.
		 * If there is no root in the whole multidomain, the function will rise an error.
		 *
		 * \param f_xr Value of the function we
		 * \param func Value of the function in the collocation points.
		 * \return
		 */
		T root_finder_increasing(T f_xr, T *func)
		{
			size_t dom_value;

			/* Check if value is in a valid range */
			if (f_xr < func[0]) {
				std::clog << "Error!!  func=" << f_xr  <<  " is not in the computing range (0)" << std::endl;
				abort();
			}
			if (f_xr > func[m_total_points - 1]) {
				std::clog << "Error!!  func=" << f_xr  <<  " is not in the computing range (n)" << std::endl;
				abort();
			}

			/* Get f = value domain */
			size_t i_minu, i_plus, i_midd;
			i_minu = 0;
			i_plus = m_domain_number;

			if (f_xr > func[get_id(i_plus, 0)])
				dom_value = i_plus;

			while(true) {
				i_midd = (i_plus + i_minu)/2;

				if ( (i_plus - i_minu) == 1) {
					dom_value = i_minu;
					break;
				}
				else if (f_xr < func[get_id(i_midd,0)])
					i_plus = i_midd;
				else
					i_minu = i_midd;
			}

			T X_root = m_lobGrid[0].root_finder(f_xr, get_pointer_dom(func, dom_value));

			return get_xp(dom_value, X_root);

		}

		/**
		 * It find the coordinate where the function f(x) have a given value.
		 * \f{eqnarray*}{
		 *     f(f_xr) - \mathrm{f_xr} = 0;
		 * \f}
		 *
		 * It requires that the function f(x) is monotonically decreasing.
		 * In this case the root, if exists is unique.
		 * If there is no root in the whole multidomain, the function will rise an error.
		 *
		 * \param f_xr Value of the function we
		 * \param func Value of the function in the collocation points.
		 * \return
		 */
		T root_finder_decreasing(T f_xr, T *func)
		{
			size_t dom_value;

			/* Check if value is in a valid range */
			if (f_xr > func[0]) {
				std::clog << "Error!!  func=" << f_xr  <<  " is not in the computing range (0)" << std::endl;
				std::clog << func[0] << std::endl;
				abort();
			}
			if (f_xr < func[m_total_points - 1]) {
				std::clog << "Error!!  func=" << f_xr  <<  " is not in the computing range (n)" << std::endl;
				std::clog << func[m_total_points - 1] << std::endl;
				abort();
			}

			/* Get f = value domain */
			size_t i_minu, i_plus, i_midd;
			i_minu = 0;
			i_plus = m_domain_number;

			while(true) {
				i_midd = (i_plus + i_minu)/2;


				if ( (i_plus - i_minu) == 1) {
					dom_value = i_minu;
					break;
				}
				else if (f_xr < func[get_id(i_midd, 0)])
					i_minu = i_midd;
				else
					i_plus = i_midd;
			}

			T X_root = m_lobGrid[0].root_finder(f_xr, get_pointer_dom(func, dom_value));

			return get_xp(dom_value, X_root);

		}

		/**
		 * It find the coordinate where the function f(x) have a given value
		 * in a given domain.
		 * \f{eqnarray*}{
		 *     f(f_xr) - \mathrm{f_xr} = 0;
		 * \f}
		 *
		 * If the function have no root in this domain, the function will return
		 * the closest value.
		 * If the function have several roots, the function will return the first
		 * it finds.
		 *
		 * \param f_xr    Value of the function we
		 * \param func    Value of the function in the collocation points.
		 * \param domain  Domain where the function look for the root.
		 * \return
		 */
		T root_finder_general(T f_xr, T *func, size_t domain)
		{
			T X_root = m_lobGrid[0].root_finder(f_xr, get_pointer_dom(func ,domain));
			return get_xp(domain, X_root);

		}


		/**
		 *
		 * This function compute the following integrals:
		 * \f{eqnarray*}{
		 *     I_L(x) = BC \; + \; \int_{x_L}^{x} \; f(x) \;  dx,
		 * \f}
		 * \f{eqnarray*}{
		 *     I_R(x) = BC \; + \; \int_{x}^{x_R} \; f(x) \;  dx
		 * \f}
		 * \f{eqnarray*}{
		 *     I_C(x) = BC \; + \; \int_{x_0}^{x} \; f(x)  \; dx
		 * \f}
		 * depending if integ_option is fixed to LEFT_BC, RIGHT_BC or CUSTOM_BC.
		 * The term BC represents the boundary condition we impose in one of the boundaries.\n
		 * For getting the integral in the physical coordinates we need to include
		 * the jacobian between the physical and pseudospectral coordinates.
		 * \param func       Array with the value of the function in collocation points
		 * \param integ      Array where it will be stored the value of the integral in the collocation points.
		 * \param in_option  Integral option: from the left, from the right or from a custom point.
		 * \param bndry_cond Value we impose in the boundary. Default value to zero.
		 * \param xp_bc      Value of x0 for CUSTOM_BC. Default value to zero.
		 *
		 */
		void comp_integral(T *func, T *integ, INTEG_OPTION in_option = INTEG_OPTION::LEFT_BC, const T bndry_cond = 0, const T &xp_bc = 0)
		{

			if (m_num_threads == 1)
				_comp_integral_seq(func, integ, in_option, bndry_cond, xp_bc);
			else {
				#if OMP == 1
					_comp_integral_omp(func, integ, in_option, bndry_cond, xp_bc);
				#else
					std::clog << "OMP is not activated. You should not be here.\n";
					abort();
				#endif
			}
		}


		/**
		 *
		 * This function compute the derivative of the whole multidomain
		 * respect to the spectral coordinates or to the multidomain ones.
		 *
		 * \param func          Array with the value of the function in collocation points
		 * \param der           Array where it will be stored the value of the derivative in the collocation points.
		 * \param deriv_option  Compute the derivative respect to multidomain coordinates or to spectral ones.
		 *                       JAC_OPTION::PHYSICAL respect to multidomain. Default.
		 *                       JAC_OPTION::SPECTRAL respect to the spectral coordinate of the domains.
		 *
		 */
		void comp_derivative(T *func, T *der, JAC_OPTION deriv_option = JAC_OPTION::PHYSICAL)
		{
			if (m_num_threads == 1)
				_comp_derivative_seq(func, der, deriv_option);
			else {
				#if OMP == 1
					_comp_derivative_omp(func, der, deriv_option);
				#else
					std::clog << "OMP is not activated. You should not be here.\n";
					abort();
				#endif
			}
		}


		/**
		 *
		 * This function compute the derivative of a single domain
		 * respect to the spectral coordinates or to the multidomain ones.
		 *
		 * \param func          Array with the value of the function in collocation points
		 * \param der           Array where it will be stored the value of the derivative in the collocation points.
		 * \param dom           Domain where to compute the derivative.
		 * \param deriv_option  Compute the derivative respect to multidomain coordinates or to spectral ones.
		 *                       JAC_OPTION::PHYSICAL respect to multidomain. Default.
		 *                       JAC_OPTION::SPECTRAL respect to the spectral coordinate of the domains.
		 *
		 */
		void comp_derivative_domain(T *func, T *der, size_t dom, JAC_OPTION deriv_option = JAC_OPTION::PHYSICAL)
		{
			m_lobGrid[0].compute_1st_der(get_pointer_dom(func, dom), get_pointer_dom(der, dom), deriv_option, \
					get_pointer_dom(m_main_dxp_dX, dom));
		}



		/**
		 * This function computes the spectral coefficients of a function and stores it in a
		 * given array
		 *
		 * \param func   array with the value of the function in collocation points.
		 * \param coeffs array for storing the spectral coefficients.
		 */
		void get_spectral_coeffs(T *func, T *coeffs) {
			for (size_t a = 0; a < m_domain_number; a++) {
				m_lobGrid[0].get_spectral_coeffs(get_pointer_dom(func, a), get_pointer_dom(coeffs, a));
			}
		}



		/**
		 * This function computes the spectral coefficients of a function in a given domain
		 * and stores it in a an array
		 *
		 * \param func   array with the value of the function in collocation points.
		 * \param coeffs array for storing the spectral coefficients. It is assumed to be
		 * 			     of the size of the spectral_doman N+1. If you want to
		 *               to stored in a multidomain array, use get_pointer_dom.
		 * \param dom    Domain where the spectral coefficients will be computed.
		 */
		void get_spectral_coeffs(T *func, T *coeffs, size_t dom) {
			m_lobGrid[0].get_spectral_coeffs(get_pointer_dom(func, dom), coeffs);
		}


		/**
		 * This function returns the pointer to one of the spectral domains used in the
		 * multidomain to perform operations with it.
		 * The multidomain class uses a spectral_domain<T> object for every thread
		 *
		 * \param num_thread   index of the spectral_domain required. Need to be less than the number of threads used by the object.
		 *                   If it is greater it will return
		 */
		spectral_domain<T>& get_spectral_domain(size_t num_thread)
		{
			if (num_thread >= m_num_threads) {
					std::clog << "You are requesting an spectral_domain greater than the number of threads used by the multidomain.\n";
					std::clog << "The argument of get_spectral_domain() needs to be less than that.\n";
					std::clog << "We return an empty spectral_domain<T>.\n";
					abort();
			}
			else {
				return m_lobGrid[num_thread];
			}

		}




	private:

		void _configure_omp(int comp_number_threads)
		{
			if (comp_number_threads == -1) {
				#if OMP == 1
					m_num_threads = omp_get_num_procs();
				#else
					m_num_threads = 1;
				#endif
			}
			else if (comp_number_threads == 1){
				m_num_threads = 1;
				m_dom_beg = new size_t[1];
				m_dom_end = new size_t[1];
				m_dom_beg[0] = 0;
				m_dom_end[0] = m_domain_number - 1;
				return;
			}
			else if (comp_number_threads > 0)
				m_num_threads = comp_number_threads;
			else {
				std::clog << " The number of threads needs to be positive or -1 (take the maximum number of processors available).\n";
				abort();
			}

// 			std::cout << "\nNumber of threads OpenMp  " << threads << std::endl;

			m_dom_beg = new size_t[m_num_threads];
			m_dom_end = new size_t[m_num_threads];

			for (size_t i = 0; i < m_num_threads; i++)
				m_dom_beg[i] = i * std::ceil( (1.*m_domain_number)/m_num_threads);

			for (size_t i = 0; i <= (m_num_threads-2); i++)
				m_dom_end[i] = m_dom_beg[i+1] - 1;

			m_dom_end[m_num_threads - 1] = m_domain_number -1;

			/* Can be a problem for high number of threads and low number of domains
			 (for example dom = 71; threads = 17)*/
			if (m_dom_end[m_num_threads - 1] < m_dom_beg[m_num_threads - 1]) {
				std::cout << "SERIOUS ERROR!!!" << std::endl;
				abort();
			}

// 			for (size_t i = 0; i < num_threads; i++)
// 				std::cout << i << "  " << dom_beg[i] << "  " << dom_end[i] <<  " " << dom_end[i] - dom_beg[i] + 1  << std::endl;

		}

		void _allocate_memory()
		{
			m_main_xp = new T[m_total_points];
			m_main_dX_dxp = new T[m_total_points];
			m_main_dxp_dX = new T[m_total_points];

			m_lobGrid = new spectral_domain<T>[m_num_threads];

			for (size_t i = 0; i < m_num_threads; i++)
				m_lobGrid[i].initialise(m_domain_points);

			m_dom_nodes = new T[m_domain_number + 1];

			/* Dual Grid variables */
			if (m_der_dual) {
				m_dual_xp     = new T[m_total_points];
				m_dual_dX_dxp = new T[m_total_points];
				m_dual_dxp_dX = new T[m_total_points];
				m_partit      = new T[m_total_points];

				m_dualDomain = new spectral_domain<T>[m_num_threads];

			for (size_t i = 0; i < m_num_threads; i++)
				m_dualDomain[i].initialise(m_domain_points);
			}

			/* Private Arrays */
			m_offset = new T[m_domain_number];
			m_aux1   = new T[m_total_points];
			m_aux2   = new T[m_total_points];
		}


		void _set_grid_uniform(T left_bndry, T right_bndry)
		{

			/*----- Create Domain nodes -----*/
			m_dom_nodes[0] = left_bndry;
			m_dom_nodes[m_domain_number] = right_bndry;

			for (size_t a = 1; a < m_domain_number; a++)
				m_dom_nodes[a] = m_dom_nodes[0] + (m_dom_nodes[m_domain_number] - m_dom_nodes[0]) * a / (T) m_domain_number;


			/*----- Construction of the Gauss-Chebyshev-Lobatto Collocation Grid -----*/
			for (size_t a = 0; a < m_domain_number; a++)
			{
				m_main_xp[get_id(a, 0)] = m_dom_nodes[a];
				m_main_xp[get_id(a, m_domain_points)] = m_dom_nodes[a+1];



				for (size_t k = 1; k < m_domain_points; k++)
					m_main_xp[get_id(a, k)] = m_half * ( (m_dom_nodes[a+1] - m_dom_nodes[a]) * m_lobGrid[0].get_X(k)
												+ (m_dom_nodes[a+1] + m_dom_nodes[a]) ) ;

			}

			T unif_dxp_dX = m_half * (m_main_xp[get_id(0, m_domain_points)] - m_main_xp[get_id(0, 0)]);
			T unif_dX_dxp = m_one / unif_dxp_dX;

			for (size_t i = 0; i < m_total_points; i++) {
				m_main_dxp_dX[i] = unif_dxp_dX;
				m_main_dX_dxp[i] = unif_dX_dxp;
			}
		}

		void _set_grid_def_bnds(T left_bndry, T right_bndry, T *boundaries)
		{
			/*----- Create Domain nodes -----*/
			m_dom_nodes[0] = left_bndry;
			m_dom_nodes[m_domain_number] = right_bndry;

			for (size_t a = 1; a < m_domain_number; a++)
				m_dom_nodes[a] = boundaries[a-1];

			/*----- Construction of the Gauss-Chebyshev-Lobatto Collocation Grid -----*/
			for (size_t a = 0; a < m_domain_number; a++) {
				m_main_xp[get_id(a, 0)] = m_dom_nodes[a];
				m_main_xp[get_id(a, m_domain_points)] = m_dom_nodes[a+1];

				for (size_t k = 1; k < m_domain_points; k++)
					m_main_xp[get_id(a, k)] = m_half * ( (m_dom_nodes[a+1] - m_dom_nodes[a]) * m_lobGrid[0].get_X(k)
												+ (m_dom_nodes[a+1] + m_dom_nodes[a]) ) ;
			}

			for (size_t i = 0; i < m_total_points; i++) {
				m_main_dxp_dX[i] = m_half * (m_main_xp[get_id(get_dom(i), m_domain_points)] - m_main_xp[get_id(get_dom(i), 0)]);
				m_main_dX_dxp[i] = m_one / m_main_dxp_dX[i];
			}
		}




		void _set_dual_grid()
		{
			T left_bndy_dom, right_bndy_dom;
			size_t a_dual;

			/*----- Construction of the Gauss-Chebyshev-Lobatto Collocation Grid -----*/
			for (size_t a = 0; a < (m_domain_number-1); a++)
			{
				left_bndy_dom  = m_half * (m_dom_nodes[a+0] + m_dom_nodes[a+1]);
				right_bndy_dom = m_half * (m_dom_nodes[a+1] + m_dom_nodes[a+2]);

				m_dual_xp[get_id(a, 0)] = left_bndy_dom;
				m_dual_xp[get_id(a, m_domain_points)] = right_bndy_dom;



				for (size_t k = 1; k < m_domain_points; k++)
					m_dual_xp[get_id(a, k)] = m_half * ( (right_bndy_dom - left_bndy_dom) * m_lobGrid[0].get_X(k)
												+ (right_bndy_dom + left_bndy_dom) ) ;

			}

			T unif_dxp_dX = m_half * (m_dual_xp[get_id(0, m_domain_points)] - m_dual_xp[get_id(0, 0)]);
			T unif_dX_dxp = m_one / unif_dxp_dX;

			for (size_t i = 0; i < m_total_points; i++) {
				m_dual_dxp_dX[i] = unif_dxp_dX;
				m_dual_dX_dxp[i] = unif_dX_dxp;
			}

			for (size_t a = 0; a < m_domain_number; a++) {
				for (size_t j = 0; j <= m_domain_points; j++) {
					T x = m_main_xp[get_id(a ,j)];
					T xl = get_l_bndry(a);
					T xr = get_r_bndry(a);

					if (get_X(m_main_xp[j]) < m_zero)
						a_dual = (a > 0)? a - 1 : a;
					else
						a_dual = a;

					T xl_d = _dual_get_l_bndry(a_dual);
					T xr_d = _dual_get_r_bndry(a_dual);

					m_partit[get_id(a ,j)] = (x - xl) * (x - xr)
						/ ((x - xl) * (x - xr) + (x - xl_d) * (x - xr_d));

				}
			}

			for (size_t j = 0; j <= m_domain_points/ 2; j++) {
				m_partit[get_id(0 ,j)] = m_one;
				m_partit[get_id(m_domain_number -1, m_domain_points - j)] = m_one;
			}

		}


		void _compute_offsets(T *integ, INTEG_OPTION integ_option)
		{

			if (integ_option == aneto::INTEG_OPTION::RIGHT_BC) {
				m_offset[m_domain_number-1] = m_zero;
				for (int a = (m_domain_number-2); a >= 0; a--)
					m_offset[a] = m_offset[a+1] + integ[get_id(a+1, 0)];
			}else {
				m_offset[0] = m_zero;
				for (size_t a = 1; a < m_domain_number; a++)
					m_offset[a] = m_offset[a-1] + integ[get_id(a-1, m_domain_points)];
			}

		}


		void _set_constants()
		{
			m_zero = 0.00;
			m_half = 0.50;
			m_one  = 1.00;
			m_two  = 2.00;
		}



		/**
		* Returns the number of domains of the dual grid.
		*
		*/
		inline size_t _dual_dom_number() const
		{
			return m_domain_number-1;
		}

		/* Returns the left boundary of a domain.
		* \param dom the domain
		*/
		inline T _dual_get_l_bndry(size_t dom) const
		{
			return m_dual_xp[get_id(dom, 0)];
		}

		/**
		* Returns the right boundary of a domain.
		* \param dom the domain
		*
		*/
		inline T _dual_get_r_bndry(size_t dom) const
		{
			return m_dual_xp[get_id(dom, m_domain_points)];
		}

		/*
		* Returns the left global boundary.
		*/
		inline T _dual_get_l_glob_bndry() const
		{
			return m_dual_xp[get_id(0, 0)];
		}

		/*
		* Returns the right global boundary.
		*/
		inline T _dual_get_r_glob_bndry() const
		{
			return m_dual_xp[get_id(_dual_dom_number() -1, m_domain_points)];
		}

		inline T _dual_get_X(T x_i) {
			size_t dom = _dual_get_dom_coord(x_i);

			return (m_two * x_i - (_dual_get_l_bndry(dom) + _dual_get_r_bndry(dom)) )/(_dual_get_r_bndry(dom) - _dual_get_l_bndry(dom));
		}

		inline T _dual_get_xp(size_t dom, T X) {
			return m_half * ((_dual_get_l_bndry(dom) + _dual_get_r_bndry(dom)) + (_dual_get_r_bndry(dom) - _dual_get_l_bndry(dom)) * X);
		}



		inline size_t _dual_get_dom_coord(T x_i) const
		{
			size_t i_minu, i_plus, i_midd;

			i_minu = 0;
			i_plus = m_domain_number - 1;

			// check v is in range
			if (x_i < _dual_get_l_glob_bndry()) {
				std::clog << "Error!!  v = " <<  x_i << " is below the computing range (x_L: " << _dual_get_l_glob_bndry() << ")" << std::endl;
				abort();
			}
			if (x_i > _dual_get_r_glob_bndry()) {
				std::clog << "Error!!  v = " <<  x_i << " is above the computing range (x_R: " << _dual_get_r_glob_bndry() << ")" << std::endl;
				abort();
			}

			while(true) {
				i_midd = (i_plus + i_minu)/2;


				if ( (i_plus - i_minu) == 1)
					return i_minu;
				else if (x_i < _dual_get_l_bndry(i_midd))
					i_plus = i_midd;
				else
					i_minu = i_midd;
			}
		}



		void _comp_integral_seq(T *y, T *integ, INTEG_OPTION integ_option = INTEG_OPTION::LEFT_BC, const T bndry_cond = 0, const T &xp_bc = 0)
		{
			INTEG_OPTION internal_int_opt;

			if (integ_option == INTEG_OPTION::RIGHT_BC)
				internal_int_opt = INTEG_OPTION::RIGHT_BC;
			else
				internal_int_opt = INTEG_OPTION::LEFT_BC;


			/* Computing partial integrals */
			for (size_t a = 0; a < m_domain_number; a++)
				m_lobGrid[0].compute_integral(get_pointer_dom(y, a), get_pointer_dom(integ, a), internal_int_opt, m_zero, JAC_OPTION::PHYSICAL, get_pointer_dom(m_main_dxp_dX, a));

			/* Computing offsets */
			_compute_offsets(integ, integ_option);

			/* Computing full integrals */
			if (integ_option == INTEG_OPTION::RIGHT_BC) {
				for (size_t i = 0; i < m_total_points; i++)
					integ[i] = bndry_cond + m_offset[get_dom(i)] + integ[i];
			}
			else if (integ_option == INTEG_OPTION::LEFT_BC) {
				for (size_t i = 0; i < m_total_points; i++)
					integ[i] = bndry_cond + m_offset[get_dom(i)] + integ[i];
			}
			else if(integ_option == INTEG_OPTION::CUSTOM_BC) {
				/* Compute bc constant in case CUSTOM_BC */
				T c_integ = - (interpolate(integ, xp_bc) + m_offset[get_dom_coord(xp_bc)]);
				for (size_t i = 0; i < m_total_points; i++)
					integ[i] = bndry_cond + c_integ + m_offset[get_dom(i)] + integ[i];
			}

		}

		void _comp_derivative_seq(T *func, T *der, JAC_OPTION deriv_option = JAC_OPTION::PHYSICAL)
		{
			/* Computing Derivative */
			for (size_t a = 0; a < m_domain_number; a++)
				m_lobGrid[0].compute_1st_der(get_pointer_dom(func, a), get_pointer_dom(der, a), deriv_option, get_pointer_dom(m_main_dxp_dX, a));

			if (m_der_dual == false)
				return;

			/* Computing with dual grid */
			for (size_t a = 0; a < _dual_dom_number(); a++) {
				/* Interpolate to dual grid */
				m_lobGrid[0].compute_spectral_func_der(get_pointer_dom(func, a), JAC_OPTION::PHYSICAL, get_pointer_dom(m_main_dxp_dX, a), COEFF_OPTION::F);
				for (size_t i = 0; i <= m_domain_points/2; i++)
					m_aux1[get_id(a, i)] = m_lobGrid[0].get_fun(get_X(m_dual_xp[get_id(a, i)]));  // idx0 = 0 does not matter when use_spec_stored == true

				m_lobGrid[0].compute_spectral_func_der(get_pointer_dom(func, a+1), JAC_OPTION::PHYSICAL, get_pointer_dom(m_main_dxp_dX, a+1), COEFF_OPTION::F);
				for (size_t i = m_domain_points/2 + 1; i <= m_domain_points; i++)
					m_aux1[get_id(a, i)] = m_lobGrid[0].get_fun(get_X(m_dual_xp[get_id(a, i)]));  // idx0 = 0 does not matter when use_spec_stored == true

				/* Middle point, if even */
				if (m_domain_points%2 == 0)
					m_aux1[get_id(a, m_domain_points/2)] = func[get_id(a + 1, 0)];

				/* Perform derivative */
				m_dualDomain[0].compute_spectral_func_der(get_pointer_dom(m_aux1, a), JAC_OPTION::PHYSICAL, get_pointer_dom(m_dual_dxp_dX, a), COEFF_OPTION::D1);

				for (size_t i = m_domain_points/2 + 1; i <= m_domain_points; i++)
					m_aux2[get_id(a, i)] = m_dualDomain[0].get_der(_dual_get_X(m_main_xp[get_id(a, i)]));

				for (size_t i = 0; i <= m_domain_points/2; i++)
					m_aux2[get_id(a+1, i)] = m_dualDomain[0].get_der(_dual_get_X(m_main_xp[get_id(a+1, i)]));

				// This point is not important because partit = 1.
				// This just assure that aux2 != nan. This can be removed if we change the previous loop.
				if (m_domain_points%2 == 0)
					m_aux2[get_id(a, m_domain_points/2)] = m_zero;
			}

			for (size_t i = 0; i < m_total_points; i++)
				der[i] = der[i] * m_partit[i] + (m_one - m_partit[i]) * m_aux2[i];

		}



#if OMP == 1

		void _comp_integral_omp(T *y, T *integ, INTEG_OPTION integ_option = INTEG_OPTION::LEFT_BC, const T bndry_cond = 0, const T &xp_bc = 0)
		{
// 			Extrae_init();
// 			Extrae_define_event_type (&type, "Kernel execution", &nvalues, values, description_values);



// 			#pragma omp parallel
			#pragma omp parallel num_threads(m_num_threads)
			{
				int idt = omp_get_thread_num();
// 				std::cout << "  #  " << idt << "   ";

// 				std::stringstream stream;
// 				stream << "#T: " << idt << "  - prec: " << mpfr::mpreal::get_default_prec() << std::endl;
// 				std::cout << stream.str();

				/* Computing partial integrals */
				for (size_t a = m_dom_beg[idt]; a <= m_dom_end[idt]; a++){
// 					Extrae_event(1, a);
					m_lobGrid[idt].compute_integral(get_pointer_dom(y, a), get_pointer_dom(integ, a), integ_option, m_zero, JAC_OPTION::PHYSICAL, get_pointer_dom(m_main_dxp_dX, a));
				}
// 			}

#pragma omp barrier

// 	Extrae_event (1000, 1);

				/* Computing offsets */
				if (idt == 0){
// 					Extrae_event(2, 0);

					_compute_offsets(integ, integ_option);
				}
// 	Extrae_event (1000, 0);

				#pragma omp barrier

// 			#pragma omp parallel num_threads(m_num_threads)
// 			{
// 				int idt = omp_get_thread_num();

				/* Computing full integrals */
				if (integ_option == INTEG_OPTION::RIGHT_BC) {
					for (size_t i = get_id(m_dom_beg[idt], 0); i <= get_id(m_dom_end[idt], dom_points()); i++)
						integ[i] = bndry_cond + integ[m_total_points-1] + m_offset[m_domain_number - 1] - m_offset[get_dom(i)] - integ[i];
				}
				else if (integ_option == INTEG_OPTION::LEFT_BC) {
					for (size_t i = get_id(m_dom_beg[idt], 0); i <= get_id(m_dom_end[idt], dom_points()); i++) {
// 						Extrae_event(3, i);

						integ[i] = bndry_cond + m_offset[get_dom(i)] + integ[i];
					}
				}
				else if(integ_option == INTEG_OPTION::CUSTOM_BC) {
					/* Compute bc constant in case CUSTOM_BC */
					T c_integ = - (interpolate(integ, xp_bc) + m_offset[get_dom_coord(xp_bc)]);
					for (size_t i = get_id(m_dom_beg[idt], 0); i <= get_id(m_dom_end[idt], dom_points()); i++)
						integ[i] = bndry_cond + c_integ + m_offset[get_dom(i)] + integ[i];
				}
			}

// 			Extrae_fini();

		}


		void _comp_derivative_omp(T *func, T *der, JAC_OPTION deriv_option = JAC_OPTION::PHYSICAL)
		{
// 			#pragma omp parallel
			#pragma omp parallel num_threads(m_num_threads) default(shared)
			{
				int idt = omp_get_thread_num();
// 				std::cout << "  #  " << idt << "   \n";

			/* Computing Derivative */
				for (size_t a = m_dom_beg[idt]; a <= m_dom_end[idt]; a++){
					m_lobGrid[idt].compute_1st_der(get_pointer_dom(func, a), get_pointer_dom(der, a) , deriv_option, get_pointer_dom(m_main_dxp_dX, a));
				}
//
// 				#pragma omp for
// 				for (size_t a = 0; a < domain_number; a++){
// 					lobGrid[idt].compute_1st_der(func, der, deriv_option, main_dxp_dX, get_id(a,0), get_id(a,0));
// 				}



			}

// 				int idt = omp_get_thread_num();
// 				std::cout << "  #  " << idt << "   \n";

		}

#endif



	};

}



#endif /* _MULTIDOMAIN_HPP */
