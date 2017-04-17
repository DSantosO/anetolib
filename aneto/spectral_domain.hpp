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


#ifndef _SPECTRAL_DOMAIN_HPP
#define _SPECTRAL_DOMAIN_HPP

#include <iostream>

#include <eigen3/unsupported/Eigen/FFT>

namespace aneto {

	/**
	 * Jacobian Option for the derivative / integral. It allows the operations to be performed in
	 * the coordinate of the spectral domain or in the one of the multidomain.
	 */
	enum class JAC_OPTION
	{
		SPECTRAL,	///< The derivative is computed respect to the spectral coordinate (X = [-1,+1]).
		PHYSICAL	///< The derivative is computed respect to the physical coordinate (multidomain).
	};

	/**
	 * Option that controls how many levels of spectral coefficients will be stored.
	 */
	enum class COEFF_OPTION
	{
		F,			///< We stored the spectral coefficients of the function.
		D1,			///< We stored the spectral coefficients of the function and the first derivative.
		D2,			///< We stored the spectral coefficients of the function and the first and second derivatives.
	};

	/**
	 * Possible options for imposing a boundary condition (BC)
	 * in the computation of the multidomain integral.
	 * Exact details can be seen in comp_integral.
	 * Be aware that in the spectral_domain function, the CUSTOM_BC is not accepted.
	 */
	enum class INTEG_OPTION
	{
		LEFT_BC,	///< Impose the BC in the left global boundary.
		RIGHT_BC,	///< Impose the BC in the right global boundary.
		CUSTOM_BC	///< Impose the BC in an arbitrary point of the multidomain.
	};

	/**
	 * Possible options for making the tranformation between the
	 * collocation points and the spectral representation.
	 * It can be done by matrix multiplication (~N^2) or by
	 * FFT (~N log(N)) optimized for number of points 2^i with
	 * i integer. The default option will select FFT if N = 2^i
	 * and matrix multiplication otherwise.
	 */
	enum class TRANSF_OPTION
	{
		MATRIX,	///< Transformation from matrix multiplication.
		FFT,	///< Transformation through FFT.
		DEFAULT	///< Select FFT option if N is a power of 2 (N = 2^i) and MATRIX otherwise.
	};



	/** \brief Template of a single spectral domain defined in X = [-1,+1]
	 *
	 * The spectral domain is generated following the Lobatto-Chebychev points:
	 * \f{eqnarray*}{
	 *     X_k = - \cos(k \pi / N);
	 * \f}
	 * Notice that the grid is defined in the range [-1, 1]. Also the number
	 * of points is N + 1 (from 0 to N) but for convention N
	 * will be denoted as grid points.
	 * \n\n
	 * The template can be used with any general type that represent real numbers where
	 * the usual operations + - * / are defined.\n
	 * Also a cos / acos functions are required for the specific type.\n
	 * \tparam T type to be used for the spectral grid
	 */
	template <class T> class spectral_domain
	{
	private:
		bool m_allocated;			// Indicate if the object is allocated.
		TRANSF_OPTION m_transf;	// Option for Chebyshev transformation
		Eigen::FFT<T> m_fft;		// Plan for compute the FFT
		std::vector<T> m_fft_func;
		std::vector<std::complex<T>> m_fft_coeffs;


		int m_spec_points;  	// Number of points of the Lobatto-Chebyshev Grid (N+1);
		T *m_X;					// Vector that stores the points of the Lobatto grid: cos( i * PI / points )
		T *m_norm		;		// Normalisation coeff: two in the boundaries, one in interior points
		T *m_weight;				// Weight: PI/(points * normalisation)
		T **m_mat_transf;					// Matrix coefficients of the transformation.

		T *m_v_in, *m_v_out;		// Auxiliary vectors for computation.
		T *m_coeffs_func;			// Aux vector that stores the spectral coeffs of the function.
		T *m_coeffs_der;			// Aux vector that stores the spectral coeffs of the derivative.
		T *m_coeffs_d2;			// Aux vector that stores the spectral coeffs of the second derivative.
		T *m_aux_der;				// Aux vector for derivative computation.

		T m_pi;			///< Pi given up to the accuracy of the created object.
		T m_zero;		///< Zero given up to the accuracy of the created object.
		T m_one;		///< One given up to the accuracy of the created object.
		T m_two;		///< Two given up to the accuracy of the created object.
		T m_half;		///< 1/2 given up to the accuracy of the created object.


	public:
		/** Void constructor.
		 *
		 * A spectral grid is created empty.
		 * It needs to be initialised in order to be used.
		 */
		spectral_domain() : m_allocated(false), m_spec_points(0) {}

		/** General Constructor.
		 *
		 * It creates a Lobatto-Chebyshev grid of po+1 points.
		 * Once created can not be modified until is clean out or deleted
		 * \param po           Number of points (po + 1)
		 * \param transf_op    Possible options for the Chebyshev tranformation.
		 *                     TRANSF_OPTION::DEFAULT is enable by default.
		 *
		 */
		spectral_domain(int po, TRANSF_OPTION transf_op = TRANSF_OPTION::DEFAULT) {
			m_allocated = false;
			initialise(po, transf_op);
		}

		/** Destructor.
		 *
		 * The used memory is freed.
		 */
		~spectral_domain()
		{
			if (m_allocated)
				_free_memory();
		}


		/**
		 * Spectral Grid Initialiser
		 *
		 * It initialises a Lobatto-Chebyshev grid of po+1 points.
		 * If the object was in use it will abort the program.
		 * \param po Number of points (po + 1)
		 * \param transf_op    Possible options for the Chebyshev tranformation.
		 *                     TRANSF_OPTION::DEFAULT is enable by default.
		 */
		void initialise(int po, TRANSF_OPTION transf_op = TRANSF_OPTION::DEFAULT) {

			if (m_allocated == true) {
				std::clog << "Error! You are trying to initialise an existing object." << std::endl;
				std::abort();
			}

			if (po >= 3) {
				m_spec_points = po;
			}
			else {
				std::clog << "Please introduced a number of points than three." << std::endl;
				abort();
			}

 			/* If TRANSF_OPTION is DEFAULT, we select FFT or MATRIX */
			if (transf_op == TRANSF_OPTION::DEFAULT) {
				/* Cheking if power of two */
				if ((po & (po - 1)) == 0) {
					m_transf = TRANSF_OPTION::FFT;
				}
				else {
					m_transf = TRANSF_OPTION::MATRIX;
				}
			}
			else {
				m_transf = transf_op;
			}

			_allocate_memory();
			_set_constants();

			_generate_points();
			_compute_transf_matrix();

			m_allocated = true;
		}

		/**
		 * Cleaning the instance.
		 *
		 * It frees the memory of the spectral grid.
		 * Once it is clean, can be initialised again.
		 */
		void clean()
		{
			if (m_allocated == false)
				return;

			_free_memory();
			m_allocated = false;
		}


		/**
		 * Return the number of points of the spectral grid.
		 * Notice that according with our convention, the grid
		 * have N+1 points and this function returns N.
		*/
		int get_points()
		{
			return m_spec_points;
		}

		/**
		 * Return the spectral coordinate of point k
		 * Be careful! The index is not checked.
		 * \param k point of the grid.
		*/
		T get_X (int k)
		{
			return m_X[k];
		}

		/** Chebyshev Transformation.
		 *
		 * A Real Fourier Transformation (Discrete Cosine) that gives us the spectral coefficients of the function.
		 * The exact normalisation is:
		 * \f{eqnarray*}{
		 *     \hat f_i = \sum_{j=0}^{N} \frac{2}{N \; m_j} \cos(i j \pi / N)
		 * \f}
		 * \param in   array with collocation points of the function.
		 * \param out  array with the spectral coefficients.
		*/
		void cheb_transf(T *in, T *out)
		{
			if (m_transf == TRANSF_OPTION::FFT)
				_fft_transf(in, out);
			else if (m_transf == TRANSF_OPTION::MATRIX)
				_matrix_transf(in, out);
			else {
				std::clog << "Error! Invalid transformation option.\n" << std::endl;
				abort();
			}

		}

		/**
		 *
		 * This function compute the following expression:
		 * \f{eqnarray*}{
		 *     I_L(xp) = BC \; + \; \int_{-1}^{X} \; f(X)  \;  \frac{ d xp}{dX}  \;  dX
		 * \f}
		 * or
		 * \f{eqnarray*}{
		 *     I_R(xp) = BC \; + \; \int_{X}^{+1} \; f(X)  \;  \frac{ d xp}{dX}  \;  dX
		 * \f}
		 * depending if i_option is fixed to INTEG_OPTION::LEFT_BC or INTEG_OPTION::RIGHT_BC
		 * The term BC represents the boundary condition we impose in one of the boundaries.
		 * Notice that the option INTEG_OPTION::CUSTOM_BC is not allowed in this function.\n
		 * For getting the integral in the physical coordinates we need to include
		 * the Jacobian between the physical and pseudospectral coordinates.
		 * \param func       array with the value of the function in collocation points
		 * \param Integ      array where it will be stored the value of the integral in the collocation points.
		 * \param i_option   Integral option from the left or from the right
		 * \param dxp_dX     array with the Jacobian values in the collocation points
		 * \param bndry_cond value we impose in the boundary. Default value to zero.
		 * \param jac_option Use the spectral coordinate or including the Jacobian.
		 * \param dxp_dX     Jacobian of the transformation between the spectral coordinate and a physical one.
		 *                   If it is no needed, NULL can be provided as parameter.
		 *
		 */
		void compute_integral(T *func, T *Integ, INTEG_OPTION i_option = INTEG_OPTION::LEFT_BC, T bndry_cond = 0, JAC_OPTION jac_option = JAC_OPTION::SPECTRAL, T *dxp_dX = NULL)
		{
			if (jac_option == JAC_OPTION::PHYSICAL and dxp_dX == NULL) {
				std::clog << "Error! For physical integral, you need to introduce the Jacobian.\n" << std::endl;
				abort();
			}

			T* aux_in;

			if (jac_option == JAC_OPTION::PHYSICAL) {
				for (int k = 0; k <= m_spec_points; k++)
					m_v_in[k] =  dxp_dX[k] * func[k];
				aux_in = m_v_in;
			}
			else {
				aux_in = func;
			}



			if (i_option == INTEG_OPTION::LEFT_BC) {
				_compute_integral_left(aux_in, Integ, bndry_cond);
			}
			else if (i_option == INTEG_OPTION::RIGHT_BC) {
				_compute_integral_right(aux_in, Integ, bndry_cond);
			}
			else {
				std::clog << "In spectral domain only INTEG_OPTION::LEFT_BC and INTEG_OPTION::RIGHT_BC are valid as an integral options2.\n";
				std::clog << "Please check that and rerun ANETO.\n";
				abort();
			}

		}


		/**
		 * Compute the first derivative.
		 *
		 * This function compute spectral derivative of a function.\n
		 * By default compute the derivative respect to the spectral coordinate X.
		 * If the physical option is enable, the Jacobian of the transformation needs to be included.\n
		 * If the collocation points of the function or the derivative are stored inside a bigger array
		 * the index to the first element can be included. This is convenient in the multidomain case.
		 *
		 * \param func  array with the value of the function in collocation points.
		 * \param der   array where it will be stored the value of the integral in the collocation points.
		 * \param der_option  option for spectral or physical derivative.
		 * \param dxp_dX    array with the Jacobian values in the collocation points.
		 */
		void compute_1st_der(T *func, T *der, JAC_OPTION der_option = JAC_OPTION::SPECTRAL, T *dxp_dX = NULL)
		{
			if (der_option == JAC_OPTION::PHYSICAL and dxp_dX == NULL) {
				std::clog << "Error! For physical derivative, you need to introduce the Jacobian.\n" << std::endl;
				abort();
			}

			cheb_transf(func, m_v_out);

			m_v_in[m_spec_points] = m_zero;
			m_v_in[m_spec_points - 1] = -m_half * m_v_out[m_spec_points];

			for (int k = (m_spec_points - 2); k >= 0; k--)
				m_v_in[k] = m_v_in[k + 2] - (k+1) * m_v_out[k+1] / (T) m_spec_points;

			cheb_transf(m_v_in, der);


			if (der_option == JAC_OPTION::SPECTRAL)
				return;

			/* For physical derivative, we need to apply the Jacobian */
			for (int k = 0; k <= m_spec_points; k++)
				der[k] /= dxp_dX[k];
		}


		/**
		 *
		 * This function computes the value of a function in any point of the domain.
		 * through the spectral representation. The spectral coefficients computed here
		 * are stored and can be used later calling the function get_X.
		 *
		 * \param X spectral coordinate where the function will be evaluated.
		 * \param func  array with the value of the function in collocation points.
		 */
		T interpolate(T X, T *func)
		{
			_compute_spectral_coeffs(func, m_coeffs_func);

			T interpolated_value = 0;
			for (int k = 0; k <= m_spec_points; k++)
				interpolated_value += m_coeffs_func[k] * _chebyshev_polynomial(k, X);

			return interpolated_value;
		}

		/**
		 * This function computes the spectral coefficients of a function and stores it in a
		 * given array. It does not stored them inside the spectral_domain.
		 *
		 * \param func   array with the value of the function in collocation points.
		 * \param coeffs array for storing the spectral coefficients.
		 */
		void get_spectral_coeffs(T *func, T *coeffs) {
			_compute_spectral_coeffs(func, coeffs);
		}



		/**
		 *
		 * This function compute the value of a function in any point of the domain.
		 * It uses a previously computed spectral coefficients that needs to
		 * be passed as parameter.
		 *
		 * \param X spectral coordinate where the function will be evaluated.
		 * \param coeffs spectral coefficients of the function already computed.
		 */
		T interpolate_with_spec_coeffs(T X, T *coeffs)
		{
			/*----- Interpolation -----*/
			T interpolated_value = 0;
			for(int k = 0; k <= m_spec_points; k++)
				interpolated_value += coeffs[k] * _chebyshev_polynomial(k, X);

			return interpolated_value;
		}


		/**
		 *
		 * This function computes the spectral coefficients and stores them internally
		 * to be used when the functions get_fun(), get_der(), get_d2() are called.
		 * They are also used if interpolate function is called with use_spec_stored
		 * enabled.
		 *
		 * An option can be included to select if you want to stored just
		 * the coeffs of the function (F), the function and the first derivative (D1)
		 * or the function and the two first derivatives (D2).
		 *
		 * \param func  array with the value of the function in collocation points.
		 * \param der_option  option for spectral or physical derivative.
		 * \param dxp_dX    array with the Jacobian values in the collocation points.
		 * \param coeff_op  option for selecting which coeffs to stored.
		 */
		void compute_spectral_func_der(T *func, JAC_OPTION der_option = JAC_OPTION::SPECTRAL, T *dxp_dX = NULL, COEFF_OPTION coeff_op = COEFF_OPTION::D1)
		{
			_compute_spectral_coeffs(func, m_coeffs_func);
			if (coeff_op == COEFF_OPTION::F)
				return;

			compute_1st_der(func, m_aux_der, der_option, dxp_dX);
			_compute_spectral_coeffs(m_aux_der, m_coeffs_der);

			if (coeff_op == COEFF_OPTION::D1)
				return;

			compute_1st_der(m_aux_der, m_aux_der, der_option, dxp_dX);
			_compute_spectral_coeffs(m_aux_der, m_coeffs_d2);

			if (coeff_op == COEFF_OPTION::D2)
				return;
		}

		/**
		 * Compute the value of the function at any point of the domain using they
		 * stored spectral coefficients.
		 * \param X coordinate.
		 * \return interpolated value of the function.
		 */
		T get_fun(T X)
		{
			/*----- Interpolation -----*/
			T interpolated_value = m_zero;
			for(int k = 0; k <= m_spec_points; k++)
				interpolated_value += m_coeffs_func[k] * _chebyshev_polynomial(k, X);

			return interpolated_value;
		}

		/**
		 * Compute the value of the derivative function at any point
		 * of the domain using they
		 * stored spectral coefficients.
		 * \param X coordinate.
		 * \return interpolated value of the derivative.
		 */
		T get_der(T X)
		{
			/*----- Interpolation -----*/
			T interpolated_value = m_zero;
			for(int k = 0; k <= m_spec_points; k++)
				interpolated_value += m_coeffs_der[k] * _chebyshev_polynomial(k, X);

			return interpolated_value;
		}

		/**
		 * Compute the value of the second derivative function at any point
		 * of the domain using they
		 * stored spectral coefficients.
		 * \param X coordinate.
		 * \return interpolated value of the second derivative.
		 */
		T get_d2(T X)
		{
			/*----- Interpolation -----*/
			T interpolated_value = m_zero;
			for(int k = 0; k <= m_spec_points; k++)
				interpolated_value += m_coeffs_d2[k] * _chebyshev_polynomial(k, X);

			return interpolated_value;
		}


		/**
		 * It find the coordinate where the function f(X) have a given value.
		 * \f{eqnarray*}{
		 *     f(f_Xr) - \mathrm{f_Xr} = 0;
		 * \f}
		 *
		 * If the function have no root in the grid, the function will return
		 * the closest value.
		 * If the function have several roots, the function will return the first
		 * it finds.
		 *
		 * \param f_Xr Value of the function we
		 * \param func Value of the function in the collocation points.
		 * \return
		 */
		T root_finder(T f_Xr, T *func)
		{
			/* Find the root using Newton's method */
			T X_root, delta_X, delta_X_old;
			delta_X_old = 1e100;

			compute_spectral_func_der(func, JAC_OPTION::SPECTRAL, NULL);

			X_root  = m_zero;
			while(true) {
				delta_X = - (get_fun(X_root) - f_Xr)  / get_der(X_root);

				if (delta_X == m_zero) {
					return X_root;
				}
				else if (_sd_abs(delta_X_old/delta_X) < 1.0001) {
					return X_root;
				}

				/* Next iteration */
				X_root = X_root + delta_X;
				delta_X_old = delta_X;

				/* Check in case Newton's method goes away spectral domain */
				if (X_root < -m_one)
					X_root = -m_one;
				if (X_root > m_one)
					X_root = m_one;

			}

			std::cout <<  "Warning, Root finding failed!!"  << std::endl;
			abort();

		}

		/**
		 * Returns Pi to the accuracy used by the object.
		 */
		T pi()
		{
			return m_pi;
		}



		private:

		void _allocate_memory()
		{

			m_X      = new T[m_spec_points+1];
			m_weight = new T[m_spec_points+1];
			m_norm = new T[m_spec_points+1];

			if (m_transf == TRANSF_OPTION::MATRIX) {
				m_mat_transf = new T *[m_spec_points+1];
				for(int i = 0; i <= m_spec_points; ++i)
					m_mat_transf[i] = new T[m_spec_points+1];
			}
			else {
				m_fft.SetFlag(m_fft.HalfSpectrum);

				m_fft_func.resize(2 * m_spec_points);
				m_fft_coeffs.resize(2 * m_spec_points);
			}

			// Private
			m_v_in = new T[m_spec_points+1];
			m_v_out = new T[m_spec_points+1];
			m_coeffs_func = new T[m_spec_points+1];
			m_coeffs_der = new T[m_spec_points+1];
			m_coeffs_d2 = new T[m_spec_points+1];
			m_aux_der = new T[m_spec_points+1];

		}

		void _free_memory()
		{
			delete [] m_X;
			delete [] m_weight;
			delete [] m_norm;

			if (m_transf == TRANSF_OPTION::MATRIX) {
				for(int i = 0; i <= m_spec_points; ++i)
					delete [] m_mat_transf[i];
				delete [] m_mat_transf;
			}
			else {


			}

			delete [] m_v_in;
			delete [] m_v_out;
			delete [] m_coeffs_func;
			delete [] m_coeffs_der;
			delete [] m_coeffs_d2;
			delete [] m_aux_der;

			m_spec_points = 0;
		}


		T _chebyshev_polynomial(int k, T X)
		{
			return cos((T)k * acos(X));
		}



		void _compute_transf_matrix()
		{

			if (m_transf != TRANSF_OPTION::MATRIX)
				return;

			for (int i = 0; i <= m_spec_points;  i++)
				for (int j = 0; j <= m_spec_points;  j++)
					m_mat_transf[i][j]    =  cos((T)(i*j) * m_pi/(T)m_spec_points)* m_two /m_norm[j]/m_spec_points;
		}

		void _generate_points()
		{
			if (m_spec_points == 0) {
				std::clog << "Warning! the points attribute in the lobatto_domain struct is zero. Cannot create matrix transformation" << std::endl;
				return;
			}

			m_X[0]      = - m_one;
			m_X[m_spec_points] = m_one;
			m_norm[0]   = m_two;
			m_norm[m_spec_points] = m_two;
			m_weight[0]   = m_pi / ((T)m_spec_points * m_norm[0]);
			m_weight[m_spec_points] = m_weight[0];


			for (int k = 1; k <= (m_spec_points-1)/2; k++)
			{
				m_X[k] = - cos((T) k * m_pi/(T)m_spec_points);
				m_X[m_spec_points - k] = - m_X[k];
			}

			if (m_spec_points%2 == 0)
				m_X[m_spec_points/2] = 0.000;


			T w = m_pi / ((T)m_spec_points);

			for (int k = 1; k < m_spec_points; k++) {
				m_norm[k] = m_one;
				m_weight[k] = w;
			}

		}

		void _compute_spectral_coeffs(T *in, T *coeffs)
		{
			T d_N = (T) m_spec_points;

			cheb_transf(in, m_v_out);

			/**----- Computing the spectral coefficients, {GF_a_spectral[J,A]_n} -----**/
			for (int n = 0; n <= m_spec_points; n++){
				if ( (n % 2) == 0 )
					coeffs[n] =  m_v_out[n]/(d_N * m_norm[n]);
				else
					coeffs[n] = - m_v_out[n]/(d_N * m_norm[n]);
			}

		}


		void _compute_integral_left (T *v_integrand, T *Integ, T left_bndry_cond)
		{
			T dN = (T) m_spec_points;
			T coeff = m_half / (T) m_spec_points;

			cheb_transf(v_integrand, m_v_out);

			m_v_in[m_spec_points] = - (coeff * m_v_out[m_spec_points - 1])/dN;

			for (int k = 1; k < m_spec_points; k++)
				m_v_in[k] = m_half * coeff * (m_v_out[k+1]/(m_norm[k+1]) - m_v_out[k-1])/((T) k);

			m_v_in[0] = m_zero;

			for (int n = 1; n <= m_spec_points; n++) {
				m_v_in[0] += m_v_in[n];
			}

			m_v_in[0] = left_bndry_cond - m_v_in[0] * m_two;

			cheb_transf(m_v_in, Integ);

			Integ[0] = left_bndry_cond;

		}

		void _compute_integral_right (T *v_integrand, T *Integ, T right_bndry_cond)
		{
			T dN = (T) m_spec_points;
			T coeff = m_half / (T) m_spec_points;

			cheb_transf(v_integrand, m_v_out);


			m_v_in[m_spec_points] = - (coeff * m_v_out[m_spec_points - 1])/dN;

			for (int k = 1; k < m_spec_points; k++)
				m_v_in[k] = m_half * coeff * (m_v_out[k+1]/(m_norm[k+1]) - m_v_out[k-1])/((T) k);

			m_v_in[0] = m_zero;
			for (int n = 1; n <= m_spec_points; n++) {
				m_v_in[0] += m_v_in[n] * (n%2==0?m_one:-m_one);
			}
			m_v_in[0] = right_bndry_cond - m_v_in[0] * m_two;

			cheb_transf(m_v_in, m_v_out);

			for (int k = 0; k <= m_spec_points; k++)
				Integ[k] = -m_v_out[k];

			Integ[m_spec_points] = right_bndry_cond;

		}

		void _matrix_transf(T *in, T *out)
		{
			T sum = m_zero;

			for (int i = 0; i <= m_spec_points;  i++) {
				sum = m_zero;
				for (int j = 0; j <= m_spec_points;  j++)
					sum += m_mat_transf[i][j] * in[j];
				out[i] = sum * (T)m_spec_points;
			}
		}

		void _fft_transf(T *in, T *out)
		{
			for (int k = 0; k <= m_spec_points; k++)
				m_fft_func[k] = in[k];
			for (int k = m_spec_points+1; k < 2*m_spec_points; k++)
				m_fft_func[k] = in[2*m_spec_points - k];

			m_fft.fwd(m_fft_coeffs, m_fft_func);

			for (int i= 0; i <= m_spec_points; i++)
				out[i] = m_fft_coeffs[i].real();

		}



		void _set_constants()
		{
			m_zero = 0.00;
			m_half = 0.50;
			m_one  = 1.00;
			m_two  = 2.00;
			m_pi   = acos(-m_one);
		}

		T _sd_abs(T x) {
				return (x > 0)? x : -x;
		}

	};




}// end namespace



#endif /* _SPECTRAL_DOMAIN_HPP */
