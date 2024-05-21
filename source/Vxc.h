#pragma once
//#define _USE_MATH_DEFINES

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

/*************************************
LDA  and LSDA
J. P. Perdew, Alex Zunger, Phys. Rev. B 23, 5048, (1981)
*************************************/

//#include "Vxc.h"


inline void calc_LSDA_ge1(double& ec, double& n_dec, const double rs, const double gamma, const double b1, const double b2) {
	const double sqrtrs = std::sqrt(rs);
	ec = gamma / (1.0 + b1 * sqrtrs + b2 * rs);
	n_dec = ec * ec * (b1 *0.5*sqrtrs + b2*rs) / (3.0*gamma);
}

inline void calc_LSDA_lt1(double& ec, double& n_dec, const double rs, const double A, const double B, const double C, const double D) {
	const double logrs = log(rs);
	ec = B + A * logrs + D * rs + C * rs * logrs;
	n_dec = -(A + (D + C) * rs + C * rs * logrs) / 3.0;
}


struct Return_XC
{
	double E_den_x; //exchange energy density, which is not energy//
	double E_den_c; //correlation energy density, which is not energy//
	double V_x; //exchange energy density, which is not energy//
	double V_c; //correlation energy density, which is not energy//

};

/*
Vxcと共にExcを返す
を返す
*/
inline Return_XC Calc_XC_LDA(double rho) {

	static const double c2 = std::cbrt(3.0 / (4.0*M_PI));
	static const double c3 = -(3.0 / (4.0*M_PI))* std::cbrt(9.0*M_PI / 4.0);

	double E_den_x = 0.0;
	double E_den_c = 0.0;
	double Vx = 0.0;
	double Vc = 0.0;

	if (rho < 1.0e-8) {
		return { E_den_x,E_den_c,Vx,Vc };
	}

	const double rs = c2 / std::cbrt(rho);
	const double ex_0 = c3 / rs;
	const double n_dex_0 = ex_0 / 3.0;


	double ec_0;
	double n_dec_0;
	
	if (rs >= 1.0) {
		calc_LSDA_ge1(ec_0, n_dec_0, rs, -0.1423, 1.0529, 0.3334);

	} else {

		calc_LSDA_lt1(ec_0, n_dec_0, rs, 0.0311, -0.0480, 0.0020, -0.0116);
	}


	return { ex_0, ec_0, ex_0 + n_dex_0, ec_0 + n_dec_0 };

}


struct Return_XC_pol
{
	double E_den_x; //exchange energy density, which is not energy//
	double E_den_c; //correlation energy density, which is not energy//
	double V_x_up; //exchange energy density, which is not energy//
	double V_x_down; //exchange energy density, which is not energy//
	double V_c_up; //correlation energy density, which is not energy//
	double V_c_down; //correlation energy density, which is not energy//

};

/*
Vxcと共にExcを返す
を返す
*/
inline Return_XC_pol Calc_XC_LSDA(double rho, double zeta) {

	if (rho < 1.0e-8) {
		return { 0.0 };
	}

	static const double c2 = std::cbrt(3.0 / (4.0*M_PI));
	static const double c3 = -(3.0 / (4.0*M_PI))* std::cbrt(9.0*M_PI / 4.0);



	const double one_p_zeta = 1.0 + zeta;
	const double one_m_zeta = 1.0 - zeta;

	const double cbrt_one_p_zeta = std::cbrt(one_p_zeta);
	const double cbrt_one_m_zeta = std::cbrt(one_m_zeta);
	const double cbrt_2_m_1 = std::cbrt(2.0) - 1.0;
	const double f_zeta = (one_p_zeta * cbrt_one_p_zeta + one_m_zeta * cbrt_one_m_zeta - 2.0) / (2.0 * cbrt_2_m_1);
	const double df_dzeta = (2.0 / 3.0)*(cbrt_one_p_zeta - cbrt_one_m_zeta) / (cbrt_2_m_1);

	//exchange//

	const double rs = c2 / std::cbrt(rho);
	const double ex_0 = c3 / rs;
	const double ex_1 = std::cbrt(2.0) * ex_0;
	const double n_dex_0 = ex_0 / 3.0;
	const double n_dex_1 = ex_1 / 3.0;

	//correlation//

	double ec_0;
	double n_dec_0;
	double ec_1;
	double n_dec_1;

	if (rs >= 1.0) {
		calc_LSDA_ge1(ec_0, n_dec_0, rs, -0.1423, 1.0529, 0.3334);
		calc_LSDA_ge1(ec_1, n_dec_1, rs, -0.0843, 1.3981, 0.2611);

	} else {

		calc_LSDA_lt1(ec_0, n_dec_0, rs, 0.0311, -0.0480, 0.0020, -0.0116);
		calc_LSDA_lt1(ec_1, n_dec_1, rs, 0.01555, -0.0269, 0.0007, -0.0048);
	}
	/*
	const double diff_1_0 = ((ex_1 + ec_1) - (ex_0 + ec_0));
	const double exc_zeta = (ex_0 + ec_0) + diff_1_0 * f_zeta;
	const double dExc_dx = exc_zeta + (n_dex_0 + n_dec_0) + ((n_dex_1 + n_dec_1) - (n_dex_0 + n_dec_0)) * f_zeta;


	*Vxc_up = dExc_dx + diff_1_0 * df_dzeta * one_m_zeta;
	*Vxc_down = dExc_dx - diff_1_0 * df_dzeta * one_p_zeta;
	return exc_zeta;
	*/


	const double diff_x_1_0 = ex_1 - ex_0;
	const double diff_c_1_0 = ec_1 - ec_0;
	const double ex_zeta = ex_0 + diff_x_1_0 *f_zeta;
	const double ec_zeta = ec_0 + diff_c_1_0 *f_zeta;
	const double dEx_dx = ex_zeta + n_dex_0 + (n_dex_1 - n_dex_0) * f_zeta;
	const double dEc_dx = ec_zeta + n_dec_0 + (n_dec_1 - n_dec_0) * f_zeta;


	const double Vx_up = dEx_dx + diff_x_1_0 * df_dzeta * one_m_zeta;
	const double Vx_down = dEx_dx - diff_x_1_0 * df_dzeta * one_p_zeta;
	const double Vc_up = dEc_dx + diff_c_1_0 * df_dzeta * one_m_zeta;
	const double Vc_down = dEc_dx - diff_c_1_0 * df_dzeta * one_p_zeta;

	Return_XC_pol res;
	res.E_den_x = ex_zeta;
	res.E_den_c = ec_zeta;
	res.V_x_up = Vx_up;
	res.V_x_down = Vx_down;
	res.V_c_up = Vc_up;
	res.V_c_down = Vc_down;
	return res;
}

