#pragma once

#ifndef RealSphericalHarmonics_H
#define RealSphericalHarmonics_H

//#define _USE_MATH_DEFINES
#include <cmath>

#include "Legendre.h"
#include "Trigonometric.h"


//functor: real spherical harmonics //


class Y_s {
public:
    double operator()(const double cos_theta, const double sin_theta, const double cos_phi, const double sin_phi) {
        return sqrt(1.0/(4.0 * M_PI));
    }
	double operator()(const double ex, const double ey, const double ez) {
        return sqrt(1.0/(4.0 * M_PI));
    }
};

class Y_pz {
public:
    double operator()(const double cos_theta, const double sin_theta, const double cos_phi, const double sin_phi) {
        return sqrt(3.0/(4.0 * M_PI)) * cos_theta;
    }
	double operator()(const double ex, const double ey, const double ez) {
        return sqrt(3.0/(4.0 * M_PI)) * ez;
    }
};

class Y_px {
public:
    double operator()(const double cos_theta, const double sin_theta, const double cos_phi, const double sin_phi) {
        return sqrt(3.0/(4.0 * M_PI)) * sin_theta * cos_phi;
    }
	double operator()(const double ex, const double ey, const double ez) {
        return sqrt(3.0/(4.0 * M_PI)) * ex;
    }
};

class Y_py {
public:
    double operator()(const double cos_theta, const double sin_theta, const double cos_phi, const double sin_phi) {
        return sqrt(3.0/(4.0 * M_PI)) * sin_theta * sin_phi;
    }
	double operator()(const double ex, const double ey, const double ez) {
        return sqrt(3.0/(4.0 * M_PI)) * ey;
    }
};


class Y_dz2 {
public:
    double operator()(const double cos_theta, const double sin_theta, const double cos_phi, const double sin_phi) {
        return sqrt(5.0/(16.0 * M_PI)) * (3.0*cos_theta*cos_theta - 1.0);
    }
	double operator()(const double ex, const double ey, const double ez) {
        return sqrt(5.0/(16.0 * M_PI)) * (3.0*ez*ez - 1.0);
    }
};




//enum STATE_L_M {STATE_S, STATE_PX, STATE_PY, STATE_PZ, STATE_DZ2};


/*
template <int N>
struct Factorial
{
    enum { value = N * Factorial<N - 1>::value };
};
 
template <>
struct Factorial<0>
{
    enum { value = 1 };
};
*/


class YlmSimple{
	const int m;
	const double coef;
	BarPlm P_lm;

	double z_prev;
	double coef_P_lm_z;

	double Fact(const int n){
		double val = 1.0;
		for(int s = 2; s <= n; s++){
			val *= (double)s;
		}
		return val;
	};
	double Coefficient(const int L, const int M){
		return (M != 0) ? 
			(M > 0) ? 
				  sqrt(((double)(2*L+1)*Fact(L-M))/( 2.0 * Fact(L+M) * M_PI))
				: sqrt(((double)(2*L+1)*Fact(L+M))/( 2.0 * Fact(L-M) * M_PI))
				: sqrt((double)(2*L+1)/( 4.0 * M_PI));
	};

public:
	YlmSimple(const int L, const int M) :
		m(M), 
		coef(Coefficient(L,M) ),
		P_lm(L, std::abs(M)),
		z_prev(2.0),
		coef_P_lm_z(coef)
	{
	};

    double operator()(const double x, const double y, const double z) {
		if(z_prev == z){
			return coef_P_lm_z * BarCosSin2(m, x,y);
		}else{
			z_prev = z;
			coef_P_lm_z = coef * P_lm(z);
			return coef_P_lm_z * BarCosSin2(m, x,y);
		}
    };
};

#if 1
template <int N>
inline double fFactorial(){
	return ((double)N) * fFactorial<N-1>();
}

template <>
inline double fFactorial<0>(){
	return 1.0;
}

class BaseYlm{
public:
	virtual double operator()(const double x, const double y, const double z) = 0;
};

template <int L, int M> 
class Ylm : public BaseYlm{
public:
    double operator()(const double x, const double y, const double z) {
		const double coef = (M != 0) ? 
			(M > 0) ? sqrt((double)(2*L+1)*(double)(fFactorial<L-M>())
				/( 2.0 * (double)(fFactorial<L+M>()) * M_PI))
				: sqrt((double)(2*L+1)*(double)(fFactorial<L+M>())
				/( 2.0 * (double)(fFactorial<L-M>()) * M_PI))
				: sqrt((double)(2*L+1)/( 4.0 * M_PI));
		Plm<L,(M>=0?M:-M)> P_lm;
		CosSin<M> Tri_m;
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    };
};

BaseYlm& ReferYlm(const int L, const int M);

#endif
/*
template <> class Ylm<0, 0> {
public:
    double operator()(const double x, const double y, const double z) {
		const double coef = sqrt(1.0/(4.0 * M_PI));
		return coef;
    }
};


template <> class Ylm<1, 0> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(3.0/(4.0 * M_PI));
		return coef * z;
    }
};

template <> class Ylm<1, 1> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(3.0/(4.0 * M_PI));
		return coef * x;
    }
};

template <> class Ylm<1,-1> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(3.0/(4.0 * M_PI));
		return coef * y;
    }
};


template <> class Ylm<2, 0> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(10.0/(8.0 * M_PI));
		return coef * 0.5*(3.0*z*z -1);
    }
};

template <> class Ylm<2, 1> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(5.0/(12.0 * M_PI));
		return coef * 3.0*z*x;
    }
};

template <> class Ylm<2,-1> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(5.0/(12.0 * M_PI));
		return coef * 3.0*z*y;
    }
};

template <> class Ylm<2, 2> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(5.0/(48.0 * M_PI));
		return coef * 3.0*(x*x - y*y);
    }
};

template <> class Ylm<2,-2> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(5.0/(48.0 * M_PI));
		return coef * 3.0*(2.0*x*y);
    }
};

template <> class Ylm<3, 0> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(42.0/(24.0 * M_PI));
		return coef * (0.5*(5.0*z*z-3.0)*z);
    }
};

template <> class Ylm<3, 1> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(14.0/(48.0 * M_PI));
		return coef * (1.5*(5.0*z*z-1.0)) * (x);
    }
};

template <> class Ylm<3,-1> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(14.0/(48.0 * M_PI));
		return coef * (1.5*(5.0*z*z-1.0)) * (y);
    }
};

template <> class Ylm<3, 2> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(7.0/(240.0 * M_PI));
		return coef * (15.0*z)*(x*x - y*y);
    }
};

template <> class Ylm<3,-2> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(7.0/(240.0 * M_PI));
		return coef *(15.0*z)*(2.0*x*y);
    }
};

template <> class Ylm<3, 3> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(7.0/(1440.0 * M_PI));
		return coef * (15.0)*(x*x*x - 3.0*x*y*y);
    }
};

template <> class Ylm<3,-3> {
public:
    double operator()(const double x, const double y, const double z) {
    	const double coef = sqrt(7.0/(1440.0 * M_PI));
		return coef * (15.0)*(3.0*x*x*y - y*y*y);
    }
};


template <> class Ylm<4, 0> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_0 P_lm;
    	const double coef = sqrt(216.0/(96.0 * M_PI));
		return coef * P_lm(z, 1.0);
    }
};


template <> class Ylm<4, 1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_1 P_lm;
		CosSin<1> Tri_m;
    	const double coef = sqrt(54.0/(240.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<4,-1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_1 P_lm;
		Sin1 Tri_m;
    	const double coef = sqrt(54.0/(240.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};


template <> class Ylm<4, 2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_2 P_lm;
		CosSin<2> Tri_m;
    	const double coef = sqrt(18.0/(1440.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<4,-2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_2 P_lm;
		Sin2 Tri_m;
    	const double coef = sqrt(18.0/(1440.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<4, 3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_3 P_lm;
		CosSin<3> Tri_m;
    	const double coef = sqrt(9.0/(10080.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<4,-3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_3 P_lm;
		Sin3 Tri_m;
    	const double coef = sqrt(9.0/(10080.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<4, 4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_4 P_lm;
		CosSin<4> Tri_m;
    	const double coef = sqrt(9.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<4,-4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_4_4 P_lm;
		Sin4 Tri_m;
    	const double coef = sqrt(9.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};



template <> class Ylm<5, 0> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_0 P_lm;
    	const double coef = sqrt(1320.0/(480.0 * M_PI));
		return coef * P_lm(z, 1.0);
    }
};


template <> class Ylm<5, 1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_1 P_lm;
		CosSin<1> Tri_m;
    	const double coef = sqrt(264.0/(1440.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<5,-1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_1 P_lm;
		Sin1 Tri_m;
    	const double coef = sqrt(264.0/(1440.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};


template <> class Ylm<5, 2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_2 P_lm;
		CosSin<2> Tri_m;
    	const double coef = sqrt(66.0/(10080.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<5,-2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_2 P_lm;
		Sin2 Tri_m;
    	const double coef = sqrt(66.0/(10080.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<5, 3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_3 P_lm;
		CosSin<3> Tri_m;
    	const double coef = sqrt(22.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<5,-3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_3 P_lm;
		Sin3 Tri_m;
    	const double coef = sqrt(22.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<5, 4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_4 P_lm;
		CosSin<4> Tri_m;
    	const double coef = sqrt(11.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<5,-4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_4 P_lm;
		Sin4 Tri_m;
    	const double coef = sqrt(11.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<5, 5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_5 P_lm;
		CosSin<5> Tri_m;
    	const double coef = sqrt(11.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<5,-5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_5_5 P_lm;
		Sin5 Tri_m;
    	const double coef = sqrt(11.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<6, 0> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_0 P_lm;
    	const double coef = sqrt(9360.0/(2880.0 * M_PI));
		return coef * P_lm(z, 1.0);
    }
};


template <> class Ylm<6, 1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_1 P_lm;
		CosSin<1> Tri_m;
    	const double coef = sqrt(1560.0/(10080.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<6,-1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_1 P_lm;
		Sin1 Tri_m;
    	const double coef = sqrt(1560.0/(10080.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<6, 2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_2 P_lm;
		CosSin<2> Tri_m;
    	const double coef = sqrt(312.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6,-2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_2 P_lm;
		Sin2 Tri_m;
    	const double coef = sqrt(312.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<6, 3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_3 P_lm;
		CosSin<3> Tri_m;
    	const double coef = sqrt(78.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6,-3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_3 P_lm;
		Sin3 Tri_m;
    	const double coef = sqrt(78.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<6, 4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_4 P_lm;
		CosSin<4> Tri_m;
    	const double coef = sqrt(26.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6,-4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_4 P_lm;
		Sin4 Tri_m;
    	const double coef = sqrt(26.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6, 5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_5 P_lm;
		CosSin<5> Tri_m;
    	const double coef = sqrt(13.0/(79833600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6,-5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_5 P_lm;
		Sin5 Tri_m;
    	const double coef = sqrt(13.0/(79833600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6, 6> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_6 P_lm;
		CosSin<6> Tri_m;
    	const double coef = sqrt(13.0/(958003200.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<6,-6> {
public:
    double operator()(const double x, const double y, const double z) {
		P_6_6 P_lm;
		Sin6 Tri_m;
    	const double coef = sqrt(13.0/(958003200.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<7, 0> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_0 P_lm;
    	const double coef = sqrt(75600.0/(20160.0 * M_PI));
		return coef * P_lm(z, 1.0);
    }
};


template <> class Ylm<7, 1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_1 P_lm;
		CosSin<1> Tri_m;
    	const double coef = sqrt(10800.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<7,-1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_1 P_lm;
		Sin1 Tri_m;
    	const double coef = sqrt(10800.0/(80640.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<7, 2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_2 P_lm;
		CosSin<2> Tri_m;
    	const double coef = sqrt(1800.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7,-2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_2 P_lm;
		Sin2 Tri_m;
    	const double coef = sqrt(1800.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<7, 3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_3 P_lm;
		CosSin<3> Tri_m;
    	const double coef = sqrt(360.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7,-3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_3 P_lm;
		Sin3 Tri_m;
    	const double coef = sqrt(360.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<7, 4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_4 P_lm;
		CosSin<4> Tri_m;
    	const double coef = sqrt(90.0/(79833600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7,-4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_4 P_lm;
		Sin4 Tri_m;
    	const double coef = sqrt(90.0/(79833600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7, 5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_5 P_lm;
		CosSin<5> Tri_m;
    	const double coef = sqrt(30.0/(958003200.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7,-5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_5 P_lm;
		Sin5 Tri_m;
    	const double coef = sqrt(30.0/(958003200.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7, 6> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_6 P_lm;
		CosSin<6> Tri_m;
    	const double coef = sqrt(15.0/(12454041600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7,-6> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_6 P_lm;
		Sin6 Tri_m;
    	const double coef = sqrt(15.0/(12454041600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7, 7> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_7 P_lm;
		Cos7 Tri_m;
    	const double coef = sqrt(15.0/(174356582400.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<7,-7> {
public:
    double operator()(const double x, const double y, const double z) {
		P_7_7 P_lm;
		Sin7 Tri_m;
    	const double coef = sqrt(15.0/(174356582400.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<8, 0> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_0 P_lm;
    	const double coef = sqrt(685440.0/(161280.0 * M_PI));
		return coef * P_lm(z, 1.0);
    }
};


template <> class Ylm<8, 1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_1 P_lm;
		CosSin<1> Tri_m;
    	const double coef = sqrt(85680.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<8,-1> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_1 P_lm;
		Sin1 Tri_m;
    	const double coef = sqrt(85680.0/(725760.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x,y);
    }
};

template <> class Ylm<8, 2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_2 P_lm;
		CosSin<2> Tri_m;
    	const double coef = sqrt(12240.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-2> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_2 P_lm;
		Sin2 Tri_m;
    	const double coef = sqrt(12240.0/(7257600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<8, 3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_3 P_lm;
		CosSin<3> Tri_m;
    	const double coef = sqrt(2040.0/(79833600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-3> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_3 P_lm;
		Sin3 Tri_m;
    	const double coef = sqrt(2040.0/(79833600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};


template <> class Ylm<8, 4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_4 P_lm;
		CosSin<4> Tri_m;
    	const double coef = sqrt(408.0/(958003200.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-4> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_4 P_lm;
		Sin4 Tri_m;
    	const double coef = sqrt(408.0/(958003200.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8, 5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_5 P_lm;
		CosSin<5> Tri_m;
    	const double coef = sqrt(102.0/(12454041600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-5> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_5 P_lm;
		Sin5 Tri_m;
    	const double coef = sqrt(102.0/(12454041600.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8, 6> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_6 P_lm;
		CosSin<6> Tri_m;
    	const double coef = sqrt(34.0/(174356582400.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-6> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_6 P_lm;
		Sin6 Tri_m;
    	const double coef = sqrt(34.0/(174356582400.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8, 7> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_7 P_lm;
		Cos7 Tri_m;
    	const double coef = sqrt(17.0/(2615348736000.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-7> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_7 P_lm;
		Sin7 Tri_m;
    	const double coef = sqrt(17.0/(2615348736000.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8, 8> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_8 P_lm;
		Cos8 Tri_m;
    	const double coef = sqrt(17.0/(41845579776000.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};

template <> class Ylm<8,-8> {
public:
    double operator()(const double x, const double y, const double z) {
		P_8_8 P_lm;
		Sin8 Tri_m;
    	const double coef = sqrt(17.0/(41845579776000.0 * M_PI));
		return coef * P_lm(z, 1.0) * Tri_m(x, y);
    }
};
*/


/*==========================================
 to check
===========================================*/

template <typename SphericalHarmonics>
double IntegralSphericalHarmonics(SphericalHarmonics Y_l_m, const int split_theta, const int split_phi){
	const double delta_theta = M_PI / ((double)split_theta);
	const double delta_phi = 2.0*M_PI / ((double)split_phi);
	double sum = 0.0;
	for(int t = 0; t < split_theta; t++){
		const double theta = delta_theta * ((double)t);
		const double z = cos(theta);
		const double sin_theta = sin(theta);
		
		for(int i = 0; i < split_phi; i++){
			const double phi = delta_phi * ((double)i);
			const double cos_phi = cos(phi);
			const double sin_phi = sin(phi);
			const double x = sin_theta * cos_phi;
			const double y = sin_theta * sin_phi;
			const double val = Y_l_m(x, y, z);
			sum += val * sin_theta * delta_theta * delta_phi;
		}
	}
	return sum;
}

/*
if result is 1.0, it is success.
*/
template <typename SphericalHarmonics>
double IntegralSphericalHarmonicsSquare(SphericalHarmonics Y_l_m, const int split_theta, const int split_phi){
	const double delta_theta = M_PI / ((double)split_theta);
	const double delta_phi = 2.0*M_PI / ((double)split_phi);
	double sum = 0.0;
	for(int t = 0; t < split_theta; t++){
		const double theta = delta_theta * ((double)t);
		const double z = cos(theta);
		const double sin_theta = sin(theta);
			
		for(int i = 0; i < split_phi; i++){
			const double phi = delta_phi * ((double)i);
			const double cos_phi = cos(phi);
			const double sin_phi = sin(phi);
			const double x = sin_theta * cos_phi;
			const double y = sin_theta * sin_phi;
			const double val = Y_l_m(x, y, z);
			sum += val * val* sin_theta * delta_theta * delta_phi;
		}
	}
	return sum;
}



#endif	//! RealSphericalHarmonics_H
