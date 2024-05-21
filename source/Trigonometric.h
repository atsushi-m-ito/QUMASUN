#ifndef Trigonometric_H
#define Trigonometric_H

//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>

/*
 cos(m phi) if m >= 0,
 sin(-m phi) if m < 0.
*/
inline double BarCosSin(const int m, const double x, const double y){
	double re = 1.0;
	double im = 0.0;

	const int s_end = std::abs(m);
	for(int s = 0; s < s_end; s++){
		double re_next = re * x - im * y;
		double im_next = re * y + im * x;
		re = re_next;
		im = im_next;
	}

	return (m >= 0) ? re : im;
}

/*
 When \cos\theta = x, \sin\theta = y are given by input argument,
 the rotated value
 cos(m phi) if m >= 0,
 sin(-m phi) if m < 0,
 is estimated as the function of x and y.
*/
inline double BarCosSin2(const int m, const double x, const double y){
	if(m==0){
		return 1.0;
	}else{
		double c0 = 1.0;
		double c_next = x;

		if(m > 0){  //cos(m phi)//
			c0 = 1.0;
			c_next = x;
		}else{		//sin(m phi)//
			c0 = 0.0;
			c_next = y;
		}

		const double xx_yy= x*x + y*y;
		const double x2 = 2.0*x;


		const int s_end = std::abs(m);
		for(int s = 2; s <= s_end; s++){
			double c_prev = c0;
			c0 = c_next;
			c_next = x2 * c0 - xx_yy * c_prev;
		}
		return c_next;
	}
}




// x = cos_theta
// y = sin_theta


template <int M> 
class CosSin {
public:
    double operator()(const double x, const double y){
        return BarCosSin2(M,x,y);
    }
};

template <> class CosSin<0> {
public:
    double operator()(const double x, const double y){
        return 1.0;
    }
};

template <> class CosSin<1> {
public:
    double operator()(const double x, const double y){
        return x;
    }
};

template <> class CosSin<-1> {
public:
    double operator()(const double x, const double y){
        return y;
    }
};


template <> class CosSin<2> {
public:
    double operator()(const double x, const double y){
        return x*x - y*y;
    }
};

template <> class CosSin<-2> {
public:
    double operator()(const double x, const double y){
        return 2.0*x*y;
    }
};

template <> class CosSin<3> {
public:
    double operator()(const double x, const double y){
        return x*x*x - 3.0*x*y*y;
    }
};

template <> class CosSin<-3> {
public:
    double operator()(const double x, const double y){
        return 3.0*x*x*y - y*y*y;
    }
};

template <> class CosSin<4> {
public:
    double operator()(const double x, const double y){
        return x*x*x*x - 6.0*x*x*y*y + y*y*y*y;
    }
};

template <> class CosSin<-4> {
public:
    double operator()(const double x, const double y){
        return 4.0*(x*x*x*y - x*y*y*y);
    }
};


template <> class CosSin<5> {
public:
    double operator()(const double x, const double y){
        return x*x*x*x*x - 10.0*x*x*x*y*y + 5.0*x*y*y*y*y;
    }
};

template <> class CosSin<-5> {
public:
    double operator()(const double x, const double y){
        return 5.0*x*x*x*x*y - 10.0*x*x*y*y*y + y*y*y*y*y;
    }
};


template <> class CosSin<6> {
public:
    double operator()(const double x, const double y){
		const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
        return x6 - 15.0*x4*y2 + 15.0*x2*y4 - y6;
    }
};

template <> class CosSin<-6> {
public:
    double operator()(const double x, const double y){
        const double x3 = x*x*x;
		const double x5 = x3*x*x;
		const double y3 = y*y*y;
		const double y5 = y3*y*y;
        return 6.0*x5*y - 20.0*x3*y3 + 6.0*x*y5;
    }
};


template <> class CosSin<7> {
public:
    double operator()(const double x, const double y){
		const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
        return x*(x6 - 21.0*x4*y2 + 35.0*x2*y4 - 7.0*y6);
    }
};

template <> class CosSin<-7> {
public:
    double operator()(const double x, const double y){
        const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
        return (7.0*x6 - 35.0*x4*y2 + 21.0*x2*y4 - y6)*y;
    }
};

template <> class CosSin<8> {
public:
    double operator()(const double x, const double y){
		const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
		const double x8 = x4*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
		const double y8 = y4*y4;
        return x8 - 28.0*x6*y2 + 70.0*x4*y4 - 28.0*x2*y6 + y8;
    }
};

template <> class CosSin<-8> {
public:
    double operator()(const double x, const double y){
        const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
        return x*y*(8.0*x6 - 56.0*x4*y2 + 56.0*x2*y4 - 8.0*y6);
    }
};

template <> class CosSin<9> {
public:
    double operator()(const double x, const double y){
		const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
		const double x8 = x4*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
		const double y8 = y4*y4;
        return x*(x8 - 36.0*x6*y2 + 126.0*x4*y4 - 84.0*x2*y6 + 9.0*y8);
    }
};

template <> class CosSin<-9> {
public:
    double operator()(const double x, const double y){
        const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
        const double x8 = x4*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
        const double y8 = y4*y4;
        return y*(9.0*x8 - 84.0*x6*y2 + 126.0*x4*y4 - 36.0*x2*y6 + y8);
    }
};


template <> class CosSin<10> {
public:
    double operator()(const double x, const double y){
		const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
		const double x8 = x4*x4;
        const double x10 = x4*x6;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
		const double y8 = y4*y4;
		const double y10 = y4*y6;
        return x10 - 45.0*x8*y2 + 210.0*x6*y4 - 210.0*x4*y6 + 45.0*x2*y8 - y10;
    }
};

template <> class CosSin<-10> {
public:
    double operator()(const double x, const double y){
        const double x2 = x*x;
		const double x4 = x2*x2;
		const double x6 = x2*x4;
        const double x8 = x4*x4;
        const double y2 = y*y;
		const double y4 = y2*y2;
		const double y6 = y2*y4;
        const double y8 = y4*y4;
        return x*y*(10.0*x8 - 120.0*x6*y2 + 252.0*x4*y4 - 120.0*x2*y6 + 10.0*y8);
    }
};

/*==========================================
 to check
===========================================*/
/*
if result is 0.0, it is success.
*/
template <typename TrigonometricFunc>
inline double IntegralTrigonometric(TrigonometricFunc func, const int split){
	const double delta = 2.0*M_PI / ((double)split);
	double sum = 0.0;
	for(int i = 0; i < split; i++){
		const double angle = delta * ((double)i);
		const double cos_phi = cos(angle);
		const double sin_phi = sin(angle);
		const double val = func(cos_phi, sin_phi);
		sum += val * delta;
	}
	return sum;
}

/*
if result is M_PI, it is success.
*/
template <typename TrigonometricFunc>
inline double IntegralTrigonometricSquare(TrigonometricFunc func, const int split){
	const double delta = 2.0*M_PI / ((double)split);
	double sum = 0.0;
	for(int i = 0; i < split; i++){
		const double angle = delta * ((double)i);
		const double cos_phi = cos(angle);
		const double sin_phi = sin(angle);
		const double val = func(cos_phi, sin_phi);
		sum += val * val * delta;
	}
	return sum;
}

#endif	//!Trigonometric_H

