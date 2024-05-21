#ifndef Legendre_H
#define Legendre_H


#include <cmath>


//functor: Associated Legendre polynomials function //

/*===========================================================================
	ルジャンドル陪関数は次の関係式を満たす。
	(1)  P_{ll}(x) =  (2l-1)!! (1-x^2)^{l/2}
	     ただし (～)!! は一つとばしの階乗
	(2)  P_{l+1,l}(x) = (2l+1) x P_{ll}(x)
	(3)  (l-m+1)P_{l+1,m}(x) = (2l+1) x P_{lm}(x) - (l+m)P_{l-1,m}(x)

	そこで、P_{lm}の計算の実行では、最初に(1)でP_{mm}を求め、次に(2)でP_{m+1,m}を求める。
	そこから漸化式(3)を用いて逐次的にP_{lm}まで計算する。

	また、実際に球面調和関数で欲しい値は、
	(4)  P_{lm}(x) = (1-x^2)^{m/2} P~_{lm}(x)
	を満たす多項式P~_{lm}である。
    式(2),(3)は左辺右辺で全てmが揃ている為、P~lmも式(2),(3)と同様の式に従う。
	よって、初期値P~_{mm}を(1)式の係数[(2m-1)!!]だけにすることで((1-x^2)^{m/2}を除くことで)、
    P~_{lm}も(2)式から逐次的に求められる。
    残った(1-x^2)^{m/2}=\sin^m(\theta)の部分は、球面調和関数では\exp(i m \phi)と合わせて
    \sin^m(\theta)\exp(i m \phi) = (\cos\phi + i \sin\phi)^m = (x+iy)^m
    where x and y is normalized, となる。
    つまり、デカルト座標(x,y,z)を引数にとって使うには、
    P_{lm}よりもP~_{lm}の方が使いやすい。

    以下の実装ではP~_{lm}をclass BarPlmとして実装する。

以下実装=================================================================*/


class BarPlm{
private:
	const int l;
	const int m;
	const double P_mm;
	/*
	double Fact2(const int n){
		double val = 1.0;
		for(int s = n; s > 1; s-=2){
			val *= (double)s;
		}
		return val;
	}
	*/
	
	double Fact2(const int n){
		double val = 1.0;
		for(int s = n; s > 1; s-=2){
			val *= (double)s;
		}
		return val;
	}
public:
	BarPlm(const int L, const int M) :
		l(L), m(M), P_mm(Fact2(2*M-1))
	{
		if(l==m){
			int a=m;
		}
	};

	double operator()(const double x){
		if(l==m){ return P_mm;}
		/*
		double p1 = P_mm;
		double p0 = P_mm * x * (double)(2*m+1);
		if(l==m+1){ return p0;}
		*/
		double p1 = 0.0;
		double p0 = P_mm;

		for(int s = m; s < l; s++){
			double p2 = p1;
			p1 = p0;
			p0 = ( + ((double)(2*s+1)) * x * p1 - ((double)(s+m)) * p2) / ((double)(s-m+1));
		}
	
		return p0;

	};
};





//以下はテンプレートの特殊化で直接書き下したもの======================================//
// z = cos(theta)
// s = sin(theta)
template <int L, int M>
class Plm {
public:
    double operator()(const double z, const double s){
        return 1.0;
    }
};

template <> class Plm<0, 0> {
public:
    double operator()(const double z, const double s){
        return 1.0;
    }
};

template <> class Plm<1, 0> {
public:
    double operator()(const double z, const double s){
        return z;
    }
};

template <> class Plm<1, 1> {
public:
    double operator()(const double z, const double s){
        return s;
    }
};

template <> class Plm<2,0> {
public:
    double operator()(const double z, const double s){
        return 0.5*(3.0*z*z -1);
    }
};

template <> class Plm<2,1> {
public:
    double operator()(const double z, const double s){
        return 3.0*s*z;
    }
};

template <> class Plm<2,2> {
public:
    double operator()(const double z, const double s){
        return 3.0*s*s;
    }
};

template <> class Plm<3,0> {
public:
    double operator()(const double z, const double s){
        return 0.5*(5.0*z*z-3.0)*z;
    }
};

template <> class Plm<3,1> {
public:
    double operator()(const double z, const double s){
        return 1.5*s*(5.0*z*z-1.0);
    }
};
template <> class Plm<3,2> {
public:
    double operator()(const double z, const double s){
        return 15.0*s*s*z;
    }
};
template <> class Plm<3,3> {
public:
    double operator()(const double z, const double s){
        return 15.0*s*s*s;
    }
};

template <> class Plm<4,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		return 0.125*(35.0*z2*z2-30.0*z2+3.0);
    }
};

template <> class Plm<4,1> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		return 2.5*s*z*(7.0*z2-3.0);
    }
};
template <> class Plm<4,2> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		return 7.5*s*s*(7.0*z2-1.0);
    }
};
template <> class Plm<4,3> {
public:
    double operator()(const double z, const double s){
		return 105.0*s*s*s*z;
    }
};
template <> class Plm<4,4> {
public:
    double operator()(const double z, const double s){
		return 105.0*s*s*s*s;
    }
};

template <> class Plm<5,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return 0.125*(63.0*z2*z2-70.0*z2+15.0)*z;
    }
};


template <> class Plm<5,1> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (15.0/8.0)*s*(21.0*z2*z2-14.0*z2+1.0);
    }
};
template <> class Plm<5,2> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (105.0/2.0)*s*s*z*(3.0*z2-1.0);
    }
};
template <> class Plm<5,3> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (105.0/2.0)*s*s*s*(9.0*z2-1.0);
    }
};
template <> class Plm<5,4> {
public:
    double operator()(const double z, const double s){	
        return 945.0*s*s*s*s*z;
    }
};
template <> class Plm<5,5> {
public:
    double operator()(const double z, const double s){	
        return 945.0*s*s*s*s*s;
    }
};

template <> class Plm<6,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return 0.0625*(231.0*z4*z2-315.0*z4+105.0*z2-5.0);
    }
};

template <> class Plm<6,1> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (21.0/8.0)*s*z*(33.0*z4-30.0*z2+5.0);
    }
};
template <> class Plm<6,2> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (105.0/8.0)*s*s*(33.0*z4-18.0*z2+1.0);
    }
};
template <> class Plm<6,3> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		return (315.0/2.0)*s*s*s*z*(11.0*z2-3.0);
    }
};
template <> class Plm<6,4> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		return (945/2.0)*s*s*s*s*(11.0*z2-1.0);
    }
};
template <> class Plm<6,5> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		return 10395.0*s*s*s*s*s*z;
    }
};
template <> class Plm<6,6> {
public:
    double operator()(const double z, const double s){
		return 10395.0*s*s*s*s*s*s;
    }
};


template <> class Plm<7,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return 0.0625*(429.0*z4*z2-693.0*z4+315.0*z2-35.0)*z;
    }
};

template <> class Plm<7,1> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (7.0/16.0)*s*(429.0*z4*z2-495.0*z4+135.0*z2-5.0);
    }
};
template <> class Plm<7,2> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (63.0/8.0)*s*s*z*(143.0*z4-110.0*z2+15.0);
    }
};
template <> class Plm<7,3> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (315.0/8.0)*s*s*s*(143.0*z4-66.0*z2+3.0);
    }
};
template <> class Plm<7,4> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (3465.0/2.0)*s*s*s*s*z*(13.0*z2-3.0);
    }
};
template <> class Plm<7,5> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (10395.0/2.0)*s*s*s*s*s*(13.0*z2-1.0);
    }
};
template <> class Plm<7,6> {
public:
    double operator()(const double z, const double s){
		return 135135.0*s*s*s*s*s*s*z;
    }
};
template <> class Plm<7,7> {
public:
    double operator()(const double z, const double s){
		return 135135.0*s*s*s*s*s*s*s;
    }
};

template <> class Plm<8,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        const double z6 = z4*z2;
		const double z8 = z4*z4;
        return 0.0078125*(6435.0*z8-12012*z6+6930.0*z4-1260.0*z2+35.0);
    }
};

template <> class Plm<8,1> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        const double z6 = z4*z2;
        return (9.0/16.0)*s*z*(715.0*z6-1001.0*z4+385.0*z2-35.0);
    }
};
template <> class Plm<8,2> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        const double z6 = z4*z2;
        return (315.0/16.0)*s*s*(143.0*z6-143.0*z4+33.0*z2-1.0);
    }
};
template <> class Plm<8,3> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (3465.0/8.0)*s*s*s*z*(39.0*z4-26.0*z2+3.0);
    }
};
template <> class Plm<8,4> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        return (10395.0/8.0)*s*s*s*s*(65.0*z4-26.0*z2+1.0);
    }
};
template <> class Plm<8,5> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (135135.0/2.0)*s*s*s*s*s*z*(5.0*z2-1.0);
    }
};
template <> class Plm<8,6> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
        return (135135.0/2.0)*s*s*s*s*s*s*(15.0*z2-1.0);
    }
};
template <> class Plm<8,7> {
public:
    double operator()(const double z, const double s){
		return (2027025.0)*s*s*s*s*s*s*s*z;
    }
};
template <> class Plm<8,8> {
public:
    double operator()(const double z, const double s){
		return (2027025.0)*s*s*s*s*s*s*s*s;
    }
};


template <> class Plm<9,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        const double z6 = z4*z2;
		const double z8 = z4*z4;
        return 0.0078125*(12155.0*z8-25740*z6+18018*z4-4620*z2+315.0)*z;
    }
};

template <> class Plm<10,0> {
public:
    double operator()(const double z, const double s){
		const double z2 = z*z;
		const double z4 = z2*z2;
        const double z6 = z4*z2;
		const double z8 = z4*z4;
        return 0.00390625*(46189.0*z8*z2-109395.0*z8 + 90090.0*z6-30030*z4+3465*z2-63)*z;
    }
};


#endif	//! Legendre_H
