#pragma once


struct SoAComplex {
	double* re;
	double* im;
};

struct OneComplex {
	double r;
	double i;
};

double real(const OneComplex& a) {
    return a.r;
}

double imag(const OneComplex& a) {
    return a.i;
}

namespace SoAC {
	inline
	double Norm(const SoAComplex& v1, size_t N) {
		double xp_re = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1.re[i] * v1.re[i]
				   + v1.im[i] * v1.im[i];
		}

		return xp_re;
	}


	inline
	void MulC(SoAComplex& v1, double c, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1.re[i] *= c;
			v1.im[i] *= c;
		}
	}

	inline
	OneComplex InnerProd(const SoAComplex& v1, const SoAComplex& v2, size_t N) {
		double xp_re = 0.0;
		double xp_im = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1.re[i] * v2.re[i] + v1.im[i] * v2.im[i];
			xp_im += v1.re[i] * v2.im[i] - v1.im[i] * v2.re[i];
		}
		return { xp_re, xp_im };
	}


	inline
	double InnerProdReal(const SoAComplex& v1, const SoAComplex& v2, size_t N) {
		double xp_re = 0.0;
		//double xp_im = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1.re[i] * v2.re[i] + v1.im[i] * v2.im[i];
			//xp_im += v1.re[i] * v2.im[i] - v1.im[i] * v2.re[i];
		}
		return xp_re;
	}


	inline
	void SetZero(SoAComplex& v1, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1.re[i] = 0.0;
		}
		for (size_t i = 0; i < N; ++i) {
			v1.im[i] = 0.0;
		}
	}

	inline
	void AddV(SoAComplex& v1, const SoAComplex& v2, const OneComplex c, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1.re[i] += v2.re[i] * c.r - v2.im[i] * c.i;
			v1.im[i] += v2.re[i] * c.i + v2.im[i] * c.r;
		}
	}

	inline
	void Copy(SoAComplex& v1, const SoAComplex& v2, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1.re[i] = v2.re[i];
		}
		for (size_t i = 0; i < N; ++i) {
			v1.im[i] = v2.im[i];
		}
	}
}
