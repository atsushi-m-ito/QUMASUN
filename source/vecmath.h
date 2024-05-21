#pragma once
#ifndef NO_COMPLEX
#include <complex>
#endif

namespace vecmath {

	template<class T>
	inline
	T InnerProd(const T* v1, const T* v2, size_t N) {
		double xp_re = 0.0;
		
		for (size_t i = 0; i < N; ++i) {
			xp_re += v1[i] * v2[i];
		}
		return xp_re;
	}

	
	template<class T>
	inline
	double Norm(const T* v1, size_t N) {
		double xp_re = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1[i] * v1[i];
		}

		return xp_re;
	}

#ifndef NO_COMPLEX

	template<>
	inline
	std::complex<double> InnerProd<std::complex<double>>(const std::complex<double>* v1, const std::complex<double>* v2, size_t N) {
		double xp_re = 0.0;
		double xp_im = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1[i].real() * v2[i].real() + v1[i].imag() * v2[i].imag();
			xp_im += v1[i].real() * v2[i].imag() - v1[i].imag() * v2[i].real();
		}
		return { xp_re, xp_im };
	}

	inline
	std::complex<double> Angle(const std::complex<double>* v1, const std::complex<double>* v2, size_t N) {
		double xp_re = 0.0;
		double xp_im = 0.0;
		double norm1 = 0.0;
		double norm2 = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1[i].real() * v2[i].real() + v1[i].imag() * v2[i].imag();
			xp_im += v1[i].real() * v2[i].imag() - v1[i].imag() * v2[i].real();

			norm1 += v1[i].real() * v1[i].real() + v1[i].imag() * v1[i].imag();
			norm2 += v2[i].real() * v2[i].real() + v2[i].imag() * v2[i].imag();
		}
		return { xp_re / sqrt(norm1 * norm2), xp_im / sqrt(norm1 * norm2) };
	}

	template<>
	inline
	double Norm<std::complex<double>>(const std::complex<double>* v1, size_t N) {
		double xp_re = 0.0;

		for (size_t i = 0; i < N; ++i) {
			xp_re += v1[i].real() * v1[i].real() + v1[i].imag() * v1[i].imag();
		}

		return xp_re;
	}
#endif


	template <class T>
	inline
	void SetZero(T* v1, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1[i] = 0.0;
		}
	}


	template <class T>
	inline
	void Copy(T* v1, const T* v2, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1[i] = v2[i];
		}
	}

	template <class T>
	inline
	void AddNorm(double* dst, const T* v, const double coef, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			dst[i] += (v[i]* v[i]) * coef;
		}
	}

#ifndef NO_COMPLEX

	template <>
	inline
	void AddNorm< std::complex<double> > (double* dst, const std::complex<double>* v, const double coef, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			dst[i] += std::norm(v[i]) * coef;
		}
	}

#endif



	template <class T>
	inline
	void MulC(T* v1, double c, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1[i] *= c;
		}
	}


	template <class T>
	inline
	double Normalize(T* v1, size_t N, double dVol) {
		const double xp_re = Norm(v1, N) * dVol;
		const double c = 1.0 / sqrt(xp_re);
		MulC(v1, c, N);
		return xp_re;
	}

	template <class T>
	inline
	void AddV(T* v1, const T* v2, double c, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1[i] += v2[i]*c;
		}
	}

	/*
	template <class T>
	inline
	void AddTo(T* v1, const T* v2, const T* v3, size_t N) {
		for (size_t i = 0; i < N; ++i) {
			v1[i] = v2[i] + v3[i];
		}
	}
	*/
}
