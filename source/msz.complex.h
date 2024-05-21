#pragma once
#include <complex>

namespace msz {
	namespace complex {
		
		inline
		std::complex<double> inner(std::complex<double>& a, std::complex<double>& b) {
			return { a.real() * b.real() + a.imag() * b.imag(),
				-a.imag() * b.real() + a.real() * b.imag() };
		}


		inline
		void omp_atomic_add(std::complex<double>& a, const std::complex<double>& b) {
			struct d2
			{
				double re; double im;
			};
			d2& ad = *(d2*)&a;

#pragma omp atomic
			ad.re += b.real();

#pragma omp atomic
			ad.im += b.imag();
		}
	}
}
