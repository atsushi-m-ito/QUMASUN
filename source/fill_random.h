#pragma once
#include <complex>
#include <random>
#include "mycomplex.h"

inline
void FillRandom(int N, double* out, std::mt19937& mt) {
    std::uniform_real_distribution<> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        out[i] = distribution(mt);
    }
}

inline
void FillRandom(int N, std::complex<double>* out, std::mt19937& mt) {
    std::uniform_real_distribution<> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        out[i].real( distribution(mt)) ;
        out[i].imag(distribution(mt));
    }
}


inline
void FillRandom(int N, MyComplex* out, std::mt19937& mt) {
    std::uniform_real_distribution<> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        out[i].re=(distribution(mt));
        out[i].im=(distribution(mt));
    }
}


inline
void FillRandom(int N, SoAComplex& out, std::mt19937& mt) {
    std::uniform_real_distribution<> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        out.re[i] = (distribution(mt));
        out.im[i] = (distribution(mt));
    }
}

