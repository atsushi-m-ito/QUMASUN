#pragma once

struct MyComplex {
	double re;
	double im;
};

inline 
MyComplex operator+(const MyComplex& a, const MyComplex& b) {
	MyComplex res;
	res.re = a.re + b.re;
	res.im = a.im + b.im;
	return res;
}

inline
MyComplex operator-(const MyComplex& a, const MyComplex& b) {
	MyComplex res;
	res.re = a.re - b.re;
	res.im = a.im - b.im;
	return res;
}


inline
MyComplex operator*(const MyComplex& a, const MyComplex& b) {
	MyComplex res;
	res.re = a.re * b.re - a.im * b.im;
	res.im = a.re * b.im + a.im * b.re;
	return res;
}

inline
MyComplex inner(const MyComplex& a, const MyComplex& b) {
	MyComplex res;
	res.re = a.re * b.re + a.im * b.im;
	res.im = a.re * b.im - a.im * b.re;
	return res;
}

struct SoAComplex {
	double* re;
	double* im;
};
