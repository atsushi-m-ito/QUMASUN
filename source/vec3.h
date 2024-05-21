//公開Level1//

#pragma once
#ifndef VEC3_H
#define VEC3_H

#include <cmath>


template<typename T>
struct vec3 {
public:
	T x;
	T y;
	T z;

	
	vec3() = default;  //for trivial class and POD//
	vec3(const T& x0, const T& y0, const T& z0) : x(x0), y(y0), z(z0) {};
	/*vec3(const T& x0) : x(x0), y(x0), z(x0) {};
	This should not be defined, because even if vec3<T> += double, there will automatically be a cast from double to vec3<T> and it will pass.*/

	inline void Set(const T& d1, const T& d2, const T& d3) {
		x = d1;
		y = d2;
		z = d3;
	}
	
	inline void Clear() {
		x = 0.0;
		y = 0.0;
		z = 0.0;

	}

	template<typename U>
	inline operator vec3<U>() {
		return vec3<U>(static_cast<U>(x), static_cast<U>(y), static_cast<U>(z));
	}

	inline vec3 Normalize() {
		const T ni = 1.0 / std::sqrt(x*x + y*y + z*z);
		*this *= ni;
		return *this;
	}


	inline vec3 operator +() const {
		return *this;
	}

	inline vec3 operator -() const {
		return vec3(-x, -y, -z);
	}


	inline vec3& operator +=(const vec3& v2) {
		x += v2.x;
		y += v2.y;
		z += v2.z;
		return *this;
	}


	inline vec3& operator -=(const vec3& v2) {
		x -= v2.x;
		y -= v2.y;
		z -= v2.z;
		return *this;
	}

	inline vec3& operator *=(const T& d) {
		x *= d;
		y *= d;
		z *= d;
		return *this;
	}

	inline vec3& operator /=(const T& d) {
		T di = 1.0 / d;
		x *= di;
		y *= di;
		z *= di;
		return *this;
	}
};


template<typename T>
inline T operator *(const vec3<T>& v1, const vec3<T>& v2) {
	return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

template<typename T>
inline vec3<T> operator +(const vec3<T>& v1, const vec3<T>& v2) {
	//return vec3<T>(v1) += v2;//FUJITSU cannot optimize//
	return vec3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template<typename T>
inline vec3<T> operator -(const vec3<T>& v1, const vec3<T>& v2) {
	//return vec3<T>(v1) -= v2;//FUJITSU cannot optimize//
	return vec3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

template<typename T>
inline vec3<T> operator *(const vec3<T>& v, const T& d) {
	//return vec3<T>(v) *= d;//FUJITSU cannot optimize//
	return vec3<T>(v.x * d, v.y * d, v.z * d);
}

template<typename T>
inline vec3<T> operator *(const T& d, const vec3<T>& v) {
	//return vec3<T>(v) *= d;//FUJITSU cannot optimize//
	return vec3<T>(v.x * d, v.y * d, v.z * d);
}

template<typename T>
inline vec3<T> operator /(const vec3<T>& v, const T& d) {
	//return vec3<T>(v) /= d;//FUJITSU cannot optimize//
	const T di = 1.0 / d;
	return vec3<T>(v.x * di, v.y * di, v.z * di);
}


template<typename T>
inline vec3<T> Unit(const vec3<T>& v) {
	const T ni = 1.0 / std::sqrt(v*v);
	return vec3<T>(v * ni);
}

template<typename T>
inline vec3<T> Cross(const vec3<T>& v1, const vec3<T>& v2) {
	return vec3<T>(
		v1.y * v2.z - v1.z * v2.y,
		v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x
		);
}

template<class T>
inline T Abs(const vec3<T>& v) {
	return sqrt(v*v);
}


typedef vec3<double> vec3d;
typedef vec3<float> vec3f;

#endif // VEC3_H

