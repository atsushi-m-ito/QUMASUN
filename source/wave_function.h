#pragma once
#include <cstdint>
#include <complex>


template <typename T>
class GridFunction {
private:
	T* data;
	size_t size;
public:
	GridFunction() :
		size(0),
		data(nullptr) {
	}

	GridFunction(size_t N) :
		size(N),
		data(new T[N])
	{
	}

	GridFunction(T* reference) :
		size(0),
		data(reference)
	{
	}

	GridFunction(const GridFunction& src) = delete;
	GridFunction& operator=(const GridFunction& src) = delete;

	virtual ~GridFunction<T>() {
		if (size > 0) {
			delete[] data;
		}
	}

	template<typename I>
	const T& operator[](I index) const {
		return data[index];
	}

	template<typename I>
	T& operator[](I index) {
		return data[index];
	}

	const T* Pointer() const {
		return data;
	}

	T* Pointer() {
		return data;
	}

	void Renew(size_t N) {
		if (size > 0) {
			delete[] data;
		}
		size = N;
		data = new T[N];		
	}

	void Refer(T* p) {
		if (size > 0) {
			delete[] data;
		}
		size = 0;
		data = p;
	}

	operator T* () {
		return data;
	}

	operator const T* () const {
		return data;
	}
};



struct GridInfo {
	size_t size_3d;
	int size_x;
	int size_y;
	int size_z;
	/*double dx;
	double dy;
	double dz;

	double BoxX() { return dx * (double)size_x; }
	double BoxY() { return dy * (double)size_y; }
	double BoxZ() { return dz * (double)size_z; }
	*/
};

struct SubgridInfo {
	size_t size_3d;
	int size_x;
	int size_y;
	int size_z;
	int begin_x;
	int begin_y;
	int begin_z;

};

template<typename T>
class RspaceFunc : public GridFunction< T >
{
public:
	RspaceFunc() : GridFunction<T>() {}
	RspaceFunc(size_t N) : GridFunction<T>(N){	}
	RspaceFunc(T* reference) : GridFunction<T>(reference) {}
};

template<typename T>
class KspaceFunc : public GridFunction< T >
{
public:
	KspaceFunc() : GridFunction<T>() {}
	KspaceFunc(size_t N) : GridFunction<T>(N) {}
	KspaceFunc(T* reference) : GridFunction<T>(reference) {}
};



