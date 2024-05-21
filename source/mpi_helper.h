#pragma once
#include <mpi.h>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>

//MPI_Datatype is estimated from template parameter
template <class T>
inline
MPI_Datatype GetDatatype() {
	if constexpr (sizeof(T) == 16) {
		return MPI_DOUBLE_COMPLEX;
	} else {
		return MPI_DATATYPE_NULL;
	}

}

template <>
inline
MPI_Datatype GetDatatype<int>() {return MPI_INT;}

template <>
inline
MPI_Datatype GetDatatype<int64_t>() { return MPI_INT64_T; }

template <>
inline
MPI_Datatype GetDatatype<uint64_t>() { return MPI_UINT64_T; }

template <>
inline
MPI_Datatype GetDatatype<double>() { return MPI_DOUBLE; }

template <>
inline
MPI_Datatype GetDatatype<std::complex<double> >() { return MPI_C_DOUBLE_COMPLEX; }

inline
int GetProcessID(const MPI_Comm& mpi_comm) {
	int proc_id;
	MPI_Comm_rank(mpi_comm, &proc_id);
	return proc_id;
}

inline
int GetNumProcess(const MPI_Comm& mpi_comm) {
	int num_procs;
	MPI_Comm_size(mpi_comm, &num_procs);
	return num_procs;
}

//check root or not
inline 
bool IsRoot(const MPI_Comm& mpi_comm) {
	constexpr int root_id = 0;
	int proc_id = GetProcessID(mpi_comm);
	return (proc_id == root_id);
}
/*
template<class FORMAT, class ...Arg>
void IfRootPrint(const MPI_Comm& mpi_comm, FORMAT format, Arg ...arg) {
	if (IsRoot(mpi_comm)) {
		printf(format, arg...);
		fflush(stdout);
	}
}
*/

template<class HEAD, class ...Arg>
constexpr size_t SizeofAny() {
	if constexpr (sizeof...(Arg) == 0) {
		return sizeof(HEAD);
	} else {
		return sizeof(HEAD) + SizeofAny<Arg...>();
	}
}


template<class HEAD, class ...Arg>
constexpr void CopyFromAny(uint8_t* data, const HEAD& head, const Arg& ...arg) {
	*(HEAD*)data = head;

	if constexpr (sizeof...(Arg) > 0) {
		CopyFromAny(data + sizeof(HEAD), arg...);
	}
}

template<class HEAD, class ...Arg>
constexpr void CopyToAny(uint8_t* data, HEAD& head, Arg& ...arg) {
	head = *(HEAD*)data;
	if constexpr (sizeof...(Arg) > 0) {
		CopyToAny(data + sizeof(HEAD), arg...);
	}
}




template<class HEAD, class ...Arg>
constexpr void TestAny(uint8_t* data, const HEAD& head, const Arg& ...arg) {
	if constexpr (std::is_same_v<HEAD, double>) {
		printf("%f\n", *(HEAD*)data);
	} else if constexpr (std::is_same_v<HEAD, int>) {
		printf("%d\n", *(HEAD*)data);
	} else if constexpr (std::is_same_v<HEAD, int64_t>) {
		printf("%zd\n", *(HEAD*)data);
	}
	if constexpr (sizeof...(Arg) > 0) {
		TestAny(data + sizeof(HEAD), arg...);
	}
}



#if 1

template<class ...Arg>
void BroadcastAny(const MPI_Comm& mpi_comm, int root_id, Arg& ...arg) {
	constexpr size_t SIZE = SizeofAny<Arg...>();
	const int proc_id = GetProcessID(mpi_comm);
	//printf("[%d] size=%zd\n", proc_id, SIZE); fflush(stdout);
	uint8_t data[SIZE];
	
	if (proc_id == root_id) {
		CopyFromAny(data, arg...);
	}
	
	MPI_Bcast(data, SIZE, MPI_BYTE, root_id, mpi_comm);
	
	if (proc_id != root_id) {
		CopyToAny(data, arg...);
	}
	
}

#else
//with tuple, but data size becomes larger than the above template meta-programing code.

template<class ...Arg>
void BroadcastAny(const MPI_Comm& mpi_comm, int root_id, Arg& ...arg) {
	using mytuple = std::tuple<Arg...>;
	mytuple t{ arg... };
	constexpr size_t SIZE = sizeof(mytuple);
	const int proc_id = GetProcessID(mpi_comm);
	printf("[%d] size=%zd\n", proc_id, SIZE); fflush(stdout);
	//uint8_t data[SIZE];
	/*
	if (proc_id == root_id) {
		CopyFromAny(data, arg...);
	}
	*/

	MPI_Bcast(&t, SIZE, MPI_BYTE, root_id, mpi_comm);

	if (proc_id != root_id) {
		std::tie(arg...) = t;
	}

}

#endif

