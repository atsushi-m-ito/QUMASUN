
# for mdgw3 icpx (with MPI)

#CXX = mpiicpx -DUSE_MPI -D_USE_SCALAPACK -O3 -ipo -std=c++17 -qmkl -w2 -wd3180 -wd161 -qopt-report=max -qopt-report-file=opt-report.txt 
CXX = mpiicpx -DUSE_MPI -DUSE_SCALAPACK -std=c++17 -restrict -O3 -ip -qmkl -xCORE-AVX512 -qopt-zmm-usage=high -w2 -wd3180 -wd161 -qopt-report=max -qopt-report-file=opt-report.txt 
#CXX = icpx -D_USE_MPI -std=c++17 -restrict -O3 -ip -qmkl -xCORE-AVX512 -qopt-zmm-usage=high -w2 -wd3180 -wd161 -qopt-report=max -qopt-report-file=opt-report.txt 
#CXX = icpc -D_USE_MPI -std=c++17 -restrict -O3 -ip -qmkl -no-prec-div -xCORE-AVX512 -qopt-report=5 -qopt-zmm-usage=high -w2 -wd3180 -wd161 -static-libstdc++ -qopt-report-file=opt-report.txt

DEST          = ../
CXXFLAGS = -I$(MKLROOT)/include/fftw
LDFLAGS       = 
LDFLAGSMPI    = 
LIBS          = -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lifcore

PROGRAM = qumasun.exe
SRC = main_mpi.cpp $(wildcard qumasun_*.cpp)
#PROGRAM = test.exe
#SRC = test_mpi_difference.cpp
#SRC = test_wave_bundle.cpp
#SRCMPI = main_mpi.cpp
OBJS = $(SRC:.cpp=.o)
HEADERS=$(wildcard *.h)

export CXX



all:	$(PROGRAM)

$(PROGRAM):	$(OBJS) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)



.SUFFIXES:	.cpp .o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -f *.o $(PROGRAM)

install:	$(PROGRAM)
		cp -p $(PROGRAM) $(DEST)


