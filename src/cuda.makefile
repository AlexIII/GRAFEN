CUDA_PATH?= /usr/local/cuda
NVCC = $(CUDA_PATH)/bin/nvcc

CXX = g++
CXXFLAGS = -O2 --std=c++17
GNU_CXXFLAGS = -fopenmp -fPIE
NVCC_FLAGS = --extended-lambda --expt-relaxed-constexpr -Xcompiler -fopenmp -Wno-deprecated-gpu-targets -arch compute_35

LDFLAGS = -lgomp

TARGET = grafen_cuda
OBJS_CPU = elFieldCU.o mobj.o Grid/Grid.o
OBJS_CUDA = calcField.o

#------------- Specify dependency paths ------------------

#GeographicLib
GEOLIB_PATH = ../deps/GeographicLib_installed
GEOLIB_INCLUDE_PATH = $(GEOLIB_PATH)/include/
GEOLIB_LIB_PATH = $(GEOLIB_PATH)/lib/

#boost
BOOST_PATH = ../deps/boost_1_73_0
BOOST_INCLUDE_PATH = $(BOOST_PATH)

#MPICH
MPI_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/mpi
MPI_LIB_PATH = /usr/lib/x86_64-linux-gnu/

#---------------------------------------------------------

all: $(TARGET)

# compile CUDA
$(OBJS_CUDA): %.o: %.cu
	$(NVCC) -o $@ $< -c $(CXXFLAGS) $(NVCC_FLAGS)

#comple CPU
$(OBJS_CPU): %.o: %.cpp
	$(CXX) -o $@ $< -c $(CXXFLAGS) $(GNU_CXXFLAGS) -I$(GEOLIB_INCLUDE_PATH) -I$(BOOST_INCLUDE_PATH) -I$(MPI_INCLUDE_PATH)

$(TARGET): $(OBJS_CPU) $(OBJS_CUDA)
	$(NVCC) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -L$(GEOLIB_LIB_PATH) -lGeographic -lrt -L$(MPI_LIB_PATH) -lmpich

clean:
	rm -rf *.o $(TARGET)
