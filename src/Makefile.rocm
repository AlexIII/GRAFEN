HIP_PATH?= $(wildcard /opt/rocm/hip)
HIPCC = $(HIP_PATH)/bin/hipcc
HIPCC_INCLUDE_PATH = /opt/rocm/include/
HIPCC_FLAGS = -Wno-unused-result

CXX = g++
CXXFLAGS = -O2 --std=c++17 -fPIE
GNU_CXXFLAGS = -fopenmp

LDFLAGS = -fopenmp

TARGET = grafen_rocm
OBJS_CPU = elFieldCU.o mobj.o Grid/Grid.o
OBJS_HIP = calcField.hip.o

#------------- Specify dependency paths ------------------

#GeographicLib
GEOLIB_PATH = ../deps/GeographicLib_installed
GEOLIB_INCLUDE_PATH = $(GEOLIB_PATH)/include/
GEOLIB_LIB_PATH = $(GEOLIB_PATH)/lib/

#boost
BOOST_PATH = ../deps/boost_1_73_0
BOOST_INCLUDE_PATH = $(BOOST_PATH)

#MPICH
MPI_INCLUDE_PATH = /usr/lib/x86_64-linux-gnu/mpich/include/
MPI_LIB_PATH = /usr/lib/x86_64-linux-gnu/

#---------------------------------------------------------


all: $(TARGET)

# convert CUDA to HIP
calcField.hip.cpp: calcField.cu cuVar.h
	/opt/rocm/bin/hipify-perl calcField.cu > calcField.hip.cpp 
	sed -i 's/#include \"cuVar\.h\"/#include \"cuVar\.hip\.h\"/g' calcField.hip.cpp
	/opt/rocm/bin/hipify-perl cuVar.h > cuVar.hip.h 

# compile HIP
$(OBJS_HIP): %.o: %.cpp
	$(HIPCC) -o $@ $< -c $(CXXFLAGS) $(HIPCC_FLAGS) -I$(HIPCC_INCLUDE_PATH)

#comple CPU
$(OBJS_CPU): %.o: %.cpp
	$(CXX) -o $@ $< -c $(CXXFLAGS) $(GNU_CXXFLAGS) -I$(GEOLIB_INCLUDE_PATH) -I$(BOOST_INCLUDE_PATH) -I$(MPI_INCLUDE_PATH)

$(TARGET): $(OBJS_CPU) $(OBJS_HIP)
	$(HIPCC) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -L$(GEOLIB_LIB_PATH) -lGeographic -lrt -L$(MPI_LIB_PATH) -lmpich

clean:
	rm -rf *.o $(TARGET) *.hip.*
