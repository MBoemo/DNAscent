CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS =
CXXFLAGS = -Wall -O2 -fopenmp -std=c++14 $(DEBUG)
CFLAGS = -Wall -std=c99 -O2 $(DEBUG)

#hdf5
H5_LIB = ../../../hdf5-1.8.14/hdf5/lib/libhdf5.a
H5_INCLUDE = -I./../../../hdf5-1.8.14/hdf5/include
LIBFLAGS += -Wl,-rpath,../../../hdf5-1.8.14/hdf5/lib -L ../../../hdf5-1.8.14/hdf5/lib -lhdf5

#hts
HTS_LIB = ../../../htslib/libhts.a
HTS_INCLUDE = -I./../../../htslib
LIBFLAGS += -Wl,-rpath,../../../htslib -L ../../../htslib/ -lhts

#fast5
FAST5_INCLUDE = -I./../../../fast5/include

#add include flags for each library
CXXFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE)

MAIN_EXECUTABLE = test_htslib_interface
EXC_SRC = test_htslib_interface.cpp

#compile the test executable
$(MAIN_EXECUTABLE): $(HTS_LIB) $(H5_LIB) $(EXC_SRC)
	$(CXX) -o $@ $(CXXFLAGS) $(EXC_SRC) ../../../src/data_IO.o ../../../src/pfasta/pfasta.o $(LIBFLAGS)
	
clean:
	rm $(MAIN_EXECUTABLE)