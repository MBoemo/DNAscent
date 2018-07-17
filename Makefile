OBJS = event_detection.o event_handling.o common.o data_IO.o Osiris_detect.o Osiris.o
CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS = -pg -L Penthus/ -l Penthus -L htslib/ -l hts -L hdf5-1.8.14/hdf5/lib/ -l hdf5 -fopenmp
CXXFLAGS = -Wall -c -O2 -std=c++11 $(DEBUG)
CFLAGS = -Wall -c -std=c99 -O2 $(DEBUG)
LFLAGS = -Wall -O2 $(DEBUG)

MAIN_EXECUTABLE = bin/Osiris

#libraries
HTS_LIB=./htslib/libhts.a
HTS_INCLUDE=-I./htslib

HDF5_LIB=-./hdf5-1.8.14/hdf5/lib/libhdf5.a
HDF5_INCLUDE=-I./hdf5-1.8.14/hdf5/include

FAST5_INCLUDE=-I./fast5/include

PENTHUS_LIB=./Penthus/libPenthus.a


all: $(MAIN_EXECUTABLE)

$(HTS_LIB):
	cd htslib && make || exit 255
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./htslib && export LD_LIBRARY_PATH

$(HDF5_LIB):
	if [ ! -e hdf5-1.8.14.tar.gz ]; then wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz; fi
	tar -xzf hdf5-1.8.14.tar.gz || exit 255
	cd hdf5-1.8.14 && ./configure --enable-threadsafe && make && make install
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./hdf5-1.8.14/hdf5/lib && export LD_LIBRARY_PATH

$(PENTHUS_LIB):
	cd Penthus && make || exit 255

event_detection.o : src/scrappie/event_detection.h src/scrappie/event_detection.c
	$(CC) $(CFLAGS) src/scrappie/event_detection.c

event_handling.o : src/scrappie/event_detection.h src/common.h src/poreModels.h src/data_IO.h src/event_handling.h src/event_handling.cpp
	$(CXX) $(CXXFLAGS) src/event_handling.cpp $(LIBFLAGS)

common.o : src/error_handling.h src/common.h src/common.cpp
	$(CXX) $(CXXFLAGS) src/common.cpp $(LIBFLAGS)

data_IO.o : src/error_handling.h src/data_IO.h src/data_IO.cpp
	$(CXX) $(CXXFLAGS) src/data_IO.cpp $(LIBFLAGS)

Osiris_detect.o : src/poreSpecificParameters.h src/common.h src/build_model.h src/data_IO.h src/error_handling.h src/event_handling.h src/poreModels.h src/Osiris_detect.h src/Osiris_detect.cpp $(HTS_LIB) $(HDF5_LIB) $(PENTHUS_LIB)
	$(CXX) $(CXXFLAGS) $(HTS_INCLUDE) $(HDF5_INCLUDE) $(FAST5_INCLUDE) src/Osiris_detect.cpp $(LIBFLAGS)

Osiris.o : src/Osiris.cpp src/Osiris_detect.h src/data_IO.h src/build_model.h src/event_handling.h
	$(CXX) $(CXXFLAGS) src/Osiris.cpp $(LIBFLAGS)

$(MAIN_EXECUTABLE) : $(OBJS) $(HTS_LIB) $(HDF5_LIB) $(PENTHUS_LIB)
	$(CXX) $(LFLAGS) $(HTS_INCLUDE) $(HDF5_INCLUDE) $(FAST5_INCLUDE) $(OBJS) -o bin/Osiris $(LIBFLAGS)

clean:
	rm $(OBJS) $(MAIN_EXECUTABLE)
