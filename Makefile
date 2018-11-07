CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS =
CXXFLAGS = -Wall -O2 -fopenmp -std=c++11 $(DEBUG)
CFLAGS = -Wall -std=c99 -O2 $(DEBUG)

#hdf5
H5_LIB = ./hdf5-1.8.14/hdf5/lib/libhdf5.a
H5_INCLUDE = -I./hdf5-1.8.14/hdf5/include
LIBFLAGS += -Wl,-rpath,$(dir $(abspath $(lastword $(MAKEFILE_LIST))))hdf5-1.8.14/hdf5/lib -L hdf5-1.8.14/hdf5/lib -lhdf5

#hts
HTS_LIB = ./htslib/libhts.a
HTS_INCLUDE = -I./htslib
LIBFLAGS += -Wl,-rpath,$(dir $(abspath $(lastword $(MAKEFILE_LIST))))htslib -L htslib/ -lhts

#fast5
FAST5_INCLUDE = -I./fast5/include

#add include flags for each library
CXXFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE)

MAIN_EXECUTABLE = bin/DNAscent

all: depend $(MAIN_EXECUTABLE)

#all each library if they're not already built
htslib/libhts.a:
	cd htslib && make && cd .. || exit 255

hdf5-1.8.14/hdf5/lib/libhdf5.a:
	if [ ! -e hdf5-1.8.14/hdf5/lib/libhdf5.a ]; then \
		wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz; \
		tar -xzf hdf5-1.8.14.tar.gz || exit 255; \
		cd hdf5-1.8.14 && \
			./configure --enable-threadsafe && \
			make && make install; \
	fi 
	
SUBDIRS = src src/scrappie src/pfasta
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))
EXE_SRC = src/DNAscent.cpp

#generate object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

depend: .depend

.depend: $(CPP_SRC) $(C_SRC) $(EXE_SRC) $(H5_LIB)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

#compile each object
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(H5_INCLUDE) -fPIC $<

#compile the main executable
$(MAIN_EXECUTABLE): src/DNAscent.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB)
	$(CXX) -o $@ $(CXXFLAGS) -fPIC $(CPP_OBJ) $(C_OBJ) $(LIBFLAGS)

clean:
	rm -f $(MAIN_EXECUTABLE) $(CPP_OBJ) $(C_OBJ) src/DNAscent.o
