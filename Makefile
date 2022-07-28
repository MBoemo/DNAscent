CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS = -lrt 
LDFLAGS ?= -ldl -llzma -lbz2 -lm -lz
CXXFLAGS = -Wall -O2 -fopenmp -std=c++14
CFLAGS = -Wall -std=c99 -O2

SPACE:= ;
SPACE+=;
null :=
space := ${null} ${null}
${space} := ${space}

CURRENT_PATH := $(subst $(lastword $(notdir $(MAKEFILE_LIST))),,$(subst $(SPACE),\$(SPACE),$(shell realpath '$(strip $(MAKEFILE_LIST))')))
PATH_SPACEFIX := $(subst ${space},\${space},${CURRENT_PATH})

ifeq ($(zstd),1)
	LDFLAGS += -lzstd
endif

#hdf5
H5_LIB = ./hdf5-1.8.14/hdf5/lib/libhdf5.a
H5_INCLUDE = -I./hdf5-1.8.14/hdf5/include

#hts
HTS_LIB = ./htslib/libhts.a
HTS_INCLUDE = -I./htslib

#tensorflow
TENS_DEPEND = tensorflow/include/tensorflow/c/c_api.h
TENS_LIB = -Wl,-rpath,${PATH_SPACEFIX}tensorflow/lib -L tensorflow/lib
TENS_INCLUDE = -I./tensorflow/include
LIBFLAGS = -ltensorflow

#fast5
FAST5_INCLUDE = -I./fast5/include

#add include flags for each library
CPPFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE) $(TENS_INCLUDE)

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

tensorflow/include/tensorflow/c/c_api.h:
	if [ ! -e tensorflow/include/tensorflow/c/c_api.h ]; then \
		mkdir tensorflow; \
		cd tensorflow; \
		wget https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-gpu-linux-x86_64-2.4.1.tar.gz; \
		tar -xzf libtensorflow-gpu-linux-x86_64-2.4.1.tar.gz || exit 255; \
		cd ..; \
	fi 
	
SUBDIRS = src src/scrappie src/pfasta src/sgsmooth
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))
EXE_SRC = src/DNAscent.cpp

#log the commit 
src/gitcommit.h: .git/HEAD .git/index
	echo "const char *gitcommit = \"$(shell git rev-parse HEAD)\";" > $@

#log the software path
src/softwarepath.h: 
	echo "const char *executablePath = \"${PATH_SPACEFIX}\";" > $@

#generate object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

depend: .depend

.depend: $(CPP_SRC) $(C_SRC) $(EXE_SRC) $(H5_LIB) $(TENS_DEPEND) src/gitcommit.h src/softwarepath.h
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

#compile each object
.cpp.o: src/gitcommit.h src/softwarepath.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS) $(H5_INCLUDE) -fPIC $<
	

#compile the main executable
$(MAIN_EXECUTABLE): src/DNAscent.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB) $(TENS_DEPEND) src/gitcommit.h src/softwarepath.h
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB) $(TENS_LIB) $(LIBFLAGS) $(LDFLAGS)

clean:
	rm -f $(MAIN_EXECUTABLE) $(CPP_OBJ) $(C_OBJ) src/DNAscent.o src/gitcommit.h src/softwarepath.h
