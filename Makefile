CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS =
CXXFLAGS = -Wall -O2 -fopenmp -std=c++11 $(DEBUG)
CFLAGS = -Wall -std=c99 -O2 $(DEBUG)

#Penthus
PENTHUS_LIB = ./Penthus/libPenthus.a
PENTHUS_INCLUDE = -I./Penthus
LIBFLAGS += -L Penthus/ -lPenthus

#hdf5
H5_LIB = ./hdf5-1.8.14/hdf5/lib/libhdf5.a
H5_INCLUDE = -I./hdf5-1.8.14/hdf5/include
LIBFLAGS += -Wl,-rpath,hdf5-1.8.14/hdf5/lib -L hdf5-1.8.14/hdf5/lib -lhdf5

#hts
HTS_LIB = ./htslib/libhts.a
HTS_INCLUDE = -I./htslib
LIBFLAGS += -L htslib/ -lhts

#fast5
FAST5_INCLUDE = -I./fast5/include

#add include flags for each library dependency
CXXFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE) $(PENTHUS_INCLUDE)

MAIN_EXECUTABLE = bin/Osiris

all: depend $(MAIN_EXECUTABLE)

htslib/libhts.a:
	cd htslib && make || exit 255

Penthus/lPenthus.a:
	cd Penthus && make || exit 255

#
# If this library is a dependency the user wants HDF5 to be downloaded and built.
#
lib/libhdf5.a:
#	if [ ! -e hdf5-1.8.14.tar.gz ]; then \
#		wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz; \
#	fi
#	tar -xzf hdf5-1.8.14.tar.gz || exit 255
#	cd hdf5-1.8.14 && \
#		./configure --enable-threadsafe --prefix=`pwd`/.. && \
#		make && make install

SUBDIRS = src src/scrappie
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))
EXE_SRC = src/Osiris.cpp

# Automatically generated object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

depend: .depend

.depend: $(CPP_SRC) $(C_SRC) $(EXE_SRC) $(H5_LIB)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

# Compile objects
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(H5_INCLUDE) -fPIC $<

# Link main executable
$(MAIN_EXECUTABLE): src/Osiris.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB) $(PENTHUS_LIB)
	$(CXX) -o $@ $(CXXFLAGS) -fPIC $(CPP_OBJ) $(C_OBJ) $(LIBFLAGS)

#$(HTS_LIB) $(H5_LIB) $(PENTHUS_LIB)
clean:
	rm -f $(MAIN_EXECUTABLE) $(CPP_OBJ) $(C_OBJ) src/Osiris.o
