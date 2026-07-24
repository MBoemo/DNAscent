CC = gcc
CXX = g++
DEBUG = -g
LIBFLAGS = -lrt 
LDFLAGS ?= -ldl -llzma -lbz2 -lm -lz
CXXFLAGS = -Wall -O2 -fopenmp -std=c++17
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

#libtorch (PyTorch C++ / LibTorch)
#The detection model is a TorchScript module loaded via torch::jit::load.
#
#By default we link against the LibTorch that ships inside the installed Python
#'torch' package: this is the exact build that scripts model.pt (via
#utils/export_torchscript.py), so ABI/opset compatibility is guaranteed and no
#multi-GB download is needed.  Override TORCH_HOME to point at a standalone
#libtorch/ directory (e.g. an official release) if you prefer.
TORCH_HOME ?= $(shell python3 -c 'import torch, os; print(os.path.dirname(torch.__file__))' 2>/dev/null)
ifeq ($(strip $(TORCH_HOME)),)
TORCH_HOME := $(CURRENT_PATH)libtorch
endif
TORCH_DEPEND = $(TORCH_HOME)/include/torch/script.h
TORCH_INCLUDE = -I$(TORCH_HOME)/include -I$(TORCH_HOME)/include/torch/csrc/api/include
TORCH_LIB = -Wl,-rpath,$(TORCH_HOME)/lib -L $(TORCH_HOME)/lib -ltorch -ltorch_cpu -lc10
#Match the ABI the LibTorch libraries were built with (CXX11 ABI for modern
#PyTorch 2.x Linux builds).
TORCH_DEFINES = -D_GLIBCXX_USE_CXX11_ABI=1
#link the CUDA runtime libraries only if this LibTorch build provides them, and
#tell the source that the CUDA stream API is available so it can overlap
#per-thread GPU inference (see runCNN in src/detect.cpp).
ifneq ($(wildcard $(TORCH_HOME)/lib/libtorch_cuda.so),)
TORCH_LIB += -ltorch_cuda -lc10_cuda
TORCH_DEFINES += -DDNASCENT_CUDA=1
#the CUDA runtime (libcudart) is needed for the stream API; conda/pip PyTorch
#ships it either in torch/lib or in the sibling nvidia/cuda_runtime wheel.
CUDART_DIR := $(dir $(firstword $(wildcard $(TORCH_HOME)/lib/libcudart.so*) $(wildcard $(TORCH_HOME)/../nvidia/cuda_runtime/lib/libcudart.so*)))
ifneq ($(strip $(CUDART_DIR)),)
TORCH_LIB += -Wl,-rpath,$(CUDART_DIR) -L $(CUDART_DIR) -l:libcudart.so.12
endif
#Build the fused selective-scan CUDA op (src/ssm_op.cu) and link it into the
#binary so the scripted model's dnascent::selective_scan op is registered.
HAVE_CUDA := 1
SSM_OBJ := src/ssm_op.o
#nvcc: the conda toolkit can be incomplete (missing nvvm/cicc); prefer a full
#system CUDA toolkit if one is present.  Override any of these on the make line.
NVCC ?= $(firstword $(wildcard /usr/local/cuda/bin/nvcc) $(shell command -v nvcc 2>/dev/null))
#nvcc rejects host gcc > 12; use g++-12 when available.
NVCC_CCBIN ?= $(firstword $(wildcard /usr/bin/g++-12) g++)
#SM architecture of the target GPU (89 = Ada; override for other cards).
CUDA_ARCH ?= 89
endif
#Fallback download (used only if TORCH_HOME resolves to ./libtorch and it is
#missing).  Set TORCH_VARIANT to match your CUDA (cpu, cu121, cu124, cu128 ...).
TORCH_VERSION ?= 2.9.1
TORCH_VARIANT ?= cu128
TORCH_ZIP ?= libtorch-cxx11-abi-shared-with-deps-$(TORCH_VERSION)%2B$(TORCH_VARIANT).zip
TORCH_URL ?= https://download.pytorch.org/libtorch/$(TORCH_VARIANT)/$(TORCH_ZIP)

#fast5
FAST5_INCLUDE = -I./fast5/include

#pod5
POD5_DEPEND = pod5-file-format/build/Release/lib/libpod5_format.a
POD5_INCLUDE = -I./pod5-file-format/build/c++
POD5_LIB = -L${PATH_SPACEFIX}pod5-file-format/build/Release/lib -lpod5_format 
POD5_LIB += -L${PATH_SPACEFIX}pod5-file-format/build/third_party/libs -larrow -ljemalloc_pic -lzstd

#add include flags for each library
CPPFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE) $(TORCH_INCLUDE) $(POD5_INCLUDE) $(TORCH_DEFINES)

DNASCENT_EXECUTABLE = bin/DNAscent

all: depend $(DNASCENT_EXECUTABLE)

#all each library if they're not already built
htslib/libhts.a:
	cd htslib; \
	make || exit 255; \
	cd ..;

hdf5-1.8.14/hdf5/lib/libhdf5.a:
	if [ ! -e hdf5-1.8.14/hdf5/lib/libhdf5.a ]; then \
		wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz; \
		tar -xzf hdf5-1.8.14.tar.gz || exit 255; \
		cd hdf5-1.8.14 && \
			./configure --enable-threadsafe && \
			make && make install; \
	fi 

$(TORCH_HOME)/include/torch/script.h:
	if [ ! -e "$(TORCH_HOME)/include/torch/script.h" ]; then \
		wget -O libtorch.zip "$(TORCH_URL)"; \
		unzip -q libtorch.zip || exit 255; \
		rm -f libtorch.zip; \
	fi
	
pod5-file-format/build/Release/lib/libpod5_format.a:
	if [ ! -e pod5-file-format/build/Release/lib/libpod5_format.a ]; then \
		pip3 install "conan<2" build; \
		conan --version; \
		cd pod5-file-format; \
		git submodule update --init --recursive; \
		pip3 install setuptools_scm==7.1.0; \
		python3 -m setuptools_scm; \
		python3 -m pod5_make_version; \
		mkdir build; \
		cd build; \
		conan install --build=missing -s build_type=Release .. && cmake -DENABLE_CONAN=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake .. && make -j; \
		cd ../..; \
	fi
	
SUBDIRS = src src/scrappie src/pfasta src/sgsmooth
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))
DNA_EXE_SRC = src/main/DNAscent.cpp

#log the commit 
src/gitcommit.h: .git/HEAD .git/index
	echo "const char *gitcommit = \"$(shell git rev-parse HEAD)\";" > $@

#log the software path
src/softwarepath.h: 
	echo "const char *executablePath = \"${PATH_SPACEFIX}\";" > $@

#generate object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

DNASCENT_OBJ = $(DNA_EXE_SRC:..cpp=.0)

depend: .depend

.depend: $(CPP_SRC) $(C_SRC) $(H5_LIB) $(TORCH_DEPEND) $(POD5_DEPEND) src/gitcommit.h src/softwarepath.h
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

#compile the fused selective-scan CUDA op with nvcc (only when CUDA is present)
src/ssm_op.o: src/ssm_op.cu
	$(NVCC) -ccbin $(NVCC_CCBIN) -O3 -std=c++17 --compiler-options -fPIC \
		-gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH) \
		$(TORCH_DEFINES) $(TORCH_INCLUDE) -c $< -o $@

#standalone shared library of the op for the Python export / validation scripts
#(torch.ops.load_library).  The C++ binary links src/ssm_op.o directly instead.
libdnascent_ssm.so: src/ssm_op.cu
	$(NVCC) -ccbin $(NVCC_CCBIN) -O3 -std=c++17 --compiler-options -fPIC -shared \
		-gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH) \
		$(TORCH_DEFINES) $(TORCH_INCLUDE) $< -o $@ \
		-Xlinker -rpath -Xlinker $(TORCH_HOME)/lib -L $(TORCH_HOME)/lib \
		-ltorch -ltorch_cpu -lc10 -ltorch_cuda -lc10_cuda

.PHONY: ssm_so
ssm_so: libdnascent_ssm.so

#compile each object
.cpp.o: src/gitcommit.h src/softwarepath.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS) $(H5_INCLUDE) -fPIC $<
	
src/main/DNAscent.o: src/gitcommit.h src/softwarepath.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

#compile the main executables
$(DNASCENT_EXECUTABLE): src/main/DNAscent.o $(CPP_OBJ) $(C_OBJ) $(SSM_OBJ) $(HTS_LIB) $(H5_LIB) $(TORCH_DEPEND) $(POD5_DEPEND) src/gitcommit.h src/softwarepath.h
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $(DNASCENT_OBJ) $(CPP_OBJ) $(C_OBJ) $(SSM_OBJ) $(HTS_LIB) $(H5_LIB) $(TORCH_LIB) $(POD5_LIB) $(LIBFLAGS) $(LDFLAGS)

clean:
	rm -f $(DNASCENT_EXECUTABLE) $(CPP_OBJ) $(C_OBJ) $(SSM_OBJ) libdnascent_ssm.so src/main/DNAscent.o src/gitcommit.h
