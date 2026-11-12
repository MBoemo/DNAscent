# DNAscent container image.
#
# This mirrors DNAscent.def and produces an OCI/Docker image that can be run
# directly with Docker or pulled by Apptainer/Singularity, e.g.:
#   singularity pull DNAscent.sif docker://ghcr.io/mboemo/dnascent:4.2.1
#
# DNAscent is software for detecting regions of BrdU and EdU incorporation in
# Oxford Nanopore reads.
#   Source: https://github.com/MBoemo/DNAscent
#   Documentation: https://dnascent.readthedocs.io
#   Web: https://www.boemogroup.org/
FROM nvidia/cuda:11.0.3-cudnn8-devel-ubuntu20.04

# The DNAscent git tag/branch to build (e.g. 4.2.1). Passed in by CI on release.
ARG DNASCENT_VERSION=master

LABEL org.opencontainers.image.title="DNAscent" \
      org.opencontainers.image.description="Software for detecting regions of BrdU and EdU incorporation in Oxford Nanopore reads." \
      org.opencontainers.image.source="https://github.com/MBoemo/DNAscent" \
      org.opencontainers.image.documentation="https://dnascent.readthedocs.io" \
      org.opencontainers.image.url="https://www.boemogroup.org/" \
      org.opencontainers.image.version="${DNASCENT_VERSION}"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Europe/London

# Install system packages
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone \
    && apt-get update \
    && apt-get install -y tzdata \
    && apt-get install -y \
        lzma-dev \
        build-essential \
        liblzma-dev \
        libbz2-dev \
        libbsd-dev \
        wget \
        git \
        zlib1g-dev \
        libncurses5-dev \
        libgdbm-dev \
        libnss3-dev \
        libssl-dev \
        libreadline-dev \
        libffi-dev \
        libsqlite3-dev \
        python3 \
        python3-pip \
        libboost-all-dev \
        cmake \
    && rm -rf /var/lib/apt/lists/*

# Install conan
RUN pip3 install "conan<2" build

# Clone and compile DNAscent
RUN mkdir -p /app \
    && cd /app \
    && git clone --recursive --branch "${DNASCENT_VERSION}" https://github.com/MBoemo/DNAscent.git \
    && cd DNAscent \
    && make

# Install vbz plugin
RUN mkdir -p /plugin \
    && cd /plugin \
    && wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz \
    && tar -xf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz \
    && rm ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz

# cuda paths and vbz plugin path
ENV CUDA_HOME=/usr/local/cuda \
    CUDA_PATH=/usr/local/cuda \
    CPATH=/usr/local/cuda/include \
    LIBRARY_PATH=/usr/local/cuda/lib64 \
    LD_LIBRARY_PATH=/usr/local/cuda/lib64 \
    PATH=/usr/local/cuda/bin:/app/DNAscent/bin:$PATH \
    HDF5_PLUGIN_PATH=/plugin/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

ENTRYPOINT ["/app/DNAscent/bin/DNAscent"]
