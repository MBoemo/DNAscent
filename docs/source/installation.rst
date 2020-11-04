.. _installation:

Download & Installation
===============================

Development was done using gcc 9.3.0 on an Ubuntu 16.04 platform. While installation on other platforms is possible, Ubuntu is the platform that is recommended and supported.

Clone the DNAscent repository with the recursive flag so that the dependencies are cloned as well.

.. code-block:: console

   git clone --recursive https://github.com/MBoemo/DNAscent.git

The DNAscent directory will appear in your current directory. Switch to the latest release and compile the software by running:

.. code-block:: console

   cd DNAscent
   git checkout 2.0.0
   make

This will put the DNAscent executable into the DNAscent/bin directory. Compilation requires a version of gcc that supports C++14, and a typical compile time for DNAscent and all of its dependencies is 5-7 minutes.

Cloning the repository recursively (see above) will provide all the required dependencies so you don't need to find them yourself. For completeness, however, they are listed here:

* pfasta (https://github.com/kloetzl/pfasta)
* fast5 (https://github.com/mateidavid/fast5.git)
* htslib (https://github.com/samtools/htslib.git)
* hdf5lib (https://support.hdfgroup.org/HDF5/)
* tinydir (https://github.com/cxong/tinydir.git)
* TensorFlow (https://www.tensorflow.org/install/lang_c)

Please note that the high throughput sequencing library (htslib) requires bzlib and lzma for compression. While these are common on most systems, if you don't have these, apt-get lzma-dev, liblzma-dev, and libbz2-dev. In addition, pfasta requires libbsd on Linux.

GPU Use
-------

The ``DNAscent detect`` executable can make use of a GPU, although this is optional (see :ref:`detect`).  DNAscent requires CUDA 10.0 and cuDNN 7.5, and information about these can be found at the following links:

* cuDNN: https://developer.nvidia.com/cudnn
* CUDA: https://developer.nvidia.com/cuda-10.0-download-archive

Always discuss any installation or version changes with your system administrator.
