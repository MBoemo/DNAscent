.. _installation:

Download & Installation
===============================

Clone the DNAscent repository with the recursive flag so that the dependencies are cloned as well.

.. code-block:: console

   git clone --recursive https://github.com/MBoemo/DNAscent.git

The DNAscent directory will appear in your current directory. Switch to the latest release and compile the software by running:

.. code-block:: console

   cd DNAscent
   git checkout 4.0.1
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

VBZ Fast5 Compression
---------------------

In new versions of MinKNOW, the fast5 files are compressed with VBZ Compression (see https://github.com/nanoporetech/vbz_compression).  To use DNAscent on these compressed fast5 files, do the following (N.B., we're assuming you don't have root permissions):

#. Go to https://github.com/nanoporetech/vbz_compression/releases and download the plugin appropriate for your processor architecture.  In this example, we'll use ``ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz``.

#. Download and unpack the plugin:

   .. code-block:: console

      wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
      tar -xf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz

#. Add the plugin to your path:

   .. code-block:: console

      export HDF5_PLUGIN_PATH=/full/path/to/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

#. Run ``DNAscent detect`` as normal.

GPU Use
-------

The ``DNAscent detect`` executable can make use of a GPU, although this is optional (see :ref:`detect`).  DNAscent requires CUDA 11.8 and cuDNN 8.9. Information about these can be found at the following links:

* cuDNN: https://developer.nvidia.com/cudnn
* CUDA: https://developer.nvidia.com/cuda-11.0-download-archive

Always discuss any installation or version changes with your system administrator.
