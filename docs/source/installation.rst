.. _installation:

Getting Started
===============================


Singularity
---------------------

We recommend running DNAscent using one of our supported Singularity images. These images contain all necessary dependencies including TensorFlow, CUDA, CuDNN, and compression plugins so that your system only needs a valid NVIDIA driver for GPU usage. If your system does not have Singularity installed, instructions are available `here <https://docs.sylabs.io/guides/3.0/user-guide/installation.html>`_.

.. code-block:: console

   singularity pull DNAscent.sif library://mboemo/dnascent/dnascent:4.2.2
   
You can run DNAscent from the image by passing the desired executable and arguments. The following example shows how to run DNAscent :ref:`detect`:

.. code-block:: console

   singularity run --nv DNAscent.sif detect -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.bam


Building from Source
---------------------

Clone the DNAscent repository with the recursive flag so that the dependencies are cloned as well.

.. code-block:: console

   git clone --recursive https://github.com/MBoemo/DNAscent.git

The DNAscent directory will appear in your current directory. Switch to the latest release and compile the software by running:

.. code-block:: console

   cd DNAscent
   git checkout 4.2.2
   make

By default, ``make`` downloads ``libtensorflow-gpu-linux-x86_64-2.4.1.tar.gz`` because this is the runtime used for the current DNAscent release models. If you need a different TensorFlow C runtime (for example to match newer CUDA stacks on newer GPUs), you can override the download settings at build time:

.. code-block:: console

   make TF_VERSION=<version> TF_VARIANT=<gpu|cpu>

Some newer TensorFlow C runtime archives are published under a versioned URL layout, for example ``https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-gpu-linux-x86_64.tar.gz``. For this layout, use:

.. code-block:: console

   make TF_URL_LAYOUT=versions TF_VERSION=<version> TF_VARIANT=<gpu|cpu>

For full control over the archive name or URL:

.. code-block:: console

   make TF_TARBALL=<archive-name.tar.gz> TF_URL=<https://.../archive-name.tar.gz>

This will put the DNAscent executable into the DNAscent/bin directory. Compilation requires a version of gcc that supports C++14, and a typical compile time for DNAscent and all of its dependencies is 5-7 minutes.

Cloning the repository recursively (see above) will provide all the required dependencies so you don't need to find them yourself. For completeness, however, they are listed here:

* pfasta (https://github.com/kloetzl/pfasta)
* fast5 (https://github.com/mateidavid/fast5.git)
* pod5 (https://github.com/nanoporetech/pod5-file-format)
* htslib (https://github.com/samtools/htslib.git)
* hdf5lib (https://support.hdfgroup.org/HDF5/)
* tinydir (https://github.com/cxong/tinydir.git)
* TensorFlow (https://www.tensorflow.org/install/lang_c)

Please note that the high throughput sequencing library (htslib) requires bzlib and lzma for compression. While these are common on most systems, if you don't have these, apt-get lzma-dev, liblzma-dev, and libbz2-dev. In addition, pfasta requires libbsd.

FAST5 files are compressed with VBZ Compression (see https://github.com/nanoporetech/vbz_compression).  To use DNAscent on these compressed FAST5 files, do the following (N.B., we're assuming you don't have root permissions):

#. Go to https://github.com/nanoporetech/vbz_compression/releases and download the plugin appropriate for your processor architecture.  In this example, we'll use ``ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz``.

#. Download and unpack the plugin:

   .. code-block:: console

      wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
      tar -xf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz

#. Add the plugin to your path:

   .. code-block:: console

      export HDF5_PLUGIN_PATH=/full/path/to/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

#. Run ``DNAscent detect`` as normal.

The ``DNAscent detect`` executable can make use of a GPU, although this is optional (see :ref:`detect`). With the default TensorFlow runtime (2.4.1), DNAscent requires a CUDA 11.x and cuDNN 8.x compatible stack. If you build with a different TensorFlow runtime, the required CUDA and cuDNN versions are determined by that TensorFlow build.

* cuDNN: https://developer.nvidia.com/cudnn
* CUDA: https://developer.nvidia.com/cuda-11.0-download-archive

