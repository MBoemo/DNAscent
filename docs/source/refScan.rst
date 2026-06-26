.. _refScan:

refScan
===============================

``DNAscent refScan`` is a ``DNAscent`` subprogram that scores each read's pore current signal against a genome-wide reference signal profile produced by ``DNAscent refSignal``.  For every position in a read's event alignment, a rolling k-mer window is evaluated: the observed signal at each position in the window is compared to the Gaussian distribution for that position stored in the reference profile, and a mean per-position log-likelihood is reported.  Positions where this score is more negative than the expected baseline are consistent with a base modification that shifts the pore current away from its unmodified value.

Usage
-----

.. code-block:: console

   To run DNAscent refScan, do:
      DNAscent refScan -b /path/to/alignment.bam -r /path/to/reference.fasta \
                       -i /path/to/index.dnascent -s /path/to/refSignal \
                       -o /path/to/output
   Required arguments are:
     -b,--bam                  path to alignment BAM file,
     -r,--reference            path to genome reference in fasta format,
     -i,--index                path to DNAscent index,
     -s,--refSignal            path to reference signal file produced by DNAscent refSignal,
     -o,--output               path to output file.
   Optional arguments are:
     -t,--threads              number of threads (default is 1 thread),
     -q,--quality              minimum mapping quality (default is 20),
     -l,--length               minimum read length in bp (default is 1000),
     -w,--minWindowCov         minimum number of positions in the k-mer window that
                               must have both read signal and reference signal for a
                               window to be scored (default is 5).
                               Maximum is the k-mer length (9 for R10.4.1, 6 for R9).

The main inputs of ``DNAscent refScan`` are the alignment BAM file for the experimental reads and the reference signal profile produced by ``DNAscent refSignal`` (passed with the ``-s`` flag).  The reference genome passed with ``-r`` and the index passed with ``-i`` must be the same files used to produce the input BAM.

The number of threads is specified with ``-t``.  ``DNAscent refScan`` parallelises by processing a batch of reads simultaneously on separate threads, so multithreading is recommended.

The ``-w`` flag controls how many of the k-mer window positions must have both an observed signal and a valid reference Gaussian (mean and non-zero standard deviation) before a window score is emitted.  Lowering this threshold increases the number of positions scored at the cost of less reliable scores near read ends and deletion-rich regions.

Method
------

For each read that passes quality filters, ``DNAscent refScan`` fetches the raw pore current from the FAST5 or POD5 file, normalises the signal, and performs the same event alignment used by ``DNAscent align``.  This maps every segment of the signal to the reference position it most likely corresponds to.

Once the event alignment is complete, ``DNAscent refScan`` slides a k-mer window across every position in the read for which signal was aligned.  The window half-widths are derived from the pore model k-mer length (``Pore_Substrate_Config.kmer_len``), so the window automatically spans 9 positions for R10.4.1 data and 6 positions for R9.4.1 data:

* **R10.4.1** — window of 9 positions: [center − 4, center + 4]
* **R9.4.1** — window of 6 positions: [center − 3, center + 2]

For each position in the window that has both an observed signal mean (from the current read) and a reference Gaussian (from the ``refSignal`` file), the per-position log-likelihood is computed as:

.. math::

   \ell_i = \log \mathcal{N}\!\left(\hat{x}_i;\, \mu_i^{\text{ref}},\, \sigma_i^{\text{ref}}\right)
          = -\frac{1}{2}\!\left[\log(2\pi) + 2\log\!\left(\sigma_i^{\text{ref}}\right) + \left(\frac{\hat{x}_i - \mu_i^{\text{ref}}}{\sigma_i^{\text{ref}}}\right)^{\!2}\right]

where :math:`\hat{x}_i` is the mean of the raw signal values at position :math:`i` in the current read, and :math:`\mu_i^{\text{ref}}` and :math:`\sigma_i^{\text{ref}}` are the reference mean and standard deviation from ``DNAscent refSignal``.

The window score reported for the centre position is the mean of :math:`\ell_i` over the :math:`n_{\text{valid}}` positions that met the validity criteria:

.. math::

   \text{score}(c) = \frac{1}{n_{\text{valid}}} \sum_{i} \ell_i

Under the null hypothesis that no modification is present, the observed signal at each position is drawn from the reference Gaussian and the expected score is approximately :math:`-\tfrac{1}{2}[\log(2\pi) + \langle 2\log\sigma^{\text{ref}}\rangle + 1]`.  A modification that shifts the pore current away from its unmodified value will produce :math:`z`-scores greater than one, driving the window score to more negative values.  Positions where the score is substantially more negative than the baseline are therefore candidate modification sites.

Output
------

``DNAscent refScan`` writes a single output file to the path specified with ``-o``.  Like other ``DNAscent`` executables, the file begins with a short header in which every line starts with a hash (``#``) character:

.. code-block:: console

   #Alignment /path/to/alignment.bam
   #Genome /path/to/reference.fasta
   #Index /path/to/index.dnascent
   #RefSignal /path/to/refSignal
   #Threads 4
   #MappingQuality 20
   #MappingLength 1000
   #MinWindowCov 5
   #Version 4.2.2
   #Commit 4cf80a7b89bdf510a91b54572f8f94d3daf9b167

You can access the header with ``grep '#' /path/to/output``.

Below the header, the file is organised into per-read blocks.  Each block begins with a read header line that starts with a greater-than (``>``) character:

.. code-block:: console

   >readID contig mappingStart mappingEnd strand

The fields are the same as those used by ``DNAscent detect`` and ``DNAscent align``:

* ``readID`` is the unique hexadecimal identifier assigned to each read by the Oxford Nanopore software,
* the read mapped between ``mappingStart`` and ``mappingEnd`` on ``contig``,
* ``strand`` is ``fwd`` if the read mapped to the forward strand, or ``rev`` if it mapped to the reverse complement.

You can count the number of reads in the output with ``grep '>' /path/to/output | wc -l``.

Below each read header, each line corresponds to one window-centre position and contains three tab-separated columns:

* the 0-based coordinate on the reference genome,
* the mean per-position log-likelihood for the k-mer window centred at that coordinate (lower values indicate greater deviation from the reference signal),
* the number of positions in the window that contributed to the score (between ``-w`` and the full k-mer length).

An example excerpt is:

.. code-block:: console

   >a4ea2872-9cb6-4218-afad-905f79204eb1 1 992440 1002318 fwd
   992444  -1.843271  9
   992445  -1.912034  9
   992446  -2.107658  9
   992447  -3.841920  9
   992448  -5.293417  9
   992449  -5.108834  9
   992450  -2.014503  9
   992451  -1.799612  9

In this example, positions 992447–992449 show notably more negative scores, suggesting a signal deviation at or near those positions that may be consistent with a base modification.

A log file is also written alongside the output.  Its name is derived from the output filename with a ``.refScan.log`` suffix.  Any reads that could not be located in the DNAscent index are reported there.
