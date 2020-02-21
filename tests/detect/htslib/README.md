# DNAscent Test - htslib Interface

This test checks DNAscent's interface with htslib in order to parse and iterate on bam files, make sure we're getting genome coordinates right, and make sure we're handling reverse complements correctly.

## Files

`reference.fasta` contains some simple sequences and `reads.fasta` shows different insertions, deletions, mismatches, and inversions of these sequences. `alignments.sorted.bam` is the result of aligning these reads to the reference with minimap2 and then sorting the bam file.

## Running

`make` will compile the test program `test_htslib_interface`.  Running the test program will show the mapping details for each test read, as well as details of the reference-to-query map.
