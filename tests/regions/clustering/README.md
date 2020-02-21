# DNAscent Test - Clustering

This test checks the ability of DNAscent regions to find the T-to-B substitution rate in BrdU-positive regions.

## Files

`testClustering.py`

## Running

Set #define TEST_CLUSTERING 1 in DNAscent regions and recompile.  Run DNAscent regions on any run that has undergone a BrdU pulse and redirect stderr to file.txt.  Run `python testClustering.py file.txt`.
