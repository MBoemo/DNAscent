# DNAscent Test - Event Alignment

This test checks the quality of DNAscent's adaptive banded event alignment and makes sure that there are no features of the alignment that could indicate a whole read should be thrown out or a call should be avoided at a particular position.

## Files

`plotEventAlignments.py`

## Running

Set #define TEST_ALIGNMENT 1 in detect.cpp and recompile.  Run DNAscent detect on 2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full barcode08 and barcode11, redirecting stderr to a file.  For each of these, run `python plotEventAlignments.py stderr.out DNAscentDetect.out`.
