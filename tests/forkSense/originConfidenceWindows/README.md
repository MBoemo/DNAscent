# DNAscent Test - Origin Confidence

DNAscent forkSense can output a bed file with replication origin locations, and specifies confidence windows of where these origins are located.  This test checks the width of these confidence windows in early S-phase cells where the confidence windows should be small.

## Files

`originWindowWidth.py`

## Running

Run DNAscent detect followd by forkSense with --markOrigins, ideally on the 1x BrdU cell cycle data released with the Nature Methods paper.  Usage for the script is:
`python originWindowWidth.py /path/to/origins_DNAscent_forkSense.bed`
This will plot a distribution of confidence window sizes.
