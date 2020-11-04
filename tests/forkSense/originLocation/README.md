# DNAscent Test - Origin Calling Accuracy

This test checks the accuracy of forkSense origin calls by making sure they are generally close to known origins in S. cerevisiae.

## Files

`distanceFromNearestOrigin.py`

This will also require a bed file of all origins from oridb (http://cerevisiae.oridb.org/).

## Running

Run DNAscent detect followd by forkSense with --markOrigins, ideally on the 1x BrdU cell cycle data released with the Nature Methods paper.  Usage for the script is:
`python distanceFromNearestOrigin.py /path/to/origins_DNAscent_forkSense.bed`
This will show a distribution of distances from DNAscent origin calls to the nearest known origin.
