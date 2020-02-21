# DNAscent Test - 6mer Bias

This test checks whether there is a bias in the 6mers that DNAscent detects.  We want to make sure that DNAscent is using context effects correctly so that even if a particular BrdU-containing 6mer doesn't shift the signal much from the thymidine-only case, using the flanking sequence still means that this 6mer can be detected.

## Files

`callsAgainstKL.py`

## Running

Run DNAscent detect, ideally on the 0-40-60-80-100% BrdU data released with the Nature Methods paper.  Usage for the script is:
`python callsAgainstKL.py /path/to/DNAscentDetect.out /path/to/DNAscent/pore_models/BrdU_full_noThreshold.model`
This will produce a scatter plot, where each point is a different 6mer.  The x-axis is the fraction of times the 6mer was called as BrdU, and the y-axis is the KL-divergence of the BrdU-containing 6mer against the thymidine-only pore model.  If well-behaved, the scatter plot should look somewhat like a vertical bar, which is to say regardless of the KL-divergence, we're still positively identifying each 6mer with about the same frequency. 
