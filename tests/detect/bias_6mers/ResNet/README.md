# DNAscent Test - 6mer Bias

This test checks whether there is a bias in the 6mers that DNAscent detects.  We want to make sure that DNAscent doesn't aggressively call BrdU just because a 6mer has a high number of thymidines in it.

## Files

`callsAgainstKL.py`

## Running

Run DNAscent detect, ideally on the 0-40-60-80-100% BrdU data released with the Nature Methods paper.  Usage for the script is:
`python callsAgainstTcontent.py /path/to/DNAscentDetect.out`
This will produce a bar plot, where each bar plot shows the average BrdU probability for 6mers with a given number of thymidines.
