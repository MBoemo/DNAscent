# DNAscent Test - Test HMM Forward Algorithm

This test crosschecks the implementation of the HMM forward algorithm in DNAscent detect against pomegranate (https://github.com/jmschrei/pomegranate) to make sure they agree.

## Files

`testHMMForward.py`
`testHMMProbeViterbi.py`
Note that pomegranate (https://github.com/jmschrei/pomegranate) also needs to be installed.

## Running

Run DNAscent detect on any sample with #define TEST_HMM 1.  Redirect stderr to file.txt and then run `python testHMMForward.py file.txt' to check agreement.  From left to right, the columns indicate:
* whether BrdU was used in the HMM (0 for thymidine only, 1 for BrdU),
* the log probability of the events from pomegranate,
* the log probability of the events from DNAscent detect.
On that same output file.txt, if it was run on 0% BrdU reads, run `python testHMMProbeViterbi.py` to see plots of Viterbi insertions and deletions and how they correlate with the log likelihood ratio.  This is meant to look for any meaningful QCs that could be added.
