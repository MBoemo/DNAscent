# DNAscent Test - Test HMM Positive Calls

This test looks at events that were called as BrdU and checks them against the span on the query sequence and the number of events passed to the HMM to determine whether we can add any meaningful QCs.

## Files

`testPositiveCalls.py`

## Running

Run DNAscent detect with #define TEST_LL 1 and redirect stderr to file.txt.  Then run `python testPositiveCalls.py file.txt` to produce the plots.
