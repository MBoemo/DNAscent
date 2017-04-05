#! /usr/bin/env python
import pysam
import sys

sam_file = pysam.Samfile(sys.argv[1])
out_files = list()

# open an output file for each reference sequence
for x in sam_file.references:
	print x
	out_files.append(pysam.Samfile(x + ".bam", "wb", template=sam_file))

for record in sam_file:
	ref_length = sam_file.lengths[record.reference_id]

	if record.aend is None or record.query_alignment_length is None:
		continue

	ref_cover = float(record.aend - record.pos) / ref_length
	query_cover = float(record.query_alignment_length) / record.query_length

	if ref_cover > 0.8 and query_cover > 0.8 and record.is_reverse == False:
		out_files[record.reference_id].write(record)
