#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

import sys

f = open(sys.argv[1],'r')
fout = open('fixedFoh.foh','w')

first = True

for line in f:
	if line[0] == '>':
		if not first:
			fout.write('<\n')
			fout.write(line)
		else:
			fout.write(line)
			first = False
	else:
		fout.write(line)

f.close()
fout.close()
