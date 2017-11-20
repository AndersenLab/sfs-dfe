

import sys

input_file="operon_F"
if len(sys.argv) > 1:
	input_file = sys.argv[1]

# Open MAXL file
with open("{}.sfs.MAXL.out".format(input_file), 'r') as f:
	maxl = f.read().strip()
maxl_out = dict([x.split(":") for x in maxl.split("\t")])

print(maxl_out)