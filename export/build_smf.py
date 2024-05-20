#!/usr/bin/env python

#quick script to build an SMF file for x variables

import sys

#take the input argument as the number of dimensions
x = sys.argv[1].strip()

#cast to an int
count = int(x)

#we will append the entire file to this array, then writelines the file at the end.
lines = ["Phase\n", "START\n", "Dimension\n"]

lines.append(x)
lines.append("\nNumber of LHS divisions\n")
lines.append("96\n")
lines.append("Upper bound\n")
lines.append("1.00000000 " * count)
lines.append("\nLower bound\n")
lines.append("0.00000000 " * count)
lines.append("\nKriging theta upper bound\n")
lines.append("400.00000000 " * count)
lines.append("\nKriging theta lower bound\n")
lines.append("0.00010000 " * count)
lines.append("\nSpacing\n")
lines.append("40 " * count)
lines.append("\n")
#write to the file
filename = "smf.inp"
with open(filename, 'w') as f:
    f.writelines(lines)




