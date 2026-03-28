import os
import sys

lines = open(sys.argv[1]).readlines()
for line in lines:
  if len(line.split()) != 0:
    if (not line.startswith('#')):
      smt = line.split()[0]  
      if not smt.startswith("[$("):
        idx_start = line.index(smt)
        idx_end = idx_start + len(smt)
        smt = "[$(" + smt + ")]"
        line = smt + '          ' + line[idx_end:-1] 
        print(line)
      else:
        print(line[:-1])
    else:
      print(line[:-1])
  else:
    print(line[:-1])