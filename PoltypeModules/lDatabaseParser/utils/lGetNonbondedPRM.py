import os
import sys



usage = 'Usage: python lGetNonbondedPRM.py xyzfile prmfile'

if not len(sys.argv) == 3:
  sys.exit(usage)

xyzfile = sys.argv[1]
prmfile = sys.argv[2]


