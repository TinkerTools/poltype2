import os
import numpy as np

files = os.listdir()
with open('filelist', 'w') as fo:
  for f in files:
    if f.endswith(".txyz"):
      atoms = list(np.loadtxt(f, usecols=(1), unpack=True, dtype = 'str', skiprows=1))
      atoms = sorted(atoms)
      fo.write(''.join(atoms) + '  ' + f + '\n')