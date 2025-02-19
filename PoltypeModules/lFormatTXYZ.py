#!/usr/bin/env python

# this is to format tinker xyz file
# in order to let it be read by pymol

# Author: Chengwen Liu

import argparse

parser = argparse.ArgumentParser(description='Input the tinker xyz file(s)')
parser.add_argument('txyzs', type=str, nargs='+', help='tinker xyz files')
args = parser.parse_args()

for txyz in args.txyzs:
  lines = open(txyz).readlines()
  with open(txyz, 'w') as f:
    f.write(lines[0])
    for i in range(1,len(lines)):
      s = lines[i].split()
      atomid = s[0]
      element = s[1]
      x = float(s[2])
      y = float(s[3])
      z = float(s[4])
      tailstr = ''.join([f'{x:>6s}' for x in s[5:]])
      linestr = f"{atomid:>6s}{element:>3s}{x:12.6f}{y:12.6f}{z:12.6f}{tailstr}\n" 
      f.write(linestr)