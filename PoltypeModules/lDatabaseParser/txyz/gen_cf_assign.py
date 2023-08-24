import os, sys

f = open("neutral.txt")
lines = f.readlines()
f.close()

for line in lines:
  terms = line.split()
  #print(f"echo {terms[0]}")
  #print(f"python ../lAssignAMOEBAplusPRM.py -xyz {terms[0]}.txyz -key tinker_final.prm -potent CF -fitting NO -new_para DATABASE -konly YES")
  print(f"cat {terms[0]}.type.cf >> total")
