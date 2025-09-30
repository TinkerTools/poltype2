
lines = open('amoebaVdwType.dat').readlines()
type_in_dat = []
for line in lines:
  if '_Class_' in line:
    s = line.split()
    typ = s[1]
    if typ not in type_in_dat:
      type_in_dat.append(typ)

lines = open('../prm/amoebaVdw.prm').readlines()
type_in_prm = []
for line in lines:
  if '_Class_' in line:
    s = line.split()
    typ = s[0]
    if typ not in type_in_prm:
      type_in_prm.append(typ)

for typ in type_in_dat:
  if typ not in type_in_prm:
    print(typ)

for typ in type_in_prm:
  if typ not in type_in_dat:
    print(typ)