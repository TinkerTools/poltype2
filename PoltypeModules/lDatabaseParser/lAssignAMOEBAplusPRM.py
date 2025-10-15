
#===================================
#        Chengwen Liu              #
#        Xudong Yang               #
#      liuchw2010@gmail.com        #
#      yang2076@utexas.edu         #    
#   University of Texas at Austin  #
#===================================

import sys
import argparse
import numpy as np
from openbabel import pybel
import networkx as nx
from valenceModule.fitting import *
from valenceModule.typing_tree import *
from valenceModule.typing_tree_assign import *
from valenceModule.modified_Seminario import *
from valenceModule.valence_utils import *

# color
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
ENDC = '\033[0m'

def write_tmp_key(txyz):
  atom2type = {}
  g = nx.Graph()
  nodes = []
  edges = []
  lines = open(txyz).readlines()
  if len(lines[0].split()) == 1:
    os.system(f"sed '1 s/$/ xxx/g' -i {txyz}")
  for line in lines[1:]:
    d = line.split()
    if int(d[0]) not in atom2type:
      atom2type[int(d[0])] = int(d[5])
    if d[0] not in nodes: 
      nodes.append(int(d[0]))
    for c in d[6:]:
      s = sorted([int(d[0]), int(c)])
      if s not in edges:
        edges.append(s)
  g.add_nodes_from(nodes)
  g.add_edges_from(edges)
 
  # find bonds and angles
  bonds = {} 
  angles = {} 
  for edge in g.edges:
    (n1, n2) = edge
    t1, t2 = atom2type[n1], atom2type[n2]
    comb = '-'.join([str(s) for s in sorted([t1,t2])])
    b = 1.0 
    if comb not in bonds:
      bonds[comb] = [b]
    else:
      bonds[comb] += [b]
  
  for node in g.nodes:
    adjs = list(g.adj[node])
    if len(adjs) > 1:
      for i in range(len(adjs)-1):
        for j in range(i+1, len(adjs)):
          n1, n2, n3 = adjs[i], node, adjs[j]
          t1, t2, t3 = atom2type[n1], atom2type[n2], atom2type[n3]
          if t1 > t3:
            t1, t3 = t3, t1
            n1, n3 = n3, n1
          a = 90.0 
          comb = '-'.join([str(s) for s in [t1,t2,t3]])
          if comb not in angles:
            angles[comb] = [a]
          else:
            angles[comb] += [a]
  
  # find torsions
  torsions = []
  for edge in g.edges:
    [a1, a2] = edge
    if (g.degree[a1] >=2) and (g.degree[a2] >= 2):
      a1_adjs = list(g.adj[a1])
      a2_adjs = list(g.adj[a2])
      for x in a1_adjs:
        for y in a2_adjs:
          if (x != a2) and (y != a1):
            if ([atom2type[x], atom2type[a1], atom2type[a2], atom2type[y]] not in torsions) \
            and ([atom2type[y], atom2type[a2], atom2type[a1], atom2type[x]] not in torsions) :
              torsions.append([atom2type[x], atom2type[a1], atom2type[a2], atom2type[y]])
  
  # find the tri- center
  tricentertypes = []
  for mol in pybel.readfile("txyz", txyz):
    natoms = len(mol.atoms)
    for i in range(natoms):
      hyb = mol.atoms[i].hyb
      if (hyb == 2) and (g.degree[i+1] == 3):
        if atom2type[i+1] not in tricentertypes:
          tricentertypes.append(str(atom2type[i+1]))
  # write
  with open(txyz.split(".")[0] + '.key', 'w') as f:
    tmpv = []
    for _, v in atom2type.items():
      if v not in tmpv:
        f.write(f'atom {v} {v} \n')
        tmpv.append(v)
      
    for bond in bonds:
      b = bond.split('-') 
      comb = '-'.join([str(s) for s in b])
      val = np.array(bonds[comb]).mean()
      f.write(f"bond   {b[0]}  {b[1]}  0.00 {val:10.4f}\n") 
    for angle in angles:
      a = angle.split('-')
      comb = '-'.join([str(s) for s in a])
      val = np.array(angles[comb]).mean()
      if a[1] in tricentertypes:
        f.write(f"anglep {a[0]}  {a[1]}  {a[2]}  0.00 {val:10.4f}\n")
      else:
        f.write(f"angle  {a[0]}  {a[1]}  {a[2]}  0.00 {val:10.4f}\n")
    
    for angle in angles:
      a = angle.split('-')
      f.write(f"strbnd {a[0]}  {a[1]}  {a[2]}  0.00 0.00\n") 
    tmp = []
    for bond in bonds:
      b = bond.split('-')
      if (b[0] in tricentertypes) and ([b[1], b[0]] not in tmp):
        f.write(f"opbend {b[1]}  {b[0]}  0  0  0.00  0.00\n")
        tmp.append([b[1], b[0]])
      if b[1] in tricentertypes and ([b[0], b[1]] not in tmp):
        f.write(f"opbend {b[0]}  {b[1]}  0  0  0.00  0.00\n")
        tmp.append([b[0], b[1]])
    for tor in torsions:
      [a,b,c,d] = tor
      f.write(f"torsion {a}  {b}  {c}  {d}  0.00        0.0   1   0.00        180.0   2    0.00   0.0   3\n")
  return

def genAtomType(txyz, key, potent):
  fname, _ = os.path.splitext(txyz)
  lines = open(txyz).readlines()
  if len(lines[0].split()) == 1:
    with open(txyz, "w") as f:
      f.write(lines[0].split("\n")[0] + " comments\n")
      for i in range(1,len(lines)):
        f.write(lines[i])
  atomnumbers, types = np.loadtxt(txyz, usecols=(0, 5,), unpack=True, dtype="str", skiprows=1)
  type_class_dict = {}
  atom_class_dict = {}
  for line in open(key).readlines():
    if "atom " in line:
      d = line.split()
      if d[1] in types:
        type_class_dict[d[1]] = d[2]
  for a,t in zip(atomnumbers, types):
    atom_class_dict[a] = type_class_dict[t]
  if sdf != None:
    inpfile = sdf
    if inpfile.endswith('mol2'):
        inpformat = 'mol2'
    else:
        inpformat = 'sdf'
    print(f"{YELLOW}Using SDF file to determine SMARTS types {ENDC}")
  else:
    inpfile = txyz
    inpformat = 'txyz'
    print(f"{YELLOW}Using tinker XYZ file to determine SMARTS types {ENDC}")
  for mol in pybel.readfile(inpformat,inpfile):
    matchDict = {}
    matchList = []
    commentsDict = {} 
    classesDict = {} 
    natoms = len(mol.atoms)
    for i in range(1, natoms+1, 1):
      matchList.append(i)
    matchDict = dict.fromkeys(matchList, 0)
    commentsDict = dict.fromkeys(matchList, 0)
    classesDict = dict.fromkeys(matchList, 0)
    potent_typefile_dict = {    "CF": ["amoebaplusCFluxType.dat", "amoebaplusCFluxType_general.dat"], \
                               "VDW": ["amoebaVdwType.dat"],  \
                             "POLAR": ["amoebaplusPolarType.dat"],  \
                            "BONDED": ["amoebaplusBondedType.dat"], \
                         "NONBONDED": ["amoebaplusNonbondedType.dat"]}
    type_files = potent_typefile_dict[potent]
    for type_file in type_files:
      lines = open(os.path.join(datfiledir, type_file)).readlines()
      for line in lines:
        if ("#" not in line[0]) and (len(line) > 10):
          data = line.split()
          myStr = data[0]
          classNum = data[2]
          className = line.split("# ")[0].split()[-1] 
          comment = line.split("# ")[-1][0:-1]
          smarts = pybel.Smarts(myStr)
          match = smarts.findall(mol)
          if match:
            for i in range(len(match)):
              matchDict[match[i][0]] = className
              commentsDict[match[i][0]] = comment
              classesDict[match[i][0]] = classNum
      out_type_file = f"{fname}.type.{potent.lower()}"
      if type_file == "amoebaplusCFluxType_general.dat":
        out_type_file = f"{fname}.type.{potent.lower()}_gen"
        
      with open(out_type_file, "w") as f:
        for atom in range(1, natoms+1, 1):
          atomtype = types[atom-1]
          atomclass = type_class_dict[atomtype]
          f.write("%5s %5s %5s %5s %5s %s\n"%(atom, atomtype, atomclass, matchDict[atom], classesDict[atom], commentsDict[atom]))
  return atom_class_dict

def assignPolar():
  genAtomType(xyz, key, 'POLAR')
  types, polars = [], []
  lines = open(os.path.join(prmfiledir,"polarize.prm")).readlines()
  for line in lines:
    dd = line.split()
    types.append(dd[1])
    polars.append(dd[2])
  smartspolarDict = dict(zip(types, polars))
  
  ttypes, stypes = np.loadtxt(f"{fname}.type.polar", usecols=(1,3), unpack=True, dtype="str")
  tinkerpolarDict = {}
  for t,s in zip(ttypes, stypes): 
    if t not in tinkerpolarDict:
      tinkerpolarDict[t] = smartspolarDict[s]
  lines = open(key).readlines()
  with open(key + "_polar", "w") as f:
    for line in lines:
      if "polarize " in line:
        dd = line.split()
        if dd[1] in ttypes:
          dd[2] = tinkerpolarDict[dd[1]]
          newline = "    ".join(dd) + "\n"
          oldline = "#" + line
          print(GREEN + "polarizability parameter found for %s"%dd[1] + ENDC)
          f.write(oldline)
          f.write(newline)
  
  
  # Match special polpair values
  # these are usually tough cases
  atomnumbers, types = np.loadtxt(xyz, usecols=(0, 5,), unpack=True, dtype="str", skiprows=1)
  atomnumbers = [int(a) for a in atomnumbers]
  atom_type_dict = dict(zip(atomnumbers, types)) 
  prmlines = open(os.path.join(prmfiledir,"amoebaPolarPair.prm")).readlines()
  # SDF is more info-rich 
  if sdf != None:
    inpfile = sdf
    inpformat = 'sdf'
  else:
    inpfile = xyz
    inpformat = 'txyz'
  
  # try to match line-by-line
  matched_polpairs = []
  comments = []
  tmp = []
  for mol in pybel.readfile(inpformat,inpfile):
    for line in prmlines:
      if ("#" not in line[0]) and (len(line) > 10):
        s = line.split()
        smt = s[0]
        polpair_params = '   '.join(s[2:4])
        polpair_comment = ' '.join(s[5:])
        smarts = pybel.Smarts(smt)
        matches = smarts.findall(mol)
        if matches != []:
          for match in matches:
            a = match[0] 
            polpair_type = atom_type_dict[a] 
            if (polpair_type not in tmp):
              tmp.append(polpair_type)
              polpair_prm_str = f"polpair {polpair_type} {polpair_params}"
              matched_polpairs.append(polpair_prm_str)
              comments.append(polpair_comment)
  
  # write the matched polpair parameters 
  with open(key + "_polar", "a") as f:
    if matched_polpairs != []:
      for c, p in zip(comments, matched_polpairs):
        print(GREEN + f"POLPAIR: {c}" + ENDC)
        f.write(f"# {c}\n")
        f.write(f"{p}\n")

  return True

def assignVdwAMOEBA():
  genAtomType(xyz, key, 'VDW')
  types = np.loadtxt(os.path.join(prmfiledir,"amoebaVdw.prm"), usecols=(0), unpack=True, dtype="str")
  
  ttypes, classs, stypes = np.loadtxt(f"{fname}.type.vdw", usecols=(1,2,3), unpack=True, dtype="str")
  ttype2class = dict(zip(ttypes, classs))
  
  vdws = []
  for line in open(os.path.join(prmfiledir,"amoebaVdw.prm")).readlines():
    ss = line.split()
    vdws.append('  '.join(ss[1:]))
  smartsvdwDict = dict(zip(types, vdws))
  tinkervdwDict = {}
  for t,s in zip(ttypes, stypes):
    if t not in tinkervdwDict:
      tinkervdwDict[t] = smartsvdwDict[s]
  lines = open(key).readlines()
  tmp = []
  with open(key + "_vdw", "w") as f:
    for ttype in ttypes:
      tclass = ttype2class[ttype]
      if tclass not in tmp:
        tmp.append(tclass)
        prmstr = tinkervdwDict[ttype]
        print(GREEN + "VdW parameters found for %s"%ttype + ENDC)
        newline = f"vdw {ttype2class[ttype]} {prmstr}\n"
        f.write(newline)
  
  # Match special vdwpair values
  # these are usually tough cases (Zn,Mg)
  atomnumbers, types = np.loadtxt(xyz, usecols=(0, 5,), unpack=True, dtype="str", skiprows=1)
  atomnumbers = [int(a) for a in atomnumbers]
  atom_type_dict = dict(zip(atomnumbers, types)) 
  prmlines = open(os.path.join(prmfiledir,"amoebaVdwPair.prm")).readlines()
  # SDF is more info-rich 
  if sdf != None:
    inpfile = sdf
    inpformat = 'sdf'
  else:
    inpfile = xyz
    inpformat = 'txyz'
  
  # try to match line-by-line
  matched_vdwpairs = []
  comments = []
  tmp = []
  for mol in pybel.readfile(inpformat,inpfile):
    for line in prmlines:
      if ("#" not in line[0]) and (len(line) > 10):
        s = line.split()
        smt = s[0]
        vdwpair_params = '   '.join(s[2:5])
        vdwpair_comment = ' '.join(s[6:])
        smarts = pybel.Smarts(smt)
        matches = smarts.findall(mol)
        if matches != []:
          for match in matches:
            a = match[0] 
            vdwpair_type = atom_type_dict[a] 
            if (vdwpair_type not in tmp):
              tmp.append(vdwpair_type)
              vdwpair_prm_str = f"vdwpair {vdwpair_type} {vdwpair_params}"
              matched_vdwpairs.append(vdwpair_prm_str)
              comments.append(vdwpair_comment)
  
  # write the matched vdwpair parameters 
  with open(key + "_vdw", "a") as f:
    if matched_vdwpairs != []:
      for c, p in zip(comments, matched_vdwpairs):
        print(GREEN + f"VdWPAIR: {c}" + ENDC)
        f.write(f"# {c}\n")
        f.write(f"{p}\n")
  return True


def assignGKAMOEBA():
  database_prm = os.path.join(prmfiledir,"amoebaGK.prm")
  prmlines = open(database_prm).readlines()
  
  # deal with txyz file so openbabel can read
  lines = open(xyz).readlines()
  if len(lines[0].split()) == 1:
    with open(xyz, "w") as f:
      f.write(lines[0].split("\n")[0] + " comments\n")
      for i in range(1,len(lines)):
        f.write(lines[i])
  atomnumbers, types = np.loadtxt(xyz, usecols=(0, 5,), unpack=True, dtype="str", skiprows=1)
  atom2types = dict(zip(atomnumbers, types))
  
  # SDF is more info-rich 
  if sdf != None:
    inpfile = sdf
    inpformat = 'sdf'
    print(f"{YELLOW}Using SDF file to determine SMARTS types {ENDC}")
  else:
    inpfile = xyz
    inpformat = 'txyz'
    print(f"{YELLOW}Using tinker XYZ file to determine SMARTS types {ENDC}")
  
  # try to match line-by-line
  type2prms = {} 
  type2comments = {} 
  tmp = []
  for mol in pybel.readfile(inpformat,inpfile):
    for line in prmlines:
      if ("#" not in line[0]) and (len(line) > 10):
        s = line.split()
        smt = s[0]
        prm = '   '.join(s[1:])

        smarts = pybel.Smarts(smt)
        matches = smarts.findall(mol)
        if matches != []:
          for match in matches:
            idx = str(match[0])
            typ = atom2types[idx]
            type2prms[typ] = prm
            type2comments[typ] = s[0]
  
  with open(key + '_gk', 'w') as f:
    for t,p in type2prms.items():
      comment = type2comments[t]
      f.write(f"# AMOEBA GK parameter matched from {comment}\n")
      f.write(f'solute {t} {p}\n')
      print(f"# AMOEBA GK parameter matched from {comment}")
      print(GREEN + f'solute {t} {p}\n' + ENDC)
  return

def assignNonbondedAMOEBAplus():
  genAtomType(xyz, key, 'NONBONDED')

  chgtrn_params = {}
  chgpen_params = {}
  vdw_params = {}

  lines = open(os.path.join(prmfiledir,"amoebaplusNonbonded.prm")).readlines()
  for line in lines:
    dd = line.split()
    if (len(dd) > 0) and ("CHGPEN" == dd[0].upper()):
      smt_type = dd[1]
      prm_str = '  '.join(dd[2:])
      if smt_type not in chgpen_params.keys():
        chgpen_params[smt_type] = prm_str
      else:
        print("Warning: there are two sets of chgpen parameters in the database")
        print(f"One is {prm_str}, and the other is {chgpen_params[smt_type]}")
    
    if (len(dd) > 0) and ("CHGTRN" == dd[0].upper()):
      smt_type = dd[1]
      prm_str = '  '.join(dd[2:])
      if smt_type not in chgtrn_params.keys():
        chgtrn_params[smt_type] = prm_str
      else:
        print("Warning: there are two sets of chgtrn parameters in the database")
        print(f"One is {prm_str}, and the other is {chgtrn_params[smt_type]}")
    
    if (len(dd) > 0) and ("VDW" == dd[0].upper()):
      smt_type = dd[1]
      prm_str = '  '.join(dd[2:])
      if smt_type not in vdw_params.keys():
        vdw_params[smt_type] = prm_str
      else:
        print("Warning: there are two sets of vdw parameters in the database")
        print(f"One is {prm_str}, and the other is {vdw_params[smt_type]}")

  ttypes, stypes = np.loadtxt(f"{fname}.type.nonbonded", usecols=(2,3), unpack=True, dtype="str")
  tinkerCPDict = {}
  tinkerCTDict = {}
  tinkerVDWDict = {}
  for t,s in zip(ttypes, stypes): 
    if t not in tinkerCPDict:
      tinkerCPDict[t] = chgpen_params[s]
      tinkerCTDict[t] = chgtrn_params[s]
      tinkerVDWDict[t] = vdw_params[s]
  lines = open(key).readlines()
  with open(key + "_nonbonded", "w") as f:
    for t in tinkerCPDict:
      line = "chgpen  %5s %s\n"%(t, tinkerCPDict[t])
      f.write(line)
    print(GREEN + "charge penetration parameters assigned from database"+ ENDC)
    for t in tinkerCTDict:
      line = "chgtrn  %5s %s\n"%(t, tinkerCTDict[t])
      f.write(line)
    print(GREEN + "charge transfer parameters assigned from database" + ENDC)
    for t in tinkerVDWDict:
      line = "vdw  %5s %s\n"%(t, tinkerVDWDict[t])
      f.write(line)
    print(GREEN+"van der Waals parameters assigned from database"+ENDC)
  return True

def assignCFlux_general():
  genAtomType(xyz, key, 'CF')
  # read in atom classes and short types (stype) 
  classs, stypes = np.loadtxt(f"{fname}.type.cf_gen", usecols=(2,3), unpack=True, dtype="str")
  class2stype = dict(zip(classs, stypes))
  
  # read in the database CFlux parameters
  # store two sets of parameters
  if cf4amoebaplus:
    lines = open(os.path.join(prmfiledir,"cflux2024_general.prm")).readlines()
  else:
    lines = open(os.path.join(prmfiledir,"cflux2023_general.prm")).readlines()
  stype2param = {}
  for line in lines:
    s = line.split()
    #bndcflux
    if len(s) == 2: 
      key0, value0 = s[0], f"{float(s[1]):10.4f}"
      k = key0.split("_")
      key1 = '_'.join([k[1], k[0]])
      value1 = "%10.4f"%(-float(s[1]))
      stype2param[key0] = value0
      stype2param[key1] = value1
    #angcflux
    if len(s) == 5: 
      key0, value0 = s[0], f"{float(s[1]):10.4f}{float(s[2]):10.4f}{float(s[3]):10.4f}{float(s[4]):10.4f}"
      k = key0.split("_")
      key1 = '_'.join([k[2], k[1], k[0]])
      value1 = f"{float(s[2]):10.4f}{float(s[1]):10.4f}{float(s[4]):10.4f}{float(s[3]):10.4f}"
      stype2param[key0] = value0
      stype2param[key1] = value1
  
  # assign parameters 
  lines = open(key).readlines()
  with open(key + "_cf","w") as f:
    for line in lines:
      if "bond " in line:
        d = line.split()
        if set(d[1:3]).issubset(set(classs)):
          s1 = class2stype[d[1]]
          s2 = class2stype[d[2]]
          comb = s1 + "_" + s2
          if d[1] == d[2]:
            stype2param[comb] = 0.0
          if comb in stype2param:
            f.write("bndcflux %s %s %s\n"%(d[1], d[2], stype2param[comb]))
            print(GREEN + "CFlux parameter assigned for bond %s-%s"%(d[1], d[2]) + ENDC)
          else:
            print(RED + f"CFlux parameter NOT found for bond %s-%s with type {comb}"%(d[1], d[2]) + ENDC)
            
            # do an average 
            jb = []
            s1_ = s1[0] + 'other'
            s2_ = s2[0] + 'other'
            if ("Cl" in s1[:2]) or ("Br" in s1[:2]):
              s1_ = s1[0:2] + 'other'
            if ("Cl" in s2[:2]) or ("Br" in s2[:2]):
              s2_ = s2[0:2] + 'other'
            for n1 in [s1, s1_]:
              for n2 in [s2, s2_]:
                combn = f"{n1}_{n2}"
                if (combn in stype2param):
                  p = stype2param[combn].split()
                  jb.append(float(p[0]))
            jb = "%10.4f"%np.array(jb).mean()
            if 'nan' in jb: 
              jb = '    0.0000'
            f.write("bndcflux %s %s %s\n"%(d[1], d[2], jb))
            print(YELLOW + f"CFlux parameter GENERATED for bond %s-%s"%(d[1], d[2]) + ENDC)
      
      if ("angle " in line) or ("anglep " in line):
        d = line.split()
        if set(d[1:4]).issubset(set(classs)):
          s1 = class2stype[d[1]]
          s2 = class2stype[d[2]]
          s3 = class2stype[d[3]]
          if int(d[1]) > int(d[3]):
            s1, s3 = s3, s1
            d[1], d[3] = d[3], d[1]
          comb = f"{s1}_{s2}_{s3}"
          if comb in stype2param:
            f.write("angcflux %s %s %s %s\n"%(d[1], d[2], d[3], stype2param[comb]))
            print(GREEN + "CFlux parameter assigned for angle %s-%s-%s"%(d[1], d[2], d[3]) + ENDC)
          else:
            print(RED + f"CFlux parameter NOT found for angle %s-%s-%s with type {comb}"%(d[1], d[2], d[3]) + ENDC)
            # do an average 
            jt1, jt2, jb1, jb2 = [], [], [], []
            s1_ = s1[0] + 'other'
            s2_ = s2[0] + 'other'
            s3_ = s3[0] + 'other'
            
            if ("Cl" in s1[:2]) or ("Br" in s1[:2]):
              s1_ = s1[0:2] + 'other'
            if ("Cl" in s2[:2]) or ("Br" in s2[:2]):
              s2_ = s2[0:2] + 'other'
            if ("Cl" in s3[:2]) or ("Br" in s3[:2]):
              s3_ = s3[0:2] + 'other'

            for n1 in [s1, s1_]:
              for n2 in [s2, s2_]:
                for n3 in [s3, s3_]:
                  combn = f"{n1}_{n2}_{n3}"
                  if combn in stype2param:
                    p = stype2param[combn].split()
                    jt1.append(float(p[0]))
                    jt2.append(float(p[1]))
                    jb1.append(float(p[2]))
                    jb2.append(float(p[3]))
            if jt1 != []:
              jt1 = "%10.4f"%np.array(jt1).mean()
              jt2 = "%10.4f"%np.array(jt2).mean()
              jb1 = "%10.4f"%np.array(jb1).mean()
              jb2 = "%10.4f"%np.array(jb2).mean()
            else:
              print(RED + "No similar angle found in database. Please assign by yourself" + ENDC)
              f.write("# Please assign the following ZEROS by yourself \n")
              jt1 = '0.0 '
              jt2 = '0.0 '
              jb1 = '0.0 '
              jb2 = '0.0 '
            f.write("angcflux %s %s %s %s\n"%(d[1], d[2], d[3], ''.join([jt1, jt2, jb1, jb2])))
            print(YELLOW + f"CFlux parameter GENERATED for angle %s-%s-%s"%(d[1], d[2], d[3]) + ENDC)
  return True

def assignCFlux(fname, key):
  genAtomType(xyz, key, 'CF')
  # read in atom classes and short types (stype) 
  atomclasses, databaseclasses = np.loadtxt(f"{fname}.type.cf", usecols=(2,4), unpack=True, dtype="str")
  class2cftype = dict(zip(atomclasses, databaseclasses))
  
  # read in the database CFlux parameters
  # store two sets of parameters
  lines = open(os.path.join(prmfiledir,"cflux2022.prm")).readlines()
  cftype2param = {}
  for line in lines:
    s = line.split()
    #bndcflux
    if len(s) == 4: 
      value0 = f"{float(s[3]):10.4f}"
      k0 = '_'.join([s[1], s[2]])
      value1 = "%10.4f"%(-float(s[3]))
      k1 = '_'.join([s[2], s[1]])
      cftype2param[k0] = value0
      cftype2param[k1] = value1
    #angcflux
    if len(s) == 8: 
      value0 = f"{float(s[4]):10.4f}{float(s[5]):10.4f}{float(s[6]):10.4f}{float(s[7]):10.4f}"
      k0 = '_'.join([s[1], s[2], s[3]])
      value1 = f"{float(s[5]):10.4f}{float(s[4]):10.4f}{float(s[7]):10.4f}{float(s[6]):10.4f}"
      k1 = '_'.join([s[3], s[2], s[1]])
      cftype2param[k0] = value0
      cftype2param[k1] = value1
  
  # assign parameters 
  lines = open(key).readlines()
  if_error = False
  with open(key + "_cf","w") as f:
    for line in lines:
      if "bond " in line:
        d = line.split()
        if set(d[1:3]).issubset(set(atomclasses)):
          s1 = class2cftype[d[1]]
          s2 = class2cftype[d[2]]
          comb = s1 + "_" + s2
          if d[1] == d[2]:
            cftype2param[comb] = 0.0
          if comb in cftype2param:
            f.write("bndcflux %s %s %s\n"%(d[1], d[2], cftype2param[comb]))
            print(GREEN + "CFlux parameter assigned for bond %s-%s"%(d[1], d[2]) + ENDC)
          else:
            print(RED + "CFlux parameter NOT found for bond %s-%s"%(d[1], d[2]) + ENDC)
            print(RED + "Try to find CF parameters for this molecule in general CF database")
            if_error = True
            break
      if ("angle " in line) or ("anglep " in line):
        d = line.split()
        if set(d[1:4]).issubset(set(atomclasses)):
          s1 = class2cftype[d[1]]
          s2 = class2cftype[d[2]]
          s3 = class2cftype[d[3]]
          if int(d[1]) > int(d[3]):
            s1, s3 = s3, s1
            d[1], d[3] = d[3], d[1]
          comb = f"{s1}_{s2}_{s3}"
          if comb in cftype2param:
            f.write("angcflux %s %s %s %s\n"%(d[1], d[2], d[3], cftype2param[comb]))
            print(GREEN + "CFlux parameter assigned for angle %s-%s-%s"%(d[1], d[2], d[3]) + ENDC)
          else:
            print(RED + "CFlux parameter NOT found for angle %s-%s-%s"%(d[1], d[2], d[3]) + ENDC)
            print(RED + "Try to find CF parameters for this molecule in general CF database")
            if_error = True
            break
  if(if_error):
    assignCFlux_general(fname, key)
  return True

def assignBonded(new_para_method, fitting = "NO"):
  atom2class = genAtomType(xyz, key, 'BONDED')
  # 2 methods to generate the new parameters: modified Seminario (Hessian); average by ranking tree (DATABASE)
  # you can choose to fit the new parameters to the given frequencies to improve the accuracy
  atomclasses, databaseclasses = np.loadtxt(f"{fname}.type.bonded", usecols=(2,4), unpack=True, dtype="str")
  tinker2database = dict(zip(atomclasses, databaseclasses))
  fitting_list = []
  para_strings_k = []
  para_strings_kbt = []
  hessian_mat = []
  coords = []

  if(new_para_method == "HESSIAN"):
    bond_list, angle_list, coords, hessian_mat = fchk_info(fname, len(atomclasses))
    eigenvalues, eigenvectors = eigen(len(atomclasses), hessian_mat)
    k_b_dict = bond_projection(bond_list, coords, eigenvalues, eigenvectors, 0.943)
    k_a_dict = bond_projection(angle_list, coords, eigenvalues, eigenvectors, 0.943)
  else:
    tree_1 = typing_tree()
    tree_1.read_ranking_file('typing_tree.log','amoebaplusBondedType.dat')
    tree_1.sorting_tree()

  # bond stretching
  # only assign force constant parameter, since equilibrium length will be from QM
  class1, class2 = np.loadtxt(os.path.join(prmfiledir,"bond.prm"), usecols=(1,2), unpack=True, dtype="str",skiprows=1)
  bondKs, bondLs = np.loadtxt(os.path.join(prmfiledir,"bond.prm"), usecols=(3,4), unpack=True, dtype="float",skiprows=1)
  classes = []
  for c1, c2 in zip(class1, class2):
    classes.append(c1 + "_" + c2)
  classBondParameterDict = dict(zip(classes, zip(bondKs, bondLs)))
  lines = open(key).readlines()
  idx = 0
  with open(key) as f:
    for line in lines[idx:]:
      d = line.split()
      if ("bond " in line):
        if (d[1] in atomclasses) and (d[2] in atomclasses):
          c1 = tinker2database[d[1]]
          c2 = tinker2database[d[2]]
          comb1 = c1 + "_" + c2
          comb2 = c2 + "_" + c1
          if comb1 in classBondParameterDict:
            para_strings_k.append("bond %s %s %10.4f %s\n"%(d[1], d[2], classBondParameterDict[comb1][0], d[4]))
            para_strings_kbt.append("bond %s %s %10.4f %10.4f\n"%(d[1], d[2], classBondParameterDict[comb1][0], classBondParameterDict[comb1][1]))
            print(GREEN + "BOND stretching parameter assigned for bond %s-%s"%(d[1], d[2]) + ENDC)
          elif comb2 in classBondParameterDict:
            para_strings_k.append("bond %s %s %10.4f %s\n"%(d[1], d[2], classBondParameterDict[comb2][0], d[4]))
            para_strings_kbt.append("bond %s %s %10.4f %10.4f\n"%(d[1], d[2], classBondParameterDict[comb2][0], classBondParameterDict[comb2][1]))
            print(GREEN + "BOND stretching parameter assigned for bond %s-%s"%(d[1], d[2]) + ENDC)
          else:
            is_wild = False
            if(new_para_method == 'DATABASE'):
              _, para, is_wild = typing_tree_assign(tree_1, 'b', comb1, classBondParameterDict)
            else: #HESSIAN
              para = modified_Seminario('b', d[1] + "_" + d[2], atomclasses, k_b_dict)

            if(fitting == 'NO'):
              if(is_wild): 
                para_strings_k.append("# Warning! This bond involves wild type. Please check your structure\n")
              para_strings_k.append("bond %s %s %10.4f %s\n"%(d[1], d[2], para[0], d[4])) # para[0]: k;  para[1]: b 
              para_strings_kbt.append("bond %s %s %10.4f %10.4f\n"%(d[1], d[2], para[0], para[1]))
              print(GREEN + "BOND stretching parameter (newly generated) assigned for bond %s-%s"%(d[1], d[2]) + ENDC)
            else:
              para_strings_k.append("bond %s %s PRM%d_ %s\n"%(d[1], d[2], len(fitting_list), d[4]))
              fitting_list.append(para)

  #angle bending

   ## find sp2 atoms
  sp2AtomTypes = findSp2AtomTypes(xyz, sdf) 
  class1, class2, class3  = np.loadtxt(os.path.join(prmfiledir, "angle.prm"), usecols=(1, 2, 3), unpack=True, dtype="str",skiprows=1)
  angleKs, angleTs  = np.loadtxt(os.path.join(prmfiledir, "angle.prm"), usecols=(4,5), unpack=True, dtype="float",skiprows=1)
  classes = []
  angKconsts = []
  angThetas = []
  '''store two sets of parameters since angle indices are interchangable''' 
  for c1, c2, c3 in zip(class1, class2, class3):
    classes.append(c1 + "_" + c2 + "_" + c3)
    classes.append(c3 + "_" + c2 + "_" + c1)
  for k,t in zip(angleKs, angleTs):
    angKconsts.append(k)
    angKconsts.append(k)
    angThetas.append(t)
    angThetas.append(t)

  classAngleParameterDict = dict(zip(classes, zip(angKconsts, angThetas)))
  
  with open(key) as f:
    for line in lines[idx:]:
      d = line.split()
      if ("angle " in line or "anglep " in line) and (len(d) > 4):
        angletype1 = d[1]
        angletype2 = d[2]
        angletype3 = d[3]
        if (angletype1 in atomclasses) and (angletype2 in atomclasses) and (angletype3 in atomclasses):
          c1 = tinker2database[angletype1]
          c2 = tinker2database[angletype2]
          c3 = tinker2database[angletype3]
          comb1 = c1 + "_" + c2 + "_" + c3
          comb2 = c3 + "_" + c2 + "_" + c1
          # change to anglep if angletype2 is sp2
          if angletype2 in sp2AtomTypes:
            d[0] = 'anglep'
          if (comb1 in classAngleParameterDict):
            para_strings_k.append("%s %s %s %s %10.4f %s\n"%(d[0], angletype1, angletype2, angletype3, classAngleParameterDict[comb1][0], d[5]))
            para_strings_kbt.append("%s %s %s %s %10.4f %10.4f\n"%(d[0], angletype1, angletype2, angletype3, classAngleParameterDict[comb1][0], classAngleParameterDict[comb1][1]))
            print(GREEN + "ANGLE bending parameter found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          elif (comb2 in classAngleParameterDict):
            para_strings_k.append("%s %s %s %s %10.4f %s\n"%(d[0], angletype3, angletype2, angletype1, classAngleParameterDict[comb2][0], d[5]))
            para_strings_kbt.append("%s %s %s %s %10.4f %10.4f\n"%(d[0],angletype3, angletype2, angletype1, classAngleParameterDict[comb2][0], classAngleParameterDict[comb2][1]))
            print(GREEN + "ANGLE bending parameter found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          else: 
            is_wild = False
            if(new_para_method == 'DATABASE'):
              _, para, is_wild = typing_tree_assign(tree_1, 'a', comb1, classAngleParameterDict)
            else: #'HESSIAN'
              para = modified_Seminario('a', d[1] + "_" + d[2] + "_" + d[3], atomclasses, k_a_dict)

            if(fitting == 'NO'):
              if(is_wild):
                para_strings_k.append("# Warning! This angle involves wild type. Please check your structure\n")
              para_strings_k.append("%s %s %s %s %10.4f %s\n"%(d[0], angletype1, angletype2, angletype3, para[0], d[5]))
              para_strings_kbt.append("%s %s %s %s %10.4f %s\n"%(d[0], angletype1, angletype2, angletype3, para[0], para[1]))
              print(GREEN + "ANGLE bending parameter (newly generated) assigned for angle %s-%s-%s"%(d[1], d[2], d[3]) + ENDC)
            else:
              para_strings_k.append("%s %s %s %s PRM%d_ %s\n"%(d[0], d[1], d[2], d[3], len(fitting_list), d[5]))
              fitting_list.append(para)

  #bond-angle coupling (strbnd term)
  class1, class2, class3  = np.loadtxt(os.path.join(prmfiledir, "strbnd.prm"), usecols=(1, 2, 3), unpack=True, dtype="str",skiprows=1)
  strbndK1, strbndK2  = np.loadtxt(os.path.join(prmfiledir, "strbnd.prm"), usecols=(4,5), unpack=True, dtype="float",skiprows=1)
  classes = []
  strbndKs = []
  '''store two sets of parameters since strbnd is asymetric''' 
  for c1, c2, c3 in zip(class1, class2, class3):
    classes.append(c1 + "_" + c2 + "_" + c3)
    classes.append(c3 + "_" + c2 + "_" + c1)
  for k1, k2 in zip(strbndK1, strbndK2):
    strbndKs.append([k1,k2])
    strbndKs.append([k2,k1])
  classStrbndKconstantDict = dict(zip(classes, strbndKs))
  
  with open(key) as f:
    for line in lines[idx:]:
      #if "strbnd " in line:
      # try to find strbnd parameters for every angle
      d = line.split()
      if ("strbnd " in line) and (len(d) > 4):
        angletype1 = d[1]
        angletype2 = d[2]
        angletype3 = d[3]
        if (angletype1 in atomclasses) and (angletype2 in atomclasses) and (angletype3 in atomclasses):
          c1 = tinker2database[angletype1]
          c2 = tinker2database[angletype2]
          c3 = tinker2database[angletype3]
          comb1 = c1 + "_" + c2 + "_" + c3
          comb2 = c3 + "_" + c2 + "_" + c1
          if (comb1 in classStrbndKconstantDict):
            tmp = "%10.4f%10.4f"%(classStrbndKconstantDict[comb1][0], classStrbndKconstantDict[comb1][1])
            para_strings_k.append("strbnd %s %s %s %s\n"%(angletype1, angletype2, angletype3, tmp))
            para_strings_kbt.append("strbnd %s %s %s %s\n"%(angletype1, angletype2, angletype3, tmp))
            print(GREEN + "STRBND coupling parameter found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          elif (comb2 in classStrbndKconstantDict):
            tmp = "%10.4f%10.4f"%(classStrbndKconstantDict[comb2][0], classStrbndKconstantDict[comb2][1])
            para_strings_k.append("strbnd %s %s %s %s\n"%(angletype3, angletype2, angletype1, tmp))
            para_strings_kbt.append("strbnd %s %s %s %s\n"%(angletype3, angletype2, angletype1, tmp))
            print(GREEN + "STRBND coupling parameter found for angle %s-%s-%s"%(angletype1, angletype2, angletype3) + ENDC)
          else: 
            _, para, is_wild = typing_tree_assign(tree_1, 'ba', comb1, classStrbndKconstantDict)
            if(fitting == 'NO'):
              if(is_wild):
                para_strings_k.append("# Warning! This strbnd involves wild type. Please check your structure\n")
              p0 = "%10.4f"%para[0] 
              p1 = "%10.4f"%para[1] 
              para_strings_k.append("strbnd %s %s %s %s %s\n"%(angletype1, angletype2, angletype3, p0, p1))
              para_strings_kbt.append("strbnd %s %s %s %s %s\n"%(angletype1, angletype2, angletype3, p0, p1))
              print(GREEN + "STRBND coupling parameter (newly generated) assigned for angle %s-%s-%s"%(d[1], d[2], d[3]) + ENDC)
            else:
              para_strings_k.append("strbnd %s %s %s PRM%d_ PRM%d_\n"%(d[1], d[2], d[3], len(fitting_list), len(fitting_list)))
              fitting_list.append(para)

  # out-of-plane bending 
  class1, class2, kopbends = np.loadtxt(os.path.join(prmfiledir,"opbend.prm"), usecols=(1,2,5), unpack=True, dtype="str",skiprows=1)
  classes = []
  for c1, c2 in zip(class1, class2):
    classes.append(c1 + "_" + c2)
  classOpbendKconstantDict = dict(zip(classes, kopbends))
  lines = open(key).readlines()
  with open(key) as f:
    for line in lines[idx:]:
      d = line.split()
      if ("opbend " in line) and (d[2] in sp2AtomTypes):
        if (d[1] in atomclasses) and (d[2] in atomclasses):
          c1 = tinker2database[d[1]]
          c2 = tinker2database[d[2]]
          comb = c1 + "_" + c2
          if comb in classOpbendKconstantDict:
            para_strings_k.append("opbend %s %s    0    0 %s\n"%(d[1], d[2], float(classOpbendKconstantDict[comb])))
            para_strings_kbt.append("opbend %s %s    0    0 %s\n"%(d[1], d[2], float(classOpbendKconstantDict[comb])))
            print(GREEN + "OPBEND parameter assigned for bond %s-%s-0-0"%(d[1], d[2]) + ENDC)
          else:
            _, para, is_wild = typing_tree_assign(tree_1, 'o', comb, classOpbendKconstantDict)
            if(is_wild):
              para_strings_k.append("# Warning! This opbend involves wild type. Please check your structure\n")
            para_strings_k.append("opbend %s %s    0    0 %s\n"%(d[1], d[2], para))
            para_strings_kbt.append("opbend %s %s    0    0 %s\n"%(d[1], d[2], para))
            print(GREEN + "OPBEND parameter (newly generated) assigned for bond %s-%s-0-0"%(d[1], d[2]) + ENDC)

  with open(key + "_bonded", "w") as f:
    if konly == "YES":
      for line in para_strings_k:
        f.write(line)
    else:
      for line in para_strings_kbt:
        f.write(line)
      
  if(fitting == 'YES'):
    with open("p0.txt", "w") as f:
      for line in fitting_list:
        f.write(line)
    fitting(fname)

  if os.path.isfile(key + "_bonded"):
    assignOpbendA09(atom2class, key+"_bonded")
    writeExplicitOpbend(key + "_bonded", xyz)
  return True

""" Assign amoeba18 opbend as a patch """
def assignOpbendA09(atom2class, keyfile):
  database_prm = os.path.join(prmfiledir,"opbend_a18.prm")
  prmlines = open(database_prm).readlines()
  
  # SDF is more info-rich 
  if sdf != None:
    inpfile = sdf
    inpformat = 'sdf'
    print(f"{YELLOW}Using SDF file to determine SMARTS types {ENDC}")
  else:
    inpfile = xyz
    inpformat = 'txyz'
    print(f"{YELLOW}Using tinker XYZ file to determine SMARTS types {ENDC}")
  
  # try to match line-by-line
  matched_opbs = {} 
  tmp = []
  for mol in pybel.readfile(inpformat,inpfile):
    for line in prmlines:
      if ("#" not in line[0]) and (len(line) > 10):
        s = line.split()
        smt = s[0]
        opb_params = s[1]
        opb_comment = line.split('%')[-1]
        idxs = line.split('%')[-2].split(',')
        smarts = pybel.Smarts(smt)
        matches = smarts.findall(mol)
        if matches != []:
          for match in matches:
            a1 = str(match[int(idxs[0])])
            a2 = str(match[int(idxs[1])])
            opb_classes = ' '.join([atom2class[a] for a in [a1, a2]])
            if (opb_classes not in tmp):
              tmp.append(opb_classes)
              opb_prm_str = f"opbend {opb_classes} 0  0  {opb_params}"
              matched_opbs[opb_classes] = opb_prm_str
  
  # replace opbend parameters if they 
  # have been matched from amoebabio18.prm
  if matched_opbs != {}:
    lines = open(keyfile).readlines()
    with open(keyfile, 'w') as f:
      for line in lines:
        if not 'opbend ' in line:
          f.write(line)
        else:
          s = line.split()
          # Make sure that line is not a comment 
          if (s[0] != '#') and (len(s) == 6):
            k = ' '.join(s[1:3])
            if k not in matched_opbs.keys():
              f.write(line)
            else:
              f.write(matched_opbs[k] + '\n')
  return

""" A helper function to write explicit opbend """
def writeExplicitOpbend(key, xyz):
  opbend_center = []
  opbend_params = {}
  
  prmlines = open(key).readlines()
  for line in prmlines:
    ss = line.split()
    if len(ss) == 6:
      if ("opbend " in line) and (ss[3] == '0') and (ss[4] == '0'):
        opbend_center.append(ss[2])
        opbend_params[f"{ss[1]}-{ss[2]}"] = f"{float(ss[5]):10.4f}"
  
  lines = open(xyz).readlines()
  atom2type = {}
  for line in lines[1:]:
    ss = line.split()
    atom2type[ss[0]] = ss[5]
  
  new_opbend_prmstrings = []
  for line in lines[1:]:
    ss = line.split()
    if (ss[5] in opbend_center) and (len(ss) == 9):
      atype1 = atom2type[ss[6]] 
      atype2 = atom2type[ss[7]] 
      atype3 = atom2type[ss[8]] 
      comb = atype1 + '-' + ss[5]
      prmstr = f"opbend {atype1} {ss[5]} {atype2} {atype3} {opbend_params[comb]}"
      if prmstr not in new_opbend_prmstrings:
        new_opbend_prmstrings.append(prmstr)
      comb = atype2 + '-' + ss[5]
      prmstr = f"opbend {atype2} {ss[5]} {atype1} {atype3} {opbend_params[comb]}"
      if prmstr not in new_opbend_prmstrings:
        new_opbend_prmstrings.append(prmstr)
      comb = atype3 + '-' + ss[5]
      prmstr = f"opbend {atype3} {ss[5]} {atype1} {atype2} {opbend_params[comb]}"
      if prmstr not in new_opbend_prmstrings:
        new_opbend_prmstrings.append(prmstr)
  
  # this will be in-place replacement
  with open(key, 'w') as f:
    for line in prmlines:
      if "opbend " not in line:
        f.write(line)
    for line in new_opbend_prmstrings:
      f.write(line + "\n")
  return

def assignTorsionAMOEBA():
  database_prm = os.path.join(prmfiledir,"amoebaTorsion.prm")
  prmlines = open(database_prm).readlines()
  
  # deal with txyz file so openbabel can read
  lines = open(xyz).readlines()
  if len(lines[0].split()) == 1:
    with open(xyz, "w") as f:
      f.write(lines[0].split("\n")[0] + " comments\n")
      for i in range(1,len(lines)):
        f.write(lines[i])
  atomnumbers, types = np.loadtxt(xyz, usecols=(0, 5,), unpack=True, dtype="str", skiprows=1)

  # construct atom2class dictionary
  atom_class_dict = {}
  type_class_dict = {}
  torsion_class = []
  for line in open(key).readlines():
    if "atom " in line:
      d = line.split()
      if d[1] in types:
        type_class_dict[d[1]] = d[2]
    
    if (len(line.split()) > 10) and ('torsion' in line.split()[0].strip()):
      d = line.split()
      torsion_class.append(' '.join(d[1:5]))
  for a,t in zip(atomnumbers, types):
    atom_class_dict[int(a)] = type_class_dict[t]
  
  # SDF is more info-rich 
  if sdf != None:
    inpfile = sdf
    inpformat = 'sdf'
    print(f"{YELLOW}Using SDF file to determine SMARTS types {ENDC}")
  else:
    inpfile = xyz
    inpformat = 'txyz'
    print(f"{YELLOW}Using tinker XYZ file to determine SMARTS types {ENDC}")
  
  # try to match line-by-line
  matched_torsions = {} 
  comments = {}
  for mol in pybel.readfile(inpformat,inpfile):
    for line in prmlines:
      if ("#" not in line[0]) and (len(line) > 10):
        s = line.split('%')
        smt = s[0]
        tor_indices = s[1].split(',')
        tor_params = s[2].strip().split(',')
        tor_comment = s[3].strip()

        smarts = pybel.Smarts(smt)
        matches = smarts.findall(mol)
        if matches != []:
          for match in matches:
            a1 = match[int(tor_indices[0]) - 1] 
            a2 = match[int(tor_indices[1]) - 1] 
            a3 = match[int(tor_indices[2]) - 1] 
            a4 = match[int(tor_indices[3]) - 1]
            tor_types = ' '.join([atom_class_dict[a] for a in [a1, a2, a3, a4]])
            tor_types_r = ' '.join([atom_class_dict[a] for a in [a4, a3, a2, a1]])
            tor_prm_str = f"torsion {tor_types} {tor_params[0]} 0.0 1 {tor_params[1]} 180.0 2 {tor_params[2]} 0.0 3"
            tor_prm_str_r = f"torsion {tor_types_r} {tor_params[0]} 0.0 1 {tor_params[1]} 180.0 2 {tor_params[2]} 0.0 3"
            if tor_types not in matched_torsions.keys():
              matched_torsions[tor_types] = [tor_prm_str]
              comments[tor_types] = [tor_comment]
            else:
              matched_torsions[tor_types] += [tor_prm_str]
              comments[tor_types] += [tor_comment]
            
            # make a copy for reverse order
            if tor_types_r not in matched_torsions.keys():
              matched_torsions[tor_types_r] = [tor_prm_str_r]
              comments[tor_types_r] = [tor_comment]
            else:
              matched_torsions[tor_types_r] += [tor_prm_str_r]
              comments[tor_types] += [tor_comment]
  
  # write the matched torsion
  with open(key + "_torsion", 'w') as f:
    for tor in torsion_class: 
      if tor in matched_torsions.keys():
        # use the last matched torsion
        prm = matched_torsions[tor][-1]
        cmt = comments[tor][-1]
        f.write(f"# Matched by lAssignAMOEBAplusPRM.py: {cmt}\n")
        f.write(f"{prm}\n")
  return

def assignTorsionAMOEBAplus():
  database_prm = os.path.join(prmfiledir,"amoebaplusTorsion.prm")
  prmlines = open(database_prm).readlines()
 
  stype2prm = {}
  stype2comment = {}
  for prmline in prmlines:
    if (prmline[0] != '#') and ('%' in prmline):
      s = prmline.split('%')
      stype = '-'.join(s[0].strip().split())
      prm = s[1].strip().split(',')
      comment = s[2].strip()
      stype2prm[stype] = prm
      stype2comment[stype] = comment 
  
  # read the key file to find atom types of torsion
  torsion_types = []
  keylines = open(key).readlines()
  for keyline in keylines:
    if ('torsion ' in keyline) and (keyline.split()[0].lower() == 'torsion'):
      s = keyline.split()
      if s[1:5] not in torsion_types:
        torsion_types.append(s[1:5])

  # read the non-bonded type file
  # this is pre-generated when running Poltype,
  # to assign nonbonded parameters
  typefile = xyz.split('.xyz')[0] + '.type.nonbonded'
  if not os.path.isfile(typefile):
    sys.exit(f'Error: {typefile} does not exist in {os.getcwd()}!')
  
  ttype2stype = {}
  lines = open(typefile).readlines()
  for line in lines:
    s = line.split()
    ttype = s[1]
    stype = s[3]
    if ttype not in ttype2stype.keys():
      ttype2stype[ttype] = stype

  with open(key + "_torsion", 'w') as f:
    for torsion_type in torsion_types:
      t1, t2, t3, t4 = torsion_type
      st1 = ttype2stype[t1]
      st2 = ttype2stype[t2]
      st3 = ttype2stype[t3]
      st4 = ttype2stype[t4]
      comb = '-'.join([st1, st2, st3, st4])
      comb_r = '-'.join([st4, st3, st2, st1])

      if comb in stype2prm.keys():
        prm = stype2prm[comb]
        comment = stype2comment[comb]
        f.write(f'# {comment}\n')
        f.write(f'torsion {t1} {t2} {t3} {t4} {prm[0]} 0.0 1 {prm[1]} 180.0 2 {prm[2]} 0.0 3\n')
      elif comb_r in stype2prm.keys():
        comment = stype2comment[comb_r]
        f.write(f'# {comment}\n')
        prm = stype2prm[comb_r]
        f.write(f'torsion {t1} {t2} {t3} {t4} {prm[0]} 0.0 1 {prm[1]} 180.0 2 {prm[2]} 0.0 3\n')
  return

if __name__ == "__main__":
  if len(sys.argv) == 1:
    sys.exit(RED + "please use '-h' option to see usage" + ENDC)
  parser = argparse.ArgumentParser()
  parser.add_argument('-xyz', dest = 'xyz', help = "tinker xyz file", required=True)  
  parser.add_argument('-key', dest = 'key', help = "tinker prm file", default=None)  
  parser.add_argument('-sdf', dest = 'sdf', help = "sdf file [optional]", default=None)  
  parser.add_argument('-potent', dest = 'potent', help = "potential energy terms", nargs='+', type=str.upper, default=[])  
  parser.add_argument('-fitting', dest = 'fitting', help = "fit the frequencies if new parameters are needed", default="NO", type=str.upper, choices=["YES", "NO"])
  parser.add_argument('-new_para', dest = 'new_para', help = "method to generate the new valence parameters", default="DATABASE", type=str.upper, choices=["HESSIAN", "DATABASE"])  
  parser.add_argument('-konly', dest = 'konly', help = "assign force constant only for valence parameters", default="YES", type=str.upper, choices=["YES", "NO"])  
  parser.add_argument('-cf4amoebaplus', dest = 'cf4amoebaplus', help = "apply special treatment of CF parms. for amoebaplus", default=False, type=bool)  
  args = vars(parser.parse_args())
  global xyz, key, sdf, cf4amoebaplus 
  xyz = args["xyz"]
  key = args["key"]
  sdf = args["sdf"]
  cf4amoebaplus = args["cf4amoebaplus"]
  if key == None:
    write_tmp_key(xyz)
    key = xyz.split('.')[0] + '.key'
  potent = args["potent"]
  new_para = args["new_para"]
  fitting = args["fitting"]
  global konly
  konly = args["konly"]
  global prmfiledir, datfiledir
  rootdir = os.path.join(os.path.split(__file__)[0])
  prmfiledir = os.path.join(rootdir, 'prm')
  datfiledir = os.path.join(rootdir, 'dat')

  fname, _ = os.path.splitext(xyz)
 
  if potent == []:
    sys.exit(RED + "potent argument must be speficied!" + ENDC)
  for p in potent:
    if (p == "POLAR"):
      assignPolar()
    elif (p == "CF"):
      assignCFlux_general() 
    elif (p == "VDW"):
      assignVdwAMOEBA() 
    elif (p == "GK"):
      assignGKAMOEBA() 
    elif (p == "BONDED"):
      assignBonded(new_para, fitting)
    elif (p == "NONBONDED"):
      assignNonbondedAMOEBAplus() 
    elif (p == "TORSION"):
      assignTorsionAMOEBA() 
    elif (p == "TORSION+"):
      assignTorsionAMOEBAplus() 
    else:
      print(RED + f"{p} term not supported!" + ENDC)
