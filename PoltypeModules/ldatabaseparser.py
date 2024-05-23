

""" 
  This program is used to assign valence and non-bonded interactions for AMOEBA and AMOEBA+ model
  The torsion parameters, are assigned by another program called torsiondatabaseparser.py
  In the long run, we will merge torsion assignment into this program too.
  - Chengwen Liu
  - Feb 2024
"""
import os
import shutil
import numpy as np
import networkx as nx
from rdkit import Chem
from openbabel import pybel 

def assign_chgpen_params(poltype):
  # same as assigning polarizability parameters,
  # charge penetration parameters are assigned on-the-fly
  # because the parameters are used for next-step ESP fitting 
  # Chengwen Liu
  # Feb 2024
  xyzfile = poltype.xyzoutfile
  sdffile = poltype.molstructfname
  datfile = poltype.amoebaplussmallmoleculenonbonded_dat
  prmfile = poltype.amoebaplussmallmoleculenonbonded_prm
  
  rawsmarts2chgpen = {}
  rawsmarts2name = {}
  name2chgpen = {}
  
  ordered_smarts = []
  lines = open(datfile).readlines()
  for line in lines:
    if (line[0] != '#') and len(line.strip()) > 10:
      d = line.split()
      rawsmarts2name[d[0]] = d[1]
      ordered_smarts.append(d[0])
  
  lines = open(prmfile).readlines()
  for line in lines:
    if (line[0] != '#'):
      d = line.split()
      if d[0].upper() == 'CHGPEN':
        name2chgpen[d[1]] = '  '.join(d[2:4])
  
  for key, name in rawsmarts2name.items():
    # it is possible that we have more smarts
    # defined but the parameters are not ready yet
    if name in name2chgpen.keys():
      rawsmarts2chgpen[key] = name2chgpen[name]
  atom2chgpen = {}
  rdkitmol = Chem.MolFromMolFile(sdffile,removeHs=False)
  atoms = range(1, len(rdkitmol.GetAtoms()) + 1)
  
  for smt in ordered_smarts:
    if smt in rawsmarts2chgpen.keys():
      chgpen = rawsmarts2chgpen[smt]
      pattern = Chem.MolFromSmarts(smt)
      match = rdkitmol.GetSubstructMatches(pattern)
      if match:
        for i in range(len(match)):	
          atom2chgpen[match[i][0]+1] = chgpen 

  lines = open(xyzfile).readlines()
  atom2type = {}
  for line in lines[1:]:
    dd = line.split()
    atom2type[dd[0]] = dd[5]
  
  types = []
  chgpen_params = []
  for atom, atype in atom2type.items():
    if atype not in types:
      chgpen_params.append(f'chgpen {atype} {atom2chgpen[int(atom)]}')
      types.append(atype)
  
  key2file = poltype.key2fnamefromavg
  with open(key2file, 'a') as f:
    f.write("#\n")
    f.write("# Charge penetration parameters assigned from database\n")
    f.write("#\n")
    f.write("FIX-CHGPEN \n")
    for prm in chgpen_params:
      f.write(prm + '\n')
  return 

# Print the zero parameters for bond/angle/strbnd/opbend
# Given the tinker xyz file
# Assumption: atom type == atom class
# helper function
def calc_bond_angle(coords):
  
  # bond
  # coords = [[x1,y1,z1], [x2,y2,z2]]
  result = 0
  ndata = len(coords)
  if ndata == 2:
    coord1 = np.array(coords[0])
    coord2 = np.array(coords[1])
    result = np.sqrt(np.square(coord1-coord2).sum())
  
  # angle
  # coords = [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
  if ndata == 3:
    coord1 = np.array(coords[0])
    coord2 = np.array(coords[1])
    coord3 = np.array(coords[2])
    vec21 = coord1 - coord2 
    vec23 = coord3 - coord2
    dot = np.dot(vec21, vec23)
    vec21norm = np.linalg.norm(vec21) 
    vec23norm = np.linalg.norm(vec23) 
    result = 180.0/np.pi * (np.arccos(dot/(vec21norm*vec23norm)))
  return result

def write_initial_parameters(sdffile, txyz):
  atom2type = {}
  g = nx.Graph()
  nodes = []
  edges = []
  lines = open(txyz).readlines()
  coords = []
  if len(lines[0].split()) == 1:
    os.system(f"sed '1 s/$/ xxx/g' -i {txyz}")
  for line in lines[1:]:
    d = line.split()
    coords.append([float(d[2]), float(d[3]), float(d[4])])
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
    b = calc_bond_angle([coords[n1-1], coords[n2-1]])
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
          a = calc_bond_angle([coords[n1-1], coords[n2-1], coords[n3-1]])
          comb = '-'.join([str(s) for s in [t1,t2,t3]])
          if comb not in angles:
            angles[comb] = [a]
          else:
            angles[comb] += [a]
  
  # find the tri- center
  tricentertypes = []
  # rdkit is more robust
  aniline_nitrogens = []
  rdkitmol = Chem.MolFromMolFile(sdffile,removeHs=False)
  pattern = Chem.MolFromSmarts('[NX3][a]')
  matches = rdkitmol.GetSubstructMatches(pattern)
  for match in matches:
    nitrogen = str(match[0] + 1)
    aniline_nitrogens.append(nitrogen)
  
  amide_nitrogens = []
  pattern = Chem.MolFromSmarts('[NX3][C]=[O]')
  matches = rdkitmol.GetSubstructMatches(pattern)
  for match in matches:
    nitrogen = str(match[0] + 1)
    amide_nitrogens.append(nitrogen)
  
  #!! dont forget to add this to valence_utils.py 
  
  excluded_nitrogens = aniline_nitrogens + amide_nitrogens
  num_atoms = rdkitmol.GetNumAtoms() 
  for i in range(num_atoms):
    atom = rdkitmol.GetAtomWithIdx(i)
    if (atom.GetHybridization() == Chem.HybridizationType.SP2) and (g.degree[i+1] == 3) and (str(i+1) not in excluded_nitrogens):
      if str(atom2type[i+1]) not in tricentertypes: 
        tricentertypes.append(str(atom2type[i+1]))
  
  with open("valence_init.dat", 'w') as f:
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
  return

def zero_special_torsions(poltype):
  xyzfile = poltype.xyzoutfile
  sdffile = poltype.molstructfname
  rdkitmol = Chem.MolFromMolFile(sdffile,removeHs=False)
  
  atom2type = {}
  lines = open(xyzfile).readlines()[1:]
  for line in lines:
    s = line.split()
    atom2type[s[0]] = s[5]
  
  zero_out_torsions = []
  smt = "[*]#[*]~[*]~[*]"
  pattern = Chem.MolFromSmarts(smt)
  matches = rdkitmol.GetSubstructMatches(pattern)
  for match in matches:
    a, b, c, d = match
    a = str(a + 1)
    b = str(b + 1)
    c = str(c + 1)
    d = str(d + 1)
  
    a_type = atom2type[a]
    b_type = atom2type[b]
    c_type = atom2type[c]
    d_type = atom2type[d]
  
    zero_out_torsions.append('-'.join([a_type,b_type,c_type,d_type]))
    zero_out_torsions.append('-'.join([d_type,c_type,b_type,a_type]))

  smt = "[*]~[*]#[*]~[*]"
  pattern = Chem.MolFromSmarts(smt)
  matches = rdkitmol.GetSubstructMatches(pattern)
  for match in matches:
    a, b, c, d = match
    a = str(a + 1)
    b = str(b + 1)
    c = str(c + 1)
    d = str(d + 1)
  
    a_type = atom2type[a]
    b_type = atom2type[b]
    c_type = atom2type[c]
    d_type = atom2type[d]
  
    zero_out_torsions.append('-'.join([a_type,b_type,c_type,d_type]))
    zero_out_torsions.append('-'.join([d_type,c_type,b_type,a_type]))
    
  if zero_out_torsions != []:  
    poltype.WriteToLog("Zeroing Torsion Parameters Involving Triple Bond")
    # write ZERO parameters in key4fname
    lines = open(poltype.key4fname).readlines()
    with open(poltype.key4fname, 'w') as f:
      for line in lines:
        if ('torsion ' in line) and (len(line.split()) > 5):
          s = line.split()
          if '-'.join(s[1:5]) in zero_out_torsions:
            line = f'torsion {" ".join(s[1:5])} 0.000 0.0 1 0.000 180.0 2 0.000 0.0 3\n' 
        f.write(line)
    
  return 

def assign_bonded_params(poltype):
    tmpxyz = poltype.xyzoutfile
    sdffile = poltype.molstructfname
    write_initial_parameters(sdffile, tmpxyz)
    tmpkey = 'tmpbonded.key'
    sdffile = poltype.molstructfname
    shutil.copy(poltype.key4fname, tmpkey)
    catcmd = f"cat valence_init.dat >> {tmpkey}"
    poltype.call_subsystem([catcmd], True)
    cmd = f'python {poltype.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {poltype.molstructfname} -potent BONDED'
    poltype.call_subsystem([cmd], True)
    
    tmpkey_b =  tmpkey + '_bonded'
    bondparams = {}
    angleparams = {}
    strbndparams = {}
    opbendparams = {}
    
    for line in open(tmpkey_b).readlines():
      if "#" not in line[0:10]:
        ss = line.split()
        if ss[0].lower() == 'bond':
          comb = '-'.join(ss[1:3])
          bondparams[comb] = '   '.join(ss)
        if ss[0].lower() in ['angle', 'anglep']:
          comb = '-'.join(ss[1:4])
          angleparams[comb] = '   '.join(ss)
        if ss[0].lower() == 'strbnd':
          comb = '-'.join(ss[1:4])
          strbndparams[comb] = '   '.join(ss)
        if ss[0].lower() == 'opbend':
          comb = '-'.join(ss[1:5])
          opbendparams[comb] = '   '.join(ss)
    
    lines = open(tmpkey).readlines()
    tmpkey_2 = 'tmpbonded.key_2'
    with open(tmpkey_2, 'w') as f:
      for line in lines:
        ss = line.split()
        if len(ss) > 3:
          if (ss[0].lower() == 'bond'):
            comb = '-'.join(ss[1:3])
            if comb in bondparams.keys():
              line = bondparams[comb] + '\n'
          if (ss[0].lower() in ['angle', 'anglep']):
            comb = '-'.join(ss[1:4])
            if comb in angleparams.keys():
              line = angleparams[comb] + '\n'
          if (ss[0].lower() == 'strbnd'):
            comb = '-'.join(ss[1:4])
            if comb in strbndparams.keys():
              line = strbndparams[comb] + '\n'
        # opbend x y 0 0
        if 'opbend ' not in line:
          f.write(line)
      # opbend x y z w
      for opbprm in opbendparams.values():
        f.write(opbprm + '\n')
    shutil.copy(tmpkey_2,poltype.key4fname)
    return 

def assign_nonbonded_params(poltype):
    tmpxyz = poltype.xyzoutfile 
    tmpkey = 'tmpnonbonded.key'
    sdffile = poltype.molstructfname
    shutil.copy(poltype.key4fname, tmpkey)
    if poltype.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
      cmd = f'python {poltype.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {poltype.molstructfname} -potent NONBONDED CF'
    else:
      cmd = f'python {poltype.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {poltype.molstructfname} -potent VDW'
    poltype.call_subsystem([cmd], True)
    
    tmpkey_v = f'{tmpkey}_vdw'
    vdwparams = {}
    
    # vdw parameters of AMOEBAplus model is in nonbonded
    if poltype.forcefield.upper() not in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
      for line in open(tmpkey_v).readlines():
        if "#" not in line[0:10]:
          ss = line.split()
          vdwparams[ss[1]] = '   '.join(ss)
    
    lines = open(tmpkey).readlines()
    tmpkey_2 = f'{tmpkey}_2'
    with open(tmpkey_2, 'w') as f:
      for line in lines:
        ss = line.split()
        if len(ss) > 3:
          if (ss[0].lower() == 'vdw') and (ss[1] in vdwparams.keys()):
            # vdw parameters of AMOEBAplus model is in nonbonded
            if poltype.forcefield.upper() not in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
              line = vdwparams[ss[1]] + '\n'
            else:
              line = '\n'
        f.write(line)
    
    if poltype.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
      tmpkey_cflux = f'{tmpkey}_cf'
      tmpkey_nonbonded = f'{tmpkey}_nonbonded'
      cflux_lines = open(tmpkey_cflux).readlines()
      nonbonded_lines = open(tmpkey_nonbonded).readlines()
      with open(tmpkey_2, 'a') as f:
        for line in cflux_lines: 
          f.write(line)
        for line in nonbonded_lines:
          # since chgpen parameters have been assigned to key_2
          if "chgpen " not in line:
            f.write(line)
    
    # for AMOEBA model, only vdw term exists
    else:
      tmpkey_vdw = f'{tmpkey}_vdw'
      vdw_lines = open(tmpkey_vdw).readlines()
      with open(tmpkey_2, 'a') as f:
        for line in vdw_lines: 
          f.write(line)
    
    # replace key4 with the file we wanted
    shutil.copy(tmpkey_2,poltype.key4fname)
    return 
