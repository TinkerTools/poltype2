import shutil
from rdkit import Chem

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

def assign_nonbonded_and_bonded(poltype):
    # Here we simply replace the existing bonded and nonbonded 
    # parameters in the key_4 file by using assignment script
    # Chengwen Liu
    # Aug. 2023

    tmpxyz = 'tmpfilename.xyz'
    tmpkey = 'tmpfilename.key'
    shutil.copy(poltype.key4fname, tmpkey)
    shutil.copy(poltype.xyzoutfile, tmpxyz)
    if poltype.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
      cmd = f'python {poltype.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {poltype.molstructfname} -potent BONDED NONBONDED CF'
    else:
      cmd = f'python {poltype.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {poltype.molstructfname} -potent BONDED VDW'
    poltype.call_subsystem([cmd], True)
    
    tmpkey_b = 'tmpfilename.key_bonded'
    tmpkey_v = 'tmpfilename.key_vdw'
    vdwparams = {}
    bondparams = {}
    angleparams = {}
    strbndparams = {}
    opbendparams = {}
    
    # vdw parameters of AMOEBAplus model is in nonbonded
    if poltype.forcefield.upper() not in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
      for line in open(tmpkey_v).readlines():
        if "#" not in line[0:10]:
          ss = line.split()
          vdwparams[ss[1]] = '   '.join(ss)
    
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
    tmpkey_2 = 'tmpfilename.key_2'
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
          if (ss[0].lower() == 'opbend'):
            comb = '-'.join(ss[1:5])
            if comb in opbendparams.keys():
              line = opbendparams[comb] + '\n'
        
        f.write(line)
    
    if poltype.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
      lines = open(tmpkey_2).readlines()
      # vdws that already exist are commented out
      with open(tmpkey_2, 'w') as f:
        for line in lines:
          if (len(line.split()) > 3) and (line.split()[0].lower() == 'vdw'):
            line = '# ' + line
          f.write(line)

      tmpkey_cflux = 'tmpfilename.key_cf'
      tmpkey_nonbonded = 'tmpfilename.key_nonbonded'
      cflux_lines = open(tmpkey_cflux).readlines()
      nonbonded_lines = open(tmpkey_nonbonded).readlines()
      with open(tmpkey_2, 'a') as f:
        for line in cflux_lines: 
          f.write(line)
        for line in nonbonded_lines:
          # since chgpen parameters have been assigned to key_2
          if "chgpen " not in line:
            f.write(line)
    # replace key4 with the file we wanted
    shutil.copy(tmpkey_2,poltype.key4fname)
    return 
