import shutil

def assign_vdw_and_bonded(poltype):
    
    ## here we simply hack poltype to use BONDED and VDW parameters from DatabaseParser
    ## use lDatabaseParser to assign BONDED and VDW parameters
    ## Chengwen Liu
    ## Aug. 2023

    tmpxyz = 'tmpfilename.xyz'
    tmpkey = 'tmpfilename.key'
    shutil.copy(poltype.key4fname, tmpkey)
    shutil.copy(poltype.xyzoutfile, tmpxyz)
    cmd = f'python {poltype.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {poltype.molstructfname} -potent BONDED VDW'
    poltype.call_subsystem([cmd], True)
    
    tmpkey_b = 'tmpfilename.key_bonded'
    tmpkey_v = 'tmpfilename.key_vdw'
    vdwparams = {}
    bondparams = {}
    angleparams = {}
    strbndparams = {}
    opbendparams = {}
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
            line = vdwparams[ss[1]] + '\n'
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
   
    # replace key4 with the file we wanted
    # key4fname will be changed globally 
    # so we are good to go forward
    shutil.copy(tmpkey_2,poltype.key4fname)
    return 
