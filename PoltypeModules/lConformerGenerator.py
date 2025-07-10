
# This program 
# - takes the .sdf input structure
# - uses RDKit to generate 500 conformers
# - optimizes the conformers with distance constraints with MMFF94 in RDKit
# - selects the extended conformer with Rg, SASA, number of Intramolecular HB, and steric interaction criteria
# - further optimize the selected extended conformer with xtb to deal with MMFF95 deficiency (puckered Nitrogen etc)
# - generates "conftest.mol" which poltype wants

# Author: Chengwen Liu
# Date: Mar. 2023
#       Feb. 2025
import os
import sys
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors3D,rdFreeSASA,rdmolfiles,ChemicalForceFields


def find_intramolecular_hbonds(mol, confId = -1, donorAtoms = [7,8,9], distTol = 2.5):
  res = []
  conf = mol.GetConformer(confId)
  for i in range(mol.GetNumAtoms()):
    atomi = mol.GetAtomWithIdx(i)
    if atomi.GetAtomicNum()==1:
      if atomi.GetDegree() != 1:
        continue
      nbr = atomi.GetNeighbors()[0]
      if nbr.GetAtomicNum() not in donorAtoms:
        continue
      for j in range(mol.GetNumAtoms()):
        if j==i:
          continue
        atomj = mol.GetAtomWithIdx(j)
        if atomj.GetAtomicNum() not in donorAtoms or mol.GetBondBetweenAtoms(i,j):
          continue
        dist = (conf.GetAtomPosition(i)- conf.GetAtomPosition(j)).Length()
        if dist<distTol:
          res.append((i,j,dist))
  return res

def get_pairAtoms(mol):
  pair_1234 = []
  for bond in mol.GetBonds():
    t1 = bond.GetBeginAtomIdx()
    t2 = bond.GetEndAtomIdx()
    comb = [t1, t2]
    # 1-2 connected
    if sorted(comb) not in pair_1234:
      pair_1234.append(sorted(comb))
    
    for neigh in mol.GetAtomWithIdx(t1).GetNeighbors():
      neig_idx = neigh.GetIdx()
      if neig_idx != t2:
        comb = [neig_idx, t2]
        # 1-3 connected
        if sorted(comb) not in pair_1234:
          pair_1234.append(sorted(comb))
        for neigh_neigh in mol.GetAtomWithIdx(neig_idx).GetNeighbors():
          neig_neig_idx = neigh_neigh.GetIdx()
          comb = [neig_neig_idx, t2]
          # 1-4 connected
          if sorted(comb) not in pair_1234:
            pair_1234.append(sorted(comb))
           
    for neigh in mol.GetAtomWithIdx(t2).GetNeighbors():
      neig_idx = neigh.GetIdx()
      if neig_idx != t1:
        comb = [neig_idx, t1]
        # 1-3 connected
        if sorted(comb) not in pair_1234:
          pair_1234.append(sorted(comb))
        for neigh_neigh in mol.GetAtomWithIdx(neig_idx).GetNeighbors():
          neig_neig_idx = neigh_neigh.GetIdx()
          # 1-4 connected
          comb = [neig_neig_idx, t1]
          if sorted(comb) not in pair_1234:
            pair_1234.append(sorted(comb))
  
  # get 1-5 and beyond connected
  pair_15beyond = []
  natoms = len(mol.GetAtoms())
  for n in range(0, natoms-1):
    for m in range(n+1, natoms):
      if [n, m] not in pair_1234:
        pair_15beyond.append([n, m])
  return pair_15beyond

# change to 2.5 so this is consistent with HB definition above
def no_close_contact(conf, atom_pairs, threshold=2.5):
  res = True
  for pair in atom_pairs:
    i, j = pair
    dist = (conf.GetAtomPosition(i)- conf.GetAtomPosition(j)).Length()
    if dist < threshold:
      res = False
      break
  return res

if __name__ == "__main__":
 
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'inputfile', help = "Input file name", required=True)  
  parser.add_argument('-p', dest = 'xtbpath', help = "Absolute path of xtb executable", required=True)  
  parser.add_argument('-o', dest = 'outputfile', help = "Output file name. Default: conftest.mol", default='conftest.mol')  
  parser.add_argument('-f', dest = 'inputformat', help = "Input file format. Default: SDF", choices = ['SDF', 'MOL2'], type=str.upper, default='SDF')  
  
  args = vars(parser.parse_args())
  
  inputfile = args['inputfile']
  outputfile = args['outputfile']
  inputformat = args['inputformat']
  xtbpath = args['xtbpath']

  if inputformat == 'MOL2':
    m1 = Chem.MolFromMol2File(inputfile,removeHs=False)
  else: 
    m1 = Chem.MolFromMolFile(inputfile,removeHs=False)
  
  # check if molecule is 2D
  coords_z = []
  for atom in m1.GetAtoms():
    atom_idx = atom.GetIdx()
    coord = m1.GetConformer().GetAtomPosition(atom_idx) 
    coords_z.append(coord.z)
  if all(abs(z) < 1e-6 for z in coords_z):
    os.system(f"obabel {inputfile} -O {inputfile} --gen3d")
    print(f" Converting {inputfile} to 3D")

  m2 = AllChem.EmbedMultipleConfs(m1, numConfs=500, useExpTorsionAnglePrefs=True,useBasicKnowledge=True, randomSeed=123456789)
  
  # find all ihb
  ihbs = []
  for i in range(m1.GetNumConformers()):
    ihb = find_intramolecular_hbonds(m1, confId=i)
    ihbs.append(ihb)

  ihb_dist ={}
  for ihb in ihbs:
    for ih in ihb:
      idx1, idx2, dist = ih
      u = f"{str(idx1)}-{str(idx2)}"
      if u not in ihb_dist:
        ihb_dist[u] = dist
      else:
        if ihb_dist[u] > dist:
          ihb_dist[u] = dist
        else:
          continue
  
  # match Amide atoms
  pattern = Chem.MolFromSmarts('[C](=[O])[N]([H])[H]')
  matches = m1.GetSubstructMatches(pattern)
  amide_torsions = [] 
  for match in matches:
    c, o, n, h1, h2 = match
    amide_torsions.append([n,c,h1,h2])

  ffps = ChemicalForceFields.MMFFGetMoleculeProperties(m1)
  converged = []
  for i in range(m1.GetNumConformers()):
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m1, ffps, confId=i, ignoreInterfragInteractions=False)
    for k,v in ihb_dist.items():
      idx1, idx2 = k.split('-')
      # atom1, atom2, relative=False, minDist=3.0, maxDist=10.0, forceConstant=50.0 
      ff.MMFFAddDistanceConstraint(int(idx1), int(idx2), False, 3.0, 10.0, 50.0)
    # add torsion restraint for Amide Nitrogen
    if amide_torsions != []:
      for amide_torsion in amide_torsions:
        t1,t2,t3,t4 = amide_torsion
        ff.MMFFAddTorsionConstraint(t1, t2, t3, t4, True, 0, 0, 1000.0)

    res=ff.Minimize(maxIts=500)
    converged.append(res) 
  
  rgs = []
  sasas = []
  num_ihbs = []
  no_steric = []
  pair_15beyond = get_pairAtoms(m1) 
  ptable = Chem.GetPeriodicTable()
  radii = [float(ptable.GetRvdw(atom.GetAtomicNum())) for atom in m1.GetAtoms()]
  for i in range(m1.GetNumConformers()):
    rg = 0.0
    ihb = range(999)
    sasa = 0.0
    conf = m1.GetConformer(i)
    no_steric.append(no_close_contact(conf, pair_15beyond))
    if converged[i] == 0:
      rg = Descriptors3D.RadiusOfGyration(m1, confId=i)
      ihb = find_intramolecular_hbonds(m1, confId=i)
      sasa = Chem.rdFreeSASA.CalcSASA(m1, radii, confIdx=i)
    rgs.append(rg)
    sasas.append(sasa)
    num_ihbs.append(len(ihb))

  ## 1. find the conformers with the fewest HB
  min_num_ihb = min(num_ihbs)
  num_ihbs_deg_idx = []
  for i in range(len(num_ihbs)):
    if num_ihbs[i] == min_num_ihb:
      num_ihbs_deg_idx.append(i)
  
  # 2. keep top 50% of 1. based on Rg
  fifty_per = int(0.5*len(num_ihbs_deg_idx))
  if fifty_per < 1:
    fifty_per = 1
  rgs_sel = [rgs[idx] for idx in num_ihbs_deg_idx] 
  top_fifty_per = list(np.argsort(rgs_sel)[-fifty_per:])

  # 3. filter out conformers with steric interactions
  for idx in top_fifty_per:
    if not no_steric[idx]:
      top_fifty_per.remove(idx)
  
  # 4. use the one with the highest SASA 
  highest_sasa = 0.0
  extended_conf = -1
  for i in top_fifty_per:
    idx = num_ihbs_deg_idx[i]
    sasa = sasas[idx]  
    if sasa > highest_sasa:
      highest_sasa = sasa
      extended_conf = idx
  
  # 5. optimize the geometry with xtb
  extended_conformer = m1.GetConformer(extended_conf)
  currdir = os.getcwd()
  confgendir = os.path.join(currdir, 'ConfGen')
  if not os.path.isdir(confgendir):
    os.mkdir(confgendir)
  conf_fname = 'extended_conformer'  
  with open(os.path.join(confgendir, f'{conf_fname}.xyz'), 'w') as f:
    f.write(f'{m1.GetNumAtoms()}\n\n') 
    for atom_idx in range(m1.GetNumAtoms()):
      atom = m1.GetAtomWithIdx(atom_idx)
      pos = extended_conformer.GetAtomPosition(atom_idx)
      f.write(f'{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n')
    
  # Write distance constraint file used with xtb
  with open(os.path.join(confgendir, 'constraint.inp'), 'w') as f:
    f.write("$constrain\n")
    f.write("  force constant=0.5\n")
    for k,v in ihb_dist.items():
      idx1, idx2 = k.split('-')
      # here we set distance involving IHB
      f.write(f"  distance: {int(idx1)+1}, {int(idx2)+1}, auto\n")
    f.write("$end\n")
  
  # Get total charge and number of unpaired electrons for later use
  total_charge = Chem.GetFormalCharge(m1)
  spin_multiplicity = m1.GetProp('spinMultiplicity') if m1.HasProp('spinMultiplicity') else 1
  unpaired_electrons = spin_multiplicity - 1  
  if not m1.HasProp('spinMultiplicity'):
    unpaired_electrons = 0
  
  # Do xtb optimization for each conformer
  os.chdir(confgendir)
  print('Running XTB optimization for extended conformer...')
  
  # Use implicit solvent for optimization if there are more than one charged atoms
  # We use crude level + 20 cycles so that the extended geometry is almost kept
  charged_atoms = [a for a in m1.GetAtoms() if a.GetFormalCharge() != 0]
  if len(charged_atoms) > 1:
    cmdstr = f'{xtbpath} {conf_fname}.xyz --opt crude --input constraint.inp --chrg {total_charge} --uhf {unpaired_electrons} --alpb water --cycles 20 > xtb_run.log '
  else:
    cmdstr = f'{xtbpath} {conf_fname}.xyz --opt crude --input constraint.inp --chrg {total_charge} --uhf {unpaired_electrons} --cycles 20 > xtb_run.log '
  os.system(cmdstr)
  
  # Replace coordinates with those in optimized xyz file
  os.chdir(currdir)
  lines = open(os.path.join(confgendir, 'xtbopt.xyz')).readlines()[2:]
  opted_coords = []
  for line in lines:
    ss = line.split()
    opted_coords.append([float(ss[1]), float(ss[2]), float(ss[3])])
  
  for atom_idx in range(extended_conformer.GetNumAtoms()):
    extended_conformer.SetAtomPosition(atom_idx, opted_coords[atom_idx])
 
  # Make a copy of the m1 object
  extended_conf_mol = Chem.Mol(m1)  
  extended_conf_mol.RemoveAllConformers()
  extended_conf_mol.AddConformer(extended_conformer, assignId=True)

  Chem.MolToMolFile(extended_conf_mol, outputfile)

# check if molecule is 2D
if inputformat == "SDF":
  coords_z = []
  lines = open(outputfile).readlines()
  for line in lines:
    ss = line.split()
    if len(ss) == 16:
      coords_z.append(float(ss[2]))
  if (all(coords_z) == 0):
    os.system(f"obabel {outputfile} -O {outputfile} --gen3d")
    print(f" Converting {outputfile} to 3D")
