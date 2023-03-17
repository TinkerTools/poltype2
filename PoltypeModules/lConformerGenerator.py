
# This program 
# - takes the .sdf input structure
# - uses RDKit to generate 500 conformers
# - optimizes the conformers with distance constraints
# - selects the extended conformer with Rg, SASA and number of Intramolecular HB criteria
# - generates "testconf.mol" which poltype wants

# Author: Chengwen Liu
# Date: Mar. 2023

import sys
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


if __name__ == "__main__":
  
  sdffile = sys.argv[1]

  m1 = Chem.MolFromMolFile(sdffile,removeHs=False)
  #AllChem.EmbedMolecule(m1)
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

  ffps = ChemicalForceFields.MMFFGetMoleculeProperties(m1)
  converged = []
  for i in range(m1.GetNumConformers()):
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m1, ffps, confId=i, ignoreInterfragInteractions=False)
    for k,v in ihb_dist.items():
      idx1, idx2 = k.split('-')
      # atom1, atom2, relative=False, minDist=3.0, maxDist=10.0, forceConstant=50.0 
      ff.MMFFAddDistanceConstraint(int(idx1), int(idx2), False, 3.0, 10.0, 50.0)
    res=ff.Minimize(maxIts=500)
    converged.append(res) 

  rgs = []
  sasas = []
  num_ihbs = []

  radii = rdFreeSASA.classifyAtoms(m1)
  for i in range(m1.GetNumConformers()):
    rg = 0.0
    ihb = range(999)
    sasa = 0.0
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
  
  # 3. use the one with the highest SASA 
  highest_sasa = 0.0
  extended_conf = -1
  for i in top_fifty_per:
    idx = num_ihbs_deg_idx[i]
    sasa = sasas[idx]  
    if sasa > highest_sasa:
      highest_sasa = sasa
      extended_conf = idx
  rdmolfiles.MolToMolFile(m1, 'conftest.mol', confId=extended_conf) 