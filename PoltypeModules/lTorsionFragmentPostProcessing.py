
# This program is used by fragmenter.py 
# - takes the {fname}.sdf file of poltype torsion fragment job 
# - further trimming on the fragments based on several predefined chemical rules
# - generates a {fname}.mol file for using with afterwards torsion fitting

# NOTE: 
# 1. the precessed file (.mol) can be the same as the original file (.sdf)

# KEY RULES IMPLEMENTED IN THIS PROGRAM:
# 1. Bond breaking only happens at single bonds
# 2. Special single bond (aromatic nitrogen) is not cut
# 3. Fused rings are treated as single rings if the fused two atoms are carbon
# 4. Long alkane chain such as propyl and ethyl are replaced with CH3-
# 5. Charged fragments (COO-) is neutralized
# 6. Substituents on a ring
#   -- 1,3- and 1,4- substituents are removed (meta- and para- )
#   -- 1,2- substituents are kept if they are polar groups
#   -- 1,2- substituents are removed if they are alkane


# Author: Chengwen Liu
# Date: Jun. 2023

import sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem,rdmolfiles,EditableMol,RingInfo


sdffile = sys.argv[1]
fname = sdffile.split('.sdf')[0]

mol = Chem.MolFromMolFile(sdffile,removeHs=False)

# all single bonds are eligible to cut,
single_bonds = []
pattern = Chem.MolFromSmarts('[*;!#1][*;!#1]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  single_bonds.append(list(match))

# except aromatic nitrogen (n-*)
aromatic_nitrogen = []
pattern = Chem.MolFromSmarts('[n][*]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  aromatic_nitrogen.append([match[0], match[1]])
  aromatic_nitrogen.append([match[1], match[0]])

# do not cut two groups on ortho- sites (1-4)
# unless the connected atom is carbon

pattern = Chem.MolFromSmarts('[*;!R;!#1][R][R][*;!R;!#1]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  mat1 = [match[0], match[1]] 
  mat1_r = [match[1], match[0]] 
  mat2 = [match[2], match[3]] 
  mat2_r = [match[3], match[2]] 
  if set(mat1) == set([1,2]): # poltype always uses 1,2
    if mat2 in single_bonds:
      single_bonds.remove(mat2)
    if mat2_r in single_bonds:
      single_bonds.remove(mat2_r)
  if set(mat2) == set([1,2]): # poltype always uses 1,2
    if mat1 in single_bonds:
      single_bonds.remove(mat1)
    if mat1_r in single_bonds:
      single_bonds.remove(mat1_r)

# get the atom indices of torsion
torsion_idx_list = []
for idx in range(len(mol.GetAtoms())):
  idx1, idx2 = [1, 2] # poltype always uses 1,2
  t1 = mol.GetAtomWithIdx(idx1)
  t2 = mol.GetAtomWithIdx(idx2)
  for neig in t1.GetNeighbors():
    neig_idx = neig.GetIdx()
    if neig_idx not in torsion_idx_list:
      torsion_idx_list.append(neig_idx)
  for neig in t2.GetNeighbors():
    neig_idx = neig.GetIdx()
    if neig_idx not in torsion_idx_list:
      torsion_idx_list.append(neig_idx)

# construct a graph based on connectivity 
g = nx.Graph()
nodes = []
edges = []
for bond in mol.GetBonds():
  idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
  s = sorted([idx1, idx2])
  if idx1 not in nodes:
    nodes.append(idx1)
  if idx2 not in nodes:
    nodes.append(idx2)
  if s not in edges:
    edges.append(s)

g.add_nodes_from(nodes)
g.add_edges_from(edges)

# detect substituents 
sub_bonds = []
for idx in range(len(mol.GetAtoms())):
  a = mol.GetAtomWithIdx(idx)
  inRing = a.IsInRing()
  if inRing:
    for neig in a.GetNeighbors():
      neig_idx = neig.GetIdx()
      sameRing = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), idx, neig_idx)
      atomNum = neig.GetAtomicNum()
      if (not sameRing) and (atomNum != 1):
        pair = [idx, neig_idx]
        pair_r = [neig_idx, idx]
        if not set(pair).issubset(set(torsion_idx_list)) and (pair in single_bonds or pair_r in single_bonds) and (pair not in aromatic_nitrogen): 
          spair = sorted([idx, neig_idx])
          if spair not in sub_bonds:
            sub_bonds.append(spair)

for sub_bond in sub_bonds:
  a1, a2 = sub_bond
  g.remove_edge(a1, a2)  

# clean up the isolated fragments
atomsToRemove = []
frags = nx.connected_components(g)
for m in frags:
  if not set(torsion_idx_list).issubset(list(m)):
    atomsToRemove += list(m)

# Replace CH3 with H 
pattern = Chem.MolFromSmarts('[*][CH3]([H])([H])[H]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  x, c, h1, h2, h3 = match
  if c not in torsion_idx_list:
    atomsToRemove += [c, h1, h2, h3] 
    xatom = mol.GetAtomWithIdx(x)
    atomic = xatom.GetAtomicNum()
    if atomic == 7:
      xatom.SetNumExplicitHs(xatom.GetTotalNumHs() + 1)


# Replace CH2CH3 with H 
pattern = Chem.MolFromSmarts('[*][CH2]([H])([H])[CH3]([H])([H])[H]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  x, c1, h11, h12, c2, h21, h22, h23 = match
  if c1 not in torsion_idx_list:
    atomsToRemove += [c1, h11, h12, c2, h21, h22, h23] 
    xatom = mol.GetAtomWithIdx(x)
    atomic = xatom.GetAtomicNum()
    if atomic == 7:
      xatom.SetNumExplicitHs(xatom.GetTotalNumHs() + 1)

# Replace NH3+ with H 
pattern = Chem.MolFromSmarts('[*][N+]([H])([H])[H]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  x, n, h1, h2, h3 = match
  if n not in torsion_idx_list:
    atomsToRemove += [n, h1, h2, h3] 

# Replace NH2 with H 
pattern = Chem.MolFromSmarts('[*][N]([H])[H]')
matches = mol.GetSubstructMatches(pattern)
for match in matches:
  x, n, h1, h2 = match
  if n not in torsion_idx_list:
    atomsToRemove += [n, h1, h2] 

# Remove atoms on distant fused rings
ring_info = mol.GetRingInfo()
atomsInRing = []
atomsInFusedRing = []
for idx in range(len(mol.GetAtoms())):
  n = RingInfo.NumAtomRings(ring_info, idx)
  if n != 0:
    atomsInRing.append(idx)
  if n > 1:
    atomsInFusedRing.append(idx)

# only keep C-C bond in fused ring as candidate to cut 
bondsInFusedRing = []
for bond in mol.GetBonds():
  idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
  if set([idx1, idx2]).issubset(set(atomsInFusedRing)):
    atomic1 = mol.GetAtomWithIdx(idx1).GetAtomicNum()
    atomic2 = mol.GetAtomWithIdx(idx2).GetAtomicNum()
    if (atomic1 == 6) and (atomic2 == 6):
      bondsInFusedRing.append([idx1, idx2])
    
allAtomsInFusedRing = []
for idx in atomsInRing:
  for bond in bondsInFusedRing:
    for fusedIdx in bond: 
      sameRing = RingInfo.AreAtomsInSameRing(ring_info, idx, fusedIdx) 
      if sameRing:
        if idx not in allAtomsInFusedRing:
          allAtomsInFusedRing.append(idx)

for idx in allAtomsInFusedRing:
  for prob_idx in torsion_idx_list:
    if prob_idx in allAtomsInFusedRing:
      sameRing = RingInfo.AreAtomsInSameRing(ring_info, idx, prob_idx) 
      if not sameRing:
        atomsToRemove.append(idx)
        for neigh in mol.GetAtomWithIdx(idx).GetNeighbors():
          if neigh.GetAtomicNum() == 1:
            neig_h = neigh.GetIdx()
            atomsToRemove.append(neig_h) 

atomsToRemove = list(set(atomsToRemove))

# clean up the isolated atoms
for atom in atomsToRemove:
  if atom in g.nodes:
    g.remove_node(atom)

frags = nx.connected_components(g)
for m in frags:
  if not set(torsion_idx_list).issubset(list(m)):
    atomsToRemove += list(m)

emol = Chem.EditableMol(mol)
atomsToRemove.sort(reverse=True)
for atom in atomsToRemove:
    emol.RemoveAtom(atom)

mol2 = emol.GetMol()
Chem.SanitizeMol(mol2)
mol2h = Chem.AddHs(mol2,addCoords=True)
Chem.Kekulize(mol2h)
Chem.SanitizeMol(mol2h)

# neutralize COO- to COOH
pattern = Chem.MolFromSmarts('[O-][C](=O)')
matches = mol2.GetSubstructMatches(pattern)
updated_mol = mol2
if len(matches) != 0:
  for match in matches:
    idx = match[0]
    mol2.GetAtomWithIdx(idx).SetFormalCharge(0)
    Chem.SanitizeMol(mol2)
    updated_mol = Chem.AddHs(updated_mol, addCoords=True)
  rdmolfiles.MolToMolFile(updated_mol, f'{fname}.mol')
else:
  rdmolfiles.MolToMolFile(mol2h, f'{fname}.mol')