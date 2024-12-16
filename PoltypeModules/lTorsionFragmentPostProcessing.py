
# This program is used by fragmenter.py 
# - takes the {fname}.sdf file of poltype torsion fragment job
# - write a new file {fname}-post.mol to fix possible broken aromaticity
# - further trimming on the fragments based on several predefined chemical rules
# - generates a {fname}.mol file for using with afterwards torsion fitting

# KEY RULES IMPLEMENTED IN THIS PROGRAM:
# 1. Bond breaking only happens at single bonds
# 2. Special bond involving nitrogen is not cut
# 3. Fused rings are treated as single rings if the fused two atoms are carbon
# 4. Long alkane chain such as propyl and ethyl are replaced with CH3-
# 5. Charged fragments (COO-) is neutralized
# 6. Substituents on a ring (OLD decision)
#   -- 1,3- and 1,4- substituents are removed (meta- and para- )
#   -- 1,2- substituents are kept if they are polar groups
#   -- 1,2- substituents are removed if they are alkane
# 6. Substituents on a ring are removed to avoid potential HB or steric effect


# Author: Chengwen Liu
# Date: Jun. 2023

import os
import sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem,rdmolfiles,EditableMol,RingInfo,Descriptors

# helper function to save required fragment
def saveFragment(mol, atomsToKeep,output, xyz2index):

  atomsToRemove = []
  for idx in range(len(mol.GetAtoms())):
    if idx not in atomsToKeep:
      atomsToRemove.append(idx)
  
  emol = Chem.EditableMol(mol)
  atomsToRemove.sort(reverse=True)
  for atom in atomsToRemove:
    emol.RemoveAtom(atom)
  
  mol1 = emol.GetMol()
  rdmolfiles.MolToMolFile(mol1, 'no_h.mol')
  
  # so far, mol2 only has the correct coordinate of required atoms
  # next we need put the correct bond type for each bond

  mol2 = Chem.MolFromMolFile('no_h.mol',removeHs=False)
 
  for bond in mol2.GetBonds():
    idx1 = bond.GetBeginAtomIdx()
    idx2 = bond.GetEndAtomIdx()
    
    coord = mol2.GetConformer().GetAtomPosition(idx1) 
    x, y, z = f"{coord.x:10.4f}", f"{coord.y:10.4f}", f"{coord.z:10.4f}"
    xyzstr = ''.join([x,y,z])
    idx1_in_par = par_xyz2index[xyzstr]
    
    coord = mol2.GetConformer().GetAtomPosition(idx2) 
    x, y, z = f"{coord.x:10.4f}", f"{coord.y:10.4f}", f"{coord.z:10.4f}"
    xyzstr = ''.join([x,y,z])
    idx2_in_par = par_xyz2index[xyzstr]
    
    if mol.GetAtomWithIdx(idx1_in_par).GetIsAromatic() and mol.GetAtomWithIdx(idx2_in_par).GetIsAromatic():
      bond.SetBondType(Chem.BondType.AROMATIC)
      print(f'Setting bond type between {idx1} and {idx2} to AROMATIC')

  Chem.SanitizeMol(mol2)
  mol2h = Chem.AddHs(mol2,addCoords=True)
  rdmolfiles.MolToMolFile(mol2h, output)
  return

if __name__ == '__main__':
  sdffile = sys.argv[1]
  parent_mol_file = sys.argv[2]
 
  fname = sdffile.split('.sdf')[0]
  fname = fname.split('.mol')[0]
  
  mol = Chem.MolFromMolFile(sdffile,removeHs=False)
  
  # check if molecule is 2D
  coords_z = []
  for atom in mol.GetAtoms():
    atom_idx = atom.GetIdx()
    coord = mol.GetConformer().GetAtomPosition(atom_idx) 
    coords_z.append(coord.z)
  if all(abs(z) < 1e-6 for z in coords_z):
    os.system(f"obabel {sdffile} -O {sdffile} --gen3d")
    print(f" Converting {sdffile} to 3D")
  
  # here we try to match the coorinates of fragment and parent molecules
  # if they match (for heavy atoms)
  # re-generate the fragment from the parent molecule
  # this is to fix the broken aromaticity poltype has 
  # Chengwen Liu
  # May 2024

  par_mol = Chem.MolFromMolFile(parent_mol_file, removeHs=False)
  
  atom_indices = []
  par_xyz2index = {}
  for atom in par_mol.GetAtoms():
    atom_idx = atom.GetIdx()
    coord = par_mol.GetConformer().GetAtomPosition(atom_idx) 
    x, y, z = f"{coord.x:10.4f}", f"{coord.y:10.4f}", f"{coord.z:10.4f}"
    xyzstr = ''.join([x,y,z])
    par_xyz2index[xyzstr] = atom_idx
  
  n_heavy_atom_match = 0
  rotbond_coords = [] 
  nitrogenMissingHs = {} 
  for i in range(len(mol.GetAtoms())):
    atom = mol.GetAtoms()[i]
    atom_idx = atom.GetIdx()
    atomNum = atom.GetAtomicNum()
    
    coord = mol.GetConformer().GetAtomPosition(atom_idx) 
    x, y, z = f"{coord.x:10.4f}", f"{coord.y:10.4f}", f"{coord.z:10.4f}"
    xyzstr = ''.join([x,y,z])
    if atom_idx == 1:
      rotbond_coords.append(xyzstr)
    if atom_idx == 2:
      rotbond_coords.append(xyzstr)

    if (xyzstr in par_xyz2index.keys()):
      atom_indices.append(par_xyz2index[xyzstr])
      if (atomNum != 1):
        n_heavy_atom_match += 1
    
    if (atom.GetAtomicNum() == 7):
      frag_neighbors = [x.GetSymbol() for x in atom.GetNeighbors()]
      par_atom_idx = par_xyz2index[xyzstr]
      par_neighbors = [x.GetSymbol() for x in par_mol.GetAtomWithIdx(par_atom_idx).GetNeighbors()]

      # Chengwen Liu
      # Oct 2024
      # check the case when parent/frag mol has N+ but different hydrogens
      if frag_neighbors.count('H') >= 2 and (atom.GetFormalCharge() == 1):
        # here we set to 2 because tetra-valent nitrogen is not allowed in RDKit
        # since the formal charge is 1, (as set above), so in the end it will 
        # have the correct number of hydrogen atoms
        nitrogenMissingHs[par_atom_idx] = 2
      # check the case when parent/frag mol has aromatic N but different hydrogens
      if (frag_neighbors.count('H') == 1) and (par_neighbors.count('H') == 0):
        nitrogenMissingHs[par_atom_idx] = 1

  # Set the correct number of hydrogens 
  for idx,numOfH in nitrogenMissingHs.items():
    atom = par_mol.GetAtomWithIdx(idx)
    atom.SetNumExplicitHs(numOfH)
  
  n_heavy_atom_total = Descriptors.HeavyAtomCount(mol)
  if n_heavy_atom_total == n_heavy_atom_match:
    try:
      saveFragment(par_mol, atom_indices, f"{fname}_tmp.mol", par_xyz2index)
    except:
      pass

  # use the {fname}_post.mol as the input file
  if os.path.isfile(f"{fname}_tmp.mol"):
    # re-order the atoms so that 1-2 is the rotbond
    tmp_mol = Chem.MolFromMolFile(f"{fname}_tmp.mol", removeHs=False)
    for atom in tmp_mol.GetAtoms():
      atom_idx = atom.GetIdx()
      coord = tmp_mol.GetConformer().GetAtomPosition(atom_idx) 
      x, y, z = f"{coord.x:10.4f}", f"{coord.y:10.4f}", f"{coord.z:10.4f}"
      xyzstr = ''.join([x,y,z])
      if xyzstr == rotbond_coords[0]:
        first_idx = atom_idx
      if xyzstr == rotbond_coords[1]:
        second_idx = atom_idx
        
    new_order = list(range(len(tmp_mol.GetAtoms())))
    new_order.remove(first_idx)
    new_order.remove(second_idx)
    new_order.insert(1, first_idx)
    new_order.insert(2, second_idx)
    tmp_mol = Chem.RenumberAtoms(tmp_mol, new_order) 
    rdmolfiles.MolToMolFile(tmp_mol, f'{fname}_post.mol')
    sdffile = f"{fname}_post.mol" 
    mol = Chem.MolFromMolFile(sdffile,removeHs=False)

  # we have the "correct" fragment to work on
  # below we focus on trimming the fragment

  # all single bonds are eligible to cut,
  single_bonds = []
  pattern = Chem.MolFromSmarts('[*;!#1][*;!#1]')
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
  for match in matches:
    single_bonds.append(list(match))
  
  # except non-trivalence nitrogen 
  special_nitrogen = []
  pattern = Chem.MolFromSmarts('[#7X2][*]')
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
  for match in matches:
    special_nitrogen.append([match[0], match[1]])
    special_nitrogen.append([match[1], match[0]])
  
  
  # record the aromatic nitrogen (n-*),
  aromatic_nitrogen = {}
  pattern = Chem.MolFromSmarts('[n][*]')
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
  for match in matches:
    n, x = match
    if n not in aromatic_nitrogen.keys():
      aromatic_nitrogen[n] = [x]
    else:
      aromatic_nitrogen[n] += [x]
 
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
  
  # OLD decision 
  # Move the polar group from ortho to para position
  # if there are two polar groups on two ortho- site,
  # only the first-appearing one (smaller index) will
  # be moved, and the second one will be deleted
  # Chengwen Liu
  # Oct 2024
  
  #atomsInSixMemRing = []
  #pattern = Chem.MolFromSmarts('[a;r6]')
  #matches = mol.GetSubstructMatches(pattern)
  #for match in matches:
  #  if len(match) == 1:
  #    atomsInSixMemRing.append(match[0])
  #
  #ortho_atom = []
  #paraH = []
  #pattern = Chem.MolFromSmarts('[*;!#1][R][R][R][R][H]')
  #matches = mol.GetSubstructMatches(pattern, uniquify=False)
  #for match in matches:
  #  mat = [match[0], match[1]] 
  #  if mat == [1,2]: # poltype always uses 1,2
  #    if (match[2] not in ortho_atom) and (match[2] in atomsInSixMemRing):
  #      ortho_atom.append(match[2])
  #      paraH.append(match[5])
  #  if mat == [2,1]: # poltype always uses 1,2
  #    if (match[2] not in ortho_atom) and (match[2] in atomsInSixMemRing):
  #      ortho_atom.append(match[2])
  #      paraH.append(match[5])
  #
  #ortho2paraH = dict(zip(ortho_atom, paraH))
  #
  ## match -OH,-NH2,-SH,-F,-Cl,-Br,-I
  #pattern = Chem.MolFromSmarts('[a][O,N,S,F,Cl,Br,I]')
  #matches = mol.GetSubstructMatches(pattern, uniquify=False)
  #mol_ed = Chem.EditableMol(mol)
  #tmp_add = []
  #for match in matches:
  #  if match[0] in ortho_atom:
  #    paraH_idx = ortho2paraH[match[0]]
  #    if paraH_idx not in tmp_add:
  #      atomicNum = mol.GetAtomWithIdx(match[1]).GetAtomicNum()
  #      add_atom = Chem.Atom(atomicNum)
  #      mol_ed.ReplaceAtom(paraH_idx, add_atom)
  #      tmp_add.append(paraH_idx)
  #mol = mol_ed.GetMol()
  #Chem.SanitizeMol(mol)

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
          if not set(pair).issubset(set(torsion_idx_list)) and (pair in single_bonds or pair_r in single_bonds) and (pair not in special_nitrogen): 
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
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
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
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
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
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
  for match in matches:
    x, n, h1, h2, h3 = match
    if n not in torsion_idx_list:
      atomsToRemove += [n, h1, h2, h3] 
  
  # Replace NH2 with H 
  pattern = Chem.MolFromSmarts('[*][N]([H])[H]')
  matches = mol.GetSubstructMatches(pattern, uniquify=False)
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
  
  # final safety check
  for idx in torsion_idx_list:
    if idx in atomsToRemove:
      atomsToRemove = []
      break
  
  # if aromatic nitrogen misses connections
  # after cut, add one hydrogen atom
  # this should be done in mol (before cut)
  # since the atom order will be different after cut
  # Chengwen Liu
  # Oct 2024

  for idx, connected_atoms in aromatic_nitrogen.items():
    for atm in connected_atoms:
      if (atm in atomsToRemove) and (idx not in atomsToRemove):
        nitrogen = mol.GetAtomWithIdx(idx)
        print(f'Setting number of attached hydrogen to be 1 on atom {idx}')
        nitrogen.SetNumExplicitHs(1)
  
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
  matches_1 = mol2.GetSubstructMatches(pattern, uniquify=False)
  # neutralize NH- to NH2
  pattern = Chem.MolFromSmarts('[N-][H]')
  matches_2 = mol2.GetSubstructMatches(pattern, uniquify=False)
  updated_mol = mol2
  if len(matches_1) != 0:
    print('COO- group detected. Neutralizing!')
    for match in matches_1:
      idx = match[0]
      oxygen = mol2.GetAtomWithIdx(idx)
      oxygen.SetFormalCharge(0)
      oxygen.SetNumExplicitHs(1)
      oxygen.UpdatePropertyCache()
      updated_mol = Chem.AddHs(updated_mol, addCoords=True)
    rdmolfiles.MolToMolFile(updated_mol, f'{fname}.mol')
  elif len(matches_2) != 0:
    print('NH- group detected. Neutralizing!')
    for match in matches_2:
      idx = match[0]
      notrogen = mol2.GetAtomWithIdx(idx)
      nitrogen.SetFormalCharge(0)
      oxygen.SetNumExplicitHs(2)
      nitrogen.UpdatePropertyCache()
      updated_mol = Chem.AddHs(updated_mol, addCoords=True)
    rdmolfiles.MolToMolFile(updated_mol, f'{fname}.mol')
  else:
    rdmolfiles.MolToMolFile(mol2h, f'{fname}.mol')
    
  # re-order the atoms so that 1-2 is the rotbond
  # which is required by poltype torsion fragment job
  final_mol = Chem.MolFromMolFile(f'{fname}.mol',removeHs=False)
  for atom in final_mol.GetAtoms():
    atom_idx = atom.GetIdx()
    coord = final_mol.GetConformer().GetAtomPosition(atom_idx) 
    x, y, z = f"{coord.x:10.4f}", f"{coord.y:10.4f}", f"{coord.z:10.4f}"
    xyzstr = ''.join([x,y,z])
    if xyzstr == rotbond_coords[0]:
      first_idx = atom_idx
    if xyzstr == rotbond_coords[1]:
      second_idx = atom_idx
  new_order = list(range(len(final_mol.GetAtoms())))
  new_order.remove(first_idx)
  new_order.remove(second_idx)
  new_order.insert(1, first_idx)
  new_order.insert(2, second_idx)
  final_mol = Chem.RenumberAtoms(final_mol, new_order) 
  rdmolfiles.MolToMolFile(final_mol, f'{fname}.mol')