
# This program:
# 1. Grow a fragment containing the atom of interest
# 2. Derive multipole parameters for the fragment 

# Author: Chengwen Liu
# Date: June 2024

import os
import sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem,rdmolfiles,EditableMol,RingInfo,ChemicalForceFields

# helper function to save required fragment
def saveFragment(mol, atomsToKeep, fragname):
  # add attached Hydrogen atoms
  for idx in atomsToKeep:
    a = mol.GetAtomWithIdx(idx)
    for neig in a.GetNeighbors():
      atomNum = neig.GetAtomicNum()
      if atomNum == 1:
        neig_idx = neig.GetIdx()
        atomsToKeep.append(neig_idx)
  
  atomsToRemove = []
  for idx in range(len(mol.GetAtoms())):
    if idx not in atomsToKeep:
      atomsToRemove.append(idx)
  
  emol = Chem.EditableMol(mol)
  atomsToRemove.sort(reverse=True)
  for atom in atomsToRemove:
      emol.RemoveAtom(atom)
  
  mol2 = emol.GetMol()
  Chem.SanitizeMol(mol2)
  mol2h = Chem.AddHs(mol2,addCoords=True)
  Chem.Kekulize(mol2h)
  Chem.SanitizeMol(mol2h)
  # using MMFF to optimize the structure
  ffps = ChemicalForceFields.MMFFGetMoleculeProperties(mol2h)
  ff = ChemicalForceFields.MMFFGetMoleculeForceField(mol2h, ffps, confId=-1, ignoreInterfragInteractions=False)
  ff.Minimize(maxIts=500)
  rdmolfiles.MolToMolFile(mol2h, fragname)
  return

def growFragment(atomidx, sdffile):
  mol = Chem.MolFromMolFile(sdffile,removeHs=False)
  
  # 1. Atom in three rings
  # 2. Atom in two rings
  # 3. Atom in one ring
  # 4. Atom not in ring
  
  atomsInThreeRing = []
  pattern = Chem.MolFromSmarts('[*;R3]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInThreeRing.append(match[0])
  
  atomsInTwoRing = []
  pattern = Chem.MolFromSmarts('[*;R2]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInTwoRing.append(match[0])
  
  atomsInOneRing = []
  pattern = Chem.MolFromSmarts('[*;R1]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInOneRing.append(match[0])
  
  atomsInRing = []
  pattern = Chem.MolFromSmarts('[*;R]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInRing.append(match[0])
  
  atomsInFiveMemRing = []
  pattern = Chem.MolFromSmarts('[*;r5]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInFiveMemRing.append(match[0])
  
  atomsInSixMemRing = []
  pattern = Chem.MolFromSmarts('[*;r6]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInSixMemRing.append(match[0])

  # Case 1: the atom is in three rings
  if (atomidx - 1) in atomsInThreeRing:
    atomsToKeep = []
    t1 = mol.GetAtomWithIdx(atomidx-1)
    for neig in t1.GetNeighbors():
      neig_idx = neig.GetIdx()
      if neig_idx in atomsInSixMemRing: 
        atomsToKeep.append(neig_idx)
        break
   
    atomsToAdd = []
    bonded_atom = atomsToKeep[0]
    for ring_atom in atomsInRing:
      sameRing1 = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), atomidx-1, ring_atom)
      sameRing2 = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), bonded_atom, ring_atom)
      if sameRing1 and sameRing2:
        atomsToAdd.append(ring_atom)

    atomsOfTwoRings = list(set(atomsToAdd + atomsToKeep))
    saveFragment(mol, atomsOfTwoRings, f"Frag_Atom{atomidx:03d}.mol")
  
  # Case 2: the atom is in two rings
  ## in this case it is a conjugated aromatic system
  ## need to break the conjugation and only have two rings
  if (atomidx - 1) in atomsInTwoRing:
    atomsToKeep = []
    t1 = mol.GetAtomWithIdx(atomidx-1)
    for neig in t1.GetNeighbors():
      neig_idx = neig.GetIdx()
      atomsToKeep.append(neig_idx)
    
    atomsToAdd = []
    for i in range(len(atomsToKeep)):
      bonded_atom = atomsToKeep[i]
      for ring_atom in atomsInRing:
        sameRing1 = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), atomidx-1, ring_atom)
        sameRing2 = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), bonded_atom, ring_atom)
        if sameRing1 and sameRing2:
          atomsToAdd.append(ring_atom)
    atomsOfTwoRings = list(set(atomsToAdd + atomsToKeep))
    saveFragment(mol, atomsOfTwoRings, f"Frag_Atom{atomidx:03d}.mol")
  
  ## Case: the atom is in one ring 
  ### in this case it is a conjugated aromatic system
  ### need to break the conjugation and only have the ring
  if (atomidx - 1) in atomsInOneRing:
    atomsToAdd = []
    for ring_atom in atomsInRing:
      sameRing = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), atomidx-1, ring_atom)
      if sameRing: 
        atomsToAdd.append(ring_atom)
    saveFragment(mol, atomsToAdd, f"Frag_Atom{atomidx:03d}.mol")
  
  return 

# need a helper function to determine where to cut
# need another function to match the 6+5 special case to cut the bond

if __name__ == '__main__':
  
  global sdffile
  sdffile = sys.argv[1]
  atom_id = int(sys.argv[2]) 
  growFragment(atom_id, sdffile)
