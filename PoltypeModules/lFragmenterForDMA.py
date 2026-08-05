
# This program:
# 1. Grow a fragment containing the atom of interest
# 2. Derive multipole parameters for the fragment 

# Author: Chengwen Liu
# Date: June 2024

import os
import sys
import shutil
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem,rdmolfiles,EditableMol,RingInfo,ChemicalForceFields

# helper function to save required fragment
def saveFragment(mol, atomsToKeep, fragname, bondsToDelete=[]):
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
  
  if bondsToDelete != []:
    for i in range(0, len(bondsToDelete), 2):
      idx1, idx2 = bondsToDelete[i], bondsToDelete[i+1]
      bond_index = mol.GetBondBetweenAtoms(idx1, idx2).GetIdx()
      emol.RemoveBond(idx1, idx2)
  
  atomsToRemove.sort(reverse=True)
  for atom in atomsToRemove:
      emol.RemoveAtom(atom)
  
  mol2 = emol.GetMol()
  Chem.SanitizeMol(mol2)
  mol2h = Chem.AddHs(mol2,addCoords=True)
  Chem.Kekulize(mol2h)
  Chem.SanitizeMol(mol2h)
  # using MMFF to optimize the structure
  #ffps = ChemicalForceFields.MMFFGetMoleculeProperties(mol2h)
  #ff = ChemicalForceFields.MMFFGetMoleculeForceField(mol2h, ffps, confId=-1, ignoreInterfragInteractions=False)
  #ff.Minimize(maxIts=500)
  rdmolfiles.MolToMolFile(mol2h, fragname)
  return

def growFragment(atomidx, sdffile):
  mol = Chem.MolFromMolFile(sdffile,removeHs=False)
  print(mol) 
  # 1. Atom in three rings
  # 2. Atom in two rings
  # 3. Atom in one ring
  # 4. Atom not in ring: not considered here
  
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

  atomsInTripleBond = []
  pattern = Chem.MolFromSmarts('[$(*#*)]')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      atomsInTripleBond.append(match[0])

  print('Atoms in triple bond',atomsInTripleBond)


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
    connectedAtoms = findConnectedAtoms(mol, atomsOfTwoRings, atomsInRing)
    try:
      saveFragment(mol, atomsOfTwoRings, f"Frag_Atom{atomidx:03d}_0.mol")
    except:
      pass

    atomsOfTwoRings += connectedAtoms
    saveFragment(mol, atomsOfTwoRings, f"Frag_Atom{atomidx:03d}.mol")
    
    if os.path.isfile(f"Frag_Atom{atomidx:03d}_0.mol"):
      testmol = Chem.MolFromMolFile(f"Frag_Atom{atomidx:03d}.mol",removeHs=False)
      if len(testmol.GetAtoms()) == len(mol.GetAtoms()):
        shutil.move(f"Frag_Atom{atomidx:03d}_0.mol", f"Frag_Atom{atomidx:03d}.mol")
      else:
        os.system(f"rm -f Frag_Atom{atomidx:03d}_0.mol")
    
    ## special case detection
    specialCaseSixFusedFiveMemRing(f"Frag_Atom{atomidx:03d}.mol")
    
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
    connectedAtoms = findConnectedAtoms(mol, atomsOfTwoRings, atomsInRing)
    
    # rear cases that the kekulization fail 
    try:
      saveFragment(mol, atomsOfTwoRings, f"Frag_Atom{atomidx:03d}_0.mol")
    except:
      pass

    atomsOfTwoRings += connectedAtoms
    saveFragment(mol, atomsOfTwoRings, f"Frag_Atom{atomidx:03d}.mol")
    
    if os.path.isfile(f"Frag_Atom{atomidx:03d}_0.mol"):
      testmol = Chem.MolFromMolFile(f"Frag_Atom{atomidx:03d}.mol",removeHs=False)
      if len(testmol.GetAtoms()) == len(mol.GetAtoms()):
        shutil.move(f"Frag_Atom{atomidx:03d}_0.mol", f"Frag_Atom{atomidx:03d}.mol")
      else:
        os.system(f"rm -f Frag_Atom{atomidx:03d}_0.mol")
    
    ## special case detection
    specialCaseSixFusedFiveMemRing(f"Frag_Atom{atomidx:03d}.mol")
      
  ## Case: the atom is in one ring 
  ### in this case it is a conjugated aromatic system
  ### need to break the conjugation and only have the ring
  if (atomidx - 1) in atomsInOneRing:
    atomsToAdd = []
    for ring_atom in atomsInRing:
      sameRing = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), atomidx-1, ring_atom)
      if sameRing: 
        atomsToAdd.append(ring_atom)
    
    connectedAtoms = findConnectedAtoms(mol, atomsToAdd, atomsInRing)
    atomsToAdd += connectedAtoms
    try:
      saveFragment(mol, atomsToAdd, f"Frag_Atom{atomidx:03d}_0.mol")
    except:
      pass

    atomsToAdd += connectedAtoms
    saveFragment(mol, atomsToAdd, f"Frag_Atom{atomidx:03d}.mol")
    
    if os.path.isfile(f"Frag_Atom{atomidx:03d}_0.mol"):
      testmol = Chem.MolFromMolFile(f"Frag_Atom{atomidx:03d}.mol",removeHs=False)
      if len(testmol.GetAtoms()) == len(mol.GetAtoms()):
        shutil.move(f"Frag_Atom{atomidx:03d}_0.mol", f"Frag_Atom{atomidx:03d}.mol")
      else:
        os.system(f"rm -f Frag_Atom{atomidx:03d}_0.mol")

  # Case when atom is part of triple bond
  if (atomidx -1) in atomsInTripleBond:
    atomsToAdd = []
    atomsToKeep = [atomidx-1]
    t1 = mol.GetAtomWithIdx(atomidx-1)
    for neig in t1.GetNeighbors():
      neig_idx = neig.GetIdx()
      atomsToAdd.append(neig_idx)
      atomsToKeep.append(neig_idx)
    print('Atoms to add ',atomsToAdd)
    print('Atoms in ring',atomsInRing) 
    # Loop through the atoms to add and check:
    # if the atom is in an aromatic ring
    # if not, then try to find coonected atoms with specific rules
    for neig in atomsToAdd:

      if neig in atomsInRing:
        print('Found atom in ring as first neighboor')
        connectedAtoms = findConnectedAtoms(mol, atomsToAdd, atomsInRing)
        atomsToKeep += connectedAtoms

      elif neig in atomsInTripleBond:
        print(f'Neig {neig} is in atomsinTripleBond')
        atomsToKeep += [neig]
      else:
        print(f'Neig {neig} is neither in atomsinTripleBOnd nor in atomsinRing')
        connectedAtoms = findConnectedAtomsGeneral(mol, atomsToAdd, atomidx, neig)
        atomsToKeep += connectedAtoms

    print('Final atoms to keep',atomsToKeep)
    saveFragment(mol, atomsToKeep, f"Frag_Atom{atomidx:03d}.mol")
    if os.path.isfile(f"Frag_Atom{atomidx:03d}_0.mol"):
      testmol = Chem.MolFromMolFile(f"Frag_Atom{atomidx:03d}.mol",removeHs=False)
      if len(testmol.GetAtoms()) == len(mol.GetAtoms()):
        shutil.move(f"Frag_Atom{atomidx:03d}_0.mol", f"Frag_Atom{atomidx:03d}.mol")
      else:
        os.system(f"rm -f Frag_Atom{atomidx:03d}_0.mol")

  return 


def find_patern(mol,atm_to_check,pat_to_check):

  """
  This function find the atoms present in a specific pattern
  and check if a given atom is part of the atom present in the 
  pattern.

  Inputs:
    -   mol: rdkit mol object
    -   atm_to_check: rdkit atom id to check
    -   pat_to_check: SMART string to check as pattern

  Outputs:
    -   is_present: if atom_to_check exist or not in the SMART pattern (bool)
    -   Search_pattern: list of rdkit atom ID present in the pattern
  """


  Search_pattern = []
  pattern = Chem.MolFromSmarts(pat_to_check)
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    if len(match) == 1:
      Search_pattern.append(match[0])

  if atm_to_check in Search_pattern:
    is_present = True
  else:
    is_present = False
  
  return is_present, Search_pattern


def findConnectedAtomsGeneral(mol, atomsToadd, atomidx, neig): 

  """
  This function attempt to find atoms connected to triple bond.
  By default, it will cut just after an atom that does not have triple bond:
  if the atom is part of triple bond, then stop after max two iterations.

  TODO:
    - If atoms neig as more than 1 neighboor, then need to add more options

  Inputs:
    -   mol: RDKIT mol object
    -   atomsToadd: list of atom neighboor to main atom to check
    -   atomidx: main atom to check ID
    -   neig: neighboor atom of the main atom to check

  Outputs:
    -   Co_atoms: list of connected atoms
  """

  Break = False
  prev_atm_idx = atomidx - 1
  current_atm_idx = neig
  Co_atoms = []
  cc = 0

  # While loop until the check for connected atoms stops
  while not Break:

    # Grab the neighboor atoms of neig
    Neig = get_neigh_atoms(mol, current_atm_idx, prev_atm_idx)
    print(f'Neig of {current_atm_idx}: {Neig}')
    if len(Neig) == 1 and list(Neig.keys())[0] not in Co_atoms:
        prev_atm_idx = current_atm_idx
        current_atm_idx = list(Neig.keys())[0]

        #if find_patern(mol,current_atm_idx,'[*;R]')[0]:
        #  atomsToAdd = list(Neig.keys())
        #  connectedAtoms = findConnectedAtoms(mol, atomsToAdd, find_patern(mol,current_atm_idx,'[*;R]')[1])
        #  atomsToAdd += connectedAtoms
        #  Co_atoms += atomsToAdd
        #  Break = True
        #else:
        Co_atoms += list(Neig.keys())
        Break = False

    else:
        Co_atoms += list(Neig.keys())
        Break = True

    cc += 1
    if cc == 1:
        Break = True

  return Co_atoms


def get_neigh_atoms(mol, Current_Idx, Prev_Idx):

  """
  This function create a dictionnary such as:
  {Atom ID: Atomic Num}
  Make sure that none of previous atom are selected

  Inputs:
    -   mol: RDKIT mol object
    -   Current_Idx: Atom id to check
    -   Prev_Idx: Atom id of the previous atom (to skip)

  Outputs:
    -   Neighboor: dictionnary with the neighboored atoms 
  
  """
  C_atm = mol.GetAtomWithIdx(Current_Idx)

  Neighboor = {}
  for neib in C_atm.GetNeighbors():
    if neib.GetIdx() != Prev_Idx:
      IDx = mol.GetAtomWithIdx(neib.GetIdx())
      Neighboor[neib.GetIdx()] = IDx.GetAtomicNum()

  return Neighboor


# helper function to determine where to cut
def findConnectedAtoms(rdkitmol, atomList, atomsInRing):
  # C-C single bond 
  eligibleBondsToCut = []
  connectedAtoms = []
  pattern = Chem.MolFromSmarts('[#6]-[#6]')
  matches = rdkitmol.GetSubstructMatches(pattern)
  for match in matches:
    eligibleBondsToCut.append(list(match))
  
  for idx in atomList:
    a = rdkitmol.GetAtomWithIdx(idx)
    for neig in a.GetNeighbors():
      atomNum = neig.GetAtomicNum()
      neig_idx = neig.GetIdx()
      if ([idx, neig_idx] not in eligibleBondsToCut) and ([neig_idx, idx] not in eligibleBondsToCut):
        if (neig_idx not in connectedAtoms) and (neig_idx not in atomList):
          connectedAtoms.append(neig_idx)
  while True:
    preLength = len(connectedAtoms)
    for idx in connectedAtoms:
      a = rdkitmol.GetAtomWithIdx(idx)
      for neig in a.GetNeighbors():
        atomNum = neig.GetAtomicNum()
        neig_idx = neig.GetIdx()
        if ([idx, neig_idx] not in eligibleBondsToCut) and ([neig_idx, idx] not in eligibleBondsToCut):
          if (neig_idx not in connectedAtoms) and (neig_idx not in atomList):
            connectedAtoms.append(neig_idx)
    
    postLength = len(connectedAtoms)
    if preLength == postLength:
      break

  return connectedAtoms

#
# helper function to deal with special case: 6+5 member rings
# if match, break the 5 member ring
#
def specialCaseSixFusedFiveMemRing(fragmol):
  bondToDelete = []
  mol = Chem.MolFromMolFile(fragmol,removeHs=False)
  pattern = Chem.MolFromSmarts('c2ccc1OCCc1c2')
  matches = mol.GetSubstructMatches(pattern)
  for match in matches:
    bondToDelete = [match[-3], match[-4]]
  
  atomsToKeep = list(range(len(mol.GetAtoms())))
  if bondToDelete != []:
    saveFragment(mol, atomsToKeep, fragmol, bondToDelete)
  return 


if __name__ == '__main__':
  
  global sdffile
  sdffile = sys.argv[1]
  atom_id = int(sys.argv[2])
  print(sdffile)
  print(atom_id) 
  growFragment(atom_id, sdffile)
