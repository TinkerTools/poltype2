import sys
import os
from openbabel import openbabel
from rdkit import Chem
import numpy as np
from openbabel import pybel
import databaseparser


def gen_canonicallabels(poltype,mol,rdkitmol=None,usesym=True,isparent=False):
    """
    Intent: Using ring membership, and graph distance to all other elements compute symmetry type.
    Input: Openbabel MOL object.
    Output: Dictionary of atom index to symmetry type.
    Referenced By: GenerateParameters in poltype.py 
    Description: First find matching indices via symmetry type, then sort atom types so that groups of heavier atoms get lower type numbers. Optionally read in custom index to type number file if user provides input. 
    """
    if rdkitmol==None:
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('mol')
        obConversion.WriteFile(mol,'symm.mol')
        temptotalcharge=poltype.totalcharge
        rdkitmol=Chem.MolFromMolFile('symm.mol',removeHs=False,sanitize=False)
        poltype.totalcharge=None
        rdkitmol,atomindextoformalcharge=poltype.CheckInputCharge(rdkitmol)
        poltype.totalcharge=temptotalcharge
        Chem.SanitizeMol(rdkitmol)
    distmat=Chem.rdmolops.GetDistanceMatrix(rdkitmol)
    indextomatchingindices=ComputeSymmetryTypes(poltype,distmat,rdkitmol,mol,usesym)

    groups=[] 
    grouptoheavy={}
    for index,matchingindices in indextomatchingindices.items():
        atom=rdkitmol.GetAtomWithIdx(index)
        atomnum=atom.GetAtomicNum()
        heavy=False
        if atomnum!=1:
            heavy=True
        if set(matchingindices) not in groups:
            groups.append(set(matchingindices))
            grouptoheavy[tuple(set(matchingindices))]=heavy
        
    sortedgroups=[]
    for group,heavy in grouptoheavy.items():
        if heavy==True:
            if group not in sortedgroups:
                sortedgroups.append(group)
    for group,heavy in grouptoheavy.items():
        if heavy==False:
            if group not in sortedgroups:
                sortedgroups.append(group)
    symclasstogrp={}
    idxtosymclass={}
    symclass=poltype.prmstartidx
    for grp in sortedgroups:
        symclasstogrp[symclass]=grp
        symclass+=1
    for symclass,grp in symclasstogrp.items():
        for index in grp:
            idxtosymclass[index+1]=symclass
    if poltype.indextotypefile!=None and isparent==True:
        idxtosymclass=ReadCustomIndexToTypeFiles(poltype,poltype.indextotypefile)
    symmetryclass=idxtosymclass.values()

    return idxtosymclass,symmetryclass


def ComputeSymmetryTypes(poltype,distmat,rdkitmol,mol,usesym):
    """
    Intent: Define symmetry type by a graph invarient vector. If two atoms have the same GI vector than they have the same type.
    Input: Pairwise topological distance matrix, MOL object 
    Output: Dictionary of atom index to matching atom indices
    Referenced By: gen_canonicallabels
    Description: 
    1. Compute GI vector for every atom
    2. Iterate over GI vectors for each atom and check if their is GI vector match, if so add to same group
    3. If usesym==False then use atom index instead of type numbers.
    """
    fprints = fingerprint(mol)
    indextomatchingindices={}
    indextoGI={}
    if usesym==True:
        atomindices=databaseparser.RingAtomicIndices(poltype,mol)
        # STEP 1
        for atom in rdkitmol.GetAtoms():
            atomidx=atom.GetIdx()
            GI=fprints[atomidx]
            indextoGI[atomidx]=GI
        # STEP 2
        for index,GI in indextoGI.items():
            for oindex,oGI in indextoGI.items():
                if GI==oGI:
                    if index not in indextomatchingindices.keys():
                        indextomatchingindices[index]=[]
                    if oindex not in indextomatchingindices[index]:
                        indextomatchingindices[index].append(oindex)
    else:
        # STEP 3
        for atom in rdkitmol.GetAtoms():
            atomidx=atom.GetIdx()
            indextomatchingindices[atomidx]=atomidx

    return indextomatchingindices


def ReadCustomIndexToTypeFiles(poltype,indextotypefile):
    """
    Intent: If user provides file mapping atom index to type, read that in as symmetry types.
    Input: indextotypefile
    Output:
    Referenced By: gen_canonicallabels
    Description: indextomatchingindices
    """
    indextomatchingindices={}
    temp=open(indextotypefile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            index=int(linesplit[0])
            typenumber=int(linesplit[1])
            indextomatchingindices[index]=typenumber
    return indextomatchingindices


def fingerprint(mol):
  fprints = []
  iteratombab = openbabel.OBMolAtomIter(mol)
  atoms = []
  elements = []
  for atm in iteratombab:
      atoms.append(str(atm.GetIdx()))
      elements.append(str(atm.GetAtomicNum()))
  
  connections = []
  for atomidx in atoms:
    atomidx=int(atomidx)
    conns=[]
    iterbond = openbabel.OBMolBondIter(mol) 
    for bond in iterbond:
        bgnidx=bond.GetBeginAtomIdx()
        endidx=bond.GetEndAtomIdx()
        if bgnidx==atomidx:
            if endidx not in conns:
                conns.append(str(endidx))
        elif endidx==atomidx:
            if bgnidx not in conns:
                conns.append(str(bgnidx))
    connections.append(conns)
  

  
  
  atom_ele_dict = dict(zip(atoms, elements))
  atom_con_dict = {}
  for atom, con in zip(atoms,connections):
    con_ele = [atom_ele_dict[c] for c in con] 
    constr = ''.join(sorted(con_ele)) 
    atom_con_dict[atom] = constr

  level = 5 
  if level > 1:
    atom_con_dict2 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict[c] for c in cons])
      atom_con_dict2[atom] = ''.join(sorted(newstr))

  # level 3 is good for chain molecules 
  if level > 2:
    atom_con_dict3 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict2[c] for c in cons])
      atom_con_dict3[atom] = ''.join(sorted(newstr))

  # level 4 is needed for ring molecules 
  if level > 3:
    atom_con_dict4 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict3[c] for c in cons])
      atom_con_dict4[atom] = ''.join(sorted(newstr))
  
  if level > 4:
    atom_con_dict5 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict4[c] for c in cons])
      atom_con_dict5[atom] = ''.join(sorted(newstr))
  
  for atom in atoms:
    fprints.append(atom_ele_dict[atom] + '-' + str(''.join(sorted(atom_con_dict[atom] + atom_con_dict2[atom] + atom_con_dict3[atom] + atom_con_dict4[atom] + atom_con_dict5[atom]))))
  
  return fprints