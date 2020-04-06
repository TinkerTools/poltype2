import sys
import os
import openbabel
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolfiles

def CalculateSymmetry(poltype,m):
    poltype.idxtosymclass={}
    smarts=rdmolfiles.MolToSmarts(m)
    m=mol_with_atom_index(poltype,m)
    smirks=rdmolfiles.MolToSmarts(m)
    correctidxarray=GrabAtomOrder(poltype,smirks)
    newsmarts=smarts.replace('=','-') # we want atoms with resonance to have same types
    rdkitmol=rdmolfiles.MolFromSmarts(newsmarts)
    m=mol_with_atom_index(poltype,rdkitmol)
    smirks=rdmolfiles.MolToSmarts(m)
    incorrectidxarray=GrabAtomOrder(poltype,smirks)
    incorrectidxtocorrectidx=dict(zip(incorrectidxarray,correctidxarray))
    matches = rdkitmol.GetSubstructMatches(rdkitmol,uniquify=False)
    tempidxtosymclass={}
    for match in matches:
        for idx1,idx2 in enumerate(match):
            incorrectatomindex1=incorrectidxarray[idx1] 
            correctatomindex1=incorrectidxtocorrectidx[incorrectatomindex1]
            incorrectatomindex2=incorrectidxarray[idx2] 
            correctatomindex2=incorrectidxtocorrectidx[incorrectatomindex2]
            if (correctatomindex1) not in tempidxtosymclass.keys():
                tempidxtosymclass[correctatomindex1]=[]
            tempidxtosymclass[correctatomindex1].append(correctatomindex2)

    symclass=poltype.prmstartidx
    for key,keylist in tempidxtosymclass.items():
        uniquekeys=list(set(keylist))
        for uniquekey in uniquekeys:
            poltype.idxtosymclass[uniquekey]=symclass
        symclass+=1

def GrabAtomOrder(poltype,smirks):
    atomorder=[]
    for i in range(len(smirks)):
        e=smirks[i]
        prevchar=smirks[i-1]
        try:
            nextchar=smirks[i+1]
        except:
            break
        if prevchar==':' and e.isdigit() and nextchar!='-' and nextchar!=')' and nextchar!=':' and nextchar!='=':
            atomindex=GrabAtomIndex(poltype,i,smirks)
            atomorder.append(atomindex)
    return atomorder

def mol_with_atom_index(poltype,mol):
    atoms = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        atomnum=atom.GetAtomicNum()
        rdkitatomidx=atom.GetIdx()
        OBatomidx=rdkitatomidx+1
        atom.SetProp( 'molAtomMapNumber', str(OBatomidx ))
    return mol

def GrabAtomIndex(poltype,i,smirks):
    num=[]
    for j in range(i,len(smirks)):
        char=smirks[j]
        if char.isdigit():
            num.append(char)
        if char==']':
            break
    atomindex=int(''.join(num))
    return atomindex 


