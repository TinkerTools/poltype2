import sys
import os
import openbabel
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolfiles

def CalculateSymmetry(poltype,m):
    poltype.idxtosymclass={}
    smarts=rdmolfiles.MolToSmarts(m)
    correctm=mol_with_atom_index(poltype,m)
    smirks=rdmolfiles.MolToSmarts(correctm)
    correctidxarray=GrabAtomOrder(poltype,smirks)
    newsmarts=ParseSMARTSString(poltype,smarts)
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
    # need to stay consistent with old poltype symmetry assisngment, so sort indexes to assign higher symmetry numbers to more massive atoms
    atomindextomass={}
    for idx in correctidxarray:
        atom=correctm.GetAtomWithIdx(idx-1)
        atommass=atom.GetMass()
        atomindextomass[idx]=atommass
    sortedatomindextomass=sorted(atomindextomass.items(), key=lambda x: x[1])
    symclass=poltype.prmstartidx
    uniquekeysalreadyused=[]
    for item in sortedatomindextomass:
        index=item[0]
        keylist=tempidxtosymclass[index]
   
        uniquekeys=list(set(keylist))
        if set(keylist) not in uniquekeysalreadyused:
            uniquekeysalreadyused.append(set(keylist))
        else: # already assigned
            continue
        for uniquekey in uniquekeys:
            poltype.idxtosymclass[uniquekey]=symclass
        symclass+=1

def ParseSMARTSString(poltype,smarts):
    newsmarts=smarts.replace('=','-').replace('@','') # we want atoms with resonance to have same types, remove @ which designates chiral center (if have resonsnance there is no chiral center0
    # need to remove charges on atoms for example charged oxygen so you can have resonance with other oxygen
    string=''
    idxlist=[]
    for eidx in range(len(newsmarts)):
        e=newsmarts[eidx]
        if e.isdigit():
            eprev=newsmarts[eidx-1]
            if eprev=='#':
                enext=newsmarts[eidx+1]
                if enext=='-':
                    idxlist.append(eidx+1)
    for eidx in range(len(newsmarts)):
        e=newsmarts[eidx]
        if eidx in idxlist:
            continue
        else:
            string+=e
    return string
      
def GrabAtomOrder(poltype,smirks):
    atomorder=[]
    for i in range(len(smirks)):
        e=smirks[i]
        prevchar=smirks[i-1]
        try:
            nextchar=smirks[i+1]
        except:
            break
        if prevchar==':' and e.isdigit() and nextchar!='/' and nextchar!='-' and nextchar!=')' and nextchar!=':' and nextchar!='=':
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


