import sys
import os
import openbabel
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import numpy as np


def CalculateSymmetry(poltype,m):
    m.UpdatePropertyCache(strict=False)
    symm=list(rdmolfiles.CanonicalRankAtoms(m, breakTies=False))
    poltype.idxtosymclass={}
    correctidxarray= np.linspace(1,m.GetNumAtoms(),num=m.GetNumAtoms(),dtype=int)
    indextosymtype=dict(zip(correctidxarray,symm))
    tempidxtosymclass={}
    for refindex,refsymtype in indextosymtype.items():
        if refindex not in tempidxtosymclass.keys():
            tempidxtosymclass[refindex]=[] 
        for index,symtype in indextosymtype.items():
            if refsymtype==symtype:
                tempidxtosymclass[refindex].append(index)
    # need to stay consistent with old poltype symmetry assisngment, so sort indexes to assign higher symmetry numbers to more massive atoms
    atomindextomass={}
    for atom in m.GetAtoms():
        idx=atom.GetIdx()
        atommass=atom.GetMass()
        atomindextomass[idx+1]=atommass
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


