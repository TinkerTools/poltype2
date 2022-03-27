import sys
import os
import openbabel
from rdkit import Chem
import numpy as np

def gen_canonicallabels(poltype,mol,rdkitmol=None):
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
    smi=Chem.MolToSmiles(rdkitmol)
    pat = Chem.MolFromSmarts(smi)
    matches=rdkitmol.GetSubstructMatches(pat, uniquify=False, maxMatches=100)
    indices=[atom.GetIdx() for atom in rdkitmol.GetAtoms()]
    matches=np.transpose(np.array(matches))
    groups=[]
    for i in range(len(matches)):
        grp=matches[i]
        if set(grp) not in groups:
            groups.append(set(grp))
    symclasstogrp={}
    idxtosymclass={}
    symclass=poltype.prmstartidx
    for grp in groups:
        symclasstogrp[symclass]=grp
        symclass+=1
    for symclass,grp in symclasstogrp.items():
        for index in grp:
            idxtosymclass[index+1]=symclass
    symmetryclass=idxtosymclass.values()
    return idxtosymclass,symmetryclass

