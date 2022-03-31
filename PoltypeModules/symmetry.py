import sys
import os
import openbabel
from rdkit import Chem
import numpy as np
from openbabel import pybel

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
    smilesorder=rdkitmol.GetProp('_smilesAtomOutputOrder')
    smilesorder=smilesorder[1:-2]
    smilesorder=smilesorder.split(',')
    smilesorder=[int(i) for i in smilesorder]
    indices=[i.GetIdx() for i in rdkitmol.GetAtoms()]
    indextotrueindex=dict(zip(indices,smilesorder))
    pbmol=pybel.readstring("smi", smi)
    mappings = pybel.ob.vvpairUIntUInt()
    success = pybel.ob.FindAutomorphisms(pbmol.OBMol, mappings)
    indextomatchingindices={}
    for maplist in mappings:
        for mapitem in maplist:
            index=indextotrueindex[mapitem[0]]
            matchingindex=indextotrueindex[mapitem[1]]

            if index not in indextomatchingindices.keys():
                indextomatchingindices[index]=[]
            if matchingindex not in indextomatchingindices[index]:
                indextomatchingindices[index].append(matchingindex)
            if index not in indextomatchingindices[index]:
                indextomatchingindices[index].append(index)

           
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
    symmetryclass=idxtosymclass.values()
    return idxtosymclass,symmetryclass

