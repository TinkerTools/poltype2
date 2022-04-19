import sys
import os
from openbabel import openbabel
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
    distmat=Chem.rdmolops.GetDistanceMatrix(rdkitmol)
    indextomatchingindices=ComputeSymmetryTypes(poltype,distmat,rdkitmol,mol)

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


def ComputeSymmetryTypes(poltype,distmat,rdkitmol,mol):
    indextomatchingindices={}
    indextoGI={}
    for atom in rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()
        GI=ComputeGIVector(poltype,atom,rdkitmol,distmat,mol)
        indextoGI[atomidx]=GI

    for index,GI in indextoGI.items():
        for oindex,oGI in indextoGI.items():
            if GI==oGI:
                if index not in indextomatchingindices.keys():
                    indextomatchingindices[index]=[]
                if oindex not in indextomatchingindices[index]:
                    indextomatchingindices[index].append(oindex)

    return indextomatchingindices


def ComputeGIVector(poltype,atom,rdkitmol,distmat,mol):
    GI=[]
    atomidx=atom.GetIdx()
    distances=distmat[atomidx]
    maxdist=max(distances)
    GI.append(maxdist)
    numneighbs=len([natom for natom in atom.GetNeighbors()])
    neighbatmnums=0
    for natom in atom.GetNeighbors():
        natomicnum=natom.GetAtomicNum()
        natomidx=natom.GetIdx()
        if natomidx!=atomidx:
            neighbatmnums+=natomicnum
            for nnatom in natom.GetNeighbors():
                nnatomidx=nnatom.GetIdx()
                nnatomicnum=nnatom.GetAtomicNum()
                if nnatomidx!=natomidx and nnatomidx!=atomidx:
                    neighbatmnums+=nnatomicnum
    GI.append(numneighbs)
    isaro=atom.GetIsAromatic()
    isinring=mol.GetAtom(atomidx+1).IsInRing()
    GI.append(isinring)
    atomicnum=atom.GetAtomicNum()
    GI.append(atomicnum)
    GI.append(neighbatmnums)

    return GI
