import sys
import os
from openbabel import openbabel
from rdkit import Chem
import numpy as np
from openbabel import pybel
import databaseparser


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
    atomindices=databaseparser.RingAtomicIndices(poltype,mol)
    for atom in rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()
        GI=ComputeGIVector(poltype,atom,rdkitmol,distmat,mol,atomindices)
        indextoGI[atomidx]=GI

    for index,GI in indextoGI.items():
        for oindex,oGI in indextoGI.items():
            if GI==oGI:
                if index not in indextomatchingindices.keys():
                    indextomatchingindices[index]=[]
                if oindex not in indextomatchingindices[index]:
                    indextomatchingindices[index].append(oindex)

    return indextomatchingindices


def ComputeGIVector(poltype,atom,rdkitmol,distmat,mol,atomindices):
    GI=[]
    atomidx=atom.GetIdx()
    distances=distmat[atomidx]
    indices=np.array(list(range(len(distances))))
    indextodist=dict(zip(indices,distances))
    newindextodist={}
    for index,dist in indextodist.items():
        oatom=rdkitmol.GetAtomWithIdx(int(index))
        atomicnum=oatom.GetAtomicNum()
        if atomicnum!=1:
            newindextodist[index]=dist

    newdistances=list(newindextodist.values())
    dtocount={}
    for d in newdistances:
        count=list(newdistances).count(d)
        dtocount[d]=count

    sorteddtocount=dict(sorted(dtocount.items()))
    distinv=[]
    for d,count in sorteddtocount.items():
        distinv.append(d*count)
    maxdist=max(newdistances)
    GI.append(maxdist)
    GI.extend(distinv)
    neighbs=[natom for natom in atom.GetNeighbors()]
    numneighbs=len(neighbs)
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
    if numneighbs==1:
        neighb=neighbs[0]
        natomidx=neighb.GetIdx()
        rings=databaseparser.GrabAllRingsContainingIndices(poltype,atomindices,[natomidx+1])
    else:
        rings=databaseparser.GrabAllRingsContainingIndices(poltype,atomindices,[atomidx+1])
    ringlens=[]
    if len(rings)>0:
        for ring in rings:
            ringlens.append(len(ring))
    GI.extend(ringlens)
    return GI
