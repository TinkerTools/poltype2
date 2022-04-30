import sys
import os
from openbabel import openbabel
from rdkit import Chem
import numpy as np
from openbabel import pybel
import databaseparser


def gen_canonicallabels(poltype,mol,rdkitmol=None,usesym=True):
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
    symmetryclass=idxtosymclass.values()
    return idxtosymclass,symmetryclass


def ComputeSymmetryTypes(poltype,distmat,rdkitmol,mol,usesym):
    indextomatchingindices={}
    indextoGI={}
    if usesym==True:
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
    else:
        for atom in rdkitmol.GetAtoms():
            atomidx=atom.GetIdx()
            indextomatchingindices[atomidx]=atomidx

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
    neighbatmnumtocount={}
    nneighbatmnumtocount={}
    if len(neighbs)==1:
        theatom=neighbs[0]
    else:
        theatom=atom
    for natom in theatom.GetNeighbors():
        natomicnum=natom.GetAtomicNum()
        if natomicnum not in neighbatmnumtocount.keys():
            neighbatmnumtocount[natomicnum]=0
        neighbatmnumtocount[natomicnum]+=1
        for nnatom in natom.GetNeighbors():
            nnatomicnum=nnatom.GetAtomicNum()
            if nnatomicnum not in nneighbatmnumtocount.keys():
                nneighbatmnumtocount[nnatomicnum]=0
            nneighbatmnumtocount[nnatomicnum]+=1

    ls=[]
    sortedneighbatmnumtocount=dict(sorted(neighbatmnumtocount.items()))
    for atmnum,count in sortedneighbatmnumtocount.items():
        ls.append(atmnum*count)
    GI.extend(ls)
    nls=[]
    sortednneighbatmnumtocount=dict(sorted(nneighbatmnumtocount.items()))
    for atmnum,count in sortednneighbatmnumtocount.items():
        nls.append(atmnum*count)
    GI.extend(nls)
    isaro=atom.GetIsAromatic()
    isinring=mol.GetAtom(atomidx+1).IsInRing()
    atomicnum=atom.GetAtomicNum()
    GI.append(atomicnum)
    atmnumtocount={}
    if numneighbs==1:
       neighb=neighbs[0]
       natomidx=neighb.GetIdx()
       rings=databaseparser.GrabAllRingsContainingIndices(poltype,atomindices,[natomidx+1])
    else:
       rings=databaseparser.GrabAllRingsContainingIndices(poltype,atomindices,[atomidx+1])
    if len(rings)>0:
        for ring in rings:
            if len(ring)<=7: # openbabel doesnt count all ring memership right someitmes when whole molecule is giant ring
                for atomindex in ring:
                    atom=rdkitmol.GetAtomWithIdx(atomindex-1)
                    atomicnum=atom.GetAtomicNum()
                    if atomicnum not in atmnumtocount.keys():
                        atmnumtocount[atomicnum]=0
                    atmnumtocount[atomicnum]+=1
        ls=[]
        sortedatmnumtocount=dict(sorted(atmnumtocount.items()))
        for atmnum,count in sortedatmnumtocount.items():
            ls.append(atmnum*count)
        GI.extend(ls)


    
    return GI
