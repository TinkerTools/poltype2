import sys
import os
from openbabel import openbabel
from rdkit import Chem
import numpy as np
from openbabel import pybel
import databaseparser
import networkx as nx
import networkx.algorithms.isomorphism as iso


def gen_canonical_labels_old(poltype,mol,rdkitmol=None,usesym=True,isparent=False):
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
        symm_mol_str = obConversion.WriteString(mol)
        temptotalcharge=poltype.totalcharge
        rdkitmol=Chem.MolFromMolBlock(symm_mol_str,removeHs=False,sanitize=False)
        poltype.totalcharge=None
        rdkitmol,atomindextoformalcharge=poltype.CheckInputCharge(rdkitmol)
        poltype.totalcharge=temptotalcharge
        Chem.SanitizeMol(rdkitmol)
    distmat=Chem.rdmolops.GetDistanceMatrix(rdkitmol)
    indextomatchingindices=ComputeSymmetryTypes_v1(poltype,distmat,rdkitmol,mol,usesym)

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

def gen_canonicallabels(poltype,mol,rdkitmol=None,usesym=True,isparent=False, debug=False):
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
        symm_mol_str = obConversion.WriteString(mol)
        temptotalcharge=poltype.totalcharge
        rdkitmol=Chem.MolFromMolBlock(symm_mol_str,removeHs=False,sanitize=False)
        poltype.totalcharge=None
        rdkitmol,atomindextoformalcharge=poltype.CheckInputCharge(rdkitmol)
        poltype.totalcharge=temptotalcharge
        Chem.SanitizeMol(rdkitmol)
    distmat=Chem.rdmolops.GetDistanceMatrix(rdkitmol)
    indextomatchingindices=ComputeSymmetryTypes(mol,usesym)

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

    if debug:
        res_old = gen_canonical_labels_old(poltype, mol, rdkitmol, usesym, isparent)
        if idxtosymclass != res_old[0]:
            obConversion = openbabel.OBConversion()
            obConversion.SetOutFormat('mol')
            symm_mol_str = obConversion.WriteString(mol)
            poltype.WriteToLog('Warning: the new algorithm generated different canonical labels from the old algorithm for the fragment:\n Old labels %s\nNew labels %s\n%s\n'%(res_old[0], idxtosymclass, symm_mol_str))
    return idxtosymclass,symmetryclass

def ComputeSymmetryTypes_v1(poltype,distmat,rdkitmol,mol,usesym):
    """
    Intent: Define symmetry type by a graph invarient vector. If two atoms have the same GI vector than they have the same type.
    Input: Pairwise topological distance matrix, MOL object 
    Output: Dictionary of atom index to matching atom indices
    Referenced By: gen_canonicallabels
    Description: 
    1. Iterate over all nodes twice, generate two copies of the original graph. For each copy, delete the node currently being iterated on. If the two graphs are isomorphic then they must be the same "type". 
    """
    indextomatchingindices={}
    if usesym==True:
        # STEP 1
        G = nx.Graph()
        iteratombab = openbabel.OBMolAtomIter(mol)
        for atm in iteratombab:
            G.add_nodes_from([(atm.GetIdx(), {"atomicnum": atm.GetAtomicNum()})])
        
        iterbond = openbabel.OBMolBondIter(mol) 
        for bond in iterbond:
            bgnidx=bond.GetBeginAtomIdx()
            endidx=bond.GetEndAtomIdx()
            G.add_edge(bgnidx, endidx)

        atomic_nums = nx.get_node_attributes(G, "atomicnum")
        center=nx.center(G)
        pairs=[]
        for node in G.nodes():
            if node-1 not in indextomatchingindices.keys():
                indextomatchingindices[node-1]=[]
            if node-1 not in indextomatchingindices[node-1]:
                indextomatchingindices[node-1].append(node-1)
            for onode in G.nodes():
                if onode-1 not in indextomatchingindices.keys():
                    indextomatchingindices[onode-1]=[]
                if onode-1 not in indextomatchingindices[onode-1]:
                    indextomatchingindices[onode-1].append(onode-1)
                node_deg=G.degree(node)
                onode_deg=G.degree(onode)
                if node_deg==onode_deg:
                    node_atomic_num=atomic_nums[node]
                    onode_atomic_num=atomic_nums[onode]
                    if node_atomic_num==onode_atomic_num:
                        pair=set([node,onode])
                        if node!=onode and pair not in pairs:
                            node_dist=[nx.shortest_path_length(G, source=node, target=i) for i in center].sort()
                            onode_dist=[nx.shortest_path_length(G, source=onode, target=i) for i in center].sort()
                            if node_dist==onode_dist:
                                G1=G.copy()
                                G2=G.copy()
                                G1.remove_node(node)
                                G2.remove_node(onode)
                                isiso=nx.is_isomorphic(G1,G2) # makes database searching a bit slower cause isomorphism checking is slow, but this is general solution
                                pairs.append(pair)
                                if isiso==True:
                                    if onode-1 not in indextomatchingindices[node-1]:
                                        indextomatchingindices[node-1].append(onode-1)
                                    if node-1 not in indextomatchingindices[onode-1]:
                                        indextomatchingindices[onode-1].append(node-1)
    else:
        # STEP 2
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




def ComputeGIVector(poltype,atom,rdkitmol,distmat,mol,atomindices):
    """
    Intent:
    Input: Atom of interest object, rdkit MOL object, pairwise topological distance matrix, openbabel MOL object, list of list of ring atomic indices.
    Output: Graph Invariant Vector
    Referenced By: ComputeSymmetryTypes
    Description:
    1. Grab topological distances to every other atom in graph.  
    2. Determine unique topological distances from atom.
    3. For each unique topological distance from atom, find all other atoms at same distance away and save the atomic number of that atom. Also save the neighboring environment (all first neighbors) of that atomic number at distance d away. 
    4. Then count how many atoms of same atomic number and neiboring enviorment are at that distance away from original atom.
    5. Now for every distance d, away from atom of interest a, you have a map of atomic number to number of occurances of atomic number (with neiboring environment) at that distance. This is in essence what the symmetry type number is.
    6. Now for each d, 0-max(d), sort the map of atomic number to occurances and add to one large array. Need to sort since the maps may have different order but same information (could be same type but have maps in different order).
    7. If you are a singly valent atom such as hydrogen check for ring membership of neighbor else check ring membership of atom of interest and add the ring atoms it belongs to to the GI array. 
    """
    GI=[]
    atomidx=atom.GetIdx()
    # STEP 1
    distances=distmat[atomidx]
    indices=np.array(list(range(len(distances))))
    indextodist=dict(zip(indices,distances))
    distancetoatomicnumneighbstocount={}
    # STEP 2
    uniquedistances=set(distances)
    # STEP 3
    for d in uniquedistances:
        if d not in distancetoatomicnumneighbstocount.keys():
            distancetoatomicnumneighbstocount[d]={}
            atomicnumneighbstocount={}
            for index,dist in indextodist.items():
                oatom=rdkitmol.GetAtomWithIdx(int(index))
                atomicnum=oatom.GetAtomicNum()
                if dist==d:
                    natomicnumtocount={}
                    natomicnums=[]
                    for natom in oatom.GetNeighbors():
                        natomicnum=natom.GetAtomicNum()
                        natomicnums.append(natomicnum)
                    uniquenatomicnums=set(natomicnums)
                    for natomnum in uniquenatomicnums:
                        count=natomicnums.count(natomnum)
                        natomicnumtocount[natomnum]=count
                    sortednatomicnumtocount=dict(sorted(natomicnumtocount.items()))
                    ls=[atomicnum]
                    for atmnum,count in sortednatomicnumtocount.items():
                        ls.append(atmnum*count)
                    ls=tuple(ls)
                    if ls not in atomicnumneighbstocount.keys():
                        atomicnumneighbstocount[ls]=0
                    atomicnumneighbstocount[ls]+=1
                    
            # STEP 4
            distancetoatomicnumneighbstocount[d]=atomicnumneighbstocount
    # STEP 5
    sorteddistancetoatomicnumcounts=dict(sorted(distancetoatomicnumneighbstocount.items()))

    # STEP 6 
    for d,dic in sorteddistancetoatomicnumcounts.items():
        newls=[]
        for neighbnum,count in dic.items():
            newls.extend(list(neighbnum)+[count])

        newls.sort()
        GI.extend(newls)
    

    neighbs=[natom for natom in atom.GetNeighbors()] # now need to include bond order information to atomic element for neighhbors since found special case that didnt work
    neighbindices=[natom.GetIdx() for natom in neighbs]
    

    # STEP 7
    numneighbs=len(neighbs)
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
