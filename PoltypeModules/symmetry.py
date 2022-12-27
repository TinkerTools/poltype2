import sys
import os
from openbabel import openbabel
from rdkit import Chem
import numpy as np
from openbabel import pybel
import databaseparser
import networkx as nx
import networkx.algorithms.isomorphism as iso


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

