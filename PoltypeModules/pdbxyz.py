import os
import submitjobs as submit
from openbabel import openbabel
import re
import shutil
import warnings
import time
import sys
from PyAstronomy import pyasl
import numpy as np
from pathlib import Path
import json
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdmolfiles


def GenerateProteinTinkerXYZFile(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    proteinindextocoordinates=GenerateUncomplexedProteinPDBFromComplexedPDB(poltype)
    poltype.uncomplexedxyzname=poltype.uncomplexedproteinpdbname.replace('.pdb','.xyz')
    poltype.complexedxyzname=poltype.uncomplexedxyzname.replace('.xyz','_comp.xyz')
    poltype.receptorligandxyzfilename=poltype.complexedxyzname
    poltype.ReadReceptorCharge()
    chainnum=DetectNumberOfChains(poltype,poltype.uncomplexedproteinpdbname)
    if chainnum==1:

        cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+poltype.prmfilepath
    else:
        cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+'ALL'+' '+poltype.prmfilepath

    submit.call_subsystem(poltype,cmdstr,wait=True)   
    
    uncomplexedatomnum=len(proteinindextocoordinates.keys()) # just make sure to move HETATMS between chains after protein xyz
    proteinpdbindextoxyzindex=MapXYZToPDB(poltype,poltype.uncomplexedproteinpdbname,poltype.uncomplexedxyzname,uncomplexedatomnum)
    atoms,coord,order,types,connections=readTXYZ(poltype,poltype.uncomplexedxyzname)
    newuncomplexedatomnum=len(atoms)
    poltype.totalreceptornumber=len(atoms)
    poltype.receptorindices=list(range(1,len(atoms)+1))
    shift= newuncomplexedatomnum- uncomplexedatomnum
    if shift!=0:
        string='WARNING! Missing atoms from original PDB have been added by PDBXYZ. Number of atoms added = '+str(shift)
        warnings.warn(string)
        poltype.WriteToLog(string)
    
    indextocoordinates,annihilateindextocoordinates,lastligindex,allligands,ligands,listofindextopdbindex,diff,ligandpdbindextoxyzindex=GrabLigandCoordinates(poltype,newuncomplexedatomnum)
    poltype.ligandindices[0]=list(annihilateindextocoordinates.keys())
    poltype.allligandindices[0]=list(indextocoordinates.keys())
    poltype.allligands=allligands
    poltype.ligands=ligands
    nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,hetatmpdbindextoxyzindex=GrabNonLigandHETATMInfo(poltype,poltype.complexedproteinpdbname,indextocoordinates,proteinindextocoordinates,lastligindex) 
    poltype.hetatmindices=list(nonlighetatmindextocoordinates.keys())
    FindWatersIonsInPocketToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,indextocoordinates)
    GenerateComplexedTinkerXYZFile(poltype,poltype.uncomplexedxyzname,indextocoordinates,newuncomplexedatomnum,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,listofindextopdbindex,diff,ligandpdbindextoxyzindex)
    poltype.pdbindextoxyzindex=CombinePDBMaps(poltype,proteinpdbindextoxyzindex,ligandpdbindextoxyzindex,hetatmpdbindextoxyzindex)
    if poltype.makexyzonly==True:
        sys.exit()

def CombinePDBMaps(poltype,proteinpdbindextoxyzindex,ligandpdbindextoxyzindex,hetatmpdbindextoxyzindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    pdbindextoxyzindex={}
    pdbindextoxyzindex.update(proteinpdbindextoxyzindex)
    pdbindextoxyzindex.update(ligandpdbindextoxyzindex)
    pdbindextoxyzindex.update(hetatmpdbindextoxyzindex)
    return pdbindextoxyzindex


def MapXYZToPDB(poltype,uncomplexedproteinpdbname,uncomplexedxyzname,uncomplexedatomnum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    pdbrdkitmol=Chem.MolFromPDBFile(uncomplexedproteinpdbname,removeHs=False,sanitize=False)   
    pdbindextocoords=poltype.GrabIndexToCoordinatesPymol(pdbrdkitmol)
    xyzfilename=poltype.ConvertTinkerXYZToCartesianXYZ(uncomplexedxyzname)
    xyzmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('xyz')
    obConversion.ReadFile(xyzmol,xyzfilename)
    obConversion.SetOutFormat('mol')
    tempname=xyzfilename.replace('.xyz','.mol')
    obConversion.WriteFile(xyzmol,tempname)
    xyzrdkitmol=Chem.MolFromMolFile(tempname,removeHs=False,sanitize=False)
    xyzindextocoords=poltype.GrabIndexToCoordinatesPymol(xyzrdkitmol)
    proteinpdbindextoxyzindex=MapXYZCoordstoPDBCoords(poltype,pdbindextocoords,xyzindextocoords,uncomplexedatomnum)

    return proteinpdbindextoxyzindex


def MapXYZCoordstoPDBCoords(poltype,pdbindextocoords,xyzindextocoords,uncomplexedatomnum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    proteinpdbindextoxyzindex={}
    for pdbindex,pdbcoords in pdbindextocoords.items():
        for xyzindex,xyzcoords in xyzindextocoords.items():
            if pdbcoords==xyzcoords:
                proteinpdbindextoxyzindex[pdbindex]=xyzindex
    indices=list(range(1,uncomplexedatomnum+1))
    for index in indices:
        if index not in proteinpdbindextoxyzindex.keys():
            proteinpdbindextoxyzindex[index]=index

    return proteinpdbindextoxyzindex







def ConvertTinkerXYZToPDB(poltype,xyzfile): 
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    cmdstr=poltype.xyzpdbpath+' '+xyzfile+' '+poltype.prmfilepath 
    submit.call_subsystem(poltype,cmdstr,wait=True)
    pdbfile=xyzfile.replace('.xyz','.pdb')
    return pdbfile


def GrabNonLigandHETATMInfo(poltype,complexedproteinpdbname,indextocoordinates,proteinindextocoordinates,lastligindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    nonlighetatmindextocoordinates={}
    nonlighetatmindextotypenum={}
    nonlighetatmindextoconnectivity={}
    nonlighetatmindextoelement={}
    pdbindextoxyzindex={}
    temp=open(complexedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,complexedproteinpdbname)
    iteratom = openbabel.OBMolAtomIter(pdbmol)
    an = pyasl.AtomicNo()
    firstindex=lastligindex+1
    count=0
    nonionelements=['H','S','C','O','Cl','I','Br','F','N']
    for atom in iteratom:
        index=atom.GetIdx()
        atomicnum=atom.GetAtomicNum()
        coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
        element=an.getElSymbol(atomicnum)
        if coords not in proteinindextocoordinates.values() and coords not in indextocoordinates.values():
            newindex=firstindex+count
            nonlighetatmindextocoordinates[newindex]=coords
            typenum=GrabNonLigandHETATMType(poltype,element)
            if typenum==None:
                continue
            nonlighetatmindextotypenum[newindex]=typenum
            shift=newindex-index 
            if element in nonionelements:
                connectivity,atomicnums,oldconn=GrabConnectivity(poltype,pdbmol,index,shift)
            else:
                connectivity=[]
            if element=='O': # if like some buffer molecule user forgot to add to ligandxyzfilenamelist, make sure O only for water
                watero=True
                for oatomicnum in atomicnums:
                    if oatomicnum!=1:
                        watero=False
                if watero==False:
                    continue
            if element=='H': # could be H on non ligand molecule (besides water)
                neighb=oldconn[0]
                natomicnum=atomicnums[0]
                if natomicnum!=8:
                    continue
                oconnectivity,oatomicnums,ooldconn=GrabConnectivity(poltype,pdbmol,neighb,shift)
                watero=True
                for oatomicnum in oatomicnums:
                    if oatomicnum!=1:
                        watero=False
                if watero==False:
                    continue


                     
            nonlighetatmindextoconnectivity[newindex]=connectivity
            nonlighetatmindextoelement[newindex]=element
            count+=1
            pdbindextoxyzindex[index]=newindex
    return nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,pdbindextoxyzindex


def GrabConnectivity(poltype,pdbmol,index,shift):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    conn=[]
    oldconn=[]
    els=[]
    bonditer=openbabel.OBMolBondIter(pdbmol)
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        obgnatom=pdbmol.GetAtom(obgnidx)
        oendatom=pdbmol.GetAtom(oendidx)
        obgnatomicnum=obgnatom.GetAtomicNum()
        oendatomicnum=oendatom.GetAtomicNum()

        bnd=[oendidx,obgnidx]
        if oendidx==index:
            newbgn=obgnidx+shift
            if newbgn not in conn:
                conn.append(newbgn)
                els.append(obgnatomicnum)
                oldconn.append(obgnidx)
        elif obgnidx==index:
            newend=oendidx+shift
            if newend not in conn:
                conn.append(newend)
                els.append(oendatomicnum)
                oldconn.append(oendidx)
    return conn,els,oldconn




def GrabNonLigandHETATMType(poltype,element): # assumes besides ligand and protein, only water and ions
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    typenum=None # if user has buffer molecule not in ligandxyzfilenamelist etc, then skip
    if element=='O':
        typenum=poltype.waterOtypenum
    elif element=='H':
        typenum=poltype.waterHtypenum
    else:
        if element in poltype.elementsymtotinktype.keys():
            typenum=poltype.elementsymtotinktype[element]

    return typenum


def readTXYZ(poltype,TXYZ):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(TXYZ,'r')
    lines = temp.readlines()[1:] #TINKER coordinate starts from second line
    atoms=[];coord=[]
    order=[];types=[];connections=[]
    for line in lines:
        data=line.split()
        order.append(data[0])
        types.append(data[5])
        connections.append(data[6:])
        atoms.append(data[1])
        coord.append([float(data[2]), float(data[3]), float(data[4])])
    return atoms,coord,order, types, connections




def SMILESMatches(poltype,ligandsmiles,rdkitmol,pdbmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ligandsmiles=ligandsmiles.replace('=','~').replace('-','~').replace('/','') # sometimes bond order can be different from PDB so remove bond order for detection :)
    p = Chem.MolFromSmarts(ligandsmiles)
    matches=rdkitmol.GetSubstructMatches(p)
    newmatches=[]
    for match in matches:
        newmatch=[i+1 for i in match]
        newmatches.append(newmatch)
    if len(newmatches)==0:
        raise ValueError(' ligandsmiles not matching anything in PDB! '+ligandsmiles)
    return newmatches,p




def GrabLigandCoordinates(poltype,proteinatomnum): # appending after protein, before pocket water and ions
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ligandpdbindextoxyzindex={}
    firstatomidx=proteinatomnum+1
    indextocoordinates={}
    annihilateindextocoordinates={}
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,poltype.complexedproteinpdbname)
    rdkitmol=Chem.MolFromPDBFile(poltype.complexedproteinpdbname,removeHs=False,sanitize=False)
    allmatches=[]
    allligands=[]
    ligands=[]
    oldindextonewindex={}
    listofindextopdbindex=[]
    smilestoindicesalreadyused={}
    for ligandsmilesidx in range(len(poltype.ligandsmileslist)):
        ligandsmiles=poltype.ligandsmileslist[ligandsmilesidx].replace('=','~').replace('-','~').replace('/','')
        xyz=poltype.ligandxyzfilenamelist[ligandsmilesidx]
        molname=xyz.replace('.xyz','.mol')
        ligm=Chem.MolFromMolFile(molname,removeHs=False,sanitize=False)
        matches,p=SMILESMatches(poltype,ligandsmiles,rdkitmol,pdbmol)
        temp=[]
        firstindices=[]
        if ligandsmiles not in smilestoindicesalreadyused.keys():
            smilestoindicesalreadyused[ligandsmiles]=[]
        for match in matches:
            match=list(match)
            xyzmatches,p=SMILESMatches(poltype,ligandsmiles,ligm,pdbmol)
            indices=xyzmatches[0]
            indextopdbindex=dict(zip(indices,match)) # in smiles basis set, need from xyz basis set
            match.sort()
            if match[0] not in smilestoindicesalreadyused[ligandsmiles]:
                firstindices.append(match[0])
                temp.append(indextopdbindex)
            if match not in allligands:
                allligands.append(match)
            for i in match:
                if i not in allmatches:
                    allmatches.append(i)
        minvalue=min(firstindices)
        minindex=firstindices.index(minvalue)
        indextopdbindex=temp[minindex]
        listofindextopdbindex.append(indextopdbindex)
        smilestoindicesalreadyused[ligandsmiles].append(minvalue)
    for xyz,smiles in poltype.ligandxyztosmiles.items():
        if xyz in poltype.annihilateligandxyzfilenamelist:
            ligxyzfnamelist=[]
            annligxyzfnamelist=[]
            for oxyz,osmiles in poltype.ligandxyztosmiles.items():
                if smiles==osmiles:
                    ligxyzfnamelist.append(oxyz)
                    if oxyz in poltype.annihilateligandxyzfilenamelist:
                        annligxyzfnamelist.append(oxyz)
            acceptablematchindices=[]
            for oxyz in annligxyzfnamelist:
                idx=ligxyzfnamelist.index(oxyz)
                acceptablematchindices.append(idx)
            matches,p=SMILESMatches(poltype,smiles,rdkitmol,pdbmol)
            firstindextomatch={}
            for match in matches:
                match=list(match)
                match.sort()
                firstindex=match[0]
                firstindextomatch[firstindex]=match
            firstindices=list(firstindextomatch.keys())
            firstindices.sort()
            truematches=[]
            for idx in range(len(firstindices)):
                firstindex=firstindices[idx]
                if idx in acceptablematchindices:
                    match=firstindextomatch[firstindex]
                    if match not in ligands:
                        ligands.append(match)
                    truematches.extend(match)
    ligands=SortLigandIndices(poltype,ligands)
    allligands=SortLigandIndices(poltype,allligands)
    firstligindex=allligands[0][0]
    diff=firstligindex-firstatomidx # difference accounting between last protein atom and first ligand atom
    firstatomidx+=diff
    atomiter=openbabel.OBMolAtomIter(pdbmol)
    count=0
    for atom in atomiter:
        atomidx=atom.GetIdx()
        res=atom.GetResidue()
        resname=res.GetName()
        if atomidx in allmatches:
            coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
            newindex=firstatomidx+count-diff
            indextocoordinates[newindex]=coords
            oldindextonewindex[atomidx]=newindex
            ligandpdbindextoxyzindex[atomidx]=newindex
            if atomidx in truematches:
                annihilateindextocoordinates[newindex]=coords
            count+=1
    if len(indextocoordinates.keys())==0:
        raise ValueError('Complexed PDB missing atoms')
    newallligands=[]
    for match in allligands:
        newallmatch=[]
        for index in match:
            nnewindex=oldindextonewindex[index]
            newallmatch.append(nnewindex)
        newallligands.append(newallmatch)
    newligands=[]
    for match in ligands:
        newmatch=[]
        for index in match:
            nnewindex=oldindextonewindex[index]
            newmatch.append(nnewindex)
        newligands.append(newmatch)

    return indextocoordinates,annihilateindextocoordinates,newindex,newallligands,newligands,listofindextopdbindex,diff ,ligandpdbindextoxyzindex


def SortLigandIndices(poltype,array):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newarray=[]
    firstindicestoarray={}
    for ls in array:
        firstindex=ls[0]
        firstindicestoarray[firstindex]=ls
    for key in sorted(firstindicestoarray):
        newarray.append(firstindicestoarray[key])


    return newarray


def GenerateComplexedTinkerXYZFile(poltype,uncomplexedxyzname,indextocoordinates,uncomplexedatomnum,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,listofindextopdbindex,diff,ligandpdbindextoxyzindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(uncomplexedxyzname,'r')
    results=temp.readlines()
    temp.close() 
    temp=open(poltype.complexedxyzname,'w')
    newatomnum=uncomplexedatomnum+len(indextocoordinates)+len(nonlighetatmindextocoordinates.keys())
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            newline=str(newatomnum)+'\n' 
        else:
            newline=line
        temp.write(newline)
    lastindex=uncomplexedatomnum
    for xyzidx in range(len(poltype.ligandxyzfilenamelist)):
        xyz=poltype.ligandxyzfilenamelist[xyzidx]
        indextopdbindex=listofindextopdbindex[xyzidx]
        atoms,coord,order,types,connections=readTXYZ(poltype,xyz)
        lastindex=WriteLigandToXYZ(poltype,temp,atoms,order,lastindex,indextocoordinates,types,connections,indextopdbindex,diff,ligandpdbindextoxyzindex)
    WriteHETATMToXYZ(poltype,temp,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement)
    temp.close()

def WriteLigandToXYZ(poltype,temp,atoms,order,lastindex,indextocoordinates,types,connections,indextopdbindex,diff,ligandpdbindextoxyzindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    shift=lastindex
    for idx in range(len(atoms)):
        element=atoms[idx]
        oldindex=int(order[idx])
        oldpdbindex=indextopdbindex[oldindex]
        pdbindex=oldpdbindex-diff
        newindex=oldindex+shift
        xyzindex=ligandpdbindextoxyzindex[oldpdbindex]
        coords=indextocoordinates[xyzindex]
        x=coords[0]
        y=coords[1]
        z=coords[2]
        typenum=types[idx]
        conns=connections[idx]
        conns=[int(i) for i in conns]
        conns=[i+shift for i in conns]
        newline='    '+str(newindex)+'  '+element+'     '+str(x)+'   '+str(y)+'   '+str(z)+'    '+str(typenum)+'     '
        for con in conns:
            newline+=str(con)+'     '
        newline+='\n'
        temp.write(newline)
    return newindex

def WriteHETATMToXYZ(poltype,temp,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for index,coords in nonlighetatmindextocoordinates.items():
        typenum=nonlighetatmindextotypenum[index]
        conns=nonlighetatmindextoconnectivity[index]
        element=nonlighetatmindextoelement[index]
        x=coords[0]
        y=coords[1]
        z=coords[2]
        newline='    '+str(index)+'  '+element+'     '+str(x)+'   '+str(y)+'   '+str(z)+'    '+str(typenum)+'     '
        for con in conns:
            newline+=str(con)+'     '
        newline+='\n'
        temp.write(newline)


def GenerateUncomplexedProteinPDBFromComplexedPDB(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    uncomplexedatomindices=[]
    indextocoordinates={}
    temp=open(poltype.complexedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            resname=line[17:21].strip()
            linesplit=line.split()
            index=int(line[6:11].strip())
            if resname in poltype.knownresiduesymbs:
                uncomplexedatomindices.append(index) 
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,poltype.complexedproteinpdbname)
    totalatoms=pdbmol.NumAtoms()
    atomiter=openbabel.OBMolAtomIter(pdbmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if atomidx in uncomplexedatomindices:
            coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
            indextocoordinates[atomidx]=coords

    
    indexestodelete=[]
    for i in range(1,totalatoms+1):
        if i not in uncomplexedatomindices:
            indexestodelete.append(i)
    
    poltype.uncomplexedproteinpdbname='uncomplexed.pdb'
    temp=open(poltype.uncomplexedproteinpdbname,'w')
    for line in results:
       if 'ATOM' in line:
           index=int(line[6:11].strip())
           if index in indexestodelete:
               pass
           else:
               temp.write(line)
       else:
           if 'HETATM' not in line:
               temp.write(line)


    temp.close()

    return indextocoordinates


def GeneratePDBFileFromXYZ(poltype,xyzfile,ligandindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    cmdstr=poltype.xyzpdbpath+' '+xyzfile+' '+poltype.prmfilepath 
    submit.call_subsystem(poltype,cmdstr,wait=True)
    newpdb=xyzfile.replace('.xyz','.pdb')
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,newpdb)
    obConversion.SetOutFormat('xyz')
    tempname='temp.xyz'
    obConversion.WriteFile(pdbmol,tempname)
    obConversion.SetInFormat('xyz')
    obConversion.SetOutFormat('pdb')
    finalname='topologyguess.pdb'
    newpdbmol=openbabel.OBMol()
    obConversion.ReadFile(newpdbmol,tempname)
    obConversion.WriteFile(newpdbmol,finalname)
    tempname=finalname.replace('.pdb','_TEMP.pdb')
    temp=open(finalname,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempname,'w')
    for line in results:
        linesplit=re.split(r'(\s+)', line)
        if 'HETATM' in line:
            othersplit=line.split()
            index=int(othersplit[1])
            if index in ligandindices:
                if 'UNL' in linesplit:
                    idx=linesplit.index('UNL')
                    linesplit[idx]='LIG'
                    line=''.join(linesplit)
        temp.write(line)
    temp.close()
    os.rename(tempname,finalname)
    if not os.path.isdir(poltype.visfolder):
        os.mkdir(poltype.visfolder)
    newpath=os.path.join(poltype.visfolder,finalname)
    shutil.copy(finalname,newpath)


def DetectNumberOfChains(poltype,pdbfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(pdbfile,'r')
    results=temp.readlines()
    temp.close()
    chainls=[]
    for line in results:
        if 'ATOM' in line and 'REMARK' not in line:
            linesplit=line.split()
            chain=line[21]
            if chain not in chainls:
                chainls.append(chain)
    chainnum=len(chainls)
    return chainnum


def WriteSEQFile(poltype,filename,code):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    
    e = Environ()
    m = Model(e, file=filename)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    seqfilename=filename.replace('.pdb','.seq')
    aln.write(file=filename.replace('.pdb','.seq'))
    return seqfilename


def FindCurrentResidueArray(poltype,filename,chaintomissingresnums={}):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    chaintocurrentresnumtothreeletter={}
    chainlettertonum={'A':0,'B':1,'C':2,'D':3}
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    threeletterseq=[]
    chainseq=[]
    resnumseq=[]
    for line in results:
        linesplit=line.split()
        if 'ATOM' in line and line[:4]=='ATOM':
            resnum=int(line[22:26])
            chain=line[21]
            chainnum=chainlettertonum[chain]
            resname=line[17:21].strip()
            foundmissing=False
            if chainnum in chaintomissingresnums.keys():
                missingresnums=chaintomissingresnums[chainnum]
                if resnum in missingresnums:
                    foundmissing=True
            if foundmissing==False:
                if chainnum not in chaintocurrentresnumtothreeletter.keys():
                    chaintocurrentresnumtothreeletter[chainnum]={}
                if resnum not in chaintocurrentresnumtothreeletter[chainnum].keys():
                    chaintocurrentresnumtothreeletter[chainnum][resnum]=resname
                    threeletterseq.append(resname)
                    chainseq.append(chainnum)
                    resnumseq.append(resnum)
    return chaintocurrentresnumtothreeletter,threeletterseq,chainseq,resnumseq


def GrabLetterSequence(poltype,seqfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    letterindextolineindex={}
    lettersequence=[]
    listoflettersequence=[]
    temp=open(seqfilename,'r')
    results=temp.readlines()
    temp.close()
    lineidx=0
    lettercount=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            if len(linesplit)==1 and ';' not in line and '>' not in line:
                letters=linesplit[0]
                array=[]
                for e in letters:
                    if e.isalpha():
                        array.append(e)
                lettersequence.extend(array)
                listoflettersequence.append(array)
                for e in array:
                    letterindextolineindex[lettercount]=lineidx
                    lettercount+=1
                lineidx+=1

    return lettersequence,letterindextolineindex,listoflettersequence


def GrabMissingLetterSequence(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    chainlettertonum={'A':0,'B':1,'C':2,'D':3}
    chaintoresnumtomissingthreeletter={}
    chaintomissingresnums={}
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    foundmissresline=False
    for line in results:
        linesplit=line.split()
        if 'M RES C SSSEQI' in line and len(linesplit)==6:
            foundmissresline=True
        if foundmissresline==True and len(linesplit)==5:
            threeletercode=linesplit[2]
            chain=linesplit[-2]
            chainnum=chainlettertonum[chain]
            resnum=int(linesplit[-1])
            if chainnum not in chaintoresnumtomissingthreeletter.keys():
                chaintoresnumtomissingthreeletter[chainnum]={}
                chaintomissingresnums[chainnum]=[]
            chaintoresnumtomissingthreeletter[chainnum][resnum]=threeletercode
            chaintomissingresnums[chainnum].append(resnum)
        elif foundmissresline==True and len(linesplit)!=5:
            if 'M RES C SSSEQI' not in line:
                foundmissresline==False
                break


    return chaintoresnumtomissingthreeletter,chaintomissingresnums


def ConvertThreeLetterToSingleLetter(poltype,chaintoresnumtomissingthreeletter,threelettercodetosinglelettercode):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    chaintoresnumtomissingsingleletter={}
    for chain,resnumtomissingthreeletter in chaintoresnumtomissingthreeletter.items():
        chaintoresnumtomissingsingleletter[chain]={}
        for resnum,threeletter in resnumtomissingthreeletter.items():
            singleletter=threelettercodetosinglelettercode[threeletter]
            chaintoresnumtomissingsingleletter[chain][resnum]=singleletter

    return chaintoresnumtomissingsingleletter



def CombineData(poltype,chaintoresnumtomissingsingleletter,chaintocurrentresnumtosingleletter):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    chaintoresnumtomissing={}
    chaintoresnumtosingleletter={}
    for chainnum,resnumtomissingsingleletter in chaintoresnumtomissingsingleletter.items():
        for resnum,missingsingleletter in resnumtomissingsingleletter.items():
            if chainnum not in chaintoresnumtomissing.keys():
                chaintoresnumtomissing[chainnum]={}
            chaintoresnumtomissing[chainnum][resnum]=True
            if chainnum not in chaintoresnumtosingleletter.keys():
                chaintoresnumtosingleletter[chainnum]={}
            chaintoresnumtosingleletter[chainnum][resnum]=missingsingleletter

    for chainnum,resnumtocurrentsingleletter in chaintocurrentresnumtosingleletter.items():
        for resnum,currentsingleletter in resnumtocurrentsingleletter.items():
            if chainnum not in chaintoresnumtomissing.keys():
                chaintoresnumtomissing[chainnum]={}
            chaintoresnumtomissing[chainnum][resnum]=False
            if chainnum not in chaintoresnumtosingleletter.keys():
                chaintoresnumtosingleletter[chainnum]={}
            chaintoresnumtosingleletter[chainnum][resnum]=currentsingleletter


    return chaintoresnumtomissing,chaintoresnumtosingleletter


def SortByKey(dic):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newdic={}
    for key in sorted(dic):
        newdic[key]=dic[key]
    return newdic



def GenerateGapAndFilledArrays(poltype,chaintoresnumtomissing,chaintoresnumtosingleletter,chainnumtoresnumtolineindex,listoflineindextochainnumresnum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gapresiduearrays=[]
    filledresiduearrays=[]
    lineindextogapresarray={}
    lineindextofilledresarray={}
    for chainnum, resnumtosingleletter in chaintoresnumtosingleletter.items():
        resnums=list(resnumtosingleletter.keys())
        resnumtomissing=chaintoresnumtomissing[chainnum]
        sortedresnums=sorted(resnums)
        gapresarray=[]
        filledresarray=[]
        for resnum in sortedresnums:
            lineindex=chainnumtoresnumtolineindex[chainnum][resnum]
            missing=resnumtomissing[resnum]
            singleletter=resnumtosingleletter[resnum]
            if missing==True:
                gaplet='-'
                filledlet=singleletter
            else:
                gaplet=singleletter
                filledlet=singleletter
            if lineindex not in lineindextogapresarray.keys():
                lineindextogapresarray[lineindex]=[]
                lineindextofilledresarray[lineindex]=[]

            lineindextogapresarray[lineindex].append(gaplet)
            lineindextofilledresarray[lineindex].append(filledlet)
    for lineindex,gapresarray in lineindextogapresarray.items():
        filledresarray=lineindextofilledresarray[lineindex]
        filledresiduearrays.append(filledresarray)
        gapresiduearrays.append(gapresarray)
    return gapresiduearrays,filledresiduearrays


def GenerateAlignmentFile(poltype,chaintoresnumtomissing,chaintoresnumtosingleletter,seqfilename,code,chainnumtoresnumtolineindex,listoflineindextochainnumresnum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gapresiduearrays,filledresiduearrays=GenerateGapAndFilledArrays(poltype,chaintoresnumtomissing,chaintoresnumtosingleletter,chainnumtoresnumtolineindex,listoflineindextochainnumresnum)
    alignmentfilename='alignment.ali'
    temp=open(seqfilename,'r')
    results=temp.readlines()
    temp.close()
    temp=open(alignmentfilename,'w')
    newcode=code+'_filled'
    WriteAlignmentFileComponent(poltype,temp,results,gapresiduearrays,filledresiduearrays,code,newcode,filled=False)
    temp.flush()
    os.fsync(temp.fileno())
    WriteAlignmentFileComponent(poltype,temp,results,gapresiduearrays,filledresiduearrays,code,newcode,filled=True)
    temp.close()
    return alignmentfilename,newcode

def WriteAlignmentFileComponent(poltype,temp,results,gapresiduearrays,filledresiduearrays,code,newcode,filled=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            if len(linesplit)==1 and '>' not in line:
                letters=linesplit[0]
                array=[]
                gapresarray=gapresiduearrays[count]
                filledresarray=filledresiduearrays[count]
                gapstring=''.join(gapresarray)
                filledstring=''.join(filledresarray)
                if count==len(gapresiduearrays)-1:
                    gapstring+='*'
                    filledstring+='*'
                gapstring+='\n'
                filledstring+='\n'
                if filled==False:
                    temp.write(gapstring) 
                else:
                    temp.write(filledstring)
                count+=1
            else:
                if filled==True:
                    if code in line:
                        line=line.replace(code,newcode)

                temp.write(line)


def GenerateLoops(poltype,alignmentfilename,newcode,code):

    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    # Loa
    log.verbose()
    env = Environ()
    env.io.atom_files_directory = ['.', '../atom_files']
    a = LoopModel(env, alnfile = alignmentfilename,knowns = code, sequence = newcode)
    a.starting_model= 1
    a.ending_model  = 1
    
    a.loop.starting_model = 1
    a.loop.ending_model   = 2
    a.loop.md_level       = refine.fast
    
    a.make()


def GrabLastGeneratedPDB(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    files=os.listdir()
    timetofilename={}
    for f in files:
        if '.pdb' in f:
            Ftime=os.path.getmtime(f)
            reltime=time.time()-Ftime
            timetofilename[reltime]=f
    mintime=min(timetofilename.keys())
    minfilename=timetofilename[mintime]
    return minfilename


def FillInMissingResidues(poltype,code):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    # Load the AutoModel class
    from Bio.PDB.PDBList import PDBList

    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(code, file_format='pdb')
    if '.ent' in filename:
        newfilename=filename.replace('.ent','.pdb')
        os.rename(filename,filename.replace('.ent','.pdb'))
        filename=newfilename
    threelettercodetosinglelettercode= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    seqfilename=WriteSEQFile(poltype,filename,code)
    lettersequence,letterindextolineindex,listoflettersequence=GrabLetterSequence(poltype,seqfilename)
    chaintoresnumtomissingthreeletter,chaintomissingresnums=GrabMissingLetterSequence(poltype,filename)
    chaintocurrentresnumtothreeletter,threeletterseq,chainseq,resnumseq=FindCurrentResidueArray(poltype,filename,chaintomissingresnums)

    chaintoresnumtomissingsingleletter=ConvertThreeLetterToSingleLetter(poltype,chaintoresnumtomissingthreeletter,threelettercodetosinglelettercode)
    chaintocurrentresnumtosingleletter=ConvertThreeLetterToSingleLetter(poltype,chaintocurrentresnumtothreeletter,threelettercodetosinglelettercode)
    singleletterseqfrompdb=[threelettercodetosinglelettercode[i] for i in threeletterseq]

    chaintoresnumtomissing,chaintoresnumtosingleletter=CombineData(poltype,chaintoresnumtomissingsingleletter,chaintocurrentresnumtosingleletter)
    chainnumtoresnumtolineindex,listoflineindextochainnumresnum=MapSeqFileToChainRes(poltype,listoflettersequence,chainseq,resnumseq,chaintoresnumtomissing)
    alignmentfilename,newcode=GenerateAlignmentFile(poltype,chaintoresnumtomissing,chaintoresnumtosingleletter,seqfilename,code,chainnumtoresnumtolineindex,listoflineindextochainnumresnum)
    GenerateLoops(poltype,alignmentfilename,newcode,code)
    finalpdb=GrabLastGeneratedPDB(poltype)


def MapSeqFileToChainRes(poltype,listoflettersequence,chainseq,resnumseq,chaintoresnumtomissing):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    listoflineindextochainnumresnum=[]
    chainnumtoresnumtolineindex={}
    count=0
    for lsidx in range(len(listoflettersequence)):
        ls=listoflettersequence[lsidx]
        lineindextochainnumres={}
        for i in range(len(ls)):
            letter=ls[i]
            chainnum=chainseq[count]
            resnum=resnumseq[count]
            if chainnum not in chainnumtoresnumtolineindex.keys():
                chainnumtoresnumtolineindex[chainnum]={}
            chainnumtoresnumtolineindex[chainnum][resnum]=lsidx
            lineindextochainnumres[i]=[chainnum,resnum]
            count+=1
        listoflineindextochainnumresnum.append(lineindextochainnumres)
    for chain,resnumtomissing in chaintoresnumtomissing.items():
        for resnum,missing in resnumtomissing.items():
            if missing==True:
                if resnum not in chainnumtoresnumtolineindex[chain].keys():
                    closestnonmissingresnum=FindClosestResnum(poltype,chainnumtoresnumtolineindex[chain],resnum)
                    lineindex=chainnumtoresnumtolineindex[chain][closestnonmissingresnum]
                    chainnumtoresnumtolineindex[chain][resnum]=lineindex
    return chainnumtoresnumtolineindex,listoflineindextochainnumresnum



def FindClosestResnum(poltype,resnumtolineindex,resnum):
    difftoresnum={}
    for oresnum,lineindex in resnumtolineindex.items():
        diff=oresnum-resnum
        difftoresnum[diff]=oresnum
    mindiff=min(difftoresnum.keys())
    closestnonmissingresnum=difftoresnum[mindiff]

    return closestnonmissingresnum



def CallPDB2PQR(poltype,pdbfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    outputfile=pdbfilename.replace('.pdb','.pqr')
    cmdstr='pdb2pqr30'+' '+pdbfilename+' '+outputfile+' '+'--titration-state-method=propka'
    os.system(cmdstr)
    finaloutputfile=outputfile.replace('.pqr','_final.pdb')
    ConvertPQRToPDB(poltype,outputfile,finaloutputfile)


def ConvertPQRToPDB(poltype,outputfile,finaloutputfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(outputfile)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(mol, outputfile)
    obConversion.SetOutFormat('pdb')
    obConversion.WriteFile(mol,finaloutputfile)



def CallFpocket(poltype,pythonpath,binpath,complexedproteinpdbname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    cmdstr=pythonpath+' '+poltype.pocketscript+' '+'--bindir='+binpath+' '+'--pdbfile='+complexedproteinpdbname
    os.system(cmdstr)


def ReadDictionaryFromFile(jsonfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    with open(jsonfile, 'r') as j:
        dic = json.loads(j.read())
    return dic


def ReadPocketGrids(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    pocketnumtocenter=ReadDictionaryFromFile('pocketnumtocenter.json')
    pocketnumtosize=ReadDictionaryFromFile('pocketnumtosize.json')
    return pocketnumtocenter,pocketnumtosize



def GrabCorrectPocketGrid(ligcenter,pocketnumtocenter,pocketnumtosize):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    difftopocketnum={}
    for pocketnum,center in pocketnumtocenter.items():
        diff=np.linalg.norm(center-ligcenter) 
        difftopocketnum[diff]=pocketnum
    mindiff=min(difftopocketnum.keys())
    minpocketnum=difftopocketnum[mindiff]
    center=np.array(pocketnumtocenter[minpocketnum])
    size=np.array(pocketnumtosize[minpocketnum])
    return center,size



def Grid(center,size):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    half=size/2
    minsize=center-half
    maxsize=center+half
    minx=minsize[0]
    miny=minsize[1]
    minz=minsize[2]
    maxx=maxsize[0]
    maxy=maxsize[1]
    maxz=maxsize[2]

    return minx,maxx,miny,maxy,minz,maxz


def DeterminePocketGrid(poltype,indextocoordinates):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    #pythonpath=poltype.which('python')
    #head,tail=os.path.split(pythonpath)
    #pythonpath=Path(head) 
    #envdir=pythonpath.parent.absolute()
    #envpath=Path(envdir)
    #allenvs=envpath.parent.absolute()
    #oldenvpath=os.path.join(allenvs,poltype.pymolenvname)
    #binpath=os.path.join(oldenvpath,'bin')
    #pythonpath=os.path.join(binpath,'python')
    #CallFpocket(poltype,pythonpath,binpath,poltype.complexedproteinpdbname)
    #pocketnumtocenter,pocketnumtosize=ReadPocketGrids(poltype)
    ligandpdbfilename,receptorpdbfilename=poltype.ExtractLigand(poltype.complexedproteinpdbname,list(indextocoordinates.values()))
    center=poltype.GrabLigandCentroid(ligandpdbfilename)
    #center,size=GrabCorrectPocketGrid(ligcenter,pocketnumtocenter,pocketnumtosize)
    size=ComputeGridAroundLigand(poltype,indextocoordinates,center)
    minx,maxx,miny,maxy,minz,maxz=Grid(center,size)
    return minx,maxx,miny,maxy,minz,maxz


def ComputeGridAroundLigand(poltype,indextocoordinates,size):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    xcol=[]
    ycol=[]
    zcol=[]
    for index,coords in indextocoordinates.items():
        xcol.append(coords[0])
        ycol.append(coords[1])
        zcol.append(coords[2])
    minx=min(xcol)
    maxx=max(xcol)
    miny=min(ycol)
    maxy=max(ycol)
    minz=min(zcol)
    maxz=max(zcol)
    extend=8
    X=maxx-minx+extend
    Y=maxy-miny+extend       
    Z=maxz-minz+extend
    size=np.array([X,Y,Z])

    return size





def SearchInGridIndicesToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextoelement,minx,maxx,miny,maxy,minz,maxz):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for index,coords in nonlighetatmindextocoordinates.items():
        element=nonlighetatmindextoelement[index]
        if element!='H': # for water only restrain oxygens
            if minx<coords[0] and coords[0]<maxx and miny<coords[1] and coords[1]<maxy and minz<coords[2] and coords[2]<maxz:
                poltype.indicestorestrain.append(index)


def FindWatersIonsInPocketToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,indextocoordinates):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.indicestorestrain=[]
    indices=list(indextocoordinates.keys())
    for indices in poltype.allligands:
        newindextocoordinates={}
        for idx in indices:
            coords=indextocoordinates[idx]
            newindextocoordinates[idx]=coords
        minx,maxx,miny,maxy,minz,maxz=DeterminePocketGrid(poltype,newindextocoordinates)
        SearchInGridIndicesToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextoelement,minx,maxx,miny,maxy,minz,maxz)


