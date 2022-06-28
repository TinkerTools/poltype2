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
    proteinindextocoordinates=GenerateUncomplexedProteinPDBFromComplexedPDB(poltype)
    poltype.uncomplexedxyzname=poltype.uncomplexedproteinpdbname.replace('.pdb','.xyz')
    poltype.complexedxyzname=poltype.uncomplexedxyzname.replace('.xyz','_comp.xyz')
    poltype.receptorligandxyzfilename=poltype.complexedxyzname
    poltype.ReadReceptorCharge()
    chainnum=DetectNumberOfChains(poltype,poltype.uncomplexedproteinpdbname)
    resarray=FindCurrentResidueArray(poltype,poltype.uncomplexedproteinpdbname)
    missingresidues=FindMissingResidues(poltype,resarray)

    if not os.path.isfile(poltype.complexedxyzname):
        if chainnum==1:

            cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+poltype.prmfilepath
        else:
            cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+'ALL'+' '+poltype.prmfilepath

        submit.call_subsystem(poltype,cmdstr,wait=True)    
    atoms,coord,order,types,connections=readTXYZ(poltype,poltype.uncomplexedxyzname)
    uncomplexedatomnum=len(proteinindextocoordinates.keys()) # just make sure to move HETATMS between chains after protein xyz
    newuncomplexedatomnum=len(atoms)
    poltype.totalproteinnumber=len(atoms)
    poltype.proteinindices=list(range(1,len(atoms)+1))
    shift= newuncomplexedatomnum- uncomplexedatomnum
    if shift!=0:
        string='WARNING! Missing atoms from original PDB have been added by PDBXYZ. Number of atoms added = '+str(shift)
        warnings.warn(string)
        poltype.WriteToLog(string)
    if len(missingresidues)!=0:
        string="WARNING! Residues are missing "+str(missingresidues)
        warnings.warn(string)
        poltype.WriteToLog(string)

    indextocoordinates,annihilateindextocoordinates,lastligindex,allligands,ligands,listofindextopdbindex=GrabLigandCoordinates(poltype,newuncomplexedatomnum)
    poltype.ligandindices[0]=list(annihilateindextocoordinates.keys())
    poltype.allligandindices[0]=list(indextocoordinates.keys())
    poltype.allligands=allligands
    poltype.ligands=ligands
    nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement=GrabNonLigandHETATMInfo(poltype,poltype.complexedproteinpdbname,indextocoordinates,proteinindextocoordinates,lastligindex) 
    poltype.hetatmindices=list(nonlighetatmindextocoordinates.keys())
    FindWatersIonsInPocketToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,indextocoordinates)
    GenerateComplexedTinkerXYZFile(poltype,poltype.uncomplexedxyzname,indextocoordinates,newuncomplexedatomnum,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,listofindextopdbindex)
    #ligandindices=poltype.ligandindices[0]
    #GeneratePDBFileFromXYZ(poltype,poltype.complexedxyzname,ligandindices)


def GrabNonLigandHETATMInfo(poltype,complexedproteinpdbname,indextocoordinates,proteinindextocoordinates,lastligindex):
    nonlighetatmindextocoordinates={}
    nonlighetatmindextotypenum={}
    nonlighetatmindextoconnectivity={}
    nonlighetatmindextoelement={}
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
    for atom in iteratom:
        index=atom.GetIdx()
        atomicnum=atom.GetAtomicNum()
        coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
        element=an.getElSymbol(atomicnum)
        if coords not in proteinindextocoordinates.values() and coords not in indextocoordinates.values():
            newindex=firstindex+count
            nonlighetatmindextocoordinates[newindex]=coords
            typenum=GrabNonLigandHETATMType(poltype,element)
            nonlighetatmindextotypenum[newindex]=typenum
            shift=newindex-index # check will this work 
            connectivity=GrabConnectivity(poltype,pdbmol,index,shift)
            nonlighetatmindextoconnectivity[newindex]=connectivity
            nonlighetatmindextoelement[newindex]=element
            count+=1
    return nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement


def GrabConnectivity(poltype,pdbmol,index,shift):
    conn=[]
    bonditer=openbabel.OBMolBondIter(pdbmol)
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        bnd=[oendidx,obgnidx]
        if oendidx==index:
            newbgn=obgnidx+shift
            if newbgn not in conn:
                conn.append(newbgn)
        elif obgnidx==index:
            newend=oendidx+shift
            if newend not in conn:
                conn.append(newend)
    return conn




def GrabNonLigandHETATMType(poltype,element): # assumes besides ligand and protein, only water and ions
    if element=='O':
        typenum=poltype.waterOtypenum
    elif element=='H':
        typenum=poltype.waterHtypenum
    else:
        typenum=poltype.elementsymtotinktype[element]

    return typenum


def readTXYZ(poltype,TXYZ):
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
    ligandsmiles=ligandsmiles.replace('=','~').replace('-','~') # sometimes bond order can be different from PDB so remove bond order for detection :)
    p = Chem.MolFromSmarts(ligandsmiles)
    matches=rdkitmol.GetSubstructMatches(p)
    newmatches=[]
    for match in matches:
        newmatch=[i+1 for i in match]
        newmatches.append(newmatch)

    return newmatches,p




def GrabLigandCoordinates(poltype,proteinatomnum): # appending after protein, before pocket water and ions
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
        ligandsmiles=poltype.ligandsmileslist[ligandsmilesidx]
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
    atomiter=openbabel.OBMolAtomIter(pdbmol)
    count=0
    for atom in atomiter:
        atomidx=atom.GetIdx()
        res=atom.GetResidue()
        resname=res.GetName()
        if atomidx in allmatches:
            coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
            newindex=firstatomidx+count
            indextocoordinates[newindex]=coords
            oldindextonewindex[atomidx]=newindex
            if atomidx in truematches:
                annihilateindextocoordinates[newindex]=coords
            count+=1
    if len(indextocoordinates.keys())==0:
        raise ValueError('Complexed PDB missing atoms')
    newallligands=[]
    for match in allligands:
        newallmatch=[]
        for index in match:
            newindex=oldindextonewindex[index]
            newallmatch.append(newindex)
        newallligands.append(newallmatch)

    newligands=[]
    for match in ligands:
        newmatch=[]
        for index in match:
            newindex=oldindextonewindex[index]
            newmatch.append(newindex)
        newligands.append(newmatch)

    return indextocoordinates,annihilateindextocoordinates,newindex,newallligands,newligands,listofindextopdbindex   


def GenerateComplexedTinkerXYZFile(poltype,uncomplexedxyzname,indextocoordinates,uncomplexedatomnum,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,listofindextopdbindex):
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
        lastindex=WriteLigandToXYZ(poltype,temp,atoms,order,lastindex,indextocoordinates,types,connections,indextopdbindex)
    WriteHETATMToXYZ(poltype,temp,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement)
    temp.close()

def WriteLigandToXYZ(poltype,temp,atoms,order,lastindex,indextocoordinates,types,connections,indextopdbindex):
    shift=lastindex
    for idx in range(len(atoms)):
        element=atoms[idx]
        oldindex=int(order[idx])
        pdbindex=indextopdbindex[oldindex]
        newindex=oldindex+shift
        coords=indextocoordinates[pdbindex]
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
    uncomplexedatomindices=[]
    indextocoordinates={}
    temp=open(poltype.complexedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            resname=line[17:21].strip()
            linesplit=line.split()
            index=int(linesplit[1])
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
    indexestodelete.sort(reverse=True)
    for idx in indexestodelete:
        atom=pdbmol.GetAtom(idx)
        pdbmol.DeleteAtom(atom)
    
    obConversion.SetOutFormat('pdb')
    poltype.uncomplexedproteinpdbname='uncomplexed.pdb'
    obConversion.WriteFile(pdbmol,poltype.uncomplexedproteinpdbname)
    return indextocoordinates


def GeneratePDBFileFromXYZ(poltype,xyzfile,ligandindices):
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
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    # Loa
    e = Environ()
    m = Model(e, file=filename)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    seqfilename=filename.replace('.pdb','.seq')
    aln.write(file=filename.replace('.pdb','.seq'))
    return seqfilename


def FindCurrentResidueArray(poltype,filename):
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    resarray=[]
    for line in results:
        linesplit=line.split()
        if 'ATOM' in line and line[:4]=='ATOM':
            resnum=int(line[22:26])
            if resnum not in resarray:
                resarray.append(resnum)
    return resarray


def FindMissingResidues(poltype,resarray):
    missingresidues=[]
    firstres=1 # dont start at first residue found, start at 1
    lastres=resarray[-1]
    allres=list(range(firstres,lastres+1))
    for res in allres:
        if res not in resarray:
            missingresidues.append(res)



    return missingresidues



def GrabLetterSequence(poltype,seqfilename):
    lettersequence=[]
    chainsequence=[]
    temp=open(seqfilename,'r')
    results=temp.readlines()
    temp.close()
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            if len(linesplit)==1:
                letters=linesplit[0]
                array=[]
                allalpha=True
                for e in letters:
                    if e.isalpha():
                        array.append(e)
                    else:
                        if e!='*':
                            allalpha=False
                            break
                if allalpha==True:
                    lettersequence.extend(array)
                    chainseq=[count]*len(array)
                    chainsequence.extend(chainseq)
                    count+=1



    return lettersequence,chainsequence


def GrabMissingLetterSequence(poltype,filename):
    chainlettertonum={'A':0,'B':1,'C':2,'D':3}
    resnumtomissingthreeletter={}
    resnumtomissingchainseq={}
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    foundmissresline=False
    for line in results:
        linesplit=line.split()
        if 'M RES C SSSEQI' in line:
            foundmissresline=True
        if foundmissresline==True and len(linesplit)==5:
            threeletercode=linesplit[2]
            chain=linesplit[-2]
            chainnum=chainlettertonum[chain]
            resnum=int(linesplit[-1])
            resnumtomissingthreeletter[resnum]=threeletercode
            resnumtomissingchainseq[resnum]=chainnum
        elif foundmissresline==True and len(linesplit)!=5:
            if 'M RES C SSSEQI' not in line:
                foundmissresline==False
                break


    return resnumtomissingthreeletter,resnumtomissingchainseq


def ConvertThreeLetterToSingleLetter(poltype,resnumtomissingthreeletter,threelettercodetosinglelettercode):
    resnumtomissingsingleletter={}
    for resnum,threeletter in resnumtomissingthreeletter.items():
        singleletter=threelettercodetosinglelettercode[threeletter]
        resnumtomissingsingleletter[resnum]=singleletter
    return resnumtomissingsingleletter



def CombineData(poltype,resnumtomissingsingleletter,resnumtocurrentsingleletter,resnumtomissingchainseq,resnumtocurrentchainseq):
    resnumtochainseq={}
    resnumtomissing={}
    resnumtosingleletter={}
    for resnum,missingsingleletter in resnumtomissingsingleletter.items():
        resnumtomissing[resnum]=True
        resnumtosingleletter[resnum]=missingsingleletter
        resnumtochainseq[resnum]=resnumtomissingchainseq[resnum]

    for resnum,currentsingleletter in resnumtocurrentsingleletter.items():
        resnumtomissing[resnum]=False
        resnumtosingleletter[resnum]=currentsingleletter
        chainseq=resnumtocurrentchainseq[resnum]
        resnumtochainseq[resnum]=chainseq
            


    return SortByKey(resnumtochainseq),SortByKey(resnumtomissing),SortByKey(resnumtosingleletter)


def SortByKey(dic):
    newdic={}
    for key in sorted(dic):
        newdic[key]=dic[key]
    return newdic



def GenerateGapAndFilledArrays(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter):
    gapresiduearrays=[]
    filledresiduearrays=[]
    chainseqtoresnums={}
    for resnum,chainseq in resnumtochainseq.items():
        if chainseq not in chainseqtoresnums.keys():
            chainseqtoresnums[chainseq]=[]
        if resnum not in chainseqtoresnums[chainseq]:
            chainseqtoresnums[chainseq].append(resnum)

    for chainseq,resnums in chainseqtoresnums.items():
        sortedresnums=sorted(resnums)
        gapresarray=[]
        filledresarray=[]
        for resnum in sortedresnums:
            missing=resnumtomissing[resnum]
            singleletter=resnumtosingleletter[resnum]
            if missing==True:
                gaplet='-'
                filledlet=singleletter
            else:
                gaplet=singleletter
                filledlet=singleletter
            gapresarray.append(gaplet)
            filledresarray.append(filledlet)
        gapresiduearrays.append(gapresarray)
        filledresiduearrays.append(filledresarray)
            

    return gapresiduearrays,filledresiduearrays


def GenerateAlignmentFile(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter,seqfilename,code):
    gapresiduearrays,filledresiduearrays=GenerateGapAndFilledArrays(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter)
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
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            if len(linesplit)==1:
                letters=linesplit[0]
                array=[]
                allalpha=True
                for e in letters:
                    if e.isalpha():
                        array.append(e)
                    else:
                        if e!='*':
                            allalpha=False
                            break
                if allalpha==True:
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
            else:
                temp.write(line)


def GenerateLoops(poltype,alignmentfilename,newcode,code):
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
    resarray=FindCurrentResidueArray(poltype,filename)
    missingresidues=FindMissingResidues(poltype,resarray)
    lettersequence,chainsequence=GrabLetterSequence(poltype,seqfilename)
    resnumtomissingthreeletter,resnumtomissingchainseq=GrabMissingLetterSequence(poltype,filename)
    resnumtocurrentsingleletter=dict(zip(resarray,lettersequence))
    resnumtocurrentchainseq=dict(zip(resarray,chainsequence))
    resnumtomissingsingleletter=ConvertThreeLetterToSingleLetter(poltype,resnumtomissingthreeletter,threelettercodetosinglelettercode)
    resnumtochainseq,resnumtomissing,resnumtosingleletter=CombineData(poltype,resnumtomissingsingleletter,resnumtocurrentsingleletter,resnumtomissingchainseq,resnumtocurrentchainseq)
    alignmentfilename,newcode=GenerateAlignmentFile(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter,seqfilename,code)
    GenerateLoops(poltype,alignmentfilename,newcode,code)
    finalpdb=GrabLastGeneratedPDB(poltype)



def CallPDB2PQR(poltype,pdbfilename):
    outputfile=pdbfilename.replace('.pdb','.pqr')
    cmdstr='pdb2pqr30'+' '+pdbfilename+' '+outputfile+' '+'--titration-state-method=propka'
    os.system(cmdstr)
    finaloutputfile=outputfile.replace('.pqr','_final.pdb')
    ConvertPQRToPDB(poltype,outputfile,finaloutputfile)


def ConvertPQRToPDB(poltype,outputfile,finaloutputfile):
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(outputfile)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(mol, outputfile)
    obConversion.SetOutFormat('pdb')
    obConversion.WriteFile(mol,finaloutputfile)



def CallFpocket(poltype,pythonpath,binpath,complexedproteinpdbname):
    cmdstr=pythonpath+' '+poltype.pocketscript+' '+'--bindir='+binpath+' '+'--pdbfile='+complexedproteinpdbname
    os.system(cmdstr)


def ReadDictionaryFromFile(jsonfile):
    with open(jsonfile, 'r') as j:
        dic = json.loads(j.read())
    return dic


def ReadPocketGrids(poltype):
    pocketnumtocenter=ReadDictionaryFromFile('pocketnumtocenter.json')
    pocketnumtosize=ReadDictionaryFromFile('pocketnumtosize.json')
    return pocketnumtocenter,pocketnumtosize



def GrabCorrectPocketGrid(ligcenter,pocketnumtocenter,pocketnumtosize):
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
    for index,coords in nonlighetatmindextocoordinates.items():
        element=nonlighetatmindextoelement[index]
        if element!='H': # for water only restrain oxygens
            if minx<coords[0] and coords[0]<maxx and miny<coords[1] and coords[1]<maxy and minz<coords[2] and coords[2]<maxz:
                poltype.indicestorestrain.append(index)


def FindWatersIonsInPocketToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextotypenum,nonlighetatmindextoconnectivity,nonlighetatmindextoelement,indextocoordinates):
    poltype.indicestorestrain=[]
    indices=list(indextocoordinates.keys())
    for indices in poltype.allligands:
        newindextocoordinates={}
        for idx in indices:
            coords=indextocoordinates[idx]
            newindextocoordinates[idx]=coords
        minx,maxx,miny,maxy,minz,maxz=DeterminePocketGrid(poltype,newindextocoordinates)
        SearchInGridIndicesToRestrain(poltype,nonlighetatmindextocoordinates,nonlighetatmindextoelement,minx,maxx,miny,maxy,minz,maxz)


