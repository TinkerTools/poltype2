import electrostaticpotential as esp
import torsiongenerator as torgen
import optimization as opt
import apicall as call
from socket import gethostname

import os
import sys
import numpy
import time
import openbabel
from rdkit import Chem
from rdkit.Chem import rdmolfiles,rdMolDescriptors
import shutil
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import svgutils.transform as sg
from cairosvg import svg2png
import copy
from collections import defaultdict
from rdkit.Chem import rdDepictor
import matplotlib.pyplot as plt
from os.path import dirname, abspath      
from itertools import combinations
from rdkit.Chem import rdFMCS
import json
from collections import Counter
from itertools import combinations
import re
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDistGeom


def GrabTorsionParametersFromFragments(poltype,torlist,rotbndindextofragmentfilepath):
    valenceprmlist=[]
    symmtorlist=[]
    for tor in torlist:
        symmtorlist.append(torgen.get_class_key(poltype,tor[0],tor[1],tor[2],tor[3]))
    torprmdic={}
    curdir=os.getcwd()
    for rotbndindex,fragmentfilepath in rotbndindextofragmentfilepath.items():
        path,filename=os.path.split(fragmentfilepath)
        os.chdir(path)
        filelist=os.listdir(os.getcwd())
        for ff in filelist:
            if '.prm' in ff:
                temp=open(ff,'r')
                results=temp.readlines()
                temp.close()
                for line in results:
                    valenceprmlist.append(line)
            if '.key_6' in ff:
                temp=open(ff,'r')
                results=temp.readlines()
                temp.close()
                for line in results:
                    if 'torsion' in line:
                        linesplit=line.split()
                        typea=int(linesplit[1])
                        typeb=int(linesplit[2])
                        typec=int(linesplit[3])
                        typed=int(linesplit[4])
                        tor=[typea,typeb,typec,typed]
                        torkey='%d %d %d %d' % (typea, typeb, typec, typed)
                        if torkey in symmtorlist:
                            torprmdic[torkey]=line
                
    os.chdir(curdir)
    temp=open(poltype.key4fname,'r')
    results=temp.readlines()
    temp.close()
    temp=open(poltype.key5fname,'w')
    for line in results:
        if 'torsion' in line:
            linesplit=line.split()
            typea=int(linesplit[1])
            typeb=int(linesplit[2])
            typec=int(linesplit[3])
            typed=int(linesplit[4])
            tor=[typea,typeb,typec,typed]
            torkey='%d %d %d %d' % (typea, typeb, typec, typed)
            if torkey in torprmdic.keys():
                torline=torprmdic[torkey]
                temp.write(torline)
            else:
                temp.write(line)
        else:
            temp.write(line)
    temp.close()
    temp=open('valence.prms','w')
    for line in valenceprmlist:
        temp.write(line)
    temp.close()

def GrabWBOMatrixGaussian(poltype,outputlog,mol):
    try:
        WBOmatrix=numpy.empty((mol.GetNumAtoms(),mol.GetNumAtoms()))
    except:
        WBOmatrix=numpy.empty((mol.NumAtoms(),mol.NumAtoms()))
    temp=open(outputlog,'r')
    results=temp.readlines()
    temp.close()
    juststartWBOmatrix=False
    currentcolnum=0
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if 'Wiberg bond index matrix' in line:
            juststartWBOmatrix=True    
        elif 'Atom' in line and juststartWBOmatrix==True:
            matcols=len(linesplit)-1
        elif 'Wiberg bond index, Totals by atom' in line and juststartWBOmatrix==True:
            return WBOmatrix
        elif line=='\n' and juststartWBOmatrix==True:
            if 'Wiberg bond index matrix' not in results[lineidx-1]:
                currentcolnum+=matcols
        elif juststartWBOmatrix==True and 'Atom' not in line and line!='\n' and '--' not in line:
            rownum=int(linesplit[0].replace('.',''))
            ele=linesplit[1]
            wborowvalues=linesplit[2:]
            wborowvalues=[float(i) for i in wborowvalues]
            for i in range(len(wborowvalues)):
                colnum=i+1+currentcolnum
                value=wborowvalues[i]	
                WBOmatrix[rownum-1,colnum-1]=float(value)
    return WBOmatrix
                
def GrabWBOMatrixPsi4(poltype,outputlog,molecule):
    try:
        WBOmatrix=numpy.empty((molecule.GetNumAtoms(),molecule.GetNumAtoms()))
    except:
        WBOmatrix=numpy.empty((molecule.NumAtoms(),molecule.NumAtoms()))

    temp=open(outputlog,'r')
    results=temp.readlines()
    temp.close()
    juststartWBOmatrix=False
    currentcolnum=0
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if 'Wiberg Bond Indices' in line:
            juststartWBOmatrix=True    
        elif 'Atomic Valences:' in line and juststartWBOmatrix==True:
            return WBOmatrix
        elif AllIntegers(poltype,line.split())==True and juststartWBOmatrix==True and line!='\n':
            colrowindex=lineidx
        elif juststartWBOmatrix==True and 'Irrep:' not in line and line!='\n' and AllIntegers(poltype,line.split())==False:
            row=line.split()[1:]
            colindexrow=results[colrowindex].split()
            rowindex=int(line.split()[0])
            for i in range(len(row)):
                value=float(row[i])
                colindex=int(colindexrow[i])
                WBOmatrix[rowindex-1,colindex-1]=value
    return WBOmatrix
                           

 
def AllIntegers(poltype,testlist):
    allintegers=True
    for value in testlist:
        if not value.isdigit():
            allintegers=False
    return allintegers

def FindEquivalentFragments(poltype,fragmentarray):
    equivalentfragmentsarray=[]
    smartsarray=[rdmolfiles.MolToSmarts(m) for m in fragmentarray]
    repeatedvalues=(Counter(smartsarray) - Counter(set(smartsarray))).keys()
    for smart in repeatedvalues:
        repeatedindexes=getIndexPositions(poltype,smartsarray,smarts)
        temp=[]
        for idx in repeatedindexes:
            temp.append(fragmentarray[idx])
        equivalentfragmentsarray.append(temp)
    return equivalentfragmentsarray

def getIndexPositions(poltype,listOfElements, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    indexPosList = []
    indexPos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            indexPos = listOfElements.index(element, indexPos)
            # Add the index position in list
            indexPosList.append(indexPos)
            indexPos += 1
        except ValueError as e:
            break
 
    return indexPosList
         
def FindEquivalentRotatableBonds(poltype,equivalentfragmentsarray,rotbndindextofragment):
    equivalentrotbndindexarrays=[]
    
    for array in equivalentfragmentsarray:
        temp=[]
        for fragmol in array:
            rotbndindex=FindRotatableBond(poltype,fragmol,rotbndindextofragment,temp)
            if rotbndindex not in temp:
                temp.append(rotbndindex)
        if len(temp)!=0:
            equivalentrotbndindexarrays.append(temp)
    return equivalentrotbndindexarrays

def FindRotatableBond(poltype,fragmol,rotbndindextofragment,temp):
    for rotbndindex in rotbndindextofragment.keys():
        m=rotbndindextofragment[rotbndindex]
        if len(m.GetAtoms())==len(fragmol.GetAtoms()) and rotbndindex not in temp:
            return rotbndindex

def FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog):
    poltypeinput={'suppressdipoleerr':'True','optmethod':poltype.optmethod,'toroptmethod':poltype.toroptmethod,'espmethod':poltype.espmethod,'torspmethod':poltype.torspmethod,'dmamethod':poltype.dmamethod,'torspbasisset':poltype.torspbasisset,'espbasisset':poltype.espbasisset,'dmabasisset':poltype.dmabasisset,'toroptbasisset':poltype.toroptbasisset,'optbasisset':poltype.optbasisset,'onlyrotbndslist':strfragrotbndindexes,'bashrcpath':poltype.bashrcpath,'externalapi':poltype.externalapi,'use_gaus':poltype.use_gaus,'use_gausoptonly':poltype.use_gausoptonly,'isfragjob':True,'poltypepath':poltype.poltypepath,'structure':tail,'numproc':poltype.numproc,'maxmem':poltype.maxmem,'maxdisk':poltype.maxdisk,'printoutput':True}
    inifilepath=poltype.WritePoltypeInitializationFile(poltypeinput)
    cmdstr='nohup'+' '+'python'+' '+poltype.poltypepath+r'/'+'poltype.py'+' '+'&'
    cmdstr='cd '+os.getcwd()+' && '+cmdstr
    molecprefix =  os.path.splitext(tail)[0]
    logname = molecprefix+ "-poltype.log"
    listofjobs.append(cmdstr)
    logpath=os.getcwd()+r'/'+logname
    if os.path.isfile(logpath): # make sure to remove logfile if exists, dont want WaitForTermination to catch previous errors before job is resubmitted
        os.remove(logpath)
    jobtooutputlog[cmdstr]=logpath    
    return listofjobs,jobtooutputlog,logpath

def SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog):
    if poltype.externalapi!=None:
        finishedjobs,errorjobs=poltype.CallJobsLocalHost(jobtooutputlog,True)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False)

    return finishedjobs,errorjobs

 
def SpawnPoltypeJobsForFragments(poltype,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,torlist,equivalentfragmentsarray,equivalentrotbndindexarrays):
    parentdir=dirname(abspath(os.getcwd()))
    listofjobs=[]
    jobtooutputlog={}
    logtoconvertidxs={}
    if equivalentfragmentsarray==[]:
        for rotbndindex in rotbndindextofragment.keys():
            fragmol=rotbndindextofragment[rotbndindex]
            fragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]
            head,tail=os.path.split(fragmentfilepath)
            os.chdir(head)
            parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
            rotbndindexes=rotbndindex.split('_')
            parentrotbndindexes=[int(i) for i in rotbndindexes]
            rotbndindexes=[int(i)-1 for i in parentrotbndindexes]
            fragrotbndindexes=[parentindextofragindex[i]+1 for i in rotbndindexes]
            wholexyz=parentdir+r'/'+poltype.xyzoutfile
            wholemol=parentdir+r'/'+poltype.molstructfname,
            strfragrotbndindexes=str(fragrotbndindexes[0])+' '+str(fragrotbndindexes[1])
            strparentrotbndindexes=str(parentrotbndindexes[0])+' '+str(parentrotbndindexes[1])
            MakeTorsionFileName(poltype,strparentrotbndindexes)
            WriteDictionaryToFile(poltype,parentindextofragindex)
            listofjobs,jobtooutputlog,newlog=FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog)
            logtoconvertidxs[newlog]=[fragmol,tail,wholexyz,wholemol,head]
    else:
        nonequivalentrotbndidxs=[] 
        for array in equivalentrotbndindexarrays:
            strfragrotbndindexes=''
            strparentrotbndindexes=''
            for i in range(len(array)):
                rotbndindex=array[i]
                if i==0:
                    equivalentrotbndindex=rotbndindex
                else:
                    nonequivalentrotbndidxs.append(rotbndindex)
                parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
                rotbndindexes=rotbndindex.split('_')
                parentrotbndindexes=[int(i) for i in rotbndindexes]
                rotbndindexes=[int(i)-1 for i in parentrotbndindexes]
                fragrotbndindexes=[parentindextofragindex[i]+1 for i in rotbndindexes]
                strfragrotbndindexes+=str(fragrotbndindexes[0])+' '+str(fragrotbndindexes[1])+','
                strparentrotbndindexes+=str(parentrotbndindexes[0])+' '+str(parentrotbndindexes[1])+','
            strfragrotbndindexes=strfragrotbndindexes[:-1]
            strparentrotbndindexes=strparentrotbndindexes[:-1]
            fragmol=rotbndindextofragment[equivalentrotbndindex]
            fragmentfilepath=rotbndindextofragmentfilepath[equivalentrotbndindex]
            head,tail=os.path.split(fragmentfilepath)
            os.chdir(head)
            MakeTorsionFileName(poltype,strparentrotbndindexes)
            parentindextofragindex=rotbndindextoparentindextofragindex[equivalentrotbndindex]
            WriteDictionaryToFile(poltype,parentindextofragindex)
            wholexyz=parentdir+r'/'+poltype.xyzoutfile
            wholemol=parentdir+r'/'+poltype.molstructfname
            listofjobs,jobtooutputlog,newlog=FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog)
            logtoconvertidxs[newlog]=[fragmol,tail,wholexyz,wholemol,head]
        for rotbndindex in rotbndindextofragment.keys():
            if rotbndindex not in nonequivalentrotbndidxs:
                fragmol=rotbndindextofragment[rotbndindex]
                fragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]
                head,tail=os.path.split(fragmentfilepath)
                os.chdir(head)
                parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
                rotbndindexes=rotbndindex.split('_')
                rotbndindexes=[int(i)-1 for i in rotbndindexes]
                fragrotbndindexes=[parentindextofragindex[i]+1 for i in rotbndindexes]
                wholexyz=parentdir+r'/'+poltype.xyzoutfile
                wholemol=parentdir+r'/'+poltype.molstructfname,
                strfragrotbndindexes=str(fragrotbndindexes[0])+' '+str(fragrotbndindexes[1])
                listofjobs,jobtooutputlog,newlog=FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog)
                logtoconvertidxs[newlog]=[fragmol,tail,wholexyz,wholemol,head]
    os.chdir(parentdir)
    finishedjobs,errorjobs=SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog)
    for log in finishedjobs:
        [fragmol,tail,wholexyz,wholemol,filepath]=logtoconvertidxs[log]
        ConvertFragIdxToWholeIdx(poltype,torlist,rotbndindextoparentindextofragindex,fragmol,tail,wholexyz,wholemol,filepath)
    
def MakeTorsionFileName(poltype,string):
    temp=open('torsions.txt','w')
    temp.write(string+'\n')
    temp.close()


def WriteDictionaryToFile(poltype,dictionary):
    newdic={}
    for key,value in dictionary.items():
        newdic[int(key)]=int(value)
    json.dump(newdic, open("parentindextofragindex.txt",'w'))


def DeleteEquivalentFragments(poltype,equivalentfragmentstodelete):
    curdir=os.getcwd()
    os.chdir('..')
    for fold in equivalentfragmentstodelete:
        if os.path.isdir(fold):
            shutil.rmtree(fold)
    os.chdir(curdir)

def get_class_key(poltype,a, b, c, d,fragidxtotypeidx):
    """
    Intent: Given a set of atom idx's, return the class key for the set (the class numbers of the atoms appended together)
    """
    cla = fragidxtotypeidx[a]
    clb = fragidxtotypeidx[b]
    clc = fragidxtotypeidx[c]
    cld = fragidxtotypeidx[d]

    if ((clb > clc) or (clb == clc and cla > cld)):
        return '%d %d %d %d' % (cld, clc, clb, cla)
    return '%d %d %d %d' % (cla, clb, clc, cld)



def ConvertFragIdxToWholeIdx(poltype,torlist,rotbndindextoparentindextofragindex,fragmol,fragmentfilename,wholexyz,wholemol,filepath):
    fragmentfileprefix=fragmentfilename.replace('.sdf','')
    if not os.path.isfile(filepath+r'/'+fragmentfileprefix+'.key_5'):
        return # if POLTYPE job failed, just try to submit the other fragment jobs instead of killing parent job
    temp=open(filepath+r'/'+fragmentfileprefix+'.key_5','r')
    fragkeyresults=temp.readlines()
    temp.close()
    temp=open(wholexyz,'r')
    wholetttxyzresults=temp.readlines()
    temp.close()
    temp=open(filepath+r'/'+'ttt.xyz','r')
    fragtttxyzresults=temp.readlines()
    temp.close()

    fragidxtotypeidx={}
    for lineidx in range(len(fragtttxyzresults)):
        line=fragtttxyzresults[lineidx]
        linesplit=line.split()
        if lineidx!=0:
            idx=int(linesplit[0])
            typeidx=int(linesplit[5])
            fragidxtotypeidx[idx]=typeidx


    wholeidxtotypeidx={}
    for lineidx in range(len(wholetttxyzresults)):
        line=wholetttxyzresults[lineidx]
        linesplit=line.split()
        if lineidx!=0:
            idx=int(linesplit[0])
            typeidx=int(linesplit[5])
            wholeidxtotypeidx[idx]=typeidx

    wholemolidxtofragidx=json.load(open(filepath+r'/'+"parentindextofragindex.txt"))
    temp={} # convert rdkit to babel 
    for key,value in wholemolidxtofragidx.items():
        temp[int(key)+1]=int(value)+1
    fragidxtowholemolidx={v: k for k, v in temp.items()}
    fragsmarts=rdmolfiles.MolToSmarts(fragmol)
    m=mol_with_atom_index(poltype,fragmol)
    fragsmirks=rdmolfiles.MolToSmarts(m) 
    fragidxarray=GrabAtomOrder(poltype,fragsmirks)
    classkeytosmilesposarray={}
    for tor in torlist:
        rotbndindex=str(tor[1])+'_'+str(tor[2])
        parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
        smilesposarray=[]
        fragtor=[]
        for index in tor:
            fragindex=parentindextofragindex[index-1]+1
            fragtor.append(fragindex)
            fragidxarraypos=fragidxarray.index(fragindex)
            smilespos=fragidxarraypos+1
            smilesposarray.append(smilespos)

        classkey=get_class_key(poltype,fragtor[0],fragtor[1],fragtor[2],fragtor[3],fragidxtotypeidx)
        classkeytosmilesposarray[classkey]=smilesposarray
    p = Chem.MolFromSmarts(fragsmarts)
    diditmatchrdkit=fragmol.HasSubstructMatch(p)

    fragtypeidxtowholemoltypeidx={}
    for fragidx in fragidxtowholemolidx.keys():
        wholemolidx=fragidxtowholemolidx[fragidx]
        fragtypeidx=fragidxtotypeidx[fragidx]
        wholemoltypeidx=wholeidxtotypeidx[wholemolidx]
        fragtypeidxtowholemoltypeidx[fragtypeidx]=wholemoltypeidx


    newtemp=open(filepath+r'/'+'valence.prms','w')
    temp=open(filepath+r'/'+fragmentfileprefix+'.key_6','w')
    for line in fragkeyresults:
        if 'torsion' in line:
            linesplit=line.split()
            typea=int(linesplit[1])
            typeb=int(linesplit[2])
            typec=int(linesplit[3])
            typed=int(linesplit[4])
            tor=[typea,typeb,typec,typed]
            torkey='%d %d %d %d' % (typea, typeb, typec, typed)
            if torkey in classkeytosmilesposarray.keys():
                smilesposarray=classkeytosmilesposarray[torkey]
                wholetypea=fragtypeidxtowholemoltypeidx[typea]
                wholetypeb=fragtypeidxtowholemoltypeidx[typeb]
                wholetypec=fragtypeidxtowholemoltypeidx[typec]
                wholetyped=fragtypeidxtowholemoltypeidx[typed]
                wholemoltypestring=str(wholetypea)+' '+str(wholetypeb)+' '+str(wholetypec)+' '+str(wholetyped)
                torprmstring=' '.join(linesplit[5:])
                torstring='torsion '+wholemoltypestring+' '+torprmstring
                temp.write(torstring+'\n')
                valencestring='"'+fragsmarts+'"'+' '+':'+' '+'['+str(smilesposarray[0])+','+str(smilesposarray[1])+','+str(smilesposarray[2])+','+str(smilesposarray[3])+','
                torprmlist=linesplit[5:]
                prms=torprmlist[0::3]
                for prm in prms:
                    valencestring+=prm+','
                valencestring=valencestring[:-1]
                valencestring+=']'+','+' '+r'\\'
                newtemp.write(valencestring+'\n')
        else:
            temp.write(line)    
    temp.close()
    newtemp.close()
    return

def GrabAtomOrder(poltype,smirks):
    atomorder=[]
    for i in range(len(smirks)):
        e=smirks[i]
        prevchar=smirks[i-1]
        try:
            nextchar=smirks[i+1]
        except:
            break
        if prevchar==':' and e.isdigit() and nextchar!='-' and nextchar!=')' and nextchar!=':' and nextchar!='=':
            atomindex=GrabAtomIndex(poltype,i,smirks)
            atomorder.append(atomindex)
    return atomorder


def GrabAtomIndex(poltype,i,smirks):
    num=[]
    for j in range(i,len(smirks)):
        char=smirks[j]
        if char.isdigit():
            num.append(char)
        if char==']':
            break
    atomindex=int(''.join(num))
    return atomindex 


def GenerateFrag(poltype,molindexlist,valtopo,mol):
    em = Chem.EditableMol(Chem.Mol())
    oldindextonewindex={}
    for i,idx in enumerate(molindexlist):
        em.AddAtom(poltype.rdkitmol.GetAtomWithIdx(idx))
        oldindextonewindex[idx]=i
    for bond in poltype.rdkitmol.GetBonds():
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        if oendidx in oldindextonewindex.keys():
            endidx=oldindextonewindex[oendidx]
        else:
            continue
        if obgnidx in oldindextonewindex.keys():
            bgnidx=oldindextonewindex[obgnidx]
        else:
            continue
        em.AddBond(bgnidx,endidx,bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    #Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap=oldindextonewindex
    newindextooldindex={v: k for k, v in oldindextonewindex.items()}
    fragidxswithchangedval=[]
    fragidxtooriginalformalcharge={}
    for atom in res.GetAtoms():
        atomidx=atom.GetIdx()
        oldindex=newindextooldindex[atomidx]
        babelidx=oldindex+1
        parentatom=poltype.rdkitmol.GetAtomWithIdx(oldindex)
        parentatomidx=parentatom.GetIdx()
        diff=0
        for nparentatom in parentatom.GetNeighbors():
            nparentidx=nparentatom.GetIdx()
            if nparentidx in oldindextonewindex.keys(): # then skip this (because it has not been removed)
                continue
            else:
                bond=poltype.rdkitmol.GetBondBetweenAtoms(parentatomidx,nparentidx)
                bondorder=bond.GetBondTypeAsDouble()
                diff-=bondorder
        parentatomformalcharge=parentatom.GetFormalCharge()
        newformalcharge=int(parentatomformalcharge+diff)
        fragatomformalcharge=atom.GetFormalCharge()
        atom.SetFormalCharge(newformalcharge)
        fragidxtooriginalformalcharge[atomidx]=parentatomformalcharge
        if diff!=0:
            fragidxswithchangedval.append(atomidx)
            #if atom.GetIsAromatic()==True:
            #    atom.SetIsAromatic(False)

    return res,oldindextonewindex,fragidxswithchangedval,fragidxtooriginalformalcharge

def GenerateBondTopology(poltype,mol):
    bondtopology=[]
    bondsettobondorder={}
    iterbond=openbabel.OBMolBondIter(mol) # iterator for all bond objects in the molecule
    for bond in iterbond:
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        aidx=a.GetIdx()
        bidx=b.GetIdx()
        bondset=set([aidx,bidx])
        bondtopology.append(bondset)
        bondorder=bond.GetBondOrder()
        bondsettobondorder[tuple(bondset)]=bondorder

    return bondtopology,bondsettobondorder

def GenerateValenceTopology(poltype,mol):
    valtopo={}
    iteratom = openbabel.OBMolAtomIter(mol)
    for atom in iteratom:
        val=atom.GetValence()
        valtopo[atom.GetIdx()]=val
    return valtopo

 
def HydrateMolecule(poltype,mol,basename,bondtopo,bondsettobondorder,parentidxtofragidx,parentmol,fragidxswithchangedval,fragidxtooriginalformalcharge):
    for atom in mol.GetAtoms():
        atomidx=atom.GetIdx()
        originalcharge=fragidxtooriginalformalcharge[atomidx]
        atom.SetFormalCharge(originalcharge) 
    hydmol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(hydmol)
    return hydmol

def ConvertRdkitMolToOBMol(poltype,mol):
    outputname='intermediate.mol'
    WriteRdkitMolToMolFile(poltype,mol,outputname)
    OBmol=ReadMolFileToOBMol(poltype,outputname)
    return OBmol

def WriteOBMolToSDF(poltype,mol,outputname):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('sdf')
    tmpconv.WriteFile(mol,outputname)


def WriteOBMolToXYZ(poltype,mol,outputname):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('xyz')
    tmpconv.WriteFile(mol,outputname)


def WriteOBMolToMol(poltype,mol,outputname):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('mol')
    tmpconv.WriteFile(mol,outputname)

def WriteRdkitMolToMolFile(poltype,mol,outputname):
    rdmolfiles.MolToMolFile(mol,outputname,kekulize=True)

def ReadRdkitMolFromMolFile(poltype,inputname):
    rdkitmol=rdmolfiles.MolFromMolFile(inputname,sanitize=False)
    return rdkitmol

def ReadMolFileToOBMol(poltype,filename):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('mol')
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,filename)
    return fragmolbabel

def ReadToOBMol(poltype,filename):
    tmpconv = openbabel.OBConversion()
    inFormat = tmpconv.FormatFromExt(filename)
    tmpconv.SetInFormat(inFormat)
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,filename)
    return fragmolbabel

     

def mol_with_atom_index(poltype,mol):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx()+1 ) )
    return mol

def mol_with_atom_index_removed(poltype,mol):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        atom=mol.GetAtomWithIdx(idx)
        atom.ClearProp('molAtomMapNumber')
    return mol



def GenerateWBOMatrix(poltype,molecule,structfname,fragmoltocharge):
    error=False
    WBOmatrix=None
    curespmethod=poltype.espmethod
    curspbasisset=poltype.espbasisset
    poltype.espmethod='HF'
    poltype.espbasisset='MINIX'
    charge=fragmoltocharge[molecule]
    inputname,outputname=esp.CreatePsi4ESPInputFile(poltype,structfname,poltype.comespfname.replace('.com','_frag.com'),molecule,poltype.maxdisk,poltype.maxmem,poltype.numproc,charge,False)
    finished,error=poltype.CheckNormalTermination(outputname)
    if finished==False and error==False:
        cmdstr='psi4 '+inputname+' '+outputname
        try:
             poltype.call_subsystem(cmdstr,True)
        except:
             error=True
    if error==False:
        WBOmatrix=GrabWBOMatrixPsi4(poltype,outputname,molecule)
    poltype.espmethod=curespmethod
    poltype.espbasisset=curspbasisset

    return WBOmatrix,outputname,error

def GenerateFragments(poltype,mol,torlist,parentWBOmatrix):
    newdir='Fragments'
    if not os.path.isdir(newdir):
        os.mkdir(newdir)
    os.chdir(newdir)
    bondtopo,bondsettobondorder=GenerateBondTopology(poltype,mol)
    valtopo=GenerateValenceTopology(poltype,mol)
    fragspath=os.getcwd()
    rotbndindextoparentindextofragindex={}
    rotbndindextofragment={}
    rotbndindextofragmentfilepath={}
    rotbndindextofragWBOmatrix={}
    rotbndindextofragfoldername={}
    rotbndindextoWBOdifference={}
    
    for tor in torlist: 
        WBOdifferencetofragWBOmatrix={}
        WBOdifferencetofoldername={}
        WBOdifferencetofragmol={}
        WBOdifferencetostructfname={}
        highlightbonds=[] 
        mols=[]
        indexes=FirstPassAtomIndexes(poltype,tor)

        fragfoldername=str(tor[1])+'_'+str(tor[2])+'_Hydrated'
        if not os.path.isdir(fragfoldername):
            os.mkdir(fragfoldername)
        os.chdir(fragfoldername)
        origfragmol,parentindextofragindex,fragidxswithchangedval,fragidxtooriginalformalcharge=GenerateFrag(poltype,indexes,valtopo,mol)
        mols.append(origfragmol)
        growfragments=[]
        fragments=[]
        charge=Chem.rdmolops.GetFormalCharge(origfragmol)
        fragmoltocharge={}
        fragmol=HydrateMolecule(poltype,origfragmol,fragfoldername,bondtopo,bondsettobondorder,parentindextofragindex,mol,fragidxswithchangedval,fragidxtooriginalformalcharge)
        charge=Chem.rdmolops.GetFormalCharge(fragmol)
        fragmoltocharge[fragmol]=charge
        mols.append(fragmol)
        filename=fragfoldername+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        hydindexes=FindAddedHydrogenIndexesOnChargedAtoms(poltype,mols,parentindextofragindex)
        combindexlist=CombinationsHydIndexes(poltype,hydindexes,fragmol)
        fragments.append(fragmol) # append case where all H no CH3 added
        fragments,fragmoltocharge=ChargedCombinations(poltype,combindexlist,filename,fragments,fragmoltocharge)
        os.chdir('..')
        fragmoltoWBOmatrices={}
        fragmoltofragfoldername={}
        fragmoltobondindexlist={}
        for i in range(len(fragments)):
            fragmol=fragments[i]
            fragfoldername=str(tor[1])+'_'+str(tor[2])+'_Index'+'_'+str(i)

            fragmoltofragfoldername[fragmol]=fragfoldername
            if not os.path.isdir(fragfoldername):
                os.mkdir(fragfoldername)
            os.chdir(fragfoldername)
            rotbndidx=str(tor[1])+'_'+str(tor[2])
            filename=fragfoldername+'.mol'

            WriteRdkitMolToMolFile(poltype,fragmol,filename)
            fragmolbabel=ReadMolFileToOBMol(poltype,filename)
            WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
            WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
            structfname=filename.replace('.mol','.sdf')
            fragWBOmatrix,outputname,error=GenerateWBOMatrix(poltype,fragmol,filename.replace('.mol','_xyzformat.xyz'),fragmoltocharge)
            if error==True:
                os.chdir('..')
                continue
            fragmentWBOvalue=fragWBOmatrix[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]] # rdkit is 0 index based so need to subtract 1, babel is 1 indexbased
            parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1] # Matrix has 0,0 so need to subtract 1 from babel index
            WBOdifference=round(numpy.abs(fragmentWBOvalue-parentWBOvalue),3)
            WBOdifferencetofragmol[WBOdifference]=fragmol
            WBOdifferencetostructfname[WBOdifference]=structfname
            rotbndindextoWBOdifference[rotbndidx]=WBOdifference
            fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist,fragmoltocharge)

            os.chdir('..')

            WBOdifferencetofragWBOmatrix[WBOdifference]=fragWBOmatrix
            WBOdifferencetofoldername[WBOdifference]=fragfoldername
        WBOdifference=min(list(WBOdifferencetofragWBOmatrix))
        fragmol=WBOdifferencetofragmol[WBOdifference]
        structfname=WBOdifferencetostructfname[WBOdifference]
        fragWBOmatrix=WBOdifferencetofragWBOmatrix[WBOdifference]
        fragfoldername=WBOdifferencetofoldername[WBOdifference]
        rotbndindextofragfoldername[rotbndidx]=fragfoldername
        os.chdir(fragfoldername)

        fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        highlightbonds.append(fragrotbndidx)
        fragpath=os.getcwd() 
        grow=False
        growfragments.append(fragmol)
        fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist,fragmoltocharge)
        curdir=os.getcwd()
        os.chdir('..') 
        growfragmoltoWBOmatrices=fragmoltoWBOmatrices.copy()
        growfragmoltofragfoldername=fragmoltofragfoldername.copy()
        growfragmoltobondindexlist=fragmoltobondindexlist.copy()


        Draw2DMoleculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,tor,'CombinationsWithIndex')
        sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in fragments]
        Draw2DMoleculesWithWBO(poltype,sanitizedfragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,tor,'CombinationsWithoutIndex')

        os.chdir(curdir)
        if WBOdifference<=poltype.WBOtol: # then we consider the fragment good enough to transfer torsion parameters, so make this fragment into .sdf file
            pass
        else:
            grow=True
            fragmol,newindexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex,fragpath,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist=GrowFragmentOut(poltype,mol,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername,bondtopo,bondsettobondorder,valtopo,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,fragspath,fragmoltocharge)
            fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist,fragmoltocharge)
        curdir=os.getcwd()
        os.chdir(fragspath)
        Draw2DMoleculesWithWBO(poltype,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,tor,'FragmentGrowthWithIndex')
        sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in growfragments]
        Draw2DMoleculesWithWBO(poltype,sanitizedfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,tor,'FragmentGrowthWithoutIndex')

        os.chdir(curdir)


        structfname=structfname.replace('_xyzformat.xyz','.sdf')
        rotbndindextofragment[rotbndidx]=fragmol
        rotbndindextofragmentfilepath[rotbndidx]=fragpath+r'/'+structfname
        rotbndindextoparentindextofragindex[rotbndidx]=parentindextofragindex
        rotbndindextofragWBOmatrix[rotbndidx]=fragWBOmatrix
        rotbndindextofragfoldername[rotbndidx]=fragfoldername
        os.chdir(fragspath)
    # now remove all folders with Hydrated in them, that was just temp storage for producing other folders
    RemoveTempFolders(poltype)
    poltype.rotbndindextofragmentfilepath=rotbndindextofragmentfilepath
    fragmentarray=[]
    for rotbndindex in rotbndindextofragment.keys():
        fragment=rotbndindextofragment[rotbndindex]
        fragmentarray.append(fragment)

    equivalentfragmentsarray=FindEquivalentFragments(poltype,fragmentarray)
    equivalentrotbndindexarrays=FindEquivalentRotatableBonds(poltype,equivalentfragmentsarray,rotbndindextofragment)
    # now we need to redraw the 2Dimages for any fragments that are equivalent (get multiple torsions from different rotatable bonds around same fragment)
    curdir=os.getcwd()
    equivalentrotbndindexes=GrabEquivalentRotBndIndexes(poltype,equivalentrotbndindexarrays,rotbndindex)
    if len(equivalentrotbndindexes)!=0:
        for bndindexes in equivalentrotbndindexes:
            highlightbonds=[]
            for bndindex in bndindexes:
                parentindextofragindex=rotbndindextoparentindextofragindex[bndindex]
                indexes=bndindex.split('_')
                indexes=[int(i) for i in indexes]
                fragrotbndidx=[parentindextofragindex[indexes[0]-1],parentindextofragindex[indexes[1]-1]]
                if fragrotbndidx not in highlightbonds:
                    highlightbonds.append(fragrotbndidx)
            for bndindex in bndindexes:
                fragmol=rotbndindextofragment[bndindex]
                fragmentfilepath=rotbndindextofragmentfilepath[bndindex]
                head,tail=os.path.split(fragmentfilepath)
                WBOdifference=rotbndindextoWBOdifference[bndindex]
                parentindextofragindex=rotbndindextoparentindextofragindex[bndindex]
                fragWBOmatrix=rotbndindextofragWBOmatrix[bndindex]
                indexes=bndindex.split('_')
                indexes=[int(i) for i in indexes]
                fragfoldername=rotbndindextofragfoldername[bndindex]
                os.chdir(head)
                basename=fragfoldername+'_WBO_'+str(round(WBOdifference,3))
                fragrotbndidx=[parentindextofragindex[indexes[0]-1],parentindextofragindex[indexes[1]-1]]
                reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
                relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
                m=mol_with_atom_index(poltype,fragmol)
                fragsmirks=rdmolfiles.MolToSmarts(m)
                structfnamemol=fragfoldername+'.mol'
                Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',m,bondindexlist=highlightbonds,smirks=fragsmirks)
                Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',m,bondindexlist=highlightbonds,smirks=fragsmirks)
        os.chdir(curdir)
    sys.exit()
    return rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentfragmentsarray,equivalentrotbndindexarrays


def GrabEquivalentRotBndIndexes(poltype,equivalentrotbndindexarrays,rotbndindex):
    equivalentrotbndindexes=[]
    for ls in equivalentrotbndindexarrays:
        if rotbndindex in ls:
            temp=[]
            for el in ls:
                if el!=rotbndindex:
                    temp.append(el)
            equivalentrotbndindexes.append(temp)
    return equivalentrotbndindexes

def RemoveTempFolders(poltype):
    foldstoremove=[]
    folds=os.listdir()
    for f in folds:
        if os.path.isdir(f) and 'Hydrated' in f:
            foldstoremove.append(f)
    for f in foldstoremove:
        shutil.rmtree(f)
    


def ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix):
    reducedparentWBOmatrix=numpy.copy(fragWBOmatrix)
    fragindextoparentindex={v: k for k, v in parentindextofragindex.items()}
    for i in range(len(fragWBOmatrix)):
        for j in range(len(fragWBOmatrix[0])):
            fragrowindex=i
            fragcolindex=j
            if fragrowindex in fragindextoparentindex.keys() and fragcolindex in fragindextoparentindex.keys():
                parentrowindex=fragindextoparentindex[fragrowindex]
                parentcolindex=fragindextoparentindex[fragcolindex]
                parentvalue=parentWBOmatrix[parentrowindex,parentcolindex]
            else:
                parentvalue=0
            reducedparentWBOmatrix[i,j]=parentvalue
    return reducedparentWBOmatrix

def CombinationsHydIndexes(poltype,hydindexes,fragmol): # only keep combinations of hydrogens that are attachated to heavy atoms that have a formal charge in the parent molecule
       
    combindexlist=[]
    for i in range(len(hydindexes)):
        comb = combinations(hydindexes, i+1)
        combindexlist.append(comb)
    return combindexlist 


def ChargedCombinations(poltype,combindexlist,molfilename,fragments,fragmoltocharge): # only worry about charged combinations where you add hydrogen on atoms that were charged in parent molecule
    fragmentsmarts=[]
    j=0
    for comb in combindexlist:
        i=0
        for idxlist in comb:
            idxlist=list(idxlist)
            fragmolcharged=ReadRdkitMolFromMolFile(poltype,molfilename) 
            fragmolchargededit=Chem.rdchem.EditableMol(fragmolcharged)
            idxlist.sort(reverse=True)
            origcharge=Chem.rdmolops.GetFormalCharge(fragmolcharged)
            for idx in idxlist:
                fragmolchargededit.RemoveAtom(idx)
                 
            newfragmolcharged=fragmolchargededit.GetMol()
            fragsmartscharged=rdmolfiles.MolToSmarts(newfragmolcharged)
            if fragsmartscharged not in fragmentsmarts:
                charge=origcharge-len(idxlist)
                fragmoltocharge[newfragmolcharged]=charge
                fragmentsmarts.append(fragsmartscharged)
                fragments.append(newfragmolcharged)
            i+=1
        j+=1
    return fragments,fragmoltocharge
                

def FindAddedHydrogenIndexesOnChargedAtoms(poltype,mols,parentindextofragindex):
    fragindextoparentindex=dict(map(reversed, parentindextofragindex.items()))
    hydindexes=[]
    hydratedmol=mols[1]
    originalmol=mols[0]
    smarts=rdmolfiles.MolToSmarts(originalmol)
    matches = hydratedmol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    firstmatch=matches[0]
    for atom in hydratedmol.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx not in firstmatch: # if its an H
            for natom in atom.GetNeighbors():
                natomidx=natom.GetIdx()
                parentindex=fragindextoparentindex[natomidx]
                parentatom=poltype.rdkitmol.GetAtomWithIdx(parentindex)
                formalcharge=parentatom.GetFormalCharge()
                if formalcharge!=0: # only worry about hydrogen combinations on atoms that have formal charge in the parent molecule
                    hydindexes.append(atomidx)
    return hydindexes

def FindAddedHydrogenIndexes(poltype,mols):
    hydindexes=[]
    hydratedmol=mols[1]
    originalmol=mols[0]
    smarts=rdmolfiles.MolToSmarts(originalmol)
    matches = hydratedmol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    firstmatch=matches[0]
    for atom in hydratedmol.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx not in firstmatch: # if its an H
            hydindexes.append(atomidx)
    return hydindexes



def GrowFragmentOut(poltype,mol,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername,bondtopo,bondsettobondorder,valtopo,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,fragspath,fragmoltocharge):
    fragfoldernamepath=os.getcwd()
    fragmentsforcomb=[]
    while not WBOdifference<=poltype.WBOtol:
        WBOdiffarray=[]
        molarray=[]
        fragmolidxtoparentindextofragindex={}
        fragmentidxtostructfname={}
        fragmolidxtofoldername={}
        fragmolidxtofragmol={}
        parentindextofragindexlist=[]

        fragments=[]
        fragmentlist=[]
        possiblefragatmidxs=GrowPossibleFragmentAtomIndexes(poltype,poltype.rdkitmol,indexes)
        for fragmolidx in range(len(possiblefragatmidxs)):
            indexlist=possiblefragatmidxs[fragmolidx]
            basename=fragfoldername+'_GrowFragment_Hydrated'+'_'+str(fragmolidx) 
            if not os.path .isdir(basename):
                os.mkdir(basename)
            os.chdir(basename)
            mols=[]
            origfragmol,parentindextofragindex,fragidxswithchangedval,fragidxtooriginalformalcharge=GenerateFrag(poltype,indexlist,valtopo,mol)
            charge=Chem.rdmolops.GetFormalCharge(origfragmol)
            mols.append(origfragmol)
            fragmol=HydrateMolecule(poltype,origfragmol,basename,bondtopo,bondsettobondorder,parentindextofragindex,mol,fragidxswithchangedval,fragidxtooriginalformalcharge)
            mols.append(fragmol)
            charge=Chem.rdmolops.GetFormalCharge(fragmol)
            fragmoltocharge[fragmol]=charge
            filename=basename+'.mol'
            WriteRdkitMolToMolFile(poltype,fragmol,filename)
            fragmolbabel=ReadMolFileToOBMol(poltype,filename)
            WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
            WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
            hydindexes=FindAddedHydrogenIndexesOnChargedAtoms(poltype,mols,parentindextofragindex)
            combindexlist=CombinationsHydIndexes(poltype,hydindexes,fragmol)
            
            fragments.append(fragmol) # include the case where all H and no H converted to CH3
            fragments,fragmoltocharge=ChargedCombinations(poltype,combindexlist,filename,fragments,fragmoltocharge)
            os.chdir('..')
            for i in range(len(fragments)):
                frag=fragments[i]
                fragmentlist.append(frag)
                parentindextofragindexlist.append(parentindextofragindex)
        for fragmolidx in range(len(fragmentlist)):
            fragmol=fragmentlist[fragmolidx]
            if fragmol not in fragmentsforcomb:
                fragmentsforcomb.append(fragmol)
            basename=fragfoldername+'_GrowFragment_'+str(fragmolidx)
            if not os.path .isdir(basename):
                os.mkdir(basename)
            os.chdir(basename)
            growfragmoltofragfoldername[fragmol]=basename
            filename=basename+'.mol'
            WriteRdkitMolToMolFile(poltype,fragmol,filename)
            parentindextofragindex=parentindextofragindexlist[i]
            fragmolbabel=ReadMolFileToOBMol(poltype,filename)
            WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
            WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
            os.chdir('..')
            fragmolidxtofragmol[fragmolidx]=fragmol
            fragmolidxtofoldername[fragmolidx]=basename
            fragmolidxtoparentindextofragindex[fragmolidx]=parentindextofragindex
            fragmentidxtostructfname[fragmolidx]=filename.replace('.mol','_xyzformat.xyz')
        WBOdifftoindexlist={}
        WBOdifftofragmol={}
        WBOdifftofragWBOmatrix={}
        WBOdifftofolder={}
        WBOdifftostructfname={}
        WBOdifftoparentindextofragindex={}
        for fragmolidx in fragmolidxtofragmol.keys():
            fragmol=fragmolidxtofragmol[fragmolidx]
            foldername=fragmolidxtofoldername[fragmolidx]
            parentindextofragindex=fragmolidxtoparentindextofragindex[fragmolidx]
            structfname=fragmentidxtostructfname[fragmolidx]
            os.chdir(foldername)

            fragWBOmatrix,outputname,error=GenerateWBOMatrix(poltype,fragmol,structfname,fragmoltocharge)
            if error==True:
                os.chdir('..')
                continue
            reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
            relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
            fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
            fragmentWBOvalue=fragWBOmatrix[fragrotbndidx[0],fragrotbndidx[1]]
            parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1]
            WBOdifference=round(numpy.abs(fragmentWBOvalue-parentWBOvalue),3)
            growfragmoltoWBOmatrices,growfragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,foldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,growfragmoltoWBOmatrices,growfragmoltobondindexlist,fragmoltocharge)

            m=mol_with_atom_index(poltype,fragmol)
            os.chdir('..')
            indexlist=list(parentindextofragindex.keys())
            WBOdifftoparentindextofragindex[WBOdifference]=parentindextofragindex
            WBOdifftoindexlist[WBOdifference]=indexlist
            WBOdifftofragmol[WBOdifference]=fragmol
            WBOdifftofragWBOmatrix[WBOdifference]=fragWBOmatrix
            WBOdifftofolder[WBOdifference]=foldername
            WBOdifftostructfname[WBOdifference]=structfname
            molarray.append(fragmol)
            WBOdiffarray.append(WBOdifference)
        WBOdifference=min(WBOdifftoindexlist.keys())
        parentindextofragindex=WBOdifftoparentindextofragindex[WBOdifference]
        indexes=WBOdifftoindexlist[WBOdifference]
        foldername=WBOdifftofolder[WBOdifference]
        structfname=WBOdifftostructfname[WBOdifference]
        RemoveTempFolders(poltype)
        os.chdir(foldername)

        fragmol=WBOdifftofragmol[WBOdifference]
        growfragments.append(fragmol)
        fragWBOmatrix=WBOdifftofragWBOmatrix[WBOdifference]
        fragpath=os.getcwd()

    curdir=os.getcwd()
    os.chdir('..')
    Draw2DMoleculesWithWBO(poltype,fragmentsforcomb,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,tor,'CombinationsWithIndex')
    sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in fragmentsforcomb]
    Draw2DMoleculesWithWBO(poltype,sanitizedfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,tor,'CombinationsWithoutIndex')
    os.chdir(curdir)
    os.chdir(fragfoldernamepath)
    PlotFragmenterResults(poltype,WBOdiffarray,molarray)
    os.chdir(curdir)
    
    return fragmol,indexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex,fragpath,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist


def GrowPossibleFragmentAtomIndexes(poltype,rdkitmol,indexes):
    possiblefragatmidxs=[]
    comblist=[]
    for bond in rdkitmol.GetBonds():
        aidx=bond.GetBeginAtomIdx()
        bidx=bond.GetEndAtomIdx()
        aatom=rdkitmol.GetAtomWithIdx(aidx)
        batom=rdkitmol.GetAtomWithIdx(bidx)
        aatomicnum=aatom.GetAtomicNum()
        batomicnum=batom.GetAtomicNum()
        bondorder=bond.GetBondTypeAsDouble()
        if bondorder>1:
            continue
        if aatomicnum==1 or batomicnum==1:
            continue
        if (aidx in indexes and bidx not in indexes): # then this means the bond is not already in the fragment but this is one of the bonds just outside of the fragment
            idx=bidx
        elif (aidx not in indexes and bidx in indexes):
            idx=aidx
        else:
            continue
        comblist.append(idx)
    combinationslist=[]
    length=len(comblist)
    for i in range(length):
        comb=combinations(comblist,i+1)
        for ls in comb:
            combinationslist.append(ls)
    for comb in combinationslist:
        indexlist=indexes.copy()
        for idx in comb:
           aromaticindexes=GrabAromaticAtoms(poltype,rdkitmol.GetAtomWithIdx(idx))
           newindexes=aromaticindexes
           for atmidx in newindexes:
               if atmidx not in indexlist:
                   indexlist.append(atmidx)
        temp=[]
        for idx in indexlist:
           neighbatom=poltype.rdkitmol.GetAtomWithIdx(idx)
           for neighbneighbatom in neighbatom.GetNeighbors():
               atomicnum=neighbneighbatom.GetAtomicNum()
               neighbneighbatomidx=neighbneighbatom.GetIdx()
               if atomicnum==1 and neighbneighbatomidx not in indexlist:
                   temp.append(neighbneighbatomidx)
               bond=poltype.rdkitmol.GetBondBetweenAtoms(neighbneighbatomidx,idx)
               bondorder=bond.GetBondTypeAsDouble()
               if bondorder>1 and neighbneighbatomidx not in indexlist:
                   temp.append(neighbneighbatomidx)
        for idx in temp:
            indexlist.append(idx)

        if indexlist not in possiblefragatmidxs:
           possiblefragatmidxs.append(indexlist)

    return possiblefragatmidxs
    

def WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist,fragmoltocharge):
    highlightbonds=[]
    structfnamemol=fragfoldername+'.mol'
    charge=fragmoltocharge[fragmol]
    # just put all the charge on the first atom for quantum
    firstatom=fragmol.GetAtomWithIdx(0)
    firstatomcurcharge=firstatom.GetFormalCharge()
    firstatom.SetFormalCharge(firstatomcurcharge+charge)
    print(Chem.MolToMolBlock(fragmol,kekulize=True),file=open(structfnamemol,'w+')) 
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('mol')
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,structfnamemol)
    tmpconv.SetOutFormat('sdf')
    structfname=fragfoldername+'.sdf'
    tmpconv.WriteFile(fragmolbabel,structfname)
    basename=fragfoldername+'_WBO_'+str(round(WBOdifference,3))
    fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
    highlightbonds.append(fragrotbndidx)
    reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
    relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
    m=mol_with_atom_index(poltype,fragmol)
    fragsmirks=rdmolfiles.MolToSmarts(m)
    Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',m,bondindexlist=highlightbonds,smirks=fragsmirks)
    Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',m,bondindexlist=highlightbonds,smirks=fragsmirks)
    temp=[relativematrix,fragWBOmatrix]
    fragmoltoWBOmatrices[fragmol]=temp 
    fragmoltobondindexlist[fragmol]=highlightbonds
    return fragmoltoWBOmatrices,fragmoltobondindexlist

def CheckIfIndexInMatches(poltype,index,indexlist):
    match=True
    for ls in indexlist:
        if index not in ls:
            match=False
            return match
    return match

def MapSMILESToParent(poltype,mol,smi,temptorlist):
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smi)
    match=sp.Match(mol)
    if match==False:
        return None,None
    indexlist=[]
    for indexls in sp.GetMapList():
        indexlist.append(indexls)
    natoms=len(indexlist[0])	
    for tor in temptorlist:
        foundall=True
        for index in tor:
            match=CheckIfIndexInMatches(index,indexlist)
            if match==False:
                foundall=False
        if foundall==True:
            return str(tor[1])+'_'+str(tor[2]),natoms
     
    return None,None

def FirstPassAtomIndexes(poltype,tor):
   molindexlist=[]
   for atom in poltype.rdkitmol.GetAtoms():
       atomindex=atom.GetIdx()
       babelatomindex=atomindex+1
       if babelatomindex in tor:
           if atomindex not in molindexlist:
               molindexlist.append(atomindex)
           for neighbatom in atom.GetNeighbors():
               neighbatomindex=neighbatom.GetIdx()
               if neighbatomindex not in molindexlist:
                   molindexlist.append(neighbatomindex)
                   if neighbatom.GetIsAromatic()==True:
                       aromaticindexes=GrabAromaticAtoms(poltype,neighbatom)
                       newindexes=aromaticindexes
                       for atmidx in newindexes:
                           if atmidx not in molindexlist:
                               molindexlist.append(atmidx)
   temp=[]
   for index in molindexlist:
       atom=poltype.rdkitmol.GetAtomWithIdx(index)
       for neighbneighbatom in atom.GetNeighbors():
           atomicnum=neighbneighbatom.GetAtomicNum()
           neighbneighbatomidx=neighbneighbatom.GetIdx()
           if atomicnum==1 and neighbneighbatomidx not in molindexlist:
               temp.append(neighbneighbatomidx)
           bond=poltype.rdkitmol.GetBondBetweenAtoms(neighbneighbatomidx,index)
           bondorder=bond.GetBondTypeAsDouble()
           if bondorder>1 and neighbneighbatomidx not in molindexlist:
               temp.append(neighbneighbatomidx)
   for idx in temp:
       molindexlist.append(idx)
   return molindexlist

def Chunks(poltype,lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
    
def ChunksList(poltype,gen):
    newlst=[]
    for item in gen:
        newlst.append(item)
    return newlst

def MatchSMARTSOB(poltype,smarts,filename):
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(filename)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(mol,filename)
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    diditmatch=sp.Match(mol)
    return diditmatch




def Draw2DMoleculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,tor,basestr):

    bondlistoflists=[]
    for frag in fragments:
        bondindexlist=fragmoltobondindexlist[frag]
        bondlist=[]
        for bondindexes in bondindexlist:
            bond=frag.GetBondBetweenAtoms(bondindexes[0],bondindexes[1])
            bondidx=bond.GetIdx()
            bondlist.append(bondidx) 
        bondlistoflists.append(bondlist)
    legendslist=[fragmoltofragfoldername[frag] for frag in fragments]
    molsperrow=3
    molsPerImage=molsperrow**2
    imagesize=400
    for i in range(len(fragments)):
        frag=fragments[i]
        rdDepictor.Compute2DCoords(frag)

    #tmats=[]
    #tmats.append(1)
    if len(fragments)>1:
        firstmol=fragments[0]
        editmol=Chem.rdchem.EditableMol(firstmol)
        firstmolcopy=editmol.GetMol()
        newmol=mol_with_atom_index_removed(poltype,firstmolcopy)
        newermol = Chem.RemoveHs(newmol)
        smarts=rdmolfiles.MolToSmarts(newermol)
        smarts=smarts.replace('@','').replace('H3','').replace('H2','').replace('H','')
        patt = Chem.MolFromSmarts(smarts)
        
        print('smarts',smarts)
        Chem.SanitizeMol(patt)
        print('fragmoltofragfoldername',fragmoltofragfoldername)
        for i in range(1,len(fragments)):
            frag=fragments[i]
            frag=mol_with_atom_index_removed(poltype,frag)
            fragsmarts=rdmolfiles.MolToSmarts(frag)
            fragpatt = Chem.MolFromSmarts(smarts)

            selfoverlap = frag.GetSubstructMatch(fragpatt) 
            editmol=Chem.rdchem.EditableMol(frag)
            fragcopy=editmol.GetMol()
            newfrag=mol_with_atom_index_removed(poltype,fragcopy)
            newerfrag = Chem.RemoveHs(newfrag)
            fragsmartsnoh=rdmolfiles.MolToSmarts(newerfrag)
            fragsmartsnoh=fragsmartsnoh.replace('@','').replace('H3','').replace('H2','').replace('H','')
            fragpattnoh = Chem.MolFromSmarts(fragsmartsnoh)
            selfoverlapnoh = frag.GetSubstructMatch(fragpattnoh)
            print('fragsmarts',fragsmarts)
            print('selfoverlap',selfoverlap)
            print('selfoverlapnoh',selfoverlapnoh)
            print('fragsmartsnoh',fragsmartsnoh)
            print('fragfoldername',fragmoltofragfoldername[frag])
            print(Chem.MolToMolBlock(frag,kekulize=True),file=open('checkfrag.mol','w+'))
            Chem.SanitizeMol(fragpattnoh)

            AllChem.EmbedMolecule(fragpattnoh)

            print(Chem.MolToMolBlock(fragpattnoh,kekulize=True),file=open('checkfragnoh.mol','w+')) 

            diditmatch=MatchSMARTSOB(poltype,smarts,'checkfrag.mol')
            print('diditmatch',diditmatch)
            overlap = frag.GetSubstructMatch(patt) # indexes of fragpatt corresponding to patt SMARTS but need the actual indexes of frag 
            atomMap = [(paid,raid) for raid,paid in enumerate(overlap)]
            AllChem.EmbedMolecule(frag)
            AllChem.EmbedMolecule(patt)
            print(Chem.MolToMolBlock(patt,kekulize=True),file=open('checkpatt.mol','w+')) 
            print('overlap',overlap)
            print('atomMap',atomMap)
            AllChem.AlignMol(frag,firstmol,atomMap=atomMap)
            #tuple=Chem.rdMolAlign.GetAlignmentTransform(frag,firstmol,atomMap=atomMap)
            #tmat=tuple[1]
            #tmats.append(tmat)
    fragmentchunks=ChunksList(poltype,Chunks(poltype,fragments,molsPerImage))
    legendschunks=ChunksList(poltype,Chunks(poltype,legendslist,molsPerImage))
    bondlistoflistschunks=ChunksList(poltype,Chunks(poltype,bondlistoflists,molsPerImage))
    #tmatchunks=ChunksList(poltype,Chunks(poltype,tmats,molsPerImage))
    #print('tmatchunks',tmatchunks)
    for i in range(len(fragmentchunks)):
        fragmentsublist=fragmentchunks[i]
        legendssublist=legendschunks[i]
        bondlistoflistssublist=bondlistoflistschunks[i]
        #tmatsublist=tmatchunks[i]
        svg=Chem.Draw.MolsToGridImage(fragmentsublist,molsPerRow=molsperrow,subImgSize=(imagesize,imagesize),legends=legendssublist,highlightBondLists=bondlistoflistssublist,useSVG=True)
        fig = sg.fromstring(svg)
        ls=range(len(fragmentsublist))
        chunks=ChunksList(poltype,Chunks(poltype,ls,molsperrow))
        indextorow={}
        for rowidx in range(len(chunks)):
            row=chunks[rowidx]
            for j in row:
                indextorow[j]=rowidx
        for j in range(len(fragmentsublist)):
            frag=fragmentsublist[j]
            bondlist=bondlistoflistssublist[j]
            legend=legendssublist[j]
            #tmat=tmatsublist[j]
            #if i==0 and j==0:
            #    tmat=1
            #print('tmat',tmat,flush=True)
            drawer=rdMolDraw2D.MolDraw2DSVG(imagesize,imagesize)
            rdDepictor.Compute2DCoords(frag)
            drawer.DrawMolecule(frag,highlightAtoms=[],highlightBonds=bondlist)
            drawer.FinishDrawing()
            newsvg = drawer.GetDrawingText().replace('svg:','')
            newfig = sg.fromstring(newsvg)
            basename=legend
            newfig.save(basename+'.svg')
            newsvg_code=newfig.to_str()
            svg2png(bytestring=newsvg_code,write_to=basename+'.png')

            atomidxtodrawcoords={}
            for bond in frag.GetBonds():
                bondidx=bond.GetIdx()
                if bondidx in bondlist:
                    begidx=bond.GetBeginAtomIdx()
                    endidx=bond.GetEndAtomIdx()
                    begatomdrawcoords=numpy.array(drawer.GetDrawCoords(begidx))
                    endatomdrawcoords=numpy.array(drawer.GetDrawCoords(endidx))
                    atomidxtodrawcoords[begidx]=begatomdrawcoords
                    atomidxtodrawcoords[endidx]=endatomdrawcoords
            
            WBOmatrixlist=fragmoltoWBOmatrices[frag]
            WBOmatrix=WBOmatrixlist[0]
            row=indextorow[j]
            x=(j-molsperrow*(row))*imagesize
            y=(row)*imagesize
            shift=numpy.array([x,y])
            for bond in frag.GetBonds():
                bondidx=bond.GetIdx()
                if bondidx in bondlist:
                    begidx=bond.GetBeginAtomIdx()
                    endidx=bond.GetEndAtomIdx()
                    begatomdrawcoords=atomidxtodrawcoords[begidx]+shift
                    endatomdrawcoords=atomidxtodrawcoords[endidx]+shift
                    bondcoords=(begatomdrawcoords+endatomdrawcoords)/2
                    WBOval=numpy.abs(WBOmatrix[begidx,endidx])
                    if WBOval==0:
                        continue
                    wbo=str(round(WBOval,4))
                    label = sg.TextElement(bondcoords[0],bondcoords[1], wbo, size=12, weight="bold")
                    array=endatomdrawcoords-begatomdrawcoords
                    if array[1]>=0:
                        pass
                    else:
                        array=-1*array
                    norm = numpy.linalg.norm(array)
                    normarray=array/norm
                    angle=numpy.abs(numpy.degrees(numpy.arccos(normarray[1])))
                    if angle>90:
                        angle=angle-90
                    if normarray[1]>=0 and normarray[0]>=0:
                        sign=-1
                    elif normarray[1]<=0 and normarray[0]<=0:
                        sign=-1
                    else:
                        sign=1
                    label.rotate(sign*angle,bondcoords[0],bondcoords[1])

                    fig.append(label)
        
        basename=basestr+'_'+'Bnd_'+str(tor[1])+'-'+str(tor[2])+'_Index_'+str(i)
        fig.save(basename+'.svg')
        svg_code=fig.to_str()
        svg2png(bytestring=svg_code,write_to=basename+'.png')


def Draw2DMoleculeWithWBO(poltype,WBOmatrix,basename,mol,bondindexlist=None,smirks=None):
    mol=mol_with_atom_index(poltype,mol)
    rdDepictor.Compute2DCoords(mol)
    drawer=rdMolDraw2D.MolDraw2DSVG(500,500)
    bondlist=[]
    if bondindexlist!=None:
        for bondindexes in bondindexlist:
            bond=mol.GetBondBetweenAtoms(bondindexes[0],bondindexes[1])
            bondidx=bond.GetIdx()
            bondlist.append(bondidx) 
    drawer.DrawMolecule(mol,highlightAtoms=[],highlightBonds=bondlist)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    fig = sg.fromstring(svg)
    for bond in mol.GetBonds():
        bondidx=bond.GetIdx()
        if bondidx in bondlist:
            begidx=bond.GetBeginAtomIdx()
            endidx=bond.GetEndAtomIdx()
            begatomdrawcoords=numpy.array(drawer.GetDrawCoords(begidx))
            endatomdrawcoords=numpy.array(drawer.GetDrawCoords(endidx))
            bondcoords=(begatomdrawcoords+endatomdrawcoords)/2
            WBOval=numpy.abs(WBOmatrix[begidx,endidx])
            if WBOval==0:
                continue
            wbo=str(round(WBOval,4))
            label = sg.TextElement(bondcoords[0],bondcoords[1], wbo, size=12, weight="bold")
            array=endatomdrawcoords-begatomdrawcoords 
            if array[1]>=0:
                pass
            else:
                array=-1*array
            norm = numpy.linalg.norm(array)
            normarray=array/norm
            angle=numpy.abs(numpy.degrees(numpy.arccos(normarray[1])))
            if angle>90:
                angle=angle-90
            if normarray[1]>=0 and normarray[0]>=0:
                sign=-1
            elif normarray[1]<=0 and normarray[0]<=0:
                sign=-1
            else:
                sign=1
            label.rotate(sign*angle,bondcoords[0],bondcoords[1])
            fig.append(label)
    if smirks!=None:
        label = sg.TextElement(25,490, smirks, size=12, weight="bold")
        fig.append(label)
    fig.save(basename+'.svg')
    svg_code=fig.to_str()
    svg2png(bytestring=svg_code,write_to=basename+'.png')

def GrabAromaticAtoms(poltype,neighbatom):
    aromaticindexes=[]
    prevringidxlen=len(aromaticindexes)
    aromaticindexes.append(neighbatom.GetIdx())
    ringidxlen=len(aromaticindexes)
    while prevringidxlen!=ringidxlen:
        for atmindex in aromaticindexes:
            atm=poltype.rdkitmol.GetAtomWithIdx(atmindex)
            if (atm.GetIsAromatic()==True) and atmindex not in aromaticindexes:
                aromaticindexes.append(atmindex)
            for natm in atm.GetNeighbors():
                if (natm.GetIsAromatic()==True) and natm.GetIdx() not in aromaticindexes:
                    aromaticindexes.append(natm.GetIdx())
        prevringidxlen=ringidxlen
        ringidxlen=len(aromaticindexes)

    return aromaticindexes


def PlotFragmenterResults(poltype,WBOdiffarray,molarray):
    fig=plt.figure()
    basename='NumberofAtomsVSWBODifference'
    plt.plot(WBOdiffarray,[m.GetNumAtoms() for m in molarray],'.') 
    plt.xlabel('WBO Difference')
    plt.ylabel('Number of atoms in fragment')
    fig.savefig(basename+'.png')


