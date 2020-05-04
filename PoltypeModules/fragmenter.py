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
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap=oldindextonewindex
    newindextooldindex={v: k for k, v in oldindextonewindex.items()}
    fragidxswithchangedval=[]
    fragidxtooriginalformalcharge={}
    for atom in res.GetAtoms():
        atomidx=atom.GetIdx()
        oldindex=newindextooldindex[atomidx]
        babelidx=oldindex+1
        parentval=valtopo[babelidx]
        fragval=len(atom.GetNeighbors())
        diff=fragval-parentval
        parentatom=poltype.rdkitmol.GetAtomWithIdx(oldindex)
        parentatomformalcharge=parentatom.GetFormalCharge()
        newformalcharge=parentatomformalcharge+diff
        atom.SetFormalCharge(newformalcharge)
        fragidxtooriginalformalcharge[atomidx]=parentatomformalcharge
        if diff!=0:
            fragidxswithchangedval.append(atomidx)

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
'''
def PreserveBondTopology(poltype,newOBmol,bondtopo,bondsettobondorder,parentidxtofragidx,parentmol,fragidxswithchangedval):
    iterbond=openbabel.OBMolBondIter(parentmol) # iterator for all bond objects in the molecule
    for bond in iterbond:
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        aidx=a.GetIdx()
        bidx=b.GetIdx()
        if (aidx-1) in parentidxtofragidx:
            fragaidx=parentidxtofragidx[aidx-1]
        else:
            continue
        if (bidx-1) in parentidxtofragidx:
            fragbidx=parentidxtofragidx[bidx-1]
        else:
            continue
        bondlist=[aidx,bidx]
        bondset=set(bondlist)
        bondorder=bondsettobondorder[tuple(bondset)]
        fraga=newOBmol.GetAtom(fragaidx+1)
        fragb=newOBmol.GetAtom(fragbidx+1)
        fragbond=newOBmol.GetBond(fraga,fragb)

        if fragbond==None: # then hydrogen may have been too far from atom and babel thinks there is no bond
            newOBmol.AddBond(fragaidx+1,fragbidx+1,bondorder)
        else:
            if fragaidx in fragidxswithchangedval or fragbidx in fragidxswithchangedval:
                fragbond.SetBondOrder(1)
            else:
                fragbond.SetBondOrder(bondorder)

    return newOBmol 
''' 


 
def HydrateMolecule(poltype,mol,basename,bondtopo,bondsettobondorder,parentidxtofragidx,parentmol,fragidxswithchangedval,fragidxtooriginalformalcharge): # need to use openbabel to do this and convert back to rdkit since rdkit has issues with Sanitization and cannot be turned off for EmbedMolecule function. Rkdit assigns in correct bond order on rings in this case (4). So convert to xyz to remove bad bond order information, then add back parent bond topology, remove any radicals OB assigns before adding hydrogens.
    for atom in mol.GetAtoms():
        atomidx=atom.GetIdx()
        originalcharge=fragidxtooriginalformalcharge[atomidx]
        atom.SetFormalCharge(originalcharge) 
    hydmol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(hydmol)
    #OBmol=ConvertRdkitMolToOBMol(poltype,mol)
    #outputname='intermediate2.xyz'
    #WriteOBMolToXYZ(poltype,OBmol,outputname)
    #newOBmol=ReadToOBMol(poltype,outputname)
    #newOBmol=PreserveBondTopology(poltype,newOBmol,bondtopo,bondsettobondorder,parentidxtofragidx,parentmol,fragidxswithchangedval)
    #newname='intermediate3.mol'
    #WriteOBMolToMol(poltype,newOBmol,newname)
    #tempname=RemoveRadicals(poltype,newname)
    #SanitizeAtomBlock(poltype,tempname)
    #finalname='final.mol'
    #cmdstr='babel'+' '+'-imol'+' '+tempname+' '+'-h'+' '+'--gen3D'+' '+'-omol'+' '+'-O'+' '+finalname
    #poltype.call_subsystem(cmdstr,True)
    #hydmol=ReadRdkitMolFromMolFile(poltype,finalname)
    return hydmol

'''
def SanitizeAtomBlock(poltype,filename):
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    tempname=filename.replace('.mol','_TEMP2.mol')
    temp=open(tempname,'w')
    for line in results:
        linesplit=line.split()
        if len(linesplit)==16: # this is the atom block
            fullsplit=re.split(r'(\s+)', line)
            idxstochange=[]
            for i in range(len(fullsplit)):
                if i>8 and fullsplit[i].isdigit():
                    idxstochange.append(i)
            for idx in idxstochange:
                fullsplit[idx]='0'
            newline=''.join(fullsplit)
            temp.write(newline)
        else:
            temp.write(line)
    temp.close()
    os.remove(filename)
    os.rename(tempname,filename)

def RemoveRadicals(poltype,filename):
    while not os.path.isfile(filename):
        time.sleep(1)
    if os.path.isfile(filename):
        filenamesplit=filename.split('.')
        ext=filenamesplit[1]
        tempname=filename.replace('.'+ext,'_TEMP.'+ext)
        temp=open(filename,'r')
        results=temp.readlines()
        temp.close()
        temp=open(tempname,'w')
        for line in results:
            if 'RAD' in line:
                pass
            else:
                temp.write(line)
        temp.close()
        return tempname

def ConvertOBMolToRdkitMol(poltype,OBmol):
    outputname='intermediate3.mol'
    WriteOBMolToMol(poltype,OBmol,outputname)
    rdkitmol=ReadRdkitMolFromMolFile(poltype,outputname)
    return rdkitmol
'''

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
    rdmolfiles.MolToMolFile(mol,outputname,kekulize=False)

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

def GenerateWBOMatrix(poltype,molecule,structfname):
    curespmethod=poltype.espmethod
    curspbasisset=poltype.espbasisset
    poltype.espmethod='HF'
    poltype.espbasisset='MINIX'
    charge=Chem.rdmolops.GetFormalCharge(molecule)
    print('molecule formal charge',charge)
    for atom in molecule.GetAtoms():
        atomidx=atom.GetIdx()
        atomcharge=atom.GetFormalCharge()
        print('atomidx ',atomidx,'formal charge',atomcharge)
    inputname,outputname=esp.CreatePsi4ESPInputFile(poltype,structfname,poltype.comespfname.replace('.com','_frag.com'),molecule,poltype.maxdisk,poltype.maxmem,poltype.numproc,charge,False)
    finished,error=poltype.CheckNormalTermination(outputname)
    if finished==False and error==False:
        cmdstr='psi4 '+inputname+' '+outputname
        poltype.call_subsystem(cmdstr,True)
    WBOmatrix=GrabWBOMatrixPsi4(poltype,outputname,molecule)
    poltype.espmethod=curespmethod
    poltype.espbasisset=curspbasisset

    return WBOmatrix,outputname

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
        highlightbonds=[] 
        mols=[]
        indexes=FirstPassAtomIndexes(poltype,tor)
        fragfoldername=str(tor[1])+'_'+str(tor[2])+'_Hydrated'
        if not os.path.isdir(fragfoldername):
            os.mkdir(fragfoldername)
        os.chdir(fragfoldername)
        
        origfragmol,parentindextofragindex,fragidxswithchangedval,fragidxtooriginalformalcharge=GenerateFrag(poltype,indexes,valtopo,mol)
        mols.append(origfragmol)
        fragments=[]
        fragments.append(origfragmol)
        print('fragfoldername',fragfoldername)
        charge=Chem.rdmolops.GetFormalCharge(molecule)
        print('molecule formal charge before hydrate',charge)

        for atom in origfragmol.GetAtoms():
            atomidx=atom.GetIdx()
            atomcharge=atom.GetFormalCharge()
            print('atom idx before hydrate ',atomidx,'charge ',atomcharge)
        fragmol=HydrateMolecule(poltype,origfragmol,fragfoldername,bondtopo,bondsettobondorder,parentindextofragindex,mol,fragidxswithchangedval,fragidxtooriginalformalcharge)


        mols.append(fragmol)
        filename=fragfoldername+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        hydindexes=FindAddedHydrogenIndexes(poltype,mols)
        print('hydindexes',hydindexes)
        combindexlist=CombinationsHydIndexes(poltype,hydindexes)
        fragments.append(fragmol) # append case where all H no CH3 added
        fragments=MethylateCombinations(poltype,combindexlist,filename,fragments)
        os.chdir('..')
        fragmoltoWBOmatrices={}
        fragmoltofragfoldername={}
        fragmoltobondindexlist={}
        for i in range(len(fragments)):
            fragmol=fragments[i]
            fragfoldername=str(tor[1])+'_'+str(tor[2])+'_Index'+'_'+str(i)
            print('fragfoldername',fragfoldername)

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
            fragWBOmatrix,outputname=GenerateWBOMatrix(poltype,fragmol,filename.replace('.mol','_xyzformat.xyz'))
            fragmentWBOvalue=fragWBOmatrix[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]] # rdkit is 0 index based so need to subtract 1, babel is 1 indexbased
            parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1] # Matrix has 0,0 so need to subtract 1 from babel index
            WBOdifference=round(numpy.abs(fragmentWBOvalue-parentWBOvalue),3)
            rotbndindextoWBOdifference[rotbndidx]=WBOdifference
            fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist)

            os.chdir('..')
            WBOdifferencetofragWBOmatrix[WBOdifference]=fragWBOmatrix
            WBOdifferencetofoldername[WBOdifference]=fragfoldername
        WBOdifference=min(list(WBOdifferencetofragWBOmatrix))
        fragWBOmatrix=WBOdifferencetofragWBOmatrix[WBOdifference]
        fragfoldername=WBOdifferencetofoldername[WBOdifference]
        rotbndindextofragfoldername[rotbndidx]=fragfoldername
        os.chdir(fragfoldername)
        fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        highlightbonds.append(fragrotbndidx)
        fragpath=os.getcwd() 
        grow=False
        fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist)
        Draw2DMolculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist)
        if WBOdifference<=poltype.WBOtol: # then we consider the fragment good enough to transfer torsion parameters, so make this fragment into .sdf file
            pass
        else:
            grow=True
            fragmol,newindexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex,fragpath=GrowFragmentOut(poltype,mol,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername,bondtopo,bondsettobondorder,valtopo)
            fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist)
        

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
                Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',structfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)
                Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',structfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)
        os.chdir(curdir)
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

def CombinationsHydIndexes(poltype,hydindexes):
    combindexlist=[]
    for i in range(len(hydindexes)):
        comb = combinations(hydindexes, i+1)
        combindexlist.append(comb)
    return combindexlist 


def MethylateCombinations(poltype,combindexlist,molfilename,fragments):
    fragmentsmarts=[]
    if combindexlist==[]:
        return fragments
    j=0
    for comb in combindexlist:
        i=0
        for idxlist in comb:
            fragmolbabel=ReadMolFileToOBMol(poltype,molfilename)
            for idx in idxlist:
                babelidx=idx+1
                atom=fragmolbabel.GetAtom(babelidx)
                convert=atom.HtoMethyl()
                if convert==False:
                    raise ValueError('Error! could not convert H to methyl '+str(babelidx)+' '+molfilename)
                    
            outputname=molfilename.replace('.mol','_frombabel_'+str(j)+'_'+str(i)+'.mol')
            WriteOBMolToMol(poltype,fragmolbabel,outputname)
            fragmolrdkit=rdmolfiles.MolFromMolFile(outputname,removeHs=False)
            fragsmarts=rdmolfiles.MolToSmarts(fragmolrdkit)
            if fragsmarts not in fragmentsmarts:
                fragmentsmarts.append(fragsmarts)
                fragments.append(fragmolrdkit)
            i+=1
        j+=1
    return fragments
                

def FindAddedHydrogenIndexes(poltype,mols):
    hydindexes=[]
    hydratedmol=mols[1]
    originalmol=mols[0]
    smarts=rdmolfiles.MolToSmarts(originalmol)
    matches = hydratedmol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    firstmatch=matches[0]
    for atom in hydratedmol.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx not in firstmatch:
            hydindexes.append(atomidx)
    if len(hydindexes)==0:
        print('smarts',smarts,'firstmatch',firstmatch)
    return hydindexes

def GrowFragmentOut(poltype,mol,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername,bondtopo,bondsettobondorder,valtopo):
    fragfoldernamepath=os.getcwd()
    while not WBOdifference<=poltype.WBOtol:
        WBOdiffarray=[]
        molarray=[]
        fragmolidxtoparentindextofragindex={}
        fragmentidxtostructfname={}
        fragmolidxtofoldername={}
        fragmolidxtofragmol={}
        parentindextofragindexlist=[]
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
            fragments=[]
            fragments.append(origfragmol)
            mols.append(origfragmol)
            fragmol=HydrateMolecule(poltype,origfragmol,basename,bondtopo,bondsettobondorder,parentindextofragindex,mol,fragidxswithchangedval,fragidxtooriginalformalcharge)
            mols.append(fragmol)
            filename=basename+'.mol'
            WriteRdkitMolToMolFile(poltype,fragmol,filename)
            fragmolbabel=ReadMolFileToOBMol(poltype,filename)
            WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
            WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
            hydindexes=FindAddedHydrogenIndexes(poltype,mols)
            combindexlist=CombinationsHydIndexes(poltype,hydindexes)
            
            fragments.append(fragmol) # include the case where all H and no H converted to CH3
            fragments=MethylateCombinations(poltype,combindexlist,filename,fragments)
            os.chdir('..')
            for i in range(len(fragments)):
                frag=fragments[i]
                fragmentlist.append(frag)
                parentindextofragindexlist.append(parentindextofragindex)
        fragmoltoWBOmatrices={}
        fragmoltofragfoldername={}
        fragmoltobondindexlist={}
        for fragmolidx in range(len(fragmentlist)):
            fragmol=fragmentlist[fragmolidx]
            basename=fragfoldername+'_GrowFragment_'+str(fragmolidx)
            if not os.path .isdir(basename):
                os.mkdir(basename)
            os.chdir(basename)
            fragmoltofragfoldername[fragmol]=basname
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

            fragWBOmatrix,outputname=GenerateWBOMatrix(poltype,fragmol,structfname)
            reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
            relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
            fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
            fragmentWBOvalue=fragWBOmatrix[fragrotbndidx[0],fragrotbndidx[1]]
            parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1]
            WBOdifference=round(numpy.abs(fragmentWBOvalue-parentWBOvalue),3)
            fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist)

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
        Draw2DMolculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist)
        WBOdifference=min(WBOdifftoindexlist.keys())
        parentindextofragindex=WBOdifftoparentindextofragindex[WBOdifference]
        indexes=WBOdifftoindexlist[WBOdifference]
        foldername=WBOdifftofolder[WBOdifference]
        structfname=WBOdifftostructfname[WBOdifference]
        RemoveTempFolders(poltype)
        os.chdir(foldername)

        fragmol=WBOdifftofragmol[WBOdifference]
        fragWBOmatrix=WBOdifftofragWBOmatrix[WBOdifference]
        fragpath=os.getcwd()
    curdir=os.getcwd()
    os.chdir(fragfoldernamepath)
    PlotFragmenterResults(poltype,WBOdiffarray,molarray)
    os.chdir(curdir)
    return fragmol,indexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex,fragpath


def GrowPossibleFragmentAtomIndexes(poltype,rdkitmol,indexes):
    possiblefragatmidxs=[]
    comblist=[]
    for bond in rdkitmol.GetBonds():
        aidx=bond.GetBeginAtomIdx()
        bidx=bond.GetEndAtomIdx()
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
           aromaticappendages=GrabAromaticAppendages(poltype,aromaticindexes,indexlist)
           newindexes=aromaticindexes+aromaticappendages
           for atmidx in newindexes:
               if atmidx not in indexlist:
                   indexlist.append(atmidx)
        if indexlist not in possiblefragatmidxs:
           possiblefragatmidxs.append(indexlist)

    return possiblefragatmidxs
    

def WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor,fragmoltoWBOmatrices,fragmoltobondindexlist):
    highlightbonds=[]
    structfnamemol=fragfoldername+'.mol'
    print(Chem.MolToMolBlock(fragmol,kekulize=False),file=open(structfnamemol,'w+')) 
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
    Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',structfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)
    Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',structfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)
    temp=[reducedparentWBOmatrix,fragWBOmatrix]
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

def RemoveRadicalsSDF(poltype,structfname):
    temp=open(structfname,'r')
    results=temp.readlines()
    temp.close()
    tempname=structfname.replace('.sdf','_temp.sdf')
    temp=open(tempname,'w')
    for line in results:
        if 'RAD' not in line:
            temp.write(line)
    temp.close()
    os.remove(structfname)
    os.rename(tempname,structfname)

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
                       aromaticappendages=GrabAromaticAppendages(poltype,aromaticindexes,[i-1 for i in tor])
                       newindexes=aromaticindexes+aromaticappendages
                       for atmidx in newindexes:
                           if atmidx not in molindexlist:
                               molindexlist.append(atmidx)
   return molindexlist

def Draw2DMolculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist):
    img=Chem.Draw.MolsToGridImage(fragments,molsPerRow=4,subImgSize=(200,200),legends=[fragmoltofragfoldername[frag] for frag in fragments],highlightBondLists=[fragmoltobondindexlist[frag] for frag in fragments])
    img.save('Fragments.svg')

def Draw2DMoleculeWithWBO(poltype,WBOmatrix,basename,structure,bondindexlist=None,smirks=None):
    mol=rdmolfiles.MolFromMolFile(structure,sanitize=True,removeHs=False)
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
            WBOval=WBOmatrix[begidx,endidx]
            if WBOval==0:
                continue
            wbo=str(round(WBOval,4))
            label = sg.TextElement(bondcoords[0],bondcoords[1], wbo, size=12, weight="bold")
            if endatomdrawcoords[1]>begatomdrawcoords[1]:
                array=endatomdrawcoords-begatomdrawcoords  
            else:
                array=begatomdrawcoords-endatomdrawcoords
            norm = numpy.linalg.norm(array)
            normarray=array/norm
            angle=numpy.degrees(numpy.arccos(normarray[1]))
            label.rotate(angle,bondcoords[0],bondcoords[1])
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

def GrabAromaticAppendages(poltype,aromaticindexes,indexlist):
    aromaticandappendages=[]
    for atomidx in aromaticindexes:
        appendageindexes=CheckForNeighboringFunctionalGroups(poltype,atomidx,indexlist,aromaticindexes)
        for index in appendageindexes:
            if index not in aromaticandappendages:
                aromaticandappendages.append(index)
    return aromaticandappendages



def PlotFragmenterResults(poltype,WBOdiffarray,molarray):
    fig=plt.figure()
    basename='NumberofAtomsVSWBODifference'
    plt.plot(WBOdiffarray,[m.GetNumAtoms() for m in molarray],'.') 
    plt.xlabel('WBO Difference')
    plt.ylabel('Number of atoms in fragment')
    fig.savefig(basename+'.png')


def CheckForNeighboringFunctionalGroups(poltype,atomidx,indexlist,aromaticindexes):
    appendage=[]
    maxcount=3
    count=0
    while count!=maxcount:
        prevlength=len(appendage)
        atom=poltype.rdkitmol.GetAtomWithIdx(atomidx)
        for natom in atom.GetNeighbors():
            nidx=natom.GetIdx()
            nval=natom.GetExplicitValence()
            if nidx not in indexlist and nidx not in aromaticindexes and nidx not in appendage:
                appendage.append(nidx)   
                if nval!=1:
                    atomidx=nidx
        count+=1
        length=len(appendage)
        if length==prevlength:
            break
    return appendage   
