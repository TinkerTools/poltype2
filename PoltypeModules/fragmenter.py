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
    print('cwd',os.getcwd(),flush=True)
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
                
def GrabWBOMatrixPsi4(poltype,outputlog,mol):
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
                           

 
def AllIntegers(poltype,testlist):
    allintegers=True
    for value in testlist:
        if not value.isdigit():
            allintegers=False
    return allintegers

def FindEquivalentFragments(poltype,fragmentarray):
    equivalentfragmentsarray=[]
    atomnumberarray=[len(m.GetAtoms()) for m in fragmentarray]
    D = defaultdict(list)
    for i,item in enumerate(atomnumberarray):
        D[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1} 
    for atomnum in D.keys():
        indices = D[atomnum]
        temp=[]
        for index in indices:
            fragmol=fragmentarray[index]
            temp.append(fragmol)
        equivalentfragmentsarray.append(temp)
    return equivalentfragmentsarray

         
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
    poltypeinput={'suppressdipoleerr':'True','optmethod':poltype.optmethod,'toroptmethod':poltype.toroptmethod,'espmethod':poltype.espmethod,'torspmethod':poltype.torspmethod,'dmamethod':poltype.dmamethod,'torspbasisset':poltype.torspbasisset,'espbasisset':poltype.espbasisset,'dmabasisset':poltype.dmabasisset,'toroptbasisset':poltype.toroptbasisset,'optbasisset':poltype.optbasisset,'fitrotbndslist':strfragrotbndindexes,'bashrcpath':poltype.bashrcpath,'externalapi':poltype.externalapi,'use_gaus':poltype.use_gaus,'use_gausoptonly':poltype.use_gausoptonly,'isfragjob':True,'poltypepath':poltype.poltypepath,'structure':tail,'numproc':poltype.numproc,'maxmem':poltype.maxmem,'maxdisk':poltype.maxdisk,'printoutput':True}
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
        finishedjobs,errorjobs=poltype.CallJobsLocalHost(listofjobs,jobtooutputlog)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(listofjobs,jobtooutputlog)

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
            logtoconvertidxs[newlog]=[fragmol,tail,wholexyz,wholemol]
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
            logtoconvertidxs[newlog]=[fragmol,tail,wholexyz,wholemol]
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
                logtoconvertidxs[newlog]=[fragmol,tail,wholexyz,wholemol]
    os.chdir(parentdir)
    finishedjobs,errorjobs=SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog)
    for log in finishedjobs:
        [fragmol,tail,wholexyz,wholemol]=logtoconvertidxs[log]
        ConvertFragIdxToWholeIdx(poltype,torlist,rotbndindextoparentindextofragindex,fragmol,tail,wholexyz,wholemol)
    
def MakeTorsionFileName(poltype,string):
    temp=open('torsions.txt','w')
    temp.write(string+'\n')
    temp.close()


def WriteDictionaryToFile(poltype,dictionary):
    json.dump(dictionary, open("parentindextofragindex.txt",'w'))


def DeleteEquivalentFragments(poltype,equivalentfragmentstodelete):
    curdir=os.getcwd()
    os.chdir('..')
    for fold in equivalentfragmentstodelete:
        if os.path.isdir(fold):
            shutil.rmtree(fold)
    os.chdir(curdir)

def ConvertFragIdxToWholeIdx(poltype,torlist,rotbndindextoparentindextofragindex,fragmol,fragmentfilename,wholexyz,wholemol):
    symmtorlist=[]
    for tor in torlist:
        classkey=torgen.get_class_key(poltype,tor[0],tor[1],tor[2],tor[3])
        if classkey not in symmtorlist:
            symmtorlist.append(classkey)
    print('torlist',torlist,flush=True)
    print('symmtorlist',symmtorlist,flush=True)

    fragmentfileprefix=fragmentfilename.replace('.sdf','')
    if not os.path.isfile(os.getcwd()+r'/'+fragmentfileprefix+'.key_5'):
        return # if POLTYPE job failed, just try to submit the other fragment jobs instead of killing parent job
    temp=open(os.getcwd()+r'/'+fragmentfileprefix+'.key_5','r')
    fragkeyresults=temp.readlines()
    temp.close()
    temp=open(wholexyz,'r')
    wholetttxyzresults=temp.readlines()
    temp.close()
    temp=open(os.getcwd()+r'/'+'ttt.xyz','r')
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

    m=mol_with_atom_index(poltype,fragmol)
    fragsmirks=rdmolfiles.MolToSmarts(m) 
    fragidxarray=GrabAtomOrder(poltype,fragsmirks)
    classkeytosmilesposarray={}
    for tor in torlist:
        rotbndindex=str(tor[1])+'_'+str(tor[2])
        parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
        smilesposarray=[]
        classkey=torgen.get_class_key(poltype,tor[0],tor[1],tor[2],tor[3])
        for index in tor:
            fragindex=parentindextofragindex[index]
            fragidxarraypos=fragidxarray.index(fragindex)
            smilespos=fragidxarraypos+1
            smilesposarray.append(smilespos)
        classkeytosmilesposarray[classkey]=smilesposarray
    fragsmarts=rdmolfiles.MolToSmarts(fragmol)
    p = Chem.MolFromSmarts(fragsmarts)
    diditmatchrdkit=fragmol.HasSubstructMatch(p)
    

    indexlistoflists=[]
    for indexls in fragmol.GetSubstructMatches(p):
        indexlistoflists.append(indexls)
    wholemolindexlist=[]
    for index in indexlistoflists[0]:
         babindex=index+1
         wholemolindexlist.append(babindex)
    fragidxtowholemolidx=dict(zip(fragidxarray, wholemolindexlist))
    fragtypeidxtowholemoltypeidx={}
    for fragidx in fragidxtowholemolidx.keys():
        wholemolidx=fragidxtowholemolidx[fragidx]
        fragtypeidx=fragidxtotypeidx[fragidx]
        wholemoltypeidx=wholeidxtotypeidx[wholemolidx]
        fragtypeidxtowholemoltypeidx[fragtypeidx]=wholemoltypeidx


    newtemp=open('valence.prms','w')
    temp=open(fragmentfileprefix+'.key_6','w')
    for line in fragkeyresults:
        if 'torsion' in line:
            linesplit=line.split()
            typea=int(linesplit[1])
            typeb=int(linesplit[2])
            typec=int(linesplit[3])
            typed=int(linesplit[4])
            tor=[typea,typeb,typec,typed]
            torkey='%d %d %d %d' % (typea, typeb, typec, typed)
            print('torkey',torkey,flush=True)
            if torkey in symmtorlist:
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


def GenerateFrag(poltype,molindexlist):
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
    return res,oldindextonewindex

 
def HydrateMolecule(poltype,mol,basename):
    hydmol=Chem.rdmolops.AddHs(mol)
    AllChem.EmbedMolecule(hydmol)
    return hydmol

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

def ReadMolFileToOBMol(poltype,filename):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('mol')
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
    inputname,outputname=esp.CreatePsi4ESPInputFile(poltype,structfname,poltype.comespfname.replace('.com','_frag.com'),molecule,poltype.maxdisk,poltype.maxmem,poltype.numproc,False)
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

    fragpath=os.getcwd()
    rotbndindextoparentindextofragindex={}
    rotbndindextofragment={}
    rotbndindextofragmentfilepath={}
    rotbndindextofragWBOmatrix={}
    rotbndindextofragfoldername={}
    for tor in torlist: 
        highlightbonds=[] 
        mols=[]
        indexes=FirstPassAtomIndexes(poltype,tor)
        fragfoldername=str(tor[1])+'_'+str(tor[2])+'_Hydrated'
        
        if not os.path.isdir(fragfoldername):
            os.mkdir(fragfoldername)
        os.chdir(fragfoldername)
        fragmol,parentindextofragindex=GenerateFrag(poltype,indexes)
        mols.append(fragmol)
        fragmol=HydrateMolecule(poltype,fragmol,fragfoldername)
        mols.append(fragmol)
        filename=fragfoldername+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        hydindexes=FindAddedHydrogenIndexes(poltype,mols)
        combindexlist=CombinationsHydIndexes(poltype,hydindexes)
        fragments=MethylateCombinations(poltype,combindexlist,filename)
        fragments.append(fragmol) # append case where all H no CH3 added
        WBOdifferencetofragWBOmatrix={}
        WBOdifferencetofoldername={}
        os.chdir('..')
        for i in range(len(fragments)):
            fragmol=fragments[i]
            fragfoldername=str(tor[1])+'_'+str(tor[2])+'_Index'+'_'+str(i)
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
            WBOdifference=numpy.abs(fragmentWBOvalue-parentWBOvalue)
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
        if WBOdifference<=poltype.WBOtol: # then we consider the fragment good enough to transfer torsion parameters, so make this fragment into .sdf file
            pass
        else:
            fragmol,newindexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex=GrowFragmentOut(poltype,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername)
        WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor)
        structfname=structfname.replace('_xyzformat.xyz','.sdf')
        rotbndindextofragment[rotbndidx]=fragmol
        rotbndindextofragmentfilepath[rotbndidx]=os.getcwd()+r'/'+structfname
        rotbndindextoparentindextofragindex[rotbndidx]=parentindextofragindex
        rotbndindextofragWBOmatrix[rotbndidx]=fragWBOmatrix
        os.chdir(fragpath)
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
    for rotbndindex,fragmol in rotbndindextofragment.items():
        highlightbonds=[]
        fragWBOmatrix=rotbndindextofragWBOmatrix[rotbndindex]
        fragfilepath=rotbndindextofragmentfilepath[rotbndindex]
        head,tail=os.path.split(fragfilepath)
        os.chdir(head)
        fragfoldername=rotbndindextofragfoldername[rotbndindex]
        structfnamemol=fragfoldername+'.mol'
        parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
        equivalentrotbndindexes=GrabEquivalentRotBndIndexes(poltype,equivalentrotbndindexarrays,rotbndindex)
        if len(equivalentrotbndindexes)!=0:
            for bndindexes in equivalentrotbndindexes:
                for bndindex in bndindexes:
                    indexes=bndindex.split('_')
                    indexes=[int(i) for i in indexes]
                    fragmentWBOvalue=fragWBOmatrix[parentindextofragindex[indexes[0]-1],parentindextofragindex[indexes[1]-1]] # rdkit is 0 index based so need to subtract 1, babel is 1 indexbased
                    parentWBOvalue=parentWBOmatrix[indexes[0]-1,indexes[1]-1] # Matrix has 0,0 so need to subtract 1 from babel index
                    WBOdifference=numpy.abs(fragmentWBOvalue-parentWBOvalue)
                    basename=fragfoldername+'_WBO_'+str(round(WBOdifference,2))
                    fragrotbndidx=[parentindextofragindex[indexes[0]-1],parentindextofragindex[indexes[1]-1]]
                    highlightbonds.append(fragrotbndidx)
            reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
            relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
            m=mol_with_atom_index(poltype,fragmol)
            fragsmirks=rdmolfiles.MolToSmarts(m) 
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
    if hydindexes==None:
        return combindexlist
    for i in range(len(hydindexes)):
        comb = combinations(hydindexes, i+1)
        combindexlist.append(comb)
    return combindexlist 


def MethylateCombinations(poltype,combindexlist,molfilename):
    fragments=[]
    if combindexlist==[]:
        return fragments
    for comb in combindexlist:
        fragmolbabel=ReadMolFileToOBMol(poltype,molfilename)
        for idxlist in comb:
            for idx in idxlist:
                babelidx=idx+1
                atom=fragmolbabel.GetAtom(babelidx)
                convert=atom.HtoMethyl()
                if convert==False:
                    raise ValueError('Error! could not convert H to methyl '+str(babelidx)+' '+molfilename)
                    
        outputname=molfilename.replace('.mol','_frombabel.mol')
        WriteOBMolToMol(poltype,fragmolbabel,outputname)
        fragmolrdkit=rdmolfiles.MolFromMolFile(outputname)
        fragments.append(fragmolrdkit)
    return fragments
                

def FindAddedHydrogenIndexes(poltype,mols):
    hydindexes=[]
    res=rdFMCS.FindMCS(mols)
    smarts=res.smartsString    
    hydratedmol=mols[1]
    originalmol=mols[0]
    matches = hydratedmol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    firstmatch=matches[0]
    for atom in hydratedmol.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx not in firstmatch:
            hydindexes.append(atomidx)

def GrowFragmentOut(poltype,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername):
    fragfoldernamepath=os.getcwd()
    WBOdiffarray=[]
    molarray=[]
    fragmolidxtoparentindextofragindex={}
    fragmentidxtostructfname={}
    fragmolidxtofoldername={}
    fragmolidxtofragmol={}
    possiblefragatmidxs=GrowPossibleFragmentAtomIndexes(poltype,poltype.rdkitmol,indexes)
    fragmentlist=[]
    parentindextofragindexlist=[]
    for fragmolidx in range(len(possiblefragatmidxs)):
        indexlist=possiblefragatmidxs[fragmolidx]
        basename=fragfoldername+'_GrowFragment_Hydrated'+'_'+str(fragmolidx) 
        if not os.path .isdir(basename):
            os.mkdir(basename)
        os.chdir(basename)
        mols=[]
        fragmol,parentindextofragindex=GenerateFrag(poltype,indexlist)
        mols.append(fragmol)
        fragmol=HydrateMolecule(poltype,fragmol,basename)
        mols.append(fragmol)
        filename=basename+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        fragmolbabel=ReadMolFileToOBMol(poltype,filename)
        WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
        WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))

        hydindexes=FindAddedHydrogenIndexes(poltype,mols)
        combindexlist=CombinationsHydIndexes(poltype,hydindexes)
        fragments=MethylateCombinations(poltype,combindexlist,filename)
        fragments.append(fragmol) # include the case where all H are not converted to CH3
        os.chdir('..')
        for i in range(len(fragments)):
            frag=fragments[i]
            fragmentlist.append(frag)
            parentindextofragindexlist.append(parentindextofragindex)
    for fragmolidx in range(len(fragmentlist)):
        fragmol=fragmentlist[fragmolidx]
        basename=fragfoldername+'_GrowFragment_'+str(fragmolidx)
        if not os.path .isdir(basename):
            os.mkdir(basename)
        os.chdir(basename)
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
        highlightbonds=[]
        fragmol=fragmolidxtofragmol[fragmolidx]
        foldername=fragmolidxtofoldername[fragmolidx]
        parentindextofragindex=fragmolidxtoparentindextofragindex[fragmolidx]
        structfname=fragmentidxtostructfname[fragmolidx]
        os.chdir(foldername)
        fragWBOmatrix,outputname=GenerateWBOMatrix(poltype,fragmol,structfname)
        reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
        relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
        fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        highlightbonds.append(fragrotbndidx)
        fragmentWBOvalue=fragWBOmatrix[fragrotbndidx[0],fragrotbndidx[1]]
        parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1]
        WBOdifference=numpy.abs(fragmentWBOvalue-parentWBOvalue)
        basename=fragfoldername+'_WBO_'+str(round(WBOdifference,2))
        m=mol_with_atom_index(poltype,fragmol)
        fragsmirks=rdmolfiles.MolToSmarts(m) 
        Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',structfname.replace('_xyzformat.xyz','.mol'),bondindexlist=highlightbonds,smirks=fragsmirks)
        Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',structfname.replace('_xyzformat.xyz','.mol'),bondindexlist=highlightbonds,smirks=fragsmirks)
        structfnamemol=basename+'.mol'
        print(Chem.MolToMolBlock(fragmol,kekulize=False),file=open(structfnamemol,'w+'))
        os.chdir('..')
        indexlist=parentindextofragindex.values()
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
    indexlist=WBOdifftoindexlist[WBOdifference]
    foldername=WBOdifftofolder[WBOdifference]
    structfname=WBOdifftostructfname[WBOdifference]
    os.chdir(foldername)
    fragmol=WBOdifftofragmol[WBOdifference]
    fragWBOmatrix=WBOdifftofragWBOmatrix[WBOdifference]
    curdir=os.getcwd()
    os.chdir(fragfoldernamepath)
    PlotFragmenterResults(poltype,WBOdiffarray,molarray)
    RemoveTempFolders(poltype)
    os.chdir(curdir)
    return fragmol,indexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex


def GrowPossibleFragmentAtomIndexes(poltype,rdkitmol,indexes):
    possiblefragatmidxs=[]
    for bond in rdkitmol.GetBonds():
        aidx=bond.GetBeginAtomIdx()
        bidx=bond.GetEndAtomIdx()
        if (aidx in indexes and bidx not in indexes): # then this means the bond is not already in the fragment but this is one of the bonds just outside of the fragment
            indexlist=indexes.copy()
            aromaticindexes=GrabRingAtoms(poltype,rdkitmol.GetAtomWithIdx(aidx))
            for atmidx in aromaticindexes:
                if atmidx not in indexlist:
                    indexlist.append(atmidx)
            if indexlist not in possiblefragatmidxs:
                possiblefragatmidxs.append(indexlist)

        elif (aidx not in indexes and bidx in indexes):
            indexlist=indexes.copy()
            aromaticindexes=GrabRingAtoms(poltype,rdkitmol.GetAtomWithIdx(bidx))
            for atmidx in aromaticindexes:
                if atmidx not in indexlist:
                    indexlist.append(atmidx)
            if indexlist not in possiblefragatmidxs:
                possiblefragatmidxs.append(indexlist)

    return possiblefragatmidxs
    

def WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,tor):
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
    basename=fragfoldername+'_WBO_'+str(round(WBOdifference,2))
    fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
    highlightbonds.append(fragrotbndidx)
    reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
    relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
    m=mol_with_atom_index(poltype,fragmol)
    fragsmirks=rdmolfiles.MolToSmarts(m) 
    Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',structfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)
    Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',structfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)
        

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

def RemoveRadicals(poltype,structfname):
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
                   if neighbatom.GetIsAromatic()==True or IsRingAtom(poltype,neighbatom)==True:
                       ringindexes=GrabRingAtoms(poltype,neighbatom)
                       ringandappendages=GrabRingAppendages(poltype,ringindexes)
                       for atmidx in ringandappendages:
                           if atmidx not in molindexlist:
                               molindexlist.append(atmidx)
   return molindexlist

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
            wbo=str(round(WBOval,3))
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
        label = sg.TextElement(250,375, smirks, size=12, weight="bold")
        fig.append(label)
    fig.save(basename+'.svg')
    svg_code=fig.to_str()
    svg2png(bytestring=svg_code,write_to=basename+'.png')

def GrabRingAtoms(poltype,neighbatom):
    ringindexes=[]
    prevringidxlen=len(ringindexes)
    ringindexes.append(neighbatom.GetIdx())
    ringidxlen=len(ringindexes)
    while prevringidxlen!=ringidxlen:
        for atmindex in ringindexes:
            atm=poltype.rdkitmol.GetAtomWithIdx(atmindex)
            if (atm.GetIsAromatic()==True or IsRingAtom(poltype,atm)==True) and atmindex not in ringindexes:
                ringindexes.append(atmindex)
            for natm in atm.GetNeighbors():
                if (natm.GetIsAromatic()==True or IsRingAtom(poltype,natm)==True) and natm.GetIdx() not in ringindexes:
                    ringindexes.append(natm.GetIdx())
        prevringidxlen=ringidxlen
        ringidxlen=len(ringindexes)

    return ringindexes

def IsRingAtom(poltype,atm):
    ringinfo=poltype.rdkitmol.GetRingInfo()
    ringsize=[3,4,5,6,7]
    isringatom=False
    for size in ringsize:
        isit=ringinfo.IsAtomInRingOfSize(atm.GetIdx(),3)
        if isit==True:
            isringatom=True
    return isringatom

def GrabRingAppendages(poltype,ringindexes):
    ringandappendages=[]
    for atomidx in ringindexes:
        appendageindexes=CheckForNeighboringFunctionalGroups(poltype,atomidx)   
        for index in appendageindexes:
            if index not in ringandappendages:
                ringandappendages.append(index)
    return ringandappendages


def CheckForNeighboringFunctionalGroups(poltype,atomidx):
    appendage=[atomidx]
    functionalgroups={'Acetylenic Carbon':'[$([CX2]#C)]','Alkyl Carbon':'[CX4]','Allenic Carbon':'[CX2](=C)=C','Carbonyl Group':'[CX3]=[OX1]','Carbonyl Group with Carbon':'[CX3](=[OX1])C','Carbonyl with Nitrogen':'[OX1]=CN','Carbonyl with Oxygen':'[CX3](=[OX1])O','Acyl Halide':'[CX3](=[OX1])[F,Cl,Br,I]','Aldehyde':'[CX3H1](=O)[#6]','Anhydride':'[CX3](=[OX1])[OX2][CX3](=[OX1])','Amide':'[NX3][CX3](=[OX1])[#6]','Amidinium':'[NX3][CX3]=[NX3+]','Carbamate':'[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]','Carbamic ester':'[NX3][CX3](=[OX1])[OX2H0]','Carbamic acid':'[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]','Carboxylate Ion':'[CX3](=O)[O-]','Carbonic Acid or Carbonic Ester':'[CX3](=[OX1])(O)O','Carbonic Acid or Carbonic Acid-Ester':'[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]','Carbonic Ester (carbonic acid diester)':'C[OX2][CX3](=[OX1])[OX2]C','Carboxylic acid':'[CX3](=O)[OX2H1]','Carboxylic acid or conjugate base':'[CX3](=O)[OX1H0-,OX2H1]','Cyanamide':'[NX3][CX2]#[NX1]','Ester Also hits anhydrides':'[#6][CX3](=O)[OX2H0][#6]','Ketone':'[#6][CX3](=O)[#6]','Ether':'[OD2]([#6])[#6]','Enamine':'[NX3][CX3]=[CX3]','Amine':'[NX3;H2,H1;!$(NC=O)]','Nitrogen':'[#7]','Phosphoric Ester':'[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]','Phosphoric Acid':'[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]','Azo Nitrogen':'[NX2]=N','Azo Nitrogen':'[NX2]=[NX2]','Nitrile':'[NX1]#[CX2]','Isonitrile':'[CX1-]#[NX2+]','Nitroso-group':'[NX2]=[OX1]','Hydroxyl':'[OX2H]','Hydroxyl in Alcohol':'[#6][OX2H]',' Hydroxyl in Carboxylic Acid':'[OX2H][CX3]=[OX1]','Hydroxyl in H-O-P-':'[OX2H]P','Enol':'[OX2H][#6X3]=[#6]','Phenol':'[OX2H][cX3]:[c]','Enol or Phenol':'[OX2H][$(C=C),$(cc)]','Hydroxyl_acidic':'[$([OH]-*=[!#6])]','Peroxide groups':'[OX2,OX1-][OX2,OX1-]','Carbo-Thiocarboxylate':'[S-][CX3](=S)[#6]','Carbo-Thioester':'S([#6])[CX3](=O)[#6]','Thiol, Sulfide or Disulfide Sulfur':'[SX2]','Thiol':'[#16X2H]','Sulfur with at-least one hydrogen':'[#16!H0]','Thioamide':'[NX3][CX3]=[SX1]','Sulfide':'[#16X2H0]','Mono-sulfide':'[#16X2H0][!#16]','Di-sulfide':'[#16X2H0][#16X2H0]','Any carbon attached to any halogen':'[#6][F,Cl,Br,I]','Halogen':'[F,Cl,Br,I]','Acyl Halide':'[CX3](=[OX1])[F,Cl,Br,I]'} # just include all functional groups, some functional groups are contained within others but most specific cases will always take precedence over less specific cases
    for funcname in functionalgroups.keys():
        SMART=functionalgroups[funcname]
        patt = Chem.MolFromSmarts(SMART)
        match=poltype.rdkitmol.HasSubstructMatch(patt)
        if match:
            matches=poltype.rdkitmol.GetSubstructMatch(patt)
            if atomidx in matches:
                for index in matches:
                    if index not in appendage:
                        appendage.append(index)
    atom=poltype.rdkitmol.GetAtomWithIdx(atomidx)
    for natom in atom.GetNeighbors():
        atomicnum=natom.GetAtomicNum()
        if atomicnum==1:
            natomidx=natom.GetIdx()
            if natomidx not in appendage:
                appendage.append(natomidx) 
    
    return appendage   


def PlotFragmenterResults(poltype,WBOdiffarray,molarray):
    fig=plt.figure()
    basename='NumberofAtomsVSWBODifference'
    plt.plot(WBOdiffarray,[m.GetNumAtoms() for m in molarray],'.') 
    plt.xlabel('WBO Difference')
    plt.ylabel('Number of atoms in fragment')
    fig.savefig(basename+'.png')

    fig=plt.figure()
    basename='NumberofRotatableBondsVSWBODifference'
    plt.plot(WBOdiffarray,[rdMolDescriptors.CalcNumRotatableBonds(m) for m in molarray],'.')
    plt.xlabel('WBO Difference')
    plt.ylabel('Number of rotatable bonds in fragment')

    fig.savefig(basename+'.png')

