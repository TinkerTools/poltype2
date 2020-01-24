from PoltypeModules import electrostaticpotential as esp
from PoltypeModules import torsiongenerator as torgen
from poltype import PolarizableTyper
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


def GrabTorsionParametersFromFragments(poltype,torlist):
    valenceprmlist=[]
    symmtorlist=[]
    for tor in torlist:
        symmtorlist.append(torgen.get_class_key(poltype,tor[0],tor[1],tor[2],tor[3]))
   
    os.chdir('Fragments')
    torprmdic={}
    foldlist=os.listdir(os.getcwd())
    for f in foldlist:
        if os.path.isdir(f):
            os.chdir(f)
            filelist=os.listdir(os.getcwd())
            for ff in filelist:
                if '.prm' in ff:
                    temp=open(ff,'r')
                    results=temp.readlines()
                    temp.close()
                    for line in results:
                        valenceprmlist.append(line)
                if '.key_5' in ff:
                    fragname=f.split('.')[0]
                    key6=fragname+'.key_6'
                    while not os.path.isfile(key6):
                        poltype.logfh.write('Waiting for final torsion parameters for '+key6+' in '+os.getcwd()+' waiting 5 min'+'\n')
                        time.sleep(5)
                    temp=open(key6,'r')
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
                            
            os.chdir('..')
    os.chdir('..')
    print('torprmdic ',torprmdic)
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
    WBOmatrix=numpy.empty((mol.GetNumAtoms(),mol.GetNumAtoms()))
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
        indices = [i for i, x in enumerate(atomnumberarray) if x == atomnum]
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
            rotbndindex=FindRotatableBond(poltype,fragmol,rotbndindextofragment)
            temp.append(rotbndindex)
        equivalentrotbndindexarrays.append(temp)

    return equivalentrotbndindexarrays

def FindRotatableBond(poltype,fragmol,rotbndindextofragment):
    for rotbndindex in rotbndindextofragment.keys():
        m=rotbndindextofragment[rotbndindex]
        if len(m.GetAtoms())==len(fragmol.GetAtoms()):
            return rotbndindex

 
def SpawnPoltypeJobsForFragments(poltype,rotbndindextofragindexmap,rotbndindextofragment):
    parentdir=dirname(dirname(abspath(os.getcwd())))
    print('parentdir ',parentdir)
    files=os.listdir(os.getcwd())
    fragmentarray=[]
    for f in files:
        if os.path.isdir(f):
            os.chdir(f)
            # need to grab SDF file
            fragfiles=os.listdir(os.getcwd())
            for fragfile in fragfiles:
                if '.sdf' in fragfile:
                    fraginput=fragfile

            rotbndindex=f
            fragment=rotbndindextofragment[rotbndindex]
            fragmentarray.append(fragment)
            os.chdir('..')
            break
    
    equivalentfragmentsarray=FindEquivalentFragments(poltype,fragmentarray)
    equivalentrotbndindexarrays=FindEquivalentRotatableBonds(poltype,equivalentfragmentsarray,rotbndindextofragment) 
    if equivalentfragmentsarray==[]: 
        os.chdir(rotbndindex)
        parentindextofragindex=rotbndindextofragindexmap[rotbndindex]
        rotbndindexes=rotbndindex.split('_')
        rotbndindexes=[int(i)-1 for i in rotbndindexes]
        fragrotbndindexes=[parentindextofragindex[i]+1 for i in rotbndindexes]
        strfragrotbndindexes=str(fragrotbndindexes[0])+','+str(fragrotbndindexes[1])
        PolarizableTyper(optmethod=poltype.optmethod,toroptmethod=poltype.toroptmethod,espmethod=poltype.espmethod,torspmethod=poltype.torspmethod,dmamethod=poltype.dmamethod,torspbasisset=poltype.torspbasisset,espbasisset=poltype.espbasisset,dmabasisset=poltype.dmabasisset,toroptbasisset=poltype.toroptbasisset,optbasisset=poltype.optbasisset,onlyrotbndlist=strfragrotbndindexes,bashrcpath=poltype.bashrcpath,externalapi=poltype.externalapi,use_gaus=poltype.use_gaus,use_gausoptonly=poltype.use_gausoptonly,isfragjob=True,poltypepath=poltype.poltypepath,structure=fraginput,numproc=poltype.numproc,maxmem=poltype.maxmem,maxdisk=poltype.maxdisk,wholexyz=parentdir+r'/'+poltype.xyzoutfile,wholemol=parentdir+r'/'+poltype.molstructfname,poltypeini=False,printoutput=True)
    else:
        for array in equivalentrotbndindexarrays:
            equivalentfragmentstodelete=[]
            strfragrotbndindexes=''
            for i in range(len(array)):
                rotbndindex=array[i]
                if i!=0: # just delete all except first
                    equivalentfragmentstodelete.append(rotbndindex)
                else:
                    equivalentrotbndindex=rotbndindex
                parentindextofragindex=rotbndindextofragindexmap[rotbndindex]
                rotbndindexes=rotbndindex.split('_')
                rotbndindexes=[int(i)-1 for i in rotbndindexes]
                fragrotbndindexes=[parentindextofragindex[i]+1 for i in rotbndindexes]
                strfragrotbndindexes+=str(fragrotbndindexes[0])+' '+str(fragrotbndindexes[1])+','
            strfragrotbndindexes=strfragrotbndindexes[:-1]
            DeleteEquivalentFragments(poltype,equivalentfragmentstodelete)
            os.chdir(equivalentrotbndindex)
            PolarizableTyper(optmethod=poltype.optmethod,toroptmethod=poltype.toroptmethod,espmethod=poltype.espmethod,torspmethod=poltype.torspmethod,dmamethod=poltype.dmamethod,torspbasisset=poltype.torspbasisset,espbasisset=poltype.espbasisset,dmabasisset=poltype.dmabasisset,toroptbasisset=poltype.toroptbasisset,optbasisset=poltype.optbasisset,fitrotbndslist=strfragrotbndindexes,bashrcpath=poltype.bashrcpath,externalapi=poltype.externalapi,use_gaus=poltype.use_gaus,use_gausoptonly=poltype.use_gausoptonly,isfragjob=True,poltypepath=poltype.poltypepath,structure=fraginput,numproc=poltype.numproc,maxmem=poltype.maxmem,maxdisk=poltype.maxdisk,wholexyz=parentdir+r'/'+poltype.xyzoutfile,wholemol=parentdir+r'/'+poltype.molstructfname,poltypeini=False,printoutput=True)
    os.chdir('..')



def DeleteEquivalentFragments(poltype,equivalentfragmentstodelete):
    curdir=os.getcwd()
    os.chdir('..')
    for fold in equivalentfragmentstodelete:
        if os.path.isdir(fold):
            shutil.rmtree(fold)
    os.chdir(curdir)

def ConvertFragIdxToWholeIdx(poltype,molstructfname,torlist,rotbndindextofragindexmap,parentWBOmatrix):
    symmtorlist=[]
    highlightbonds=[] 
    for tor in torlist:
        rotbndindex=str(tor[1])+'_'+str(tor[2])
        parentindextofragindex=rotbndindextofragindexmap[rotbndindex]
        rotbnd=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        highlightbonds.append(rotbnd)
        classkey=torgen.get_class_key(poltype,tor[0],tor[1],tor[2],tor[3])
        if classkey not in symmtorlist:
            symmtorlist.append(classkey)
    while not os.path.isfile(poltype.key5fname):
        logfh.write('Waiting for key 5 from fragment'+'\n')
        time.sleep(5)
    temp=open(poltype.key5fname,'r')
    fragkeyresults=temp.readlines()
    temp.close()
    while not os.path.isfile(poltype.wholexyz):
        logfh.write('Waiting for parent molecule .xyz_2 file'+'\n')
        time.sleep(5)
    temp=open(poltype.wholexyz,'r')
    wholetttxyzresults=temp.readlines()
    temp.close()
    temp=open(poltype.xyzoutfile,'r')
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

    
    fragmol=rdmolfiles.MolFromMolFile(poltype.molstructfnamemol,removeHs=False)
    fragsmarts=rdmolfiles.MolToSmarts(fragmol) # we can use this to generate the parmaters that need to go into valence.py
    m=mol_with_atom_index(poltype,fragmol)
    fragsmirks=rdmolfiles.MolToSmarts(m) # first parse the smarts to get an array corresponding to atom indexes in fragment, then remake SMARTS without indexes, match to parent
    string='_frag_psi4.log'
    files=os.listdir(os.getcwd())
    for f in files:
        if string in f:
            fragfolder=os.getcwd()
            parentindextofragindex=rotbndindextofragindexmap[fragfolder]
            fragfoldersplit=fragfolder.split()
            rotbnd=[int(i) for i in fragfoldersplit]
            fragWBOmatrix=GrabWBOMatrixPsi4(poltype,f,fragmol)
            fragmentWBOvalue=fragWBOmatrix[parentindextofragindex[rotbnd[0]-1],parentindextofragindex[rotbnd[1]-1]]
            parentWBOvalue=parentWBOmatrix[rotbnd[0]-1,rotbnd[1]-1]
            WBOdifference=numpy.abs(fragmentWBOvalue-parentWBOvalue)
            basename=poltype.molstructfname+'_'+str(round(WBOdifference,2))
            Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,os.getcwd(),poltype.molstructfnamemol,bondindexlist=highlightbonds,smirks=fragsmirks)

    fragidxarray=[]
    for eidx in range(len(fragsmirks)):
         e=fragsmirks[eidx]
         prev=fragsmirks[eidx-1]
         if e.isdigit() and prev==':':
             fragidxarray.append(int(e))
    classkeytosmilesposarray={}
    for tor in torlist:
        smilesposarray=[]
        classkey=torgen.get_class_key(poltype,tor[0],tor[1],tor[2],tor[3])
        for index in tor:
            fragidxarraypos=fragidxarray.index(index)
            smilespos=fragidxarraypos+1
            smilesposarray.append(smilespos)
        classkeytosmilesposarray[classkey]=smilesposarray
    parentmol = openbabel.OBMol()
    obConversion=openbabel.OBConversion()
    obConversion.ReadFile(parentmol, poltype.wholemol)
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,fragsmarts)
    match=sp.Match(parentmol)
    indexlistoflists=[]
    for indexls in sp.GetMapList():
        indexlistoflists.append(indexls)
    wholemolindexlist=[]
    for index in indexlistoflists[0]:
         wholemolindexlist.append(index)
    fragidxtowholemolidx=dict(zip(fragidxarray, wholemolindexlist))
    fragtypeidxtowholemoltypeidx={}
    for fragidx in fragidxtowholemolidx.keys():
        wholemolidx=fragidxtowholemolidx[fragidx]
        fragtypeidx=fragidxtotypeidx[fragidx]
        wholemoltypeidx=wholeidxtotypeidx[wholemolidx]
        fragtypeidxtowholemoltypeidx[fragtypeidx]=wholemoltypeidx
    newtemp=open('valence.prms','w')
    temp=open(poltype.key6fname,'w')
    for line in fragkeyresults:
        if 'torsion' in line:
            linesplit=line.split()
            typea=int(linesplit[1])
            typeb=int(linesplit[2])
            typec=int(linesplit[3])
            typed=int(linesplit[4])
            tor=[typea,typeb,typec,typed]
            torkey='%d %d %d %d' % (typea, typeb, typec, typed)
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
def GenerateFrag(poltype,molindexlist,basename,makefolder):
    if makefolder==True:
        if not os.path.isdir(basename):
            os.mkdir(basename)
        os.chdir(basename)
    conf = poltype.rdkitmol.GetConformer()
    atomidxstoremove=[]
    for atom in poltype.rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx not in molindexlist:
            atomidxstoremove.append(atomidx)
    fragmol=copy.copy(poltype.rdkitmol)
    rwmol = Chem.RWMol(fragmol)
    sortedatomidxstoremove=sorted(atomidxstoremove, reverse=True)
    for index in sortedatomidxstoremove:
        rwmol.RemoveAtom(index)
    parentindextofragindex={}
    positiontoparentindex={}
    for index in molindexlist:
        pos=conf.GetAtomPosition(index)
        positiontoparentindex[tuple(pos)]=index
    fragconf = rwmol.GetConformer()
    for atom in rwmol.GetAtoms():
        atomindex=atom.GetIdx()
        pos=tuple(fragconf.GetAtomPosition(atomindex)) 
        if pos in positiontoparentindex.keys():
            parentindex=positiontoparentindex[pos]
            parentindextofragindex[parentindex]=atomindex


    filename=basename+'intermediate.mol'
    hydmol=Chem.AddHs(rwmol)
    rdmolfiles.MolToMolFile(hydmol,filename)
    m=rdmolfiles.MolFromMolFile(filename)
    hydmol=Chem.AddHs(m)
    AllChem.EmbedMolecule(hydmol)
    filename=basename+'intermediate.mol'
    rdmolfiles.MolToMolFile(hydmol,filename)
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('mol')
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,filename)
    tmpconv.SetOutFormat('xyz')
    structfname=basename+'intermediate.xyz'
    tmpconv.WriteFile(fragmolbabel,structfname)
 
    return hydmol,parentindextofragindex,structfname

      

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
    rotbndindextofragindexmap={}
    rotbndindextofragment={}
    for tor in torlist: 
        indexes=FirstPassAtomIndexes(poltype,tor)
        fragfoldername=str(tor[1])+'_'+str(tor[2])
        fragmol,parentindextofragindex,structfname=GenerateFrag(poltype,indexes,fragfoldername,True)
        fragWBOmatrix,outputname=GenerateWBOMatrix(poltype,fragmol,structfname)
        fragmentWBOvalue=fragWBOmatrix[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1]
        WBOdifference=numpy.abs(fragmentWBOvalue-parentWBOvalue)
        rotbndidx=tuple([tor[1],tor[2]])
        fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        rotbndindextofragindexmap[fragfoldername]=parentindextofragindex
        sameasparent=False
        if WBOdifference<=poltype.WBOtol: # then we consider the fragment good enough to transfer torsion parameters, so make this fragment into .sdf file
            pass
        else:
            fragmol,newindexes,fragWBOmatrix=GrowFragmentOut(poltype,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername)
            if set(newindexes)==set(indexes):
                sameasparent=True
        WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,[fragrotbndidx],sameasparent)
        os.chdir('..')
        rotbndindextofragment[fragfoldername]=fragmol
    return rotbndindextofragindexmap,rotbndindextofragment

def GrowFragmentOut(poltype,parentWBOmatrix,indexes,WBOdifference,tor,fragfoldername):
    WBOdiffarray=[]
    molarray=[]
    while WBOdifference>poltype.WBOtol:
        fragmentlist=[]
        fragmolidxtoindexlist={}
        fragmentidxtostructfname={}
        possiblefragatmidxs=GrowPossibleFragmentAtomIndexes(poltype,poltype.rdkitmol,indexes)
        for fragmolidx in range(len(possiblefragatmidxs)):
            indexlist=possiblefragatmidxs[fragmolidx]
            basename=fragfoldername+'_'+str(round(WBOdifference,2))
            fragmol,parentindextofragindex,structfname=GenerateFrag(poltype,indexlist,basename,False)
            if set(indexlist)==set(indexes): # then fragment is the parent! copy all files so you dont rerun QM
                return fragmol,indexes,parentWBOmatrix
            fragmentlist.append(fragmol)
            fragmolidxtoindexlist[fragmolidx]=parentindextofragindex.values()
            fragmentidxtostructfname[fragmolidx]=structfname
        WBOdifftoindexlist={}
        WBOdifftofragmol={}
        WBOdifftofragWBOmatrix={}
        for fragmoldx in range(len(fragmentlist)):
           fragmol=fragmentlist[fragmolidx]
           structfname=fragmenttostructfname[fragmolidx]
           fragWBOmatrix,outputname=GenerateWBOMatrix(poltype,fragmol,structfname)
           Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename,structfnamemol)
           fragmentWBOvalue=fragWBOmatrix[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
           parentWBOvalue=parentWBOmatrix[tor[1]-1,tor[2]-1]
           WBOdifference=numpy.abs(fragmentWBOvalue-parentWBOvalue)
           basename=fragfoldername+'_'+str(round(WBOdifference,2))
           structfnamemol=basename+'.mol'
           print(Chem.MolToMolBlock(fragmol),file=open(structfnamemol,'w+')) 
           WBOdifftoindexlist[WBOdifference]=fragmolidxtoindexlist[fragmolidx]
           WBOdifftofragmol[WBOdifference]=fragmol
           WBOdifftofragWBOmatrix[WBOdifference]=fragWBOmatrix
        WBOdifference=min(WBOdifftofragmentlist.keys())
        WBOdiffarray.append(WBOdifference)
        indexes=WBOdifftoindexlist[WBOdifference]
        fragmol=WBOdifftofragmol[WBOdifference]
        molarray.append(fragmol)
        fragWBOmatrix=WBOdifftofragWBOmatrix[WBOdifference]
    PlotFragmenterResults(poltype,WBOdiffarray,molarray)
    return fragmol,indexes,fragWBOmatrix


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
            possiblefragatmidxs.append(indexlist)

        elif (aidx not in indexes and bidx in indexes):
            indexlist=indexes.copy()
            aromaticindexes=GrabRingAtoms(poltype,rdkitmol.GetAtomWithIdx(bidx))
            for atmidx in aromaticindexes:
                if atmidx not in indexlist:
                    indexlist.append(atmidx)
            possiblefragatmidxs.append(indexlist)

    return possiblefragatmidxs
    

def WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,bondindexlist,sameasparent):
    if sameasparent==False:
        structfnamemol=fragfoldername+'.mol'
        print(Chem.MolToMolBlock(fragmol),file=open(structfnamemol,'w+')) 
        tmpconv = openbabel.OBConversion()
        tmpconv.SetInFormat('mol')
        fragmolbabel=openbabel.OBMol()
        tmpconv.ReadFile(fragmolbabel,structfnamemol)
        tmpconv.SetOutFormat('sdf')
        structfname=fragfoldername+'.sdf'
        tmpconv.WriteFile(fragmolbabel,structfname)
        Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,fragfoldername,structfnamemol,bondindexlist)
    else:
        parentdir=dirname(dirname(abspath(os.getcwd())))
        files=os.listdir(parentdir)
        for f in files:
            if os.path.isfile(f):
                shutil.copy(parentdir+r'/'+f,f.replace(poltype.molecprefix,fragfoldername))
    

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
            wbo=str(round(WBOmatrix[begidx,endidx],3))
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
    fig.plot([m.GetNumAtoms() for m in molarray],WBOdiffarray) 
    fig.save(basename+'.png')


    fig=plt.figure()
    basename='NumberofRotatableBondsVSWBODifference'
    fig.plot([rdMolDescriptors.CalcNumRotatableBonds(m) for m in molarray],WBOdiffarray) 
    fig.save(basename+'.png')

