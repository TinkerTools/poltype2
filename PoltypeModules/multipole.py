import electrostaticpotential as esp
import time
import os
import sys
from openbabel import openbabel
from openbabel import pybel
import shutil
import re
from collections import deque
from itertools import product,permutations
import numpy as np
import re
from rdkit import Chem

def AddPolarizeCommentsToKey(poltype,keyfilename,polartypetotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    tempname=keyfilename.replace('.key','_temp.key')
    time.sleep(1)
    temp=open(tempname,'w')
    for line in results:
        if 'polarize' in line and '#' not in line:
            linesplit=line.split()
            typenum=linesplit[1]
            if typenum in polartypetotransferinfo.keys():
                comments=polartypetotransferinfo[typenum]
                line=comments+line
        temp.write(line)
        temp.flush()
        sys.stdout.flush()
    temp.close()
    os.remove(keyfilename)
    os.rename(tempname,keyfilename)
    time.sleep(1)


def SanitizeMultipoleFrames(poltype,keyfilename): # pearl script for averging only understands 0 and not empty spaces for parsing
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    tempname=keyfilename.replace('.key','_temp.key')
    temp=open(tempname,'w')
    extraspace='     '
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            realsplit=re.split(r'(\s+)', line)
            if len(linesplit)==3:
                realsplit=realsplit[:3]+[extraspace]+['0']+[extraspace]+['0']+realsplit[3:]
            elif len(linesplit)==4:
                realsplit=realsplit[:5]+[extraspace]+['0']+realsplit[5:]
            line=''.join(realsplit)
        temp.write(line)
    temp.close()
    os.remove(keyfilename)
    os.rename(tempname,keyfilename)

def CheckIfAllAtomsSameClass(poltype,classlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    allsymm=True
    if len(classlist)>=1:
        firstsymm=poltype.idxtosymclass[classlist[0].GetIdx()]
        for atom in classlist:
            atomidx=atom.GetIdx()
            symmetry=poltype.idxtosymclass[atomidx]
            if symmetry!=firstsymm:
                allsymm=False
    return allsymm

def RemoveFromList(poltype,atomlist,atm):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newatmlist=[]
    idx=atm.GetIdx()
    for newatm in atomlist:
        if idx!=newatm.GetIdx():
            newatmlist.append(newatm)
    return newatmlist


def AtLeastOneHeavyLf1NeighbNotAtom(poltype,lf1atom,atom):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    foundatleastoneheavy=False
    checkneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(lf1atom)]
    for neighb in checkneighbs:
        if neighb.GetIdx()!=atom.GetIdx() and neighb.GetAtomicNum()!=1:
            foundatleastoneheavy=True
    return foundatleastoneheavy


def AtLeastOneHeavyNeighb(poltype,atom):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    foundatleastoneheavy=False
    checkneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(atom)]
    for neighb in checkneighbs:
        if neighb.GetAtomicNum()!=1:
            foundatleastoneheavy=True
    return foundatleastoneheavy


def GrabHeavyAtomIdx(poltype,lf1atom,atom):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    checkneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(lf1atom)]
    for neighb in checkneighbs:
        if neighb.GetIdx()!=atom.GetIdx() and neighb.GetAtomicNum()!=1:
            return neighb.GetIdx()

def FindUniqueNonRepeatingNeighbors(poltype,nlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    d={}
    for b in nlist:
        typenum= poltype.idxtosymclass[b.GetIdx()]
        if typenum not in d.keys():
            d[typenum]=0
        d[typenum]+=1
    uniquetypeneighbsnorepeat=[]
    for typenum in d.keys():
        repeat=d[typenum]
        if repeat==1:
            idx=GrabIndexFromUniqueTypeNumber(poltype,nlist,typenum)
            uniquetypeneighbsnorepeat.append(idx)
    sorteduniquetypeneighbsnorepeat=sorted(uniquetypeneighbsnorepeat,reverse=True)
    return sorteduniquetypeneighbsnorepeat

def GrabIndexFromUniqueTypeNumber(poltype,nlist,typenum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for b in nlist:
        idx=b.GetIdx()
        typenumber= poltype.idxtosymclass[idx]
        if typenumber==typenum:
            return idx

def GrabIndexesFromUniqueTypeNumber(poltype,nlist,typenum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    idxlist=[]
    for b in nlist:
        idx=b.GetIdx()
        typenumber= poltype.idxtosymclass[idx]
        if typenumber==typenum:
            if idx not in idxlist:
                idxlist.append(idx)
    return idxlist




def CheckIfNeighbHasSameType(poltype,a,neighbs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    check=False
    reftype=poltype.idxtosymclass[a.GetIdx()]
    neighbtypes=[poltype.idxtosymclass[b.GetIdx()] for b in neighbs]
    if reftype in neighbtypes:
        check=True
    return check

def GrabNumberOfConnectedHydrogens(poltype,highestsymneighbnorepeat):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(highestsymneighbnorepeat)]
    hydnum=0
    for atom in atomneighbs:
        if atom.GetAtomicNum()==1:
            hydnum+=1
    return hydnum

   
def gen_peditinfile(poltype,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog('Local Frame Symmetry Detection')
    poltype.WriteToLog('Define Polarization Groups and Polarization Parameters')
    lfzerox = [ False ] * mol.NumAtoms()
    if poltype.usepoleditframes==False:
        atomindextoremovedipquad={} # for methane need to make dipole and quadupole on the carbon zeroed out, will return this for post proccesing the keyfile after poledit is run
        atomindextoremovedipquadcross={}
        atomtypetospecialtrace={} # for H on CH4 need to make sure Qxx=Qyy=-1/2*Qzz
        idxtobisecthenzbool={}
        idxtobisectidxs={}
        idxtotrisecbool={}
        idxtotrisectidxs={}
        # Assign local frame for atom a based on the symmetry classes of the atoms it is bound to
        atomiter=openbabel.OBMolAtomIter(mol)
        poltype.localframe2 = list(poltype.localframe2)
        for a in atomiter:
            idxtobisecthenzbool[a.GetIdx()]=False
            idxtotrisecbool[a.GetIdx()]=False
        for atom in openbabel.OBMolAtomIter(mol):
            atomidx=atom.GetIdx()
            hyb=atom.GetHyb()
            atomicnum=atom.GetAtomicNum()
            atomneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(atom)]
            val=len(atomneighbs)
            neighbtypes=list([poltype.idxtosymclass[b.GetIdx()] for b in atomneighbs])
            uniqueneighbtypes=list(set([poltype.idxtosymclass[b.GetIdx()] for b in atomneighbs]))
            sorteduniquetypeneighbsnorepeat=FindUniqueNonRepeatingNeighbors(poltype,atomneighbs)
            foundcase=False
            numhyds=GrabNumberOfConnectedHydrogens(poltype,atom)
            atomisinring=atom.IsInRing()
            if len(sorteduniquetypeneighbsnorepeat)!=0:
                highestsymneighbnorepeatidx=sorteduniquetypeneighbsnorepeat[0]
                highestsymneighbnorepeat=mol.GetAtom(highestsymneighbnorepeatidx)
                highestsymneighbnorepeathyb=highestsymneighbnorepeat.GetHyb()
                highestsymneighbnorepeatatomicnum=highestsymneighbnorepeat.GetAtomicNum()
                numhydsneighb=GrabNumberOfConnectedHydrogens(poltype,highestsymneighbnorepeat)
                neighbsofneighb=[neighb for neighb in openbabel.OBAtomAtomIter(highestsymneighbnorepeat)]
                highestsymneighbnorepeatval=len(neighbsofneighb)
                uniqueneighbtypesofhighestsymneighbnorepeat=list(set([poltype.idxtosymclass[b.GetIdx()] for b in neighbsofneighb]))
                neighbsofneighbwithoutatom=RemoveFromList(poltype,neighbsofneighb,atom)
                neighbswithoutatom=RemoveFromList(poltype,atomneighbs,atom)
                uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom=list(set([poltype.idxtosymclass[b.GetIdx()] for b in neighbsofneighbwithoutatom]))
                neighbindices=list([b.GetIdx() for b in atomneighbs])
                if atomicnum==7 and val==1:
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    foundcase=True

                elif highestsymneighbnorepeatval==3 and CheckIfAllAtomsSameClass(poltype,[neighb for neighb in openbabel.OBAtomAtomIter(highestsymneighbnorepeat)]) and numhydsneighb==3: # then this is like the H on Ammonia and we can use z-then bisector
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    idxtobisecthenzbool[atomidx]=True
                    bisectidxs=[atm.GetIdx() for atm in neighbsofneighbwithoutatom]
                    idxtobisectidxs[atomidx]=bisectidxs
                    foundcase=True
                elif ((len(uniqueneighbtypes)==2 and val==4)) and val==highestsymneighbnorepeatval and len(uniqueneighbtypesofhighestsymneighbnorepeat)==2 and CheckIfNeighbHasSameType(poltype,atom,neighbsofneighbwithoutatom)==False: # then this is like CH3PO3, we want z-onlytrisector would also work
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    foundcase=True
                elif ((val==4 and len(uniqueneighbtypes)==2 and highestsymneighbnorepeatval==3 and highestsymneighbnorepeathyb!=2) or (val==3 and hyb!=2 and len(uniqueneighbtypes)==2 and (highestsymneighbnorepeatval==3 or highestsymneighbnorepeatval==4))) and len(uniqueneighbtypesofhighestsymneighbnorepeat)==2 and len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==1 and (numhyds==2 or numhydsneighb==2):  # then this is like methyl-amine and we can use the two atoms with same symmetry class to do a z-then-bisector, need len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==1 otherwise hits dimethylamine. Need to specify numhyds otherwise will hit nitrobenzene N/C
                    idxtobisecthenzbool[atomidx]=True
                    if atomicnum==6:
                        bisectidxs=[atm.GetIdx() for atm in neighbsofneighbwithoutatom]
                    else:
                        neighbswithoutatom=RemoveFromList(poltype,atomneighbs,highestsymneighbnorepeat)
                        bisectidxs=[atm.GetIdx() for atm in neighbswithoutatom]
                       
                    idxtobisectidxs[atomidx]=bisectidxs
                    poltype.localframe1[atomidx - 1] = highestsymneighbnorepeatidx
                    foundcase=True
                elif (val==3 and len(uniqueneighbtypes)==2 and (highestsymneighbnorepeatval==4) and hyb!=2) and len(uniqueneighbtypesofhighestsymneighbnorepeat)==3 and len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==2: # ethylamine
                        neighbswithoutatom=RemoveFromList(poltype,atomneighbs,highestsymneighbnorepeat)
                        bisectidxs=[atm.GetIdx() for atm in neighbswithoutatom]
                        idxtobisectidxs[atomidx]=bisectidxs
                        poltype.localframe1[atomidx - 1] = highestsymneighbnorepeatidx
                        foundcase=True
                        idxtobisecthenzbool[atomidx]=True
                    # now make sure neighboring atom (lf1) also is using z-then-bisector
                elif val==1 and highestsymneighbnorepeatval==3 and numhydsneighb==1 and len(uniqueneighbtypes)==1 and len(uniqueneighbtypesofhighestsymneighbnorepeat)==2 and highestsymneighbnorepeatatomicnum==7: #dimethylamine
                    idxtobisecthenzbool[atomidx]=True
                    bisectidxs=[atm.GetIdx() for atm in neighbsofneighbwithoutatom]
                    idxtobisectidxs[atomidx]=bisectidxs
                    poltype.localframe1[atomidx - 1] = highestsymneighbnorepeatidx
                    foundcase=True
                elif val==3 and highestsymneighbnorepeatval==1 and numhydsneighb==0 and len(uniqueneighbtypes)==2 and len(uniqueneighbtypesofhighestsymneighbnorepeat)==1 and atomicnum==7: #dimethylamine
                    idxtobisecthenzbool[atomidx]=True

                    neighbswithoutatom=RemoveFromList(poltype,atomneighbs,highestsymneighbnorepeat)
                    bisectidxs=[atm.GetIdx() for atm in neighbswithoutatom]
                    idxtobisectidxs[atomidx]=bisectidxs
                    poltype.localframe1[atomidx - 1] = highestsymneighbnorepeatidx
                    foundcase=True
                elif val==1 and highestsymneighbnorepeatval==3 and len(uniqueneighbtypes)==1 and len(uniqueneighbtypesofhighestsymneighbnorepeat)==2 and len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==1: # H on benzene
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    foundcase=True
                elif val==1 and highestsymneighbnorepeatval==1 and len(uniqueneighbtypes)==1: # H on Iodine
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    foundcase=True
                elif val==4 and len(uniqueneighbtypes)==2 and (len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==0 or len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==1): # N(CH3)(CH3)(CH3)H
                    poltype.localframe1[atomidx-1]=highestsymneighbnorepeat.GetIdx()
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    foundcase=True


            if foundcase==False:
                firstneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(atomneighbs[0])]
                if val==1 and CheckIfAllAtomsSameClass(poltype,[neighb for neighb in openbabel.OBAtomAtomIter(atomneighbs[0])]) and len(firstneighbs)==4: # then this is like H in Methane, we want Z-only
                 
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    zonly=True
                    atomtypetospecialtrace[atomidx]=True
                    atomindextoremovedipquadcross[atomidx]=True
                elif CheckIfAllAtomsSameClass(poltype,atomneighbs) and not AtLeastOneHeavyNeighb(poltype,atom) and val==4: # then this is like carbon in Methane, we want Z-only
                    idxlist=GrabIndexesFromUniqueTypeNumber(poltype,atomneighbs,uniqueneighbtypes[0])
                    poltype.localframe1[atomidx-1]=idxlist[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                    atomindextoremovedipquad[atomidx]=True
                elif CheckIfAllAtomsSameClass(poltype,atomneighbs) and val==3 and numhyds==3: # then this is like Ammonia and we can use a trisector here which behaves like Z-only
                    idxtotrisecbool[atomidx]=True
                    trisectidxs=[atm.GetIdx() for atm in atomneighbs]
                    idxtotrisectidxs[atomidx]=trisectidxs
                    lfzerox[atomidx - 1]=True # need to zero out the x components just like for z-only case
                elif (val==2 and CheckIfAllAtomsSameClass(poltype,atomneighbs)) or (val==4 and len(uniqueneighbtypes)==2) and len(sorteduniquetypeneighbsnorepeat)==0:  # then this is like middle propane carbon or oxygen in water
                    idxlist=GrabIndexesFromUniqueTypeNumber(poltype,atomneighbs,uniqueneighbtypes[0])
                    poltype.localframe1[atomidx-1]=-1*idxlist[0]
                    poltype.localframe2[atomidx-1]=-1*idxlist[1]
                elif len(uniqueneighbtypes)==2 and CheckIfNeighbHasSameType(poltype,atom,atomneighbs) and atomisinring==False: # this handles, ethane,ethene...., z-only
                    poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1] = 0
                    lfzerox[atomidx - 1]=True
                elif len(uniqueneighbtypes)==2 and val==3: # benzene, one carbon in aniline
                    typenumtocnt={}
                    neighbtypes=[poltype.idxtosymclass[b.GetIdx()] for b in atomneighbs]
                    for typenum in neighbtypes:
                        cnt=neighbtypes.count(typenum)
                        typenumtocnt[typenum]=cnt
                    maxcnt=max(typenumtocnt.values())
                    for typenum in uniqueneighbtypes:
                        cnt=typenumtocnt[typenum]
                        if cnt==maxcnt:
                            break 
                    idxlist=GrabIndexesFromUniqueTypeNumber(poltype,atomneighbs,typenum)
                    poltype.localframe1[atomidx-1]=-1*idxlist[0]
                    poltype.localframe2[atomidx-1]=-1*idxlist[1]

                else:
                    if len(sorteduniquetypeneighbsnorepeat)==1:
                        neighboffirstneighbs=[]
                        for n in openbabel.OBAtomAtomIter(highestsymneighbnorepeat):
                            neighboffirstneighbs.append(n)
                        newneighbs=RemoveFromList(poltype,neighboffirstneighbs,atom)
                        newsorteduniquetypeneighbsnorepeat=FindUniqueNonRepeatingNeighbors(poltype,newneighbs)
                        indices=[atom.GetIdx() for atom in neighboffirstneighbs]
                        indices=[atom.GetIdx() for atom in newneighbs]

                        sorteduniquetypeneighbsnorepeat+=newsorteduniquetypeneighbsnorepeat
                    if len(sorteduniquetypeneighbsnorepeat)>=2:
                        aatomidx=atomidx
                        batomidx=sorteduniquetypeneighbsnorepeat[0]
                        catomidx=sorteduniquetypeneighbsnorepeat[1]
                        theindices=[aatomidx,batomidx,catomidx]
                        linear=False
                        combs=permutations(theindices)
                        for comb in combs:
                            atoms=[mol.GetAtom(k) for k in comb]
                            aatom=atoms[0]
                            batom=atoms[1]
                            catom=atoms[2] 
                            angle=mol.GetAngle(aatom,batom,catom)
                            if angle<0:
                                angle=angle+360
                            angletol=3.5
                            if np.abs(angle)<=angletol:
                                linear=True
                        if linear==True:
                            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
                            poltype.localframe2[atomidx - 1] = 0
                            lfzerox[atomidx - 1]=True
                            foundcase=True
                            continue
                    else:
                        neighbs=[]
                        for n in openbabel.OBAtomAtomIter(atom):
                            neighbs.append(n)
                        neighbindices=[k.GetIdx() for k in neighbs]
                        poltype.localframe1[atomidx-1]=neighbindices[0]
                        poltype.localframe2[atomidx - 1] = 0
                        lfzerox[atomidx - 1]=True
                        foundcase=True
                        continue



                    sorteduniquetypeneighbsnorepeattypes=[poltype.idxtosymclass[i] for i in sorteduniquetypeneighbsnorepeat]
                    indicesofsametypetomove=[]
                    atomtype=poltype.idxtosymclass[atomidx]
                    for i in range(len(sorteduniquetypeneighbsnorepeat)):
                        index=sorteduniquetypeneighbsnorepeat[i]
                        typenum=sorteduniquetypeneighbsnorepeattypes[i]
                        if typenum==atomtype:
                            indicesofsametypetomove.append(index)
                    newsorteduniquetypeneighbsnorepeat=[]
                    for index in sorteduniquetypeneighbsnorepeat:  
                        if index not in indicesofsametypetomove:
                            newsorteduniquetypeneighbsnorepeat.append(index)
                    for index in indicesofsametypetomove:
                        newsorteduniquetypeneighbsnorepeat.append(index)
                    symtypetocount={}
                    for typenum in sorteduniquetypeneighbsnorepeattypes:
                        cnt=sorteduniquetypeneighbsnorepeattypes.count(typenum)
                        symtypetocount[typenum]=cnt
                    mincnt=min(symtypetocount.values())
                    if mincnt==1:
                        for symtype,cnt in symtypetocount.items():
                            if cnt==mincnt:
                                break
                        for i in range(len(sorteduniquetypeneighbsnorepeat)):
                            index=sorteduniquetypeneighbsnorepeat[i]
                            typenum=sorteduniquetypeneighbsnorepeattypes[i]
                            if typenum==symtype:
                                break
                        specialindex=index
                        thenewsorteduniquetypeneighbsnorepeat=[specialindex]
                        for i in range(len(sorteduniquetypeneighbsnorepeat)):
                            index=sorteduniquetypeneighbsnorepeat[i]
                            if index!=specialindex:
                                if i==0:
                                    thenewsorteduniquetypeneighbsnorepeat.insert(0,index)
                                else:
                                    thenewsorteduniquetypeneighbsnorepeat.append(index)
                        newsorteduniquetypeneighbsnorepeat=thenewsorteduniquetypeneighbsnorepeat[:]
                    poltype.localframe1[atomidx - 1]=newsorteduniquetypeneighbsnorepeat[0]
                    poltype.localframe2[atomidx - 1]=newsorteduniquetypeneighbsnorepeat[1] 
    # write out the local frames
    iteratom = openbabel.OBMolAtomIter(mol)
    f = open (poltype.peditinfile, 'w')


    rdkitmol = Chem.MolFromMolFile(poltype.molstructfname,removeHs=False)
    if poltype.indextompoleframefile==None:
        if poltype.usepoleditframes==False:
            for a in iteratom:
                if not idxtobisecthenzbool[a.GetIdx()] and not idxtotrisecbool[a.GetIdx()]:
                    f.write(str(a.GetIdx()) + " " + str(poltype.localframe1[a.GetIdx() - 1]) + " " + str(poltype.localframe2[a.GetIdx() - 1]) + "\n")
                elif idxtobisecthenzbool[a.GetIdx()] and not idxtotrisecbool[a.GetIdx()]:
                    bisectidxs=idxtobisectidxs[a.GetIdx()]
                    f.write(str(a.GetIdx()) + " " + str(poltype.localframe1[a.GetIdx() - 1]) + " -" + str(bisectidxs[0])+ " -" + str(bisectidxs[1]) + "\n")
                else:
                    trisecidxs=idxtotrisectidxs[a.GetIdx()]
                    f.write(str(a.GetIdx()) + " -" + str(trisecidxs[0])+ " -" + str(trisecidxs[1]) + " -" + str(trisecidxs[2])+ "\n")
        else:
            # Tinker poledit frame is mostly fine
            # here we add a patch to deal with some molecules
            # Chengwen Liu
            frames = []
            # NH3
            pattern = Chem.MolFromSmarts('[NH3]([H])([H])[H]')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,h1,h2,h3 = match
              frames.append(f'{n+1} -{h1+1} -{h2+1} -{h3+1}')
            # NC3
            pattern = Chem.MolFromSmarts('[NX3]([CX4])([CX4])[CX4]')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,c1,c2,c3 = match
              frames.append(f'{n+1} -{c1+1} -{c2+1} -{c3+1}')
            # NHC2
            pattern = Chem.MolFromSmarts('[NX3H1]([H])([CX4])[CX4]')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,h,c1,c2 = match
              frames.append(f'{n+1} {h+1} -{c1+1} -{c2+1}')
            # nHc2
            pattern = Chem.MolFromSmarts('[#7X3H1;R]([H])([#6;R])[#6;R]')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,h,c1,_ = match
              frames.append(f'{h+1} {n+1} {c1+1}')
            # NH2-C
            pattern = Chem.MolFromSmarts('[NH2]([H])([H])[CX4]')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,h1,h2,c = match
              frames.append(f'{n+1} {c+1} -{h1+1} -{h2+1}')
              frames.append(f'{c+1} {n+1} -{h1+1} -{h2+1}')
            # NH2-C=O
            pattern = Chem.MolFromSmarts('[NH2]([H])([H])[CX3]=O')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,h1,h2,c,_ = match
              # z-then-x
              frames.append(f'{h1+1} {n+1} {c+1}')
              frames.append(f'{h2+1} {n+1} {c+1}')
              # bisect 
              frames.append(f'{n+1} -{h1+1} -{h2+1}')
            # NH2-c 
            # this will give better quality of the ESP
            # comparing to poledit.x frame (Z-bisec)
            pattern = Chem.MolFromSmarts('[NH2]([H])([H])[c]')
            matches = rdkitmol.GetSubstructMatches(pattern)
            for match in matches:
              n,h1,h2,_ = match
              frames.append(f'{h1+1} {n+1} {h2+1}')
              frames.append(f'{h2+1} {n+1} {h1+1}')
              frames.append(f'{n+1} -{h1+1} -{h2+1}')
            
            if frames != []:
              for frame in frames:
                f.write(f'{frame}\n')

    f.write("\n")
    if poltype.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
       f.write('P\n')
    else:
       f.write('A\n')

    # pybel is not robust, it will add hydrogens for some charged/zwitter-ion molecules
    # use rdkit instead!
    # Chengwen Liu
    # July 2023

    atoms = range(1, len(rdkitmol.GetAtoms()) + 1)
    polarDict = dict.fromkeys(atoms, 0)
    lines = open(poltype.updatedsmallmoleculepolarizeprmlib).readlines()
    for line in lines:
      d = line.split()
      smt = d[0]
      pol = d[2]
      pattern = Chem.MolFromSmarts(smt)
      match = rdkitmol.GetSubstructMatches(pattern)
      if match:
        for i in range(len(match)):	
          polarDict[match[i][0]+1] = pol
    for k,v in polarDict.items():  
      f.write(f"{k} {v}\n")
    
    f.write("\n")
    f.flush()
    os.fsync(f.fileno())
    f.write("2\n")
    f.write("N\n")
    f.write("Y\n")


    f.flush()
    os.fsync(f.fileno())

    f.close()
    return lfzerox
    

def post_proc_localframes(poltype,keyfilename, lfzerox):
    """
    Intent: This method runs after the tinker tool Poledit has run and created an
    initial *.key file. The local frames for each multipole are "post processed".
    Zeroed out x-components of the local frame are set back to their original values
    If certain multipole values were not zeroed out by poltype this method zeroes them out
    Input:
       keyfilename: string containing file name of *.key file
       lfzerox: array containing a boolean about whether the x-component of the local frame
                for a given atom should be zeroed or not. Filled in method 'gen_peditin'.
                'lfzerox' is true for atoms that are only bound to one other atom (valence = 1)
                that have more than one possible choice for the x-component of their local frame
    Output: *.key file is edited
    Referenced By: main
    Description:
    1. Move the original *.key file to *.keyb
    2. Create new empty file *.key
    3. Iterate through the lines of *.keyb.
        a. If poledit wrote out the local frame with an x-component missing or as 0
        Then rewrite it with the original x-component (lf2) found in gen_peditin
        b. If poledit did not zero out the local frame x-component for an atom but
        lfzerox is true, zero out the necessary multipole components manually
    4.  Write the edited multipole information to *.key
    """
    # mv *.key to *.keyb
    tmpfname = keyfilename + "b"
    shutil.move(keyfilename, tmpfname)
    keyfh = open(tmpfname)
    lines = keyfh.readlines()
    newlines = []
    # open *.key which is currently empty
    tmpfh = open(keyfilename, "w")
    mpolelines = 0
    # iterate over lines in the *.key file

    for ln1 in range(len(lines)):
        # skip over the lines containing multipole values
        if mpolelines > 0:
            mpolelines -= 1
            continue
        elif 'multipole' in lines[ln1]:
            # Check what poledit wrote as localframe2
            tmplst = lines[ln1].split()
            if len(tmplst) == 5:
                (keywd,atmidx,lf1,lf2,chg) = tmplst
            elif len(tmplst) == 6:
                (keywd,atmidx,lf1,lf2,lf3,chg) = tmplst

            elif len(tmplst) == 4:
                (keywd,atmidx,lf1,chg) = tmplst
                lf2 = '0'
            # manually zero out components of the multipole if they were not done by poledit
            if lfzerox[int(atmidx) - 1]:
                tmpmp = list(map(float, lines[ln1+1].split()))
                tmpmp[0] = 0.
                lines[ln1+1] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0.
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)


            if len(tmplst) == 4:
                linesplit=re.split(r'(\s+)', lines[ln1])
                newtmplist=linesplit[:len(linesplit)-4]
                newtmplist.append('    0')
                newtmplist.append(linesplit[-4])
                newtmplist.append(linesplit[-3])
                newtmplist.append(linesplit[-2])
                templine=''.join(newtmplist)
                newlines.extend(templine)
                newlines.extend(lines[ln1+1:ln1+5])
            else:
                newlines.extend(lines[ln1:ln1+5])
            mpolelines = 4
        # append any lines unrelated to multipoles as is
        else:
            newlines.append(lines[ln1])
    # write out the new lines to *.key
    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()

def rm_esp_terms_keyfile(poltype,keyfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpfname = keyfilename + "_tmp"
    temp=open(tmpfname,'w')
    anothertemp=open(keyfilename,'r')
    results=anothertemp.readlines()
    anothertemp.close()
    passedatomblock=False
    foundatomblock=False
    foundpolarize=False
    for line in results:
        if ('potential-offset' in line or 'none' in line or 'fix-monopole' in line) and 'atom' not in line:
            pass
        else:
            linesplit=line.split()
            if len(linesplit)>0:
                if linesplit[0]=='atom':
                    if foundatomblock==False:
                        foundatomblock=True
                    temp.write(line)

                else:
                    if foundatomblock==True:
                        passedatomblock=True
                    if linesplit[0]=='polarize':
                        foundpolarize=True
                    if passedatomblock==True and foundpolarize==False and linesplit[0]!='#': # then old multipoles
                        pass
                    else:
                        temp.write(line)
            else:
                temp.write(line)

    temp.close()
    shutil.move(tmpfname, keyfilename)
   
    
def prepend_keyfile(poltype,keyfilename,optmol,dipole=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    while not os.path.isfile(keyfilename):
        time.sleep(5)
        poltype.WriteToLog('Waiting for '+keyfilename)
    tmpfname = keyfilename + "_tmp"
    tmpfh = open(tmpfname, "w")
    keyfh = open(keyfilename, "r")
    if len(poltype.paramhead.split()) > 1:
        tmpfh.write("parameters \"" + poltype.paramhead + "\"\n")
    else:
        tmpfh.write("parameters " + poltype.paramhead + "\n")
    tmpfh.write('OPENMP-THREADS '+str(poltype.numproc)+'\n')
    tmpfh.write("bondterm none\n")
    tmpfh.write("angleterm none\n")
    tmpfh.write("torsionterm none\n")
    tmpfh.write("vdwterm none\n")
    tmpfh.write("fix-monopole\n")
    tmpfh.write("digits 8\n")
    tmpfh.write(f"potential-offset {poltype.potentialoffset}\n")
    tmpfh.write("RESP-WEIGHT "+str(poltype.esprestweight)+"\n\n")


    logname=poltype.logespfname
    if poltype.fitqmdipole==True and os.path.exists(logname):
        qmdipole=esp.GrabQMDipoles(poltype,optmol,logname)
        tmpfh.write('TARGET-DIPOLE'+' '+str(qmdipole[0])+' '+str(qmdipole[1])+' '+str(qmdipole[2])+'\n')

    for line in keyfh:
        tmpfh.write(line)
    shutil.move(tmpfname, keyfilename)


def CheckFileForString(poltype,filetocheck):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    temp=open(filetocheck,'r')
    results=temp.readlines()
    temp.close()
    usemp2=False
    for line in results:
        if 'MP2' in line and 'Density' in line:
            usemp2=True
    if usemp2==False:
        string='CC'
    else:
        string='MP2'
    return string


def gen_gdmain(poltype,gdmainfname,molecprefix,fname,dmamethod,new_gdma):
    """
    Intent: Generate GDMA input file for the molecule
    Input:
        gdmainfname: file name for the gdma input file
        molecprefix: name of the molecule
        fname: file name for the *-dma.fchk file
    Output: *.gdmain file is created
    Referenced By: run_gdma
    Description:
        1. create pointer file dma.fchk to *-dma.fchk
        2. create *.gdmain file and write in all the necessary information for the gdma run
    """
    fnamesym = "dma.fchk"
    if not os.path.isfile(fnamesym):
        try:
            os.symlink(fname,fnamesym)
        except Exception as e:
            print(e)

    gdmainp = poltype.gdmainp
    #punfname = os.path.splitext(fname)[0] + ".punch"
    punfname = "dma.punch"

    if poltype.dmamethod=='MP2':
        densitystring=CheckFileForString(poltype,fnamesym)
    else:
        densitystring='SCF'
    gdmainp.update(Title="%s gdmain"%molecprefix, File="%s density %s"%(fnamesym, densitystring), Punch=punfname)
    gdmainp.write_file(gdmainfname)
    if not os.path.isfile(gdmainfname):
        raise ValueError("Unable to write GDMA input file")

def run_gdma(poltype):
    """
    Intent: Runs GDMA to find multipole information
    The GDMA program carries out distributed multipole analysis of the wavefunctions
    calculated by Gaussian (using the *-dma.fchk file)
    Input:
    Output: *.gdmaout is created; this is used as input by poledit, a tinker program
    Referenced By: main
    Description:
    1. Generates the gdma input file by calling 'gen_gdmain'
    2. Runs the following command: gdma < *.gdmain > *.gdmaout
    """
    poltype.WriteToLog("Gaussian Distributed Multipole Analysis (GDMA)")

    if not os.path.isfile(poltype.fckdmafname):
        poltype.fckdmafname = os.path.splitext(poltype.fckdmafname)[0]

    assert os.path.isfile(poltype.fckdmafname), "Error: " + poltype.fckdmafname + " does not exist."+' '+os.getcwd()

    poltype.gdmainfname = poltype.assign_filenames ( "gdmainfname" , ".gdmain")
    gen_gdmain(poltype,poltype.gdmainfname,poltype.molecprefix,poltype.fckdmafname,poltype.dmamethod,poltype.new_gdma)

    cmdstr = poltype.gdmaexe + " < " + poltype.gdmainfname + " > " + poltype.gdmafname
    poltype.call_subsystem([cmdstr],True)

    assert os.path.getsize(poltype.gdmafname) > 0, "Error: " + os.getcwd() +' '+os.path.basename(poltype.gdmaexe) + " cannot create .gdmaout file."
   
def AverageMultipoles(poltype,optmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gen_avgmpole_groups_file(poltype)
    avgmpolecmdstr = '"' + poltype.avgmpolesexe + "\" " + poltype.keyfname + " " + poltype.xyzfname + " " + poltype.grpfname + " " + poltype.key2fnamefromavg + " " + poltype.xyzoutfile + " " + str(poltype.prmstartidx)
    poltype.call_subsystem([avgmpolecmdstr],True)    
    prepend_keyfile(poltype,poltype.key2fnamefromavg,optmol,True)
    statexyzatominfo,oldindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=poltype.GrabXYZInfo(poltype.xyzoutfile)
    babelindextoneighbs={}
    atomiter=openbabel.OBMolAtomIter(optmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        iteratomatom = openbabel.OBAtomAtomIter(atom)
        for natom in iteratomatom:
            natomidx=natom.GetIdx()
            if atomidx not in babelindextoneighbs.keys():
                babelindextoneighbs[atomidx]=[]
            babelindextoneighbs[atomidx].append(natomidx)
    for atomidx in babelindextoneighbs.keys():
        babelneighbs=set(babelindextoneighbs[atomidx])
        tinkneighbs=set(indextoneighbs[atomidx])
        if babelneighbs!=tinkneighbs:
            string='Neighbors for atom index '+str(atomidx)+' not consistent between tinker = '+str(tinkneighbs)+' and babel = '+str(babelneighbs)
            raise ValueError(string) 

def gen_avgmpole_groups_file(poltype):
    """
    Intent: Print out *-groups.txt which is a map from symm class to idx
    Also, symm class labels are altered from 1, 2, 3, ... to 401, 402, ...
    (or some other labeling system dependent on 'prmstartidx')
    Input:
    Output:
        *-groups.txt: map from symm group to atom idx
    Referenced By: main
    Description: -
    """
    symclasstoidxlist={}
    for idx,symclass in poltype.idxtosymclass.items():
        if symclass not in symclasstoidxlist.keys():
            symclasstoidxlist[symclass]=[]
        symclasstoidxlist[symclass].append(idx)
    f = open(poltype.grpfname,"w")
    for symclass,idxlist in symclasstoidxlist.items():
        string=str(symclass)+' '
        for idx in idxlist:
            string+=str(idx)+' '
        string+='\n'
        f.write(string)
    f.close()
