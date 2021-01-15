import math
import sys
import os
from packaging import version
import rdkit
from rdkit.Chem import rdFMCS
import openbabel
from rdkit.Chem import rdmolfiles
import itertools
import re
from rdkit import Chem
import copy
import symmetry as symm
from rdkit.Chem.Lipinski import RotatableBondSmarts
import numpy as np
import json
import torsionfit as torfit

   
def appendtofile(poltype, vf,newname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo):
    temp=open(vf,'r')
    results=temp.readlines()
    temp.close()
    foundatomblock=False
    atomline=False
    wroteout=False
    tempname=vf.replace('.key','_temp.key')
    f=open(tempname,'w')
    for line in results:
        if 'atom' in line:
            atomline=True
            if foundatomblock==False:
                foundatomblock=True
        else:
            atomline=False
        if foundatomblock==True and atomline==False and wroteout==False:
            wroteout=True
            f.write('\n')
            for line,transferinfo in vdwprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
            f.write('\n')
            for line,transferinfo in bondprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
            f.write('\n')
            for line,transferinfo in angleprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
            f.write('\n')
            for line,transferinfo in strbndprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
            f.write('\n')
            for line,transferinfo in opbendprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
            f.write('\n')
            for line,transferinfo in torsionprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
            f.write('\n')
                                    
        else:
            f.write(line)
    f.close()
    os.rename(tempname,newname)
     

def ReadSmallMoleculeLib(poltype,filepath):
    temp=open(filepath,'r')
    results=temp.readlines()
    temp.close()
    smartsatomordertoelementtinkerdescrip={}
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==0:
            continue
        elementsymb=linesplit[0]
        newline=' '.join(linesplit[1:])
        newsplit=newline.split('%')
        tinkerdescrip=newsplit[0].lstrip().rstrip()
        smarts=newsplit[1].lstrip().rstrip()
        atomindices=newsplit[2].lstrip().rstrip()
        atomorderlist=atomindices.split()
        atomorderlist=tuple([int(i) for i in atomorderlist])
        ls=[smarts,atomorderlist]
        newls=[elementsymb,tinkerdescrip]
        smartsatomordertoelementtinkerdescrip[tuple(ls)]=tuple(newls)
    return smartsatomordertoelementtinkerdescrip

def GrabParameters(poltype,fname):
    atomdefs=[]
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    opbendprms=[]
    polarizeprms=[]
    vdwprms=[]
    mpoleprms=[]
    temp=open(fname,'r')
    results=temp.readlines()
    temp.close()
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'atom' in line:
            atomdefs.append(line)
        elif 'bond' in line:
            bondprms.append(line)
        elif 'angle' in line or 'anglep' in line:
            angleprms.append(line)
        elif 'torsion' in line:
            torsionprms.append(line) 
        elif 'strbnd' in line:
            strbndprms.append(line)
        elif 'opbend' in line:
            opbendprms.append(line)
        elif 'polarize' in line:
            polarizeprms.append(line)
        elif 'vdw' in line: 
            vdwprms.append(line)
        elif 'multipole' in line:
            mpolelist=[line,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
            for mpoleline in mpolelist:
                mpoleprms.append(mpoleline)

    return atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms
 
def ShiftPoltypeNumbers(poltype,filename,keyfilename):
    oldtypetonewtype={}
    maxnumberfromprm=GrabMaxTypeNumber(poltype,filename)
    maxnumberfromkey=GrabMaxTypeNumber(poltype,keyfilename)
    minnumberfromkey=GrabMinTypeNumber(poltype,keyfilename)
    typenumbers=list(range(minnumberfromkey,maxnumberfromkey+1))
    newmaxnumber=maxnumberfromprm+1
    shift=minnumberfromkey-newmaxnumber
    oldtypetonewtype={}
    for typenum in typenumbers:
        newtypenum=typenum-shift
        oldtypetonewtype[typenum]=newtypenum

    return oldtypetonewtype,shift

def GrabMaxTypeNumber(poltype,parameterfile):
    maxnumberfromprm=1
    temp=open(parameterfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            atomtype=int(linesplit[1])
            if atomtype>maxnumberfromprm:
                maxnumberfromprm=atomtype
    return maxnumberfromprm

def GrabMinTypeNumber(poltype,parameterfile):
    minnumberfromprm=10000
    temp=open(parameterfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            atomtype=int(linesplit[1])
            if atomtype<minnumberfromprm:
                minnumberfromprm=atomtype
    return minnumberfromprm

def ShiftParameterDefintions(poltype,parameterarray,oldtypetonewtype):
    newparameterarray=[]
    for array in parameterarray:
        newarray=[]
        for line in array:
            try:
                linesplitall=re.split(r'(\s+)', line)
            except:
                print('Error ',line)
            for i in range(len(linesplitall)):
                element=linesplitall[i]
                if RepresentsInt(poltype,element):
                    oldtypenum=np.abs(int(element))
                    if oldtypenum in oldtypetonewtype.keys():
                        newtypenum=oldtypetonewtype[oldtypenum]
                        typenum=newtypenum
                    else:
                        typenum=oldtypenum
                    if '-' in element:
                        typenum=-typenum
                    linesplitall[i]=str(typenum)
            newline=''.join(linesplitall)
            newarray.append(newline)
        newparameterarray.append(newarray)
    return newparameterarray

def WriteToPrmFile(poltype,atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms,prmpath):
    # assumes that there is white space at end of every parameter block, just add a line to key_5
    temp=open(prmpath,'r')
    results=temp.readlines()
    temp.close()
    lastidx=len(results)-1
    newprmname=prmpath.replace('.prm','_temp.prm')
    temp=open(newprmname,'w')
    linelist=[]
    firstpassatomdef=False
    firstpassbond=False
    firstpassangle=False
    firstpasstorsion=False
    firstpassstrbnd=False
    firstpasspolarize=False
    firstpassvdw=False
    firstpassmpole=False
    firstpassopbend=False
    passedatomdefblock=False
    passedbondblock=False
    passedangleblock=False
    passedtorsionblock=False
    passedstrbndblock=False
    passedpolarizeblock=False
    passedvdwblock=False
    passedmpoleblock=False
    passedopbendblock=False

    
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'type' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line and 'unit' not in line and 'scale' not in line:
            prewhitespacelinesplit=re.split(r'(\s+)', line)
            if 'atom' in line and firstpassatomdef==False:
                firstpassatomdef=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'atom' not in line and firstpassatomdef==True and passedatomdefblock==False:
                passedatomdefblock=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                for atomdef in atomdefs:
                    temp.write(atomdef)
                temp.write(line)
            elif 'vdw' in line and firstpassvdw==False:
                firstpassvdw=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'vdw' not in line and firstpassvdw==True and passedvdwblock==False:
                passedvdwblock=True
                for vdwline in vdwprms:
                    temp.write(vdwline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'bond' in line and firstpassbond==False:
                firstpassbond=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'bond' not in line and firstpassbond==True and passedbondblock==False:
                passedbondblock=True
                for bondline in bondprms:
                    temp.write(bondline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'angle' in line and firstpassangle==False:
                firstpassangle=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'angle' not in line and firstpassangle==True and passedangleblock==False:
                passedangleblock=True
                for angleline in angleprms:
                    temp.write(angleline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'torsion' in line and firstpasstorsion==False:
                firstpasstorsion=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'torsion' not in line and firstpasstorsion==True and passedtorsionblock==False:
                passedtorsionblock=True
                for torsionline in torsionprms:
                    temp.write(torsionline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'strbnd' in line and firstpassstrbnd==False:
                firstpassstrbnd=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'strbnd' not in line and firstpassstrbnd==True and passedstrbndblock==False:
                passedstrbndblock=True
                for strbndline in strbndprms:
                    temp.write(strbndline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'polarize' in line and firstpasspolarize==False:
                firstpasspolarize=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'polarize' not in line and firstpasspolarize==True and passedpolarizeblock==False:
                passedpolarizeblock=True
                for polarizeline in polarizeprms:
                    temp.write(polarizeline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'opbend' in line and firstpassopbend==False:
                firstpassopbend=True
                temp.write(line)
            elif 'opbend' not in line and firstpassopbend==True and passedopbendblock==False:
                passedopbendblock=True
                for opbendline in opbendprms:
                    temp.write(opbendline)
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'multipole' in line and firstpassmpole==False:
                firstpassmpole=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'multipole' not in line and firstpassmpole==True and passedmpoleblock==False and line=='\n':
                passedmpoleblock=True
                for mpoleline in mpoleprms:
                    temp.write(mpoleline)
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()

    temp.flush()
    os.fsync(temp.fileno())
    sys.stdout.flush()
    temp.close()
    os.remove(prmpath)
    os.rename(newprmname,prmpath)

def GrabParametersFromPrmFile(poltype,bondtinkerclassestopoltypeclasses,opbendtinkerclassestopoltypeclasses,opbendtinkerclassestotrigonalcenterbools,angletinkerclassestopoltypeclasses,torsiontinkerclassestopoltypeclasses,poltypetoprmtype,atomtinkerclasstopoltypeclass,typestoframedefforprmfile,fname,skipmultipole=False):
    temp=open(fname,'r') 
    results=temp.readlines()
    temp.close()
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    mpoleprms=[]
    opbendprms=[]
    polarizeprms=[]
    vdwprms=[]
    torsiontopitor={}
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        linesplitall=re.split(r'(\s+)', line)
        if '#' in line:
            continue
        if 'bond' in line and 'cubic' not in line and 'quartic' not in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundbond=False
            if tuple(bondclasslist) in bondtinkerclassestopoltypeclasses.keys():
                bondtup=tuple(bondclasslist)
                foundbond=True
            elif tuple(bondclasslist[::-1]) in bondtinkerclassestopoltypeclasses.keys():
                bondtup=tuple(bondclasslist[::-1])
                foundbond=True
            if foundbond==True:
                classes=bondtinkerclassestopoltypeclasses[bondtup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    newline=''.join(linesplitall)
                    bondprms.append(newline)           
        elif 'angle-cubic' not in line and 'angle-quartic' not in line and 'pentic' not in line and 'sextic' not in line and ('angle' in line or 'anglep' in line) :
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            foundangle=False
            if tuple(angleclasslist) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist)
                foundangle=True
            elif tuple(angleclasslist[::-1]) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundangle=True
            if foundangle==True:
                classes=angletinkerclassestopoltypeclasses[angletup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    newline=''.join(linesplitall)
                    angleprms.append(newline) 
        
        elif 'torsion' in line and 'torsionunit' not in line:
            torsionclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
            foundtorsion=False
            if tuple(torsionclasslist) in torsiontinkerclassestopoltypeclasses.keys():
                torsiontup=tuple(torsionclasslist)
                foundtorsion=True
            elif tuple(torsionclasslist[::-1]) in torsiontinkerclassestopoltypeclasses.keys():
                torsiontup=tuple(torsionclasslist[::-1])
                foundtorsion=True
            if foundtorsion==True:   
                classes=torsiontinkerclassestopoltypeclasses[torsiontup]
                for boundcls in classes:
                    if linesplitall[0]=='': # sometimes when torsion is added back to the key file, it has a space in front of it
                        linesplitall[4]=str(boundcls[0])    
                        linesplitall[6]=str(boundcls[1])  
                        linesplitall[8]=str(boundcls[2])
                        linesplitall[10]=str(boundcls[3])
                    else:
                        linesplitall[2]=str(boundcls[0])    
                        linesplitall[4]=str(boundcls[1])  
                        linesplitall[6]=str(boundcls[2])
                        linesplitall[8]=str(boundcls[3])

                    newline=''.join(linesplitall)
                    torsionprms.append(newline)

        elif 'strbnd' in line:
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            foundstrbnd=False
            if tuple(angleclasslist) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist)
                foundstrbnd=True
            elif tuple(angleclasslist[::-1]) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundstrbnd=True
            if foundstrbnd==True:
                classes=angletinkerclassestopoltypeclasses[angletup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    newline=''.join(linesplitall)
                    strbndprms.append(newline) 

               
        
        elif 'opbend' in line and 'opbendtype' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundopbend=False

            if tuple(bondclasslist) in opbendtinkerclassestopoltypeclasses.keys():
                bondtup=tuple(bondclasslist)
                foundopbend=True

            if foundopbend==True:
                classes=opbendtinkerclassestopoltypeclasses[bondtup]
                boolarray=opbendtinkerclassestotrigonalcenterbools[bondtup]

                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    newlinesplitall=linesplitall[:]
                    newlinesplitall[4]=str(boundcls[0])    
                    newlinesplitall[2]=str(boundcls[1])
                    newline=''.join(linesplitall)
                    if boolarray[1]==True: 
                        opbendprms.append(newline)

        elif 'multipole' in line and skipmultipole==False:
            newlinesplit=linesplit[1:-1]
            frames=[int(i) for i in newlinesplit]
            grabit=False
            if tuple(frames) in typestoframedefforprmfile.keys():
                theframe=tuple(frames)
                grabit=True
            elif tuple(frames[::-1]) in typestoframedefforprmfile.keys():
                theframe=tuple(frames[::-1])
                grabit=True
            if grabit==True:
                chgpartofline=linesplitall[-2:] # include space
                phrasepartofline=linesplitall[:2] # include space
                spacetoadd='   '
                newlist=[]
                newlist.extend(phrasepartofline)
                framedef=typestoframedefforprmfile[atomtype]
                for typenum in framedef:
                    newlist.append(str(typenum))
                    newlist.append(spacetoadd)
                newlist=newlist[:-1] # remove last space then add space before chg
                newlist.extend(chgpartofline)
                newline=''.join(newlist)
                mpolelist=[newline,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
                for mpoleline in mpolelist:
                    mpoleprms.append(mpoleline)

        elif 'polarize' in line: 
            atomtype=int(linesplit[1])
            if atomtype in poltypetoprmtype.keys():
                prmtype=poltypetoprmtype[atomtype]
                newline=line.replace('\n','')+' '+str(prmtype)+'\n'
                polarizeprms.append(newline)
         
        elif 'vdw' in line and 'type' not in line and 'scale' not in line: 
            atomclass=int(linesplit[1])
            atomclasslist=tuple([atomclass])
            if atomclasslist in atomtinkerclasstopoltypeclass.keys():
                prmclasses=atomtinkerclasstopoltypeclass[atomclasslist]
                for prmclasslist in prmclasses:
                    for prmclass in prmclasslist:
                        linesplit[1]=str(prmclass)
                        newline=' '.join(linesplit)+'\n'
                        vdwprms.append(newline)

        elif 'pitor' in line:
            linesplit=line.split()
            b=linesplit[1]
            c=linesplit[2]
            for tinkerclasses,poltypeclasses in torsiontinkerclassestopoltypeclasses.items():
                tb=str(tinkerclasses[1])
                tc=str(tinkerclasses[2])
                if (b==tb and c==tc) or (b==tc and c==tb):
                    for cls in poltypeclasses:
                        torsiontopitor[tuple(cls)]=line
                        torsiontopitor[tuple(cls[::-1])]=line

    
    return bondprms,angleprms,torsionprms,strbndprms,mpoleprms,opbendprms,polarizeprms,vdwprms,torsiontopitor


def GrabTypeAndClassNumbers(poltype,prmfile):
    temp=open(prmfile,'r')
    results=temp.readlines()
    temp.close()
    elementtinkerdescriptotinkertype={}
    tinkertypetoclass={}
    for line in results:
        if 'atom' in line:
            linesplit=line.split()    
            newlinesplit=linesplit[1:-3]
            tinkertype=newlinesplit[0]
            classtype=newlinesplit[1]
            element=newlinesplit[2]
            tinkerdescrip=' '.join(newlinesplit[3:])
            ls=[element,tinkerdescrip]
            elementtinkerdescriptotinkertype[tuple(ls)]=tinkertype
            tinkertypetoclass[tinkertype]=classtype
    return elementtinkerdescriptotinkertype,tinkertypetoclass



def GrabAtomsForParameters(poltype,mol):
    # now we can define arrays to collect bonds, angles and torsions
    listoftorsionsforprm=[]
    listofbondsforprm=[]
    listofanglesforprm=[]
    listofatomsforprm=[]
    for atom in openbabel.OBMolAtomIter(mol):
        atomidx=atom.GetIdx()-1
        listofatomsforprm.append([atomidx])
        neighbs=[natom for natom in openbabel.OBAtomAtomIter(atom)]
        for natom in neighbs:
            nidx=natom.GetIdx()-1
            bondset=[nidx,atomidx]
            if bondset not in listofbondsforprm and bondset[::-1] not in listofbondsforprm:
                listofbondsforprm.append(bondset)
            nextneighbs=[nextatom for nextatom in openbabel.OBAtomAtomIter(natom)]
            for nextneighb in nextneighbs:
                nextneighbidx=nextneighb.GetIdx()-1
                if nextneighbidx!=atomidx:
                    angleset=[nextneighbidx,nidx,atomidx]
                    if angleset not in listofanglesforprm and angleset[::-1] not in listofanglesforprm:
                        listofanglesforprm.append(angleset)
                    nextnextneighbs=[nextnextatom for nextnextatom in openbabel.OBAtomAtomIter(nextneighb)]
                    for nextnextatom in nextnextneighbs:
                        nextnextatomidx=nextnextatom.GetIdx()-1
                        if nextnextatomidx!=nidx:                
                            torsionset=[nextnextatomidx,nextneighbidx,nidx,atomidx]
                            if torsionset not in listoftorsionsforprm and torsionset[::-1] not in listoftorsionsforprm:
                                listoftorsionsforprm.append(torsionset)
    return listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm


def GenerateSMARTSListFromAtomList(poltype,listforprm,rdkitmol,mol,maxatomsize):
    listforprmtosmartslist={}
    atomnum=rdkitmol.GetNumAtoms()
    for ls in listforprm:
        smartslist=[]
        atomindiceslist=GenerateAllPossibleFragmentIndices(poltype,ls,rdkitmol,maxatomsize)
        anyringatoms,ringindices=LookForRingAtoms(poltype,mol,ls)
        lengthstoatomindices={}
        for atomindices in atomindiceslist:
            lengthstoatomindices[len(atomindices)]=atomindices
        atomindiceslist=list(lengthstoatomindices.values()) 

        if anyringatoms==True:
            ringsize=FindRingSize(poltype,mol,ringindices)
            if ringsize!=0:
                singlepathindices=GenerateAllPossibleSinglePathFragmentIndices(poltype,ls,rdkitmol,ringsize)
                singlepathindices.sort(key=len)
                lengthtosinglepathindices={}
                for indices in singlepathindices:
                    lengthtosinglepathindices[len(indices)]=indices
                singlepathindices=list(lengthtosinglepathindices.values()) 

                atomindiceslist.extend(singlepathindices)
        for subls in atomindiceslist:
            fragsmarts,smartsfortransfer=GenerateFragmentSMARTS(poltype,rdkitmol,mol,subls)
            fragsmarts=FilterSMARTS(poltype,fragsmarts)
            smartsfortransfer=FilterSMARTS(poltype,smartsfortransfer)

            smartslist.append([fragsmarts,smartsfortransfer])
        listforprmtosmartslist[tuple(ls)]=smartslist
    return listforprmtosmartslist

def FilterSMARTS(poltype,smarts):
    smarts=smarts.replace('-]',']')
    return smarts

def FindRingSize(poltype,mol,ringindices):
    ringsizes=[]
    for index in ringindices:
        atom=mol.GetAtom(index)
        for i in range(3,9+1):
            isinringofsize=atom.IsInRingSize(i)
            if isinringofsize==True:
                ringsizes.append(i)
    if len(ringsizes)!=0:
        maxringsize=max(ringsizes)
    else:
        maxringsize=0

    return maxringsize


def LookForRingAtoms(poltype,mol,ls):
    anyringatoms=False
    babelindices=[i+1 for i in ls]
    ringindices=[]
    for index in babelindices:
        atom=mol.GetAtom(index)
        isinring=atom.IsInRing()
        if isinring==True:
            anyringatoms=True
            ringindices.append(index)
    return anyringatoms,ringindices



def GenerateAllPossibleSinglePathFragmentIndices(poltype,ls,rdkitmol,ringsize): # for rings
    atomindiceslist=[]
    for n in range(1,ringsize+1):
        paths=Chem.rdmolops.FindAllPathsOfLengthN(rdkitmol,n,useBonds=False,useHs=True,rootedAtAtom=ls[0])
        for path in paths:
            allin=True
            for index in ls:
                if index not in path:
                    allin=False
            if allin==True:
                path=list(dict.fromkeys(path)) # sometimes atoms repeated, need to throw out for SMARTS generation
                atomindiceslist.append(list(path))
    return atomindiceslist


def GenerateAllPossibleFragmentIndices(poltype,ls,rdkitmol,maxatomsize):
    atomindiceslist=[copy.deepcopy(ls)]
    oldindexlist=copy.deepcopy(ls)
    indexlist=[]
    count=0
    neighbcount=0
    while set(oldindexlist)!=set(indexlist):
        if count!=0:
            oldindexlist=copy.deepcopy(indexlist)
        if len(oldindexlist)>maxatomsize:
            break
        if neighbcount>=2:
            break
        neighborindexes=GrabNeighboringIndexes(poltype,oldindexlist,rdkitmol)
        neighbcount+=1
       
        for i in range(len(neighborindexes)):
            combs=list(itertools.combinations(neighborindexes, i+1))
            for comb in combs:
                newcomb=copy.deepcopy(oldindexlist)
                newcomb.extend(comb) 
                if newcomb not in atomindiceslist:
                    atomindiceslist.append(newcomb)
        count+=1
        newindexlist=copy.deepcopy(oldindexlist)
        for index in neighborindexes:
            if index not in newindexlist:
                newindexlist.append(index)
        indexlist=newindexlist
  
    return atomindiceslist


def GrabNeighboringIndexes(poltype,indexlist,rdkitmol):
    neighborindexes=[]
    for index in indexlist:
        atom=rdkitmol.GetAtomWithIdx(index)
        for natom in atom.GetNeighbors():
            natomidx=natom.GetIdx()
            if natomidx not in neighborindexes and natomidx not in indexlist:
                neighborindexes.append(natomidx)
    return neighborindexes
            
             
def GenerateFragmentSMARTS(poltype,rdkitmol,mol,ls):
    rdkitsmarts=rdmolfiles.MolToSmarts(rdkitmol)
    newmol = Chem.Mol()
    mw = Chem.RWMol(newmol)
    # need to treat ring bonds as aromatic since all transferred parameters from amoeba09 are aromatic rings
    oldindextonewindex={}
    aromaticindices=[]
    for i,idx in enumerate(ls):
        oldatom=rdkitmol.GetAtomWithIdx(idx)
        mw.AddAtom(oldatom)
        oldindextonewindex[idx]=i
        oldatombabel=mol.GetAtom(idx+1)
        isinring=oldatombabel.IsInRing()
        if isinring==True:
            aromaticindices.append(i)

    

    atomswithcutbonds=[]
    aromaticbonds=[]
    bonditer=rdkitmol.GetBonds()
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        babelbond=mol.GetBond(oendidx+1,obgnidx+1)
        isinring=babelbond.IsInRing()
        if oendidx in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            if oldindextonewindex[oendidx] not in atomswithcutbonds:
                atomswithcutbonds.append(oldindextonewindex[oendidx])
            continue
        if oendidx not in oldindextonewindex.keys() and obgnidx in oldindextonewindex.keys():
            if oldindextonewindex[obgnidx] not in atomswithcutbonds:
                atomswithcutbonds.append(oldindextonewindex[obgnidx])
            continue
        if oendidx not in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            continue
        endidx=oldindextonewindex[oendidx]
        bgnidx=oldindextonewindex[obgnidx]
        if isinring==True:
            aromaticbonds.append([endidx,bgnidx]) 
        bondorder=bond.GetBondType()
        mw.AddBond(bgnidx,endidx,bondorder)



    smarts=rdmolfiles.MolToSmarts(mw)

    for atom in mw.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx in aromaticindices:
            atom.SetIsAromatic(True)


    for bond in mw.GetBonds():
        endidx = bond.GetEndAtomIdx()
        bgnidx = bond.GetBeginAtomIdx()
        temp=[endidx,bgnidx]
        if temp in aromaticbonds or temp[::-1] in aromaticbonds:
            bond.SetIsAromatic(True)

    smartsfortransfer=rdmolfiles.MolToSmarts(mw)

    return smarts,smartsfortransfer


def MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listforprmtosmartslist,parametersmartslist,mol):
    listforprmtoparametersmarts={}
    listforprmtosmarts={}
    for ls,smartslist in listforprmtosmartslist.items(): 
        parametersmartstomatchlen={}
        parametersmartstosmartslist={}
        parametersmartstomatchlen,parametersmartstosmartslist=MatchAllPossibleSMARTSToParameterSMARTS(poltype,smartslist,parametersmartslist,parametersmartstomatchlen,parametersmartstosmartslist)
        if len(parametersmartstomatchlen.keys())==0:
            smartslist=ReplaceSMARTSBondsWithGenericBonds(poltype,smartslist,ls,mol)
            parametersmartstomatchlen,parametersmartstosmartslist=MatchAllPossibleSMARTSToParameterSMARTS(poltype,smartslist,parametersmartslist,parametersmartstomatchlen,parametersmartstosmartslist)
            if len(parametersmartstomatchlen.keys())==0:
                smartslist=ReplaceEachAtomIdentityOnEnd(poltype,smartslist[0]) 
                parametersmartstomatchlen,parametersmartstosmartslist=MatchAllPossibleSMARTSToParameterSMARTS(poltype,smartslist,parametersmartslist,parametersmartstomatchlen,parametersmartstosmartslist)
                if len(parametersmartstomatchlen.keys())==0:
                    smartslist=ReplaceAtomIdentitiesOnEnd(poltype,smartslist[0])
                    parametersmartstomatchlen,parametersmartstosmartslist=MatchAllPossibleSMARTSToParameterSMARTS(poltype,smartslist,parametersmartslist,parametersmartstomatchlen,parametersmartstosmartslist)
                    if len(parametersmartstomatchlen.keys())==0:
                        smartslist=ReplaceAllAtomIdentities(poltype,smartslist[0])
                        parametersmartstomatchlen,parametersmartstosmartslist=MatchAllPossibleSMARTSToParameterSMARTS(poltype,smartslist,parametersmartslist,parametersmartstomatchlen,parametersmartstosmartslist)
                        listforprmtoparametersmarts,listforprmtosmarts=GrabBestMatch(poltype,parametersmartstomatchlen,parametersmartstosmartslist,listforprmtoparametersmarts,listforprmtosmarts,ls)
                        
                    else:
                        listforprmtoparametersmarts,listforprmtosmarts=GrabBestMatch(poltype,parametersmartstomatchlen,parametersmartstosmartslist,listforprmtoparametersmarts,listforprmtosmarts,ls)

                else:
                    listforprmtoparametersmarts,listforprmtosmarts=GrabBestMatch(poltype,parametersmartstomatchlen,parametersmartstosmartslist,listforprmtoparametersmarts,listforprmtosmarts,ls)
                    
            else:
                listforprmtoparametersmarts,listforprmtosmarts=GrabBestMatch(poltype,parametersmartstomatchlen,parametersmartstosmartslist,listforprmtoparametersmarts,listforprmtosmarts,ls)
        else:
            listforprmtoparametersmarts,listforprmtosmarts=GrabBestMatch(poltype,parametersmartstomatchlen,parametersmartstosmartslist,listforprmtoparametersmarts,listforprmtosmarts,ls)
    return listforprmtoparametersmarts,listforprmtosmarts

def ReplaceEachAtomIdentityOnEnd(poltype,smartsls):
    newsmartslist=[]
    templeft=[]
    tempright=[]
    for smarts in smartsls:
        templeft.append(ReplaceLeftAtomIdentity(poltype,smarts))
        tempright.append(ReplaceRightAtomIdentity(poltype,smarts))
    newsmartslist=[templeft,tempright]
    return newsmartslist

def ReplaceLeftAtomIdentity(poltype,smarts):
    smartsplit=smarts.split('~') 
    tempsplit=smartsplit[:]
    addp=False
    if tempsplit[0][-1]=='(':
        addp=True
    tempsplit[0]='[*]'
    if addp==True:
       tempsplit[0]+='('
    newsmarts='~'.join(tempsplit)
    return newsmarts

def ReplaceRightAtomIdentity(poltype,smarts):
    smartsplit=smarts.split('~') 
    tempsplit=smartsplit[:]
    tempsplit[-1]='[*]'
    newsmarts='~'.join(tempsplit)
    return newsmarts

def ReplaceAtomIdentitiesOnEnd(poltype,smartsls):
    newsmartslist=[]
    for smarts in smartsls:
        smartsplit=smarts.split('~') 
        tempsplit=copy.deepcopy(smartsplit)
        tempsplit[0]='[*]'
        tempsplit[-1]='[*]'
        newsmarts='~'.join(tempsplit)
        newsmartslist.append(newsmarts)
    newsmartslist=[newsmartslist]
    return newsmartslist


def ReplaceAllAtomIdentities(poltype,smartsls):
    newsmartslist=[]
    for smarts in smartsls:
        smartsplit=smarts.split('~') 
        tempsplit=copy.deepcopy(smartsplit)
        for idx in range(len(tempsplit)):
            tempsplit[idx]='[*]'
        newsmarts='~'.join(tempsplit)
        newsmartslist.append(newsmarts) 
    newsmartslist=[newsmartslist]
    return newsmartslist


def GenerateElementCountsDictionary(poltype,structure):
    countdic={}
    for atom in structure.GetAtoms():
        atomicnum=atom.GetAtomicNum()
        if atomicnum not in countdic.keys():
            countdic[atomicnum]=0
        countdic[atomicnum]+=1
    return countdic 



def GrabBestMatch(poltype,parametersmartstomatchlen,parametersmartstosmartslist,listforprmtoparametersmarts,listforprmtosmarts,ls):
    maxparametersmartsmatchlength=max(parametersmartstomatchlen.values())
    parametersmartslist=GrabKeysFromValue(poltype,parametersmartstomatchlen,maxparametersmartsmatchlength)
    prmsmartstodiff={}
    for prmsmarts in parametersmartslist:
        smartsls=parametersmartstosmartslist[prmsmarts]
        smartsfortransfer=smartsls[1]
        substructure = Chem.MolFromSmarts(smartsfortransfer)
        substructurenumatoms=substructure.GetNumAtoms()
        structure = Chem.MolFromSmarts(prmsmarts)
        mols = [structure,substructure]
        res=rdFMCS.FindMCS(mols)
        structurenumatoms=structure.GetNumAtoms()
        atomnum=res.numAtoms
        diff=np.abs(atomnum-structurenumatoms)
        prmsmartstodiff[prmsmarts]=diff
    mindiff=min(prmsmartstodiff.values())
    parametersmartslist=GrabKeysFromValue(poltype,prmsmartstodiff,mindiff)
    smartstocountdic={}
    atomicnums=[]
    for prmsmarts in parametersmartslist:
        countdic=GenerateElementCountsDictionary(poltype,Chem.MolFromSmarts(prmsmarts))
        for atomicnum in countdic.keys():
            if atomicnum not in atomicnums:
                atomicnums.append(atomicnum)
        smartstocountdic[prmsmarts]=countdic
    for smarts,countdic in smartstocountdic.items():
        for atomicnum in atomicnums:
            if atomicnum not in countdic.keys():
                countdic[atomicnum]=0
        smartstocountdic[smarts]=countdic

    mastercountdic=GenerateElementCountsDictionary(poltype,poltype.rdkitmol)
    for atomicnum in atomicnums:
        if atomicnum not in mastercountdic.keys():
            mastercountdic[atomicnum]=0
    difftosmarts={}
    for smarts,countdic in smartstocountdic.items():
        totdiff=0
        for atomicnum,counts in countdic.items():
            mastercounts=mastercountdic[atomicnum]
            diff=np.abs(mastercounts-counts)
            totdiff+=diff 
        difftosmarts[totdiff]=smarts
    
    
    mindiff=min(difftosmarts.keys())
    parametersmarts=difftosmarts[mindiff]
    smartsls=parametersmartstosmartslist[parametersmarts]
    listforprmtoparametersmarts[ls]=parametersmarts 
    listforprmtosmarts[ls]=smartsls
    return listforprmtoparametersmarts,listforprmtosmarts


def RemoveHydrogens(poltype,mol):
    em = Chem.EditableMol(mol)
    indexestoremove=[]
    for atom in mol.GetAtoms():
        atomicnum=atom.GetAtomicNum()
        atomidx=atom.GetIdx()
        if atomicnum==1:
            indexestoremove.append(atomidx)
    indexestoremove.sort(reverse=True)
    for idx in indexestoremove:
        em.RemoveAtom(idx)
    newmol = em.GetMol()
    return newmol

def ReplaceSMARTSBondsWithGenericBonds(poltype,smartslist,atomindices,mol):
    newsmartslist=[]
    for smartsls in smartslist:
        temp=[]
        for smarts in smartsls:
            newsmarts=smarts.replace('-','~').replace('=','~').replace(':','~')
            temp.append(newsmarts)
        newsmartslist.append(temp)
    return newsmartslist
         

def MatchAllPossibleSMARTSToParameterSMARTS(poltype,smartslist,parametersmartslist,parametersmartstomatchlen,parametersmartstosmartslist):
    for smartls in smartslist:
        smartsfortransfer=smartls[1]
        smarts=smartls[0]
        for parametersmarts in parametersmartslist:
            substructure = Chem.MolFromSmarts(smartsfortransfer)
            substructurenumatoms=substructure.GetNumAtoms()
            structure = Chem.MolFromSmarts(parametersmarts)
            structurenumatoms=structure.GetNumAtoms()
            if structurenumatoms>=substructurenumatoms:
                diditmatch=structure.HasSubstructMatch(substructure)
                if diditmatch==True:
                    matches=structure.GetSubstructMatches(substructure)
                    firstmatch=matches[0]
                    parametersmartstomatchlen[parametersmarts]=len(firstmatch)
                    parametersmartstosmartslist[parametersmarts]=smartls
    return parametersmartstomatchlen,parametersmartstosmartslist


def CheckIfOrientationIsSameAsPrevious(poltype,currentorientation,smartsindextoparametersmartsindex,moleculeindextosmartsindices,atomindices):
    orientation=[]
    keep=True
    allpossiblesmartsindices=[moleculeindextosmartsindices[i] for i in atomindices]
    allpossibleprmsmartsindices=[]
    for smartsindices in allpossiblesmartsindices:
        prmindices=[smartsindextoparametersmartsindex[i] for i in smartsindices]
        allpossibleprmsmartsindices.append(prmindices)

    if len(allpossibleprmsmartsindices)==3:
        orientation=[set(allpossibleprmsmartsindices[1])]
    elif len(allpossibleprmsmartsindices)==4:
        orientation=[set(allpossibleprmsmartsindices[1]),set(allpossibleprmsmartsindices[2])]
    if orientation not in currentorientation and len(currentorientation)!=0:
        keep=False

    return orientation,keep

def GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol):
    atomindicestotinkertypes={}
    atomindicestotinkerclasses={}
    atomindicestoparametersmartsatomorders={}
    atomindicestoelementtinkerdescrips={}
    atomindicestosmartsatomorders={}
    for atomindices,parametersmarts in atomindicesforprmtoparametersmarts.items():
        smartsls=atomindicesforprmtosmarts[atomindices]
        smarts=smartsls[0]
        smartsfortransfer=smartsls[1]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)
        moleculeindextosmartsindices={}
        smartsindextomoleculeindices={} 
        currentorientation=[]
        orient=[]
        for match in matches:
            allin=True
            for idx in atomindices:
                if idx not in match:
                    allin=False
            if allin==True:
                indices=list(range(len(match)))
                smartsindextomoleculeindex=dict(zip(indices,match)) 
                moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}
                smartsindices=[moleculeindextosmartsindex[i] for i in atomindices]
                revsmartsindices=smartsindices[::-1]
                if revsmartsindices in currentorientation:
                    continue
                if len(smartsindices)==3:
                    orientation=[smartsindices[1]]
                elif len(smartsindices)==4:
                    orientation=[smartsindices[1],smartsindices[2]]
                else:
                    orientation=[]
                if len(orientation)!=0 and len(orient)!=0 and orientation not in orient:
                    continue
                if len(orientation)!=0 and orientation not in orient:
                    orient.append(orientation)
               
                if smartsindices not in currentorientation:
                    currentorientation.append(smartsindices)
                
                for smartsindex,moleculeindex in smartsindextomoleculeindex.items():
                    if smartsindex not in smartsindextomoleculeindices.keys(): 
                        smartsindextomoleculeindices[smartsindex]=[]
                    if moleculeindex not in smartsindextomoleculeindices[smartsindex]:
                        smartsindextomoleculeindices[smartsindex].append(moleculeindex) 
        moleculeindextosmartsindices=ReverseSpecialDictionary(poltype,smartsindextomoleculeindices)
        structure = Chem.MolFromSmarts(parametersmarts)
        substructure = Chem.MolFromSmarts(smartsfortransfer)
        matches=structure.GetSubstructMatches(substructure)
        trueparametersmartsindextosmartsindices={}
        currentorientation=[]
        for match in matches:
            indices=list(range(len(match)))
            smartsindextoparametersmartsindex=dict(zip(indices,match)) 

            orientation,keep=CheckIfOrientationIsSameAsPrevious(poltype,currentorientation,smartsindextoparametersmartsindex,moleculeindextosmartsindices,atomindices)
            if orientation not in currentorientation:
                currentorientation.append(orientation)
            if keep==False:
                continue
            parametersmartsindextosmartsindex={v: k for k, v in smartsindextoparametersmartsindex.items()}
            for parametersmartsindex,smartsindex in parametersmartsindextosmartsindex.items():
                if parametersmartsindex not in trueparametersmartsindextosmartsindices.keys():
                    trueparametersmartsindextosmartsindices[parametersmartsindex]=[]
                if smartsindex not in trueparametersmartsindextosmartsindices[parametersmartsindex]:
                    trueparametersmartsindextosmartsindices[parametersmartsindex].append(smartsindex)

        smartsindextoparametersmartsindices=ReverseSpecialDictionary(poltype,trueparametersmartsindextosmartsindices)
        allpossiblebabelindicesprmsmarts=GrabSymmetryIdenticalIndicesBabel(poltype,structure,moleculeindextosmartsindices,smartsindextoparametersmartsindices,atomindices)
        for i in range(len(allpossiblebabelindicesprmsmarts)):
            babelindicessmarts=allpossiblebabelindicesprmsmarts[i]
            rdkitindicessmarts=[j-1 for j in babelindicessmarts]
            indicesindic=[]
            indicesnotindic=[]
            for idx in rdkitindicessmarts: # sometimes type classes for three oxygen are same but two are double bonded so matches will not match single bond
                if idx in trueparametersmartsindextosmartsindices.keys():
                    indicesindic.append(idx)
                else:
                    indicesnotindic.append(idx)
            for idx in indicesnotindic:
                trueparametersmartsindextosmartsindices[idx]=trueparametersmartsindextosmartsindices[indicesindic[0]]
        smartsindextoparametersmartsindices=ReverseSpecialDictionary(poltype,trueparametersmartsindextosmartsindices)

                
        
         
        allpossiblerdkitindices=ReturnSymmetryIdenticalIndicesSMARTSAndParameterSMARTS(poltype,allpossiblebabelindicesprmsmarts,smartsindextomoleculeindices,trueparametersmartsindextosmartsindices)
        wildcards=CountWildCards(poltype,smartsfortransfer)
        parametersmartsatomorderlists=[]
        masterdic={}
        for parametersmartsatomorderlist,elementtinkerdescrip in smartsatomordertoelementtinkerdescrip.items():
            prmsmarts=parametersmartsatomorderlist[0]
            atomorderlist=parametersmartsatomorderlist[1]
            if prmsmarts==parametersmarts:
                parametersmartsatomorderlists.append(parametersmartsatomorderlist)
        for rdkitindex in range(len(allpossiblerdkitindices)):
            rdkitindexlist=allpossiblerdkitindices[rdkitindex]
            for parametersmartsatomorderlist in parametersmartsatomorderlists:
                atomorderlist=parametersmartsatomorderlist[1]
                elementtinkerdescrip=smartsatomordertoelementtinkerdescrip[parametersmartsatomorderlist]
                for atomorder in atomorderlist:
                    for idx in rdkitindexlist:
                        smartsindices=moleculeindextosmartsindices[idx]
                        for smartsindex in smartsindices:
                            smartsorder=smartsindex+1
                            parametersmartsindices=smartsindextoparametersmartsindices[smartsindex]
                            for parametersmartsatomindex in parametersmartsindices: 
                                parametersmartsatomorder=parametersmartsatomindex+1
                                if parametersmartsatomorder==atomorder:
                                    if rdkitindex not in masterdic.keys():
                                        masterdic[rdkitindex]={}
                                    molindextodescrips=masterdic[rdkitindex]
                                    if idx not in molindextodescrips.keys():
                                         molindextodescrips[idx]={}
                                    dic=molindextodescrips[idx]
                                    if elementtinkerdescrip not in dic.keys():
                                        dic[elementtinkerdescrip]={}
                                    innerdic=dic[elementtinkerdescrip]

                                    atomtinkertype=elementtinkerdescriptotinkertype[elementtinkerdescrip]
                                    atomtinkerclass=tinkertypetoclass[atomtinkertype]
                                    innerdic['prmsmartsatomorder']=parametersmartsatomorder
                                    innerdic['smartsatomorder']=smartsorder
                                    innerdic['tinkertype']=int(atomtinkertype)
                                    innerdic['tinkerclass']=int(atomtinkerclass)
                                    dic[elementtinkerdescrip]=innerdic
                                    molindextodescrips[idx]=dic
                                    masterdic[rdkitindex]=molindextodescrips
                                if len(atomindices)==wildcards:
                                    break
        arraydics=list(masterdic.values())
        allpossiblemolidxs=[]
        for dic in arraydics:
            allpossiblemolidxs.append(list(set(list(dic.keys()))))
        allpossibledescrips=[]
        firstmolidxs=GrabFirstInstanceAndCheckIfNotAlreadyUsed(poltype,allpossiblemolidxs)
        
        nextdics=[]
        for dicidx in range(len(arraydics)):
            dic=arraydics[dicidx]
            molidx=firstmolidxs[dicidx]
            nextdic=dic[molidx]
            nextdics.append(nextdic)
            allpossibledescrips.append(list(set(list(nextdic.keys()))))
             
        firstdescrips=GrabFirstInstanceAndCheckIfNotAlreadyUsed(poltype,allpossibledescrips)
        finaldics=[]
        finaldescrips=[]
        for dicidx in range(len(nextdics)):
            nextdic=nextdics[dicidx]
            descrip=firstdescrips[dicidx]
            finaldescrips.append(descrip)
            innerdic=nextdic[descrip]
            finaldics.append(innerdic)
        if len(finaldescrips)!=len(atomindices):
            raise ValueError('Missing SMARTS-Atom Order -> Tinker Description, Element entry for smarts '+smartsfortransfer+' and matching to parameter smarts '+parametersmarts+' for indices '+str(atomindices))
        atomindicestotinkertypes[atomindices]=[dic['tinkertype'] for dic in finaldics]
        atomindicestotinkerclasses[atomindices]=[dic['tinkerclass'] for dic in finaldics]
        atomindicestoparametersmartsatomorders[atomindices]=[parametersmarts,[dic['prmsmartsatomorder'] for dic in finaldics]]
        atomindicestoelementtinkerdescrips[atomindices]=finaldescrips
        atomindicestosmartsatomorders[atomindices]=[smarts,[dic['smartsatomorder'] for dic in finaldics]]
    return atomindicestotinkertypes,atomindicestotinkerclasses,atomindicestoparametersmartsatomorders,atomindicestoelementtinkerdescrips,atomindicestosmartsatomorders

def GrabFirstInstanceAndCheckIfNotAlreadyUsed(poltype,allpossible):
    alreadyused=[]
    firstinstance=[]
    for ls in allpossible:
        first=ls[0]
        if first not in alreadyused:
            alreadyused.append(first)
            firstinstance.append(first)
        else:
            if len(ls)>1:
                second=ls[1]
                alreadyused.append(second)
                firstinstance.append(second)
            else:
                alreadyused.append(first)
                firstinstance.append(first)
    return firstinstance


def ReverseSpecialDictionary(poltype,keytoarray):
    reversedkeytoarray={} 
    smartsindextoparametersmartsindices={}
    for key,array in keytoarray.items():
        for index in array:
            if index not in reversedkeytoarray.keys():
                reversedkeytoarray[index]=[]
            if key not in reversedkeytoarray[index]:
                reversedkeytoarray[index].append(key)
    return reversedkeytoarray


def CountWildCards(poltype,string):
    count=0
    for e in string:
        if e=='*':
            count+=1
    return count

def ReturnSymmetryIdenticalIndicesSMARTSAndParameterSMARTS(poltype,allpossiblebabelindicessmarts,smartsindextomoleculeindices,trueparametersmartsindextosmartsindices):
    allpossiblerdkitindices=[]
    for i in range(len(allpossiblebabelindicessmarts)):
        babelindicessmarts=allpossiblebabelindicessmarts[i]
        rdkitindicessmarts=[j-1 for j in babelindicessmarts]

        indices=[]
        allpossiblerdkitindicessmarts=[trueparametersmartsindextosmartsindices[j] for j in rdkitindicessmarts]
        for rdkitindicessmarts in allpossiblerdkitindicessmarts:
            allpossiblerdkitindicesmolsmarts=[smartsindextomoleculeindices[j] for j in rdkitindicessmarts]
            for rdkitindicesmolsmarts in allpossiblerdkitindicesmolsmarts:
               for index in rdkitindicesmolsmarts:
                   if index not in indices:
                       indices.append(index)
        allpossiblerdkitindices.append(indices)
    return allpossiblerdkitindices


def GrabSymmetryIdenticalIndicesBabel(poltype,mol,moleculeindextosmartsindices,smartsindextoparametersmartsindices,atomindices,show=None):
    show=True
    fragmentfilepath='fragment.mol'
    if os.path.isfile(fragmentfilepath):
        os.remove(fragmentfilepath)
    rdmolfiles.MolToMolFile(mol,fragmentfilepath)
    obConversion = openbabel.OBConversion()
    fragbabelmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(fragmentfilepath)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(fragbabelmol, fragmentfilepath)
    fragidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,fragbabelmol)
    allpossiblefragsmartsidxs=[moleculeindextosmartsindices[i] for i in atomindices]
    allpossiblebabelindicestotal=[]
    for smartsindices in allpossiblefragsmartsidxs: 
        allpossibleparametersmartsrdkitindices=[smartsindextoparametersmartsindices[i] for i in smartsindices]
        flatlist=[]
            
        for prmsmartsrdkitindices in allpossibleparametersmartsrdkitindices: 
            babelindices=[i+1 for i in prmsmartsrdkitindices]
            symclasses=[fragidxtosymclass[babelindex] for babelindex in babelindices]
            allpossiblebabelindices=[GrabKeysFromValue(poltype,fragidxtosymclass,symclass) for symclass in symclasses]
            for ls in allpossiblebabelindices:
                for idx in ls:
                    flatlist.append(idx)

        allpossiblebabelindicestotal.append(flatlist)
    return allpossiblebabelindicestotal

def GrabSMARTSList(poltype,smartsatomordertoelementtinkerdescrip):
    smartslist=[]
    for smartsatomorder in smartsatomordertoelementtinkerdescrip.keys():
        smarts=smartsatomorder[0]
        if smarts not in smartslist:
            smartslist.append(smarts)
    return smartslist 

def CheckForPlanerAngles(poltype,listofanglesforprm,mol):
    listofanglesthatneedplanarkeyword=[]
    shoulduseanglep = version.parse(poltype.versionnum) >= version.parse("8.7")
    for ls in listofanglesforprm:
        a = mol.GetAtom(ls[0]+1)
        b = mol.GetAtom(ls[1]+1)
        c = mol.GetAtom(ls[2]+1)
        anglep=False
        if b.GetHyb()==2 and shoulduseanglep==True: # only for SP2 hyb middle atoms use angp
            neighbs=list(openbabel.OBAtomAtomIter(b))
            if len(neighbs)==3:
                anglep=True
        if anglep==True:
            listofanglesthatneedplanarkeyword.append(ls)
    return listofanglesthatneedplanarkeyword

def ModifyAngleKeywords(poltype,angleprms,listofanglesthatneedplanarkeywordtinkerclassestopoltypeclasses):
    newangleprms=[]
    for line in angleprms:
        found=False
        linesplit=line.split()
        temp=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
        for ls,polclassesls in listofanglesthatneedplanarkeywordtinkerclassestopoltypeclasses.items():
            inline=True
            for i in temp:
                if i not in polclassesls:
                    inline=False
            if inline==True:
                found=True
                if 'anglep' not in line:
                    newline=line.replace('angle','anglep')
                else:
                    newline=line
                break
        if found==False:
            newline=line.replace('anglep','angle')

        newangleprms.append(newline)
    return newangleprms

def FilterList(poltype,allitems,listbabel):
    newallitems=[]
    for ls in allitems:
        revls=ls[::-1]
        if list(ls) in listbabel or list(revls) in listbabel:
            newallitems.append(ls)
    return newallitems

def AddOptimizedBondLengths(poltype,mol,bondprms,bondlistbabel):
    newbondprms=[]
    for line in bondprms:
        linesplit=line.split()
        bondtypes=[int(linesplit[1]),int(linesplit[2])]
        bondindices=[]
        for prmtype in bondtypes:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,prmtype)
            bondindices.append(keylist)
        allbonds = list(itertools.product(bondindices[0], bondindices[1]))
        allbonds=[x for x in allbonds if len(x) == len(set(x))]
        allbonds=FilterList(poltype,allbonds,bondlistbabel)
        tot=0

        for bond in allbonds:
            try:
                blen = mol.GetBond(bond[0],bond[1]).GetLength()
                tot+=blen
            except:
                pass
        avgbondlength=round(tot/len(allbonds),2)
        linesplit=re.split(r'(\s+)', line)   
        linesplit[8]=str(avgbondlength)
        line=''.join(linesplit)
        newbondprms.append(line)
    return newbondprms
      

def AddOptimizedAngleLengths(poltype,mol,angleprms,anglelistbabel):
    newangleprms=[]
    for line in angleprms:
        linesplit=line.split()
        angletypes=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
        angleindices=[]
        for prmtype in angletypes:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,prmtype)
            angleindices.append(keylist)
        allangles = list(itertools.product(angleindices[0], angleindices[1],angleindices[2]))
        allangles=[x for x in allangles if len(x) == len(set(x))]
        allangles=FilterList(poltype,allangles,anglelistbabel)

        tot=0
        for angle in allangles:
            a = mol.GetAtom(angle[0])
            b = mol.GetAtom(angle[1])
            c = mol.GetAtom(angle[2])
            anglelen = mol.GetAngle(a,b,c)
            tot+=anglelen
        avganglelength=round(tot/len(allangles),2)
        linesplit=re.split(r'(\s+)', line)
        linesplit=linesplit[:11]
        linesplit.append('\n')
        linesplit[10]=str(avganglelength)
        line=''.join(linesplit)
        newangleprms.append(line)
    return newangleprms
 


def GrabKeysFromValue(poltype,dic,thevalue):
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist
          

def CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol):
    allhydrogentor=True
    atomobjects=[mol.GetAtom(i) for i in babelindices]
    a,b,c,d=atomobjects[:]
    aidx,bidx,cidx,didx=babelindices[:]
    aatomicnum=a.GetAtomicNum()
    datomicnum=d.GetAtomicNum()
    if aatomicnum!=1 or datomicnum!=1:
        allhydrogentor=False
    else:
        torlist=[]
        iteratomatom = openbabel.OBAtomAtomIter(b)
        for iaa in iteratomatom:
            iteratomatom2 = openbabel.OBAtomAtomIter(c)
            for iaa2 in iteratomatom2:
                ta = iaa.GetIdx()
                tb = bidx
                tc = cidx
                td = iaa2.GetIdx()
                if ((ta != tc and td != tb) and not (ta == aidx and td == didx)):
                    torlist.append([ta,tb,tc,td])
        for tor in torlist:
            atoms=[mol.GetAtom(i) for i in tor]
            aatomnum=atoms[0].GetAtomicNum()
            datomnum=atoms[3].GetAtomicNum()
            if aatomnum!=1 or datomnum!=1:
                allhydrogentor=False
    return allhydrogentor    




def FindMissingTorsions(poltype,torsionindicestoparametersmartsenv,rdkitmol,mol,indextoneighbidxs):
    torsionsmissing=[]
    for torsionindices,smartsenv in torsionindicestoparametersmartsenv.items():
        aidx,bidx,cidx,didx=torsionindices[:]
        babelindices=[i+1 for i in torsionindices]
        abidx,bbidx,cbidx,dbidx=babelindices[:]
        if len(poltype.onlyrotbndslist)!=0:
            if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                if torsionindices not in torsionsmissing:
                    torsionsmissing.append(torsionindices)
                    continue
        bond=mol.GetBond(bbidx,cbidx)
        bondorder=bond.GetBondOrder()
        if bondorder!=1: # then dont zero out
            continue 
        babelatoms=[mol.GetAtom(i) for i in babelindices]
        ringbools=[a.IsInRing() for a in babelatoms]
        contin=False
        bnd=[babelindices[1],babelindices[2]]
        if bnd in poltype.partialdoublebonds or bnd[::-1] in poltype.partialdoublebonds:
            continue 
        ringb=ringbools[1]
        ringc=ringbools[2]
        if ringb==True and ringc==True:
            contin=True

        if contin==True:
            continue
        allhydrogentor=CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol)
        firstneighborindexes=indextoneighbidxs[bidx]
        secondneighborindexes=indextoneighbidxs[cidx]
        neighborindexes=firstneighborindexes+secondneighborindexes
        smarts=smartsenv[0]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)

        matcharray=[]
        for match in matches:
            for idx in match:
                if idx not in matcharray:
                    matcharray.append(idx)
            
        check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
        atoma=rdkitmol.GetAtomWithIdx(aidx)
        atomd=rdkitmol.GetAtomWithIdx(didx)
        atomicnumatoma=atoma.GetAtomicNum()
        atomicnumatomd=atomd.GetAtomicNum()
        if poltype.transferanyhydrogentor==True: 
            if atomicnumatoma==1 or atomicnumatomd==1:
                if allhydrogentor==True:
                    if torsionindices not in torsionsmissing:
                        torsionsmissing.append(torsionindices)
                else:
                    continue
        if '~' in smarts or '*' in smarts:
            if torsionindices not in torsionsmissing:
                torsionsmissing.append(torsionindices)
            continue

        if check==False:
            if torsionindices not in torsionsmissing:
                torsionsmissing.append(torsionindices)
    return torsionsmissing 


def FindAllNeighborIndexes(poltype,rdkitmol):
    indextoneighbidxs={}
    for atm in rdkitmol.GetAtoms():
        atmidx=atm.GetIdx()
        if atmidx not in indextoneighbidxs.keys():
            indextoneighbidxs[atmidx]=[]
        for neighbatm in atm.GetNeighbors():
            neighbatmidx=neighbatm.GetIdx()
            if neighbatmidx not in indextoneighbidxs[atmidx]:
                indextoneighbidxs[atmidx].append(neighbatmidx)
            

    return indextoneighbidxs

def CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,smartsindexes):
    check=True
    for idx in neighborindexes:
        if idx not in smartsindexes:
            check=False
    return check

def ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms):
    newtorsionprms=[]
    for line in torsionprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])])
        for sublist in torsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                splitafter=newlinesplit[12:]
                splitafter[0]='0'
                splitafter[4]='0'
                splitafter[10]='0'
                newlinesplit=newlinesplit[:12]+splitafter
                line=''.join(newlinesplit)
        newtorsionprms.append(line)

    return newtorsionprms


def GrabTorsionParameterCoefficients(poltype,torsionprms):
    torsionkeystringtoparameters={}
    for line in torsionprms:
        linesplit=line.split() 
        key = '%s %s %s %s' % (linesplit[1], linesplit[2], linesplit[3], linesplit[4])
        parameters=[float(linesplit[5]),float(linesplit[8]),float(linesplit[11])]
        torsionkeystringtoparameters[key]=parameters
       
    return torsionkeystringtoparameters

def PruneDictionary(poltype,keysubset,dic):
    newdic={}
    for key,value in dic.items():
        if key in keysubset:
            newdic[key]=value
    return newdic


def TinkerClassesToPoltypeClasses(poltype,indicestotinkerclasses):
    tinkerclassestopoltypeclasses={}
    for indices,tinkerclasses in indicestotinkerclasses.items():
        babelindices=[i+1 for i in indices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        if len(indices)>1:
            first=int(tinkerclasses[0])
            second=int(tinkerclasses[-1])
            if first>second:
                pass
            elif first<second:
                tinkerclasses=tinkerclasses[::-1]
                poltypeclasses=poltypeclasses[::-1]
            else:
                if len(indices)==4:
                    first=int(tinkerclasses[0])
                    second=int(tinkerclasses[1])
                    third=int(tinkerclasses[2])
                    fourth=int(tinkerclasses[3])
                    firstsum=first+second
                    secondsum=third+fourth
                    if firstsum>secondsum:
                        pass
                    elif firstsum<secondsum:
                        tinkerclasses=tinkerclasses[::-1]
                        poltypeclasses=poltypeclasses[::-1]

                     
        if tuple(tinkerclasses) not in tinkerclassestopoltypeclasses.keys():
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)]=[]
        if tuple(poltypeclasses) not in tinkerclassestopoltypeclasses[tuple(tinkerclasses)]: 
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)].append(tuple(poltypeclasses))
    return tinkerclassestopoltypeclasses


def ConvertIndicesDictionaryToPoltypeClasses(poltype,indicestovalue,indicestotinkerclasses,tinkerclassestopoltypeclasses):
    poltypeclassestovalue={}
    for indices,value in indicestovalue.items():
        tinkerclasses=tuple(indicestotinkerclasses[indices])
        if tinkerclasses in tinkerclassestopoltypeclasses.keys():
            poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses])
        else:
            poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses[::-1]])

        poltypeclassestovalue[poltypeclasses]=value
    return poltypeclassestovalue 


def SearchForPoltypeClasses(poltype,prmline,poltypeclasseslist):
    for listofpoltypeclasses in poltypeclasseslist:
        for poltypeclasses in listofpoltypeclasses:
            allin=True
            for poltypeclass in poltypeclasses:
                if str(poltypeclass) not in prmline:
                    allin=False 
            if allin==True:
                return listofpoltypeclasses

def MapParameterLineToTransferInfo(poltype,prms,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,poltypeclassestosmartsatomordersext,defaultvalues=None):
    prmstotransferinfo={}
    for line in prms:
        poltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestoparametersmartsatomorders.keys())
        if poltypeclasses==None:
            poltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())
            smartsatomorders=poltypeclassestosmartsatomordersext[poltypeclasses]
            transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' from external database'+'\n'
        else:
            parametersmartsatomorders=poltypeclassestoparametersmartsatomorders[poltypeclasses]
            smartsatomorders=poltypeclassestosmartsatomorders[poltypeclasses]
            elementtinkerdescrips=poltypeclassestoelementtinkerdescrips[poltypeclasses]
            transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' '+'to SMARTS from parameter file'+' '+str(parametersmartsatomorders)+' '+'with tinker type descriptions '+str(elementtinkerdescrips)+'\n'
            warn=False
            if defaultvalues!=None:
                for value in defaultvalues:
                    if str(value) in line:
                        warn=True
            if warn==True:
                transferinfoline+='# '+'WARNING DEFAULT MM3 OPBEND VALUES USED '+'\n'
            smarts=smartsatomorders[0]
            if '~' in smarts or '*' in smarts:
                transferinfoline+='# '+'WARNING WILDCARDS USED IN SMARTS PARAMETER MATCHING'+'\n'
                
        prmstotransferinfo[line]=transferinfoline

    return prmstotransferinfo

def FindMaximumCommonSubstructures(poltype,parametersmartslist,rdkitmol):
    parametersmartstomaxcommonsubstructure={}
    maxatomsize=0
    for parametersmarts in parametersmartslist:
        mols = [rdkitmol,Chem.MolFromSmarts(parametersmarts)]
        res=rdFMCS.FindMCS(mols)
        atomnum=res.numAtoms
        smartsmcs=res.smartsString
        if atomnum>0:
            if atomnum>maxatomsize:
                maxatomsize=atomnum
            parametersmartstomaxcommonsubstructure[parametersmarts]=smartsmcs
    return parametersmartstomaxcommonsubstructure,maxatomsize 

def GrabPlanarBonds(poltype,listofbondsforprm,mol): # used for checking missing opbend parameters later
    planarbonds=[]
    for bond in listofbondsforprm:
        revbond=bond[::-1]
        totalbonds=[bond,revbond]
        for lsidx in range(len(totalbonds)):
            bond=totalbonds[lsidx]
            a = mol.GetAtom(bond[0]+1)
            b = mol.GetAtom(bond[1]+1)
            ainring=a.IsInRing()
            binring=b.IsInRing()
            aisaromatic=a.IsAromatic()
            bisaromatic=b.IsAromatic()
            if ainring==True and binring==True:
                if aisaromatic==True or bisaromatic==True:
                    if bond not in planarbonds:
                        planarbonds.append(bond)
                else:
                    if b.GetHyb()==2 and len(list(openbabel.OBAtomAtomIter(b)))==3:
                        if bond not in planarbonds:
                            planarbonds.append(bond)

            else:
                if b.GetHyb()==2 and len(list(openbabel.OBAtomAtomIter(b)))==3:
                    if bond not in planarbonds:
                        planarbonds.append(bond)
    return planarbonds 


def FindPotentialMissingParameterTypes(poltype,prms,tinkerclassestopoltypeclasses):
    poltypeclasseslist=list(tinkerclassestopoltypeclasses.values())
    newpoltypeclasseslist=[]
    for poltypeclassesls in poltypeclasseslist:
        newpoltypeclassesls=[]
        for poltypeclasses in poltypeclassesls:
            revpoltypeclasses=poltypeclasses[::-1]
            newpoltypeclassesls.append(poltypeclasses)
            newpoltypeclassesls.append(revpoltypeclasses)
        newpoltypeclasseslist.append(newpoltypeclassesls)

    missingprms=[]
    foundprms=[]
    for line in prms:
        poltypeclassesls=SearchForPoltypeClasses(poltype,line,newpoltypeclasseslist)
        linesplit=line.split()
        prms=[int(linesplit[1]),int(linesplit[2])]
        if poltypeclassesls!=None:
            for poltypeclasses in poltypeclassesls:
                poltypeclasses=tuple([int(i) for i in poltypeclasses])
                if poltypeclasses[0]==prms[0] and poltypeclasses[1]==prms[1]:
                    foundprms.append(poltypeclasses) 
    for poltypeclassesls in newpoltypeclasseslist:
        for poltypeclasses in poltypeclassesls:
            if poltypeclasses not in foundprms and poltypeclasses not in missingprms:
                missingprms.append(poltypeclasses)
 
    return missingprms 


def ConvertPoltypeClassesToIndices(poltype,missingprmtypes):
    missingprmindices=[]
    for poltypeclasses in missingprmtypes:
        indices=[]
        for poltypeclass in poltypeclasses:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
            keylist=[i-1 for i in keylist]
            indices.append(keylist)
        missingprmindices.append(indices)
    finallist=[]
    for indices in missingprmindices:
        combs = list(itertools.product(*indices))
        for comb in combs:
            finallist.append(comb)
    return finallist

def FilterIndices(poltype,potentialmissingprmindices,indices):
    missingprmindices=[]
    for prmindices in potentialmissingprmindices:
        if list(prmindices) in indices:
            missingprmindices.append(prmindices)
    return missingprmindices
        

def GrabTypeAndDescriptions(poltype,prmfile):
    elementmm3descriptomm3type={}
    temp=open(prmfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line and 'piatom' not in line:
            linesplit=line.split()    
            newlinesplit=linesplit[1:-3]
            mm3type=newlinesplit[0]
            element=newlinesplit[1]
            mm3descrip=' '.join(newlinesplit[2:])
            ls=[element,mm3descrip]
            elementmm3descriptomm3type[tuple(ls)]=mm3type
    return elementmm3descriptomm3type 


def DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools):
    newopbendprms=[]
    defaultvalues=[]
    for opbendprmindices in missingopbendprmindices:
        babelindices=[i+1 for i in opbendprmindices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        for atomindex in babelindices:
            atom=mol.GetAtom(atomindex)
            neighbs=list(openbabel.OBAtomAtomIter(atom))
            if len(neighbs)==3:
                atomnum=atom.GetAtomicNum()
                if atomnum==6:
                    match=False
                    for natom in neighbs:
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==8:
                            bnd=mol.GetBond(atom,natom)
                            BO=bnd.GetBondOrder()
                            if BO==2:
                                match=True
                    if match==True:
                        opbendvalue=round(1*71.94,2)
                    else:
                        opbendvalue=round(.2*71.94,2)
                elif atomnum==7:
                    count=0
                    for natom in neighbs:
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==8:
                           count+=1
                    if count==2:
                        opbendvalue=round(1.5*71.94,2)
                    else:
                        opbendvalue=round(.05*71.94,2)
                else:
                    opbendvalue=round(poltype.defopbendval*71.94,2)

                defaultvalues.append(opbendvalue)
                boolarray=opbendbondindicestotrigonalcenterbools[opbendprmindices]
                firstprmline='opbend'+' '+str(poltypeclasses[0])+' '+str(poltypeclasses[1])+' '+'0'+' '+'0'+' '+str(opbendvalue)+'\n'

                if boolarray[1]==True: 
                    newopbendprms.append(firstprmline)

    return newopbendprms,defaultvalues


def WriteOutList(poltype,ls,filename):
    np.savetxt(filename,ls)

def ReadList(poltype,filename):
    newls=[]
    if os.stat(filename).st_size != 0:
        ls=np.loadtxt(filename)
        if ls.ndim==1: # read only in to list, but need list of list
            ls=[ls.tolist()]       
        for subls in ls:
            newsubls=[int(i) for i in subls]
            newls.append(newsubls)
    return newls

def CheckIfParametersExist(poltype,potentialmissingindices,prms):
    missingprmindices=[]
    for indices in potentialmissingindices:
        babelindices=[i+1 for i in indices]
        symtypes=[poltype.idxtosymclass[i] for i in babelindices]
        found=False
        symtypes=[str(i) for i in symtypes]
        for prmline in prms:
            linesplit=prmline.split()
            parms=[linesplit[1],linesplit[2]]
            if parms[0]==symtypes[0] and parms[1]==symtypes[1]:
               found=True
        if found==False:
            missingprmindices.append(indices)
    return missingprmindices
                
def WriteDictionaryToFile(poltype,dic,filename):
    json.dump(dic, open(filename,'w'))        

def ReadDictionaryFromFile(poltype,filename):
    return json.load(open(filename))

def GenerateSMARTSMatchLine(poltype,rdkitmol,rdkitindex):
    smarts=rdmolfiles.MolToSmarts(rdkitmol)
    smartsmol=Chem.MolFromSmarts(smarts)
    matches=rdkitmol.GetSubstructMatches(smartsmol)
    for match in matches:
        lastmatch=match
    lastindex=lastmatch.index(rdkitindex)
    atomorder=lastindex+1
    string='%'+' '+smarts+' % '+str(atomorder)
    return string

def GenerateAtomSMARTSMap(poltype,rdkitmol):
    lines=[]
    descrip='"%s"'% poltype.molecprefix
    for atom in rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()
        symbol=atom.GetSymbol()
        string=GenerateSMARTSMatchLine(poltype,rdkitmol,atomidx)
        line=' '+symbol+' '+descrip+' '+string+'\n'
        lines.append(line)
    return lines        

def AppendToSMARTSMapFile(poltype,lines,filename):
    temp=open(filename,'a')
    for line in lines:
        temp.write(line)
    temp.close()


def RepresentsInt(poltype,s):
    try: 
        int(s)
        return True
    except ValueError:
        return False



def AddKeyFileParametersToParameterFile(poltype,rdkitmol):   
    atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms=GrabParameters(poltype,poltype.keyfiletoaddtodatabase)
    oldtypetonewtype,shift=ShiftPoltypeNumbers(poltype,poltype.smallmoleculeprmlib,poltype.keyfiletoaddtodatabase)
    result=ShiftParameterDefintions(poltype,[atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms],oldtypetonewtype)
    atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms=result[:] 
    WriteToPrmFile(poltype,atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms,poltype.smallmoleculeprmlib)
    lines=GenerateAtomSMARTSMap(poltype,rdkitmol)
    AppendToSMARTSMapFile(poltype,lines,poltype.smallmoleculesmartstotinkerdescrip)

def ReadExternalDatabase(poltype):
    temp=open(poltype.externalparameterdatabase,'r')
    results=temp.readlines()
    temp.close()
    bondsmartsatomordertoparameters={} 
    anglesmartsatomordertoparameters={}
    strbndsmartsatomordertoparameters={}
    torsionsmartsatomordertoparameters={}
    opbendsmartsatomordertoparameters={}
    vdwsmartsatomordertoparameters={}
    for line in results:
        linesplit=line.split()
        if linesplit[0]=='#':
            continue
        keyword=linesplit[0]
        linesplit=linesplit[1:]
        newline=' '.join(linesplit)
        linesplit=newline.split('%')
        smarts=linesplit[1].lstrip().rstrip()
        atomorderstring=linesplit[2]
        atomorderstringlist=atomorderstring.split(',')
        atomorderlist=tuple([int(i) for i in atomorderstringlist])
        prmstring=linesplit[3]
        prmstringlist=prmstring.split(',')
        prmlist=[float(i) for i in prmstringlist] 
        smartsatomorder=tuple([smarts,atomorderlist])
        if keyword=='bond':
            bondsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='angle':
            anglesmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='strbnd':
            strbndsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='torsion':
            torsionsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='opbend':
            opbendsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='vdw':
            vdwsmartsatomordertoparameters[smartsatomorder]=prmlist
    return bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters


def ConvertToPoltypeClasses(poltype,torsionsmissing):
    newtorsionsmissing=[]
    for sublist in torsionsmissing:
        newsublist=[i+1 for i in sublist]
        a,b,c,d=newsublist[:]
        sorttor=torfit.sorttorsion(poltype,[poltype.idxtosymclass[a],poltype.idxtosymclass[b],poltype.idxtosymclass[c],poltype.idxtosymclass[d]])
        if sorttor not in newtorsionsmissing:
            newtorsionsmissing.append(sorttor)
    return newtorsionsmissing


def MatchExternalSMARTSToMolecule(poltype,rdkitmol,smartsatomordertoparameters):
    indicestoextsmartsmatchlength={}
    indicestoextsmarts={}
    for smartsatomorder,parameters in smartsatomordertoparameters.items():
        smarts=smartsatomorder[0]
        atomorderlist=smartsatomorder[1]
        substructure = Chem.MolFromSmarts(smarts)
        diditmatch=rdkitmol.HasSubstructMatch(substructure)
        if diditmatch==True:
            matches=rdkitmol.GetSubstructMatches(substructure)
            for match in matches:
                indices=list(range(len(match)))
                smartsindextomoleculeindex=dict(zip(indices,match)) 
                smartsindexlist=[i-1 for i in atomorderlist] 
                moleculeindices=tuple([smartsindextomoleculeindex[i] for i in smartsindexlist])
                if moleculeindices not in indicestoextsmartsmatchlength.keys():
                    indicestoextsmartsmatchlength[moleculeindices]=len(match)
                    indicestoextsmarts[moleculeindices]=smarts
    return indicestoextsmartsmatchlength,indicestoextsmarts

def CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,indicestoextsmartsmatchlength,indicesforprmtoparametersmarts,indicesforprmtosmarts,indicestoextsmarts):
    newindicestoextsmarts={}
    for indices,extsmartsmatchlength in indicestoextsmartsmatchlength.items():
        if indices in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices]
            smarts=smartsls[0]
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms:
                del indicesforprmtosmarts[indices]
                del indicesforprmtoparametersmarts[indices]
                newindicestoextsmarts[indices]=indicestoextsmarts[indices]
        elif indices[::-1] in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices[::-1]]
            smarts=smartsls[0]
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms:
                del indicesforprmtosmarts[indices[::-1]]
                del indicesforprmtoparametersmarts[indices[::-1]]
                newindicestoextsmarts[indices[::-1]]=indicestoextsmarts[indices[::-1]]

    return indicesforprmtoparametersmarts,indicesforprmtosmarts,newindicestoextsmarts

def AddExternalDatabaseSMARTSMatchParameters(poltype,prms,indicestoextsmarts,smartsatomordertoparameters,keyword):
    poltypeclassestosmartsatomordersext={}
    for indices,extsmarts in indicestoextsmarts.items():
        for smartsatomorder,parameters in smartsatomordertoparameters.items():
            smarts=smartsatomorder[0]
            if smarts==extsmarts:
                line=keyword+' '
                babelindices=[i+1 for i in indices]
                poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
                for poltypeclass in poltypeclasses:
                    line+=str(poltypeclass)+' '
                if keyword=='opbend':
                    line+='0 0 '
                for prm in parameters:
                    line+=str(prm)+' '
                line+='\n' 
                prms.append(line)
                poltypeclassestosmartsatomordersext[tuple([tuple(poltypeclasses)])]=smartsatomorder
    return prms,poltypeclassestosmartsatomordersext


def FilterBondSMARTSEnviorment(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses):
    opbendbondindicestosmartsatomorders={}
    opbendbondindicestotinkerclasses={}
    for bondindices,smartsatomorders in bondindicestosmartsatomorders.items(): 
        smarts=smartsatomorders[0]
        brackcount=CountBrackets(poltype,smarts)
        revbondindices=bondindices[::-1]
        if brackcount==2:
            pass
        else:
            temp=bondindicestotinkerclasses[bondindices]
            revtemp=temp[::-1]
            opbendbondindicestotinkerclasses[bondindices]=temp
            opbendbondindicestotinkerclasses[revbondindices]=revtemp
            opbendbondindicestosmartsatomorders[bondindices]=smartsatomorders 
            opbendbondindicestosmartsatomorders[revbondindices]=smartsatomorders 

    return opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders

def CountBrackets(poltype,string):
    count=0
    for e in string:
        if e=='[':
            count+=1
    return count

def CheckTrigonalCenters(poltype,listofbondsforprm,mol):
    opbendbondindicestotrigonalcenterbools={}
    for bondindices in listofbondsforprm:
        bondindices=tuple(bondindices)
        boolarray=[]
        babelindices=[i+1 for i in bondindices]
        atoms=[mol.GetAtom(i) for i in babelindices]     
        for a in atoms:
            hyb=a.GetHyb()
            neighbs=list(openbabel.OBAtomAtomIter(a))
            if len(neighbs)==3 and hyb==2:
                boolarray.append(True)
            else:
                boolarray.append(False)
        opbendbondindicestotrigonalcenterbools[bondindices]=boolarray
        revboolarray=boolarray[::-1]
        revbondindices=bondindices[::-1]
        opbendbondindicestotrigonalcenterbools[revbondindices]=revboolarray
    return opbendbondindicestotrigonalcenterbools


def CorrectPitorEnergy(poltype,torsionprms,torsiontopitor):
    newtorsionprms=[]
    for torline in torsionprms:
        torlinesplit=torline.split()
        tor=tuple([int(torlinesplit[1]),int(torlinesplit[2]),int(torlinesplit[3]),int(torlinesplit[4])])
        if tor in torsiontopitor:
            pitorline=torsiontopitor[tor]
            pitorlinesplit=pitorline.split()
            prm=float(pitorlinesplit[3])
            torprm=float(torlinesplit[8])
            newtorprm=prm+torprm
            torlinesplit[8]=str(newtorprm)
        torline=' '.join(torlinesplit)+'\n'
        newtorsionprms.append(torline)
    return newtorsionprms


def FindMissingParameters(poltype,indicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs):
    missing=[]
    for indices,smartsatomorders in indicestosmartsatomorders.items():
        nindexes=[]
        for idx in indices: 
            neighborindexes=indextoneighbidxs[idx]
            for nidx in neighborindexes:
                if nidx not in nindexes:
                    nindexes.append(nidx)
        smarts=smartsatomorders[0]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)
        matcharray=[]
        for match in matches:
            for idx in match:
                if idx not in matcharray:
                    matcharray.append(idx)
            
        check=CheckIfNeighborsExistInSMARTMatch(poltype,nindexes,matcharray)
               
        if check==False or '*' in smarts or '~' in smarts:
            missing.append(indices)
    return missing

def ReduceMissingVdwByTypes(poltype,vdwmissing):
    reducedvdwmissing=[]
    typesfound=[]
    for vdwarray in vdwmissing:
        vdwatomidx=vdwarray[0]
        vdwbabelidx=vdwatomidx+1
        vdwtype=poltype.idxtosymclass[vdwbabelidx]
        if vdwtype in typesfound:
            continue
        else:
            typesfound.append(vdwtype)
            reducedvdwmissing.append(vdwbabelidx)
    return reducedvdwmissing



def ConvertToBabelList(poltype,listforprm):
    babellist=[]
    for ls in listforprm:
        babells=[i+1 for i in ls]
        babellist.append(babells)
    return babellist



def AddReverseKeys(poltype,tinkerclassestopoltypeclasses):
    newtinkerclassestopoltypeclasses={}
    for tinkerclasses,poltypeclasses in tinkerclassestopoltypeclasses.items():
        revtinkerclasses=tinkerclasses[::-1]
        revpoltypeclasses=[]
        for ls in poltypeclasses:
            revls=ls[::-1]
            revpoltypeclasses.append(revls) 
        newtinkerclassestopoltypeclasses[tinkerclasses]=poltypeclasses
        newtinkerclassestopoltypeclasses[revtinkerclasses]=revpoltypeclasses
    return newtinkerclassestopoltypeclasses


def TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools):
    opbendtinkerclassestotrigonalcenterbools={}
    for bondindices,tinkerclasses in opbendbondindicestotinkerclasses.items():
        boolarray=opbendbondindicestotrigonalcenterbools[bondindices]
        opbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]=boolarray

    return opbendtinkerclassestotrigonalcenterbools



def FilterDictionaries(poltype,dics,ls):
    newdics=[]
    for dic in dics:
        newdic={}
        for key,value in dic.items():
            if key in ls:
                newdic[key]=value
        newdics.append(newdic)
    return newdics

def ConvertListOfListToListOfTuples(poltype,listoflist):
    listoftuples=[]
    for item in listoflist:
        tup=tuple(item)
        listoftuples.append(tup)
    return listoftuples


def GrabSmallMoleculeAMOEBAParameters(poltype,optmol,mol,rdkitmol):
    bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters=ReadExternalDatabase(poltype)
    smartsatomordertoelementtinkerdescrip=ReadSmallMoleculeLib(poltype,poltype.smallmoleculesmartstotinkerdescrip)
    elementtinkerdescriptotinkertype,tinkertypetoclass=GrabTypeAndClassNumbers(poltype,poltype.smallmoleculeprmlib)
    listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForParameters(poltype,mol)
    planarbonds=GrabPlanarBonds(poltype,listofbondsforprm,mol)
    listofanglesthatneedplanarkeyword=CheckForPlanerAngles(poltype,listofanglesforprm,mol)
    parametersmartslist=GrabSMARTSList(poltype,smartsatomordertoelementtinkerdescrip)
    parametersmartstomaxcommonsubstructuresmarts,maxatomsize=FindMaximumCommonSubstructures(poltype,parametersmartslist,rdkitmol)
    if len(list(parametersmartstomaxcommonsubstructuresmarts.keys()))!=0:
         parametersmartslist=list(parametersmartstomaxcommonsubstructuresmarts.keys())
    atomindicesforprmtosmartslist=GenerateSMARTSListFromAtomList(poltype,listofatomsforprm,rdkitmol,mol,maxatomsize)
    bondsforprmtosmartslist=GenerateSMARTSListFromAtomList(poltype,listofbondsforprm,rdkitmol,mol,maxatomsize)
    planarbondsforprmtosmartslist=GenerateSMARTSListFromAtomList(poltype,listofbondsforprm,rdkitmol,mol,maxatomsize)
    anglesforprmtosmartslist=GenerateSMARTSListFromAtomList(poltype,listofanglesforprm,rdkitmol,mol,maxatomsize)
    planaranglesforprmtosmartslist=GenerateSMARTSListFromAtomList(poltype,listofanglesthatneedplanarkeyword,rdkitmol,mol,maxatomsize)
    torsionsforprmtosmartslist=GenerateSMARTSListFromAtomList(poltype,listoftorsionsforprm,rdkitmol,mol,maxatomsize)
    bondindicestoextsmartsmatchlength,bondindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,bondsmartsatomordertoparameters)
    angleindicestoextsmartsmatchlength,angleindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,anglesmartsatomordertoparameters)
    strbndindicestoextsmartsmatchlength,strbndindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,strbndsmartsatomordertoparameters)
    torsionindicestoextsmartsmatchlength,torsionindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,torsionsmartsatomordertoparameters)

    opbendindicestoextsmartsmatchlength,opbendindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,opbendsmartsatomordertoparameters)
    vdwindicestoextsmartsmatchlength,vdwindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,vdwsmartsatomordertoparameters)

    atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,atomindicesforprmtosmartslist,parametersmartslist,mol)
    bondsforprmtoparametersmarts,bondsforprmtosmarts=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,bondsforprmtosmartslist,parametersmartslist,mol)
    planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,planarbondsforprmtosmartslist,parametersmartslist,mol)
    anglesforprmtoparametersmarts,anglesforprmtosmarts=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,anglesforprmtosmartslist,parametersmartslist,mol)
    planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,planaranglesforprmtosmartslist,parametersmartslist,mol)
    torsionsforprmtoparametersmarts,torsionsforprmtosmarts=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,torsionsforprmtosmartslist,parametersmartslist,mol)
    atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,vdwindicestoextsmarts=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,vdwindicestoextsmartsmatchlength,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,vdwindicestoextsmarts)
    bondsforprmtoparametersmarts,bondsforprmtosmarts,bondindicestoextsmarts=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,bondindicestoextsmartsmatchlength,bondsforprmtoparametersmarts,bondsforprmtosmarts,bondindicestoextsmarts)

    planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,bondindicestoextsmarts=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,bondindicestoextsmartsmatchlength,planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,bondindicestoextsmarts)

    anglesforprmtoparametersmarts,anglesforprmtosmarts,angleindicestoextsmarts=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,angleindicestoextsmartsmatchlength,anglesforprmtoparametersmarts,anglesforprmtosmarts,angleindicestoextsmarts)

    planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,angleindicestoextsmarts=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,angleindicestoextsmartsmatchlength,planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,angleindicestoextsmarts)

    torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionindicestoextsmarts=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,torsionindicestoextsmartsmatchlength,torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionindicestoextsmarts)

    atomindextotinkertype,atomindextotinkerclass,atomindextoparametersmartsatomorder,atomindextoelementtinkerdescrip,atomindextosmartsatomorder=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

    bondindicestotinkertypes,bondindicestotinkerclasses,bondindicestoparametersmartsatomorders,bondindicestoelementtinkerdescrips,bondindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,bondsforprmtoparametersmarts,bondsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

    opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders=FilterBondSMARTSEnviorment(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses)
    newplanarbonds=ConvertListOfListToListOfTuples(poltype,planarbonds)
    newdics=FilterDictionaries(poltype,[opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders],newplanarbonds)
    opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders=newdics[:]
    planarbondindicestotinkertypes,planarbondindicestotinkerclasses,planarbondindicestoparametersmartsatomorders,planarbondindicestoelementtinkerdescrips,planarbondindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

    angleindicestotinkertypes,angleindicestotinkerclasses,angleindicestoparametersmartsatomorders,angleindicestoelementtinkerdescrips,angleindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,anglesforprmtoparametersmarts,anglesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

    planarangleindicestotinkertypes,planarangleindicestotinkerclasses,planarangleindicestoparametersmartsatomorders,planarangleindicestoelementtinkerdescrips,planarangleindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

    torsionindicestotinkertypes,torsionindicestotinkerclasses,torsionindicestoparametersmartsatomorders,torsionindicestoelementtinkerdescrips,torsionindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,torsionsforprmtoparametersmarts,torsionsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

    indextoneighbidxs=FindAllNeighborIndexes(poltype,rdkitmol)
    torsionsmissing=FindMissingTorsions(poltype,torsionindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
    vdwmissing=FindMissingParameters(poltype,atomindextosmartsatomorder,rdkitmol,mol,indextoneighbidxs)
    missingvdwatomindices=ReduceMissingVdwByTypes(poltype,vdwmissing)
    bondmissing=FindMissingParameters(poltype,bondindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
    anglemissing=FindMissingParameters(poltype,angleindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
    torsionsmissingindicestotinkerclasses=PruneDictionary(poltype,torsionsmissing,torsionindicestotinkerclasses)
    atomtinkerclasstopoltypeclass=TinkerClassesToPoltypeClasses(poltype,atomindextotinkerclass)
    bondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondindicestotinkerclasses)
   
    planarbondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarbondindicestotinkerclasses)
    opbendtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,opbendbondindicestotinkerclasses)
    opbendtinkerclassestopoltypeclasses=AddReverseKeys(poltype,opbendtinkerclassestopoltypeclasses)
    opbendbondindicestotrigonalcenterbools=CheckTrigonalCenters(poltype,listofbondsforprm,mol)
    opbendtinkerclassestotrigonalcenterbools=TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools)

    angletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,angleindicestotinkerclasses)
    planarangletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarangleindicestotinkerclasses)
    torsiontinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionindicestotinkerclasses)
    torsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionsmissingindicestotinkerclasses)
    torsionpoltypeclassestoparametersmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,torsionindicestoparametersmartsatomorders,torsionindicestotinkerclasses,torsiontinkerclassestopoltypeclasses)
    torsionpoltypeclassestosmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,torsionindicestosmartsatomorders,torsionindicestotinkerclasses,torsiontinkerclassestopoltypeclasses)
    torsionpoltypeclassestoelementtinkerdescrips=ConvertIndicesDictionaryToPoltypeClasses(poltype,torsionindicestoelementtinkerdescrips,torsionindicestotinkerclasses,torsiontinkerclassestopoltypeclasses)
    anglepoltypeclassestoparametersmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,angleindicestoparametersmartsatomorders,angleindicestotinkerclasses,angletinkerclassestopoltypeclasses)

    anglepoltypeclassestosmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,angleindicestosmartsatomorders,angleindicestotinkerclasses,angletinkerclassestopoltypeclasses)

    anglepoltypeclassestoelementtinkerdescrips=ConvertIndicesDictionaryToPoltypeClasses(poltype,angleindicestoelementtinkerdescrips,angleindicestotinkerclasses,angletinkerclassestopoltypeclasses)

    bondpoltypeclassestoparametersmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,bondindicestoparametersmartsatomorders,bondindicestotinkerclasses,bondtinkerclassestopoltypeclasses)

    bondpoltypeclassestosmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses,bondtinkerclassestopoltypeclasses)

    bondpoltypeclassestoelementtinkerdescrips=ConvertIndicesDictionaryToPoltypeClasses(poltype,bondindicestoelementtinkerdescrips,bondindicestotinkerclasses,bondtinkerclassestopoltypeclasses)

    atompoltypeclasstoparametersmartsatomorder=ConvertIndicesDictionaryToPoltypeClasses(poltype,atomindextoparametersmartsatomorder,atomindextotinkerclass,atomtinkerclasstopoltypeclass)

    atompoltypeclasstosmartsatomorder=ConvertIndicesDictionaryToPoltypeClasses(poltype,atomindextosmartsatomorder,atomindextotinkerclass,atomtinkerclasstopoltypeclass)

    atompoltypeclassestoelementtinkerdescrip=ConvertIndicesDictionaryToPoltypeClasses(poltype,atomindextoelementtinkerdescrip,atomindextotinkerclass,atomtinkerclasstopoltypeclass)
    poltypetoprmtype={} # dont need polarize parameters 
    typestoframedefforprmfile={} # dont need multipole parameters
    fname=poltype.smallmoleculeprmlib
    bondprms,angleprms,torsionprms,strbndprms,mpoleprms,opbendprms,polarizeprms,vdwprms,torsiontopitor=GrabParametersFromPrmFile(poltype,bondtinkerclassestopoltypeclasses,opbendtinkerclassestopoltypeclasses,opbendtinkerclassestotrigonalcenterbools,angletinkerclassestopoltypeclasses,torsiontinkerclassestopoltypeclasses,poltypetoprmtype,atomtinkerclasstopoltypeclass,typestoframedefforprmfile,fname,True)
    torsionprms=CorrectPitorEnergy(poltype,torsionprms,torsiontopitor)
    bondprms,bondpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,bondprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'bond')
    angleprms,anglepoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,angleprms,angleindicestoextsmarts,anglesmartsatomordertoparameters,'angle')
    
    strbndprms,strbndpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,strbndprms,angleindicestoextsmarts,strbndsmartsatomordertoparameters,'strbnd')
    torsionprms,torsionpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,torsionprms,torsionindicestoextsmarts,torsionsmartsatomordertoparameters,'torsion')
    opbendprms,opbendpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,opbendprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'opbend')
    vdwprms,vdwpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,vdwprms,vdwindicestoextsmarts,vdwsmartsatomordertoparameters,'vdw')
    angleprms=ModifyAngleKeywords(poltype,angleprms,planarangletinkerclassestopoltypeclasses)
    bondlistbabel=ConvertToBabelList(poltype,listofbondsforprm)
    anglelistbabel=ConvertToBabelList(poltype,listofanglesforprm)
    bondprms=AddOptimizedBondLengths(poltype,optmol,bondprms,bondlistbabel)
    angleprms=AddOptimizedAngleLengths(poltype,optmol,angleprms,anglelistbabel)
    torsionkeystringtoparameters=GrabTorsionParameterCoefficients(poltype,torsionprms)
    torsionprms=ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms)
    potentialmissingopbendprmtypes=FindPotentialMissingParameterTypes(poltype,opbendprms,planarbondtinkerclassestopoltypeclasses)
    potentialmissingopbendprmindices=ConvertPoltypeClassesToIndices(poltype,potentialmissingopbendprmtypes)
    potentialmissingopbendprmindices=FilterIndices(poltype,potentialmissingopbendprmindices,planarbonds)
    missingopbendprmindices=CheckIfParametersExist(poltype,potentialmissingopbendprmindices,opbendprms)
    defaultvalues=None
    if len(missingopbendprmindices)!=0:
        newopbendprms,defaultvalues=DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools)
        opbendprms.extend(newopbendprms)
    vdwprmstotransferinfo=MapParameterLineToTransferInfo(poltype,vdwprms,atompoltypeclasstoparametersmartsatomorder,atompoltypeclasstosmartsatomorder,atompoltypeclassestoelementtinkerdescrip,vdwpoltypeclassestosmartsatomordersext)
    bondprmstotransferinfo=MapParameterLineToTransferInfo(poltype,bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips,bondpoltypeclassestosmartsatomordersext)
    opbendprmstotransferinfo=MapParameterLineToTransferInfo(poltype,opbendprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips,opbendpoltypeclassestosmartsatomordersext,defaultvalues)
    angleprmstotransferinfo=MapParameterLineToTransferInfo(poltype,angleprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,anglepoltypeclassestosmartsatomordersext)
    strbndprmstotransferinfo=MapParameterLineToTransferInfo(poltype,strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,strbndpoltypeclassestosmartsatomordersext)
    torsionprmstotransferinfo=MapParameterLineToTransferInfo(poltype,torsionprms,torsionpoltypeclassestoparametersmartsatomorders,torsionpoltypeclassestosmartsatomorders,torsionpoltypeclassestoelementtinkerdescrips,torsionpoltypeclassestosmartsatomordersext)
    torsionsmissing=ConvertToPoltypeClasses(poltype,torsionsmissing)
    WriteOutList(poltype,torsionsmissing,poltype.torsionsmissingfilename)
    WriteDictionaryToFile(poltype,torsionkeystringtoparameters,poltype.torsionprmguessfilename)
    WriteOutList(poltype,missingvdwatomindices,poltype.vdwmissingfilename)

    return bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,torsionsmissing,torsionkeystringtoparameters,missingvdwatomindices
