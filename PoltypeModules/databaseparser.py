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
from rdkit import DataStructs
import torsiongenerator as torgen
from itertools import combinations
  
def appendtofile(poltype, vf,newname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo):
    temp=open(vf,'r')
    results=temp.readlines()
    temp.close()
    foundatomblock=False
    atomline=False
    wroteout=False
    tempname=vf.replace('.key','_temp.key')
    f=open(tempname,'w')
    linestoskip=[]
    for line in results:
        if 'atom' in line:
            atomline=True
            if foundatomblock==False:
                foundatomblock=True
        elif 'polarize' in line:
            atomline=False
            linesplit=line.split()
            atomtype=linesplit[1]
            resid=linesplit[3:]
            residstring=' '.join(resid)
            for polarprmsline,transferinfo in polarprmstotransferinfo.items(): 
                polarlinesplit=polarprmsline.split()
                polartype=polarlinesplit[1]
                if atomtype==polartype:
                    f.write(transferinfo)
                    f.write(polarprmsline.replace('\n','')+' '+residstring+'\n')
                    linestoskip.append(line)
                    break
                    
        else:
            atomline=False
        if foundatomblock==True and atomline==False and wroteout==False:
            wroteout=True
            f.write('\n')
            if poltype.forcefield=='AMOEBA+':
                for line,transferinfo in amoebaplusvdwprmstotransferinfo.items():
                    f.write(transferinfo)
                    f.write(line)

            else:
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
            for line in soluteprms:
                f.write(line)
            f.write('\n')
            for line,transferinfo in tortorprmstotransferinfo.items():
                if 'tortors' in line:
                    f.write(transferinfo)
                f.write(line)
            f.write('\n')

            if poltype.forcefield=='AMOEBA+':
                for line,transferinfo in ctprmstotransferinfo.items():
                    f.write(transferinfo)
                    f.write(line)
                f.write('\n')
                for line,transferinfo in cpprmstotransferinfo.items():
                    f.write(transferinfo)
                    f.write(line)
                f.write('\n')
                for line,transferinfo in bondcfprmstotransferinfo.items():
                    f.write(transferinfo)
                    f.write(line)
                f.write('\n')
                for line,transferinfo in anglecfprmstotransferinfo.items():
                    f.write(transferinfo)
                    f.write(line)
                f.write('\n')


                                    
        else:
            if line not in linestoskip:
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
                if RepresentsInt(element):
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
            rev=False
            prm1=linesplitall[8]
            prm2=linesplitall[10]
            if tuple(angleclasslist) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist)
                foundstrbnd=True
            elif tuple(angleclasslist[::-1]) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundstrbnd=True
                rev=True
            if foundstrbnd==True:
                classes=angletinkerclassestopoltypeclasses[angletup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    if rev==True:
                        linesplitall[8]=prm2
                        linesplitall[10]=prm1
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
                        if nextnextatomidx!=nidx and nextnextatomidx!=atomidx:                
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
    brackindex=tempsplit[0].index(']')
    extra=''
    for i in range(len(tempsplit[0][brackindex:])):
        char=tempsplit[0][brackindex+i]
        extra+=char
    if len(extra)>1:
        extra=extra[1:]
    else:
        extra=''
    tempsplit[0]='[*]'
    tempsplit[0]+=extra
    newsmarts='~'.join(tempsplit)
    return newsmarts

def ReplaceRightAtomIdentity(poltype,smarts):
    smartsplit=smarts.split('~') 
    tempsplit=smartsplit[:]
    if ']' in tempsplit[-1]:
        idx=-1
    else:
        idx=-2
    brackindex=tempsplit[idx].index(']')
    extra=''
    for i in range(len(tempsplit[idx][brackindex:])):
        char=tempsplit[idx][brackindex+i]
        extra+=char
    if len(extra)>1:
        extra=extra[1:]
    else:
        extra=''

    tempsplit[idx]='[*]'

    tempsplit[idx]+=extra
    newsmarts='~'.join(tempsplit)
    return newsmarts

def ReplaceAtomIdentitiesOnEnd(poltype,smartsls):
    newsmartslist=[]
    for smarts in smartsls:
        smartsplit=smarts.split('~') 
        tempsplit=copy.deepcopy(smartsplit)
        newsmarts=ReplaceLeftAtomIdentity(poltype,smarts)
        newsmarts=ReplaceRightAtomIdentity(poltype,newsmarts)

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
    prmsmartstomolnum={}
    for prmsmarts in parametersmartslist:
        structure = Chem.MolFromSmarts(prmsmarts)
        structurenumatoms=structure.GetNumAtoms()
        prmsmartstomolnum[prmsmarts]=structurenumatoms
    minfragmentsize=min(prmsmartstomolnum.values())
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
        sizediff=structurenumatoms-minfragmentsize # some fragments have greater absolute size so without adding this than we are biasing smaller fragments that match to be a "closer match"
        diff=np.abs(atomnum-structurenumatoms)-sizediff
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
            tmpsmarts=''
            for eidx in range(len(newsmarts)):
                e=newsmarts[eidx]
                replace=False
                if eidx>0:
                    if e=='#':
                        if newsmarts[eidx-1]=='[':
                            pass
                        else:
                            replace=True
                if replace==True:
                    tmpsmarts+='~'
                else:
                    tmpsmarts+=e
            newsmarts=tmpsmarts
            temp.append(newsmarts)
        newsmartslist.append(temp)
    return newsmartslist
         

def CountRings(poltype,smartsfortransfer):
    counts=0
    for idx in range(len(smartsfortransfer)):
        e=smartsfortransfer[idx]
        if e=='1':
            preve=smartsfortransfer[idx-1]
            if preve!='#':
                counts+=1
    return counts



def CheckMatch(poltype,match,atomindices,smartsfortransfer):
    allin=True
    for idx in atomindices:
        if idx not in match:
            allin=False

    validmatch=True
    if allin==True:
        branchcounts=smartsfortransfer.count('(')
        ringcounts=CountRings(poltype,smartsfortransfer)
        if branchcounts==0 and len(atomindices)>1 and ringcounts==0:
            indices=[]
            for idx in atomindices:
                matchidx=match.index(idx)
                indices.append(matchidx)
            checkconsecleft=checkConsecutive(indices)
            revindices=indices[::-1]
            checkconsecright=checkConsecutive(revindices)
            if checkconsecleft==False and checkconsecright==False:
                validmatch=False
    else:
        validmatch=False

    return validmatch 

def checkConsecutive(l): 
    consec=True
    previdx=None
    for i in range(len(l)):
        idx=l[i]
        if previdx!=None:
            if idx<previdx:
                consec=False
        previdx=idx
    return consec
        



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
        matches=list(rdkitmol.GetSubstructMatches(substructure,maxMatches=1000))
        sp=openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,smarts)
        diditmatch=sp.Match(poltype.mol)
        babelmatches=sp.GetMapList()
        newbabelmatches=[]
        for mtch in babelmatches:
            newmtch=[i-1 for i in mtch]
            newbabelmatches.append(newmtch)
        for newmtch in newbabelmatches:
            if newmtch not in matches:
                matches.append(newmtch)

        for match in matches:
            validmatch=CheckMatch(poltype,match,atomindices,smartsfortransfer)
            if validmatch==True:
               indices=list(range(len(match)))
               smartsindextomoleculeindex=dict(zip(indices,match)) 
               moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}
            
        structure = Chem.MolFromSmarts(parametersmarts)
        fragmentfilepath='fragment.mol'
        if os.path.isfile(fragmentfilepath):
            os.remove(fragmentfilepath)
        rdmolfiles.MolToMolFile(structure,fragmentfilepath)
        obConversion = openbabel.OBConversion()
        fragbabelmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(fragmentfilepath)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(fragbabelmol, fragmentfilepath)
        fragidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,fragbabelmol)

        substructure = Chem.MolFromSmarts(smartsfortransfer)
        matches=structure.GetSubstructMatches(substructure)
        firstmatch=matches[0]

        indices=list(range(len(firstmatch)))
        smartsindextoparametersmartsindex=dict(zip(indices,firstmatch)) 
        parametersmartsordertoelementtinkerdescrip={}
        for parametersmartsatomorderlist,elementtinkerdescrip in smartsatomordertoelementtinkerdescrip.items():
            prmsmarts=parametersmartsatomorderlist[0]
            atomorderlist=parametersmartsatomorderlist[1]
            if prmsmarts==parametersmarts:
                atomorderlist=parametersmartsatomorderlist[1]

                for atomorder in atomorderlist:
                    
                    parametersmartsordertoelementtinkerdescrip[atomorder]=elementtinkerdescrip 

        for fragidx,symclass in fragidxtosymclass.items():
            indexes=GrabKeysFromValue(poltype,fragidxtosymclass,symclass)
            specialindex=None
            for index in indexes:
                if index in parametersmartsordertoelementtinkerdescrip.keys():
                    specialindex=index
            if specialindex!=None:
                elementtinkerdescrip=parametersmartsordertoelementtinkerdescrip[specialindex]
                for index in indexes:
                    parametersmartsordertoelementtinkerdescrip[index]=elementtinkerdescrip
        smartindices=[moleculeindextosmartsindex[i] for i in atomindices]
        parametersmartindices=[smartsindextoparametersmartsindex[i] for i in smartindices]
        parametersmartsorders=[i+1 for i in parametersmartindices]
        elementtinkerdescrips=[parametersmartsordertoelementtinkerdescrip[i] for i in parametersmartsorders]
        tinkertypes=[elementtinkerdescriptotinkertype[i] for i in elementtinkerdescrips]
        tinkerclasses=[tinkertypetoclass[i] for i in tinkertypes]
        tinkerclasses=[int(i) for i in tinkerclasses]
        tinkertypes=[int(i) for i in tinkertypes]
        smartsorders=[i+1 for i in smartindices]
        atomindicestotinkertypes[atomindices]=tinkertypes
        atomindicestotinkerclasses[atomindices]=tinkerclasses
        atomindicestoparametersmartsatomorders[atomindices]=[parametersmarts,parametersmartsorders]
        atomindicestoelementtinkerdescrips[atomindices]=elementtinkerdescrips             
        atomindicestosmartsatomorders[atomindices]=[smarts,smartsorders]


    return atomindicestotinkertypes,atomindicestotinkerclasses,atomindicestoparametersmartsatomorders,atomindicestoelementtinkerdescrips,atomindicestosmartsatomorders


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
        newline=line.replace('anglef','angle')

        if found==False:
            newline=newline.replace('anglep','angle')

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
        if len(allbonds)==0:
            pass 
        else:
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
        if len(allangles)==0:
            pass
        else:
            avganglelength=round(tot/len(allangles),2)
            linesplit=re.split(r'(\s+)', line)
            linesplit=linesplit[:11]
            linesplit.append('\n')
            linesplit[10]=str(avganglelength)
            line=''.join(linesplit)
        line+='\n'
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

def CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol):
    allhydrogentoroneside=True
    atomobjects=[mol.GetAtom(i) for i in babelindices]
    a,b,c,d=atomobjects[:]
    aidx,bidx,cidx,didx=babelindices[:]
    aatomicnum=a.GetAtomicNum()
    datomicnum=d.GetAtomicNum()
    if aatomicnum!=1 and datomicnum!=1:
        allhydrogentoroneside=False
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
            if aatomnum!=1 and datomnum!=1:
                allhydrogentoroneside=False
    return allhydrogentoroneside


def FindAllConsecutiveRotatableBonds(poltype,mol,listofbondsforprm):
    totalbondscollector=[]
    newrotbnds=[]
    for rotbnd in listofbondsforprm:
        b=rotbnd[0]+1
        c=rotbnd[1]+1
        batom=poltype.mol.GetAtom(b)
        catom=poltype.mol.GetAtom(c)
        bval=batom.GetValence()
        cval=catom.GetValence()
        if bval<2 or cval<2:
            continue
        newrotbnds.append(rotbnd)
    combs=list(combinations(newrotbnds,2)) 
    for comb in combs:
        firstbnd=comb[0]
        secondbnd=comb[1]
        total=firstbnd[:]+secondbnd[:]
        totalset=set(total)
        if len(totalset)==3:
            if firstbnd[-1]!=secondbnd[-1] and firstbnd[0]!=secondbnd[-1] and firstbnd[0]!=secondbnd[0]:
                secondbnd=secondbnd[::-1]
            elif firstbnd[0]==secondbnd[-1] and firstbnd[-1]!=secondbnd[-1] and firstbnd[0]!=secondbnd[0]:
                firstbnd=firstbnd[::-1]
            elif firstbnd[0]!=secondbnd[-1] and firstbnd[-1]!=secondbnd[-1] and firstbnd[0]==secondbnd[0]:
                firstbnd=firstbnd[::-1]
                secondbnd=secondbnd[::-1]


            catomindex=firstbnd[1]+1
            catom=poltype.mol.GetAtom(catomindex)
            isinring=catom.IsInRing()
            if isinring==True:
                continue
           
            totalbondscollector.append([firstbnd,secondbnd])
    return totalbondscollector

def FindAdjacentMissingTorsionsForTorTor(poltype,torsionsmissing,totalbondscollector,tortorsmissing):
    for bndlist in totalbondscollector:
        first=bndlist[0]
        second=bndlist[1]
        foundfirst=CheckIfRotatableBondInMissingTorsions(poltype,first,torsionsmissing) 
        foundsecond=CheckIfRotatableBondInMissingTorsions(poltype,second,torsionsmissing) 
        foundtortormissing=CheckIfRotatableBondInMissingTorTors(poltype,first,second,tortorsmissing)
        if foundtortormissing==False:
            continue
        if (foundfirst==False and foundsecond==True and poltype.tortor==True):
            b=first[0]+1
            c=first[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tor)
        elif (foundfirst==True and foundsecond==False and poltype.tortor==True):
            b=second[0]+1
            c=second[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tor)
        elif (foundfirst==False and foundsecond==False and poltype.tortor==True):
            b=first[0]+1
            c=first[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tor)
            b=second[0]+1
            c=second[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tor)

    return torsionsmissing

def CheckIfRotatableBondInMissingTorsions(poltype,rotbnd,torsionsmissing):
    found=False
    for tor in torsionsmissing:
        bnd=[tor[1],tor[2]]
        if bnd==rotbnd or bnd[::-1]==rotbnd:
            found=True
            break
    return found

def CheckIfRotatableBondInMissingTorTors(poltype,firstbnd,secondbnd,tortormissing):
    foundtortormissing=False
    for tortor in tortormissing:
        first=[tortor[1],tortor[2]]
        second=[tortor[2],tortor[3]]
        if (first==firstbnd or first[::-1]==firstbnd) and (second==secondbnd or second[::-1]==secondbnd):
            foundtortormissing=True
        elif (first==secondbnd or first[::-1]==secondbnd) and (second==firstbnd or second[::-1]==firstbnd):
            foundtortormissing=True

    return foundtortormissing

def FindMissingTorTors(poltype,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,rdkitmol,mol,indextoneighbidxs,totalbondscollector):
    tortorsmissing=[]
    tortorsfound=[]
    for tortorindices,extsmarts in tortorindicestoextsmarts.items(): #only for ones that have match
        for tortorsmartsatomorder,parameters in tortorsmartsatomordertoparameters.items():
            smarts=tortorsmartsatomorder[0]
            if smarts==extsmarts:
                aidx,bidx,cidx,didx,eidx=tortorindices[:]
                firstneighborindexes=indextoneighbidxs[aidx]
                secondneighborindexes=indextoneighbidxs[bidx]
                thirdneighborindexes=indextoneighbidxs[cidx]
                fourthneighborindexes=indextoneighbidxs[didx]
                fifthneighborindexes=indextoneighbidxs[eidx]
                neighborindexes=firstneighborindexes+secondneighborindexes+thirdneighborindexes+fourthneighborindexes+fifthneighborindexes
                substructure = Chem.MolFromSmarts(smarts)
                matches=rdkitmol.GetSubstructMatches(substructure)
                matcharray=[]
                for match in matches:
                    for idx in match:
                        if idx not in matcharray:
                            matcharray.append(idx)
                check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
                if check==False:
                    tortorsmissing.append(tortorindices)
                else:
                    tortorsfound.append(tortorindices)
    for bndlist in totalbondscollector:
        first=bndlist[0]
        second=bndlist[1]
        b,c=first[:]
        d=second[0]
        aatom,dnewatom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b+1),poltype.mol.GetAtom(c+1))
        bnewatom,eatom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(c+1),poltype.mol.GetAtom(d+1))
        a=aatom.GetIdx()
        dnew=dnewatom.GetIdx()
        bnew=bnewatom.GetIdx()
        e=eatom.GetIdx()
        indices=[a-1,b,c,d,e-1]
        inonlyrotbnds=CheckIfRotatableBondsInOnlyRotBnds(poltype,first,second)
        foundtortormissing=CheckIfRotatableBondInMissingTorTors(poltype,[bnew-1,c],[c,dnew-1],tortorsmissing)
        foundtortor=CheckIfRotatableBondInMissingTorTors(poltype,[bnew-1,c],[c,dnew-1],tortorsfound)
        if foundtortormissing==False and foundtortor==False and inonlyrotbnds==True:
            tortorsmissing.append(indices)

    return tortorsmissing

def CheckIfRotatableBondsInOnlyRotBnds(poltype,first,second):
    inonlyrotbnds=True
    first=[i+1 for i in first]
    second=[i+1 for i in second]
    if (first in poltype.onlyrotbndslist or first[::-1] in poltype.onlyrotbndslist) and (second in poltype.onlyrotbndslist or second[::-1] in poltype.onlyrotbndslist):
        pass
    else:
        inonlyrotbnds=False

    return inonlyrotbnds 


def FindMissingTorsions(poltype,torsionindicestoparametersmartsenv,rdkitmol,mol,indextoneighbidxs):
    torsionsmissing=[]
    poormatchingaromatictorsions=[]
    for torsionindices,smartsenv in torsionindicestoparametersmartsenv.items():
        aidx,bidx,cidx,didx=torsionindices[:]
        babelindices=[i+1 for i in torsionindices]
        allhydrogentor=CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol)
        allhydrogentoroneside=CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol)
        atoma=rdkitmol.GetAtomWithIdx(aidx)
        atomd=rdkitmol.GetAtomWithIdx(didx)
        atomicnumatoma=atoma.GetAtomicNum()
        atomicnumatomd=atomd.GetAtomicNum()
        abidx,bbidx,cbidx,dbidx=babelindices[:]
        if len(poltype.onlyrotbndslist)!=0:
            if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                if atomicnumatoma==1 or atomicnumatomd==1:
                    if allhydrogentor==False and allhydrogentoroneside==False:
                       continue
                else:

                    if torsionindices not in torsionsmissing:
                        torsionsmissing.append(torsionindices)
                        continue
            else:
                continue
        bond=mol.GetBond(bbidx,cbidx)
        bondorder=bond.GetBondOrder()
        if bondorder!=1: # then dont zero out
            continue 
        
        babelatoms=[mol.GetAtom(i) for i in babelindices]
        atomvals=[a.GetValence() for a in babelatoms]
        atomnums=[a.GetAtomicNum() for a in babelatoms]
        batomnum=atomnums[1]
        catomnum=atomnums[2]
           
        ringbools=[a.IsInRing() for a in babelatoms]
        arobools=[a.IsAromatic() for a in babelatoms]
        contin=False
        bnd=[babelindices[1],babelindices[2]]
        if (bnd in poltype.partialdoublebonds or bnd[::-1] in poltype.partialdoublebonds) and poltype.rotalltors==False:
            continue 
        ringb=ringbools[1]
        ringc=ringbools[2]
        aroa=arobools[0]
        arob=arobools[1]
        aroc=arobools[2]
        arod=arobools[3]
        hybs=[a.GetHyb() for a in babelatoms]
        hybb=hybs[1]
        hybc=hybs[2]
        anyarot2=CheckForAnyAromaticsInRing(poltype,babelindices[1])
        anyarot3=CheckForAnyAromaticsInRing(poltype,babelindices[2])

        if contin==True:
            continue

        firstneighborindexes=indextoneighbidxs[aidx]
        secondneighborindexes=indextoneighbidxs[bidx]
        thirdneighborindexes=indextoneighbidxs[cidx]
        fourthneighborindexes=indextoneighbidxs[didx]
        neighborindexes=firstneighborindexes+secondneighborindexes+thirdneighborindexes+fourthneighborindexes
        smarts=smartsenv[0]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)

        matcharray=[]
        for match in matches:
            for idx in match:
                if idx not in matcharray:
                    matcharray.append(idx)
        check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
        
        if poltype.transferanyhydrogentor==True and (anyarot2==False and anyarot3==False): 
            if atomicnumatoma==1 or atomicnumatomd==1:
                if allhydrogentor==False and allhydrogentoroneside==False:
                    continue
                else:
                    if torsionindices not in torsionsmissing:
                        torsionsmissing.append(torsionindices)
        if '~' in smarts or '*' in smarts:
            if (anyarot2==False and anyarot3==False):
                if torsionindices not in torsionsmissing:
                    torsionsmissing.append(torsionindices)
            else:
                if torsionindices not in poormatchingaromatictorsions:
                    poormatchingaromatictorsions.append(torsionindices)
            continue

        if check==False:
            if (anyarot2==False and anyarot3==False):
                if torsionindices not in torsionsmissing:
                    torsionsmissing.append(torsionindices)
            else:
                if torsionindices not in poormatchingaromatictorsions:
                    poormatchingaromatictorsions.append(torsionindices)
        else:
            if poltype.rotalltors==True:
                if torsionindices not in torsionsmissing:
                    torsionsmissing.append(torsionindices)

    return torsionsmissing,poormatchingaromatictorsions 

def GrabRingAtoms(poltype,neighbatom):
    ri = poltype.rdkitmol.GetRingInfo()
    ls=ri.AtomRings()
    atmidx=neighbatom.GetIdx()
    for grp in ls:
        if atmidx in grp:
            return grp
    grp=[]
    return grp
    
def CheckForAnyAromaticsInRing(poltype,babelindex):
    rdkitindex=babelindex-1
    rdkitatom=poltype.rdkitmol.GetAtomWithIdx(rdkitindex)
    ringindexes=GrabRingAtoms(poltype,rdkitatom)
    anyaro=False
    for index in ringindexes:
        atm=poltype.rdkitmol.GetAtomWithIdx(index)
        if atm.GetIsAromatic()==True:
            anyaro=True
    return anyaro


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

def ZeroOutMissingStrbnd(poltype,anglemissingtinkerclassestopoltypeclasses,strbndprms):
    newstrbndprms=[]
    for line in strbndprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
        for sublist in anglemissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[8]='0'
                newlinesplit[10]='0'
                line=''.join(newlinesplit)
        newstrbndprms.append(line)

    return newstrbndprms

def AssignAngleGuessParameters(poltype,angletinkerclassestoexampleindices,anglemissingtinkerclassestopoltypeclasses,angleprms):
    newangleprms=[]
    for line in angleprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
        for tinkerclasses,sublist in anglemissingtinkerclassestopoltypeclasses.items():
            found=False
            if classes in sublist:
                found=True
            elif classes[::-1] in sublist: 
                found=True
            if found==True:
                if tinkerclasses in angletinkerclassestoexampleindices.keys(): 
                    exampleindices=angletinkerclassestoexampleindices[tinkerclasses]
                elif tinkerclasses[::-1] in angletinkerclassestoexampleindices.keys(): 
                    exampleindices=angletinkerclassestoexampleindices[tinkerclasses[::-1]]

                atoms=[poltype.rdkitmol.GetAtomWithIdx(k) for k in exampleindices]
                atomicnums=[a.GetAtomicNum() for a in atoms]
                atomicval=[a.GetExplicitValence() for a in atoms]

                angleguess=AngleGuess(poltype,atomicnums[0],atomicnums[1],atomicnums[2],atomicval[0],atomicval[1],atomicval[2])
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[8]=str(angleguess)
                line=''.join(newlinesplit)
        newangleprms.append(line)

    return newangleprms

def AssignBondGuessParameters(poltype,bondtinkerclassestoexampleindices,bondmissingtinkerclassestopoltypeclasses,bondprms):
    newbondprms=[]
    for line in bondprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2])])
        for tinkerclasses,sublist in bondmissingtinkerclassestopoltypeclasses.items():
            found=False
            if classes in sublist:
                found=True
            elif classes[::-1] in sublist: 
                found=True
            if found==True:
                if tinkerclasses in bondtinkerclassestoexampleindices.keys(): 
                    exampleindices=bondtinkerclassestoexampleindices[tinkerclasses]
                elif tinkerclasses[::-1] in bondtinkerclassestoexampleindices.keys(): 
                    exampleindices=bondtinkerclassestoexampleindices[tinkerclasses[::-1]]

                atoms=[poltype.rdkitmol.GetAtomWithIdx(k) for k in exampleindices]
                atomicnums=[a.GetAtomicNum() for a in atoms]
                atomicval=[a.GetExplicitValence() for a in atoms]

                bondguess=BondGuess(poltype,atomicnums[0],atomicnums[1],atomicval[0],atomicval[1])
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[6]=str(bondguess)
                line=''.join(newlinesplit)
        newbondprms.append(line)

    return newbondprms




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


def DefaultAromaticMissingTorsions(poltype,arotorsionsmissingtinkerclassestopoltypeclasses,torsionprms): # transfer bezene aromatic torsion paramerers from amoeba09
    newtorsionprms=[]
    for line in torsionprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])])
        for sublist in arotorsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                splitafter=newlinesplit[10:]
                splitafter[0]='-.670'
                splitafter[6]='10.684' # add Pitor energy back
                splitafter[12]='0'
                newlinesplit=newlinesplit[:10]+splitafter
                line=''.join(newlinesplit)
        newtorsionprms.append(line)

    return newtorsionprms


def GrabTorsionParameterCoefficients(poltype,torsionprms):
    torsionkeystringtoparameters={}
    for line in torsionprms:
        linesplit=line.split() 
        key = '%s %s %s %s' % (linesplit[1], linesplit[2], linesplit[3], linesplit[4])
        parameters=[float(linesplit[5]),float(linesplit[8]),float(linesplit[11])]
        allzero=True
        for prm in parameters:
            if prm!=0:
                allzero=False
        if allzero==False:
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
    poltypeclassesalreadyassigned=[]
    for indices,tinkerclasses in indicestotinkerclasses.items():
        babelindices=[i+1 for i in indices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        revpoltypeclasses=poltypeclasses[::-1]
        revtinkerclasses=tinkerclasses[::-1]
        if poltypeclasses in poltypeclassesalreadyassigned or revpoltypeclasses in poltypeclassesalreadyassigned:
            continue
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

        poltypeclassesalreadyassigned.append(poltypeclasses)                
        if tuple(tinkerclasses) not in tinkerclassestopoltypeclasses.keys():
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)]=[]
        if tuple(poltypeclasses) not in tinkerclassestopoltypeclasses[tuple(tinkerclasses)]: 
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)].append(tuple(poltypeclasses))

    return tinkerclassestopoltypeclasses


def ConvertIndicesDictionaryToPoltypeClasses(poltype,indicestovalue,indicestotinkerclasses,tinkerclassestopoltypeclasses):
    poltypeclassestovalue={}
    for indices,value in indicestovalue.items():
        if indices in indicestotinkerclasses.keys():
            tinkerclasses=tuple(indicestotinkerclasses[indices])
            if tinkerclasses in tinkerclassestopoltypeclasses.keys():
                poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses])
            elif tinkerclasses[::-1] in tinkerclassestopoltypeclasses.keys():
                poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses[::-1]])
            else:
                continue

            poltypeclassestovalue[poltypeclasses]=value
    return poltypeclassestovalue 

def CheckIfStringInStringList(poltype,string,stringlist):
    found=False
    for e in stringlist:
        if e==string:
            found=True
    return found

def GrabTypesFromPrmLine(poltype,ls):
    typenums=[]
    for e in ls:
        isint=RepresentsInt(e)  
        if isint==True:
            typenums.append(e) 
    if len(typenums)>=4:
        if typenums[-2]=='0' and typenums[-1]=='0':
            typenums=typenums[:2]
    if len(typenums)>4:
        typenums=typenums[:4]     
    return typenums


def SearchForPoltypeClasses(poltype,prmline,poltypeclasseslist):
    listofpoltypeclasses=None
    poltypeclasses=None
    allin=None
    prmlinesplit=prmline.split()
    typenums=GrabTypesFromPrmLine(poltype,prmlinesplit)
    for listofpoltypeclasses in poltypeclasseslist:
        for poltypeclasses in listofpoltypeclasses:
            allin=True
            if int(typenums[0])!=poltypeclasses[0]: # then try flipping
                poltypeclasses=poltypeclasses[::-1]
            for i in range(len(typenums)):
                poltypeclass=int(typenums[i])
                otherpoltypeclass=poltypeclasses[i]
                if poltypeclass!=otherpoltypeclass:
                    allin=False 
            if allin==True:
                return listofpoltypeclasses,poltypeclasses
    if allin!=None:
        if allin==False:
            listofpoltypeclasses=None
            poltypeclasses=None

    return listofpoltypeclasses,poltypeclasses


def MapParameterLineToTransferInfo(poltype,prms,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,poltypeclassestosmartsatomordersext,newpoltypeclassestocomments,newpoltypeclassestosmartslist,defaultvalues=None,keyword=None):
    prmstotransferinfo={}
    for line in prms:
        if keyword!=None:
            if keyword not in line:
                transferinfo='# blank'
            else:
                poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())
                smartsatomorders=poltypeclassestosmartsatomordersext[poltypeclasses]
                transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' from external database'+'\n'

        else:
            poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestoparametersmartsatomorders.keys())
            if poltypeclasses==None:
                     
                poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())
                if poltypeclasses==None:
                    poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,newpoltypeclassestocomments.keys())
                    comments=newpoltypeclassestocomments[poltypeclasses]
                    smartslist=newpoltypeclassestosmartslist[poltypeclasses]
                    commentstring=' '.join(comments)
                    smartsliststring=' '.join(smartslist)
                    transferinfoline='#'+' '+'updated valence parameter matched comments='+commentstring+' '+'SMARTS match = '+smartsliststring+'\n'
                else:
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
        poltypeclassesls,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,newpoltypeclasseslist)
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
    with open(filename, 'w') as filehandle:
        for listitem in ls:
            filehandle.write('%s\n' % listitem)

def ReadTorsionList(poltype,filename):
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                current = line[:-1]
                current=current.replace('[','').replace(']','')
                current=current.split(',')
                current=[int(i) for i in current]
                newls.append(current)

    return newls


def ReadTorTorList(poltype,filename):
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                current = line[:-1]
                current=current.replace('[','').replace(']','')
                current=current.split(',')
                current=[int(i) for i in current]
                newls.append(current)

    return newls



def ReadVdwList(poltype,filename):
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                current = line[:-1]
                newls.append(int(current))

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
    tortorsmartsatomordertoparameters={}
    tortorsmartsatomordertogrid={}

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
        elif keyword=='tortors':
            gridstring=linesplit[3].lstrip().rstrip()
            grid=gridstring.split()
            prmstring=linesplit[4]
            prmstringlist=prmstring.split(',')
            prmstringsplits=[s.lstrip().rstrip().split() for s in prmstringlist]
            tortorsmartsatomordertoparameters[smartsatomorder]=prmstringsplits
            tortorsmartsatomordertogrid[smartsatomorder]=grid
    return bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters,tortorsmartsatomordertoparameters,tortorsmartsatomordertogrid


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
            if '~' in smarts or '*' in smarts:
                poormatch=True
            else:
                poormatch=False
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms or poormatch==True:
                del indicesforprmtosmarts[indices]
                del indicesforprmtoparametersmarts[indices]
                newindicestoextsmarts[indices]=indicestoextsmarts[indices]
        elif indices[::-1] in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices[::-1]]
            smarts=smartsls[0]
            if '~' in smarts or '*' in smarts:
                poormatch=True
            else:
                poormatch=False
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms:
                del indicesforprmtosmarts[indices[::-1]]
                del indicesforprmtoparametersmarts[indices[::-1]]
                newindicestoextsmarts[indices[::-1]]=indicestoextsmarts[indices[::-1]]

    return indicesforprmtoparametersmarts,indicesforprmtosmarts,newindicestoextsmarts

def AddExternalDatabaseSMARTSMatchParameters(poltype,prms,indicestoextsmarts,smartsatomordertoparameters,keyword,smartsatomordertogrid=None):
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
                if keyword=='tortors':
                    grid=smartsatomordertogrid[smartsatomorder]
                    gridstring=' '.join(grid)
                    line+=gridstring
                    line+='\n'
                    prms.append(line)
                    for point in parameters:
                        string=' '.join(point)
                        string+='\n'
                        prms.append(string)
                    pass
                else:
               
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
    torsiontocount={}
    middletocount={}
     
    for torsion in torsiontopitor.keys():
        middle=tuple([torsion[1],torsion[2]])
        if middle not in middletocount.keys():
            middletocount[middle]=0
        middletocount[middle]+=1
    for torsion in torsiontopitor.keys():
        middle=tuple([torsion[1],torsion[2]])
        count=middletocount[middle]
        torsiontocount[torsion]=count
    for torline in torsionprms:
        torlinesplit=torline.split()
        tor=tuple([int(torlinesplit[1]),int(torlinesplit[2]),int(torlinesplit[3]),int(torlinesplit[4])])
        if tor in torsiontopitor:
            if tor in torsiontopitor.keys():
                pitorline=torsiontopitor[tor]
            else:
                pitorline=torsiontopitor[tor[::-1]]

            count=torsiontocount[tor]
            pitorlinesplit=pitorline.split()
            prm=float(pitorlinesplit[3])/count
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



def ConvertToBabelList(poltype,listforprm,indicestoclasses):
    babellist=[]
    for ls in listforprm:
        add=True
        if poltype.transferqmequilvalues==False:
            fwd=tuple(ls)
            rev=fwd[::-1]
            if fwd in indicestoclasses.keys() or rev in indicestoclasses.keys():
                add=False
        if add==True:
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
        revtinkerclasses=tinkerclasses[::-1]
        revboolarray=boolarray[::-1]
        opbendtinkerclassestotrigonalcenterbools[tuple(revtinkerclasses)]=revboolarray

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


def AddExternalDatabaseMatches(poltype, indicestosmartsatomorder,extindicestoextsmarts,smartsatomordertoparameters):
    newindicestosmartsatomorder=indicestosmartsatomorder.copy()
    for smartsatomorder in smartsatomordertoparameters.keys():
        smarts=smartsatomorder[0]
        for indices,extsmarts in extindicestoextsmarts.items():
            if smarts==extsmarts:
                newindicestosmartsatomorder[indices]=list(smartsatomorder)         
        
    return newindicestosmartsatomorder   
   
def AngleGuess(poltype,ita,itb,itc,iva,ivb,ivc):
    radian=57.29577951
    angunit = 1.0 / radian**2
    if (itb==6):
       if (ita==1):
          if (ivb==4):
             if (itc==1):
                angguess = 34.50
             elif (itc==6):
                angguess = 38.0
             elif (itc==7):
                angguess = 50.60
             elif (itc==8):
                angguess = 51.50
             elif (itc==9):
                angguess = 50.0
             else:
                angguess = 35.0
          elif (ivb==3):
             angguess = 32.00
          else:
             angguess = 32.00
       elif (ita==6):
          if (ivb==4):
             if (itc==6):
                angguess = 60.00
             elif (itc==7):
                angguess = 80.00
             elif (itc==8):
                angguess = 88.00
             elif (itc==9):
                angguess = 89.00
             elif (itc==14):
                angguess = 65.00
             elif (itc==15):
                angguess = 60.00
             elif (itc==16):
                angguess = 53.20
             elif(itc==17):
                angguess = 55.00
             else:
                angguess = 50.00
          elif (ivb==3):
             angguess = 60.00
          else:
             angguess = 60.00
       elif (ita==8):
          if (ivb==4):
             if (itc==8):
                angguess = 65.00
             elif (itc==9):
                angguess = 65.00
             elif (itc==15):
                angguess = 60.00
             elif (itc==16):
                angguess = 65.00
             else:
                angguess = 65.00
             
          elif (ivb==3):
             angguess = 50.00
          else:
             angguess = 60.00
          
       else:
          angguess = 60.00
       
    elif (itb==8):
       if (ita==1):
          if (itc==1):
             angguess = 34.05
          elif (itc==6):
             angguess = 65.00
          else:
             angguess = 60.00
          
       elif (ita==6):
          if (itc==6):
             angguess = 88.50
          elif (itc==8):
             if (iva==1 or ivc==1):
                angguess = 122.30
             else:
                angguess = 85.00
             
          elif (itc==15):
             angguess = 80.30
          else:
             angguess = 80.0
          
       else:
          angguess = 80.0
       
    elif (itb==15):
       if (ita==1):
          angguess = 30.0
       elif (ita==6):
          if (itc==6):
             angguess = 75.00
          elif (itc==8):
             angguess = 80.00
          else:
             angguess = 75.00
         
       elif (ita==8):
          if (itc==8):
             if (iva==1 and ivc==1):
                angguess = 89.88
             elif (iva==1 or ivc==1):
                angguess = 75.86
             else:
                angguess = 65.58
             
          else:
             angguess = 70.00
          
       else:
          angguess = 75.00
       
    elif (itb==16):
       if (ita==1):
          angguess = 30.00
       elif (ita==6):
          if (itc==16):
             angguess = 72.00
          else:
             angguess = 80.00
          
       elif (ita==8):
          if (itc==8):
             if (iva==1 and ivc==1):
                angguess = 168.00
             elif (iva==1 or ivc==1):
                angguess = 85.00
             else:
                angguess = 80.00
             
          elif (itc==16):
             angguess = 75.00
          else:
             angguess = 75.00
          
       else:
          angguess = 75.00
       
    elif (ita==1):
       angguess = 35.00
    else:
       angguess = 65.00
    
    angguess = angguess / (angunit*radian**2)
    return angguess
 
 
def BondGuess(poltype,ita,itb,iva,ivb):
    bndunit=1
    if (ita==1):
         if (itb==6):
            if (ivb==3):
               bndguess = 410.0
            elif (ivb==4):
               bndguess = 400.0
            else:
               bndguess = 400.0
            
         elif (itb==7):
            bndguess = 520.0
         elif (itb==8):
            bndguess = 560.0
         elif (itb==9):
            bndguess = 500.0
         elif (itb==14):
            bndguess = 200.0
         elif (itb==15):
            bndguess = 230.0
         elif (itb==16):
            bndguess = 260.0
         else:
            bndguess = 300.0
         
    elif (ita==6):
       if (itb==6):
          if (iva==3 and ivb==3):
             bndguess = 680.0
          elif (iva==4 or ivb==4):
             bndguess = 385.0
          else:
             bndguess = 350.0
          
       elif (itb==7):
          if (iva==3 and ivb==2):
             bndguess = 435.0
          elif (iva==3 and ivb==3):
             bndguess = 250.0
          elif (iva==4):
             bndguess = 400.0
          else:
             bndguess = 450.0
          
       elif (itb==8):
          if (ivb==1):
             bndguess = 680.0
          elif (ivb==2):
             bndguess = 465.0
          else:
             bndguess = 465.0
          
       elif (itb==9):
          bndguess = 350.0
       elif (itb==14):
          bndguess = 350.0
       elif (itb==15):
          bndguess = 350.0
       elif (itb==16):
          bndguess = 216.0
       elif (itb==17):
          bndguess = 350.0
       else:
          bndguess = 450.0
       
    elif (ita==7):
       if (itb==7):
          if (iva==1):
             bndguess = 1613.0
          elif (iva==2 and ivb==2):
             bndguess = 950.0
          else:
             bndguess = 850.0
          
       elif (itb==8):
          if (ivb==1):
             bndguess = 900.0
          else:
             bndguess = 750.0
          
       elif (itb==14):
          bndguess = 450.0
       elif (itb==15):
          bndguess = 500.0
       elif (itb==16):
          bndguess = 550.0
       else:
          bndguess = 600.0
       
    elif (ita==8):
       if (itb==8):
          bndguess = 750.0
       elif (itb==14):
          bndguess = 500.0
       elif (itb==15):
          if (iva==2):
             bndguess = 450.0
          elif (iva==1):
             bndguess = 775.0
          else:
             bndguess = 450.0
          
       elif (itb==16):
          bndguess = 606.0
       elif (itb==17):
          bndguess = 500.0
       else:
          bndguess = 600.0
       
    elif (ita==14):
       if (itb==14):
          bndguess = 400.0
       elif (itb==15):
          bndguess = 450.0
       elif (itb==16):
          bndguess = 500.0
       elif (itb==17):
          bndguess = 650.0
       else:
          bndguess = 450.0
      
    elif (ita==16):
       if (itb==16):
          bndguess = 188.0
       else:
          bndguess = 250.0
       
    elif (ita==17):
       bndguess = 300.0
    else:
       bndguess = 350.0
      
    bndguess = bndguess / bndunit
    return bndguess 

def ReverseDictionaryValueList(poltype,keytovalues):
    valuetokey={}
    for key,value in keytovalues.items():
        for ls in value:
            valuetokey[tuple(ls)]=list(key)
    return valuetokey

def ReverseDictionary(poltype,keytovalues):
    valuetokey={}
    for key,value in keytovalues.items():
        valuetokey[tuple(value)]=list(key)
    return valuetokey

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def ReadDatabaseSmartsMap(poltype,databasepath):
    smartstoatomclass={}
    atomclasstoclassname={}
    atomclasstocomment={}
    temp=open(databasepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            if '#' not in first:
                smarts=linesplit[0]
                if RepresentsInt(linesplit[1])==True:
                    tinkerclass=int(linesplit[1])
                else:
                    tinkerclass=linesplit[1]
                tinkername=linesplit[2]
                comment=' '.join(linesplit[3:])
                comment=comment.replace('\n','').replace('#','').lstrip().rstrip()
                smartstoatomclass[smarts]=tinkerclass
                atomclasstocomment[tinkerclass]=comment
                atomclasstoclassname[tinkerclass]=tinkername
             

    return smartstoatomclass, atomclasstoclassname, atomclasstocomment
 

def MatchAllSmartsToAtomIndices(poltype,smartstoatomclass): #rdkit 0 index based
    atomindextoallsmarts={}
    atomindextoallsmartsmatches={}
    for smarts in smartstoatomclass.keys():
        substructure = Chem.MolFromSmarts(smarts)
        diditmatch=poltype.rdkitmol.HasSubstructMatch(substructure)
        if diditmatch==True:
            matches=list(poltype.rdkitmol.GetSubstructMatches(substructure))
            for match in matches:
                atomindex=match[0]
                if atomindex not in atomindextoallsmarts.keys():
                    atomindextoallsmarts[atomindex]=[]
                    atomindextoallsmartsmatches[atomindex]=[] 
                atomindextoallsmarts[atomindex].append(smarts) 
                if match not in atomindextoallsmartsmatches[atomindex]:
                    atomindextoallsmartsmatches[atomindex].append(match)   
               
    return atomindextoallsmarts,atomindextoallsmartsmatches


def MapIndicesToCommentsAtom(poltype,atomindextoallsmarts,smartstocomment,listofatomsforprm):
    atomcommentstolistofsmartslist={}
    atomindicestolistofatomcomments={}
    for atoms in listofatomsforprm:
        aindex=atoms[0] 
        if aindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            combs = list(itertools.product(asmartslist)) 
            for comb in combs:
                asmarts=comb[0]
                acomment=smartstocomment[asmarts]
                comments=tuple([acomment])
                smartslist=[asmarts]
                if comments not in atomcommentstolistofsmartslist.keys():
                    atomcommentstolistofsmartslist[comments]=[]
                if tuple(atoms) not in atomindicestolistofatomcomments.keys(): 
                    atomindicestolistofatomcomments[tuple(atoms)]=[]
                if smartslist not in atomcommentstolistofsmartslist[comments]: 
                    atomcommentstolistofsmartslist[comments].append(smartslist)   
                if comments not in atomindicestolistofatomcomments[tuple(atoms)]: 
                    atomindicestolistofatomcomments[tuple(atoms)].append(comments)    
    return atomcommentstolistofsmartslist,atomindicestolistofatomcomments


def MapIndicesToCommentsBondAngle(poltype,atomindextoallsmarts,smartstocomment,listofbondsforprm,listofanglesforprm):
    bondcommentstolistofsmartslist={}
    bondindicestolistofbondcomments={}
    for bond in listofbondsforprm:
        aindex=bond[0] 
        bindex=bond[1] 
        asmartslist=atomindextoallsmarts[aindex]
        bsmartslist=atomindextoallsmarts[bindex]
        combs = list(itertools.product(asmartslist,bsmartslist)) 
        for comb in combs:
            asmarts=comb[0]
            bsmarts=comb[1]
            acomment=smartstocomment[asmarts]
            bcomment=smartstocomment[bsmarts]
            comments=tuple([acomment,bcomment])
            smartslist=[asmarts,bsmarts]
            if comments not in bondcommentstolistofsmartslist.keys():
                bondcommentstolistofsmartslist[comments]=[]
            if tuple(bond) not in bondindicestolistofbondcomments.keys(): 
                bondindicestolistofbondcomments[tuple(bond)]=[]
            if smartslist not in bondcommentstolistofsmartslist[comments]: 
                bondcommentstolistofsmartslist[comments].append(smartslist)   
            if comments not in bondindicestolistofbondcomments[tuple(bond)]: 
                bondindicestolistofbondcomments[tuple(bond)].append(comments)    

    anglecommentstolistofsmartslist={}
    angleindicestolistofanglecomments={}
    for angle in listofanglesforprm:
        aindex=angle[0] 
        bindex=angle[1] 
        cindex=angle[2] 
        asmartslist=atomindextoallsmarts[aindex]
        bsmartslist=atomindextoallsmarts[bindex]
        csmartslist=atomindextoallsmarts[cindex]
        combs = list(itertools.product(asmartslist,bsmartslist,csmartslist)) 
        for comb in combs:
            asmarts=comb[0]
            bsmarts=comb[1]
            csmarts=comb[2]
            acomment=smartstocomment[asmarts]
            bcomment=smartstocomment[bsmarts]
            ccomment=smartstocomment[csmarts]
            comments=tuple([acomment,bcomment,ccomment])
            smartslist=[asmarts,bsmarts,csmarts]
            if comments not in anglecommentstolistofsmartslist.keys():
                anglecommentstolistofsmartslist[comments]=[]
            if tuple(angle) not in angleindicestolistofanglecomments.keys(): 
                angleindicestolistofanglecomments[tuple(angle)]=[]
            if smartslist not in anglecommentstolistofsmartslist[comments]: 
                anglecommentstolistofsmartslist[comments].append(smartslist)   
            if comments not in angleindicestolistofanglecomments[tuple(angle)]: 
                angleindicestolistofanglecomments[tuple(angle)].append(comments)    

    return bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments




def MapIndicesToClasses(poltype,atomindextoallsmarts,smartstoatomclass,listofbondsforprm,listofanglesforprm,planarbonds):
    bondclassestolistofsmartslist={}
    angleclassestolistofsmartslist={}
    strbndclassestolistofsmartslist={}
    opbendclassestolistofsmartslist={}
    bondindicestolistofbondclasses={}
    angleindicestolistofangleclasses={}
    strbndindicestolistofstrbndclasses={}
    opbendindicestolistofopbendclasses={}
    for bond in listofbondsforprm:
        aindex=bond[0]
        bindex=bond[1] 
        if aindex in atomindextoallsmarts.keys() and bindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            bsmartslist=atomindextoallsmarts[bindex]   
            combs = list(itertools.product(asmartslist, bsmartslist)) 
            for comb in combs:
                asmarts=comb[0]
                bsmarts=comb[1]
                aclass=smartstoatomclass[asmarts]
                bclass=smartstoatomclass[bsmarts]
                bondclasses=tuple([aclass,bclass])
                smartslist=[asmarts,bsmarts]
                if bondclasses not in bondclassestolistofsmartslist.keys():
                    bondclassestolistofsmartslist[bondclasses]=[]
                if tuple(bond) not in bondindicestolistofbondclasses.keys(): 
                    bondindicestolistofbondclasses[tuple(bond)]=[]
                if smartslist not in bondclassestolistofsmartslist[bondclasses]: 
                    bondclassestolistofsmartslist[bondclasses].append(smartslist)   
                 
                if bondclasses not in bondindicestolistofbondclasses[tuple(bond)]: 
                    bondindicestolistofbondclasses[tuple(bond)].append(bondclasses)    

                if bond in planarbonds or bond[::-1] in planarbonds:
                    if bondclasses not in opbendclassestolistofsmartslist.keys():
                        opbendclassestolistofsmartslist[bondclasses]=[]
                    if tuple(bond) not in opbendindicestolistofopbendclasses.keys():
                        opbendindicestolistofopbendclasses[tuple(bond)]=[]
                    if smartslist not in opbendclassestolistofsmartslist[bondclasses]: 
                        opbendclassestolistofsmartslist[bondclasses].append(smartslist)   
                    if bondclasses not in opbendindicestolistofopbendclasses[tuple(bond)]:
                        opbendindicestolistofopbendclasses[tuple(bond)].append(bondclasses)    


    for angle in listofanglesforprm:
        aindex=angle[0]
        bindex=angle[1] 
        cindex=angle[2] 
        if aindex in atomindextoallsmarts.keys() and bindex in atomindextoallsmarts.keys() and cindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            bsmartslist=atomindextoallsmarts[bindex]   
            csmartslist=atomindextoallsmarts[cindex]   
            combs = list(itertools.product(asmartslist,bsmartslist,csmartslist)) 
            for comb in combs:
                asmarts=comb[0]
                bsmarts=comb[1]
                csmarts=comb[2]
                aclass=smartstoatomclass[asmarts]
                bclass=smartstoatomclass[bsmarts]
                cclass=smartstoatomclass[csmarts]
                angleclasses=tuple([aclass,bclass,cclass])
                smartslist=[asmarts,bsmarts,csmarts]
                if angleclasses not in angleclassestolistofsmartslist.keys():
                    angleclassestolistofsmartslist[angleclasses]=[]
                    strbndclassestolistofsmartslist[angleclasses]=[]
                if tuple(angle) not in angleindicestolistofangleclasses.keys():
                    angleindicestolistofangleclasses[tuple(angle)]=[]
                if tuple(angle) not in strbndindicestolistofstrbndclasses.keys():
                    strbndindicestolistofstrbndclasses[tuple(angle)]=[]
                if smartslist not in angleclassestolistofsmartslist[angleclasses]: 
                    angleclassestolistofsmartslist[angleclasses].append(smartslist)     
                if smartslist not in strbndclassestolistofsmartslist[angleclasses]:
                    strbndclassestolistofsmartslist[angleclasses].append(smartslist)     
                if angleclasses not in angleindicestolistofangleclasses[tuple(angle)]: 
                    angleindicestolistofangleclasses[tuple(angle)].append(angleclasses)     
                if angleclasses not in strbndindicestolistofstrbndclasses[tuple(angle)]: 
                    strbndindicestolistofstrbndclasses[tuple(angle)].append(angleclasses)     

    return bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses

def SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments):
    atomcommentstoparameters={}
    temp=open(poltype.latestsmallmoleculepolarizeprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if '#' not in line:
            comment=linesplit[0]
            prm=linesplit[1]
            atomcommentstoparameters[tuple([comment])]=[prm]

    atomindicestolistofatomcomments=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist,atomcommentstoparameters)
    return atomindicestolistofatomcomments,atomcommentstoparameters


def SearchForParametersViaCommentsNonBondedAMOEBAPlus(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments):
    atomcommentstocpparameters={}
    atomcommentstoctparameters={}
    atomcommentstovdwparameters={}
    atomindicestolistofatomcommentscp=atomindicestolistofatomcomments.copy()
    atomindicestolistofatomcommentsct=atomindicestolistofatomcomments.copy()
    atomindicestolistofatomcommentsvdw=atomindicestolistofatomcomments.copy()

    temp=open(poltype.amoebaplusnonbondedprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        first=linesplit[0]
        if '#' not in first:
            comment=linesplit[0]
            cpprms=linesplit[1:2+1]
            linesplit[3]=str(float(linesplit[3])*3.01147) # conversion factor for ct
            ctprms=linesplit[3:4+1]
            vdwprms=linesplit[5:7+1]
            atomcommentstocpparameters[tuple([comment])]=cpprms
            atomcommentstoctparameters[tuple([comment])]=ctprms
            atomcommentstovdwparameters[tuple([comment])]=vdwprms


    atomindicestolistofatomcommentscp=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcommentscp,atomcommentstolistofsmartslist,atomcommentstocpparameters)
    atomindicestolistofatomcommentsct=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcommentsct,atomcommentstolistofsmartslist,atomcommentstoctparameters)
    atomindicestolistofatomcommentsvdw=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcommentsvdw,atomcommentstolistofsmartslist,atomcommentstovdwparameters)



    return atomindicestolistofatomcommentscp,atomcommentstocpparameters,atomindicestolistofatomcommentsct,atomcommentstoctparameters,atomindicestolistofatomcommentsvdw,atomcommentstovdwparameters


def SearchForParametersViaCommentsChargeFlux(poltype,bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments):
    bondcommentstocfparameters={}
    anglecommentstocfparameters={}
    temp=open(poltype.amoebapluscfprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            if '#' not in first:
                if 'cflux-b' in line:
                    comments=linesplit[2:]
                    prms=[linesplit[1]] 
                    bondcommentstocfparameters[tuple(comments)]=prms
                elif 'cflux-a' in line:
                    prms=linesplit[1:4+1]
                    comments=linesplit[5:]
                    anglecommentstocfparameters[tuple(comments)]=prms

    bondindicestolistofbondcomments=RemoveIndicesThatDontHaveParameters(poltype,bondindicestolistofbondcomments,bondcommentstolistofsmartslist,bondcommentstocfparameters)
    angleindicestolistofanglecomments=RemoveIndicesThatDontHaveParameters(poltype,angleindicestolistofanglecomments,anglecommentstolistofsmartslist,anglecommentstocfparameters)

    return bondindicestolistofbondcomments,bondcommentstocfparameters,angleindicestolistofanglecomments,anglecommentstocfparameters




def SearchForParameters(poltype,bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses):
    bondclassestoparameters={}
    angleclassestoparameters={}
    strbndclassestoparameters={}
    opbendclassestoparameters={}
    temp=open(poltype.latestsmallmoleculeprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'bond' in line:
            bondclasses=tuple([int(linesplit[1]),int(linesplit[2])]) 
            prms=[float(linesplit[3]),float(linesplit[4])]
            bondclassestoparameters[bondclasses]=prms
        elif 'angle' in line:
            angleclasses=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]) 
            prms=[float(linesplit[4]),float(linesplit[5])]
            angleclassestoparameters[angleclasses]=prms
        elif 'strbnd' in line:
            angleclasses=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]) 
            prms=[float(linesplit[4]),float(linesplit[5])]
            strbndclassestoparameters[angleclasses]=prms
        elif 'opbend' in line:
            bondclasses=tuple([int(linesplit[1]),int(linesplit[2])]) 
            prms=[float(linesplit[5])]
            opbendclassestoparameters[bondclasses]=prms

    bondindicestolistofbondclasses=RemoveIndicesThatDontHaveParameters(poltype,bondindicestolistofbondclasses,bondclassestolistofsmartslist,bondclassestoparameters)
    angleindicestolistofangleclasses=RemoveIndicesThatDontHaveParameters(poltype,angleindicestolistofangleclasses,angleclassestolistofsmartslist,angleclassestoparameters)
    strbndindicestolistofstrbndclasses=RemoveIndicesThatDontHaveParameters(poltype,strbndindicestolistofstrbndclasses,strbndclassestolistofsmartslist,strbndclassestoparameters)
    opbendindicestolistofopbendclasses=RemoveIndicesThatDontHaveParameters(poltype,opbendindicestolistofopbendclasses,opbendclassestolistofsmartslist,opbendclassestoparameters)

 
    return bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses,bondclassestoparameters,angleclassestoparameters,strbndclassestoparameters,opbendclassestoparameters

def RemoveIndicesThatDontHaveParameters(poltype,indicestolistofclasses,classestolistofsmartslist,classestoparameters):
    indicestodelete=[]
    for i in range(len(indicestolistofclasses.keys())):
        indices=list(indicestolistofclasses.keys())[i]
        listofclasses=indicestolistofclasses[indices]
        classesmissing=[]
        classesexisting=[]
        for classes in listofclasses:
            if classes not in classestoparameters.keys() or classes[::-1] not in classestoparameters.keys():
                classesmissing.append(classes)
            else:
                classesexisting.append(classes)
        if len(classesmissing)==len(listofclasses):
            indicestodelete.append(indices)
        indicestolistofclasses[indices]=classesexisting

    for indices in indicestodelete:
        del indicestolistofclasses[indices]
    return indicestolistofclasses

def FindBestSMARTSMatch(poltype,indicestolistofclasses,classestolistofsmartslist):
    indicestoclasses={}
    indicestosmartslist={}
    for indices,listofclasses in indicestolistofclasses.items():
        fragmentsarray=[]
        classesarray=[]
        smartslistarray=[] 
        for classes in listofclasses:
            listofsmartslist=classestolistofsmartslist[classes]
            for smartslist in listofsmartslist:
                matchedindices=MatchAllSMARTS(poltype,smartslist,indices) 
                fragmentmol=CreateFragment(poltype,matchedindices)
                fragmentsarray.append(fragmentmol)
                classesarray.append(classes)
                smartslistarray.append(smartslist)
        classes,smartslist=DetermineBestSMARTSMatch(poltype,fragmentsarray,classesarray,smartslistarray)         
        indicestoclasses[indices]=classes
        indicestosmartslist[indices]=smartslist

    return indicestoclasses,indicestosmartslist

def DetermineBestSMARTSMatch(poltype,listoffragments,listofclasses,listofsmartslist):  
    tanimotoarray=[]     
    for i in range(len(listoffragments)):
        fragment=listoffragments[i]
        classes=listofclasses[i]
        smartslist=listofsmartslist[i]
        ms=[poltype.rdkitmol,fragment]
        fps = [Chem.RDKFingerprint(x) for x in ms]
        tanimoto=DataStructs.FingerprintSimilarity(fps[0],fps[1], metric=DataStructs.DiceSimilarity) 
        tanimotoarray.append(tanimoto)
    maxvalue=max(tanimotoarray)
    nums=tanimotoarray.count(maxvalue)
    if nums>1:
        smartslisttocompare=[]
        for i in range(len(listofsmartslist)):
            smartslist=listofsmartslist[i]
            tanimoto=tanimotoarray[i]
            if tanimoto==maxvalue:
                smartslisttocompare.append(smartslist)
        smartslist=ChooseMostDescriptiveSMARTSList(poltype,smartslisttocompare) 
        index=listofsmartslist.index(smartslist)
        classes=listofclasses[index]
    else:
        maxindex=tanimotoarray.index(maxvalue)
        classes=listofclasses[maxindex]
        smartslist=listofsmartslist[maxindex]
    return classes,smartslist 


def ChooseMostDescriptiveSMARTSList(poltype,smartslisttocompare):
    lengtharray=[]
    for smartslist in smartslisttocompare:
        lengths=[len(e) for e in smartslist]
        total=sum(lengths)
        lengtharray.append(total)
    maxlength=max(lengtharray)
    maxidx=lengtharray.index(maxlength)
    smartslist=smartslisttocompare[maxidx]
    return smartslist


def CreateFragment(poltype,matchedindices):
    m = Chem.Mol()
    em = Chem.EditableMol(m)
    oldindextonewindex={}
    for i,idx in enumerate(matchedindices):
        oldatom=poltype.rdkitmol.GetAtomWithIdx(idx)
        newindex=em.AddAtom(oldatom)
        oldindextonewindex[idx]=newindex

    for bond in poltype.rdkitmol.GetBonds():
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        if oendidx not in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            continue
        elif oendidx in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            continue
        elif oendidx not in oldindextonewindex.keys() and obgnidx in oldindextonewindex.keys():
            continue

        endidx=oldindextonewindex[oendidx]
        bgnidx=oldindextonewindex[obgnidx]
        bondorder=bond.GetBondType()
        em.AddBond(bgnidx,endidx,bondorder)

    return em.GetMol()



def MatchAllSMARTS(poltype,smartslist,indices):
    matchedindices=[]
    for idx in range(len(smartslist)):
        smarts=smartslist[idx]
        index=indices[idx]
        substructure = Chem.MolFromSmarts(smarts)
        diditmatch=poltype.rdkitmol.HasSubstructMatch(substructure)
        if diditmatch==True:
            matches=poltype.rdkitmol.GetSubstructMatches(substructure)
            for match in matches:
                if index in match:
                    for atomindex in match:
                        if atomindex not in matchedindices:
                            matchedindices.append(atomindex)
                    break
        else:
            raise ValueError('SMARTS match does not exist')


    return matchedindices


def GrabNewParameters(poltype,indicestoclasses,classestoparameters,keyword,indicestosmartslist,atomclasstocomment):
    indicestopoltypeclasses={}
    newpoltypeclassestocomments={}
    newpoltypeclassestosmartslist={}

    newprms=[]
    for indices,classes in indicestoclasses.items():
        smartslist=indicestosmartslist[indices]
        comments=[atomclasstocomment[k] for k in classes]
        babelindices=[k+1 for k in indices]
        poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
        indicestopoltypeclasses[indices]=poltypeclasses
        newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
        newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
        poltypeclasses=[str(k) for k in poltypeclasses]
        if keyword=='opbend':
            poltypeclasses.append('0')
            poltypeclasses.append('0')
        parameters=classestoparameters[classes]
        parameters=[str(k) for k in parameters]
        poltypeclassesstring=' '.join(poltypeclasses)
        parametersstring=' '.join(parameters)
        newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
        newprms.append(newprm)
        

    return indicestopoltypeclasses,newprms,newpoltypeclassestocomments,newpoltypeclassestosmartslist


def RemoveOldParametersKeepNewParameters(poltype,prms,newindicestopoltypeclasses,keyword,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,removefromdic=True):
    newprms=[]
    dellist=[]
    appendlist=[]
    for prm in prms:
        prmsplit=prm.split()
        if keyword=='bond' or keyword=='opbend':
            endindex=2
        elif keyword=='angle' or keyword=='strbnd':
            endindex=3
        prmsplit=prmsplit[1:endindex+1]
        classes=[int(k) for k in prmsplit]
        if classes in newindicestopoltypeclasses.values() or classes[::-1] in newindicestopoltypeclasses.values():
            pass
            if removefromdic==True:
                poltypeclasseslist,poltypeclasses=SearchForPoltypeClasses(poltype,prm,poltypeclassestoparametersmartsatomorders.keys())
                if poltypeclasseslist!=None:
                    newpoltypeclasseslist=[]
                    for cls in poltypeclasseslist:
                        if cls==poltypeclasses:
                            pass
                        else:
                            newpoltypeclasseslist.append(cls)
                    newpoltypeclasseslist=tuple(newpoltypeclasseslist)
                    if poltypeclasseslist not in dellist:
                        dellist.append(poltypeclasseslist)
                        appendlist.append(newpoltypeclasseslist)
        else:
            newprms.append(prm)
    for i in range(len(dellist)):
        item=dellist[i]
        newitem=appendlist[i]
        value=poltypeclassestoparametersmartsatomorders[item]
        if item in poltypeclassestoparametersmartsatomorders.keys():
            del poltypeclassestoparametersmartsatomorders[item]
        if item in poltypeclassestosmartsatomorders.keys():
            del poltypeclassestosmartsatomorders[item]
        if item in poltypeclassestoelementtinkerdescrips.items():
            del poltypeclassestoelementtinkerdescrips[item]
        poltypeclassestoparametersmartsatomorders[newitem]=value 
        poltypeclassestosmartsatomorders[newitem]=value
        poltypeclassestoelementtinkerdescrips[newitem]=value


    return newprms,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips

def MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstocommentpolar):
    smartstocomment={}
    for smarts,atomclass in smartstoatomclasspolar.items():
        comment=atomclasstocommentpolar[atomclass]
        smartstocomment[smarts]=comment

    return smartstocomment

def GrabSmartsToSoluteRadiiMap(poltype):
    smartstosoluteradiiprms={}
    temp=open(poltype.smartstosoluteradiimap,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        smarts=linesplit[0]
        prms=linesplit[1:]
        smartstosoluteradiiprms[smarts]=prms
    return smartstosoluteradiiprms


def MakeListOfListValues(poltype,keytovalue):
    newkeytovalue={}
    for key,value in keytovalue.items():
        newvalues=[]
        for item in value:
            newitem=[item]
            newvalues.append(newitem)
        newkeytovalue[key]=newvalues    

    return newkeytovalue

def GenerateSoluteParameters(poltype,atomindicestosmartslist,smartstosoluteradiiprms):
    soluteprms=[]
    for atomindices,smartslist in atomindicestosmartslist.items():
        atomindex=atomindices[0]
        smarts=smartslist[0]
        babelatomindex=atomindex+1
        poltypeclass=poltype.idxtosymclass[babelatomindex]
        prms=smartstosoluteradiiprms[smarts]
        prmsline=' '.join(prms)
        line='SOLUTE'+' '+str(poltypeclass)+' '+smarts+' '+prmsline+'\n'
        if line not in soluteprms:
            soluteprms.append(line)

    return soluteprms

def RemoveIndicesMatchedFromNewDatabase(poltype,indicestosmartsatomorders,indicestopoltypeclasses):
    removelist=[]
    for indices,smartsatomorders in indicestosmartsatomorders.items():
        if indices in indicestopoltypeclasses.keys() or indices[::-1] in indicestopoltypeclasses.keys():
            removelist.append(indices)
    for indices in removelist:
        del indicestosmartsatomorders[indices]


    return indicestosmartsatomorders



def GrabSmallMoleculeAMOEBAParameters(poltype,optmol,mol,rdkitmol):
    
    listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForParameters(poltype,mol)
    if poltype.forcefield=='AMOEBA+':
        smartstoatomclasscf, atomclasstoclassnamecf, atomclasstocommentcf=ReadDatabaseSmartsMap(poltype,poltype.amoebapluscfsmartstocommentmap) 
        atomindextoallsmartscf,atomindextoallsmartsmatchescf=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasscf)
        smartstocommentcf=MapSMARTSToComments(poltype,smartstoatomclasscf,atomclasstoclassnamecf)
        bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments=MapIndicesToCommentsBondAngle(poltype,atomindextoallsmartscf,smartstocommentcf,listofbondsforprm,listofanglesforprm)
        bondindicestolistofbondcomments,bondcommentstocfparameters,angleindicestolistofanglecomments,anglecommentstocfparameters=SearchForParametersViaCommentsChargeFlux(poltype,bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments)

        bondindicestocommentscf,bondindicestosmartslistcf=FindBestSMARTSMatch(poltype,bondindicestolistofbondcomments,bondcommentstolistofsmartslist)
        angleindicestocommentscf,angleindicestosmartslistcf=FindBestSMARTSMatch(poltype,angleindicestolistofanglecomments,anglecommentstolistofsmartslist)
        commentlistcf=list(atomclasstoclassnamecf.values())
        atomcommentstocommentcf=dict(zip(commentlistcf,commentlistcf))
        bondcfindicestopoltypeclasses,bondcfprms,bondcfpoltypecommentstocomments,bondcfpoltypecommentstosmartslist=GrabNewParameters(poltype,bondindicestocommentscf,bondcommentstocfparameters,'cflux-b',bondindicestosmartslistcf,atomcommentstocommentcf) 
        anglecfindicestopoltypeclasses,anglecfprms,anglecfpoltypecommentstocomments,anglecfpoltypecommentstosmartslist=GrabNewParameters(poltype,angleindicestocommentscf,anglecommentstocfparameters,'cflux-a',angleindicestosmartslistcf,atomcommentstocommentcf) 
        smartstoatomclassnonbondedplus, atomclasstoclassnamenonbondedplus, atomclasstocommentnonbondedplus=ReadDatabaseSmartsMap(poltype,poltype.amoebaplusnonbondedsmartstocommentmap) 
        atomindextoallsmartsnonbondedplus,atomindextoallsmartsmatchesnonbondedplus=MatchAllSmartsToAtomIndices(poltype,smartstoatomclassnonbondedplus)
        smartstocommentnonbondedplus=MapSMARTSToComments(poltype,smartstoatomclassnonbondedplus,atomclasstoclassnamenonbondedplus)
        atomcommentstolistofsmartslistnonbondedplus,atomindicestolistofatomcommentsnonbondedplus=MapIndicesToCommentsAtom(poltype,atomindextoallsmartsnonbondedplus,smartstocommentnonbondedplus,listofatomsforprm)

        atomindicestolistofatomcommentscp,atomcommentstocpparameters,atomindicestolistofatomcommentsct,atomcommentstoctparameters,atomindicestolistofatomcommentsvdw,atomcommentstovdwparameters=SearchForParametersViaCommentsNonBondedAMOEBAPlus(poltype,atomcommentstolistofsmartslistnonbondedplus,atomindicestolistofatomcommentsnonbondedplus)
        atomindicestocommentscp,atomindicestosmartslistcp=FindBestSMARTSMatch(poltype,atomindicestolistofatomcommentscp,atomcommentstolistofsmartslistnonbondedplus)
        atomindicestocommentsct,atomindicestosmartslistct=FindBestSMARTSMatch(poltype,atomindicestolistofatomcommentsct,atomcommentstolistofsmartslistnonbondedplus)
        atomindicestocommentsvdw,atomindicestosmartslistvdw=FindBestSMARTSMatch(poltype,atomindicestolistofatomcommentsvdw,atomcommentstolistofsmartslistnonbondedplus)

        commentlistnonbondedplus=list(atomclasstoclassnamenonbondedplus.values())
        atomcommentstocommentnonbondedplus=dict(zip(commentlistnonbondedplus,commentlistnonbondedplus))

        cpindicestopoltypeclasses,cpprms,cppoltypecommentstocomments,cppoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocommentscp,atomcommentstocpparameters,'chgpen',atomindicestosmartslistcp,atomcommentstocommentnonbondedplus) 
        ctindicestopoltypeclasses,ctprms,ctpoltypecommentstocomments,ctpoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocommentsct,atomcommentstoctparameters,'chgtrn',atomindicestosmartslistct,atomcommentstocommentnonbondedplus) 
        newvdwindicestopoltypeclasses,newvdwprms,newvdwpoltypecommentstocomments,newvdwpoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocommentsvdw,atomcommentstovdwparameters,'vdw',atomindicestosmartslistvdw,atomcommentstocommentnonbondedplus) 


    smartstosoluteradiiprms=GrabSmartsToSoluteRadiiMap(poltype)   
    atomindextoallsmartssolute,atomindextoallsmartsmatchessolute=MatchAllSmartsToAtomIndices(poltype,smartstosoluteradiiprms)
    atomindices=list(atomindextoallsmartssolute.keys())
    atomindices=[tuple([k]) for k in atomindices]
    atomindicestolistofatomindices=dict(zip(atomindices,atomindices))
    atomindextoallsmartssolute=MakeListOfListValues(poltype,atomindextoallsmartssolute)
    atomindicestoindices,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomindices,atomindextoallsmartssolute)
    soluteprms=GenerateSoluteParameters(poltype,atomindicestosmartslist,smartstosoluteradiiprms)
    smartstoatomclasspolar, atomclasstoclassnamepolar, atomclasstocommentpolar=ReadDatabaseSmartsMap(poltype,poltype.latestsmallmoleculesmartstotypespolarize) 
    atomindextoallsmartspolar,atomindextoallsmartsmatchespolar=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasspolar)
    smartstocomment=MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstoclassnamepolar)
    atomcommentstolistofsmartslist,atomindicestolistofatomcomments=MapIndicesToCommentsAtom(poltype,atomindextoallsmartspolar,smartstocomment,listofatomsforprm)
    atomindicestolistofatomcomments,atomcommentstoparameters=SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments)
    atomindicestocomments,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist)
    commentlist=list(atomclasstoclassnamepolar.values())
    atomcommentstocomment=dict(zip(commentlist,commentlist))
    newpolarindicestopoltypeclasses,newpolarprms,newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocomments,atomcommentstoparameters,'polarize',atomindicestosmartslist,atomcommentstocomment) 

    smartstoatomclass, atomclasstoclassname, atomclasstocomment=ReadDatabaseSmartsMap(poltype,poltype.latestsmallmoleculesmartstotinkerclass) 
    atomindextoallsmarts,atomindextoallsmartsmatches=MatchAllSmartsToAtomIndices(poltype,smartstoatomclass)
    bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters,tortorsmartsatomordertoparameters,tortorsmartsatomordertogrid=ReadExternalDatabase(poltype)
    smartsatomordertoelementtinkerdescrip=ReadSmallMoleculeLib(poltype,poltype.smallmoleculesmartstotinkerdescrip)
    elementtinkerdescriptotinkertype,tinkertypetoclass=GrabTypeAndClassNumbers(poltype,poltype.smallmoleculeprmlib)
    planarbonds=GrabPlanarBonds(poltype,listofbondsforprm,mol)
    bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses=MapIndicesToClasses(poltype,atomindextoallsmarts,smartstoatomclass,listofbondsforprm,listofanglesforprm,planarbonds)
    bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses,bondclassestoparameters,angleclassestoparameters,strbndclassestoparameters,opbendclassestoparameters=SearchForParameters(poltype,bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses)
    bondindicestoclasses,bondindicestosmartslist=FindBestSMARTSMatch(poltype,bondindicestolistofbondclasses,bondclassestolistofsmartslist)
    angleindicestoclasses,angleindicestosmartslist=FindBestSMARTSMatch(poltype,angleindicestolistofangleclasses,angleclassestolistofsmartslist)
    strbndindicestoclasses,strbndindicestosmartslist=FindBestSMARTSMatch(poltype,strbndindicestolistofstrbndclasses,strbndclassestolistofsmartslist)
    opbendindicestoclasses,opbendindicestosmartslist=FindBestSMARTSMatch(poltype,opbendindicestolistofopbendclasses,opbendclassestolistofsmartslist)
    newangleindicestopoltypeclasses,newangleprms,newanglepoltypeclassestocomments,newanglepoltypeclassestosmartslist=GrabNewParameters(poltype,angleindicestoclasses,angleclassestoparameters,'angle',angleindicestosmartslist,atomclasstocomment) 
    newbondindicestopoltypeclasses,newbondprms,newbondpoltypeclassestocomments,newbondpoltypeclassestosmartslist=GrabNewParameters(poltype,bondindicestoclasses,bondclassestoparameters,'bond',bondindicestosmartslist,atomclasstocomment) 
    newstrbndindicestopoltypeclasses,newstrbndprms,newstrbndpoltypeclassestocomments,newstrbndpoltypeclassestosmartslist=GrabNewParameters(poltype,strbndindicestoclasses,strbndclassestoparameters,'strbnd',strbndindicestosmartslist,atomclasstocomment) 
    newopbendindicestopoltypeclasses,newopbendprms,newopbendpoltypeclassestocomments,newopbendpoltypeclassestosmartslist=GrabNewParameters(poltype,opbendindicestoclasses,opbendclassestoparameters,'opbend',opbendindicestosmartslist,atomclasstocomment) 
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
    tortorindicestoextsmartsmatchlength,tortorindicestoextsmarts=MatchExternalSMARTSToMolecule(poltype,rdkitmol,tortorsmartsatomordertoparameters)

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
 
    originalbondindicestosmartsatomorders=bondindicestosmartsatomorders.copy()
    originalangleindicestosmartsatomorders=angleindicestosmartsatomorders.copy()

    formissingangleindicestosmartsatomorders=RemoveIndicesMatchedFromNewDatabase(poltype,angleindicestosmartsatomorders,newangleindicestopoltypeclasses) 

    formissingbondindicestosmartsatomorders=RemoveIndicesMatchedFromNewDatabase(poltype,bondindicestosmartsatomorders,newbondindicestopoltypeclasses) 
    indextoneighbidxs=FindAllNeighborIndexes(poltype,rdkitmol)
    bondindicestosmartsatomorders=originalbondindicestosmartsatomorders.copy()
    angleindicestosmartsatomorders=originalangleindicestosmartsatomorders.copy()

    totalbondscollector=FindAllConsecutiveRotatableBonds(poltype,mol,listofbondsforprm)
    tortorsmissing=FindMissingTorTors(poltype,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,rdkitmol,mol,indextoneighbidxs,totalbondscollector)
    torsionsmissing,poormatchingaromatictorsions=FindMissingTorsions(poltype,torsionindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
    torsionsmissing=FindAdjacentMissingTorsionsForTorTor(poltype,torsionsmissing,totalbondscollector,tortorsmissing)

    atomindextosmartsatomorder=AddExternalDatabaseMatches(poltype, atomindextosmartsatomorder,vdwindicestoextsmarts,vdwsmartsatomordertoparameters)
    vdwmissing=FindMissingParameters(poltype,atomindextosmartsatomorder,rdkitmol,mol,indextoneighbidxs)
    missingvdwatomindices=ReduceMissingVdwByTypes(poltype,vdwmissing)
    bondmissing=FindMissingParameters(poltype,formissingbondindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
    anglemissing=FindMissingParameters(poltype,formissingangleindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
    anglemissingindicestotinkerclasses=PruneDictionary(poltype,anglemissing,angleindicestotinkerclasses)
    angletinkerclassestoexampleindices=ReverseDictionary(poltype,angleindicestotinkerclasses)
    bondmissingindicestotinkerclasses=PruneDictionary(poltype,bondmissing,bondindicestotinkerclasses)
    bondtinkerclassestoexampleindices=ReverseDictionary(poltype,bondmissingindicestotinkerclasses)

    torsionsmissingindicestotinkerclasses=PruneDictionary(poltype,torsionsmissing,torsionindicestotinkerclasses)
    atomtinkerclasstopoltypeclass=TinkerClassesToPoltypeClasses(poltype,atomindextotinkerclass)
    bondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondindicestotinkerclasses)
    
    arotorsionsmissingindicestotinkerclasses=PruneDictionary(poltype,poormatchingaromatictorsions,torsionindicestotinkerclasses)
    planarbondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarbondindicestotinkerclasses)
    opbendtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,opbendbondindicestotinkerclasses)
    opbendtinkerclassestopoltypeclasses=AddReverseKeys(poltype,opbendtinkerclassestopoltypeclasses)
    opbendbondindicestotrigonalcenterbools=CheckTrigonalCenters(poltype,listofbondsforprm,mol)
    opbendtinkerclassestotrigonalcenterbools=TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools)
    angletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,angleindicestotinkerclasses)
    anglepoltypeclassestotinkerclasses=ReverseDictionaryValueList(poltype,angletinkerclassestopoltypeclasses)
    planarangletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarangleindicestotinkerclasses)
    torsiontinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionindicestotinkerclasses)
    torsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionsmissingindicestotinkerclasses)
    arotorsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,arotorsionsmissingindicestotinkerclasses)
    anglemissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,anglemissingindicestotinkerclasses)
    bondmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondmissingindicestotinkerclasses)
    
    bondpoltypeclassestotinkerclasses=ReverseDictionaryValueList(poltype,bondtinkerclassestopoltypeclasses)
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
    angleprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,angleprms,newangleindicestopoltypeclasses,'angle',anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips)
    opbendprms,blankbondpoltypeclassestoparametersmartsatomorders,blankbondpoltypeclassestosmartsatomorders,blankbondpoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,opbendprms,newopbendindicestopoltypeclasses,'opbend',bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips,removefromdic=False)

    opbendpoltypeclassestoparametersmartsatomorders=bondpoltypeclassestoparametersmartsatomorders.copy()
    opbendpoltypeclassestosmartsatomorders=bondpoltypeclassestosmartsatomorders.copy()
    opbendpoltypeclassestoelementtinkerdescrips=bondpoltypeclassestoelementtinkerdescrips.copy()
    bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,bondprms,newbondindicestopoltypeclasses,'bond',bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips)

    strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,strbndprms,newstrbndindicestopoltypeclasses,'strbnd',anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips)
    bondprms,bondpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,bondprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'bond')
    angleprms,anglepoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,angleprms,angleindicestoextsmarts,anglesmartsatomordertoparameters,'angle')
    strbndprms,strbndpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,strbndprms,angleindicestoextsmarts,strbndsmartsatomordertoparameters,'strbnd')
    torsionprms,torsionpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,torsionprms,torsionindicestoextsmarts,torsionsmartsatomordertoparameters,'torsion')
    opbendprms,opbendpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,opbendprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'opbend')
    vdwprms,vdwpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,vdwprms,vdwindicestoextsmarts,vdwsmartsatomordertoparameters,'vdw')
    tortorprms=[]
    tortorprms,tortorpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,tortorprms,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,'tortors',tortorsmartsatomordertogrid)
    

    strbndprms=ZeroOutMissingStrbnd(poltype,anglemissingtinkerclassestopoltypeclasses,strbndprms)
    angleprms=AssignAngleGuessParameters(poltype,angletinkerclassestoexampleindices,anglemissingtinkerclassestopoltypeclasses,angleprms)
    bondprms=AssignBondGuessParameters(poltype,bondtinkerclassestoexampleindices,bondmissingtinkerclassestopoltypeclasses,bondprms)
    angleprms.extend(newangleprms)
    bondprms.extend(newbondprms)
    strbndprms.extend(newstrbndprms)
    opbendprms.extend(newopbendprms)
    angleprms=ModifyAngleKeywords(poltype,angleprms,planarangletinkerclassestopoltypeclasses)
    bondlistbabel=ConvertToBabelList(poltype,listofbondsforprm,bondindicestoclasses)
    anglelistbabel=ConvertToBabelList(poltype,listofanglesforprm,angleindicestoclasses)
    bondprms=AddOptimizedBondLengths(poltype,optmol,bondprms,bondlistbabel)
    angleprms=AddOptimizedAngleLengths(poltype,optmol,angleprms,anglelistbabel)
    torsionprms=ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms)
    torsionprms=DefaultAromaticMissingTorsions(poltype,arotorsionsmissingtinkerclassestopoltypeclasses,torsionprms)
    torsionkeystringtoparameters=GrabTorsionParameterCoefficients(poltype,torsionprms)
    potentialmissingopbendprmtypes=FindPotentialMissingParameterTypes(poltype,opbendprms,planarbondtinkerclassestopoltypeclasses)
    potentialmissingopbendprmindices=ConvertPoltypeClassesToIndices(poltype,potentialmissingopbendprmtypes)
    potentialmissingopbendprmindices=FilterIndices(poltype,potentialmissingopbendprmindices,planarbonds)
    missingopbendprmindices=CheckIfParametersExist(poltype,potentialmissingopbendprmindices,opbendprms)
    defaultvalues=None
    if len(missingopbendprmindices)!=0:
        newopbendprms,defaultvalues=DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools)
        opbendprms.extend(newopbendprms)
    polarprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newpolarprms,{},{},{},{},newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist)

    vdwprmstotransferinfo=MapParameterLineToTransferInfo(poltype,vdwprms,atompoltypeclasstoparametersmartsatomorder,atompoltypeclasstosmartsatomorder,atompoltypeclassestoelementtinkerdescrip,vdwpoltypeclassestosmartsatomordersext,{},{})
    if poltype.forcefield=='AMOEBA+':
        amoebaplusvdwprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newvdwprms,{},{},{},{},newvdwpoltypecommentstocomments,newvdwpoltypecommentstosmartslist)
        ctprmstotransferinfo=MapParameterLineToTransferInfo(poltype,ctprms,{},{},{},{},ctpoltypecommentstocomments,ctpoltypecommentstosmartslist)
        cpprmstotransferinfo=MapParameterLineToTransferInfo(poltype,cpprms,{},{},{},{},cppoltypecommentstocomments,cppoltypecommentstosmartslist)

        bondcfprmstotransferinfo=MapParameterLineToTransferInfo(poltype,bondcfprms,{},{},{},{},bondcfpoltypecommentstocomments,bondcfpoltypecommentstosmartslist)
        anglecfprmstotransferinfo=MapParameterLineToTransferInfo(poltype,anglecfprms,{},{},{},{},anglecfpoltypecommentstocomments,anglecfpoltypecommentstosmartslist)


    else:
        amoebaplusvdwprmstotransferinfo={}
        ctprmstotransferinfo={}
        cpprmstotransferinfo={}
        bondcfprmstotransferinfo={}
        anglecfprmstotransferinfo={}

    tortorprmstotransferinfo=MapParameterLineToTransferInfo(poltype,tortorprms,{},{},{},tortorpoltypeclassestosmartsatomordersext,{},{},defaultvalues=None,keyword='tortors')
    
    bondprmstotransferinfo=MapParameterLineToTransferInfo(poltype,bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips,bondpoltypeclassestosmartsatomordersext,newbondpoltypeclassestocomments,newbondpoltypeclassestosmartslist)
    opbendprmstotransferinfo=MapParameterLineToTransferInfo(poltype,opbendprms,opbendpoltypeclassestoparametersmartsatomorders,opbendpoltypeclassestosmartsatomorders,opbendpoltypeclassestoelementtinkerdescrips,opbendpoltypeclassestosmartsatomordersext,newopbendpoltypeclassestocomments,newopbendpoltypeclassestosmartslist,defaultvalues)
    angleprmstotransferinfo=MapParameterLineToTransferInfo(poltype,angleprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,anglepoltypeclassestosmartsatomordersext,newanglepoltypeclassestocomments,newanglepoltypeclassestosmartslist)
    strbndprmstotransferinfo=MapParameterLineToTransferInfo(poltype,strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,strbndpoltypeclassestosmartsatomordersext,newstrbndpoltypeclassestocomments,newstrbndpoltypeclassestosmartslist)
    torsionprmstotransferinfo=MapParameterLineToTransferInfo(poltype,torsionprms,torsionpoltypeclassestoparametersmartsatomorders,torsionpoltypeclassestosmartsatomorders,torsionpoltypeclassestoelementtinkerdescrips,torsionpoltypeclassestosmartsatomordersext,{},{})
    torsionsmissing=ConvertToPoltypeClasses(poltype,torsionsmissing)
    WriteOutList(poltype,torsionsmissing,poltype.torsionsmissingfilename)
    WriteDictionaryToFile(poltype,torsionkeystringtoparameters,poltype.torsionprmguessfilename)
    WriteOutList(poltype,missingvdwatomindices,poltype.vdwmissingfilename)
    WriteOutList(poltype,tortorsmissing,poltype.tortormissingfilename)
    return bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,torsionsmissing,torsionkeystringtoparameters,missingvdwatomindices,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo,tortorsmissing
