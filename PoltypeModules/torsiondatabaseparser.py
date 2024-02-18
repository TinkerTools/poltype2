from tqdm import tqdm
import warnings
import math
import sys
import os
from packaging import version
import rdkit
from rdkit.Chem import rdFMCS
from openbabel import openbabel
from rdkit.Chem import rdmolfiles
import itertools
import re
from rdkit import Chem
import copy
import symmetry as symm
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdMolDescriptors import CalcNumRings
import numpy as np
import json
import torsionfit as torfit
from rdkit import DataStructs
import torsiongenerator as torgen
from itertools import combinations
import shutil
 

def CheckIfStringIsFloat(string):
    """
    Intent: Check if string is float value
    Input: String
    Output: boolean indicating if string is a float or not
    Referenced By: appendtofile
    Description: 
    1. Assume is not float
    2. Try to convert to float, if works set boolean to True
    3. If fail return boolean
    """
    # STEP 1
    isfloat=False
    # STEP 2
    try:
        float(string)
        isfloat=True
    except:
        # STEP 3
        pass
    return isfloat


def appendtofile(poltype, vf,newname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,tortorprmstotransferinfo):
    """
    Intent: Append the parameters found from database search to key file
    Input: Original key file, new name of updated key file, various dictionaries of parameters to comments about 
    Output: -  
    Referenced By: poltype.py - GenerateParameters 
    Description: 
    1. Iterate over lines of key files 
    2. Detect when find atom block and when after atom block
    3. If polarize or multipole in line, skip this line (only adding valence and vdw parameters for AMOEBA from database here).
    4. Write out AMOEBA vdw parameters 
    5. Write out bond parameters
    6. Write out angle parameters
    7. Write out stretch bend parameters
    8. Write out op-bend parameters
    9. Write out torsion parameters
    10.Write out solute parameters
    11.Write out tor-tor parameters
    12.If user gives key file parameters as input to write in, then write those to key file
    """
    tempname=vf.replace('.key','_temp.key')
    f=open(tempname,'w')
    if poltype.writeoutpolarize==True and poltype.writeoutmultipole==True: 
        temp=open(vf,'r')
        results=temp.readlines()
        temp.close()
        foundatomblock=False
        atomline=False
        wroteout=False
        linestoskip=[]
        # STEP 1
        for theline in results:
            linesplit=theline.split()
            # STEP 2
            if 'atom' in theline:
                atomline=True
                if foundatomblock==False:
                    foundatomblock=True

                        
            else:
                atomline=False
            # STEP 3
            if 'polarize' in theline and poltype.writeoutpolarize==False:
                linestoskip.append(theline)
            if 'multipole' in theline and poltype.writeoutmultipole==False:
                linestoskip.append(theline)
            if len(linesplit)>0:
                if CheckIfStringIsFloat(linesplit[0])==True and poltype.writeoutmultipole==False:
                    linestoskip.append(theline)


            if foundatomblock==True and atomline==False and wroteout==False:
                wroteout=True
                f.write('\n')
                # STEP 4
                if (poltype.writeoutvdw==True) and (poltype.forcefield.upper() == 'AMOEBA'):
                    for line in vdwprmstotransferinfo.keys():
                         f.write(line)
                         f.write('\n')
                
                # STEP 5
                if poltype.writeoutbond==True:
                    for line in bondprmstotransferinfo.keys():
                        f.write(line)
                        f.write('\n')
                # STEP 6 
                if poltype.writeoutangle==True:
                    for line in angleprmstotransferinfo.keys():
                        f.write(line)
                        f.write('\n')
                # STEP 7
                if poltype.writeoutstrbnd==True:
                    for line in strbndprmstotransferinfo.keys():
                        f.write(line)
                        f.write('\n')
                # STEP 8
                if poltype.writeoutopbend==True:
                    for line in opbendprmstotransferinfo.keys():
                        f.write(line)
                        f.write('\n')
                # STEP 9
                if poltype.writeouttorsion==True:
                    for line,transferinfo in torsionprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                # STEP 10
                if poltype.forcefield.upper() == 'AMOEBA':
                    for line in soluteprms:
                        f.write(line)
                        f.write('\n')
                # STEP 11
                for line,transferinfo in tortorprmstotransferinfo.items():
                    if 'tortors' in line:
                        f.write(transferinfo)
                    f.write(line)
                    f.write('\n')
                # STEP 12
                                        
            else:
                if theline not in linestoskip:
                    f.write(theline)
    
    # STEP 13
    if poltype.inputkeyfile!=None:
        if poltype.writeouttorsion==True:
            for line,transferinfo in torsionprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
                f.write('\n')
        temp=open(poltype.inputkeyfile,'r')
        results=temp.readlines()
        temp.close()
        f.write("parameters " + poltype.paramhead + "\n")
        for line in results: # handle case where user gives giant prm file as input
            linesplit=line.split()
            if len(linesplit)>1:
                if linesplit[0]=='torsion':
                    if poltype.writeouttorsion==True:
                        for theline,transferinfo in torsionprmstotransferinfo.items():
                            thelinesplit=theline.split()
                            tors=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
                            thetors=[int(thelinesplit[1]),int(thelinesplit[2]),int(thelinesplit[3]),int(thelinesplit[4])]
                            if tors==thetors or tors==thetors[::-1]:
                                line=theline
                                break
                    
            f.write(line)

    f.close()
    os.rename(tempname,newname)

def ReadSmallMoleculeLib(poltype,filepath):
    """
    Intent: Read amoeba09 SMARTS, atom order -> amoeba09 tinker description file
    Input: Filepath
    Output: Dictionary of SMARTS, atom order -> amoeba09 tinker description
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines of file  
    2. Grab symbol
    3. Grab tinker description
    4. Grab SMARTS string
    5. Grab atom order list
    6. Put items in dictionary 
    """
    temp=open(filepath,'r')
    results=temp.readlines()
    temp.close()
    smartsatomordertoelementtinkerdescrip={}
    # STEP 1
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==0:
            continue
        # STEP 2
        elementsymb=linesplit[0]
        newline=' '.join(linesplit[1:])
        newsplit=newline.split('%')
        # STEP 3
        tinkerdescrip=newsplit[0].lstrip().rstrip()
        # STEP 4
        smarts=newsplit[1].lstrip().rstrip()
        # STEP 5
        atomindices=newsplit[2].lstrip().rstrip()
        atomorderlist=atomindices.split()
        atomorderlist=tuple([int(i) for i in atomorderlist])
        ls=[smarts,atomorderlist]
        newls=[elementsymb,tinkerdescrip]
        # STEP 6
        smartsatomordertoelementtinkerdescrip[tuple(ls)]=tuple(newls)
    return smartsatomordertoelementtinkerdescrip

def GrabParameters(poltype,fname):
    """
    Intent: Grab parameters from parameter file
    Input:
    Output:
    Referenced By: mutation.py GrabBgnToEndPrms
    Description: 
    1. Iterate over lines of key file
    2. If atom keyword, save atom line
    3. If bond keyword, save bond parameter line
    4. If angle keyword, save angle parameter line
    5. If torsion keyword, save torsion parameter line
    6. If strbnd keyword, save stretch bend parameter line
    7. If opbend keyword, save opbend parameter line
    8. If polarize keyword, save polarize parameter line
    9. If vdw keyword, save vdw parameter line
    10. If multipole keyword, save multipole parameter lines 
    """
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
    # STEP 1
    for lineidx in range(len(results)):
        line=results[lineidx]
        # STEP 2
        if 'atom' in line:
            atomdefs.append(line)
        # STEP 3
        elif 'bond' in line:
            bondprms.append(line)
        # STEP 4
        elif 'angle' in line or 'anglep' in line:
            angleprms.append(line)
        # STEP 5
        elif 'torsion' in line:
            torsionprms.append(line) 
        # STEP 6
        elif 'strbnd' in line:
            strbndprms.append(line)
        # STEP 7
        elif 'opbend' in line:
            opbendprms.append(line)
        # STEP 8
        elif 'polarize' in line:
            polarizeprms.append(line)
        # STEP 9
        elif 'vdw' in line: 
            vdwprms.append(line)
        # STEP 10
        elif 'multipole' in line:
            mpolelist=[line,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
            for mpoleline in mpolelist:
                mpoleprms.append(mpoleline)

    return atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms
 
def GrabParametersFromPrmFile(poltype,bondtinkerclassestopoltypeclasses,opbendtinkerclassestopoltypeclasses,opbendtinkerclassestotrigonalcenterbools,angletinkerclassestopoltypeclasses,torsiontinkerclassestopoltypeclasses,poltypetoprmtype,atomtinkerclasstopoltypeclass,typestoframedefforprmfile,fname,skipmultipole=False):
    """
    Intent: Grab parameters from amoeba09 parameter file.
    Input: dictionaries of tinker class -> poltype class numbers for each parameter type
    Output: arrays of parameter lines for each parameter type in amoeba09
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines of parameter file
    2. If bond keyword in line, check if bond classes are in input dictionaries and if so, save parameter line 
    3. If angle keyword in line, check if angle classes are in input dictionaries and if so, save parameter line 
    4. If torsion keyword in line, check if torsion classes are in input dictionaries and if so, save parameter line 
    5. If strbnd keyword in line, check if strbnd classes are in input dictionaries and if so, save parameter line 
    6. If opbend keyword in line, check if opbend classes are in input dictionaries and if so, save parameter line 
    7. If multipole keyword in line, check if multipole classes are in input dictionaries and if so, save parameter line 
    8. If polarize keyword in line, check if polarize classes are in input dictionaries and if so, save parameter line 
    9. If vdw keyword in line, check if vdw classes are in input dictionaries and if so, save parameter line 
    10. If pitor keyword in line, check if middle two classes are in torsion input dictionaries, if so map parameter line to torsion parameters (later on will modify torsion parameters to account for not adding pitor parameters to keyfile)
    """
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
    # STEP 1
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        linesplitall=re.split(r'(\s+)', line)
        if '#' in line:
            continue
        # STEP 2
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
        # STEP 3       
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
        # STEP 4 
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
        # STEP 5
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

               
        # STEP 6 
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
        # STEP 7
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
        # STEP 8 
        elif 'polarize' in line: 
            atomtype=int(linesplit[1])
            if atomtype in poltypetoprmtype.keys():
                prmtype=poltypetoprmtype[atomtype]
                newline=line.replace('\n','')+' '+str(prmtype)+'\n'
                polarizeprms.append(newline)
        # STEP 9 
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
        # STEP 10
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
    """
    Intent: Grab dictionaries of tinker element description -> tinker type and tinker type -> tinker class from amoeba09 parameter file. 
    Input: amoeba09 parameterfile
    Output: dictionaries of tinker element description -> tinker type and tinker type -> tinker class 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines of parameter file, if atom keyword is in line
    2. Grab tinker type
    3. Grab class type
    4. Grab element
    5. Grab tinker description
    6. Put everything into dictionaries 
    """
    temp=open(prmfile,'r')
    results=temp.readlines()
    temp.close()
    elementtinkerdescriptotinkertype={}
    tinkertypetoclass={}
    # STEP 1
    for line in results:
        if 'atom' in line:
            linesplit=line.split()    
            newlinesplit=linesplit[1:-3]
            # STEP 2
            tinkertype=newlinesplit[0]
            # STEP 3
            classtype=newlinesplit[1]
            # STEP 4
            element=newlinesplit[2]
            # STEP 5
            tinkerdescrip=' '.join(newlinesplit[3:])
            ls=[element,tinkerdescrip]
            # STEP 6
            elementtinkerdescriptotinkertype[tuple(ls)]=tinkertype
            tinkertypetoclass[tinkertype]=classtype
    return elementtinkerdescriptotinkertype,tinkertypetoclass



def GrabAtomsForParameters(poltype,mol):
    """
    Intent: Need a list of atoms, bonds, angles and torsions for input molecule to process and parse database later.
    Input: Openbabel mol object
    Output: list of atoms, bonds, angles and torsions
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over atoms append atom indices to atom array
    2. Iterate over neighbor of atoms and if bond (atom,natom) hasnt been found append to array
    3. Iterate over neighbor of neighbor and if first atom and neighbor of neighbor are not same atoms, append angle to array
    4. Iterate over neighbors again to find torsion, make sure no indices are repeating, and append torsion to array
    """
    # now we can define arrays to collect bonds, angles and torsions
    listoftorsionsforprm=[]
    listofbondsforprm=[]
    listofanglesforprm=[]
    listofatomsforprm=[]
    # STEP 1
    for atom in openbabel.OBMolAtomIter(mol):
        atomidx=atom.GetIdx()-1
        listofatomsforprm.append([atomidx])
        neighbs=[natom for natom in openbabel.OBAtomAtomIter(atom)]
        # STEP 2
        for natom in neighbs:
            nidx=natom.GetIdx()-1
            bondset=[nidx,atomidx]
            if bondset not in listofbondsforprm and bondset[::-1] not in listofbondsforprm:
                listofbondsforprm.append(bondset)
            # STEP 3
            nextneighbs=[nextatom for nextatom in openbabel.OBAtomAtomIter(natom)]
            for nextneighb in nextneighbs:
                nextneighbidx=nextneighb.GetIdx()-1
                if nextneighbidx!=atomidx:
                    angleset=[nextneighbidx,nidx,atomidx]
                    if angleset not in listofanglesforprm and angleset[::-1] not in listofanglesforprm:
                        listofanglesforprm.append(angleset)
                    nextnextneighbs=[nextnextatom for nextnextatom in openbabel.OBAtomAtomIter(nextneighb)]
                    # STEP 4
                    for nextnextatom in nextnextneighbs:
                        nextnextatomidx=nextnextatom.GetIdx()-1
                        if nextnextatomidx!=nidx and nextnextatomidx!=atomidx:                
                            torsionset=[nextnextatomidx,nextneighbidx,nidx,atomidx]
                            if torsionset not in listoftorsionsforprm and torsionset[::-1] not in listoftorsionsforprm:
                                listoftorsionsforprm.append(torsionset)
    return listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm



def ExtendByNeighbors(poltype,ls):    
    """
    Intent: Add neighbors to input array, when wanting to check if neighbors exist in SMARTS match
    Input: Array of atom indices
    Output: Array of atom indices plus neighbors
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. If input array is of length 1 (vdw atom matches) want neighbor of neighbors also (so match neighbor of neighbor on hydrogen for example to get enough information if match is correct)
    2. Iterate over atom indices in input array
    3. Iterate over neighboring atom indices of atom index and append to new array
    4. If input array is of length 1 also iterate over neighbor of neighbors and append indices to new array 
    """
    # STEP 1
    extendneighbofneighb=False
    if len(ls)==1:
        extendneighbofneighb=True 
    newls=[]
    # STEP 2
    for atomidx in ls:
        if atomidx not in newls:
            newls.append(atomidx)
        atom=poltype.rdkitmol.GetAtomWithIdx(atomidx)
        # STEP 3
        for natom in atom.GetNeighbors():
            natomidx=natom.GetIdx()
            if natomidx not in newls:
                newls.append(natomidx)     
            # STEP 4 
            if extendneighbofneighb==True and len(atom.GetNeighbors())==1:
                for nnatom in natom.GetNeighbors():
                    nnatomidx=nnatom.GetIdx()
                    if nnatomidx not in newls:
                        newls.append(nnatomidx)


    return newls
             


def GenerateFragmentSMARTS(poltype,ls):
    """
    Intent: When matching input molecule to database need a way to take part of input molecule, convert to SMARTS to attempt SMARTS matching. This functions generate SMARTS from input array of atom indices in molecule. 
    Input: Array of atom indices from input molecule
    Output: SMART string, second SMARTS that is modified so that any atoms that were aromatic in original mol are aromatic in SMARTS
    Referenced By: GenerateFragmentSMARTSList
    Description:
    1. Create new mol object
    2. Iterate over input array of atom indices
    3. Grab the atom object from input molecule and add atom to newly created mol object
    4. Keep track of old atom index and new atom index in dicttoinary
    5. Keep track if atom is aromatic atom or not in original molecule
    6. Iterate over bonds
    7. If one atom is in new molecule and but the other atom is not, then keep track of this bond as a cut bond. If both are in new molecule then add bond to new molecule.
    8. Generate SMARTS from new mol  
    9. Iterate back over atoms in new mol and if they were aromatic atoms in original molecule, make them aromatic in new mol
    10.Iterate back over bonds in new mol and if they were aromatic bonds in original molecule, make them aromatic in new mol 
    11.Generate new SMARTS from modified mol (containing aromaticity information)
    """
    # STEP 1
    newmol = Chem.Mol()
    mw = Chem.RWMol(newmol)
    # need to treat ring bonds as aromatic since all transferred parameters from amoeba09 are aromatic rings
    oldindextonewindex={}
    aromaticindices=[]
    # STEP 2
    for i,idx in enumerate(ls):
        # STEP 3
        oldatom=poltype.rdkitmol.GetAtomWithIdx(idx)
        mw.AddAtom(oldatom)
        # STEP 4
        oldindextonewindex[idx]=i
        oldatombabel=poltype.mol.GetAtom(idx+1)
        # STEP 5
        isaro=oldatombabel.IsAromatic()
        isinring=oldatombabel.IsInRing()
        hyb=oldatombabel.GetHyb()
        if isaro==True and isinring==True and hyb==2:
            aromaticindices.append(i)

    atomswithcutbonds=[]
    aromaticbonds=[]
    bonditer=poltype.rdkitmol.GetBonds()
    # STEP 6
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        babelbond=poltype.mol.GetBond(oendidx+1,obgnidx+1)
        isinring=babelbond.IsInRing()
        # STEP 7
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
        if isinring==True and endidx in aromaticindices and bgnidx in aromaticindices:
            aromaticbonds.append([endidx,bgnidx]) 
        bondorder=bond.GetBondType()
        mw.AddBond(bgnidx,endidx,bondorder)
    # STEP 8
    smarts=rdmolfiles.MolToSmarts(mw)
    # STEP 9
    for atom in mw.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx in aromaticindices:
            atom.SetIsAromatic(True)
    # STEP 10
    for bond in mw.GetBonds():
        endidx = bond.GetEndAtomIdx()
        bgnidx = bond.GetBeginAtomIdx()
        temp=[endidx,bgnidx]
        if temp in aromaticbonds or temp[::-1] in aromaticbonds:
            bond.SetIsAromatic(True)
    # STEP 11
    smartsfortransfer=rdmolfiles.MolToSmarts(mw)
    
    return smarts,smartsfortransfer


def GenerateFragmentSMARTSList(poltype,ls):
    """
    Intent: When wanting to find matches to amoeba09 database, need a way to generate many SMARTS of varying lengths centered around inpyt atom indices from input molecule. This is sort of like maximum common substructure search but only centered on atom desired so need custom code. 
    Input: List of atom indices from input molecule want to generate SMARTS strings for
    Output: List of SMARTS strings containing matches to part of molecule with input atom indices
    Referenced By: MatchAtomIndicesSMARTSToParameterSMARTS
    Description: 
    1. Generate all possible fragment indices containing input atom indices up to length 14 (4 SP3 atoms for torsion have max 14 atoms)
    2. Iterate over each list of atom indices just generated
    3. Generate a SMARTS string for each one and append to array 
    """
    fragsmartslist=[]
    # STEP 1
    atomindiceslist=GenerateAllPossibleFragmentIndices(poltype,ls,poltype.rdkitmol,14)
    # STEP 2
    for thels in atomindiceslist:
        # STEP 3
        smarts,smartsfortransfer=GenerateFragmentSMARTS(poltype,thels)
        if smartsfortransfer not in fragsmartslist:
            fragsmartslist.append(smartsfortransfer)


    return fragsmartslist


def FindHowManyBondTypes(poltype,trymol,atomidx):
    """
    Intent: Want to define a "type" that does not take into account global symmetry but only surrounding envioronment. Define it with atomicnumber, number of single,aromatic, double and triple bonds.
    Input: Mol object, atom index that want a type for
    Output: The atom type (array of atomic number, number of bond types for each bond type for that atom)
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. Grab the atom of interest from input atom index and mol object
    2. Iterate over neighbors of atom of interest
    3. Grab bond type and increase count of respective bond type (single, aromatic, double, triple)
    4. Construct atom type using information of atomic number as well as how many single, aromatic, double and triple bonds.
    """
    numsinglebonds=0
    numaromaticbonds=0
    numdoublebonds=0
    numtriplebonds=0
    # STEP 1
    atom=trymol.GetAtomWithIdx(atomidx)
    # STEP 2
    for natom in atom.GetNeighbors():
        natomidx=natom.GetIdx()
        bond=trymol.GetBondBetweenAtoms(atomidx,natomidx)
        # STEP 3
        bondtype=bond.GetBondTypeAsDouble() 
        if bondtype==1:
            numsinglebonds+=1
        elif bondtype==1.5:
            numaromaticbonds+=1
        elif bondtype==2:
            numdoublebonds+=1
        elif bondtype==3:
            numtriplebonds+=1
    atomicnum=atom.GetAtomicNum()
    # STEP 4
    atomtype=tuple([atomicnum,numsinglebonds,numaromaticbonds,numdoublebonds,numtriplebonds])

    return atomtype



def MatchAllPossibleSMARTSToParameterSMARTS(poltype,parametersmartslist,parametersmartstosmartslist,ls,mol,parametersmartstordkitmol,parametersmartstosmartsmatchingtoindices,parametersmartstomolsmatchingtoindices,parametersmartstordkitmolmatchingindices,parametersmartstoprmmolmatchingindices):
    """
    Intent: For a given input array of atom indices (one for atom, two for bond, three for angle etc) then try to match SMARTS containing those indices to amoeba09 database in order to find parameters suitable to the given environment.  
    Input: List of SMARTS strings from amoeba09, dictionary mapping amoeba09 SMARTS to final matched SMARTS from input molecule, input list of atom indices from molecule trying to find SMARTS matches to amoeba09 for, mol object, dictionary of amoeba09 SMARTS -> mol object, dictionary of amoeba09 SMARTS ->  matching indices of input molecule,  dictionary of amoeba09 SMARTS ->  matching indices of SMARTS generated from input molecule
    Output: dictionary mapping amoeba09 SMARTS to final matched SMARTS from input molecule, dictionary mapping amoeba09 SMARTS to boolean specifying if the SMARTS matches all neighbors of input list of atom indices in input molecule
    Referenced By: MatchAtomIndicesSMARTSToParameterSMARTS
    Description: 
    1. Iterate over amoeba09 SMARTS
    2. Grab amoeba09 SMARTS
    3. Grab mol object for amoeba09 SMARTS
    4. Grab list of SMARTS that matched to amoeba09 SMARTS
    5. Grab list of mols from SMARTS that matched to amoeba09 SMARTS
    6. Iterate over list of SMARTS that matched to amoeba09 SMARTS
    7. Grab SMARTS that matched to amoeba09 SMARTS
    8. Grab corresponding mol object for that SMARTS 
    9. Assume the SMARTS match is a "good" match for now
    10.Grab the first match of input molecule indices that matched from SMARTS 
    11.If the number of atoms in match is greater than or equal to number of atoms in input array of atom indices
    12.Iterate over all matches
    13.If there is an atom index in ls (input list of atom indices wish to find amoeba09 parameters for) that does not exist in match, then this is not a "good" match and skip that one. 
    14.Check if the matches are consecutively connected in the input molecule (otherwise SMARTS match wouldnt make sense)
    15.Define new array that includes neighbors of (ls) input array of atom indices. Determine if the match contains all the neighbors or not.
    16. Sanity check to ensure the number of rings in amoeba09 SMARTS doesnt exceed number of rings in SMARTS match from input mol (redundant now?)
    17. Grab the hybridizations of each atom in input molecule that matched to SMARTS, later will use to filter if match has same hybridzation as amoeba09 SMARTS
    18. Count the number of times each "type" defined by FindHowManyBondTypes, occurs for array of input atom indices + neighbors. Will later use as filtration step for potential SMARTS matches to amoeba09 parameter SMARTS
    19. Match generated SMARTS to input molecule
    20. Iterate over matches and if all indices in input array exist in match and match indices are consecutive
    21. Generate dictionary of input molecule index to SMARTS index (the atom matched in input molecule to index of atom in SMARTS)
    22. Grab indices that map current SMARTS to current amoeba09 SMARTS
    23. SMARTS matching from current SMARTS->amoeba09 SMARTS returns many matches, try and pick one that includes most neighbors from input array of atom indices from input molecule
    24. Skip any matches that dont have correct hybridization between input molecule atoms and atoms matched with amoeba09 SMARTS
    25. Compute various "scores" for filtering out possible SMARTS that match to amoeba09 SMARTS to find the "best" match. The first score is how many atoms from input array + neighbors are included in match (want to maximize this), the second score is the difference between number of local "types" in input molecule vs "types" in amoeba09 SMARTS (want to minimize this, a difference of 0 means molecule must be the same), the third score is the difference between number of each atomic number occuring in input molecule vs amoeba09 SMARTS (want to minimize).
    26. Apply filters sequentially.
    27. If no good matches are found in case of vdw matches, just try and match element only and keep as match (better than wild card match). 
    """
    smartsmcstomol={}
    prmsmartstomcssmarts={}
    parametersmartstoscore={}
    parametersmartstonumcommon={}
    parametersmartstootherscore={}
    parametersmartstothefinalscore={}
    parametersmartstofoundallneighbs={}
    newls=copy.deepcopy(ls)
    # STEP 1
    for parametersmartsidx in tqdm(range(len(parametersmartslist)),desc='amoeba09 search for '+str(ls)):
        # STEP 2
        parametersmarts=parametersmartslist[parametersmartsidx]
        # STEP 3
        prmmol=parametersmartstordkitmol[parametersmarts]
        # STEP 4
        smartsmatchingtoindices=parametersmartstosmartsmatchingtoindices[parametersmarts]
        # STEP 5
        molsmatchingtoindices=parametersmartstomolsmatchingtoindices[parametersmarts]
        prmsmartsatomnum=prmmol.GetNumAtoms()
        thesmartstonumneighbs={}
        thesmartstotypescore={}
        thesmartstoelementscore={}
        thesmartstothemol={}
        # STEP 6
        for theidx in range(len(smartsmatchingtoindices)):
            # STEP 7
            thesmarts=smartsmatchingtoindices[theidx]
            # STEP 8
            themol=molsmatchingtoindices[theidx]
            thesmartstothemol[thesmarts]=themol
            # STEP 9
            diditmatchprmmol=True
            diditmatch=True
            if diditmatch==True and diditmatchprmmol==True:
                matches=parametersmartstordkitmolmatchingindices[parametersmarts][theidx]
                # STEP 10
                firstmatch=matches[0]
                # STEP 11 
                if len(firstmatch)>=len(ls):
                    # STEP 12
                    for match in matches:
                        goodmatch=True
                        matchidxs=[]
                        # STEP 13
                        for idx in ls:
                            if idx not in match:
                                goodmatch=False
                                break
                            else:
                                matchidx=match.index(idx)
                                matchidxs.append(matchidx)
                        # STEP 14
                        if len(ls)>1 and goodmatch==True:
                            goodmatch=CheckIfConsecutivelyConnected(poltype,matchidxs,themol)
                        # STEP 15
                        newls=ExtendByNeighbors(poltype,ls) 
                        score=0 
                        allneighbsin=True
                        matchidxs=[]
                        for idx in newls:
                            if idx in match:
                                score+=1
                                matchidx=match.index(idx)
                                matchidxs.append(matchidx)
                            else:
                                allneighbsin=False
                        if goodmatch==True:
                             
                            break
                    # STEP 16 
                    prmmolnumrings=CountRingsInSMARTS(poltype,parametersmarts)
                    molnumrings=CalcNumRings(poltype.rdkitmol)
                    if prmmolnumrings>molnumrings:
                        goodmatch=False


                    if goodmatch==True:
                        # STEP 17
                        rdkitatoms=[poltype.rdkitmol.GetAtomWithIdx(r) for r in firstmatch]
                        rdkithybs=[r.GetHybridization() for r in rdkitatoms]
                        rdkitmatch=copy.deepcopy(firstmatch)
                        idx=ls[0]
                            


                        # STEP 18
                        rdkitatomictypetotype={}
                        rdkitatomicnumtonum={}
                        for atom in poltype.rdkitmol.GetAtoms():
                            atomicnum=atom.GetAtomicNum()
                            atomidx=atom.GetIdx()
                            if atomicnum not in rdkitatomicnumtonum.keys():
                                rdkitatomicnumtonum[atomicnum]=0
                            rdkitatomicnumtonum[atomicnum]+=1
                            if atomidx in newls:
                                atomtype=FindHowManyBondTypes(poltype,poltype.rdkitmol,atomidx)
                                
                                if atomtype not in rdkitatomictypetotype.keys():
                                    rdkitatomictypetotype[atomtype]=0
                                rdkitatomictypetotype[atomtype]+=1
                        # STEP 19 
                        sp=openbabel.OBSmartsPattern()
                        openbabel.OBSmartsPattern.Init(sp,thesmarts)
                        diditmatch=sp.Match(poltype.mol)
                        babelmatches=sp.GetMapList()
                        newbabelmatches=[]
                        for mtch in babelmatches:
                            newmtch=[i-1 for i in mtch]
                            newbabelmatches.append(newmtch)
                        for newmtch in newbabelmatches:
                            if newmtch not in matches:
                                matches.append(newmtch)
                        # STEP 20
                        for match in matches:
                            validmatch=CheckMatch(poltype,match,ls,thesmarts,themol)
                            if validmatch==True:
                               indices=list(range(len(match)))
                               smartsindextomoleculeindex=dict(zip(indices,match)) 
                               moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}


                        # STEP 21
                        smartindices=[moleculeindextosmartsindex[i] for i in ls]
                        # STEP 22
                        prmmatches=parametersmartstoprmmolmatchingindices[parametersmarts][theidx]
                        # STEP 23 
                        firstmatch=TryAndPickMatchWithNeighbors(poltype,prmmatches,smartindices,themol)
                        # STEP 24
                        indices=list(range(len(firstmatch)))
                        smartsindextoparametersmartsindex=dict(zip(indices,firstmatch)) 
                        prmsmartsindices=[smartsindextoparametersmartsindex[i] for i in smartindices]
                        prmatoms=[prmmol.GetAtomWithIdx(r) for r in firstmatch]
                        prmhybs=[r.GetHybridization() for r in prmatoms]
                        allhybssame=True
                        for r in range(len(rdkithybs)):
                            rdkithyb=rdkithybs[r]
                            prmhyb=prmhybs[r]
                            if rdkithyb!=prmhyb:
                                allhybssame=False
                        if allhybssame==False:
                            continue
                        # STEP 25 
                        prmsmartsatomicnumtonum={}
                        for atomicnum,num in rdkitatomicnumtonum.items():
                            if atomicnum not in prmsmartsatomicnumtonum.keys():
                                prmsmartsatomicnumtonum[atomicnum]=0
                        prmsmartsatomictypetotype={}
                        for atomicnum,num in rdkitatomictypetotype.items():
                            prmsmartsatomictypetotype[atomicnum]=0
                        extra=prmsmartsindices[:]
                        for atom in prmmol.GetAtoms():
                            atomicnum=atom.GetAtomicNum()
                            atomidx=atom.GetIdx()
                            if atomicnum not in prmsmartsatomicnumtonum.keys():
                                prmsmartsatomicnumtonum[atomicnum]=0
                            prmsmartsatomicnumtonum[atomicnum]+=1
                            if atomidx in prmsmartsindices:
                                atomtype=FindHowManyBondTypes(poltype,prmmol,atomidx)

                                if atomtype not in prmsmartsatomictypetotype.keys():
                                    prmsmartsatomictypetotype[atomtype]=0
                                prmsmartsatomictypetotype[atomtype]+=1
                                for natom in atom.GetNeighbors():
                                    natomidx=natom.GetIdx()
                                    atomicnum=natom.GetAtomicNum()

                                    if natomidx not in extra:
                                        atomtype=FindHowManyBondTypes(poltype,prmmol,natomidx)
                                        if atomtype not in prmsmartsatomictypetotype.keys():
                                            prmsmartsatomictypetotype[atomtype]=0
                                        prmsmartsatomictypetotype[atomtype]+=1
                                        extra.append(natomidx)
                                    if len(ls)==1 and len(atom.GetNeighbors())==1:
                                        for nnatom in natom.GetNeighbors():
                                            nnatomidx=nnatom.GetIdx()
                                            atomicnum=nnatom.GetAtomicNum()

                                            if nnatomidx not in extra:
                                                atomtype=FindHowManyBondTypes(poltype,prmmol,nnatomidx)
                                                if atomtype not in prmsmartsatomictypetotype.keys():
                                                    prmsmartsatomictypetotype[atomtype]=0
                                                prmsmartsatomictypetotype[atomtype]+=1
                                                extra.append(natomidx)




                        otherscore=0
                        for atomicnum,num in rdkitatomictypetotype.items():
                            prmnum=prmsmartsatomictypetotype[atomicnum]
                            diff=np.abs(prmnum-num)

                            otherscore+=diff

                        globalscore=0
                        for atomicnum,num in rdkitatomicnumtonum.items():
                            prmnum=prmsmartsatomicnumtonum[atomicnum]
                            diff=np.abs(prmnum-num)

                            globalscore+=diff

                        thesmartstonumneighbs[thesmarts]=score
                        thesmartstotypescore[thesmarts]=otherscore
                        thesmartstoelementscore[thesmarts]=globalscore
    

        if len(thesmartstonumneighbs.keys())>0:
            maxscore=max(thesmartstonumneighbs.values())            
            for thesmarts,score in thesmartstonumneighbs.items():
                themol=thesmartstothemol[thesmarts]
                otherscore=thesmartstotypescore[thesmarts]
                if score==maxscore:

                    prmsmartstomcssmarts[parametersmarts]=thesmarts
                    parametersmartstoscore[parametersmarts]=score
                    parametersmartstootherscore[parametersmarts]=otherscore
                    smartsmcstomol[thesmarts]=themol
                    parametersmartstofoundallneighbs[parametersmarts]=allneighbsin
                    parametersmartstothefinalscore[parametersmarts]=thesmartstoelementscore[thesmarts]
                    break
    foundmin=False
    # STEP 26
    parametersmartstofinalscore={}
    parametersmartstolastscore={}
    if len(parametersmartstoscore.keys())>0:
        minscore=max(parametersmartstoscore.values())
        for parametersmarts in parametersmartstoscore.keys():
            score=parametersmartstoscore[parametersmarts]
            otherscore=parametersmartstootherscore[parametersmarts]
            smartsmcs=prmsmartstomcssmarts[parametersmarts]
            lastscore=parametersmartstothefinalscore[parametersmarts]
            if score==minscore:
                parametersmartstolastscore[parametersmarts]=lastscore          
                parametersmartstofinalscore[parametersmarts]=otherscore
        minscore=min(parametersmartstofinalscore.values())
        finalprmsmartstoscore={}
        for parametersmarts in parametersmartstofinalscore.keys():
            score=parametersmartstofinalscore[parametersmarts]
            smartsmcs=prmsmartstomcssmarts[parametersmarts]
            lastscore=parametersmartstolastscore[parametersmarts]
            if score==minscore:
                finalprmsmartstoscore[parametersmarts]=lastscore
                
        minscore=min(finalprmsmartstoscore.values())
        for parametersmarts,score in finalprmsmartstoscore.items():
            if score==minscore:
                smartsmcs=prmsmartstomcssmarts[parametersmarts] 
                mcsmol=smartsmcstomol[smartsmcs]
                foundmin=True
                break

        if foundmin==True:
            smartls=[smartsmcs,smartsmcs]
            if parametersmarts not in parametersmartstosmartslist.keys():
                parametersmartstosmartslist[parametersmarts]=smartls
    else:
        # STEP 27
        if len(ls)==1:
            atom=poltype.rdkitmol.GetAtomWithIdx(ls[0])
            atomicnum=atom.GetAtomicNum()
            string='[#'+str(atomicnum)+']'
            othermol=Chem.MolFromSmarts(string)
            for parametersmarts in parametersmartslist:
                prmmol=Chem.MolFromSmarts(parametersmarts)
                mols = [othermol,prmmol]
                diditmatch=mols[1].HasSubstructMatch(mols[0])
                if diditmatch==True:
                    matches=mols[1].GetSubstructMatches(mols[0])
                    firstmatch=matches[0]
                    smartls=[string,string]
                    if parametersmarts not in parametersmartstosmartslist.keys():
                        parametersmartstosmartslist[parametersmarts]=smartls
                        parametersmartstofoundallneighbs[parametersmarts]=False
                    break

         

    return parametersmartstosmartslist,parametersmartstofoundallneighbs


def GenerateAllPossibleFragmentIndices(poltype,ls,rdkitmol,maxatomsize):
    """
    Intent: When generating SMARTS from input molecule to match to amoeba09 SMARTS, need to generate many atomic indices first for each possible SMARTS.
    Input: Input array of atom indices from input molecule, mol object, maximum number of atoms to search out in molecule.
    Output: List of list of atomic indices each containing input array of atom indices
    Referenced By: GenerateFragmentSMARTSList
    Description: 
    1. Start with list of current atom indices
    2. Define a second list that will constantly be updated with more atoms and stop when the lists are the same. After first iteration set oldindexlist=indexlist, if at end of iteration, no new neighbors are added, then will exit while loop.
    3. If length of list is greater than maxatomsize or if searched further than two neighbors away, then quit
    4. Grab neighboring indices and for each neighbor, generate a list of original atom indices+some of neigbors, append to list of atom indices to generate SMARTS from
    5. Add neighboring indices to indexlist 
    """
    # STEP 1
    atomindiceslist=[copy.deepcopy(ls)]
    oldindexlist=copy.deepcopy(ls)
    indexlist=[]
    count=0
    neighbcount=0
    # STEP 2
    while set(oldindexlist)!=set(indexlist):
        if count!=0:
            oldindexlist=copy.deepcopy(indexlist)
        # STEP 3
        if len(oldindexlist)>maxatomsize:
            break
        if neighbcount>=2:
            break
        # STEP 4
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
        # STEP 5
        newindexlist=copy.deepcopy(oldindexlist)
        for index in neighborindexes:
            if index not in newindexlist:
                newindexlist.append(index)
        indexlist=newindexlist
  
    return atomindiceslist


def GrabNeighboringIndexes(poltype,indexlist,rdkitmol):
    """
    Intent: Grab neighboring indices for when generating all possible SMARTS to match to a trial amoeba09 SMARTS.
    Input: Array of atom indices, mol object
    Output: Array of neighboring atom indices
    Referenced By: GenerateAllPossibleFragmentIndices
    Description: 
    1. Iterate over atom indices in input list
    2. Grab neighbors of each atom
    3. Append to array of neighbors 
    """
    neighborindexes=[]
    # STEP 1
    for index in indexlist:
        atom=rdkitmol.GetAtomWithIdx(index)
        # STEP 2
        for natom in atom.GetNeighbors():
            natomidx=natom.GetIdx()
            # STEP 3
            if natomidx not in neighborindexes and natomidx not in indexlist:
                neighborindexes.append(natomidx)
    return neighborindexes


def CountRingsInSMARTS(poltype,parametersmarts):
    """
    Intent: Sanity check to ensure number of rings in amoeba09 SMARTS, doesnt exceed number of rings in input molecule 
    Input: SMARTS strings
    Output: Number of rings
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. Iterate over SMARTS characters
    2. Search for when after ] bracket (after atom) and if there is a digit after that (specifies ring number), append to array of ring numbers
    3. Return length of ring numbers array 
    """
    closedbrack=True
    numbers=[]
    # STEP 1
    for e in parametersmarts:
        # STEP 2
        if e=='[':
            closedbrack=False
        elif e==']':
            closedbrack=True
        elif e.isdigit():
            if closedbrack==True:
                if e not in numbers:
                    numbers.append(e)
    # STEP 3
    rings=len(numbers)
    return rings

def CheckIfConsecutivelyConnected(poltype,matchidxs,mcsmol):
    """
    Intent: Check if indices in match are consecutively connected (needed for adequate SMARTS match to amoeba09 SMARTS).
    Input: Array of indices in SMARTS that match to input molecule, mol from SMARTS 
    Output: Boolean specifying if all indices are consecutive or not.
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. Iterate over atom indices (except last one)
    2. Grab atom from SMARTS mol that corresponds to index
    3. Grab next index
    4. Grab neighboring indices of current index
    5. If the next atom index in sequence does not exist in neighboring indices of current index, then the match is not consecutively connected. 
    """
    goodmatch=True
    # STEP 1
    for i in range(len(matchidxs)-1):
        matchidx=matchidxs[i]
        # STEP 2
        atom=mcsmol.GetAtomWithIdx(matchidx)
        # STEP 3
        nextmatchidx=matchidxs[i+1]
        natoms=atom.GetNeighbors()
        # STEP 4
        natomidxs=[a.GetIdx() for a in natoms]
        # STEP 5
        if nextmatchidx not in natomidxs:
            goodmatch=False


    return goodmatch





def MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listforprm,parametersmartslist,mol,parametersmartstordkitmol):
    """
    Intent: For each list of atom indices (all atoms, all bonds, all angles etc), find the "best" amoeba09 SMARTS match. 
    Input: List of list of atom indices, list of all amoeba09 SMARTS, mol object, dictionary of amoeba09 SMARTS to corresponding mol object.
    Output: Dictionary of atom indices (vdw, bond, angle, torsion) -> best amoeba09 SMARTS match, Dictionary of atom indices (vdw, bond, angle, torsion) -> boolean if amoeba09 SMARTS matches to all neighbors of atom indices 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each list of atom indices
    2. For each list of atom indices, generate a list of potential SMARTS to match to any amoeba09 SMARTS
    3. Save results in list for later 
    4. Iterate over list of amoeba09 SMARTS
    5. Iterate over list of all possible SMARTS to match to current amoeba09 SMARTS
    6. If the current SMARTS matches to both input molecule and to amoeba09 SMARTS molecule, then will save for later filtering
    7. Save SMARTS match results in dictionaries for later processing and filtering
    8. Iterate over each list of atom indices again
    9. Call MatchAllPossibleSMARTSToParameterSMARTS using current list of atom indices and the saved list of amoeba09 SMARTS and potential SMARTS matches from input molecule. This will output the "best" SMARTS match for each amoeba09 SMARTS.
    10. Choose the amoeba09 SMARTS, that has the largest SMARTS matching to it (most information in input molecule environment)
    11. If no match was found, assign a wildcard SMARTS (worst case)  
    """
    listforprmtoparametersmarts={}
    listforprmtosmarts={}
    listforprmtomatchallneighbs={}
    allfragsmartslist=[]
    # STEP 1
    for ls in listforprm:
        # STEP 2
        fragsmartslist=GenerateFragmentSMARTSList(poltype,ls)
        # STEP 3
        for fragsmarts in fragsmartslist:
            if fragsmarts not in allfragsmartslist:
                allfragsmartslist.append(fragsmarts)

    parametersmartstosmartsmatchingtoindices={}
    parametersmartstomolsmatchingtoindices={}
    parametersmartstordkitmolmatchingindices={}
    parametersmartstoprmmolmatchingindices={}
    # STEP 4
    for parametersmarts in parametersmartslist:
        prmmol=parametersmartstordkitmol[parametersmarts]
        finalfragsmartslist=[]
        smartsmatchingtoindices=[]
        molsmatchingtoindices=[]
        rdkitmolmatchingindices=[]
        prmmolmatchingindices=[]
        # STEP 5
        for fragsmarts in allfragsmartslist:
            fragmol=Chem.MolFromSmarts(fragsmarts)
            diditmatch=poltype.rdkitmol.HasSubstructMatch(fragmol)
            # STEP 6
            if diditmatch==True:
                diditmatch=prmmol.HasSubstructMatch(fragmol)
                if diditmatch==True:
                    matches=list(poltype.rdkitmol.GetSubstructMatches(fragmol,maxMatches=10000))
                    prmmatches=prmmol.GetSubstructMatches(fragmol)
                    molsmatchingtoindices.append(fragmol)
                    smartsmatchingtoindices.append(fragsmarts)
                    rdkitmolmatchingindices.append(matches)
                    prmmolmatchingindices.append(prmmatches)
        # STEP 7
        parametersmartstosmartsmatchingtoindices[parametersmarts]=smartsmatchingtoindices
        parametersmartstomolsmatchingtoindices[parametersmarts]=molsmatchingtoindices
        parametersmartstordkitmolmatchingindices[parametersmarts]=rdkitmolmatchingindices
        parametersmartstoprmmolmatchingindices[parametersmarts]=prmmolmatchingindices
    # STEP 8
    for ls in listforprm:
        # STEP 9 
        parametersmartstomatchlen={}
        parametersmartstosmartslist={}
        parametersmartstosmartslist,parametersmartstofoundallneighbs=MatchAllPossibleSMARTSToParameterSMARTS(poltype,parametersmartslist,parametersmartstosmartslist,ls,mol,parametersmartstordkitmol,parametersmartstosmartsmatchingtoindices,parametersmartstomolsmatchingtoindices,parametersmartstordkitmolmatchingindices,parametersmartstoprmmolmatchingindices)
        # STEP 10
        if len(parametersmartstosmartslist.keys())!=0:
            parametersmartstosmartslen={}
            for prmsmarts,newls in parametersmartstosmartslist.items():
                smarts=newls[0]
                smartslen=len(smarts)
                parametersmartstosmartslen[prmsmarts]=smartslen
            
            valuelist=list(parametersmartstosmartslen.values())
            maxvalue=max(valuelist)
            for prmsmarts,smartslen in parametersmartstosmartslen.items():
                if smartslen==maxvalue:
                    maxprmsmarts=prmsmarts
                    break 
            maxsmartsls=parametersmartstosmartslist[maxprmsmarts]
            matchallneighbs=parametersmartstofoundallneighbs[maxprmsmarts]
        else:
            # STEP 11
            matchallneighbs=False
            wild='[*]'
            smarts=''
            for i in range(len(ls)):
                smarts+=wild+'~'
            smarts=smarts[:-1]
            maxsmartsls=[smarts,smarts]
            maxprmsmarts='[#6](-[H])(-[H])(-[H])-[#6](-[H])(-[H])(-[H])'
        
        listforprmtoparametersmarts[tuple(ls)]=maxprmsmarts
        listforprmtosmarts[tuple(ls)]=maxsmartsls
        listforprmtomatchallneighbs[tuple(ls)]=matchallneighbs
    return listforprmtoparametersmarts,listforprmtosmarts,listforprmtomatchallneighbs



def CheckMatch(poltype,match,atomindices,smarts,substructure):
    """
    Intent: Want to ensure that when matching SMARTS to input molecule, the atom indices list (vdw, bond, angle etc) are all consecutively matched in the SMARTS match. 
    Input: Array of indices in input molecule that correspond to the SMARTS match, array of indices of interest from input molecule, SMARTS string, mol for SMARTS string 
    Output: Boolean specifying if match is consecutive or not
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS , GenerateAtomIndexToAtomTypeAndClassForAtomList
    Description: 
    1. If all indices in molecule of interest are not in match array, then its not a valid match to begin with
    2. Assume indices are consecutive
    3. For each index in indices of interest, find location of that index in match array
    4. Save each location found in new array
    5. Now using the mol and location indices array, check if each atom that should be consecutive is consecutive in that match array, this will return boolean if consecutive or not
    """
    # STEP 1
    allin=True
    for idx in atomindices:
        if idx not in match:
            allin=False
    # STEP 2
    validmatch=True
    if allin==True:
        indices=[]
        # STEP 3
        for idx in atomindices:
            matchidx=match.index(idx)
            # STEP 4
            indices.append(matchidx)
        # STEP 5
        checkconsec=CheckConsecutiveTorsion(poltype,indices,substructure)
        if checkconsec==False:
            validmatch=False
    else:
        validmatch=False
    return validmatch 
       
def CheckConsecutiveTorsion(poltype,indices,substructure):
    """
    Intent: Check for each parameter type, bond, angle, torsion, if the match to SMARTS is consecutive. 
    Input: Array of atom indices of interest in molecule sorted by location found in SMARTS match, mol for SMARTS
    Output: Boolean if consecutive or not
    Referenced By: CheckMatch
    Description: 
    1. For bond, just check if there exists bond between indices or not
    2. For angle, just check if there exists bonds between a,b and b,c or not
    3. For torsion, just check if there exists bonds between a,b and b,c and c,d or not
    """
    consec=True
    # STEP 1
    if len(indices)==2:
        a,b=indices[:]
        bond=substructure.GetBondBetweenAtoms(a,b)
        if bond==None:
            consec=False
    # STEP 2
    elif len(indices)==3:
        a,b,c=indices[:]
        bond=substructure.GetBondBetweenAtoms(a,b)
        if bond==None:
            consec=False
        bond=substructure.GetBondBetweenAtoms(b,c)
        if bond==None:
            consec=False
    # STEP 3
    elif len(indices)==4: 
        a,b,c,d=indices[:]
        bond=substructure.GetBondBetweenAtoms(a,b)
        if bond==None:
            consec=False
        bond=substructure.GetBondBetweenAtoms(b,c)
        if bond==None:
            consec=False
        bond=substructure.GetBondBetweenAtoms(c,d)
        if bond==None:
            consec=False


    return consec


def GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol):
    """
    Intent: For each list of atom indices with the "best" amoeba09 SMARTS match, determine the tinker type and class number and store in dictionaries for later use.
    Input: Dictionary of atom indices -> "best" amoeba09 SMARTS ,  Dictionary of atom indices -> "best" SMARTS that matches to amoeba09 SMARTS and input molecule, dictionary of amoeba09 SMARTS, element -> tinker type description, dictionary of element tinker description -> tinker type, dictionary of tinker type -> tinker class, input mol object
    Output: Dictionary of atom indices -> tinker types, Dictionary of atom indices -> tinker classes, Dictionary of atom indices -> amoeba09 SMARTS + atom order, Dictionary of atom indices -> tinker description + element, Dictionary of atom indices -> SMARTS + atom order
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of list of atom indices to best amoeba09 SMARTS
    2. Grab corresponding best SMARTS match to input molecule 
    3. Generate mol object from SMARTS
    4. Generate matches between SMARTS and input molecule
    5. Iterate over matches
    6. Check if match is valid (consecutively matches atom indices in molecule).
    7. If match is valid, then generate a dictionary of input molecule index -> index of match in SMARTS 
    8. Generate rdkit mol object from amoeba09 SMARTS
    9. Output .mol file and have babel read as input to generate babel mol object
    10. Compute dictionary of atom index (in amoeba09 SMARTS) -> type number, will later use to help matching indices of same atom type 
    11. Generate matches for SMARTS to amoeba09 SMARTS
    12. Try and pick match from 11. that includes neighbors of atom indices of interest in SMARTS
    13. Generate dictionary of SMARTS atom index -> amoeba09 SMARTS atom index
    14. Generate a dictionary of amoeba09 SMARTS and atom order (location of atoms of interest in amoeba09 SMARTS) -> element + tinker description using input dictionaries and current amoeba09 SMARTS
    15. Using the type dictionary generated earlier, iterate over amoeba09 molecule and for atoms of same type, fill in dictionary from 14 with information of atoms having same type
    16. Using dictionary from 15 and input dictionaries, fill in desired dictionaries mapping atom indices -> tinker types/classes/amoeba09 SMARTS/SMARTS  
    """
    atomindicestotinkertypes={}
    atomindicestotinkerclasses={}
    atomindicestoparametersmartsatomorders={}
    atomindicestoelementtinkerdescrips={}
    atomindicestosmartsatomorders={}
    # STEP 1
    for atomindices,parametersmarts in atomindicesforprmtoparametersmarts.items():
        # STEP 2
        smartsls=atomindicesforprmtosmarts[atomindices]
        smarts=smartsls[0]
        smartsfortransfer=smartsls[1]
        # STEP 3
        substructure = Chem.MolFromSmarts(smarts)
        # STEP 4
        matches=list(rdkitmol.GetSubstructMatches(substructure,maxMatches=10000))
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
        # STEP 5
        for match in matches:
            # STEP 6
            validmatch=CheckMatch(poltype,match,atomindices,smarts,substructure)
            # STEP 7
            if validmatch==True:
               indices=list(range(len(match)))
               smartsindextomoleculeindex=dict(zip(indices,match)) 
               moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}
        # STEP 8
        structure = Chem.MolFromSmarts(parametersmarts)
        fragmentfilepath='fragment.mol'
        if os.path.isfile(fragmentfilepath):
            os.remove(fragmentfilepath)
        # STEP 9
        rdmolfiles.MolToMolFile(structure,fragmentfilepath)
        obConversion = openbabel.OBConversion()
        fragbabelmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(fragmentfilepath)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(fragbabelmol, fragmentfilepath)
        # STEP 10
        fragidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,fragbabelmol,rdkitmol=structure)
        smartindices=[moleculeindextosmartsindex[i] for i in atomindices]
        # STEP 11
        substructure = Chem.MolFromSmarts(smartsfortransfer)
        matches=structure.GetSubstructMatches(substructure)
        # STEP 12
        firstmatch=TryAndPickMatchWithNeighbors(poltype,matches,smartindices,substructure)
        indices=list(range(len(firstmatch)))
        # STEP 13
        smartsindextoparametersmartsindex=dict(zip(indices,firstmatch)) 
        # STEP 14
        parametersmartsordertoelementtinkerdescrip={}
        for parametersmartsatomorderlist,elementtinkerdescrip in smartsatomordertoelementtinkerdescrip.items():
            prmsmarts=parametersmartsatomorderlist[0]
            atomorderlist=parametersmartsatomorderlist[1]
            if prmsmarts==parametersmarts:
                atomorderlist=parametersmartsatomorderlist[1]
                for atomorder in atomorderlist:
                    parametersmartsordertoelementtinkerdescrip[atomorder]=elementtinkerdescrip 
        # STEP 15
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
        # STEP 16
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

def TryAndPickMatchWithNeighbors(poltype,matches,smartindices,substructure):
    """
    Intent: Trying to pick a match from SMARTS -> amoeba09 SMARTS that includes more neighbors (corresponds to including more neighbors of input molecule around atom indices of interest)
    Input: List of list of atom indices that match to amoeba09 SMARTS from SMARTS, the indices in SMARTS that correspond to same indices of interst in input molecule 
    Output: Match that attempt to include neighbors of SMARTS mol matches to indices of interest in input molecule 
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS,GenerateAtomIndexToAtomTypeAndClassForAtomList
    Description:
    1. Generate array of atom indices from SMARTS mol that include neighbors of input smartindices
    2. Iterate over each match and check if all indices from step 1. are included in match. If so, then use this match. 
    """
    # STEP 1
    extendedsmartindices=[]
    for idx in smartindices:
        atom=substructure.GetAtomWithIdx(idx)
        extendedsmartindices.append(idx)
        for natm in atom.GetNeighbors():
            natmidx=natm.GetIdx()
            if natmidx not in extendedsmartindices:
                extendedsmartindices.append(natmidx)
    # STEP 2
    for match in matches:
        indices=list(range(len(match)))
        smartsindextoparametersmartsindex=dict(zip(indices,match)) 
        foundmap=True
        for idx in extendedsmartindices:
            if idx not in smartsindextoparametersmartsindex:
                foundmap=False
        if foundmap==True:
            return match
    firstmatch=matches[0]
    return firstmatch
   
        

def GrabSMARTSList(poltype,smartsatomordertoelementtinkerdescrip):
    """
    Intent: Just extract SMARTS list from dictionary containing SMARTS + atom order -> element + tinker description
    Input: dictionary containing SMARTS + atom order -> element + tinker description
    Output: list of SMARTS
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over input dictionary
    2. Extract SMARTS from key
    3. Save to list 
    """
    smartslist=[]
    # STEP 1
    for smartsatomorder in smartsatomordertoelementtinkerdescrip.keys():
        # STEP 2
        smarts=smartsatomorder[0]
        # STEP 3
        if smarts not in smartslist:
            smartslist.append(smarts)
    return smartslist 

def CheckForPlanerAngles(poltype,listofanglesforprm,mol):
    """
    Intent: Tinker requires angle parameters that should be planer to be anglep rather than just angle, so need to detect this.
    Input: List of all angles in input molecule, mol object
    Output: List of angles in input molecule that are planar
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each angle in molecule
    2. Check if the hybridization of the middle atom is SP2, if so and that atom has 3 neighbors, then this is a planar angle. 
    """
    listofanglesthatneedplanarkeyword=[]
    # STEP 1
    for ls in listofanglesforprm:
        a = mol.GetAtom(ls[0]+1)
        b = mol.GetAtom(ls[1]+1)
        c = mol.GetAtom(ls[2]+1)
        # STEP 2
        anglep=False
        if b.GetHyb()==2: # only for SP2 hyb middle atoms use angp
            neighbs=list(openbabel.OBAtomAtomIter(b))
            if len(neighbs)==3:
                anglep=True
        if anglep==True:
            listofanglesthatneedplanarkeyword.append(ls)
    return listofanglesthatneedplanarkeyword

def ModifyAngleKeywords(poltype,angleprms,listofanglesthatneedplanarkeywordtinkerclassestopoltypeclasses):
    """
    Intent: Modify angle parameters taken from prm file to have anglep keyword if the angle is planer in input molecule
    Input: Array of angle parameter lines, dictionary of planar angles tinker classes -> poltype classes 
    Output: Modified array of angle parameter lines
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each angle parameter line
    2. Grab array of tinker classes from parameter line
    3. Iterate over input dictionary
    4. If there is a match between tinker classes in parameter line and current iteration over dictionary, then replace angle with anglep 
    5. Append modified parameter line to array
    """
    newangleprms=[]
    # STEP 1
    for line in angleprms:
        found=False
        linesplit=line.split()
        # STEP 2
        temp=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
        # STEP 3
        for ls,polclassesls in listofanglesthatneedplanarkeywordtinkerclassestopoltypeclasses.items():
            # STEP 4
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
        # STEP 5
        newangleprms.append(newline)
    return newangleprms

def FilterList(poltype,allitems,listbabel):
    """
    Intent: When wanting to replace bond length/angle lengths with QM optimized bond length/angle lengths, need a way to grab all bonds of same type but also ensure they are physically connected ( you can enumerate many combinations of indices for each bond type but not all are physically connected). 
    Input: All possible bonds/angles, list of physical bonds/angles
    Output: List of bonds/angles with same types that are also in input molecule
    Referenced By: AddOptimizedBondLengths , AddOptimizedAngleLengths
    Description: 
    1. Iterate over each bond/angle in input list
    2. Check if the bond/angle exists in molecule, if so append to array 
    """
    newallitems=[]
    # STEP 1
    for ls in allitems:
        revls=ls[::-1]
        # STEP 2
        if list(ls) in listbabel or list(revls) in listbabel: # need to check reverse too
            newallitems.append(ls)
    return newallitems

def AddOptimizedBondLengths(poltype,optmol,bondprms,bondlistbabel):
    """
    Intent: Modify bond parameters to use QM geometry optimized bond length values 
    Input: QM optimized mol object, list of bond parameter lines, list of bonds in in put molecule
    Output: Modified list of bond parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over list of bond parameter lines
    2. Grab poltype classes from parameter line
    3. For each poltype class, generate list of atom indices that have same type
    4. Now generate all combinations of indices that share same type
    5. Filter indices that arent consective in the molecule
    6. Grab the bond lengths for each bond in list of bonds with same types
    7. Average the bond lengths
    8. Replace bond length in parameter line with the average bond length 
    """
    newbondprms=[]
    # STEP 1
    for line in bondprms:
        linesplit=line.split()
        # STEP 2
        bondtypes=[int(linesplit[1]),int(linesplit[2])]
        bondindices=[]
        # STEP 3
        for prmtype in bondtypes:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,prmtype)
            bondindices.append(keylist)
        # STEP 4
        allbonds = list(itertools.product(bondindices[0], bondindices[1]))
        allbonds=[x for x in allbonds if len(x) == len(set(x))]
        # STEP 5
        allbonds=FilterList(poltype,allbonds,bondlistbabel)
        tot=0
        # STEP 6
        for bond in allbonds:
            blen = optmol.GetBond(int(bond[0]),int(bond[1])).GetLength()
            tot+=blen
        if len(allbonds)==0:
            pass 
        else:
            # STEP 7
            avgbondlength=round(tot/len(allbonds),2)
            linesplit=re.split(r'(\s+)', line)   
            # STEP 8
            linesplit[8]=str(avgbondlength)
            line=''.join(linesplit)
        newbondprms.append(line)
    return newbondprms
      

def AddOptimizedAngleLengths(poltype,optmol,angleprms,anglelistbabel):
    """
    Intent: Modify angle parameters to use QM geometry optimized angle length values 
    Input: QM optimized mol object, list of angle parameter lines, list of angles in in put molecule
    Output: Modified list of angle parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over list of angle parameter lines
    2. Grab poltype classes from parameter line
    3. For each poltype class, generate list of atom indices that have same type
    4. Now generate all combinations of indices that share same type
    5. Filter indices that arent consective in the molecule
    6. Grab the angle lengths for each angle in list of angles with same types
    7. Average the angle lengths
    8. Replace angle length in parameter line with the average angle length 
    """
    newangleprms=[]
    # STEP 1
    for line in angleprms:
        linesplit=line.split()
        # STEP 2
        angletypes=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
        angleindices=[]
        # STEP 3
        for prmtype in angletypes:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,prmtype)
            angleindices.append(keylist)
        # STEP 4
        allangles = list(itertools.product(angleindices[0], angleindices[1],angleindices[2]))
        allangles=[x for x in allangles if len(x) == len(set(x))]
        # STEP 5
        allangles=FilterList(poltype,allangles,anglelistbabel)
        tot=0
        # STEP 6
        for angle in allangles:
            angle=[int(i) for i in angle]
            a = optmol.GetAtom(angle[0])
            b = optmol.GetAtom(angle[1])
            c = optmol.GetAtom(angle[2])
            anglelen = optmol.GetAngle(a,b,c)
            tot+=anglelen
        if len(allangles)==0:
            pass
        else:
            # STEP 7
            avganglelength=round(tot/len(allangles),2)
            linesplit=re.split(r'(\s+)', line)
            linesplit=linesplit[:11]
            linesplit.append('\n')
            transferredangle=float(linesplit[10])
            linesplit[10]=str(avganglelength)
            line=''.join(linesplit)
        line+='\n'
        # STEP 8
        newangleprms.append(line)
    return newangleprms
 


         

def CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol):
    """
    Intent: For determining the correct restraints (usually only want heavy atoms unless special case like hydrogens on Nitrogen atoms), also for determining in some case for torsion missing from database and on a ring that has heavy atoms (not all hydrogen torsion) then can transfer something simple depending on hybridization of atoms on ring.  
    Input: Atom indices for torsion, mol object
    Output: Boolean specifying if all torsion around rotatable bond are hydrogen torsions.
    Referenced By: FindMissingTorsion, torsiongenerator.py - get_torlist , RotatableBondRestraints ,  FrozenBondRestraints 
    Description:
    1. Assume all torsions are hydrogen torsions
    2. Iterate over neighbors of b in a-b-c-d
    3. Iterate over neighbors of c in a-b-c-d
    4. Find new a' and d' such that a is not a' and d is not d' and append torsion to list
    5. For each torsion in list grab atomic numbers of a and d
    6. If there exists a torsion with atomic numbers on a/d not a hydrogen then not all torsions are hydrogen torsions 
    """
    # STEP 1
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
        # STEP 2
        iteratomatom = openbabel.OBAtomAtomIter(b)
        for iaa in iteratomatom:
            # STEP 3
            iteratomatom2 = openbabel.OBAtomAtomIter(c)
            for iaa2 in iteratomatom2:
                ta = iaa.GetIdx()
                tb = bidx
                tc = cidx
                td = iaa2.GetIdx()
                # STEP 4
                if ((ta != tc and td != tb) and not (ta == aidx and td == didx)):
                    torlist.append([ta,tb,tc,td])
        # STEP 5
        for tor in torlist:
            atoms=[mol.GetAtom(i) for i in tor]
            aatomnum=atoms[0].GetAtomicNum()
            datomnum=atoms[3].GetAtomicNum()
            # STEP 6
            if aatomnum!=1 or datomnum!=1:
                allhydrogentor=False
    return allhydrogentor    

def CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol):
    """
    Intent: For determining the correct restraints (usually only want heavy atoms unless special case like hydrogens on Nitrogen atoms), also for determining in some case for torsion missing from database and on a ring that has heavy atoms (not all hydrogen torsion) then can transfer something simple depending on hybridization of atoms on ring.
    Input: Atom indices for torsion, mol object  
    Output: Boolean specifying if all a (or d) in a-b-c-d have hydrogen torsion only on one side.
    Referenced By: FindMissingTorsion, torsiongenerator.py - get_torlist , RotatableBondRestraints ,  FrozenBondRestraints 
    Description:
    1. Assume there is only hydrogen torsion on one side
    2. Iterate over neighbors of b in a-b-c-d
    3. Iterate over neighbors of c in a-b-c-d
    4. Find new a' and d' such that a is not a' and d is not d' and append torsion to list
    5. For each torsion in list grab atomic numbers of a and d
    6. If both a and d are not hydrogens, then there doesnt exist all hydrogen torsion on one side   
    """
    # STEP 1
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
        # STEP 2 
        iteratomatom = openbabel.OBAtomAtomIter(b)
        for iaa in iteratomatom:
            # STEP 3
            iteratomatom2 = openbabel.OBAtomAtomIter(c)
            for iaa2 in iteratomatom2:
                ta = iaa.GetIdx()
                tb = bidx
                tc = cidx
                td = iaa2.GetIdx()
                # STEP 4
                if ((ta != tc and td != tb) and not (ta == aidx and td == didx)):
                    torlist.append([ta,tb,tc,td])
        # STEP 5
        for tor in torlist:
            atoms=[mol.GetAtom(i) for i in tor]
            aatomnum=atoms[0].GetAtomicNum()
            datomnum=atoms[3].GetAtomicNum()
            # STEP 6
            if aatomnum!=1 and datomnum!=1:
                allhydrogentoroneside=False
    

    return allhydrogentoroneside


def FindAllConsecutiveRotatableBonds(poltype,mol,listofbondsforprm):
    """
    Intent: Important for determining missing tor-tor parameters if user decided they wanted tor-tor parameterization
    Input: Mol object, list of bonds in molecule
    Output: Array of consecutive rotatable bonds [[b1,c1],[b2,c2],...] 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over list of bonds
    2. Grab the valence for b and c, if either valence is 1, then it must not be rotatable
    3. Append rotatable bond to new array
    4. Generate all pair combinations of rotatable bonds from previous array
    5. If the current pair of rotatable bonds shares an atom index, then they are consecutive, 
    6. If no atoms in the rotatable bond pair are in a ring, then append to array of consecutive rotatable bonds. 
    """
    totalbondscollector=[]
    newrotbnds=[]
    # STEP 1
    for rotbnd in listofbondsforprm:
        b=rotbnd[0]+1
        c=rotbnd[1]+1
        batom=poltype.mol.GetAtom(b)
        catom=poltype.mol.GetAtom(c)
        # STEP 2
        bval=len([neighb for neighb in openbabel.OBAtomAtomIter(batom)])
        cval=len([neighb for neighb in openbabel.OBAtomAtomIter(catom)])
        if bval<2 or cval<2:
            continue
        # STEP 3
        newrotbnds.append(rotbnd)
    # STEP 4
    combs=list(combinations(newrotbnds,2)) 
    for comb in combs:
        firstbnd=comb[0]
        secondbnd=comb[1]
        total=firstbnd[:]+secondbnd[:]
        totalset=set(total)
        # STEP 5
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
            # STEP 6
            if isinring==True:
                continue
           
            totalbondscollector.append([firstbnd,secondbnd])
    return totalbondscollector

def FindAdjacentMissingTorsionsForTorTor(poltype,torsionsmissing,totalbondscollector,tortorsmissing):
    """
    Intent: If one torsion is missing but an adjacent torsion is not, then consider that tor-tor to be missing. Since that tor-tor is missing, then need to also scan the torsion that was determined to not be missing at first.   
    Input: Array of missing torsions, array of all consecutive rotatable bonds, array of missing tor-tors.
    Output: Updated missing torsions array
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over list of consecutive rotatable bonds
    2. For each bond, check if they are in torsionsmissing array
    3. Also check if the tor-tor is missing, if not then can skip this tor-tor
    4. If either the first bond is missing and second is not or second is missing and first is not then add the missing torsion to missing torsion array
    """
    # STEP 1
    for bndlist in totalbondscollector:
        first=bndlist[0]
        second=bndlist[1]
        # STEP 2
        foundfirst=CheckIfRotatableBondInMissingTorsions(poltype,first,torsionsmissing) 
        foundsecond=CheckIfRotatableBondInMissingTorsions(poltype,second,torsionsmissing)
        # STEP 3
        foundtortormissing=CheckIfRotatableBondInMissingTorTors(poltype,first,second,tortorsmissing)
        if foundtortormissing==False:
            continue
        # STEP 4
        if (foundfirst==False and foundsecond==True and poltype.tortor==True):
            b=second[0]+1
            c=second[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tuple(tor))
        elif (foundfirst==True and foundsecond==False and poltype.tortor==True):
            b=first[0]+1
            c=first[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tuple(tor))
        

    return torsionsmissing

def CheckIfRotatableBondInMissingTorsions(poltype,rotbnd,torsionsmissing):
    """
    Intent: Check if input rotatable bond exists in any torsions in missing torsion array
    Input: Rotatable bond array, array of missing torsions
    Output: Boolean specifying if found rotatable bond in array of missing torsions
    Referenced By: FindAdjacentMissingTorsionsForTorTor
    Description: 
    1. Assume it doesnt exist in missing torsion array
    2. Iterate over array of missing torsions
    3. Grab the bond from current torsion
    4. If the input rotatable bond is same (or reverse) of current bond then it was found 
    """
    # STEP 1
    found=False
    # STEP 2
    for tor in torsionsmissing:
        # STEP 3
        bnd=[tor[1],tor[2]]
        # STEP 4
        if bnd==rotbnd or bnd[::-1]==rotbnd:
            found=True
            break
    return found

def CheckIfRotatableBondInMissingTorTors(poltype,firstbnd,secondbnd,tortormissing):
    """
    Intent: Check if rotatable bonds are in missing tor-tor array. 
    Input: First rotatable bond in tor-tor, second rotatable bond in tor-tor
    Output: Boolean specifying if the rotatable bonds were found in missing tor-tor array or not.
    Referenced By: FindAdjacentMissingTorsionsForTorTor
    Description: 
    1. Assume the rotatable bond pair does not exist in missing tor-tor array
    2. Iterate over missing tor-tors
    3. Extract current first and second rotatable bonds
    4. If the current first rotatable bond matches input first rotatable bond and current second bond matches input second bond (or current first matches second input and current second matches first input), then determine the tor-tor to be found in missing tor-tor array 
    """
    # STEP 1
    foundtortormissing=False
    # STEP 2
    for tortor in tortormissing:
        # STEP 3
        first=[tortor[1],tortor[2]]
        second=[tortor[2],tortor[3]]
        # STEP 4
        if (first==firstbnd or first[::-1]==firstbnd) and (second==secondbnd or second[::-1]==secondbnd):
            foundtortormissing=True
        elif (first==secondbnd or first[::-1]==secondbnd) and (second==firstbnd or second[::-1]==firstbnd):
            foundtortormissing=True

    return foundtortormissing

def FindMissingTorTors(poltype,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,rdkitmol,mol,indextoneighbidxs,totalbondscollector):
    """
    Intent: Determine which tor-tors are missing from database (which determines which ones are needed to be parameterized).  
    Input: Dictionary of tor-tor atom indices to SMARTS from database, dictionary of tor-tor SMARTS + atom order -> array of tor-tor parameter lines, rkdit mol object, openbabel mol object, dictionary of atom index -> neighboring indices, array of consecutive rotatable bonds 
    Output: Array of missing tor-tors
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of tor-tor atom indices -> SMARTS
    2. Iterate over dictionary of tor-tor SMARTS + atom order -> tor-tor parameters
    3. If the SMARTS from both dictoinaries match then
    4. Determine the neighboring atom indices of all atom indices from the current tor-tor
    5. Match the SMARTS to input molecule
    6. Check if the matched atom indices from SMARTS contain all of the neighbors found
    7. If they contain all neighbors then tor-tor is determined not to be missing, if they do then its determined to be missing
    8. For tor-tors not in database, need to iterate over array of consecutive rotatable bonds to determine if need to add tor-tor to array ofmissing tor-tors
    9. Grab the first and second rotatatable bonds for current tor-tor
    10. If tor-tor was not already found in database earlier and also not already appended to array of missing tor-tors, then append the tor-tor to array of missing tor-tors
    11. If tor-tor is in onlyrottortorlist array then add to array of missing tor-tors 
    """
    tortorsmissing=[]
    tortorsfound=[]
    # STEP 1
    for tortorindices,extsmarts in tortorindicestoextsmarts.items(): #only for ones that have match
        # STEP 2
        for tortorsmartsatomorder,parameters in tortorsmartsatomordertoparameters.items():
            smarts=tortorsmartsatomorder[0]
            # STEP 3
            if smarts==extsmarts:
                # STEP 4
                aidx,bidx,cidx,didx,eidx=tortorindices[:]
                firstneighborindexes=indextoneighbidxs[aidx]
                secondneighborindexes=indextoneighbidxs[bidx]
                thirdneighborindexes=indextoneighbidxs[cidx]
                fourthneighborindexes=indextoneighbidxs[didx]
                fifthneighborindexes=indextoneighbidxs[eidx]
                neighborindexes=firstneighborindexes+secondneighborindexes+thirdneighborindexes+fourthneighborindexes+fifthneighborindexes
                substructure = Chem.MolFromSmarts(smarts)
                # STEP 5
                matches=rdkitmol.GetSubstructMatches(substructure)
                matcharray=[]
                for match in matches:
                    for idx in match:
                        if idx not in matcharray:
                            matcharray.append(idx)
                # STEP 6
                check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
                # STEP 7
                if check==False:
                    tortorsmissing.append(tortorindices)
                else:
                    tortorsfound.append(tortorindices)
    if poltype.onlyrotbndslist==None:
        poltype.onlyrotbndslist=[]
    # STEP 8
    for bndlist in totalbondscollector:
        # STEP 9
        first=bndlist[0]
        second=bndlist[1]
        babelfirst=[i+1 for i in first]
        babelsecond=[i+1 for i in second]
        b,c=first[:]
        d=second[0]
        aatom,dnewatom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b+1),poltype.mol.GetAtom(c+1))
        bnewatom,eatom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(c+1),poltype.mol.GetAtom(d+1))
        a=aatom.GetIdx()
        dnew=dnewatom.GetIdx()
        bnew=bnewatom.GetIdx()
        e=eatom.GetIdx()
        indices=[a-1,b,c,d,e-1]
        # STEP 10
        foundtortormissing=CheckIfRotatableBondInMissingTorTors(poltype,[bnew-1,c],[c,dnew-1],tortorsmissing)
        foundtortor=CheckIfRotatableBondInMissingTorTors(poltype,[bnew-1,c],[c,dnew-1],tortorsfound)
        if foundtortormissing==False and foundtortor==False:
            if len(poltype.onlyrottortorlist)==0:
                tortorsmissing.append(indices)
            else:
                ls=[b+1,c+1,d+1]
                if ls in poltype.onlyrottortorlist or ls[::-1] in poltype.onlyrottortorlist:
                    tortorsmissing.append(indices)

        else:
            # STEP 11
            ls=[b+1,c+1,d+1]
            if ls in poltype.onlyrottortorlist or ls[::-1] in poltype.onlyrottortorlist:
                tortorsmissing.append(indices)

    return tortorsmissing


def RingAtomicIndices(poltype,mol):
    """
    Intent: Grab all rings in input molecule
    Input: mol object
    Output: List of list of ring atomic indices
    Referenced By: FindMissingTorsions , DefaultAromaticMissingTorsions
    Description: 
    1. Call SSSR method
    2. Iterate over ring objects
    3. Grab atomic indices for current ring
    4. Append ring atomic indices to array 
    """
    # STEP 1
    sssr = mol.GetSSSR()
    atomindices=[]
    # STEP 2
    for ring in sssr:
        # STEP 3
        ringatomindices=GrabRingAtomIndices(poltype,mol,ring)
        # STEP 4
        atomindices.append(ringatomindices)        
    return atomindices

def GrabRingAtomIndices(poltype,mol,ring):
    """
    Intent: For a given ring object, extract atomic indices
    Input: Mol object, ring object
    Output: Array of ring atomic indices
    Referenced By: RingAtomicIndices
    Description: 
    1. Iterate over mol object
    2. If atom is in ring object, then append to array
    """
    ringatomindices=[]
    atomiter=openbabel.OBMolAtomIter(mol)
    # STEP 1
    for atom in atomiter:
        atomidx=atom.GetIdx()
        # STEP 2
        if ring.IsInRing(atomidx)==True:
            ringatomindices.append(atomidx)
    return ringatomindices

def GrabRingAtomIndicesFromInputIndex(poltype,atomindex,atomindices):
    """
    Intent: Only return ring that contains the input atom index
    Input: Atom index of interest, array of rings
    Output: Ring of interest
    Referenced By: FindMissingTorsions
    Description: 
    1. Iterate over rings in array of rings
    2. If atom index of interest is in ring, then return ring 
    """
    # STEP 1
    for ring in atomindices:
        # STEP 2
        if atomindex in ring:
            return ring

    ring=None
    return ring

def GrabIndicesInRing(poltype,babelindices,ring):
    """
    Intent: Grab only the indices from input array that are in input ring 
    Input: Array of atom indices, ring indices
    Output: Array of indices from babelindices that are in ring
    Referenced By: FindMissingTorsions
    Description: 
    1. Iterate over indices in ring
    2. If index is also in input array of atom indices, then save in new array
    """
    ringtorindices=[]
    # STEP 1
    for index in ring:
        # STEP 2
        if index in babelindices:
            ringtorindices.append(index)

    return ringtorindices


def FindMissingTorsions(poltype,torsionindicestoparametersmartsenv,rdkitmol,mol,indextoneighbidxs):
    """
    Intent: Determine which torsions are missing and hence need to be parameterized. 
    Input: Dictionary of torsion atomic indices -> amoeba09 SMARTS + atom order, rdkit mol object, openbabel mol object, dictionary of atom index -> neighboring indices 
    Output: Array of missing torsion indices
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over torsion atomic indices and amoeba09 SMARTS
    2. If user specifies to not parameterize torsion, then dont add any torsions to array of missing torsions.
    3. If detect linear angle (for example nitrile group) then dont want to parameterize torsion. 
    4. Check if the amoeba09 SMARTS match contains all neighbors of the torsion, if not the match is "bad"
    5. If user specifies to parameterize torsions about rotatable bond in onlyrotbndslist, then add torsion to array of missing torsions
    6. If there exists a * or ~ then the SMARTS transfer is determined to be "bad"
    7. If user specifies to not parameterize torsion in dontrotbndslist, then transfer the amoeba09 SMARTS match
    8. Check for special case when whole molecule is a giant ring (then set ringbond=False)
    9. If SMARTS is "bad" transfer
    10.If bond is ring bond and the torsion contains no SP2 atoms, then want to transfer non-aromatic torsion parameters. Transfer hydrogen torsion from alkane. If user specifies to scan ring bond, can add to list and poltype will fragment molecule, otherwise will transfer heavy torsion from alkane also. 
    11. If ring bond and also middle two atoms are SP2, then transfer torsion from benzene
    12. If one of the middle atoms is SP2, or both are not SP2 then transfer torsion from alkane. 
    13. If middle bond is not in a ring, and torsion is a hydrogen torsion, then if middle two atoms are both SP2 transfer from benzene, else transfer from alkane
    14. For heavy torsion not in ring, append to array of missing torsions (if not double/triple bond).
    15. Otherwise if the SMARTS transfer is "good", check if user specified to parameterize all torsions or this specific torsion, then append torsion to array of missing torsions in this case 
    """
    torsionsmissing=[]
    poormatchingaromatictorsions=[]
    poormatchingpartialaromatictorsions=[]
    torsionstozerooutduetocolinear=[]
    # STEP 1
    for torsionindices,smartsenv in torsionindicestoparametersmartsenv.items():
        # STEP 2
        if poltype.dontdotor==True:
            continue
        aidx,bidx,cidx,didx=torsionindices[:]
        babelindices=[i+1 for i in torsionindices]

        allhydrogentor=CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol)
        allhydrogentoroneside=CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol)
        atoma=rdkitmol.GetAtomWithIdx(aidx)
        atomd=rdkitmol.GetAtomWithIdx(didx)
        atomicnumatoma=atoma.GetAtomicNum()
        atomicnumatomd=atomd.GetAtomicNum()
        abidx,bbidx,cbidx,dbidx=babelindices[:]
        
        bond=mol.GetBond(bbidx,cbidx)
        bondorder=bond.GetBondOrder()
        babelatoms=[mol.GetAtom(i) for i in babelindices]
        aatom,batom,catom,datom=babelatoms[:]

        middlebond=mol.GetBond(babelindices[1],babelindices[2])
        ringbond=middlebond.IsInRing()
        # STEP 3
        firstangle=mol.GetAngle(aatom,batom,catom)
        secondangle=mol.GetAngle(batom,catom,datom)
        if firstangle<0:
            firstangle=firstangle+360
        if secondangle<0:
            secondangle=secondangle+360
        angletol=3.5
        if np.abs(180-firstangle)<=angletol or np.abs(180-secondangle)<=angletol:
            torsionstozerooutduetocolinear.append(torsionindices)
            continue
         
        
        atomvals=[len([neighb for neighb in openbabel.OBAtomAtomIter(a)]) for a in babelatoms]
        atomnums=[a.GetAtomicNum() for a in babelatoms]
        batomnum=atomnums[1]
        catomnum=atomnums[2]
           
        ringbools=[a.IsInRing() for a in babelatoms]
        arobools=[a.IsAromatic() for a in babelatoms]
        bnd=[babelindices[1],babelindices[2]]
         
        ringb=ringbools[1]
        ringc=ringbools[2]
        aroa=arobools[0]
        arob=arobools[1]
        aroc=arobools[2]
        arod=arobools[3]
        hybs=[a.GetHyb() for a in babelatoms]
        hybb=hybs[1]
        hybc=hybs[2]
        # STEP 4
        firstneighborindexes=indextoneighbidxs[aidx]
        secondneighborindexes=indextoneighbidxs[bidx]
        thirdneighborindexes=indextoneighbidxs[cidx]
        fourthneighborindexes=indextoneighbidxs[didx]
        neighborindexes=firstneighborindexes+secondneighborindexes+thirdneighborindexes+fourthneighborindexes
        smarts=smartsenv[0]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)
        for match in matches:
            allin=True
            for idx in torsionindices:
                if idx not in match:
                    allin=False
            if allin==True:
                matcharray=match
                break

        check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
        # STEP 5
        if len(poltype.onlyrotbndslist)!=0:
            if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                check=False 
        # STEP 6
        if '~' in smarts or '*' in smarts:
            check=False
        # STEP 7
        if [bbidx,cbidx] in poltype.dontrotbndslist or [cbidx,bbidx] in poltype.dontrotbndslist:
            check=True
        # STEP 8
        if ringbond==True:
            atomindices=RingAtomicIndices(poltype,mol)
            therings=torgen.GrabAllRingsContainingMostIndices(poltype,atomindices,babelindices,2)
            if (len(therings)>0) and poltype.dontfrag==False:
                if len(therings[0])>7: # special case where whole molecule is a ring then dont consider ring bond
                    if hybs[1]!=2 and hybs[2]!=2:
                        ringbond=False

            ring=GrabRingAtomIndicesFromInputIndex(poltype,babelindices[1],atomindices)
            ringtorindices=GrabIndicesInRing(poltype,babelindices,ring)
        if (bnd in poltype.partialdoublebonds or bnd[::-1] in poltype.partialdoublebonds) and poltype.rotalltors==False and ([bbidx,cbidx] not in poltype.onlyrotbndslist and [cbidx,bbidx] not in poltype.onlyrotbndslist) and check==True:
            continue
        # STEP 9
        if check==False:
            if ringbond==True:
                # STEP 10    
                if (2 not in hybs): # non-aromatic torsion want parameters for 
                    if poltype.transferanyhydrogentor==True and (atomicnumatoma==1 or atomicnumatomd==1) and (allhydrogentor==False and allhydrogentoroneside==False): # then here transfer torsion because can pick up most QM-MM on heavy atoms, less parameters to fit
                        poormatchingpartialaromatictorsions.append(torsionindices)
                    else: # if dont have heavy atoms on either side then just fit the hydrogen torsion
                        if poltype.nonaroringtor1Dscan==True or poltype.refinenonaroringtors==True: 
                            found=False
                            if len(poltype.onlyrotbndslist)!=0:
                                if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                                    found=True
                            else:
                                found=True
                            if found==True:
                                if len(ring)>3:
                                    if torsionindices not in torsionsmissing and poltype.dontfrag==False: # make sure fragmenter is on (wont work for < 25 atoms by default)
                                        torsionsmissing.append(torsionindices)
                            else:
                                poormatchingpartialaromatictorsions.append(torsionindices)
                        else:
                            poormatchingpartialaromatictorsions.append(torsionindices)

                # STEP 11        
                elif hybs[1]==2 and hybs[2]==2:
                    if torsionindices not in poormatchingaromatictorsions:
                        poormatchingaromatictorsions.append(torsionindices)
                # STEP 12
                elif (hybs[1]==2 and hybs[2]!=2) or (hybs[1]!=2 and hybs[2]==2):
                    if torsionindices not in poormatchingpartialaromatictorsions:
                        poormatchingpartialaromatictorsions.append(torsionindices)
                elif (hybs[1]!=2 and hybs[2]!=2) and (hybs[0]==2 or hybs[3]==2):
                    if torsionindices not in poormatchingpartialaromatictorsions:
                        poormatchingpartialaromatictorsions.append(torsionindices)

            else:
                # STEP 13
                if poltype.transferanyhydrogentor==True and poltype.rotalltors==False and (atomicnumatoma==1 or atomicnumatomd==1) and (allhydrogentor==False and allhydrogentoroneside==False): # then here transfer torsion because can pick up most QM-MM on heavy atoms, less parameters to fit
                    if hybs[1]==2 and hybs[2]==2:
                        if torsionindices not in poormatchingaromatictorsions:
                            poormatchingaromatictorsions.append(torsionindices)
                    elif (hybs[1]==2 and hybs[2]!=2) or (hybs[1]!=2 and hybs[2]==2):
                        if torsionindices not in poormatchingpartialaromatictorsions:
                            poormatchingpartialaromatictorsions.append(torsionindices)
                    elif (hybs[1]!=2 and hybs[2]!=2) and (hybs[0]==2 or hybs[3]==2):
                        if torsionindices not in poormatchingpartialaromatictorsions:
                            poormatchingpartialaromatictorsions.append(torsionindices)
                    else:
                        if torsionindices not in poormatchingpartialaromatictorsions:
                            poormatchingpartialaromatictorsions.append(torsionindices)


                else:
                    # STEP 14
                    if len(poltype.onlyrotbndslist)!=0:
                        if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                            if bondorder!=1: 
                                continue

                            if torsionindices not in torsionsmissing:
                                torsionsmissing.append(torsionindices)
                    else:
                        if bondorder!=1: # for double/triple bond transfer something stiff (benzene heavy torsion) if missing parameters
                            poormatchingaromatictorsions.append(torsionindices)
                        if torsionindices not in torsionsmissing:
                            torsionsmissing.append(torsionindices)

                    continue
        # STEP 15
        else:
            if poltype.rotalltors==True and ringbond==False:
                if bondorder<2:
                    if torsionindices not in torsionsmissing:
                        torsionsmissing.append(torsionindices)
            if len(poltype.onlyrotbndslist)!=0:
                if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                    if bondorder<2:
                        if torsionindices not in torsionsmissing:
                            torsionsmissing.append(torsionindices)
    return torsionsmissing,poormatchingaromatictorsions,poormatchingpartialaromatictorsions,torsionstozerooutduetocolinear


def GrabAllRingsContainingIndices(poltype,atomindices,babelindices):
    """
    Intent: Needed to check if whole molecule is one big ring or not.   
    Input: Array of rings, array of atomic indices corresponding to torsion 
    Output: Array of rings containing atomic indices from torsion
    Referenced By: DefaultAromaticMissingTorsions, FindMissingTorsions
    Description:
    1. Iterate over rings in array of rings
    2. Assume all atomic indices exist in ring
    3. Iterate over atomic indices in array of atomic indices, if atomic index doesnt exist in ring, then not all indices are in ring
    4. If all atomic indices are in ring, then append ring to array of rings 
    """
    rings=[]
    # STEP 1
    for ring in atomindices:
        # STEP 2
        allin=True
        # STEP 3
        for atomindex in babelindices:
            if atomindex not in ring:
                allin=False
        # STEP 4
        if allin==True:
            rings.append(ring)
    return rings


def FindAllNeighborIndexes(poltype,rdkitmol):
    """
    Intent: Used for checking if SMARTS match contains all neighbors of atom indices attempting to match parameters for
    Input: Rdkit mol object
    Output: Dictionary of atomic indices -> neighboring atomic indices
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over atoms in mol object
    2. Iterate over neighboring atoms for current atom
    3. Append all neighbors to dictionary for current atom index entry 
    """
    indextoneighbidxs={}
    # STEP 1
    for atm in rdkitmol.GetAtoms():
        atmidx=atm.GetIdx()
        if atmidx not in indextoneighbidxs.keys():
            indextoneighbidxs[atmidx]=[]
        # STEP 2
        for neighbatm in atm.GetNeighbors():
            neighbatmidx=neighbatm.GetIdx()
            # STEP 3
            if neighbatmidx not in indextoneighbidxs[atmidx]:
                indextoneighbidxs[atmidx].append(neighbatmidx)
            

    return indextoneighbidxs

def CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,smartsindexes):
    """
    Intent: Check to see if atomic indices from SMARTS match contain neighboring indices for parameters being assigned
    Input: Array of neighboring indices, array of indices from SMARTS match
    Output: Boolean specifying if all neighbors exist in SMARTS match or notn
    Referenced By: FindMissingTorsions , MatchExternalSMARTSToMolecule , FindMissingParameters
    Description: 
    1. Assume all neighbors exist in SMARTS match
    2. Iterate over neighboring indices
    3. If index not in neighboring indices, then not all neighbors exist in SMARTS match 
    """
    # STEP 1
    check=True
    # STEP 2
    for idx in neighborindexes:
        # STEP 3
        if idx not in smartsindexes:
            check=False
    return check

def ZeroOutMissingStrbnd(poltype,anglemissingtinkerclassestopoltypeclasses,strbndprms):
    """
    Intent: When SMARTS matches for strbnd are not "good" then instead of assigning generic parameters, zero this term out 
    Input: Dictionary of missing angle indices -> corresponding poltype type numbers , array of strbnd parameter lines 
    Output: Modifed array of strbnd parameter lines
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over strbnd parameter lines
    2. Grab the classes from parameter line
    3. Iterate over keys from dictionary of missing angle indices
    4. If there is a match, then modify parameters in line to be 0
    5. Append to array of modified strbnd parameter lines
    """
    newstrbndprms=[]
    # STEP 1
    for line in strbndprms:
        linesplit=line.split()
        # STEP 2
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
        # STEP 3
        for sublist in anglemissingtinkerclassestopoltypeclasses.values():
            # STEP 4
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[8]='0'
                newlinesplit[10]='0'
                line=''.join(newlinesplit)
        # STEP 5
        newstrbndprms.append(line)

    return newstrbndprms

def AssignAngleGuessParameters(poltype,anglemissingtinkerclassestopoltypeclasses,angleprms,indextoneighbidxs):
    """
    Intent: For angle parameters that have a "bad" SMARTS match, just assign generic angle parameters based on element and valance 
    Input: Dictionary of missing atomic angle indices -> poltype classes, array of angle parameter lines, dictionary of atomic index -> neighboring atomic indices
    Output: Modified array of angle parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of angle parameters lines
    2. Extract classes from the line 
    3. Iterate over dictionary of missing angle indices
    4. If the classes exist in current dictionary item
    5. Grab all atomic indices that have the same atom types
    6. Extract atomic numbers and explict valencies for each atomic index
    7. Generate angle guess parameters based on atomic numbers and valencies 
    8. Update the line with parameter guess
    9. Append line to modified angle parameter array 
    """
    newangleprms=[]
    # STEP 1
    for line in angleprms:
        linesplit=line.split()
        # STEP 2
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
        # STEP 3
        for tinkerclasses,sublist in anglemissingtinkerclassestopoltypeclasses.items():
            # STEP 4
            found=False
            if classes in sublist:
                found=True
            elif classes[::-1] in sublist: 
                found=True
            if found==True:
                # STEP 5
                indices=[]
                for poltypeclass in classes:
                    keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
                    indices.append(keylist)
                exampleindices=None
                combs = list(itertools.product(*indices))
                for comb in combs:
                    poscomb=[np.abs(i) for i in comb]
                    poscomb=[k-1 for k in poscomb]
                    checkconsec=CheckIfAtomsConnected(poltype,poscomb,indextoneighbidxs)
                    if checkconsec==True:
                        exampleindices=[k-1 for k in comb]
                        break
                # STEP 6
                exampleindices=[int(i) for i in exampleindices]
                atoms=[poltype.rdkitmol.GetAtomWithIdx(k) for k in exampleindices]
                atomicnums=[a.GetAtomicNum() for a in atoms]
                atomicval=[a.GetExplicitValence() for a in atoms]
                # STEP 7
                angleguess=AngleGuess(poltype,atomicnums[0],atomicnums[1],atomicnums[2],atomicval[0],atomicval[1],atomicval[2])
                # STEP 8
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[8]=str(angleguess)
                line=''.join(newlinesplit)
        # STEP 9
        newangleprms.append(line)

    return newangleprms


def CheckIfAtomsConnected(poltype,poscomb,endindextoneighbs):
    """
    Intent: When generating all combinations of atomic indices that map to same type numbers, only pick tuples of indices that are conseuctively connected (otherwise doesnt exist) 
    Input: Array of atom indices, dictionary of atomic index -> neighboring indices
    Output: Boolean that specifies if all atoms are consecutively connected or not
    Referenced By: AssignBondGuessParameters , AssignAngleGuessParameters
    Description:
    1. Assume that atoms are consecutively connected
    2. Iterate over indices in array of atomic indices
    3. Grab neighboring indices
    4. If the next index in array of atomic indices is not in neighbors, then indices are not consecutive 
    """
    # STEP 1
    checkconsec=True
    if len(poscomb)>1:
        # STEP 2
        for i in range(len(poscomb)-1):
            index=poscomb[i]
            nextindex=poscomb[i+1]
            # STEP 3
            indexneighbs=endindextoneighbs[index]
            # STEP 4
            if nextindex not in indexneighbs:
                checkconsec=False



    return checkconsec

def AssignBondGuessParameters(poltype,bondmissingtinkerclassestopoltypeclasses,bondprms,indextoneighbidxs):
    """
    Intent: For bond parameters that have a "bad" SMARTS match, just assign generic bond parameters based on element and valance 
    Input: Dictionary of missing atomic bond indices -> poltype classes, array of bond parameter lines, dictionary of atomic index -> neighboring atomic indices
    Output: Modified array of bond parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of bond parameters lines
    2. Extract classes from the line 
    3. Iterate over dictionary of missing angle indices
    4. If the classes exist in current dictionary item
    5. Grab all atomic indices that have the same atom types
    6. Extract atomic numbers and explict valencies for each atomic index
    7. Generate bond guess parameters based on atomic numbers and valencies 
    8. Update the line with parameter guess
    9. Append line to modified bond parameter array 
    """
    newbondprms=[]
    # STEP 1
    for line in bondprms:
        linesplit=line.split()
        # STEP 2
        classes=tuple([int(linesplit[1]),int(linesplit[2])])
        # STEP 3
        for tinkerclasses,sublist in bondmissingtinkerclassestopoltypeclasses.items():
            # STEP 4
            found=False
            if classes in sublist:
                found=True
            elif classes[::-1] in sublist: 
                found=True
            if found==True:
                # STEP 5
                indices=[]
                for poltypeclass in classes:
                    keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
                    indices.append(keylist)
                exampleindices=None
                combs = list(itertools.product(*indices))
                for comb in combs:
                    poscomb=[np.abs(i) for i in comb]
                    poscomb=[k-1 for k in poscomb]
                    checkconsec=CheckIfAtomsConnected(poltype,poscomb,indextoneighbidxs)
                    if checkconsec==True:
                        exampleindices=[k-1 for k in comb]
                        break 
                # STEP 6
                exampleindices=[int(i) for i in exampleindices]
                atoms=[poltype.rdkitmol.GetAtomWithIdx(k) for k in exampleindices]
                atomicnums=[a.GetAtomicNum() for a in atoms]
                atomicval=[a.GetExplicitValence() for a in atoms]
                # STEP 7
                bondguess=BondGuess(poltype,atomicnums[0],atomicnums[1],atomicval[0],atomicval[1])
                # STEP 8
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[6]=str(bondguess)
                line=''.join(newlinesplit)
        # STEP 9
        newbondprms.append(line)

    return newbondprms




def ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms):
    """
    Intent: For torsions that need to be fit, need to zero out the torsion parameters prior to fitting. 
    Input: Dictionary of missing torsion tinker classes -> poltype types, array of torsion parameter lines
    Output: Modifed array of torsion parameter lines
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over array of torsion parameter lines
    2. Extract classes from current parameter line
    3. Iterate over dictionary items
    4. If current classes exist in current dictionary item, then zero out the parameters in the line
    5. Append modified line to array of modified torsion parameter lines 
    """
    newtorsionprms=[]
    # STEP 1
    for line in torsionprms:
        linesplit=line.split()
        # STEP 2
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])])
        # STEP 3
        for sublist in torsionsmissingtinkerclassestopoltypeclasses.values():
            # STEP 4
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                splitafter=newlinesplit[10:]
                splitafter[0]='0'
                splitafter[6]='0'
                splitafter[12]='0'
                newlinesplit=newlinesplit[:10]+splitafter
                line=''.join(newlinesplit)
        # STEP 5
        newtorsionprms.append(line)

    return newtorsionprms



def DefaultAromaticMissingTorsions(poltype,arotorsionsmissingtinkerclassestopoltypeclasses,partialarotorsionsmissingtinkerclassestopoltypeclasses,torsionprms,mol): # transfer bezene aromatic torsion paramerers from amoeba09
    """
    Intent: Some torsions parameters need to transfer from alkane or from benezne, this function hardcodes that torsion transfer
    Input: Dictionary of torsions that need benzene transfer tinker classes -> poltype types, dictionary of torsion tinker classes that need alkane transfer -> poltype types , array of torsion parameter lines, mol object
    Output: Modified torsion parameter line array
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Define hardcoded parameters and tinker type descriptions (for adding helpful comments to parameter transfer in key file)
    2. Find all hydrogens types there are in the molecule, will help with counting how many hydrogen classes in torsion
    3. Iterate over lines of torsion parameter array
    4. Extract torsion classes from parameter line
    5. Determine possible corresponding atomic indices from classes
    6. Determine hybridization of atoms in torsion and if a ring bond or not
    7. Count how many hydrogens in current torsion 
    8. Check if curren torsion is in either dictionary for benzene transfers or alkane transfers, if so determine which parameters to transfer and which comments to add based on number of hydrogens
    9. Update current torsion line with new parameters
    10. If torsion is on ring bond and depending on bond order (if double/triple etc), then add extra comment about bond order to key file 
    11. Append modified torsion to array of modified torsion parameter lines
    """
    # STEP 1
    newtorsionprms=[]
    poorarohydcounttoprms={0:[-.67,6.287,0],1:[.55,6.187,-.55],2:[0,6.355,0]}
    poorpartialarohydcounttoprms={0:[.854,-.374,.108],1:[0,0,.108],2:[0,0,.299]}
    poorarohydcounttodescrips={0:["Benzene C","Benzene C","Benzene C","Benzene C"],1:["Benzene HC","Benzene C","Benzene C","Benzene C"],2:["Benzene HC","Benzene C","Benzene C","Benzene HC"]}
    poorpartialarohydcounttodescrips={0:["Alkane -CH2-","Alkane -CH2-","Alkane -CH2-","Alkane -CH2-"],1:["Alkane -H2C-","Alkane -CH2-","Alkane -CH2-","Alkane -CH2-"],2:["Alkane -H2C-","Alkane -CH2-","Alkane -CH2-","Alkane -H2C-"]}
    ringextra='# Ring bond detected for missing torsion'+'\n'
    doubleextra='# Double bond detected for missing torsion'+'\n'
    singleextra='# Single bond detected for missing torsion'+'\n'
    hydextra='# Transferring hydrogen torsion to reduce fit parameters'+'\n'
    arotorsionlinetodescrips={} 
    # STEP 2
    hydclasses=[]
    for atom in poltype.rdkitmol.GetAtoms():
        atomicnum=atom.GetAtomicNum()
        babelindex=atom.GetIdx()+1
        symclass=poltype.idxtosymclass[babelindex]
        if atomicnum==1:
            hydclasses.append(symclass)
    # STEP 3 
    for line in torsionprms:
        # STEP 4
        linesplit=line.split()
        a=int(linesplit[1])
        b=int(linesplit[2])
        c=int(linesplit[3])
        d=int(linesplit[4])
        allcombs=[]
        ls=[a,b,c,d]
        # STEP 5
        for typenum in ls: 
            indices=GrabKeysFromValue(poltype,poltype.idxtosymclass,typenum)
            allcombs.append(indices)
        combs=list(itertools.product(*allcombs))
        for comb in combs:
            aindex=comb[0]
            bindex=comb[1]
            cindex=comb[2]
            dindex=comb[3]
            firstbond=mol.GetBond(aindex,bindex)
            midbond=mol.GetBond(bindex,cindex)
            lastbond=mol.GetBond(cindex,dindex)
            if midbond!=None and firstbond!=None and lastbond!=None:
                break 
        # STEP 6
        babelindices=list(comb)
        babelatoms=[mol.GetAtom(i) for i in babelindices]
        hybs=[a.GetHyb() for a in babelatoms]
        ringbond=midbond.IsInRing()
        if ringbond==True:
            atomindices=RingAtomicIndices(poltype,mol)
            therings=torgen.GrabAllRingsContainingMostIndices(poltype,atomindices,babelindices,2)
            if (len(therings)>0) and poltype.dontfrag==False:
                if len(therings[0])>7: # special case where whole molecule is a ring then dont consider ring bond
                    if hybs[1]!=2 and hybs[2]!=2:
                        ringbond=False

        bondorder=midbond.GetBondOrder()      
        classes=tuple([a,b,c,d])
        # STEP 7
        hydcount=0
        if a in hydclasses:
            hydcount+=1
        if d in hydclasses:
            hydcount+=1
        # STEP 8
        found=False
        for sublist in arotorsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                found=True
                prms=poorarohydcounttoprms[hydcount]
                descrips=poorarohydcounttodescrips[hydcount]
                break
        for sublist in partialarotorsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                found=True
                prms=poorpartialarohydcounttoprms[hydcount]
                descrips=poorpartialarohydcounttodescrips[hydcount]
                break
        # STEP 9
        if found==True:
            newlinesplit=re.split(r'(\s+)', line)
            splitafter=newlinesplit[10:]
            splitafter[0]=str(prms[0])
            splitafter[6]=str(prms[1])
            splitafter[12]=str(prms[2])
            newlinesplit=newlinesplit[:10]+splitafter
            line=''.join(newlinesplit)
            # STEP 10 
            if ringbond==True:
                extra=ringextra
            elif ringbond==False and bondorder!=1:
                extra=doubleextra
            else:
                extra=singleextra
            if hydcount>0 and ringbond==False:
                extra+=hydextra
            arotorsionlinetodescrips[line]=[extra,descrips]
        # STEP 11
        newtorsionprms.append(line)

    return newtorsionprms,arotorsionlinetodescrips


def GrabTorsionParameterCoefficients(poltype,torsionprms):
    """
    Intent: When want to refine non-aromatic ring torsion parameters
    Input: Array of torsion parameter lines
    Output: Dictionary of torsion classes -> parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over array of torsion parameters
    2. Grab torsion classes from current parameter line and also parameters
    3. Append classes and parameters to dictionary 
    """
    torsionkeystringtoparameters={}
    # STEP 1
    for line in torsionprms:
        linesplit=line.split()
        # STEP 2
        key = '%s %s %s %s' % (linesplit[1], linesplit[2], linesplit[3], linesplit[4])
        parameters=[float(linesplit[5]),float(linesplit[8]),float(linesplit[11])]
        allzero=True
        for prm in parameters:
            if prm!=0:
                allzero=False
        # STEP 3
        if allzero==False:
            torsionkeystringtoparameters[key]=parameters
       
    return torsionkeystringtoparameters

def PruneDictionary(poltype,keysubset,dic):
    """
    Intent: Want to make a dictionary that is a subset of input dictionary with only keys given as input, for example only want dictionary of missing bonds/angles from dictionary of all bonds/angles. 
    Input: Array of keys that is subset of keys in input dictionary, input dictionary
    Output: Subset dictionary of input dictionary, array of keys not included in subet
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of keys
    2. If key (or the reverse) is in the dictionary then
    3. Append key and value to new dictionary
    4. Append removed key to array of removed keys
    """
    newdic={}
    removeditems=[]
    # STEP 1
    for key in keysubset:
        # STEP 2 & 3
        if key in dic.keys():
            value=dic[key]
            newdic[key]=value

        elif key[::-1] in dic.keys():
            value=dic[key[::-1]]
            newdic[key]=value
        else:
            # STEP 4
            removeditems.append(key)

    return newdic,removeditems


def TinkerClassesToPoltypeClasses(poltype,indicestotinkerclasses,formatorder=True):
    """
    Intent: Make a dictionary of tinker classes (for vdw/bond/angle/torsion) -> list of all poltype types that match in molecule
    Input: Dictionary of atomic indices -> tinker classes, boolean that gives canonical ordering to tinker classes tuple and poltype types tuple
    Output: Dictionary of tinker classes (for vdw/bond/angle/torsion) -> list of all poltype types that match in molecule
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of atomic indices to tinker classes
    2. Generate poltype types from atomic indices and internal dictionary of atomic index -> poltype type
    3. Give canoncal ordering to tinker classes and poltype types
    4. Append poltype types array to dictionary
    """
    tinkerclassestopoltypeclasses={}
    poltypeclassesalreadyassigned=[]
    # STEP 1
    for indices,tinkerclasses in indicestotinkerclasses.items():
        # STEP 2
        babelindices=[i+1 for i in indices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        revtinkerclasses=tinkerclasses[::-1]
        if poltypeclasses in poltypeclassesalreadyassigned:
            continue
        # STEP 3
        if formatorder==True:
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
        # STEP 4
        if tuple(tinkerclasses) not in tinkerclassestopoltypeclasses.keys():
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)]=[]
        if tuple(poltypeclasses) not in tinkerclassestopoltypeclasses[tuple(tinkerclasses)]: 
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)].append(tuple(poltypeclasses))

    return tinkerclassestopoltypeclasses


def ConvertIndicesDictionaryToPoltypeClasses(poltype,indicestovalue,indicestotinkerclasses,tinkerclassestopoltypeclasses):
    """
    Intent: Convert dictionary of atomic indices -> values to poltype types -> value
    Input: Dictionary of atomic indices -> value, dictionary of atomic indices -> tinker classes, dictionary of tinker classes -> poltype types
    Output: Dictionary of poltype classes -> value
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of atomic indices -> value
    2. If indices are in dictionary of indices -> tinker classes then
    3. Extract list of poltype type arrays from tinker classes
    4. Compute the poltype types for the current atomic indices
    5. If there is match between current poltype types and the ones from indices-> tinker classes dictionary then append to new dictionary 
    """
    poltypeclassestovalue={}
    # STEP 1
    for indices,value in indicestovalue.items():
        # STEP 2
        if indices in indicestotinkerclasses.keys():
            tinkerclasses=tuple(indicestotinkerclasses[indices])
            # STEP 3
            if tinkerclasses in tinkerclassestopoltypeclasses.keys():
                poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses])
            elif tinkerclasses[::-1] in tinkerclassestopoltypeclasses.keys():
                poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses[::-1]])
            else:
                continue
            # STEP 4
            babelindices=[i+1 for i in indices]
            symclasses=tuple([poltype.idxtosymclass[i] for i in babelindices])
            # STEP 5
            for ls in poltypeclasses:
                if ls==symclasses or ls==symclasses[::-1]:
                    ls=tuple([ls])
                    if ls not in poltypeclassestovalue.keys():
                        poltypeclassestovalue[ls]=[]
                        if value not in poltypeclassestovalue[ls]:
                            poltypeclassestovalue[ls].append(value)
    return poltypeclassestovalue 

def GrabTypesFromPrmLine(poltype,ls):
    """
    Intent: Used for when searching for poltype types dictionary items that contain types from parameter line
    Input: Array of strings from parameter line
    Output: The type numbers from parameter line
    Referenced By: SearchForPoltypeClasses
    Description:
    1. Iterate over strings in array
    2. Keep only strings that are integers (types/classes are integers)
    3. Keep only types up to certain length in array (integers can also be in parameter line like torsion has integers for phase etc) 
    """
    typenums=[]
    # STEP 1
    for e in ls:
        isint=RepresentsInt(e)  
        # STEP 2
        if isint==True:
            typenums.append(e)
    # STEP 3
    if len(typenums)>=4:
        if typenums[-2]=='0' and typenums[-1]=='0':
            typenums=typenums[:2]
    if len(typenums)>4:
        typenums=typenums[:4]     
    return typenums


def SearchForPoltypeClasses(poltype,prmline,poltypeclasseslist):
    """
    Intent: Given a parameter line, extract poltype types from list of poltype types. 
    Input: Parameter line, list of list of poltype types (for vdw or bond or angle or torsion)
    Output:
    Referenced By: FindPotentialMissingParameterTypes , MapParameterLineToTransferInfo
    Description:
    1. Extract type numbers from input parameter line
    2. Iterate over array of poltype types array
    3. Iterate over each poltype types array in list
    4. Check for a match between type numbers from parameter line and current types in iteration
    5. If there is a match return current poltype types array and list of poltype types array
    6. If no matches return None
    """
    listofpoltypeclasses=None
    poltypeclasses=None
    allin=None
    prmlinesplit=prmline.split()
    # STEP 1
    typenums=GrabTypesFromPrmLine(poltype,prmlinesplit)
    # STEP 2
    for listofpoltypeclasses in poltypeclasseslist:
        # STEP 3
        for poltypeclasses in listofpoltypeclasses:
            # STEP 4
            allin=True
            if int(typenums[0])!=poltypeclasses[0]: # then try flipping
                poltypeclasses=poltypeclasses[::-1]
            for i in range(len(typenums)):
                poltypeclass=int(typenums[i])
                otherpoltypeclass=poltypeclasses[i]
                if poltypeclass!=otherpoltypeclass:
                    allin=False 
            # STEP 5
            if allin==True:
                return listofpoltypeclasses,poltypeclasses
    # STEP 6
    if allin!=None:
        if allin==False:
            listofpoltypeclasses=None
            poltypeclasses=None

    return listofpoltypeclasses,poltypeclasses


def MapParameterLineToTransferInfo(poltype,prms,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,poltypeclassestosmartsatomordersext,newpoltypeclassestocomments,newpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses,defaultvalues=None,keyword=None):
    """
    Intent: Want to add comments to key file when adding parameters from database so user can read and understand where parameters came from.  
    Input: Array of parameter lines, dictionary of poltypes types -> amoeba09 SMARTS + atom orders, dictionary of poltypes types -> SMARTS (that matches input molecule) + atom orders, dictionary of poltype types -> element + amoeba09 tinker description, dictionary of poltype types -> SMARTS (from external database of SMARTS) + atom order, dictionary of poltype types -> comments (from amoeba21 database), dictionary of poltypes types -> SMARTS (from amoeba21 database), dictionary of torsion parameter line -> comment (from missing torsions transfering from alkane/benzene), array of missing vdw types, array of missing torsions, array of missing bonds types, array of missing angle types, keyword boolean to specify if using default values, keyword boolean to specify if transfering comments for tor-tor 
    Output: Dictionary of paramteter line to add to key -> array of comment lines to add above parameter line
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of parameter lines
    2. Check if parameter line is in dictionary of torsion parameter line -> comment (from missing torsions transfering from alkane/benzene), if so extract the comments from dictionary
    3. If not, then check if parameter line is in amoeba09 SMARTS dictionary
    4. If not, then check if parameter line is in external SMARTS dictionary
    5. If not, then check if parameter line is in amoeba21 dictionary
    6. Construct comments appropriate to which database the parameters came from.
    7. If default values are used for bond/angle/opbend, then add comment for warning in key file
    8. Check for missing vdw/bond/angle/strbnd/torsion parameters and add comment in key file that parameters are missing from database 
    """
    prmstotransferinfo={}
    # STEP 1
    for line in prms:
        warn=False
        # STEP 2
        if keyword!=None:
           if keyword not in line:
               transferinfo='# blank'
           else:
               poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())
               smartsatomorders=poltypeclassestosmartsatomordersext[poltypeclasses]
               transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' from external database'+'\n'
        else:
           if line in arotorsionlinetodescrips.keys():
               tup=arotorsionlinetodescrips[line]
               extra=tup[0]
               descrips=tup[1]
               transferinfoline=extra+'# Transferring from '+str(descrips)+'\n' 
           else:
               # STEP 3 
               poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestoparametersmartsatomorders.keys())
               if poltypeclasses==None:
                   # STEP 4    
                   poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())
                   # STEP 5
                   if poltypeclasses==None:
                       poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,newpoltypeclassestocomments.keys())
                       # STEP 6
                       comments=newpoltypeclassestocomments[poltypeclasses]
                       smartslist=newpoltypeclassestosmartslist[poltypeclasses]
                       commentstring=' '.join(comments)
                       smartsliststring=' '.join(smartslist)
                       transferinfoline='# amoeba21'+' '+'comments='+commentstring+' '+'SMARTS match = '+smartsliststring+'\n'
                   else:
                       smartsatomorders=poltypeclassestosmartsatomordersext[poltypeclasses]
                       transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' from external database'+'\n'

               else:
                   # STEP 6
                   parametersmartsatomorders=poltypeclassestoparametersmartsatomorders[poltypeclasses]
                   smartsatomorders=poltypeclassestosmartsatomorders[poltypeclasses]
                   elementtinkerdescrips=poltypeclassestoelementtinkerdescrips[poltypeclasses]
                   transferinfoline='# amoeba09'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' '+'to SMARTS from parameter file'+' '+str(parametersmartsatomorders)+' '+'with tinker type descriptions '+str(elementtinkerdescrips)+'\n'
                   # STEP 7
                   if defaultvalues!=None:
                       for value in defaultvalues:
                           if str(value) in line:
                               warn=True
                   if warn==True:
                       transferinfoline+='# '+'WARNING DEFAULT MM3 OPBEND VALUES USED '+'\n'
                   smarts=smartsatomorders[0]
                   if '~' in smarts or '*' in smarts:
                       transferinfoline+='# '+'WARNING WILDCARDS USED IN SMARTS PARAMETER MATCHING'+'\n'
                       warn=True
        # STEP 8
        showtransferinfo=True
        extraline=''
        if 'vdw' in line:
            linesplit=line.split()
            vdwtype=int(linesplit[1])
            if vdwtype in missingvdwtypes:
                warn=True
                extraline+='# Missing vdw parameters'+'\n'

        if 'torsion' in line:
            linesplit=line.split()
            a=int(linesplit[1])
            b=int(linesplit[2])
            c=int(linesplit[3])
            d=int(linesplit[4])
            ls=[a,b,c,d]
            if ls in torsionsmissing or ls[::-1] in torsionsmissing:
                extraline+='# Missing torsion parameters, will attempt to fit parameters'+'\n'
                showtransferinfo=True


        transferinfoline+=extraline
        if showtransferinfo==False:
            transferinfoline=extraline
        if warn==True:
            poltype.WriteToLog(transferinfoline)
            warnings.warn(transferinfoline)
        prmstotransferinfo[line]=transferinfoline

    return prmstotransferinfo

def FindMaximumCommonSubstructures(poltype,parametersmartslist,rdkitmol):
    """
    Intent: Reduce the amount of amoeba09 SMARTS matching, by only trying amoeba09 SMARTS that have common substructure with input molecule. 
    Input: List of amoeba09 SMARTS, rdkit mol object
    Output: Dictionary of amoeba09 SMARTS -> maximum common substructure, maximum number of atoms from common substructure 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over amoeba09 SMARTS array
    2. Compute max common substructure between mol and amoeba09 SMARTS
    3. Compute atom size and extract the SMARTS from max common substructure
    4. If there is a common substructure, save the substructure SMARTS and amoeba09 SMARTS in a dictionary
    """
    parametersmartstomaxcommonsubstructure={}
    maxatomsize=0
    # STEP 1
    for parametersmarts in parametersmartslist:
        # STEP 2
        mols = [rdkitmol,Chem.MolFromSmarts(parametersmarts)]
        res=rdFMCS.FindMCS(mols)
        # STEP 3
        atomnum=res.numAtoms
        smartsmcs=res.smartsString
        if atomnum>0:
            if atomnum>maxatomsize:
                maxatomsize=atomnum
            # STEP 4
            parametersmartstomaxcommonsubstructure[parametersmarts]=smartsmcs
    return parametersmartstomaxcommonsubstructure,maxatomsize 

def GrabPlanarBonds(poltype,listofbondsforprm,mol): # used for checking missing opbend parameters later
    """
    Intent: Need to know which bonds are planar for finding opbend parameters
    Input: Array of bonds in mol, mol object
    Output: Array of planar bonds
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of bonds in mol
    2. Check if the bond (or its reverse) has last atom as SP2 and connected to 3 atoms (or if in aromatic ring)
    3. If so then append to array of planar bonds
    """
    planarbonds=[]
    # STEP 1
    for bond in listofbondsforprm:
        revbond=bond[::-1]
        totalbonds=[bond,revbond]
        # STEP 2
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
                        # STEP 3
                        planarbonds.append(bond)
                else:
                    if b.GetHyb()==2 and len(list(openbabel.OBAtomAtomIter(b)))==3:
                        if bond not in planarbonds:
                            # STEP 3
                            planarbonds.append(bond)

            else:
                if b.GetHyb()==2 and len(list(openbabel.OBAtomAtomIter(b)))==3:
                    if bond not in planarbonds:
                        # STEP 3
                        planarbonds.append(bond)
    return planarbonds 


def FindPotentialMissingParameterTypes(poltype,prms,tinkerclassestopoltypeclasses):
    """
    Intent: Want to find missing opbend parameter types by first searching for opbend parameter types in dictionary of planar bonds -> poltype types
    Input: Array of opbend parameter lines, dictionary of tinker classes -> poltype types for planar bonds
    Output: Array of potential missing opbend parameters (need further filtering)
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Get array of every poltype types array from array poltype types (flattening list)
    2. Iterate over array of opbend parameter lines
    3. Search for poltype types array in and see if there is match to parameter line
    4. If there is a match append to array of parameters found
    5. Iterate over poltypes types arrays if they dont exist in found parameter array, then append to missing parameter array 
    """
    # STEP 1
    poltypeclasseslist=list(tinkerclassestopoltypeclasses.values())
    newpoltypeclasseslist=[]
    for poltypeclassesls in poltypeclasseslist:
        newpoltypeclassesls=[]
        for poltypeclasses in poltypeclassesls:
            newpoltypeclassesls.append(poltypeclasses)
        newpoltypeclasseslist.append(newpoltypeclassesls)
    missingprms=[]
    foundprms=[]
    # STEP 2
    for line in prms:
        # STEP 3
        poltypeclassesls,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,newpoltypeclasseslist)
        linesplit=line.split()
        prms=[int(linesplit[1]),int(linesplit[2])]
        if poltypeclassesls!=None:
            # STEP 4
            for poltypeclasses in poltypeclassesls:
                poltypeclasses=tuple([int(i) for i in poltypeclasses])
                if poltypeclasses[0]==prms[0] and poltypeclasses[1]==prms[1]:
                    foundprms.append(poltypeclasses) 
    # STEP 5
    for poltypeclassesls in newpoltypeclasseslist:
        for poltypeclasses in poltypeclassesls:
            if poltypeclasses not in foundprms and poltypeclasses not in missingprms:
                missingprms.append(poltypeclasses)
 
    return missingprms 

def ConvertPoltypeClassesToIndices(poltype,missingprmtypes):
    """
    Intent: Function to convert poltypes types -> back into indices 
    Input: Array of missing parameter types
    Output: Array of indices matching poltype types  
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of poltype type arrays
    2. Iterate over array of poltype types
    3. Grab all possible atomic indices from poltype type number
    4. Append indices to array of indices 
    5. Generate all combinations of indices (multiple occurances of types in different parts of molecule)
    """
    missingprmindices=[]
    # STEP 1
    for poltypeclasses in missingprmtypes:
        indices=[]
        # STEP 2
        for poltypeclass in poltypeclasses:
            # STEP 3
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
            keylist=[i-1 for i in keylist]
            indices.append(keylist)
        # STEP 4
        missingprmindices.append(indices)
    # STEP 5
    finallist=[]
    for indices in missingprmindices:
        combs = list(itertools.product(*indices))
        for comb in combs:
            finallist.append(comb)
    return finallist

def FilterIndices(poltype,potentialmissingprmindices,indices):
    """
    Intent: Remove potential missing opbend parameters unless they have a planar bond
    Input: Array of bond indices for opbend, array of planar bond indices 
    Output: Modified array of bond indices for opbend
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of opbend bond indices
    2. Check if indices are in planar indices, if so append to new array
    """
    missingprmindices=[]
    # STEP 1
    for prmindices in potentialmissingprmindices:
        # STEP 2
        if list(prmindices) in indices:
            missingprmindices.append(prmindices)
    return missingprmindices
        

def DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools):
    """
    Intent: For missing opbend parameters, assign default MM3 values
    Input: Array of missing opbend parameter indices, mol object, dictionary of bond indices -> boolean if its trigonal center or not
    Output: Array of opbend parameters, default values assigned used later for adding comments in key 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of missing bond indices
    2. Compute poltype types from current bond indices
    3. Check if the b in a-b is trigonal center
    4. If so, then if check if 3 neighbors (might be redundant) 
    5. If carbon bonded to double bonded oxygen (carbonyl), then assign corresponding MM3 parameters
    6. If nitro group, assign corresponding MM3 parameters
    7. Otherwise assign very generic MM3 parameters for all other cases
    """
    newopbendprms=[]
    defaultvalues=[]
    # STEP 1
    for opbendprmindices in missingopbendprmindices:
        # STEP 2
        babelindices=[i+1 for i in opbendprmindices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        boolarray=opbendbondindicestotrigonalcenterbools[opbendprmindices]
        # STEP 3
        if boolarray[1]==True:
            atomindex=int(babelindices[1])
            atom=mol.GetAtom(atomindex)
            neighbs=list(openbabel.OBAtomAtomIter(atom))
            # STEP 4
            if len(neighbs)==3:
                atomnum=atom.GetAtomicNum()
                # STEP 5
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
                # STEP 6
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
                # STEP 7
                else:
                    opbendvalue=round(poltype.defopbendval*71.94,2)

                defaultvalues.append(opbendvalue)
                firstprmline='opbend'+' '+str(poltypeclasses[0])+' '+str(poltypeclasses[1])+' '+'0'+' '+'0'+' '+str(opbendvalue)+'\n'
                newopbendprms.append(firstprmline)

    return newopbendprms,defaultvalues


def WriteOutList(poltype,ls,filename):
    """
    Intent: Write out array to output file
    Input: Array, output filename
    Output: -
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each item in list
    2. Write out item to filename as new line
    """
    with open(filename, 'w') as filehandle:
        # STEP 1
        for listitem in ls:
            # STEP 2
            filehandle.write('%s\n' % listitem)

def ReadTorsionList(poltype,filename):
    """
    Intent: Read back torsion list from output file
    Input: filename containing missing torsions 
    Output: Array of missing torsoins
    Referenced By: poltype.py - GenerateParameters
    Description:
    1. Iterate over lines of file
    2. Remove brackets and split string by comma into array
    3. Convert items in array to integer and append list to list of list
    """
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            # STEP 1
            for line in filehandle:
                # STEP 2
                current = line[:-1]
                current=current.replace('[','').replace(']','')
                current=current.split(',')
                # STEP 3
                current=[int(i) for i in current]
                newls.append(current)

    return newls


def ReadTorTorList(poltype,filename):
    """
    Intent: Read tor-tor from output filename into array
    Input: Output filename
    Output: Array of missing tor-tor 
    Referenced By: poltype.py - GenerateParameters
    Description:
    1. Iterate over lines of file
    2. Remove brackets and split string by comma into array
    3. Convert items in array to integer and append list to list of list
    """
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
    """
    Intent: Read missing vdw from filename
    Input: Output filename of missing vdw atom indices
    Output: Array of missing vdw atom indices
    Referenced By: poltype.py - GenerateParameters
    Description:
    1. Iterate over lines in file
    2. Convert string to integer and append to array
    """
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            # STEP 1
            for line in filehandle:
                current = line[:-1]
                # STEP 2
                newls.append(int(current))

    return newls


def CheckIfParametersExist(poltype,potentialmissingindices,prms):
    """
    Intent: If a given opbend indices are not in given opbend parameters then it must be missing.
    Input: Array of potential missing opbend indices, array of parameter lines
    Output: Array of missing opbend parameter indices  
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of missing opbend indices
    2. Convert indices -> poltype types
    3. Iterate over lines in input parameter lines
    4. Grab the types from parameter lines
    5. If there is a match between types, then the parameters are determined to be found
    6. If not found, append indices to array of missing opbend parameter indices
    """
    missingprmindices=[]
    # STEP 1
    for indices in potentialmissingindices:
        # STEP 2
        babelindices=[i+1 for i in indices]
        symtypes=[poltype.idxtosymclass[i] for i in babelindices]
        found=False
        symtypes=[str(i) for i in symtypes]
        # STEP 3
        for prmline in prms:
            linesplit=prmline.split()
            # STEP 4
            parms=[linesplit[1],linesplit[2]]
            # STEP 5
            if parms[0]==symtypes[0] and parms[1]==symtypes[1]:
               found=True
        # STEP 6
        if found==False:
            missingprmindices.append(indices)
    return missingprmindices
                
def WriteDictionaryToFile(poltype,dic,filename):
    """
    Intent: Write contents of dictionary (such as torsion types -> parameters for non-aromatic torsion ring refinement) to disk
    Input: Dictoinary, filename
    Output: - 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Dump dictionary to filename
    """
    # STEP 1
    json.dump(dic, open(filename,'w'))        

def ReadDictionaryFromFile(poltype,filename):
    """
    Intent: Read dictionary of torsion -> torsion parameter guess for non-aromatic ring refinement
    Input: Filename to read in dictionary values
    Output: Dictionary of torsion -> parameter transfer from database
    Referenced By: poltype.py - GenerateParameters
    Description: 
    1. Load dictionary from file
    """
    return json.load(open(filename))



def ReadExternalDatabase(poltype):
    """
    Intent: Read database of SMARTS for all parameters (supersedes all other databases)
    Input: -
    Output: For each parameter type, a dictionary of SMARTS + atom order in SMARTS -> parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over lines of database file
    2. Determine if in first part of database (logical SMARTS only) or second part SMILES + atom order only
    3. Determine the keyword (which parameter type) and also the SMARTS + atom order and the parameters
    4. Based on keyword put SMARTS + atom order and parameters in appropriate dictionary, tor-tor requires extra formatting due to the nature of how syntax is when you input it into the keyfile
    """
    temp=open(poltype.externalparameterdatabase,'r')
    results=temp.readlines()
    temp.close()
    bondsmartsatomordertoparameters={} 
    anglesmartsatomordertoparameters={}
    strbndsmartsatomordertoparameters={}
    torsionsmartsatomordertoparameters={}
    tortorsmartsatomordertoparameters={}
    tortorsmartsatomordertogrid={}
    smartsatomordertotorvdwdb={}
    opbendsmartsatomordertoparameters={}
    vdwsmartsatomordertoparameters={}
    founddelim=False
    # STEP 1
    for line in results:
        linesplit=line.split()
        if len(linesplit)==0:
            continue
        # STEP 2
        if 'DELIM' in line:
            founddelim=True
        if linesplit[0]=='#':
            continue
        # STEP 3
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
        # STEP 4
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
        if founddelim==True:
            smartsatomordertotorvdwdb[smartsatomorder]=True
        else:
            smartsatomordertotorvdwdb[smartsatomorder]=False
    return bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters,tortorsmartsatomordertoparameters,tortorsmartsatomordertogrid,smartsatomordertotorvdwdb


def ConvertToPoltypeClasses(poltype,torsionsmissing):
    """
    Intent: Convert missing torsion in atom indices to atom types instead 
    Input: Array of missing torsions (in atomic indices)
    Output: Array of missing torsions (in atomic types) 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over torsions in array of missing torsions
    2. Use canonical algorithm for determining a and d for representing the torsions around bond b-c
    3. Convert atomic indices to atomic types and append to array of missing torsion types 
    """
    newtorsionsmissing=[]
    # STEP 1
    for sublist in torsionsmissing:
        newsublist=[i+1 for i in sublist]
        a,b,c,d=newsublist[:]
        batom=poltype.mol.GetAtom(b)
        catom=poltype.mol.GetAtom(c)
        # STEP 2
        aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,batom,catom)
        newa=aatom.GetIdx()
        newd=datom.GetIdx()
        # STEP 3
        sorttor=torfit.sorttorsion(poltype,[poltype.idxtosymclass[a],poltype.idxtosymclass[b],poltype.idxtosymclass[c],poltype.idxtosymclass[d]])
        newsorttor=torfit.sorttorsion(poltype,[poltype.idxtosymclass[newa],poltype.idxtosymclass[b],poltype.idxtosymclass[c],poltype.idxtosymclass[newd]])

        if sorttor not in newtorsionsmissing:
            newtorsionsmissing.append(sorttor)
        if newsorttor not in newtorsionsmissing:
            newtorsionsmissing.append(newsorttor)

    return newtorsionsmissing



def GrabKeysFromValue(poltype,dic,thevalue):
    """
    Intent: Given a dictionary value such as type number, find all corresponding keys or atom indices that map to that type number.
    Input: Dictionary, value of interest
    Output: List of keys corresponding to value of interest 
    Referenced By: Many functions
    Description:
    1. Iterate over dictionary 
    2. If current value equals value of interest, then append current key to list
    """
    keylist=[]
    # STEP 1
    for key,value in dic.items():
        # STEP 2
        if value==thevalue:
            keylist.append(key)
    return keylist



def GrabAllPossibleMoleculeIndices(poltype,typeindices):
    """
    Intent: Grab all possible atomic indices that map to input type numbers
    Input: Array of type numbers
    Output: Array of atomic indices that map to input type numbers
    Referenced By: MatchExternalSMARTSToMolecule
    Description:
    1. Iterate over each type number in input type array
    2. Grab all atomic indices that correspond to the type number
    3. Append to array
    4. Generate all combinations of indices that map to type numbers
    """
    thelist=[]
    # STEP 1
    for typenum in typeindices:
        # STEP 2
        allindices=GrabKeysFromValue(poltype,poltype.idxtosymclass,typenum)
        # STEP 3
        thelist.append(allindices)
    # STEP 4
    listofmoleculeindices=list(itertools.product(*thelist))


    return listofmoleculeindices


def MatchExternalSMARTSToMolecule(poltype,rdkitmol,smartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb,torsionindicestoextsmarts=None):
    """
    Intent: Take database of external SMARTS for all parameter types and attempt to match to input molecule
    Input: Rdkit mol object, dictionary of SMARTS + atom order -> parameters, dictionary of atom index -> neighboring atom indices, dictionary of SMARTS + atom order -> boolean if came from SMILES database or logical SMARTS database (in same file),  
    Output: Dictionary of atomic indices (atom/bond/angle/torsion) -> SMARTS length , Dictionary of atomic indices -> SMARTS, Dictionary of atomic indices -> SMARTS + atom order , Dictionary of SMARTS + atom order -> parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. If already matched torsion external SMARTS and trying to match vdw SMARTS, only match vdw SMARTS consistent with vdw from torsion SMARTS (make an array of SMARTS that only has matches from torsion parameters)
    2. Iterate over dictionary of SMARTS+atom orders
    3. Try to match only vdw from torsion matches in DrugBank database set
    4. If logical SMARTS (above the database from DrugBank in file), then if not exact match, skip (where as from DrugBank molecules dont need exact match, only torsion+neighbors etc)
    5. Try to match the maximum common substructure between DrugBank SMILES and input molecule to input molecule and generate map from index in max common substructure SMARTS -> atom index in input molecule
    6. Based on that map and input atom order asscoiated with DrugBank smiles, determine if all atom indices in input molecule (for torsion or vdw) exist in that match
    7. Now determine if all neighbors exist in the SMARTS match also, otherwise its not "transferable" enough
    8. If SMARTS (logical smarts exact match) or SMILES from DrugBank matches input indices+neighbors, then append match and parameters and SMARTS to dictionaries 
    9. Organize SMARTS in dictionaries around rotatable bond (if multiple matches were found), then want to choose matches that are largest == more transferable
    10. Choose best smarts match by longest length
    11. Input the best SMARTS match into dictionaries, if another atomic indices pair has same type, then also add matches for those indices too  
    """
    indicestoextsmartsmatchlength={}
    indicestoextsmarts={}
    indicestoextsmartsatomorder={}
    indicestosmartslist={}
    indicestomatchlist={}
    indicestoatomorderlist={}
    indicestosmartsatomorder={}
    indicestoprmlist={}
    restrictedsmarts=[]
    failedmatches=[]
    if poltype.quickdatabasesearch==False: # debugging keyword 
        # STEP 1
        if torsionindicestoextsmarts!=None:
            for torsionindices,extsmarts in torsionindicestoextsmarts.items():
                if extsmarts not in restrictedsmarts:
                    restrictedsmarts.append(extsmarts)
        listsmartsatomordertoparameters=list(smartsatomordertoparameters.keys())
        # STEP 2
        for k in tqdm(range(len(listsmartsatomordertoparameters)),desc='Searching SMARTS database'):
            smartsatomorder=listsmartsatomordertoparameters[k]
            parameters=smartsatomordertoparameters[smartsatomorder]
            smarts=smartsatomorder[0]
            fromtorvdwdb=smartsatomordertotorvdwdb[smartsatomorder]
            atomorderlist=smartsatomorder[1]
            # STEP 3
            if len(restrictedsmarts)!=0:
                if smarts not in restrictedsmarts and fromtorvdwdb==True: # try to keep vdw matches consistent with torsion matches in DrugBank testing set
                    continue
            if (len(restrictedsmarts)==0 and len(atomorderlist)==1) and fromtorvdwdb==True:
                continue 
            substructure = Chem.MolFromSmarts(smarts)
            mols = [rdkitmol,substructure]
            res=rdFMCS.FindMCS(mols)
            atomnum=res.numAtoms
            smartsmcs=res.smartsString
            diditmatch=False
            diditmatchexactly=rdkitmol.HasSubstructMatch(substructure)
            # STEP 4
            if fromtorvdwdb==False and diditmatchexactly==False:
                continue
            if atomnum>=len(atomorderlist):
                # STEP 5
                mcssubstructure = Chem.MolFromSmarts(smartsmcs)
                mcsmatches=substructure.GetSubstructMatches(mcssubstructure)
                if len(mcsmatches)>0:
                    # STEP 6
                    mcsmatch=mcsmatches[0]
                    mcsindices=list(range(len(mcsmatch)))
                    mcssmartsindextosmartsmolindex=dict(zip(mcsindices,mcsmatch)) 
                    smartsmolindexlist=[i-1 for i in atomorderlist]
                    foundindices=True
                    
                    for i in smartsmolindexlist:
                        if i not in mcssmartsindextosmartsmolindex.values():
                            foundindices=False
                    if foundindices==True:
                        # STEP 7
                        matches=rdkitmol.GetSubstructMatches(mcssubstructure)
                        if len(matches)>0:
                            thematch=matches[0]
                            theindices=list(range(len(thematch)))
                            mcssmartsindextomolindex=dict(zip(theindices,thematch)) 
                            smartsmolindextomcssmartsindex={v: k for k, v in mcssmartsindextosmartsmolindex.items()} 
                            mcssmartsindices=[smartsmolindextomcssmartsindex[i] for i in smartsmolindexlist]
                            indices=[mcssmartsindextomolindex[i] for i in mcssmartsindices]
                            mcsorder=[i+1 for i in indices]
                            
                            for match in matches:
                                allin=True
                                for idx in indices:
                                    if idx not in match:
                                        allin=False
                                if allin==True:
                                    matcharray=match
                                    break
                            nindexes=[]
                            for idx in indices: 
                                neighborindexes=indextoneighbidxs[idx]
                                for nidx in neighborindexes:
                                    if nidx not in nindexes:
                                        nindexes.append(nidx)
                                    if len(indices)==1 and len(neighborindexes)==1:
                                        nneighborindexes=indextoneighbidxs[nidx]
                                        for nnidx in nneighborindexes:
                                            if nnidx not in nindexes:
                                                nindexes.append(nnidx)
            
                            diditmatch=CheckIfNeighborsExistInSMARTMatch(poltype,nindexes,matcharray)  
                            
                            if fromtorvdwdb==False:
                                diditmatch=True # only care about neighbor matching for large SMARTS from part of DB genreated by fragmenter, top of DB can have shorter SMARTS 
                            if diditmatch==True:
                                moleculeindices=tuple(indices)
                                if len(indices)==4: # then torsion matches need to be consistent for all torsion around bond
                                    b=indices[1]
                                    c=indices[2]
                                    batombabel=poltype.mol.GetAtom(b+1)
                                    catombabel=poltype.mol.GetAtom(c+1)
                                    bisinring=batombabel.IsInRing()
                                    cisinring=catombabel.IsInRing()
                                    bhyb=batombabel.GetHyb()
                                    chyb=catombabel.GetHyb()
                                    if (bisinring==True and cisinring==True) and (bhyb==2 and chyb==2): # dont transfer to SP2 torsion from database, let this case be handled by benzene transfer (easier to understand in output file also)
                                        continue
                                # STEP 8
                                if moleculeindices not in indicestosmartslist.keys():
                                    indicestosmartslist[moleculeindices]=[]
                                    indicestomatchlist[moleculeindices]=[]
                                    indicestoatomorderlist[moleculeindices]=[]
                                    indicestoprmlist[moleculeindices]=[]
                                    indicestosmartsatomorder[moleculeindices]=[]

                                if smartsmcs not in indicestosmartslist[moleculeindices]:
                                    mcssmartsatomorder=tuple([smartsmcs,tuple(mcsorder)])
                                    indicestosmartslist[moleculeindices].append(smartsmcs)
                                    indicestomatchlist[moleculeindices].append(matcharray)
                                    indicestoatomorderlist[moleculeindices].append(mcssmartsatomorder)
                                    indicestosmartsatomorder[moleculeindices].append(smartsatomorder)
                                    indicestoprmlist[moleculeindices].append(parameters)
            else:
                failedmatches.append(smarts)
        # STEP 9
        rotatablebondtosmartslist={}
        for moleculeindices,smartslist in indicestosmartslist.items():
            if len(moleculeindices)==4: 
                b=moleculeindices[1]
                c=moleculeindices[2]
                if b<c:
                    rotbnd=[b,c]
                elif b>c:
                    rotbnd=[c,b]
                else:
                    rotbnd=[b,c] 
                rotbnd=tuple(rotbnd)
                if rotbnd not in rotatablebondtosmartslist.keys():
                    rotatablebondtosmartslist[rotbnd]=[]
                for smarts in smartslist:
                    if smarts not in rotatablebondtosmartslist[rotbnd]:
                        rotatablebondtosmartslist[rotbnd].append(smarts)
        rotatablebondtotruesmartslist={}
        for rotbnd,allsmartslist in rotatablebondtosmartslist.items(): # this section appears to be redundant
            for smarts in allsmartslist:
                allin=True
                for moleculeindices,smartslist in indicestosmartslist.items():
                    if len(moleculeindices)==4: 
                       b=moleculeindices[1]
                       c=moleculeindices[2]
                       if b<c:
                           otherrotbnd=[b,c]
                       elif b>c:
                           otherrotbnd=[c,b]
                       else:
                           otherrotbnd=[b,c] 
                       otherrotbnd=tuple(otherrotbnd)
                       if rotbnd==otherrotbnd: # if same rotbond and found other match not in  
                           if smarts not in smartslist:
                               allin=False
                if allin==True: 
                    if rotbnd not in rotatablebondtotruesmartslist.keys():
                        rotatablebondtotruesmartslist[rotbnd]=[]
                    if smarts not in rotatablebondtotruesmartslist[rotbnd]:
                        rotatablebondtotruesmartslist[rotbnd].append(smarts)
        # STEP 10
        rotatablebondtotruesmarts={}
        for rotbnd,smartslist in rotatablebondtotruesmartslist.items():
            smartslistlen=[len(i) for i in smartslist]
            maxlen=max(smartslistlen)
            maxidx=smartslistlen.index(maxlen)
            maxsmarts=smartslist[maxidx]
            rotatablebondtotruesmarts[rotbnd]=maxsmarts

        for moleculeindices,smartslist in indicestosmartslist.items():
           if len(moleculeindices)==4: 
               b=moleculeindices[1]
               c=moleculeindices[2]
               if b<c:
                   rotbnd=[b,c]
               elif b>c:
                   rotbnd=[c,b]
               else:
                   rotbnd=[b,c] 
               rotbnd=tuple(rotbnd)
               if rotbnd in rotatablebondtotruesmarts.keys():
                   maxsmarts=rotatablebondtotruesmarts[rotbnd]
               else:
                   smartslistlen=[len(i) for i in smartslist]
                   maxlen=max(smartslistlen)
                   maxidx=smartslistlen.index(maxlen)
                   maxsmarts=smartslist[maxidx]
          
        
           else:
               smartslistlen=[len(i) for i in smartslist]
               maxlen=max(smartslistlen)
               maxidx=smartslistlen.index(maxlen)
               maxsmarts=smartslist[maxidx]
           for j in range(len(indicestosmartsatomorder[moleculeindices])):
               originalsmartsatomorder=indicestosmartsatomorder[moleculeindices][j]
               thesmarts=indicestosmartslist[moleculeindices][j]
               if thesmarts==maxsmarts: 
                   break

           # STEP 11
           smarts=maxsmarts
           matcharray=indicestomatchlist[moleculeindices][j]
           smartsatomorder=indicestoatomorderlist[moleculeindices][j]
           prms=indicestoprmlist[moleculeindices][j]
           babelindices=[i+1 for i in moleculeindices]
           typeindices=[poltype.idxtosymclass[i] for i in babelindices]
           listofmoleculeindices=GrabAllPossibleMoleculeIndices(poltype,typeindices) 
           for babelindices in listofmoleculeindices:
               moleculeindices=tuple([i-1 for i in babelindices])
               if moleculeindices not in indicestoextsmarts.keys() and moleculeindices[::-1] not in indicestoextsmarts.keys():
                   indicestoextsmartsmatchlength[moleculeindices]=len(matcharray)
                   indicestoextsmarts[moleculeindices]=smarts
                   indicestoextsmartsatomorder[moleculeindices]=smartsatomorder
                   smartsatomordertoparameters[smartsatomorder]=prms
    return indicestoextsmartsmatchlength,indicestoextsmarts,indicestoextsmartsatomorder,smartsatomordertoparameters



def CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,indicestoextsmartsmatchlength,indicesforprmtoparametersmarts,indicesforprmtosmarts,indicestoextsmarts,indicesforprmtomatchallneighbs,indicestoextsmartsatomorders):
    """
    Intent: If have amoeba09 match and match from external smarts, then just use match from external database
    Input: dictionary of atomic indices -> SMARTS match length, dictionary of atomic indices -> amoeba09 SMILES , dictionary of atomic indices -> SMARTS that matches to input molecule and amoeba09 SMILES , dicionary of atomic indices -> SMARTS matches from external database, dictionary of atomic indices for amoeba09 matches -> boolean if it matched all neighbors , dictionary of atomic indices -> external SMARTS database matches + atom order
    Output: modified dictionary of atomic indices -> amoeba09 SMILES  , modified dictionary of atomic indices -> SMARTS that matches to input molecule and amoeba09 SMILES , modified dicionary of atomic indices -> SMARTS matches from external database , modified dictionary of atomic indices -> external SMARTS database matches + atom order
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dicionary of atomic indices -> external SMARTS match length
    2. If there is also a match in amoeba09 database for those atomic indices
    3. Remove indices from amoeba09 matches dictionary 
    """
    newindicestoextsmarts={}
    newindicestoextsmartsatomorder={}
    # STEP 1
    for indices,extsmartsmatchlength in indicestoextsmartsmatchlength.items():
        # STEP 2
        if indices in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices]
            smarts=smartsls[0]
            if indices in indicestoextsmartsatomorders.keys():
                smartsatomorder=indicestoextsmartsatomorders[indices]
                extsmarts=indicestoextsmarts[indices]
            else:
                smartsatomorder=indicestoextsmartsatomorders[indices[::-1]]
                extsmarts=indicestoextsmarts[indices[::-1]]
            poormatch=True
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            # STEP 3
            if extsmartsmatchlength>substructurenumatoms or poormatch==True:
                del indicesforprmtosmarts[indices]
                del indicesforprmtoparametersmarts[indices]
                newindicestoextsmarts[indices]=extsmarts
                newindicestoextsmartsatomorder[indices]=smartsatomorder
        elif indices[::-1] in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices[::-1]]
            smarts=smartsls[0]
            # STEP 3
            if indices[::-1] in indicestoextsmartsatomorders.keys():
                smartsatomorder=indicestoextsmartsatomorders[indices[::-1]]
                extsmarts=indicestoextsmarts[indices[::-1]]
            else:
                smartsatomorder=indicestoextsmartsatomorders[indices]
                extsmarts=indicestoextsmarts[indices]

            poormatch=True
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms:
                del indicesforprmtosmarts[indices[::-1]]
                del indicesforprmtoparametersmarts[indices[::-1]]
                 
                newindicestoextsmarts[indices[::-1]]=extsmarts
                newindicestoextsmartsatomorder[indices[::-1]]=smartsatomorder
    return indicesforprmtoparametersmarts,indicesforprmtosmarts,newindicestoextsmarts,newindicestoextsmartsatomorder

def AddExternalDatabaseSMARTSMatchParameters(poltype,prms,indicestoextsmarts,smartsatomordertoparameters,keyword,indicestoextsmartsatomorder,smartsatomordertogrid=None):
    """
    Intent: If match exists for atomic indices from external SMARTS database, then take those parameters and update them in array of parameter lines (previous array of parameter lines would have been from amoeba09)
    Input: Array of parameter lines, dictionary of atomic indices -> SMARTS from external database, dictionary of SMARTS + atom order -> parameters , keyword (need to add extra info depending on if opbend, torsion or tor-tors matches), dictionary of atom ic indices to SMARTS + atom order from external database, dictionary of SMARTS + atom order -> dihedral angle energy grid (for tor-tor) 
    Output: Modified array of parameter lines , dictionary of poltype types -> SMARTS + atom order from external database
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of indices -> SMARTS matches
    2. Iterate over dictionary of SMARTS + atom order -> parameters
    3. If same SMARTS and atoms matching in SMARTS then 
    4. Construct parameter line based on keyword (tinker has specific format for each keyword)
    5. Add modified parameter line to array of parameter lines, also append match information to dictionary of poltype types -> SMARTS + atom order from external database
    """
    poltypeclassestosmartsatomordersext={}
    # STEP 1
    for indices,extsmarts in indicestoextsmarts.items():
        thesmartsatomorder=indicestoextsmartsatomorder[indices]
        theatomorder=thesmartsatomorder[1]
        # STEP 2
        for smartsatomorder,parameters in smartsatomordertoparameters.items():
            smarts=smartsatomorder[0]
            atomorder=smartsatomorder[1]
            # STEP 3       
            if smarts==extsmarts and atomorder==theatomorder:
                # STEP 4
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
                elif keyword=='torsion':
                    for prmidx in range(len(parameters)):
                        prm=parameters[prmidx]
                        if prmidx==1:
                            phase=180
                        else:
                            phase=0
                        fold=prmidx+1
                        line+=str(prm)+' '+str(phase)+' '+str(fold)+' ' 
                    line+='\n' 

                else:
               
                    for prm in parameters:
                        line+=str(prm)+' '
                    line+='\n'
                # STEP 5 
                prms.append(line)

                poltypeclassestosmartsatomordersext[tuple([tuple(poltypeclasses)])]=smartsatomorder
    return prms,poltypeclassestosmartsatomordersext


def ExtractOpbendFromBond(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses):
    """
    Intent: Convert dictionaries of bond indices -> SMARTS atom order to dictionaries of opbend indices -> SMARTS atom order and opbend indices -> tinker classes
    Input: Dictionaries of bond indices -> SMARTS atom order, dictionary of bond indices -> tinker classes
    Output: opbebd indices -> tinker classes, opbend indices -> SMARTS + atom orders 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over dictionary of bond indices -> SMARTS + atom orders
    2. Extract tinker classes to add to dictionary of indices -> tinker classes 
    3. Opbend is directionally important depending on which atom is trigonal center , so add bond indices + reverse and will check later which one is needed for opbend (second atom should be trigonal center)
    """
    opbendbondindicestosmartsatomorders={}
    opbendbondindicestotinkerclasses={}
    # STEP 1
    for bondindices,smartsatomorders in bondindicestosmartsatomorders.items(): 
        smarts=smartsatomorders[0]
        revbondindices=bondindices[::-1]
        # STEP 2
        temp=bondindicestotinkerclasses[bondindices]
        revtemp=temp[::-1]
        # STEP 3
        opbendbondindicestotinkerclasses[bondindices]=temp
        opbendbondindicestotinkerclasses[revbondindices]=revtemp
        opbendbondindicestosmartsatomorders[bondindices]=smartsatomorders 
        opbendbondindicestosmartsatomorders[revbondindices]=smartsatomorders 

    return opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders


def CheckTrigonalCenters(poltype,listofbondsforprm,mol):
    """
    Intent: Determine which of the possible opbend bond indices orientiation are trigonal centers ( need to know to put in correct format in the keyfile) 
    Input: List of bonds in input molecule, mol object
    Output: opbend bond indices -> boolean of which atoms in the bond are trigonal centers
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over all bonds in input molecule 
    2. Iterate over all atoms in bond
    3. If SP2 atom and also have 3 neighbors then it is a trigonal center (except if amine which may accidentally have hydrogens in plane and babel might think its SP2).
    4. Append array of booleans to dictionary 
    """
    opbendbondindicestotrigonalcenterbools={}
    # STEP 1
    for bondindices in listofbondsforprm:
        bondindices=tuple(bondindices)
        boolarray=[]
        babelindices=[i+1 for i in bondindices]
        atoms=[mol.GetAtom(i) for i in babelindices]    
        # STEP 2 
        for a in atoms:
            # STEP 3
            hyb=a.GetHyb()
            numhyds=0
            neighbs=list(openbabel.OBAtomAtomIter(a))
            for natom in neighbs:
                natomicnum=natom.GetAtomicNum()
                if natomicnum==1:
                    numhyds+=1
            idx=a.GetIdx()
            atomicnum=a.GetAtomicNum()
            if len(neighbs)==3 and hyb==2:
                if atomicnum==7 and numhyds==2:
                    boolarray.append(False)

                else:
                    boolarray.append(True)
            else:
                boolarray.append(False)
        # STEP 4
        opbendbondindicestotrigonalcenterbools[bondindices]=boolarray
        revboolarray=boolarray[::-1]
        revbondindices=bondindices[::-1]
        opbendbondindicestotrigonalcenterbools[revbondindices]=revboolarray
    return opbendbondindicestotrigonalcenterbools


def CorrectPitorEnergy(poltype,torsionprms,torsiontopitor):
    """
    Intent: If input torsion parameters also have pitor from amoeba09 (term no longer used for aromatic torsion), then want to add contribution of pitor back to the second torsion term instead of including pitor in the key file
    Input: Array of torsion parameters lines, dictionary of torsion parameter line -> pitor line (if there was match)
    Output: Modified array of torsion parameter lines
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Determine how many times the middle bond in torsion occurs out of all pitor matches, will use this to determine how many torsion will partition the pitor parameter into and add evenly to each torsion around bond.
    2. Iterate over each line in torsion parameter array
    3. Grab pitorsion line that corresponds to match in dictionary (if that torsion has a match)
    4. Take pitor parameter and divide it by the number of torsions that occur around bond
    5. Modify the torsion parameter line second parameter (dominant in SP2 torsions) with addition of pitor 
     
    """
    newtorsionprms=[]
    torsiontocount={}
    middletocount={}
    # STEP 1 
    for torsion in torsiontopitor.keys():
        middle=tuple([torsion[1],torsion[2]])
        if middle not in middletocount.keys():
            middletocount[middle]=0
        middletocount[middle]+=1
    for torsion in torsiontopitor.keys():
        middle=tuple([torsion[1],torsion[2]])
        count=middletocount[middle]
        torsiontocount[torsion]=count
    # STEP 2
    for torline in torsionprms:
        torlinesplit=torline.split()
        tor=tuple([int(torlinesplit[1]),int(torlinesplit[2]),int(torlinesplit[3]),int(torlinesplit[4])])
        if tor in torsiontopitor:
            # STEP 3
            if tor in torsiontopitor.keys():
                pitorline=torsiontopitor[tor]
            else:
                pitorline=torsiontopitor[tor[::-1]]

            count=torsiontocount[tor]
            pitorlinesplit=pitorline.split()
            # STEP 4
            prm=float(pitorlinesplit[3])/count
            torprm=float(torlinesplit[8])
            newtorprm=prm+torprm
            # STEP 5
            torlinesplit[8]=str(newtorprm)
        torline=' '.join(torlinesplit)+'\n'
        newtorsionprms.append(torline)
    return newtorsionprms


def FindMissingParameters(poltype,indicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs):
    """
    Intent: Determine if parameters are missing for bond, angle or vdw via if SMARTS matches neighbors or not
    Input: Dictionary of atomic indices -> SMARTS + atom order, rdkit mol object, openbabel mol object, dictionary of atomic index -> neighboring atomic indices 
    Output: Array of missing atomic indices for either, bond, angle or vdw
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of atomic indices -> SMARTS + atom order
    2. If vdw and valance ==1 (like hydrogen), then grab neighbors of neighbors, else just grab neighbors of indices (otherwise not enough information to check if hydrogen match is transferable just by neighboring atom)
    3. Match the SMARTS to input molecule and check if at a minimum the atomic indices exist in match
    4. If so, also check if all neighboring indices exist in the match.
    5. If not then add to array of missing indices (special case if user specifies to only determine certain atoms as missing vdw) (fragmenter uses this when want fragment job for only vdw dimer QM of certain atom)  
    """
    missing=[]
    # STEP 1
    for indices,smartsatomorders in indicestosmartsatomorders.items():
        # STEP 2
        nindexes=[]
        for idx in indices: 
            neighborindexes=indextoneighbidxs[idx]
            for nidx in neighborindexes:
                if nidx not in nindexes:
                    nindexes.append(nidx)
                if len(indices)==1 and len(neighborindexes)==1:
                    nneighborindexes=indextoneighbidxs[nidx]
                    for nnidx in nneighborindexes:
                        if nnidx not in nindexes:
                            nindexes.append(nnidx) 
        smarts=smartsatomorders[0]
        substructure = Chem.MolFromSmarts(smarts)
        # STEP 3
        matches=rdkitmol.GetSubstructMatches(substructure)
        for match in matches:
            allin=True
            for idx in indices:
                if idx not in match:
                    allin=False
            if allin==True:
                matcharray=match
                break

        # STEP 4 
        check=CheckIfNeighborsExistInSMARTMatch(poltype,nindexes,matcharray)
        # STEP 5
        if check==False or '*' in smarts or '~' in smarts:
            if len(indices)==1: # vdw
                index=indices[0]
                if poltype.onlyvdwatomindex!=None:
                    idx=indices[0]+1
                    if idx==poltype.onlyvdwatomindex:
                        missing.append(indices)
                else:
                    missing.append(indices)
            else:
                missing.append(indices)
        else:
            if len(indices)==1: # vdw
                index=indices[0]
                if poltype.onlyvdwatomindex==index:
                    missing.append(indices)


    return missing

def ReduceMissingVdwByTypes(poltype,vdwmissing):
    """
    Intent: If vdw missing for atoms of same type, just choose one to keep as missing (will use later in fragmenter if user specifies to generate dimer QM)
    Input: Array of missing vdw atom indices 
    Output: Modified array of missing vdw atom indices 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over array of missing vdw atom indices
    2. For each atom index, determine type number and append to an array
    3. If the same type has already occured in array, then skip atom index, else append atom index to new array of missing vdw atom indices 
    """
    reducedvdwmissing=[]
    typesfound=[]
    # STEP 1
    for vdwarray in vdwmissing:
        vdwatomidx=vdwarray[0]
        vdwbabelidx=vdwatomidx+1
        # STEP 2
        vdwtype=poltype.idxtosymclass[vdwbabelidx]
        # STEP 3
        if vdwtype in typesfound:
            continue
        else:
            typesfound.append(vdwtype)
            reducedvdwmissing.append(vdwbabelidx)
    return reducedvdwmissing



def ConvertToBabelList(poltype,listforprm):
    """
    Intent: just convert rdkit indices -> babel indices for bond and angle list
    Input: Array of bond or angles in molecule
    Output: Array of bond or angles with shifted indices (add +1)
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each bond or angle in list of bonds/angles
    2. Shift indices by +1 and append to new array 
    """
    babellist=[]
    # STEP 1
    for ls in listforprm:
        # STEP 2
        babells=[i+1 for i in ls]
        babellist.append(babells)

    return babellist



def AddReverseKeys(poltype,tinkerclassestopoltypeclasses):
    """
    Intent: For opbend tinker classes -> poltype classes , add reverse keys until determine which orientation is correct (opbend needs second class to be trigonal center) 
    Input: opbend tinker classes -> poltype classes
    Output: New dictionary of tinker classes -> poltype classes 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over dictionary of opbend tinker classes -> poltype classes
    2. Take the reverse of the tinker classes key and append both to new dictionary 
    """
    newtinkerclassestopoltypeclasses={}
    # STEP 1
    for tinkerclasses,poltypeclasses in tinkerclassestopoltypeclasses.items():
        # STEP 2
        revtinkerclasses=tinkerclasses[::-1]
        revpoltypeclasses=[]
        for ls in poltypeclasses:
            revls=ls[::-1]
            revpoltypeclasses.append(revls) 
        newtinkerclassestopoltypeclasses[tinkerclasses]=poltypeclasses
        newtinkerclassestopoltypeclasses[revtinkerclasses]=revpoltypeclasses
    return newtinkerclassestopoltypeclasses


def TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools):
    """
    Intent: Convert dictionary of opbend bond indices -> trigonal center to dictionary of opbend tinker classes -> trigional centers  
    Input: Dictionary of opbend bond indices -> trigonal center
    Output: Dictionary of opbend tinker classes -> trigional centers
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of opbend bond indices -> tinker classes
    2. Determine the array of trigonal center booleans from input dictionary and bond indices
    3. Save boolean array in new dictionary with tinker classes as the key 
    """
    opbendtinkerclassestotrigonalcenterbools={}
    # STEP 1
    for bondindices,tinkerclasses in opbendbondindicestotinkerclasses.items():
        # STEP 2
        boolarray=opbendbondindicestotrigonalcenterbools[bondindices]
        # STEP 3
        append=True
        if tuple(tinkerclasses) in opbendtinkerclassestotrigonalcenterbools.keys(): # sometimes different parts of molecule match to same amoeba09 opbend but one may be trigonal center and another may not be (if ring can pucker etc), so when mapping tinker classes -> list of poltype classes, need to enfore keeping any dictionary entry that says tinkerclass is not trigonal center. For the trigonal center that also matched to same amoeba09 tinker classes, those parameters will be "determined missing" then defaults added for that one later.
            prevboolarray=opbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]
            prevfalsecount=prevboolarray.count(False)
            falsecount=boolarray.count(False)
            if prevfalsecount>falsecount:
                append=False
        if append==True:
            opbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]=boolarray
            temp=opbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]
    newopbendtinkerclassestotrigonalcenterbools={}
    for tinkerclasses in opbendtinkerclassestotrigonalcenterbools.keys():
        boolarray=opbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]
        revtinkerclasses=tinkerclasses[::-1]
        revboolarray=boolarray[::-1]
        newopbendtinkerclassestotrigonalcenterbools[tuple(revtinkerclasses)]=revboolarray
        newopbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]=boolarray

    return newopbendtinkerclassestotrigonalcenterbools


def FilterDictionaries(poltype,dics,ls):
    """
    Intent: Remove keys in input list of dictionaries that are not in input list (want only planar bonds in dictionary for opbend)
    Input: List of dictionaries for opbend, list of planar bond indices 
    Output: Modified array of dictionaries
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over array of dictionaries  
    2. For each key in dictionary if key is not in input list of keys then dont add to new dictionary
    3. Append new dictionary to new array of dictionaries 
    """
    newdics=[]
    # STEP 1
    for dic in dics:
        newdic={}
        # STEP 2
        for key,value in dic.items():
            if key in ls:
                newdic[key]=value
        # STEP 3
        newdics.append(newdic)
    return newdics

def ConvertListOfListToListOfTuples(poltype,listoflist):
    """
    Intent: Sometimes need tuple for adding to dictionary keys, do this for planar bonds 
    Input: List of planar bonds ( a list) 
    Output: List of tuples
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over list of list
    2. Convert list to tuple and append to new array 
    """
    listoftuples=[]
    # STEP 1
    for item in listoflist:
        # STEP 2
        tup=tuple(item)
        listoftuples.append(tup)
    return listoftuples


def AddExternalDatabaseMatches(poltype, indicestosmartsatomorder,extindicestoextsmarts,smartsatomordertoparameters):
    """
    Intent: Add matches from external SMARTS database for vdw to dictionary of atom index > SMARTS + atom order (which originally only had amoeba09 SMARTS matches). 
    Input: dictionary of atom index > SMARTS + atom order , dictionary of atom indices -> SMARTS from external database match, dictionary of SMARTS + atom order -> parameters from external datbase match
    Output: modified dictionary of atom index > SMARTS + atom order
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of SMARTS + atom order -> parameters from external database match
    2. If find match in atom indices -> SMARTS from external database, then add atom indices and SMARTS to dictionary of atom index > SMARTS + atom order
    """
    newindicestosmartsatomorder=indicestosmartsatomorder.copy()
    # STEP 1
    for smartsatomorder in smartsatomordertoparameters.keys():
        smarts=smartsatomorder[0]
        for indices,extsmarts in extindicestoextsmarts.items():
            # STEP 2
            if smarts==extsmarts:
                newindicestosmartsatomorder[indices]=list(smartsatomorder)         
        
    return newindicestosmartsatomorder   
   
def AngleGuess(poltype,ita,itb,itc,iva,ivb,ivc):
    """
    Intent: If missing angle parameters, want to assign generic values based on element and valence (from tinker program valence)
    Input: Atomic numbers for three atoms, valencies for three atoms 
    Output: Angle force constant parameter guess
    Referenced By: AssignAngleGuessParameters
    Description:
    1. Copied special cases from valence program in tinker 
    """
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
    """
    Intent: If bond parameters are missing, then assign generic defaults based on atomic number and valencies
    Input: Atomic numbers and valencies
    Output: Bond force constant
    Referenced By: AssignBondGuessParameters
    Description:
    1. Copied special cases from valence program in tinker 
    """
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

def RepresentsInt(s):
    """
    Intent: Various functions when parsing parmaters and want to check if string is an integer or not
    Input: String
    Output: Boolean specifying if string is an integer or not
    Referenced By: Many functions
    Description: 
    1. Try to convert string to integer, if successful return True
    2. If fail, return False
    """
    # STEP 1
    try: 
        int(s)
        return True
    # STEP 2
    except ValueError:
        return False

def ReadDatabaseSmartsMap(poltype,databasepath): # smartsString classNumber className # comments
    """
    Intent: For AMOEBA+ parameters and latest bond/angle/strbnd/opbend AMOEBA parameters read in map of SMARTS -> atom class and atom class -> description and atom class -> comment
    Input: Database file
    Output: map of SMARTS -> atom class and atom class -> description and atom class -> comment
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over lines in database file
    2. Grab the SMARTS, tinker name/description and comment and put them into dictionaries 
    """
    smartstoatomclass={}
    atomclasstoclassname={}
    atomclasstocomment={}
    temp=open(databasepath,'r')
    results=temp.readlines()
    temp.close()
    # STEP 1
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            # STEP 2
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
 
def ReadDatabaseSmartsMapPolarize(poltype,databasepath): # smartsString className # comments
    """
    Intent: Read polarize parameter database
    Input: Database file
    Output: SMARTS -> atom class , atomclasstocomment , smartstoactualcomment
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines in database file
    2. Grab the SMARTS, tinker name/description and comment and put them into dictionaries 
    """
    smartstoatomclass={}
    atomclasstocomment={}
    smartstoactualcomment={}
    temp=open(databasepath,'r')
    results=temp.readlines()
    temp.close()
    # STEP 1
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            # STEP 2
            if '# ' in first or first=='#':
                continue
            else:
                smarts=linesplit[0]
                if RepresentsInt(linesplit[1])==True:
                    tinkerclass=int(linesplit[1])
                else:
                    tinkerclass=linesplit[1]
                comment=' '.join(linesplit[3:])
                comment=comment.replace('\n','').replace('#','').lstrip().rstrip()
                smartstoatomclass[smarts]=tinkerclass
                atomclasstocomment[tinkerclass]=comment
                smartstoactualcomment[smarts]=comment

    return smartstoatomclass, atomclasstocomment,smartstoactualcomment
 
def MatchAllSmartsToAtomIndices(poltype,smartstoatomclass): #rdkit 0 index based
    """
    Intent: Given input SMARTS from database, attempt to match to input molecule
    Input: Dictionary of SMARTS -> atom class
    Output: Dictionary of atom index -> SMARTS that matched, dictionary of atom index -> list of SMARTS matches (array of indices in input molecule corresponding to SMARTS that matched)
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over SMARTS in dictionary of SMARTS -> atom class
    2. Check if the SMARTS matches to the input molecule
    3. Find other atoms with same type and want to ensure matches for all atoms of same type are included (where as rdkit only return match for one atom out of group of same type).
    """
    atomindextoallsmarts={}
    atomindextoallsmartsmatches={}
    # STEP 1
    for smarts in smartstoatomclass.keys():
        substructure = Chem.MolFromSmarts(smarts)
        atomclass=smartstoatomclass[smarts]
        diditmatch=poltype.rdkitmol.HasSubstructMatch(substructure)
        # STEP 2
        if diditmatch==True:
            matches=list(poltype.rdkitmol.GetSubstructMatches(substructure))
            newmatches=[]
            # STEP 3
            for match in matches:
                newmatches.append(match)
                lastatom=match[-1]
                firstatom=match[0]
                firsttype=poltype.idxtosymclass[firstatom+1]
                secondtype=poltype.idxtosymclass[lastatom+1]
                if firsttype==secondtype: # special case where both ends of match are same type but the match only returns one way so need to add reverse match also
                    newmatches.append(match[::-1])
                otheratomindices=[firstatom]
                for idx,typenum in poltype.idxtosymclass.items():
                    rdkitidx=idx-1
                    if typenum==firsttype and rdkitidx!=firstatom:
                        otheratomindices.append(rdkitidx)

            for match in newmatches:
                for atomindex in otheratomindices:
                    if atomindex not in atomindextoallsmarts.keys():
                        atomindextoallsmarts[atomindex]=[]
                        atomindextoallsmartsmatches[atomindex]=[] 
                    if smarts not in atomindextoallsmarts[atomindex]:
                        atomindextoallsmarts[atomindex].append(smarts) 
                    if match not in atomindextoallsmartsmatches[atomindex]:
                        atomindextoallsmartsmatches[atomindex].append(match)   
               
    return atomindextoallsmarts,atomindextoallsmartsmatches


def MapIndicesToCommentsAtom(poltype,atomindextoallsmarts,smartstocomment,listofatomsforprm):
    """
    Intent: From SMARTS matches found, construct maps of comments -> all SMARTS with same comment and map of atom indices -> list of comments (that correpsond to SMARTS that were matched) 
    Input: Dictionary of atom index -> SMARTS list, dictionary of SMARTS -> comment, list of all atom indexes in molecule 
    Output: maps of comments -> all SMARTS with same comment and map of atom indices -> list of comments (that correpsond to SMARTS that were matched) 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over atoms in list of atom indices
    2. If atom index in dictionary of atom index -> SMARTS list
    3. Using map of SMARTS -> comments, add items to dictionaries 
    """
    atomcommentstolistofsmartslist={}
    atomindicestolistofatomcomments={}
    # STEP 1
    for atoms in listofatomsforprm:
        aindex=atoms[0]
        # STEP 2
        if aindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            combs = list(itertools.product(asmartslist)) 
            # STEP 3
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
    """
    Intent: Using dictionaries that map to only each atom for SMARTS etc, construct dictionaries for bond and angle matches   
    Input: Dictionary of atom index -> list of SMARTS, dicionary of SMARTS -> comment, list of all bonds in molecule, list of all angles in molecule 
    Output: Map of array of comments for each atom in bond -> list of SMARTS , map of bond indices -> list of comments , Map of array of comments for each atom in bond -> list of SMARTS , map of bond indices -> list of comments
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over list bonds in molecule
    2. Grab the SMARTS list for each atom in the bond
    3. Generate combination of all SMARTS for each list
    4. For each combination of SMARTS, grab the corresponding comments and add information into dictionaries
    5. Iterate over list angles in molecule
    6. Grab the SMARTS list for each atom in the angle
    7. Generate combination of all SMARTS for each list
    8. For each combination of SMARTS, grab the corresponding comments and add information into dictionaries
    """
    bondcommentstolistofsmartslist={}
    bondindicestolistofbondcomments={}
    # STEP 1
    for bond in listofbondsforprm:
        aindex=bond[0] 
        bindex=bond[1] 
        # STEP 2
        asmartslist=atomindextoallsmarts[aindex]
        bsmartslist=atomindextoallsmarts[bindex]
        # STEP 3
        combs = list(itertools.product(asmartslist,bsmartslist)) 
        # STEP 4
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
    # STEP 5
    for angle in listofanglesforprm:
        aindex=angle[0] 
        bindex=angle[1] 
        cindex=angle[2] 
        # STEP 6
        asmartslist=atomindextoallsmarts[aindex]
        bsmartslist=atomindextoallsmarts[bindex]
        csmartslist=atomindextoallsmarts[cindex]
        # STEP 7
        combs = list(itertools.product(asmartslist,bsmartslist,csmartslist)) 
        # STEP 8
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
    """
    Intent: Convert dictionaries that have atomic indices as keys to having class numbers as keys instead (parameters map from class numbers)  
    Input: Dictionary of atom index -> SMARTS list, dictionary of SMARTS -> atom class , list of bonds in molecule, list of angles in molecule, list of planar bonds in molecule 
    Output: Dictionaries of classes -> SMARTS list for bond/angle/strbnd/opbend and dictionaries of indices -> list of all classes (for each SMARTS that matched) for bond/angle/strbnd/opbend  
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over all bonds in molecule
    2. If bond exists in dictionary of atomindex -> SMARTS list
    3. Generate all combinations of SMARTS for each atom in bond
    4. For each combination find the class number and then store in dictionary 
    5. Iterate over all angles in molecule
    6. If angle exists in dictionary of atomindex -> SMARTS list
    7. Generate all combinations of SMARTS for each atom in angle
    8. For each combination find the class number and then store in dictionary 
    """
    bondclassestolistofsmartslist={}
    angleclassestolistofsmartslist={}
    strbndclassestolistofsmartslist={}
    opbendclassestolistofsmartslist={}
    bondindicestolistofbondclasses={}
    angleindicestolistofangleclasses={}
    strbndindicestolistofstrbndclasses={}
    opbendindicestolistofopbendclasses={}
    # STEP 1
    for bond in listofbondsforprm:
        aindex=bond[0]
        bindex=bond[1] 
        # STEP 2 
        if aindex in atomindextoallsmarts.keys() and bindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            bsmartslist=atomindextoallsmarts[bindex]
            # STEP 3
            combs = list(itertools.product(asmartslist, bsmartslist)) 
            # STEP 4
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

    # STEP 5
    for angle in listofanglesforprm:
        aindex=angle[0]
        bindex=angle[1] 
        cindex=angle[2] 
        # STEP 6
        if aindex in atomindextoallsmarts.keys() and bindex in atomindextoallsmarts.keys() and cindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            bsmartslist=atomindextoallsmarts[bindex]   
            csmartslist=atomindextoallsmarts[cindex]   
            # STEP 7
            combs = list(itertools.product(asmartslist,bsmartslist,csmartslist)) 
            # STEP 8
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
    """
    Intent: The naming is odd here, comments in this function are really like a pseudo SMARTS string that tells what the type is. These "pseudo smarts" (just a name) are associated with the parameters in the parameter file. 
    Input: Dictionary of atom comments -> list of SMARTS , dictionary of atom indices -> list of comments  
    Output: Modified dictionary of atom indices -> atom comments , dictionary of comment -> parameters  
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Generate dictionary of atom indices -> list of comments by reading database file
    2. Modify dictionary of atom indices -> atom comments to only include entries that have parameters associated with comments
    """
    # STEP 1
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
    # STEP 2
    atomindicestolistofatomcomments=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist,atomcommentstoparameters,prin=True)
    return atomindicestolistofatomcomments,atomcommentstoparameters

def SearchForParameters(poltype,bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses):
    """
    Intent: Generate dictionaries of classes -> parameters for bond/angle/strbnd/opbend
    Input: Dictionaries of classes to list of SMARTS for bond/angle/strbnd/opbend and dictionaries of indices -> classes for bond/angle/strbnd/opbend 
    Output: Modified dictionaries of indices -> classes,  dictionaries of classes -> parameters for bond/angle/strbnd/opbend
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Generate dictionaries of classes -> parameters for bond/angle/strbnd/opbend by reading in parameterfile
    2. Modify dictionaries of atom indices -> atom comments to only include entries that have parameters associated with comments
    """
    bondclassestoparameters={}
    angleclassestoparameters={}
    strbndclassestoparameters={}
    opbendclassestoparameters={}
    temp=open(poltype.latestsmallmoleculeprmlib,'r')
    results=temp.readlines()
    temp.close()
    # STEP 1
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
    # STEP 2
    bondindicestolistofbondclasses=RemoveIndicesThatDontHaveParameters(poltype,bondindicestolistofbondclasses,bondclassestolistofsmartslist,bondclassestoparameters)
    angleindicestolistofangleclasses=RemoveIndicesThatDontHaveParameters(poltype,angleindicestolistofangleclasses,angleclassestolistofsmartslist,angleclassestoparameters)
    strbndindicestolistofstrbndclasses=RemoveIndicesThatDontHaveParameters(poltype,strbndindicestolistofstrbndclasses,strbndclassestolistofsmartslist,strbndclassestoparameters)
    opbendindicestolistofopbendclasses=RemoveIndicesThatDontHaveParameters(poltype,opbendindicestolistofopbendclasses,opbendclassestolistofsmartslist,opbendclassestoparameters)

    return bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses,bondclassestoparameters,angleclassestoparameters,strbndclassestoparameters,opbendclassestoparameters

def RemoveIndicesThatDontHaveParameters(poltype,indicestolistofclasses,classestolistofsmartslist,classestoparameters,prin=False):
    """
    Intent: Just check to ensure that assigned class has parameters in parameter file (like a sanity check)
    Input: Dictionary of indices -> list of classes, dictionary of classes -> list of SMARTS, dictionary of classes -> parmaters
    Output: Modified indices -> list of classes
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of indices -> list of classes
    2. For each classes in list of classes, check if there exists parameters in dictionary of classes -> parameters
    3. If parameters dont exist, then want to remove classes from final dictionary
    """
    indicestodelete=[]
    # STEP 1
    for i in range(len(indicestolistofclasses.keys())):
        indices=list(indicestolistofclasses.keys())[i]
        listofclasses=indicestolistofclasses[indices]
        classesmissing=[]
        classesexisting=[]
        # STEP 2
        for classes in listofclasses:
            # STEP 3
            if classes not in classestoparameters.keys() and classes[::-1] not in classestoparameters.keys():
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
    """
    Intent: When have several possible SMARTS matches to choose from, want to choose the one with the largest SMARTS string which corresponds to more specificity about envioroment and is more transferable. 
    Input: Dictionary of indices -> list of classes, dictionary of classes -> list of SMARTS
    Output: Dictionary of indices -> classes, dictionary of indices -> SMARTS
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of indices -> list of classes
    2. For each classes in list of classes, find list of SMARTS that all map to those classes
    3. Iterate over SMARTS in list of SMARTS
    4. Find all atom indices that are matched from list of SMARTS
    5. Create a fragment molecule from input molecule with only those indices that were matched
    6. Save the molecule and SMARTS list in array
    7. For each SMARTS list and fragment generated, determine the best SMARTS match via tanimoto similarity between input molecule and generated fragment molecule
    8. Save best SMARTS list (like for bond there are two SMARTS per atom) and assicated classes in dicitonaries 
    """
    indicestoclasses={}
    indicestosmartslist={}
    # STEP 1
    for indices,listofclasses in indicestolistofclasses.items():
        fragmentsarray=[]
        classesarray=[]
        smartslistarray=[] 
        # STEP 2
        for classes in listofclasses:
            listofsmartslist=classestolistofsmartslist[classes]
            # STEP 3
            for smartslist in listofsmartslist:
                # STEP 4
                matchedindices=MatchAllSMARTS(poltype,smartslist,indices) 
                # STEP 5
                fragmentmol=CreateFragment(poltype,matchedindices)
                # STEP 6
                fragmentsarray.append(fragmentmol)
                classesarray.append(classes)
                smartslistarray.append(smartslist)
        # STEP 7
        classes,smartslist=DetermineBestSMARTSMatch(poltype,fragmentsarray,classesarray,smartslistarray)        
        # STEP 8
        indicestoclasses[indices]=classes
        indicestosmartslist[indices]=smartslist

    return indicestoclasses,indicestosmartslist

def DetermineBestSMARTSMatch(poltype,listoffragments,listofclasses,listofsmartslist):  
    """
    Intent: Want to find array of SMARTS (one for each atom in bond or angle etc) that has closest match to input molecule. 
    Input: List of fragments, list of classes, list of SMARTS lists
    Output: Best classes and best smarts list
    Referenced By: FindBestSMARTSMatch
    Description:
    1. Iterate over list of fragment molecules
    2. Grab list of classes and SMARTS that correspond
    3. Compute fingerprints via rdkit for input molecule and fragment molecule
    4. Compute tanimoto similairty between input molecule and fragment then save in array for finding closest match
    5. Grab max tanimoto value
    6. If two matches have same tanimoto simliairty value, then choose the SMARTS that is largest (count length of all SMARTS in list)  
    7. Otherwise choose SMARTS list with highest tanimoto value 
    """
    tanimotoarray=[]    
    # STEP 1
    for i in range(len(listoffragments)):
        fragment=listoffragments[i]
        # STEP 2
        classes=listofclasses[i]
        smartslist=listofsmartslist[i]
        ms=[poltype.rdkitmol,fragment]
        # STEP 3
        fps = [Chem.RDKFingerprint(x) for x in ms]
        # STEP 4
        tanimoto=DataStructs.FingerprintSimilarity(fps[0],fps[1], metric=DataStructs.DiceSimilarity) 
        tanimotoarray.append(tanimoto)
    # STEP 5
    maxvalue=max(tanimotoarray)
    nums=tanimotoarray.count(maxvalue)
    # STEP 6
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
        # STEP 7
        maxindex=tanimotoarray.index(maxvalue)
        classes=listofclasses[maxindex]
        smartslist=listofsmartslist[maxindex]
    return classes,smartslist 


def ChooseMostDescriptiveSMARTSList(poltype,smartslisttocompare):
    """
    Intent: If two molecules have same tanimoto similarity to input molecule then count the length of all SMARTS in SMARTS list and choose the ones with largest length (more descriptive envioronment)
    Input: list of SMARTS list
    Output: SMARTS list that is most descriptive 
    Referenced By: DetermineBestSMARTSMatch
    Description:
    1. Iterate over SMARTS lists
    2. Compute the total length and save in array
    3. Find the max length and then grab the SMARTS list with max length
    """
    lengtharray=[]
    # STEP 1
    for smartslist in smartslisttocompare:
        # STEP 2
        lengths=[len(e) for e in smartslist]
        total=sum(lengths)
        lengtharray.append(total)
    # STEP 3
    maxlength=max(lengtharray)
    maxidx=lengtharray.index(maxlength)
    smartslist=smartslisttocompare[maxidx]
    return smartslist


def CreateFragment(poltype,matchedindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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


def GrabNewParameters(poltype,indicestoclasses,classestoparameters,keyword,indicestosmartslist,atomclasstocomment,opbendbondindicestotrigonalcenterbools=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestopoltypeclasses={}
    newpoltypeclassestocomments={}
    newpoltypeclassestosmartslist={}

    newprms=[]
    for indices,classes in indicestoclasses.items():
        smartslist=indicestosmartslist[indices]
        comments=[atomclasstocomment[k] for k in classes]
        babelindices=[k+1 for k in indices]
        if opbendbondindicestotrigonalcenterbools!=None:
            if indices in opbendbondindicestotrigonalcenterbools.keys():
                tribools=opbendbondindicestotrigonalcenterbools[indices]
            elif indices[::-1] in opbendbondindicestotrigonalcenterbools.keys():
                tribools=opbendbondindicestotrigonalcenterbools[indices[::-1]]
            if tribools[1]==False: # need second one to be trigonal center:
                babelindices=babelindices[::-1]
                indices=indices[::-1]

 
        poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
        indicestopoltypeclasses[indices]=poltypeclasses
        newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
        newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
        poltypeclasses=[str(k) for k in poltypeclasses]
        if keyword=='opbend':
            poltypeclasses.append('0')
            poltypeclasses.append('0')
        if classes in classestoparameters.keys():
            parameters=classestoparameters[classes]
        elif classes[::-1] in classestoparameters.keys():
            parameters=classestoparameters[classes[::-1]]

        parameters=[str(k) for k in parameters]
        poltypeclassesstring=' '.join(poltypeclasses)
        parametersstring=' '.join(parameters)
        newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
        newprms.append(newprm)
        if opbendbondindicestotrigonalcenterbools!=None:
            if tribools[0]==True and tribools[1]==True:
                babelindices=babelindices[::-1]
                indices=indices[::-1]
                poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
                indicestopoltypeclasses[indices]=poltypeclasses
                newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
                newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
                poltypeclasses=[str(k) for k in poltypeclasses]
                if keyword=='opbend':
                    poltypeclasses.append('0')
                    poltypeclasses.append('0')
                parameters=[str(k) for k in parameters]
                poltypeclassesstring=' '.join(poltypeclasses)
                parametersstring=' '.join(parameters)
                newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
                newprms.append(newprm)

        

    return indicestopoltypeclasses,newprms,newpoltypeclassestocomments,newpoltypeclassestosmartslist

def GrabNewParametersPolarize(poltype,indicestoclasses,classestoparameters,keyword,indicestosmartslist,smartstoactualcomment):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestopoltypeclasses={}
    newpoltypeclassestocomments={}
    newpoltypeclassestosmartslist={}

    newprms=[]
    for indices,classes in indicestoclasses.items():
        smartslist=indicestosmartslist[indices]
        comments=[smartstoactualcomment[k] for k in smartslist]
        babelindices=[k+1 for k in indices]
        poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
        indicestopoltypeclasses[indices]=poltypeclasses
        newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
        newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
        poltypeclasses=[str(k) for k in poltypeclasses]
        if keyword=='opbend':
            poltypeclasses.append('0')
            poltypeclasses.append('0')
        if classes in classestoparameters.keys():
            parameters=classestoparameters[classes]
        elif classes[::-1] in classestoparameters.keys():
            parameters=classestoparameters[classes[::-1]]

        parameters=[str(k) for k in parameters]
        poltypeclassesstring=' '.join(poltypeclasses)
        parametersstring=' '.join(parameters)
        newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
        newprms.append(newprm)
        

    return indicestopoltypeclasses,newprms,newpoltypeclassestocomments,newpoltypeclassestosmartslist




def RemoveOldParametersKeepNewParameters(poltype,prms,newindicestopoltypeclasses,keyword,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,removefromdic=True):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newpoltypeclassestoparametersmartsatomorders=copy.deepcopy(poltypeclassestoparametersmartsatomorders)
    newpoltypeclassestosmartsatomorders=copy.deepcopy(poltypeclassestosmartsatomorders)
    newpoltypeclassestoelementtinkerdescrips=copy.deepcopy(poltypeclassestoelementtinkerdescrips) 
    newprms=[]
    dellist=[]
    appendlist=[]
    for poltypeclassesls in newpoltypeclassestoparametersmartsatomorders.keys():
        newpoltypeclasseslist=[]
        for poltypeclasses in poltypeclassesls:
            poltypeclasses=list(poltypeclasses)
            if poltypeclasses in newindicestopoltypeclasses.values() or poltypeclasses[::-1] in newindicestopoltypeclasses.values():
                pass
            else:
                newpoltypeclasseslist.append(tuple(poltypeclasses))
        newpoltypeclasseslist=tuple(newpoltypeclasseslist)
        if poltypeclassesls not in dellist:
            dellist.append(poltypeclassesls)
            if len(newpoltypeclasseslist)!=0:
                appendlist.append(newpoltypeclasseslist)


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
        else:
            newprms.append(prm)
    for i in range(len(dellist)):
        item=dellist[i]
        parametersmartatomordervalue=newpoltypeclassestoparametersmartsatomorders[item]
        smartatomordervalue=newpoltypeclassestosmartsatomorders[item]
        descrips=newpoltypeclassestoelementtinkerdescrips[item]
        if item in newpoltypeclassestoparametersmartsatomorders.keys():
            del newpoltypeclassestoparametersmartsatomorders[item]
        if item in newpoltypeclassestosmartsatomorders.keys():
            del newpoltypeclassestosmartsatomorders[item]
        if item in newpoltypeclassestoelementtinkerdescrips.keys():
            del newpoltypeclassestoelementtinkerdescrips[item]
    for newitem in appendlist:
        if len(newitem)!=0:
            newpoltypeclassestoparametersmartsatomorders[newitem]=parametersmartatomordervalue
            newpoltypeclassestosmartsatomorders[newitem]=smartatomordervalue
            newpoltypeclassestoelementtinkerdescrips[newitem]=descrips
     
    return newprms,newpoltypeclassestoparametersmartsatomorders,newpoltypeclassestosmartsatomorders,newpoltypeclassestoelementtinkerdescrips

def MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstocommentpolar):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    smartstocomment={}
    for smarts,atomclass in smartstoatomclasspolar.items():
        comment=atomclasstocommentpolar[atomclass]
        smartstocomment[smarts]=comment

    return smartstocomment

def GrabSmartsToSoluteRadiiMap(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newkeytovalue={}
    for key,value in keytovalue.items():
        newvalues=[]
        for item in value:
            newitem=[item]
            newvalues.append(newitem)
        newkeytovalue[key]=newvalues    

    return newkeytovalue

def GenerateSoluteParameters(poltype,atomindicestosmartslist,smartstosoluteradiiprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    soluteprms=[]
    for atomindices,smartslist in atomindicestosmartslist.items():
        atomindex=atomindices[0]
        smarts=smartslist[0]
        babelatomindex=atomindex+1
        poltypeclass=poltype.idxtosymclass[babelatomindex]
        prms=smartstosoluteradiiprms[smarts]
        prmsline=' '.join(prms)
        commentline='#SOLUTE-SMARTS'+' '+str(poltypeclass)+' '+smarts+'\n'
        line='SOLUTE'+' '+str(poltypeclass)+' '+prmsline+'\n'
        if line not in soluteprms:
            soluteprms.append(commentline)
            soluteprms.append(line)

    return soluteprms

def RemoveIndicesMatchedFromNewDatabase(poltype,indicestosmartsatomorders,indicestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    removelist=[]
    for indices,smartsatomorders in indicestosmartsatomorders.items():
        if indices in indicestopoltypeclasses.keys() or indices[::-1] in indicestopoltypeclasses.keys():
            removelist.append(indices)
    for indices in removelist:
        del indicestosmartsatomorders[indices]


    return indicestosmartsatomorders


def TrimDictionary(poltype,dictotrim,dicref):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keylist=[]
    for key,value in dictotrim.items():
        if key not in dicref.keys() and key[::-1] not in dicref.keys():
            keylist.append(key)
    for key in keylist:
        del dictotrim[key]
    return dictotrim

def AddDictionaryItems(poltype,dictoaddto,dicref):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for key,value in dicref.items():
        dictoaddto[key]=value
    return dictoaddto


def RemovePoltypeClassesFromNewMatches(poltype,missingtinkerclassestopoltypeclasses,poltypeclassestoparametersmartsatomorders):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    newmissingtinkerclassestopoltypeclasses={}
    for missingtinkerclasses,poltypeclasses in missingtinkerclassestopoltypeclasses.items():
        for polcls in poltypeclasses:
            tuppolcls=tuple([tuple(polcls)])
            if tuppolcls in poltypeclassestoparametersmartsatomorders.keys() or tuppolcls[::-1] in poltypeclassestoparametersmartsatomorders.keys():
                if missingtinkerclasses not in newmissingtinkerclassestopoltypeclasses.keys():
                    newmissingtinkerclassestopoltypeclasses[missingtinkerclasses]=[]
                newmissingtinkerclassestopoltypeclasses[missingtinkerclasses].append(polcls)
  

    return newmissingtinkerclassestopoltypeclasses

def ExtractMissingPoltypeClasses(poltype,missingtinkerclassestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    missingpoltypeclasses=[]
    for tinkerclasses,poltypeclassesls in missingtinkerclassestopoltypeclasses.items():
        for cls in poltypeclassesls:
            missingpoltypeclasses.append(list(cls))
    return missingpoltypeclasses 


def MergeDictionaries(poltype,add,addto):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for key,value in add.items():
        addto[key]=value

    return addto


def GetPolarIndexToPolarizePrm(poltype,polarprmstotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    polarindextopolarizeprm={}
    for line in polarprmstotransferinfo.keys():
        linesplit=line.split()
        typenum=int(linesplit[1])
        prm=float(linesplit[2])
        for index in poltype.idxtosymclass:
            othertypenum=poltype.idxtosymclass[index]
            if typenum==othertypenum:
                polarindextopolarizeprm[index]=prm   
         

    return polarindextopolarizeprm


def TestBondAngleEquilValues(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpkeyfile='testbondangleequilvalues.key'
    tmpxyzfile='testbondangleequilvalues.xyz'
    alzout='testbondangleequilvaluesalz.out'
    shutil.copy(poltype.key4fname,tmpkeyfile)
    shutil.copy(poltype.xyzoutfile,tmpxyzfile)
    cmd = poltype.minimizeexe+' '+tmpxyzfile+' -k '+tmpkeyfile+' 0.1 > testbondangleequilvalues.out'
    poltype.call_subsystem([cmd], True)

    cmd = poltype.analyzeexe+' '+tmpxyzfile+'_2'+' -k '+tmpkeyfile+' '+' d > '+alzout
    poltype.call_subsystem([cmd], True)

    bondindicestonewbondequilvalues,angleindicestonewbondequilvalues=CheckBondAngleDeviationsFromQM(poltype,alzout)
   
    bondtypeindicestonewbondequilvalues=ConvertIndicesToTypeIndices(poltype,bondindicestonewbondequilvalues)
    angletypeindicestonewbondequilvalues=ConvertIndicesToTypeIndices(poltype,angleindicestonewbondequilvalues)
    ChangeBondAngleEquilValues(poltype,bondtypeindicestonewbondequilvalues,angletypeindicestonewbondequilvalues)

def ChangeBondAngleEquilValues(poltype,bondtypeindicestonewbondequilvalues,angletypeindicestonewbondequilvalues): 
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.key4fname,'r')
    results=temp.readlines()
    temp.close()
    tempname=poltype.key4fname.replace('.key','_TEMP.key')
    temp=open(tempname,'w')
    for line in results:
        if '#' not in line:
            found=False
            linesplit=line.split()
            if 'angle' in line and poltype.writeoutangle==True:
                typeindices=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
                if typeindices in angletypeindicestonewbondequilvalues.keys():
                    equilvalue=angletypeindicestonewbondequilvalues[typeindices]
                    found=True
                elif typeindices[::-1] in angletypeindicestonewbondequilvalues.keys():
                    equilvalue=angletypeindicestonewbondequilvalues[typeindices[::-1]]
                    found=True
                if found==True:
                    linesplit=re.split(r'(\s+)', line)
                    linesplit=linesplit[:11]
                    linesplit.append('\n')
                    linesplit[10]=str(equilvalue)
                    line=''.join(linesplit)
                
            if 'bond' in line and poltype.writeoutbond==True:
                typeindices=tuple([int(linesplit[1]),int(linesplit[2])])
                if typeindices in bondtypeindicestonewbondequilvalues.keys():
                    equilvalue=bondtypeindicestonewbondequilvalues[typeindices]
                    found=True
                elif typeindices[::-1] in bondtypeindicestonewbondequilvalues.keys():
                    equilvalue=bondtypeindicestonewbondequilvalues[typeindices[::-1]]
                    found=True
                if found==True:
                    linesplit=re.split(r'(\s+)', line)
                    linesplit=linesplit[:11]
                    linesplit.append('\n')
                    linesplit[8]=str(equilvalue)
                    line=''.join(linesplit)
        temp.write(line)
    temp.close()
    os.remove(poltype.key4fname)
    os.rename(tempname,poltype.key4fname)
             


def ConvertIndicesToTypeIndices(poltype,indicestonewequilvalues):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    typeindicestonewequilvalues={}
    for indices,newequilvalues in indicestonewequilvalues.items():
        typeindices=[]
        for index in indices:
            symclass=poltype.idxtosymclass[index]
            typeindices.append(symclass)
        typeindicestonewequilvalues[tuple(typeindices)]=newequilvalues


    return typeindicestonewequilvalues


def CheckBondAngleDeviationsFromQM(poltype,alzout):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bondindicestonewbondequilvalues={}
    angleindicestonewbondequilvalues={}
    temp=open(alzout,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Additional' in line or 'Classes' in line or 'Individual' in line or 'Atom' in line or 'Stretch' in line or 'Bending' in line or 'Torsional' in line:
            continue
        linesplit=line.split()
        tol=None
        keyword=None
        if 'Angle' in line:
            indexlist=[1,2,3]
            equilindices=[4,5]
            indices,qmequil,currentequil=GrabIndicesAndEquilValues(linesplit,indexlist,equilindices)
            tol=1
            keyword='angle'
    
        if 'Bond' in line:
            indexlist=[1,2]
            equilindices=[3,4]
            indices,qmequil,currentequil=GrabIndicesAndEquilValues(linesplit,indexlist,equilindices)
            tol=2
            keyword='bond'
            
        if tol!=None:
            deviation=(np.abs(qmequil-currentequil)*100)/qmequil
            if deviation>=tol:
                shift=qmequil-(currentequil-qmequil)
                if keyword=='bond':
                    bondindicestonewbondequilvalues[tuple(indices)]=shift
                elif keyword=='angle':
                    angleindicestonewbondequilvalues[tuple(indices)]=shift



    return bondindicestonewbondequilvalues,angleindicestonewbondequilvalues

def GrabIndicesAndEquilValues(linesplit,indexlist,equilindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indices=[]
    for index in indexlist:
        indexele=linesplit[index] 
        indexelesplit=indexele.split('-')
        theindex=int(indexelesplit[0])
        indices.append(theindex)
    qmequil=float(linesplit[equilindices[0]])
    currentequil=float(linesplit[equilindices[1]])

    return indices,qmequil,currentequil

def RemoveIndicesMatchedFromExternalDatabase(poltype,indicestosmartsatomorder,indicestoextsmarts):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    missingindicestosmartsatomorders={}
    for indices,smartsatomorder in indicestosmartsatomorder.items():
        if indices in indicestoextsmarts.keys() or indices[::-1] in indicestoextsmarts.keys():
            pass
        else:
            missingindicestosmartsatomorders[indices]=smartsatomorder
    return missingindicestosmartsatomorders


def ExtractTransferInfo(poltype,polarprmstotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    polartypetotransferinfo={}
    for polarprms,transferinfo in polarprmstotransferinfo.items():
        linesplit=polarprms.split()
        polartype=linesplit[1]
        polartypetotransferinfo[polartype]=transferinfo

    return polartypetotransferinfo


def StiffenZThenBisectorAngleConstants(poltype,keyfilename):
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
    tempname=keyfilename.replace('.key','TEMP.key')
    temp=open(tempname,'w')
    angles=[]
    for line in results:
        if 'multipole ' in line and '#' not in line:
            linesplit=line.split()
            if len(linesplit)==6:
                newlinesplit=linesplit[1:-2]
                newlinesplit=[int(i) for i in newlinesplit]
                newlinesplit=[int(np.abs(i)) for i in newlinesplit]
                a=newlinesplit[0]
                b=newlinesplit[1]
                c=newlinesplit[2]
                angle=tuple([b,a,c])
                angles.append(angle)
    for line in results:
        if 'angle ' in line and '#' not in line:
            linesplit=line.split()
            indices=[linesplit[1],linesplit[2],linesplit[3]]
            indices=[int(i) for i in indices]
            indices=tuple(indices)
            for angle in angles:
                if angle[1]==indices[1]:
                    value=float(linesplit[4])
                    value=str(2*value)
                    linesplit[4]=value
                    string='# Doubling angle force constant because z-then-bisector frame '+'\n'
                    line=string+' '.join(linesplit)
        temp.write(line)
                    
                    


    temp.close()
    os.remove(keyfilename)
    os.rename(tempname,keyfilename)


def GenerateRdkitMolObjectsParameterSMARTS(poltype,parametersmartslist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    parametersmartstordkitmol={}
    temptotalcharge=poltype.totalcharge
    for parametersmarts in parametersmartslist:
        prmmol=Chem.MolFromSmarts(parametersmarts)
        prmmol.UpdatePropertyCache()
        poltype.totalcharge=None
        prmmol,atomindextoformalcharge=poltype.CheckInputCharge(prmmol)
        Chem.SanitizeMol(prmmol)
        poltype.totalcharge=temptotalcharge
        parametersmartstordkitmol[parametersmarts]=prmmol


    return parametersmartstordkitmol


def GrabVdwParametersFromParent(poltype,oldvdwprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    currentdir=os.getcwd()
    parenttypestofragtypes=json.load(open("parentsymclasstofragsymclasses.txt"))
    vdwkey='parentvdw.key'
    parenttypes=list(parenttypestofragtypes.keys())
    parenttypes=[str(i) for i in parenttypes]
    fragtypes=[list(parenttypestofragtypes.values())]
    new=[]
    for fragtype in fragtypes:
        for actualtype in fragtype:
            new.append(str(actualtype[0]))
    parentvdwtransferinfo={}
    fragtypes=new[:]
    vdwprms=GrabVdwParametersFromKeyFile(poltype,vdwkey,parenttypes)
    os.chdir(currentdir)
    newvdwprms=[]
    for oldvdwprmline in oldvdwprms:
        linesplit=oldvdwprmline.split()
        typenum=linesplit[1]
        if typenum not in fragtypes:
            newvdwprms.append(oldvdwprmline)
    for vdwprmline in vdwprms:
        linesplit=vdwprmline.split()
        typenum=linesplit[1]
        fragtype=str(parenttypestofragtypes[typenum][0])
        linesplit[1]=fragtype
        otherline='# Transferring from parent vdw type '+typenum+'\n'
        parentvdwtransferinfo[fragtype]=otherline
        newline=' '.join(linesplit)+'\n'
        newvdwprms.append(newline)

    return newvdwprms,parentvdwtransferinfo


def GrabVdwParametersFromKeyFile(poltype,vdwkey,parenttypes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    vdwprms=[]
    temp=open(vdwkey,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'vdw ' in line:
            typenum=linesplit[1]
            for parenttype in parenttypes:
                if str(parenttype)==typenum:
                    vdwprms.append(line)

    return vdwprms



def AddParentVdwTransferInfo(poltype,parentvdwtransferinfo,vdwprmstotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for line in vdwprmstotransferinfo.keys():
        linesplit=line.split()
        typenum=linesplit[1]
        if str(typenum) in parentvdwtransferinfo.keys():
            otherline=parentvdwtransferinfo[typenum]
            vdwprmstotransferinfo[line]+=otherline



    return vdwprmstotransferinfo

def AddTrivialOPBendForAmine(poltype,opbendprms,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomiter=openbabel.OBMolAtomIter(mol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        atomicnum=atom.GetAtomicNum()
        if atomicnum==7:
            hcount=0
            neighbs=[natom for natom in openbabel.OBAtomAtomIter(atom)]
            atomidxtoatomicnum={}
            for neighb in neighbs:
                natomidx=neighb.GetIdx()
                natomicnum=neighb.GetAtomicNum()
                atomidxtoatomicnum[natomidx]=natomicnum
                if natomicnum==1:
                    hcount+=1
            if hcount==2:
                for natomidx,atomicnum in atomidxtoatomicnum.items():
                    if atomicnum!=1:
                        theidx=natomidx
                    elif atomicnum==1:
                        oidx=natomidx
                indices=[theidx,atomidx]
                symclasses=[poltype.idxtosymclass[i] for i in indices]
                symclasses=[str(i) for i in symclasses]
                string=' '.join(symclasses)
                string='opbend '+string+' 0 0 0'+'\n'
                opbendprms.append(string)
                indices=[oidx,atomidx]
                symclasses=[poltype.idxtosymclass[i] for i in indices]
                symclasses=[str(i) for i in symclasses]
                string=' '.join(symclasses)
                string='opbend '+string+' 0 0 0'+'\n'
                opbendprms.append(string)



    return opbendprms

def GenerateTinkerClassesToPoltypeClasses(poltype,indices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tinkerclassestopoltypeclasses={}
    for indexls in indices:
        indexls=[i+1 for i in indexls]
        poltypeclasses=[poltype.idxtosymclass[i] for i in indexls]
        tup=tuple(poltypeclasses)
        ls=[tup]
        tinkerclassestopoltypeclasses[tup]=ls

    return tinkerclassestopoltypeclasses



def ForceSameResonanceTypesSameMatches(poltype,atomindicestocomments,atomindicestosmartslist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    for atomindices,comments in atomindicestocomments.items():
        smartslist=atomindicestosmartslist[atomindices]
        atomindex=atomindices[0]
        atomtype=poltype.idxtosymclass[atomindex+1]
        equivalentindices=[]
        for idx,typenum in poltype.idxtosymclass.items():
            if typenum==atomtype:
                equivalentindices.append(idx-1)
        equivalentindextosmartslist={}
        for idx in equivalentindices:
            oatomindices=tuple([idx])
            osmartslist=atomindicestosmartslist[oatomindices]
            equivalentindextosmartslist[oatomindices]=osmartslist
        lentoidx={}
        for eidx,smartslist in equivalentindextosmartslist.items():
            length=len(smartslist[0])
            lentoidx[length]=eidx
        maxlen=max(lentoidx.keys())
        maxidx=lentoidx[maxlen]
        maxsmartslist=atomindicestosmartslist[maxidx]
        maxcomments=atomindicestocomments[maxidx]
        atomindicestocomments[atomindices]=maxcomments
        atomindicestosmartslist[atomindices]=maxsmartslist


    return atomindicestocomments,atomindicestosmartslist


def GrabSmallMoleculeAMOEBAParameters(poltype,optmol,mol,rdkitmol,polarize=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    if polarize==True:
        listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForParameters(poltype,mol)

        smartstoatomclasspolar, atomclasstocommentpolar,smartstoactualcomment=ReadDatabaseSmartsMapPolarize(poltype,poltype.latestsmallmoleculesmartstotypespolarize) 
        atomclasstoatomclass=dict(zip(atomclasstocommentpolar.keys(),atomclasstocommentpolar.keys()))
       
        atomindextoallsmartspolar,atomindextoallsmartsmatchespolar=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasspolar)
        smartstocomment=MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstoatomclass)
        
        atomcommentstolistofsmartslist,atomindicestolistofatomcomments=MapIndicesToCommentsAtom(poltype,atomindextoallsmartspolar,smartstocomment,listofatomsforprm)
        atomindicestolistofatomcomments,atomcommentstoparameters=SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments)
        atomindicestocomments,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist)

        atomindicestocomments,atomindicestosmartslist=ForceSameResonanceTypesSameMatches(poltype,atomindicestocomments,atomindicestosmartslist)




        newpolarindicestopoltypeclasses,newpolarprms,newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist=GrabNewParametersPolarize(poltype,atomindicestocomments,atomcommentstoparameters,'polarize',atomindicestosmartslist,smartstoactualcomment) 
        
        polarprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newpolarprms,{},{},{},{},newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist,{},[],[],[],[])
        polarindextopolarizeprm=GetPolarIndexToPolarizePrm(poltype,polarprmstotransferinfo)
        polartypetotransferinfo=ExtractTransferInfo(poltype,polarprmstotransferinfo)
        return polarindextopolarizeprm,polartypetotransferinfo

    else:

    
        listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForParameters(poltype,mol)
        smartstosoluteradiiprms=GrabSmartsToSoluteRadiiMap(poltype)   
        atomindextoallsmartssolute,atomindextoallsmartsmatchessolute=MatchAllSmartsToAtomIndices(poltype,smartstosoluteradiiprms)
        atomindices=list(atomindextoallsmartssolute.keys())
        atomindices=[tuple([k]) for k in atomindices]
        atomindicestolistofatomindices=dict(zip(atomindices,atomindices))
        atomindextoallsmartssolute=MakeListOfListValues(poltype,atomindextoallsmartssolute)
        atomindicestoindices,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomindices,atomindextoallsmartssolute)
        soluteprms=GenerateSoluteParameters(poltype,atomindicestosmartslist,smartstosoluteradiiprms)
        smartstoatomclasspolar, atomclasstocommentpolar,smartstoactualcomment=ReadDatabaseSmartsMapPolarize(poltype,poltype.latestsmallmoleculesmartstotypespolarize) 
        atomclasstoatomclass=dict(zip(atomclasstocommentpolar.keys(),atomclasstocommentpolar.keys()))
        atomindextoallsmartspolar,atomindextoallsmartsmatchespolar=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasspolar)
        smartstocomment=MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstoatomclass)
        atomcommentstolistofsmartslist,atomindicestolistofatomcomments=MapIndicesToCommentsAtom(poltype,atomindextoallsmartspolar,smartstocomment,listofatomsforprm)
        atomindicestolistofatomcomments,atomcommentstoparameters=SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments)
        atomindicestocomments,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist)
        newpolarindicestopoltypeclasses,newpolarprms,newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist=GrabNewParametersPolarize(poltype,atomindicestocomments,atomcommentstoparameters,'polarize',atomindicestosmartslist,smartstoactualcomment) 
        smartstoatomclass, atomclasstoclassname, atomclasstocomment=ReadDatabaseSmartsMap(poltype,poltype.latestsmallmoleculesmartstotinkerclass) 
        atomindextoallsmarts,atomindextoallsmartsmatches=MatchAllSmartsToAtomIndices(poltype,smartstoatomclass)
        bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters,tortorsmartsatomordertoparameters,tortorsmartsatomordertogrid,smartsatomordertotorvdwdb=ReadExternalDatabase(poltype)
        smartsatomordertoelementtinkerdescrip=ReadSmallMoleculeLib(poltype,poltype.smallmoleculesmartstotinkerdescrip)
        elementtinkerdescriptotinkertype,tinkertypetoclass=GrabTypeAndClassNumbers(poltype,poltype.smallmoleculeprmlib)
        planarbonds=GrabPlanarBonds(poltype,listofbondsforprm,mol)
        bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses=MapIndicesToClasses(poltype,atomindextoallsmarts,smartstoatomclass,listofbondsforprm,listofanglesforprm,planarbonds)
        bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses,bondclassestoparameters,angleclassestoparameters,strbndclassestoparameters,opbendclassestoparameters=SearchForParameters(poltype,bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses)
        bondindicestoclasses,bondindicestosmartslist=FindBestSMARTSMatch(poltype,bondindicestolistofbondclasses,bondclassestolistofsmartslist)

 
        angleindicestoclasses,angleindicestosmartslist=FindBestSMARTSMatch(poltype,angleindicestolistofangleclasses,angleclassestolistofsmartslist)
 
        strbndindicestoclasses,strbndindicestosmartslist=FindBestSMARTSMatch(poltype,strbndindicestolistofstrbndclasses,strbndclassestolistofsmartslist)
        opbendindicestoclasses,opbendindicestosmartslist=FindBestSMARTSMatch(poltype,opbendindicestolistofopbendclasses,opbendclassestolistofsmartslist)
        opbendbondindicestotrigonalcenterbools=CheckTrigonalCenters(poltype,listofbondsforprm,mol)
        newangleindicestopoltypeclasses,newangleprms,newanglepoltypeclassestocomments,newanglepoltypeclassestosmartslist=GrabNewParameters(poltype,angleindicestoclasses,angleclassestoparameters,'angle',angleindicestosmartslist,atomclasstocomment) 
        newbondindicestopoltypeclasses,newbondprms,newbondpoltypeclassestocomments,newbondpoltypeclassestosmartslist=GrabNewParameters(poltype,bondindicestoclasses,bondclassestoparameters,'bond',bondindicestosmartslist,atomclasstocomment) 
        newstrbndindicestopoltypeclasses,newstrbndprms,newstrbndpoltypeclassestocomments,newstrbndpoltypeclassestosmartslist=GrabNewParameters(poltype,strbndindicestoclasses,strbndclassestoparameters,'strbnd',strbndindicestosmartslist,atomclasstocomment) 
        newopbendindicestopoltypeclasses,newopbendprms,newopbendpoltypeclassestocomments,newopbendpoltypeclassestosmartslist=GrabNewParameters(poltype,opbendindicestoclasses,opbendclassestoparameters,'opbend',opbendindicestosmartslist,atomclasstocomment,opbendbondindicestotrigonalcenterbools) 
 
        
        listofanglesthatneedplanarkeyword=CheckForPlanerAngles(poltype,listofanglesforprm,mol)
        parametersmartslist=GrabSMARTSList(poltype,smartsatomordertoelementtinkerdescrip)
        parametersmartstomaxcommonsubstructuresmarts,maxatomsize=FindMaximumCommonSubstructures(poltype,parametersmartslist,rdkitmol)
        if len(list(parametersmartstomaxcommonsubstructuresmarts.keys()))!=0:
             parametersmartslist=list(parametersmartstomaxcommonsubstructuresmarts.keys())
        indextoneighbidxs=FindAllNeighborIndexes(poltype,rdkitmol)
        bondindicestoextsmartsmatchlength,bondindicestoextsmarts,bondindicestoextsmartsatomorder,bondsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,bondsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        angleindicestoextsmartsmatchlength,angleindicestoextsmarts,angleindicestoextsmartsatomorder,anglesmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,anglesmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        strbndindicestoextsmartsmatchlength,strbndindicestoextsmarts,strbndindicestoextsmartsatomorder,strbndsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,strbndsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        torsionindicestoextsmartsmatchlength,torsionindicestoextsmarts,torsionindicestoextsmartsatomorders,torsionsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,torsionsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        opbendindicestoextsmartsmatchlength,opbendindicestoextsmarts,opbendindicestoextsmartsatomorder,opbendsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,opbendsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        vdwindicestoextsmartsmatchlength,vdwindicestoextsmarts,vdwindicestoextsmartsatomorder,vdwsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,vdwsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb,torsionindicestoextsmarts)
        tortorindicestoextsmartsmatchlength,tortorindicestoextsmarts,tortorindicestoextsmartsatomorders,tortorsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,tortorsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        parametersmartstordkitmol=GenerateRdkitMolObjectsParameterSMARTS(poltype,parametersmartslist)
        atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,atomindicesforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofatomsforprm,parametersmartslist,mol,parametersmartstordkitmol)
        bondsforprmtoparametersmarts,bondsforprmtosmarts,bondsforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofbondsforprm,parametersmartslist,mol,parametersmartstordkitmol)
        planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,planarbondsforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,planarbonds,parametersmartslist,mol,parametersmartstordkitmol)
        anglesforprmtoparametersmarts,anglesforprmtosmarts,anglesforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofanglesforprm,parametersmartslist,mol,parametersmartstordkitmol)
        planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,planaranglesforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofanglesforprm,parametersmartslist,mol,parametersmartstordkitmol)
        torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionsforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listoftorsionsforprm,parametersmartslist,mol,parametersmartstordkitmol)
        atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,vdwindicestoextsmarts,vdwindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,vdwindicestoextsmartsmatchlength,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,vdwindicestoextsmarts,atomindicesforprmtomatchallneighbs,vdwindicestoextsmartsatomorder)
        bondsforprmtoparametersmarts,bondsforprmtosmarts,bondindicestoextsmarts,bondindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,bondindicestoextsmartsmatchlength,bondsforprmtoparametersmarts,bondsforprmtosmarts,bondindicestoextsmarts,bondsforprmtomatchallneighbs,bondindicestoextsmartsatomorder)
        planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,planarbondindicestoextsmarts,planarbondindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,bondindicestoextsmartsmatchlength,planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,bondindicestoextsmarts,planarbondsforprmtomatchallneighbs,bondindicestoextsmartsatomorder)
        anglesforprmtoparametersmarts,anglesforprmtosmarts,angleindicestoextsmarts,angleindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,angleindicestoextsmartsmatchlength,anglesforprmtoparametersmarts,anglesforprmtosmarts,angleindicestoextsmarts,anglesforprmtomatchallneighbs,angleindicestoextsmartsatomorder)
        planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,planarangleindicestoextsmarts,planarangleindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,angleindicestoextsmartsmatchlength,planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,angleindicestoextsmarts,planaranglesforprmtomatchallneighbs,angleindicestoextsmartsatomorder)
        torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionindicestoextsmarts,torsionindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,torsionindicestoextsmartsmatchlength,torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionindicestoextsmarts,torsionsforprmtomatchallneighbs,torsionindicestoextsmartsatomorders)
        torsionindicestoextsmartsatomorders=TrimDictionary(poltype,torsionindicestoextsmartsatomorders,torsionindicestoextsmarts)
        atomindextotinkertype,atomindextotinkerclass,atomindextoparametersmartsatomorder,atomindextoelementtinkerdescrip,atomindextosmartsatomorder=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)
 
        bondindicestotinkertypes,bondindicestotinkerclasses,bondindicestoparametersmartsatomorders,bondindicestoelementtinkerdescrips,bondindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,bondsforprmtoparametersmarts,bondsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

 
        opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders=ExtractOpbendFromBond(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses)

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
        formissingvdwindicestosmartsatomorders=RemoveIndicesMatchedFromExternalDatabase(poltype,atomindextosmartsatomorder,vdwindicestoextsmarts) 
        bondindicestosmartsatomorders=originalbondindicestosmartsatomorders.copy()
        angleindicestosmartsatomorders=originalangleindicestosmartsatomorders.copy()

        totalbondscollector=FindAllConsecutiveRotatableBonds(poltype,mol,listofbondsforprm)
        tortorsmissing=FindMissingTorTors(poltype,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,rdkitmol,mol,indextoneighbidxs,totalbondscollector)
        torsionindicestosmartsatomorders=AddDictionaryItems(poltype,torsionindicestosmartsatomorders,torsionindicestoextsmartsatomorders)
        torsionsmissing,poormatchingaromatictorsions,poormatchingpartialaromatictorsions,torsionstozerooutduetocolinear=FindMissingTorsions(poltype,torsionindicestosmartsatomorders,rdkitmol,optmol,indextoneighbidxs)
        torsionsmissing=FindAdjacentMissingTorsionsForTorTor(poltype,torsionsmissing,totalbondscollector,tortorsmissing)
        atomindextosmartsatomorder=AddExternalDatabaseMatches(poltype, atomindextosmartsatomorder,vdwindicestoextsmarts,vdwsmartsatomordertoparameters)
        vdwmissing=FindMissingParameters(poltype,formissingvdwindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
        missingvdwatomindices=ReduceMissingVdwByTypes(poltype,vdwmissing)
        bondmissing=FindMissingParameters(poltype,formissingbondindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
        anglemissing=FindMissingParameters(poltype,formissingangleindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
        anglemissingindicestotinkerclasses,removedangleindices=PruneDictionary(poltype,anglemissing,angleindicestotinkerclasses)
        bondmissingindicestotinkerclasses,removedbondindices=PruneDictionary(poltype,bondmissing,bondindicestotinkerclasses)
        torsionsmissingindicestotinkerclasses,removedtorsionindices=PruneDictionary(poltype,torsionsmissing,torsionindicestotinkerclasses)
        torsionszerooutindicestotinkerclasses,removedzerotorsionindices=PruneDictionary(poltype,torsionstozerooutduetocolinear,torsionindicestotinkerclasses)
        atomtinkerclasstopoltypeclass=TinkerClassesToPoltypeClasses(poltype,atomindextotinkerclass)
        bondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondindicestotinkerclasses)
        arotorsionsmissingindicestotinkerclasses,removedaromissingtorsionindices=PruneDictionary(poltype,poormatchingaromatictorsions,torsionindicestotinkerclasses)
        partialarotorsionsmissingindicestotinkerclasses,removedpartialarotorsionmpartialarotorsionmissinggindices=PruneDictionary(poltype,poormatchingpartialaromatictorsions,torsionindicestotinkerclasses)
        planarbondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarbondindicestotinkerclasses,False)
        opbendtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,opbendbondindicestotinkerclasses)
        opbendtinkerclassestopoltypeclasses=AddReverseKeys(poltype,opbendtinkerclassestopoltypeclasses)
        opbendtinkerclassestotrigonalcenterbools=TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools)

        angletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,angleindicestotinkerclasses)

        planarangletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarangleindicestotinkerclasses)
        torsiontinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionindicestotinkerclasses)
        torsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionsmissingindicestotinkerclasses)
        torsionszeroouttinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionszerooutindicestotinkerclasses)

        arotorsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,arotorsionsmissingindicestotinkerclasses)
        partialarotorsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,partialarotorsionsmissingindicestotinkerclasses)
        anglemissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,anglemissingindicestotinkerclasses)
        bondmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondmissingindicestotinkerclasses)
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
        opbendprms,blankbondpoltypeclassestoparametersmartsatomorders,blankbondpoltypeclassestosmartsatomorders,blankbondpoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,opbendprms,newopbendindicestopoltypeclasses,'opbend',bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips)
        
        opbendpoltypeclassestosmartsatomorders=bondpoltypeclassestosmartsatomorders.copy()
        opbendpoltypeclassestoelementtinkerdescrips=bondpoltypeclassestoelementtinkerdescrips.copy()
        bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,bondprms,newbondindicestopoltypeclasses,'bond',bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips)
        strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,strbndprms,newstrbndindicestopoltypeclasses,'strbnd',anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips)
        bondprms,bondpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,bondprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'bond',bondindicestoextsmartsatomorder)
        angleprms,anglepoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,angleprms,angleindicestoextsmarts,anglesmartsatomordertoparameters,'angle',angleindicestoextsmartsatomorder)
        strbndprms,strbndpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,strbndprms,angleindicestoextsmarts,strbndsmartsatomordertoparameters,'strbnd',angleindicestoextsmartsatomorder)
        torsionprms,torsionpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,torsionprms,torsionindicestoextsmarts,torsionsmartsatomordertoparameters,'torsion',torsionindicestoextsmartsatomorder)
        opbendprms,opbendpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,opbendprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'opbend',bondindicestoextsmartsatomorder)
        vdwprms,vdwpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,vdwprms,vdwindicestoextsmarts,vdwsmartsatomordertoparameters,'vdw',vdwindicestoextsmartsatomorder)
        tortorprms=[]
        tortorprms,tortorpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,tortorprms,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,'tortors',tortorindicestoextsmartsatomorders,tortorsmartsatomordertogrid)
        anglemissingtinkerclassestopoltypeclasses=RemovePoltypeClassesFromNewMatches(poltype,anglemissingtinkerclassestopoltypeclasses,anglepoltypeclassestoparametersmartsatomorders)
        bondmissingtinkerclassestopoltypeclasses=RemovePoltypeClassesFromNewMatches(poltype,bondmissingtinkerclassestopoltypeclasses,bondpoltypeclassestoparametersmartsatomorders)
        missinganglepoltypeclasses=ExtractMissingPoltypeClasses(poltype,anglemissingtinkerclassestopoltypeclasses)
        missingbondpoltypeclasses=ExtractMissingPoltypeClasses(poltype,bondmissingtinkerclassestopoltypeclasses)
        strbndprms=ZeroOutMissingStrbnd(poltype,anglemissingtinkerclassestopoltypeclasses,strbndprms)
        angleprms=AssignAngleGuessParameters(poltype,anglemissingtinkerclassestopoltypeclasses,angleprms,indextoneighbidxs)
        bondprms=AssignBondGuessParameters(poltype,bondmissingtinkerclassestopoltypeclasses,bondprms,indextoneighbidxs)
        angleprms.extend(newangleprms)
        bondprms.extend(newbondprms)
        strbndprms.extend(newstrbndprms)
        opbendprms.extend(newopbendprms)
        angleprms=ModifyAngleKeywords(poltype,angleprms,planarangletinkerclassestopoltypeclasses)
        bondlistbabel=ConvertToBabelList(poltype,listofbondsforprm)
        anglelistbabel=ConvertToBabelList(poltype,listofanglesforprm)
        bondprms=AddOptimizedBondLengths(poltype,optmol,bondprms,bondlistbabel)
        angleprms=AddOptimizedAngleLengths(poltype,optmol,angleprms,anglelistbabel)
        
        torsionsmissingtinkerclassestopoltypeclasses=MergeDictionaries(poltype,torsionszeroouttinkerclassestopoltypeclasses,torsionsmissingtinkerclassestopoltypeclasses)
        extratormissingtinkerclassestopoltypeclasses=GenerateTinkerClassesToPoltypeClasses(poltype,removedtorsionindices) 
        torsionsmissingtinkerclassestopoltypeclasses=MergeDictionaries(poltype,extratormissingtinkerclassestopoltypeclasses,torsionsmissingtinkerclassestopoltypeclasses)

        torsionprms=ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms)
        torsionprms,arotorsionlinetodescrips=DefaultAromaticMissingTorsions(poltype,arotorsionsmissingtinkerclassestopoltypeclasses,partialarotorsionsmissingtinkerclassestopoltypeclasses,torsionprms,mol)
        torsionkeystringtoparameters=GrabTorsionParameterCoefficients(poltype,torsionprms)
        potentialmissingopbendprmtypes=FindPotentialMissingParameterTypes(poltype,opbendprms,planarbondtinkerclassestopoltypeclasses)
        potentialmissingopbendprmindices=ConvertPoltypeClassesToIndices(poltype,potentialmissingopbendprmtypes)
        potentialmissingopbendprmindices=FilterIndices(poltype,potentialmissingopbendprmindices,planarbonds)
        missingopbendprmindices=CheckIfParametersExist(poltype,potentialmissingopbendprmindices,opbendprms)
        torsionsmissing=ConvertToPoltypeClasses(poltype,torsionsmissing)
        missingvdwtypes=[poltype.idxtosymclass[i] for i in missingvdwatomindices]
        defaultvalues=None
        if len(missingopbendprmindices)!=0:
            newopbendprms,defaultvalues=DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools)
            opbendprms.extend(newopbendprms)
        opbendprms=AddTrivialOPBendForAmine(poltype,opbendprms,mol)
        if poltype.isfragjob==True and len(poltype.onlyrotbndslist)!=0:
            vdwprms,parentvdwtransferinfo=GrabVdwParametersFromParent(poltype,vdwprms)
        polarprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newpolarprms,{},{},{},{},newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        vdwprmstotransferinfo=MapParameterLineToTransferInfo(poltype,vdwprms,atompoltypeclasstoparametersmartsatomorder,atompoltypeclasstosmartsatomorder,atompoltypeclassestoelementtinkerdescrip,vdwpoltypeclassestosmartsatomordersext,{},{},arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        tortorprmstotransferinfo=MapParameterLineToTransferInfo(poltype,tortorprms,{},{},{},tortorpoltypeclassestosmartsatomordersext,{},{},arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses,defaultvalues=None,keyword='tortors')
        bondprmstotransferinfo=MapParameterLineToTransferInfo(poltype,bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips,bondpoltypeclassestosmartsatomordersext,newbondpoltypeclassestocomments,newbondpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        opbendprmstotransferinfo=MapParameterLineToTransferInfo(poltype,opbendprms,blankbondpoltypeclassestoparametersmartsatomorders,opbendpoltypeclassestosmartsatomorders,opbendpoltypeclassestoelementtinkerdescrips,opbendpoltypeclassestosmartsatomordersext,newopbendpoltypeclassestocomments,newopbendpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses,defaultvalues)
        
        angleprmstotransferinfo=MapParameterLineToTransferInfo(poltype,angleprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,anglepoltypeclassestosmartsatomordersext,newanglepoltypeclassestocomments,newanglepoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        strbndprmstotransferinfo=MapParameterLineToTransferInfo(poltype,strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,strbndpoltypeclassestosmartsatomordersext,newstrbndpoltypeclassestocomments,newstrbndpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        torsionprmstotransferinfo=MapParameterLineToTransferInfo(poltype,torsionprms,torsionpoltypeclassestoparametersmartsatomorders,torsionpoltypeclassestosmartsatomorders,torsionpoltypeclassestoelementtinkerdescrips,torsionpoltypeclassestosmartsatomordersext,{},{},arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        if poltype.isfragjob==True and len(poltype.onlyrotbndslist)!=0:
            vdwprmstotransferinfo=AddParentVdwTransferInfo(poltype,parentvdwtransferinfo,vdwprmstotransferinfo)
        WriteOutList(poltype,torsionsmissing,poltype.torsionsmissingfilename)
        WriteDictionaryToFile(poltype,torsionkeystringtoparameters,poltype.torsionprmguessfilename)
        WriteOutList(poltype,missingvdwatomindices,poltype.vdwmissingfilename)
        WriteOutList(poltype,tortorsmissing,poltype.tortormissingfilename)
        return bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,torsionsmissing,torsionkeystringtoparameters,missingvdwatomindices,soluteprms,tortorprmstotransferinfo,tortorsmissing
