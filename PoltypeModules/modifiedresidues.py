import os
import sys
import openbabel
import time
from rdkit import Chem
from MDAnalysis import Universe, Merge
from MDAnalysis.analysis.align import alignto
import re
import numpy
import shutil
import copy
from rdkit.Chem import rdFMCS
import databaseparser as db

def GenIndexToTypeIndexDic(poltype,tinkxyzfile):
    temp=open(tinkxyzfile,'r')
    xyzresults=temp.readlines()
    temp.close()
    idxtotypeidx={}
    for lineidx in range(len(xyzresults)):
        line=xyzresults[lineidx]
        linesplit=line.split()
        if lineidx!=0:
            idx=int(linesplit[0])
            typeidx=int(linesplit[5])
            idxtotypeidx[idx]=typeidx
    return idxtotypeidx

def IsItABoundaryAtom(poltype,neighbidxs,modproidxs,proOBmol):
    foundbound=False
    for neighbidx in neighbidxs:
        if neighbidx not in modproidxs: # then this is an atom still inside MCS on ligand side and we define this as boundary atom
            foundbound=True
    return foundbound

def GrabModifiedResidueProteinIndexes(poltype,modifiedproteinpdbname,knownresiduesymbs,totalatomnumber):
    temp=open(modifiedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    modproidxs=[]
    for line in results:
        if 'ATOM' in line:
            linesplit=line.split()
            PDBcode=linesplit[3]
            atomindex=int(linesplit[1])
            try:
                resnumber=linesplit[5]
                int(resnumber)
            except:
                resnumber=linesplit[4]
            if PDBcode not in knownresiduesymbs and atomindex!=totalatomnumber:
                proatomidx=int(linesplit[1])
                modproidxs.append(proatomidx)
                modresnumber=int(resnumber)
                modresiduelabel=PDBcode
    return modproidxs,modresnumber,modresiduelabel

def GrabPDBInfo(poltype,pdbfilename):
    # first grab the first residue number, then grab the last residue number
    temp=open(pdbfilename,'r')
    results=temp.readlines()
    temp.close()
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'ATOM' in line:
            atomidx=int(line[7:11+1])
            resnum=int(line[23:26+1])
            if atomidx==1:
                firstresnum=resnum
            lastresnum=resnum
    proidxtothreelettercode={}
    proidxtoatomlabel={}
    proidxtoresnum={}
    for line in results:
        if 'ATOM' in line:
            atomidx=int(line[6:10+1].lstrip().rstrip())
            atomlabel=line[12:15+1].lstrip().rstrip()
            resname=line[17:19+1].lstrip().rstrip()
            resnum=int(line[22:25+1].lstrip().rstrip())
            proidxtothreelettercode[atomidx]=resname
            proidxtoatomlabel[atomidx]=atomlabel
            proidxtoresnum[atomidx]=resnum
    return proidxtothreelettercode,proidxtoatomlabel,proidxtoresnum,firstresnum,lastresnum


def PreservePDBAtomLabels(poltype,pdbfilename,idxtoatomlabel):
    tempname=pdbfilename.replace('.pdb','_temp.pdb')
    ftemp=open(tempname,'w')
    temp=open(pdbfilename,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=re.split(r'(\s+)',line)
        if 'ATOM' in line:
            atomidx=int(line[6:10+1].lstrip().rstrip())
            oldatomlabel=line[12:15+1].lstrip().rstrip()
            if atomidx in idxtoatomlabel.keys():
                atomlabel=idxtoatomlabel[atomidx]
                oldlen=len(oldatomlabel)
                newlen=len(atomlabel)
                diff=newlen-oldlen
                magdiff=numpy.abs(diff)
                if diff<0:
                    atomlabel+=magdiff*' '
                elif diff>0:
                    space=linesplit[5]
                    spacelength=len(space)
                    newlength=spacelength-magdiff
                    linesplit[5]=space[:newlength]
                linesplit[4]=atomlabel
                newline=''.join(linesplit)
                ftemp.write(newline)
            else:
                ftemp.write(line)
        else:
            ftemp.write(line)
    temp.close()
    ftemp.close()
    os.remove(pdbfilename)
    os.rename(tempname,pdbfilename)
            

def GrabProteinTypeNumbers(poltype,pdbfilename,knownresiduesymbs,libpath,modproidxtotypenumber):

    proidxtothreelettercode,proidxtoatomlabel,proidxtoresnum,firstresnum,lastresnum=GrabPDBInfo(poltype,pdbfilename)

    proidxtotypeidx={}
    proidxtoboundtypeidx={}
 
    for proidx in proidxtoatomlabel.keys():
        atomlabel=proidxtoatomlabel[proidx]
        reslabel=proidxtothreelettercode[proidx]
        resnum=proidxtoresnum[proidx]
        if reslabel in knownresiduesymbs:
            if resnum==firstresnum:
                reslabel='N'+reslabel
            elif resnum==lastresnum:
                reslabel='C'+reslabel
            else:
                pass


        atomlabeltoprotypeidx=GrabTypeLibIdxs(poltype,reslabel,libpath)
        if atomlabel in atomlabeltoprotypeidx.keys():
            protypeidx=atomlabeltoprotypeidx[atomlabel]
        else:
            protypeidx=0
        if proidx in modproidxtotypenumber.keys():
            protypeidx=modproidxtotypenumber[proidx]
        if protypeidx==0:
            raise ValueError("Type number for index "+str(proidx)+' '+atomlabel+' was not found')
     
            
        proidxtotypeidx[proidx]=int(protypeidx)


    return proidxtotypeidx,firstresnum,lastresnum


def GrabTypeLibIdxs(poltype,reslabel,libpath):
    atomlabeltoprotypeidx={}
    temp=open(libpath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line:
            libreslabel=line.split()[4]
            if libreslabel==reslabel:
                linesplit=line.split()
                typenumber=linesplit[1]
                atomlabel=linesplit[3]
                atomlabeltoprotypeidx[atomlabel]=typenumber

    return atomlabeltoprotypeidx



def GenerateXYZFile(poltype,molstruct, structfname,idxtoatomlabel): # dont use idxtoatomlabel, just use element symbol
    tmpconv = openbabel.OBConversion()
    tmpfh = open(structfname, "w")
    iteratom = openbabel.OBMolAtomIter(molstruct)
    tmpfh.write('%6d\n' % (molstruct.NumAtoms()))
    etab = openbabel.OBElementTable()
    for ia in iteratom:
        linestring='%6d %2s %13.6f %11.6f %11.6f %5s' % (ia.GetIdx(), etab.GetSymbol(ia.GetAtomicNum()), ia.x(), ia.y(), ia.z(), '0')
        iteratomatom = openbabel.OBAtomAtomIter(ia)
        neighbors = []
        for iaa in iteratomatom:
            neighbors.append(iaa.GetIdx())
        neighbors = sorted(neighbors)
        neighbstring=''
        for iaa in neighbors:
            neighbstring+='%5d' % iaa+' '
        linestring+=' '
        linestring+=neighbstring
        tmpfh.write(linestring+'\n')
    tmpfh.close()

    return 


def GenerateIdxToAtomLabel(poltype,pdbname):
    idxtoatomlabel={}
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            idxtoatomlabel[atomindex]=atomlabel
    
    return idxtoatomlabel
                    
def GenerateProteinTinkerXYZFile(poltype,modifiedproteinpdbname,modproidxtotypenumber,amoebabioprmpath,proidxtoprotype,knownresiduesymbs):
    if not os.path.isfile(modifiedproteinpdbname):
        newpdbname=GeneratePDBTopologyFile(poltype,modifiedproteinpdbname)
    else:
        newpdbname=modifiedproteinpdbname 
    # now read pdb into openbabel structure and then also generate idxtoatomlabel dictionary
    obConversion = openbabel.OBConversion()
    molstruct = openbabel.OBMol()
    obConversion.SetInFormat('.pdb')
    obConversion.ReadFile(molstruct, newpdbname) 
    # now create idxtoatomlabel
    idxtoatomlabel=GenerateIdxToAtomLabel(poltype,newpdbname)
    structfname=newpdbname.replace('.pdb','.xyz')
    GenerateXYZFile(poltype,molstruct, structfname,idxtoatomlabel)
    proidxtoprotype,firstresnum,lastresnum=GrabProteinTypeNumbers(poltype,modifiedproteinpdbname,knownresiduesymbs,poltype.libpath,modproidxtotypenumber)
    temp=open(structfname,'r')
    results=temp.readlines()
    temp.close()
    protinkxyz='temp.xyz'
    temp=open(protinkxyz,'w')
    for line in results:
        nowhitespacelinesplit=line.split()
        if len(nowhitespacelinesplit)>1 and 'pdb' not in line: # then I must be in the atom block
            prewhitespacelinesplit=re.split(r'(\s+)', line) # I want to preserve white space when editing lines in tinker xyz file
            proidx=int(prewhitespacelinesplit[2])
            protype=proidxtoprotype[proidx]
            if len(str(protype))==2:
                string=' '+str(protype)
            elif len(str(protype))==1:
                string='  '+str(protype)
            else:
                string=str(protype)
            prewhitespacelinesplit[12]=string
            temp.write(''.join(prewhitespacelinesplit))
        else:
            temp.write(line)
    temp.close()
    os.remove(structfname)
    os.rename(protinkxyz,structfname)
    return structfname,firstresnum,lastresnum

def GrabProteinParameterTypeNumbers(poltype,amoebabioprmpath):
    proteinprmtypenumbers=[]
    temp=open(amoebabioprmpath,'r')
    results=temp.readline()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            typenumber=int(linesplit[1])
            proteinprmtypenumbers.append(typenumber)
    return proteinprmtypenumbers


def GrabLigandParameterTypeNumbers(poltype,key5fname):
    ligandtypenumbers=[]
    temp=open(key5fname,'r')
    results=temp.readline()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            typenumber=int(linesplit[1])
            ligandtypenumbers.append(typenumber)
         
    return ligandtypenumbers


def AssignNewTypeNumbers(poltype,modproidxs,ligidxtotypeidx,proidxtoligidx):
    modproidxtotypenumber={}
    for proidx in modproidxs:
        ligidx=proidxtoligidx[proidx]
        ligtype=ligidxtotypeidx[ligidx]
        modproidxtotypenumber[proidx]=ligtype

    return modproidxtotypenumber


def GrabTinkerClassNumber(poltype,boundtype,prmfile):
    temp=open(prmfile,'r')
    results=temp.readlines()
    temp.close()
    theatomclass=None
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            atomtype=linesplit[1]
            atomclass=linesplit[2]
            if str(boundtype)==atomtype:
                theatomclass=int(atomclass)
    return theatomclass


def GrabLigandTypeToBoundryTypeForKeyFile(poltype,listofvalterms,proidxtoligidx,ligidxtotypeidx,proidxtotypeidx,prmfile,modproidxs,mincorenumber,maxcorenumber): # be careful because need to convert type number to class number for valence only if it is a protein type number (poltype type numbers are the same as class numbers) terms!
    valligtypestoboundtypes={}
    for val in listofvalterms:
        valtype=[]
        for proidx in val:
            ligidx=proidxtoligidx[proidx]
            ligtype=ligidxtotypeidx[ligidx]
            valtype.append(ligtype)
        if tuple(valtype) not in valligtypestoboundtypes.keys() and tuple(valtype[::-1]) not in valligtypestoboundtypes.keys():
            valtup=tuple(valtype)
            boundtypes=[]
            for proidx in val:
                if proidx not in modproidxs: # this is normal part of protein
                    boundtype=proidxtotypeidx[proidx]
                else:
                    ligidx=proidxtoligidx[proidx]
                    boundtype=ligidxtotypeidx[ligidx]
                        
                if boundtype<mincorenumber or boundtype>maxcorenumber:
                    boundclass=GrabTinkerClassNumber(poltype,boundtype,prmfile)
                    boundtypes.append(boundclass)
                else:
                    boundtypes.append(boundtype)
            if valtype!=boundtypes:
                valligtypestoboundtypes[valtup]=boundtypes
            
    return valligtypestoboundtypes


def GrabPrmClassToBoundryClassForPrmFile(poltype,listofvalterms,proidxtoligidx,ligidxtotypeidx,proidxtotypeidx,prmfile,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber): # be careful because need to convert type number to class number for valence only if it is a protein type number (poltype type numbers are the same as class numbers) terms!
    prmclassestoboundclasses={}
    for val in listofvalterms:
        valtype=[]
        prmclasses=[]
        boundclasses=[]
        skip=False
        for proidx in val:
            protype=proidxtotypeidx[proidx]
            if protype<mincorenumber or protype>maxcorenumber: # if not poltype type number
                valtype.append(protype)
                prmclasses.append(prmtypetoprmclass[protype])
                boundclasses.append(prmtypetoprmclass[protype])
            else:
                if proidx in proboundidxtoprotype.keys():
                    valtype.append(proboundidxtoprotype[proidx])
                    prmclasses.append(prmtypetoprmclass[proboundidxtoprotype[proidx]])
                    boundclasses.append(protype)
                else:
                    skip=True
                    valtype.append(protype)
                    prmclasses.append(protype)
                    boundclasses.append(protype)
        if skip==True:
            continue
        if tuple(prmclasses) in prmclassestoboundclasses.keys() and tuple(prmclasses[::-1]) not in prmclassestoboundclasses.keys(): # then we need to make a list
            boundclassesprev=prmclassestoboundclasses[tuple(prmclasses)]
            if boundclasses not in boundclassesprev:
                prmclassestoboundclasses[tuple(prmclasses)].append(boundclasses)
        elif tuple(prmclasses[::-1]) in prmclassestoboundclasses.keys() and tuple(prmclasses) not in prmclassestoboundclasses.keys(): # then we need to make a list
            boundclassesprev=prmclassestoboundclasses[tuple(prmclasses[::-1])]
            if boundclasses not in boundclassesprev: 
                prmclassestoboundclasses[tuple(prmclasses[::-1])].append(boundclasses)
        else: 
            prmclassestoboundclasses[tuple(prmclasses)]=[]
            prmclassestoboundclasses[tuple(prmclasses)].append(boundclasses)
    return prmclassestoboundclasses
            


def GrabParametersFromKeyFile(poltype,key,torsionligtypestoboundtypesforkey,ligboundaryatomtypes,poltypetoprmtype,atomtypetoframedef,mincorenumber,maxcorenumber):
    
    temp=open(key,'r')
    results=temp.readlines()
    temp.close()
    torsionprms=[]
    polarizeprms=[]
    mpoleprms=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        linesplitall=re.split(r'(\s+)', line)

        
        if 'torsion' in line and 'unit' not in line and '#' not in line:
            torsiontypelist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
            foundtorsion=False
               
            if tuple(torsiontypelist) in torsionligtypestoboundtypesforkey.keys():
                torsiontup=tuple(torsiontypelist)
                foundtorsion=True
            elif tuple(torsiontypelist[::-1]) in torsionligtypestoboundtypesforkey.keys():
                torsiontup=tuple(torsiontypelist[::-1])
                foundtorsion=True
            if foundtorsion==True:
                allgreaterthanpoltypenums=True
                boundtypes=torsionligtypestoboundtypesforkey[torsiontup]
                allgreaterthanpoltypenums=True # if all numbers are >=401 then will already be in poltype key file, no need to add again
                for tortype in boundtypes:
                    if tortype<mincorenumber or tortype>maxcorenumber:
                        allgreaterthanpoltypenums=False
                if allgreaterthanpoltypenums==True:
                    foundtorsion=False
            if foundtorsion==True:
                             
                if linesplitall[0]=='': # sometimes when torsion is added back to the key file, it has a space in front of it
                    linesplitall[4]=str(boundtypes[0])    
                    linesplitall[6]=str(boundtypes[1])  
                    linesplitall[8]=str(boundtypes[2])
                    linesplitall[10]=str(boundtypes[3])
                else:
                    linesplitall[2]=str(boundtypes[0])    
                    linesplitall[4]=str(boundtypes[1])  
                    linesplitall[6]=str(boundtypes[2])
                    linesplitall[8]=str(boundtypes[3])


                newline=''.join(linesplitall)
                torsionprms.append(newline)
 
 
        elif 'polarize' in line: # only want to add neighbors to polline that are already in original poltype definition
            atomtype=int(linesplit[1])
            if atomtype in poltypetoprmtype.keys():
                prmtype=poltypetoprmtype[atomtype]
                newline=line.replace('\n','')+' '+str(prmtype)+'\n'
                polarizeprms.append(newline)
                
 
        elif 'multipole' in line:
            atomtype=int(linesplit[1])
            if atomtype in atomtypetoframedef.keys():
                chgpartofline=linesplitall[-4:] # include space
                phrasepartofline=linesplitall[:2] # include space
                spacetoadd='   '
                newlist=[]
                newlist.extend(phrasepartofline)
                framedef=atomtypetoframedef[atomtype]
                for typenum in framedef:
                    newlist.append(str(typenum))
                    newlist.append(spacetoadd)
                newlist=newlist[:-1] # remove last space then add space before chg
                newlist.extend(chgpartofline)
                newline=''.join(newlist)
                mpolelist=[newline,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
                for mpoleline in mpolelist:
                    mpoleprms.append(mpoleline)


    
    return torsionprms,polarizeprms,mpoleprms


def CountNumberPoltypeNums(poltype,typelist,mincorenumber,maxcorenumber):
    count=0
    for typenum in typelist:
        if typenum>=mincorenumber and typenum<=maxcorenumber:
            count+=1
    return count

def CountNumberOfNonPoltypeNums(poltype,typelist,mincorenumber,maxcorenumber):
    count=0
    for typenum in typelist:
        if typenum<mincorenumber or typenum>maxcorenumber:
            count+=1
    return count


def GrabIdxsFromType(poltype,ligtypes,ligidxtotypeidx):
    ligidxs=[]
    typenums=[]
    for typenum in ligtypes:
        for ligidx in ligidxtotypeidx.keys():
            ligtype=ligidxtotypeidx[ligidx]
            if typenum==ligtype and typenum not in typenums:
                ligidxs.append(ligidx)
                typenums.append(typenum)
                
    return ligidxs



def GrabMultipoleFrameDefintions(poltype,key,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,modproidxs,ligidxtoproidx,prosideboundligidx,proboundidxtoprotype,mincorenumber,maxcorenumber,modligidxs):
    boundaryatomtypeidxs=[ligidxtotypeidx[i] for i in boundaryatomidxs]
    temp=open(key,'r')
    results=temp.readlines()
    temp.close()
    atomtypetoframedef={}
    proteintypestoframedefforprmfile={}
    for line in results:
        linesplit=line.split()
        if 'multipole' in line:
            atomidx=int(linesplit[1])
            removedchg=linesplit[:-1]
            removedstring=removedchg[1:]
            frames=[int(i.replace('-','')) for i in removedstring]
            
            numberofboundatoms=CountNumberOfAtomsInBoundryList(poltype,frames,boundaryatomtypeidxs)
            if numberofboundatoms==1:
                protypes=[]
                ligidxs=GrabIdxsFromType(poltype,frames,ligidxtotypeidx)
                if ligidxs[0] in modligidxs:
                    firstprotype=proidxtoprotype[ligidxtoproidx[ligidxs[0]]]
                    if firstprotype>=mincorenumber and firstprotype<=maxcorenumber:
                        for i in range(len(ligidxs)):
                            ligidx=ligidxs[i]
                            proidx=ligidxtoproidx[ligidx]
                            string=removedstring[i]
                            protype=proidxtoprotype[proidx]
                            if '-' in string:
                                protype=-protype
                            protypes.append(protype)
                        numberofpoltypetypes=CountNumberPoltypeNums(poltype,protypes,mincorenumber,maxcorenumber)
                        if numberofpoltypetypes!=len(protypes):
                            atomtypetoframedef[firstprotype]=protypes
                    else:
                        allprotypelist=[]
                        for i in range(len(ligidxs)):
                            ligidx=ligidxs[i]
                            proidx=ligidxtoproidx[ligidx]
                            string=removedstring[i]
                            if proidx in proboundidxtoprotype.keys():
                                protype=proboundidxtoprotype[proidx]
                            else:
                                protype=proidxtoprotype[proidx]
                            if '-' in string:
                                protype=-protype
                            protypes.append(protype)
                            allprotypelist.append(proidxtoprotype[proidx])
                            numberofnonpoltypetypes=CountNumberOfNonPoltypeNums(poltype,protypes,mincorenumber,maxcorenumber)
                            if numberofnonpoltypetypes!=len(protypes):
                                proteintypestoframedefforprmfile[tuple(allprotypelist)]=protypes
    return atomtypetoframedef,proteintypestoframedefforprmfile


def WriteNewKeyFile(poltype,atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms,writekey,amoebabioprmpath):
    # assumes that there is white space at end of every parameter block, just add a line to key_5
    temp=open(writekey,'a')
    temp.write('\n\n\n\n\n\n\n\n\n\n')
    temp.close()
    temp=open(writekey,'r')
    results=temp.readlines()
    temp.close()
    lastidx=len(results)-1
    temp=open(writekey,'w')
    linelist=[]
    firstpassatomdef=True
    firstpassbond=True
    firstpassangle=True
    firstpasstorsion=True
    firstpassstrbnd=True
    firstpasspolarize=True
    firstpassvdw=True
    firstpassvdw=True
    firstpassmpole=True
    firstpassopbend=True
    linelist.append('parameters '+poltype.ModifiedResiduePrmPath+'\n')
    passedatomdefblock=False
    passedbondblock=False
    passedangleblock=False
    passedtorsionblock=False
    passedstrbndblock=False
    passedpolarizeblock=False
    passedvdwblock=False
    passedvdwblock=False
    passedmpoleblock=False
    passedopbendblock=False

    
    for lineidx in range(len(results)):
        line=results[lineidx]
        prewhitespacelinesplit=re.split(r'(\s+)', line)
        if 'atom' in line and firstpassatomdef==False:
            firstpassatomdef=True
            if line not in linelist:
                linelist.append(line)
        elif 'atom' not in line and firstpassatomdef==True and passedatomdefblock==False:
            passedatomdefblock=True
            for atomdefline in atomdefs:
                linelist.append(atomdefline)
            if line not in linelist:
                linelist.append(line)
        elif 'bond' in line and firstpassbond==False:
            firstpassbond=True
            if line not in linelist:
                linelist.append(line)
        elif 'bond' not in line and firstpassbond==True and passedbondblock==False:
            passedbondblock=True
            for bondline in bondprms:
                linelist.append(bondline)
            if line not in linelist:
                linelist.append(line)
        elif 'angle' in line and firstpassangle==False:
            firstpassangle=True
            if line not in linelist:
                linelist.append(line)
        elif 'angle' not in line and firstpassangle==True and passedangleblock==False:
            passedangleblock=True
            for angleline in angleprms:
                linelist.append(angleline)
            if line not in linelist:
                linelist.append(line)
        elif 'torsion' in line and firstpasstorsion==False:
            firstpasstorsion=True
            if line not in linelist:
                linelist.append(line)
        elif 'torsion' not in line and firstpasstorsion==True and passedtorsionblock==False:
            passedtorsionblock=True
            for torsionline in torsionprms:
                linelist.append(torsionline) 
            if line not in linelist:
                linelist.append(line)           
        elif 'strbnd' in line and firstpassstrbnd==False:
            firstpassstrbnd=True
            if line not in linelist:
                linelist.append(line)
        elif 'strbnd' not in line and firstpassstrbnd==True and passedstrbndblock==False:
            passedstrbndblock=True
            for strbndline in strbndprms:
                linelist.append(strbndline)
            if line not in linelist:
                linelist.append(line)
        elif 'polarize' in line and firstpasspolarize==False:
            firstpasspolarize=True
            if line not in linelist:
                linelist.append(line)
        elif 'polarize' not in line and firstpasspolarize==True and passedpolarizeblock==False:
            passedpolarizeblock=True
            for polarizeline in polarizeprms:
                linelist.append(polarizeline)
            if line not in linelist:
                linelist.append(line)
        elif 'opbend' in line and firstpassopbend==False:
            firstpassopbend=True
            if line not in linelist:
                linelist.append(line)
        elif 'opbend' not in line and firstpassopbend==True and passedopbendblock==False:
            passedopbendblock=True
            for opbendline in opbendprms:
                linelist.append(opbendline)
            if line not in linelist:
                linelist.append(line)
        elif 'vdw' in line and firstpassvdw==False:
            firstpassvdw=True
            if line not in linelist:
                linelist.append(line)
        elif 'vdw' not in line and firstpassvdw==True and passedvdwblock==False:
            passedvdwblock=True
            for vdwline in vdwprms:
                linelist.append(vdwline)
            if line not in linelist:
                linelist.append(line)
        elif 'multipole' in line and firstpassmpole==False:
            firstpassmpole=True
            linelist.append(line)
        elif 'multipole' not in line and firstpassmpole==True and passedmpoleblock==False and line=='\n':
            passedmpoleblock=True
            for mpoleline in mpoleprms:
                linelist.append(mpoleline)
            if line not in linelist:
                linelist.append(line)

        else:
            if line not in linelist:
                linelist.append(line)

    for line in linelist:
        temp.write(line)
        temp.flush()
        os.fsync(temp.fileno())
    temp.close()



def GenerateFragBabel(poltype,molindexlist,pdbname):
    atomidxstoremove=[]
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,pdbname)

    atomiter=openbabel.OBMolAtomIter(pdbmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if atomidx not in molindexlist:
            atomidxstoremove.append(atomidx)
        
    sortedatomidxstoremove=sorted(atomidxstoremove, reverse=True)
    for index in sortedatomidxstoremove:
        atom=pdbmol.GetAtom(index)
        pdbmol.DeleteAtom(atom)
    return pdbmol

def GenerateFragRdkit(poltype,molindexlist,rdkitmol,pdbname):
    conf = rdkitmol.GetConformer()
    atomidxstoremove=[]
    for atom in rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx not in molindexlist:
            atomidxstoremove.append(atomidx)
    fragmol=Chem.rdmolfiles.MolFromPDBFile(pdbname,removeHs=False,sanitize=False)
    rwmol = Chem.RWMol(fragmol)
    sortedatomidxstoremove=sorted(atomidxstoremove,reverse=True)
    for index in sortedatomidxstoremove:
        rwmol.RemoveAtom(index)
    return rwmol


def OldIndexToNewIndexBabel(poltype,pdbname,fragmol):
    parentindextofragindex={}
    positiontoparentindex={}
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,pdbname)

    atomiter=openbabel.OBMolAtomIter(pdbmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        pos=[atom.GetX(),atom.GetY(),atom.GetZ()]
        positiontoparentindex[tuple(pos)]=atomidx

    fragiter=openbabel.OBMolAtomIter(fragmol)
    for atom in fragiter:
        atomindex=atom.GetIdx()
        pos=tuple([atom.GetX(),atom.GetY(),atom.GetZ()])
        if pos in positiontoparentindex.keys():
            parentindex=positiontoparentindex[pos]
            parentindextofragindex[parentindex]=atomindex
    return parentindextofragindex

def OldIndexToNewIndexRdkit(poltype,pdbname,fragmol):
    parentindextofragindex={}
    positiontoparentindex={}
    pdbmol=Chem.rdmolfiles.MolFromPDBFile(pdbname,removeHs=False,sanitize=False)
    pdbconf = pdbmol.GetConformer()
    for atom in pdbmol.GetAtoms():
        atomidx=atom.GetIdx()
        position = pdbconf.GetAtomPosition(atomidx)
        pos=[position.x,position.y,position.z]
        positiontoparentindex[tuple(pos)]=atomidx
    fragconf = fragmol.GetConformer()
    for atom in fragmol.GetAtoms():
        atomindex=atom.GetIdx()
        position = fragconf.GetAtomPosition(atomindex)
        pos=tuple([position.x,position.y,position.z])
        if pos in positiontoparentindex.keys():
            parentindex=positiontoparentindex[pos]
            parentindextofragindex[parentindex]=atomindex
    return parentindextofragindex



def NeighboringModResidueAtomIndexes(poltype,modresnumber,modifiedproteinpdbname,modproidxs,proOBmol,proboundidxs):
    neighbidxs=[]
    temp=open(modifiedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    for lineindex in range(len(results)):
        line=results[lineindex]
        if 'ATOM' in line:
            linesplit=line.split()
            resnumber=int(line[23:26+1])
            proatomidx=int(linesplit[1])
            if proatomidx in proboundidxs:
                proatom=proOBmol.GetAtom(proatomidx)
                atomatomiter=openbabel.OBAtomAtomIter(proatom)
                for natom in atomatomiter:
                    natomidx=natom.GetIdx()
                    if natomidx not in modproidxs and natomidx not in neighbidxs:
                        neighbidxs.append(natomidx)
                    newatomatomiter=openbabel.OBAtomAtomIter(natom)
                    for nnatom in newatomatomiter:
                        nnatomidx=nnatom.GetIdx()
                        if nnatomidx not in modproidxs and nnatomidx not in neighbidxs and nnatom.GetAtomicNum()==8: #dont chop the carbonyl oxygen
                            neighbidxs.append(nnatomidx)

    return neighbidxs           


def FindBoundaryAtomIdxs(poltype,proOBmol,modproidxs):
    proboundaryatomidxs=[]
    for atom in openbabel.OBMolAtomIter(proOBmol):
        atomidx=atom.GetIdx()
        if atomidx in modproidxs:
            proatom=proOBmol.GetAtom(atomidx)
            neighbs=[natom for natom in openbabel.OBAtomAtomIter(proatom)]
            neighbidxs=[i.GetIdx() for i in neighbs]
            foundbound=IsItABoundaryAtom(poltype,neighbidxs,modproidxs,proOBmol)
            if foundbound==True:
                proboundaryatomidxs.append(atomidx)
    return proboundaryatomidxs

def CountNumberOfAtomsInModifiedResidue(poltype,atomidxlist,modproidxs):
    count=0
    for idx in atomidxlist:
        if idx in modproidxs:
            count+=1
    return count

def CountNumberOfAtomsInBoundryList(poltype,atomidxlist,proboundidxs):
    count=0
    for idx in atomidxlist:
        if idx in proboundidxs:
            count+=1
    return count


def GrabAtomsForValenceTermsAcrossBoundary(poltype,ligOBmol,proOBmol,ligidxtoproidx,proboundidxs,boundaryatomidxs,modproidxs):
    # now we can define arrays to collect bonds, angles and torsions
    listoftorsionsforkey=[]
    listofbondsforprm=[]
    listofanglesforprm=[]
    listoftorsionsforprm=[]
    for atom in openbabel.OBMolAtomIter(ligOBmol):
        atomidx=atom.GetIdx()
        if atomidx in ligidxtoproidx.keys():
            proidx=ligidxtoproidx[atomidx]
            proatom=proOBmol.GetAtom(proidx)
            neighbs=[natom for natom in openbabel.OBAtomAtomIter(proatom)]
            for natom in neighbs:
                nidx=natom.GetIdx()
                bondset=[nidx,proidx]
                numberofmodatoms=CountNumberOfAtomsInModifiedResidue(poltype,bondset,modproidxs) # redundant for bond
                numberofboundatoms=CountNumberOfAtomsInBoundryList(poltype,bondset,proboundidxs)
                if numberofmodatoms==1 and numberofboundatoms==1: 
                    if bondset not in listofbondsforprm and bondset[::-1] not in listofbondsforprm:
                        listofbondsforprm.append(bondset)
                nextneighbs=[nextatom for nextatom in openbabel.OBAtomAtomIter(natom)]
                for nextneighb in nextneighbs:
                    nextneighbidx=nextneighb.GetIdx()
                    if nextneighbidx!=proidx:
                        angleset=[nextneighbidx,nidx,proidx]
                        numberofboundatoms=CountNumberOfAtomsInBoundryList(poltype,angleset,proboundidxs)
                        numberofmodatoms=CountNumberOfAtomsInModifiedResidue(poltype,angleset,modproidxs)
                        if numberofboundatoms==1 and angleset not in listofanglesforprm and angleset[::-1] not in listofanglesforprm:
                            listofanglesforprm.append(angleset)
                        nextnextneighbs=[nextnextatom for nextnextatom in openbabel.OBAtomAtomIter(nextneighb)]
                        for nextnextatom in nextnextneighbs:
                            nextnextatomidx=nextnextatom.GetIdx()
                            if nextnextatomidx!=nidx:                
                                torsionset=[nextnextatomidx,nextneighbidx,nidx,proidx]
                                numberofboundatoms=CountNumberOfAtomsInBoundryList(poltype,torsionset,proboundidxs)
                                numberofmodatoms=CountNumberOfAtomsInModifiedResidue(poltype,torsionset,modproidxs)
                                if numberofboundatoms>=1 and torsionset not in listoftorsionsforkey and torsionset[::-1] not in listoftorsionsforkey and torsionset not in listoftorsionsforprm and torsionset[::-1] not in listoftorsionsforprm:
                                    if torsionset[0] in proboundidxs or torsionset[-1] in proboundidxs:
                                        if numberofmodatoms==4:
                                            listoftorsionsforkey.append(torsionset)
                                        else: # then number of mod atoms is 1 since boundary atom is in modproidxs
                                            listoftorsionsforprm.append(torsionset)
                                    else: # then the boundry atom must be in the middle two
                                        listoftorsionsforprm.append(torsionset)
    return listoftorsionsforkey,listofbondsforprm,listofanglesforprm,listoftorsionsforprm


def GrabPolarizeBoundaryLigTypeToProType(poltype,ligOBmol,proOBmol,modproidxs,boundaryatomidxs,prosideboundligidx):
    ligtypes=[] # need to know list of indexes that polarize terms needs to be changed
    # be careful we also want to use protein types instead of ligand types for the atom type being polarized if across the boundary on the protein side
    ligtypetoprotype={} # use this for converting ligtypes that should be protypes to protypes
    for atom in openbabel.OBMolAtomIter(ligOBmol):
        atomidx=atom.GetIdx()
        if atomidx in boundaryatomidxs or atomidx in prosideboundligidx:
            proidx=ligidxtoproidx[atomidx]
            proatom=proOBmol.GetAtom(proidx)
            proneighbs=[natom for natom in openbabel.OBAtomAtomIter(proatom)]
            proneighbidxs=[i.GetIdx() for i in proneighbs]
            ligtypenum=ligidxtotypeidx[atomidx]
            protypenum=proidxtoprotype[proidx]
            if proidx not in modproidxs:
                ligtypetoprotype[ligtypenum]=protypenum

            ligtypes.append(ligtypenum)
            idxarray=[]
            for proneighbidx in proneighbidxs:
                if proneighbidx in modproidxs: # then use ligand type number
                    idx=proidxtoligidx[proneighbidx]
                    neighbtypenum=ligidxtotypeidx[idx]
                    idxarray.append(idx)
                else:
                    neighbtypenum=proidxtoprotype[proneighbidx]
                    ligtypetoprotype[ligtypenum]=neighbtypenum
                    idxarray.append(proneighbidx)
            


    return ligtypes,ligtypetoprotype


def GrabKnownResidueSymbs(poltype): # if adding parameters to the library then need to 
    knownresiduesymbs=['A', 'G', 'U', 'A3', 'A5','DA', 'DC', 'DG', 'DT', 'G3', 'G5', 'U3', 'U5','DA3', 'DA5', 'DC3', 'DC5', 'DG3', 'DG5', 'DT3', 'DT5','ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIS','ILE', 'LEU', 'LYD', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYD', 'TYR', 'VAL', 'CYD', 'CYS']
    
    return knownresiduesymbs

def MatchCorrectProteinAtomsToCorrespondingHydrogen(poltype,proboundidxs,proOBmol,polOBmol,proidxtoligidx,ligidxtoproidx):
    atomiter=openbabel.OBMolAtomIter(proOBmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        atomatomiter=openbabel.OBAtomAtomIter(atom)
        for natom in atomatomiter:
            natomidx=natom.GetIdx()
    addedproatoms=[]
    for proidx in proboundidxs:
        proatom=proOBmol.GetAtom(proidx)
        atomatomiter=openbabel.OBAtomAtomIter(proatom)
        for proneighbatom in atomatomiter:
            proneighbatomidx=proneighbatom.GetIdx()
            if proneighbatomidx not in proidxtoligidx.keys(): # then this must be the atom corresponding to one of the Hydrogens
                ligidx=proidxtoligidx[proidx]
                ligatom=polOBmol.GetAtom(ligidx)
                atomatomiterlig=openbabel.OBAtomAtomIter(ligatom)
                for ligneighb in atomatomiterlig:
                    ligneighbidx=ligneighb.GetIdx()
                    if ligneighbidx not in proidxtoligidx.values(): # then this is the corresponding added hydrogen
                        proidxtoligidx[proneighbatomidx]=ligneighbidx
                        addedproatoms.append(proneighbatomidx)
    return proidxtoligidx,addedproatoms


def GrabCarbonylCarbonAndAmideNitrogen(poltype,modresnumber,pdbname):
    backboneprobound=[]
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())    
            if atomindex!=1:
                resnumber=int(line[23:26+1])
                nextline=results[lineidx+1]
                nextatomlabel=nextline[12:16].rstrip().lstrip()
                prevline=results[lineidx-1]
                prevatomlabel=prevline[12:16].rstrip().lstrip()

                if resnumber==modresnumber:
                    if atomlabel=='C' and prevatomlabel=='CA':
                        backboneprobound.append(atomindex)
                    if atomlabel=='N' and nextatomlabel=='CA':
                        backboneprobound.append(atomindex)
                if len(backboneprobound)==2:
                    break

    return backboneprobound



def CheckIfSMARTSInLibrary(poltype,SMARTSToTypelibpath,smarts):
    check=False
    if os.path.isfile(SMARTSToTypelibpath):
        temp=open(SMARTSToTypelibpath,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if smarts in line:
                check=True
    return check

def GrabSideChainIndexes(poltype,residuenumber,pdbname):
    sidechainindexes=[]
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            if resnumber==residuenumber:
               if atomlabel!='N' and atomlabel!='HN' and atomlabel!='CA' and atomlabel!='C' and atomlabel!='O' and atomlabel!='H' and atomlabel!='HA':
                   sidechainindexes.append(atomindex)
    return sidechainindexes

def GrabPDBIndexes(poltype,residuenumber,pdbname):
    indexes=[]
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            if resnumber==residuenumber:
               indexes.append(atomindex)
    return indexes

def FindAtomIndexesToKeepAfterMutation(poltype,residuenumber,atomindexlist,pdbname):
    indexestokeep=[]
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            if atomindex not in atomindexlist: 
                indexestokeep.append(atomindex)
    return indexestokeep

def ModifyAlignmentPDB(poltype,pdbname):
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    tempname=pdbname.replace('.pdb','_temp.pdb')
    temp=open(tempname,'w')
    for line in results:
        if 'HETATM' in line:
            linesplitall=re.split(r'(\s+)', line)
            linesplitall[6]=poltype.modifiedresiduepdbcode
            linesplitall[8]=poltype.mutatedresiduenumber
            if len(poltype.mutatedresiduenumber)==1:
                pass
            elif len(poltype.mutatedresiduenumber)==2:
                spaces=linesplitall[7]
                linesplitall[7]=spaces[:-1]
                otherspaces=linesplitall[9]
                linesplitall[9]=otherspaces[:-1]
            elif len(poltype.mutatedresiduenumber)==3:
                spaces=linesplitall[7]
                linesplitall[7]=spaces[:-2]
                otherspaces=linesplitall[9]
                linesplitall[9]=otherspaces[:-2]
            newline=''.join(linesplitall)
            temp.write(newline)
        else:
            temp.write(line)
    temp.close()
    os.remove(pdbname)
    os.rename(tempname,pdbname)

def FindCorrectReferenceIndexes(poltype,allindexes,pdbindexes):
    backboneindexesreference=[]
    for match in allindexes:
        foundmatch=True
        for index in match:
            if index not in pdbindexes:
                foundmatch=False
        if foundmatch==True:
            return match
    if foundmatch==False:
        return None

def GrabModifiedResidueLines(poltype,pdbname,resid):
    lines=[]
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    passedfirstres=False
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            reslabel=line[17:20]
            chainid=line[21]
            if resid==resnumber and passedfirstres==False and reslabel==poltype.modifiedresiduepdbcode:
                lines.append(line)
            elif resnumber==1:
                passedfirstres=True
    return lines

def GrabChainId(poltype,resid,pdbname):
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    passedfirstres=False
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            chainid=line[21]
            if resnumber==resid:
                return chainid

            
def ModifyChainIdentifier(poltype,lines,chainid):
    newlines=[]
    for line in lines:
        newline=line[:21]+chainid+line[22:]
        newlines.append(newline)
    return newlines
            

def MoveModifiedLines(poltype,lines,resid,pdbname,firstresnumber):
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    passedfirstres=False
    tempname=pdbname.replace('.pdb','_temp.pdb')
    temp=open(tempname,'w')
    wrotelines=False
    for line in results:
        if 'ATOM' in line and 'UNL' not in line: # sometimes UNL atoms from next amino acid (like an alpha carbon that does not continue on chain at end of file) needs to be removed
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            chainid=line[21]
            reslabel=line[17:20]
            if resnumber==firstresnumber:
                passedfirstres=True
            if resnumber==resid+1 and wrotelines==False:
                wrotelines=True
                for l in lines:
                    temp.write(l)
            if resnumber==resid and passedfirstres==True:
                newline=line[:17]+poltype.modifiedresiduepdbcode+line[20:]
                temp.write(newline)
            if passedfirstres==True and resnumber!=resid:
                temp.write(line)
        elif 'ATOM' in line and 'UNL' in line:
            pass    
        else:
            temp.write(line)
    temp.close()
    os.remove(pdbname)
    os.rename(tempname,pdbname)


def DeleteHydrogensAttachedToBackbone(poltype,molobj,indexlist):
    dellist=[]
    for index in indexlist:
        atom=molobj.GetAtom(index)
        atomatomiter=openbabel.OBAtomAtomIter(atom)
        for natom in atomatomiter:
            if natom.GetAtomicNum()==1:
                dellist.append(natom.GetIdx())
    dellist.sort(reverse=True)
    for index in dellist:
        atom=molobj.GetAtom(index)
        molobj.DeleteAtom(atom)
    return molobj 


def SwapDictionaryKeys(poltype,refdic,swapdic):
    newdic={}
    for key in list(refdic.keys()):
        if key in swapdic.keys():
            newkey=swapdic[key]
            value=refdic[key]
            newdic[newkey]=value
    return newdic

def MissingIndexes(poltype,proidxtoatomlabel):
    missingproidxs=[]
    minproidx=min(list(proidxtoatomlabel.keys()))
    maxproidx=max(list(proidxtoatomlabel.keys()))
    for i in range(minproidx,maxproidx+1):
        if i not in proidxtoatomlabel.keys():
            missingproidxs.append(i)
    return missingproidxs

    


def GenerateModifiedProteinPDB(poltype):      
    sidechainindexes=GrabSideChainIndexes(poltype,int(poltype.mutatedresiduenumber),poltype.unmodifiedproteinpdbname)
    atomindexestokeepafter=FindAtomIndexesToKeepAfterMutation(poltype,int(poltype.mutatedresiduenumber),sidechainindexes,poltype.unmodifiedproteinpdbname) # only do this for indexes after residue
    unmodpdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.SetOutFormat('pdb')
    obConversion.ReadFile(unmodpdbmol,poltype.unmodifiedproteinpdbname)
    proidxtothreelettercode,proidxtoatomlabel,proidxtoresnum,firstresnum,lastresnum=GrabPDBInfo(poltype,poltype.unmodifiedproteinpdbname)
    fragmolafter=GenerateFragBabel(poltype,atomindexestokeepafter,poltype.unmodifiedproteinpdbname)
    unmodidxtonewidx=OldIndexToNewIndexBabel(poltype,poltype.unmodifiedproteinpdbname,fragmolafter)
    refname='RefAfter.pdb'
    obConversion.WriteFile(fragmolafter,refname)

    proidxtoatomlabel=SwapDictionaryKeys(poltype,proidxtoatomlabel,unmodidxtonewidx)

    missingproidxs=MissingIndexes(poltype,proidxtoatomlabel)
    PreservePDBAtomLabels(poltype,refname,proidxtoatomlabel)
    inFormat = obConversion.FormatFromExt(poltype.mutatedsidechain)
    obConversion.SetInFormat(inFormat)   

     
    sidechainmol=openbabel.OBMol()
    obConversion.ReadFile(sidechainmol,poltype.mutatedsidechain)
    charge=sidechainmol.GetTotalCharge()

    backbonesmiles='C(C=O)N'
        
    backboneindexessidechain=MatchSMARTSToOBMol(poltype,backbonesmiles,sidechainmol)
    sidechainmol=DeleteHydrogensAttachedToBackbone(poltype,sidechainmol,backboneindexessidechain) # we will use H on backbone from referernce pdb not mutated sidechain pdb
    backboneindexessidechain=MatchSMARTSToOBMol(poltype,backbonesmiles,sidechainmol)
    refmol=openbabel.OBMol()
    obConversion.SetInFormat('pdb')
    
    obConversion.ReadFile(refmol,refname)
    pdbname='MutatedSideChain.pdb'
    obConversion.WriteFile(sidechainmol,pdbname)
    
    ModifyAlignmentPDB(poltype,pdbname)

    allindexes=MatchSMARTSToOBMolGrabAllMatches(poltype,backbonesmiles,refmol) 
    pdbindexes=GrabPDBIndexes(poltype,int(poltype.mutatedresiduenumber),refname)
    backboneindexesreference=FindCorrectReferenceIndexes(poltype,allindexes,pdbindexes)

    selectstringsidechainlist=[]
    for i in range(len(backboneindexessidechain)):
        index=backboneindexessidechain[i]
        selectstringsidechainlist.append('bynum '+str(index))
    
    selectstringreferencelist=[]
    for i in range(len(backboneindexesreference)):
        index=backboneindexesreference[i]
        selectstringreferencelist.append('bynum '+str(index))
    mergestring='not ('
    for i in range(len(backboneindexessidechain)):
        index=backboneindexessidechain[i]
        mergestring+='bynum'+' '+str(index)+' '+'or'+' '
    mergestring=mergestring[:-4]
    mergestring+=')'
    ref = Universe(refname)
    sc = Universe(pdbname)
    alignto(sc,ref,select={"mobile": selectstringsidechainlist,"reference": selectstringreferencelist})
    u = Merge(sc.select_atoms(mergestring),ref.select_atoms('all'))
    output=poltype.unmodifiedproteinpdbname.replace('.pdb','_Mutated.pdb')
    u.atoms.write(output) # this puts modified res on top of PDB but need to move to correct spot
    lines=GrabModifiedResidueLines(poltype,output,int(poltype.mutatedresiduenumber)) # MDAnalysis puts merged lines at top of pdb, need to move them and change chain back to original chain

    firstresnumber=GrabFirstResidueNumber(poltype,poltype.unmodifiedproteinpdbname) 
    chainid=GrabChainId(poltype,firstresnumber,output)

    lines=ModifyChainIdentifier(poltype,lines,chainid)
    # FIX ME, need to connect alpha carbon to side chain, PDB topology generated from matching SMARTS in library to this file, so need to connect here
    
    MoveModifiedLines(poltype,lines,int(poltype.mutatedresiduenumber),output,firstresnumber)
    newpdbfilename,connectedatomidx=ConnectSideChainToBackbone(poltype,output)
    GeneratePDBTopologyFile(poltype,newpdbfilename)
    

    return poltype.modifiedproteinpdbname,charge,connectedatomidx,backboneindexesreference


def ConnectSideChainToBackbone(poltype,pdbfilename):
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
   
    obConversion.ReadFile(pdbmol,pdbfilename)
    
    temp=open(pdbfilename,'r')
    results=temp.readlines()
    temp.close()
    atomindextoposition={}
    for line in results:
        if 'ATOM' in line and 'UNL' not in line: # sometimes UNL atoms from next amino acid (like an alpha carbon that does not continue on chain at end of file) needs to be removed
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            linesplit=line.split()
            position=[float(linesplit[6]),float(linesplit[7]),float(linesplit[8])]
            if resnumber==int(poltype.mutatedresiduenumber):
                if atomlabel=='CA':
                    refidx=atomindex
                    refpos=position
    bonditer=openbabel.OBMolBondIter(pdbmol)
    connectedidxs=[]
    for bond in bonditer:
        bgnatomidx=bond.GetBeginAtomIdx()
        endatomidx=bond.GetEndAtomIdx()
        if bgnatomidx==refidx:
            connectedidxs.append(endatomidx)
        elif endatomidx==refidx:
            connectedidxs.append(bgnatomidx)
    for line in results:
        if 'ATOM' in line and 'UNL' not in line: # sometimes UNL atoms from next amino acid (like an alpha carbon that does not continue on chain at end of file) needs to be removed
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            linesplit=line.split()
            position=[float(linesplit[6]),float(linesplit[7]),float(linesplit[8])]
            if resnumber==int(poltype.mutatedresiduenumber): # quick hack not general...
                if atomindex not in connectedidxs:
                    atomindextoposition[atomindex]=position
    atomindextodistance={}
    minatomindextodistance={}
    for atomindex,position in atomindextoposition.items():
        distance=numpy.linalg.norm(numpy.array(position)-numpy.array(refpos))
        atomindextodistance[atomindex]=distance
        if atomindex<10:
            minatomindextodistance[atomindex]=distance

    mindistance=min(minatomindextodistance.values())
    distancetoatomindex={v: k for k, v in atomindextodistance.items()}
    minatomindex=distancetoatomindex[mindistance]
    tempfilename=pdbfilename.replace('.pdb','_temp.pdb')
    newpdbfilename=pdbfilename.replace('.pdb','_mutsidechainconnected.pdb')
    obConversion.SetOutFormat('pdb')
    obConversion.WriteFile(pdbmol,tempfilename)
    newpdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(newpdbmol,tempfilename)
    

    searchpos=atomindextoposition[minatomindex]
    searchposref=atomindextoposition[refidx]
    temp=open(tempfilename,'r')
    results=temp.readlines()
    temp.close()
    atomindextoposition={}
    for line in results:
        if 'ATOM' in line and 'UNL' not in line: # sometimes UNL atoms from next amino acid (like an alpha carbon that does not continue on chain at end of file) needs to be removed
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            linesplit=line.split()
            position=[float(linesplit[6]),float(linesplit[7]),float(linesplit[8])]
            if position==searchpos:
                minatomindex=atomindex
            elif position==searchposref:
                refidx=atomindex

    newpdbmol.AddBond(refidx,minatomindex,1)
    obConversion.WriteFile(newpdbmol,newpdbfilename)

    return newpdbfilename,minatomindex
    
                
 

def GeneratePDBTopologyFile(poltype,filename):
    finalmol=openbabel.OBMol() # need to pass through babel again to change atom indexes
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(finalmol,filename)
    unmodidxtonewidx=OldIndexToNewIndexBabel(poltype,poltype.unmodifiedproteinpdbname,finalmol)
    # first add modified residue bond connection topology to library
    # first need to reindex the atoms then need to delete all bonds, then use our library residue_connections.txt to rebuild bond topology 
    lines=LibAddition(poltype,'MutatedSideChain.pdb',poltype.modifiedresiduepdbcode) 
    check=IsStringInFile(poltype,poltype.modifiedresiduepdbcode,poltype.topologylibpath)
    if check==False:
        AppendLinesToFile(poltype,lines,poltype.topologylibpath)
    finalmol=GeneratePDBTopologyObject(poltype,finalmol,filename)
    finalname=filename.replace('.pdb','_Final.pdb')
    obConversion.WriteFile(finalmol,finalname)
    poltype.modifiedproteinpdbname=finalname
    return poltype.modifiedproteinpdbname

def GeneratePDBTopologyObject(poltype,molecule,filename):
    molecule.BeginModify()
    molecule=RemoveAllBonds(poltype,molecule)
    smartstoatomorder,smartstobondtopo,smartstolabel,canonicallabels=ParseTopologyLib(poltype)
    molecule=BuildPDBTopology(poltype,molecule,smartstoatomorder,smartstobondtopo,filename)
    return molecule


def BuildPDBTopology(poltype,finalmol,smartstoatomorder,smartstobondtopo,filename):
    tempmol=openbabel.OBMol() # need to pass through babel again to change atom indexes
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(tempmol,filename)

    for smarts in smartstoatomorder.keys():
        atomorderlist=smartstoatomorder[smarts]
        bondtopo=smartstobondtopo[smarts]
        sp = openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,smarts)
        diditmatchOB=sp.Match(tempmol)

        if diditmatchOB==True:
            listofmatchlists=sp.GetMapList()
            residtomatch={}
            for match in listofmatchlists: # dont want to rebuild bonds more than once for same residue
                firstatomindex=match[0]
                firstatom=finalmol.GetAtom(firstatomindex)
                atomres=firstatom.GetResidue()
                atomresnum=atomres.GetNum()
                residtomatch[atomresnum]=match
            for resid in residtomatch.keys():
                match=residtomatch[resid]
                atomordertomatchedindex={}
                for i in range(len(match)):
                    matchedindex=match[i]
                    atomorder=atomorderlist[i]
                    atomordertomatchedindex[atomorder]=matchedindex
                for bondlist in bondtopo:
                    firstatomorder=bondlist[0]
                    secondatomorder=bondlist[1]
                    bondorder=bondlist[2]
                    firstatomindex=atomordertomatchedindex[firstatomorder]
                    secondatomindex=atomordertomatchedindex[secondatomorder]
                    finalmol.AddBond(firstatomindex,secondatomindex,bondorder)
    return finalmol
                    
                    
            
        
 

def AppendLinesToFile(poltype,lines,filename):
    temp=open(filename,'a')
    temp.write('\n')
    for line in lines:
        temp.write(line)
    temp.close()

def IsStringInFile(poltype,string,filename):
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    check=False
    for line in results:
        if string in line:
            check=True
    return check

def ParseTopologyLib(poltype):
    smartstoatomorder={}
    smartstobondtopo={}
    smartstolabel={}
    NAcanonical=['A', 'G', 'U', 'A3', 'A5','DA', 'DC', 'DG', 'DT', 'G3', 'G5', 'U3', 'U5','DA3', 'DA5', 'DC3', 'DC5', 'DG3', 'DG5', 'DT3', 'DT5']
    AAcanonicalcore=['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIS','ILE', 'LEU', 'LYD', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYD', 'TYR', 'VAL', 'CYD', 'CYS']
    AAcanonicalN=['N'+lab for lab in AAcanonicalcore]
    AAcanonicalC=["C"+lab for lab in AAcanonicalcore]
    canonicallabels=NAcanonical+AAcanonicalcore+AAcanonicalN+AAcanonicalC
    temp=open(poltype.topologylibpath,'r')
    results=temp.readlines()
    temp.close()
    passedconnect=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==0:
            if passedconnect==True:
                passedconnect=False
            continue
        if line[0]=='#' and 'CONNECT' not in line: # then this is comment with label
            label=linesplit[0].replace('#','')
            smarts=results[lineidx+1].replace('\n','')
            atomorderlist=[int(i) for i in results[lineidx+2].replace('\n','').split()]
            smartstoatomorder[smarts]=atomorderlist
            smartstolabel[smarts]=label
             
        elif line[0]=='#' and 'CONNECT' in line:
            passedconnect=True
        if len(linesplit)==3 and 'CHG' not in line and passedconnect==True:
            if smarts not in smartstobondtopo.keys():
                smartstobondtopo[smarts]=[]
            bondline=[int(i) for i in linesplit]
            smartstobondtopo[smarts].append(bondline)
    return smartstoatomorder,smartstobondtopo,smartstolabel,canonicallabels

def RemoveAllBonds(poltype,mol):
    bonditer=openbabel.OBMolBondIter(mol)
    bondlist=[]
    for bond in bonditer:
        bgnidx=bond.GetBeginAtomIdx()
        endidx=bond.GetEndAtomIdx()
        bgnatom=mol.GetAtom(bgnidx)
        endatom=mol.GetAtom(endidx)
        bgnatomres=bgnatom.GetResidue()
        bgnatomresnum=bgnatomres.GetNum()
        endatomres=endatom.GetResidue()
        endatomresnum=endatomres.GetNum()
        bgnatomlabel=GetAtomLabel(poltype,bgnatomres,bgnidx)
        endatomlabel=GetAtomLabel(poltype,endatomres,endidx)
        if bgnatomresnum!=endatomresnum:
            if (bgnatomlabel=='N' and endatomlabel=='C') or (bgnatomlabel=='C' and endatomlabel=='N'):
                continue
        bondlist.append(bond)
    for bond in bondlist:
        mol.DeleteBond(bond)
    return mol

def GetAtomLabel(poltype,residue,atomidx):
     for atom in openbabel.OBResidueAtomIter(residue):
         atName = residue.GetAtomID(atom)
         if atom.GetIdx()==atomidx:
             return atName.lstrip().rstrip()
 


def mol_with_atom_index(poltype,mol):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx()+1 ) )
    return mol

def ConvertSDFToMOL(poltype,filename):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('sdf')
    molbabel=openbabel.OBMol()
    tmpconv.ReadFile(molbabel,filename)
    tmpconv.SetOutFormat('mol')
    structfname=filename.replace('.sdf','.mol')
    tmpconv.WriteFile(molbabel,structfname)
    return structfname

def ConvertToSDF(poltype,filename):
    tmpconv = openbabel.OBConversion()
    inFormat = tmpconv.FormatFromExt(filename)
    tmpconv.SetInFormat(inFormat)
    molbabel=openbabel.OBMol()
    tmpconv.ReadFile(molbabel,filename)
    tmpconv.SetOutFormat('sdf')
    namesplit=filename.split('.')
    structfname=namesplit[0]+'.sdf'
    tmpconv.WriteFile(molbabel,structfname)
    return structfname



def GrabConnectRecord(poltype,molfile):
    if '.sdf' not in molfile:
        molfile=ConvertToSDF(poltype,molfile)

    temp=open(molfile,'r')
    results=temp.readlines()
    temp.close()
    lines=[]
    for line in results:
        linesplit=line.split()
        if len(linesplit)==7 and 'CHG' not in line:
            connectsplit=linesplit[:3] 
            connectline=' '.join(connectsplit)
            lines.append(connectline)
    return lines

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
    return ' '.join(atomorder)


def GrabAtomIndex(poltype,i,smirks):
    num=[]
    for j in range(i,len(smirks)):
        char=smirks[j]
        if char.isdigit():
            num.append(char)
        if char==']':
            break
    atomindex=''.join(num)
    return atomindex 

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

def MatchSMARTSrdkit(poltype,smarts,filename):
    mol=Chem.rdmolfiles.MolFromPDBFile(filename,removeHs=False,sanitize=False)
    p = Chem.MolFromSmarts(smarts)
    diditmatchrdkit=mol.HasSubstructMatch(p)
    return diditmatchrdkit

def ChiralFilter(poltype,smarts):
    return smarts.replace('@','')

def LibAddition(poltype,filename,label): 
    lines=[]
    mol=Chem.rdmolfiles.MolFromPDBFile(filename,removeHs=False)
    smarts=Chem.rdmolfiles.MolToSmarts(mol)
    smarts=ChiralFilter(poltype,smarts)
    diditmatchOB=MatchSMARTSOB(poltype,smarts,filename)
    diditmatchrdkit=MatchSMARTSrdkit(poltype,smarts,filename)
    if diditmatchrdkit==False:
        raise ValueError('rdkit SMARTS DOES NOT MATCH ',smarts)

    #if diditmatchOB==False:
    #     raise ValueError('OB SMARTS DOES NOT MATCH ',smarts)
        
    m=mol_with_atom_index(poltype,mol)
    smirks=Chem.rdmolfiles.MolToSmarts(m)
    atomorder=GrabAtomOrder(poltype,smirks)
    connect=GrabConnectRecord(poltype,filename)
    lines.append('#'+label+' '+filename.replace('.sdf','')+'\n')
    lines.append(smarts+'\n')
    lines.append(atomorder+'\n')
    lines.append('#CONNECT'+'\n')
    for line in connect:
        lines.append(line+'\n')
    lines.append('\n')
    return lines


def GrabFirstResidueNumber(poltype,pdbname):
    temp=open(pdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            resnumber=int(line[23:26+1])
            return resnumber


def GenerateModifiedProteinPoltypeInput(poltype):

    if poltype.unmodifiedproteinpdbname!=None:
        poltype.modifiedproteinpdbname,charge,connectedatomidx,backboneindexesreference=GenerateModifiedProteinPDB(poltype)
    knownresiduesymbs=GrabKnownResidueSymbs(poltype)

    # first need to generate protein XYZ file from PDB
    # so we have proidxtotypeidx, however since we will be transfering some of the protein parameters near boundary of poltype job and new residue in protein, we need to no which protein indexes correspond to type numbers in original protein parameter file vs what protein indexes will have new type numbers for the new residue

    obConversion = openbabel.OBConversion()
    proOBmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(poltype.modifiedproteinpdbname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(proOBmol,poltype.modifiedproteinpdbname) 
    totalatomnumber=proOBmol.NumAtoms()
    # have to regenerate PDB topology since when reading in the PDB, babel likes to guess new bonds by distance...
    proOBmol=GeneratePDBTopologyObject(poltype,proOBmol,poltype.modifiedproteinpdbname)
    
    modproidxs,modresnumber,modresiduelabel=GrabModifiedResidueProteinIndexes(poltype,poltype.modifiedproteinpdbname,knownresiduesymbs,totalatomnumber) # these indexes need to be excluded when calling pdbxyz.x then added back

    proboundidxs=FindBoundaryAtomIdxs(poltype,proOBmol,modproidxs)
    # now to reduce the number of atoms needed for QM and paramterization, just only include up to three atoms away starting from boundary atom on modified protein side


    modneighbidxs=NeighboringModResidueAtomIndexes(poltype,modresnumber,poltype.modifiedproteinpdbname,modproidxs,proOBmol,proboundidxs)
    # FIX ME
    # now need to grab residue number of modified residue, then use that to grab atom indexes of neighboring residues, then add those indexes to modified residue atom indexes and make fragment molecule with those atoms for poltype

    # so first step is to identify boundary atoms in a list iterate through all ligand atoms and if it has a neighbor in unmatched list, then it is a boundary atom
    choppedfragidxs=modproidxs+modneighbidxs
    chopmodmol=GenerateFragBabel(poltype,choppedfragidxs,poltype.modifiedproteinpdbname)
    # add hydrogens but only want to add to backbone, so remove any additional added hydrogens not on the backbone afterwords
    
    obConversion.SetOutFormat("pdb")
    refname='ModifiedRes.pdb'
    obConversion.WriteFile(chopmodmol,refname)
    newchopmodmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(refname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(newchopmodmol,refname)
    newchopmodmol.AddHydrogens() 
    proidxtoligidx=OldIndexToNewIndexBabel(poltype,poltype.modifiedproteinpdbname,newchopmodmol)
    modligidxs=[proidxtoligidx[j] for j in modproidxs]
    obConversion.SetOutFormat("pdb")
    refhydname='ModifiedResHydrated.pdb'
    obConversion.WriteFile(newchopmodmol,refhydname)
    positions,backboneidxs=GrabPositions(poltype,refname)
    # look for C or N on residue not same as modified residue, only keep additional hydrogens on those atoms
    idxstodelete=GrabAtomIndexesToDelete(poltype,refhydname,positions,newchopmodmol)
    modmol=openbabel.OBMol()
    obConversion.ReadFile(modmol,refhydname)
    modmol=DeleteAtomIndexes(poltype,idxstodelete,modmol)
    obConversion.SetOutFormat("mol")
    corename='ModifiedRes.mol'
    obConversion.WriteFile(modmol,corename) 
    time.sleep(2)
    m = Chem.MolFromMolFile(corename,removeHs=False)
    smarts=Chem.MolToSmarts(m)
    check=CheckIfSMARTSInLibrary(poltype,poltype.SMARTSToTypelibpath,smarts)
    boundaryatomidxs=[proidxtoligidx[i] for i in proboundidxs]
    backboneprobound=GrabCarbonylCarbonAndAmideNitrogen(poltype,modresnumber,poltype.modifiedproteinpdbname)
    backboneligbound=[proidxtoligidx[i] for i in backboneprobound]
    # okay now need to add hydrogens and convert to sdf
    polOBmol = openbabel.OBMol()
    # just take first atom and set total charge as formal charge for the first atom to make sure quantum has correct total charge
    inFormat = obConversion.FormatFromExt(corename)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(polOBmol,corename)

    polOBmol=AddTotalCharge(poltype,charge,polOBmol)
    molname=corename.replace('.mol','.sdf')
    obConversion.WriteFile(polOBmol,molname)
    
    polrdkitmol=Chem.rdmolfiles.MolFromPDBFile(refname)
    # first using the carbonyl carbon and backbone N of modprotein indexes indexes iterate over the neighbors and if the neighbor is not already in the dictionary then that must be the atom index that matches to the hydrogen indexes

    ligidxtoproidx = dict([v,k] for k,v in proidxtoligidx.items()) 
    proidxtoligidx,addedproatoms=MatchCorrectProteinAtomsToCorrespondingHydrogen(poltype,backboneprobound,proOBmol,polOBmol,proidxtoligidx,ligidxtoproidx)
    ligidxtoproidx = dict([v,k] for k,v in proidxtoligidx.items()) 
    return knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference,modligidxs

def AddTotalCharge(poltype,charge,molecule):
    atomiter=openbabel.OBMolAtomIter(molecule)
    for atom in atomiter:
        if atom.GetIdx()==1:
            atom.SetFormalCharge(charge)
    return molecule



def DeleteAtomIndexes(poltype,idxstodelete,molecule):
    sortedatomidxstoremove=sorted(idxstodelete, reverse=True)
    for index in sortedatomidxstoremove:
        atom=molecule.GetAtom(index)
        molecule.DeleteAtom(atom)
    return molecule

def GrabAtomIndexesToDelete(poltype,pdbfilename,positions,molecule):
    newpositions,backboneidxs=GrabPositions(poltype,pdbfilename)
    idxstodelete=[]
    temp=open(pdbfilename,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line and 'UNL' not in line: # sometimes UNL atoms from next amino acid (like an alpha carbon that does not continue on chain at end of file) needs to be removed
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            linesplit=line.split()
            position=tuple([float(linesplit[6]),float(linesplit[7]),float(linesplit[8])])
            check=IsAtomConnectedToBackboneAtom(poltype,position,backboneidxs,molecule)
  
            if check==False and position in newpositions and position not in positions:
                idxstodelete.append(atomindex)
    return idxstodelete

def IsAtomConnectedToBackboneAtom(poltype,position,backboneidxs,molecule):
    atomiter=openbabel.OBMolAtomIter(molecule)
    check=False
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if atomidx in backboneidxs:
            niter=openbabel.OBAtomAtomIter(atom)
            for natom in niter:
                pos=tuple([round(natom.GetX(),3),round(natom.GetY(),3),round(natom.GetZ(),3)])
                if pos==position:
                    check=True
    return check
    


def GrabPositions(poltype,pdbfilename):
    positions=[]
    backboneidxs=[]
    temp=open(pdbfilename,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line and 'UNL' not in line: # sometimes UNL atoms from next amino acid (like an alpha carbon that does not continue on chain at end of file) needs to be removed
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            resnumber=int(line[23:26+1])
            linesplit=line.split()
            position=tuple([float(linesplit[6]),float(linesplit[7]),float(linesplit[8])])
            if resnumber!=int(poltype.mutatedresiduenumber) and (atomlabel=='C' or atomlabel=='N'):
                backboneidxs.append(atomindex)
            positions.append(position)

    return positions,backboneidxs



def GrabLibraryInfo(poltype,proidxtoprotype,modresiduelabel,proOBmol):
    modresiduedic={}
    temp=open(poltype.modifiedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    etab = openbabel.OBElementTable()
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())
            residuelabel=line[17:20].rstrip().lstrip()
            atomtype=proidxtoprotype[atomindex]
            atomclass=atomtype # poltype assigns same class and type numbers
            proatom=proOBmol.GetAtom(atomindex)
            atomicnum=proatom.GetAtomicNum()
            atomicmass=etab.GetMass(atomicnum)
            if residuelabel==modresiduelabel:
               if modresiduelabel not in modresiduedic.keys():
                   modresiduedic[modresiduelabel]={}
                   modresiduedic[modresiduelabel]['atomindex']={}
                   modresiduedic[modresiduelabel]['atomtype']={}
                   modresiduedic[modresiduelabel]['atomclass']={}
                   modresiduedic[modresiduelabel]['atomlabel']={}
                   modresiduedic[modresiduelabel]['atomicnum']={}
                   modresiduedic[modresiduelabel]['atomicmass']={}
               if atomindex not in modresiduedic[modresiduelabel]['atomindex']: # then assume that this atom is not in any of the lists
                   modresiduedic[modresiduelabel]['atomindex'][atomindex]=atomindex
                   modresiduedic[modresiduelabel]['atomtype'][atomindex]=atomtype
                   modresiduedic[modresiduelabel]['atomclass'][atomindex]=atomclass
                   modresiduedic[modresiduelabel]['atomlabel'][atomindex]=atomlabel
                   modresiduedic[modresiduelabel]['atomicnum'][atomindex]=atomicnum
                   modresiduedic[modresiduelabel]['atomicmass'][atomindex]=atomicmass

    return modresiduedic

def GenerateLibFileToAdd(poltype,modresiduedic,modresiduelabel):
    tempname='ModifiedResidueLibAddition.txt'
    temp=open(tempname,'w')
    temp2=open(poltype.libpath,'a')
    temp2.write('\n')
    temp.write('# '+'RESIDUE '+modresiduelabel+'\n')
    temp2.write('# '+'RESIDUE '+modresiduelabel+'\n')
    dic=modresiduedic[modresiduelabel]
    for atomindex in dic['atomindex']:
        atomtype=str(dic['atomtype'][atomindex])
        atomclass=str(dic['atomclass'][atomindex])
        atomlabel=str(dic['atomlabel'][atomindex])
        atomicnum=str(dic['atomicnum'][atomindex])
        atomicmass=str(dic['atomicmass'][atomindex])
        atomicmassdigits=atomicmass.split('.')[0]
        if len(atomtype)==3:
            firstspace='  '
        elif len(atomtype)==2:
            firstspace='   '
        elif len(atomtype)==1:
            firstspace='    '
        if len(atomicmassdigits)==2:
            massspace=' '
        elif len(atomicmassdigits)==1:
            massspace='  '
        if len(atomlabel)==1:
            resspace='     '
        elif len(atomlabel)==2:
            resspace='    '
        elif len(atomlabel)==3:
            resspace='   '
        elif len(atomlabel)==4:
            resspace='  '
        temp.write('atom'+firstspace+atomtype+firstspace+atomclass+' '+atomlabel+resspace+modresiduelabel+'                    '+atomicnum+massspace+atomicmass+'  '+'0'+'\n')
        temp2.write('atom'+firstspace+atomtype+firstspace+atomclass+' '+atomlabel+resspace+modresiduelabel+'                    '+atomicnum+massspace+atomicmass+'  '+'0'+'\n')
    temp.close()
    temp2.close()
           


def GrabProSideBoundIdxs(poltype,proOBmol,boundaryatomidxs,ligidxtoproidx,modproidxs,proidxtoligidx):
    prosideboundligidx=[]
    for index in boundaryatomidxs:
        proindex=ligidxtoproidx[index]
        proatom=proOBmol.GetAtom(proindex)
        proneighbs=[natom for natom in openbabel.OBAtomAtomIter(proatom)]
        proneighbidxs=[i.GetIdx() for i in proneighbs]
        for proneighbidx in proneighbidxs:
            if proneighbidx not in modproidxs: # this means that you are boundary atom on the protein side, now convert back to ligand index
                ligidx=proidxtoligidx[proneighbidx]
                prosideboundligidx.append(ligidx)
    return prosideboundligidx


def GrabPolarizeNeighbors(poltype,ligtype,keyfile):
    temp=open(keyfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'polarize' in line:
            linesplit=line.split()
            poltype=int(linesplit[1])
            if ligtype==poltype:
                stringneighbs=linesplit[4:]
                ligneighbs=[int(i) for i in stringneighbs]
    return ligneighbs

def GrabResidueTypesToPolarizeTypesAcrossBoundary(poltype,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,ligidxtoproidx,key,mincorenumber,maxcorenumber):    
    prmtypetopoltype={} # this dictionary has keys as prm type numbers that need poltype type numbers appended to their definition in the new key file to handle boundaries
    for ligidx in boundaryatomidxs:
        ligtype=ligidxtotypeidx[ligidx]
        # now we need to grab the polarization neighbors for this from the first key file ( need to convert indexes to protein indexes)
        ligneighbs=GrabPolarizeNeighbors(poltype,ligtype,key)
        ligidxs=GrabIdxsFromType(poltype,ligneighbs,ligidxtotypeidx)
        for ligidx in ligidxs:
            proidx=ligidxtoproidx[ligidx]
            protype=proidxtoprotype[proidx]
            if protype<mincorenumber or protype>maxcorenumber: # then this is a type number from the prm file and we need to append ligtype to the definition of protype that is in the prm file and add that to the key file
                prmtypetopoltype[protype]=ligtype
    return prmtypetopoltype

def GrabPolarizeDefinitionLinesFromPrmFile(poltype,prmtypetopoltype,proprmfile):
    prmtypetoprmdeflines={}
    temp=open(proprmfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'polarize' in line:
            linesplit=line.split()
            prmtype=int(linesplit[1])
            if prmtype in prmtypetopoltype.keys():
                prmtypetoprmdeflines[prmtype]=line.replace('\n','')
    return prmtypetoprmdeflines

def AppendPoltypeTypeNumbers(poltype,prmtypetoprmdeflines,prmtypetopoltype):
    newprmdeflines=[]
    for prmtype in prmtypetoprmdeflines.keys():
        prmdefline=prmtypetoprmdeflines[prmtype]
        poltype=prmtypetopoltype[prmtype]
        newprmdefline=prmdefline+' '+str(poltype)+'\n'
        newprmdeflines.append(newprmdefline)
    return newprmdeflines


def GrabClassNumbersForProtein(poltype,proidxtoprotype,check):
    prmtypetoprmclass={}
    if check==False:
        temp=open(poltype.amoebabioprmpath,'r')
    else:
        temp=open(poltype.ModifiedResiduePrmPath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'atom' in line:
            typenum=int(linesplit[1])
            classnum=int(linesplit[2])
            if typenum in proidxtoprotype.values():
                prmtypetoprmclass[typenum]=classnum
        



    return prmtypetoprmclass

def GrabProBoundIdxToProType(poltype,proboundidxs,proOBmol,modproidxs,proidxtotypeidx,transferableidxs,transferableidxsformatching,firstresnum,lastresnum):
    
    proboundidxtoprotype={}
    rdkittransferableidxs=[i-1 for i in transferableidxsformatching]
    
    rdkitmol=Chem.rdmolfiles.MolFromPDBFile(poltype.modifiedproteinpdbname,removeHs=False,sanitize=False)
    mol=GenerateFragRdkit(poltype,rdkittransferableidxs,rdkitmol,poltype.modifiedproteinpdbname)
    oldindextonewindex=OldIndexToNewIndexRdkit(poltype,poltype.modifiedproteinpdbname,mol)
    newindextooldindex={v: k for k, v in oldindextonewindex.items()}
    m=mol_with_atom_index(poltype,mol)
    fragsmirks=Chem.rdmolfiles.MolToSmarts(m)

    fragsmarts=Chem.rdmolfiles.MolToSmarts(mol)
    p = Chem.MolFromSmarts(fragsmarts)
    matches=mol.GetSubstructMatches(p) 
    firstmatchmol=matches[0]
    matches=rdkitmol.GetSubstructMatches(p) # just grab a match from somewhere else in protein without modproidxs
    for match in matches:
        allnotmodandtermres=True
        for idx in match:
            babelidx=idx+1
            atom=proOBmol.GetAtom(babelidx)
            atomres=atom.GetResidue()
            atomresnum=atomres.GetNum()
            if atomresnum==int(firstresnum) or atomresnum==int(lastresnum):
                allnotmodandtermres=False
            if babelidx in modproidxs:
                allnotmodandtermres=False
        if allnotmodandtermres==True: # then just use these to get type numbers, need to make sure you do not match C or N terminal!! only middle
            firstmatch=match
    pdbmatchidxtofragmolmatchidx=dict(zip(firstmatch, firstmatchmol))
    firstmatchfragmol=[pdbmatchidxtofragmolmatchidx[i] for i in firstmatch]
    oldboundidxsrdkit=[newindextooldindex[i] for i in firstmatchfragmol] 
    oldboundidxsbabel=[i+1 for i in oldboundidxsrdkit]
    for i in range(len(firstmatch)):
        newproidxbabel=firstmatch[i]+1
        newprotype=proidxtotypeidx[newproidxbabel]
        oldidx=oldboundidxsbabel[i]
        proboundidxtoprotype[oldidx]=newprotype
    for i in range(len(transferableidxs)):
        boundidx=transferableidxs[i]
        if boundidx not in proboundidxtoprotype.keys():
            protype=proidxtotypeidx[boundidx]
            proboundidxtoprotype[boundidx]=protype
    return proboundidxtoprotype
                 


def ReadSMARTSToTypeLib(poltype,smarts,modmol):
    ligidxtotypeidx={}
    temp=open(poltype.SMARTSToTypelibpath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>1:
            smartstring=linesplit[0]
            atompos=int(linesplit[1])
            typeidx=int(linesplit[2])
            if smarts==smartstring:
                grabmatch=MatchSMARTSToOBMol(poltype,smarts,modmol)
                ligidx=grabmatch[atompos-1]
                ligidxtotypeidx[ligidx]=typeidx
    return ligidxtotypeidx

       
def GenerateModifiedProteinXYZAndKey(poltype,knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference,modligidxs):

    obConversion = openbabel.OBConversion()
    ligOBmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(molname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(ligOBmol,molname) 
    ligandindexes=set(range(1,ligOBmol.NumAtoms()))

    
    if check==False:
        ligidxtotypeidx=GenIndexToTypeIndexDic(poltype,poltype.tmpxyzfile)
    else:
        ligidxtotypeidx=ReadSMARTSToTypeLib(poltype,smarts,modmol)
    modproidxtotypenumber=AssignNewTypeNumbers(poltype,modproidxs,ligidxtotypeidx,proidxtoligidx) # will try to keep modified residue type numbers as type numbers in poltype job but then for boudnary parts need to assign new type numbers 
    proidxtoprotype,firstresnum,lastresnum=GrabProteinTypeNumbers(poltype,poltype.modifiedproteinpdbname,knownresiduesymbs,poltype.libpath,modproidxtotypenumber)
    protinkxyz,firstresnum,lastresnum=GenerateProteinTinkerXYZFile(poltype,poltype.modifiedproteinpdbname,modproidxtotypenumber,poltype.amoebabioprmpath,proidxtoprotype,knownresiduesymbs)

    listoftorsionsforkey,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForValenceTermsAcrossBoundary(poltype,ligOBmol,proOBmol,ligidxtoproidx,proboundidxs,boundaryatomidxs,modproidxs) # execlude all cases where only 1 atom is in modproidxs and the others are in protein, seperate those into different lists since those can be transfered from prm file and the others from the key file

    # for stretch bend (strbnd) we can use the angle types to grab the strbnd parameters

    # for  polarize just use list of single atom types, (the list of boundary atom types), do not need to transfer vdw, that does not cross boundary polarize neighbors do cross the boundary



    transferableidxs=GetIdxsFromListOfTorsionsForPrm(poltype,listoftorsionsforprm) # now take connected atomidx (usually CB) and add any hydrogens connected
    transferableidxsformatching=AddHydrogensAndBackBoneOnConnectedAtomIdxToGetRightClass(poltype,connectedatomidx,proOBmol,backboneindexesreference)
    proboundidxtoprotype=GrabProBoundIdxToProType(poltype,proboundidxs,proOBmol,modproidxs,proidxtoprotype,transferableidxs,transferableidxsformatching,firstresnum,lastresnum)
     
    # now convert protein index numbers to ligand index, then we will convert ligand index to type numbers to grab parameters. We also need dictionaries that map each type in ligand type numbers to what will be represented in final key file (so for a torsion across a boundary we want the ligand type numbers to grab parmeters but we want to have protein type numbers for part of torsion on the protein side in our key file)

    # these parameters will be taken from the key file so need to no what to convert the indexes to once they are taken from the key file
    if check==False:
        maxcorenumber=db.GrabMaxTypeNumber(poltype,poltype.key5fname)
        mincorenumber=db.GrabMinTypeNumber(poltype,poltype.key5fname)
        key=poltype.key5fname
        if not os.path.isfile(poltype.ModifiedResiduePrmPath):
            shutil.copy(poltype.amoebabioprmpath,poltype.ModifiedResiduePrmPath)
        if not os.path.isfile(poltype.SMARTSToTypelibpath):
            temp=open(poltype.SMARTSToTypelibpath,'a')
            temp.write('\n')
            temp.close()

    else:
        maxcorenumber=db.GrabMaxTypeNumberModifiedCore(poltype,poltype.ModifiedResiduePrmPath)
        mincorenumber=db.GrabMinTypeNumberModifiedCore(poltype,poltype.ModifiedResiduePrmPath)
        key=poltype.ModifiedResiduePrmPath

    writekey=poltype.key5fname.replace('.key_5','.key_6')
    torsionligtypestoboundtypesforkey=GrabLigandTypeToBoundryTypeForKeyFile(poltype,listoftorsionsforkey,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,mincorenumber,maxcorenumber)
    prmtypetoprmclass=GrabClassNumbersForProtein(poltype,proidxtoprotype,check) # need this for grepping torsion parameters from prm file
    bondprmclassestoboundclasses=GrabPrmClassToBoundryClassForPrmFile(poltype,listofbondsforprm,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber)
    
    angleprmclassestoboundclasses=GrabPrmClassToBoundryClassForPrmFile(poltype,listofanglesforprm,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber)
    torsionprmclassestoboundclasses=GrabPrmClassToBoundryClassForPrmFile(poltype,listoftorsionsforprm,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber)
    ligboundaryatomtypes=[ligidxtotypeidx[i] for i in boundaryatomidxs]
    # for the polarize and multipole we need to not only have list of atomindexes (ligand indexes) that are boundary atoms on the ligand side, but also the protein side (mulitpole and polarize use neighbors that may go across the boundaries)
       
    prosideboundligidx=GrabProSideBoundIdxs(poltype,proOBmol,boundaryatomidxs,ligidxtoproidx,modproidxs,proidxtoligidx)
    prmtypetopoltype=GrabResidueTypesToPolarizeTypesAcrossBoundary(poltype,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,ligidxtoproidx,key,mincorenumber,maxcorenumber)
    prmtypetoprmdeflines=GrabPolarizeDefinitionLinesFromPrmFile(poltype,prmtypetopoltype,poltype.amoebabioprmpath)
    prmdeflines=AppendPoltypeTypeNumbers(poltype,prmtypetoprmdeflines,prmtypetopoltype) # now also grab them from the key file and then append the prm types then combine with prmdeflines here
    poltypetoprmtype = dict([v,k] for k,v in prmtypetopoltype.items())
    
    atomtypetoframedef,proteintypestoframedefforprmfile=GrabMultipoleFrameDefintions(poltype,key,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,modproidxs,ligidxtoproidx,prosideboundligidx,proboundidxtoprotype,mincorenumber,maxcorenumber,modligidxs) # for multipole frames on the boundary (on ligand side), we need to grab the original index frame definitions then convert to protein type number when appropriate, also need to consider atoms on other side of boundary (protein side)

    # need to add charge
    stitchtorsionprms,stitchpolarizeprms,stitchmpoleprms=GrabParametersFromKeyFile(poltype,key,torsionligtypestoboundtypesforkey,ligboundaryatomtypes,poltypetoprmtype,atomtypetoframedef,mincorenumber,maxcorenumber)
    if check==False:
        fname=poltype.amoebabioprmpath
    else:
        fname=poltype.ModifiedResiduePrmPath
    polclasstoprmclass={}
    stitchbondprms,stitchangleprms,stitchnewtorsionprms,stitchstrbndprms,stitchnewmpoleprms,stitchopbendprms,uneededpolarizeprms,unneededvdwprms=db.GrabParametersFromPrmFile(poltype,bondprmclassestoboundclasses,angleprmclassestoboundclasses,torsionprmclassestoboundclasses,poltypetoprmtype,polclasstoprmclass,proteintypestoframedefforprmfile,fname)
    stitchtorsionprms=stitchtorsionprms+stitchnewtorsionprms
    stitchmpoleprms=stitchnewmpoleprms+stitchmpoleprms
    stitchpolarizeprms=prmdeflines+stitchpolarizeprms
    stitchatomdefs=[]
    stitchvdwprms=[]
    WriteNewKeyFile(poltype,stitchatomdefs,stitchbondprms,stitchangleprms,stitchtorsionprms,stitchstrbndprms,stitchopbendprms,stitchpolarizeprms,stitchvdwprms,stitchmpoleprms,writekey,poltype.amoebabioprmpath)
    if check==False:
        coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms=db.GrabParameters(poltype,poltype.key5fname)
        oldtypetonewtype,shift=db.ShiftPoltypeNumbers(poltype,poltype.ModifiedResiduePrmPath,poltype.key5fname)
        mincorenumber=mincorenumber-shift
        maxcorenumber=maxcorenumber-shift 
        # now I need to convert all arrays
        coreresult=db.ShiftParameterDefintions(poltype,[coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms],oldtypetonewtype)
        coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms=coreresult[:] 
        stitchresult=db.ShiftParameterDefintions(poltype,[stitchatomdefs,stitchbondprms,stitchangleprms,stitchtorsionprms,stitchstrbndprms,stitchopbendprms,stitchpolarizeprms,stitchvdwprms,stitchmpoleprms],oldtypetonewtype)
        stitchatomdefs,stitchbondprms,stitchangleprms,stitchtorsionprms,stitchstrbndprms,stitchopbendprms,stitchpolarizeprms,stitchvdwprms,stitchmpoleprms=stitchresult[:]
        writekey=poltype.key5fname.replace('.key_5','.key_7')
        WriteNewKeyFile(poltype,coreatomdefs+stitchatomdefs,corebondprms+stitchbondprms,coreangleprms+stitchangleprms,coretorsionprms+stitchtorsionprms,corestrbndprms+stitchstrbndprms,coreopbendprms+stitchopbendprms,corepolarizeprms+stitchpolarizeprms,corevdwprms+stitchvdwprms,corempoleprms+stitchmpoleprms,writekey,poltype.amoebabioprmpath)       
        db.WriteToPrmFile(poltype,coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms,poltype.ModifiedResiduePrmPath)
        proidxtoprotype=ShiftDictionaryValueTypes(poltype,proidxtoprotype,oldtypetonewtype)
        ligidxtotypeidx=ShiftDictionaryValueTypes(poltype,ligidxtotypeidx,oldtypetonewtype)
        modproidxtotypenumber=ShiftDictionaryValueTypes(poltype,modproidxtotypenumber,oldtypetonewtype)
        # now I need to make a text file that takes parameters and can add to library
        modresiduedic=GrabLibraryInfo(poltype,proidxtoprotype,modresiduelabel,proOBmol)
        protinkxyz,firstresnum,lastresnum=GenerateProteinTinkerXYZFile(poltype,poltype.modifiedproteinpdbname,modproidxtotypenumber,poltype.amoebabioprmpath,proidxtoprotype,knownresiduesymbs)

        cmdstr=CallAnalyze(poltype,protinkxyz,writekey,'e','alz.out')
        poltype.call_subsystem(cmdstr,True)
        error=ReadAnalyzeOutput(poltype)
        if error==False:
            GenerateLibFileToAdd(poltype,modresiduedic,modresiduelabel)
            GenerateSMARTSToTypeFileToAdd(poltype,modmol,ligidxtotypeidx,smarts)
        else:
            raise ValueError(' alz.out for PDB XYZ file and .key_7 have errors')
        CorrectTotalCharge(poltype,poltype.unmodifiedproteinpdbname,writekey,protinkxyz,mincorenumber,maxcorenumber)

def AddHydrogensAndBackBoneOnConnectedAtomIdxToGetRightClass(poltype,connectedatomidx,proOBmol,backboneindexesreference):
    transferableidxsformatching=[connectedatomidx]
    transferableidxsformatching.extend(backboneindexesreference)
    hydidxs=[]
    for atomidx in transferableidxsformatching:
        atom=proOBmol.GetAtom(atomidx)
        atomatomiter=openbabel.OBAtomAtomIter(atom)
        for atom in atomatomiter:
            if atom.GetAtomicNum()==1:
                atomidx=atom.GetIdx()
                if atomidx not in hydidxs:
                    hydidxs.append(atomidx)
    transferableidxsformatching.extend(hydidxs)
    return transferableidxsformatching
   

def GetIdxsFromListOfTorsionsForPrm(poltype,listoftorsionsforprm):
    transferableidxs=[]
    for torsion in listoftorsionsforprm:
        for idx in torsion:
            if idx not in transferableidxs:
                transferableidxs.append(idx)
    return transferableidxs


def CallAnalyze(poltype,xyzfile,keyfile,option,outfile):
    cmdstr=poltype.analyzeexe+' '+xyzfile+' '+'-k'+' '+keyfile+' '+option+' > '+outfile
    return cmdstr

def ReadAnalyzeOutput(poltype):
    error=False
    temp=open('alz.out','r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Tinker is Unable to Continue' in line or 'error' in line or 'Error' in line:
            error=True
    return error

def ShiftDictionaryValueTypes(poltype,dictionary,oldtypetonewtype):
    for key,value in dictionary.items():
        if value in oldtypetonewtype.keys():
            newvalue=oldtypetonewtype[value]
            dictionary[key]=newvalue
    return dictionary



                   
    

def GrabMaxTypeNumberModifiedCore(poltype,parameterfile):
    maxnumberfromprm=1
    temp=open(parameterfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line and 'Modified' in line:
            linesplit=line.split()
            atomtype=int(linesplit[1])
            if atomtype>maxnumberfromprm:
                maxnumberfromprm=atomtype
    return maxnumberfromprm

def GrabMinTypeNumberModifiedCore(poltype,parameterfile):
    minnumberfromprm=10000
    temp=open(parameterfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line and 'Modified' in line:
            linesplit=line.split()
            atomtype=int(linesplit[1])
            if atomtype<minnumberfromprm:
                minnumberfromprm=atomtype
    return minnumberfromprm



def GenerateSMARTSToTypeFileToAdd(poltype,modmol,ligidxtotypeidx,smarts):
    libname='SMARTSToType.txt'
    grabmatch=MatchSMARTSToOBMol(poltype,smarts,modmol)
    atompostotype={}
    for atomposidx in range(len(grabmatch)):
        ligidx=grabmatch[atomposidx]
        ligtypeidx=ligidxtotypeidx[ligidx]
        atompostotype[atomposidx]=ligtypeidx
    temp=open(libname,'w')
    temp2=open(poltype.SMARTSToTypelibpath,'a')
    for atomposidx in atompostotype.keys():
        atomtype=atompostotype[atomposidx]
        temp.write(smarts+' '+str(atomposidx+1)+' '+ str(atomtype)+'\n')
        temp2.write(smarts+' '+str(atomposidx+1)+' '+ str(atomtype)+'\n')
    temp2.write('\n')
    temp.close()
    temp2.close()

def MatchSMARTSToOBMol(poltype,smarts,modmol):
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    diditmatch=sp.Match(modmol)
    listofmatchlists=sp.GetMapList()
    for matchlist in listofmatchlists:
        grabmatch=matchlist
    if diditmatch==True:
        return grabmatch
    else:
        print('Houston, we have a problem')
        return None

def MatchSMARTSToPDB(poltype,smarts,modmol):
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    diditmatch=sp.Match(modmol)
    listofmatchlists=sp.GetMapList()
    for matchlist in listofmatchlists:
        firstatomidx=matchlist[0]
        firstatom=modmol.GetAtom(firstatomidx)
        atomres=firstatom.GetResidue()
        atomresnum=atomres.GetNum()
        if atomresnum==int(poltype.mutatedresiduenumber):
            grabmatch=matchlist
            break
    if diditmatch==True:
        return grabmatch
    else:
        print('Houston, we have a problem')
        return None


def MatchSMARTSToOBMolGrabAllMatches(poltype,smarts,modmol):
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    diditmatch=sp.Match(modmol)
    listofmatchlists=sp.GetMapList()
    matches=[]
    for matchlist in listofmatchlists:
        grabmatch=matchlist
        matches.append(grabmatch)
    if diditmatch==True:
        return matches
    else:
        print('Houston, we have a problem')
        return None 

def BabelGenerateMol2File(poltype,pdbfile):
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('pdb')
    molbabel=openbabel.OBMol()
    tmpconv.ReadFile(molbabel,pdbfile)
    tmpconv.SetOutFormat('mol2')
    structfname=pdbfile.replace('.pdb','.mol2')
    tmpconv.WriteFile(molbabel,structfname)
    return structfname 


def ReadTotalProteinChargeBeforeMerge(poltype,mol2file):
    temp=open(mol2file,'r')
    results=temp.readlines()
    temp.close()
    chg=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)==9:
            chg+=float(linesplit[-1])
    totalchg=round(chg)
    return totalchg

def ReadTotalProteinChargeAfterMerge(poltype,alzfile):
    temp=open(alzfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Total Electric Charge : ' in line:
            linesplit=line.split()
            chg=float(linesplit[-2])
    return chg

def ChargeToRedistribute(poltype,proteinchargebeforemerge,ligandcharge,proteinchargeaftermerge):
    return proteinchargeaftermerge-(proteinchargebeforemerge+ligandcharge)

def ReadMonopolesFromMultipolesAcrossBoundary(poltype,keyfile,mincorenumber,maxcorenumber):
    temp=open(keyfile,'r')
    results=temp.readlines()
    temp.close()
    monopoles=[]
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            monopole=float(linesplit[-1])
            frames=linesplit[1:-1]
            frames=[int(i) for i in frames]
            polnums=CountNumberPoltypeNums(poltype,frames,mincorenumber,maxcorenumber)
            if polnums!=len(frames):
                monopoles.append(monopole)
    return monopoles

def ModifyMonopoles(poltype,monopoles,chgtoredistribute):
    nummonopoles=len(monopoles)
    chgpermonopole=chgtoredistribute/nummonopoles 
    modmonopoles=[]
    for monopole in monopoles:
        if monopole>0 and chgpermonopole>0:
            modmonopole=monopole-chgpermonopole
        elif monopole<0 and chgpermonopole>0:
            modmonopole=monopole-chgpermonopole
        elif monopole>0 and chgpermonopole<0:
            modmonopole=monopole+chgpermonopole
        elif monopole<0 and chgpermonopole<0:
            modmonopole=monopole+chgpermonopole
        modmonopoles.append(modmonopole)
    return modmonopoles


def ModifyMonopolesFromMultipolesAcrossBoundary(poltype,keyfile,modmonopoles,mincorenumber,maxcorenumber):
    temp=open(keyfile,'r')
    results=temp.readlines()
    temp.close()
    tempname='temp.txt'
    temp=open(tempname,'w')
    modcount=0
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            frames=linesplit[1:-1]
            frames=[int(i) for i in frames]
            polnums=CountNumberPoltypeNums(poltype,frames,mincorenumber,maxcorenumber)
            if polnums!=len(frames):
                linesplit=re.split(r'(\s+)', line)
                linesplit[-3]=str(modmonopoles[modcount])
                modcount+=1
                newline=''.join(linesplit)
                temp.write(newline)
        else:
            temp.write(line)
    temp.close()
    os.remove(keyfile)
    os.rename(tempname,keyfile)

def CorrectTotalCharge(poltype,pdbfile,writekey,protinkxyz,mincorenumber,maxcorenumber):
    mol2file=BabelGenerateMol2File(poltype,pdbfile)
    proteinchargebeforemerge=ReadTotalProteinChargeBeforeMerge(poltype,mol2file)
    alzfile='alz_charge.out'
    cmdstr=CallAnalyze(poltype,protinkxyz,writekey,'m',alzfile)
    poltype.call_subsystem(cmdstr,True)
    proteinchargeaftermerge=ReadTotalProteinChargeAfterMerge(poltype,alzfile)
    ligandcharge=poltype.totalcharge
    chgtoredistribute=ChargeToRedistribute(poltype,proteinchargebeforemerge,ligandcharge,proteinchargeaftermerge)
    monopoles=ReadMonopolesFromMultipolesAcrossBoundary(poltype,writekey,mincorenumber,maxcorenumber)
    modmonopoles=ModifyMonopoles(poltype,monopoles,chgtoredistribute)
    totnewchg=sum(modmonopoles)
    
    ModifyMonopolesFromMultipolesAcrossBoundary(poltype,writekey,modmonopoles,mincorenumber,maxcorenumber)
