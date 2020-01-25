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

def GrabModifiedResidueProteinIndexes(poltype,modifiedproteinpdbname,knownresiduesymbs):
    temp=open(modifiedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    modproidxs=[]
    
    for line in results:
        if 'ATOM' in line:
            linesplit=line.split()
            PDBcode=line[17:19+1]
            resnumber=line[23:26+1]
            if PDBcode not in knownresiduesymbs:
                proatomidx=int(linesplit[1])
                modproidxs.append(proatomidx)
                modresnumber=int(resnumber)
                modresiduelabel=PDBcode
    return modproidxs,modresnumber,modresiduelabel


def GrabProteinTypeNumbers(poltype,pdbfilename,knownresiduesymbs,libpath,modproidxtotypenumber):
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


    proidxtotypeidx={}
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
        if protypeidx==0 and proidx in modproidxtotypenumber.keys():
            protypeidx=modproidxtotypenumber[proidx]
     
            
        proidxtotypeidx[proidx]=int(protypeidx)


    return proidxtotypeidx


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
    
    # first take original PDB, convert to mol2 to preserve the atom label and residue names, then need to convert back to PDB to generate any connectivity via open babel that does not exist
    mol2fname=modifiedproteinpdbname.replace('.pdb','.mol2')
    cmdstr='babel'+' '+'-ipdb'+' '+modifiedproteinpdbname+' '+'-omol2'+' '+'-O'+' '+mol2fname
    print('Calling: '+cmdstr)
    os.system(cmdstr)
    newpdbname=mol2fname.replace('.mol2','_conect.pdb')
    cmdstr='babel'+' '+'-imol2'+' '+mol2fname+' '+'-opdb'+' '+'-O'+' '+newpdbname
    print('Calling: '+cmdstr)
    os.system(cmdstr)
    # now read pdb into openbabel structure and then also generate idxtoatomlabel dictionary
    obConversion = openbabel.OBConversion()
    molstruct = openbabel.OBMol()
    obConversion.SetInFormat('.pdb')
    obConversion.ReadFile(molstruct, newpdbname) 
    # now create idxtoatomlabel
    idxtoatomlabel=GenerateIdxToAtomLabel(poltype,newpdbname)
    structfname=newpdbname.replace('.pdb','.xyz')
    GenerateXYZFile(poltype,molstruct, structfname,idxtoatomlabel)
    proidxtoprotype=GrabProteinTypeNumbers(poltype,modifiedproteinpdbname,knownresiduesymbs,poltype.libpath,modproidxtotypenumber)
    temp=open(structfname,'r')
    results=temp.readlines()
    temp.close()
    protinkxyz=structfname.replace('.xyz','_Final.xyz')
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
                    onebool=True
                else:
                    ligidx=proidxtoligidx[proidx]
                    boundtype=ligidxtotypeidx[ligidx]
                        
                if boundtype<mincorenumber or boundtype>maxcorenumber:
                    boundclass=GrabTinkerClassNumber(poltype,boundtype,prmfile)
                    boundtypes.append(boundclass)
                else:
                    boundtypes.append(boundtype)
            valligtypestoboundtypes[valtup]=boundtypes
            
    return valligtypestoboundtypes


def GrabPrmClassToBoundryClassForPrmFile(poltype,listofvalterms,proidxtoligidx,ligidxtotypeidx,proidxtotypeidx,prmfile,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber): # be careful because need to convert type number to class number for valence only if it is a protein type number (poltype type numbers are the same as class numbers) terms!
    prmclassestoboundclasses={}
    for val in listofvalterms:
        valtype=[]
        boundclasses=[]
        for proidx in val:
            protype=proidxtotypeidx[proidx]
            if protype<mincorenumber or protype>maxcorenumber:
                valtype.append(protype)
                boundclasses.append(prmtypetoprmclass[protype])
            else:
                valtype.append(proboundidxtoprotype[proidx])
                boundclasses.append(protype)
        valclasses=[prmtypetoprmclass[i] for i in valtype]
        if tuple(valclasses) in prmclassestoboundclasses.keys() and tuple(valclasses[::-1]) not in prmclassestoboundclasses.keys(): # then we need to make a list
            prmclassestoboundclasses[tuple(valclasses)].append(boundclasses)
        elif tuple(valclasses[::-1]) in prmclassestoboundclasses.keys() and tuple(valclasses) not in prmclassestoboundclasses.keys(): # then we need to make a list
            prmclassestoboundclasses[tuple(valclasses[::-1])].append(boundclasses)
        else:
            prmclassestoboundclasses[tuple(valclasses)]=[]
            prmclassestoboundclasses[tuple(valclasses)].append(boundclasses)
    
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

        
        if 'torsion' in line and 'unit' not in line:
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

def GrabParametersFromPrmFile(poltype,bondprmclassestoboundclasses,angleprmclassestoboundclasses,torsionprmclassestoboundclasses,ligboundaryatomtypes,poltypetoprmtype,proteintypestoframedefforprmfile,check):
    if check==False:
        temp=open(poltype.amoebabioprmpath,'r')
    else:
        temp=open(poltype.ModifiedResiduePrmPath,'r')
    results=temp.readlines()
    temp.close()
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    pitorsprms=[]
    mpoleprms=[]
    opbendprms=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        linesplitall=re.split(r'(\s+)', line)
        if 'bond' in line and 'cubic' not in line and 'quartic' not in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundbond=False
            if tuple(bondclasslist) in bondprmclassestoboundclasses.keys():
                bondtup=tuple(bondclasslist)
                foundbond=True
            elif tuple(bondclasslist[::-1]) in bondprmclassestoboundclasses.keys():
                bondtup=tuple(bondclasslist[::-1])
                foundbond=True
            if foundbond==True:
                boundclasses=bondprmclassestoboundclasses[bondtup]
                for boundcls in boundclasses:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    newline=''.join(linesplitall)
                    bondprms.append(newline)           
        elif 'angle-cubic' not in line and 'angle-quartic' not in line and 'pentic' not in line and 'sextic' not in line and ('angle' in line or 'anglep' in line) :
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            foundangle=False
            if tuple(angleclasslist) in angleprmclassestoboundclasses.keys():
                angletup=tuple(angleclasslist)
                foundangle=True
            elif tuple(angleclasslist[::-1]) in angleprmclassestoboundclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundangle=True
            if foundangle==True:
                boundclasses=angleprmclassestoboundclasses[angletup]
                for boundcls in boundclasses:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    newline=''.join(linesplitall)
                    angleprms.append(newline) 
        
        elif 'torsion' in line and 'torsionunit' not in line:
            torsionclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
            foundtorsion=False
            if tuple(torsionclasslist) in torsionprmclassestoboundclasses.keys():
                torsiontup=tuple(torsionclasslist)
                foundtorsion=True
            elif tuple(torsionclasslist[::-1]) in torsionprmclassestoboundclasses.keys():
                torsiontup=tuple(torsionclasslist[::-1])
                foundtorsion=True
            if foundtorsion==True:   
                boundclasses=torsionprmclassestoboundclasses[torsiontup]
                for boundcls in boundclasses:
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
            if tuple(angleclasslist) in angleprmclassestoboundclasses.keys():
                angletup=tuple(angleclasslist)
                foundstrbnd=True
            elif tuple(angleclasslist[::-1]) in angleprmclassestoboundclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundstrbnd=True
            if foundstrbnd==True:
                boundclasses=angleprmclassestoboundclasses[angletup]
                for boundcls in boundclasses:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    newline=''.join(linesplitall)
                    strbndprms.append(newline) 

               
        elif 'pitors' in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundpitors=False
            if tuple(bondclasslist) in bondprmclassestoboundclasses.keys():
                bondtup=tuple(bondclasslist)
                foundpitors=True
            elif tuple(bondclasslist[::-1]) in bondprmclassestoboundclasses.keys():
                bondtup=tuple(bondclasslist[::-1])
                foundpitors=True
            if foundpitors==True:
                boundclasses=bondprmclassestoboundclasses[bondtup]
                for boundcls in boundclasses:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    newline=''.join(linesplitall)
                    pitorsprms.append(newline)

        elif 'opbend' in line and 'opbendtype' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundopbend=False
            reversedopbend=False
            if tuple(bondclasslist) in bondprmclassestoboundclasses.keys():
                bondtup=tuple(bondclasslist)
                foundopbend=True
            elif tuple(bondclasslist[::-1]) in bondprmclassestoboundclasses.keys():
                bondtup=tuple(bondclasslist[::-1])
                foundopbend=True
                reversedopbend=True # need to be careful because for example 3 409 0 0 will be different than 409 3 0 0, so need to swap the boundindexes if the bondtypelist is reversed

            if foundopbend==True:
                boundclasses=bondprmclassestoboundclasses[bondtup]
                for boundcls in boundclasses:
                    if reversedopbend==True:
                        linesplitall[4]=str(boundcls[0])    
                        linesplitall[2]=str(boundcls[1]) 
                    else:
                        linesplitall[2]=str(boundcls[0])    
                        linesplitall[4]=str(boundcls[1])  
                    newline=''.join(linesplitall)
                    opbendprms.append(newline)
 
        elif 'multipole' in line:
            newlinesplit=linesplit[1:-1]
            frames=[int(i) for i in newlinesplit]
            grabit=False
            if tuple(frames) in proteintypestoframedefforprmfile.keys():
                theframe=tuple(frames)
                grabit=True
            elif tuple(frames[::-1]) in proteintypestoframedefforprmfile.keys():
                theframe=tuple(frames[::-1])
                grabit=True
            if grabit==True:
                chgpartofline=linesplitall[-2:] # include space
                phrasepartofline=linesplitall[:2] # include space
                spacetoadd='   '
                newlist=[]
                newlist.extend(phrasepartofline)
                framedef=proteintypestoframedefforprmfile[atomtype]
                for typenum in framedef:
                    newlist.append(str(typenum))
                    newlist.append(spacetoadd)
                newlist=newlist[:-1] # remove last space then add space before chg
                newlist.extend(chgpartofline)
                newline=''.join(newlist)
                mpolelist=[newline,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
                for mpoleline in mpolelist:
                    mpoleprms.append(mpoleline)


    
    return bondprms,angleprms,torsionprms,strbndprms,pitorsprms,mpoleprms,opbendprms

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



def GrabMultipoleFrameDefintions(poltype,key,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,modproidxs,ligidxtoproidx,prosideboundligidx,proboundidxtoprotype,mincorenumber,maxcorenumber):
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


def WriteNewKeyFile(poltype,atomdefs,bondprms,angleprms,torsionprms,strbndprms,pitorprms,opbendprms,polarizeprms,vdwprms,mpoleprms,writekey,amoebabioprmpath):
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
    firstpasspitors=True
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
    passedpitorsblock=False
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
        elif 'pitors' in line and firstpasspitors==False:
            firstpasspitors=True
            if line not in linelist:
                linelist.append(line)
        elif 'pitors' not in line and firstpasspitors==True and passedpitorsblock==False:
            passedpitorsblock=True
            for pitorsline in pitorprms:
                linelist.append(pitorsline)
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


def WriteToPrmFile(poltype,atomdefs,bondprms,angleprms,torsionprms,strbndprms,pitorprms,opbendprms,polarizeprms,vdwprms,mpoleprms,prmpath):
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
    firstpasspitors=False
    firstpassmpole=False
    firstpassopbend=False
    passedatomdefblock=False
    passedbondblock=False
    passedangleblock=False
    passedtorsionblock=False
    passedstrbndblock=False
    passedpolarizeblock=False
    passedvdwblock=False
    passedpitorsblock=False
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
            elif 'pitors' in line and firstpasspitors==False:
                firstpasspitors=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            elif 'pitors' not in line and firstpasspitors==True and passedpitorsblock==False:
                passedpitorsblock=True
                for pitorsline in pitorprms:
                    temp.write(pitorsline)
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

def GenerateFragBabel(poltype,molindexlist,mol):
    molidxtonewmolidx={}
    newmol=openbabel.OBMol() # define new OBMol object for the fragment
    atomlist=[] # list of atom objects from mol object
    newatomlist=[] # list of blank atom objects for fragment mol object
    count=1
    for index in molindexlist: # iterate over indexes in torsion
        atom=mol.GetAtom(index) # grab the atom object via index number
        molidx=atom.GetIdx()
        atomlist.append(atom) # append the atom object to list
        newatom=newmol.NewAtom()
        newatom=newatom.Duplicate(atom)
        newatomlist.append(newatom) # just put into blank atom objects
        molidxtonewmolidx[molidx]=count
        count+=1

    bondorderidxdic={} # key=(atom index1 of bond,atom index2 of bond), value=bondorder
    iterbond = openbabel.OBMolBondIter(mol) # iterator for all bond objects in the molecule
    for bond in iterbond:
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        aidx=a.GetIdx()
        bidx=b.GetIdx()
        if aidx in molindexlist and bidx in molindexlist: # check to make sure we want these atoms
            newaidx=molindexlist.index(aidx)+1 # new atom indexes are always in the order of atoms added to molecule via newatomlist above, +1 is because python starts at 0, atom indexes start at 1
            newbidx=molindexlist.index(bidx)+1
            bondorder=bond.GetBondOrder()
            bondorderidxdic[(newaidx,newbidx)]=bondorder

        else:
            continue

    for key in bondorderidxdic: # add back the bond between atoms in original mol object to the fragment mol object
        key=list(key)
        newaidx=key[0]
        newbidx=key[1]
        bondorder=bondorderidxdic[(newaidx,newbidx)]
        newmol.AddBond(newaidx,newbidx,bondorder)
    
    return newmol,molidxtonewmolidx


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
                numberofmodatoms=CountNumberOfAtomsInModifiedResidue(poltype,bondset,modproidxs)
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
                            if angleset[1] in proboundidxs:
                                listofanglesforprm.append(angleset)
                            else: # then boundary is start or end of angle
                                if numberofmodatoms==1:
                                    listofanglesforprm.append(angleset)

            
                        nextnextneighbs=[nextnextatom for nextnextatom in openbabel.OBAtomAtomIter(nextneighb)]
                        for nextnextatom in nextnextneighbs:
                            nextnextatomidx=nextnextatom.GetIdx()
                            if nextnextatomidx!=nidx:                
                                torsionset=[nextnextatomidx,nextneighbidx,nidx,proidx]
                                numberofboundatoms=CountNumberOfAtomsInBoundryList(poltype,torsionset,proboundidxs)
                                numberofmodatoms=CountNumberOfAtomsInModifiedResidue(poltype,torsionset,modproidxs)
                                if numberofboundatoms==1 and torsionset not in listoftorsionsforkey and torsionset[::-1] not in listoftorsionsforkey and torsionset not in listoftorsionsforprm and torsionset[::-1] not in listoftorsionsforprm:
                                    if torsionset[0] in proboundidxs or torsionset[-1] in proboundidxs:
                                        if numberofmodatoms==4:
                                            listoftorsionsforkey.append(torsionset)
                                        else: # then number of mod atoms is 1 since boundary atom is in modproidxs
                                            listoftorsionsforprm.append(torsionset)
                                    else: # then the boundry atom must be in the middle two
                                        if numberofmodatoms==2:
                                            listoftorsionsforprm.append(torsionset)
                                        else:
                                            listoftorsionsforkey.append(torsionset)
    return listoftorsionsforkey,listofbondsforprm,listofanglesforprm,listoftorsionsforprm


def GrabProSideBoundIdxs(poltype,proOBmol,boundaryatomidxs):
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
    knownresiduesymbs=[]
    temp=open(poltype.libpath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            resname=linesplit[4]
            if resname[0]=='N' or resname[0]=='C':
                resname=resname[1:]
            if resname not in knownresiduesymbs:
                knownresiduesymbs.append(resname)

    return knownresiduesymbs

def MatchCorrectProteinAtomsToCorrespondingHydrogen(poltype,proboundidxs,proOBmol,polOBmol,proidxtoligidx,ligidxtoproidx):
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
    for line in results:
        if 'ATOM' in line:
            atomlabel=line[12:16].rstrip().lstrip()
            atomindex=int(line[6:11].rstrip().lstrip())    
            resnumber=int(line[23:26+1])
            if resnumber==modresnumber-1:
                if atomlabel=='C':
                    backboneprobound.append(atomindex)
            elif resnumber==modresnumber+1:
                if atomlabel=='N':
                    backboneprobound.append(atomindex)

    return backboneprobound

def RemoveUnwantedHydrogens(poltype,backboneligbound,polOBmol,proidxtoligidx,refname,obConversion):
    iteratom = openbabel.OBMolAtomIter(polOBmol)
    hydtokeep=[]
    for atom in iteratom:
        atomidx=atom.GetIdx()
        if atomidx in backboneligbound:
            iteratomatom=openbabel.OBAtomAtomIter(atom)
            for natom in iteratomatom:
                if natom.GetAtomicNum()==1:
                    hydtokeep.append(natom.GetIdx())
    iteratom = openbabel.OBMolAtomIter(polOBmol)
    deletehyds=[]
    for atom in iteratom:
        if atom.GetIdx() not in proidxtoligidx.values() and atom.GetIdx() not in hydtokeep:
            deletehyds.append(atom.GetIdx())
    for atomidx in deletehyds:
        atom=polOBmol.GetAtom(atomidx)
        polOBmol.DeleteAtom(atom)
    molname=refname.replace('.pdb','.sdf')
    obConversion.WriteFile(polOBmol,molname)
    newproidxtoligidx={}
    if len(deletehyds)!=0:
        for atomidx in deletehyds:
            for proidx in proidxtoligidx.keys():
                ligidx=proidxtoligidx[proidx]
                if ligidx>atomidx:
                    newproidxtoligidx[proidx]=ligidx-1
                else:
                    newproidxtoligidx[proidx]=ligidx
        proidxtoligidx=newproidxtoligidx
    return proidxtoligidx 


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
    match=None
    for match in allindexes:
        foundmatch=True
        for index in match:
            if index not in pdbindexes:
                foundmatch=False
        if foundmatch==True:
            return match
    return match

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
            chainid=line[21]
            if resid==resnumber and passedfirstres==False:
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
            

def MoveModifiedLines(poltype,lines,resid,pdbname):
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
            if resnumber==1:
                passedfirstres=True
            if resnumber==resid+1 and wrotelines==False:
                wrotelines=True
                for l in lines:
                    temp.write(l)
            if passedfirstres==False and resnumber==resid:
                pass
            elif passedfirstres==True and resnumber==resid: # make sure that pdbcode is not GLY but instead modified PDB code
                newline=line[:17]+poltype.modifiedresiduepdbcode+line[20:]
                temp.write(newline)
            elif passedfirstres==True and resnumber!=resid:
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


def GenerateModifiedProteinPDB(poltype):      
    sidechainindexes=GrabSideChainIndexes(poltype,int(poltype.mutatedresiduenumber),poltype.unmodifiedproteinpdbname)
    atomindexestokeepafter=FindAtomIndexesToKeepAfterMutation(poltype,int(poltype.mutatedresiduenumber),sidechainindexes,poltype.unmodifiedproteinpdbname) # only do this for indexes after residue
    unmodpdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.SetOutFormat('pdb')
    obConversion.ReadFile(unmodpdbmol,poltype.unmodifiedproteinpdbname)
    fragmolafter,unmodidxtonewidx=GenerateFragBabel(poltype,atomindexestokeepafter,unmodpdbmol)
    refname='RefAfter.pdb'
    obConversion.WriteFile(fragmolafter,refname)
    chainid=GrabChainId(poltype,int(poltype.mutatedresiduenumber),refname)


    inFormat = obConversion.FormatFromExt(poltype.mutatedsidechain)
    obConversion.SetInFormat(inFormat)        
    sidechainmol=openbabel.OBMol()
    obConversion.ReadFile(sidechainmol,poltype.mutatedsidechain)
    sidechainmol.AddHydrogens()
    backbonesmiles='C(C=O)N'
    backboneindexessidechain=MatchSMARTSToOBMol(poltype,backbonesmiles,sidechainmol)
    sidechainmol=DeleteHydrogensAttachedToBackbone(poltype,sidechainmol,backboneindexessidechain) # we will use H on backbone from referernce pdb not mutated sidechain pdb
    refmol=openbabel.OBMol()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(refmol,refname)
    
    pdbname='MutatedSideChain.pdb'
    obConversion.WriteFile(sidechainmol,pdbname)
    ModifyAlignmentPDB(poltype,pdbname)
    allindexes=MatchSMARTSToOBMolGrabAllMatches(poltype,backbonesmiles,refmol) 
    pdbindexes=GrabPDBIndexes(poltype,int(poltype.mutatedresiduenumber),poltype.unmodifiedproteinpdbname)
    backboneindexesreference=FindCorrectReferenceIndexes(poltype,allindexes,pdbindexes)

    selectstringsidechainlist=[]
    for index in backboneindexessidechain:
        selectstringsidechainlist.append('bynum '+str(index))
    
    selectstringreferencelist=[]
    for index in backboneindexesreference:
        selectstringreferencelist.append('bynum '+str(index))
    mergestring='not ('
    for index in backboneindexessidechain:
        mergestring+='bynum'+' '+str(index)+' '+'or'+' '
    mergestring=mergestring[:-4]
    mergestring+=')'
    ref = Universe(refname)
    sc = Universe(pdbname)
    alignto(sc,ref,select={"mobile": selectstringsidechainlist,"reference": selectstringreferencelist})
    u = Merge(sc.select_atoms(mergestring),ref.select_atoms('all'))
    u.trajectory.ts.dimensions = ref.trajectory.ts.dimensions
    output=poltype.unmodifiedproteinpdbname.replace('.pdb','_Mutated.pdb')
    u.atoms.write(output)
    lines=GrabModifiedResidueLines(poltype,output,int(poltype.mutatedresiduenumber)) # MDAnalysis puts merged lines at top of pdb, need to move them and change chain back to original chain
    lines=ModifyChainIdentifier(poltype,lines,chainid)
    MoveModifiedLines(poltype,lines,int(poltype.mutatedresiduenumber),output)
    finalmol=openbabel.OBMol() # need to pass through babel again to change atom indexes
    obConversion.ReadFile(finalmol,output)
    finalname=output.replace('.pdb','_Final.pdb')
    obConversion.WriteFile(finalmol,finalname)
    poltype.modifiedproteinpdbname=finalname
    return poltype.modifiedproteinpdbname


def GenerateModifiedProteinPoltypeInput(poltype):

    if poltype.unmodifiedproteinpdbname!=None:
        poltype.modifiedproteinpdbname=GenerateModifiedProteinPDB(poltype)
    knownresiduesymbs=GrabKnownResidueSymbs(poltype)

    # first need to generate protein XYZ file from PDB
    # so we have proidxtotypeidx, however since we will be transfering some of the protein parameters near boundary of poltype job and new residue in protein, we need to no which protein indexes correspond to type numbers in original protein parameter file vs what protein indexes will have new type numbers for the new residue

    modproidxs,modresnumber,modresiduelabel=GrabModifiedResidueProteinIndexes(poltype,poltype.modifiedproteinpdbname,knownresiduesymbs) # these indexes need to be excluded when calling pdbxyz.x then added back
    obConversion = openbabel.OBConversion()
    proOBmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(poltype.modifiedproteinpdbname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(proOBmol,poltype.modifiedproteinpdbname) 
    
    proboundidxs=FindBoundaryAtomIdxs(poltype,proOBmol,modproidxs)
    # now to reduce the number of atoms needed for QM and paramterization, just only include up to three atoms away starting from boundary atom on modified protein side



    modneighbidxs=NeighboringModResidueAtomIndexes(poltype,modresnumber,poltype.modifiedproteinpdbname,modproidxs,proOBmol,proboundidxs)

    # now need to grab residue number of modified residue, then use that to grab atom indexes of neighboring residues, then add those indexes to modified residue atom indexes and make fragment molecule with those atoms for poltype

    # so first step is to identify boundary atoms in a list iterate through all ligand atoms and if it has a neighbor in unmatched list, then it is a boundary atom

    choppedfragidxs=modproidxs+modneighbidxs
    chopmodmol,proidxtoligidx=GenerateFragBabel(poltype,choppedfragidxs,proOBmol) 
    obConversion.SetOutFormat("pdb")
    refname='ModifiedRes.pdb'
    obConversion.WriteFile(chopmodmol,refname)
    modmol=openbabel.OBMol()
    obConversion.ReadFile(modmol,refname)
    obConversion.SetOutFormat("mol")
    corename='ModifiedResNoHyd.mol'
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
    inFormat = obConversion.FormatFromExt(refname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(polOBmol,refname) 
    polOBmol.AddHydrogens()
    molname=refname.replace('.pdb','.sdf')
    obConversion.WriteFile(polOBmol,molname)
    proidxtoligidx= RemoveUnwantedHydrogens(poltype,backboneligbound,polOBmol,proidxtoligidx,refname,obConversion)
    # for the extra hydrogen on the backbone N and the backbone carbonyl carbon we need to add a map from those hydrogen indexes to the corresponding atom that should be in the protein
    
    # first using the carbonyl carbon and backbone N of modprotein indexes indexes iterate over the neighbors and if the neighbor is not already in the dictionary then that must be the atom index that matches to the hydrogen indexes
    ligidxtoproidx = dict([v,k] for k,v in proidxtoligidx.items())  
    proidxtoligidx,addedproatoms=MatchCorrectProteinAtomsToCorrespondingHydrogen(poltype,backboneprobound,proOBmol,polOBmol,proidxtoligidx,ligidxtoproidx)
    ligidxtoproidx = dict([v,k] for k,v in proidxtoligidx.items()) 

    return knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check



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

def GrabProBoundIdxToProType(poltype,proboundidxs,proOBmol,modproidxs):
    proboundidxtoprotype={}
    for proboundidx in proboundidxs:
        proatom=proOBmol.GetAtom(proboundidx)
        if proatom.GetAtomicNum()==6:
            proboundidxtoprotype[proboundidx]=9
        elif proatom.GetAtomicNum()==7:
            proboundidxtoprotype[proboundidx]=7
        atomatomiter=openbabel.OBAtomAtomIter(proatom)
        for natom in atomatomiter:
            natomidx=natom.GetIdx()
            if natomidx in modproidxs:
                if proatom.GetAtomicNum()==7:
                    if natom.GetAtomicNum()==1:
                        proboundidxtoprotype[natomidx]=10
                    elif natom.GetAtomicNum()==6:
                        proboundidxtoprotype[natomidx]=8
                elif proatom.GetAtomicNum()==6:
                    if natom.GetAtomicNum()==8:
                        proboundidxtoprotype[natomidx]=11
                    elif natom.GetAtomicNum()==6:
                        proboundidxtoprotype[natomidx]=8

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

def GrabCoreParameters(poltype,key5fname):
    atomdefs=[]
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    pitorprms=[]
    opbendprms=[]
    polarizeprms=[]
    vdwprms=[]
    mpoleprms=[]
    temp=open(key5fname,'r')
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
        elif 'pitor' in line:
            pitorprms.append(line)
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

    return atomdefs,bondprms,angleprms,torsionprms,strbndprms,pitorprms,opbendprms,polarizeprms,vdwprms,mpoleprms
        
def GenerateModifiedProteinXYZAndKey(poltype,knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check):

    
    # WARNING, may need to append the bio18 prm file to the end of final key file, for some reason cannot refererence the parameter file with key word parameters
    obConversion = openbabel.OBConversion()
    ligOBmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(molname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(ligOBmol,molname) 

    ligandindexes=set(range(1,ligOBmol.NumAtoms()))

    
    proboundidxtoprotype=GrabProBoundIdxToProType(poltype,proboundidxs,proOBmol,modproidxs)
    if check==False:
        ligidxtotypeidx=GenIndexToTypeIndexDic(poltype,poltype.tmpxyzfile)
    else:
        ligidxtotypeidx=ReadSMARTSToTypeLib(poltype,smarts,modmol)

    modproidxtotypenumber=AssignNewTypeNumbers(poltype,modproidxs,ligidxtotypeidx,proidxtoligidx) # will try to keep modified residue type numbers as type numbers in poltype job but then for boudnary parts need to assign new type numbers 
    proidxtoprotype=GrabProteinTypeNumbers(poltype,poltype.modifiedproteinpdbname,knownresiduesymbs,poltype.libpath,modproidxtotypenumber)
    GenerateProteinTinkerXYZFile(poltype,poltype.modifiedproteinpdbname,modproidxtotypenumber,poltype.amoebabioprmpath,proidxtoprotype,knownresiduesymbs)

    listoftorsionsforkey,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForValenceTermsAcrossBoundary(poltype,ligOBmol,proOBmol,ligidxtoproidx,proboundidxs,boundaryatomidxs,modproidxs) # execlude all cases where only 1 atom is in modproidxs and the others are in protein, seperate those into different lists since those can be transfered from prm file and the others from the key file

    # for stretch bend (strbnd) we can use the angle types to grab the strbnd parameters

    # for  polarize just use list of single atom types, (the list of boundary atom types), do not need to transfer vdw, that does not cross boundary polarize neighbors do cross the boundary

    # for pitor use bond list (just two consecutive atoms)


    print('listofbondsforprm ',listofbondsforprm)
    print('listofanglesforprm ',listofanglesforprm)
    print('listoftorsionsforkey ',listoftorsionsforkey)
    print('listoftorsionsforprm ',listoftorsionsforprm)

    
    
    # now convert protein index numbers to ligand index, then we will convert ligand index to type numbers to grab parameters. We also need dictionaries that map each type in ligand type numbers to what will be represented in final key file (so for a torsion across a boundary we want the ligand type numbers to grab parmeters but we want to have protein type numbers for part of torsion on the protein side in our key file)

    # these parameters will be taken from the key file so need to no what to convert the indexes to once they are taken from the key file
    if check==False:
        maxcorenumber=GrabMaxTypeNumber(poltype,poltype.key5fname)
        mincorenumber=GrabMinTypeNumber(poltype,poltype.key5fname)
        key=poltype.key5fname
        writekey=poltype.key5fname.replace('.key_5','.key_6')
        if not os.path.isfile(poltype.ModifiedResiduePrmPath):
            shutil.copy(poltype.amoebabioprmpath,poltype.ModifiedResiduePrmPath)
        if not os.path.isfile(poltype.SMARTSToTypelibpath):
            temp=open(poltype.SMARTSToTypelibpath,'a')
            temp.write('\n')
            temp.close()

    else:
        maxcorenumber=GrabMaxTypeNumberModifiedCore(poltype,poltype.ModifiedResiduePrmPath)
        mincorenumber=GrabMinTypeNumberModifiedCore(poltype,poltype.ModifiedResiduePrmPath)
        key=poltype.ModifiedResiduePrmPath
        writekey=poltype.key2fname

    torsionligtypestoboundtypesforkey=GrabLigandTypeToBoundryTypeForKeyFile(poltype,listoftorsionsforkey,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,mincorenumber,maxcorenumber)
    prmtypetoprmclass=GrabClassNumbersForProtein(poltype,proidxtoprotype,check) # need this for grepping torsion parameters from prm file
    bondprmclassestoboundclasses=GrabPrmClassToBoundryClassForPrmFile(poltype,listofbondsforprm,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber)
    angleprmclassestoboundclasses=GrabPrmClassToBoundryClassForPrmFile(poltype,listofanglesforprm,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber)
    torsionprmclassestoboundclasses=GrabPrmClassToBoundryClassForPrmFile(poltype,listoftorsionsforprm,proidxtoligidx,ligidxtotypeidx,proidxtoprotype,poltype.amoebabioprmpath,modproidxs,proboundidxtoprotype,prmtypetoprmclass,mincorenumber,maxcorenumber)


    print('bondprmclassestoboundclasses ',bondprmclassestoboundclasses)
    print('angleprmclassestoboundclasses ',angleprmclassestoboundclasses)
    print('torsionligtypestoboundtypesforkey ',torsionligtypestoboundtypesforkey)
    print('torsionprmclassestoboundclasses ',torsionprmclassestoboundclasses)
    ligboundaryatomtypes=[ligidxtotypeidx[i] for i in boundaryatomidxs]

    # for the polarize and multipole we need to not only have list of atomindexes (ligand indexes) that are boundary atoms on the ligand side, but also the protein side (mulitpole and polarize use neighbors that may go across the boundaries)
       
    prosideboundligidx=GrabProSideBoundIdxs(poltype,proOBmol,boundaryatomidxs,ligidxtoproidx,modproidxs,proidxtoligidx)
    prmtypetopoltype=GrabResidueTypesToPolarizeTypesAcrossBoundary(poltype,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,ligidxtoproidx,key,mincorenumber,maxcorenumber)
    prmtypetoprmdeflines=GrabPolarizeDefinitionLinesFromPrmFile(poltype,prmtypetopoltype,poltype.amoebabioprmpath)
    prmdeflines=AppendPoltypeTypeNumbers(poltype,prmtypetoprmdeflines,prmtypetopoltype) # now also grab them from the key file and then append the prm types then combine with prmdeflines here
    poltypetoprmtype = dict([v,k] for k,v in prmtypetopoltype.items())
        
    atomtypetoframedef,proteintypestoframedefforprmfile=GrabMultipoleFrameDefintions(poltype,key,boundaryatomidxs,ligidxtotypeidx,proidxtoprotype,modproidxs,ligidxtoproidx,prosideboundligidx,proboundidxtoprotype,mincorenumber,maxcorenumber) # for multipole frames on the boundary (on ligand side), we need to grab the original index frame definitions then convert to protein type number when appropriate, also need to consider atoms on other side of boundary (protein side)

    # need to add charge
    stitchtorsionprms,stitchpolarizeprms,stitchmpoleprms=GrabParametersFromKeyFile(poltype,key,torsionligtypestoboundtypesforkey,ligboundaryatomtypes,poltypetoprmtype,atomtypetoframedef,mincorenumber,maxcorenumber)
    stitchbondprms,stitchangleprms,stitchnewtorsionprms,stitchstrbndprms,stitchpitorprms,stitchnewmpoleprms,stitchopbendprms=GrabParametersFromPrmFile(poltype,bondprmclassestoboundclasses,angleprmclassestoboundclasses,torsionprmclassestoboundclasses,ligboundaryatomtypes,poltypetoprmtype,proteintypestoframedefforprmfile,check)
    stitchtorsionprms=stitchtorsionprms+stitchnewtorsionprms
    stitchmpoleprms=stitchnewmpoleprms+stitchmpoleprms
    stitchpolarizeprms=prmdeflines+stitchpolarizeprms
    stitchatomdefs=[]
    stitchvdwprms=[]
    WriteNewKeyFile(poltype,stitchatomdefs,stitchbondprms,stitchangleprms,stitchtorsionprms,stitchstrbndprms,stitchpitorprms,stitchopbendprms,stitchpolarizeprms,stitchvdwprms,stitchmpoleprms,writekey,poltype.amoebabioprmpath)
    if check==False:
        coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,corepitorprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms=GrabCoreParameters(poltype,poltype.key5fname)
        oldtypetonewtype=ShiftPoltypeNumbers(poltype)
        # now I need to convert all arrays
        coreresult=ShiftParameterDefintions(poltype,[coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,corepitorprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms],oldtypetonewtype)
        coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,corepitorprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms=coreresult[:] 
        stitchresult=ShiftParameterDefintions(poltype,[stitchatomdefs,stitchbondprms,stitchangleprms,stitchtorsionprms,stitchstrbndprms,stitchpitorprms,stitchopbendprms,stitchpolarizeprms,stitchvdwprms,stitchmpoleprms],oldtypetonewtype)
        stitchatomdefs,stitchbondprms,stitchangleprms,stitchtorsionprms,stitchstrbndprms,stitchpitorprms,stitchopbendprms,stitchpolarizeprms,stitchvdwprms,stitchmpoleprms=stitchresult[:]
        writekey=poltype.key5fname.replace('.key_5','.key_7')
        WriteNewKeyFile(poltype,coreatomdefs+stitchatomdefs,corebondprms+stitchbondprms,coreangleprms+stitchangleprms,coretorsionprms+stitchtorsionprms,corestrbndprms+stitchstrbndprms,corepitorprms+stitchpitorprms,coreopbendprms+stitchopbendprms,corepolarizeprms+stitchpolarizeprms,corevdwprms+stitchvdwprms,corempoleprms+stitchmpoleprms,writekey,poltype.amoebabioprmpath)       
        WriteToPrmFile(poltype,coreatomdefs,corebondprms,coreangleprms,coretorsionprms,corestrbndprms,corepitorprms,coreopbendprms,corepolarizeprms,corevdwprms,corempoleprms,poltype.ModifiedResiduePrmPath)
        #print('proidxtoprotype ',proidxtoprotype)
        proidxtoprotype=ShiftDictionaryValueTypes(poltype,proidxtoprotype,oldtypetonewtype)
        ligidxtotypeidx=ShiftDictionaryValueTypes(poltype,ligidxtotypeidx,oldtypetonewtype)
        # now I need to make a text file that takes parameters and can add to library
        modresiduedic=GrabLibraryInfo(poltype,proidxtoprotype,modresiduelabel,proOBmol)
        GenerateLibFileToAdd(poltype,modresiduedic,modresiduelabel)
        GenerateSMARTSToTypeFileToAdd(poltype,modmol,ligidxtotypeidx,smarts)

def ShiftDictionaryValueTypes(poltype,dictionary,oldtypetonewtype):
    for key,value in dictionary.items():
        if value in oldtypetonewtype.keys():
            newvalue=oldtypetonewtype[value]
            dictionary[key]=newvalue
    return dictionary



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
                    oldtypenum=numpy.abs(int(element))
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
                    
    
def RepresentsInt(poltype,s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def ShiftPoltypeNumbers(poltype):
    oldtypetonewtype={}
    maxnumberfromprm=GrabMaxTypeNumber(poltype,poltype.ModifiedResiduePrmPath)
    maxnumberfromkey=GrabMaxTypeNumber(poltype,poltype.key5fname)
    minnumberfromkey=GrabMinTypeNumber(poltype,poltype.key5fname)
    typenumbers=list(range(minnumberfromkey,maxnumberfromkey+1))
    newmaxnumber=maxnumberfromprm+1
    shift=minnumberfromkey-newmaxnumber
    oldtypetonewtype={}
    for typenum in typenumbers:
        newtypenum=typenum-shift
        oldtypetonewtype[typenum]=newtypenum
    return oldtypetonewtype

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


