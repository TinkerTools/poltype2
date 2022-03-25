import os
import sys
import keyfilemodifications as keymods
import openbabel
from rdkit import Chem
import re
import itertools
import numpy as np
import shutil
import submitjobs as submit
from rdkit.Chem import rdMolTransforms
import pylab as plt
from scipy.interpolate import interp1d
from scipy.optimize import fmin

def BgnTypeToNewType(poltype,bgnstatexyz,endstatexyz):
    bgnstatexyzatominfo,bgnstateindextotypeindex,bgnstateatomnum,bgnindextocoords,bgnindextoneighbs=GrabXYZInfo(poltype,bgnstatexyz)
    endstatexyzatominfo,endstateindextotypeindex,endstateatomnum,endindextocoords,endindextoneighbs=GrabXYZInfo(poltype,endstatexyz)
    bgnstateindextoendstateindex=BgnIndexToNewIndex(poltype,bgnstatexyz,endstatexyz)
    bgnstatetypeindextoendstatetypeindex={}
    bgnindextoendcoords={}
    typestocopy=[]
    for bgnstateindex in bgnstateindextotypeindex.keys():
        endstateindex=bgnstateindextoendstateindex[bgnstateindex]
        endcoords=endindextocoords[endstateindex]
        
        bgnindextoendcoords[bgnstateindex]=endcoords
        bgnstatetypeindex=int(bgnstateindextotypeindex[bgnstateindex])
        endstatetypeindex=int(endstateindextotypeindex[endstateindex])
        if bgnstatetypeindex not in bgnstatetypeindextoendstatetypeindex.keys():
            bgnstatetypeindextoendstatetypeindex[bgnstatetypeindex]=endstatetypeindex
        else:
            typestocopy.append(bgnstatetypeindex)
    return bgnstatetypeindextoendstatetypeindex,bgnstateindextotypeindex,endstateindextotypeindex,bgnindextoendcoords,typestocopy,endindextoneighbs


def GrabXYZInfo(poltype,xyzfile):
    temp=open(xyzfile,'r')
    xyzfileresults=temp.readlines()
    temp.close()
    xyzatominfo=[]
    indextotypeindex={}
    indextocoords={}
    indextoneighbs={}
    for lineidx in range(len(xyzfileresults)):
        line=xyzfileresults[lineidx]
        linesplit=line.split()
        if lineidx==0:
            xyzatomnum=int(linesplit[0])
        else:
            if len(linesplit)>1:
                index=int(linesplit[0])
                typenum=int(linesplit[5])
                coords=[linesplit[2],linesplit[3],linesplit[4]]
                indextotypeindex[index]=typenum
                xyzatominfo.append(linesplit[2:])
                indextocoords[index]=coords
                if len(linesplit)>=6:
                    neighbs=linesplit[6:]
                    neighbs=[int(i) for i in neighbs]
                else:
                    neighbs=[]
                indextoneighbs[index]=neighbs
    return xyzatominfo,indextotypeindex,xyzatomnum,indextocoords,indextoneighbs

def CorrectBondOrder(poltype,themol,structurefile):
    if structurefile!=None:
        obConversion = openbabel.OBConversion()
        mol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(structurefile)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(mol, structurefile)
        iterbond = openbabel.OBMolBondIter(mol)
        for bond in iterbond:
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            BO=bond.GetBondOrder()
            aidx=a.GetIdx()
            bidx=b.GetIdx()
            thebond=themol.GetBond(aidx,bidx)
            theBO=thebond.GetBondOrder()
            thebond.SetBondOrder(BO)



    return themol


def BgnIndexToNewIndex(poltype,bgnstatexyz,endstatexyz):
    bgncartxyz=poltype.ConvertTinktoXYZ(bgnstatexyz,bgnstatexyz.replace('.xyz','_cart.xyz'))
    endcartxyz=poltype.ConvertTinktoXYZ(endstatexyz,endstatexyz.replace('.xyz','_cart.xyz'))
    obConversion = openbabel.OBConversion()
    bgnmol = openbabel.OBMol()
    endmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(bgncartxyz)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(bgnmol, bgncartxyz)
    obConversion.ReadFile(endmol, endcartxyz)
    obConversion.SetOutFormat('mol')
    bgnmolfile=bgncartxyz.replace('_cart.xyz','.mol')
    endmolfile=endcartxyz.replace('_cart.xyz','.mol')
    CorrectBondOrder(poltype,bgnmol,poltype.bgnstatestructurefile)
    CorrectBondOrder(poltype,endmol,poltype.endstatestructurefile)
    obConversion.WriteFile(bgnmol,bgnmolfile)
    obConversion.WriteFile(endmol,endmolfile)
    bgnmol=Chem.MolFromMolFile(bgnmolfile,removeHs=False,sanitize=False)
    endmol=Chem.MolFromMolFile(endmolfile,removeHs=False,sanitize=False)
    bgnmol,atomindextoformalcharge=poltype.CheckLigandInputCharge(bgnmol)
    endmol,atomindextoformalcharge=poltype.CheckLigandInputCharge(endmol)
    Chem.SanitizeMol(bgnmol)
    Chem.SanitizeMol(endmol)
    matches=endmol.GetSubstructMatches(bgnmol)
    firstmatch=matches[0]
    firstmatch=[i+1 for i in firstmatch]
    indices=list(range(len(firstmatch)))
    indices=[i+1 for i in indices]
    bgnstateindextoendstateindex=dict(zip(indices,firstmatch)) 
    return bgnstateindextoendstateindex

def UpdateXYZAndKeyFiles(poltype,bgnstatetypeindextoendstatetypeindex,bgnstatexyz,bgnstatekey,bgnindextoendcoords):
    
    tempbgnstatexyz=bgnstatexyz.replace('.xyz','_TEMP.xyz')
    ReplaceXYZCoords(poltype,bgnstatexyz,bgnindextoendcoords,tempbgnstatexyz)
    newbgnxyzfile=bgnstatexyz.replace('.xyz','_new.xyz')
    newbgnkeyfile=bgnstatekey.replace('.key','_new.key')
    ReplaceKeywordsInFile(poltype,bgnstatexyz,newbgnxyzfile,bgnstatetypeindextoendstatetypeindex)
    ReplaceKeywordsInFile(poltype,bgnstatekey,newbgnkeyfile,bgnstatetypeindextoendstatetypeindex)
    return newbgnxyzfile,newbgnkeyfile

def ReplaceXYZCoords(poltype,bgnstatexyz,bgnindextoendcoords,tempname,replace=True):
    temp=open(bgnstatexyz,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempname,'w')
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        prettylinesplit=re.split(r'(\s+)', line)
        if lineidx!=0 and len(linesplit)>1:
            index=int(linesplit[0])
            coords=bgnindextoendcoords[index]
            prettylinesplit[6]=str(coords[0])
            prettylinesplit[8]=str(coords[1])
            prettylinesplit[10]=str(coords[2])
        line=''.join(prettylinesplit)
        temp.write(line)
    temp.close()
    if replace==True:
        os.remove(bgnstatexyz)
        os.rename(tempname,bgnstatexyz)



def ReplaceKeywordsInFile(poltype,bgnstate,newbgnfile,bgnstatetypeindextoendstatetypeindex):
    temp=open(bgnstate,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newbgnfile,'w')
    for line in results:
        linesplit=re.split(r'(\s+)', line)
        for eidx in range(len(linesplit)):
            e=linesplit[eidx]
            if e.isdigit():
                digit=int(e)
                if digit in bgnstatetypeindextoendstatetypeindex.keys():
                    newdigit=bgnstatetypeindextoendstatetypeindex[digit]
                    linesplit[eidx]=str(newdigit)
        line=''.join(linesplit)
        temp.write(line)
    temp.close() 


def CompareParameters(poltype,endstatexyz,endstatekey,newbgnxyzfile,newbgnkeyfile,bgnstatetypeindextoendstatetypeindex,bgnstateindextotypeindex,endstateindextotypeindex,analyzepath,energycutoff,bgnstatekey,potentialexe,poleditpath,endindextoneighbs):
    endalzout='endkeyalz.out'
    poltype.CallAnalyze(endstatexyz,endstatekey,endalzout,analyzepath,'d')
    endvdwindicestoenergy,endbondindicestoenergy,endangleindicestoenergy,endopbendindicestoenergy,endtorsionindicestoenergy=GrabOutput(poltype,endalzout) 
    endvdwtypeindicestoenergy=ConvertIndicesToTypes(poltype,endvdwindicestoenergy,endstateindextotypeindex)
    endbondtypeindicestoenergy=ConvertIndicesToTypes(poltype,endbondindicestoenergy,endstateindextotypeindex)
    endangletypeindicestoenergy=ConvertIndicesToTypes(poltype,endangleindicestoenergy,endstateindextotypeindex)
    endopbendtypeindicestoenergy=ConvertIndicesToTypes(poltype,endopbendindicestoenergy,endstateindextotypeindex)
    endtorsiontypeindicestoenergy=ConvertIndicesToTypes(poltype,endtorsionindicestoenergy,endstateindextotypeindex)
    endprmtypes=list(bgnstatetypeindextoendstatetypeindex.values())
    endvdwtypetoprmline=GrabParameterLines(poltype,endstatekey,'vdw',endprmtypes)
    endbondtypetoprmline=GrabParameterLines(poltype,endstatekey,'bond',endprmtypes)
    endangletypetoprmline=GrabParameterLines(poltype,endstatekey,'angle',endprmtypes)
    endopbendtypetoprmline=GrabParameterLines(poltype,endstatekey,'opbend',endprmtypes)
    endtorsiontypetoprmline=GrabParameterLines(poltype,endstatekey,'torsion',endprmtypes)
    endtypeindicestoenergylist=[endvdwtypeindicestoenergy,endbondtypeindicestoenergy,endangletypeindicestoenergy,endopbendtypeindicestoenergy,endtorsiontypeindicestoenergy]
    stringlist=['vdw','bond','angle','opbend','torsion']
    newbgnalzout='bgnkeyalz.out'
    poltype.CallAnalyze(newbgnxyzfile,newbgnkeyfile,newbgnalzout,analyzepath,'d')
    newbgnvdwindicestoenergy,newbgnbondindicestoenergy,newbgnangleindicestoenergy,newbgnopbendindicestoenergy,newbgntorsionindicestoenergy=GrabOutput(poltype,newbgnalzout) 
    newbgnvdwtypeindicestoenergy=ConvertIndicesToTypes(poltype,newbgnvdwindicestoenergy,bgnstateindextotypeindex)
    newbgnbondtypeindicestoenergy=ConvertIndicesToTypes(poltype,newbgnbondindicestoenergy,bgnstateindextotypeindex)
    newbgnangletypeindicestoenergy=ConvertIndicesToTypes(poltype,newbgnangleindicestoenergy,bgnstateindextotypeindex)
    newbgnopbendtypeindicestoenergy=ConvertIndicesToTypes(poltype,newbgnopbendindicestoenergy,bgnstateindextotypeindex)
    newbgntorsiontypeindicestoenergy=ConvertIndicesToTypes(poltype,newbgntorsionindicestoenergy,bgnstateindextotypeindex)
    bgnprmtypes=list(bgnstatetypeindextoendstatetypeindex.keys())
    bgnvdwtypetoprmline=GrabParameterLines(poltype,bgnstatekey,'vdw',bgnprmtypes)
    bgnbondtypetoprmline=GrabParameterLines(poltype,bgnstatekey,'bond',bgnprmtypes)
    bgnangletypetoprmline=GrabParameterLines(poltype,bgnstatekey,'angle',bgnprmtypes)
    bgnopbendtypetoprmline=GrabParameterLines(poltype,bgnstatekey,'opbend',bgnprmtypes)
    bgntorsiontypetoprmline=GrabParameterLines(poltype,bgnstatekey,'torsion',bgnprmtypes)
    bgntypeindicestoenergylist=[newbgnvdwtypeindicestoenergy,newbgnbondtypeindicestoenergy,newbgnangletypeindicestoenergy,newbgnopbendtypeindicestoenergy,newbgntorsiontypeindicestoenergy]
    bgntypeindicestocommentslist=FindComments(poltype,bgntypeindicestoenergylist,bgnstatekey,bgnprmtypes)
    endtypeindicestocommentslist=FindComments(poltype,endtypeindicestoenergylist,endstatekey,endprmtypes)
    vdwbgntypeindicestocomments=bgntypeindicestocommentslist[0]
    vdwendtypeindicestocomments=endtypeindicestocommentslist[0]
    endtypetoprmlinelist=[endvdwtypetoprmline,endbondtypetoprmline,endangletypetoprmline,endopbendtypetoprmline,endtorsiontypetoprmline]
    bgntypetoprmlinelist=[bgnvdwtypetoprmline,bgnbondtypetoprmline,bgnangletypetoprmline,bgnopbendtypetoprmline,bgntorsiontypetoprmline]
    endtypeprmtoamoeba09matchlist=FindAMOEBA09Matches(poltype,endtypetoprmlinelist) 
    bgntypeprmtoamoeba09matchlist=FindAMOEBA09Matches(poltype,bgntypetoprmlinelist)
    vdwendtypeprmtoamoeba09match=endtypeprmtoamoeba09matchlist[0]
    vdwbgntypeprmtoamoeba09match=bgntypeprmtoamoeba09matchlist[0]

    CompareVDW(poltype,bgnvdwtypetoprmline,endvdwtypetoprmline,bgnstatetypeindextoendstatetypeindex,bgnstateindextotypeindex,endstateindextotypeindex,vdwbgntypeindicestocomments,vdwendtypeindicestocomments,vdwbgntypeprmtoamoeba09match,vdwendtypeprmtoamoeba09match)
    endvdwtypetoprmline=AddPairVdwParameterLines(poltype,endvdwtypetoprmline,endstateindextotypeindex)
    bgnvdwtypetoprmline=AddPairVdwParameterLines(poltype,bgnvdwtypetoprmline,bgnstateindextotypeindex)
    endtypetoprmlinelist=[endvdwtypetoprmline,endbondtypetoprmline,endangletypetoprmline,endopbendtypetoprmline,endtorsiontypetoprmline]
    bgntypetoprmlinelist=[bgnvdwtypetoprmline,bgnbondtypetoprmline,bgnangletypetoprmline,bgnopbendtypetoprmline,bgntorsiontypetoprmline]
    
    CompareEnergies(poltype,endtypeindicestoenergylist,bgntypeindicestoenergylist,bgnstatetypeindextoendstatetypeindex,energycutoff,stringlist,endtypetoprmlinelist,bgntypetoprmlinelist,bgnstateindextotypeindex,endstateindextotypeindex,bgntypeindicestocommentslist,endtypeindicestocommentslist,bgntypeprmtoamoeba09matchlist,endtypeprmtoamoeba09matchlist)
    newbgnkeyfilename=CompareESP(poltype,newbgnxyzfile,newbgnkeyfile,endstatexyz,endstatekey,potentialexe,analyzepath,bgnstatetypeindextoendstatetypeindex,endstateindextotypeindex,poleditpath,endindextoneighbs)
    if poltype.perturbedkeyfilename!=None:
        shutil.copy(newbgnkeyfilename,poltype.perturbedkeyfilename)

    tortypes=list(endtorsiontypetoprmline.keys())
    bgntorsionindicestocommentslist=bgntypeindicestocommentslist[-1] 
    endtorsionindicestocommentslist=endtypeindicestocommentslist[-1]
    bgntransfertors=FindTransferTorsions(poltype,bgntorsionindicestocommentslist)
    endtransfertors=FindTransferTorsions(poltype,endtorsionindicestocommentslist)
    CompareDihedralScans(poltype,tortypes,analyzepath,newbgnxyzfile,newbgnkeyfile,endstatexyz,endstatekey)

def FindTransferTorsions(poltype,torsionindicestocommentslist):
    transfertors=[]
    fittors=[]
    for torsionindices,commentslist in torsionindicestocommentslist.items():
        transfer=True
        for comment in commentslist:
            if 'Fitted' in comment:
                transfer=False
        if len(commentslist)==0:
            transfer=False
        if transfer==True:
            if torsionindices not in transfertors:
                if len(torsionindices)==4:
                    transfertors.append(torsionindices)
        else:
            if len(torsionindices)==4 and len(commentslist)!=0:
                fittors.append(torsionindices)
    if len(transfertors)!=0:
        print('Transferrered torsions')
        for tor in transfertors:
            print(tor)
    if len(fittors)!=0:
        print('Fit torsions')
        for tor in fittors:
            print(tor)

    return transfertors

def GrabTorsionIndices(poltype,stateindextotypeindex,indextoneighbs,tortype):
    allindices=[]
    for typenum in tortype:
        indexes=GrabKeysFromValue(poltype,stateindextotypeindex,typenum)
        allindices.append(indexes)
    combs = list(itertools.product(*allindices))
    for comb in combs:
        poscomb=[np.abs(i) for i in comb]
        checkconsec=CheckIfAtomsConnected(poltype,poscomb,indextoneighbs)
        if checkconsec==True:
            torsion=comb[:]
        break
    return torsion


def CompareDihedralScans(poltype,tortypes,analyzepath,newbgnxyzfile,newbgnkeyfile,endstatexyz,endstatekey):
    xyzfilename=poltype.ConvertTinkerXYZToCartesianXYZ(newbgnxyzfile)
    bgnmol=poltype.ReadLigandRdkitMol(xyzfilename)
    bgnmol,atomindextoformalcharge=poltype.CheckLigandInputCharge(bgnmol)
    Chem.SanitizeMol(bgnmol)      
    xyzfilename=poltype.ConvertTinkerXYZToCartesianXYZ(endstatexyz)
    endmol=poltype.ReadLigandRdkitMol(xyzfilename)
    endmol,atomindextoformalcharge=poltype.CheckLigandInputCharge(endmol)
    Chem.SanitizeMol(endmol)     

    bgnstatexyzatominfo,bgnstateindextotypeindex,bgnstateatomnum,bgnindextocoords,bgnindextoneighbs=GrabXYZInfo(poltype,newbgnxyzfile)
    endstatexyzatominfo,endstateindextotypeindex,endstateatomnum,endindextocoords,endindextoneighbs=GrabXYZInfo(poltype,endstatexyz)
    middletypes=[]
    cutoff=1
    for tortype in tortypes:
        btype=int(tortype[1])
        ctype=int(tortype[2])
        ls=[btype,ctype]
        ls.sort()
        if ls in middletypes:
            continue
        middletypes.append(ls)
        torsion=GrabTorsionIndices(poltype,bgnstateindextotypeindex,bgnindextoneighbs,tortype)    
        b=torsion[1]
        c=torsion[2]
        batom=bgnmol.GetAtomWithIdx(b-1)
        catom=bgnmol.GetAtomWithIdx(c-1)
        isinringb=batom.IsInRing()
        isinringc=catom.IsInRing()
        if isinringb==True and isinringc==True:
            continue
        dihedral_energies_bgn,dihedral_degrees=DihedralAngleScan(poltype,newbgnxyzfile,newbgnkeyfile,torsion,bgnmol,analyzepath)
        torsion=GrabTorsionIndices(poltype,endstateindextotypeindex,endindextoneighbs,tortype)    
        b=torsion[1]
        c=torsion[2]
        batom=endmol.GetAtomWithIdx(b)
        catom=endmol.GetAtomWithIdx(c)
        isinringb=batom.IsInRing()
        isinringc=catom.IsInRing()

        dihedral_energies_end,dihedral_degrees=DihedralAngleScan(poltype,endstatexyz,endstatekey,torsion,endmol,analyzepath)
        rmse=PlotConformationalEnergy(poltype,dihedral_energies_bgn,dihedral_energies_end,tortype,dihedral_degrees)
        if rmse>cutoff:
            print('Torsion '+str(tortype)+' has RMSE(Bgn,End) of '+str(rmse))


def PlotConformationalEnergy(poltype,dihedral_energies_bgn,dihedral_energies_end,tortype,dihedral_degrees):
    newbgn=[]
    newend=[]
    newdeg=[]
    for i in range(len(dihedral_energies_bgn)):
        bgne=dihedral_energies_bgn[i]
        ende=dihedral_energies_end[i]
        ang=dihedral_degrees[i]
        if bgne>50 or ende>50:
            continue
        else:
            newbgn.append(bgne)
            newend.append(ende)
            newdeg.append(ang)
    minRMSD=0
    if len(newbgn)>3:
        a,b,c,d=tortype[:]
        figfname='EndType_'+str(a)+'-'+str(b)+'-'+str(c)+'-'+str(d)+'.png'
        plt.style.use('ggplot')
        xpoints=np.array([newdeg[i] for i in range(len(newdeg))])
        x_new = np.linspace(xpoints.min(),xpoints.max(),500)
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        l1, = ax.plot(newdeg,newbgn,'ro',color='red',label='Begin')
        f = interp1d(xpoints,newbgn, kind='quadratic')
        y_smooth=f(x_new)
        ax.plot(x_new,y_smooth,color='red')
        l2, = ax.plot(newdeg,newend,'bo',color='blue',label='End')
        f = interp1d(xpoints,newend, kind='quadratic')
        y_smooth=f(x_new)
        ax.plot(x_new,y_smooth,color='blue')
        plt.legend(handles=[l1,l2],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
        ax.set_xlabel('Dihedral Angle Degree')
        ax.set_ylabel('AMOEBA Total Energy kcal/mol')

        def RMSD(c):
            return np.sqrt(np.mean(np.square(np.add(np.subtract(newbgn,newend),c))))
        result=fmin(RMSD,.5)
        minRMSD=RMSD(result[0])
        string='RMSD(Bgn,End)=%s'%(round(minRMSD,2))
        ax.set_title(string)
        fig.savefig(figfname)

    return minRMSD




def FindAMOEBA09Matches(poltype,typetoprmlinelist):
    typeprmtoamoeba09matchlist=[]
    temp=open(poltype.amoeba09prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for typetoprmline in typetoprmlinelist:
        typeprmtoamoeba09match={}
        for typeindices,prmline in typetoprmline.items():
            typeprmtoamoeba09match[typeindices]=[]
            typeprmtoamoeba09match[typeindices[::-1]]=[]

        atomclasstocomment={}
        for lineidx in range(len(results)):
            line=results[lineidx]
            indices=None
            if '#' not in line and 'BUFF' not in line and 'scale' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line and 'unit' not in line and 'ALLINGER' not in line and 'cutoff' not in line:
                if 'atom' in line:
                    linesplit=line.split()
                    linesplit=linesplit[:-3]
                    atomclass=linesplit[2]
                    comment=linesplit[4:]
                    if atomclass not in atomclasstocomment.keys():
                        atomclasstocomment[atomclass]=[]
                    atomclasstocomment[atomclass].append(comment)

                if 'vdw' in line:
                    indices=[1]
                if 'bond ' in line or 'opbend ' in line:
                    indices=[1,2]
                elif 'angle ' in line or 'strbnd ' in line:
                    indices=[1,2,3]
                elif 'torsion ' in line:
                    indices=[1,2,3,4]
                if indices!=None:
                    linesplit=line.split()
                    typenums=[linesplit[i] for i in indices]
                    prms=linesplit[indices[-1]+1:]
                    for typeindices,prmline in typetoprmline.items():
                        if 'vdw' in prmline:
                            newindices=[1]
                        if 'bond ' in prmline or 'opbend ' in prmline:
                            newindices=[1,2]
                        elif 'angle ' in prmline or 'strbnd ' in prmline:
                            newindices=[1,2,3]
                        elif 'torsion ' in prmline:
                            newindices=[1,2,3,4]
                        if newindices!=None:
                            newlinesplit=prmline.split()
                        
                            if len(newindices)==len(indices):
                                newprms=newlinesplit[newindices[-1]+1:]
                                prms=[float(i) for i in prms]
                                newprms=[float(i) for i in newprms]
                                if len(prms)==len(newprms):
                                    match=True
                                    for i in range(len(prms)):
                                        cp=prms[i]
                                        np=newprms[i]
                                        if cp!=np:
                                            match=False
                                    if match==True:
                                        comments=[atomclasstocomment[i] for i in typenums]
                                        allin=True
                                        if comments not in typeprmtoamoeba09match[typeindices]:
                                            typeprmtoamoeba09match[typeindices].append(comments)
                                            if 'opbend' not in line and comments not in typeprmtoamoeba09match[typeindices[::-1]]:
                                                typeprmtoamoeba09match[typeindices[::-1]].append(comments)
        

        typeprmtoamoeba09matchlist.append(typeprmtoamoeba09match)

    return typeprmtoamoeba09matchlist



def FindComments(poltype,typeindicestoenergylist,statekey,prmtypes):
    temp=open(statekey,'r')
    results=temp.readlines()
    temp.close()
    typeindicestocommentslist=[]
    for typeindicestoenergy in typeindicestoenergylist:
        typeindicestocomments={}
        for typeindices in typeindicestoenergy.keys():
            typeindicestocomments[typeindices]=[]

        for lineidx in range(len(results)):
            line=results[lineidx]
            indices=None
            if '#' in line or 'cutoff' in line:
                continue
            if 'vdw' in line and 'vdwtype' not in line and 'scale' not in line:
                indices=[1]
            if 'bond ' in line or 'opbend ' in line:
                indices=[1,2]
            elif 'angle ' in line or 'strbnd ' in line:
                indices=[1,2,3]
            elif 'torsion ' in line:
                indices=[1,2,3,4]
            if indices!=None:
                linesplit=line.split()
                typenums=[linesplit[i] for i in indices]
                typenums=[int(i) for i in typenums]
                typenums=tuple(typenums)
                found=False
                if 'vdw' in line and 'vdwtype' not in line:
                    indiceslist=prmtypes[:]
                    indiceslist=[tuple([i]) for i in indiceslist]
                else:
                    indiceslist=list(typeindicestoenergy.keys())
                for typeindices in indiceslist:
                    if typenums==typeindices:
                        found=True
                    else:
                        if 'opbend' not in line:
                            if typenums==typeindices[::-1]:
                                found=True
                    if found==True:
                        lines=[]
                        for i in range(1,4):
                            prevline=results[lineidx-i] 
                            if '#' in prevline and '##' not in prevline:
                                if prevline!='\n' and '] = [' not in prevline and 'WARNING' not in prevline:
                                    lines.append(prevline)
                        typeindicestocomments[typeindices]=lines
                        if 'opbend' not in line:
                            typeindicestocomments[typeindices[::-1]]=lines
                        break



        typeindicestocommentslist.append(typeindicestocomments)

    return typeindicestocommentslist


def CompareVDW(poltype,bgnvdwtypetoprmline,endvdwtypetoprmline,bgnstatetypeindextoendstatetypeindex,bgnstateindextotypeindex,endstateindextotypeindex,vdwbgntypeindicestocomments,vdwendtypeindicestocomments,vdwbgntypeprmtoamoeba09match,vdwendtypeprmtoamoeba09match):
    for bgntype,bgnprmline in bgnvdwtypetoprmline.items():
        endtype=bgnstatetypeindextoendstatetypeindex[bgntype[0]]
        endtype=tuple([endtype])
        endprmline= endvdwtypetoprmline[endtype]
        bgnsplit=bgnprmline.split()
        bgnparms=bgnsplit[2:]
        bgnparms=[float(i) for i in bgnparms]
        endsplit=endprmline.split()
        endparms=endsplit[2:]
        endparms=[float(i) for i in endparms]
        diff=False
        bgnindices=TypeIndicesToIndices(poltype,bgntype,bgnstateindextotypeindex)
        endindices=TypeIndicesToIndices(poltype,endtype,endstateindextotypeindex)
        bgncomments=vdwbgntypeindicestocomments[bgntype]
        endcomments=vdwendtypeindicestocomments[endtype]
        bgnmatches=vdwbgntypeprmtoamoeba09match[bgntype]
        endmatches=vdwendtypeprmtoamoeba09match[endtype]

        if len(bgnparms)==3 and len(endparms)==2:
            endparms.append(1)
        elif len(bgnparms)==2 and len(endparms)==3:
            bgnparms.append(1)

        
        for i in range(len(bgnparms)):
            bgn=bgnparms[i]
            end=endparms[i]
            if np.abs(bgn-end)>.00001:
                diff=True
        if diff==True:
            print('\n')
            print('Begin parameters '+bgnprmline.replace('\n','')+' '+str(bgnindices)+'\n')
            for line in bgncomments:
                print(line)
            for ls in bgnmatches:
                print('amoeba09 match '+str(ls))
            print('\n')
            print('End parameters '+endprmline.replace('\n','')+' '+str(endindices)+'\n')
            for line in endcomments:
                print(line)
            for ls in endmatches:
                print('amoeba09 match '+str(ls))

def GrabEndFrames(poltype,endstatekey,bgnstatetypeindextoendstatetypeindex):
    typenumtoframe={}
    temp=open(endstatekey,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'multipole' in line:
            firsttype=int(linesplit[1])
            if firsttype in bgnstatetypeindextoendstatetypeindex.values():
                frame=linesplit[1:-1]
                frame=[int(i) for i in frame]
                typenumtoframe[firsttype]=frame

    return typenumtoframe

def CompareESP(poltype,newbgnxyzfile,newbgnkeyfile,endstatexyz,endstatekey,potentialexe,analyzepath,bgnstatetypeindextoendstatetypeindex,endstateindextotypeindex,poleditpath,endindextoneighbs):
    newbgnxyzfilename=newbgnxyzfile.replace('.xyz','_potentialcompare.xyz')
    newbgnkeyfilename=newbgnkeyfile.replace('.key','_potentialcompare.key')
    shutil.copy(endstatexyz,newbgnxyzfilename) # need same atom order
    shutil.copy(newbgnkeyfile,newbgnkeyfilename)
    shutil.copy(newbgnkeyfile,'mutated.key')
    typenumtoframe=GrabEndFrames(poltype,endstatekey,bgnstatetypeindextoendstatetypeindex)
    atomindextoframeindices=ConvertTypesToIndices(poltype,endstateindextotypeindex,typenumtoframe,endindextoneighbs)
    if poltype.rotateframes==True:
        RotateBeginKeyFrames(poltype,poleditpath,atomindextoframeindices,newbgnxyzfilename,newbgnkeyfilename)
    endalzout='checkendalz.out' 
    poltype.CallAnalyze(endstatexyz,endstatekey,endalzout,analyzepath,'d')
    bgnalzout='checkbgnalz.out' 
    poltype.CallAnalyze(newbgnxyzfilename,newbgnkeyfilename,bgnalzout,analyzepath,'d')
    cmdstr=potentialexe+' '+'4'+' '+endstatexyz+' '+newbgnxyzfilename+' '+'Y'+' > espcompare.out'
    submit.call_subsystem(poltype,cmdstr,True,False,'espcompare.out')    
    return newbgnkeyfilename

def RotateBeginKeyFrames(poltype,poleditpath,atomindextoframeindices,bgnxyzfilename,bgnkeyfilename):
    temp=open('poledit.in','w')
    for atomindex,frameindices in atomindextoframeindices.items():
        frameindices=[str(i) for i in frameindices]
        string=' '.join(frameindices)
        temp.write(string+'\n')
    temp.write('\n') 
    temp.close()
    cmdstr=poleditpath+' '+'2'+' '+bgnxyzfilename+' '+'-k'+' '+bgnkeyfilename+' '+'<'+' '+'poledit.in'
    submit.call_subsystem(poltype,cmdstr,True,False)
    newkey=bgnkeyfilename+'_2'
    os.remove(bgnkeyfilename)
    os.rename(newkey,bgnkeyfilename)
    os.remove(bgnxyzfilename+'_2') 


def ConvertTypesToIndices(poltype,endstateindextotypeindex,typenumtoframe,endindextoneighbs):
    atomindextoframeindices={}
    for typenum,frame in typenumtoframe.items():
        allindices=[]
        neg=False
        for typenum in frame: 
            if typenum<0:
               neg=True
            else:
               neg=False
            typenum=np.abs(typenum)
             
            indexes=GrabKeysFromValue(poltype,endstateindextotypeindex,typenum)
            if neg==True:
                indexes=[-i for i in indexes]

            allindices.append(indexes)
        combs = list(itertools.product(*allindices))
        for comb in combs:
            poscomb=[np.abs(i) for i in comb]
            theset=list(set(poscomb))
            checkconsec=CheckIfAtomsConnected(poltype,poscomb,endindextoneighbs)
            if len(theset)==len(comb) and checkconsec==True:
                index=poscomb[0]
                atomindextoframeindices[index]=comb   

    return atomindextoframeindices

def CheckIfAtomsConnected(poltype,poscomb,endindextoneighbs):
    checkconsec=True
    if len(poscomb)>1:
        for i in range(len(poscomb)-1):
            index=poscomb[i]
            nextindex=poscomb[i+1]
            indexneighbs=endindextoneighbs[index]
            if nextindex not in indexneighbs:
                checkconsec=False



    return checkconsec



def GrabKeysFromValue(poltype,dic,thevalue):
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist


def AddPairVdwParameterLines(poltype,vdwtypetoprmline,stateindextotypeindex):
    newvdwtypetoprmline={}
    vdwtypes=[]
    for vdwtype,prmline in vdwtypetoprmline.items():
        vdwtype=vdwtype[0]
        vdwtypes.append(vdwtype)
    combs=itertools.combinations(vdwtypes, 2)
    for comb in combs: 
        prmlines=[]    
        for vdwtype in comb:
            vdwtypels=tuple([vdwtype])
            prmline=vdwtypetoprmline[vdwtypels]
            prmlines.append(prmline) 
        ls=tuple(list(comb))
        newprmline=''
        for prmline in prmlines:
            newprmline+=prmline+'\n'
        newvdwtypetoprmline[ls]=newprmline    

 
    return newvdwtypetoprmline

def GrabParameterLines(poltype,statekey,string,prmtypes):
    typeindicestoprmline={}
    classtotypes={}
    temp=open(statekey,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            if 'atom'==linesplit[0]:
                typenum=int(linesplit[1])
                classnum=int(linesplit[2])
                if classnum not in classtotypes.keys():
                    classtotypes[classnum]=[]
                if typenum not in classtotypes[classnum]:
                    classtotypes[classnum].append(typenum)
            if string in line and '#' not in line and 'BUFF' not in line and 'scale' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line and 'unit' not in line and 'ALLINGER' not in line and 'cutoff' not in line:
                if string=='vdw':
                    classnums=[int(linesplit[1])]
                elif string=='bond':
                    classnums=[int(linesplit[1]),int(linesplit[2])]
                elif string=='angle':
                    classnums=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
                elif string=='opbend' or string=='torsion':
                    classnums=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
                typenums=[]
                for classnum in classnums:
                    if classnum in classtotypes.keys():
                        typenum=classtotypes[classnum]
                        typenums.append(typenum)
                    else:
                        continue
                

                if len(typenums)==0:
                    continue
                cont=False
                for typenumls in typenums:
                    for typenum in typenumls:
                        if typenum not in prmtypes:
                            cont=True
                if cont==True:
                    continue
                combs=list(itertools.product(*typenums))
                for comb in combs:
                    typeindices=tuple(comb)
                    typeindicestoprmline[typeindices]=line



    return typeindicestoprmline

def TypeIndicesToIndices(poltype,typenums,idxtosymclass):
    allindices=[]
    for typenum in typenums: 
        indexes=GrabKeysFromValue(poltype,idxtosymclass,typenum)
        allindices.append(indexes)
    return allindices


def CompareEnergies(poltype,endtypeindicestoenergylist,bgntypeindicestoenergylist,bgnstatetypeindextoendstatetypeindex,energycutoff,stringlist,endtypetoprmlinelist,bgntypetoprmlinelist,bgnstateindextotypeindex,endstateindextotypeindex,bgntypeindicestocommentslist,endtypeindicestocommentslist,bgntypeprmtoamoeba09matchlist,endtypeprmtoamoeba09matchlist):
    
    for i in range(len(endtypeindicestoenergylist)):
        endtypeindicestoenergy=endtypeindicestoenergylist[i]
        bgntypeindicestoenergy=bgntypeindicestoenergylist[i]
        string=stringlist[i]
        bgntypeindicestoprmline=bgntypetoprmlinelist[i]
        endtypeindicestoprmline=endtypetoprmlinelist[i]
        bgntypeindicestocomments=bgntypeindicestocommentslist[i]
        endtypeindicestocomments=endtypeindicestocommentslist[i]
        endtypeormtoamoeba09match=endtypeprmtoamoeba09matchlist[i]
        bgntypeormtoamoeba09match=bgntypeprmtoamoeba09matchlist[i]
        for bgntypeindices,bgnenergy in bgntypeindicestoenergy.items():
            cutoff=energycutoff
            endtypeindices=[]
            if bgntypeindices in bgntypeindicestoprmline.keys():
                bgnprmline=bgntypeindicestoprmline[bgntypeindices]
            elif bgntypeindices[::-1] in bgntypeindicestoprmline.keys():
                bgnprmline=bgntypeindicestoprmline[bgntypeindices[::-1]]

            for bgntypeindex in bgntypeindices:
                endtypeindex=bgnstatetypeindextoendstatetypeindex[bgntypeindex]
                endtypeindices.append(endtypeindex)

            endtypeindices=tuple(endtypeindices)
            if endtypeindices in endtypeindicestoprmline.keys():
                endprmline=endtypeindicestoprmline[endtypeindices]
            elif endtypeindices[::-1] in endtypeindicestoprmline.keys():
                endprmline=endtypeindicestoprmline[endtypeindices[::-1]]
            if endtypeindices in endtypeindicestoenergy.keys():
                endenergy=endtypeindicestoenergy[endtypeindices]
            elif endtypeindices[::-1] in endtypeindicestoenergy.keys():
                endenergy=endtypeindicestoenergy[endtypeindices[::-1]]
            bgnindices=TypeIndicesToIndices(poltype,bgntypeindices,bgnstateindextotypeindex)
            endindices=TypeIndicesToIndices(poltype,endtypeindices,endstateindextotypeindex)
            bgncomments=[]
            if bgntypeindices in bgntypeindicestocomments.keys(): 
                bgncomments=bgntypeindicestocomments[bgntypeindices]
            endcomments=[]
            if endtypeindices in endtypeindicestocomments.keys():
                endcomments=endtypeindicestocomments[endtypeindices]
            if bgntypeindices in bgntypeormtoamoeba09match.keys():
                bgnmatches=bgntypeormtoamoeba09match[bgntypeindices]
            else:
                bgnmatches=[]
            if endtypeindices in endtypeormtoamoeba09match.keys():
                endmatches=endtypeormtoamoeba09match[endtypeindices]
            else:
                endmatches=[]
            diff=np.abs(bgnenergy-endenergy)
            if diff>cutoff:
                print('Detected energy difference of '+str(diff)+' for interaction type '+string +' with type indices '+str(bgntypeindices)+' and '+str(endtypeindices)+', with tolerance of '+str(cutoff))
                print('Begin parameters '+bgnprmline.replace('\n','')+' '+str(bgnindices)+'\n')
                for line in bgncomments:
                    print(line)
                for ls in bgnmatches:
                    print('amoeba09 match' +str(ls))
                print('End parameters '+endprmline.replace('\n','')+' '+str(endindices)+'\n')
                for line in endcomments:
                    print(line)
                for ls in endmatches:
                    print('amoeba09 match' +str(ls))


def ConvertIndicesToTypes(poltype,indicestoenergy,indextotypeindex):
    typeindicestoenergy={}
    for indices,energy in indicestoenergy.items():
        typeindices=[]
        for i in range(len(indices)):
            index=indices[i]
            typenum=indextotypeindex[index]
            typeindices.append(typenum)
        typeindices=tuple(typeindices)
        if typeindices not in typeindicestoenergy.keys():
             typeindicestoenergy[typeindices]=[]
        typeindicestoenergy[typeindices].append(energy)
    for typeindices,energyls in typeindicestoenergy.items():
        energy=sum(energyls)/len(energyls)
        typeindicestoenergy[typeindices]=energy


    return typeindicestoenergy


def GrabTinkerEnergy(poltype,alzfname):
    tot_energy = None
    if os.path.isfile(alzfname):  
        tmpfh = open(alzfname, 'r')
        for line in tmpfh:
            m = re.search(r'Potential Energy :\s+(\-*\d+\.\d+)',line)
            if not m is None:
                tot_energy = float(m.group(1))
        tmpfh.close()
    return tot_energy


def GrabOutput(poltype,alzout):
    
    vdwindicestoenergy={}
    bondindicestoenergy={}
    angleindicestoenergy={}
    opbendindicestoenergy={}
    torsionindicestoenergy={}

    temp=open(alzout,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()

        if 'Additional' in line or 'Classes' in line or 'Individual' in line or 'Atom' in line or 'Stretch' in line or 'Bending' in line or 'Torsional' in line:
            continue
        if 'VDW' in line:
            indexlist=[1,2]
            indices,energy=GrabIndicesAndEnergy(poltype,linesplit,indexlist)
            vdwindicestoenergy[indices]=energy
        if 'Bond' in line:
            indexlist=[1,2]
            indices,energy=GrabIndicesAndEnergy(poltype,linesplit,indexlist)
            bondindicestoenergy[indices]=energy
        if 'Angle' in line:
            indexlist=[1,2,3]
            indices,energy=GrabIndicesAndEnergy(poltype,linesplit,indexlist)
            angleindicestoenergy[indices]=energy
        if 'O-P-Bend' in line:
            indexlist=[1,2,3,4]
            indices,energy=GrabIndicesAndEnergy(poltype,linesplit,indexlist)
            opbendindicestoenergy[indices]=energy
        if 'Torsion' in line and 'Pi' not in line:
            indexlist=[1,2,3,4]
            indices,energy=GrabIndicesAndEnergy(poltype,linesplit,indexlist)
            torsionindicestoenergy[indices]=energy


    return vdwindicestoenergy,bondindicestoenergy,angleindicestoenergy,opbendindicestoenergy,torsionindicestoenergy

def GrabIndicesAndEnergy(poltype,linesplit,indexlist):
    indices=[]
    for index in indexlist:
        indexele=linesplit[index] 
        indexelesplit=indexele.split('-')
        theindex=int(indexelesplit[0])
        indices.append(theindex)
    energy=float(linesplit[-1])
    return tuple(indices),energy


def SanitizeKey(poltype,keyfilename):
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    newresults=[]
    oldkey=False
    for line in results:
        linesplit=line.split()
        if 'polarize' in line: 
            if '.39' not in line:
                oldkey=True
                if len(linesplit)==3:
                    newlinesplit=linesplit[0:2+1]+['.39']
                else:
                    newlinesplit=linesplit[0:2+1]+['.39']+linesplit[3:]
                line=' '.join(newlinesplit)+'\n'
        elif 'opbend' in line and 'unit' not in line and 'ALL' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line:
            if ' 0 ' not in line:
                oldkey=True
                prm=str(float(linesplit[3])*71.94)
                newlinesplit=[linesplit[0]]+[linesplit[2]]+[linesplit[1]]+['0']+['0']+[prm]
                line=' '.join(newlinesplit)+'\n'
        newresults.append(line)
    if oldkey==True:
        newline='allinger'+'\n'
        newresults.insert(1,newline)
    tempkey=keyfilename.replace('.key','_TEMP.key')
    temp=open(tempkey,'w')
    for line in newresults:
        temp.write(line)
    temp.close()
    os.remove(keyfilename)
    os.rename(tempkey,keyfilename)  
    keymods.RemoveKeyWord(poltype,keyfilename,'parameters')
    string='parameters '+poltype.prmfilepath+'\n'
    keymods.AddKeyWord(poltype,keyfilename,string)


def CompareBgnEndParameters(poltype):
    SanitizeKey(poltype,poltype.bgnstatekey)
    SanitizeKey(poltype,poltype.endstatekey)
    bgnstatetypeindextoendstatetypeindex,bgnstateindextotypeindex,endstateindextotypeindex,bgnindextoendcoords,typestocopy,endindextoneighbs=BgnTypeToNewType(poltype,poltype.bgnstatexyz,poltype.endstatexyz)

    newbgnxyzfile,newbgnkeyfile=UpdateXYZAndKeyFiles(poltype,bgnstatetypeindextoendstatetypeindex,poltype.bgnstatexyz,poltype.bgnstatekey,bgnindextoendcoords)
    
    CompareParameters(poltype,poltype.endstatexyz,poltype.endstatekey,newbgnxyzfile,newbgnkeyfile,bgnstatetypeindextoendstatetypeindex,bgnstateindextotypeindex,endstateindextotypeindex,poltype.analyzepath,poltype.energycutoff,poltype.bgnstatekey,poltype.potentialpath,poltype.poleditpath,endindextoneighbs)


def DihedralAngleScan(poltype,xyzfilename,keyfilename,torsion,mol,analyzepath):
    rdkittorsion=[i-1 for i in torsion]
    a,b,c,d=rdkittorsion[:]
    dihedral_energies = []
    dihedral_degrees = [i for i in range(0, 360, 30)]
    prevtempname=xyzfilename
    conformer = mol.GetConformer(0)
    torang = rdMolTransforms.GetDihedralDeg(conformer, a, b, c, d)
    for deg in dihedral_degrees:
        angle=round((torang+deg)%360)
        rdMolTransforms.SetDihedralDeg(conformer, a, b, c, d, angle)
        indextocoordinates={}
        for i in range(len(mol.GetAtoms())):
            pos = conformer.GetAtomPosition(i) 
            vec=np.array([float(pos.x),float(pos.y),float(pos.z)])
            indextocoordinates[i+1]=vec
        tempname=xyzfilename.replace('.xyz','_'+str(a)+'-'+str(b)+'-'+str(c)+'-'+str(d)+'_'+str(deg)+'.xyz')
        ReplaceXYZCoords(poltype,prevtempname,indextocoordinates,tempname,replace=False)
        poltype.ConvertTinkerXYZToCartesianXYZ(tempname)
        alzout=tempname.replace('.xyz','.alz')
        poltype.CallAnalyze(tempname,keyfilename,alzout,analyzepath,'d')
        prevtempname=tempname
        e=GrabTinkerEnergy(poltype,alzout)
        dihedral_energies.append(e)
    minvalue=min(dihedral_energies)
    dihedral_energies=[i - minvalue for i in dihedral_energies]

    return dihedral_energies,dihedral_degrees
