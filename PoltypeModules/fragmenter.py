import electrostaticpotential as esp
import torsiongenerator as torgen
import torsionfit
import os
import numpy
import openbabel
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import shutil
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
import svgutils.transform as sg
from cairosvg import svg2png
import matplotlib.pyplot as plt
from os.path import dirname, abspath
from itertools import combinations
import json
from collections import Counter
from rdkit.Geometry import Point3D
import sys # used for terminaing job after fragmenter finishes and troubleshooting
import symmetry as symm
import shlex


def AssignTotalCharge(poltype,molecule,babelmolecule):
    atomicnumtoformalchg={1:{2:1},5:{4:1},6:{3:-1},7:{2:-1,4:1},8:{1:-1,3:1},15:{4:1},16:{1:-1,3:1,5:-1},17:{0:-1,4:3},9:{0:-1},35:{0:-1},53:{0:-1}}
    totchg=Chem.rdmolops.GetFormalCharge(molecule)
    totchg=0
    for atom in molecule.GetAtoms():
        atomidx=atom.GetIdx()
        atomnum=atom.GetAtomicNum()
        val=atom.GetExplicitValence()
        valtochg=atomicnumtoformalchg[atomnum]
        if val not in valtochg.keys(): # then assume chg=0
            chg=0
        else:
            chg=valtochg[val]

        polneighb=False
        if atomnum==6:
            for natom in atom.GetNeighbors():
                natomicnum=natom.GetAtomicNum()
                if natomicnum==7 or natomicnum==8 or natomicnum==16:
                    polneighb=True
            if polneighb and val==3:
                chg=1
        totchg+=chg
        atom.SetFormalCharge(chg)
    return molecule


def GrabKeysFromValue(poltype,dic,thevalue):
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist


def GrabVdwAndTorsionParametersFromFragments(poltype,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor):
    valenceprmlist=[]
    parentsymmtorlist=[]
    allparenttortorskeys=[]
    for array in equivalentrotbndindexarrays:
        vdwfragment=False
        for rotbndindex in array:
            if '_' not in rotbndindex:
                vdwfragment=True
            if vdwfragment==False:
                rotkey=rotbndindex.replace('_',' ')
                tors,maintortors,tortor=GrabParentTorsions(poltype,rotbndindextoringtor,rotbndindex,rotkey)
                if len(maintortors)>0:
                    firsttor=maintortors[0]
                    secondtor=maintortors[1]
                    tortorclskey,tortoratomidxs=torsionfit.GenerateTorTorClasskey(poltype,firsttor,secondtor)
                    fwdsplit=tortorclskey.split()        
                    revsplit=fwdsplit[::-1]
                    rev='%d %d %d %d %d' % (int(revsplit[0]), int(revsplit[1]), int(revsplit[2]), int(revsplit[3]),int(revsplit[4]))
                    allparenttortorskeys.append(tortorclskey)
                    allparenttortorskeys.append(rev)
                for torsion in tors:
                    fwd=torgen.get_class_key(poltype,torsion[0],torsion[1],torsion[2],torsion[3])
                    fwdsplit=fwd.split()        
                    revsplit=fwdsplit[::-1]
                    rev='%d %d %d %d' % (int(revsplit[0]), int(revsplit[1]), int(revsplit[2]), int(revsplit[3]))
                    if fwd not in parentsymmtorlist:
                        parentsymmtorlist.append(fwd)
                    if rev not in parentsymmtorlist:
                        parentsymmtorlist.append(rev)
    classkeytoparameters={}
    classkeytofragmentfilename={}
    classkeytotorsionindexescollected={}
    classkeytoatomindexescollected={}
    classkeytosmartscollected={}
    classkeytosmartsposarraycollected={}
    tortorclasskeytogridpts={}
    tortorclasskeytogridlinesarray={}
    tortorclasskeytotorsionindexescollected={}
    tortorclasskeytosmartscollected={}
    tortorclasskeytosmartsposarraycollected={}
    curdir=os.getcwd()
    for rotbndindex,fragmentfilepath in rotbndindextofragmentfilepath.items():
        path,filename=os.path.split(fragmentfilepath)
        os.chdir(path)
        vdwfragment=False
        if '_' not in rotbndindex:
            vdwfragment=True 
         
        filelist=os.listdir(os.getcwd())
        foundkey5=False
        foundkey=False

        for ff in filelist:
            if '.key' in ff:
                foundkey=True
            if '.key_5' in ff:
                foundkey5=True
                parentindextofragindex=json.load(open("parentindextofragindex.txt"))
                parentsymclasstofragsymclass=json.load(open("parentsymclasstofragsymclass.txt"))
                classkeytosmartsposarray=json.load(open("classkeytosmartsposarray.txt"))
                classkeytosmarts=json.load(open("classkeytosmarts.txt"))
                parentclasskeytofragclasskey=json.load(open("parentclasskeytofragclasskey.txt"))
                if vdwfragment==False:
                    classkeytotorsionindexes=json.load(open("classkeytotorsionindexes.txt"))
                    parenttortorclasskeytofragtortorclasskey=json.load(open("parenttortorclasskeytofragtortorclasskey.txt"))
                    tortorclasskeytosmartsposarray=json.load(open("tortorclasskeytosmartsposarray.txt"))
                    tortorclasskeytosmarts=json.load(open("tortorclasskeytosmarts.txt"))
                    tortorclasskeytotorsionindexes=json.load(open("tortorclasskeytotorsionindexes.txt"))

                    fragsymmtorlist=[]
                    for tor in parentsymmtorlist:
                        if tor in parentclasskeytofragclasskey.keys():
                            fragclasskey=parentclasskeytofragclasskey[tor]
                            fragsymmtorlist.append(fragclasskey)
                    fragsymmtortorlist=[]
                    for tortorclskey in allparenttortorskeys:
                        if tortorclskey in parenttortorclasskeytofragtortorclasskey.keys():
                            fragtortorclskey=parenttortorclasskeytofragtortorclasskey[tortorclskey]
                            fragsymmtortorlist.append(fragtortorclskey) 
                else:
                    classkeytoatomindexes=json.load(open("classkeytoatomindexes.txt"))

                temp=open(ff,'r')
                results=temp.readlines()
                temp.close()
                for lineidx in range(len(results)):
                    line=results[lineidx]
                    newline=line.strip()
                    linesplit=newline.split()
                    if 'torsion' in line and '#' not in line and vdwfragment==False:
                        typea=int(linesplit[1])
                        typeb=int(linesplit[2])
                        typec=int(linesplit[3])
                        typed=int(linesplit[4])
                        prms=linesplit[5:]
                        tor=[typea,typeb,typec,typed]
                        torkey='%d %d %d %d' % (typea, typeb, typec, typed)
                        revtorkey='%d %d %d %d' % (typed, typec, typeb, typea)
                        if torkey in fragsymmtorlist or revtorkey in fragsymmtorlist:
                            if torkey in parentclasskeytofragclasskey.values():
                                fwdclasskeys=GrabKeysFromValue(poltype,parentclasskeytofragclasskey,torkey)
                                revclasskeys=GrabKeysFromValue(poltype,parentclasskeytofragclasskey,revtorkey)
                                classkeys=fwdclasskeys+revclasskeys
                            elif revtorkey in parentclasskeytofragclasskey.values():
                                fwdclasskeys=GrabKeysFromValue(poltype,parentclasskeytofragclasskey,torkey)
                                revclasskeys=GrabKeysFromValue(poltype,parentclasskeytofragclasskey,revtorkey)
                                classkeys=fwdclasskeys+revclasskeys
                            for classkey in classkeys:
                                smartsposarray=classkeytosmartsposarray[classkey]
                                torsionindexes=classkeytotorsionindexes[classkey]
                                smarts=classkeytosmarts[classkey]
                                classkeytotorsionindexescollected[classkey]=torsionindexes
                                classkeytosmartscollected[classkey]=smarts
                                classkeytosmartsposarraycollected[classkey]=smartsposarray
                                classkeytoparameters[classkey]=prms
                                classkeytofragmentfilename[classkey]=filename
                    elif 'tortors' in line and vdwfragment==False:
                        typea=int(linesplit[1])
                        typeb=int(linesplit[2])
                        typec=int(linesplit[3])
                        typed=int(linesplit[4])
                        typee=int(linesplit[5])
                        gridpts=[int(linesplit[6]),int(linesplit[7])]
                        tortorkey='%d %d %d %d %d' % (typea, typeb, typec, typed, typee)
                        revtortorkey='%d %d %d %d %d' % (typee, typed, typec, typeb, typea)
                        gridlinesarray=ParseGridLines(poltype,results,lineidx)
                        if tortorkey in fragsymmtortorlist or revtortorkey in fragsymmtortorlist:
                            if tortorkey in parenttortorclasskeytofragtortorclasskey.values():
                                fwdclasskeys=GrabKeysFromValue(poltype,parenttortorclasskeytofragtortorclasskey,tortorkey)
                                revclasskeys=GrabKeysFromValue(poltype,parenttortorclasskeytofragtortorclasskey,revtortorkey)
                                classkeys=fwdclasskeys+revclasskeys
                            elif revtortorkey in parenttortorclasskeytofragtortorclasskey.values():
                                fwdclasskeys=GrabKeysFromValue(poltype,parenttortorclasskeytofragtortorclasskey,tortorkey)
                                revclasskeys=GrabKeysFromValue(poltype,parenttortorclasskeytofragtortorclasskey,revtortorkey)
                                classkeys=fwdclasskeys+revclasskeys

                            for classkey in classkeys:
                                smartsposarray=tortorclasskeytosmartsposarray[classkey]
                                torsionindexes=tortorclasskeytotorsionindexes[classkey]
                                smarts=tortorclasskeytosmarts[classkey]
                                tortorclasskeytotorsionindexescollected[classkey]=torsionindexes
                                tortorclasskeytosmartscollected[classkey]=smarts
                                tortorclasskeytosmartsposarraycollected[classkey]=smartsposarray
                                tortorclasskeytogridpts[classkey]=gridpts
                                tortorclasskeytogridlinesarray[classkey]=gridlinesarray
                                classkeytofragmentfilename[classkey]=filename
                    elif 'vdw' in line and '#' not in line and vdwfragment==True:
                        fragclasskey=linesplit[1]
                        fragsymclass=int(fragclasskey)
                        parentclasskeys=GrabKeysFromValue(poltype,parentsymclasstofragsymclass,fragsymclass)
                        parentsymclass=parentclasskeys[0]
                        classkey=str(parentsymclass)
                        prms=linesplit[2:]
                        smartsposarray=classkeytosmartsposarray[classkey]
                        atomindexes=classkeytoatomindexes[classkey]
                        smarts=classkeytosmarts[classkey]
                        classkeytoatomindexescollected[classkey]=atomindexes
                        classkeytosmartscollected[classkey]=smarts
                        classkeytosmartsposarraycollected[classkey]=smartsposarray
                        classkeytoparameters[classkey]=prms
                        classkeytofragmentfilename[classkey]=filename




        if foundkey5==False and foundkey==True:
            raise ValueError('Fragment job did not finish '+filename)
    os.chdir(curdir)
    temp=open(poltype.key4fname,'r')
    results=temp.readlines()
    temp.close()
    temp=open(poltype.key5fname,'w')
    for line in results:
        fitline="# Fitted from Fragment "
        linesplit=line.split()
        if 'torsion' in line and '#' not in line:
            typea=int(linesplit[1])
            typeb=int(linesplit[2])
            typec=int(linesplit[3])
            typed=int(linesplit[4])
            torkey='%d %d %d %d' % (typea, typeb, typec, typed)
            rev='%d %d %d %d' % (typed,typec,typeb,typea)
            if torkey in classkeytoparameters.keys():
                valenceprmlist=ConstructTorsionLineFromFragment(poltype,torkey,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytotorsionindexescollected,temp,valenceprmlist,fitline)

            elif rev in classkeytoparameters.keys():
                valenceprmlist=ConstructTorsionLineFromFragment(poltype,rev,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytotorsionindexescollected,temp,valenceprmlist,fitline)

            else:
                temp.write(line)
        elif 'vdw' in line and '#' not in line:
            classkey=linesplit[1]
            if classkey in classkeytoatomindexescollected.keys():
                valenceprmlist=ConstructVdwLineFromFragment(poltype,classkey,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytoatomindexescollected,temp,valenceprmlist,fitline)


            
        else:
            temp.write(line)
    for tortorkey in tortorclasskeytogridpts.keys():
        valenceprmlist=ConstructTorsionTorsionLineFromFragment(poltype,tortorkey,classkeytofragmentfilename,tortorclasskeytogridpts,tortorclasskeytosmartsposarraycollected,tortorclasskeytosmartscollected,tortorclasskeytotorsionindexescollected,temp,valenceprmlist,fitline,tortorclasskeytogridlinesarray)
        WriteGridPoints(poltype,tortorkey,tortorclasskeytogridlinesarray,temp)



    temp.close()
    WriteOutDatabaseLines(poltype,valenceprmlist)



def ConstructVdwLineFromFragment(poltype,key,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytoatomindexescollected,temp,valenceprmlist,fitline):
    filename=classkeytofragmentfilename[key]
    prms=classkeytoparameters[key]
    parameters=' '.join(prms)
    line='vdw '+key+' '+parameters+'\n'
    smartspos=classkeytosmartsposarraycollected[key]
    smarts=classkeytosmartscollected[key]
    atomindexes=classkeytoatomindexescollected[key]
    fitline+=' SMARTS '+smarts+' vdw atom indexes = '+atomindexes+' with smarts torsion indices '+smartspos+' from fragment '+filename+"\n"
    valencestring='vdw'+' % '+smarts+' % '+smartspos+' % '
    prmsplit=prms.split()
    for prm in prmsplit:
        valencestring+=prm+','
    valencestring=valencestring[:-1]
    valencestring+='\n'
    temp.write(fitline)
    temp.write('# '+valencestring)
    temp.write(line)
    valenceprmlist.append(valencestring)
    return valenceprmlist





def WriteOutDatabaseLines(poltype,valenceprmlist):
    newtemp=open(poltype.databaseprmfilename,'w')
    for line in valenceprmlist:
        newtemp.write(line)
    newtemp.close()


def WriteGridPoints(poltype,key,tortorclasskeytogridlinesarray,temp):
    gridarray=tortorclasskeytogridlinesarray[key]
    for line in gridarray:
        temp.write(line)
        
    


def ParseGridLines(poltype,results,lineidx):
    gridlinesarray=[]
    for lineindex in range(len(results)):
        line=results[lineindex]
        if lineindex>lineidx:
            linesplit=line.split()
            if len(linesplit)==3:
                gridlinesarray.append(line)
            else:
                break
        
    return gridlinesarray

def ConstructTorsionLineFromFragment(poltype,key,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytotorsionindexescollected,temp,valenceprmlist,fitline):
    filename=classkeytofragmentfilename[key]
    prms=classkeytoparameters[key]
    parameters=' '.join(prms)
    torline='torsion '+key+' '+parameters+'\n'
    smartspos=classkeytosmartsposarraycollected[key]
    smarts=classkeytosmartscollected[key]
    torsionindexes=classkeytotorsionindexescollected[key]
    fitline+=' SMARTS '+smarts+' torsion atom indexes = '+torsionindexes+' with smarts torsion indices '+smartspos+' from fragment '+filename+"\n"
    valencestring='torsion'+' % '+smarts+' % '+smartspos+' % '
    newprms=prms[0::3]
    folds=prms[2::3]
    folds=[int(i) for i in folds]
    foldtoprms=dict(zip(folds,newprms))
    for i in range(1,4):
        if i not in foldtoprms.keys():
            foldtoprms[i]=str(0)
    for fold in sorted(foldtoprms.keys()):
        prm=foldtoprms[fold]
        valencestring+=prm+','
    valencestring=valencestring[:-1]
    valencestring+='\n'
    temp.write(fitline)
    temp.write('# '+valencestring)
    temp.write(torline)
    valenceprmlist.append(valencestring)
    return valenceprmlist


def ConstructTorsionTorsionLineFromFragment(poltype,key,classkeytofragmentfilename,tortorclasskeytogridpts,tortorclasskeytosmartsposarraycollected,tortorclasskeytosmartscollected,tortorclasskeytotorsionindexescollected,temp,valenceprmlist,fitline,tortorclasskeytogridlinesarray):
    filename=classkeytofragmentfilename[key]
    gridpts=tortorclasskeytogridpts[key]
    gridptsarray=[str(i) for i in gridpts]
    gridpts=' '.join(gridptsarray)
    gridlinesarray=tortorclasskeytogridlinesarray[key]
    torline='tortors '+key+' '+gridpts+'\n'
    smartspos=tortorclasskeytosmartsposarraycollected[key]
    smarts=tortorclasskeytosmartscollected[key]
    torsionindexes=tortorclasskeytotorsionindexescollected[key]
    fitline+=' SMARTS '+smarts+' torsion atom indexes = '+torsionindexes+' with smarts torsion indices '+smartspos+' from fragment '+filename+"\n"
    gridptscommastr=','.join(gridptsarray)
    valencestring='tortors'+' % '+smarts+' % '+smartspos+' % '+gridptscommastr+' % '
    valencestring+='\n'
    temp.write(fitline)
    temp.write('# '+valencestring)
    temp.write(torline)
    for gridline in gridlinesarray:
        gridline=gridline.replace('\n','')
        valencestring+=gridline+','
    valencestring=valencestring[:-1]
    valenceprmlist.append(valencestring)
    return valenceprmlist




def GrabWBOMatrixGaussian(poltype,outputlog,mol):
    try:
        WBOmatrix=numpy.empty((mol.GetNumAtoms(),mol.GetNumAtoms()))
    except Exception:
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
        elif 'Atom' in line and juststartWBOmatrix:
            matcols=len(linesplit)-1
        elif 'Wiberg bond index, Totals by atom' in line and juststartWBOmatrix:
            return WBOmatrix
        elif line=='\n' and juststartWBOmatrix:
            if 'Wiberg bond index matrix' not in results[lineidx-1]:
                currentcolnum+=matcols
        elif juststartWBOmatrix and 'Atom' not in line and line!='\n' and '--' not in line:
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
    except Exception:
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
        elif 'Atomic Valences:' in line and juststartWBOmatrix:
            return WBOmatrix
        elif AllIntegers(poltype,line.split()) and juststartWBOmatrix and line!='\n':
            colrowindex=lineidx
        elif juststartWBOmatrix and 'Irrep:' not in line and line!='\n' and not AllIntegers(poltype,line.split()):
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

def FindEquivalentFragments(poltype,fragmentarray,namearray):
    equivalentnamesarray=[]
    equivalentnamesarrayset=[]
    smartsarray=[rdmolfiles.MolToSmarts(m) for m in fragmentarray]
    for smartidx in range(len(smartsarray)):
        refsmarts=smartsarray[smartidx]
        refname=namearray[smartidx]
        reffragment=fragmentarray[smartidx]
        refsmartsmol=Chem.MolFromSmarts(refsmarts)
        nametemp=[]
        nametemp.append(refname)
        for anothersmartidx in range(len(smartsarray)):
            if anothersmartidx!=smartidx:
                smarts=smartsarray[anothersmartidx]
                name=namearray[anothersmartidx]
                fragment=fragmentarray[anothersmartidx]
                smartsmol=Chem.MolFromSmarts(smarts)
                match = refsmartsmol.HasSubstructMatch(smartsmol)
                if match==True:
                    nametemp.append(name)
        if set(nametemp) not in equivalentnamesarrayset:
            equivalentnamesarrayset.append(set(nametemp))
            equivalentnamesarray.append(set(nametemp))
     # need unique way to always order the same way so dont redo QM if list order is different

    newequivalentnamesarray=[]
    for array in equivalentnamesarray:
        firstrotidx=[]
        secondrotidx=[]
        newarray=[]
        for name in array:
            namesplit=name.split('_')
            namesplit=[int(i) for i in namesplit]
            if len(namesplit)>2:
                for i in range(0,len(namesplit),2):
                    firstrotidx.append(namesplit[i])
                    secondrotidx.append(namesplit[i+1])
            else:
                firstrotidx.append(namesplit[0])
                secondrotidx.append(namesplit[1])
        firstrotidx=[str(i)+'_' for i in firstrotidx]
        secondrotidx=[str(i) for i in secondrotidx]  
        for i in range(len(firstrotidx)):
            firstrot=firstrotidx[i]
            secondrot=secondrotidx[i]
            rotbnd=firstrot+secondrot
            newarray.append(rotbnd)
        if len(namesplit)>2:
            #newarray=newarray[::-1]
            newarray=['_'.join(newarray)] 
        newequivalentnamesarray.append(newarray)
    return newequivalentnamesarray

def FindRotatableBond(poltype,fragmol,rotbndindextofragment,temp):
    for rotbndindex in rotbndindextofragment.keys():
        m=rotbndindextofragment[rotbndindex]
        if len(m.GetAtoms())==len(fragmol.GetAtoms()) and rotbndindex not in temp:
            return rotbndindex

def FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog):
    poltypeinput={'dontdovdwscan':poltype.dontdovdwscan,'refinenonaroringtors':poltype.refinenonaroringtors,'tortor':poltype.tortor,'maxgrowthcycles':poltype.maxgrowthcycles,'suppressdipoleerr':'True','optmethod':poltype.optmethod,'toroptmethod':poltype.toroptmethod,'espmethod':poltype.espmethod,'torspmethod':poltype.torspmethod,'dmamethod':poltype.dmamethod,'torspbasisset':poltype.torspbasisset,'espbasisset':poltype.espbasisset,'dmabasisset':poltype.dmabasisset,'toroptbasisset':poltype.toroptbasisset,'optbasisset':poltype.optbasisset,'bashrcpath':poltype.bashrcpath,'externalapi':poltype.externalapi,'use_gaus':poltype.use_gaus,'use_gausoptonly':poltype.use_gausoptonly,'isfragjob':True,'poltypepath':poltype.poltypepath,'structure':tail,'numproc':poltype.numproc,'maxmem':poltype.maxmem,'maxdisk':poltype.maxdisk,'printoutput':True}
    if strfragrotbndindexes!=None:
        poltypeinput['onlyrotbndslist']=strfragrotbndindexes
    inifilepath=poltype.WritePoltypeInitializationFile(poltypeinput)
    cmdstr='nohup'+' '+'python'+' '+poltype.poltypepath+r'/'+'poltype.py'+' '+'&'
    cmdstr='cd '+shlex.quote(os.getcwd())+' && '+cmdstr
    molecprefix =  os.path.splitext(tail)[0]
    logname = molecprefix+ "-poltype.log"
    listofjobs.append(cmdstr)
    logpath=os.getcwd()+r'/'+logname
    if os.path.isfile(logpath): # make sure to remove logfile if exists, dont want WaitForTermination to catch previous errors before job is resubmitted
        os.remove(logpath)
    jobtooutputlog[cmdstr]=logpath
    return listofjobs,jobtooutputlog,logpath

def SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog):
    if poltype.externalapi is not None:
        finishedjobs,errorjobs=poltype.CallJobsLocalHost(jobtooutputlog,True)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False)

    return finishedjobs,errorjobs



def SpawnPoltypeJobsForFragments(poltype,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor):
    parentdir=dirname(abspath(os.getcwd()))
    listofjobs=[]
    jobtooutputlog={}
    for array in equivalentrotbndindexarrays:
        strfragrotbndindexes=''
        strparentrotbndindexes=''
        fragrotbnds=[]
        vdwfragment=False
        for i in range(len(array)):
            rotbndindex=array[i]
            fragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]
            head,tail=os.path.split(fragmentfilepath)
            os.chdir(head)

            parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
            if i==0:
                equivalentrotbndindex=rotbndindex
                equivalentparentindextofragindex=parentindextofragindex
            else:
                parentindextofragindex=equivalentparentindextofragindex[rotbndindex]    
            
            MakeFileName(poltype,equivalentrotbndindex,'equivalentfragment.txt')
            if '_' not in rotbndindex:
                vdwfragment=True
            if vdwfragment==False:
                rotbndindexes=rotbndindex.split('_')
                parentrotbndindexes=[int(i) for i in rotbndindexes]
                rotbndindexes=[int(i)-1 for i in parentrotbndindexes]
                fragrotbndindexes=[parentindextofragindex[i]+1 for i in rotbndindexes]
                for j in range(0,len(fragrotbndindexes),2):
                    fragrotbnd=str(fragrotbndindexes[j])+' '+str(fragrotbndindexes[j+1])
                    if fragrotbnd not in fragrotbnds:
                        fragrotbnds.append(fragrotbnd)
                        strfragrotbndindexes+=str(fragrotbndindexes[j])+' '+str(fragrotbndindexes[j+1])+','
                
                    strparentrotbndindexes+=str(parentrotbndindexes[j])+' '+str(parentrotbndindexes[j+1])+','
        if vdwfragment==True:
            strfragrotbndindexes=strfragrotbndindexes[:-1]

        if rotbndindex in rotbndindextoringtor.keys() or vdwfragment==True:
            strfragrotbndindexes=None
        strparentrotbndindexes=strparentrotbndindexes[:-1]
        fragmol=rotbndindextofragment[equivalentrotbndindex]
        fragmentfilepath=rotbndindextofragmentfilepath[equivalentrotbndindex]
        obConversion = openbabel.OBConversion()
        fragbabelmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(fragmentfilepath)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(fragbabelmol, fragmentfilepath)
        fragidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,fragbabelmol) 

        head,tail=os.path.split(fragmentfilepath)
        os.chdir(head)
        if vdwfragment==False:
            MakeFileName(poltype,strparentrotbndindexes,'torsions.txt')
        parentindextofragindex=rotbndindextoparentindextofragindex[equivalentrotbndindex]

        tempdic={}
        for parentindex,fragindex in parentindextofragindex.items():
            if type(parentindex)==str:
                newtemp=parentindextofragindex[parentindex]
                tempdic[parentindex]={}
                for idx,value in newtemp.items():
                    tempdic[parentindex][idx+1]=value+1
            else: 
                tempdic[parentindex+1]=fragindex+1
        parentindextofragindex=tempdic
        parentsymclasstofragsymclass={}
        
        for parentindex,fragindex in parentindextofragindex.items():
            if type(parentindex)!=str:
                parentsymclass=poltype.idxtosymclass[parentindex]
                fragsymclass=fragidxtosymclass[fragindex]
                parentsymclasstofragsymclass[parentsymclass]=fragsymclass
            else:
                newdic=parentindextofragindex[parentindex]
                parentsymclasstofragsymclass[parentindex]={}
                for idx,fragidx in newdic.items():
                    parentsymclass=poltype.idxtosymclass[idx]
                    fragsymclass=fragidxtosymclass[fragidx]
                    parentsymclasstofragsymclass[parentindex][parentsymclass]=fragsymclass
        WriteDictionaryToFile(poltype,parentsymclasstofragsymclass,"parentsymclasstofragsymclass.txt")
        WriteDictionaryToFile(poltype,parentindextofragindex,"parentindextofragindex.txt")
        tempmol=mol_with_atom_index_removed(poltype,fragmol) 
        fragsmarts=rdmolfiles.MolToSmarts(tempmol)
        m=mol_with_atom_index(poltype,fragmol)
        fragsmirks=rdmolfiles.MolToSmarts(m)
        fragidxarray=GrabAtomOrder(poltype,fragsmirks)
        classkeytosmartsposarray={}
        classkeytosmarts={}
        classkeytotorsionindexes={}
        classkeytoatomindexes={}

        tortorclasskeytosmartsposarray={}
        tortorclasskeytosmarts={}
        tortorclasskeytotorsionindexes={}

        parentclasskeytofragclasskey={}
        parenttortorclasskeytofragtortorclasskey={}
        copied=parentindextofragindex.copy()
        anothercopy=parentsymclasstofragsymclass.copy()
        if vdwfragment==False:
            for rotbndindex in array:
                if rotbndindex!=equivalentrotbndindex:
                    parentindextofragindex=copied[rotbndindex]
                    parentsymclasstofragsymclass=anothercopy[rotbndindex]
                rotkey=rotbndindex.replace('_',' ')
                tors,maintortors,tortor=GrabParentTorsions(poltype,rotbndindextoringtor,rotbndindex,rotkey)

                for torsion in tors:
                    classkey=torgen.get_class_key(poltype,torsion[0], torsion[1], torsion[2], torsion[3])
                    fragclasskey=GenerateFragmentClassKey(poltype,classkey,parentsymclasstofragsymclass)
                    smilesposstring,fragtorstring=GenerateSMARTSPositionStringAndAtomIndices(poltype,torsion,parentindextofragindex,fragidxarray)
                    parentclasskeytofragclasskey[classkey]=fragclasskey
                    classkeytosmartsposarray[classkey]=smilesposstring
                    classkeytosmarts[classkey]=fragsmarts
                    classkeytotorsionindexes[classkey]=fragtorstring
            if tortor==True:
                firsttor=maintortors[0]
                secondtor=maintortors[1]
                tortorclskey,tortoratomidxs=torsionfit.GenerateTorTorClasskey(poltype,firsttor,secondtor)
                fragclasskey=GenerateFragmentClassKey(poltype,tortorclskey,parentsymclasstofragsymclass)
                smilesposstring,fragtorstring=GenerateSMARTSPositionStringAndAtomIndices(poltype,tortoratomidxs,parentindextofragindex,fragidxarray)

                parenttortorclasskeytofragtortorclasskey[tortorclskey]=fragclasskey
                tortorclasskeytosmartsposarray[tortorclskey]=smilesposstring
                tortorclasskeytosmarts[tortorclskey]=fragsmarts
                tortorclasskeytotorsionindexes[tortorclskey]=fragtorstring


            WriteDictionaryToFile(poltype,tortorclasskeytosmartsposarray,"tortorclasskeytosmartsposarray.txt")
            WriteDictionaryToFile(poltype,tortorclasskeytosmarts,"tortorclasskeytosmarts.txt")
            WriteDictionaryToFile(poltype,tortorclasskeytotorsionindexes,"tortorclasskeytotorsionindexes.txt")

            WriteDictionaryToFile(poltype,parenttortorclasskeytofragtortorclasskey,"parenttortorclasskeytofragtortorclasskey.txt")
            WriteDictionaryToFile(poltype,parentclasskeytofragclasskey,"parentclasskeytofragclasskey.txt")
            WriteDictionaryToFile(poltype,classkeytotorsionindexes,"classkeytotorsionindexes.txt")
        else:
            rotbndindex=array[0]
            vdwatomindex=int(rotbndindex)
            classkey=str(poltype.idxtosymclass[vdwatomindex])
            ls=[vdwatomindex] 
            smilesposstring,fragatomstring=GenerateSMARTSPositionStringAndAtomIndices(poltype,ls,parentindextofragindex,fragidxarray)
            classkeytosmartsposarray[classkey]=smilesposstring
            classkeytosmarts[classkey]=fragsmarts
            classkeytoatomindexes[classkey]=fragatomstring
            WriteDictionaryToFile(poltype,classkeytoatomindexes,"classkeytoatomindexes.txt")

        WriteDictionaryToFile(poltype,classkeytosmartsposarray,"classkeytosmartsposarray.txt")
        WriteDictionaryToFile(poltype,classkeytosmarts,"classkeytosmarts.txt")

        wholexyz=parentdir+r'/'+poltype.xyzoutfile
        wholemol=parentdir+r'/'+poltype.molstructfname
        

        parentatoms=poltype.rdkitmol.GetNumAtoms()


        listofjobs,jobtooutputlog,newlog=FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog)
    os.chdir(parentdir)
    finishedjobs,errorjobs=SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog)
    return equivalentrotbndindexarrays,rotbndindextoringtor


def GenerateFragmentClassKey(poltype,clskey,parentsymclasstofragsymclass):
    classkeysplit=clskey.split()
    classkeysplit=[int(i) for i in classkeysplit]
    fragclasskeysplit=[parentsymclasstofragsymclass[i] for i in classkeysplit]
    fragclasskeysplit=[str(i) for i in fragclasskeysplit]
    fragclasskey=' '.join(fragclasskeysplit)
    return fragclasskey

def GenerateSMARTSPositionStringAndAtomIndices(poltype,torsion,parentindextofragindex,fragidxarray):
    smilesposarray=[]
    fragtor=[]
    for index in torsion:
        fragindex=parentindextofragindex[index]
        fragtor.append(fragindex)
        fragidxarraypos=fragidxarray.index(fragindex)
        smilespos=fragidxarraypos+1
        smilesposarray.append(smilespos)
    smilesposarray=[str(i) for i in smilesposarray]
    smilesposstring=','.join(smilesposarray)
    fragtor=[str(i) for i in fragtor]
    fragtorstring=[str(i) for i in fragtor]
    fragtorstring=','.join(fragtorstring)
    return smilesposstring,fragtorstring


def GrabParentTorsions(poltype,rotbndindextoringtor,rotbndindex,rotkey):
    tors=[]
    tortor=False
    maintortors=[]
    if rotbndindex in rotbndindextoringtor.keys():
        torset=rotbndindextoringtor[rotbndindex]
        for tor in torset:
            tors.append(tor)
    elif rotkey in poltype.rotbndlist.keys():
        tors=list(poltype.rotbndlist[rotkey])
    else:
        tortor=True
        rotkeysplit=rotkey.split()
        rotkeys=[]
        maintortors=[]
        for j in range(0,len(rotkeysplit),2):
            curkey=str(rotkeysplit[j])+' '+str(rotkeysplit[j+1])
            rotkeys.append(curkey)
        for key in rotkeys: 
            keytors=list(poltype.rotbndlist[key])
            maintortors.append(keytors[0])
            tors.extend(keytors)
    return tors,maintortors,tortor


def CountUnderscores(poltype,string):
    count=0
    for e in string:
        if e=='_':
            count+=1
    return count

def MakeFileName(poltype,string,filename):
    temp=open(filename,'w')
    temp.write(string+'\n')
    temp.close()


def WriteDictionaryToFile(poltype,dictionary,filename):
    with open(filename,'w') as f: 
        json.dump(dictionary, f)



def GrabAtomOrder(poltype,smirks):
    atomorder=[]
    for i in range(len(smirks)):
        e=smirks[i]
        prevchar=smirks[i-1]
        try:
            nextchar=smirks[i+1]
        except Exception:
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

def GrabIndexToCoordinates(poltype,mol):
    indextocoordinates={}
    iteratom = openbabel.OBMolAtomIter(mol)
    for atom in iteratom:
        atomidx=atom.GetIdx()
        coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
        indextocoordinates[atomidx]=coords
    return indextocoordinates

def AddInputCoordinatesAsDefaultConformer(poltype,m,indextocoordinates):
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        x,y,z = indextocoordinates[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    return m


def GenerateFrag(poltype,molindexlist,mol):
    molindexlist=[i+1 for i in molindexlist]
    em = openbabel.OBMol()
    oldindextonewindex={}
    oldindextoformalcharge={}
    for i,idx in enumerate(molindexlist):
        oldatom=poltype.mol.GetAtom(idx)
        em.AddAtom(oldatom)
        oldindextonewindex[idx]=i+1
        formalcharge=oldatom.GetFormalCharge()
        oldindextoformalcharge[idx]=formalcharge
        spinmult=oldatom.GetSpinMultiplicity()
    for oldindex,formalcharge in oldindextoformalcharge.items():
        newindex=oldindextonewindex[oldindex]
        atm=em.GetAtom(newindex) 
        spinmult=atm.GetSpinMultiplicity()
        atm.SetFormalCharge(formalcharge)

    atomiter=openbabel.OBMolAtomIter(em)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        formalcharge=atom.GetFormalCharge()
    atomswithcutbonds=[]
    bonditer=openbabel.OBMolBondIter(poltype.mol)
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
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
        bondorder=bond.GetBondOrder()
        diditwork=em.AddBond(bgnidx,endidx,bondorder)
        

    filename='frag.mol'
    WriteOBMolToMol(poltype,em,filename)
    indextocoordinates=GrabIndexToCoordinates(poltype,em) # need to convert indexes now
    nem=ReadToOBMol(poltype,filename)
    nem.AddHydrogens()

    hydindexes=[]
    atomiter=openbabel.OBMolAtomIter(nem)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        atomvec=[atom.GetX(),atom.GetY(),atom.GetZ()]
        if atomidx not in indextocoordinates.keys():
            indextocoordinates[atomidx]=atomvec
            hydindexes.append(atomidx)
    hydindexestokeep=[]
    for hydratedidx in atomswithcutbonds:
        atom=nem.GetAtom(hydratedidx)
        atomatomiter=openbabel.OBAtomAtomIter(atom)
        for natom in atomatomiter:
            natomidx=natom.GetIdx()
            if natomidx in hydindexes and natomidx not in hydindexestokeep: # then this one needs to be keeped
                hydindexestokeep.append(natomidx)
    hydindexestodelete=[]
    for hydidx in hydindexes:
        if hydidx not in hydindexestokeep:
            hydindexestodelete.append(hydidx)
    hydindexestodelete.sort(reverse=True)
    for hydidx in hydindexestodelete:
        atom=nem.GetAtom(hydidx)
        nem.DeleteAtom(atom)
        del indextocoordinates[hydidx]
    outputname='hydrated.mol'
    WriteOBMolToMol(poltype,nem,outputname)
    newmol=rdmolfiles.MolFromMolFile(outputname,removeHs=False)
    newmol.UpdatePropertyCache(strict=False)
    AllChem.EmbedMolecule(newmol)
    rdkitindextocoordinates={}
    for idx,coords in indextocoordinates.items():
        rdkitidx=idx-1
        rdkitindextocoordinates[rdkitidx]=coords
    newmol=AddInputCoordinatesAsDefaultConformer(poltype,newmol,rdkitindextocoordinates)

    rdkitoldindextonewindex={}
    for oldindex,newindex in oldindextonewindex.items():
        rdkitoldindex=oldindex-1
        rdkitnewindex=newindex-1
        rdkitoldindextonewindex[rdkitoldindex]=rdkitnewindex
    newmol=AssignTotalCharge(poltype,newmol,nem)
    return newmol,rdkitoldindextonewindex

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
    rdmolfiles.MolToMolFile(mol,outputname,kekulize=True)

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

def mol_with_atom_index_removed(poltype,mol):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        atom=mol.GetAtomWithIdx(idx)
        atom.ClearProp('molAtomMapNumber')
    return mol



def GenerateWBOMatrix(poltype,molecule,moleculebabel,structfname):
    error=False
    WBOmatrix=None
    curespmethod=poltype.espmethod
    curspbasisset=poltype.espbasisset
    poltype.espmethod='HF'
    poltype.espbasisset='MINIX'
    charge=Chem.rdmolops.GetFormalCharge(molecule)

    inputname,outputname=esp.CreatePsi4ESPInputFile(poltype,structfname,poltype.comespfname.replace('.com','_frag.com'),moleculebabel,poltype.maxdisk,poltype.maxmem,poltype.numproc,charge,False)
    finished,error=poltype.CheckNormalTermination(outputname)
    if not finished and not error:
        cmdstr='psi4 '+inputname+' '+outputname
        try:
             poltype.call_subsystem(cmdstr,True)
        except Exception:
             error=True
    if not error:
        WBOmatrix=GrabWBOMatrixPsi4(poltype,outputname,molecule)
    poltype.espmethod=curespmethod
    poltype.espbasisset=curspbasisset

    return WBOmatrix,outputname,error

def GenerateFragments(poltype,mol,torlist,parentWBOmatrix,missingvdwatomsets):

    newdir='Fragments'
    if not os.path.isdir(newdir):
        os.mkdir(newdir)
    os.chdir(newdir)
    fragspath=os.getcwd()
    rotbndindextoparentindextofragindex={}
    rotbndindextofragment={}
    rotbndindextofragmentfilepath={}
    rotbndindextofragWBOmatrix={}
    rotbndindextofragfoldername={}
    rotbndindextoWBOdifference={}
    rotbndindextoringtor={} 
    tempmaxgrowthcycles=poltype.maxgrowthcycles
    for torset in torlist:
        extendedtorindexes=[]
        for tor in torset:
            indexes=FirstPassAtomIndexes(poltype,tor)
            for index in indexes:
                if index not in extendedtorindexes:
                    extendedtorindexes.append(index)
        if torset in poltype.nonaroringtorsets or torset in missingvdwatomsets:
            poltype.maxgrowthcycles=0
        else:
            poltype.maxgrowthcycles=tempmaxgrowthcycles
        
        if torset in missingvdwatomsets:
            vdwfrag=True
        else:
            vdwfrag=False 

    
        WBOdifferencetofragWBOmatrix={}
        WBOdifferencetofoldername={}
        WBOdifferencetofragmol={}
        WBOdifferencetostructfname={}
        highlightbonds=[]
        fragfoldername=''
        for tor in torset:
            if len(tor)>1:
                fragfoldername+=str(tor[1])+'_'+str(tor[2])+'_'
            else:
                fragfoldername+=str(tor[0])+'_' # special case for vdw atom types
        fragfoldername+='Hydrated'
        if not os.path.isdir(fragfoldername):
            os.mkdir(fragfoldername)
        os.chdir(fragfoldername)
        fragmol,parentindextofragindex=GenerateFrag(poltype,extendedtorindexes,mol)
        growfragments=[]
        filename=fragfoldername+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        os.chdir('..')
        fragmoltoWBOmatrices={}
        fragmoltofragfoldername={}
        fragmoltobondindexlist={}
        fragfoldername=''
        for tor in torset:
            if len(tor)>1:
                fragfoldername+=str(tor[1])+'_'+str(tor[2])+'_'
            else:
                fragfoldername+=str(tor[0])+'_'
        fragfoldername+='Index'+'_'+str(0)
        if not os.path.isdir(fragfoldername):
            os.mkdir(fragfoldername)
        os.chdir(fragfoldername)
        rotbndidx=''
        for tor in torset:
            if len(tor)>1:
                rotbndidx+=str(tor[1])+'_'+str(tor[2])+'_'
            else:
                rotbndidx+=str(tor[0])+'_'
        rotbndidx=rotbndidx[:-1]
        if torset in poltype.nonaroringtorsets:
            rotbndindextoringtor[rotbndidx]=torset
        filename=fragfoldername+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        fragmoltofragfoldername[fragmol]=fragfoldername
        fragmolbabel=ReadMolFileToOBMol(poltype,filename)
        WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
        WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
        structfname=filename.replace('.mol','.sdf')
        fragWBOmatrix,outputname,error=GenerateWBOMatrix(poltype,fragmol,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
        if error:
            os.chdir('..')
            continue

        if torset in missingvdwatomsets:
            torset=GenerateFakeTorset(poltype,mol,parentindextofragindex)
        fragmentWBOvalues=numpy.array([round(fragWBOmatrix[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]],3) for tor in torset]) # rdkit is 0 index based so need to subtract 1, babel is 1 indexbased
        parentWBOvalues=numpy.array([round(parentWBOmatrix[tor[1]-1,tor[2]-1],3) for tor in torset]) # Matrix has 0,0 so need to subtract 1 from babel index
        WBOdifference=numpy.amax(numpy.abs(fragmentWBOvalues-parentWBOvalues))
        WBOdifferencetofragmol[WBOdifference]=fragmol
        WBOdifferencetostructfname[WBOdifference]=structfname
        rotbndindextoWBOdifference[rotbndidx]=WBOdifference
        fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,torset,fragmoltoWBOmatrices,fragmoltobondindexlist)

        os.chdir('..')

        WBOdifferencetofragWBOmatrix[WBOdifference]=fragWBOmatrix
        WBOdifferencetofoldername[WBOdifference]=fragfoldername
        WBOdifference=min(list(WBOdifferencetofragWBOmatrix))
        fragmol=WBOdifferencetofragmol[WBOdifference]
        structfname=WBOdifferencetostructfname[WBOdifference]
        fragWBOmatrix=WBOdifferencetofragWBOmatrix[WBOdifference]
        fragfoldername=WBOdifferencetofoldername[WBOdifference]
        rotbndindextofragfoldername[rotbndidx]=fragfoldername
        os.chdir(fragfoldername)
        for tor in torset:
            fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
            highlightbonds.append(fragrotbndidx)
        fragpath=os.getcwd()
        grow=False
        growfragments.append(fragmol)
        fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,torset,fragmoltoWBOmatrices,fragmoltobondindexlist)
        curdir=os.getcwd()
        os.chdir('..')
        growfragmoltoWBOmatrices=fragmoltoWBOmatrices.copy()
        growfragmoltofragfoldername=fragmoltofragfoldername.copy()
        growfragmoltobondindexlist=fragmoltobondindexlist.copy()

        fragments=[fragmol]
        Draw2DMoleculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,torset,'CombinationsWithIndex')
        sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in fragments]
        Draw2DMoleculesWithWBO(poltype,sanitizedfragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,torset,'CombinationsWithoutIndex')

        os.chdir(curdir)
        if WBOdifference<=poltype.WBOtol: # then we consider the fragment good enough to transfer torsion parameters, so make this fragment into .sdf file
            pass
        else:
            grow=True
            possiblefragatmidxs=GrowPossibleFragmentAtomIndexes(poltype,poltype.rdkitmol,extendedtorindexes)
            if len(possiblefragatmidxs)!=0 and poltype.maxgrowthcycles!=0:
                fragmol,newindexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex,fragpath,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist=GrowFragmentOut(poltype,mol,parentWBOmatrix,extendedtorindexes,WBOdifference,torset,fragfoldername,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,fragspath)
                fragmoltoWBOmatrices,fragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,torset,fragmoltoWBOmatrices,fragmoltobondindexlist)
        curdir=os.getcwd()
        os.chdir(fragspath)
        growfragments=[mol_with_atom_index(poltype,frag) for frag in growfragments]
        Draw2DMoleculesWithWBO(poltype,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,torset,'FragmentGrowthWithIndex')
        sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in growfragments]
        Draw2DMoleculesWithWBO(poltype,sanitizedfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,torset,'FragmentGrowthWithoutIndex')

        os.chdir(curdir)


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
    namearray=[]
    for rotbndindex in rotbndindextofragment.keys():
        if '_' in rotbndindex:
            fragment=rotbndindextofragment[rotbndindex]
            fragmentarray.append(fragment)
            namearray.append(rotbndindex)
    equivalentrotbndindexarrays=FindEquivalentFragments(poltype,fragmentarray,namearray)
    for rotbndindex in rotbndindextofragment.keys():
        if '_' not in rotbndindex: # for vdwfragments
            array=[rotbndindex]
            equivalentrotbndindexarrays.append(array)
 
    tempdic=rotbndindextoparentindextofragindex
    rotbndindextoparentindextofragindex=AddMatchingParentIndexesBetweenFragments(poltype,equivalentrotbndindexarrays,rotbndindextoparentindextofragindex,rotbndindextofragment)
    # now we need to redraw the 2Dimages for any fragments that are equivalent (get multiple torsions from different rotatable bonds around same fragment)
    '''
    curdir=os.getcwd()
    for rotbndindex in rotbndindextofragment.keys():
        if len(equivalentrotbndindexarrays)!=0:
            bndindextohighlightbonds={}
            for bndindexes in equivalentrotbndindexarrays:
                parenthighlightbonds=[]
     
                for bndindex in bndindexes:
                    parentindextofragindex=rotbndindextoparentindextofragindex[bndindex]
                    indexes=bndindex.split('_')
                    indexes=[int(i) for i in indexes]
                    parentrotbndidx=[indexes[0]-1,indexes[1]-1]

                    if parentrotbndidx not in parenthighlightbonds:
                        parenthighlightbonds.append(parentrotbndidx)
                for bndindex in bndindexes:
                    fraghighlightbonds=[]
                    parentindextofragindex=rotbndindextoparentindextofragindex[bndindex]
                    for parentrotbndidx in parenthighlightbonds:
                        fragrotbndidx=[]
                        for parentidx in parentrotbndidx:
                            if parentidx in parentindextofragindex.keys():
                                fragidx=parentindextofragindex[parentidx]
                                fragrotbndidx.append(fragidx)
                        if fragrotbndidx not in fraghighlightbonds and len(fragrotbndidx)==2:
                            fraghighlightbonds.append(fragrotbndidx)
                    bndindextohighlightbonds[bndindex]=fraghighlightbonds
                for bndindex in bndindexes:
                    highlightbonds=bndindextohighlightbonds[bndindex]
                    fragmol=rotbndindextofragment[bndindex]
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
                    Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',m,bondindexlist=highlightbonds,smirks=fragsmirks)
                    Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',m,bondindexlist=highlightbonds,smirks=fragsmirks)
            os.chdir(curdir)
    '''
    return rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor


def GenerateFakeTorset(poltype,mol,parentindextofragindex):
    bonditer=openbabel.OBMolBondIter(poltype.mol)
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        rdkitoendidx=oendidx-1
        rdkitobgnidx=obgnidx-1
        if rdkitoendidx in parentindextofragindex.keys() and rdkitobgnidx in parentindextofragindex.keys(): 
            tor=[1,oendidx,obgnidx,1]
            torset=[tuple(tor)]
            return torset



def AddMatchingParentIndexesBetweenFragments(poltype,equivalentrotbndindexarrays,rotbndindextoparentindextofragindex,rotbndindextofragment):
   newrotbndindextoparentindextofragindex={}
   for rotbndindexarray in equivalentrotbndindexarrays:
       for refrotbndindex in rotbndindexarray:
           reffragment=rotbndindextofragment[refrotbndindex]
           refparentindextofragindex=rotbndindextoparentindextofragindex[refrotbndindex]
           reffragindextoparentindex={value: key for key, value in refparentindextofragindex.items()} 
           newrefparentindextofragindex=refparentindextofragindex.copy()
            
           for rotbndindex in rotbndindexarray:
               if rotbndindex!=refrotbndindex:
                   if rotbndindex not in newrefparentindextofragindex.keys():
                       newrefparentindextofragindex[rotbndindex]={} 
                   parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
                   fragment=rotbndindextofragment[rotbndindex]
                   match = reffragment.HasSubstructMatch(fragment)
                   if match==True:
                       matches=reffragment.GetSubstructMatches(fragment)
                       for match in matches:
                           rotbndkey=rotbndindex.replace('_',' ')
                           tors=poltype.rotbndlist[rotbndkey]
                           for tor in tors:
                               indexes=[int(j)-1 for j in tor]
                               fragindexes=[parentindextofragindex[j] for j in indexes]
                               equivalentfragindexes=[match[j] for j in fragindexes]
                               equivalentrefindexes=[reffragindextoparentindex[j] for j in equivalentfragindexes]

                               for j in range(len(indexes)):
                                   index=indexes[j]
                                   equivalentfragindex=equivalentfragindexes[j]
                                   newrefparentindextofragindex[rotbndindex][index]=equivalentfragindex
           newrotbndindextoparentindextofragindex[refrotbndindex]=newrefparentindextofragindex
           
   return newrotbndindextoparentindextofragindex
 
  
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

def GrowFragmentOut(poltype,mol,parentWBOmatrix,indexes,WBOdifference,torset,fragfoldername,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,fragspath):
    fragfoldernamepath=os.getcwd()
    fragmentsforcomb=growfragments.copy()
    fragmentsforcombwbo=[WBOdifference]
    growthcycles=0
    WBOdiffarray=[]
    molarray=[]
    while not WBOdifference<=poltype.WBOtol:
        growthcycles+=1
        WBOdiffarray=[]
        molarray=[]
        fragmolidxtoparentindextofragindex={}
        fragmentidxtostructfname={}
        fragmolidxtofoldername={}
        fragmolidxtofragmol={}
        fragmolidxtofragmolbabel={}
        fragments=[]
        possiblefragatmidxs=GrowPossibleFragmentAtomIndexes(poltype,poltype.rdkitmol,indexes)
        if len(possiblefragatmidxs)!=0:
            for fragmolidx in range(len(possiblefragatmidxs)):
                indexlist=possiblefragatmidxs[fragmolidx]

                basename=fragfoldername+'_GrowFragment_'+str(fragmolidx)
                fragmol,parentindextofragindex=GenerateFrag(poltype,indexlist,mol)
                fragments.append(fragmol) # include the case where all H and no H converted to CH3
                if fragmol not in fragmentsforcomb:
                    fragmentsforcomb.append(fragmol)
                if not os.path .isdir(basename):
                    os.mkdir(basename)
                os.chdir(basename)
                growfragmoltofragfoldername[fragmol]=basename
                filename=basename+'.mol'
                WriteRdkitMolToMolFile(poltype,fragmol,filename)
                fragmolbabel=ReadMolFileToOBMol(poltype,filename)
                WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
                WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
                os.chdir('..')
                fragmolidxtofragmol[fragmolidx]=fragmol
                fragmolidxtofragmolbabel[fragmolidx]=fragmolbabel

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
                fragmolbabel=fragmolidxtofragmolbabel[fragmolidx]

                fragWBOmatrix,outputname,error=GenerateWBOMatrix(poltype,fragmol,fragmolbabel,structfname)
                if error:
                    os.chdir('..')
                    continue
                reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
                relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
                
                fragrotbndidxs=[[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]] for tor in torset]
                fragmentWBOvalues=numpy.array([round(fragWBOmatrix[fragrotbndidx[0],fragrotbndidx[1]],3) for fragrotbndidx in fragrotbndidxs])
                parentWBOvalues=numpy.array([round(parentWBOmatrix[tor[1]-1,tor[2]-1],3) for tor in torset])
                WBOdifference=numpy.amax(numpy.abs(fragmentWBOvalues-parentWBOvalues))
                fragmentsforcombwbo.append(WBOdifference)
                growfragmoltoWBOmatrices,growfragmoltobondindexlist=WriteOutFragmentInputs(poltype,fragmol,foldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,torset,growfragmoltoWBOmatrices,growfragmoltobondindexlist)

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
            WBOdifference=min(WBOdifftoindexlist.keys())
            parentindextofragindex=WBOdifftoparentindextofragindex[WBOdifference]
            indexes=WBOdifftoindexlist[WBOdifference]
            foldername=WBOdifftofolder[WBOdifference]
            structfname=WBOdifftostructfname[WBOdifference]
            RemoveTempFolders(poltype)
            os.chdir(foldername)

            fragmol=WBOdifftofragmol[WBOdifference]
            growfragments.append(fragmol)
            fragWBOmatrix=WBOdifftofragWBOmatrix[WBOdifference]
            fragpath=os.getcwd()
        else:
            break
        if poltype.maxgrowthcycles!=None:
            if growthcycles<=poltype.maxgrowthcycles:
                break


    curdir=os.getcwd()
    os.chdir('..')
    wbotofragments = dict(zip(fragmentsforcombwbo[1:], fragmentsforcomb[1:]))
    sortedwbotofragments={k: v for k, v in sorted(wbotofragments.items(), key=lambda item: item[0],reverse=True)}
    sorted_list=list(sortedwbotofragments.values())
    sorted_list.insert(0,fragmentsforcomb[0])
    Draw2DMoleculesWithWBO(poltype,sorted_list,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,torset,'CombinationsWithIndex')
    sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in sorted_list]
    Draw2DMoleculesWithWBO(poltype,sanitizedfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,torset,'CombinationsWithoutIndex')
    os.chdir(curdir)
    os.chdir(fragfoldernamepath)
    PlotFragmenterResults(poltype,WBOdiffarray,molarray)
    os.chdir(curdir)

    return fragmol,indexes,fragWBOmatrix,structfname,WBOdifference,parentindextofragindex,fragpath,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist


def GrowPossibleFragmentAtomIndexes(poltype,rdkitmol,indexes):
    possiblefragatmidxs=[]
    comblist=[]
    for bond in rdkitmol.GetBonds():
        aidx=bond.GetBeginAtomIdx()
        bidx=bond.GetEndAtomIdx()
        aatom=rdkitmol.GetAtomWithIdx(aidx)
        batom=rdkitmol.GetAtomWithIdx(bidx)
        aatomicnum=aatom.GetAtomicNum()
        batomicnum=batom.GetAtomicNum()
        bondorder=bond.GetBondTypeAsDouble()
        if bondorder>1:
            continue
        if aatomicnum==1 or batomicnum==1:
            continue
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
           newindexes=aromaticindexes
           for atmidx in newindexes:
               if atmidx not in indexlist:
                   indexlist.append(atmidx)
        temp=[]
        for idx in indexlist:
           neighbatom=poltype.rdkitmol.GetAtomWithIdx(idx)
           for neighbneighbatom in neighbatom.GetNeighbors():
               atomicnum=neighbneighbatom.GetAtomicNum()
               neighbneighbatomidx=neighbneighbatom.GetIdx()
               if atomicnum==1 and neighbneighbatomidx not in indexlist:
                   temp.append(neighbneighbatomidx)
               bond=poltype.rdkitmol.GetBondBetweenAtoms(neighbneighbatomidx,idx)
               bondorder=bond.GetBondTypeAsDouble()
               if bondorder>1 and neighbneighbatomidx not in indexlist:
                   temp.append(neighbneighbatomidx)
        for idx in temp:
            indexlist.append(idx)

        if indexlist not in possiblefragatmidxs:
           possiblefragatmidxs.append(indexlist)
    return possiblefragatmidxs


def WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,torset,fragmoltoWBOmatrices,fragmoltobondindexlist):
    highlightbonds=[]
    structfnamemol=fragfoldername+'.mol'
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('mol')
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,structfnamemol)
    tmpconv.SetOutFormat('sdf')
    structfname=fragfoldername+'.sdf'
    tmpconv.WriteFile(fragmolbabel,structfname)
    basename=fragfoldername+'_WBO_'+str(round(WBOdifference,3))
    for tor in torset:
        fragrotbndidx=[parentindextofragindex[tor[1]-1],parentindextofragindex[tor[2]-1]]
        highlightbonds.append(fragrotbndidx)
    reducedparentWBOmatrix=ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix)
    relativematrix=numpy.subtract(reducedparentWBOmatrix,fragWBOmatrix)
    m=mol_with_atom_index(poltype,fragmol)
    fragsmirks=rdmolfiles.MolToSmarts(m)
    Draw2DMoleculeWithWBO(poltype,fragWBOmatrix,basename+'_Absolute',m,bondindexlist=highlightbonds,smirks=fragsmirks)
    Draw2DMoleculeWithWBO(poltype,relativematrix,basename+'_Relative',m,bondindexlist=highlightbonds,smirks=fragsmirks)
    temp=[relativematrix,fragWBOmatrix]
    fragmoltoWBOmatrices[fragmol]=temp
    fragmoltobondindexlist[fragmol]=highlightbonds
    return fragmoltoWBOmatrices,fragmoltobondindexlist

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
                   if neighbatom.GetIsAromatic():
                       aromaticindexes=GrabAromaticAtoms(poltype,neighbatom)
                       newindexes=aromaticindexes
                       for atmidx in newindexes:
                           if atmidx not in molindexlist:
                               molindexlist.append(atmidx)
   temp=[]
   for index in molindexlist:
       atom=poltype.rdkitmol.GetAtomWithIdx(index)
       for neighbneighbatom in atom.GetNeighbors():
           atomicnum=neighbneighbatom.GetAtomicNum()
           neighbneighbatomidx=neighbneighbatom.GetIdx()
           if atomicnum==1 and neighbneighbatomidx not in molindexlist:
               temp.append(neighbneighbatomidx)
           bond=poltype.rdkitmol.GetBondBetweenAtoms(neighbneighbatomidx,index)
           bondorder=bond.GetBondTypeAsDouble()
           if bondorder>1 and neighbneighbatomidx not in molindexlist:
               temp.append(neighbneighbatomidx)
   for idx in temp:
       molindexlist.append(idx)
   return molindexlist

def Chunks(poltype,lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ChunksList(poltype,gen):
    newlst=[]
    for item in gen:
        newlst.append(item)
    return newlst


def Draw2DMoleculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,torset,basestr):

    bondlistoflists=[]
    for frag in fragments:
        bondindexlist=fragmoltobondindexlist[frag]
        bondlist=[]
        for bondindexes in bondindexlist:
            bond=frag.GetBondBetweenAtoms(bondindexes[0],bondindexes[1])
            bondidx=bond.GetIdx()
            bondlist.append(bondidx)
        bondlistoflists.append(bondlist)
    legendslist=[fragmoltofragfoldername[frag] for frag in fragments]
    molsperrow=3
    molsPerImage=molsperrow**2
    imagesize=400
    for i in range(len(fragments)):
        frag=fragments[i]
        rdDepictor.Compute2DCoords(frag)
    newfragments=[]
    if len(fragments)>1:
        firstmol=fragments[0]
        editmol=Chem.rdchem.EditableMol(firstmol)
        firstmolcopy=editmol.GetMol()
        newmol=mol_with_atom_index_removed(poltype,firstmolcopy)
        newermol = Chem.RemoveHs(newmol)
        smarts=rdmolfiles.MolToSmarts(newermol)
        smarts=smarts.replace('@','').replace('H3','').replace('H2','').replace('H','')
        patt = Chem.MolFromSmarts(smarts)
        newfragments.append(firstmol)

        for i in range(1,len(fragments)):
            frag=fragments[i]
            frag=mol_with_atom_index_removed(poltype,frag)
            overlap = frag.GetSubstructMatch(patt) # indexes of fragpatt corresponding to patt SMARTS but need the actual indexes of frag
            atomMap = [(paid,raid) for raid,paid in enumerate(overlap)]
            try:
                AllChem.AlignMol(frag,firstmol,atomMap=atomMap)
            except:
                return 
            newfragments.append(frag)
    fragmentchunks=ChunksList(poltype,Chunks(poltype,newfragments,molsPerImage))
    legendschunks=ChunksList(poltype,Chunks(poltype,legendslist,molsPerImage))
    bondlistoflistschunks=ChunksList(poltype,Chunks(poltype,bondlistoflists,molsPerImage))
    for i in range(len(fragmentchunks)):
        fragmentsublist=fragmentchunks[i]
        legendssublist=legendschunks[i]
        bondlistoflistssublist=bondlistoflistschunks[i]
        svg=Chem.Draw.MolsToGridImage(fragmentsublist,molsPerRow=molsperrow,subImgSize=(imagesize,imagesize),legends=legendssublist,highlightBondLists=bondlistoflistssublist,useSVG=True)
        fig = sg.fromstring(svg)
        ls=range(len(fragmentsublist))
        chunks=ChunksList(poltype,Chunks(poltype,ls,molsperrow))
        indextorow={}
        for rowidx in range(len(chunks)):
            row=chunks[rowidx]
            for j in row:
                indextorow[j]=rowidx
        for j in range(len(fragmentsublist)):
            frag=fragmentsublist[j]
            bondlist=bondlistoflistssublist[j]
            legend=legendssublist[j]
            drawer=rdMolDraw2D.MolDraw2DSVG(imagesize,imagesize)
            drawer.DrawMolecule(frag,highlightAtoms=[],highlightBonds=bondlist)

            atomidxtodrawcoords={}
            for bond in frag.GetBonds():
                bondidx=bond.GetIdx()
                if bondidx in bondlist:
                    begidx=bond.GetBeginAtomIdx()
                    endidx=bond.GetEndAtomIdx()
                    begatomdrawcoords=numpy.array(drawer.GetDrawCoords(begidx))
                    endatomdrawcoords=numpy.array(drawer.GetDrawCoords(endidx))
                    atomidxtodrawcoords[begidx]=begatomdrawcoords
                    atomidxtodrawcoords[endidx]=endatomdrawcoords

            WBOmatrixlist=fragmoltoWBOmatrices[frag]
            WBOmatrix=WBOmatrixlist[0]
            row=indextorow[j]
            x=(j-molsperrow*(row))*imagesize
            y=(row)*imagesize
            shift=numpy.array([x,y])
            for bond in frag.GetBonds():
                bondidx=bond.GetIdx()
                if bondidx in bondlist:
                    begidx=bond.GetBeginAtomIdx()
                    endidx=bond.GetEndAtomIdx()
                    begatomdrawcoords=atomidxtodrawcoords[begidx]+shift
                    endatomdrawcoords=atomidxtodrawcoords[endidx]+shift
                    bondcoords=(begatomdrawcoords+endatomdrawcoords)/2
                    WBOval=numpy.abs(WBOmatrix[begidx,endidx])
                    if WBOval==0:
                        continue
                    wbo=str(round(WBOval,4))
                    label = sg.TextElement(bondcoords[0],bondcoords[1], wbo, size=12, weight="bold")
                    array=endatomdrawcoords-begatomdrawcoords
                    if array[1]>=0:
                        pass
                    else:
                        array=-1*array
                    norm = numpy.linalg.norm(array)
                    normarray=array/norm
                    angle=numpy.abs(numpy.degrees(numpy.arccos(normarray[1])))
                    if angle>90:
                        angle=angle-90
                    if normarray[1]>=0 and normarray[0]>=0:
                        sign=-1
                    elif normarray[1]<=0 and normarray[0]<=0:
                        sign=-1
                    else:
                        sign=1
                    label.rotate(sign*angle,bondcoords[0],bondcoords[1])

                    fig.append(label)
     
        basename=basestr+'_'+'Bnd_'
        for tor in torset:
            basename+=str(tor[1])+'-'+str(tor[2])+'_'
        basename+='Index_'+str(i)
        fig.save(basename+'.svg')
        svg_code=fig.to_str()
        svg2png(bytestring=svg_code,write_to=basename+'.png')


def Draw2DMoleculeWithWBO(poltype,WBOmatrix,basename,mol,bondindexlist=None,smirks=None,imgsize=None):
    mol=mol_with_atom_index(poltype,mol)
    rdDepictor.Compute2DCoords(mol)
    if imgsize==None:
        drawer=rdMolDraw2D.MolDraw2DSVG(500,500)
    else:
        drawer=rdMolDraw2D.MolDraw2DSVG(imgsize,imgsize)
    bondlist=[]
    if bondindexlist is not None:
        for bondindexes in bondindexlist:
            bond=mol.GetBondBetweenAtoms(bondindexes[0],bondindexes[1])
            if bond!=None:
                bondidx=bond.GetIdx()
                bondlist.append(bondidx)
    try:
        mol.UpdatePropertyCache()
    except:
        return
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
            WBOval=numpy.abs(WBOmatrix[begidx,endidx])
            if WBOval==0:
                continue
            wbo=str(round(WBOval,4))
            label = sg.TextElement(bondcoords[0],bondcoords[1], wbo, size=12, weight="bold")
            array=endatomdrawcoords-begatomdrawcoords
            if array[1]>=0:
                pass
            else:
                array=-1*array
            norm = numpy.linalg.norm(array)
            normarray=array/norm
            angle=numpy.abs(numpy.degrees(numpy.arccos(normarray[1])))
            if angle>90:
                angle=angle-90
            if normarray[1]>=0 and normarray[0]>=0:
                sign=-1
            elif normarray[1]<=0 and normarray[0]<=0:
                sign=-1
            else:
                sign=1
            label.rotate(sign*angle,bondcoords[0],bondcoords[1])
            fig.append(label)
    if smirks is not None:
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
            if atm.GetIsAromatic() and atmindex not in aromaticindexes:
                aromaticindexes.append(atmindex)
            for natm in atm.GetNeighbors():
                if natm.GetIsAromatic() and natm.GetIdx() not in aromaticindexes:
                    aromaticindexes.append(natm.GetIdx())
        prevringidxlen=ringidxlen
        ringidxlen=len(aromaticindexes)

    return aromaticindexes


def PlotFragmenterResults(poltype,WBOdiffarray,molarray):
    fig=plt.figure()
    basename='NumberofAtomsVSWBODifference'
    plt.plot(WBOdiffarray,[m.GetNumAtoms() for m in molarray],'.')
    plt.xlabel('WBO Difference')
    plt.ylabel('Number of atoms in fragment')
    fig.savefig(basename+'.png')


