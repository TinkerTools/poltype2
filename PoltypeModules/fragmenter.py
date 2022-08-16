import electrostaticpotential as esp
import torsiongenerator as torgen
import torsionfit
import os
import numpy
from openbabel import openbabel
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
import itertools
import apicall as call
import math
import re

def AssignTotalCharge(poltype,molecule,babelmolecule):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist


def GrabVdwAndTorsionParametersFromFragments(poltype,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor,openkey,newkey,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    valenceprmlist={}
    parentsymmtorlist=[]
    allparenttortorskeys=[]
    curdir=os.getcwd()
    for array in equivalentrotbndindexarrays:
        vdwfragment=False
        firstfrag=array[0]
        if '_' not in firstfrag:
            vdwfragment=True
        if vdwfragment==False:
            parentrotbndindexes=rotbndindextoparentrotbndindexes[firstfrag]
            tors,maintortors,tortor,nonaroringfrag,rotbndindextotors=GrabParentTorsions(poltype,rotbndindextoringtor,array,parentrotbndindexes)
            if len(maintortors)>0:
                firsttor=maintortors[0]
                secondtor=maintortors[1]
                tortorclskey,tortoratomidxs=torsionfit.GenerateTorTorClasskey(poltype,firsttor,secondtor,poltype.idxtosymclass,poltype.rdkitmol)
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
    smartstovdwlinelist={}
    classkeytoparameters={}
    classkeytofragmentfilename={}
    classkeytotorsionindexescollected={}
    classkeytoatomindexescollected={}
    classkeytosmartscollected={}
    classkeytosmartsposarraycollected={}
    classkeytofitresults={}
    tortorclasskeytogridpts={}
    tortorclasskeytogridlinesarray={}
    tortorclasskeytotorsionindexescollected={}
    tortorclasskeytosmartscollected={}
    tortorclasskeytosmartsposarraycollected={}
    for array in equivalentrotbndindexarrays:
        rotbndindex=array[0]
        fragidxarray=rotbndindextosmartsindexarray[rotbndindex]
        fragindextosmartspos=GenerateAtomIndexToSMARTSPosition(poltype,fragidxarray)
        fragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]
        path,filename=os.path.split(fragmentfilepath)
        os.chdir(path)
        vdwfragment=False
        if '_' not in rotbndindex:
            vdwfragment=True 
        filelist=os.listdir(os.getcwd())
        foundkey5=False
        foundkey=False
        post=newkey.split('_')[-1]
        for ff in filelist:
            if '.key' in ff:
                foundkey=True
            if post in ff:
                foundkey5=True
                fragidxtosymclass=json.load(open("fragidxtosymclass.txt"))
                parentindextofragindex=json.load(open("parentindextofragindex.txt"))
                parentsymclasstofragsymclasses=json.load(open("parentsymclasstofragsymclasses.txt"))
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
                    if line.strip().startswith('torsion') and '#' not in line and 'Missing' not in line and vdwfragment==False:
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
                    elif 'RMSD(MM2,QM)' in line and vdwfragment==False:
                        typea=int(linesplit[2])
                        typeb=int(linesplit[3])
                        typec=int(linesplit[4])
                        typed=int(linesplit[5])
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
                                fragclasssplit=classkey.split()
                                linesplit[2]=fragclasssplit[0]
                                linesplit[3]=fragclasssplit[1]
                                linesplit[4]=fragclasssplit[2]
                                linesplit[5]=fragclasssplit[3]
                                classkeytofitresults[classkey]=' '.join(linesplit)+'\n'
                    
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
                    elif 'vdw' in line and '#' not in line:
                        for clskey,smrts in classkeytosmarts.items():
                            pass
                        linesplit=line.split() 
                        fragclasskey=linesplit[1]
                        for fragidx,symclass in fragidxtosymclass.items():
                            if symclass==int(fragclasskey):
                                break
                        fragsymclass=int(fragclasskey)
                        prms=linesplit[2:]
                        fragidx=int(fragidx)-1
                        if fragidx in fragindextosmartspos.keys():
                            smartspos=fragindextosmartspos[fragidx]
                            smilesposarray=[smartspos]
                            smilesposarray=[str(i) for i in smilesposarray]
                            smilespos=','.join(smilesposarray)
                            valencestring='vdw'+' % '+smrts+' % '+smilespos+' % '
                            for prm in prms:
                                valencestring+=prm+','
                            valencestring=valencestring[:-1]
                            valencestring+='\n'
                            if smrts not in smartstovdwlinelist.keys():
                                smartstovdwlinelist[smrts]=[]
                            smartstovdwlinelist[smrts].append(valencestring)
                           
                        for ls in parentsymclasstofragsymclasses.values():
                            first=ls[0]
                            if fragsymclass in ls:
                                parentclasskeys=GrabKeysFromValue(poltype,parentsymclasstofragsymclasses,ls)
                                for parentsymclass in parentclasskeys:
                                    classkey=str(parentsymclass)
                                    if vdwfragment==True:
                                        if classkey in classkeytosmartsposarray.keys():
                                            smartsposarray=classkeytosmartsposarray[classkey]
                                            atomindexes=classkeytoatomindexes[classkey]
                                            smarts=classkeytosmarts[classkey]
                                            classkeytoatomindexescollected[classkey]=atomindexes
                                            classkeytosmartscollected[classkey]=smarts
                                            classkeytosmartsposarraycollected[classkey]=smartsposarray
                                            classkeytoparameters[classkey]=prms
                                            classkeytofragmentfilename[classkey]=filename

                    elif 'RMSD(MM,QM)' in line and vdwfragment==True:
                        fragclasskey=linesplit[2]
                        fragsymclass=int(fragclasskey)
                        for ls in parentsymclasstofragsymclasses.values():
                            first=ls[0]
                            if fragsymclass in ls:
                                parentclasskeys=GrabKeysFromValue(poltype,parentsymclasstofragsymclasses,ls)
                                for parentsymclass in parentclasskeys:
                                    classkey=str(parentsymclass)
                                    if classkey in classkeytosmartsposarray.keys():
                                        fragclasssplit=classkey.split()
                                        linesplit[2]=fragclasssplit[0]
                                        classkeytofitresults[classkey]=' '.join(linesplit)+'\n'



        if foundkey5==False and foundkey==True:
            poltype.WriteToLog('Fragment job did not finish '+filename)
            raise ValueError('Fragment job did not finish '+filename)
    os.chdir(curdir)
    temp=open(openkey,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkey,'w')
    valkeytosmarts={}
    for line in results:
        fitline="# Fitted from Fragment "
        linesplit=line.split()
        if line.strip().startswith('torsion') and '#' not in line and 'Missing' not in line:
            typea=int(linesplit[1])
            typeb=int(linesplit[2])
            typec=int(linesplit[3])
            typed=int(linesplit[4])
            torkey='%d %d %d %d' % (typea, typeb, typec, typed)
            rev='%d %d %d %d' % (typed,typec,typeb,typea)
            if typeb<typec:
                valkey=tuple([typeb,typec])
            else:
                valkey=tuple([typec,typeb])
            if torkey in classkeytoparameters.keys():
                valenceprmlist,valkeytosmarts=ConstructTorsionLineFromFragment(poltype,torkey,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytotorsionindexescollected,temp,valenceprmlist,fitline,classkeytofitresults,valkey,valkeytosmarts)

            elif rev in classkeytoparameters.keys():
                valenceprmlist,valkeytosmarts=ConstructTorsionLineFromFragment(poltype,rev,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytotorsionindexescollected,temp,valenceprmlist,fitline,classkeytofitresults,valkey,valkeytosmarts)

            else:
                temp.write(line)
        elif 'vdw' in line and '#' not in line and 'Missing' not in line:
            classkey=linesplit[1]
            if classkey in classkeytoatomindexescollected.keys():
                valenceprmlist=ConstructVdwLineFromFragment(poltype,classkey,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytoatomindexescollected,temp,valenceprmlist,fitline,classkeytofitresults)

            else:
                temp.write(line)

            
        else:
            temp.write(line)
    for tortorkey in tortorclasskeytogridpts.keys():
        valenceprmlist=ConstructTorsionTorsionLineFromFragment(poltype,tortorkey,classkeytofragmentfilename,tortorclasskeytogridpts,tortorclasskeytosmartsposarraycollected,tortorclasskeytosmartscollected,tortorclasskeytotorsionindexescollected,temp,valenceprmlist,fitline,tortorclasskeytogridlinesarray)
        WriteGridPoints(poltype,tortorkey,tortorclasskeytogridlinesarray,temp)



    temp.close()
    WriteOutDatabaseLines(poltype,valenceprmlist,valkeytosmarts,smartstovdwlinelist)



def ConstructVdwLineFromFragment(poltype,key,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytoatomindexescollected,temp,valenceprmlist,fitline,classkeytofitresults,writeout=True):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    filename=classkeytofragmentfilename[key]
    prms=classkeytoparameters[key]
    parameters=' '.join(prms)
    line='vdw '+key+' '+parameters+'\n'
    smartspos=classkeytosmartsposarraycollected[key]
    smarts=classkeytosmartscollected[key]
    atomindexes=classkeytoatomindexescollected[key]
    if key in classkeytofitresults.keys():
        fitresultsline=classkeytofitresults[key]
        
    else:
        fitresultsline='' 

    fitline+=' SMARTS '+smarts+' vdw atom indexes = '+atomindexes+' with smarts torsion indices '+smartspos+' from fragment '+filename+"\n"
    valencestring='vdw'+' % '+smarts+' % '+smartspos+' % '
    for prm in prms:
        valencestring+=prm+','
    valencestring=valencestring[:-1]
    valencestring+='\n'
    if writeout==True:
        temp.write(fitline)
        if fitresultsline!='':
            temp.write(fitresultsline)
        temp.write('# '+valencestring)
        temp.write(line)
    if key not in valenceprmlist.keys():
        valenceprmlist[key]=[] 
    if valencestring not in valenceprmlist[key]: 
        if fitresultsline!='':
            valenceprmlist[key].append(fitresultsline)

        valenceprmlist[key].append(valencestring)
    return valenceprmlist





def WriteOutDatabaseLines(poltype,valenceprmlist,valkeytosmarts,smartstovdwlinelist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    newtemp=open(poltype.databaseprmfilename,'w')
    for key,ls in valenceprmlist.items():
        for line in ls:
            if '#' in line:
                 newtemp.write(line)
        for line in ls:
            if '#' not in line:
                 newtemp.write(line)
        if key in valkeytosmarts.keys():
            smarts=valkeytosmarts[key]
            vdwlinelist=smartstovdwlinelist[smarts]
            for line in vdwlinelist:
                newtemp.write(line) 
        newtemp.write('\n')
    newtemp.close()


def WriteGridPoints(poltype,key,tortorclasskeytogridlinesarray,temp):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gridarray=tortorclasskeytogridlinesarray[key]
    for line in gridarray:
        temp.write(line)
        
    


def ParseGridLines(poltype,results,lineidx):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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

def ConstructTorsionLineFromFragment(poltype,key,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarraycollected,classkeytosmartscollected,classkeytotorsionindexescollected,temp,valenceprmlist,fitline,classkeytofitresults,valkey,valkeytosmarts):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    filename=classkeytofragmentfilename[key]
    prms=classkeytoparameters[key]
    parameters=' '.join(prms)
    torline='torsion '+key+' '+parameters+'\n'
    smartspos=classkeytosmartsposarraycollected[key]
    smarts=classkeytosmartscollected[key]
    torsionindexes=classkeytotorsionindexescollected[key]
    if key in classkeytofitresults.keys():
        fitresultsline=classkeytofitresults[key]
    else:
        fitresultsline='' 
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
    if fitresultsline!='':
        temp.write(fitresultsline)
    temp.write('# '+valencestring)
    temp.write(torline)
    if valkey not in valenceprmlist.keys():
        valenceprmlist[valkey]=[] 
    if valkey not in valkeytosmarts.keys():
        valkeytosmarts[valkey]=smarts
    if valencestring not in valenceprmlist[valkey]:
        if fitresultsline!='':
            valenceprmlist[valkey].append(fitresultsline)

        valenceprmlist[valkey].append(valencestring)

    return valenceprmlist,valkeytosmarts


def ConstructTorsionTorsionLineFromFragment(poltype,key,classkeytofragmentfilename,tortorclasskeytogridpts,tortorclasskeytosmartsposarraycollected,tortorclasskeytosmartscollected,tortorclasskeytotorsionindexescollected,temp,valenceprmlist,fitline,tortorclasskeytogridlinesarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    if key not in valenceprmlist.keys():
        valenceprmlist[key]=[] 

    if valencestring not in valenceprmlist[key]: 
        valenceprmlist[key].append(valencestring)

    return valenceprmlist




def GrabWBOMatrixGaussian(poltype,outputlog,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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


def FindEquivalentFragments(poltype,fragmentarray,namearray,parentindextofragindexarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    equivalentnamesarray=[]
    equivalentnamesarrayset=[]
    nametoparentindextofragindex=dict(zip(namearray,parentindextofragindexarray))
    smartsarray=[rdmolfiles.MolToSmarts(m) for m in fragmentarray]
    for smartidx in range(len(smartsarray)):
        refsmarts=smartsarray[smartidx]
        refname=namearray[smartidx]
        reffragment=fragmentarray[smartidx]
        refsmartsmol=Chem.MolFromSmarts(refsmarts)
        refsmartsatoms=refsmartsmol.GetNumAtoms()

        nametemp=[]
        nametemp.append(refname)
        for anothersmartidx in range(len(smartsarray)):
            if anothersmartidx!=smartidx:
                smarts=smartsarray[anothersmartidx]
                name=namearray[anothersmartidx]
                fragment=fragmentarray[anothersmartidx]
                smartsmol=Chem.MolFromSmarts(smarts)
                smartsatoms=smartsmol.GetNumAtoms()
                match = refsmartsmol.HasSubstructMatch(smartsmol)
                if match==True and refsmartsatoms==smartsatoms:
                    nametemp.append(name)
        if set(nametemp) not in equivalentnamesarrayset:
            equivalentnamesarrayset.append(set(nametemp))
            equivalentnamesarray.append(set(nametemp))
    # need unique way to always order the same way so dont redo QM if list order is different
    newequivalentnamesarray=[]
    for array in equivalentnamesarray:
        Sumarray=[]
        for name in array:
            namesplit=name.split('_')
            namesplit=[int(i) for i in namesplit]
            Sum=sum(namesplit)
            Sumarray.append(Sum)
        sortedarray=[x for _, x in sorted(zip(Sumarray,array), key=lambda pair: pair[0])] 
        newequivalentnamesarray.append(sortedarray) 
    return newequivalentnamesarray


def FindRotatableBond(poltype,fragmol,rotbndindextofragment,temp):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for rotbndindex in rotbndindextofragment.keys():
        m=rotbndindextofragment[rotbndindex]
        if len(m.GetAtoms())==len(fragmol.GetAtoms()) and rotbndindex not in temp:
            return rotbndindex

def CopyAllQMDataAndRename(poltype,molecprefix,parentdir):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    os.chdir(parentdir)
    files=os.listdir()
    for f in files:
        if poltype.molecprefix in f:
            fsplit=f.split(poltype.molecprefix)
            secondpart=fsplit[1]
            newfname=molecprefix+secondpart 
            destination=curdir+r'/'+newfname
            
    os.chdir(curdir)    









def FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog,fragmol,parentdir,vdwfragment,strfragvdwatomindex,onlyfittorsions,jobtoinputfilepaths):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempmaxmem,tempmaxdisk,tempnumproc=poltype.PartitionResources()
    poltypeinput={'xtbtorresconstant':poltype.xtbtorresconstant,'deleteallnonqmfiles':poltype.deleteallnonqmfiles,'debugmode':poltype.debugmode,'atmidx':poltype.prmstartidx,'parentname':poltype.parentname,'use_gau_vdw':poltype.use_gau_vdw,'use_qmopt_vdw':poltype.use_qmopt_vdw,'onlyvdwatomindex':poltype.onlyvdwatomindex,'tordebugmode':poltype.tordebugmode,'dovdwscan':poltype.dovdwscan,'refinenonaroringtors':poltype.refinenonaroringtors,'tortor':poltype.tortor,'maxgrowthcycles':poltype.maxgrowthcycles,'suppressdipoleerr':'True','toroptmethod':poltype.toroptmethod,'espmethod':poltype.espmethod,'torspmethod':poltype.torspmethod,'dmamethod':poltype.dmamethod,'torspbasisset':poltype.torspbasisset,'espbasisset':poltype.espbasisset,'dmabasisset':poltype.dmabasisset,'toroptbasisset':poltype.toroptbasisset,'optbasisset':poltype.optbasisset,'bashrcpath':poltype.bashrcpath,'externalapi':poltype.externalapi,'use_gaus':poltype.use_gaus,'use_gausoptonly':poltype.use_gausoptonly,'isfragjob':True,'poltypepath':poltype.poltypepath,'structure':tail,'numproc':tempnumproc,'maxmem':tempmaxmem,'maxdisk':tempmaxdisk,'printoutput':True,'toroptmethodlist':','.join(poltype.toroptmethodlist),'torspmethodlist':','.join(poltype.torspmethodlist)}
    if strfragrotbndindexes!=None:
        poltypeinput['onlyrotbndslist']=strfragrotbndindexes
        rotbnds=strfragrotbndindexes.split(',')
        rotbnds.remove('') 
        if len(rotbnds)==1:
            poltypeinput['tortor']='False'
   
    if vdwfragment==True:
        poltypeinput['dontdotor']=True
        poltypeinput['onlyvdwatomindex']=strfragvdwatomindex
    if len(onlyfittorsions)!=0:
        string=''
        for tor in onlyfittorsions:
            stringtor=[str(i) for i in tor]
            stringtor=' '.join(stringtor)
            string+=stringtor+','
        string=string[:-1]
        poltypeinput['onlyfittorstogether']=string 

    inifilepath=poltype.WritePoltypeInitializationFile(poltypeinput)
    cmdstr='python'+' '+poltype.poltypepath+r'/'+'poltype.py'
    cmdstr='cd '+shlex.quote(os.getcwd())+' && '+cmdstr
    
    jobtoinputfilepaths[cmdstr]=inifilepath
    molecprefix =  os.path.splitext(tail)[0]
    logname = molecprefix+ "-poltype.log"
    listofjobs.append(cmdstr)
    logpath=os.getcwd()+r'/'+logname
    jobtooutputlog[cmdstr]=logpath
    if os.path.isfile(logpath): # make sure to remove logfile if exists, dont want WaitForTermination to catch previous errors before job is resubmitted
        os.remove(logpath)
    

    b = Chem.MolToSmiles(poltype.rdkitmol)
    a = Chem.MolToSmiles(fragmol)
    if a==b:
        CopyAllQMDataAndRename(poltype,molecprefix,parentdir)
    return listofjobs,jobtooutputlog,logpath,jobtoinputfilepaths


def SetupClusterSubmission(poltype,listofjobs,parentdir):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    jobtologlistfilenameprefix=os.path.join(parentdir,'FragmentJobToLog'+'_'+poltype.molecprefix)
    jobtooutputfiles={}
    jobtoabsolutebinpath={}
    scratchdir=poltype.scratchdir
    for cmdstr in listofjobs:
        jobtooutputfiles[cmdstr]=None
        jobtoabsolutebinpath[cmdstr]=None


    return jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix

def SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog('Submitting Fragment Jobs') 
    if poltype.fragmenterdebugmode==False:
        if poltype.externalapi is not None and poltype.fragmentjobslocal==False:
            call.CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
            finshedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog,False)

        else:
            finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,True)
    else:
        finishedjobs=list(jobtooutputlog.keys())
        errorjobs=[]

    return finishedjobs,errorjobs

def ChangeNumpyIntToIntDicKeys(poltype,dic):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newdic={}
    for key,value in dic.items():
        newdic[int(key)]=value
    return newdic


def SpawnPoltypeJobsForFragments(poltype,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog('Spawning Poltype Fragment Jobs')
    parentdir=dirname(abspath(os.getcwd()))
    listofjobs=[]
    jobtooutputlog={}
    jobtoinputfilepaths={}
    listofrotbndindexes=[]
    rotbndindextoparentrotbndindexes={}
    rotbndindextosmartsindexarray={}
    for arrayidx in range(len(equivalentrotbndindexarrays)):
        array=equivalentrotbndindexarrays[arrayidx]
        strfragrotbndindexes=''
        strparentrotbndindexes=''
        fragrotbnds=[]
        vdwfragment=False
        vdwparentindices=[]
        parentsymclasstofragsymclasses={}
        fragmol=rotbndindextofragment[array[0]]
        fragmentfilepath=rotbndindextofragmentfilepath[array[0]]

        obConversion = openbabel.OBConversion()
        fragbabelmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(fragmentfilepath)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(fragbabelmol, fragmentfilepath)
        fragidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,fragbabelmol)
        rdkitmol=rdmolfiles.MolFromMolFile(fragmentfilepath,removeHs=False)
        for i in range(len(array)):
            rotbndindex=array[i]
            if '_' not in rotbndindex:
                vdwfragment=True

            fragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]
            head,tail=os.path.split(fragmentfilepath)
            os.chdir(head)
            
            if i==0:
                equivalentrotbndindex=rotbndindex

            parentindextofragindex=rotbndindextoparentindextofragindex[equivalentrotbndindex]
            otherparentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
            otherfragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]
            otherrdkitmol=rdmolfiles.MolFromMolFile(otherfragmentfilepath,removeHs=False)
            matches=rdkitmol.GetSubstructMatches(otherrdkitmol)
            match=matches[0]
            indices=list(range(len(match)))
            othermolindextoequivmolindex=dict(zip(indices,match)) 
            if vdwfragment==False:
                MakeFileName(poltype,equivalentrotbndindex,'equivalentfragment.txt')
                rotbndindexes=rotbndindex.split('_')
                parentrotbndindexes=[int(j) for j in rotbndindexes]
                rotbndindexes=[int(j)-1 for j in parentrotbndindexes]
                rotbndindexes=ConvertSameBondTypes(poltype,rotbndindexes,parentindextofragindex,otherparentindextofragindex,othermolindextoequivmolindex)
                fragrotbndindexes=[parentindextofragindex[j] for j in rotbndindexes]
                fragrotbndindexes=[j+1 for j in fragrotbndindexes]
                for j in range(0,len(fragrotbndindexes),2):
                    fragrotbnd=str(fragrotbndindexes[j])+' '+str(fragrotbndindexes[j+1])
                    if fragrotbnd not in fragrotbnds:
                        fragrotbnds.append(fragrotbnd)
                        strfragrotbndindexes+=str(fragrotbndindexes[j])+' '+str(fragrotbndindexes[j+1])+','
                
                    strparentrotbndindexes+=str(parentrotbndindexes[j])+' '+str(parentrotbndindexes[j+1])+','
            else:
                rotbndindex=array[i]
                vdwatomindex=int(rotbndindex)
                vdwparentindices.append(vdwatomindex)
        if vdwfragment==True:
            strfragrotbndindexes=strfragrotbndindexes[:-1]

        if vdwfragment==True:
            strfragrotbndindexes=None

        strparentrotbndindexes=strparentrotbndindexes[:-1]
        rotbndindextoparentrotbndindexes[equivalentrotbndindex]=strparentrotbndindexes   
        parentindextofragindex=rotbndindextoparentindextofragindex[equivalentrotbndindex]
        fragindextoparentindex={v: k for k, v in parentindextofragindex.items()}
        for parentindex,fragindex in parentindextofragindex.items():
            parentsymclass=poltype.idxtosymclass[parentindex+1]
            fragsymclass=fragidxtosymclass[fragindex+1]
            if parentsymclass not in parentsymclasstofragsymclasses.keys():
                parentsymclasstofragsymclasses[parentsymclass]=[]
            if fragsymclass not in parentsymclasstofragsymclasses[parentsymclass]: 
                parentsymclasstofragsymclasses[parentsymclass].append(fragsymclass)
        fragmentfilepath=rotbndindextofragmentfilepath[equivalentrotbndindex]

        head,tail=os.path.split(fragmentfilepath)
        os.chdir(head)
        if vdwfragment==False:
            MakeFileName(poltype,strparentrotbndindexes,'torsions.txt')
    

        fragidxtosymclass=ChangeNumpyIntToIntDicKeys(poltype,fragidxtosymclass) 
        WriteDictionaryToFile(poltype,fragidxtosymclass,"fragidxtosymclass.txt")
        WriteDictionaryToFile(poltype,parentsymclasstofragsymclasses,"parentsymclasstofragsymclasses.txt")
        WriteDictionaryToFile(poltype,parentindextofragindex,"parentindextofragindex.txt")
        tempmol=mol_with_atom_index_removed(poltype,fragmol) 
        fragsmarts=rdmolfiles.MolToSmarts(tempmol)
        m=mol_with_atom_index(poltype,fragmol)
        fragsmirks=rdmolfiles.MolToSmarts(m)
        fragidxarray=GrabAtomOrder(poltype,fragsmirks)
        rotbndindextosmartsindexarray[equivalentrotbndindex]=fragidxarray
        classkeytosmartsposarray={}
        classkeytosmarts={}
        classkeytotorsionindexes={}
        classkeytoatomindexes={}
        tortorclasskeytosmartsposarray={}
        tortorclasskeytosmarts={}
        tortorclasskeytotorsionindexes={}
        parentclasskeytofragclasskey={}
        parenttortorclasskeytofragtortorclasskey={}
        strfragvdwatomindex=None
        fragindextoparentindex={v: k for k, v in parentindextofragindex.items()}
        onlyfittorsions=[]
        equivalentmolstruct=ReadToOBMol(poltype,fragmentfilepath)
        if vdwfragment==False:
            vdwkeypath=os.path.join(poltype.startdir,poltype.key5fname) 
            shutil.copy(vdwkeypath,'parentvdw.key')
            tors,maintortors,tortor,nonaroringfrag,rotbndindextotors=GrabParentTorsions(poltype,rotbndindextoringtor,array,strparentrotbndindexes)
            for rotbndindex,tors in rotbndindextotors.items():
                otherparentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
                otherfragmentfilepath=rotbndindextofragmentfilepath[rotbndindex]        
                molstruct=ReadToOBMol(poltype,otherfragmentfilepath)
                indextoreferenceindex=MatchOBMols(poltype,molstruct,equivalentmolstruct,otherparentindextofragindex,parentindextofragindex)              
                for torsion in tors:
                    rdkittorsion=[k-1 for k in torsion]
                    trueotherparentindextofragindex={}
                    for otherparentindex,fragindex in otherparentindextofragindex.items():
                        truefragindex=indextoreferenceindex[fragindex]
                        trueotherparentindextofragindex[otherparentindex]=truefragindex
                    fragtor=[trueotherparentindextofragindex[k] for k in rdkittorsion]
                    fragtorbabel=[k+1 for k in fragtor]
                    if nonaroringfrag==True:
                        if fragtorbabel not in onlyfittorsions:
                            onlyfittorsions.append(fragtorbabel)

                    fragclasskey=[fragidxtosymclass[k] for k in fragtorbabel]
                    fragclasskey=[str(k) for k in fragclasskey]
                    fragclasskey=' '.join(fragclasskey)
                    classkey=torgen.get_class_key(poltype,torsion[0],torsion[1],torsion[2],torsion[3])
                    smilesposstring,fragtorstring=GenerateSMARTSPositionStringAndAtomIndices(poltype,torsion,trueotherparentindextofragindex,fragidxarray)
                    parentclasskeytofragclasskey[classkey]=fragclasskey

                    classkeytosmartsposarray[classkey]=smilesposstring
                    classkeytosmarts[classkey]=fragsmarts
                    classkeytotorsionindexes[classkey]=fragtorstring
            if tortor==True:
                firsttor=maintortors[0]
                secondtor=maintortors[1]
                tortorclskey,tortoratomidxs=torsionfit.GenerateTorTorClasskey(poltype,firsttor,secondtor,poltype.idxtosymclass,poltype.rdkitmol)
                firstrdkittor=[k-1 for k in firsttor]
                secondrdkittor=[k-1 for k in secondtor]
                firstfragtor=[parentindextofragindex[k] for k in firstrdkittor]
                secondfragtor=[parentindextofragindex[k] for k in secondrdkittor]

                firstfragtorbabel=[k+1 for k in firstfragtor]
                secondfragtorbabel=[k+1 for k in secondfragtor]

                fragtortorclskey,fragtortoratomidxs=torsionfit.GenerateTorTorClasskey(poltype,firstfragtorbabel,secondfragtorbabel,fragidxtosymclass,fragmol)
                smilesposstring,fragtorstring=GenerateSMARTSPositionStringAndAtomIndices(poltype,tortoratomidxs,parentindextofragindex,fragidxarray)
                parenttortorclasskeytofragtortorclasskey[tortorclskey]=fragtortorclskey
                tortorclasskeytosmartsposarray[tortorclskey]=smilesposstring
                tortorclasskeytosmarts[tortorclskey]=fragsmarts
                tortorclasskeytotorsionindexes[tortorclskey]=fragtorstring

            WriteDictionaryToFile(poltype,tortorclasskeytosmartsposarray,"tortorclasskeytosmartsposarray.txt")
            WriteDictionaryToFile(poltype,tortorclasskeytosmarts,"tortorclasskeytosmarts.txt")
            WriteDictionaryToFile(poltype,tortorclasskeytotorsionindexes,"tortorclasskeytotorsionindexes.txt")

            WriteDictionaryToFile(poltype,parenttortorclasskeytofragtortorclasskey,"parenttortorclasskeytofragtortorclasskey.txt")
            WriteDictionaryToFile(poltype,classkeytotorsionindexes,"classkeytotorsionindexes.txt")
        else:
            vdwindextoequivindex={}
            equivvdwindex=vdwparentindices[0]
            for vdwatomindex in vdwparentindices:
                vdwindextoequivindex[vdwatomindex]=equivvdwindex
            for k in range(len(vdwparentindices)): 
                vdwatomindex=vdwparentindices[k]
                classkey=str(poltype.idxtosymclass[vdwatomindex])
                ls=[vdwindextoequivindex[vdwatomindex]] 
                smilesposstring,fragatomstring=GenerateSMARTSPositionStringAndAtomIndices(poltype,ls,parentindextofragindex,fragidxarray)
                fragvdwatomindex=parentindextofragindex[vdwindextoequivindex[vdwatomindex]-1]
                if k==0:            
                    strfragvdwatomindex=str(fragvdwatomindex+1)
                fragvdwatomindex+=1
                fragclasskey=fragidxtosymclass[fragvdwatomindex]
                fragclasskey=int(fragclasskey)
                parentclasskeytofragclasskey[classkey]=fragclasskey
                classkeytosmartsposarray[classkey]=smilesposstring
                classkeytosmarts[classkey]=fragsmarts
                classkeytoatomindexes[classkey]=fragatomstring
                WriteDictionaryToFile(poltype,classkeytoatomindexes,"classkeytoatomindexes.txt")
        WriteDictionaryToFile(poltype,classkeytosmartsposarray,"classkeytosmartsposarray.txt")
        WriteDictionaryToFile(poltype,classkeytosmarts,"classkeytosmarts.txt")
        WriteDictionaryToFile(poltype,parentclasskeytofragclasskey,"parentclasskeytofragclasskey.txt")
        wholexyz=parentdir+r'/'+poltype.xyzoutfile
        wholemol=parentdir+r'/'+poltype.molstructfname
        parentatoms=poltype.rdkitmol.GetNumAtoms()
        listofjobs,jobtooutputlog,newlog,jobtoinputfilepaths=FragmentJobSetup(poltype,strfragrotbndindexes,tail,listofjobs,jobtooutputlog,tempmol,parentdir,vdwfragment,strfragvdwatomindex,onlyfittorsions,jobtoinputfilepaths)
    os.chdir(parentdir)
    if poltype.setupfragjobsonly==True:
        sys.exit()


    jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix=SetupClusterSubmission(poltype,listofjobs,parentdir)
    finishedjobs,errorjobs=SubmitFragmentJobs(poltype,listofjobs,jobtooutputlog,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
    os.chdir(parentdir)

    return equivalentrotbndindexarrays,rotbndindextoringtor,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray


def MatchOBMols(poltype,molstruct,equivalentmolstruct,parentindextofragindex,equivalentparentindextofragindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    fragindices=list(parentindextofragindex.values())
    equivalentfragindices=list(equivalentparentindextofragindex.values())
    indextoreferenceindex={}
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('mol')
    outputname='temp.mol'
    tmpconv.WriteFile(equivalentmolstruct,outputname)
    newmol=rdmolfiles.MolFromMolFile(outputname,removeHs=False)
    smarts=rdmolfiles.MolToSmarts(newmol).replace('@','')
    tmpconv.WriteFile(molstruct,outputname)
    molstructrdkit=rdmolfiles.MolFromMolFile(outputname,removeHs=False)
    p = Chem.MolFromSmarts(smarts)
    matches=molstructrdkit.GetSubstructMatches(p)
    for othermatch in matches: 
       indices=list(range(len(othermatch)))
       smartsindextomoleculeindex=dict(zip(indices,othermatch)) 
       matches=newmol.GetSubstructMatches(p)
       for match in matches: 
           indices=list(range(len(match)))
           smartsindextoequivalentmoleculeindex=dict(zip(indices,match)) 
           moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}
           for moleculeindex,smartsindex in moleculeindextosmartsindex.items():
               refindex=smartsindextoequivalentmoleculeindex[smartsindex]
               indextoreferenceindex[moleculeindex]=refindex
   
    return indextoreferenceindex


def ConvertSameBondTypes(poltype,rotbndindexes,parentindextofragindex,otherparentindextofragindex,othermolindextoequivmolindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    foundall=True
    for index in rotbndindexes:
        if index not in parentindextofragindex.keys():
            foundall=False
    fragindextoparentindex = {v: k for k, v in parentindextofragindex.items()}
    if foundall==False:
        comb=[]
        for index in rotbndindexes:
            otherfragindex=otherparentindextofragindex[index]
            equivfragindex=othermolindextoequivmolindex[otherfragindex]
            parentindex=fragindextoparentindex[equivfragindex]
            comb.append(parentindex)

        bonded=CheckIfIndicesBonded(poltype,comb,poltype.rdkitmol)
        if bonded==True:
            return comb

    else:
        return rotbndindexes


def CheckIfIndicesBonded(poltype,comb,themol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bonded=False
    firstatomidx=comb[0]
    secondatomidx=comb[1]
    firstatom=themol.GetAtomWithIdx(firstatomidx)
    for atom in firstatom.GetNeighbors():
        atomidx=atom.GetIdx()
        if atomidx==secondatomidx:
            bonded=True

    return bonded





def GenerateAtomIndexToSMARTSPosition(poltype,fragidxarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomindextosmartspos={}
    totalatoms=len(fragidxarray)
    for i in range(totalatoms):
        atomindex=i+1
        if (atomindex+1) in fragidxarray:
            fragidxarraypos=fragidxarray.index(atomindex+1)
            smilespos=fragidxarraypos+1
            atomindextosmartspos[atomindex]=smilespos   

    return atomindextosmartspos

def GenerateSMARTSPositionStringAndAtomIndices(poltype,torsion,parentindextofragindex,fragidxarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    smilesposarray=[]
    fragtor=[]
    for index in torsion:
        fragindex=parentindextofragindex[index-1]
        fragtor.append(fragindex)
        fragidxarraypos=fragidxarray.index(fragindex+1)
        smilespos=fragidxarraypos+1
        smilesposarray.append(smilespos)
    smilesposarray=[str(i) for i in smilesposarray]
    smilesposstring=','.join(smilesposarray)
    fragtor=[str(i) for i in fragtor]
    fragtorstring=[str(i) for i in fragtor]
    fragtorstring=','.join(fragtorstring)
    return smilesposstring,fragtorstring


def GrabParentTorsions(poltype,rotbndindextoringtor,array,strparentrotbndindexes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tors=[]
    tortor=False
    nonaroringfrag=False
    maintortors=[]
    listoftors=[]
    rotbndindextotors={}
    for i in range(len(array)):
        rotbndindex=array[i]
        rotkey=rotbndindex.replace('_',' ')
        if rotbndindex in rotbndindextoringtor.keys():
            nonaroringfrag=True
            torset=rotbndindextoringtor[rotbndindex]
            for tor in torset:
                tors.append(tor)
               
        elif rotkey in poltype.rotbndlist.keys():
            tors.extend(list(poltype.rotbndlist[rotkey]))
            listoftors.append(list(poltype.rotbndlist[rotkey]))
            rotbndindextotors[rotbndindex]=list(poltype.rotbndlist[rotkey])
        else:
            tortor=True
            rotkeysplit=rotkey.split()
            rotkeys=[]
            maintortors=[]
            for j in range(0,len(rotkeysplit),2):
                curkey=str(rotkeysplit[j])+' '+str(rotkeysplit[j+1])
                rotkeys.append(curkey)
            for keyidx in range(len(rotkeys)):
                key=rotkeys[keyidx]
                keytors=list(poltype.rotbndlist[key])

                maintortors.append(keytors[0])
                tors.extend(keytors)
    
    
    return tors,maintortors,tortor,nonaroringfrag,rotbndindextotors


def CountUnderscores(poltype,string):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    count=0
    for e in string:
        if e=='_':
            count+=1
    return count

def MakeFileName(poltype,string,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(filename,'w')
    temp.write(string+'\n')
    temp.close()


def WriteDictionaryToFile(poltype,dictionary,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    with open(filename,'w') as f: 
        json.dump(dictionary, f)



def GrabAtomOrder(poltype,smirks):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indextocoordinates={}
    iteratom = openbabel.OBMolAtomIter(mol)
    for atom in iteratom:
        atomidx=atom.GetIdx()
        coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
        indextocoordinates[atomidx]=coords
    return indextocoordinates

def AddInputCoordinatesAsDefaultConformer(poltype,m,indextocoordinates):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        x,y,z = indextocoordinates[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    return m


def GrabNeigbsBabel(poltype,atomidx):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    neigbindexes=[]
    atom=poltype.mol.GetAtom(atomidx)
    atomatomiter=openbabel.OBAtomAtomIter(atom)
    for natom in atomatomiter:
        natomidx=natom.GetIdx()
        neigbindexes.append(natomidx)

    return neigbindexes


def GenerateFrag(poltype,molindexlist,mol,torset):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    vdwfrag=True
    for tor in torset:
        if len(tor)>1:
            vdwfrag=False   

 
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

    bondstoremove=[]
    
    atomswithcutbonds=[]
    if poltype.refinenonaroringtors==False:
        bonditer=openbabel.OBMolBondIter(poltype.mol)
        for bond in bonditer:
            oendidx = bond.GetEndAtomIdx()
            obgnidx = bond.GetBeginAtomIdx()
            if oendidx in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
                continue
            if oendidx not in oldindextonewindex.keys() and obgnidx in oldindextonewindex.keys():
                continue
            if oendidx not in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
                continue
            ringbond=bond.IsInRing()
            bgnatom=bond.GetBeginAtom()
            endatom=bond.GetEndAtom()
            bgnhyb=bgnatom.GetHyb()
            endhyb=endatom.GetHyb()
            if ringbond==True and bgnhyb!=2 and endhyb!=2 and vdwfrag==False:
                atomindices=RingAtomicIndices(poltype,mol)
                a=molindexlist[0]
                b=molindexlist[1]
                c=molindexlist[2]
                d=molindexlist[3]
                ring=GrabRingAtomIndicesFromInputIndex(poltype,[a,b,c,d],atomindices)
                if ring!=None:
                    if c in ring and b in ring:
                        if len(ring)==4: 
                           ls=[a,d] 
                           bondstoremove.append(ls)
                           atomswithcutbonds.append(oldindextonewindex[a])
                           atomswithcutbonds.append(oldindextonewindex[d])
                        elif len(ring)==5:
                           if a in ring:
                               idx=a
                           elif d in ring:
                               idx=d
                            
                           neigbindexes=GrabNeigbsBabel(poltype,idx)
                           for idx in neigbindexes:
                               if idx in ring and idx!=a and idx!=b and idx!=c and idx!=d:
                                   theidx=idx
                                   break
                           ls=[a,theidx] 
                           atomswithcutbonds.append(oldindextonewindex[a])
                           atomswithcutbonds.append(oldindextonewindex[theidx])
                           bondstoremove.append(ls)
                        elif len(ring)==6: 
                           ls=[]
                           for idx in ring:
                               if idx!=a and idx!=b and idx!=c and idx!=d:
                                   ls.append(idx) 
                           newls=[]
                           for idx in ls:
                               if idx in oldindextonewindex.keys():
                                   newls.append(idx)
                           bondstoremove.append(newls)
                           for idx in newls:
                               atomswithcutbonds.append(oldindextonewindex[idx])




    bonditer=openbabel.OBMolBondIter(poltype.mol)
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        ls=[oendidx,obgnidx]
        if ls in bondstoremove or ls[::-1] in bondstoremove:
            continue
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
    RemoveRadicals(poltype,filename)
    indextocoordinates=GrabIndexToCoordinates(poltype,em) 
    nem=ReadToOBMol(poltype,filename)
    nem.AddHydrogens()
    filename='fragalladdedhyd.mol'
    WriteOBMolToMol(poltype,nem,filename)

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
    try:
        status = AllChem.EmbedMolecule(newmol,maxAttempts=5000)
        indextocoordinates=UpdateAtomNumbers(poltype,indextocoordinates)
        rdkitindextocoordinates={}
        for idx,coords in indextocoordinates.items():
            rdkitidx=idx-1
            rdkitindextocoordinates[rdkitidx]=coords
        newmol=AddInputCoordinatesAsDefaultConformer(poltype,newmol,rdkitindextocoordinates)
    except:
        pass
    rdkitoldindextonewindex={}
    for oldindex,newindex in oldindextonewindex.items():
        rdkitoldindex=oldindex-1
        rdkitnewindex=newindex-1
        rdkitoldindextonewindex[rdkitoldindex]=rdkitnewindex
    newmol=AssignTotalCharge(poltype,newmol,nem)
    return newmol,rdkitoldindextonewindex


def RemoveRadicals(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    temp=open(filename,'w')
    for line in results:
        linesplit=line.split()
        if len(linesplit)==16:
            nicesplit=re.split(r'(\s+)', line)
            nicesplit[20]='0'
            line=''.join(nicesplit)
        temp.write(line)
    temp.close()

def UpdateAtomNumbers(poltype,indextocoordinates):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newindextocoordinates={}
    ls=list(indextocoordinates.keys())
    for i in range(len(ls)):
        originalbabelindex=ls[i]
        newbabelindex=i+1
        coords=indextocoordinates[originalbabelindex]
        newindextocoordinates[newbabelindex]=coords
    return newindextocoordinates


def WriteOBMolToSDF(poltype,mol,outputname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('sdf')
    atomiter=openbabel.OBMolAtomIter(mol)

    tmpconv.WriteFile(mol,outputname)


def WriteOBMolToXYZ(poltype,mol,outputname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('xyz')
    tmpconv.WriteFile(mol,outputname)


def WriteOBMolToMol(poltype,mol,outputname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('mol')
    tmpconv.WriteFile(mol,outputname)

def WriteRdkitMolToMolFile(poltype,mol,outputname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    rdmolfiles.MolToMolFile(mol,outputname,kekulize=True)
    rdmolfiles.MolToMolFile(mol,outputname,kekulize=False)


def ReadRdkitMolFromMolFile(poltype,inputname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    rdkitmol=rdmolfiles.MolFromMolFile(inputname,sanitize=False)
    return rdkitmol

def ReadMolFileToOBMol(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpconv = openbabel.OBConversion()
    tmpconv.SetInFormat('mol')
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,filename)
    return fragmolbabel

def ReadToOBMol(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpconv = openbabel.OBConversion()
    inFormat = tmpconv.FormatFromExt(filename)
    tmpconv.SetInFormat(inFormat)
    fragmolbabel=openbabel.OBMol()
    tmpconv.ReadFile(fragmolbabel,filename)
    return fragmolbabel



def mol_with_atom_index(poltype,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx()+1 ) )
    return mol

def mol_with_atom_index_removed(poltype,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        atom=mol.GetAtomWithIdx(idx)
        atom.ClearProp('molAtomMapNumber')
    return mol



def GenerateWBOMatrix(poltype,molecule,moleculebabel,structfname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    error=False
    WBOmatrix=None
    curespmethod=poltype.espmethod
    curspbasisset=poltype.espbasisset
    poltype.espmethod='HF'
    poltype.espbasisset='MINIX'
    charge=Chem.rdmolops.GetFormalCharge(molecule)
    inputname,outputname=esp.CreatePsi4ESPInputFile(poltype,structfname,poltype.comespfname.replace('.com','_frag.com'),moleculebabel,poltype.maxdisk,poltype.maxmem,poltype.numproc,charge,False)
    finished,error=poltype.CheckNormalTermination(outputname)
    cmdstr='psi4 '+inputname+' '+outputname
    if finished==False:
        try:
             poltype.call_subsystem([cmdstr],True)

             temp={cmdstr:outputname} 
             finishedjobs,errorjobs=poltype.WaitForTermination(temp,False)

        except Exception:
             error=True
    try:
        WBOmatrix=GrabWBOMatrixPsi4(poltype,outputname,molecule)
            
    except:
        cmdstr='psi4 '+inputname+' '+outputname
        temp={cmdstr:outputname}
        poltype.call_subsystem([cmdstr],True)
        finishedjobs,errorjobs=poltype.WaitForTermination(temp,False)

        WBOmatrix=GrabWBOMatrixPsi4(poltype,outputname,molecule)

    poltype.espmethod=curespmethod
    poltype.espbasisset=curspbasisset

    return WBOmatrix,outputname,error



def GenerateFragments(poltype,mol,torlist,parentWBOmatrix,missingvdwatomsets,nonaroringtorlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
        if torset in poltype.nonaroringtorsets or torset in missingvdwatomsets or torset in nonaroringtorlist:
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
        fragmol,parentindextofragindex=GenerateFrag(poltype,extendedtorindexes,mol,torset)
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
            torsbeingfit=poltype.nonarotortotorsbeingfit[torset]
            rotbndindextoringtor[rotbndidx]=torsbeingfit
        filename=fragfoldername+'.mol'
        WriteRdkitMolToMolFile(poltype,fragmol,filename)
        fragmoltofragfoldername[fragmol]=fragfoldername
        fragmolbabel=ReadMolFileToOBMol(poltype,filename)
        WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
        WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
        structfname=filename.replace('.mol','.sdf')
        try:
            fragWBOmatrix,outputname,error=GenerateWBOMatrix(poltype,fragmol,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
        except: # assume that fragmenter changed (code changed so structure changed and then need to delete and restart that fragment. Can remove this when code becomes more stable
            os.chdir('..')
            shutil.rmtree(fragfoldername)
            if not os.path.isdir(fragfoldername):
                os.mkdir(fragfoldername)
            os.chdir(fragfoldername)
            WriteRdkitMolToMolFile(poltype,fragmol,filename)
            WriteOBMolToXYZ(poltype,fragmolbabel,filename.replace('.mol','_xyzformat.xyz'))
            WriteOBMolToSDF(poltype,fragmolbabel,filename.replace('.mol','.sdf'))
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
                
        curdir=os.getcwd()
        os.chdir(fragspath)
        growfragments=[mol_with_atom_index(poltype,frag) for frag in growfragments]
        Draw2DMoleculesWithWBO(poltype,growfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,torset,'FragmentGrowthWithIndex')
        sanitizedfragments=[mol_with_atom_index_removed(poltype,frag) for frag in growfragments]
        Draw2DMoleculesWithWBO(poltype,sanitizedfragments,growfragmoltoWBOmatrices,growfragmoltofragfoldername,growfragmoltobondindexlist,torset,'FragmentGrowthWithoutIndex')

        os.chdir(curdir)


        structfname=structfname.replace('_xyzformat.xyz','.sdf')
        structfnamemol=structfname.replace('.sdf','.mol')
        rotbndindextofragment[rotbndidx]=fragmol
        rotbndindextofragmentfilepath[rotbndidx]=fragpath+r'/'+structfnamemol
        rotbndindextoparentindextofragindex[rotbndidx]=parentindextofragindex
        rotbndindextofragWBOmatrix[rotbndidx]=fragWBOmatrix
        rotbndindextofragfoldername[rotbndidx]=fragfoldername
        os.chdir(fragspath)
    # now remove all folders with Hydrated in them, that was just temp storage for producing other folders
    RemoveTempFolders(poltype)
    poltype.rotbndindextofragmentfilepath=rotbndindextofragmentfilepath
    vdwfragmentarray=[]
    vdwnamearray=[]
    torfragmentarray=[]
    tornamearray=[]
    vdwparentindextofragindexarray=[]
    torparentindextofragindexarray=[]

    for rotbndindex in rotbndindextofragment.keys():
        fragment=rotbndindextofragment[rotbndindex]
        parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]

        if '_' not in rotbndindex:
            vdwfragmentarray.append(fragment)
            vdwnamearray.append(rotbndindex)
            vdwparentindextofragindexarray.append(parentindextofragindex)
        else:
            torfragmentarray.append(fragment)
            tornamearray.append(rotbndindex)
            torparentindextofragindexarray.append(parentindextofragindex)

    vdwequivalentrotbndindexarrays=FindEquivalentFragments(poltype,vdwfragmentarray,vdwnamearray,vdwparentindextofragindexarray)
    torequivalentrotbndindexarrays=FindEquivalentFragments(poltype,torfragmentarray,tornamearray,torparentindextofragindexarray)
    equivalentrotbndindexarrays=[]
    equivalentrotbndindexarrays.extend(vdwequivalentrotbndindexarrays)
    equivalentrotbndindexarrays.extend(torequivalentrotbndindexarrays)
    tempdic={} 
    for rotbndindex in rotbndindextofragment.keys():
        if '_' not in rotbndindex:
            parentindextofragindex=rotbndindextoparentindextofragindex[rotbndindex]
            tempdic[rotbndindex]=parentindextofragindex     
    for key,value in tempdic.items():
        rotbndindextoparentindextofragindex[key]=value 

    return rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor




def GenerateFakeTorset(poltype,mol,parentindextofragindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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




  
def RemoveTempFolders(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    foldstoremove=[]
    folds=os.listdir()
    for f in folds:
        if os.path.isdir(f) and 'Hydrated' in f:
            foldstoremove.append(f)
    for f in foldstoremove:
        shutil.rmtree(f)

def ReduceParentMatrix(poltype,parentindextofragindex,fragWBOmatrix,parentWBOmatrix):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
                fragmol,parentindextofragindex=GenerateFrag(poltype,indexlist,mol,torset)
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
           aromaticindexes=GrabAromaticAtoms(poltype,idx+1)
           newindexes=aromaticindexes
           for atmidx in newindexes:
               if atmidx not in indexlist:
                   indexlist.append(atmidx)
           if idx not in indexlist:
               indexlist.append(idx)
        temp=[]
        for idx in indexlist:
           neighbatom=poltype.rdkitmol.GetAtomWithIdx(idx)
           atomicnum=neighbatom.GetAtomicNum()
           if atomicnum==15: # add all neighbors, seems to be issue with cutting bonds on P then adding H, h is added in strange way
               for neighbneighbatom in neighbatom.GetNeighbors():
                   neighbneighbatomidx=neighbneighbatom.GetIdx()
                   if neighbneighbatomidx not in indexlist and neighbneighbatomidx not in temp:
                       temp.append(neighbneighbatomidx) 

           for neighbneighbatom in neighbatom.GetNeighbors():
               atomicnum=neighbneighbatom.GetAtomicNum()
               neighbneighbatomidx=neighbneighbatom.GetIdx()
               if atomicnum==1 and neighbneighbatomidx not in indexlist:
                   temp.append(neighbneighbatomidx)
               bond=poltype.rdkitmol.GetBondBetweenAtoms(neighbneighbatomidx,idx)
               bondorder=bond.GetBondTypeAsDouble()
               babelneighbatom=poltype.mol.GetAtom(neighbneighbatomidx+1)

               babelneighbatomisinring=babelneighbatom.IsInRing()
               if bondorder>1 and neighbneighbatomidx not in indexlist and babelneighbatomisinring==False:
                   temp.append(neighbneighbatomidx)
        for idx in temp:
            if idx not in indexlist:
                indexlist.append(idx)

        if indexlist not in possiblefragatmidxs:
           possiblefragatmidxs.append(indexlist)
    return possiblefragatmidxs


def WriteOutFragmentInputs(poltype,fragmol,fragfoldername,fragWBOmatrix,parentWBOmatrix,WBOdifference,parentindextofragindex,torset,fragmoltoWBOmatrices,fragmoltobondindexlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
   """
   Intent:
   Input:
   Output:
   Referenced By: 
   Description: 
   """
   molindexlist=[i-1 for i in tor]
   if len(molindexlist)==4: # then add all a,d for a-b-c-d, this way if a is H, then neighbors of the other a still get added to fragment
       for atom in poltype.rdkitmol.GetAtoms():
           atomindex=atom.GetIdx()
           babelatomindex=atomindex+1
           if atomindex==molindexlist[1] or atomindex==molindexlist[2]:
               for neighbatom in atom.GetNeighbors():
                   neighbatomindex=neighbatom.GetIdx()
                   if neighbatomindex not in molindexlist:
                       molindexlist.append(neighbatomindex)
   originalatomindexlist=molindexlist[:]
   for atom in poltype.rdkitmol.GetAtoms():
       atomindex=atom.GetIdx()
       babelatomindex=atomindex+1
       grabneighbs=False
       if atomindex in originalatomindexlist:
           if atomindex not in molindexlist:
               molindexlist.append(atomindex)
           grabneighbs=True
           if grabneighbs==True:
               for neighbatom in atom.GetNeighbors():
                   neighbatomindex=neighbatom.GetIdx()
                   babelneighbatom=poltype.mol.GetAtom(neighbatomindex+1)
                   babelneighbatomisinring=babelneighbatom.IsInRing()
                   babelneighbatomhyb=babelneighbatom.GetHyb()
                   if neighbatomindex not in molindexlist:
                       molindexlist.append(neighbatomindex)
                       if neighbatom.GetIsAromatic() or (babelneighbatomisinring==True and babelneighbatomhyb==2):
                           aromaticindexes=GrabAromaticAtoms(poltype,neighbatomindex+1)
                           newindexes=aromaticindexes
                           for atmidx in newindexes:
                               if atmidx not in molindexlist:
                                   molindexlist.append(atmidx)

   temp=[]
   for index in molindexlist:
       atom=poltype.rdkitmol.GetAtomWithIdx(index)
       atomicnum=atom.GetAtomicNum()
       if atomicnum==15: # add all neighbors, seems to be issue with cutting bonds on P then adding H, h is added in strange way
           for neighbneighbatom in atom.GetNeighbors():
               neighbneighbatomidx=neighbneighbatom.GetIdx()
               if neighbneighbatomidx not in temp:
                   temp.append(neighbneighbatomidx) 
       for neighbneighbatom in atom.GetNeighbors():
           atomicnum=neighbneighbatom.GetAtomicNum()
           neighbneighbatomidx=neighbneighbatom.GetIdx()
           babelneighbatom=poltype.mol.GetAtom(neighbneighbatomidx+1)
           babelneighbatomisinring=babelneighbatom.IsInRing()

           if atomicnum==1 and neighbneighbatomidx not in molindexlist:
               temp.append(neighbneighbatomidx)
           bond=poltype.rdkitmol.GetBondBetweenAtoms(neighbneighbatomidx,index)
           bondorder=bond.GetBondTypeAsDouble()
           if bondorder>1 and neighbneighbatomidx not in molindexlist and babelneighbatomisinring==False:
               temp.append(neighbneighbatomidx)
   for idx in temp:
       if idx not in molindexlist:
           molindexlist.append(idx)
   return molindexlist

def Chunks(poltype,lst, n):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ChunksList(poltype,gen):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newlst=[]
    for item in gen:
        newlst.append(item)
    return newlst


def Draw2DMoleculesWithWBO(poltype,fragments,fragmoltoWBOmatrices,fragmoltofragfoldername,fragmoltobondindexlist,torset,basestr):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
        smarts=smarts.replace('@','').replace('/[H]','').replace('\[H]','').replace('H3','').replace('H2','').replace('H','')
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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


def RingAtomicIndices(poltype,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    sssr = mol.GetSSSR()
    atomindices=[]
    for ring in sssr:
        ringatomindices=GrabRingAtomIndices(poltype,mol,ring)
        atomindices.append(ringatomindices)        
    return atomindices

def GrabRingAtomIndicesFromInputIndex(poltype,atomindexlist,atomindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ringtocount={}
    for ring in atomindices:
        count=0
        for atomindex in atomindexlist:
            if atomindex in ring:
                count+=1
        ringtocount[tuple(ring)]=count
    if len(ringtocount.keys())!=0:
        maxcount=max(ringtocount.values())
        if maxcount!=0:
            for ring,count in ringtocount.items():
                if count==maxcount:
                    return ring

    ring=None
    return ring

def GrabRingAtomIndices(poltype,mol,ring):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ringatomindices=[]
    atomiter=openbabel.OBMolAtomIter(mol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if ring.IsInRing(atomidx)==True:
            ringatomindices.append(atomidx)
    return ringatomindices



def GrabAromaticAtoms(poltype,neighbatomidx):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomindices=RingAtomicIndices(poltype,poltype.mol)
    ring=GrabRingAtomIndicesFromInputIndex(poltype,[neighbatomidx],atomindices)
    aromaticindexes=[]
    if ring!=None:
        babelatoms=[poltype.mol.GetAtom(atmindex) for atmindex in ring]
        hybs=[babelatom.GetHyb() for babelatom in babelatoms]
        allflat=True
        for hyb in hybs:
            if hyb!=2:
                allflat=False
        if allflat==True:
             ring=[i-1 for i in ring]
             aromaticindexes.extend(ring) 
    
    return aromaticindexes


def PlotFragmenterResults(poltype,WBOdiffarray,molarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    fig=plt.figure()
    basename='NumberofAtomsVSWBODifference'
    plt.plot(WBOdiffarray,[m.GetNumAtoms() for m in molarray],'.')
    plt.xlabel('WBO Difference')
    plt.ylabel('Number of atoms in fragment')
    fig.savefig(basename+'.png')


