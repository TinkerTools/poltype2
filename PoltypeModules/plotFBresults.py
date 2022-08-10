import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats, optimize, interpolate
from sklearn.metrics import mean_squared_error
from math import sqrt
import csv
from scipy.interpolate import interp1d
from matplotlib.pyplot import cm
from openbabel import openbabel
from PyAstronomy import pyasl
import mdtraj as md
import csv
import itertools
from PIL import Image
import pandas as pd
import random
from random import randint
import shutil
from pathlib import Path


def GrabJobDirectories(fbdir):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    jobdirs=[]
    files=os.listdir(fbdir)
    for f in files:
        path=os.path.join(fbdir,f)
        if os.path.isdir(path):
            newfiles=os.listdir(path)
            found=False
            for newf in newfiles:
                if 'optimize.out' in newf:
                    found=True
            if found==True:
                jobdirs.append(path)
    return jobdirs


def GrabOutputFiles(jobdirs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    outputfiles=[]
    for path in jobdirs:
        os.chdir(path)
        files=os.listdir()
        for f in files:
            if 'optimize.out' in f and 'nohup' not in f:
                filepath=os.path.join(path,f)
                outputfiles.append(filepath)

    return outputfiles


def GrabCubeFiles(groupedpoltypedirs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    nametocubefiles={}
    curdir=os.getcwd()
    groupednames=[]
    for grp in groupedpoltypedirs:
        newgrp=[]
        for path in grp:
            os.chdir(path)
            files=os.listdir()
            array=[]
            for f in files:
                if f=='ESP.cube' or f=='Dt.cube':
                    filepath=os.path.join(path,f)
                    array.append(filepath)
                if '-esp.log' in f:
                    name=f.replace('-esp.log','')
                    newgrp.append(name)
            if len(array)!=0:
                nametocubefiles[name]=array
        groupednames.append(newgrp)
    

    os.chdir(curdir)
    return nametocubefiles,groupednames


def GrabPoltypeDirectories(jobdirs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltypedirs=[]
    groupedpoltypedirs=[]
    curdir=os.getcwd()
    for path in jobdirs:
        os.chdir(path)
        files=os.listdir()
        for f in files:
            if 'poltype.ini' in f:
                filepath=os.path.join(path,f)
                paths=ReadInputGrabPaths(filepath)
                poltypedirs.extend(paths)
                groupedpoltypedirs.append(paths)

    os.chdir(curdir)
    return poltypedirs,groupedpoltypedirs


def ReadInputGrabPaths(filepath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(filepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'poltypepathlist' in line:
            split=line.split('=')
            res=split[1]
            paths=res.split(',')
            paths=[a.replace('\n','') for a in paths]
            paths=[a.strip() for a in paths] 


    return paths


def GrabResultsAllMolecules(outputfiles):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tpdiclist=[]
    qmdiclist=[]
    newoutputfiles=[]
    for outputfile in outputfiles:
        try:
            tpdic,qmdic,outputfile=GrabResults(outputfile)
            tpdiclist.append(tpdic)
            qmdiclist.append(qmdic)
            newoutputfiles.append(outputfile)
        except:
            pass
    return tpdiclist,qmdiclist,newoutputfiles

def GrabResults(outputfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tpdic={}
    qmdic={}
    temp=open(outputfile,'r')
    results=temp.readlines()
    temp.close()
    foundblock=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if '#============================' in line:
            prevline=results[lineidx-1]
            if 'Temperature  Pressure' in prevline or 'Label' in prevline:
                foundblock=True
                prevprevline=results[lineidx-2]
                prevprevlinesplit=prevprevline.split()
                refkeyword='Ref'
                calckeyword='Calc'
                errorkeyword='CalcErr'
                losskeyword='Loss'
                prevprevprevline=results[lineidx-3]
                prevprevprevlinesplit=prevprevprevline.split()
                foundQM=False
                foundliq=False
                foundweight=False
                if 'Enthalpy' in prevprevline:
                   index=prevprevlinesplit.index('Enthalpy')
                   totalname=prevprevlinesplit[index-1]
                   namesplit=totalname.split('_')
                   name=namesplit[1]
                   propertyval='Enthalpy'
                   foundliq=True
                elif 'Density' in prevprevline:
                   propertyval='Density'
                   index=prevprevlinesplit.index('Density')
                   totalname=prevprevlinesplit[index-1]
                   namesplit=totalname.split('_')
                   name=namesplit[1]
                   foundliq=True
                elif 'Interaction Energies' in prevprevline:
                   propertyval='Interaction Energies'
                   index=prevprevprevlinesplit.index('Target:')
                   targetname=prevprevprevlinesplit[index+1]
                   namesplit=targetname.split('_')
                   name=namesplit[1]
                   foundQM=True
                continue
            if 'Property Name' in prevline:
                foundblock=True
                foundweight=False
                foundQM=False
                foundliq=False
                prevprevline=results[lineidx-2]
                prevprevlinesplit=prevprevline.split()
            if 'Property Name' in prevline and 'Condensed Phase' in prevprevline:
                foundweight=True
                index=prevprevlinesplit.index('Condensed')
                totalname=prevprevlinesplit[index-1]
                namesplit=totalname.split('_')
                name=namesplit[1]
                continue
                

        if '--------------------------------------------------' in line:
            foundblock=False 
        if foundblock==True:
            linesplit=line.split()
            if foundweight==True:
                linesplit=line.split()
                weightkeyword='weight'
                if 'Density' in line:
                    propertyval='Density'
                elif 'Enthalpy' in line:
                    propertyval='Enthalpy'
                else:
                    continue

                weight=float(linesplit[-2])
                tparray=list(tpdic[name].keys())
                for tp in tparray:   
                    tpdic[name][tp][propertyval][weightkeyword]=weight
            if foundliq==True:
                if name not in tpdic.keys():
                    tpdic[name]={}
                temp=float(linesplit[0])
                pressure=float(linesplit[1])
                tp=tuple([temp,pressure])
                if tp not in tpdic[name].keys():
                    tpdic[name][tp]={}
                if propertyval not in tpdic[name][tp].keys():
                    tpdic[name][tp][propertyval]={}
                if refkeyword not in tpdic[name][tp][propertyval].keys():
                    tpdic[name][tp][propertyval][refkeyword]=[]
                if calckeyword not in tpdic[name][tp][propertyval].keys():
                    tpdic[name][tp][propertyval][calckeyword]=[]
                if errorkeyword not in tpdic[name][tp][propertyval].keys():
                    tpdic[name][tp][propertyval][errorkeyword]=[]
                if losskeyword not in tpdic[name][tp][propertyval].keys():
                    tpdic[name][tp][propertyval][losskeyword]=[]

                ref=float(linesplit[3])
                calc=float(linesplit[4])
                error=float(linesplit[6])
                loss=float(linesplit[-1])
                tpdic[name][tp][propertyval][refkeyword].append(ref)
                tpdic[name][tp][propertyval][errorkeyword].append(error)
                tpdic[name][tp][propertyval][calckeyword].append(calc)
                tpdic[name][tp][propertyval][losskeyword].append(loss)
            elif foundQM==True:
                if targetname not in qmdic.keys():
                    qmdic[targetname]={}
                label=int(linesplit[0])
                ref=float(linesplit[2])
                calc=float(linesplit[1])
                if label not in qmdic[targetname].keys():
                    qmdic[targetname][label]={} 
                if refkeyword not in qmdic[targetname][label].keys():
                    qmdic[targetname][label][refkeyword]=[]
                if calckeyword not in qmdic[targetname][label].keys():
                    qmdic[targetname][label][calckeyword]=[]
                qmdic[targetname][label][refkeyword].append(ref)
                qmdic[targetname][label][calckeyword].append(calc)
 
    return tpdic,qmdic,outputfile


def PlotAllFBJobs(tpdiclist,qmtargetnamedic,nametofilenametoformula,truenametoindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    
    nametotptofinalprops={}
    nametoformulatormse={}
    nametoformulatomse={}
    listofmoleculenametoproptoimagenames=[]
    for tpdic in tpdiclist:
        nametotptofinalprops,moleculetoproptoimagenames=PlotFBLiq(tpdic,nametotptofinalprops)
        listofmoleculenametoproptoimagenames.append(moleculetoproptoimagenames)
    nametofigs={}
    nametoaxes={}
    for qmtarget,targetdic in qmtargetnamedic.items():
        namesplit=qmtarget.split('-')
        name='-'.join(namesplit[:-2])
        namesplit=name.split('-')
        if name.count('-')>=2:
            name='-'.join(namesplit[:-1])
        namesplit=name.split('_')
        name=namesplit[0]
        if name not in nametofigs.keys():
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            nametofigs[name]=fig
            nametoaxes[name]=ax1

        filenametoformula=nametofilenametoformula[name]
        n=len(list(filenametoformula.keys()))
        color = list(iter(cm.rainbow(np.linspace(0, 1, n))))
        filenametocolor=dict(zip(filenametoformula.keys(),color))
        if qmtarget in filenametoformula.keys():
            indices=truenametoindices[qmtarget] 
            namepath=os.path.join(os.getcwd(),name)
            if not os.path.isdir(namepath):
                os.mkdir(name)
            os.chdir(name)
            formula=filenametoformula[qmtarget] 
            c=filenametocolor[qmtarget]
            x,rmsvalues,msvalues=PlotFBQM(targetdic,formula,indices)
            nametoaxes[name].scatter(x,rmsvalues, s=10, c=[c],marker="s", label=formula)
            x=np.array(x)
            x_new = np.linspace(x.min(),x.max(),500)
            try:
                f = interp1d(x,rmsvalues, kind='quadratic')
                y_smooth=f(x_new)
                nametoaxes[name].plot(x_new,y_smooth,color=c)
            except:
                pass
            ykey='RMSE Interaction Energy (kcal/mol)'
            xkey='Iteration'
            nametoaxes[name].set_ylabel(ykey,fontsize=12)
            nametoaxes[name].set_xlabel(xkey,fontsize=12)
            newtitle=name+' RMSE Interaction Energy Vs Iteration'
            nametoaxes[name].set_title(newtitle)
            nametoaxes[name].legend()
            imagename=newtitle+'.png'
            nametofigs[name].savefig(imagename)
            if name not in nametoformulatormse.keys():
                nametoformulatormse[name]={}
                nametoformulatomse[name]={}
            nametoformulatormse[name][formula]=rmsvalues[-1]
            nametoformulatomse[name][formula]=msvalues[-1]
            os.chdir('..')

    return nametotptofinalprops,nametoformulatormse,nametoformulatomse,listofmoleculenametoproptoimagenames

def GrabStructure(formula):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    files=os.listdir()
    for f in files:
        if formula in f and '.xyz' in f and 'tinkermincart' not in f:
            return f

def GrabMolecule(newname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    obConversion = openbabel.OBConversion()
    newmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(newname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(newmol,newname)
    return newmol


def GrabDistanceAndElements(mol,indices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atoms=[mol.GetAtom(i) for i in indices]
    coords=[np.array([atom.GetX(),atom.GetY(),atom.GetZ()]) for atom in atoms]
    distance=np.linalg.norm(coords[0]-coords[1])
    atomicnums=[atom.GetAtomicNum() for atom in atoms]
    an = pyasl.AtomicNo()
    elements=[an.getElSymbol(atomicnum) for atomicnum in atomicnums]
    prefix='-'.join(elements)
    return distance,prefix

def PlotFBQM(targetdic,formula,indices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    distances=targetdic['Distances']
    refvalues=np.transpose(np.array(targetdic['Ref']))
    calcvalues=np.transpose(np.array(targetdic['Calc']))
    rmsvalues=[]
    msvalues=[]
    struct=GrabStructure(formula)
    mol=GrabMolecule(struct)
    distance,prefix=GrabDistanceAndElements(mol,indices)
    distances=np.array(distances)*distance
    for i in range(len(refvalues)):
        refs=np.array(refvalues[i])
        refs=refs-min(refs)
        calcs=np.array(calcvalues[i])
        calcs=calcs-min(calcs)
        rms = sqrt(mean_squared_error(calcs, refs))
        rmsvalues.append(rms)
        diff=calcs-refs
        ms=np.mean(diff)
        msvalues.append(ms)
        label1='AMOEBA'
        label2='QM'
        ykey='Interaction Energy (kcal/mol)'
        xkey=prefix+' Distance (Angstroms)'
        labels=[label1,label2]
        title=formula+' Interaction Energy Vs Distance'
        if i==len(refvalues)-1:
            if len(distances)==len(calcs):
                ScatterPlot2D(title,distances,calcs,refs,xkey,ykey,labels)
                title=formula+' AMOEBA Energy Vs QM Energy'
                ykey='AMOEBA Interaction Energy (kcal/mol)'
                xkey='QM Interaction Energy (kcal/mol)'
                labels=[label1]
                ScatterPlot1D(title,refs,calcs,xkey,ykey,labels,interpolate=False,correlation=True)
    x=list(range(len(rmsvalues)))
    label1='RMSE Interaction Energy'
    ykey='RMSE Interaction Energy (kcal/mol)'
    xkey='Iteration'
    labels=[label1]
    title=formula+' RMSE vs Iteration'
    #ScatterPlot1D(title,x,rmsvalues,xkey,ykey,labels)
    return x,rmsvalues,msvalues


def PlotFBLiq(tpdic,nametotptofinalprops):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    moleculetoproptoimagenames={}
    for name in tpdic.keys():
        dic=tpdic[name]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        n=len(list(dic.keys()))
        color = iter(cm.rainbow(np.linspace(0, 1, n)))
        proptofigs={}
        proptoaxes={}
        if name not in nametotptofinalprops.keys():
            nametotptofinalprops[name]={}
        for tp in dic.keys():
            innerdic=dic[tp]
            c = next(color)
            if tp not in nametotptofinalprops[name].keys():
                nametotptofinalprops[name][tp]={}
            for propertyval in innerdic.keys():
                relprop=propertyval+'RelativeError'
                if propertyval not in nametotptofinalprops[name][tp].keys():
                    nametotptofinalprops[name][tp][propertyval]={}
                if propertyval not in proptofigs.keys():
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    proptofigs[propertyval]=fig
                    proptoaxes[propertyval]=ax1
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    proptofigs[relprop]=fig
                    proptoaxes[relprop]=ax1


                innermostdic=innerdic[propertyval]
                refkey='Ref'
                calckey='Calc'
                errkey='CalcErr' 
                losskey='Loss'
                weightkey='weight'
                weight=innermostdic[weightkey]
                refarray=np.array(innermostdic[refkey])
                calcarray=np.array(innermostdic[calckey])
                errarray=np.array(innermostdic[errkey])
                lossarray=np.array(innermostdic[losskey])*weight
                if propertyval=='Enthalpy':
                    scale=0.239006 # convert kj to kcal
                    units='(kcal/mol)'
                else:
                    scale=1
                    units='$kg/m^3$'
                refarray=scale*refarray
                calcarray=scale*calcarray
                nametotptofinalprops[name][tp][propertyval][refkey]=refarray[-1]
                nametotptofinalprops[name][tp][propertyval][calckey]=calcarray[-1]
                nametotptofinalprops[name][tp][propertyval][errkey]=errarray[-1]
                nametotptofinalprops[name][tp][propertyval][losskey]=lossarray[-1]
                err=np.abs(refarray-calcarray)
                relerr=100*(err/refarray)
                label1=propertyval+' AMOEBA'
                label2=propertyval+' Exp'
                ykey=propertyval+units
                xkey='Iterations'
                labels=[label1,label2]
                string='T=%s,P=%s'%(str(tp[0]),str(tp[1]))
                title=name+' '+propertyval+' '+string
                x=list(range(len(refarray)))
                #ScatterPlot2D(title,x,calcarray,refarray,xkey,ykey,labels)
                proptoaxes[propertyval].scatter(x,calcarray, s=10, c=[c], marker="o", label='AMOEBA '+string)
                proptoaxes[propertyval].errorbar(x, calcarray, yerr=errarray,c=c, fmt="o")
                proptoaxes[propertyval].scatter(x,refarray, s=10, c=[c], marker="s", label='Exp '+string)              
                proptoaxes[propertyval].set_ylabel(ykey,fontsize=12)
                proptoaxes[propertyval].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+propertyval
                proptoaxes[propertyval].set_title(newtitle)
                proptoaxes[propertyval].legend()
                x=np.array(x)
                refarray=np.array(refarray)
                calcarray=np.array(calcarray)
                x_new = np.linspace(x.min(),x.max(),500)
                try:
                    f = interp1d(x,calcarray, kind='quadratic')
                    y_smooth=f(x_new)
                    proptoaxes[propertyval].plot(x_new,y_smooth,color=c)
                    f = interp1d(x,refarray, kind='quadratic')
                    y_smooth=f(x_new)
                    proptoaxes[propertyval].plot(x_new,y_smooth,color=c)
                except:
                    pass
                imagename=newtitle.replace(' ','_')+'.png'
                proptofigs[propertyval].savefig(imagename)
                if name not in moleculetoproptoimagenames.keys():
                    moleculetoproptoimagenames[name]={}
                moleculetoproptoimagenames[name][propertyval]=os.path.join(os.getcwd(),imagename)

                proptoaxes[relprop].scatter(x,relerr, s=10, c=[c], marker="s", label=string)
                try:
                    f = interp1d(x,relerr, kind='quadratic')
                    y_smooth=f(x_new)
                    proptoaxes[relprop].plot(x_new,y_smooth,color=c)
                except:
                    pass
                ykey='Relative Error '+propertyval+units
                proptoaxes[relprop].set_ylabel(ykey,fontsize=12)
                proptoaxes[relprop].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+'Relative Error '+propertyval
                proptoaxes[relprop].set_title(newtitle)
                proptoaxes[relprop].legend()
                imagename=newtitle+'.png'
                proptofigs[relprop].savefig(imagename)

 
        os.chdir('..')
    return nametotptofinalprops,moleculetoproptoimagenames

def ScatterPlot2D(title,x,y1,y2,xkey,ykey,labels):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    x=np.array(x)
    y1=np.array(y1)
    y2=np.array(y2)
    plt.figure()
    ax = plt.axes()
    calc=ax.scatter(x,y1,color='r')
    ref=ax.scatter(x,y2,color='b')
    imagename=title+'.png'
    ax.set_ylabel(ykey,fontsize=12)
    ax.set_xlabel(xkey,fontsize=12)
    plt.title(title)
    plt.legend((calc, ref),(labels[0], labels[1]),fontsize=8)
    x_new = np.linspace(x.min(),x.max(),500)
    try:
        f = interp1d(x,y1, kind='quadratic')
        y_smooth=f(x_new)
        plt.plot(x_new,y_smooth,color='red')
        f = interp1d(x,y2, kind='quadratic')
        y_smooth=f(x_new)
        plt.plot(x_new,y_smooth,color='blue')
    except:
        pass
    plt.savefig(imagename)


def ScatterPlot1D(title,x,y1,xkey,ykey,labels,interpolate=True,correlation=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    x=np.array(x)
    y1=np.array(y1)
    plt.figure()
    ax = plt.axes()
    calc=ax.scatter(x,y1,color='r')
    imagename=title+'.png'
    ax.set_ylabel(ykey,fontsize=12)
    ax.set_xlabel(xkey,fontsize=12)
    plt.title(title)
    x=np.array(x)
    x_new = np.linspace(x.min(),x.max(),500)
    if interpolate==True:
        try:
            f = interp1d(x,y1, kind='quadratic')
            y_smooth=f(x_new)
            plt.plot(x_new,y_smooth,color='red')
        except:
            pass
    if correlation==True:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y1)
        rsquare=round(r_value*r_value,2)
        ax.text(0.1, .9,'$R^2$=%s'%(rsquare),horizontalalignment="center",verticalalignment="center",transform = ax.transAxes)
        ax.plot(np.array(x), slope*np.array(x) + intercept,'k-')

    plt.savefig(imagename)


def GrabFinalNeatLiquidTrajectories(jobdirs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    moltotptoarc={}
    moltomaxiter={}
    ext='.arc'
    curdir=os.getcwd()
    for i in range(len(jobdirs)):
        jobdir=jobdirs[i]
        os.chdir(jobdir)
        if os.path.isdir('optimize.tmp'):
            os.chdir('optimize.tmp')
            files=os.listdir()
            for f in files:
                if 'Liquid' in f:
                    os.chdir(f)
                    namesplit=f.split('_')
                    name=namesplit[1]
                    if name not in moltotptoarc.keys():
                        moltotptoarc[name]={}
                    newfiles=os.listdir()
                    numtofolder={}
                    for newf in newfiles:
                        if 'iter' in newf:
                            split=newf.split('_')
                            num=int(split[1])
                            numtofolder[num]=newf
                    maxnum=max(numtofolder.keys())
                    moltomaxiter[name]=maxnum
                    maxfolder=numtofolder[maxnum]
                    os.chdir(maxfolder)
                    morefiles=os.listdir()
                    for finalf in morefiles:
                        if os.path.isdir(finalf):
                            os.chdir(finalf)
                            evenmorefiles=os.listdir()
                            for last in evenmorefiles:
                                 if ext in last and 'gas' not in last:
                                     path=os.path.join(os.getcwd(),last)
                                     moltotptoarc[name][finalf]=path
 
                            os.chdir('..')


   


                    os.chdir('..') 
              

 
                    os.chdir('..')
    os.chdir(curdir)
    return moltotptoarc,moltomaxiter


def GrabDimerDistanceInfo(qmdiclist,jobdirs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    qmtargetnamedic={}
    nametodimerstructs={}
    truenametoindices={}
    for i in range(len(jobdirs)):
        jobdir=jobdirs[i]
        os.chdir(jobdir)
        qmdic=qmdiclist[i]
        if os.path.isdir('targets'):
            os.chdir('targets')
            for name in qmdic.keys():
                namedic=qmdic[name]
                os.chdir(name)
                temp=open('all.arc','r')
                results=temp.readlines()
                temp.close()
                indextoname={}
                dimernametotinkxyzlines={}
                tinkxyzlines=[]
                dimername=None
                num=None
                for line in results:
                    linesplit=line.split()
                    if len(linesplit)==3:
                        if dimername!=None:
                            if num==1:
                                dimernametotinkxyzlines[dimername]=tinkxyzlines           
                        tinkxyzlines=[]
                        tinkxyzlines.append(line)
                        dimername=linesplit[1]
                        prefix=dimername.replace('.xyz','')
                        prefixsplit=prefix.split('_')
                        num=float(prefixsplit[-1])
                        index=linesplit[2]
                        indextoname[index]=dimername
                    else:
                        tinkxyzlines.append(line)
                if num==1:
                    dimernametotinkxyzlines[dimername]=tinkxyzlines      
                split=name.split('_')
                truename='_'.join(split[1:-1]) 
                if truename not in nametodimerstructs.keys():
                    nametodimerstructs[truename]={}
                nametodimerstructs[truename].update(dimernametotinkxyzlines)
                indextotruename={}
                indextodistance={}
                for index,name in indextoname.items():
                    split=name.split('_')
                    truename='_'.join(split[:-1]) 
                    tail=split[-1]
                    distance=tail.replace('.xyz','')
                    distance=float(distance)
                    indextotruename[index]=truename
                    indextodistance[index]=distance
                truenametoindex={}
                for index,truename in indextotruename.items():
                    if truename not in truenametoindex.keys():
                        truenametoindex[truename]=[]
                    truenametoindex[truename].append(index)
                for truename,indices in truenametoindex.items():
                    distances=[indextodistance[k] for k in indices]
                    smallindicestodistance=dict(zip(indices,distances))
                    sortedindextodistance={k: v for k, v in sorted(smallindicestodistance.items(), key=lambda item: item[1])}
                    sortedindices=list(sortedindextodistance.keys())
                    sortedindices=[int(k) for k in sortedindices]
                    sorteddistances=list(sortedindextodistance.values())
                    #try:
                    valuedic=[namedic[k] for k in sortedindices]
                    refvalues=[k['Ref'] for k in valuedic]
                    calcvalues=[k['Calc'] for k in valuedic]
                    if truename not in qmtargetnamedic.keys():
                        qmtargetnamedic[truename]={}
                        split=truename.split('_')
                        indices=[split[-3],split[-2]]
                        indices=[int(i) for i in indices]

                        truenametoindices[truename]=indices
                    if 'Distances' not in qmtargetnamedic[truename].keys():
                        qmtargetnamedic[truename]['Distances']=sorteddistances
                    if 'Ref' not in qmtargetnamedic[truename].keys():
                        qmtargetnamedic[truename]['Ref']=refvalues
                    if 'Calc' not in qmtargetnamedic[truename].keys():
                        qmtargetnamedic[truename]['Calc']=calcvalues
                    #except:
                    #    pass
            
                os.chdir('..')
    return qmtargetnamedic,nametodimerstructs,truenametoindices 


def Chunks(lst, n):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ChunksList(gen):
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


def PlotAllESPSurfaces(nametocubefiles,pymolscript,python2path):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for name in nametocubefiles.keys():
        cubefiles=nametocubefiles[name]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        CallPlotESPSurfaces(name,cubefiles,pymolscript,python2path)
        os.chdir('..')



def PlotAllDimers(nametodimerstructs,truenametoindices,pymolscript,python2path):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    nametofilenametoformula={} 
    for name in nametodimerstructs.keys():
        dic=nametodimerstructs[name]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        filenames=[]
        filenametoformula={}
        count=0
        allindices=[]
        for dimername,tinkerxyzlines in dic.items():
            split=dimername.split('_')
            truename='_'.join(split[0:-1]) 
            indices=truenametoindices[truename]
            allindices.append(indices)
            WriteFile(dimername,tinkerxyzlines)
            newfilename=dimername.replace('.xyz','cartesian.xyz')
            ConvertTinkerXYZToCartesian(dimername,newfilename)
            monomer=GrabMonomer(newfilename)
            monomerformula=GetMonomerFormula(monomer)
            dimerformula=GetDimerFormula(monomerformula,newfilename)
            os.remove(monomer)
            os.remove(dimername)
            dimerformula+='_'+str(count)
            renamedfile=dimerformula+'.xyz'
            os.rename(newfilename,renamedfile) 
            filenames.append(renamedfile) 
            dimersplit=dimername.split('_')
            dimer='_'.join(dimersplit[:-1])
            filenametoformula[dimer]=dimerformula
            count+=1
        nametofilenametoformula[name]=filenametoformula
        CallPlotDimers3D(filenames,allindices,pymolscript,python2path)
        os.chdir('..')
    return nametofilenametoformula


def CallPlotDimers3D(filenames,allindices,pymolscript,python2path):
    filenamestring=','.join(filenames)
    allindicesstring=''
    for ls in allindices:
        ls=[str(i) for i in ls]
        innerstring='_'.join(ls)
        allindicesstring+=innerstring+ ','
    allindicesstring=allindicesstring[:-1]
    cmdstr=python2path+' '+pymolscript+' '+'--filenames='+filenamestring+' '+'--allindices='+allindicesstring
    print('cmdstr',cmdstr,flush=True)
    os.system(cmdstr)


def CallPlotESPSurfaces(name,cubefiles,pymolscript,python2path):
    cubefilestring=','.join(cubefiles)
    cmdstr=python2path+' '+pymolscript+' '+'--cubefiles='+cubefilestring+' '+'--name='+name
    print('cmdstr',cmdstr,flush=True)
    os.system(cmdstr)



def GetMonomerFormula(monomer):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    monomerformula=''
    temp=open(monomer,'r')
    results=temp.readlines()
    temp.close()
    elementcounts={}
    for line in results:
        linesplit=line.split()
        if len(linesplit)>1:
            element=linesplit[0]
            if element not in elementcounts.keys():
                elementcounts[element]=1
            else:
                elementcounts[element]+=1 
    sortedelementcounts=sorted(elementcounts.keys(), key=lambda x:x.lower())
    for element in sortedelementcounts:
        counts=elementcounts[element]
        total=str(counts)
        if counts==1:
            total=''
        monomerformula+=element+total

    return monomerformula


def GetDimerFormula(monomerformula,newfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    dimerformula=''
    if 'water' in newfilename:
        dimerformula+=monomerformula+'-'+'H2O'
    else:
        dimerformula+=monomerformula+'-'+monomerformula
    return dimerformula


def GrabMonomer(filename):
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
    indextolines={}
    count=0
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx==0:
            linesplit=line.split()
            atomnum=int(linesplit[0])

            if 'water' in filename:
                newatomnum=atomnum-3
            else:
                newatomnum=int(atomnum/2)
        elif lineidx==1:
            pass
        else:
            linesplit=line.split()
            if len(linesplit)>1:
                count+=1  
                if count<=newatomnum:
                    indextolines[count]=line
    monomer=filename.replace('.xyz','monomer.xyz')
    temp=open(monomer,'w')
    temp.write(str(newatomnum)+'\n') 
    temp.write('\n')
    for index,line in indextolines.items():
        temp.write(line)
    temp.close()
    return monomer


def WriteFile(dimername,tinkerxyzlines):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(dimername,'w')
    for line in tinkerxyzlines:
        temp.write(line)
    temp.close()



def ConvertTinkerXYZToCartesian(filename,newfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(filename,'r')
    tempwrite=open(newfilename,'w')
    results=temp.readlines()
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx==0:
            linesplit=line.split()
            tempwrite.write(linesplit[0]+'\n')
            tempwrite.write('\n')
            tempwrite.flush()
            os.fsync(tempwrite.fileno())
        else:
            linesplit=line.split()
            newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
            tempwrite.write(newline)
            tempwrite.flush()
            os.fsync(tempwrite.fileno())
    temp.close()
    tempwrite.close()



def SmallestDivisor(n):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    a=[]
    for i in range(2,n+1):
        if(n%i==0):
            a.append(i)
    a.sort()
    return a[0]




def ParseXYZ(xyzpath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(xyzpath,'r')
    results=temp.readlines()
    temp.close()
    indextoconn={}
    indextoelement={}
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            totalatomnum=int(linesplit[0])
        if len(linesplit)>1 and '90.00' not in line:
            index=int(linesplit[0])
            element=linesplit[1]
            indextoelement[index]=element
            conn=linesplit[6:]
            conn=[int(i) for i in conn]
            indextoconn[index]=conn
    groups=[]
    for index,conn in indextoconn.items():
        current=[index]+conn
        found=False
        for grpidx in range(len(groups)):
            breakout=False
            grp=groups[grpidx]
            for idx in current:
                if idx in grp: 
                    breakout=True
            if breakout==True:
                found=True
                break
        if found==True:
            for idx in current:
                if idx not in grp:
                    grp.append(idx)
            groups[grpidx]=grp

        else:
            groups.append(current) 
    numatomspermol=len(groups[0])
    nummolecules=len(groups)

    return nummolecules,numatomspermol,indextoelement

def PlotAllRDFs(moltotptoarc,neatliqnametoindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for mol,tptoarc in moltotptoarc.items():
        if not os.path.isdir(mol):
            os.mkdir(mol)
        os.chdir(mol)
        if not os.path.isdir('RDF'):
            os.mkdir('RDF')
        os.chdir('RDF')
        if mol in neatliqnametoindices.keys():
            indicesarray=neatliqnametoindices[mol]
            for tp,arcpath in tptoarc.items():
                t = md.load_arc(arcpath)
                xyzpath=arcpath.replace('.arc','.xyz')
                nummolecules,numatomspermol,indextoelement=ParseXYZ(xyzpath)
                pairs,indicesmats,dimerformula=GeneratePairs(nummolecules,numatomspermol,indextoelement,indicesarray)
                for i in range(len(pairs)):
                    indicesmat=indicesmats[i]
                    pair=pairs[i] 
                    PlotRDF(mol,tp,t,indicesmat,pair,dimerformula)
        os.chdir('..')
        os.chdir('..')


def GeneratePairs(nummolecules,numatomspermol,indextoelement,indicesarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    elementcounts={}
    for index,element in indextoelement.items():
        if index>numatomspermol:
            break
        if element not in elementcounts.keys():
            elementcounts[element]=1
        else:
            elementcounts[element]+=1
        
    sortedelementcounts=sorted(elementcounts.keys(), key=lambda x:x.lower())
    monomerformula=''
    for element in sortedelementcounts:
        counts=elementcounts[element]
        total=str(counts)
        if counts==1:
            total=''
        monomerformula+=element+total
    dimerformula=monomerformula+'-'+monomerformula
    pairs=[]
    indicesmats=[]
    for indices in indicesarray:
        first=indices[0]
        last=indices[1]
        mat=[]
        for i in range(nummolecules-1):
            newlast=last+numatomspermol 
            newpair=[first,newlast]
            mat.append(newpair)
        firstelement=indextoelement[first]
        lastelement=indextoelement[last]
        pair=firstelement+'..'+lastelement
        pairs.append(pair) 
        indicesmats.append(mat) 

    return pairs,indicesmats,dimerformula


def PlotRDF(mol,tp,t,indicesmat,pair,dimerformula):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    plt.figure()
    title=mol+' '+dimerformula+' '+pair+' '+tp+' '+'RDF'
    imagename=title+'.png'
    r_max = 1
    r_min = 0.01
    rdf=md.compute_rdf(t,indicesmat,(r_min, r_max))
    plt.plot(*rdf, color='black',label="mdtraj", alpha=0.5)
    ykey='RDF'
    xkey='Distance (nm)'+' '+pair
    plt.ylabel(ykey,fontsize=12)
    plt.xlabel(xkey,fontsize=12)
    plt.title(title)
    plt.show()
    plt.savefig(imagename)


def ExtractNeatLiquidIndices(truenametoindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    neatliqnametoindices={}
    for truename,indices in truenametoindices.items():
        if 'water' in truename:
            continue
        namesplit=truename.split('-')
        name='-'.join(namesplit[:-2])
        namesplit=name.split('-')
        if name.count('-')>=2:
            name='-'.join(namesplit[:-1])
        namesplit=name.split('_')
        name=namesplit[0]
        if name not in neatliqnametoindices.keys():
             neatliqnametoindices[name]=[]
        neatliqnametoindices[name].append(indices)
    return neatliqnametoindices


def WriteOutParamTable(moltotypetoprms,moltotypetoelement,moltotypetofitrad,moltotypetofitdep,moltotypetofitred):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempname='SummaryParams.csv'
    nametoallprmlines={}
    elementtoprmlines={}
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header=['Name','Type','Element','Radius','Depth','Reduction','Fit Radius','Fit Depth','Fit Reduction']
        energy_writer.writerow(header)
        for mol,typetoprms in moltotypetoprms.items():
            typetoelement=moltotypetoelement[mol]
            if mol in moltotypetofitrad.keys():
                typetofitrad=moltotypetofitrad[mol]
                typetofitdep=moltotypetofitdep[mol]
                typetofitred=moltotypetofitred[mol]
                for typenum,prms in typetoprms.items():
                    array=[]
                    element=typetoelement[typenum]
                    radius=round(prms[0],3)
                    depth=round(prms[1],3)
                    shrink=round(prms[2],3)
                    fitrad=typetofitrad[typenum]
                    fitdep=typetofitdep[typenum]
                    fitred=typetofitred[typenum]
                    array.append(mol)
                    array.append(typenum)
                    array.append(element)
                    array.append(radius)
                    array.append(depth)
                    array.append(shrink)
                    array.append(fitrad)
                    array.append(fitdep)
                    array.append(fitred)



                    energy_writer.writerow(array) 
                    if mol not in nametoallprmlines.keys():
                        nametoallprmlines[mol]=[]
                    nametoallprmlines[mol].append(array)

                    if element not in elementtoprmlines.keys():
                        elementtoprmlines[element]=[]
                    elementtoprmlines[element].append(array)
        energy_writer.writerow([])
        for element,prmlines in elementtoprmlines.items():
            for prmline in prmlines:
                energy_writer.writerow(prmline)

    return nametoallprmlines,elementtoprmlines,tempname



def PlotLiquidPropsCorrelation(nametotptofinalprops):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    for name,tptofinalprops in nametotptofinalprops.items():
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        n=len(list(tptofinalprops.keys()))
        color = iter(cm.rainbow(np.linspace(0, 1, n)))
        proptofigs={}
        proptoaxes={}
        propnametocalcs={}
        propnametorefs={}
        propnametotitle={}
        for tp,proptodic in tptofinalprops.items():
            string='T=%s,P=%s'%(str(tp[0]),str(tp[1]))
            c = next(color)
            for propname in proptodic.keys():
                if propname not in proptofigs.keys():
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    proptofigs[propname]=fig
                    proptoaxes[propname]=ax1

                if propname=='Enthalpy':
                    units='(kcal/mol)'
                else:
                    units='$kg/m^3$'
                xkey=propname+' '+'AMOEBA'+' '+units
                ykey=propname+' '+'QM'+' '+units
                innerdic=proptodic[propname]
                calcarray=innerdic['Calc']
                refarray=innerdic['Ref']
                calcerrarray=innerdic['CalcErr']   
                proptoaxes[propname].scatter(calcarray,refarray, s=10, c=[c], marker="o", label=string)
                proptoaxes[propname].errorbar(calcarray, refarray, yerr=calcerrarray,c=c, fmt="o")
                proptoaxes[propname].set_ylabel(ykey,fontsize=12)
                proptoaxes[propname].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+propname+ ' AMOEBA Vs QM'
                proptoaxes[propname].set_title(newtitle)
                proptoaxes[propname].legend()
                if propname not in propnametocalcs.keys():
                    propnametocalcs[propname]=[]
                    propnametorefs[propname]=[]
                propnametocalcs[propname].append(calcarray)
                propnametorefs[propname].append(refarray)
                propnametotitle[propname]=newtitle
        for propname,calcarray in propnametocalcs.items():
            refarray=propnametorefs[propname]
            newtitle=propnametotitle[propname]
            slope, intercept, r_value, p_value, std_err = stats.linregress(calcarray,refarray)
            rsquare=round(r_value*r_value,2)
            imagename=newtitle+'.png'
            newtitle+=' $R^2$=%s'%(rsquare)
            proptoaxes[propname].plot(np.array(calcarray), slope*np.array(calcarray) + intercept,'k-')

            proptoaxes[propname].set_title(newtitle)
            proptofigs[propname].savefig(imagename)
            rms = round(sqrt(mean_squared_error(calcarray,refarray)),3)
            print('RMSE for '+str(name)+' '+propname+' '+str(rms))

        

    os.chdir(curdir)

def PlotLiquidPropsVsTemp(nametotptofinalprops):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    nametotemparray={}
    nametoproptokeytoarray={}
    for name,tptofinalprops in nametotptofinalprops.items():
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        temparray=[]
        proptokeytoarray={}
        for tp, finalprops in tptofinalprops.items():
            temp=tp[0]
            pressure=tp[1]
            temparray.append(temp)
            for propname,dic in finalprops.items():
                if propname not in proptokeytoarray.keys():
                    proptokeytoarray[propname]={}
                for key,value in dic.items():
                    if key not in proptokeytoarray[propname].keys():
                        proptokeytoarray[propname][key]=[]
                    proptokeytoarray[propname][key].append(value)
        
        for propname in proptokeytoarray.keys():
            if propname=='Enthalpy':
                units='(kcal/mol)'
            else:
                units='$kg/m^3$'
            xkey='Temperature (K)'
            ykey=propname+' '+units

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            innerdic=proptokeytoarray[propname]
            calcarray=innerdic['Calc']
            refarray=innerdic['Ref']
            calcerrarray=innerdic['CalcErr']   
            ax1.scatter(temparray,calcarray, s=10, c='red', marker="o", label='AMOEBA')
            ax1.errorbar(temparray, calcarray, yerr=calcerrarray,c='red', fmt="o")
            ax1.scatter(temparray,refarray, s=10, c='blue', marker="s", label='Exp')              
            ax1.set_ylabel(ykey,fontsize=12)
            ax1.set_xlabel(xkey,fontsize=12)
            newtitle=name+' '+propname+ ' Vs T'
            ax1.set_title(newtitle)
            ax1.legend()
            x=np.array(temparray)
            refarray=np.array(refarray)
            calcarray=np.array(calcarray)
            x_new = np.linspace(x.min(),x.max(),500)
            imagename=newtitle+'.png'
            fig.savefig(imagename)


        nametoproptokeytoarray[name]=proptokeytoarray
        nametotemparray[name]=temparray

        os.chdir('..')
    os.chdir(curdir)

    return nametotemparray,nametoproptokeytoarray


def WriteOutPropTable(nametotptofinalprops,moltomaxiter,targetdensityerror,targetenthalpyerror):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempname='SummaryProps.csv'
    densityerrors=[]
    densityrelerrors=[]
    enthalpyerrors=[]
    enthalpyrelerrors=[]
    nametopropavgerrors={}
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header=['Name','Temperature','Pressure','Density Ref','Density Calc','Density CalcErr','Density Error','Density RelError','Density Direction','Enthalpy Ref','Enthalpy Calc','Enthalpy CalcErr','Enthalpy Error','Enthalpy RelError','Enthalpy Direction','Max Iter','Enthalpy Loss','Density Loss','Enthalpy Reached Target?','Density Reached Target?','Enthalpy Target','Density Target']
        energy_writer.writerow(header)
        for name,tptofinalprops in nametotptofinalprops.items():
            maxiter=moltomaxiter[name]
            namedensityerrors=[]
            namedensityrelerrors=[]
            nameenthalpyerrors=[]
            nameenthalpyrelerrors=[]

            for tp, finalprops in tptofinalprops.items():
                array=[0]*len(header)
                headeridx=header.index('Max Iter')
                array[headeridx]=maxiter
                temp=tp[0]
                pressure=tp[1]
                array[0]=name
                array[1]=temp
                array[2]=pressure
                for propname,dic in finalprops.items():
                    all=[]
                    for key,value in dic.items():
                        thename=propname+' '+key
                        headeridx=header.index(thename)
                        value=round(value,3)
                        array[headeridx]=value
                        all.append(value)
                        if key=='Ref':
                            truevalue=value
                        elif key=='Calc':
                            calcvalue=value
                    direrror=all[0]-all[1]
                    error=np.abs(round(direrror,2))
                    if (truevalue-calcvalue)>0:
                        direction='Underpredicted'
                    elif (truevalue-calcvalue)<0:
                        direction='Overpredicted'
                    relerror=round((100*error)/truevalue,2)
                    thename=propname+' '+'Error'
                    headeridx=header.index(thename)
                    array[headeridx]=error
                    thename=propname+' '+'RelError'
                    headeridx=header.index(thename)
                    array[headeridx]=relerror
                    reachedtargeterror=True
                    if propname=='Density':
                        densityerrors.append(error)
                        densityrelerrors.append(relerror)
                        namedensityerrors.append(error)
                        namedensityrelerrors.append(relerror)
                        targeterror=targetdensityerror
                        if error>targetdensityerror:
                            reachedtargeterror=False
                    elif propname=='Enthalpy':
                        enthalpyerrors.append(error)
                        enthalpyrelerrors.append(relerror)
                        nameenthalpyerrors.append(error)
                        nameenthalpyrelerrors.append(relerror)
                        targeterror=targetenthalpyerror
                        if error>targetenthalpyerror:
                            reachedtargeterror=False
                    thename=propname+' '+'Reached Target?'
                    headeridx=header.index(thename)
                    array[headeridx]=reachedtargeterror
                    thename=propname+' '+'Target'
                    headeridx=header.index(thename)
                    array[headeridx]=targeterror
                    thename=propname+' '+'Direction'
                    headeridx=header.index(thename)
                    array[headeridx]=direction


                energy_writer.writerow(array)
            nameavedensityerr=np.mean(np.array(namedensityerrors))
            nameavedensityrelerr=np.mean(np.array(namedensityrelerrors))
            nameaveenthalpyerr=np.mean(np.array(nameenthalpyerrors))
            nameaveenthalpyrelerr=np.mean(np.array(nameenthalpyrelerrors))
            if name not in nametopropavgerrors.keys():
                nametopropavgerrors[name]={}
            nametopropavgerrors[name]['DensityErr']=nameavedensityerr
            nametopropavgerrors[name]['DensityRelErr']=nameavedensityrelerr
            nametopropavgerrors[name]['EnthalpyErr']=nameaveenthalpyerr
            nametopropavgerrors[name]['EnthalpyRelErr']=nameaveenthalpyrelerr

    avedensityerr=round(np.mean(np.array(densityerrors)),3)
    avedensityrelerr=round(np.mean(np.array(densityrelerrors)),3)
    aveenthalpyerr=round(np.mean(np.array(enthalpyerrors)),3)
    aveenthalpyrelerr=round(np.mean(np.array(enthalpyrelerrors)),3)
    print('Average Density Error ',avedensityerr)
    print('Average Relative Density Error ',avedensityrelerr)
    print('Average Enthalpy Error ',aveenthalpyerr)
    print('Average Relative Enthalpy Error ',aveenthalpyrelerr)
    return nametopropavgerrors,tempname


def GrabFinalParameters(jobdirs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    prmfiles=[]
    curdir=os.getcwd()
    for jobdir in jobdirs:
        os.chdir(jobdir)
        if os.path.isdir('result'):
            os.chdir('result')
            if os.path.isdir('optimize'):
                os.chdir('optimize')
                files=os.listdir()
                numtoprmfilepath={}
                for f in files: 
                    if '.prm' in f and '_' in f:
                        split=f.split('_')
                        suffix=split[1].replace('.prm','')   
                        num=int(suffix)
                        prmfilepath=os.path.join(os.getcwd(),f)
                        numtoprmfilepath[num]=prmfilepath
                if len(numtoprmfilepath.keys())>0:
                    maxnum=max(numtoprmfilepath.keys())
                    maxprmfilepath=numtoprmfilepath[maxnum]
                    prmfiles.append(maxprmfilepath) 

    os.chdir(curdir)
    return prmfiles



def GrabMolFromTypes(moltotypes,typenum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for mol,types in moltotypes.items():
        for otypenum in types:
            if typenum==otypenum:
                return mol


def GrabParameterValues(prmfiles):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    moltotypetoprms={}
    moltotypetoelement={}
    moltotypetofitrad={}
    moltotypetofitdep={}
    moltotypetofitred={}
    for prmfile in prmfiles:
        temp=open(prmfile,'r')
        results=temp.readlines()
        temp.close()
        foundatomblock=False
        typetoprms={}
        typetoelement={}
        moltotypes={}
        passedatomblock=False
        for line in results:
            linesplit=line.split()
            if len(linesplit)>1:
                first=linesplit[0]
                if first=='atom' and passedatomblock==False and 'AMOEBAWater' not in line:
                    foundatomblock=True
                    typenum=linesplit[1]
                    element=linesplit[3]
                    typetoelement[typenum]=element
                    mol=linesplit[4].replace('"','')
                    mol=mol.replace('_3D','')
                    if mol not in moltotypes.keys():
                        moltotypes[mol]=[]

                    moltotypes[mol].append(typenum)     

                if 'forcefield' in line:
                    passedatomblock=True
        foundatomblock=False
        passedatomblock=False
        for line in results:
            linesplit=line.split()
            if len(linesplit)>1:
                first=linesplit[0]
                if foundatomblock==False and first=='vdw':
                    typenum=linesplit[1]
                    mol=GrabMolFromTypes(moltotypes,typenum)
                    if mol not in moltotypetofitrad.keys():
                        moltotypetofitrad[mol]={}
                        moltotypetofitdep[mol]={}
                        moltotypetofitred[mol]={}
                    fitarray=[]
                    if 'PRM' in line:
                        foundprm=False
                        appendarray=True
                        for eidx in range(len(linesplit)):
                            e=linesplit[eidx]
                            if e=='PRM':
                                foundprm=True
                            if foundprm==True and e!='PRM':
                                if e.isdigit() and appendarray==True:
                                    fitarray.append(int(e))
                                else:
                                    appendarray=False
                    if 2 in fitarray:
                        moltotypetofitrad[mol][typenum]=True
                    else:
                        moltotypetofitrad[mol][typenum]=False

                    if 3 in fitarray:
                        moltotypetofitdep[mol][typenum]=True
                    else:
                        moltotypetofitdep[mol][typenum]=False


                    if 4 in fitarray:
                        moltotypetofitred[mol][typenum]=True
                    else:
                        moltotypetofitred[mol][typenum]=False

                if first=='vdw':
                    typenum=linesplit[1]
                    mol=GrabMolFromTypes(moltotypes,typenum)
                    radius=float(linesplit[2])
                    depth=float(linesplit[3])
                    try:
                        shift=float(linesplit[4])
                    except:
                        shift=1
                    typetoprms[typenum]=[radius,depth,shift]
                    if typenum not in moltotypetofitrad[mol].keys():
                        moltotypetofitrad[mol][typenum]=False
                    if typenum not in moltotypetofitdep[mol].keys():
                        moltotypetofitdep[mol][typenum]=False
                    if typenum not in moltotypetofitred[mol].keys():
                        moltotypetofitred[mol][typenum]=False




        for mol,types in moltotypes.items():
            if mol not in moltotypetoprms.keys():
                moltotypetoprms[mol]={}
                moltotypetoelement[mol]={}
            for typenum in types: 
                if typenum in typetoprms.keys():
                    prms=typetoprms[typenum]
                    moltotypetoprms[mol][typenum]=prms
                    moltotypetoelement[mol][typenum]=typetoelement[typenum] 

    return moltotypetoprms,moltotypetoelement,moltotypetofitrad,moltotypetofitdep,moltotypetofitred


def WriteOutQMTable(nametoformulatormse,nametoformulatomse):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempname='SummaryQM.csv'
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header=['Name','Dimer Name','RMSE','MSE']
        energy_writer.writerow(header)
        for name,formulatormse in nametoformulatormse.items():
            for formula,rmse in formulatormse.items():
                mse=nametoformulatomse[name][formula]
                array=[name,formula,round(rmse,2),round(mse,2)]
                energy_writer.writerow(array)
    return tempname

                    
def PlotLiquidPropsCorrelationForGroups(nametotptofinalprops,groupednames,nametotemparray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    grptoproptofig={}
    grptoproptoaxes={}
    grptopropnametocalc={}
    grptopropnametorefs={}
    grptopropnametotitle={}
    for grp in groupednames:
        n=0
        if len(grp)>1:
            grp=tuple(grp)
            for name in grp:
                if name in nametotptofinalprops.keys():
                    tptofinalprops=nametotptofinalprops[name]
                    n+=len(list(tptofinalprops.keys()))
            color = cm.rainbow(np.linspace(0, 1, n))
            if grp not in grptoproptofig.keys():
                grptoproptofig[grp]={}
                grptoproptoaxes[grp]={}
            totalname=''
            for name in grp:
                if name in nametotemparray.keys():
                    totalname+=name+','
            totalname=totalname[:-1]
            count=0
            for nameidx in range(len(grp)):
                name=grp[nameidx]
                string=' '+name
                if name in nametotptofinalprops.keys():
                    tptofinalprops=nametotptofinalprops[name]
                    for tp,proptodic in tptofinalprops.items():
                        newstring=string+' T=%s'%(str(tp[0]))
                        c=color[count]
                        count+=1
                        for propname in proptodic.keys():
                            if propname not in grptoproptofig[grp].keys():
                                 fig = plt.figure()
                                 ax1 = fig.add_subplot(111)
                                 grptoproptofig[grp][propname]=fig
                                 grptoproptoaxes[grp][propname]=ax1
                            if propname=='Enthalpy':
                                units='(kcal/mol)'
                            else:
                                units='$kg/m^3$'
                            xkey=propname+' '+'AMOEBA'+' '+units
                            ykey=propname+' '+'Exp'+' '+units
                            innerdic=proptodic[propname]
                            calcarray=innerdic['Calc']
                            refarray=innerdic['Ref']
                            calcerrarray=innerdic['CalcErr']   
                            grptoproptoaxes[grp][propname].scatter(calcarray,refarray, s=10, c=[c], marker="o", label=newstring)
                            grptoproptoaxes[grp][propname].errorbar(calcarray, refarray, yerr=calcerrarray,c=c, fmt="o")
                            grptoproptoaxes[grp][propname].set_ylabel(ykey,fontsize=12)
                            grptoproptoaxes[grp][propname].set_xlabel(xkey,fontsize=12)
                            newtitle=totalname+' '+propname+ ' AMOEBA Vs Exp'
                            grptoproptoaxes[grp][propname].set_title(newtitle)
                            grptoproptoaxes[grp][propname].legend()
                            if grp not in grptopropnametocalc.keys():
                                grptopropnametocalc[grp]={}
                                grptopropnametorefs[grp]={}
                                grptopropnametotitle[grp]={}
                            if propname not in grptopropnametocalc[grp].keys():
                                grptopropnametocalc[grp][propname]=[]
                                grptopropnametorefs[grp][propname]=[]
                                grptopropnametotitle[grp][propname]={}
                            grptopropnametocalc[grp][propname].append(calcarray)
                            grptopropnametorefs[grp][propname].append(refarray)
                            grptopropnametotitle[grp][propname]=newtitle
    for grp,propnametocalc in grptopropnametocalc.items():
        propnametorefs=grptopropnametorefs[grp] 
        propnametotitle=grptopropnametotitle[grp]
        for propname,calcarray in propnametocalc.items():
            refarray=propnametorefs[propname]
            title=propnametotitle[propname]
            slope, intercept, r_value, p_value, std_err = stats.linregress(calcarray,refarray)
            rsquare=round(r_value*r_value,5)
            newtitle=title+' $R^2$=%s'%(rsquare)
            grptoproptoaxes[grp][propname].plot(np.array(calcarray), slope*np.array(calcarray) + intercept,'k-')
            grptoproptoaxes[grp][propname].set_title(newtitle)
            imagename=title+'.png'
            grptoproptofig[grp][propname].savefig(imagename)
            rms = round(sqrt(mean_squared_error(calcarray,refarray)),3)
            print('RMSE for '+str(grp)+' '+propname+' '+str(rms))



 
def PlotLiquidPropsVsTempForGroups(nametotemparray,nametoproptokeytoarray,groupednames):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    grptoproptofig={}
    grptoproptoaxes={}
    for grp in groupednames:
        
        if len(grp)>1:
            grp=tuple(grp)
            n=len(grp)
            color = cm.rainbow(np.linspace(0, 1, n))
            if grp not in grptoproptofig.keys():
                grptoproptofig[grp]={}
                grptoproptoaxes[grp]={}
            totalname=''
            for name in grp:
                if name in nametotemparray.keys():
                    totalname+=name+','
            totalname=totalname[:-1]
            for nameidx in range(len(grp)):
                name=grp[nameidx]
                if name in nametotemparray.keys():
                    temparray=nametotemparray[name]
                    string=' '+name        
                    proptokeytoarray=nametoproptokeytoarray[name]
                    c=color[nameidx]
                    for propname in proptokeytoarray.keys():
                        if propname not in grptoproptofig[grp].keys():
                            fig = plt.figure()
                            ax1 = fig.add_subplot(111)
                            grptoproptofig[grp][propname]=fig
                            grptoproptoaxes[grp][propname]=ax1

                        if propname=='Enthalpy':
                            units='(kcal/mol)'
                        else:
                            units='$kg/m^3$'
                        xkey='Temperature (K)'
                        ykey=propname+' '+units

                        innerdic=proptokeytoarray[propname]
                        calcarray=innerdic['Calc']
                        refarray=innerdic['Ref']
                        calcerrarray=innerdic['CalcErr']   
                        grptoproptoaxes[grp][propname].scatter(temparray,calcarray, s=10, c=[c], marker="o", label='AMOEBA'+string)
                        grptoproptoaxes[grp][propname].errorbar(temparray, calcarray, yerr=calcerrarray,c=c, fmt="o")
                        grptoproptoaxes[grp][propname].scatter(temparray,refarray, s=10, c=[c], marker="s", label='Exp'+string)              
                        grptoproptoaxes[grp][propname].set_ylabel(ykey,fontsize=12)
                        grptoproptoaxes[grp][propname].set_xlabel(xkey,fontsize=12)
                        newtitle=totalname+' '+propname+ ' Vs T'
                        grptoproptoaxes[grp][propname].set_title(newtitle)
                        grptoproptoaxes[grp][propname].legend()
                        imagename=newtitle+'.png'
                        grptoproptofig[grp][propname].savefig(imagename)


def ComputeQMAverages(nametoformulatormse):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    rmsehomodimers=[]
    rmseheterodimers=[]
    rmsealldimers=[]
    nametodimerqmerror={}
    for name,formulatormse in nametoformulatormse.items():
        namermsehomodimers=[]
        namermseheterodimers=[]
        namermsealldimers=[]

        for formula,rmse in formulatormse.items():
            homo=True
            if 'H2O' in formula:
                homo=False
            if homo==True:
                rmsehomodimers.append(rmse)
                namermsehomodimers.append(rmse)
            else:
                rmseheterodimers.append(rmse)
                namermseheterodimers.append(rmse)
            rmsealldimers.append(rmse)
            namermsealldimers.append(rmse)
        if len(namermsehomodimers)!=0:
            namehomoaverage=np.mean(np.array(namermsehomodimers))
        else:
            namehomoaverage=0
        nameheteroaverage=np.mean(np.array(namermseheterodimers))
        nameallaverage=np.mean(np.array(namermsealldimers))
        if name not in nametodimerqmerror.keys():
            nametodimerqmerror[name]={}
        nametodimerqmerror[name]['homo']=namehomoaverage
        nametodimerqmerror[name]['hetero']=nameheteroaverage
        nametodimerqmerror[name]['all']=nameallaverage
    homoaverage=round(np.mean(np.array(rmsehomodimers)),3)
    heteroaverage=round(np.mean(np.array(rmseheterodimers)),3)
    allaverage=round(np.mean(np.array(rmsealldimers)),3)
    print('Homodimer average RMSE ',homoaverage)
    print('Heterodimer average RMSE ',heteroaverage)
    print('Dimer average RMSE ',allaverage)
    return nametodimerqmerror

def WriteOutParamAndAvgErrorTable(nametopropavgerrors,nametodimerqmerror,nametoallprmlines):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempname='SummaryParamsAvgPropErrorQMError.csv'
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header=['Name','Element','Rad','Dep','Red','Homo RMSE','Hetero RMSE','QM RMSE','Density Err','Density RelErr','Enthalpy Err','Enthalpy RelErr']
        energy_writer.writerow(header)
        for name,prmlines in nametoallprmlines.items():
            homoavg=0
            heteroavg=0
            allavg=0
            if name in nametodimerqmerror.keys():
                dimerqmerrors=nametodimerqmerror[name]
                homoavg=round(dimerqmerrors['homo'],2)
                heteroavg=round(dimerqmerrors['hetero'],2)
                allavg=round(dimerqmerrors['all'],2)
            newprmlines=[]
            for array in prmlines:
                array.append(homoavg)
                array.append(heteroavg)
                array.append(allavg)
                newprmlines.append(array)
            finalprmlines=[]
            if name in nametopropavgerrors.keys():
                propavgerrors=nametopropavgerrors[name]
                densityavgerr=round(propavgerrors['DensityErr'],2)
                densityavgrelerr=round(propavgerrors['DensityRelErr'],2)
                enthalpyavgerr=round(propavgerrors['EnthalpyErr'],2)
                enthalpyavgrelerr=round(propavgerrors['EnthalpyRelErr'],2)
                for array in newprmlines:
                    array.append(densityavgerr)
                    array.append(densityavgrelerr)
                    array.append(enthalpyavgerr)
                    array.append(enthalpyavgrelerr)
                    finalprmlines.append(array)
            else:
                for array in newprmlines:
                    finalprmlines.append(array)

            for array in finalprmlines: 
                energy_writer.writerow(array)



def FindSameParameters(prmcsvname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    data = pd.read_csv(prmcsvname, quotechar='"', skipinitialspace=True)
    radiuscol=data.loc[:,['Radius']]
    depthcol=data.loc[:,['Depth']]
    redcol=data.loc[:,['Reduction']]
    raddepred=data.loc[:,['Radius','Depth','Reduction']]
    dataslicetosameindices=FindSameParametersSingleCol(raddepred)
    raddatatosameindices=FindSameParametersSingleCol(radiuscol)
    depdatatosameindices=FindSameParametersSingleCol(depthcol)
    reddatatosameindices=FindSameParametersSingleCol(redcol)
    ls=[raddatatosameindices,depdatatosameindices,reddatatosameindices]
    for dic in ls:
        for value,loc in dic.items():
            thevalue=value[0]
            alreadyhave=False
            for ovalue,oloc in dataslicetosameindices.items():
                if thevalue in ovalue:
                    for locidx in loc:
                        if locidx in oloc:
                            alreadyhave=True
            if alreadyhave==False:
                dataslicetosameindices[value]=loc
    return data,dataslicetosameindices


def FindSameParametersSingleCol(col):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    datatosameindices={}
    values=[]
    for i in range(len(col)):
        valuedf=col.iloc[i]
        temp=[]
        for k in range(len(valuedf)):
            value=valuedf.iloc[k]
            temp.append(value)
        if temp not in values:
            values.append(temp)
    nparray=col.to_numpy()
    for value in values:
        indices=np.where((nparray==np.array(value)).all(axis=1))[0]
        if len(indices)>2: # duplicate tables at bottom but sorted differently
            datatosameindices[tuple(value)]=list(indices)


    return datatosameindices



def GenerateIndicesToColor(data,dataslicetosameindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    rowindextocolumnindices={}
    rowindextocolorindex={}
    header=list(data.head())
    radidx=header.index('Radius')
    depidx=header.index('Depth')
    redidx=header.index('Reduction')
    ls=[radidx,depidx,redidx]
    cidx=0
    for dataslice,sameindices in dataslicetosameindices.items():
        colindices=tuple(ls[:len(dataslice)])
        for rowindex in sameindices:
            if rowindex not in rowindextocolumnindices.keys():
                rowindextocolumnindices[rowindex]=colindices
            if rowindex not in rowindextocolorindex.keys():
                rowindextocolorindex[rowindex]=cidx
        cidx+=1
    return rowindextocolumnindices,rowindextocolorindex

def MakeExcelFile(propcsvname,prmcsvname,qmcsvname,rowindextocolumnindices,rowindextocolorindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    from openpyxl.styles import Color, PatternFill, Font, Border
    from openpyxl.styles import colors
    from openpyxl.cell import Cell
    spreadsheet='ParametersQMPropsTable.xlsx' 
    df = pd.read_csv(prmcsvname)
    df2=pd.read_csv(propcsvname)
    df3=pd.read_csv(qmcsvname)
    colors = (
    '00000000', '00FFFFFF', '00FF0000', '0000FF00', '000000FF', #0-4
    '00FFFF00', '00FF00FF', '0000FFFF', '00000000', '00FFFFFF', #5-9
    '00FF0000', '0000FF00', '000000FF', '00FFFF00', '00FF00FF', #10-14
    '0000FFFF', '00800000', '00008000', '00000080', '00808000', #15-19
    '00800080', '00008080', '00C0C0C0', '00808080', '009999FF', #20-24
    '00993366', '00FFFFCC', '00CCFFFF', '00660066', '00FF8080', #25-29
    '000066CC', '00CCCCFF', '00000080', '00FF00FF', '00FFFF00', #30-34
    '0000FFFF', '00800080', '00800000', '00008080', '000000FF', #35-39
    '0000CCFF', '00CCFFFF', '00CCFFCC', '00FFFF99', '0099CCFF', #40-44
    '00FF99CC', '00CC99FF', '00FFCC99', '003366FF', '0033CCCC', #45-49
    '0099CC00', '00FFCC00', '00FF9900', '00FF6600', '00666699', #50-54
    '00969696', '00003366', '00339966', '00003300', '00333300', #55-59
    '00993300', '00993366', '00333399', '00333333',  #60-63
    )
    
    no_of_colors=len(set(list(rowindextocolorindex.values())))
    colorarray=[]
    count=0
    for col in colors:
        if col!='00FFFFFF' and col!='00000000' and col!='00FFFF00' and count<no_of_colors:
            colorarray.append(col)
            count+=1

 
    with pd.ExcelWriter(spreadsheet, engine="openpyxl") as writer:
        sheet_name = "Parameters"
        df.to_excel(writer, sheet_name=sheet_name)
        sheet = writer.sheets[sheet_name]
        for rowindex,columnindices in rowindextocolumnindices.items():
            colorindex=rowindextocolorindex[rowindex]
            col=colorarray[colorindex]
            for idx in columnindices:
                alpha=chr(ord('@')+idx+2)
                cellnum=alpha+str(rowindex+1+1)
                cell=sheet[cellnum]
                cell.font = Font(color=col)   
        sheet_name = "Properties"
        df2.to_excel(writer, sheet_name=sheet_name)
        sheet_name = "QM"
        df3.to_excel(writer, sheet_name=sheet_name)



def GenerateImagePropGrid(listofmoleculenametoproptoimagenames):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for moleculenametoproptoimagename in listofmoleculenametoproptoimagenames:
        proptoimages={}
        molnames=[]
        for moleculename,proptoimage in moleculenametoproptoimagename.items():
            molnames.append(moleculename)
            for prop,image in proptoimage.items():
                if prop not in proptoimages.keys():
                    proptoimages[prop]=[]
                proptoimages[prop].append(image)
        for prop,images in proptoimages.items():
            WriteOutGridImage(molnames,images,prop)


def WriteOutGridImage(moleculenames,filenamearray,prop):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    molsPerImage=len(filenamearray)
    if (molsPerImage % 2) == 0 or (molsPerImage ** 0.5) % 1==0:
        n=molsPerImage
    else:
        n=molsPerImage+1 
    if n==2 or n==1:
        molsperrow=1
    else:
        molsperrow=SmallestDivisor(n)
    ximagesize=640
    yimagesize=480
    filenamechunks=ChunksList(Chunks(filenamearray,molsPerImage))
    prevmatslen=len(filenamechunks[0])
    for i in range(len(filenamechunks)):
        filenamesublist=filenamechunks[i]
        imagenames=[]
        for j in range(len(filenamesublist)):
            filename=filenamesublist[j]
            ls=range(len(filenamesublist))
            chunks=ChunksList(Chunks(ls,molsperrow))
            indextorow={}
            for rowidx in range(len(chunks)):
                row=chunks[rowidx]
                for j in row:
                    indextorow[j]=rowidx
            
            fileprefix=filename.split('.')[0]
            imagename=fileprefix+'.png'
            imagenames.append(imagename)
            
            
        if i>0:
            factor=1
        else:
            factor=0
        firstj=i*prevmatslen+factor
        secondj=firstj+len(filenamesublist)-factor
        prevmatslen=len(filenamesublist)
         
        basename=','.join(moleculenames)+'_'+prop
        indextoimage={}
        for index in range(len(filenamesublist)):
            imagename=imagenames[index]
            image=Image.open(imagename)
            indextoimage[index]=image
        
        cols=len(set(list(indextorow.values())))
        dest = Image.new('RGB', (ximagesize*molsperrow,yimagesize*cols))
        for j in range(len(filenamesublist)):
            row=indextorow[j]
            x=(j-molsperrow*(row))*ximagesize
            y=(row)*yimagesize
            dest.paste(indextoimage[j],(x,y))
        dest.show()
        finalname=basename+'.png'
        dest.save(finalname)
        newname=os.path.join(curdir,finalname)




def PlotForceBalanceResults(fbdir,targetdensityerror,targetenthalpyerror,pymolscript,python2path):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    jobdirs=GrabJobDirectories(fbdir)
    outputfiles=GrabOutputFiles(jobdirs)
    outputfilestojobdir=dict(zip(outputfiles,jobdirs))
    tpdiclist,qmdiclist,newoutputfiles=GrabResultsAllMolecules(outputfiles)
    jobdirs=[outputfilestojobdir[i] for i in newoutputfiles]
    qmtargetnamedic,nametodimerstructs,truenametoindices=GrabDimerDistanceInfo(qmdiclist,jobdirs)
    neatliqnametoindices=ExtractNeatLiquidIndices(truenametoindices)
    os.chdir(curdir)
    nametofilenametoformula=PlotAllDimers(nametodimerstructs,truenametoindices,pymolscript,python2path)
    nametotptofinalprops,nametoformulatormse,nametoformulatomse,listofmoleculenametoproptoimagenames=PlotAllFBJobs(tpdiclist,qmtargetnamedic,nametofilenametoformula,truenametoindices)
    moltotptoarc,moltomaxiter=GrabFinalNeatLiquidTrajectories(jobdirs)
    nametopropavgerrors,propcsvname=WriteOutPropTable(nametotptofinalprops,moltomaxiter,targetdensityerror,targetenthalpyerror)
    PlotAllRDFs(moltotptoarc,neatliqnametoindices)
    prmfiles=GrabFinalParameters(jobdirs)
    moltotypetoprms,moltotypetoelement,moltotypetofitrad,moltotypetofitdep,moltotypetofitred=GrabParameterValues(prmfiles)
    nametoallprmlines,elementtoprmlines,prmcsvname=WriteOutParamTable(moltotypetoprms,moltotypetoelement,moltotypetofitrad,moltotypetofitdep,moltotypetofitred)
    data,dataslicetosameindices=FindSameParameters(prmcsvname)
    rowindextocolumnindices,rowindextocolorindex=GenerateIndicesToColor(data,dataslicetosameindices) 
    qmcsvname=WriteOutQMTable(nametoformulatormse,nametoformulatomse)
    poltypedirs,groupedpoltypedirs=GrabPoltypeDirectories(jobdirs)
    nametocubefiles,groupednames=GrabCubeFiles(groupedpoltypedirs)
    PlotAllESPSurfaces(nametocubefiles,pymolscript,python2path)
    nametotemparray,nametoproptokeytoarray=PlotLiquidPropsVsTemp(nametotptofinalprops)
    PlotLiquidPropsCorrelation(nametotptofinalprops)
    PlotLiquidPropsCorrelationForGroups(nametotptofinalprops,groupednames,nametotemparray)
    PlotLiquidPropsVsTempForGroups(nametotemparray,nametoproptokeytoarray,groupednames)
    nametodimerqmerror=ComputeQMAverages(nametoformulatormse)
    WriteOutParamAndAvgErrorTable(nametopropavgerrors,nametodimerqmerror,nametoallprmlines)
    MakeExcelFile(propcsvname,prmcsvname,qmcsvname,rowindextocolumnindices,rowindextocolorindex)
    GenerateImagePropGrid(listofmoleculenametoproptoimagenames) 
