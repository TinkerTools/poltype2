import os
import sys
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
from pymol import cmd,preset,util
from PIL import Image
from pymol.vfont import plain
from pymol.cgo import CYLINDER,cyl_text



fbdir=sys.argv[1]


def GrabJobDirectories(fbdir):
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
                prevprevprevline=results[lineidx-3]
                prevprevprevlinesplit=prevprevprevline.split()
                foundQM=False
                foundliq=False
                if 'Enthalpy' in prevprevline:
                   index=prevprevlinesplit.index('Enthalpy')
                   totalname=prevprevlinesplit[index-1]
                   namesplit=totalname.split('_')
                   name=namesplit[1]
                   property='Enthalpy'
                   foundliq=True
                elif 'Density' in prevprevline:
                   property='Density'
                   index=prevprevlinesplit.index('Density')
                   totalname=prevprevlinesplit[index-1]
                   namesplit=totalname.split('_')
                   name=namesplit[1]
                   foundliq=True
                elif 'Interaction Energies' in prevprevline:
                   property='Interaction Energies'
                   index=prevprevprevlinesplit.index('Target:')
                   targetname=prevprevprevlinesplit[index+1]
                   namesplit=targetname.split('_')
                   name=namesplit[1]
                   foundQM=True
                continue
        if '--------------------------------------------------' in line:
            foundblock=False 
        if foundblock==True:
            linesplit=line.split()
            if foundliq==True:
                if name not in tpdic.keys():
                    tpdic[name]={}
                temp=float(linesplit[0])
                pressure=float(linesplit[1])
                tp=tuple([temp,pressure])
                if tp not in tpdic[name].keys():
                    tpdic[name][tp]={}
                if property not in tpdic[name][tp].keys():
                    tpdic[name][tp][property]={}
                if refkeyword not in tpdic[name][tp][property].keys():
                    tpdic[name][tp][property][refkeyword]=[]
                if calckeyword not in tpdic[name][tp][property].keys():
                    tpdic[name][tp][property][calckeyword]=[]
                if errorkeyword not in tpdic[name][tp][property].keys():
                    tpdic[name][tp][property][errorkeyword]=[]

                ref=float(linesplit[3])
                calc=float(linesplit[4])
                error=float(linesplit[6])
                tpdic[name][tp][property][refkeyword].append(ref)
                tpdic[name][tp][property][errorkeyword].append(error)
                tpdic[name][tp][property][calckeyword].append(calc)
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
    nametotptofinalprops={}
    nametoformulatormse={}
    nametoformulatomse={}
    for tpdic in tpdiclist:
        nametotptofinalprops=PlotFBLiq(tpdic,nametotptofinalprops)
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
            nametoaxes[name].scatter(x,rmsvalues, s=10, c=c,marker="s", label=formula)
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

    return nametotptofinalprops,nametoformulatormse,nametoformulatomse

def GrabStructure(formula):
    files=os.listdir()
    for f in files:
        if formula in f and '.xyz' in f and 'tinkermincart' not in f:
            return f

def GrabMolecule(newname):
    obConversion = openbabel.OBConversion()
    newmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(newname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(newmol,newname)
    return newmol


def GrabDistanceAndElements(mol,indices):
    atoms=[mol.GetAtom(i) for i in indices]
    coords=[np.array([atom.GetX(),atom.GetY(),atom.GetZ()]) for atom in atoms]
    distance=np.linalg.norm(coords[0]-coords[1])
    atomicnums=[atom.GetAtomicNum() for atom in atoms]
    an = pyasl.AtomicNo()
    elements=[an.getElSymbol(atomicnum) for atomicnum in atomicnums]
    prefix='-'.join(elements)
    return distance,prefix

def PlotFBQM(targetdic,formula,indices):
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
        label2='Target'
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
            for property in innerdic.keys():
                relprop=property+'RelativeError'
                if property not in nametotptofinalprops[name][tp].keys():
                    nametotptofinalprops[name][tp][property]={}
                if property not in proptofigs.keys():
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    proptofigs[property]=fig
                    proptoaxes[property]=ax1
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    proptofigs[relprop]=fig
                    proptoaxes[relprop]=ax1


                innermostdic=innerdic[property]
                refkey='Ref'
                calckey='Calc'
                errkey='CalcErr' 
                refarray=np.array(innermostdic[refkey])
                calcarray=np.array(innermostdic[calckey])
                errarray=np.array(innermostdic[errkey])
                if property=='Enthalpy':
                    scale=0.239006 # convert kj to kcal
                    units='(kcal/mol)'
                else:
                    scale=1
                    units='$kg/m^3$'
                refarray=scale*refarray
                calcarray=scale*calcarray
                nametotptofinalprops[name][tp][property][refkey]=refarray[-1]
                nametotptofinalprops[name][tp][property][calckey]=calcarray[-1]
                nametotptofinalprops[name][tp][property][errkey]=errarray[-1]
                err=np.abs(refarray-calcarray)
                relerr=100*(err/refarray)
                label1=property+' AMOEBA'
                label2=property+' Target'
                ykey=property+units
                xkey='Iterations'
                labels=[label1,label2]
                string='T=%s,P=%s'%(str(tp[0]),str(tp[1]))
                title=name+' '+property+' '+string
                x=list(range(len(refarray)))
                #ScatterPlot2D(title,x,calcarray,refarray,xkey,ykey,labels)
                proptoaxes[property].scatter(x,calcarray, s=10, c=c, marker="o", label='AMOEBA '+string)
                proptoaxes[property].errorbar(x, calcarray, yerr=errarray,c=c, fmt="o")
                proptoaxes[property].scatter(x,refarray, s=10, c=c, marker="s", label='Target '+string)              
                proptoaxes[property].set_ylabel(ykey,fontsize=12)
                proptoaxes[property].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+property
                proptoaxes[property].set_title(newtitle)
                proptoaxes[property].legend()
                x=np.array(x)
                refarray=np.array(refarray)
                calcarray=np.array(calcarray)
                x_new = np.linspace(x.min(),x.max(),500)
                try:
                    f = interp1d(x,calcarray, kind='quadratic')
                    y_smooth=f(x_new)
                    proptoaxes[property].plot(x_new,y_smooth,color=c)
                    f = interp1d(x,refarray, kind='quadratic')
                    y_smooth=f(x_new)
                    proptoaxes[property].plot(x_new,y_smooth,color=c)
                except:
                    pass
                imagename=newtitle+'.png'
                proptofigs[property].savefig(imagename)
                proptoaxes[relprop].scatter(x,relerr, s=10, c=c, marker="s", label=string)
                try:
                    f = interp1d(x,relerr, kind='quadratic')
                    y_smooth=f(x_new)
                    proptoaxes[relprop].plot(x_new,y_smooth,color=c)
                except:
                    pass
                ykey='Relative Error '+property+units
                proptoaxes[relprop].set_ylabel(ykey,fontsize=12)
                proptoaxes[relprop].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+'Relative Error '+property
                proptoaxes[relprop].set_title(newtitle)
                proptoaxes[relprop].legend()
                imagename=newtitle+'.png'
                proptofigs[relprop].savefig(imagename)

 
        os.chdir('..')
    return nametotptofinalprops

def ScatterPlot2D(title,x,y1,y2,xkey,ykey,labels):
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
                    try:
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
                    except:
                        pass
            
                os.chdir('..')

    return qmtargetnamedic,nametodimerstructs,truenametoindices 


def Chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ChunksList(gen):
    newlst=[]
    for item in gen:
        newlst.append(item)
    return newlst


def PlotAllESPSurfaces(nametocubefiles):
    for name in nametocubefiles.keys():
        cubefiles=nametocubefiles[name]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        PlotESPSurfaces(name,cubefiles)
        os.chdir('..')



def PlotAllDimers(nametodimerstructs,truenametoindices):
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
        PlotDimers3D(filenames,allindices)
        os.chdir('..')
    return nametofilenametoformula


def GetMonomerFormula(monomer):
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
    dimerformula=''
    if 'water' in newfilename:
        dimerformula+=monomerformula+'-'+'H2O'
    else:
        dimerformula+=monomerformula+'-'+monomerformula
    return dimerformula


def GrabMonomer(filename):
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
    temp=open(dimername,'w')
    for line in tinkerxyzlines:
        temp.write(line)
    temp.close()



def ConvertTinkerXYZToCartesian(filename,newfilename):
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
    a=[]
    for i in range(2,n+1):
        if(n%i==0):
            a.append(i)
    a.sort()
    return a[0]


def PlotESPSurfaces(name,cubefiles):
    imagesize=1180
    dpi=300
    imagename=name+'_ESP.'+'png'
    cmd.delete('all')
    for filename in cubefiles:
        cmd.load(filename)
    preset.ball_and_stick(selection='all', mode=1)
    cmd.bg_color("white")
    cmd.color("grey50","all")
    util.cnc("all")
    cmd.set('label_size',26) 
    cmd.set('depth_cue',0)
    cmd.set('ray_trace_fog',0) 
    cmd.isosurface('Dt2','Dt', 0.001)
    cmd.ramp_new('espcol', 'ESP', [-.01,-.005,0,.005,.01], ['red','orange', 'yellow','green', 'blue'])
    cmd.set('surface_color', 'espcol', 'Dt2')
    cmd.zoom()
    cmd.disable('espcol')
    cmd.png(imagename, imagesize,imagesize,dpi,1)
    cmd.save(name+'_ESP.'+'pse')



def PlotDimers3D(filenamearray,allindices):
    molsPerImage=len(filenamearray)
    if (molsPerImage % 2) == 0 or (molsPerImage ** 0.5) % 1==0:
        n=molsPerImage
    else:
        n=molsPerImage+1 
    molsperrow=SmallestDivisor(n)
    if n==2:
        molsperrow=1
    imagesize=1180
    dpi=300
    filenamechunks=ChunksList(Chunks(filenamearray,molsPerImage))
    indiceschunks=ChunksList(Chunks(allindices,molsPerImage))
    prevmatslen=len(filenamechunks[0])
    for i in range(len(filenamechunks)):
        filenamesublist=filenamechunks[i]
        indicessublist=indiceschunks[i]
        imagenames=[]
        for j in range(len(filenamesublist)):
            filename=filenamesublist[j]
            indices=indicessublist[j]
            ls=range(len(filenamesublist))
            chunks=ChunksList(Chunks(ls,molsperrow))
            indextorow={}
            for rowidx in range(len(chunks)):
                row=chunks[rowidx]
                for j in row:
                    indextorow[j]=rowidx
            
            fileprefix=filename.split('.')[0]
            imagename=fileprefix+'_3D.'+'png'
            imagenames.append(imagename)
            cmd.delete('all')
            cmd.load(filename)
            preset.ball_and_stick(selection='all', mode=1)
            cmd.bg_color("white")
            cmd.color("grey50","all")
            util.cnc("all")
            cmd.set('label_size',26) 
            cmd.set('depth_cue',0)
            cmd.set('ray_trace_fog',0) 
            firstidx=indices[0]
            secondidx=indices[1]
            firstlab=str(firstidx)
            secondlab=str(secondidx)
            lab=firstlab+'-'+secondlab
            cmd.distance(lab,'index '+firstlab,'index '+secondlab)
            cmd.zoom(lab)
            cmd.png(imagename, imagesize,imagesize,dpi,1)
            cmd.save(fileprefix+'_3D.'+'pse')
        if i>0:
            factor=1
        else:
            factor=0
        firstj=i*prevmatslen+factor
        secondj=firstj+len(filenamesublist)-factor
        prevmatslen=len(filenamesublist)

        basename=fileprefix+'_'+str(firstj)+'-'+str(secondj)
        indextoimage={}
        for index in range(len(filenamesublist)):
            imagename=imagenames[index]
            image=Image.open(imagename)
            indextoimage[index]=image
        cols=len(set(list(indextorow.values())))
        dest = Image.new('RGB', (imagesize*molsperrow,imagesize*cols))
        for j in range(len(filenamesublist)):
            row=indextorow[j]
            x=(j-molsperrow*(row))*imagesize
            y=(row)*imagesize
            dest.paste(indextoimage[j],(x,y))
        dest.show()
        dest.save(basename+'.png')

def ParseXYZ(xyzpath):
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


def WriteOutParamTable(moltotypetoprms,moltotypetoelement):
    tempname='SummaryParams.csv'
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header=['Name','Element','Radius','Depth','Reduction']
        energy_writer.writerow(header)
        for mol,typetoprms in moltotypetoprms.items():
            typetoelement=moltotypetoelement[mol]
            for type,prms in typetoprms.items():
                array=[]
                element=typetoelement[type]
                radius=round(prms[0],3)
                depth=round(prms[1],3)
                shrink=round(prms[2],3)
                array.append(mol)
                array.append(element)
                array.append(radius)
                array.append(depth)
                array.append(shrink)
                energy_writer.writerow(array)  


def PlotLiquidPropsVsTemp(nametotptofinalprops):
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
            ax1.scatter(temparray,refarray, s=10, c='blue', marker="s", label='Target')              
            ax1.set_ylabel(ykey,fontsize=12)
            ax1.set_xlabel(xkey,fontsize=12)
            newtitle=name+' '+propname+ ' Vs Temperature'
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



def WriteOutPropTable(nametotptofinalprops,moltomaxiter):
    tempname='SummaryProps.csv'
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header=['Name','Temperature','Pressure','Density Ref','Density Calc','Density CalcErr','Density Error','Density RelError','Enthalpy Ref','Enthalpy Calc','Enthalpy CalcErr','Enthalpy Error','Enthalpy RelError','Max Iter']
        energy_writer.writerow(header)
        for name,tptofinalprops in nametotptofinalprops.items():
            maxiter=moltomaxiter[name]

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
                    error=all[0]-all[1]
                    error=np.abs(round(error,2))
                    relerror=round((100*error)/truevalue,2)
                    thename=propname+' '+'Error'
                    headeridx=header.index(thename)
                    array[headeridx]=error
                    thename=propname+' '+'RelError'
                    headeridx=header.index(thename)
                    array[headeridx]=relerror
                energy_writer.writerow(array) 


def GrabFinalParameters(jobdirs):
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


def GrabParameterValues(prmfiles):
    moltotypetoprms={}
    moltotypetoelement={}
    for prmfile in prmfiles:
        temp=open(prmfile,'r')
        results=temp.readlines()
        temp.close()
        foundatomblock=False
        typetoprms={}
        typetoelement={}
        moltotypes={}
        for line in results:
            linesplit=line.split()
            if len(linesplit)>1:
                first=linesplit[0]
                if first=='atom':
                    foundatomblock=True
                    typenum=linesplit[1]
                    element=linesplit[3]
                    typetoelement[typenum]=element
                    mol=linesplit[4].replace('"','')
                    if mol not in moltotypes.keys():
                        moltotypes[mol]=[]
                    moltotypes[mol].append(typenum)      
                if foundatomblock==False:
                    typenum=linesplit[1]
                    radius=float(linesplit[2])
                    depth=float(linesplit[3])
                    try:
                        shift=float(linesplit[4])
                    except:
                        shift=1
                    typetoprms[typenum]=[radius,depth,shift]
        for mol,types in moltotypes.items():
            if mol not in moltotypetoprms.keys():
                moltotypetoprms[mol]={}
                moltotypetoelement[mol]={}
            for type in types: 
                if type in typetoprms.keys():
                    prms=typetoprms[type]
                    moltotypetoprms[mol][type]=prms
                    moltotypetoelement[mol][type]=typetoelement[type] 

    return moltotypetoprms,moltotypetoelement


def WriteOutQMTable(nametoformulatormse,nametoformulatomse):
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

                     
def PlotLiquidPropsVsTempForGroups(nametotemparray,nametoproptokeytoarray,groupednames):
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
                        grptoproptoaxes[grp][propname].scatter(temparray,calcarray, s=10, c=c, marker="o", label='AMOEBA'+string)
                        grptoproptoaxes[grp][propname].errorbar(temparray, calcarray, yerr=calcerrarray,c=c, fmt="o")
                        grptoproptoaxes[grp][propname].scatter(temparray,refarray, s=10, c=c, marker="s", label='Target'+string)              
                        grptoproptoaxes[grp][propname].set_ylabel(ykey,fontsize=12)
                        grptoproptoaxes[grp][propname].set_xlabel(xkey,fontsize=12)
                        newtitle=totalname+' '+propname+ ' Vs Temperature'
                        grptoproptoaxes[grp][propname].set_title(newtitle)
                        grptoproptoaxes[grp][propname].legend()
                        imagename=newtitle+'.png'
                        grptoproptofig[grp][propname].savefig(imagename)





curdir=os.getcwd()
jobdirs=GrabJobDirectories(fbdir)
outputfiles=GrabOutputFiles(jobdirs)
outputfilestojobdir=dict(zip(outputfiles,jobdirs))
tpdiclist,qmdiclist,newoutputfiles=GrabResultsAllMolecules(outputfiles)
jobdirs=[outputfilestojobdir[i] for i in newoutputfiles]
qmtargetnamedic,nametodimerstructs,truenametoindices=GrabDimerDistanceInfo(qmdiclist,jobdirs)
neatliqnametoindices=ExtractNeatLiquidIndices(truenametoindices)
os.chdir(curdir)
nametofilenametoformula=PlotAllDimers(nametodimerstructs,truenametoindices)
nametotptofinalprops,nametoformulatormse,nametoformulatomse=PlotAllFBJobs(tpdiclist,qmtargetnamedic,nametofilenametoformula,truenametoindices)
moltotptoarc,moltomaxiter=GrabFinalNeatLiquidTrajectories(jobdirs)
WriteOutPropTable(nametotptofinalprops,moltomaxiter)
PlotAllRDFs(moltotptoarc,neatliqnametoindices)
prmfiles=GrabFinalParameters(jobdirs)
moltotypetoprms,moltotypetoelement=GrabParameterValues(prmfiles)
WriteOutParamTable(moltotypetoprms,moltotypetoelement)
WriteOutQMTable(nametoformulatormse,nametoformulatomse)
poltypedirs,groupedpoltypedirs=GrabPoltypeDirectories(jobdirs)
nametocubefiles,groupednames=GrabCubeFiles(groupedpoltypedirs)
PlotAllESPSurfaces(nametocubefiles)
nametotemparray,nametoproptokeytoarray=PlotLiquidPropsVsTemp(nametotptofinalprops)
PlotLiquidPropsVsTempForGroups(nametotemparray,nametoproptokeytoarray,groupednames)

