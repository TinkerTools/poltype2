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


def GrabResultsAllMolecules(outputfiles):
    tpdiclist=[]
    qmdiclist=[]
    newoutputfiles=[]
    for outputfile in outputfiles:
        #try:
        tpdic,qmdic,outputfile=GrabResults(outputfile)
        tpdiclist.append(tpdic)
        qmdiclist.append(qmdic)
        newoutputfiles.append(outputfile)
        #except:
            #pass
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
                ref=float(linesplit[3])
                calc=float(linesplit[4])
                tpdic[name][tp][property][refkeyword].append(ref)
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


def PlotAllFBJobs(tpdiclist,qmtargetnamedic):
    for tpdic in tpdiclist:
        PlotFBLiq(tpdic)
    for qmtarget,targetdic in qmtargetnamedic.items():
        namesplit=qmtarget.split('-')
        name=namesplit[0]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        PlotFBQM(qmtarget,targetdic)
        os.chdir('..')

def PlotFBQM(qmtarget,targetdic):
    distances=targetdic['Distances']
    refvalues=np.transpose(np.array(targetdic['Ref']))
    calcvalues=np.transpose(np.array(targetdic['Calc']))
    rmsvalues=[]
    for i in range(len(refvalues)):
        refs=refvalues[i]
        calcs=calcvalues[i]
        rms = sqrt(mean_squared_error(calcs, refs))
        rmsvalues.append(rms)
        label1='AMOEBA'
        label2='Target'
        ykey='Interaction Energy'
        xkey='Relative Distance'
        labels=[label1,label2]
        title=qmtarget+' Iter='+str(i)
        if i==len(refvalues)-1:
            ScatterPlot2D(title,distances,calcs,refs,xkey,ykey,labels)
    x=list(range(len(rmsvalues)))
    label1='RMSE Interaction Energy'
    ykey='RMSE Interaction Energy'
    xkey='Iteration'
    labels=[label1]
    title=qmtarget+' RMSE vs Iteration'
    ScatterPlot1D(title,x,rmsvalues,xkey,ykey,labels)

    
         


def PlotFBLiq(tpdic):
    for name in tpdic.keys():
        dic=tpdic[name]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        n=len(list(dic.keys()))
        color = iter(cm.rainbow(np.linspace(0, 1, n)))
        proptofigs={}
        proptoaxes={}
        for tp in dic.keys():
            innerdic=dic[tp]
            c = next(color)
            for property in innerdic.keys():
                relprop=property+'RelativeError'
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
                refarray=np.array(innermostdic[refkey])
                calcarray=np.array(innermostdic[calckey])
                if property=='Enthalpy':
                    scale=0.239006 # convert kj to kcal
                else:
                    scale=1
                refarray=scale*refarray
                calcarray=scale*calcarray
                err=np.abs(refarray-calcarray)
                relerr=100*(err/refarray)
                label1=property+' AMOEBA'
                label2=property+' Target'
                ykey=property
                xkey='Iterations'
                labels=[label1,label2]
                string='T=%s,P=%s'%(str(tp[0]),str(tp[1]))
                title=name+' '+property+' '+string
                x=list(range(len(refarray)))
                ScatterPlot2D(title,x,calcarray,refarray,xkey,ykey,labels)
                proptoaxes[property].scatter(x,calcarray, s=10, c=c, marker="s", label='AMOEBA '+string)
                proptoaxes[property].scatter(x,refarray, s=10, c=c, marker="o", label='Target '+string)              
                proptoaxes[property].set_ylabel(ykey,fontsize=12)
                proptoaxes[property].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+property
                proptoaxes[property].set_title(newtitle)
                proptoaxes[property].legend()
                x=np.array(x)
                refarray=np.array(refarray)
                calcarray=np.array(calcarray)
                x_new = np.linspace(x.min(),x.max(),500)
                f = interp1d(x,calcarray, kind='quadratic')
                y_smooth=f(x_new)
                proptoaxes[property].plot(x_new,y_smooth,color=c)
                f = interp1d(x,refarray, kind='quadratic')
                y_smooth=f(x_new)
                proptoaxes[property].plot(x_new,y_smooth,color=c)
                imagename=newtitle+'.png'
                proptofigs[property].savefig(imagename)
                proptoaxes[relprop].scatter(x,relerr, s=10, c=c, marker="s", label=string)
                f = interp1d(x,relerr, kind='quadratic')
                y_smooth=f(x_new)
                proptoaxes[relprop].plot(x_new,y_smooth,color=c)
                ykey='Relative Error '+property
                proptoaxes[relprop].set_ylabel(ykey,fontsize=12)
                proptoaxes[relprop].set_xlabel(xkey,fontsize=12)
                newtitle=name+' '+'Relative Error '+property
                proptoaxes[relprop].set_title(newtitle)
                proptoaxes[relprop].legend()
                imagename=newtitle+'.png'
                proptofigs[relprop].savefig(imagename)

 
        os.chdir('..')

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


def ScatterPlot1D(title,x,y1,xkey,ykey,labels):
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
    try:
        f = interp1d(x,y1, kind='quadratic')
        y_smooth=f(x_new)
        plt.plot(x_new,y_smooth,color='red')
    except:
        pass
    plt.savefig(imagename)



def GrabDimerDistanceInfo(qmdiclist,jobdirs):
    qmtargetnamedic={}
    nametodimerstructs={}
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
                split=name.split('_')
                truename='_'.join(split[1:-1]) 
                nametodimerstructs[truename]=dimernametotinkxyzlines
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
                    valuedic=[namedic[k] for k in sortedindices]
                    refvalues=[k['Ref'] for k in valuedic]
                    calcvalues=[k['Calc'] for k in valuedic]
                    if truename not in qmtargetnamedic.keys():
                        qmtargetnamedic[truename]={}
                    if 'Distances' not in qmtargetnamedic[truename].keys():
                        qmtargetnamedic[truename]['Distances']=sorteddistances
                    if 'Ref' not in qmtargetnamedic[truename].keys():
                        qmtargetnamedic[truename]['Ref']=refvalues
                    if 'Calc' not in qmtargetnamedic[truename].keys():
                        qmtargetnamedic[truename]['Calc']=calcvalues
            
                os.chdir('..')

    return qmtargetnamedic,nametodimerstructs 


def Chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ChunksList(gen):
    newlst=[]
    for item in gen:
        newlst.append(item)
    return newlst

def PlotAllDimers(nametodimerstructs):
    for name in nametodimerstructs.keys():
        dic=nametodimerstructs[name]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        filenames=[]
        for dimername,tinkerxyzlines in dic.items():
            WriteFile(dimername,tinkerxyzlines)
            newfilename=dimername.replace('.xyz','cartesian.xyz')
            ConvertTinkerXYZToCartesian(dimername,newfilename)
            filenames.append(newfilename) 
        PlotDimers3D(filenames)
        os.chdir('..')


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


def PlotDimers3D(filenamearray):
    from pymol import cmd,preset,util
    from PIL import Image
    from pymol.vfont import plain
    from pymol.cgo import CYLINDER,cyl_text
    molsperrow=2
    molsPerImage=molsperrow**2
    imagesize=400
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
            cmd.zoom()
            cmd.ray(imagesize,imagesize)
            cmd.png(imagename, imagesize,imagesize)
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
        dest = Image.new('RGB', (imagesize*molsperrow,imagesize*molsperrow))
        for j in range(len(filenamesublist)):
            row=indextorow[j]
            x=(j-molsperrow*(row))*imagesize
            y=(row)*imagesize
            dest.paste(indextoimage[j],(x,y))
        dest.show()
        dest.save(basename+'.png')

curdir=os.getcwd()
jobdirs=GrabJobDirectories(fbdir)
outputfiles=GrabOutputFiles(jobdirs)
outputfilestojobdir=dict(zip(outputfiles,jobdirs))
tpdiclist,qmdiclist,newoutputfiles=GrabResultsAllMolecules(outputfiles)
jobdirs=[outputfilestojobdir[i] for i in newoutputfiles]
qmtargetnamedic,nametodimerstructs=GrabDimerDistanceInfo(qmdiclist,jobdirs)
os.chdir(curdir)
PlotAllFBJobs(tpdiclist,qmtargetnamedic)
#PlotAllDimers(nametodimerstructs)
