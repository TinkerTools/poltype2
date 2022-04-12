import productiondynamics as prod 
import submitjobs
import tables
import terminate as term
import time
import numpy as np
import shutil
import keyfilemodifications as keymods
import os
import restraints
import sys


def AverageBoxSizeFromNPTArc(poltype,arcpath,firstframe,lastframe,framestep,index):
    firstline=True
    framecount=0
    framestoextract=np.arange(firstframe,lastframe+1,framestep)
    framearray=[]
    extractingframe=False
    aaxisarray=[]
    baxisarray=[]
    caxisarray=[]
    with open(arcpath) as infile:
        for line in infile:
            if firstline==True:
                firstlinesplit=line.split()
                framestring=firstlinesplit[0].lstrip().rstrip()
                firstline=False
            linesplit=line.split()
            if '90.000000' in line and framecount in framestoextract:
                aaxis=float(linesplit[0])
                baxis=float(linesplit[1])
                caxis=float(linesplit[2])
                aaxisarray.append(aaxis)
                baxisarray.append(baxis)
                caxisarray.append(caxis)
            if framestring in line and (len(linesplit)==1):
                extractingframe=False
                framecount+=1
    aaxisaverage=np.mean(np.array(aaxisarray))
    baxisaverage=np.mean(np.array(baxisarray))
    caxisaverage=np.mean(np.array(caxisarray))
    if index!=None: 
        poltype.tabledict[index]['Average Box Size']=aaxisaverage
        tables.WriteTableUpdateToLog(poltype)
    return aaxisaverage,baxisaverage,caxisaverage



def EquilbriateDynamicCommand(poltype,steps,ensemble,temp,outputfilename,equilboxfilename,configkeyfilename,NPT=False):
    head,tail=os.path.split(outputfilename)
    outputfilename=tail
    if NPT==False:
        cmdstr=poltype.truedynamicpath+' '+equilboxfilename+' '+ '-k'+' '+configkeyfilename+' '+str(steps)+' '+ str(poltype.equiltimestep)+' '+ str(poltype.equilwritefreq)+' '+str(ensemble)+' '+str(temp)+' '+ ' > '+outputfilename  
    else:
        cmdstr=poltype.truedynamicpath+' '+equilboxfilename+' '+ '-k'+' '+configkeyfilename+' '+str(steps)+' '+ str(poltype.equiltimestep)+' '+ str(poltype.equilwritefreq)+' '+str(ensemble)+' '+str(temp)+' '+str(poltype.pressure)+' '+ ' > '+outputfilename  

    return cmdstr         

def ExecuteEquilibriation(poltype):
    jobtolog=[]
    jobtoidx=[]
    jobtojobpath=[]
    jobtorestconstant=[]
    jobtoinputfilepaths=[]
    jobtooutputfiles=[]
    jobtoabsolutebinpath=[]
    for i in range(len(poltype.equiloutputarray)):
        jobtolog.append({})
        jobtoidx.append({})
        jobtojobpath.append({})
        jobtorestconstant.append({})
        jobtoinputfilepaths.append({})
        jobtooutputfiles.append({})
        jobtoabsolutebinpath.append({})

    for idx in range(len(poltype.equilibriatescheme)):
        temp=poltype.equilibriatescheme[idx]
        restrainpositionconstant=poltype.equilibriaterestscheme[idx]
        for i in range(len(poltype.equiloutputarray)):
            equiloutputarray=poltype.equiloutputarray[i]
            equilboxfilename=poltype.equilboxfilename[i]
            configkeyfilename=poltype.configkeyfilename[i][0]
            if idx==len(poltype.equilibriatescheme)-1:
                steps=poltype.lastequilstepsNVT   
            else:
                steps=poltype.equilstepsNVT
            cmdstr=EquilbriateDynamicCommand(poltype,steps,poltype.NVTensem,temp,equiloutputarray[idx],equilboxfilename,configkeyfilename)
            jobtolog[i][cmdstr]=equiloutputarray[idx]
            jobtoidx[i][cmdstr]=idx
            jobtojobpath[i][cmdstr]=poltype.outputpath
            jobtorestconstant[i][cmdstr]=restrainpositionconstant
            head,tail=os.path.split(equiloutputarray[idx]) 
            arcname=equilboxfilename.replace('.xyz','.arc')
            dynname=equilboxfilename.replace('.xyz','.dyn')

            inputfilepaths=[poltype.outputpath+equilboxfilename,poltype.outputpath+configkeyfilename,poltype.outputpath+poltype.prmfilepath]
            if os.path.isfile(poltype.outputpath+arcname):
                inputfilepaths.append(poltype.outputpath+arcname)
            if os.path.isfile(poltype.outputpath+dynname):
                inputfilepaths.append(poltype.outputpath+dynname)

            outputfilenames=[tail,arcname,dynname]
            absbinpath=poltype.which(poltype.truedynamicpath)
            jobtoinputfilepaths[i][cmdstr]=inputfilepaths
            jobtooutputfiles[i][cmdstr]=outputfilenames
            jobtoabsolutebinpath[i][cmdstr]=absbinpath


    temp=poltype.equilibriatescheme[-1]
    for i in range(len(poltype.equiloutputarray)):
        equiloutputarray=poltype.equiloutputarray[i]
        equilboxfilename=poltype.equilboxfilename[i]
        configkeyfilename=poltype.configkeyfilename[i][0]
        if i==0 and poltype.solvation==True and poltype.addsolvionwindows==True:
            idx=-2
        else:
            idx=-1

        cmdstr=EquilbriateDynamicCommand(poltype,poltype.equilstepsNPT,poltype.NPTensem,temp,equiloutputarray[idx],equilboxfilename,configkeyfilename,NPT=True)
        jobtolog[i][cmdstr]=equiloutputarray[idx]
        jobtoidx[i][cmdstr]=idx
        jobtojobpath[i][cmdstr]=poltype.outputpath
        jobtorestconstant[i][cmdstr]=0
        head,tail=os.path.split(equiloutputarray[idx]) 
        outputfilenames=[tail]
        absbinpath=poltype.which(poltype.truedynamicpath)
        arcname=equilboxfilename.replace('.xyz','.arc')
        dynname=equilboxfilename.replace('.xyz','.dyn')
        inputfilepaths=[poltype.outputpath+equilboxfilename,poltype.outputpath+configkeyfilename,poltype.outputpath+poltype.prmfilepath]
        if os.path.isfile(poltype.outputpath+arcname):
            inputfilepaths.append(poltype.outputpath+arcname)
        if os.path.isfile(poltype.outputpath+dynname):
            inputfilepaths.append(poltype.outputpath+dynname)

        jobtoinputfilepaths[i][cmdstr]=inputfilepaths
        jobtooutputfiles[i][cmdstr]=outputfilenames
        jobtoabsolutebinpath[i][cmdstr]=absbinpath

        if i==0 and poltype.solvation==True and poltype.addsolvionwindows==True:
            shutil.copy(poltype.ionboxfilename,poltype.ionequilboxfilename)
            cmdstr=EquilbriateDynamicCommand(poltype,poltype.equilstepsionNPT,poltype.NPTensem,temp,equiloutputarray[-1],poltype.ionequilboxfilename,poltype.ionkeyfilename,NPT=True)
            jobtolog[i][cmdstr]=equiloutputarray[-1]
            jobtoidx[i][cmdstr]=-1
            jobtojobpath[i][cmdstr]=poltype.outputpath
            jobtorestconstant[i][cmdstr]=0
            head,tail=os.path.split(equiloutputarray[-1]) 
            outputfilenames=[tail]
            absbinpath=poltype.which(poltype.truedynamicpath)
            arcname=poltype.ionequilboxfilename.replace('.xyz','.arc')
            dynname=poltype.ionequilboxfilename.replace('.xyz','.dyn')
            inputfilepaths=[poltype.outputpath+poltype.ionequilboxfilename,poltype.outputpath+poltype.ionkeyfilename,poltype.outputpath+poltype.prmfilepath]
            if os.path.isfile(poltype.outputpath+arcname):
                inputfilepaths.append(poltype.outputpath+arcname)
            if os.path.isfile(poltype.outputpath+dynname):
                inputfilepaths.append(poltype.outputpath+dynname)

            jobtoinputfilepaths[i][cmdstr]=inputfilepaths
            jobtooutputfiles[i][cmdstr]=outputfilenames
            jobtoabsolutebinpath[i][cmdstr]=absbinpath

            

    jobtologdic=jobtolog[0]
    for k in range(len(jobtologdic)):
        newjobtolog={}
        newjobtojobpath={}
        newjobtoinputfilepaths={}
        newjobtooutputfiles={}
        newjobtoabsolutebinpath={}
        for i in range(len(poltype.equiloutputarray)):
            jobtologdic=jobtolog[i]
            jobtoidxdic=jobtoidx[i]
            jobtojobpathdic=jobtojobpath[i]
            jobtoinputfilepathsdic=jobtoinputfilepaths[i]
            jobtooutputfilesdic=jobtooutputfiles[i]
            jobtoabsolutebinpathdic=jobtoabsolutebinpath[i]
            jobtorestconstantdic=jobtorestconstant[i]
            equiloutputarray=poltype.equiloutputarray[i]
            stepsarray=poltype.equiloutputstepsarray[i]
            job=list(jobtologdic.keys())[k]
            idx=jobtoidxdic[job]
            restrainpositionconstant=jobtorestconstantdic[job]
            outarray=[equiloutputarray[idx]]
            finished=term.CheckFilesTermination(poltype,outarray,stepsarray,True)[0]
            recentlyupdated=term.CheckFilesRecentlyUpdated(poltype,outarray)
            configkeyfilename=poltype.configkeyfilename[i][0]
            ligandindices=poltype.ligandindices[i]
            if finished==False and recentlyupdated==False:
                files=os.listdir()
                for f in files:
                    if '.end' in f:
                        os.remove(f)
                newjobtolog[job]=jobtologdic[job] # only submit one equilibriate job at a time (increasing temp...)
                newjobtojobpath[job]=jobtojobpathdic[job] 
                newjobtoinputfilepaths[job]=jobtoinputfilepathsdic[job]
                newjobtooutputfiles[job]=jobtooutputfilesdic[job]
                newjobtoabsolutebinpath[job]=jobtoabsolutebinpathdic[job]
                if poltype.complexation==True and i==0:
                    keymods.RemoveKeyWords(poltype,configkeyfilename,['restrain-group','restrain-position'])
                    restraints.AddHarmonicRestrainGroupTermsToKeyFile(poltype,configkeyfilename,poltype.restraintdistance,restrainpositionconstant)
                    totalatomnumberxyzfilename=poltype.totalatomnumberxyzfilename[i]
                    if restrainpositionconstant!=0 and poltype.restrainreceptorligand==True: 
                        resposstring='restrain-position -'+str(1)+' '+str(totalatomnumberxyzfilename-len(ligandindices))+' '+str(restrainpositionconstant)+' '+str(poltype.equilrestrainsphereradius)+'\n'
                        keymods.AddKeyWord(poltype,poltype.outputpath+configkeyfilename,resposstring)
                    resposstring='restrain-position -'+str(poltype.norotpair[0])+' '+str(poltype.norotpair[0])+' '+str(poltype.restrainpositionconstant)+' '+str(poltype.norotrestrainsphereradius)+'\n'
                    keymods.AddKeyWord(poltype,poltype.outputpath+configkeyfilename,resposstring)
                    resposstring='restrain-position -'+str(poltype.norotpair[1])+' '+str(poltype.norotpair[1])+' '+str(poltype.restrainpositionconstant)+' '+str(poltype.norotrestrainsphereradius)+'\n'
                    keymods.AddKeyWord(poltype,poltype.outputpath+configkeyfilename,resposstring)
        submitjobs.SubmitJobs(poltype,newjobtolog,newjobtojobpath,newjobtoinputfilepaths,newjobtooutputfiles,newjobtoabsolutebinpath,poltype.outputpath+poltype.equiljobsfilename)
        messages=[]
        for i in range(len(poltype.equiloutputarray)):
            equiloutputarray=poltype.equiloutputarray[i]
            stepsarray=poltype.equiloutputstepsarray[i]
            jobtoidxdic=jobtoidx[i]
            jobtologdic=jobtolog[i]
            job=list(jobtologdic.keys())[k]
            idx=jobtoidxdic[job]
            outarray=[equiloutputarray[idx]]
            finished=term.CheckFilesTermination(poltype,outarray,stepsarray,True)[0]

            while finished==False:
                msg='System equilibriation is not complete '
                if msg not in messages:
                    poltype.WriteToLog(msg,prin=True)
                    messages.append(msg)
                time.sleep(poltype.waitingtime)
                finished=term.CheckFilesTermination(poltype,outarray,stepsarray,True)[0]

        poltype.WriteToLog('Single system equilibriation job is complete ',prin=True)

    poltype.WriteToLog('System equilibriation is complete ',prin=True)


def ExtractLastTinkerFrame(poltype,arcfilename,newname):
    poltype.WriteToLog('Extracting the last frame of '+arcfilename,prin=True)
    framenum=int(os.path.getsize(poltype.outputpath+arcfilename)/os.path.getsize(poltype.outputpath+arcfilename.replace('.arc','.xyz')))
    ExtractTinkerFrames(poltype,poltype.outputpath+arcfilename,framenum,framenum,1,framenum,newname)

def ExtractTinkerFrames(poltype,arcpath,firstframe,lastframe,framestep,totalnumberframes,newname):
    firstline=True
    framecount=0
    framestoextract=np.arange(firstframe,lastframe+1,framestep)
    framearray=[]
    extractingframe=False
    with open(arcpath) as infile:
        for line in infile:
            if firstline==True:
                firstlinesplit=line.split()
                framestring=firstlinesplit[0].lstrip().rstrip()
                firstline=False
           
            linesplit=line.split()
            if framestring in line and (len(linesplit)==1):
                extractingframe=False
                framecount+=1
                if len(framearray)!=0:
                    numberofzeroes=len(str(totalnumberframes))-len(str(framecount))
                    zerostring=''
                    for i in range(numberofzeroes):
                        zerostring+='0'
                    framename=arcpath.replace('.arc','.'+zerostring+str(framecount))
                    temp=open(framename,'w')
                    for saveline in framearray:
                        temp.write(saveline)
                    temp.close()
                    framearray=[]
                if framecount in framestoextract:
                    extractingframe=True
            if(extractingframe):
                framearray.append(line)
            if len(framearray)!=0: # for last frame
                numberofzeroes=len(str(totalnumberframes))-len(str(framecount))
                zerostring=''
                for i in range(numberofzeroes):
                    zerostring+='0'
                framename=arcpath.replace('.arc','.'+zerostring+str(framecount))
                temp=open(framename,'w')
                for saveline in framearray:
                    temp.write(saveline)
                temp.close()
    if os.path.exists(framename):
        os.rename(framename,newname)
    return
  
def EquilibriationProtocol(poltype):
    if poltype.equilfinished==False:
        if poltype.restrainatomgroup1==None and poltype.restrainatomgroup2==None and poltype.complexation==True and poltype.restrainreceptorligand==True:
            restraints.ComputeIdealGroupRestraints(poltype,poltype.minboxfilename[0])
        if poltype.complexation==True and poltype.restrainreceptorligand==True: 
            dist=restraints.AverageCOMGroups(poltype,poltype.minboxfilename[0])
            poltype.restraintdistance=dist

        for i in range(len(poltype.minboxfilename)):
            firstxyz=poltype.minboxfilename[i]
            secondxyz=poltype.equilboxfilename[i]
            shutil.copy(firstxyz,secondxyz)

        ExecuteEquilibriation(poltype)
        firstframe=poltype.equilframenum-poltype.equilframenumNPT
        lastframe=poltype.equilframenum
        for i in range(len(poltype.equilarcboxfilename)):
            equilarcboxfilename=poltype.equilarcboxfilename[i]
            configkeyfilenamelist=poltype.configkeyfilename[i]
            
            if poltype.productiondynamicsNVT==True:

                aaxis,baxis,caxis=AverageBoxSizeFromNPTArc(poltype,poltype.outputpath+equilarcboxfilename,firstframe,lastframe,1,i)
                for configkeyfilename in configkeyfilenamelist:
                    
                    keymods.RemoveKeyWords(poltype,configkeyfilename,['axis'])
                    keymods.AddKeyWord(poltype,configkeyfilename,'a-axis'+' '+str(aaxis)+'\n')
            proddynboxfilename=poltype.proddynboxfilename[i]
            proddynboxfilenamepymol=poltype.proddynboxfilenamepymol[i]
            if not os.path.isfile(poltype.outputpath+proddynboxfilename):
                ExtractLastTinkerFrame(poltype,equilarcboxfilename,poltype.proddynboxfilename[i])
            if i==0 and poltype.solvation==True and poltype.addsolvionwindows==True:
                if not os.path.isfile(poltype.outputpath+poltype.ionproddynboxfilename):
                    ExtractLastTinkerFrame(poltype,poltype.ionequilarcfilename,poltype.ionproddynboxfilename)

            

            for configkeyfilename in configkeyfilenamelist:
                keymods.RemoveKeyWords(poltype,configkeyfilename,['restrain','group'])

            if poltype.complexation==True and poltype.proddyngrprests==True and i==0:
                equildist=restraints.AverageCOMGroups(poltype,equilarcboxfilename)
                poltype.restraintdistance=equildist
                if poltype.restrainreceptorligand==True:
                    for configkeyfilename in configkeyfilenamelist:
                        restraints.AddHarmonicRestrainGroupTermsToKeyFile(poltype,poltype.outputpath+configkeyfilename,equildist,poltype.distancerestraintconstant)
                    restraints.GroupRestraintFreeEnergyFix(poltype)

            if not os.path.isfile(poltype.outputpath+proddynboxfilenamepymol):
                poltype.PymolReadableFile(proddynboxfilename,proddynboxfilenamepymol)

        
