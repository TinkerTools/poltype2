import os
import numpy as np
import time
import math


def CheckTrajectory(poltype,traj):
    foundatomnum=False
    badtraj=False
    poltype.WriteToLog('Checking for issues in trajectory '+traj,prin=True)
    atomindextoentrylength={} 
    with open(traj) as infile:
        for line in infile:
            linesplit=line.split()
            if foundatomnum==False and len(linesplit)==1:
                atomnum=int(linesplit[0])
                foundatomnum=True
                lastatomindex=atomnum
            boxline=False
            if '90.000' in line and len(linesplit)==6: 
                boxline=True
            if boxline==False and len(linesplit)!=1:
                if linesplit[0].isdigit():
                    atomindex=int(linesplit[0])
                    lastatomindex=atomindex
                    if len(atomindextoentrylength.keys())!=atomnum:
                        atomindextoentrylength[atomindex]=len(linesplit)
                    else:
                        if atomindex in atomindextoentrylength.keys():
                            trueentrylength=atomindextoentrylength[atomindex]
                            if len(linesplit)!=trueentrylength:
                                badtraj=True
                                break
                        else:
                            badtraj=True 
                            break
                else:
                    badtraj=True
                    break
            elif boxline==True:
                if lastatomindex!=atomnum:
                    badtraj=True
                    break
            lastline=line
 
    if badtraj==False:
        pass
    else:
        poltype.WriteToLog('Deleting trajectory '+traj,prin=True)
    return badtraj

def CheckBARFile(poltype,filepath):
    temp=open(filepath,'r')
    results=temp.readlines()
    temp.close()
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)==2:
            count+=1
    if count==2:
        deletefile=False
    else:
        deletefile=True
    return deletefile

def CheckFileTermination(poltype,f,steps=None,equil=False,firsttime=True):
   term=False
   deletefile=False
   if poltype.checktraj==False:
       firsttime=False
   if os.path.exists(f):
       temp=open(f,'r')
       results=temp.readlines()
       temp.close()
       head,tail=os.path.split(f)
       files=os.listdir(head)
       foundfilename=False
       if equil==False:
           for tf in files:
               if '.xyz' in tf and 'pymol' not in tf:
                   foundfilename=True
                   xyzfilepath=os.path.join(head,tf)
                   dynpath=xyzfilepath.replace('.xyz','.dyn')
                   filepath=xyzfilepath.replace('.xyz','.arc')

       if firsttime==True:
           if foundfilename==True:
               if os.path.isfile(filepath):
                   try:
                        deletefile=CheckTrajectory(poltype,filepath)
                   except:
                        deletefile=True
                   if deletefile==True:
                       os.remove(filepath) # remove trajectory and also outputfile for trajectory
                       os.remove(dynpath)
       for tf in files:
           if '.' in tf:
               tfsplit=tf.split('.')
               ext=tfsplit[-1]
               #if ext=='bar':
               #    filepath=os.path.join(head,tf)
               #    deletefile=CheckBARFile(poltype,filepath)
               #    if deletefile==True#:
               #        os.remove(filepath)
       error=False
       for line in results:
           if 'INDUCE  --  Warning, Induced Dipoles are not Converged' in line:
               error=True
               errorline=line
           if 'Normal Termination' in line or 'Potential Energy Values Written To :' in line or 'Free Energy Difference via BAR' in line:
               term=True
               if 'Potential Energy Values Written To :' in line:
                   foundbar=False
                   for otherf in files:
                       if '.' in tf:
                           tfsplit=otherf.split('.')
                           ext=tfsplit[-1]
                           if ext=='bar':
                               foundbar=True
                   if foundbar==False:
                       term=False

               '''
               if 'Potential Energy Values Written To :' in line:
                   linesplit=line.split()
                   lastterm=linesplit[-1]
                   head,tail=os.path.split(lastterm)
                   tailsplit=tail.split('.')
                   ext=tailsplit[-1]
                   if ext!='bar': # if its bar_2 then not terminated need to redo, program only reads .bar
                       term=False
               '''
               break
           if 'Dynamics Steps' in line:
               if equil==True:
                   stepsize=float(poltype.equiltimestep)*10**-6
                   writefreq=float(poltype.equilwritefreq)*.001             
               else:
                   stepsize=float(poltype.proddyntimestep)*10**-6
                   writefreq=float(poltype.proddynwritefreq)*.001
               totaltime=steps*stepsize 
               totalframes=math.floor(totaltime/writefreq)
               stepsperframe=int(writefreq/stepsize)
               tolratio=(stepsperframe*totalframes)/steps
               
               linesplit=line.split()
               try:
                   stepnum=float(linesplit[-3])
                   ratio=np.abs(stepnum)/steps
                   if np.abs(ratio-tolratio)<.01:       
                       term=True

                   if foundfilename==True and term==True:
                       if not os.path.isfile(filepath): 

                           term=False
               except:
                   deletefile=True 
                   if firsttime==True:
                       if foundfilename==True:
                           os.remove(filepath)
                           os.remove(dynpath)
   if error==True:
       poltype.WriteToLog('Error Detected '+f+' '+errorline)

   return term,deletefile,error

def CheckFilesTermination(poltype,outputfilepathlist,steps=None,equil=False,firsttime=False):
    finished=True
    islist=False
    if type(steps)==list:
        islist=True
    outputfilestodelete=[]
    jobsfinished=[]
    anyerrors=False
    for i in range(len(outputfilepathlist)):
        outputfilepath=outputfilepathlist[i]
        if steps!=None:
            if islist==False:
                truesteps=float(steps)
            else:
                truesteps=float(steps[i])
        else:
            truesteps=steps   
        terminate,deletefile,error=CheckFileTermination(poltype,outputfilepath,truesteps,equil,firsttime)
        if error==True:
            anyerrors=True
        if deletefile==True:
            outputfilestodelete.append(outputfilepath)
            terminate=False
        if terminate==False:
            finished=False
            poltype.nextfiletofinish=outputfilepath
        else:
            jobsfinished.append(outputfilepath)
    for filename in outputfilestodelete:
        poltype.WriteToLog('File being deleted '+filename,prin=True)
        os.remove(filename)
    if len(outputfilepathlist)==0:
        percentfinished=100
    else:
        percentfinished=round((len(jobsfinished)*100)/len(outputfilepathlist),2)
    if anyerrors==True:
        raise ValueError('Errors Detected in Output Files!')
    return finished,percentfinished


def CheckFileRecentlyUpdated(poltype,outputfilepath):
    updated=False
    if os.path.isfile(outputfilepath):
        Ftime=os.path.getmtime(outputfilepath)
        reltime=time.time()-Ftime
        htime=reltime*0.000277778
        updatetime=.25 # hours.
        if htime<updatetime:
            updated=True

    return updated


def CheckFilesRecentlyUpdated(poltype,outputfilepathlist):
    recentlyupdated=False
    for i in range(len(outputfilepathlist)):
        outputfilepath=outputfilepathlist[i]
        updated=CheckFileRecentlyUpdated(poltype,outputfilepath)
        if updated==True:
            recentlyupdated=True
            break

    return recentlyupdated
