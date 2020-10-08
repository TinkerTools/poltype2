import os
import sys
import subprocess
import time
import getopt
import re
from itertools import islice

global cpunodelistfilepath
global gpucardlistfilepath
global bashrcpath
global restrictedprogramtonumber
global currentrestrictedprogramtoprocesslist
global pidfile
global loggerfile
global jobtoinfo
global masterloghandle
global writepath
global cpuprogramexceptionlist
global gpuprogramexceptionlist
global cpuprogramlist
global gpuprogramlist

cpunodelistfilepath='/home/bdw2292/ExternalAPIRenLab/cpunodes.txt'
gpucardlistfilepath='/home/bdw2292/ExternalAPIRenLab/gpucards.txt'
bashrcpath=None
jobinfofilepath=None
pidfile='/home/bdw2292/ExternalAPIRenLab/daemon.pid'
loggerfile='/home/bdw2292/ExternalAPIRenLab/logger.txt'
jobtoinfo='/home/bdw2292/ExternalAPIRenLab/jobtoinfo.txt'
writepath='/home/bdw2292/ExternalAPIRenLab/'

masterloghandle=open(loggerfile,'a',buffering=1)
sleeptime=60
cpuprogramexceptionlist=['psi4','g09','g16',"cp2k.ssmp","mpirun_qchem","dynamic.x"]
gpuprogramexceptionlist=['dynamic_omm.x']
restrictedprogramtonumber={'bar.x':10,'bar_omm.x':10}
currentrestrictedprogramtoprocesslist={'bar.x':[],'bar_omm.x':[]}
cpuprogramlist=['psi4','g09','g16',"cp2k.ssmp","mpirun_qchem","dynamic.x"]
gpuprogramlist=['dynamic_omm.x','bar_omm.x']


opts, xargs = getopt.getopt(sys.argv[1:],'',["bashrcpath=","jobinfofilepath="])
for o, a in opts:
    if o in ("--bashrcpath"):
        bashrcpath=a
    elif o in ("--jobinfofilepath"):
        jobinfofilepath=a


def RemoveDeadAndAlreadyActiveNodes(nodelist,programexceptionlist):
    newnodelist=[]
    for node in nodelist:
        keepnode=False
        exceptionstring=''
        for exception in programexceptionlist:
            exceptionstring+=exception+'|'
        exceptionstring=exceptionstring[:-1]
        cmdstr='pgrep '+exceptionstring
        output=CheckOutputFromExternalNode(node,cmdstr)
        if output==False:
            keepnode=True
        else:
            cmdstr1='ps -p %s'%(output)
            cmdstr2='ps -p %s'%(output)+' -u'
            output1=CheckOutputFromExternalNode(node,cmdstr1)
            output2=CheckOutputFromExternalNode(node,cmdstr2)
            if type(output1)==str:
                WriteToLogFile(output1)
            if type(output2)==str:
                WriteToLogFile(output2)


        cmdstr='nohup ls > temp.out'
        process = subprocess.Popen(cmdstr,shell=True)
        nodedead=CheckForDeadNode(process,node) 
        if nodedead==True:
            WriteToLogFile('node '+node+' is dead')
            keepnode=False
        if keepnode==True:
            newnodelist.append(node)
        else:
            WriteToLogFile('cannot use node '+node+' because of program exception')
    return newnodelist
            
def CheckOutputFromExternalNode(node,cmdstr):
    output=True
    if node[-1].isdigit() and '-' in node:
        node=node[:-2]
    job='ssh %s "%s"'%(node,cmdstr)

    try: # if it has output that means this process is running
         output=subprocess.check_output(job,stderr=subprocess.STDOUT,shell=True)
         output=ConvertOutput(output)
    except: #if it fails, that means no process with the program is running or node is dead/unreachable
         output=False
    return output 

def ReadNodeList(nodelistfilepath):
    nodelist=[]
    if os.path.isfile(nodelistfilepath):
        temp=open(nodelistfilepath,'r')
        results=temp.readlines()
        for line in results:
            if '#' not in line:
                newline=line.replace('\n','')
                linesplit=newline.split()
                nodelist.append(linesplit[0])

        temp.close()
    return nodelist

def chunks(lst, n):
    ls=[]
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        ls.append(lst[i:i + n])
    return ls

def DistributeJobsToNodes(nodes,jobs,jobtoscratchspace,prevnodetojoblist):
    nodetojoblist={}
    nodetojoblist.update(prevnodetojoblist)
    if len(nodes)!=0:
        if len(jobs)>len(nodes):
            remainder=len(jobs) % len(nodes)
            wholejobs=len(jobs)-remainder
            jobspernode=int(wholejobs/len(nodes))
        else:
            jobspernode=1
            wholejobs=len(jobs)
        jobspernodelist=chunks(jobs[:wholejobs],jobspernode)
        newnodes=[]
        for i in range(len(jobspernodelist)):
           joblist=jobspernodelist[i]
           node=nodes[i]
           for job in joblist:
               scratchspaceneeded=jobtoscratchspace[job]
               if scratchspaceneeded!=None:
                   scratchavail=CheckScratchSpace(node)
                   isitenough=CheckIfEnoughScratch(scratchspaceneeded,scratchavail)
               else:
                   isitenough=True
               if isitenough==False:
                   WriteToLogFile('not enough scratch space for job = '+job+' on node '+node)
               else:
                   newnodes.append(node)
        if len(newnodes)!=0:
            if len(jobs)>len(newnodes):
                remainder=len(jobs) % len(newnodes)
                wholejobs=len(jobs)-remainder
                jobspernode=int(wholejobs/len(newnodes))
            else:
                wholejobs=len(jobs)
                jobspernode=1
            jobspernodelist=chunks(jobs[:wholejobs],jobspernode)
           
            for i in range(len(jobspernodelist)):
               joblist=jobspernodelist[i]
               node=newnodes[i]
               for job in joblist:
                   if node not in nodetojoblist.keys():
                       nodetojoblist[node]=[]
                   nodetojoblist[node].append(job)
            if len(jobs)>len(newnodes):
                res = list(islice(reversed(jobs), 0, remainder)) 
                res.reverse() 
                nodetojoblist=AssignJobForEveryNode(res,newnodes,nodetojoblist,jobtoscratchspace)
        return nodetojoblist


def AssignJobForEveryNode(jobs,nodes,nodetojoblist,jobtoscratchspace):
    newnodes=[]
    for i in range(len(jobs)):
        node=nodes[i]
        job=jobs[i]
        scratchspaceneeded=jobtoscratchspace[job]
        if scratchspaceneeded!=None:
            scratchavail=CheckScratchSpace(node)
            isitenough=CheckIfEnoughScratch(scratchspaceneeded,scratchavail)
        else:
            isitenough=True
        if isitenough==False:
            WriteToLogFile('not enough scratch space for job = '+job+' on node '+node)
            continue
        else:
            newnodes.append(node)
    for i in range(len(newnodes)):
        node=newnodes[i]
        job=jobs[i]
        if node not in nodetojoblist.keys():
            nodetojoblist[node]=[]
        nodetojoblist[node].append(job)
    return nodetojoblist


def CheckScratchSpace(node):
    cmdstr='df -h'
    job='ssh %s "%s"'%(node,cmdstr)
    p = subprocess.Popen(job, stdout=subprocess.PIPE,shell=True)
    output = p.communicate()
    lines=output[0].decode("utf-8").split('\n')[1:-1]
    d={}
    for line in lines:
        linesplit=line.split()
        if len(linesplit)==5 or len(linesplit)==6:
            avail = re.split('\s+', line)[3]
            mount = re.split('\s+', line)[5]
            d[mount] = avail
    if '/scratch' in d.keys(): 
        scratchavail=d['/scratch']
    else:
        scratchavail='0G'
        WriteToLogFile(' node '+node+' has no scratch')
    if scratchavail==False:
        cmdstr="du -h /scratch | sort -n -r | head -n 15"
        output=CheckOutputFromExternalNode(node,cmdstr)
        WriteToLogFile(output)

    return scratchavail


def CheckIfEnoughScratch(scratchinput,scratchspace):
    enoughscratchspace=False
    if scratchinput==None:
        enoughscratchspace=True
        return enoughscratchspace
    inputspace,inputunit=SplitScratch(scratchinput)
    availspace,availunit=SplitScratch(scratchspace)
    if inputunit=='G' and availunit=='M' or (inputunit=='T' and availunit=='M') or (inputunit=='T' and availunit=='G'):
        pass
    elif inputunit==availunit:
        if float(inputspace)<float(availspace):
            enoughscratchspace=True
    elif  (inputunit=='M' and availunit=='G'):
        inputspace=float(inputspace)*.001
        availspace=float(availspace)
        if float(inputspace)<float(availspace):
            enoughscratchspace=True
    elif  (inputunit=='M' and availunit=='T'):
        inputspace=float(inputspace)*.000001
        availspace=float(availspace)
        if float(inputspace)<float(availspace):
            enoughscratchspace=True
    elif  (inputunit=='G' and availunit=='T'):
        inputspace=float(inputspace)*.001
        availspace=float(availspace)
        if float(inputspace)<float(availspace):
            enoughscratchspace=True
     
       
    return enoughscratchspace

def SplitScratch(string):
    for eidx in range(len(string)):
        e=string[eidx]
        if not e.isdigit() and e!='.':
            index=eidx
            break
    space=string[:index]
    diskunit=string[index]
    return space,diskunit



def CheckIfPreviousJobsFinished(jobtoprocess,previousjobs,finishedjoblist,jobtologhandle,node,polledjobs):
    previousjobsfinished=True
    for job in previousjobs:
        if job not in jobtoprocess.keys() and job not in polledjobs:
            previousjobsfinished=False
        else:
            loghandle=jobtologhandle[job]
            finishedjoblist,term,polledjobs=PollProcess(jobtoprocess,job,finishedjoblist,loghandle,node,polledjobs)
            if term==False:
                previousjobsfinished=False
    return previousjobsfinished,finishedjoblist

def WriteToLogFile(string,loghandle=None):
    now = time.strftime("%c",time.localtime())
    if loghandle!=None:
        loghandle.write(now+' '+string+'\n')
        loghandle.flush()
        os.fsync(loghandle.fileno())
    masterloghandle.write(now+' '+string+'\n')
    masterloghandle.flush()
    os.fsync(masterloghandle.fileno())


def CallSubprocess(node,jobpath,bashrcpath,job,loghandle,wait=False):
    if node[-1].isdigit() and '-' in node:
        node=node[:-2]
        cardvalue=node[-1]
        job=SpecifyGPUCard(cardvalue,job)
    cmdstr = 'ssh %s "cd %s;source %s;%s"' %(str(node),jobpath,bashrcpath,job)
    process = subprocess.Popen(cmdstr, stdout=loghandle,stderr=loghandle,shell=True)
    WriteToLogFile('Calling: '+cmdstr,loghandle)
    nodedead=False
    if wait==True: # grab output from subprocess:
        nodedead=CheckForDeadNode(process,node)    
    return process,nodedead

def CheckForDeadNode(process,node):
    nodedead=False
    output, err = process.communicate()
    if process.returncode != 0:
        err=ConvertOutput(err)
        WriteToLogFile(err+' '+'on node '+node)
        nodedead=True    
    return nodedead

def ConvertOutput(output):
    if output!=None:
        output=output.rstrip()
        if type(output)!=str:
            output=output.decode("utf-8")
    return output

def MakeScratch(node,jobpath,bashrcpath,loghandle,scratchdir):
    cmdstr='[ -d "%s" ] && echo "Directory Exists"'%(scratchdir)
    output=CheckOutputFromExternalNode(node,cmdstr)
    while CheckOutputFromExternalNode(node,cmdstr)==False:
        mkstr='mkdir '+scratchdir
        process,nodedead=CallSubprocess(node,jobpath,bashrcpath,mkstr,loghandle)                            
        time.sleep(2)


def SubmitJob(node,jobpath,bashrcpath,job,loghandle,jobtoprocess,jobtoscratchdir,jobinfo):
    if job in jobtoscratchdir.keys():
        scratchdir=jobtoscratchdir[job]
        if scratchdir!=None:
            MakeScratch(node,jobpath,bashrcpath,loghandle,scratchdir)
        process,nodedead=CallSubprocess(node,jobpath,bashrcpath,job,loghandle)
        jobtoprocess[job]=process
        RemoveJobInfoFromQueue(jobinfo,jobtoprocess)
    return jobtoprocess

def PollProcess(jobtoprocess,job,finishedjoblist,loghandle,node,polledjobs):
    process=jobtoprocess[job]
    poll=process.poll()
    polledjobs.append(job)
    if poll!=None:
        if job not in finishedjoblist:
            finishedjoblist.append(job)
            WriteToLogFile(job+' '+'has terminated on node '+node,loghandle)
        for program in restrictedprogramtonumber.keys():
            if program in job:
                plist=currentrestrictedprogramtoprocesslist[program]
                if process in plist:
                    currentrestrictedprogramtoprocesslist[program].remove(process)
        term=True
    else:
        for program in restrictedprogramtonumber.keys():
            if program in job:
                plist=currentrestrictedprogramtoprocesslist[program]
                if process not in plist:
                    currentrestrictedprogramtoprocesslist[program].append(process) 
        WriteToLogFile(job+' '+'has not terminated on node '+node,loghandle)
        term=False
    return finishedjoblist,term,polledjobs


def SubmitJobs(cpunodetojoblist,gpunodetojoblist,bashrcpath,sleeptime,jobtologhandle,jobinfo):
    jobnumber=len(jobtologhandle.keys())
    jobtoprocess={}
    finishedjoblist=[]
    while len(finishedjoblist)!=jobnumber:
        jobinfo,jobtoprocess,finishedjoblist=SubmitJobsLoop(cpunodetojoblist,jobtologhandle,jobinfo,jobtoprocess,finishedjoblist)
        jobinfo,jobtoprocess,finishedjoblist=SubmitJobsLoop(gpunodetojoblist,jobtologhandle,jobinfo,jobtoprocess,finishedjoblist)
        time.sleep(sleeptime)
        WriteToLogFile('*************************************')
        cpunodes,gpucards=GrabCPUGPUNodes()
        jobinfo=ReadTempJobInfoFiles(jobinfo)
        AddJobInfoToDictionary(jobinfo,jobtoinfo,jobtoprocess)
        jobinfo=ReadJobInfoFromFile(jobinfo,jobtoinfo)
        cpujobs,gpujobs=PartitionJobs(jobinfo,cpuprogramlist,gpuprogramlist)
        jobtologhandle=CreateNewLogHandles(jobinfo['logname'],jobtologhandle)
        cpunodetojoblist=DistributeJobsToNodes(cpunodes,cpujobs,jobinfo['scratchspace'],cpunodetojoblist)
        gpunodetojoblist=DistributeJobsToNodes(gpucards,gpujobs,jobinfo['scratchspace'],gpunodetojoblist)
        jobnumber=len(jobtologhandle.keys())
    WriteToLogFile('All jobs have finished ')

def SubmitJobsLoop(nodetojoblist,jobtologhandle,jobinfo,jobtoprocess,finishedjoblist):
    polledjobs=[]
    for node in nodetojoblist.keys():
        joblist=nodetojoblist[node]
        for i in range(len(joblist)):
            job=joblist[i]
            loghandle=jobtologhandle[job]
            logname=jobinfo['logname'][job]
            jobpath,tail=os.path.split(logname)
            submit=True
            if job not in jobtoprocess.keys():

                previousjobs=joblist[:i]
                previousjobsfinished,finishedjoblist=CheckIfPreviousJobsFinished(jobtoprocess,previousjobs,finishedjoblist,jobtologhandle,node,polledjobs)
                if previousjobsfinished==False:
                    submit=False
                for program,number in restrictedprogramtonumber.items():
                    if program in job:
                        plist=currentrestrictedprogramtoprocesslist[program]
                        currentnumber=len(plist)
                        if currentnumber>=number:
                            submit=False
                if submit==True:
                    jobtoprocess=SubmitJob(node,jobpath,bashrcpath,job,loghandle,jobtoprocess,jobinfo['scratch'],jobinfo)
            elif job in jobtoprocess.keys() and job not in polledjobs:
                finishedjoblist,term,polledjobs=PollProcess(jobtoprocess,job,finishedjoblist,loghandle,node,polledjobs)
    return jobinfo,jobtoprocess,finishedjoblist

def GPUCardToNode(gpucards):
    gpucardtonode={}
    for gpucard in gpucards:
        gpunode=gpucard[:-2]
        gpucardtonode[gpucard]=gpunode
    return gpucardtonode

def PruneBadNodes(gpucardtonode,gpunodes):
    gpucards=[]
    for gpucard,gpunode in gpucardtonode.items():
        if gpunode in gpunodes:
            gpucards.append(gpucard)
    return gpucards

def SpecifyGPUCard(cardvalue,job):
    string='export CUDA_VISIBLE_DEVICES='+str(cardvalue)
    job=string+';'+job
    return job


def RemoveJobInfoFromQueue(jobinfo,jobtoprocess):
    newjobinfo=RemoveAlreadySubmittedJobs(jobtoprocess,jobinfo)
    WriteOutJobInfo(newjobinfo,jobtoinfo,jobtoprocess)

def RemoveAlreadySubmittedJobs(jobtoprocess,jobinfo):
    newjobinfo={}
    for key in jobinfo.keys():
        d=jobinfo[key]
        if key not in newjobinfo.keys():
            newjobinfo[key]={}
        for job in d.keys():
            if job not in jobtoprocess.keys():
                newjobinfo[key][job]=d[job]
    return newjobinfo

def WriteOutJobInfo(jobinfo,filepath,jobtoprocess):
    bufsize=1
    if os.path.isfile(filepath):
        os.remove(filepath)
    temp=open(filepath,'w',buffering=bufsize)
    jobtologname=jobinfo['logname']
    jobtoscratch=jobinfo['scratch']
    jobtoscratchspace=jobinfo['scratchspace']
    for job,log in jobtologname.items():
        if job in jobtoprocess.keys():
            continue
        scratch=jobtoscratch[job]
        scratchspace=jobtoscratchspace[job]
        temp.write('--job='+job+' '+'--outputlogpath='+log+' '+'--scratchdir='+scratch+' '+'--scratchspace='+scratchspace+'\n')
        temp.flush()
        os.fsync(temp.fileno())

    temp.close()


def AddJobInfoToDictionary(jobinfo,filepath,jobtoprocess):
    jobinfoprev=ReadJobInfoFromFile(jobinfo,jobtoinfo)
    for key in jobinfo.keys():
        prevd=jobinfo[key]
        jobinfo[key].update(prevd)
    WriteOutJobInfo(jobinfo,filepath,jobtoprocess)


def ParseJobInfo(line):
    linesplit=line.split('--')[1:]
    linesplit=[e.rstrip() for e in linesplit]
    job=None
    logname=None
    scratch=None
    scratchspace=None
    for line in linesplit:
        if "job=" in line:
            job=line.replace('job=','')
        elif "outputlogpath=" in line:
            logname=line.replace('outputlogpath=','')
        elif "scratchdir=" in line:
            scratch=line.replace('scratchdir=','')
        elif "scratchspace=" in line:
            scratchspace=line.replace('scratchspace=','')
    return job,logname,scratch,scratchspace


def ReadTempJobInfoFiles(jobinfo):
    curdir=os.getcwd()
    os.chdir(writepath)
    files=os.listdir()
    dellist=[]
    for f in files:
        if '_TEMP' in f:
            jobinfo=ReadJobInfoFromFile(jobinfo,f)
            dellist.append(f)
    for f in dellist:
        if os.path.isfile(f):
            os.remove(f)
    return jobinfo

   
def CreateNewLogHandles(jobtologname,jobtologhandle):
    createdloghandles=[]
    lognametologhandle={}
    for job,logname in jobtologname.items():
        if job in jobtologhandle.keys():
            continue
        if logname not in createdloghandles:
            createdloghandles.append(logname)
            loghandle=open(logname,'w')
            lognametologhandle[logname]=loghandle
        else:
            loghandle=lognametologhandle[logname]
        jobtologhandle[job]=loghandle
    return jobtologhandle

def ReadJobInfoFromFile(jobinfo,filename):
    if os.path.isfile(filename):
        temp=open(filename,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            job,log,scratch,scratchspace=ParseJobInfo(line)
            if job==None or log==None:
                continue
            jobinfo['logname'][job]=log
            jobinfo['scratch'][job]=scratch
            jobinfo['scratchspace'][job]=scratchspace
    return jobinfo



def WritePIDFile():
    pid=str(os.getpid())
    temp=open(pidfile, 'w',buffering=1)
    temp.write(pid+'\n')
    temp.flush()
    os.fsync(temp.fileno())
    temp.close()

def PartitionJobs(jobinfo,cpuprogramlist,gpuprogramlist):
    cpujobs=[]
    gpujobs=[]
    d=jobinfo['logname']
    for job in d.keys():
        for program in cpuprogramlist:
            if program in job:
                cpujobs.append(job)
        for program in gpuprogramlist:
            if program in job:
                gpujobs.append(job)
    return cpujobs,gpujobs

def GrabCPUGPUNodes():
    WriteToLogFile('*************************************')
    WriteToLogFile("Checking available CPU nodes")
    cpunodes=ReadNodeList(cpunodelistfilepath)
    cpunodes=RemoveDeadAndAlreadyActiveNodes(cpunodes,cpuprogramexceptionlist)
    gpucards=ReadNodeList(gpucardlistfilepath)
    gpucardtonode=GPUCardToNode(gpucards)
    gpunodes=list(set(gpucardtonode.values()))
    WriteToLogFile('*************************************')
    WriteToLogFile("Checking available GPU cards")
    gpunodes=RemoveDeadAndAlreadyActiveNodes(gpunodes,gpuprogramexceptionlist)
    gpucards=PruneBadNodes(gpucardtonode,gpunodes)
    return cpunodes,gpucards

jobinfo={}
jobinfo['logname']={}
jobinfo['scratch']={}
jobinfo['scratchspace']={}
jobinfo=ReadJobInfoFromFile(jobinfo,jobinfofilepath)
jobtoprocess={}
if os.path.isfile(pidfile):
    head,tail=os.path.split(jobinfofilepath)
    tempfilepath=writepath+tail.replace('.txt','_TEMP.txt')
    WriteOutJobInfo(jobinfo,tempfilepath,jobtoprocess)
    sys.exit()
else:
    WritePIDFile()
    AddJobInfoToDictionary(jobinfo,jobtoinfo,jobtoprocess)
    try:
        jobtologhandle={}
        cpunodetojoblist={} 
        gpunodetojoblist={}
        jobinfo=ReadJobInfoFromFile(jobinfo,jobtoinfo)
        cpujobs,gpujobs=PartitionJobs(jobinfo,cpuprogramlist,gpuprogramlist)
        jobtologhandle=CreateNewLogHandles(jobinfo['logname'],jobtologhandle)
        cpunodes,gpucards=GrabCPUGPUNodes()
        cpunodetojoblist=DistributeJobsToNodes(cpunodes,cpujobs,jobinfo['scratchspace'],cpunodetojoblist)
        gpunodetojoblist=DistributeJobsToNodes(gpucards,gpujobs,jobinfo['scratchspace'],gpunodetojoblist)
        SubmitJobs(cpunodetojoblist,gpunodetojoblist,bashrcpath,sleeptime,jobtologhandle,jobinfo)
    finally:
        if os.path.isfile(pidfile):
            os.unlink(pidfile)
    
