import os
import sys
import subprocess

def CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,jobinfofilepath,makejobfileonly,jobtooutputfilepath=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(jobinfofilepath,'w')
    for job,inputfilepaths in jobtoinputfilepaths.items():
        if jobtooutputfilepath!=None:
            outputfilepath=jobtooutputfilepath[job]
        else:
            head,tail=os.path.split(inputfilepaths[0])
            outputfilepath=head
        outputfiles=jobtooutputfiles[job]
        outputfiles=[os.path.join(outputfilepath,i) for i in outputfiles]
        inputfilestr=','.join(inputfilepaths)
        outputfilestr=','.join(outputfiles)
        binpath=jobtoabsolutebinpath[job] 
        string='--job='
        if poltype.runjobslocally==True:
            string+=' cd '+outputfilepath+' ; '
        string+=job+' '
        string+='--numproc='+str(1)+' '+'--disk=0GB'+' '+'--inputfilepaths='+inputfilestr
        if poltype.runjobslocally==False:
            string+=' '+'--outputfilepaths='+outputfilestr
        if '9 ' in job:
            string+=' '+'--gpujob'
            string+=' --ram=0GB'
        else:
            string+=' --ram=3GB'
        if poltype.email!=None:
            string+=' '+'--email='+poltype.email
        string+='\n'
        temp.write(string)

    temp.close()
    if makejobfileonly==False:
        if poltype.printjobsleft==True:
            sys.exit()

        
        if poltype.externalapi!=None and poltype.submitlocally==False:
            if poltype.bashrcpath!=None:
                cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
            else:
                cmdstr='python'+' '+poltype.externalapi+' '+'--jobinfofilepath='+jobinfofilepath

            poltype.WriteToLog('Calling external API ')
            call_subsystem(poltype,cmdstr,wait=False,skiperrors=False)



def SubmitJobs(poltype,jobtolog,jobtojobpath,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,jobinfofilepath,jobtooutputfilepath=None,makejobfileonly=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    if len(jobtolog.keys())!=0:
        CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,jobinfofilepath,makejobfileonly,jobtooutputfilepath)
    if makejobfileonly==True:
        sys.exit()
    if poltype.submitlocally==True:
        CallJobsSeriallyLocalHost(poltype,jobtolog,jobtojobpath,jobtooutputfiles)


def CallJobsSeriallyLocalHost(poltype,jobtolog,jobtojobpath,jobtooutputfiles):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    jobs=list(jobtolog.keys())
    for jobidx in range(len(jobs)):
        job=jobs[jobidx]
        remainingjobs=len(jobs)-(jobidx+1)+1
        outputfiles=jobtooutputfiles[job]
        for outputfile in outputfiles:
            if 'Solv' in outputfile:
                jobtype='Solv'
            elif 'Comp' in outputfile:
                jobtype='Comp'
            else:
                continue
            totaltime,outputfiletype,freq,serial=poltype.GrabOutputFileTypeAndTime(outputfile,jobtype,jobidx,jobs)
            ETA=totaltime
            ETA_string=poltype.ETAString(ETA)
            poltype.WriteToLog('Job: '+job+' ETA: '+ETA_string)
            totalETA=ETA*remainingjobs
            ETA_string=poltype.ETAString(totalETA)
            poltype.WriteToLog('Jobtype: '+outputfiletype+' ETA: '+ETA_string)
            

        percentcomplete=round((jobidx+1)*100/len(jobs),2)
        poltype.WriteToLog('Percent of jobs complete = '+str(percentcomplete)+'%')
        jobpath=jobtojobpath[job]
        os.chdir(jobpath)
        call_subsystem(poltype,job,wait=True)
    os.chdir(curdir)

def call_subsystem(poltype,cmdstr,wait=False,skiperrors=False,outputfilename=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog("Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
    p = subprocess.Popen(cmdstr, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if wait==True:
        saved=[]
        for line in p.stdout:
            poltype.WriteToLog(line)
            saved.append(line)
        p.wait()
        if outputfilename==None:
            check=ParseForErrorsInOutPut(poltype,saved)
        else:
            check=ParseForErrorsInOutPut(poltype,[],outputfilename)
        if (p.returncode != 0 and skiperrors==False) or check==False:
            poltype.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
            raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())

def ParseForErrorsInOutPut(poltype,saved,outputfilename=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    check=True
    if outputfilename==None:
        pass
    else:
        temp=open(outputfilename,'r')
        saved=temp.readlines()
        temp.close()
 
    for line in saved:
        if not isinstance(line, str):
            line=line.decode("utf-8")
        if 'Tinker is Unable to Continue' in line or 'Terminating with uncaught exceptio' in line or 'OLDATM  --  A PDB Atom of Biotype' in line:
            poltype.WriteToLog(line)
            print(line)
            check=False 

    return check
