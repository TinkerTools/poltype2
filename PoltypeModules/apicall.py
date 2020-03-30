import os
import sys
import subprocess

def CallExternalAPI(poltype,jobtolog,jobinfofilenameprefix,scratchdir):
    poltype.WriteToLog('Calling external API ')
    jobinfofilepath=jobinfofilenameprefix+'.txt'
    temp=open(jobinfofilepath,'w')
    for job,log in jobtolog.items():
        temp.write('--job='+job+' '+'--outputlogpath='+log+' '+'--scratchdir='+scratchdir+' '+'--scratchspace='+poltype.maxdisk+'\n')
    temp.close()
    cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
    p = subprocess.Popen(cmdstr, shell=True,stdout=poltype.logfh, stderr=poltype.logfh)
    if p.wait() != 0:
        poltype.WriteToLog("ERROR: " + cmdstr)
        sys.exit(1)

