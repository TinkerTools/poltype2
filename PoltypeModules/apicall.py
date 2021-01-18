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
    if poltype.bashrcpath!=None:
        cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
    else:
        cmdstr='python'+' '+poltype.externalapi+' '+'--jobinfofilepath='+jobinfofilepath
    print('cmdstr '+cmdstr,flush=True)
    poltype.call_subsystem(cmdstr,True)

