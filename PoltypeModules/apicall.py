import os
import sys
import subprocess

def CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobinfofilenameprefix,jobtooutputfilepath={}):
    poltype.WriteToLog('Calling external API ')
    jobinfofilepath=jobinfofilenameprefix+'.txt'
    temp=open(jobinfofilepath,'w')
    for job,inputfilepaths in jobtoinputfilepaths.items():
        if len(jobtooutputfilepath)==0:
            outputfilepath=os.path.split(inputfilepaths[0])[0]
        else:
            outputfilepath=jobtooutputfilepath[job]
        outputfiles=jobtooutputfiles[job]
        outputfiles=[os.path.join(outputfilepath,i) for i in outputfiles]
        inputfilestr=','.join(inputfilepaths)
        outputfilestr=','.join(outputfiles)
        binpath=jobtoabsolutebinpath[job]
        temp.write('--job='+job+' '+'--scratchpath='+scratchdir+' '+'--numproc='+str(poltype.numproc)+' '+'--ram='+poltype.maxmem+' '+'--disk='+poltype.maxdisk+' '+'--inputfilepaths='+inputfilestr+' '+'--outputfilepaths='+outputfilestr+' '+'--absolutepathtobin='+binpath+' '+'--username='+poltype.username+'\n')
    temp.close()
    if poltype.bashrcpath!=None:
        cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
    else:
        cmdstr='python'+' '+poltype.externalapi+' '+'--jobinfofilepath='+jobinfofilepath
    poltype.call_subsystem([cmdstr],True)

