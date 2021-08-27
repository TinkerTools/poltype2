import os
import sys
import subprocess

def CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobinfofilenameprefix):
    poltype.WriteToLog('Calling external API ')
    jobinfofilepath=jobinfofilenameprefix+'.txt'
    temp=open(jobinfofilepath,'w')
    head,tail=os.path.split(scratchdir)
    for job,inputfilepaths in jobtoinputfilepaths.items():
        outputfiles=jobtooutputfiles[job]
        inputfilestr=','.join(inputfilepaths)
        outputfilestr=','.join(outputfiles)
        binpath=jobtoabsolutebinpath[job]
        temp.write('--job='+job+' '+'--scratchdir='+tail+' '+'--scratchpath='+head+' '+'--scratchspace='+poltype.maxdisk+' '+'--numproc='+str(poltype.numproc)+' '+'--ram='+poltype.maxmem+' '+'--inputfilepaths='+inputfilestr+' '+'--outputfiles='+outputfilestr+' '+'--absolutepathtobin='+binpath+'\n')
    temp.close()
    if poltype.bashrcpath!=None:
        cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
    else:
        cmdstr='python'+' '+poltype.externalapi+' '+'--jobinfofilepath='+jobinfofilepath
    print('cmdstr '+cmdstr,flush=True)
    poltype.call_subsystem(cmdstr,True)

