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
        if 'poltype.ini' in inputfilepaths:
            inputfilestr=inputfilepaths
        else:
            inputfilestr=','.join(inputfilepaths)
        binpath=jobtoabsolutebinpath[job]
        string='--job='+job+' '+'--scratchpath='+scratchdir+' '+'--numproc='+str(poltype.numproc)+' '+'--ram='+poltype.maxmem+' '+'--disk='+poltype.maxdisk+' '+'--inputfilepaths='+inputfilestr+' '+'--username='+poltype.username
        if binpath!=None:
            string+=' '+'--absolutepathtobin='+binpath
        if outputfiles!=None:
            outputfiles=[os.path.join(outputfilepath,i) for i in outputfiles]
            outputfilestr=','.join(outputfiles)
            string+='--outputfilepaths='+outputfilestr

        temp.write(string+'\n')
    temp.close()
    if poltype.bashrcpath!=None:
        cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
    else:
        cmdstr='python'+' '+poltype.externalapi+' '+'--jobinfofilepath='+jobinfofilepath
    poltype.call_subsystem([cmdstr],True)

