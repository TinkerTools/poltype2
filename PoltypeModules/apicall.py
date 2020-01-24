import os
import sys
import subprocess

def CallExternalAPI(poltype,listofcommands,scriptname,scratchdir):
    poltype.WriteToLog('Calling external API ')
    shellname=scriptname+'.sh'
    temp=open(shellname,'w')
    temp.write('source '+poltype.bashrcpath+'\n')
    for cmd in listofcommands:
        temp.write(cmd+'\n')
    temp.close()
    os.system('chmod +x '+shellname)
    cmdstr='python'+' '+poltype.externalapi+' '+'--bashrcpath='+poltype.bashrcpath+' '+'--jobpath='+os.getcwd()+' '+'--scratchdir='+scratchdir+' '+'--shellname='+shellname
    p = subprocess.Popen(cmdstr, shell=True,stdout=poltype.logfh, stderr=poltype.logfh)
    if p.wait() != 0:
        poltype.WriteToLog("ERROR: " + cmdstr)
        sys.exit(1)




   
