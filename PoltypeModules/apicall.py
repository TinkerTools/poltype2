import os

def CallExternalAPI(poltype,listofcommands,scriptname,scratchdir):
    poltype.WriteToLog('Calling external API ')
    shellname=scriptname+'.sh'
    temp=open(shellname,'w')
    temp.write('source '+poltype.bashrcpath+'\n')
    for cmd in listofcommands:
        temp.write(cmd+'\n')
    temp.close()
    os.system('chmod +x '+shellname)
    cmdstr='python'+' '+poltype.externalapi+' '+poltype.bashrcpath+' '+os.getcwd()+' '+scratchdir+' '+shellname
    os.system(cmdstr)



   
