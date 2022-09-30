import os
import shutil
   
def CopyKeyFileAsBackup(poltype,keypath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    maxnumber=FindCurrentLastKeyFile(poltype,keypath)
    newkeyfile=keypath+'_'+str(maxnumber+1)
    if not os.path.exists('BackupKeys'):
        os.mkdir('BackupKeys')
    shutil.copy(keypath,os.path.join('BackupKeys',newkeyfile))
    head,newkeyfile=os.path.split(newkeyfile)
    return newkeyfile

def FindCurrentLastKeyFile(poltype,keypath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    files=os.listdir()
    head,tail=os.path.split(keypath)
    maxnumber=0
    for f in files:
        if tail in f:
            if '_'==f[-2] or '_'==f[-3]:
                split=f.split('_')
                number=int(split[-1])
            else:
                number=0
            if number>maxnumber:
                maxnumber=number
    return maxnumber
        
 
def AddKeyWord(poltype,keypath,string):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    #newkeyfile=CopyKeyFileAsBackup(poltype,keypath)
    poltype.WriteToLog('Adding key words to '+keypath+' '+string,prin=True)
    #poltype.WriteToLog('Making copy of old keyfile '+newkeyfile,prin=True)
    read=open(keypath,'r')
    results=read.readlines()
    read.close()
    tempkeyname=keypath.replace('.key','-t.key')
    temp=open(tempkeyname,'w')
    temp.write(string)
    for line in results:
        temp.write(line)
    temp.close()
    os.remove(keypath)
    os.rename(tempkeyname,keypath)

def CheckIfStringAlreadyInKeyfile(poltype,keypath,string):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    found=False
    temp=open(keypath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if string in line:
            found=True
    return found
    
def RemoveKeyWords(poltype,keypath,keystringlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog('Removing key words from '+keypath+' '+str(keystringlist),prin=True)
    read=open(keypath,'r')
    results=read.readlines()
    read.close()
    tempname=keypath.replace('.key','-t.key')
    temp=open(tempname,'w')
    for line in results:
        write=True
        for keystring in keystringlist:
            if keystring in line:
                write=False
        if write==True:
            temp.write(line)
    temp.close()
    os.remove(keypath)
    os.rename(tempname,keypath)  
    
def InsertKeyfileHeader(poltype,keyfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog('Adding header to key file '+keyfilename,prin=True)
    string='parameters '+poltype.prmfilepath+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='archive'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='integrator '+poltype.integrator+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
       AddKeyWord(poltype,keyfilename,string)
    string='thermostat '+poltype.thermostat+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='ewald'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='vdw-cutoff '+str(poltype.vdwcutoff)+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='ewald-cutoff '+str(poltype.ewaldcutoff)+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='polar-eps '+str(poltype.polareps)+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='polar-predict'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='barostat'+' '+poltype.barostatmethod+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='neighbor-list'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='vdw-correction'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='vdw-annihilate'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='pme-grid 64 64 64'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='pme-order 5'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='polarization MUTUAL'+'\n'
    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)
    string='heavy-hydrogen'+'\n'
    if poltype.heavyhyd==True:
       if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
           AddKeyWord(poltype,keyfilename,string)
    if poltype.externalapi==None:
        string='OPENMP-THREADS'+' '+str(poltype.numproc)+'\n'
    else:
        string='OPENMP-THREADS'+' '+str(1)+'\n'

    if not CheckIfStringAlreadyInKeyfile(poltype,keyfilename,string):
        AddKeyWord(poltype,keyfilename,string)




def AddMultipoleDefinitionsForIonIndices(poltype,ionindexes,charge,keyfilepath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for ionindex in ionindexes:
        AddIonMultipoleDefinition(poltype,ionindex,charge,keyfilepath)

def AddIonMultipoleDefinition(poltype,ionindex,charge,keyfilepath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(keyfilepath,'a')
    chgline='multipole '+'-'+str(ionindex)+'   0   0                '+str(charge)+'\n'
    temp.write(chgline)
    temp.write('                                       0.00000   0.00000  -0.00000'+'\n')
    temp.write('                                       0.00000'+'\n')
    temp.write('                                       0.00000  -0.00000'+'\n')
    temp.write('                                      -0.00000   0.00000  -0.00000'+'\n')
    temp.close()


