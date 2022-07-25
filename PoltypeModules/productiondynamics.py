import submitjobs
import keyfilemodifications as keymods
import os
import shutil
import boxsetup as box
import terminate as term
import sys
import time
import numpy as np
import mutation as mutate
import equilbriation as equil
import mdtraj as md


def ExecuteProductionDynamics(poltype):
   """
   Intent:
   Input:
   Output:
   Referenced By: 
   Description: 
   """
   jobtolog={}
   jobtojobpath={}
   jobtoinputfilepaths={}
   jobtooutputfiles={}
   jobtoabsolutebinpath={}
   for j in range(len(poltype.simfoldname)):
       simfoldname=poltype.simfoldname[j]
       os.chdir(poltype.outputpath+simfoldname)
       subfolders=os.listdir(os.getcwd())
       proddynboxfilename=poltype.proddynboxfilename[j]
       proddynboxkeyfilename=poltype.proddynboxkeyfilename[j]
       proddynarcboxfilename=poltype.proddynarcboxfilename[j]
       proddynoutfilepathlistoflist=poltype.proddynoutfilepath[j]
       for k in range(len(proddynoutfilepathlistoflist)):
           proddynoutfilepath=proddynoutfilepathlistoflist[k]
           for i in range(len(proddynoutfilepath)):
               outputfilepath=proddynoutfilepath[i]
               if 'Gas' in outputfilepath:
                   dynamicpath=poltype.dynamicpath
                   proddyntimestep=1
                   proddynsteps=str(int((poltype.proddyntime*1000000)/proddyntimestep)) # hard code .1 for Gas phase no shielding from box
                   nvt=True
                   ensemble=2
               else:
                   dynamicpath=poltype.truedynamicpath
                   proddynsteps=poltype.proddynsteps
                   proddyntimestep=poltype.proddyntimestep
                   nvt=poltype.productiondynamicsNVT
                   ensemble=poltype.proddynensem
               path,tail=os.path.split(outputfilepath)
               os.chdir(path)
               fold=os.path.basename(os.path.normpath(path))
               keyfilename=proddynboxkeyfilename.replace(poltype.foldername+'_',poltype.foldername+'_'+fold+'_').replace('0.','0-')
               boxfilename=proddynboxfilename.replace(poltype.foldername+'_',poltype.foldername+'_'+fold+'_').replace('0.','0-')

               arcfilename=proddynarcboxfilename.replace(poltype.foldername+'_',poltype.foldername+'_'+fold+'_').replace('0.','0-')


               cmdstr=ProductionDynamicsCommand(poltype,boxfilename,keyfilename,proddynsteps,ensemble,outputfilepath,dynamicpath,proddyntimestep,nvt)
               terminate,deletefile,error=term.CheckFileTermination(poltype,outputfilepath,float(poltype.proddynsteps))
               if terminate==False:
                   if os.path.exists(os.path.join(path,arcfilename)):
                       stepstaken,dynoutfiles=CheckLastNumberDynamicsStepsCompletedAllTimes(poltype)
                       newfile=GenerateBackupDynFile(poltype,dynoutfiles)
                       shutil.copy(outputfilepath,newfile)
                       newstepstotake=int(proddynsteps)-stepstaken
                       cmdstr=ProductionDynamicsCommand(poltype,arcfilename,keyfilename,newstepstotake,ensemble,outputfilepath,dynamicpath,proddyntimestep,nvt)
                   if '_gpu' in cmdstr and poltype.cpujobsonly==True:
                       pass
                   else:
                       jobtolog[cmdstr]=poltype.outputpath+poltype.proddynjobsfilename
                       jobtojobpath[cmdstr]=path
                       xyzname=arcfilename.replace('.arc','.xyz')
                       dynname=arcfilename.replace('.arc','.dyn')
                       inputfilepaths=[os.path.join(path,xyzname),os.path.join(path,keyfilename),poltype.outputpath+poltype.prmfilepath]
                       arcpath=os.path.join(path,arcfilename)
                       dynpath=os.path.join(path,dynname)
                       if os.path.isfile(arcpath):
                           inputfilepaths.append(arcpath)
                       if os.path.isfile(dynpath):
                           inputfilepaths.append(dynpath)
                       head,tail=os.path.split(outputfilepath)
                       outputfilenames=[tail,arcfilename,dynname]
                       absbinpath=poltype.which(dynamicpath)
                       jobtoinputfilepaths[cmdstr]=inputfilepaths
                       jobtooutputfiles[cmdstr]=outputfilenames
                       jobtoabsolutebinpath[cmdstr]=absbinpath
   

               os.chdir('..')
       os.chdir('..')
   if len(jobtolog.keys())!=0:
       if poltype.equilonly==True:
           submitjobs.SubmitJobs(poltype,jobtolog,jobtojobpath,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,poltype.outputpath+poltype.proddynjobsfilename,makejobfileonly=True)
       else:
           poltype.WriteToLog('Running production dynamics',prin=True)
           submitjobs.SubmitJobs(poltype,jobtolog,jobtojobpath,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,poltype.outputpath+poltype.proddynjobsfilename)

   if poltype.equilonly==True:
       sys.exit()

def ProductionDynamicsCommand(poltype,inputxyzname,keyfile,steps,ensemble,outputpath,dynamicpath,proddyntimestep,nvt):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    head,tail=os.path.split(outputpath) 
    if nvt==True:
        cmdstr=dynamicpath+' '+inputxyzname+' -k ' + keyfile + ' '+str(steps)+' '+str(proddyntimestep)+' '+str(poltype.proddynwritefreq)+' '+str(ensemble)+' '+str(poltype.equilibriatescheme[-1])+' '+' > '+tail
    else:
        cmdstr=dynamicpath+' '+ inputxyzname+' '+ '-k'+' '+ keyfile+' '+str(steps)+' '+ str(proddyntimestep)+' '+ str(poltype.proddynwritefreq)+' '+str(ensemble)+' '+str(poltype.equilibriatescheme[-1])+' '+str(poltype.pressure)+' '+' > '+tail

    return cmdstr
   
def ModifyKeyForGasPhase(poltype,keyfilepath,changeligandindices,index):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keymods.RemoveKeyWords(poltype,keyfilepath,['axis','ewald','cutoff','correction'])
    
    if poltype.changegasphaseintegrator==True:
        keymods.RemoveKeyWords(poltype,keyfilepath,['integrator'])
        string='integrator stochastic'+'\n'
        if not keymods.CheckIfStringAlreadyInKeyfile(poltype,keyfilepath,string):
            keymods.AddKeyWord(poltype,keyfilepath,string)
    else:
        string='RESPA-INNER 0.1'+'\n'
        if not keymods.CheckIfStringAlreadyInKeyfile(poltype,keyfilepath,string):
            keymods.AddKeyWord(poltype,keyfilepath,string)
    if changeligandindices==True:
        keymods.RemoveKeyWords(poltype,keyfilepath,['ligand'])
        AddLigandIndices(poltype,[poltype.ligandindices[index]],keyfilepath)


def SearchNearestNonPerturbedFolder(poltype,folderlist,index):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    fold=folderlist[index]
    invertedlist=folderlist[::-1]
    invertedindex=invertedlist.index(fold)
    if 'Gas' in fold:
        for k in range(invertedindex+1,len(invertedlist)):
            newfold=invertedlist[k]
            if 'Key' not in newfold:
                break
    else:
        for k in range(index+1,len(folderlist)):
            newfold=folderlist[k]
            if 'Key' not in newfold:
                break
    return newfold




def AppendAllXYZ(poltype,xyzfilenamelist,key):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    finalxyz='alllig.xyz'
    if len(xyzfilenamelist)==1:
        return xyzfilenamelist[0]
    firstxyz=xyzfilenamelist[0]
    for i in range(1,len(xyzfilenamelist)):
        xyz=xyzfilenamelist[i]
        firstxyz=AppendXYZ(poltype,firstxyz,xyz,key)
    shutil.copy(firstxyz,finalxyz)

    return finalxyz



def AppendXYZ(poltype,firstxyz,xyz,key):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.outputpath+'xyzedit.in','w')
    temp.write('22'+'\n')
    temp.write(str(xyz)+'\n')
    temp.write('\n')
    temp.close()
    cmdstr=poltype.xyzeditpath+' '+firstxyz+' '+'-k'+' '+key+' <'+' '+'xyzedit.in'
    submit.call_subsystem(poltype,cmdstr,wait=True)    
    if '.xyz_' in firstxyz:
        split=firstxyz.split('.xyz_')
        num=split[-1]
        num=int(num)
        newxyz=split[:-1]+str(num)
    else:
        newxyz=firstxyz+'_2'

    return newxyz



def SetupProductionDynamics(poltype,simfoldname,lambdafolderlist,index,proddynboxfilename,lambdakeyfilenamelist,xyzfilenamelist):
   """
   Intent:
   Input:
   Output:
   Referenced By: 
   Description: 
   """
   poltype.WriteToLog('Setting up dynamics for '+simfoldname,prin=True)
   if poltype.addsolvionwindows==True:
       string='ele-lambda'+'\n'
       keymods.AddKeyWord(poltype,poltype.ionkeyfilename,string)
       string='vdw-lambda'+'\n'
       keymods.AddKeyWord(poltype,poltype.ionkeyfilename,string)
       typelist=list(poltype.solviontocount.keys())
       iontypenumber=typelist[0]
       ioncount=poltype.solviontocount[iontypenumber]
       ionindexes=box.GrabIonIndexes(poltype,ioncount,poltype.ionproddynboxfilename,iontypenumber)
       string='ligand'+' '
       for idx in ionindexes:
           string+=str(idx)+','
       string=string[:-1]
       string+='\n'
       keymods.AddKeyWord(poltype,poltype.ionkeyfilename,string)
   if not os.path.isdir(poltype.outputpath+simfoldname):
       os.mkdir(poltype.outputpath+simfoldname)
   os.chdir(poltype.outputpath+simfoldname)
   if len(poltype.mutlambdascheme)==0:
       mut=False
   else:
       mut=True
   for k in range(len(lambdafolderlist)):
       folderlist=lambdafolderlist[k]
       estatlambdascheme=poltype.estatlambdascheme[k]
       vdwlambdascheme=poltype.vdwlambdascheme[k]
       restlambdascheme=poltype.restlambdascheme[k]
       foldtonextfold={}
       for i in range(len(folderlist)-1):
           fold=folderlist[i]
           if 'Gas' in fold:
               nextfold=folderlist[i-1]
           else:
               nextfold=folderlist[i+1]
           foldtonextfold[fold]=nextfold
       for i in range(len(folderlist)):
           fold=folderlist[i]
           if fold in foldtonextfold.keys():
               nextfold=foldtonextfold[fold]
           if mut==False:
               elelamb=estatlambdascheme[i]
               vdwlamb=vdwlambdascheme[i]
           else:
               mutlambda=poltype.mutlambdascheme[i]
           if poltype.complexation==True and poltype.restrainreceptorligand==True and index==0:
               reslambda=restlambdascheme[i]
           else:
               reslambda=0
           if not os.path.isdir(fold):
               os.mkdir(fold)
           os.chdir(fold)

           files=os.listdir()
           for f in files:
               if '.end' in f or '.err' in f:
                   os.remove(f) # just in case option to stop all simulations was used
               if '.out' in f and "BAR" not in f:
                   outputfile=os.path.join(os.getcwd(),f) 
                   terminate,deletefile,error=term.CheckFileTermination(poltype,outputfile,float(poltype.proddynsteps))
                   #if terminate==True:
                   #    continue
           newfoldpath=os.getcwd()+'/'
           newtempkeyfile=proddynboxfilename.replace('.xyz','.key').replace(poltype.foldername+'_',poltype.foldername+'_'+fold+'_').replace('0.','0-')

           if mut==False:
               if 'Ion' in fold:
                   shutil.copyfile(poltype.outputpath+poltype.ionkeyfilename,os.path.join(newfoldpath+newtempkeyfile))

               elif 'Key' in fold:
                   numinterpols=poltype.foldernametonuminterpols[fold]
                   interpolindex=poltype.foldernametointerpolindex[fold]
                   if numinterpols==1 or i==0 or i==len(folderlist)-1:
                       shutil.copyfile(poltype.outputpath+poltype.foldernametolambdakeyfilename[fold],os.path.join(newfoldpath+newtempkeyfile))
                   else:
                       thiskey=poltype.outputpath+poltype.foldernametolambdakeyfilename[fold]
                       nextkey=poltype.outputpath+poltype.foldernametolambdakeyfilename[nextfold]
                       poltype.bgnstatekey=thiskey
                       poltype.endstatekey=nextkey
                       poltype.bgnstatexyz=poltype.ligandxyzfilenamelist[0]
                       poltype.endstatexyz=poltype.ligandxyzfilenamelist[0]
                       mutate.SingleTopologyMutationProtocol(poltype)

                       arrayoflinearrays=mutate.MutateAllParameters(poltype,poltype.bgnlinetoendline,None,numinterpols,interpolindex)
                       mutate.GenerateKeyFile(poltype,arrayoflinearrays,os.path.join(newfoldpath+newtempkeyfile))
                       keymods.InsertKeyfileHeader(poltype,os.path.join(newfoldpath+newtempkeyfile))

               else:
                   shutil.copyfile(poltype.outputpath+lambdakeyfilenamelist[0],os.path.join(newfoldpath+newtempkeyfile))
               if 'Gas' in fold:

                   if poltype.binding==True:
                       changeligandindices=True
                   else:
                       changeligandindices=False
                   ModifyKeyForGasPhase(poltype,os.path.join(newfoldpath+newtempkeyfile),changeligandindices,index)

           else:
               arrayoflinearrays=mutate.MutateAllParameters(poltype,poltype.bgnlinetoendline,mutlambda,None,None)
               mutate.GenerateKeyFile(poltype,arrayoflinearrays,os.path.join(newfoldpath+newtempkeyfile))
               keymods.InsertKeyfileHeader(poltype,os.path.join(newfoldpath+newtempkeyfile))

           outputboxname=os.path.join(newfoldpath+proddynboxfilename.replace(poltype.foldername+'_',poltype.foldername+'_'+fold+'_').replace('0.','0-'))

           outputarcname=outputboxname.replace('.xyz','.arc')
           if 'Key' in fold:
               newfold=SearchNearestNonPerturbedFolder(poltype,folderlist,i)
               curdir=os.getcwd()
               os.chdir('..')
               os.chdir(newfold)
               otherfoldpath=os.getcwd()+'/'
               os.chdir(curdir)
               tempoutputboxname=os.path.join(otherfoldpath+proddynboxfilename.replace(poltype.foldername+'_',poltype.foldername+'_'+fold+'_').replace('0.','0-'))
               lastarcpath=tempoutputboxname.replace('.xyz','.arc')
               outputsignalfile=outputarcname.replace('.arc','.fep')
               if poltype.fep==True:
                   if not os.path.isfile(outputarcname):
                       poltype.WriteToLog('Copying arc file from  '+lastarcpath+' to '+outputarcname,prin=True)
                       shutil.copy(lastarcpath,outputarcname)
                       with open(outputsignalfile, 'w') as fp:
                           pass
               else:
                   if os.path.isfile(outputsignalfile):
                       os.remove(outputsignalfile)
                       if os.path.isfile(outputarcname):
                           os.remove(outputarcname)
           shutil.copy(os.path.join(poltype.outputpath,poltype.prmfilepath),os.getcwd())

           if 'Gas' not in fold and 'Ion' not in fold:
 
               shutil.copyfile(poltype.outputpath+proddynboxfilename,outputboxname)

           elif 'Ion' not in fold and 'Gas' in fold:
               if poltype.binding==False:
                   finalxyz=AppendAllXYZ(poltype,xyzfilenamelist,newfoldpath+newtempkeyfile)
                   shutil.copyfile(poltype.outputpath+finalxyz,outputboxname)
               else:
                   liquidfoldernames=lambdafolderlist[0]
                   liquidindex=len(liquidfoldernames)-1-i
                   liquidfolder=liquidfoldernames[liquidindex]
                   ligandindices=poltype.ligandindices[index]
                   solvligandindices=poltype.ligandindices[1]
                   ligandindicestosolvligandindices=dict(zip(ligandindices,solvligandindices))
                   ExtractLigandDynamicsFromLiquidBox(poltype,outputboxname,liquidfolder,ligandindices,ligandindicestosolvligandindices) 
           else:
               shutil.copyfile(poltype.outputpath+poltype.ionproddynboxfilename,outputboxname)

           if mut==False:
               ModifyLambdaKeywords(poltype,newfoldpath,newtempkeyfile,elelamb,vdwlamb,reslambda)
           if poltype.binding==True and 'Gas' not in fold:
               xyzpath=outputboxname
               keypath=newfoldpath+newtempkeyfile
               alzout=newfoldpath+'checknetcharge.alz'
               poltype.CheckNetChargeIsZero(xyzpath,keypath,alzout) 
           os.chdir('..')
   os.chdir('..')


def ExtractLigandDynamicsFromLiquidBox(poltype,outputboxname,liquidfolder,ligandindices,ligandindicestosolvligandindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    outputarcname=outputboxname.replace('.xyz','.arc')
    arcexists=False
    curdir=os.getcwd()
    os.chdir('..')
    os.chdir(liquidfolder)
    files=os.listdir()
    for f in files:
        if '.arc' in f:
            arcfile=f
            arcexists=True
    if arcexists==True:
        poltype.WriteToLog('Reading in ARC for Gas Phase Extraction '+arcfile,prin=True)
        indextovecs={}
        t = md.load_arc(arcfile)
        for i in range(t.n_frames):
            for index in ligandindices:
                zeroindex=index-1
                vec=t.xyz[i,zeroindex,:]
                vec=vec*10
                trueindex=ligandindicestosolvligandindices[index]
                if trueindex not in indextovecs.keys():
                    indextovecs[trueindex]=[]
                indextovecs[trueindex].append(vec)

        totalatomnum=len(ligandindices)
        WriteOutArcFile(poltype,totalatomnum,indextovecs,poltype.ligandindextoneighbslist,poltype.ligandindextosymlist,poltype.ligandindextotypenumlist,outputarcname,t.n_frames)


    os.chdir(curdir)



def MergeMaps(poltype,diclist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    dic={}
    for d in diclist:
        dic.update(d)

    return dic




def WriteOutArcFile(poltype,totalatomnum,indextovecs,ligandindextoneighbslist,ligandindextosymlist,ligandindextotypelist,outputarcname,framenum):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ligandindextoneighbs=MergeMaps(poltype,ligandindextoneighbslist)
    ligandindextosym=MergeMaps(poltype,ligandindextosymlist)
    ligandindextotype=MergeMaps(poltype,ligandindextotypelist)
    temp=open(outputarcname,'w')
    for i in range(framenum):
        temp.write(str(totalatomnum)+'\n')
        for index in indextovecs.keys():
            sym=ligandindextosym[index]
            neighbs=ligandindextoneighbs[index]
            typenum=ligandindextotype[index]
            vecs=indextovecs[index]
            vec=vecs[i]
            neighbstring=''
            for neighb in neighbs:
                neighbstring+=str(neighb)+' '
            string=str(index)+' '+sym+' '+str(vec[0])+' '+str(vec[1])+' '+str(vec[2])+' '+str(typenum)+' '+neighbstring +'\n'
            temp.write(string)

    temp.close()




def ModifyLambdaKeywords(poltype,newfoldpath,newtempkeyfile,elelamb,vdwlamb,reslambda):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(newfoldpath+newtempkeyfile,'r')
    results=temp.readlines()
    temp.close()
    newkeyfile=open(newfoldpath+newtempkeyfile,'w')
    for line in results:
        linesplit=line.split()
        if "ele-lambda" in line:
            newline=line.replace('\n'," "+str(elelamb)+'\n')
            newkeyfile.write(newline)
        elif "vdw-lambda" in line:
            newline=line.replace('\n'," "+str(vdwlamb)+'\n')
            newkeyfile.write(newline)
        elif "restrain-groups 1 2" in line:
            group1index=linesplit[1]
            group2index=linesplit[2]
            constant=linesplit[3]
            teatherdist=linesplit[-1]
            constant=str(float(constant)*float(reslambda))
            if poltype.flatbotrest==False:
                restrainstring='restrain-groups '+str(group1index)+' '+ str(group2index)+' '+str(constant)+' '+str(teatherdist)+' '+str(teatherdist)
            else:
                restrainstring='restrain-groups '+str(group1index)+' '+ str(group2index)+' '+str(constant)+' '+'0'+' '+str(teatherdist)

            
            newkeyfile.write(restrainstring+'\n')      
        elif 'restrain-distance' in line:
            index1=linesplit[1]
            index2=linesplit[2]
            constant=linesplit[3]
            constant=str(float(constant)*float(reslambda))
            lowerdist=linesplit[4]
            upperdist=linesplit[5]
            string='restrain-distance '+index1+' '+index2+' '+constant+' '+lowerdist+' '+upperdist
            newkeyfile.write(string+'\n')        

        elif 'restrain-angle' in line:
            index1=linesplit[1]
            index2=linesplit[2]
            index3=linesplit[3]
            constant=linesplit[4]
            constant=str(float(constant)*float(reslambda))
            lowerangle=linesplit[5]
            upperangle=linesplit[6]
            string='restrain-angle '+index1+' '+index2+' '+index3+' '+constant+' '+lowerangle+' '+upperangle
            newkeyfile.write(string+'\n')        

        elif 'restrain-torsion' in line:
            index1=linesplit[1]
            index2=linesplit[2]
            index3=linesplit[3]
            index4=linesplit[4]
            constant=linesplit[5]
            constant=str(float(constant)*float(reslambda))
            lowerangle=linesplit[6]
            upperangle=linesplit[7]
            string='restrain-torsion '+index1+' '+index2+' '+index3+' '+index4+' '+constant+' '+lowerangle+' '+upperangle
            newkeyfile.write(string+'\n')        

        elif 'multipole' in line:
            typeindex=linesplit[1]
            if '-' in typeindex:
                xframeindex=linesplit[2]
                yframeindex=linesplit[3]
                if int(xframeindex)==0 and int(yframeindex)==0:
                    chg=float(linesplit[4])
                    chg=chg*(float(elelamb))
                    newline='multipole '+typeindex+' '+xframeindex+' '+yframeindex+' '+str(chg)+'\n'
                    newkeyfile.write(newline)
                else:
                    newkeyfile.write(line)
            else:
                newkeyfile.write(line)

        else:
            newkeyfile.write(line)
        newkeyfile.flush()
        os.fsync(newkeyfile.fileno())
        sys.stdout.flush()

    newkeyfile.close()

 
def CheckLastNumberDynamicStepsCompleted(poltype,outputfilepath):
   """
   Intent:
   Input:
   Output:
   Referenced By: 
   Description: 
   """
   steps=None
   if os.path.isfile(outputfilepath):
       temp=open(outputfilepath,'r')
       results=temp.readlines()
       temp.close()
       for line in results:
           if 'Dynamics Steps' in line:
               linesplit=line.split()
               if linesplit[-3].isdigit():
                   steps=int(linesplit[-3])
       return steps
   return steps


def DetermineIonIndicesToModifyCharge(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    totalionindexes=[]
    totalcharge=[]
    if poltype.addsolvionwindows==False:
        index=len(poltype.ligandchargelist)
    else:
        index=len(poltype.ligandchargelist)-1
    for i in range(index):
        chg=int(poltype.ligandchargelist[i])
        proddynboxfilename=poltype.proddynboxfilename[i]
        ionindexes=[]
        charge=0
        if chg>0:
            Clions=np.absolute(chg)
            iontypenumber=poltype.elementsymtotinktype['Cl']
            ionindexes=box.GrabIonIndexes(poltype,Clions,proddynboxfilename,iontypenumber)
            charge=-1
        elif chg<0:
            Kions=np.absolute(chg)
            iontypenumber=poltype.elementsymtotinktype['K']
            ionindexes=box.GrabIonIndexes(poltype,Kions,proddynboxfilename,iontypenumber)
            charge=1
        totalionindexes.append(ionindexes)
        totalcharge.append(charge)
    return totalionindexes,totalcharge



def GrabTinkerFiles(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    arcpaths=[]
    keypaths=[]
    curdir=os.getcwd()
    for i in range(len(poltype.lambdakeyfilename)): # Comp,Solv or just NeatLiq
        simfoldname=poltype.simfoldname[i]
        os.chdir(poltype.outputpath+simfoldname)
        lambdafolderlist=poltype.lambdafolderlist[i]
        temparc=[]
        tempkey=[]
        for k in range(len(lambdafolderlist)): # Lambda Liq, Lambda Gas
            folderlist=lambdafolderlist[k]
            estatlambdascheme=poltype.estatlambdascheme[k]
            vdwlambdascheme=poltype.vdwlambdascheme[k]
            restlambdascheme=poltype.restlambdascheme[k]
            temparc2=[]
            tempkey2=[]
            for j in range(len(folderlist)): # Lambda List
                fold=folderlist[j]
                os.chdir(fold)
                files=os.listdir()
                for f in files:
                    path=os.path.join(os.getcwd(),f)
                    if '.key' in f:
                        tempkey2.append(path)
                    if '.arc' in f:
                        temparc2.append(path)
                os.chdir('..')
            temparc.append(temparc2)
            tempkey.append(tempkey2)
        arcpaths.append(temparc)
        keypaths.append(tempkey)


    os.chdir(curdir)
    return arcpaths,keypaths


def AddLigandIndices(poltype,ligandindices,keyfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    string='ligand'+' '
    firstindices=ligandindices[0]
    lastindices=ligandindices[-1]
    firstligidx=str(firstindices[0])
    lastligidx=str(lastindices[-1])
    string+='-'+firstligidx+' '+lastligidx

    string+='\n'
    keymods.AddKeyWord(poltype,keyfilename,string)



def ProductionDynamicsProtocol(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    if poltype.didinputsubmitlocally==False:
        poltype.submitlocally=False
 
    if poltype.prodmdfinished==False:
        if len(poltype.mutlambdascheme)==0 and poltype.neatliquidsim==False:
            for i in range(len(poltype.lambdakeyfilename)):
                lambdakeyfilenamelist=poltype.lambdakeyfilename[i]
                configkeyfilenamelist=poltype.configkeyfilename[i]
                for k in range(len(lambdakeyfilenamelist)):
                    lambdakeyfilename=lambdakeyfilenamelist[k]
                    configkeyfilename=configkeyfilenamelist[k]
                    shutil.copyfile(poltype.outputpath+configkeyfilename,poltype.outputpath+lambdakeyfilename)
                    AddLigandIndices(poltype,[poltype.ligandindices[i]],lambdakeyfilename)
                    string='ele-lambda'+'\n'
                    keymods.AddKeyWord(poltype,lambdakeyfilename,string)
                    string='vdw-lambda'+'\n'
                    keymods.AddKeyWord(poltype,lambdakeyfilename,string)
                    if poltype.complexation==True and i==0 and poltype.needrot==True:
                        string='restrain-position '+str(poltype.norotpair[0])+','+str(poltype.norotpair[1])+' '+str(poltype.restrainpositionconstant)+' '+str(poltype.norotrestrainsphereradius)+'\n'
                        keymods.AddKeyWord(poltype,lambdakeyfilename,string)

        if poltype.binding==True:
            totalionindexes,totalcharge=DetermineIonIndicesToModifyCharge(poltype)
        for i in range(len(poltype.lambdakeyfilename)):
            lambdakeyfilenamelist=poltype.lambdakeyfilename[i]
            if poltype.binding==True:
                ionindexes=totalionindexes[i]
                charge=totalcharge[i]
            for j in range(len(lambdakeyfilenamelist)):
                lambdakeyfilename=lambdakeyfilenamelist[j]
                if poltype.binding==True:
                    keymods.AddMultipoleDefinitionsForIonIndices(poltype,ionindexes,charge,lambdakeyfilename)
            simfoldname=poltype.simfoldname[i]
            lambdafolderlist=poltype.lambdafolderlist[i]
            index=i
            proddynboxfilename=poltype.proddynboxfilename[i]
            xyzfilenamelist=poltype.xyzfilename[i]
            proddynoutfilepathlist=poltype.proddynoutfilepath[i]
            for j in range(len(proddynoutfilepathlist)):
                proddynoutfilepath=proddynoutfilepathlist[j]
                SetupProductionDynamics(poltype,simfoldname,lambdafolderlist,index,proddynboxfilename,lambdakeyfilenamelist,xyzfilenamelist)
        
        ExecuteProductionDynamics(poltype)
        messages=[]
        for i in range(len(poltype.proddynoutfilepath)):
            proddynoutfilepathlist=poltype.proddynoutfilepath[i]
            for j in range(len(proddynoutfilepathlist)):
                proddynoutfilepath=proddynoutfilepathlist[j]
                checkfin=term.CheckFilesTermination(poltype,proddynoutfilepath,poltype.proddynsteps)
                finished=checkfin[0]
                percentfinished=checkfin[1]
                while finished==False:
                    msg='Production dynamics is not complete, '+str(percentfinished)+'% of jobs finished'
                    if msg not in messages:
                        poltype.WriteToLog(msg,prin=True)
                        messages.append(msg)
                    time.sleep(poltype.waitingtime)
                    checkfin=term.CheckFilesTermination(poltype,proddynoutfilepath,poltype.proddynsteps)
                    finished=checkfin[0]
                    percentfinished=checkfin[1]
            if poltype.generatepdbtrajs==True and poltype.complexation==True and i==0:
                arcpaths,keypaths=GrabTinkerFiles(poltype)
                liqarcpath=arcpaths[0][0][0]
                pdbfilename=liqarcpath.replace('.arc','.pdb')
                poltype.GeneratePDBFromARC(liqarcpath,pdbfilename)


        poltype.WriteToLog('System dynamics is complete ',prin=True)
    if poltype.alchemical==False:
        sys.exit()
    if poltype.neatliquidsim==True:
        arcpaths,keypaths=GrabTinkerFiles(poltype)
        liqarcpath=arcpaths[0][0][0]
        liqxyzpath=liqarcpath.replace('.arc','.xyz')
        liqkeypath=keypaths[0][0][0]
        gasarcpath=arcpaths[0][1][0]
        gasxyzpath=gasarcpath.replace('.arc','.xyz')
        gaskeypath=keypaths[0][1][0]
        aaxis,baxis,caxis,aaxisarray,baxisarray,caxisarray=equil.AverageBoxSizeFromNPTArc(poltype,liqarcpath,1,10000000,1,None)
        aaxisarray=aaxisarray*10**-10 # convert to m
        baxisarray=baxisarray*10**-10
        caxisarray=caxisarray*10**-10
        aaxis=aaxis**10**-10
        baxis=baxis**10**-10
        caxis=caxis**10**-10
        vavg=aaxis*baxis*caxis
        mass,nummols=GrabSystemMass(poltype,liqxyzpath,liqkeypath)
        avgdensity,stddev=ComputeDensity(poltype,mass,aaxisarray,baxisarray,caxisarray)
        poltype.WriteToLog('Average density is (kg/m^3) '+str(avgdensity)+' +/- '+str(stddev)) 
        avgliqenergy,stddevliqenergy=GrabEnergy(poltype,liqarcpath,liqkeypath)
        avggasenergy,stddevgasenergy=GrabEnergy(poltype,gasarcpath,gaskeypath)
        R = 1.98720425864083 * 0.001 
        T = float(poltype.equilibriatescheme[-1])
        avgliqenergy=avgliqenergy/nummols
        stddevliqenergy=stddevliqenergy/nummols
        Hvap = (-avgliqenergy + avggasenergy + R*T)
        Hvaperr=np.sqrt(np.square(stddevliqenergy)+np.square(stddevgasenergy))
        poltype.WriteToLog('Average Hvap is (kcal/mol) '+str(Hvap)+' +/- '+str(Hvaperr)) 
        sys.exit()



def ComputeDensity(poltype,mass,aaxisarray,baxisarray,caxisarray):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    densityarray=[]
    for i in range(len(aaxisarray)):
        a=aaxisarray[i]
        b=baxisarray[i]
        c=caxisarray[i]
        v=a*b*c
        d=mass/v
        densityarray.append(d)
    densityarray=np.array(densityarray)
    avgdensity=np.mean(densityarray)
    stddev=np.std(densityarray)/np.sqrt(np.size(densityarray))

    return avgdensity,stddev


def GrabSystemMass(poltype,xyzpath,keypath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    alzout='xyz.alz'
    if not os.path.isfile(alzout):
        poltype.CallAnalyze(xyzpath,keypath,alzout,poltype.analyzepath,'G')
    temp=open(alzout,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'Total System Mass' in line:
            mass=float(linesplit[-1])*.001*(1/(6.0221408*10**23)) # kg
        if 'Number of Molecules' in line:
            nummols=int(linesplit[-1])

    return mass,nummols



def GrabEnergy(poltype,arcpath,keypath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    alzout='energy.alz'
    poltype.CallAnalyze(arcpath,keypath,alzout,poltype.trueanalyzepath,'E')
    temp=open(alzout,'r')
    results=temp.readlines()
    temp.close()
    energyarray=[]
    for line in results:
        if 'Total Potential Energy : ' in line:
            linesplit=line.split()
            energy=float(linesplit[-2])
            energyarray.append(energy)
    energyarray=np.array(energyarray)
    avgliqenergy=np.mean(energyarray)
    stddevliqenergy=np.std(energyarray)/np.sqrt(np.size(energyarray))

    return avgliqenergy,stddevliqenergy

def CheckLastNumberDynamicStepsCompleted(poltype,outputfilepath):
   steps=None
   if os.path.isfile(outputfilepath):
       temp=open(outputfilepath,'r')
       results=temp.readlines()
       temp.close()
       for line in results:
           if 'Dynamics Steps' in line:
               linesplit=line.split()
               if linesplit[-3].isdigit():
                   steps=int(linesplit[-3])
       return steps
   return steps


def CheckLastNumberDynamicsStepsCompletedAllTimes(poltype):
    files=os.listdir()
    total=0
    dynoutfiles=[]
    for f in files:
        if '.out' in f:
            stepstaken=CheckLastNumberDynamicStepsCompleted(poltype,f)
            total+=stepstaken
            dynoutfiles.append(f)
    return total,dynoutfiles


def GenerateBackupDynFile(poltype,dynoutfiles):
    highestnum=0
    for f in dynoutfiles:
        split=f.split('.out')
        pre=split[0]
        if pre.isdigit():
            dig=int(pre)
            if dig>highestnum:
                highestnum=dig
    if highestnum==0:
        highestnum=1
    else:
        highestnum+=1
    newfile=str(highestnum)+'.out'
    return newfile


