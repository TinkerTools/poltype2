import keyfilemodifications as keymods
import tables
import numpy as np
import itertools
import os
import submitjobs as submit
import numpy as np
import shutil
import sys
from rdkit.Chem import rdmolfiles

def ComputeBoxSize(poltype):
    longestdimx=[]
    longestdimy=[]
    longestdimz=[]

    for i in range(len(poltype.xyzfilename)):
        xyzfilename=poltype.xyzfilename[i]
        if poltype.complexation==True and poltype.solvation==True:

            if i==0: # complexation
                x,y,z=FindDimensionsOfMoleculeTinkerXYZ(poltype,poltype.outputpath+xyzfilename)
            elif i==1: # solvation

                longestdim=FindDimensionsOfMoleculeTinker(poltype,poltype.outputpath+xyzfilename)
                x=longestdim
                y=longestdim
                z=longestdim

        elif poltype.complexation==False and (poltype.solvation==True or poltype.neatliquidsim==True):
            longestdim=FindDimensionsOfMoleculeTinker(poltype,poltype.outputpath+xyzfilename)
            x=longestdim
            y=longestdim
            z=longestdim
        longestdimx.append(x)
        longestdimy.append(y)
        longestdimz.append(z)


             
    poltype.WriteToLog('Longest Dimension '+str(longestdimx)+' '+str(longestdimy)+' '+str(longestdimz))
   

    if poltype.aaxis==None and poltype.baxis==None and poltype.caxis==None:
        poltype.aaxislist=[]
        poltype.baxislist=[]
        poltype.caxislist=[]
        for i in range(len(longestdimx)):
            xdim=longestdimx[i]
            ydim=longestdimy[i]
            zdim=longestdimz[i]
            buffervalue=poltype.bufferlen[i]
            xvalue=round(xdim+buffervalue,1)
            yvalue=round(ydim+buffervalue,1)
            zvalue=round(zdim+buffervalue,1)

            if poltype.solvation==True and poltype.salthfe==False: # ion component only, then want water box to be big.
                xvalue=max(xvalue,50)
                yvalue=max(yvalue,50)
                zvalue=max(zvalue,50)


            poltype.aaxislist.append(xvalue)
            poltype.baxislist.append(yvalue)
            poltype.caxislist.append(zvalue)


    else:
        poltype.aaxislist=[poltype.aaxis,poltype.aaxis]
        poltype.baxislist=[poltype.baxis,poltype.baxis]
        poltype.caxislist=[poltype.caxis,poltype.caxis]

    poltype.WriteToLog('Initial Box Length'+' '+str(poltype.aaxislist)+' '+str(poltype.baxislist)+' '+str(poltype.caxislist))
    for i in range(len(poltype.tabledict)):
        poltype.tabledict[i]['Initial Box Length']=str(poltype.aaxislist[i])+','+str(poltype.baxislist[i])+','+str(poltype.caxislist[i])

    tables.WriteTableUpdateToLog(poltype)
    poltype.volume=[]
    for i in range(len(poltype.configkeyfilename)):
        keylist=poltype.configkeyfilename[i]
        for j in range(len(keylist)):
            keyname=keylist[j]
            aaxis=poltype.aaxislist[i]
            keymods.AddKeyWord(poltype,poltype.outputpath+keyname,'a-axis '+str(aaxis)+'\n')
            baxis=poltype.baxislist[i]
            if baxis!=None:
                keymods.AddKeyWord(poltype,poltype.outputpath+keyname,'b-axis '+str(baxis)+'\n')
            caxis=poltype.caxislist[i]
            if caxis!=None:
                keymods.AddKeyWord(poltype,poltype.outputpath+keyname,'c-axis '+str(caxis)+'\n')
            volume=aaxis*baxis*caxis
            poltype.volume.append(volume)   

 
def ComputePhysiologicalIonNumber(poltype):
    
    if poltype.addphysioions==True:
        if poltype.binding==True:
            commasplit=poltype.listofsaltcons.split(',')
            for ioncomplex in commasplit:
                ioncomplex=ioncomplex.lstrip().rstrip()
                equalsignsplit=ioncomplex.split('=')
                ioncomplexstring=equalsignsplit[0].lstrip().rstrip()
                conc=float(equalsignsplit[1].lstrip().rstrip()) # in mM
                complexnum=[conc*6.022*10**-7*i for i in poltype.volume]
                Sum=[0,0]
                for el in poltype.elementsymtotinktype.keys():
                    if el in ioncomplexstring:
                        nextindex=ioncomplexstring.find(el)+len(el) # if K returns K index +1, if Cl returns Cl index +2, want to know if there is 2 or 3... after element, nextindex should always be defined because last elemt will be ]
                        if ioncomplexstring[nextindex].isdigit():
                            multfactor=int(ioncomplexstring[nextindex])
                        else:
                            multfactor=1
                        tinktype=poltype.elementsymtotinktype[el]
                        for i in range(len(poltype.tabledict)):
                            if el not in poltype.tabledict[i].keys():
                                poltype.tabledict[i][el]=0
                            poltype.tabledict[i][el]+=int(round(complexnum[i]*multfactor))
                            Sum[i]+=int(round(complexnum[i]*multfactor))
                            if tinktype in poltype.iontypetoionnumberphysio[i].keys():
                                poltype.iontypetoionnumberphysio[i][tinktype]+=int(round(complexnum[i]*multfactor))
                            else:
                                poltype.iontypetoionnumberphysio[i][tinktype]=int(round(complexnum[i]*multfactor))
        if poltype.solvation==True and poltype.binding==False:
            for i in range(len(poltype.tabledict)):
                poltype.tabledict[i]['Physio Counterions']=0
        else:
            for i in range(len(poltype.tabledict)):
                poltype.tabledict[i]['Physio Counterions']=Sum[i]
        tables.WriteTableUpdateToLog(poltype)

def ComputeNeutralizingIonNumber(poltype,systemcharge): # this can be just ligand charge or ligand and receptor charge depending on whether this is complexation or solvation    
    tinktypetoelementsym={v: k for k, v in poltype.elementsymtotinktype.items()}
    for i in range(len(systemcharge)):
        chg=systemcharge[i]

        if poltype.solvation==True and poltype.binding==False:
            chg=0
            if i>0:
                continue
        else:
            if chg>0:
                cltinktype=poltype.elementsymtotinktype['Cl']
                poltype.iontypetoionnumberneut[i][cltinktype]=np.abs(int(chg))
                if 'Cl' not in poltype.tabledict[i].keys():
                    poltype.tabledict[i]['Cl']=0
                poltype.tabledict[i]['Cl']+=np.abs(int(chg))

            elif chg<0:
                ktinktype=poltype.elementsymtotinktype['K']
                poltype.iontypetoionnumberneut[i][ktinktype]=np.abs(int(chg))
                if 'K' not in poltype.tabledict[i].keys():
                    poltype.tabledict[i]['K']=0
                poltype.tabledict[i]['K']+=np.abs(int(chg))

        poltype.tabledict[i]['Neut Counterions']=np.abs(int(chg))
    tables.WriteTableUpdateToLog(poltype)


def ComputeWaterNumber(poltype):
    poltype.waternum=[int(round(.0334*i)) for i in poltype.volume] # number of waters per cubic angstroms
    poltype.WriteToLog('Water Number '+str(poltype.waternum))

def RemoveBoxSizeTerms(poltype,keyfile):
    keymods.RemoveKeyWords(poltype,keyfile,['axis'])


def FindDimensionsOfMoleculeTinker(poltype,structurefilepath):
    veclist=[]
    temp=open(structurefilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            vec=np.array([float(linesplit[2]),float(linesplit[3]),float(linesplit[4])])
            veclist.append(vec)

    pairs=list(itertools.combinations(veclist, 2))
    disttodiffvec={}
    
    for pairidx in range(len(pairs)):
        pair=pairs[pairidx]
        progress=(pairidx*100)/len(pairs)
        diff=np.array(pair[0])-np.array(pair[1])
        dist=np.linalg.norm(diff)
        
        disttodiffvec[dist]=diff
    distlist=list(disttodiffvec.keys())
    if len(distlist)!=0:
        mindist=np.amax(np.array(distlist))
        diffvec=disttodiffvec[mindist]
    else:
        mindist=0
    return mindist

def FindDimensionsOfMoleculeTinkerXYZ(poltype,structurefilepath):
    veclist=[]
    temp=open(structurefilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            vec=np.array([float(linesplit[2]),float(linesplit[3]),float(linesplit[4])])
            veclist.append(vec)
    xdistlist=[]
    ydistlist=[]
    zdistlist=[]
    for i in range(len(veclist)):
        vec=veclist[i]
        x,y,z=vec[:]
        xdistlist.append(x)
        ydistlist.append(y)
        zdistlist.append(z)
    maxx=max(xdistlist)
    minx=min(xdistlist)
    x=np.abs(maxx-minx) 
    maxy=max(ydistlist)
    miny=min(ydistlist)
    y=np.abs(maxy-miny) 
    maxz=max(zdistlist)
    minz=min(zdistlist)
    z=np.abs(maxz-minz) 


    return x,y,z





def TotalAtomNumber(poltype,xyzfilename):
    atomnum=[FindNumberTinkerXYZAtoms(poltype,i) for i in xyzfilename]
    poltype.totalatomnumberxyzfilename=atomnum
    totalatomnum=[]
    for i in range(len(atomnum)):
        atoms=atomnum[i]
        waters=poltype.waternum[i]
        tot=atoms+3*waters
        totalatomnum.append(tot)
        poltype.tabledict[i]['Total Atom Number']=tot
    tables.WriteTableUpdateToLog(poltype)



def FindNumberTinkerXYZAtoms(poltype,xyzfilename):
    temp=open(poltype.outputpath+xyzfilename,'r')
    results=temp.readlines()
    temp.close()
    atomnum=int(results[0].replace('\n','').rstrip().lstrip())
    return atomnum

def GrabTypeNumber(poltype,typedescrip):
    temp=open(poltype.prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if typedescrip in line:
            linesplit=line.split()
            typenum=linesplit[1]
    return typenum
   
def CreateWaterXYZ(poltype):
    temp=open(poltype.outputpath+'water.xyz','w')
    temp.write('3'+'\n')
    temp.write('1  O      8.019814    5.892935    0.449481   '+str(poltype.waterOtypenum)+'     2     3'+'\n')
    temp.write('2  H      7.906681    4.942582    0.463215   '+str(poltype.waterHtypenum)+'     1'+'\n')
    temp.write('3  H      8.149325    6.169963   -0.457590   '+str(poltype.waterHtypenum)+'     1'+'\n')
    temp.close() 

      

def CreateSolventBox(poltype,aaxis,baxis,caxis,waternum,filename,key):
    temp=open(poltype.outputpath+'xyzedit.in','w')
    temp.write('23'+'\n')
    temp.write(str(waternum)+'\n')
    if baxis==None:
        baxis=aaxis
    if caxis==None:
        caxis=aaxis
    temp.write(str(aaxis)+','+str(baxis)+','+str(caxis)+'\n')
    temp.write('Y'+'\n')
    temp.close()
    cmdstr=poltype.xyzeditpath+' '+filename+' '+'-k'+' '+key+' <'+' '+'xyzedit.in'
    submit.call_subsystem(poltype,cmdstr,wait=True)    
    newfilename=filename+'_2'
    return newfilename
 
def SoakMoleculeInSolventBox(poltype,xyzfilename,keyfilename):
    cmdstr=poltype.xyzeditpath+' '+xyzfilename+' '+'-k'+' '+keyfilename+' '+'24'+' '+'water.xyz_2'
    submit.call_subsystem(poltype,cmdstr,wait=True)    

def AddIonToSolventBox(poltype,solutexyzfilename,keyfilename,tinktype,ionnum,count,writesolute=True):
    soluteatomnum=FindNumberTinkerXYZAtoms(poltype,solutexyzfilename)
    num=2+count
    inputfile='xyzedit_'+str(num)+'.in'
    temp=open(poltype.outputpath+inputfile,'w')
    if writesolute==True:
        temp.write('1'+' '+str(soluteatomnum)+'\n')
    else:
        temp.write('0'+'\n')
    string=''
    string+=str(tinktype)+' '+str(ionnum)+' '
    string+='\n'
    temp.write(string)
    temp.write('\n')
    temp.close()
    cmdstr=poltype.xyzeditpath+' '+solutexyzfilename+'_'+str(num)+' '+'-k'+' '+keyfilename+' '+'25'+' '+' < '+inputfile
    submit.call_subsystem(poltype,cmdstr,wait=True)    



def AddIonsToSolventBox(poltype,solutexyzfilename,keyfilename,boxfilename,prevcount,iontypetoionnumberneut,iontypetoionnumberphysio,writesolute=True):
    count=0
    for tinktype in iontypetoionnumberneut.keys():
        ionnum=iontypetoionnumberneut[tinktype]
        AddIonToSolventBox(poltype,solutexyzfilename,keyfilename,tinktype,ionnum,count,writesolute)
        count+=1
    for tinktype in iontypetoionnumberphysio.keys():
        ionnum=iontypetoionnumberphysio[tinktype]
        AddIonToSolventBox(poltype,solutexyzfilename,keyfilename,tinktype,ionnum,count,writesolute)
        count+=1
    os.rename(solutexyzfilename+'_'+str(count+prevcount),boxfilename)
    


def RemoveTempFiles(poltype):
    files=os.listdir()
    for f in files:
        if '.' in f:
            fsplit=f.split('.')
            end=fsplit[1]
            if 'xyz_' in end or 'in' in end and 'ini' not in end:
                os.remove(f) 


def GrabIonIndexes(poltype,ionnumber,boxfilename,iontypenumber):
    ionindexes=[]
    count=0     
    temp=open(boxfilename,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if count==ionnumber:
            break
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            index=int(linesplit[0])
            typenum=int(linesplit[5]) 
            if typenum==int(iontypenumber):
                ionindexes.append(index)
                count+=1 
    return ionindexes

def GrabTotalMass(xyzfilename):
    symboltomass = {'H': 1.007825, 'C': 12.01, 'O': 15.9994, 'N': 14.0067, 'S': 31.972071,'P': 30.973762,'F':18.998,'Cl':35.453,'Br':79.904,'I':126.904}
    temp=open(xyzfilename,'r')
    results=temp.readlines()
    temp.close()
    atomicsyms=[]
    for line in results:
        linesplit=line.split()
        if len(linesplit)>1:
            atomicsym=linesplit[1]
            atomicsyms.append(atomicsym)
    mass=0
    for atomicsym in atomicsyms:
        mass+=symboltomass[atomicsym]

    return mass


def BoxSetupProtocol(poltype):

    poltype.WriteToLog('Computing volume ',prin=True)
    ComputeBoxSize(poltype)

    if poltype.addsolvionwindows==True: # create box file for ions
        ionaxis=poltype.bufferlen[0]
        ionvolume=ionaxis**3
        ionwaternum=int(round(.0334*ionvolume)) 
        poltype.volume.append(ionvolume)
    for i in range(len(poltype.configkeyfilename)):
        keylist=poltype.configkeyfilename[i]
        for j in range(len(keylist)):
            key=keylist[j]
            if not os.path.isfile(poltype.outputpath+key):
                string='a-axis '+str(poltype.aaxislist[i])+'\n'
                keymods.AddKeyWord(poltype,key,string)
    if poltype.neatliquidsim==False:
        ComputeWaterNumber(poltype)
        ComputePhysiologicalIonNumber(poltype)
        ComputeNeutralizingIonNumber(poltype,poltype.systemcharge)
        TotalAtomNumber(poltype,poltype.xyzfilename)
        CreateWaterXYZ(poltype)
    for i in range(len(poltype.boxfilename)):
        boxxyzfilename=poltype.boxfilename[i]
        xyzfilename=poltype.xyzfilename[i]
        key=poltype.configkeyfilename[i][0]
        aaxis=poltype.aaxislist[i]
        baxis=poltype.baxislist[i]
        caxis=poltype.caxislist[i]
        if not os.path.isfile(poltype.outputpath+boxxyzfilename):
            if poltype.usepreequilibriatedbox==False:
                if poltype.neatliquidsim==False:
                    newfilename=CreateSolventBox(poltype,aaxis,baxis,caxis,poltype.waternum[i],'water.xyz',key)
                else:
                    mass=GrabTotalMass(xyzfilename)
                    mass=mass*1.66054*10**(-27) # convert daltons to Kg
                    ls=[aaxis,baxis,caxis]
                    ls=[axis*10**-10 for axis in ls]
                    numbermolecules=int(poltype.density*(ls[0]*ls[1]*ls[2])/mass)
                    newfilename=CreateSolventBox(poltype,aaxis,baxis,caxis,numbermolecules,xyzfilename,key)
                    os.rename(newfilename,boxxyzfilename)

            else:
                boxsize=[aaxis,baxis,caxis]
                poltype.TrimPreEquilibriatedBox(boxsize)
            if poltype.neatliquidsim==False:
                SoakMoleculeInSolventBox(poltype,xyzfilename,key)
                AddIonsToSolventBox(poltype,xyzfilename,key,boxxyzfilename,2,poltype.iontypetoionnumberneut[i],poltype.iontypetoionnumberphysio[i])
        RemoveTempFiles(poltype)
        
    poltype.xyzfilesize=[float(os.path.getsize(i)) for i in poltype.boxfilename] # in bytes
    poltype.equilarcfilesize=[i*poltype.equilframenum*10**-9 for i in poltype.xyzfilesize] # in GB
    poltype.singlepertubationfilesize=[i*int(poltype.proddynframenum)*10**-9 for i in poltype.xyzfilesize] # in GB
    poltype.totalperturbationfilesize=[len(poltype.estatlambdascheme)*i for i in poltype.singlepertubationfilesize] # in GB
    poltype.totalfilesize=[]
    for i in range(len(poltype.equilarcfilesize)):
        equilsize=poltype.equilarcfilesize[i]
        perturbsize=poltype.totalperturbationfilesize[i]
        tot=equilsize+perturbsize
        poltype.totalfilesize.append(tot)
        poltype.tabledict[i]['Prod MD Arc File Space']=tot
    solv=False
    if poltype.complexation==False and poltype.solvation==True:
        solv=True
        index=0
    elif poltype.complexation==True and poltype.solvation==True:
        solv=True
        index=1
    if solv==True:
        for i in range(poltype.totalatomnumberxyzfilename[index]):
            idx=i+1
            poltype.ligandindices[index].append(idx)
    if poltype.complexation==True:
        uncomplexedatomnum=poltype.totalatomnumberxyzfilename[0]-poltype.totalatomnumberxyzfilename[1]
        ligandnum=poltype.totalatomnumberxyzfilename[1]
        indices=np.arange(uncomplexedatomnum+1,ligandnum+uncomplexedatomnum+1,1)
        poltype.ligandindices[0]=indices


    if poltype.addsolvionwindows==True:
        shutil.copy(poltype.configkeyfilename[0][0],poltype.ionkeyfilename)
        keymods.RemoveKeyWords(poltype,poltype.ionkeyfilename,['axis'])
        string='a-axis '+str(ionaxis)+'\n'
        keymods.AddKeyWord(poltype,poltype.ionkeyfilename,string)
        newfilename=CreateSolventBox(poltype,ionaxis,ionaxis,ionaxis,ionwaternum,'water.xyz',poltype.ionkeyfilename)
        filename='water.xyz'
        AddIonsToSolventBox(poltype,filename,poltype.ionkeyfilename,poltype.ionboxfilename,2,poltype.solviontocount,poltype.iontypetoionnumberphysio[-1],False)

    

    if poltype.boxonly==True:
        sys.exit() 
        
