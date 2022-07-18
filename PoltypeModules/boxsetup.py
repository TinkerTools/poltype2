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
import parametercomparison
import time

def ComputeBoxSize(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    longestdimx=[]
    longestdimy=[]
    longestdimz=[]

    for i in range(len(poltype.xyzfilename)):
        xyzfilenamelist=poltype.xyzfilename[i]
        longestdimxlist=[]
        longestdimylist=[]
        longestdimzlist=[]
        for xyz in xyzfilenamelist:
            if poltype.complexation==True:

                if i==0: # complexation
                    x,y,z=FindDimensionsOfMoleculeTinkerXYZ(poltype,poltype.outputpath+xyz)
                elif i==1: # solvation

                    longestdim=FindDimensionsOfMoleculeTinker(poltype,poltype.outputpath+xyz)
                    x=longestdim
                    y=longestdim
                    z=longestdim

            elif poltype.complexation==False and (poltype.solvation==True or poltype.neatliquidsim==True):
                longestdim=FindDimensionsOfMoleculeTinker(poltype,poltype.outputpath+xyz)
                x=longestdim
                y=longestdim
                z=longestdim
            longestdimxlist.append(x)
            longestdimylist.append(y)
            longestdimzlist.append(z)
        longestdimx.append(longestdimxlist)
        longestdimy.append(longestdimylist)
        longestdimz.append(longestdimzlist)

    maxdirection=[]
    numberofsubboxes=[]
    subaaxislist=[]
    subbaxislist=[]
    subcaxislist=[]
    if poltype.aaxis==None and poltype.baxis==None and poltype.caxis==None:
        poltype.aaxislist=[]
        poltype.baxislist=[]
        poltype.caxislist=[]
        for i in range(len(longestdimx)):
            xdimlist=longestdimx[i]
            ydimlist=longestdimy[i]
            zdimlist=longestdimz[i]
            numberofsubboxes.append(len(xdimlist))
            maxxdim=max(xdimlist)
            maxydim=max(ydimlist)
            maxzdim=max(zdimlist)
            suballdim=[maxxdim,maxydim,maxzdim]
            alldim=suballdim.copy()
            maxdim=max(alldim)
            maxidx=alldim.index(maxdim)
            maxdirection.append(maxidx)
            maxxdim=alldim[0]
            maxydim=alldim[1]
            maxzdim=alldim[2]
            subxdim=suballdim[0]
            subydim=suballdim[1]
            subzdim=suballdim[2]

            buffervalue=poltype.bufferlen[i]
            xvalue=round(maxxdim+buffervalue,1)
            yvalue=round(maxydim+buffervalue,1)
            zvalue=round(maxzdim+buffervalue,1)
            subxvalue=round(subxdim+buffervalue,1)
            subyvalue=round(subydim+buffervalue,1)
            subzvalue=round(subzdim+buffervalue,1)
            alldim=[xvalue,yvalue,zvalue]
            alldim[maxidx]=alldim[maxidx]*len(xdimlist) # repeats for solvation box with multiple ligands
            xvalue=alldim[0]
            yvalue=alldim[1]
            zvalue=alldim[2]
            poltype.aaxislist.append(xvalue)
            poltype.baxislist.append(yvalue)
            poltype.caxislist.append(zvalue)
            subaaxislist.append(subxvalue)
            subbaxislist.append(subyvalue)
            subcaxislist.append(subzvalue)

    else:
        poltype.aaxislist=[poltype.aaxis,poltype.aaxis]
        poltype.baxislist=[poltype.baxis,poltype.baxis]
        poltype.caxislist=[poltype.caxis,poltype.caxis]
        subaaxislist=[poltype.aaxis,poltype.aaxis]
        subbaxislist=[poltype.baxis,poltype.baxis]
        subcaxislist=[poltype.caxis,poltype.caxis]
        suballdim=[poltype.aaxis,poltype.baxis,poltype.caxis]
        maxdim=max(suballdim)
        maxidx=suballdim.index(maxdim)
        maxdirection=[maxidx,maxidx]
        numberofsubboxes=[1,1]


    for i in range(len(poltype.tabledict)):
        poltype.tabledict[i]['Initial Box Length']=str(poltype.aaxislist[i])+','+str(poltype.baxislist[i])+','+str(poltype.caxislist[i])

    tables.WriteTableUpdateToLog(poltype)
    poltype.volume=[]
    subvolumelist=[]
    for i in range(len(poltype.configkeyfilename)):
        keylist=poltype.configkeyfilename[i]
        for j in range(len(keylist)):
            subaaxis=subaaxislist[i]
            subbaxis=subbaxislist[i]
            subcaxis=subcaxislist[i]

            keyname=keylist[j]
            aaxis=poltype.aaxislist[i]
            keymods.AddKeyWord(poltype,poltype.outputpath+keyname,'a-axis '+str(subaaxis)+'\n')
            baxis=poltype.baxislist[i]
            if baxis!=None:
                keymods.AddKeyWord(poltype,poltype.outputpath+keyname,'b-axis '+str(subbaxis)+'\n')
            caxis=poltype.caxislist[i]
            if caxis!=None:
                keymods.AddKeyWord(poltype,poltype.outputpath+keyname,'c-axis '+str(subcaxis)+'\n')
            volume=aaxis*baxis*caxis
            poltype.volume.append(volume)  
            subvolume=subaaxis*subbaxis*subcaxis
            subvolumelist.append(subvolume)

    return maxdirection,numberofsubboxes,subaaxislist,subbaxislist,subcaxislist,subvolumelist
 
def ComputePhysiologicalIonNumber(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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


def ComputeWaterNumber(poltype,subvolumelist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.waternum=[int(round(.0334*i)) for i in poltype.volume] # number of waters per cubic angstroms
    poltype.WriteToLog('Water Number '+str(poltype.waternum))
    subwaternum=[int(round(.0334*i)) for i in subvolumelist]
    return subwaternum

def RemoveBoxSizeTerms(poltype,keyfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keymods.RemoveKeyWords(poltype,keyfile,['axis'])


def FindDimensionsOfMoleculeTinker(poltype,structurefilepath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    veclist=[]
    temp=open(structurefilepath,'r')
    results=temp.readlines()
    temp.close()
    atomindexlist=[]
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            vec=np.array([float(linesplit[2]),float(linesplit[3]),float(linesplit[4])])
            veclist.append(vec)
            atomindex=int(linesplit[0])
            atomindexlist.append(atomindex)
    atompairs=list(itertools.combinations(atomindexlist, 2))
    disttodiffvec={}
    disttopairs={} 
    for pairidx in range(len(atompairs)):
        atompair=atompairs[pairidx]
        vec1=veclist[atompair[0]-1]
        vec2=veclist[atompair[1]-1]
        diff=np.array(vec1)-np.array(vec2)
        dist=np.linalg.norm(diff)
        disttopairs[dist]=atompair
        disttodiffvec[dist]=diff
    distlist=list(disttodiffvec.keys())
    if len(distlist)!=0:
        maxdist=np.amax(np.array(distlist))
        diffvec=disttodiffvec[maxdist]
        pair=disttopairs[maxdist]
    else:
        maxdist=0
    return maxdist

def FindDimensionsOfMoleculeTinkerXYZ(poltype,structurefilepath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    veclist=[]
    temp=open(structurefilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            vec=np.array([float(linesplit[2]),float(linesplit[3]),float(linesplit[4])])
            typenumber=int(linesplit[5])
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
    maxxatomindex=xdistlist.index(maxx)+1
    minxatomindex=xdistlist.index(minx)+1
    x=np.abs(maxx-minx) 
    maxy=max(ydistlist)
    miny=min(ydistlist)
    maxyatomindex=ydistlist.index(maxy)+1
    minyatomindex=ydistlist.index(miny)+1
    y=np.abs(maxy-miny) 
    maxz=max(zdistlist)
    minz=min(zdistlist)
    maxzatomindex=zdistlist.index(maxz)+1
    minzatomindex=zdistlist.index(minz)+1
    z=np.abs(maxz-minz) 

    return x,y,z





def TotalAtomNumber(poltype,xyzfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomnum=[]
    for subls in xyzfilename:
        tot=0
        for xyz in subls:
            atoms=FindNumberTinkerXYZAtoms(poltype,xyz)
            tot+=atoms
        atomnum.append(tot)
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.outputpath+xyzfilename,'r')
    results=temp.readlines()
    temp.close()
    atomnum=int(results[0].replace('\n','').rstrip().lstrip())
    return atomnum

def GrabTypeNumber(poltype,typedescrip):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if typedescrip in line:
            linesplit=line.split()
            typenum=linesplit[1]
    return typenum
   
def CreateWaterXYZ(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.outputpath+'water.xyz','w')
    temp.write('3'+'\n')
    temp.write('1  O      8.019814    5.892935    0.449481   '+str(poltype.waterOtypenum)+'     2     3'+'\n')
    temp.write('2  H      7.906681    4.942582    0.463215   '+str(poltype.waterHtypenum)+'     1'+'\n')
    temp.write('3  H      8.149325    6.169963   -0.457590   '+str(poltype.waterHtypenum)+'     1'+'\n')
    temp.close() 

      

def CreateSolventBox(poltype,aaxis,baxis,caxis,waternum,filename,key):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newfilename=filename+'_2'
    if not os.path.isfile(newfilename):
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
    return newfilename



def SoakMoleculeInSolventBox(poltype,xyzfilename,keyfilename,solventbox):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    cmdstr=poltype.xyzeditpath+' '+xyzfilename+' '+'-k'+' '+keyfilename+' '+'24'+' '+solventbox
    submit.call_subsystem(poltype,cmdstr,wait=True)    
    newname=NewTinkerFileName(poltype,xyzfilename)
    return newname

def AddIonToSolventBox(poltype,solutexyzfilename,keyfilename,tinktype,ionnum,soluteindices,writesolute=True):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    inputfile='xyzedit'+'.in'
    temp=open(poltype.outputpath+inputfile,'w')
    if writesolute==True:
        if len(soluteindices)==1:
            indices=soluteindices[0]
            temp.write('-1'+' '+str(indices[-1])+'\n')
        else:
            string=''
            for indices in soluteindices:
                for index in indices:
                    string+=str(index)+','
            string=string[:-1]
            temp.write(string+'\n')
    else:
        temp.write('0'+'\n')
    string=''
    string+=str(tinktype)+' '+str(ionnum)+' '
    string+='\n'
    temp.write(string)
    temp.write('\n')
    temp.close()
    cmdstr=poltype.xyzeditpath+' '+solutexyzfilename+' '+'-k'+' '+keyfilename+' '+'25'+' '+' < '+inputfile
    submit.call_subsystem(poltype,cmdstr,wait=True)    
    newsolutexyzfilename=NewTinkerFileName(poltype,solutexyzfilename)
    return newsolutexyzfilename


def AddIonsToSolventBox(poltype,solutexyzfilename,keyfilename,iontypetoionnumberneut,iontypetoionnumberphysio,soluteindices,writesolute=True):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    iontypetototalionnum=CombineDictionaries(poltype,[iontypetoionnumberneut,iontypetoionnumberphysio])
    for tinktype in iontypetototalionnum.keys():
        ionnum=iontypetototalionnum[tinktype]
        solutexyzfilename=AddIonToSolventBox(poltype,solutexyzfilename,keyfilename,tinktype,ionnum,soluteindices,writesolute)
    return solutexyzfilename



def CombineDictionaries(poltype,diclist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    combined={}
    for dic in diclist:
        for key,value in dic.items():
            if key not in combined.keys():
                combined[key]=0
            combined[key]+=value

    return combined


def RemoveTempFiles(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    files=os.listdir()
    for f in files:
        if '.' in f:
            fsplit=f.split('.')
            end=fsplit[1]
            if 'xyz_' in end or 'in' in end and 'ini' not in end:
                os.remove(f) 


def GrabIonIndexes(poltype,ionnumber,boxfilename,iontypenumber):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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




def WrapStrayAtomsBackInBox(poltype,newname,key):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.outputpath+'xyzedit.in','w')
    temp.write('18'+'\n')
    temp.write('\n')
    temp.close()
    cmdstr=poltype.xyzeditpath+' '+newname+' '+'-k'+' '+key+' <'+' '+'xyzedit.in'
    submit.call_subsystem(poltype,cmdstr,wait=True)   
    filename=NewTinkerFileName(poltype,newname)
    return filename



def ShiftSolventBoxCoordinates(poltype,newfilename,key,shiftvector):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.outputpath+'xyzedit.in','w')
    temp.write('12'+'\n')
    temp.write(str(shiftvector[0])+' '+str(shiftvector[1])+' '+str(shiftvector[2])+'\n')
    temp.write('\n')
    temp.close()
    cmdstr=poltype.xyzeditpath+' '+newfilename+' '+'-k'+' '+key+' <'+' '+'xyzedit.in'
    submit.call_subsystem(poltype,cmdstr,wait=True)   
    filename=NewTinkerFileName(poltype,newfilename)
    return filename


def NewTinkerFileName(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    time.sleep(1)
    prevcount=1
    numtofiles={}
    for f in os.listdir():
        if '.xyz_' in f:
            split=f.split('.xyz_')
            prefix=''.join(split[:-1])
            if prefix in filename:
                prevcount=int(split[-1])     
                numtofiles[prevcount]=f

    if len(numtofiles.keys())>0:
        maxcount=max(numtofiles.keys())
        maxfile=numtofiles[maxcount]
        outfilename=maxfile
    else:
        count=prevcount+1
        outfilename=filename+'_'+str(count)
    return outfilename




def AppendXYZs(poltype,xyztoappend,key):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    firstxyz=xyztoappend[0]
    for i in range(1,len(xyztoappend)):
        xyz=xyztoappend[i]
        firstxyz=AppendXYZ(poltype,firstxyz,xyz,key)

    return firstxyz


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
    filename=NewTinkerFileName(poltype,firstxyz)
    return filename


def CorrectVolume(poltype,xyz,axis):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempxyz=xyz.replace('.xyz','_TEMP.xyz')
    temp=open(xyz,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempxyz,'w')
    for line in results:
        linesplit=line.split()
        if '90.00' in line and len(linesplit)==6:
            linesplit[0]=str(axis[0])
            linesplit[1]=str(axis[1])
            linesplit[2]=str(axis[2])
            line=' '.join(linesplit)+'\n'
        temp.write(line)
    temp.close()
    os.rename(tempxyz,xyz)
    return xyz


def GrabIndicesFromType(poltype,xyz,types):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indices=[]
    statexyzatominfo,stateindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=parametercomparison.GrabXYZInfo(poltype,xyz)
    for stateindex,typeindex in stateindextotypeindex.items():
        if typeindex in types:
            indices.append(stateindex)

    return indices



def MakeLigandsConsecutive(poltype,finalxyz,allligandtypes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(finalxyz,'r')
    results=temp.readlines()
    temp.close()
    liglines=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            xyzatomnum=int(linesplit[0])
        else:
            if len(linesplit)>1 and '90.00' not in line:
                index=int(linesplit[0])
                typenum=int(linesplit[5])
                if typenum in allligandtypes:
                    liglines.append(line)
    tempname=finalxyz.replace('.xyz','_TEMP.xyz')
    temp=open(tempname,'w')
    first=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            xyzatomnum=int(linesplit[0])
            temp.write(line)
        else:
            if len(linesplit)>1 and '90.00' not in line:
                if first==False:
                    for ln in liglines:
                        temp.write(ln)
                    first=True
                if line not in liglines:
                    temp.write(line)
            else:
                temp.write(line)
    temp.close()
    os.rename(tempname,finalxyz)
    return finalxyz


def BoxSetupProtocol(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog('Computing volume ',prin=True)
    maxdirection,numberofsubboxes,subaaxislist,subbaxislist,subcaxislist,subvolumelist=ComputeBoxSize(poltype)
    ligandtypes=[]
    allligandtypes=[]
    for xyz in poltype.ligandxyzfilenamelist:
        statexyzatominfo,stateindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=parametercomparison.GrabXYZInfo(poltype,xyz)
        types=list(stateindextotypeindex.values())
        if xyz in poltype.annihilateligandxyzfilenamelist:
            ligandtypes.extend(types)
        allligandtypes.extend(types)
    if poltype.addsolvionwindows==True: # create box file for ions
        ionaxis=poltype.bufferlen[0]
        ionvolume=ionaxis**3
        ionwaternum=int(round(.0334*ionvolume)) 
        poltype.volume.append(ionvolume)

    if poltype.neatliquidsim==False:
        subwaternum=ComputeWaterNumber(poltype,subvolumelist)
        ComputePhysiologicalIonNumber(poltype)
        ComputeNeutralizingIonNumber(poltype,poltype.systemcharge)
        TotalAtomNumber(poltype,poltype.xyzfilename)
        CreateWaterXYZ(poltype)
    for i in range(len(poltype.boxfilename)): 
        boxxyzfilename=poltype.boxfilename[i]
        xyzfilenamelist=poltype.xyzfilename[i]
        key=poltype.configkeyfilename[i][0]
        aaxis=poltype.aaxislist[i]
        baxis=poltype.baxislist[i]
        caxis=poltype.caxislist[i]
        subaaxis=subaaxislist[i]
        subbaxis=subbaxislist[i]
        subcaxis=subcaxislist[i]
        if poltype.neatliquidsim==False:
            subwater=subwaternum[i]
        maxdir=maxdirection[i]
        numboxes=numberofsubboxes[i]
        if not os.path.isfile(poltype.outputpath+boxxyzfilename):
            if poltype.usepreequilibriatedbox==False:
                if poltype.neatliquidsim==False:
                    newfilename=CreateSolventBox(poltype,subaaxis,subbaxis,subcaxis,subwater,'water.xyz',key)
                else:
                    xyzfilename=xyzfilenamelist[0]
                    mass=GrabTotalMass(xyzfilename)
                    mass=mass*1.66054*10**(-27) # convert daltons to Kg
                    ls=[aaxis,baxis,caxis]
                    ls=[axis*10**-10 for axis in ls]
                    numbermolecules=int(poltype.density*(ls[0]*ls[1]*ls[2])/mass)
                    newfilename=CreateSolventBox(poltype,subaaxis,subbaxis,subcaxis,numbermolecules,xyzfilename,key)
                    os.rename(newfilename,boxxyzfilename)

            else:
                boxsize=[aaxis,baxis,caxis]
                poltype.TrimPreEquilibriatedBox(boxsize)
            if poltype.neatliquidsim==False:
                axis=[aaxis,baxis,caxis]
                subaxis=[subaaxis,subbaxis,subcaxis]
                xyztoappend=[]
                soluteindices=[]
                atmshift=0
                for idx in range(len(xyzfilenamelist)):
                    xyzfilename=xyzfilenamelist[idx]
                    soluteatomnum=FindNumberTinkerXYZAtoms(poltype,xyzfilename)
                    indices=list(range(1,soluteatomnum+1))
                    indices=[k+atmshift for k in indices]
                    soluteindices.append(indices)
                    maxlen=axis[maxdir]
                    shift=maxlen-((idx+1)/numboxes)*maxlen
                    shiftvector=axis.copy()
                    shiftvector[maxdir]=shift
                    for k in range(len(shiftvector)):
                        if k!=maxdir:
                            shiftvector[k]=0
                    newname=SoakMoleculeInSolventBox(poltype,xyzfilename,key,newfilename)
                    wrappedname=WrapStrayAtomsBackInBox(poltype,newname,key)
                    shiftedfilename=ShiftSolventBoxCoordinates(poltype,wrappedname,key,shiftvector)
                    newatomnum=FindNumberTinkerXYZAtoms(poltype,shiftedfilename)
                    atmshift=newatomnum
                    xyztoappend.append(shiftedfilename)
                finalxyz=AppendXYZs(poltype,xyztoappend,key)
                finalxyz=CorrectVolume(poltype,finalxyz,axis)
                if i==1 and poltype.binding==True:
                    finalxyz=MakeLigandsConsecutive(poltype,finalxyz,allligandtypes)
                keymods.RemoveKeyWords(poltype,key,['axis'])
                keymods.AddKeyWord(poltype,key,'a-axis '+str(aaxis)+'\n')
                keymods.AddKeyWord(poltype,key,'b-axis '+str(baxis)+'\n')
                keymods.AddKeyWord(poltype,key,'c-axis '+str(caxis)+'\n')

                newname=AddIonsToSolventBox(poltype,finalxyz,key,poltype.iontypetoionnumberneut[i],poltype.iontypetoionnumberphysio[i],soluteindices)

                shutil.copy(newname,boxxyzfilename)
                if poltype.binding==True:
                    alzout='checknetcharge.alz'
                    poltype.CheckNetChargeIsZero(boxxyzfilename,key,alzout)

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
        thexyz=poltype.boxfilename[index] 
        allligandindices=GrabIndicesFromType(poltype,thexyz,allligandtypes)
        ligandindices=GrabIndicesFromType(poltype,thexyz,ligandtypes)
        poltype.ligandindices[index]=ligandindices
        poltype.allligandindices[index]=allligandindices

    if poltype.complexation==True:
        thexyz=poltype.boxfilename[0]
        allligandindices=GrabIndicesFromType(poltype,thexyz,allligandtypes)
        ligandindices=GrabIndicesFromType(poltype,thexyz,ligandtypes)
        poltype.ligandindices[0]=ligandindices
        poltype.allligandindices[0]=allligandindices

    if poltype.addsolvionwindows==True:
        shutil.copy(poltype.configkeyfilename[0][0],poltype.ionkeyfilename)
        keymods.RemoveKeyWords(poltype,poltype.ionkeyfilename,['axis'])
        string='a-axis '+str(ionaxis)+'\n'
        keymods.AddKeyWord(poltype,poltype.ionkeyfilename,string)
        newfilename=CreateSolventBox(poltype,ionaxis,ionaxis,ionaxis,ionwaternum,'water.xyz',poltype.ionkeyfilename)
        filename='water.xyz'
        soluteindices=[[0]]
        newname=AddIonsToSolventBox(poltype,filename,poltype.ionkeyfilename,poltype.solviontocount,poltype.iontypetoionnumberphysio[-1],soluteindices,False)
        os.rename(newname,poltype.ionboxfilename)
    

    if poltype.boxonly==True:
        sys.exit() 
        
