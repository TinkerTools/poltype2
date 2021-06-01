import symmetry as symm
import optimization as opt
import electrostaticpotential as esp
import databaseparser
import torsionfit as torfit
import apicall as call
import os
import sys
import openbabel
import shutil
from socket import gethostname
import re
import time
import numpy
from rdkit.Chem.rdmolfiles import MolFromMol2File
from rdkit.Chem import rdMolTransforms as rdmt
from rdkit import Chem
from itertools import product,combinations
import shlex
from rdkit.Chem import rdmolfiles
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles

def __init__(poltype):
    PolarizableTyper.__init__(poltype)



def DefaultMaxRange(poltype,torsions):
    poltype.rotbndtomaxrange={}
    for torsionset in torsions:
        for torsion in torsionset:
            a,b,c,d=torsion[0:4]
            key=str(b)+' '+str(c)
            poltype.rotbndtomaxrange[key]=360


def RemoveCommentsFromKeyFile(poltype,keyfilename):
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    passedatoms=False
    tempname='temp.key'
    newtemp=open(tempname,'w')
    for line in results:
        if 'atom' in line:
            passedatoms=True 
        if passedatoms==True and '#' in line:
            continue
        else:
            newtemp.write(line)
    newtemp.close()
        
    shutil.move(tempname, keyfilename)

def ExecuteOptJobs(poltype,listofstructurestorunQM,phaselist,optmol,torset,variabletorlist,torsionrestraint,mol,currentopt):
    jobtooutputlog={}
    listofjobs=[]
    outputlogs=[]
    initialstructures=[]
    scratchdir=None
    for i in range(len(listofstructurestorunQM)):
        torxyzfname=listofstructurestorunQM[i]
        phaseangles=phaselist[i]
        inputname,outputlog,cmdstr,scratchdir=GenerateTorsionOptInputFile(poltype,torxyzfname,torset,phaseangles,optmol,variabletorlist,mol,currentopt)
        finished,error=poltype.CheckNormalTermination(outputlog,errormessages=None,skiperrors=True)
        if finished==True and 'opt' in outputlog:
            opt.GrabFinalXYZStructure(poltype,outputlog,outputlog.replace('.log','.xyz'),mol)
        if finished==False:
            listofjobs.append(cmdstr)
            jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputlog


        outputlogs.append(outputlog)
        initialstructures.append(torxyzfname)   
    return outputlogs,listofjobs,scratchdir,jobtooutputlog,initialstructures

def CheckIfPsi4Log(poltype,outputlog):
    check=False
    temp=open(outputlog,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Psi4' in line:
            check=True
            break 
    return check 


def ExecuteSPJobs(poltype,optoutputlogs,phaselist,optmol,torset,variabletorlist,torsionrestraint,outputlogtophaseangles,mol):
    jobtooutputlog={}
    listofjobs=[]
    outputnames=[]
    for i in range(len(optoutputlogs)):
        phaseangles=phaselist[i]
        outputlog=optoutputlogs[i]
        if not poltype.use_gaus:
            prevstrctfname=outputlog.replace('.log','.xyz')
            inputname,outputname=CreatePsi4TorESPInputFile(poltype,prevstrctfname,optmol,torset,phaseangles,mol)
            cmdstr='psi4 '+inputname+' '+outputname
        else:
            inputname,outputname=GenerateTorsionSPInputFileGaus(poltype,torset,optmol,phaseangles,outputlog,mol)
            cmdstr = 'GAUSS_SCRDIR='+poltype.scrtmpdirgau+' '+poltype.gausexe+' '+inputname
        finished,error=poltype.CheckNormalTermination(outputname,errormessages=None,skiperrors=True)
        if finished==True and error==False:
            pass
        else:
            listofjobs.append(cmdstr)
            jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputname

        outputnames.append(outputname)
        outputlogtophaseangles[outputname]=phaseangles
    if not poltype.use_gaus:
        
        return outputnames,listofjobs,poltype.scrtmpdirpsi4,jobtooutputlog,outputlogtophaseangles
    else:
        return outputnames,listofjobs,poltype.scrtmpdirgau,jobtooutputlog,outputlogtophaseangles



def CreateGausTorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist,mol,currentopt):

    prefix='%s-opt-%s_' % (poltype.molecprefix,currentopt)
    postfix='.com' 
    toroptcomfname,angles=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
    strctfname = os.path.splitext(toroptcomfname)[0] + '.log'
    gen_torcomfile(poltype,toroptcomfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,optmol,torxyzfname,mol) # prevstruct is just used for Atomic Num and charge,torxyzfname is used for xyz coordinates
    # Append restraints to *opt*.com file
    tempname=toroptcomfname.replace('.com','_temp.com')
    temp=open(toroptcomfname,'r')
    results=temp.readlines()
    temp.close()
    tmpfh = open(tempname, "w")
    foundatomblock=False
    writeres=False
    for k in range(len(results)):
        line=results[k]
        linesplit=line.split() 
        if len(linesplit)==4 and foundatomblock==False and '#' not in line:
            foundatomblock=True
        if len(linesplit)!=4 and foundatomblock==True and writeres==False:
            writeres=True
            tmpfh.write('\n')
            restlist=[]
            # Fix all torsions around the rotatable bond b-c
            count=0
            for i in range(len(torset)):
                tor=torset[i]
                a,b,c,d=tor[:]
                indices=[a,b,c,d]
                allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,indices,mol)
                allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,indices,mol)

                aatom=mol.GetAtom(a)
                datom=mol.GetAtom(d)
                aatomicnum=aatom.GetAtomicNum()
                datomicnum=datom.GetAtomicNum()
                if (aatomicnum==1 or datomicnum==1) and (allhydtors==False and allhydtorsoneside==False):
                    continue
                if (allhydtorsoneside==True or allhydtors==True) and (aatomicnum==1 or datomicnum==1):
                    if count>=1:
                        continue
                    count+=1

                babelindices=[b,c]
                string=' '.join([str(b),str(c)])
                if string in poltype.rotbndlist.keys(): 
                    tors=poltype.rotbndlist[string]
                    count=0
                    for resttors in tors:
                        rta,rtb,rtc,rtd = resttors
                        indices=[rta,rtb,rtc,rtd]
                        allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,indices,mol)
                        allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,indices,mol)
                        aatom=mol.GetAtom(rta)
                        datom=mol.GetAtom(rtd)
                        rtaatomicnum=aatom.GetAtomicNum()
                        rtdatomicnum=datom.GetAtomicNum()
                        if (rtaatomicnum==1 or rtdatomicnum==1) and (allhydtors==False and allhydtorsoneside==False):
                            continue
                        if (allhydtorsoneside==True or allhydtors==True) and (rtaatomicnum==1 or rtdatomicnum==1):

                      
                            if count>=1:
                                continue
                            count+=1

                        rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
                        if resttors not in variabletorlist and resttors not in restlist:
                            tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))
                            restlist.append(resttors)
            # Leave all torsions around other rotatable bonds fixed
            for rotkey,torsions in poltype.rotbndlist.items():
                for resttors in torsions:
                    rta,rtb,rtc,rtd = resttors
                    rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
                    if resttors not in restlist and resttors not in variabletorlist:
                        tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))
                        restlist.append(resttors)
            
            tmpfh.write("\n")
        else:
            tmpfh.write(line)
    tmpfh.write("\n")
    tmpfh.close()
    os.remove(toroptcomfname)
    shutil.copy(tempname,toroptcomfname)
    outputname= os.path.splitext(toroptcomfname)[0] + '.log'
    return toroptcomfname,outputname
    

def GenerateTorsionOptInputFile(poltype,torxyzfname,torset,phaseangles,optmol,variabletorlist,mol,currentopt):
    if  poltype.use_gaus==False and poltype.use_gausoptonly==False:
        inputname,outputname=CreatePsi4TorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist,mol,currentopt)
        cmdstr='psi4 '+inputname+' '+outputname
    else:
        inputname,outputname=CreateGausTorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist,mol,currentopt)
        cmdstr='GAUSS_SCRDIR='+poltype.scrtmpdirgau+' '+poltype.gausexe+' '+inputname
    if poltype.use_gaus==False and poltype.use_gausoptonly==False:
        return inputname,outputname,cmdstr,poltype.scrtmpdirpsi4

    else:
        return inputname,outputname,cmdstr,poltype.scrtmpdirgau
        

def GenerateTorsionSPInputFileGaus(poltype,torset,optmol,phaseangles,prevstrctfname,mol):
    prevstruct = opt.load_structfile(poltype,prevstrctfname)
    prefix='%s-sp-' % (poltype.molecprefix)
    postfix='.com' 
    torspcomfname,angles=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)

    torsplogfname = os.path.splitext(torspcomfname)[0] + '.log'
    gen_torcomfile(poltype,torspcomfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,prevstruct,None,mol)
    outputname=torspcomfname.replace('.com','.log')
    return torspcomfname,outputname


def tinker_minimize_angles(poltype,torset,optmol,variabletorlist,phaselist,prevstrctfname,torsionrestraint,bondtopology):
    tinkerstructnamelist=[]
    # load prevstruct
    firsttor=torset[0]
    a,b,c,d=firsttor[0:4]
    tor=tuple([a,b,c,d])
    newtorxyzfname=None  
    # create xyz and key and write restraint then minimize, getting .xyz_2
    count=-1
    center=prevstrctfname
    prevphaseanglelist=None
    energyarray=numpy.empty(phaselist.shape[0])
    failedgridpoints=[]

        
    for rowindex in range(len(phaselist)): # we need to send back minimized structure in XYZ (not tinker) format to load for next tinker minimization,but append the xyz_2 tinker XYZ file so that com file can be generated from that 
        phaseanglelist=phaselist[rowindex]
        count+=1
        if len(phaseanglelist)==2:
            a=phaseanglelist[1]
            idx=1
        else:
            a=phaseanglelist[0]
            idx=0
        needtogobacktocenter=False
        try:
            preva=prevphaseanglelist[idx]
            if (a>0 and preva<0) or (a<0 and preva>0):
                needtogobacktocenter=True
        except: # first time no prevphaseanglelist
            pass
        if newtorxyzfname==None or needtogobacktocenter==True:
            prevstruct = opt.load_structfile(poltype,center)
        else:
            prevstruct =opt.load_structfile(poltype,cartxyz)

        prevstruct = opt.PruneBonds(poltype,prevstruct,bondtopology) # sometimes extra bonds are made when atoms get too close during minimization
        prevstruct=opt.rebuild_bonds(poltype,prevstruct,optmol)
        
        torxyzfname,tmpkeyfname,torminlogfname=tinker_minimize_filenameprep(poltype,torset,optmol,variabletorlist,phaseanglelist,torsionrestraint,prevstruct,'_preQMOPTprefit',poltype.key4fname,'../')
        prevstrctfname,torxyzfname,newtorxyzfname,keyfname=tinker_minimize(poltype,torset,optmol,variabletorlist,phaseanglelist,torsionrestraint,prevstruct,'_preQMOPTprefit',poltype.key4fname,'../',torxyzfname,tmpkeyfname,torminlogfname)
        toralzfname = os.path.splitext(torxyzfname)[0] + '.alz'
        term=AnalyzeTerm(poltype,toralzfname)
        if term==False:
            tinker_analyze(poltype,newtorxyzfname,keyfname,toralzfname)
        if term==False:
            failedgridpoints.append(phaseanglelist)
        else:
            tot_energy,tor_energy=torfit.GrabTinkerEnergy(poltype,toralzfname)
            energyarray[rowindex]=tot_energy

        cartxyz=ConvertTinktoXYZ(poltype,newtorxyzfname,newtorxyzfname.replace('.xyz_2','_cart.xyz'))
        tinkerstructnamelist.append(newtorxyzfname)
        try:
            if len(phaseanglelist)==2: # then check for new center make sure hopping around 2D grid correctly (clockwise then counterclockwise)
                a=phaseanglelist[0]
                aprev=prevphaseanglelist[0]
                if a!=aprev:
                    center=cartxyz    
        except: # no prevphaseanglelist
            pass

        prevphaseanglelist=phaseanglelist
    return tinkerstructnamelist,failedgridpoints,energyarray



def tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,prevstrctfname,torsionrestraint,designatexyz,keybase,keybasepath,bondtopology):
    prevstruct = opt.load_structfile(poltype,prevstrctfname) 
    prevstruct=opt.rebuild_bonds(poltype,prevstruct,optmol)
    prevstruct = opt.PruneBonds(poltype,prevstruct,bondtopology)
    torxyzfname,tmpkeyfname,torminlogfname=tinker_minimize_filenameprep(poltype,torset,optmol,variabletorlist,phaseangles,torsionrestraint,prevstruct,designatexyz,keybase,keybasepath)
    cartxyz,torxyzfname,newtorxyzfname,keyfname=tinker_minimize(poltype,torset,optmol,variabletorlist,phaseangles,torsionrestraint,prevstruct,designatexyz,keybase,keybasepath,torxyzfname,tmpkeyfname,torminlogfname)
    toralzfname = os.path.splitext(torxyzfname)[0] + '.alz'
    term=AnalyzeTerm(poltype,toralzfname)
    if term==False:
        tinker_analyze(poltype,newtorxyzfname,keyfname,toralzfname)

    return cartxyz,newtorxyzfname

def AnalyzeTerm(poltype,filename):
    term=False
    if os.path.isfile(filename):
        temp=open(filename,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if 'Total Potential Energy :' in line:
                term=True
    return term

def tinker_analyze(poltype,torxyzfname,keyfname,toralzfname):
    alzcmdstr=poltype.analyzeexe+' -k '+keyfname+' '+torxyzfname+' ed > %s' % toralzfname
    poltype.call_subsystem(alzcmdstr,True)


def tinker_minimize_filenameprep(poltype,torset,optmol,variabletorlist,phaseanglelist,torsionrestraint,prevstruct,designatexyz,keybase,keybasepath):
    prefix='%s-opt-' % (poltype.molecprefix)
    postfix='%s.xyz' % (designatexyz)
    torxyzfname,angles=GenerateFilename(poltype,torset,phaseanglelist,prefix,postfix,optmol)

    prefix='tmp-'
    postfix='%s.key' % (designatexyz)
    tmpkeyfname,angles=GenerateFilename(poltype,torset,phaseanglelist,prefix,postfix,optmol)

    torminlogfname=torxyzfname.replace('.xyz','.out')
    return torxyzfname,tmpkeyfname,torminlogfname 

       
def tinker_minimize(poltype,torset,optmol,variabletorlist,phaseanglelist,torsionrestraint,prevstruct,designatexyz,keybase,keybasepath,torxyzfname,tmpkeyfname,torminlogfname):
    save_structfile(poltype,prevstruct,torxyzfname)
    shutil.copy(keybasepath+keybase, tmpkeyfname)
    tmpkeyfh = open(tmpkeyfname,'a')
    restlist=[]
    count=0
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[0:4]
        indices=[a,b,c,d]
        torang = optmol.GetTorsion(a,b,c,d)
        phaseangle=phaseanglelist[i]
        allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,indices,optmol)
        allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,indices,optmol)

        aatom=optmol.GetAtom(a)
        datom=optmol.GetAtom(d)
        aatomicnum=aatom.GetAtomicNum()
        datomicnum=datom.GetAtomicNum()
        if (aatomicnum==1 or datomicnum==1) and (allhydtors==False and allhydtorsoneside==False):
            continue
        if (allhydtorsoneside==True or allhydtors==True) and (aatomicnum==1 or datomicnum==1):

       
            if count>=1:
                continue
            count+=1

        tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (a,b,c,d,torsionrestraint,round((torang+phaseangle)%360),round((torang+phaseangle)%360)))
        for key in poltype.rotbndlist.keys():
            torlist=poltype.rotbndlist[key]
            count=0
            for res in torlist:
                resa,resb,resc,resd = res[0:4]
                indices=[resa,resb,resc,resd]
                totalatoms=optmol.NumAtoms()
                allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,indices,optmol)
                allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,indices,optmol)

                aatom=optmol.GetAtom(resa)
                datom=optmol.GetAtom(resd)
                aatomicnum=aatom.GetAtomicNum()
                datomicnum=datom.GetAtomicNum()
                if (aatomicnum==1 or datomicnum==1) and (allhydtors==False and allhydtorsoneside==False):
                    continue
                if (allhydtorsoneside==True or allhydtors==True) and (aatomicnum==1 or datomicnum==1):

               
                    if count>=1:
                        continue
                    count+=1


                if res not in restlist:
                    if (resa,resb,resc,resd) not in torset and (resd,resc,resb,resa) not in torset and (resa,resb,resc,resd) not in variabletorlist and (resd,resc,resb,resa) not in variabletorlist:
                        if (b==resb and c==resc) or (b==resc and c==resb):
                            secondang = optmol.GetTorsion(resa,resb,resc,resd)
                            tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (resa,resb,resc,resd,torsionrestraint,round((secondang+phaseangle)%360),round((secondang+phaseangle)%360)))
                        else:
                            tmpkeyfh.write('restrain-torsion %d %d %d %d %f\n' % (resa,resb,resc,resd,torsionrestraint))
                        restlist.append(res)
    tmpkeyfh.close()
    mincmdstr = poltype.minimizeexe+' '+torxyzfname+' -k '+tmpkeyfname+' 0.1'+' '+'>'+torminlogfname
    term,error=poltype.CheckNormalTermination(torminlogfname)
    if term==True and error==False:
        pass
    else:
        poltype.call_subsystem(mincmdstr,True)
    filename=torxyzfname+'_2'
    newfilename=filename.replace('.xyz_2','_xyzformat.xyz')
    cartxyz=ConvertTinktoXYZ(poltype,torxyzfname+'_2',newfilename)
    return cartxyz,torxyzfname,torxyzfname+'_2',tmpkeyfname

def GenerateFilename(poltype,torset,phaseangles,prefix,postfix,mol):
    oldtorset=torset[:]
    if torset in poltype.torsettofilenametorset.keys():
        tortorindex=poltype.torsettotortorindex[torset]
        phaseangle=phaseangles[tortorindex]
        phaseindex=phaseangles.index(phaseangle)
        torset=poltype.torsettofilenametorset[torset]
    angles=[]
    filename=prefix
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[:]
        phaseangle=phaseangles[i]
        torang = mol.GetTorsion(a,b,c,d)
        angle=round((torang+phaseangle)%360)
        if oldtorset in poltype.torsettofilenametorset.keys():
            if phaseindex==i:
                angles.append(angle)

        else:
            angles.append(angle)
        filename+='%d-%d-%03d'%(b,c,angle)
        filename+='_'
    filename=filename[:-1]
    filename+=postfix
    return filename,angles


def GeneratePhaseLists(poltype,torset,optmol):
    phaselists=[]
    currentanglelist=[]
    for tor in torset:
        a,b,c,d = tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        currentanglelist.append(round(torang,1))
        key=str(b)+' '+str(c)
        anginc=poltype.rotbndtoanginc[key]
        maxrange=poltype.rotbndtomaxrange[key]
        phaselist=range(0,maxrange,anginc)
        halfindex=int(len(phaselist)/2)
        clock=list(phaselist[:halfindex+1])
        counterclock=[-1*i for i in phaselist[1:halfindex+1]]
        full=clock+counterclock
        checkfull=[]
        for i in range(len(full)):
            value=full[i]
            if value<0:
                value=value+360
            checkfull.append(value)
        check=[]
        removeindices=[]
        for i in range(len(checkfull)):
            value=checkfull[i]
            if value in check:
                removeindices.append(i)
            else:
                check.append(value)
        for index in removeindices: # should only be 1
            del full[index]
        phaselists.append(full)
    currentanglelist=numpy.array(currentanglelist)
    return phaselists,currentanglelist


def ConvertCartesianXYZToTinkerXYZ(poltype,cartxyz,tinkerxyz):
    temp=open('../'+poltype.xyzfname,'r')
    results=temp.readlines()
    temp.close() 
    tinkerxyzstuff=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            totalatoms=linesplit[0]
        if len(linesplit)>=6:
            index=linesplit[0]
            element=linesplit[1] 
            tinkertype=linesplit[5]
            if len(linesplit)>6:
                connections=linesplit[6:] 
            else:
                connections=[]
            tinkerxyzstuff.append([index,element,tinkertype,connections])
    cartxyzstuff=[]
    temp=open(cartxyz,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)==4:
            x=linesplit[1]
            y=linesplit[2]
            z=linesplit[3]
            cartxyzstuff.append([x,y,z]) 
    temp=open(tinkerxyz,'w')
    temp.write(totalatoms+'\n') 
    for i in range(len(tinkerxyzstuff)):
        tink=tinkerxyzstuff[i]
        coords=cartxyzstuff[i]
        index=tink[0]
        element=tink[1]
        tinkertype=tink[2]
        connections=tink[3]
        x=coords[0]
        y=coords[1]
        z=coords[2]
        line=index+' '+element+' '+x+' '+y+' '+z+' '+tinkertype+' '
        for item in connections:
            line+=item+' '
        line+='\n'
        temp.write(line)
    temp.close()


def gen_torsion(poltype,optmol,torsionrestraint,mol):
    """
    Intent: For each rotatable bond, rotate the torsion about that bond about 
    30 degree intervals. At each interval use Gaussian SP to find the energy of the molecule at
    that dihedral angle. Create an energy profile for each rotatable bond: 
    "QM energy vs. dihedral angle" 
    Input:
        mol: OBMol object
    Output:
    Referenced By: main
    Description:
    1. Create and change to directory 'qm-torsion'
    2. For each torsion in torlist (essentially, for each rotatable bond)
        a. Rotate the torsion value by interval of 30 from 30 to 180, then -30 to -210
        b. Find energy using Gaussian SP
    """
    if not os.path.isdir('qm-torsion'):
        os.mkdir('qm-torsion')
    
    os.chdir('qm-torsion')
    files=os.listdir(os.getcwd())
     
    poltype.optoutputtotorsioninfo={}

    fullrange=[]
    optnumtotorsettofulloutputlogs={}
    optnumtotorsettofulloutputlogs['1']={}
    optnumtotorsettofulloutputlogs['2']={}
    optnumtofulllistofjobs={}
    optnumtofulllistofjobs['1']=[]
    optnumtofulllistofjobs['2']=[]

    fulljobtolog={}
    fulljobtooutputlog={}
    poltype.torsettophaselist={}
    torsettooptoutputlogs={}
    bondtopology=GenerateBondTopology(poltype,optmol)
    poltype.idealangletensor={}
    poltype.tensorphases={}
    torsettooutputlogtoinitialstructure={}
    torsettoinistructophaselist={}
    
    for torset in poltype.torlist:
        phaseangles=[0]*len(torset)
        if poltype.use_gaus==False and poltype.use_gausoptonly==False:
            prefix='%s-opt-' % (poltype.molecprefix)
            postfix='-opt.xyz' 
            prevstrctfname,angles=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
            tempdir=os.getcwd()
            os.chdir('..')
            shutil.copy(poltype.logoptfname.replace('.log','.xyz'),os.path.join(tempdir,prevstrctfname))
            os.chdir(tempdir)
        else:
            prefix='%s-opt-' % (poltype.molecprefix)
            postfix='.log' 
            prevstrctfname,angles=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
            # copy *-opt.log found early by Gaussian to 'prevstrctfname'
            tempdir=os.getcwd()
            os.chdir('..')
            shutil.copy(poltype.logoptfname,os.path.join(tempdir,prevstrctfname))
            os.chdir(tempdir)


        torsettooutputlogtoinitialstructure[tuple(torset)]={}
        optnumtotorsettofulloutputlogs['1'][tuple(torset)]=[]
        optnumtotorsettofulloutputlogs['2'][tuple(torset)]=[]
        variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
        phaselists,currentanglelist=GeneratePhaseLists(poltype,torset,optmol) 
        flatphaselist=numpy.array(list(product(*phaselists)))
        poltype.tensorphases[tuple(torset)]=flatphaselist
        flatphaselist=FlattenArray(poltype,flatphaselist)
        idealangletensor=flatphaselist+currentanglelist
        poltype.idealangletensor[tuple(torset)]=idealangletensor
        minstrctfname = prevstrctfname
        prevstrctfname = minstrctfname
        listoftinkertorstructures,failedgridpoints,energyarray=tinker_minimize_angles(poltype,torset,optmol,variabletorlist,flatphaselist,prevstrctfname,torsionrestraint,bondtopology)
        if len(torset)==2 and torset not in poltype.nonaroringtorsets:
            indicestokeep,firsttorindices,secondtorindices=Determine1DTorsionSlicesOnTorTorSurface(poltype,energyarray,flatphaselist,torset)
            flatphaselist,listoftinkertorstructures=RemoveTorTorPoints(poltype,energyarray,flatphaselist,indicestokeep,listoftinkertorstructures)
            for toridx in range(len(torset)):
                tor=torset[toridx]
                if toridx==0:
                    poltype.torsettophaselist[tuple([tor])]=firsttorindices
                elif toridx==1:
                    poltype.torsettophaselist[tuple([tor])]=secondtorindices

        else:
            gridspacing=int(len(flatphaselist)/poltype.defaultmaxtorsiongridpoints)
            if gridspacing>1: 
                locstoremove=[]
                for i in range(1,len(flatphaselist)+1):
                    phases=flatphaselist[i-1]
                    if i % gridspacing == 0:
                        pass
                    else:
                        indexes=numpy.where((flatphaselist == phases).all(axis=1))[0][0]
                        locstoremove.append(indexes)
                flatphaselist = numpy.delete(flatphaselist,locstoremove , axis=0)
                listoftinkertorstructures=numpy.delete(listoftinkertorstructures,locstoremove,axis=0)
        poltype.torsettophaselist[tuple(torset)]=flatphaselist

        poltype.SanitizeAllQMMethods()
        
        outputlogs,listofjobs,scratchdir,jobtooutputlog,initialstructures=ExecuteOptJobs(poltype,listoftinkertorstructures,flatphaselist,optmol,torset,variabletorlist,torsionrestraint,mol,'1')
        if poltype.tordebugmode==True:
            sys.exit() 

        torsettooptoutputlogs[tuple(torset)]=outputlogs
        dictionary = dict(zip(outputlogs,initialstructures))  
        torsettooutputlogtoinitialstructure[tuple(torset)].update(dictionary)
        inistructophaselist = dict(zip(initialstructures,flatphaselist))  
        torsettoinistructophaselist[tuple(torset)]=inistructophaselist
        lognames=[]
        for job in listofjobs:
            log=jobtooutputlog[job]
            lognames.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+r'/'+poltype.logfname)
        jobtolog=dict(zip(listofjobs, lognames)) 
        optnumtotorsettofulloutputlogs['1'][tuple(torset)].extend(outputlogs)
        optnumtofulllistofjobs['1'].extend(listofjobs)
        fulljobtolog.update(jobtolog) 
        fulljobtooutputlog.update(jobtooutputlog)

    jobtologlistfilenameprefix=os.getcwd()+r'/'+'QMOptJobToLog'+'_1'+'_'+poltype.molecprefix
    if poltype.externalapi!=None:
        if len(optnumtofulllistofjobs['1'])!=0:
            call.CallExternalAPI(poltype,fulljobtolog,jobtologlistfilenameprefix,scratchdir)
        finishedjobs,errorjobs=poltype.WaitForTermination(fulljobtooutputlog)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(fulljobtooutputlog,True)
    for torset in poltype.torlist:
        outputlogtoinitialstructure=torsettooutputlogtoinitialstructure[tuple(torset)]
        tryoutputlogs=optnumtotorsettofulloutputlogs['1'][tuple(torset)]
        inistructophaselist=torsettoinistructophaselist[tuple(torset)]
        for outputlog in tryoutputlogs:
            initialtinkerstructure=outputlogtoinitialstructure[outputlog]
            phaseangles=inistructophaselist[initialtinkerstructure]
            finished,error=poltype.CheckNormalTermination(outputlog)
            cartxyz=outputlog.replace('.log','.xyz')
            opt.GrabFinalXYZStructure(poltype,outputlog,cartxyz,mol)
            tinkerxyz=outputlog.replace('.log','_tinker.xyz')
            ConvertCartesianXYZToTinkerXYZ(poltype,cartxyz,tinkerxyz) 
            if finished==False and 'opt' in outputlog: # if there is an error, then finished=False
                bondtoposame=CheckBondTopology(poltype,outputlog,initialtinkerstructure)
                maxiter=4
                attempts=0
                while bondtoposame==False:
                    if attempts>=maxiter or finished==True:
                        break
                    inputname,outputlog,cmdstr,scratchdir=GenerateTorsionOptInputFile(poltype,tinkerxyz,torset,phaseangles,optmol,variabletorlist,mol,'1')
                    jobtolog={cmdstr:outputlog}
                    otherfinishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtolog,True)
                    opt.GrabFinalXYZStructure(poltype,outputlog,cartxyz,mol)
                    ConvertCartesianXYZToTinkerXYZ(poltype,cartxyz,tinkerxyz) 
                    bondtoposame=CheckBondTopology(poltype,outputlog,initialtinkerstructure)
                    finished,error=poltype.CheckNormalTermination(outputlog)
                    attempts+=1
                if outputlog not in finishedjobs: # assume if not converged in 5 cycles (or 20 etc) then bond lengths are okay, if not torsionfit will catch later
                    finishedjobs.append(outputlog)

            else:
                if finished==True and outputlog not in finishedjobs:
                    finishedjobs.append(outputlog)

    firstfinishedjobs=finishedjobs[:]
    for job in firstfinishedjobs:
        if job not in finishedjobs:
            finishedjobs.append(job)
    fulllistofjobs=[]
    fulljobtolog={}
    fulljobtooutputlog={}
    torsettospoutputlogs={}
    outputlogtophaseangles={}
    for torset in poltype.torlist:
        variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        outputlogs=torsettooptoutputlogs[tuple(torset)]
        finishedoutputlogs=[]
        finishedflatphaselist=[]
        for i in range(len(outputlogs)):
            outputlog=outputlogs[i]
            phaselist=flatphaselist[i]
            if outputlog in finishedjobs and outputlog not in errorjobs:
                finishedoutputlogs.append(outputlog)
                finishedflatphaselist.append(phaselist)
        outputlogs,listofjobs,scratchdir,jobtooutputlog,outputlogtophaseangles=ExecuteSPJobs(poltype,finishedoutputlogs,finishedflatphaselist,optmol,torset,variabletorlist,torsionrestraint,outputlogtophaseangles,mol)
        lognames=[]
        torsettospoutputlogs[tuple(torset)]=outputlogs
        for job in listofjobs:
           log=jobtooutputlog[job]
           lognames.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+r'/'+poltype.logfname)
        jobtolog=dict(zip(listofjobs, lognames)) 
        fulljobtooutputlog.update(jobtooutputlog)
        fulllistofjobs.extend(listofjobs)
        fulljobtolog.update(jobtolog)

    for torset in poltype.torlist:
        variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        outputlogs=torsettospoutputlogs[tuple(torset)]
        optoutputlogs=torsettooptoutputlogs[tuple(torset)]
        for i in range(len(outputlogs)):
            outputlog=outputlogs[i]
            optoutputlog=optoutputlogs[i]
            if optoutputlog in finishedjobs and outputlog not in errorjobs:
                phaseangles=outputlogtophaseangles[outputlog]
                poltype.optoutputtotorsioninfo[outputlog]= [torset,optmol,variabletorlist,phaseangles,bondtopology,optoutputlog]
    jobtologlistfilenameprefix=os.getcwd()+r'/'+'QMSPJobToLog'+'_'+poltype.molecprefix
    if poltype.externalapi!=None:
        if len(fulllistofjobs)!=0:
            call.CallExternalAPI(poltype,fulljobtolog,jobtologlistfilenameprefix,scratchdir)
        finshedjobs,errorjobs=poltype.WaitForTermination(fulljobtooutputlog)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(fulljobtooutputlog,True)
    os.chdir('..')



def CheckBondTopology(poltype,outputlog,initialtinkerstructure):
    bondtoposame=True
    if poltype.use_gaus==False and poltype.use_gausoptonly==False:
        fname=outputlog.replace('.log','.xyz')
    else:
        fname=outputlog
    optmol = opt.load_structfile(poltype,fname)
    cartxyz=ConvertTinktoXYZ(poltype,initialtinkerstructure,initialtinkerstructure.replace('.xyz','_cart.xyz'))
    inioptmol = opt.load_structfile(poltype,cartxyz)
    isnear=opt.CompareBondLengths(poltype,inioptmol,optmol,outputlog)
    if isnear==False:
        bondtoposame=False

    return bondtoposame

def GenerateBondTopology(poltype,optmol):
    bondtopology=[]
    iterbond=openbabel.OBMolBondIter(optmol) # iterator for all bond objects in the molecule
    for bond in iterbond:
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        aidx=a.GetIdx()
        bidx=b.GetIdx()
        bondtopology.append(set([aidx,bidx]))
    return bondtopology


def FindPartialDoubleBonds(poltype,rdkitmol):
    poltype.partialdoublebonds=[]
    amidesmarts='[NX3][CX3](=[OX1])'
    acidsmarts='[*][CX3](=O)[OX2]'
    smartslist=[amidesmarts,acidsmarts]
    for smartsidx in range(len(smartslist)):
        smarts=smartslist[smartsidx]
        sp = openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,smarts)
        diditmatch=sp.Match(poltype.mol)
        matches=sp.GetMapList()
        for match in matches:
            if smartsidx==0:
                babelbond=[match[0],match[1]]
            else:
                babelbond=[match[0],match[1]]
            if babelbond not in poltype.partialdoublebonds:
                poltype.partialdoublebonds.append(babelbond)

def RdkitIsInRing(poltype,atom):
    isinringofsize=False
    i=None
    for i in range(3,13+1):
        isinringofsize=atom.IsInRingSize(i)
        if isinringofsize==True:
            break
    return isinringofsize,i


def get_torlist(poltype,mol,missed_torsions):
    """
    Intent: Find unique rotatable bonds.
    Input:
        mol: An openbabel molecule structure
    Output:
        torlist: contains list of 4-ples around which torsion scans are done.
        rotbndlist: contains a hash (indexed by middle 2 atoms surrounding
            rotatable bond) of lists that contains all possible combinations
            around each rotatable bond.
    Referenced By: main
    Description:
    1. Iterate over bonds
        a. Check 'IsRotor()' (is the bond rotatable?)
        b. Find the atoms 1 and 4 (of the highest possible sym_class) of a possible torsion about atoms t2 and t3 of the rotatable bond (calls find_tor_restraint_idx)
        c. Check if this torsion is in user provided toromitlist
        d. Check if this torsion is found in the look up table
        e. If it neither c nor d are true, then append this torsion to 'rotbndlist' for future torsion scanning
        f. Find other possible torsions around the bond t2-t3 and repeat steps c through e
    """
    torlist = []
    hydtorsionlist=[]
    rotbndlist = {}
    iterbond = openbabel.OBMolBondIter(mol)
    nonaroringtorlist=[]
    for bond in iterbond:
        skiptorsion=True
        # is the bond rotatable
        t2 = bond.GetBeginAtom()
        t3 = bond.GetEndAtom()
        BO=bond.GetBondOrder()
        t2idx=t2.GetIdx()
        t3idx=t3.GetIdx()
        bnd=[t2idx,t3idx]
        t2val=t2.GetValence()
        t3val=t3.GetValence()
        if BO>1:
            continue
        ringbond=bond.IsInRing()
        if t2val<2 or t3val<2:
            continue 
        t1,t4 = find_tor_restraint_idx(poltype,mol,t2,t3)
        babelatoms=[t1,t2,t3,t4]
        aromatics=[i.IsAromatic() for i in babelatoms]
        hybs=[i.GetHyb() for i in babelatoms]
        if ringbond==True:
            if 2 in hybs or True in aromatics:
                continue
        
        sortedtor=torfit.sorttorsion(poltype,[poltype.idxtosymclass[t1.GetIdx()],poltype.idxtosymclass[t2.GetIdx()],poltype.idxtosymclass[t3.GetIdx()],poltype.idxtosymclass[t4.GetIdx()]])
        if(sortedtor in missed_torsions or sortedtor[::-1] in missed_torsions) and len(poltype.onlyrotbndslist)==0:
            skiptorsion = False
        onlyrot=False
        if [t2.GetIdx(),t3.GetIdx()] in poltype.onlyrotbndslist or [t3.GetIdx(),t2.GetIdx()] in poltype.onlyrotbndslist:
            skiptorsion = False
            onlyrot=True

        if poltype.rotalltors==True:
            skiptorsion=False
        if (bnd in poltype.partialdoublebonds or bnd[::-1] in poltype.partialdoublebonds) and poltype.rotalltors==False and onlyrot==False:
            skiptorsion=True
        babelindices=[t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx()]
        t1atomicnum=t1.GetAtomicNum()
        t4atomicnum=t4.GetAtomicNum()
        allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol)
        allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol)
        if (t1atomicnum==1 or t4atomicnum==1) and (allhydtors==False and allhydtorsoneside==False) and poltype.rotalltors==False:
            skiptorsion=True
        if (t1atomicnum==1 or t4atomicnum==1):
            hydtorsionlist.append(sortedtor)
        unq=get_uniq_rotbnd(poltype,t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx())
        if ringbond==True and poltype.allownonaromaticringscanning==False and poltype.refinenonaroringtors==False and len(poltype.onlyrotbndslist)==0 and poltype.dontfrag==False:
            nonaroringtorlist.append(unq)
            skiptorsion=False
        elif ringbond==True and poltype.refinenonaroringtors==True:
            skiptorsion=True
        rotbndkey = '%d %d' % (unq[1],unq[2])
        # store the torsion and temporary torsion value found by openbabel in torlist
        tor = mol.GetTorsion(t1,t2,t3,t4)
        if not skiptorsion:
            torlist.append(unq)
        # store torsion in rotbndlist
        rotbndlist[rotbndkey] = []
        rotbndlist[rotbndkey].append(unq)
        # write out rotatable bond to log
        #Find other possible torsions about this rotatable bond
        iteratomatom = openbabel.OBAtomAtomIter(bond.GetBeginAtom())
        for iaa in iteratomatom:
            iteratomatom2 = openbabel.OBAtomAtomIter(bond.GetEndAtom())
            for iaa2 in iteratomatom2:
                a = iaa.GetIdx()
                b = t2.GetIdx()
                c = t3.GetIdx()
                d = iaa2.GetIdx()
                if ((iaa.GetIdx() != t3.GetIdx() and iaa2.GetIdx() != t2.GetIdx()) and not (iaa.GetIdx() == t1.GetIdx() and iaa2.GetIdx() == t4.GetIdx())):
                    rotbndlist[rotbndkey].append(get_uniq_rotbnd(poltype,iaa.GetIdx(),t2.GetIdx(),t3.GetIdx(),iaa2.GetIdx()))
    return (torlist ,rotbndlist,hydtorsionlist,nonaroringtorlist)



def get_all_torsions(poltype,mol):
    poltype.alltorsionslist=[]
    iterbond = openbabel.OBMolBondIter(mol)
    for bond in iterbond:
        t2 = bond.GetBeginAtom()
        t3 = bond.GetEndAtom()
        t2idx=t2.GetIdx()
        t3idx=t3.GetIdx()
        t2val=t2.GetValence()
        t3val=t3.GetValence()
        if (t2val>=2 and t3val>=2):
            t1,t4 = find_tor_restraint_idx(poltype,mol,t2,t3)
            unq=get_uniq_rotbnd(poltype,t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx())
            poltype.alltorsionslist.append(unq)
            iteratomatom = openbabel.OBAtomAtomIter(bond.GetBeginAtom())
            for iaa in iteratomatom:
                iteratomatom2 = openbabel.OBAtomAtomIter(bond.GetEndAtom())
                for iaa2 in iteratomatom2:
                    a = iaa.GetIdx()
                    b = t2.GetIdx()
                    c = t3.GetIdx()
                    d = iaa2.GetIdx()
                    if ((iaa.GetIdx() != t3.GetIdx() and iaa2.GetIdx() != t2.GetIdx()) and not (iaa.GetIdx() == t1.GetIdx() and iaa2.GetIdx() == t4.GetIdx())):
                        poltype.alltorsionslist.append(get_uniq_rotbnd(poltype,iaa.GetIdx(),t2.GetIdx(),t3.GetIdx(),iaa2.GetIdx()))


            
    return




def get_torlist_opt_angle(poltype,optmol, torlist):
    tmplist = []
    for tor in torlist:
        a,b,c,d=obatom2idx(poltype,tor[0:4])
        e = optmol.GetTorsion(a,b,c,d)
        tmplist.append([a,b,c,d,e % 360])
    return tmplist

def DetermineAngleIncrementAndPointsNeededForEachTorsionSet(poltype,mol,rotbndlist):
    poltype.rotbndtoanginc={}
    poltype.torsionsettonumptsneeded={}
    # if here are multiple torsions to fit per rotatable bond, make sure there are enough angles for QM profile to do fitting
    for key in rotbndlist:
        keylist=rotbndlist[key]
        torsionlist=[]
        maintor=keylist[0]
        for tor in keylist:
            a2,b2,c2,d2=tor[0:4]
            obaa = mol.GetAtom(a2)
            obab = mol.GetAtom(b2)
            obac = mol.GetAtom(c2)
            obad = mol.GetAtom(d2)
            tpdkey = get_class_key(poltype,a2, b2, c2, d2)
            aatomicnum=obaa.GetAtomicNum()
            datomicnum=obad.GetAtomicNum()
            babelindices=[a2,b2,c2,d2]
            torsionsmissing=databaseparser.ReadTorsionList(poltype,os.path.join(os.path.abspath(os.getcwd()),poltype.torsionsmissingfilename))
            classes=[poltype.idxtosymclass[i] for i in babelindices]
            if classes not in torsionsmissing and classes[::-1] not in torsionsmissing: # then probably H torsions transferred
                continue
            if tpdkey not in torsionlist:
                torsionlist.append(tpdkey)
        torset=tuple([maintor])
        if poltype.fitfirsttorsionfoldphase==True:
            prmnum=len(poltype.nfoldlist)*len(torsionlist)+1+1

        else:
            prmnum=len(poltype.nfoldlist)*len(torsionlist)+1
        poltype.torsionsettonumptsneeded[torset]=prmnum
        if poltype.tordatapointsnum==None:
            ang=round(360/(prmnum)) # offset parameter is the +1
            if ang> 30:
                ang=30
        else:
            ang=round(360/poltype.tordatapointsnum)
        poltype.rotbndtoanginc[key]=ang
    return poltype.rotbndtoanginc

def find_tor_restraint_idx(poltype,mol,b1,b2):
    """
    Intent: Find the atoms 1 and 4 about which torsion angles are restrained
    Given b1, b2, finds the torsion: t1 b1 b2 t4
    Input:
        mol: OBMol object
        b1: first atom of the rotatable bond (t2 in the torsion)
        b2: second atom of the rotatable bond (t3 in the torsion)
    Output:
        t1: atom 1 in the torsion
        t4: atom 4 in the torsion
    Referenced By: get_torlist
    Description:
    1. Find the heaviest (heaviest meaning of highest sym class) 
       atom bound to atom b1 (that is not b2)
    2. Find the heaviest atom bound to atom b2 (that is not b1)
    3. These two atoms are returned as atoms 1 and 4 for the torsion
    """
    b1idx = b1.GetIdx()
    b2idx = b2.GetIdx()
    iteratomatom = openbabel.OBAtomAtomIter(b1)
    b1nbridx = list(map(lambda x: x.GetIdx(), iteratomatom))
    del b1nbridx[b1nbridx.index(b2idx)]    # Remove b2 from list
    assert(b1nbridx is not [])
    minb1class = min(b1nbridx,key=lambda x: poltype.idxtosymclass[x])

    iteratomatom = openbabel.OBAtomAtomIter(b2)
    b2nbridx = list(map(lambda x: x.GetIdx(), iteratomatom))
    del b2nbridx[b2nbridx.index(b1idx)]    # Remove b1 from list
    assert(b2nbridx is not [])
    minb2class = min(b2nbridx, key= lambda x:poltype.idxtosymclass[x])

    t1 = mol.GetAtom(minb1class)
    t4 = mol.GetAtom(minb2class)

    return t1,t4


def GrabFirstHeavyAtomIdx(poltype,indices,mol):
    atoms=[mol.GetAtom(i) for i in indices]
    atomicnums=[a.GetAtomicNum() for a in atoms]
    
    heavyidx=indices[0]
    for i in range(len(indices)):
        idx=indices[i]
        atomnum=atomicnums[i]
        if atomnum!=1:
            heavyidx=idx
    return heavyidx


def ConvertTinktoXYZ(poltype,filename,newfilename):
    temp=open(os.getcwd()+r'/'+filename,'r')
    tempwrite=open(os.getcwd()+r'/'+newfilename,'w')
    results=temp.readlines()
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx==0:
            linesplit=line.split()
            tempwrite.write(linesplit[0]+'\n')
            tempwrite.write('\n')
            tempwrite.flush()
            os.fsync(tempwrite.fileno())
        else:
            linesplit=line.split()
            newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
            tempwrite.write(newline)
            tempwrite.flush()
            os.fsync(tempwrite.fileno())
    temp.close()
    tempwrite.close()
    return filename.replace('.xyz_2','_xyzformat.xyz')


def CreatePsi4TorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist,mol,currentopt):

    prefix='%s-opt-%s_' % (poltype.molecprefix,currentopt)
    postfix='.psi4'  
    inputname,angles=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)

    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (mol.GetTotalCharge(),1))
    iteratom = openbabel.OBMolAtomIter(optmol)
    etab = openbabel.OBElementTable()
    if os.path.isfile(torxyzfname):
        xyzstr = open(torxyzfname,'r')
        xyzstrl = xyzstr.readlines()
        i = 0
        for atm in iteratom:
            i = i + 1
            ln = xyzstrl[i]
            temp.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), float(ln.split()[2]),float(ln.split()[3]),float(ln.split()[4])))
        xyzstr.close()
    temp.write('}'+'\n')

    # Fix all torsions around the rotatable bond b-c
    temp.write('set optking { '+'\n')
    temp.write('  frozen_dihedral = ("'+'\n')
    restlist=[]
    firsttor=True
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[0:4]
        key=' '.join([str(b),str(c)])
        if key in poltype.rotbndlist.keys():
            count=0
            for resttors in poltype.rotbndlist[key]:
                rta,rtb,rtc,rtd = resttors
                indices=[rta,rtb,rtc,rtd]
                allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,indices,mol)
                allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,indices,mol)

                aatom=mol.GetAtom(rta)
                datom=mol.GetAtom(rtd)
                rtaatomicnum=aatom.GetAtomicNum()
                rtdatomicnum=datom.GetAtomicNum()
                if (rtaatomicnum==1 or rtdatomicnum==1) and (allhydtors==False and allhydtorsoneside==False):
                    continue
                if (allhydtorsoneside==True or allhydtors==True) and (rtaatomicnum==1 or rtdatomicnum==1):

               
                    if count>=1:
                        continue
                    count+=1

                if resttors not in variabletorlist:
                    rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
                    if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
                       (optmol.GetAtom(rtd).GetAtomicNum() != 1):
                        restlist.append(resttors)
                        if not firsttor:
                            temp.write(', %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                        else:
                            temp.write('    %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                            firsttor=True
            
        # Leave all torsions around other rotatable bonds fixed
    for rotkey,torsions in poltype.rotbndlist.items():
        count=0
        for resttors in torsions:
            rta,rtb,rtc,rtd = resttors
            indices=[rta,rtb,rtc,rtd]
            allhydtors=databaseparser.CheckIfAllTorsionsAreHydrogen(poltype,indices,mol)
            allhydtorsoneside=databaseparser.CheckIfAllTorsionsAreHydrogenOneSide(poltype,indices,mol)

            aatom=mol.GetAtom(rta)
            datom=mol.GetAtom(rtd)
            rtaatomicnum=aatom.GetAtomicNum()
            rtdatomicnum=datom.GetAtomicNum()
            if (rtaatomicnum==1 or rtdatomicnum==1) and (allhydtors==False and allhydtorsoneside==False):
                continue
            if (allhydtorsoneside==True or allhydtors==True) and (rtaatomicnum==1 or rtdatomicnum==1):

           
                if count>=1:
                    continue
                count+=1

            rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
            if resttors not in restlist and resttors not in variabletorlist:
                restlist.append(resttors)
                if not firsttor:
                    temp.write(', %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                else:
                    temp.write('    %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                    firsttor=True
    temp.write('  ")'+'\n')
    temp.write('}'+'\n')

    if poltype.toroptpcm==True:
        temp.write('set {'+'\n')
        temp.write('  g_convergence GAU_LOOSE'+'\n')
        temp.write('  scf_type pk'+'\n')
        temp.write('  pcm true'+'\n')
        temp.write('  pcm_scf_type total '+'\n')
        temp.write('  geom_maxiter '+str(poltype.optmaxcycle)+'\n')
        temp.write('}'+'\n')
        temp.write('pcm = {'+'\n')
        temp.write('  Units = Angstrom'+'\n')
        temp.write('  Medium {'+'\n')
        temp.write('    SolverType = IEFPCM'+'\n')
        temp.write('    Solvent = Water'+'\n')
        temp.write('  }'+'\n')
        temp.write('    Cavity {'+'\n')
        temp.write('    RadiiSet = UFF'+'\n')
        temp.write('    Type = GePol'+'\n')
        temp.write('    Scaling = False'+'\n')
        temp.write('    Area = 0.3'+'\n')
        temp.write('    Mode = Implicit'+'\n')
        temp.write('  }'+'\n')
        temp.write('}'+'\n')
    else:
        temp.write('set {'+'\n')
        temp.write('  geom_maxiter '+str(poltype.optmaxcycle)+'\n')
        temp.write('  g_convergence GAU_LOOSE'+'\n')
        temp.write('  dynamic_level 1'+'\n')
        temp.write('}'+'\n')

    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('for _ in range(1):'+'\n')
    temp.write('  try:'+'\n')
    if poltype.toroptpcm==True:
        temp.write('    set opt_coordinates cartesian'+'\n')
    spacedformulastr=optmol.GetSpacedFormula()
    if ('I ' in spacedformulastr):
        temp.write('    basis {'+'\n')
        temp.write('    ['+' '+poltype.toroptbasissetfile+' '+poltype.iodinetoroptbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.toroptbasissetfile,poltype.iodinetoroptbasissetfile,'    ')
        temp.write('    }'+'\n')
        temp.write('    set opt_coordinates both'+'\n')
        temp.write("    optimize('%s')" % (poltype.toroptmethod.lower())+'\n')

    else:

        temp.write('    set opt_coordinates both'+'\n')
        temp.write("    optimize('%s/%s')" % (poltype.toroptmethod.lower(),poltype.toroptbasisset)+'\n')
    if poltype.freq:
        temp.write('    scf_e,scf_wfn=freq(%s/%s,return_wfn=True)'%(poltype.toroptmethod.lower(),poltype.toroptbasisset)+'\n')
    temp.write('    break'+'\n')
    temp.write('  except OptimizationConvergenceError:'+'\n')
    temp.write('    break'+'\n')
    temp.write('  else:'+'\n')
    temp.write('    try:'+'\n')
    temp.write('      set opt_coordinates cartesian'+'\n')
    if ('I ' in spacedformulastr):
        temp.write("      optimize('%s')" % (poltype.toroptmethod.lower())+'\n')
    else:
        temp.write("      optimize('%s/%s')" % (poltype.toroptmethod.lower(),poltype.toroptbasisset)+'\n')

    if poltype.freq:
        temp.write('      scf_e,scf_wfn=freq(%s/%s,return_wfn=True)'%(poltype.toroptmethod.lower(),poltype.toroptbasisset)+'\n')
    temp.write('      break'+'\n')
    temp.write('    except OptimizationConvergenceError:'+'\n')
    temp.write('      '+'pass'+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    outputname=inputname.replace('.psi4','.log')
    return inputname,outputname


def gen_torcomfile (poltype,comfname,numproc,maxmem,maxdisk,prevstruct,xyzf,mol):
    """
    Intent: Create *.com file for qm torsion calculations 
    Input:
        comfname: com file name
        numproc: number of processors
        maxmem: max memory size
        prevstruct: OBMol object
        xyzf: xyzfile with information to create *.com file
    Output:
        *.com is written
    Referenced By: tor_opt_sp 
    Description: -
    """
    opt.write_com_header(poltype,comfname,os.path.splitext(comfname)[0] + ".chk",maxdisk,maxmem,numproc)
    tmpfh = open(comfname, "a")

    optimizeoptlist = ["ModRedundant","maxcycle=%s"%(poltype.optmaxcycle),'Loose']

    optstr=opt.gen_opt_str(poltype,optimizeoptlist)
    if ('-opt-' in comfname):
        if ('I ' in poltype.mol.GetSpacedFormula()):
            poltype.toroptbasisset='gen'
            iodinebasissetfile=poltype.iodinetoroptbasissetfile
            basissetfile=poltype.toroptbasissetfile

        if poltype.toroptpcm==True:
            operationstr = "%s %s/%s SCRF=(PCM)" % (optstr,poltype.toroptmethod,poltype.toroptbasisset)
        else:
            operationstr = "%s %s/%s" % (optstr,poltype.toroptmethod,poltype.toroptbasisset)
        commentstr = poltype.molecprefix + " Rotatable Bond Optimization on " + gethostname()
    else:
        if ('I ' in poltype.mol.GetSpacedFormula()):
            prevbasisset=poltype.torspbasisset
            poltype.torspbasisset='gen'
            iodinebasissetfile=poltype.iodinetorspbasissetfile
            basissetfile=poltype.torspbasissetfile

        if poltype.torsppcm==True:
            operationstr = "#P %s/%s SP SCF=(qc,maxcycle=800) SCRF=(PCM) Pop=NBORead" % (poltype.torspmethod,poltype.torspbasisset)
        else:       
            operationstr = "#P %s/%s SP SCF=(qc,maxcycle=800) Pop=NBORead" % (poltype.torspmethod,poltype.torspbasisset)

        commentstr = poltype.molecprefix + " Rotatable Bond SP Calculation on " + gethostname()   
    if ('I ' in poltype.mol.GetSpacedFormula()):
        operationstr+=' pseudo=read'
    string=' MaxDisk=%s \n'%(maxdisk)
    operationstr+=string

    tmpfh.write(operationstr)
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d\n' % (mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity()))
    iteratom = openbabel.OBMolAtomIter(prevstruct)
    etab = openbabel.OBElementTable()
    if xyzf!=None:
        if os.path.isfile(xyzf):
            xyzstr = open(xyzf,'r')
            xyzstrl = xyzstr.readlines()
            i = 0
            for atm in iteratom:
                i = i + 1
                ln = xyzstrl[i]
                tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), float(ln.split()[2]),float(ln.split()[3]),float(ln.split()[4])))
            tmpfh.write('\n')
            xyzstr.close()
    else:
        for atm in iteratom:
            tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), atm.x(),atm.y(),atm.z()))
        tmpfh.write('\n')

    if ('I ' in poltype.mol.GetSpacedFormula()):
        formulalist=poltype.mol.GetSpacedFormula().lstrip().rstrip().split()
        elementtobasissetlines=GenerateElementToBasisSetLines(poltype,poltype.basissetpath+basissetfile)
        for element,basissetlines in elementtobasissetlines.items():
            if element in poltype.mol.GetSpacedFormula():
                for line in basissetlines: 
                    tmpfh.write(line)



        temp=open(poltype.basissetpath+iodinebasissetfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '!' not in line:
                tmpfh.write(line)
        tmpfh.write('\n')



        
    if 'opt' not in comfname and poltype.dontfrag==False: 
        tmpfh.write('$nbo bndidx $end'+'\n')
        tmpfh.write('\n')



    tmpfh.close()


def GenerateElementToBasisSetLines(poltype,basissetfile):
    elementtobasissetlines={}
    temp=open(basissetfile,'r')
    results=temp.readlines()
    temp.close()
    lines=[]
    element=None
    for line in results:
        linesplit=line.split()
        if len(linesplit)==2 and linesplit[0].isalpha() and linesplit[1]=='0':
            if element!=None:
                elementtobasissetlines[element]=lines
            element=linesplit[0]
            lines=[line]
        else:
            lines.append(line)
    return elementtobasissetlines
 



def save_structfile(poltype,molstruct, structfname):
    """
    Intent: Output the data in the OBMol structure to a file (such as *.xyz)
    Input:
        molstruct: OBMol structure
        structfname: output file name
    Output:
        file is output to structfname
    Referenced By: tor_opt_sp, compute_mm_tor_energy
    Description: -
    """
    strctext = os.path.splitext(structfname)[1]
    tmpconv = openbabel.OBConversion()
    if strctext in '.xyz':
        tmpfh = open(structfname, "w")
        iteratom = openbabel.OBMolAtomIter(molstruct)
        etab = openbabel.OBElementTable()
        tmpfh.write('%6d   %s\n' % (molstruct.NumAtoms(), molstruct.GetTitle()))
        for ia in iteratom:
            tmpfh.write( '%6d %2s %13.6f %11.6f %11.6f %5d' % (ia.GetIdx(), etab.GetSymbol(ia.GetAtomicNum()), ia.x(), ia.y(), ia.z(), poltype.idxtosymclass[ia.GetIdx()]))
            iteratomatom = openbabel.OBAtomAtomIter(ia)
            neighbors = []
            for iaa in iteratomatom:
                neighbors.append(iaa.GetIdx())
            neighbors = sorted(neighbors)
            for iaa in neighbors:
                tmpfh.write('%5d' % iaa)
            tmpfh.write('\n')
    else:
        inFormat = openbabel.OBConversion.FormatFromExt(structfname)
        tmpconv.SetOutFormat(inFormat)
    return tmpconv.WriteFile(molstruct, structfname)

def save_structfileXYZ(poltype,molstruct, structfname):
        
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('xyz')
    return tmpconv.WriteFile(molstruct, structfname)

def get_class_key(poltype,a, b, c, d):
    """
    Intent: Given a set of atom idx's, return the class key for the set (the class numbers of the atoms appended together)
    """
    cla = poltype.idxtosymclass[a]
    clb = poltype.idxtosymclass[b]
    clc = poltype.idxtosymclass[c]
    cld = poltype.idxtosymclass[d]
    if ((clb > clc) or (clb == clc and cla > cld)):
        return '%d %d %d %d' % (cld, clc, clb, cla)
    return '%d %d %d %d' % (cla, clb, clc, cld)

def get_uniq_rotbnd(poltype,a, b, c, d):
    """
    Intent: Return the atom idx's defining a rotatable bond in the order of the class key
    found by 'get_class_key'
    """
    cla = poltype.idxtosymclass[a]
    clb = poltype.idxtosymclass[b]
    clc = poltype.idxtosymclass[c]
    cld = poltype.idxtosymclass[d]

    tmpkey = '%d %d %d %d' % (cla,clb,clc,cld)
    if (get_class_key(poltype,a,b,c,d) == tmpkey):
        return (a, b, c, d)
    return (d, c, b, a)

def obatom2idx(poltype,obatoms):
    """
    Intent: Given a list of OBAtom objects, return a list of their corresponding idx's
    Referenced By: get_torlist_opt_angle
    """
    atmidxlist = []
    for obatm in obatoms:
        atmidxlist.append(obatm.GetIdx())
    return atmidxlist

def rads(poltype,degrees):
    """
    Intent: Convert degrees to radians
    """
    if type(degrees) == type(list):
        return [ deg * numpy.pi / 180 for deg in degrees ]
    else:
        return degrees * numpy.pi / 180

def write_arr_to_file(poltype,fname, array_list):
    """
    Intent: Write out information in array to file
    Input:
        fname: file name
        array_list: array with data to be printed
    Output:
        file is written to 'filename'
    Referenced By: fit_rot_bond_tors, eval_rot_bond_parms
    Description: - 
    """
    outfh = open(fname,'w')
    rows = zip(*array_list)
    for cols in rows:
        for ele in cols:
            outfh.write(str(ele)+' ')
        outfh.write("\n")
    outfh.close()
    return 

def CreatePsi4TorESPInputFile(poltype,prevstrctfname,optmol,torset,phaseangles,mol,makecube=None):
    prefix='%s-sp-'%(poltype.molecprefix)
    postfix='.psi4'  
    inputname,angles=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
    finalstruct= opt.load_structfile(poltype,prevstrctfname)

    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (mol.GetTotalCharge(), 1))
    iteratom = openbabel.OBMolAtomIter(finalstruct)
    etab = openbabel.OBElementTable()
    for atm in iteratom:
        temp.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()),atm.GetX(),atm.GetY(),atm.GetZ()))
    temp.write('}'+'\n')
    if poltype.torsppcm==True:
        temp.write('set {'+'\n')
        temp.write(' scf_type pk'+'\n')
        temp.write(' pcm true'+'\n')
        temp.write(' pcm_scf_type total '+'\n')
        temp.write('}'+'\n')
        temp.write('pcm = {'+'\n')
        temp.write(' Units = Angstrom'+'\n')
        temp.write(' Medium {'+'\n')
        temp.write(' SolverType = IEFPCM'+'\n')
        temp.write(' Solvent = Water'+'\n')
        temp.write(' }'+'\n')
        temp.write(' Cavity {'+'\n')
        temp.write(' RadiiSet = UFF'+'\n')
        temp.write(' Type = GePol'+'\n')
        temp.write(' Scaling = False'+'\n')
        temp.write(' Area = 0.3'+'\n')
        temp.write(' Mode = Implicit'+'\n')
        temp.write(' }'+'\n')
        temp.write('}'+'\n')
    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('set freeze_core True'+'\n')
    spacedformulastr=optmol.GetSpacedFormula()
    if ('I ' in spacedformulastr):
        temp.write('basis {'+'\n')
        temp.write('['+' '+poltype.torspbasissetfile+' '+poltype.iodinetorspbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.torspbasissetfile,poltype.iodinetorspbasissetfile,'')
        temp.write('}'+'\n')
        temp.write("E, wfn = energy('%s',return_wfn=True)" % (poltype.torspmethod.lower())+'\n')

    else:

        temp.write("E, wfn = energy('%s/%s',return_wfn=True)" % (poltype.torspmethod.lower(),poltype.torspbasisset)+'\n')
    temp.write('oeprop(wfn,"WIBERG_LOWDIN_INDICES")'+'\n')

    temp.write('clean()'+'\n')
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname


def ReadInBasisSet(poltype,tmpfh,normalelementbasissetfile,otherelementbasissetfile,space):
    newtemp=open(poltype.basissetpath+normalelementbasissetfile,'r')
    results=newtemp.readlines()
    newtemp.close()
    for line in results:
        if '!' not in line:
            tmpfh.write(space+line)


    newtemp=open(poltype.basissetpath+otherelementbasissetfile,'r')
    results=newtemp.readlines()
    newtemp.close()
    for line in results:
        if '!' not in line:
            tmpfh.write(space+line)
    return tmpfh



def RemoveDuplicateRotatableBondTypes(poltype):
    tortorotbnd={}
    for key,tors in poltype.rotbndlist.items():
        for tor in tors:
            classkeynew=get_class_key(poltype,tor[0],tor[1],tor[2],tor[3])
            tortorotbnd[tuple(tor)]=classkeynew
    listofduptors=[]
    for key,value in tortorotbnd.items():
        duptors=[k for k,v in tortorotbnd.items() if v == value]
        keylist=[]
        for tor in duptors:
            rotbndkey=str(tor[1])+' '+str(tor[2])
            if rotbndkey not in keylist:
                keylist.append(rotbndkey)
        if len(keylist)>=2 and duptors not in listofduptors:
            listofduptors.append(duptors)
    for dup in listofduptors: # doesnt matter which one is first, just remove duplicates
        for i in range(len(dup)-1):
            tor=list(dup[i])
            if tor in poltype.torlist:
                poltype.torlist.remove(tor)
    
    return poltype.torlist,poltype.rotbndlist 



def PrependStringToKeyfile(poltype,keyfilename,string):
    """
    Intent: Adds a header to the key file given by 'keyfilename'
    """
    tmpfname = keyfilename + "_tmp"
    tmpfh = open(tmpfname, "w")
    keyfh = open(keyfilename, "r")
    tmpfh.write(string+'\n')

    for line in keyfh:
        tmpfh.write(line)
    shutil.move(tmpfname, keyfilename)

def RemoveStringFromKeyfile(poltype,keyfilename,string):
    """
    Intent: Adds a header to the key file given by 'keyfilename'
    """
    tmpfname = keyfilename + "_tmp"
    tmpfh = open(tmpfname, "w")
    keyfh = open(keyfilename, "r")

    for line in keyfh:
        if string in line:
            pass
        else:
            tmpfh.write(line)
    shutil.move(tmpfname, keyfilename)





def UpdateMaxRange(poltype,torsion,maxrange):
    a,b,c,d=torsion[0:4]
    key=str(b)+' '+str(c)
    poltype.rotbndtomaxrange[key]=maxrange


    
 



def FindMainConsecutiveTorsions(poltype,tortorsmissing):
    tortors=[]
    for tortor in tortorsmissing:
        a,b,c,d,e=tortor[:]
        aidx=a+1
        bidx=b+1
        cidx=c+1
        didx=d+1
        eidx=e+1
        temp=[]
        firstrot=[bidx,cidx]
        firstrot=[str(i) for i in firstrot]
        secondrot=[cidx,didx]
        secondrot=[str(i) for i in secondrot]
        firstkey=' '.join(firstrot)
        revfirstkey=' '.join(firstrot[::-1])
        secondkey=' '.join(secondrot)
        revsecondkey=' '.join(secondrot[::-1])
        first=tuple([aidx,bidx,cidx,didx])
        second=tuple([bidx,cidx,didx,eidx])
        if firstkey in poltype.rotbndlist.keys():
            temp.append(first)
        elif revfirstkey in poltype.rotbndlist.keys():
            temp.append(first[::-1])
        if secondkey in poltype.rotbndlist.keys():
            temp.append(second)
        elif revsecondkey in poltype.rotbndlist.keys():
            temp.append(second[::-1])

        tortors.append(temp)
    return tortors

def RemoveRotatableBondFromTorlist(poltype,tor):
    indicestoremove=[]
    rotbnd=[tor[1],tor[2]]
    for grpidx in range(len(poltype.torlist)):
        grp=poltype.torlist[grpidx]
        ls=grp[0]
        other=[ls[1],ls[2]]
        if rotbnd==other or rotbnd[::-1]==other:
            indicestoremove.append(grpidx)
    for index in indicestoremove:
        del poltype.torlist[index]



def UpdateTorsionSets(poltype,tortor):
    for tor in tortor:
        torset=tuple([tuple(tor)])
        RemoveRotatableBondFromTorlist(poltype,tor)
    poltype.torsettovariabletorlist[tuple(tortor)]=[]
    poltype.torlist.append(tuple(tortor))

def FlattenArray(poltype,flatphaselist):
    flatphaselist= flatphaselist.reshape(-1, flatphaselist.shape[-1])
    return flatphaselist 

def Determine1DTorsionSlicesOnTorTorSurface(poltype,energyarray,phaseanglelist,torset):
    indicestokeep=[]
    firsttorindices=[]
    secondtorindices=[]
    arraylength=len(energyarray)
    sqrt=int(numpy.sqrt(arraylength))
    energymatrix=energyarray.reshape((sqrt,sqrt))
    minvalue=numpy.amin(energymatrix)
    minindices=numpy.transpose(numpy.where(energymatrix==minvalue))[0]
    rowindex=minindices[0]
    colindex=minindices[1]
    phaseanglematrix=phaseanglelist.reshape((sqrt,sqrt,2))
    rowslice=phaseanglematrix[rowindex,:,:]
    colslice=phaseanglematrix[:,colindex,:]
    for phase in rowslice:
        phase=phase.tolist()
        indicestokeep.append(phase)
        secondtorindices.append(phase)
    poltype.torsettotortorphaseindicestokeep[tuple([torset[1]])]=secondtorindices

    for phase in colslice:
        phase=phase.tolist()
        indicestokeep.append(phase)
        firsttorindices.append(phase)
    poltype.torsettotortorphaseindicestokeep[tuple([torset[0]])]=firsttorindices
    return indicestokeep,firsttorindices,secondtorindices

def RemoveTorTorPoints(poltype,energyarray,phaseanglelist,indicestokeep,listofstructurestorunQM):
    energyarray=numpy.array(energyarray)
    listofstructurestorunQM=numpy.array(listofstructurestorunQM)
    sortedenergyarrayindices=numpy.argsort(energyarray)
    sortedenergyarray=energyarray[sortedenergyarrayindices]
    sortedphaseanglelist=phaseanglelist[sortedenergyarrayindices]
    sortedlistofstructurestorunQM=listofstructurestorunQM[sortedenergyarrayindices]
    numpts=len(sortedenergyarray)
    newphasearray=[]
    newlistofstructurestorunQM=[]
    for i in range(len(sortedphaseanglelist)):
        struct=sortedlistofstructurestorunQM[i]
        phase=sortedphaseanglelist[i]
        phase=phase.tolist()
        found=phase in indicestokeep
        if found==True:
            newphasearray.append(phase)
            newlistofstructurestorunQM.append(struct)
    for i in range(len(sortedphaseanglelist)):
        struct=sortedlistofstructurestorunQM[i]
        phase=sortedphaseanglelist[i]
        phase=phase.tolist()
        found=phase in newphasearray
        if found==False:
            currentlength=len(newphasearray)
            if currentlength<=poltype.defaultmaxtorsiongridpoints:
                newphasearray.append(phase)
                newlistofstructurestorunQM.append(struct)

    correctorderphasearray=[]
    correctlistofstructurestorunQM=[]
    for i in range(len(phaseanglelist)):
        struct=listofstructurestorunQM[i]
        phase=phaseanglelist[i]
        phase=phase.tolist()
        found=phase in newphasearray
        if found==True:
            correctorderphasearray.append(phase)
            correctlistofstructurestorunQM.append(struct)
    return correctorderphasearray,correctlistofstructurestorunQM
            

def PrepareTorsionTorsion(poltype,optmol,mol,tortorsmissing):
    tortors=FindMainConsecutiveTorsions(poltype,tortorsmissing)
    for tortor in tortors:
        UpdateTorsionSets(poltype,tortor)

