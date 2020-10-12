import symmetry as symm
import optimization as opt
import electrostaticpotential as esp
import valence
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
from itertools import product

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

def ExecuteOptJobs(poltype,listofstructurestorunQM,phaselist,optmol,torset,variabletorlist,torsionrestraint):
    jobtooutputlog={}
    listofjobs=[]
    outputlogs=[]
    for i in range(len(listofstructurestorunQM)):
        torxyzfname=listofstructurestorunQM[i]
        phaseangles=phaselist[i]
        inputname,outputlog,cmdstr,scratchdir=GenerateTorsionOptInputFile(poltype,torxyzfname,torset,phaseangles,optmol,variabletorlist)
        finished,error=poltype.CheckNormalTermination(outputlog)
        if finished==True and 'opt' in outputlog:
            opt.GrabFinalXYZStructure(poltype,outputlog,outputlog.replace('.log','.xyz'))

        if finished==False:
            if os.path.isfile(outputlog):
                statinfo=os.stat(outputlog)
                size=statinfo.st_size
                if size!=0:
                    os.remove(outputlog)
                    listofjobs.append(cmdstr)
                    jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputlog
            else:
                listofjobs.append(cmdstr)
                jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputlog


        outputlogs.append(outputlog)
        
    return outputlogs,listofjobs,scratchdir,jobtooutputlog

def TinkerMinimizePostQMOpt(poltype,outputlogs,phaselist,optmol,torset,variabletorlist,torsionrestraint,finishedjobs,errorjobs,bondtopology):
    torxyznames=[]
    finishedoutputlogs=[]
    cartxyznames=[]
    for i in range(len(outputlogs)):
        outputlog=outputlogs[i]
        phaseangles=phaselist[i]
        if outputlog in finishedjobs and outputlog not in errorjobs:
            if not poltype.use_gaus:
                finalstruct=outputlog.replace('.log','.xyz')
                cartxyz,torxyzfname=tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,finalstruct,torsionrestraint,'_postQMOPTprefit',poltype.key4fname,'../',bondtopology)
            else:
                cartxyz,torxyzfname=tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,outputlog,torsionrestraint,'_postQMOPTprefit',poltype.key4fname,'../',bondtopology)
            torxyznames.append(torxyzfname)
            finishedoutputlogs.append(outputlog)
            cartxyznames.append(cartxyz)
    return torxyznames,finishedoutputlogs,cartxyznames

def ExecuteSPJobs(poltype,torxyznames,cartxyznames,optoutputlogs,phaselist,optmol,torset,variabletorlist,torsionrestraint,outputlogtophaseangles,outputlogtocartxyz):
    jobtooutputlog={}
    listofjobs=[]
    outputnames=[]
    for i in range(len(torxyznames)):
        torxyzfname=torxyznames[i]
        phaseangles=phaselist[i]
        outputlog=optoutputlogs[i]
        cartxyzname=cartxyznames[i]
        if not poltype.use_gaus:
            finalstruct=outputlog.replace('.log','_opt.xyz')
            inputname,outputname=CreatePsi4TorESPInputFile(poltype,finalstruct,torxyzfname,optmol,torset,phaseangles)
            cmdstr='cd '+os.getcwd()+' && '+'psi4 '+inputname+' '+outputname
        else:
            inputname,outputname=GenerateTorsionSPInputFileGaus(poltype,torxyzfname,torset,optmol,phaseangles,outputlog)
            cmdstr = 'cd '+os.getcwd()+' && '+'GAUSS_SCRDIR='+poltype.scrtmpdirgau+' '+poltype.gausexe+' '+inputname
        finished,error=poltype.CheckNormalTermination(outputname)
        if finished==True and error==False:
            pass
        else:
            if os.path.isfile(outputname):
                statinfo=os.stat(outputname)
                size=statinfo.st_size
                if size!=0:
                    os.remove(outputname)
                    listofjobs.append(cmdstr)
                    jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputname
            else:
                listofjobs.append(cmdstr)
                jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputname

        outputnames.append(outputname)
        outputlogtophaseangles[outputname]=phaseangles
        outputlogtocartxyz[outputname]=cartxyzname
    if not poltype.use_gaus:
        
        return outputnames,listofjobs,poltype.scrtmpdirpsi4,jobtooutputlog,outputlogtophaseangles,outputlogtocartxyz
    else:
        return outputnames,listofjobs,poltype.scrtmpdirgau,jobtooutputlog,outputlogtophaseangles,outputlogtocartxyz



def CreateGausTorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist):
    toroptcomfname = '%s-opt-' % (poltype.molecprefix)
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[:]
        phaseangle=phaseangles[i]
        torang = optmol.GetTorsion(a,b,c,d)
        toroptcomfname+='%d-%d-%03d' % (b,c,round((torang+phaseangle)%360))
        toroptcomfname+='_'
    toroptcomfname=toroptcomfname[:-1]
    toroptcomfname+='.com'     
    strctfname = os.path.splitext(toroptcomfname)[0] + '.log'
    gen_torcomfile(poltype,toroptcomfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,optmol,torxyzfname) # prevstruct is just used for Atomic Num and charge,torxyzfname is used for xyz coordinates
    # Append restraints to *opt*.com file
    tmpfh = open(toroptcomfname, "a")
    restlist=[]
    # Fix all torsions around the rotatable bond b-c
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[:]
        for resttors in poltype.rotbndlist[' '.join([str(b),str(c)])]:
            rta,rtb,rtc,rtd = resttors
            rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
            if resttors not in variabletorlist and resttors not in restlist:
                if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
                   (optmol.GetAtom(rtd).GetAtomicNum() != 1):
                    tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))
                    restlist.append(resttors)

    # Leave all torsions around other rotatable bonds fixed
    for rotkey,torsions in poltype.rotbndlist.items():
        for resttors in torsions:
            rta,rtb,rtc,rtd = resttors
            rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
            if resttors not in restlist and resttors not in variabletorlist:
                if (optmol.GetAtom(rta).GetAtomicNum() != 1) and (optmol.GetAtom(rtd).GetAtomicNum() != 1):
                    tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))
                    restlist.append(resttors)
    
    tmpfh.write("\n")
    tmpfh.close()
    outputname= os.path.splitext(toroptcomfname)[0] + '.log'
    return toroptcomfname,outputname
    

def GenerateTorsionOptInputFile(poltype,torxyzfname,torset,phaseangles,optmol,variabletorlist):
    if  poltype.use_gaus==False and poltype.use_gausoptonly==False:
        inputname,outputname=CreatePsi4TorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist)
        cmdstr='cd '+os.getcwd()+' && '+'psi4 '+inputname+' '+outputname
    else:
        inputname,outputname=CreateGausTorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist)
        cmdstr='cd '+os.getcwd()+' && '+'GAUSS_SCRDIR='+poltype.scrtmpdirgau+' '+poltype.gausexe+' '+inputname
    if poltype.use_gaus==False and poltype.use_gausoptonly==False:
        return inputname,outputname,cmdstr,poltype.scrtmpdirpsi4

    else:
        return inputname,outputname,cmdstr,poltype.scrtmpdirgau
        

def GenerateTorsionSPInputFileGaus(poltype,torxyzfname,torset,optmol,phaseangles,prevstrctfname):
    prevstruct = opt.load_structfile(poltype,prevstrctfname)
    torspcomfname = '%s-sp-' % (poltype.molecprefix)
    for i in range(len(torset)):
        tor=torset[i]
        phaseangle=phaseangles[i]
        a,b,c,d=tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        torspcomfname+='%d-%d-%03d' % (b,c,round((torang+phaseangle)%360))
        torspcomfname+='_'
    torspcomfname=torspcomfname[:-1]
    torspcomfname+='.com'     
    torsplogfname = os.path.splitext(torspcomfname)[0] + '.log'
    gen_torcomfile(poltype,torspcomfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,prevstruct,torxyzfname)
    outputname=torspcomfname.replace('.com','.log')
    return torspcomfname,outputname


def tinker_minimize_angles(poltype,torset,optmol,variabletorlist,phaselist,prevstrctfname,torsionrestraint,bondtopology):
    tinkerstructnamelist=[]
    # load prevstruct
    firsttor=torset[0]
    a,b,c,d=firsttor[0:4]
    tor=tuple([a,b,c,d])
    if tor in poltype.nonaroringtors:
        nonaroringtor=True
    else:
        nonaroringtor=False
    # create xyz and key and write restraint then minimize, getting .xyz_2
    for phaseanglelist in phaselist: # we need to send back minimized structure in XYZ (not tinker) format to load for next tinker minimization,but append the xyz_2 tinker XYZ file so that com file can be generated from that 
        prevstruct = opt.load_structfile(poltype,prevstrctfname)
        prevstruct = opt.PruneBonds(poltype,prevstruct,bondtopology) # sometimes extra bonds are made when atoms get too close during minimization
        prevstruct=opt.rebuild_bonds(poltype,prevstruct,optmol)
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat('mol')
        obConversion.SetOutFormat('mol2')
        obConversion.WriteFile(prevstruct, 'temp.mol2')
        try:
            rdmol=MolFromMol2File('temp.mol2',False,False)   
        except:
            rdmol=MolFromMol2File('temp.mol2',True,False)
        if nonaroringtor==False:
            conf = rdmol.GetConformer()
            for i in range(len(torset)):
                tor=torset[i]
                phaseangle=phaseanglelist[i]
                a,b,c,d=tor[0:4]
                dihedral = optmol.GetTorsion(a,b,c,d)
                newdihedral=round((dihedral+phaseangle)%360)
                rdmt.SetDihedralDeg(conf, a-1, b-1, c-1, d-1, newdihedral)
                rdmol.UpdatePropertyCache(strict=False)

        try:
            print(Chem.MolToMolBlock(rdmol,kekulize=True),file=open('tempout.mol','w+'))
        except:
            print(Chem.MolToMolBlock(rdmol,kekulize=False),file=open('tempout.mol','w+'))

        obConversion.ReadFile(prevstruct, 'tempout.mol')
        

        prevstrctfname,torxyzfname,newtorxyzfname,keyfname=tinker_minimize(poltype,torset,optmol,variabletorlist,phaseanglelist,torsionrestraint,prevstruct,'_preQMOPTprefit',poltype.key4fname,'../')
        tinkerstructnamelist.append(newtorxyzfname)
    return tinkerstructnamelist


def tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,prevstrctfname,torsionrestraint,designatexyz,keybase,keybasepath,bondtopology):
    for i in range(len(phaseangles)):
        phaseangle=phaseangles[i]
        tor=tuple(torset[i])
        a,b,c,d=tor[0:4] 
        if tor in poltype.nonaroringtors:
            nonaroringtor=True
        else:
            nonaroringtor=False
        if 'post' not in designatexyz and nonaroringtor==False:    
            dihedral = optmol.GetTorsion(a,b,c,d)
            dihedral=round((dihedral+phaseangle)%360)
            currentdihedral=prevstruct.GetTorsion(a,b,c,d)
            currentdihedral=round((currentdihedral)%360)
            diff= numpy.abs(currentdihedral-dihedral)
            tol=.01
            if diff>tol and diff!=360:
                raise ValueError('Difference of '+str(diff)+' is greater than '+str(tol)+' for target dihedral of '+str(dihedral)+' and current dihedral of '+str(currentdihedral)+' '+os.getcwd())

    prevstruct = opt.load_structfile(poltype,prevstrctfname) # this should be a logfile
    prevstruct=opt.rebuild_bonds(poltype,prevstruct,optmol)
    prevstruct = opt.PruneBonds(poltype,prevstruct,bondtopology)
    cartxyz,torxyzfname,newtorxyzfname,keyfname=tinker_minimize(poltype,torset,optmol,variabletorlist,phaseangles,torsionrestraint,prevstruct,designatexyz,keybase,keybasepath)
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

       
def tinker_minimize(poltype,torset,optmol,variabletorlist,phaseanglelist,torsionrestraint,prevstruct,designatexyz,keybase,keybasepath):

    torxyzfname = '%s-opt-'%(poltype.molecprefix)
    for i in range(len(torset)):
        tor=torset[i]
        phaseangle=phaseanglelist[i]
        a,b,c,d=tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        torxyzfname+='%d-%d-%03d' % (b,c,round((torang+phaseangle)%360))
        torxyzfname+='_'
    torxyzfname=torxyzfname[:-1]
    torxyzfname+='%s.xyz' % (designatexyz)

    tmpkeyfname = 'tmp-'
    for i in range(len(torset)):
        tor=torset[i]
        phaseangle=phaseanglelist[i]
        a,b,c,d=tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        tmpkeyfname+='%d-%d-%03d' % (b,c,round((torang+phaseangle)%360))
        tmpkeyfname+='_'
    tmpkeyfname=tmpkeyfname[:-1]
    tmpkeyfname+='%s.key' % (designatexyz)

    torminlogfname=torxyzfname.replace('.xyz','.out')
    save_structfile(poltype,prevstruct,torxyzfname)
    shutil.copy(keybasepath+keybase, tmpkeyfname)
    tmpkeyfh = open(tmpkeyfname,'a')
    restlist=[]
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        phaseangle=phaseanglelist[i]
        tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (a,b,c,d,torsionrestraint,round((torang+phaseangle)%360),round((torang+phaseangle)%360)))
        for key in poltype.rotbndlist.keys():
            torlist=poltype.rotbndlist[key]
            for res in torlist:
                resa,resb,resc,resd = res[0:4]
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
    
    cartxyz=ConvertTinktoXYZ(poltype,torxyzfname+'_2')
    return cartxyz,torxyzfname,torxyzfname+'_2',tmpkeyfname

def GenerateFilename(poltype,torset,phaseangles,prefix,postfix,mol):
    filename=prefix
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[:]
        phaseangle=phaseangles[i] 
        torang = mol.GetTorsion(a,b,c,d)
        filename+='%d-%d-%03d'%(b,c,round((torang+phaseangle)%360))
        filename+='_'
    filename=filename[:-1]
    filename+=postfix
    return filename

def gen_torsion(poltype,optmol,torsionrestraint):
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

    listofstructurestorunQM=[]
    fullrange=[]
    fulloutputlogs=[]
    fulllistofjobs=[]
    fulljobtolog={}
    fulljobtooutputlog={}
    poltype.torsettophaselist={}
    torsettooptoutputlogs={}
    bondtopology=GenerateBondTopology(poltype,optmol)
    poltype.idealangletensor={}
    poltype.tensorphases={}
    for torset in poltype.torlist:
        variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
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
            clock=list(phaselist[:int(len(phaselist)/2)+1])
            counterclock=[-1*i for i in phaselist[1:int(len(phaselist)/2)+1]][::-1]
            phaselists.append(counterclock+clock)
        currentanglelist=numpy.array(currentanglelist)
        flatphaselist=numpy.array(list(product(*phaselists)))
        origshape=flatphaselist.shape
        prms=poltype.torsionsettonumptsneeded[tuple(torset)]
        prmstoremove=len(flatphaselist)-prms
        lastdim=len(phaselists)
        shape=[]
        for i in range(len(phaselists)):
            phaselist=phaselists[i]
            shape.append(len(phaselist))
        shape.append(lastdim)
        shape=tuple(shape)
        flatphaselist=flatphaselist.reshape(shape) # now we can remove diagonal elements
        idealangletensor=flatphaselist+currentanglelist
        poltype.idealangletensor[tuple(torset)]=idealangletensor
        poltype.tensorphases[tuple(torset)]=flatphaselist
        ax2idx=len(shape)-1-1
        ax1idx=ax2idx-1

        #diags = numpy.transpose(numpy.diagonal(flatphaselist, axis1=ax1idx, axis2=ax2idx))
        flatphaselist=flatphaselist.reshape(origshape)
        print('reshaped ',flatphaselist)
        '''
        for diag in diags:
            indexes=numpy.where((flatphaselist == diag).all(axis=1))
            flatphaselist = numpy.delete(flatphaselist,indexes , axis=0)
        ''' 
        # we could remove more points if we wanted....
        phaseangles=[0]*len(torset)
        if poltype.use_gaus==False and poltype.use_gausoptonly==False:
            prefix='%s-opt-' % (poltype.molecprefix)
            postfix='-opt.xyz' 
            prevstrctfname=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
            cmd = 'cp ../%s %s' % (poltype.logoptfname.replace('.log','.xyz'),prevstrctfname)
            poltype.call_subsystem(cmd,True)

        else:
            prefix='%s-opt-' % (poltype.molecprefix)
            postfix='.log' 
            prevstrctfname=GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
            # copy *-opt.log found early by Gaussian to 'prevstrctfname'
            cmd = 'cp ../%s %s' % (poltype.logoptfname,prevstrctfname)
            poltype.call_subsystem(cmd,True)


        minstrctfname = prevstrctfname
        prevstrctfname = minstrctfname
        listoftinkertorstructures=tinker_minimize_angles(poltype,torset,optmol,variabletorlist,flatphaselist,prevstrctfname,torsionrestraint,bondtopology)
        listofstructurestorunQM.extend(listoftinkertorstructures)
        outputlogs,listofjobs,scratchdir,jobtooutputlog=ExecuteOptJobs(poltype,listoftinkertorstructures,flatphaselist,optmol,torset,variabletorlist,torsionrestraint)
        poltype.torsettophaselist[tuple(torset)]=flatphaselist
        torsettooptoutputlogs[tuple(torset)]=outputlogs


        lognames=[]
        for job in listofjobs:
            log=jobtooutputlog[job]
            lognames.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+r'/'+poltype.logfname)
        jobtolog=dict(zip(listofjobs, lognames)) 
        fulloutputlogs.extend(outputlogs)
        fulllistofjobs.extend(listofjobs)
        fulljobtolog.update(jobtolog) 
        fulljobtooutputlog.update(jobtooutputlog)

    jobtologlistfilenameprefix=os.getcwd()+r'/'+'QMOptJobToLog'+'_'+poltype.molecprefix
    if poltype.externalapi!=None:
        if len(fulllistofjobs)!=0:
            call.CallExternalAPI(poltype,fulljobtolog,jobtologlistfilenameprefix,scratchdir)
        finishedjobs,errorjobs=poltype.WaitForTermination(fulljobtooutputlog)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(fulljobtooutputlog,True)
    for outputlog in fulloutputlogs:
        finished,error=poltype.CheckNormalTermination(outputlog)
        if finished==True and 'opt' in outputlog:
            opt.GrabFinalXYZStructure(poltype,outputlog,outputlog.replace('.log','.xyz'))

        if finished==True and outputlog not in finishedjobs:
            finishedjobs.append(outputlog)
    fulltorxyznames=[]
    fullfinishedoutputlogsSP=[]
    fullcartxyznames=[]
    fulloutputlogsSP=[]
    fulllistofjobs=[]
    fulljobtolog={}
    fulljobtooutputlog={}
    torsettospoutputlogs={}
    torsettocartxyz={}
    outputlogtophaseangles={}
    outputlogtocartxyz={}
    for torset in poltype.torlist:
        variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        outputlogs=torsettooptoutputlogs[tuple(torset)]
        torxyznames,finishedoutputlogs,cartxyznames=TinkerMinimizePostQMOpt(poltype,outputlogs,flatphaselist,optmol,torset,variabletorlist,torsionrestraint,finishedjobs,errorjobs,bondtopology)
        fulltorxyznames.extend(torxyznames)
        fullfinishedoutputlogsSP.extend(finishedoutputlogs)
        fullcartxyznames.extend(cartxyznames)
        torsettocartxyz[tuple(torset)]=cartxyznames

        outputlogs,listofjobs,scratchdir,jobtooutputlog,outputlogtophaseangles,outputlogtocartxyz=ExecuteSPJobs(poltype,torxyznames,cartxyznames,finishedoutputlogs,flatphaselist,optmol,torset,variabletorlist,torsionrestraint,outputlogtophaseangles,outputlogtocartxyz)
        lognames=[]
        torsettospoutputlogs[tuple(torset)]=outputlogs
        for job in listofjobs:
           log=jobtooutputlog[job]
           lognames.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+r'/'+poltype.logfname)
        jobtolog=dict(zip(listofjobs, lognames)) 
        fulljobtooutputlog.update(jobtooutputlog)
        fulloutputlogsSP.extend(outputlogs)
        fulllistofjobs.extend(listofjobs)
        fulljobtolog.update(jobtolog)

    for torset in poltype.torlist:
        variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        outputlogs=torsettospoutputlogs[tuple(torset)]
        for i in range(len(outputlogs)):
            outputlog=outputlogs[i]
            cartxyzname=outputlogtocartxyz[outputlog] 
            phaseangles=outputlogtophaseangles[outputlog]
            poltype.optoutputtotorsioninfo[outputlog]= [torset,optmol,variabletorlist,phaseangles,cartxyzname,bondtopology]
    jobtologlistfilenameprefix=os.getcwd()+r'/'+'QMSPJobToLog'+'_'+poltype.molecprefix
    if poltype.externalapi!=None:
        if len(fulllistofjobs)!=0:
            call.CallExternalAPI(poltype,fulljobtolog,jobtologlistfilenameprefix,scratchdir)
        finshedjobs,errorjobs=poltype.WaitForTermination(fulljobtooutputlog)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(fulljobtooutputlog,True)
    os.chdir('..')

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




def get_torlist(poltype,mol):
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
    rotbndlist = {}
    iterbond = openbabel.OBMolBondIter(mol)
    v1 = valence.Valence(poltype.versionnum,poltype.logfname,poltype.dontfrag,poltype.isfragjob,poltype.allownonaromaticringscanning)
    v1.setidxtoclass(poltype.idxtosymclass)
    v1.torguess(mol,False,[])
    missed_torsions = v1.get_mt()
    for bond in iterbond:
        skiptorsion=True
        # is the bond rotatable
        t2 = bond.GetBeginAtom()
        t3 = bond.GetEndAtom()
        t2idx=t2.GetIdx()
        t3idx=t3.GetIdx()
        
        t2val=t2.GetValence()
        t3val=t3.GetValence()
        ringbond=False
        if t2.IsInRing()==True and t3.IsInRing()==True:
            ringbond=True

        arobond=False
        if t2.IsAromatic()==True and t3.IsAromatic()==True:
            arobond=True
        if arobond==True:
            continue
        if ringbond==True and poltype.allownonaromaticringscanning==False:
            continue
        if ((bond.IsRotor()) or [t2.GetIdx(),t3.GetIdx()] in poltype.fitrotbndslist or [t3.GetIdx(),t2.GetIdx()] in poltype.fitrotbndslist or [t2.GetIdx(),t3.GetIdx()] in poltype.onlyrotbndslist or [t3.GetIdx(),t2.GetIdx()] in poltype.onlyrotbndslist or (poltype.rotalltors and t2val>=2 and t3val>=2)):
            t1,t4 = find_tor_restraint_idx(poltype,mol,t2,t3)
            # is the torsion in toromitlist
            if(torfit.sorttorsion(poltype,[poltype.idxtosymclass[t1.GetIdx()],poltype.idxtosymclass[t2.GetIdx()],poltype.idxtosymclass[t3.GetIdx()],poltype.idxtosymclass[t4.GetIdx()]]) in missed_torsions) and len(poltype.onlyrotbndslist)==0:
                skiptorsion = False
            if [t2.GetIdx(),t3.GetIdx()] in poltype.fitrotbndslist or [t3.GetIdx(),t2.GetIdx()] in poltype.fitrotbndslist or [t2.GetIdx(),t3.GetIdx()] in poltype.onlyrotbndslist or [t3.GetIdx(),t2.GetIdx()] in poltype.onlyrotbndslist:
                skiptorsion = False # override previous conditions if in list
            if poltype.rotalltors==True:
                skiptorsion=False

            unq=get_uniq_rotbnd(poltype,t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx())
            rotbndkey = '%d %d' % (unq[1],unq[2])
            if (not skiptorsion):
                # store the torsion and temporary torsion value found by openbabel in torlist
                tor = mol.GetTorsion(t1,t2,t3,t4)
                if not skiptorsion:
                    torlist.append(unq)
                # store torsion in rotbndlist
                rotbndlist[rotbndkey] = []
                rotbndlist[rotbndkey].append(unq)
                # write out rotatable bond to log
                poltype.logfh.write('Rotatable bond found about %s\n' %
                str(rotbndlist[rotbndkey][0]))


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

                            if not skiptorsion:    
                                rotbndlist[rotbndkey].append(get_uniq_rotbnd(poltype,iaa.GetIdx(),t2.GetIdx(),t3.GetIdx(),iaa2.GetIdx()))

            else:
                continue
            
    return (torlist ,rotbndlist)



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
            if tpdkey not in torsionlist:
                torsionlist.append(tpdkey)
        
        if poltype.tordatapointsnum==None:
            prmnum=len(poltype.nfoldlist)*(len(torsionlist))+1
            torset=tuple([maintor])
            poltype.torsionsettonumptsneeded[torset]=prmnum
            ang=round(360/(prmnum)) # offset parameter is the +1
            if ang> 30:
                ang=30
           
        else:
            ang=round(360/poltype.tordatapointsnum)
        ratio=round(360/ang)
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
    maxb1class = max(b1nbridx,key=lambda x: poltype.idxtosymclass[x])

    iteratomatom = openbabel.OBAtomAtomIter(b2)
    b2nbridx = list(map(lambda x: x.GetIdx(), iteratomatom))
    del b2nbridx[b2nbridx.index(b1idx)]    # Remove b1 from list
    assert(b2nbridx is not [])
    maxb2class = max(b2nbridx, key= lambda x:poltype.idxtosymclass[x])


    t1 = mol.GetAtom(maxb1class)
    t4 = mol.GetAtom(maxb2class)


    return t1,t4

def ConvertTinktoXYZ(poltype,filename):
    temp=open(os.getcwd()+r'/'+filename,'r')
    tempwrite=open(os.getcwd()+r'/'+filename.replace('.xyz_2','_xyzformat.xyz'),'w')
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


def CreatePsi4TorOPTInputFile(poltype,torset,phaseangles,optmol,torxyzfname,variabletorlist):
    inputname = '%s-opt-'%(poltype.molecprefix)
    for i in range(len(torset)):
        tor=torset[i]
        a,b,c,d=tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        phaseangle=phaseangles[i]
        inputname+='%d-%d-%03d'%(b,c,round((torang+phaseangle)%360))
        inputname+='_'
    inputname=inputname[:-1]
    inputname+='.psi4'     
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (optmol.GetTotalCharge(),1))
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
            for resttors in poltype.rotbndlist[key]:
                rta,rtb,rtc,rtd = resttors
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
        for resttors in torsions:
            rta,rtb,rtc,rtd = resttors
            rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
            if resttors not in restlist and resttors not in variabletorlist:
                if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
                   (optmol.GetAtom(rtd).GetAtomicNum() != 1):
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
        temp.write('  basis '+poltype.toroptbasisset+'\n')
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
        temp.write('}'+'\n')

    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('for _ in range(1):'+'\n')
    temp.write('  try:'+'\n')
    if poltype.toroptpcm==True:
        temp.write('    set opt_coordinates cartesian'+'\n')
    temp.write("    optimize('%s/%s')" % (poltype.toroptmethod.lower(),poltype.toroptbasisset)+'\n')
    if poltype.freq:
        temp.write('    scf_e,scf_wfn=freq(%s/%s,return_wfn=True)'%(poltype.toroptmethod.lower(),poltype.toroptbasisset)+'\n')
    temp.write('    break'+'\n')
    temp.write('  except OptimizationConvergenceError:'+'\n')
    temp.write('    try:'+'\n')
    temp.write('      set opt_coordinates cartesian'+'\n')
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

def gen_torcomfile (poltype,comfname,numproc,maxmem,maxdisk,prevstruct,xyzf):
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
    optimizeoptlist = ["maxcycle=%s"%(poltype.optmaxcycle),'Loose']

    optstr=opt.gen_opt_str(poltype,optimizeoptlist)

    if ('-opt-' in comfname):
        if poltype.toroptpcm==True:
            operationstr = "%s %s/%s MaxDisk=%s SCRF=(PCM)\n" % (optstr,poltype.toroptmethod,poltype.toroptbasisset, maxdisk)
        else:
            operationstr = "%s %s/%s MaxDisk=%s\n" % (optstr,poltype.toroptmethod,poltype.toroptbasisset, maxdisk)
        commentstr = poltype.molecprefix + " Rotatable Bond Optimization on " + gethostname()
    else:
#        operationstr = "#m06L/%s SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=%s\n" % (torspbasisset, maxdisk)
        if poltype.torsppcm==True:
            operationstr = "#P %s/%s SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=%s SCRF=(PCM) Pop=NBORead\n" % (poltype.torspmethod,poltype.torspbasisset, maxdisk)
        else:       
            operationstr = "#P %s/%s SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=%s Pop=NBORead\n" % (poltype.torspmethod,poltype.torspbasisset, maxdisk)

        commentstr = poltype.molecprefix + " Rotatable Bond SP Calculation on " + gethostname()   

    bset=re.search('6-31\S+',operationstr)
    if ('I ' in poltype.mol.GetSpacedFormula()):
        operationstr=re.sub(r'6-31\S+',r'Gen',operationstr)
    tmpfh.write(operationstr)
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d\n' % (prevstruct.GetTotalCharge(), prevstruct.GetTotalSpinMultiplicity()))
    iteratom = openbabel.OBMolAtomIter(prevstruct)
    etab = openbabel.OBElementTable()
    if os.path.isfile(xyzf):
        xyzstr = open(xyzf,'r')
        xyzstrl = xyzstr.readlines()
        i = 0
        for atm in iteratom:
            i = i + 1
            ln = xyzstrl[i]
            tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), float(ln.split()[2]),float(ln.split()[3]),float(ln.split()[4])))
        tmpfh.write('\n')
        tmpfh.close()
        xyzstr.close()
    else:
        for atm in iteratom:
            tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), atm.x(),atm.y(),atm.z()))
        tmpfh.write('\n')

        tmpfh.write('$nbo bndidx $end'+'\n')
        tmpfh.write('\n')

        tmpfh.close()

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
            outfh.write(str(ele))
        outfh.write("\n")


def CreatePsi4TorESPInputFile(poltype,finalstruct,torxyzfname,optmol,torset,phaseangles,makecube=None):
    inputname= '%s-sp-'%(poltype.molecprefix)
    for i in range(len(torset)):
        tor=torset[i]
        phaseangle=phaseangles[i]
        a,b,c,d=tor[0:4]
        torang = optmol.GetTorsion(a,b,c,d)
        inputname+='%d-%d-%03d' % (b,c,round((torang+phaseangle)%360))
        inputname+='_'
    inputname=inputname[:-1]
    inputname+='.psi4'    
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (optmol.GetTotalCharge(), 1))
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
    if poltype.torsppcm==True:
        temp.write('set {'+'\n')
        temp.write(' e_convergence 10 '+'\n')
        temp.write(' d_convergence 10 '+'\n')
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
    else:
        temp.write('set {'+'\n')
        temp.write(' e_convergence 10 '+'\n')
        temp.write(' d_convergence 10 '+'\n')
        temp.write('}'+'\n')
    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write("E, wfn = energy('%s/%s',return_wfn=True)" % (poltype.torspmethod.lower(),poltype.torspbasisset)+'\n')
    temp.write('oeprop(wfn,"WIBERG_LOWDIN_INDICES")'+'\n')

    temp.write('clean()'+'\n')
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname

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
    keystokeep=[]
    for i in range(len(poltype.torlist)):
        tor=poltype.torlist[i]
        rotbndkey=str(tor[1])+' '+str(tor[2])
        if rotbndkey not in keystokeep:
            keystokeep.append(rotbndkey)
    keystoremove=[]
    for key,value in poltype.rotbndlist.items():
        if key not in keystokeep:
            keystoremove.append(key) 
    for key in keystoremove:
        del poltype.rotbndlist[key]
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

