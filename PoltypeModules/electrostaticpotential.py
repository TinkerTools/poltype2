import optimization as opt
import electrostaticpotential as esp
import torsiongenerator as torgen
import apicall as call
import os
import sys
from socket import gethostname
import re
import shutil
import time
import numpy as np
from openbabel import openbabel
import shlex
import warnings
from scipy.optimize import fmin
from PyAstronomy import pyasl

def CheckIfLogFileUsingGaussian(poltype,f):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    use_gaus=False
    temp=open(f,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Entering Gaussian System' in line:
            use_gaus=True
            break 
    return use_gaus


def gen_esp_grid(poltype,mol,gridnamelist,espnamelist,fchknamelist,cubenamelist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    """
    Intent: Find the QM Electrostatic Potential Grid which can be used for multipole fitting
    Input:
    Output: 
            *.grid: CUBEGEN input file; written out by tinker's potential utility
            *.cube: CUBEGEN output file
            *.cube_2: *.cube cleaned up and with QM potential values for each point of the grid 
    Referenced By: main
    Description:
    1. Run tinker's potential utility with option 1 which says:
       "(1) Create an Input File for Gaussian CUBEGEN"
    2. Move the output to *.grid
    3. Run Gaussian CUBEGEN. Outputs *.cube
    4. Run tinker's potential utility with option 2 on *.cube. Option 2 reads:
       "(2) Get QM Potential from a Gaussian Cube File"
       Outputs *.cube_2
    """
    potnamelist=[]
    use_gaus=CheckIfLogFileUsingGaussian(poltype,poltype.logespfname)
    if use_gaus==False:
        for i in range(len(gridnamelist)):
            gridname=gridnamelist[i]
            espname=espnamelist[i]
            fchkname=fchknamelist[i]
            cubename=cubenamelist[i]
            potfile=cubename.replace('.cube','.pot')
            Vvals,gridpts=GrabGridData(poltype,gridname,espname)

            # Generate a "cube" file.  I have no idea what the format should be (it's not a
            # regular cube file) so I reverse engineered one by looking at TINKER source.
            with open(cubename, 'w') as fp:
                fp.write(" %s/%s potential calculation\n\n" % (poltype.espmethod,poltype.espbasisset))
                fp.write("%5d\n%5d\n\n\n" % (0,len(Vvals)))
                for xyz,v in zip(gridpts, Vvals):
                    fp.write("%s %s\n" % (xyz.rstrip(), v))
            if not os.path.isfile(potfile):
                genqmpotcmd = poltype.potentialexe + " 2 " + cubename
                poltype.call_subsystem([genqmpotcmd],True)
            potnamelist.append(potfile)

    for i in range(len(gridnamelist)):
        gridname=gridnamelist[i]
        espname=espnamelist[i]
        fchkname=fchknamelist[i]
        cubename=cubenamelist[i]
        potfile=cubename.replace('.cube','.pot') 
        if not os.path.isfile(cubename):
            fckfname = fchkname
            if not poltype.espfit:
                fckfname = poltype.fckdmafname

            if not os.path.isfile(fckfname):
                fckfname = os.path.splitext(fckfname)[0]

            assert os.path.isfile(fckfname), "Error: " + fckfname + " does not exist."+' '+os.getcwd()
            if poltype.espmethod.upper() =='MP2':
                densitystring='MP2'
            else:
                densitystring='SCF'
            gencubecmd = poltype.cubegenexe + " 0 potential=%s "%(densitystring) + fckfname + " " + cubename + " -5 h < " + gridname
            poltype.call_subsystem([gencubecmd],True)
    # Run potential
            if not os.path.isfile(potfile):
                genqmpotcmd = poltype.potentialexe + " 2 " + cubename
                poltype.call_subsystem([genqmpotcmd],True)
            potnamelist.append(potfile)
    return potnamelist
       
def GrabGridData(poltype,gridname,espname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(espname,'r')
    results=temp.readlines()
    temp.close()
    Vvals=[]
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            Vvals.append(linesplit[0])
        
    poltype.WriteToLog("Calling: " + "Generating CUBE File from PSI4")
    with open(gridname, 'r') as fp:
        gridpts = fp.readlines()
    return Vvals,gridpts


      
def CreatePsi4ESPInputFile(poltype,comfilecoords,comfilename,mol,maxdisk,maxmem,numproc,charge,makecube=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    inputname=comfilename.replace('.com','.psi4')
    temp=open(inputname,'w')
    if not poltype.sameleveldmaesp:
      temp.write('molecule { '+'\n')
      temp.write('%d %d\n' % (charge, mol.GetTotalSpinMultiplicity()))
      for lineidx in range(len(results)):
          line=results[lineidx]
          linesplit=line.split()
          if len(linesplit)==4 and '#' not in line:
              temp.write(line)
      temp.write('}'+'\n')
    temp.write('memory '+maxmem+'\n')
    temp.write('set_num_threads(%s)'%(numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('set maxiter '+str(poltype.scfmaxiter)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write('set PROPERTIES_ORIGIN ["COM"]'+'\n')
    temp.write("set cubeprop_tasks ['esp']"+'\n')
    if not poltype.sameleveldmaesp:
      temp.write('set basis %s '%(poltype.espbasisset)+'\n')
      if poltype.allowradicals==True:
          temp.write('set reference uhf '+'\n')
          temp.write("G, wfn = gradient('%s', return_wfn=True)" % (poltype.espmethod.lower())+'\n')
      else:
          spacedformulastr=mol.GetSpacedFormula()
          if ('I ' in spacedformulastr):
              temp.write('basis {'+'\n')
              temp.write('assign '+poltype.espbasisset+'\n')
              temp.write('assign I '+poltype.iodineespbasisset+'\n')
              temp.write('}'+'\n')
          if makecube==True:
              temp.write("E, wfn = properties('%s',properties=['dipole','GRID_ESP','WIBERG_LOWDIN_INDICES','MULLIKEN_CHARGES'],return_wfn=True)" % (poltype.espmethod.lower())+'\n')
          else:
              temp.write("E, wfn = properties('%s',properties=['dipole','WIBERG_LOWDIN_INDICES','MULLIKEN_CHARGES'],return_wfn=True)" % (poltype.espmethod.lower())+'\n')

      temp.write('cubeprop(wfn)'+'\n')
    temp.write("wfn = psi4.core.Wavefunction.from_file('dma.wfn.npy')\n")
    temp.write('oeprop(wfn,"GRID_ESP","WIBERG_LOWDIN_INDICES","MULLIKEN_CHARGES")\n')
    if not poltype.sameleveldmaesp:
      temp.write('fchk(wfn, "%s.fchk")'%(comfilename.replace('.com',''))+'\n')

    temp.write('clean()'+'\n')
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname


def CreatePsi4DMAInputFile(poltype,comfilecoords,comfilename,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    inputname=comfilename.replace('.com','.psi4')
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (poltype.totalcharge, mol.GetTotalSpinMultiplicity()))
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('no_reorient'+'\n')
    temp.write('}'+'\n')
    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write('set PROPERTIES_ORIGIN ["COM"]'+'\n')
    temp.write("set cubeprop_tasks ['esp']"+'\n')
    temp.write('set basis %s '%(poltype.dmabasisset)+'\n')
    if poltype.allowradicals==True:
        temp.write('set reference uhf '+'\n')

        temp.write("G, wfn = gradient('%s', return_wfn=True)" % (poltype.dmamethod.lower())+'\n')
    else:
        spacedformulastr=mol.GetSpacedFormula()
        if ('I ' in spacedformulastr):
            temp.write('basis {'+'\n')
            temp.write('assign '+poltype.dmabasisset+'\n')
            temp.write('assign I '+poltype.iodinedmabasisset+'\n')
            temp.write('}'+'\n')
            temp.write("E, wfn = properties('%s',properties=['dipole'],return_wfn=True)" % (poltype.dmamethod.lower())+'\n')
        
        else:
        
            temp.write("E, wfn = properties('%s',properties=['dipole'],return_wfn=True)" % (poltype.dmamethod.lower())+'\n')
              
    temp.write("wfn.to_file('dma.wfn')\n")
    temp.write('cubeprop(wfn)'+'\n')
    temp.write('fchk(wfn, "%s.fchk")'%(comfilename.replace('.com',''))+'\n')
    if poltype.sameleveldmaesp:
        temp.write('fchk(wfn, "%s.fchk")'%(comfilename.replace('.com','').replace('-dma', '-esp'))+'\n')
      
    header=poltype.molstructfname.split('.')[0]
    temp.write('clean()'+'\n')
    temp.close()
    return inputname






def GrabFinalPsi4Energy(poltype,logname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    energy=None
    temp=open(logname,'r')
    results=temp.readlines()
    temp.close()
    foundfinalenergies=False
    for line in results:
        linesplit=line.split()
        if 'Energies <====================' in line or '=> Energetics <=' in line:
            foundfinalenergies=True
        if foundfinalenergies==True and 'Total Energy' in line and 'SCS' not in line:
            energy=float(linesplit[3])
        if foundfinalenergies==True and 'SCS Total Energy' in line:
            energy=float(linesplit[4])
    return energy
            
def CheckRMSPD(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    rmspdexists=False
    poltype.failedrmspd=False
    if os.path.isfile('RMSPD.txt'):
        temp=open('RMSPD.txt','r')
        for line in temp.readlines():
            if 'Root Mean Square Potential Difference :' in line:
                RMSPD=line.split(':')[1].strip()
                rmspdexists=True
        temp.close()
        if rmspdexists==True:
            relRMSPD=ComputeRelativeRMSPD(poltype)
            reltol=3
            string1='Warning: RMSPD of QM and MM optimized structures is high, RMSPD = '+ RMSPD+' Absolute tolerance is '+str(poltype.maxRMSPD)+' kcal/mol '+ 'and relative RMSPD='+str(relRMSPD)+'% relative tolerance is '+str(reltol)+'%'
            string2='RMSPD = '+ RMSPD+' Absolute tolerance is '+str(poltype.maxRMSPD)+' kcal/mol '+ 'and relative RMSPD='+str(relRMSPD)+'% relative tolerance is '+str(reltol)+'%'

            if float(RMSPD)>poltype.maxRMSPD and relRMSPD>reltol:
                poltype.failedrmspd=True
                poltype.WriteToLog(string1)
                warnings.warn(string1)
                if poltype.issane==True and poltype.skipespfiterror==False:
                    raise ValueError(os.getcwd()+string1)
            else:
                poltype.WriteToLog(string2)

    return rmspdexists

def IsFloat(poltype,string):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    isfloat=False
    try:
        float(string)
        isfloat=True
    except:
        pass
    return isfloat

def ComputeRelativeRMSPD(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    target=[]
    model=[]
    if os.path.isfile('RMSPD.txt'):
       temp=open('RMSPD.txt','r')
       for line in temp.readlines():
           linesplit=line.split()
           if len(linesplit)==6:
               allfloats=True
               for e in linesplit:
                   isfloat=IsFloat(poltype,e)
                   if isfloat==False:
                       allfloats=False
               if allfloats==True:
                   targetvalue=float(linesplit[-2])
                   modelvalue=float(linesplit[-3])
                   target.append(targetvalue)
                   model.append(modelvalue)
    model=np.array(model)
    model=np.add(1,model)
    target=np.array(target)
    target=np.add(1,target)
    def RMSDRel(c):
        return np.sqrt(np.mean(np.square(np.add(np.divide(model-target,target),c))))
    resultRel=fmin(RMSDRel,.5)
    minRMSDRel=round(RMSDRel(resultRel[0])*100,2)
    return minRMSDRel


def gen_comfile(poltype,comfname,numproc,maxmem,maxdisk,chkname,tailfname,mol):
    """
    Intent: Create *.com file for qm dma and sp
    Input:
        comfname: com file name
        numproc: number of processors
        maxmem: max memory size
        chkname: chk file name
        tailfname: tmp com file 
        mol: OBMol object
    Output:
        *.com is written
    Referenced By: run_gaussian
    Description: -
    """

    opt.write_com_header(poltype,comfname,chkname,maxdisk,maxmem,numproc)
    tmpfh = open(comfname, "a")
    #NOTE: Need to pass parameter to specify basis set
    if ('dma' in comfname):
        if ('I ' in poltype.mol.GetSpacedFormula()):
            iodinebasissetfile=poltype.iodinedmabasissetfile 
            basissetfile=poltype.dmabasissetfile 
            #poltype.dmamethod='wB97XD'
        if poltype.dmamethod.upper() == 'MP2':
            densitystring='MP2'
        else:
            densitystring='SCF'


        opstr="#P %s/%s Sp Density=%s" % (poltype.dmamethod,poltype.dmabasisset, densitystring)
    elif ('pop' in comfname):
        opstr="#P HF/%s MaxDisk=%s Pop=SaveMixed" % (poltype.popbasisset)
    else:
        if ('I ' in poltype.mol.GetSpacedFormula()):
            poltype.espbasisset='gen'
            iodinebasissetfile=poltype.iodineespbasissetfile 
            basissetfile=poltype.espbasissetfile 
            #poltype.espmethod='wB97XD'


        if poltype.espmethod.upper() =='MP2':
            densitystring='MP2'
        else:
            densitystring='SCF'
        
        if poltype.dontfrag==False: 
            opstr="#P %s/%s Sp Density=%s SCF=Save Pop=NBORead" % (poltype.espmethod,poltype.espbasisset, densitystring)
        else:
            opstr="#P %s/%s Sp Density=%s SCF=Save" % (poltype.espmethod,poltype.espbasisset, densitystring)


    if ('I ' in poltype.mol.GetSpacedFormula()):
        opstr+=' pseudo=read'
    string=' MaxDisk=%s \n'%(maxdisk)
    opstr+=string


    tmpfh.write(opstr)
    commentstr = poltype.molecprefix + " Gaussian SP Calculation on " + gethostname()
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d\n' % (mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity()))
    tmpfh.close()

    iteratombab = openbabel.OBMolAtomIter(mol)
    tmpfh = open(comfname, "a")
    for atm in iteratombab:
        symb=poltype.indextoatomicsymbol[atm.GetIdx()]
        tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (symb, atm.x(), atm.y(), atm.z()))


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
    if not 'dma' in comfname and poltype.dontfrag==False:  
        tmpfh.write('\n')
        tmpfh.write('$nbo bndidx $end'+'\n')
        tmpfh.write('\n')
    else:
        tmpfh.write('\n')
        tmpfh.write('\n')


    tmpfh.close()


def GenerateElementToBasisSetLines(poltype,basissetfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
 
def CombineFiles(poltype,namelist,name):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    totalresults=[]
    for fname in namelist:
        temp=open(fname,'r')
        results=temp.readlines()
        temp.close()
        totalresults.extend(results)
    temp=open(name,'w')
    for line in totalresults:
        temp.write(line)
    temp.close()
    return name



def ElectrostaticPotentialFitting(poltype,xyzfnamelist,keyfnamelist,potnamelist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    combinedxyz=CombineFiles(poltype,xyzfnamelist,'combined.xyz')
    combinedpot=CombineFiles(poltype,potnamelist,'combined.pot')
    optmpolecmd = poltype.potentialexe + " 6 " + combinedxyz + " -k " + poltype.key2fnamefromavg + " " + combinedpot + " N "+str(poltype.espgrad)
    if poltype.deletedfiles==True:
        poltype.call_subsystem([optmpolecmd],True)
    else:
        poltype.call_subsystem([optmpolecmd],True)
    shutil.copy(combinedxyz.replace('.xyz','.key'),poltype.key3fnamefrompot)
    return combinedxyz,combinedpot



def ElectrostaticPotentialComparison(poltype,combinedxyz,combinedpot):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    if poltype.deleteallnonqmfiles==True: 
        poltype.WriteToLog("")
        poltype.WriteToLog("=========================================================")
        poltype.WriteToLog("Electrostatic Potential Comparison\n")
        cmd=poltype.potentialexe + ' 5 ' + combinedxyz + ' ' + '-k'+' '+ poltype.key3fname+' '+ combinedpot + ' N > RMSPD.txt'
        poltype.call_subsystem([cmd],True)
        rmspdexists=CheckRMSPD(poltype)

def SPForDMA(poltype,optmol,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    if poltype.use_gaus==False or poltype.use_gausoptonly==True:

        gen_comfile(poltype,poltype.comdmafname.replace('.com','_temp.com'),poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkdmafname,poltype.comtmp,optmol)
        term,error=poltype.CheckNormalTermination(poltype.logdmafname,errormessages=None,skiperrors=True)
        inputname=CreatePsi4DMAInputFile(poltype,poltype.logoptfname.replace('.log','.xyz'),poltype.comdmafname,mol)
        if term==False:
            poltype.WriteToLog("Gas Phase Single Point for GDMA")
            cmdstr=' '.join(['psi4 ',poltype.psi4_args,inputname,poltype.logdmafname])
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logdmafname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdirpsi4
            inputfilepath=os.path.join(os.getcwd(),inputname)
            jobtoinputfilepaths={cmdstr:[inputfilepath]}
            jobtooutputfiles={cmdstr:[poltype.logdmafname,poltype.fckdmafname]}
            jobtoabsolutebinpath={cmdstr:poltype.which('psi4')}
            jobtologlistfilenameprefix=os.getcwd()+r'/'+'dma_jobtolog_'+poltype.molecprefix


            if poltype.externalapi!=None:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog,False)
            else:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False)

            term,error=poltype.CheckNormalTermination(poltype.logdmafname)
            if error:
                poltype.RaiseOutputFileError(poltype.logdmafname) 

        
    else:
        term,error=poltype.CheckNormalTermination(poltype.logdmafname,errormessages=None,skiperrors=True)
        if not term:
            poltype.WriteToLog("Gas Phase Single Point for GDMA")
            gen_comfile(poltype,poltype.comdmafname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkdmafname,poltype.comtmp,optmol)
            cmdstr = poltype.gausexe + " " + poltype.comdmafname
            
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logdmafname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdirgau
            jobtologlistfilenameprefix=os.getcwd()+r'/'+'dma_jobtolog_'+poltype.molecprefix
            inputfilepath=os.path.join(os.getcwd(),poltype.comdmafname)
            jobtoinputfilepaths={cmdstr:[inputfilepath]}
            jobtooutputfiles={cmdstr:[poltype.logdmafname,poltype.chkdmafname]}
            jobtoabsolutebinpath={cmdstr:poltype.which(poltype.gausexe)}

            if poltype.externalapi!=None:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog,False)
            else:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False)

            poltype.call_subsystem([cmdstr],True)
            cmdstr = poltype.formchkexe + " " + poltype.chkdmafname
            poltype.call_subsystem([cmdstr],True)
            term,error=poltype.CheckNormalTermination(poltype.logdmafname)
            if error:
                poltype.RaiseOutputFileError(poltype.logdmafname) 


def SPForESP(poltype,optmolist,molist,xyzfnamelist,keyfnamelist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gridnamelist=[]
    espnamelist=[] 
    fchknamelist=[]
    cubenamelist=[]
    for i in range(len(optmolist)):
        optmol=optmolist[i]
        mol=molist[i]
        xyzfname=xyzfnamelist[i]
        keyfname=keyfnamelist[i]
        if i==0:
            suffix=''
           
        else:
            suffix='_'+str(i+1)
        optsuffix='_'+str(i+1)

        logoptfname=poltype.logoptfname.replace('_1.log',optsuffix+'.log')
        comespfname=poltype.comespfname.replace('.com',suffix+'.com')
        fckespfname=poltype.fckespfname.replace('.fchk',suffix+'.fchk')
        chkespfname=poltype.chkespfname.replace('.chk',suffix+'.chk')
        logespfname=poltype.logespfname.replace('.log',suffix+'.log')
        cubename=poltype.qmespfname.replace('.cube',suffix+'.cube')
        cubenamelist.append(cubename)
        fchknamelist.append(fckespfname)
        espgrdfname=xyzfname.replace('.xyz','.grid')
        gridnamelist.append(espgrdfname)
        espname=espgrdfname.replace('.grid','esp.dat')
        espnamelist.append(espname)
        if not os.path.isfile(espgrdfname):
            gengridcmd = poltype.potentialexe + " 1 " + xyzfname+' -k '+keyfname
            poltype.call_subsystem([gengridcmd],True)
        if poltype.use_gaus==False or poltype.use_gausoptonly==True:
            shutil.copy(espgrdfname, 'grid.dat') 
            inputname,outputname=CreatePsi4ESPInputFile(poltype,logoptfname.replace('.log','.xyz'),comespfname,mol,poltype.maxdisk,poltype.maxmem,poltype.numproc,poltype.totalcharge,True)
            term,error=poltype.CheckNormalTermination(outputname,errormessages=None,skiperrors=True)
            if term==False:
                poltype.WriteToLog("Gas Phase High Level Single Point for Multipole Refinement")
                cmdstr=' '.join(['psi4 ',poltype.psi4_args,inputname,outputname])
                jobtooutputlog={cmdstr:os.getcwd()+r'/'+outputname}
                jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
                scratchdir=poltype.scrtmpdirpsi4
                jobtologlistfilenameprefix=os.getcwd()+r'/'+'esp_jobtolog_'+poltype.molecprefix 
                inputfilepath=os.path.join(os.getcwd(),inputname)
                anotherinputfile=os.path.join(os.getcwd(),'grid.dat')
                jobtoinputfilepaths={cmdstr:[inputfilepath,anotherinputfile]}
                jobtooutputfiles={cmdstr:[outputname,fckespfname,'grid_esp.dat']}
                jobtoabsolutebinpath={cmdstr:poltype.which('psi4')}


                if poltype.externalapi!=None:
                    if len(jobtooutputlog.keys())!=0:
                        call.CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
                    finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog,False)
                else:
                    finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False,wait=True)

                term,error=poltype.CheckNormalTermination(outputname)
                if error:
                    poltype.RaiseOutputFileError(outputname)
                 
                shutil.copy('grid_esp.dat',espname) 


        else:
            term,error=poltype.CheckNormalTermination(logespfname,errormessages=None,skiperrors=True)
            if poltype.espfit and not term:
                poltype.WriteToLog("Gas Phase High Level Single Point for Multipole Refinement")
                gen_comfile(poltype,comespfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,chkespfname,poltype.comtmp,optmol)
                cmdstr = poltype.gausexe + " " + comespfname
                if poltype.sameleveldmaesp:
                  comdmafname = comespfname.replace('-esp', '-dma')
                  logdmafname = logespfname.replace('-esp', '-dma')
                  fckdmafname = fckespfname.replace('-esp', '-dma')
                  chkdmafname = chkespfname.replace('-esp', '-dma')
                  cmdstr = f"ln -sf {comdmafname}  {comespfname}; ln -sf {logdmafname}  {logespfname}; ln -sf {fckdmafname} {fckespfname}; ln -sf {chkdmafname} {chkespfname}"
                jobtooutputlog={cmdstr:os.getcwd()+r'/'+logespfname}
                jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
                scratchdir=poltype.scrtmpdirgau
                jobtologlistfilenameprefix=os.getcwd()+r'/'+'esp_jobtolog_'+poltype.molecprefix
                inputfilepath=os.path.join(os.getcwd(),comespfname)
                jobtoinputfilepaths={cmdstr:[inputfilepath]}
                jobtooutputfiles={cmdstr:[logespfname,chkespfname]}
                jobtoabsolutebinpath={cmdstr:poltype.which(poltype.gausexe)}


                if poltype.externalapi!=None:
                    if len(jobtooutputlog.keys())!=0:
                        call.CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
                    finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog,False)
                else:
                    finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False,wait=True)

                cmdstr = poltype.formchkexe + " " + chkespfname
                poltype.call_subsystem([cmdstr],True)
                term,error=poltype.CheckNormalTermination(logespfname)
                if error:
                    poltype.RaiseOutputFileError(logespfname)

    sys.stdout.flush()
    return gridnamelist,espnamelist,fchknamelist,cubenamelist


def CheckIfLogFileIsGaussian(poltype,logname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gaus=False
    temp=open(logname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Gaussian System' in line:
            gaus=True
            break

    return gaus


def GrabQMDipoles(poltype,optmol,logname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    gaus=False
    gaus=CheckIfLogFileIsGaussian(poltype,logname)
    if gaus==False:
        temp=open(logname,'r')
        results=temp.readlines()
        temp.close()
        firstway=False
        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Dipole Moment: [D]' in line:
                firstway=True
                nextline=results[lineidx+1]
                nextlinesplit=nextline.split()
                dipole=np.array([float(nextlinesplit[1]),float(nextlinesplit[3]),float(nextlinesplit[5])])
        if firstway==False:
            dipole=[]
            for lineidx in range(len(results)):
                line=results[lineidx]
                linesplit=line.split()
                if 'Dipole' in line:
                    dipole.append(float(linesplit[-1])*2.5417464519) # lee pings branch affected output somehow sometimes?
            dipole=np.array(dipole)
                
    else:
        temp=open(logname,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            linesplit=line.split()
            
            if 'Tot' in line and 'X' in line and 'Y' in line and 'Z' in line:
                dipole=np.array([float(linesplit[1]),float(linesplit[3]),float(linesplit[5])])
                dipole=ConvertDipoleToCOMFrame(poltype,dipole,optmol)
    dipole=np.array([round(dipole[0],3),round(dipole[1],3),round(dipole[2],3)])
    return dipole


def CheckDipoleMoments(poltype,optmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    logname=poltype.logespfname
    if poltype.sameleveldmaesp:
      logname = poltype.logdmafname
    if not os.path.exists(logname):
        return
    dipole=GrabQMDipoles(poltype,optmol,logname)
    qmdipole=round(np.linalg.norm(dipole),3)
    while not os.path.exists(poltype.tmpkeyfile):
        time.sleep(1)
    torgen.RemoveStringFromKeyfile(poltype,poltype.tmpkeyfile,'solvate')
    cmd=poltype.analyzeexe + ' ' + poltype.xyzoutfile+' '+'-k'+' '+poltype.tmpkeyfile + ' em | grep -A11 Charge'+'>'+'MMDipole.txt'
    poltype.call_subsystem([cmd],True)

    while not os.path.isfile('MMDipole.txt'):
        time.sleep(1)
    temp=open('MMDipole.txt','r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Dipole Moment' in line:
            linesplit=line.split()
            mmdipole=float(linesplit[-2])
            diff=qmdipole-mmdipole
            if qmdipole!=0:
                ratio=np.abs(diff/qmdipole)
            else:
                ratio=0

            if ratio>poltype.dipoletol and poltype.suppressdipoleerr==False and np.abs(diff)>poltype.absdipoletol:
                if poltype.esprestweight==.1:
                    string='Relative error of '+str(ratio)+' for QMDipole '+str(qmdipole)+' and '+str(mmdipole)+' for MMDipole '+'is bigger than '+str(poltype.dipoletol)+' '+os.getcwd()
                    raise ValueError(string)

            else:
                 string='Relative error of '+str(ratio)+' for QMDipole '+str(qmdipole)+' and '+str(mmdipole)+' for MMDipole '+' tolerance = '+str(poltype.dipoletol)+' '+os.getcwd()
                 poltype.WriteToLog(string)


def ConvertDipoleToCOMFrame(poltype,dipole,optmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    nucsum = 0.0
    masssum=0.0
    atomiter=openbabel.OBMolAtomIter(optmol)
    for atom in atomiter:
        atomicnum=atom.GetAtomicNum()
        nucsum += atomicnum
        masssum+=atom.GetAtomicMass()
    #COC
    coc_x = 0.0
    coc_y = 0.0
    coc_z = 0.0

    com_x = 0.0
    com_y = 0.0
    com_z = 0.0

    newatomiter=openbabel.OBMolAtomIter(optmol)
    for atom in newatomiter:
        atomicnum=atom.GetAtomicNum()
        atomicmass=atom.GetAtomicMass()
        x=atom.GetX()
        y=atom.GetY()
        z=atom.GetZ()
        coc_x = coc_x + x*atomicnum/nucsum
        coc_y = coc_y + y*atomicnum/nucsum
        coc_z = coc_z + z*atomicnum/nucsum
        com_x = com_x + x*atomicmass/masssum
        com_y = com_y + y*atomicmass/masssum
        com_z = com_z + z*atomicmass/masssum

    COM=np.array([com_x,com_y,com_z])
    COC=np.array([coc_x,coc_y,coc_z])
    debye2eA= 0.20819434
    # 1 debye = 0.393456 ebohr 
    bohr_e2debye = 0.393456
    bohr2A = 0.529177210
    debye2eA= 0.20819434
    vec_gaus = poltype.totalcharge*(COC-COM)/bohr2A/bohr_e2debye
    dipole=dipole-vec_gaus
    return dipole
