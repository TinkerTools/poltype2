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
import openbabel

def gen_esp_grid(poltype,mol):
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
    if poltype.use_gaus==False or poltype.use_gausoptonly==True:
        Vvals,gridpts=GrabGridData(poltype)

        while len(gridpts) != len(Vvals):
            time.sleep(60)
            sys.stdout.flush()
            poltype.WriteToLog('Waiting to flush system stdout for grid.dat and grid_esp.dat')
            Vvals,gridpts=GrabGridData(poltype)
   
        # Generate a "cube" file.  I have no idea what the format should be (it's not a
        # regular cube file) so I reverse engineered one by looking at TINKER source.
        with open(poltype.qmespfname, 'w') as fp:
            fp.write(" %s/%s potential calculation\n\n" % (poltype.espmethod,poltype.espbasisset))
            fp.write("%5d\n%5d\n\n\n" % (0,len(Vvals)))
            for xyz,v in zip(gridpts, Vvals):
                fp.write("%s %s\n" % (xyz.rstrip(), v))
        
    if not os.path.isfile(poltype.qmespfname):
        fckfname = poltype.fckespfname
        if not poltype.espfit:
            fckfname = poltype.fckdmafname

        if not os.path.isfile(fckfname):
            fckfname = os.path.splitext(fckfname)[0]

        assert os.path.isfile(fckfname), "Error: " + fckfname + " does not exist."+' '+os.getcwd()
        if poltype.espmethod=='MP2':
            densitystring='MP2'
        else:
            densitystring='SCF'
        gencubecmd = poltype.cubegenexe + " 0 potential=%s "%(densitystring) + fckfname + " " + poltype.qmespfname + " -5 h < " + poltype.espgrdfname
        poltype.call_subsystem(gencubecmd,True)
    # Run potential
    if not os.path.isfile(poltype.qmesp2fname):
        genqmpotcmd = poltype.potentialexe + " 2 " + poltype.qmespfname
        poltype.call_subsystem(genqmpotcmd,True)
       
def GrabGridData(poltype):
    temp=open('grid_esp.dat','r')
    results=temp.readlines()
    temp.close()
    Vvals=[]
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            Vvals.append(linesplit[0])
        
    poltype.WriteToLog("Calling: " + "Generating CUBE File from PSI4")
    with open('grid.dat', 'r') as fp:
        gridpts = fp.readlines()
    return Vvals,gridpts


def CreatePsi4ESPInputFile(poltype,comfilecoords,comfilename,mol,maxdisk,maxmem,numproc,makecube=None):
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    inputname=comfilename.replace('.com','_psi4.dat')
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (poltype.totalcharge, 1))
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')
    temp.write('memory '+maxmem+'\n')
    temp.write('set_num_threads(%s)'%(numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scratchdir)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write('set PROPERTIES_ORIGIN ["COM"]'+'\n')
    temp.write("E, wfn = properties('%s/%s',properties=['dipole'],return_wfn=True)" % (poltype.espmethod.lower(),poltype.espbasisset)+'\n')
    if makecube==True:
       temp.write('oeprop(wfn,"GRID_ESP","WIBERG_LOWDIN_INDICES")'+'\n')
    else:
       temp.write('oeprop(wfn,"WIBERG_LOWDIN_INDICES")'+'\n')

    temp.write('clean()'+'\n')
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname

def CreatePsi4DMAInputFile(poltype,comfilecoords,comfilename,mol):
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    inputname=comfilename.replace('.com','_psi4.dat')
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (poltype.totalcharge, 1))
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')
    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scratchdir)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write('set PROPERTIES_ORIGIN ["COM"]'+'\n')
    temp.write("E, wfn = energy('%s/%s',properties=['dipole'],return_wfn=True)" % (poltype.dmamethod.lower(),poltype.dmabasisset)+'\n')
    temp.write('fchk(wfn, "%s.fchk")'%(comfilename.replace('.com',''))+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    return inputname

def GrabFinalPsi4Energy(poltype,logname):
    energy=None
    temp=open(logname,'r')
    results=temp.readlines()
    temp.close()
    foundfinalenergies=False
    for line in results:
        linesplit=line.split()
        if 'Energies <====================' in line:
            foundfinalenergies=True
        if foundfinalenergies==True and 'Total Energy' in line and 'SCS' not in line:
            energy=float(linesplit[3])
        if foundfinalenergies==True and 'SCS Total Energy' in line:
            energy=float(linesplit[4])
    return energy
            
def CheckRMSPD(poltype):
    rmspdexists=False
    if os.path.isfile('RMSPD.txt'):
        temp=open('RMSPD.txt','r')
        for line in temp.readlines():
            if 'Root Mean Square Potential Difference :' in line:
                RMSPD=line.split(':')[1].strip()
                rmspdexists=True
        temp.close()
        if rmspdexists==True:
            if float(RMSPD)>poltype.maxRMSPD:
                poltype.WriteToLog('Warning: RMSPD of QM and MM optimized structures is high, RMSPD = '+ RMSPD+' Tolerance is '+str(poltype.maxRMSPD)+' kcal/mol ')
            
                raise ValueError(os.getcwd()+' '+'Warning: RMSPD of QM and MM optimized structures is high, RMSPD = ',RMSPD)
    return rmspdexists

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
        if poltype.dmamethod=='MP2':
            densitystring='MP2'
        else:
            densitystring='SCF'
        opstr="#P %s/%s Sp Density=%s MaxDisk=%s\n" % (poltype.dmamethod,poltype.dmabasisset, densitystring,maxdisk)
    elif ('pop' in comfname):
        opstr="#P HF/%s MaxDisk=%s Pop=SaveMixed\n" % (poltype.popbasisset, maxdisk)
    else:
        if poltype.espmethod=='MP2':
            densitystring='MP2'
        else:
            densitystring='SCF'
        
        opstr="#P %s/%s Sp Density=%s SCF=Save Guess=Huckel MaxDisk=%s Pop=NBORead\n" % (poltype.espmethod,poltype.espbasisset, densitystring,maxdisk)


    bset=re.search('(?i)(6-31|aug-cc)\S+',opstr)
    if ('I ' in mol.GetSpacedFormula()):
        opstr=re.sub(r'(?i)(6-31|aug-cc)\S+',r'Gen',opstr)
    tmpfh.write(opstr)
    commentstr = poltype.molecprefix + " Gaussian SP Calculation on " + gethostname()
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d\n' % (mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity()))
    tmpfh.close()

    iteratombab = openbabel.OBMolAtomIter(mol)
    tmpfh = open(comfname, "a")
    etab = openbabel.OBElementTable()
    for atm in iteratombab:
        tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), atm.x(), atm.y(), atm.z()))


    tmpfh.write('\n')
    tmpfh.write('$nbo bndidx $end'+'\n')
    tmpfh.write('\n')
    tmpfh.close()

def ElectrostaticPotentialFitting(poltype):
    optmpolecmd = poltype.potentialexe + " 6 " + poltype.xyzoutfile + " -k " + poltype.key2fname + " " + poltype.qmesp2fname + " N 0.1"
    poltype.call_subsystem(optmpolecmd,True)

def ElectrostaticPotentialComparison(poltype):
    rmspdexists=CheckRMSPD(poltype)
    if rmspdexists==False:
        poltype.WriteToLog("")
        poltype.WriteToLog("=========================================================")
        poltype.WriteToLog("Electrostatic Potential Comparison\n")
        cmd=poltype.potentialexe + ' 5 ' + poltype.xyzoutfile + ' ' + '-k'+' '+ poltype.key3fname+' '+ poltype.qmesp2fname + ' N > RMSPD.txt'
        poltype.call_subsystem(cmd,True)
    rmspdexists=CheckRMSPD(poltype)

def SPForDMA(poltype,optmol,mol):
    if poltype.use_gaus==False or poltype.use_gausoptonly==True:
        if os.path.isfile(poltype.chkdmafname):
            os.remove(poltype.chkdmafname)

        gen_comfile(poltype,poltype.comdmafname.replace('.com','_temp.com'),poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkdmafname,poltype.comtmp,optmol)
        poltype.WriteToLog("Calling: " + "Psi4 Gradient for DMA")
        term,error=poltype.CheckNormalTermination(poltype.logdmafname)
        inputname=CreatePsi4DMAInputFile(poltype,poltype.logoptfname.replace('.log','.xyz'),poltype.comdmafname,mol)
        if term==False:
            cmdstr='cd '+os.getcwd()+' && '+'psi4 '+inputname+' '+poltype.logdmafname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logdmafname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scratchdir
            jobtologlistfilenameprefix=os.getcwd()+r'/'+'dma_jobtolog_'+poltype.molecprefix
            if os.path.isfile(poltype.logdmafname):
                os.remove(poltype.logdmafname) # if chk point exists just remove logfile, there could be error in it and we dont want WaitForTermination to catch error before job is resubmitted by daemon 


            if poltype.externalapi!=None:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilenameprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)
            else:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog)

            term,error=poltype.CheckNormalTermination(poltype.logdmafname)
            if error:
                poltype.RaiseOutputFileError(poltype.logdmafname) 

        
    else:
        term,error=poltype.CheckNormalTermination(poltype.logdmafname)
        if not term:
            if os.path.isfile(poltype.chkdmafname):
                os.remove(poltype.chkdmafname)
            gen_comfile(poltype,poltype.comdmafname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkdmafname,poltype.comtmp,optmol)
            cmdstr = 'cd '+os.getcwd()+' && '+'GAUSS_SCRDIR=' + poltype.scrtmpdir + ' ' + poltype.gausexe + " " + poltype.comdmafname
            
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logdmafname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdir
            jobtologlistfilenameprefix=os.getcwd()+r'/'+'dma_jobtolog_'+poltype.molecprefix
            if os.path.isfile(poltype.logdmafname):
                os.remove(poltype.logdmafname)
            if poltype.externalapi!=None:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilenameprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)
            else:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog)

            poltype.call_subsystem(cmdstr,True)
            cmdstr = poltype.formchkexe + " " + poltype.chkdmafname
            poltype.call_subsystem(cmdstr,True)
            term,error=poltype.CheckNormalTermination(poltype.logdmafname)
            if error:
                poltype.RaiseOutputFileError(poltype.logdmafname) 


def SPForESP(poltype,optmol,mol):
    if not os.path.isfile(poltype.espgrdfname):
        gengridcmd = poltype.potentialexe + " 1 " + poltype.xyzfname+' -k '+poltype.keyfname
        poltype.call_subsystem(gengridcmd,True)
    if poltype.use_gaus==False or poltype.use_gausoptonly==True:
        shutil.copy(poltype.espgrdfname, 'grid.dat') 
        inputname,outputname=CreatePsi4ESPInputFile(poltype,poltype.logoptfname.replace('.log','.xyz'),poltype.comespfname,mol,poltype.maxdisk,poltype.maxmem,poltype.numproc,True)
        term,error=poltype.CheckNormalTermination(outputname)
        if term==False:
            poltype.WriteToLog("Calling: " + "Psi4 Gradient for ESP")
            cmdstr='cd '+os.getcwd()+' && '+'psi4 '+inputname+' '+outputname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+outputname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scratchdir
            jobtologlistfilenameprefix=os.getcwd()+r'/'+'esp_jobtolog_'+poltype.molecprefix 
            if os.path.isfile(outputname):
                os.remove(outputname)

            if poltype.externalapi!=None:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilenameprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)
            else:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog)

            term,error=poltype.CheckNormalTermination(outputname)
            if error:
                poltype.RaiseOutputFileError(outputname) 


    else:
        term,error=poltype.CheckNormalTermination(poltype.logespfname)
        if poltype.espfit and not term:
            if os.path.isfile(poltype.chkespfname):
                os.remove(poltype.chkespfname)
            gen_comfile(poltype,poltype.comespfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkespfname,poltype.comtmp,optmol)
            cmdstr = 'cd '+os.getcwd()+' && '+'GAUSS_SCRDIR=' + poltype.scrtmpdir + ' ' + poltype.gausexe + " " + poltype.comespfname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logespfname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdir
            jobtologlistfilenameprefix=os.getcwd()+r'/'+'esp_jobtolog_'+poltype.molecprefix
            if os.path.isfile(poltype.logespfname):
                os.remove(poltype.logespfname)

            if poltype.externalapi!=None:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilenameprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)
            else:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog)

            cmdstr = poltype.formchkexe + " " + poltype.chkespfname
            poltype.call_subsystem(cmdstr,True)
            term,error=poltype.CheckNormalTermination(poltype.logespfname)
            if error:
                poltype.RaiseOutputFileError(poltype.logespfname)

    sys.stdout.flush()



def GrabQMDipoles(poltype,optmol,logname):
    poltype.WriteToLog("")
    poltype.WriteToLog("=========================================================")
    poltype.WriteToLog("QM Dipole moment\n")
    if poltype.use_gaus==False or poltype.use_gausoptonly==True:
        temp=open(logname,'r')
        results=temp.readlines()
        temp.close()
        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Dipole Moment: [D]' in line:
                nextline=results[lineidx+1]
                nextlinesplit=nextline.split()
                dipole=np.array([float(nextlinesplit[1]),float(nextlinesplit[3]),float(nextlinesplit[5])])

                
    else:
        grepcmd = 'grep -A7 "Dipole moment" ' + logname+'>'+'QMDipole.txt'
        poltype.call_subsystem(grepcmd,True)
        temp=open('QMDipole.txt','r')
        results=temp.readlines()
        temp.close()
        for line in results:
            linesplit=line.split()
            if 'Tot' in line:
                dipole=np.array([float(linesplit[1]),float(linesplit[3]),float(linesplit[5])])
                dipole=ConvertDipoleToCOMFrame(poltype,dipole,optmol)
    dipole=np.array([round(dipole[0],3),round(dipole[1],3),round(dipole[2],3)])
    return dipole


def CheckDipoleMoments(poltype,optmol):
    if poltype.use_gaus:
        logname=poltype.logespfname
    else:
        logname=poltype.logespfname.replace('.log','_psi4.log')
    dipole=GrabQMDipoles(poltype,optmol,logname)
    qmdipole=round(np.linalg.norm(dipole),3)
    poltype.WriteToLog('QM Dipole moment = '+str(qmdipole))    
    poltype.WriteToLog("")
    poltype.WriteToLog("=========================================================")
    poltype.WriteToLog("MM Dipole moment\n")
    torgen.RemoveStringFromKeyfile(poltype,poltype.tmpkeyfile,'solvate')
    cmd=poltype.analyzeexe + ' ' +  poltype.xyzoutfile+' '+'-k'+' '+poltype.tmpkeyfile + ' em | grep -A11 Charge'+'>'+'MMDipole.txt'
    poltype.call_subsystem(cmd,True)
    while not os.path.isfile('MMDipole.txt'):
        time.sleep(1)
    temp=open('MMDipole.txt','r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Dipole Moment' in line:
            linesplit=line.split()
            mmdipole=float(linesplit[-2])
            poltype.WriteToLog('MM Dipole moment = '+str(mmdipole))
            diff=qmdipole-mmdipole
            ratio=np.abs(diff/qmdipole)
            if ratio>poltype.dipoletol and poltype.surpressdipoleerr==False:
                raise ValueError('Relative error of '+str(ratio)+' for QMDipole '+str(qmdipole)+' and '+str(mmdipole)+' for MMDipole '+'is bigger than '+str(poltype.dipoletol)+' '+os.getcwd()) 

def ConvertDipoleToCOMFrame(poltype,dipole,optmol):
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
