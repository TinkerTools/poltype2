from . import optimization as opt
from . import electrostaticpotential as esp
import os
import sys
from socket import gethostname
import re

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
    # Create a *.grid file which is an input file for Gaussian CUBEGEN
    if not os.path.isfile(poltype.espgrdfname):
        gengridcmd = poltype.potentialexe + " 1 " + poltype.xyzfname+' -k '+poltype.keyfname
        poltype.call_subsystem(gengridcmd,True)
        
    #    shutil.move(xyzoutfile,espgrdfname)
    # Run CUBEGEN
    if poltype.use_psi4 or poltype.use_psi4SPonly:
        temp=open('grid_esp.dat','r')
        results=temp.readlines()
        temp.close()
        Vvals=[]
        for line in results:
            Vvals.append(line.split()[0])
            
        poltype.WriteToLog("Calling: " + "Generating CUBE File from PSI4")
        with open('grid.dat', 'r') as fp:
            gridpts = fp.readlines()
        if len(gridpts) != len(Vvals):
            raise Exception('Dimension error in potential calculation!')
        # Generate a "cube" file.  I have no idea what the format should be (it's not a
        # regular cube file) so I reverse engineered one by looking at TINKER source.
        with open(poltype.qmespfname, 'w') as fp:
            fp.write(" %s/%s potential calculation\n\n" % (poltype.espmethod,poltype.espbasisset))
            fp.write("%5d\n%5d\n\n\n" % (0,len(Vvals)))
            for xyz,v in zip(gridpts, Vvals):
                fp.write("%s %s\n" % (xyz.rstrip(), v))
        
    elif not os.path.isfile(poltype.qmespfname):
        fckfname = poltype.fckespfname
        if not poltype.espfit:
            fckfname = poltype.fckdmafname

        if not os.path.isfile(fckfname):
            fckfname = os.path.splitext(fckfname)[0]

        assert os.path.isfile(fckfname), "Error: " + fckfname + " does not exist."
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
        

def CreatePsi4ESPInputFile(poltype,comfilecoords,mol,molecprefix,a,b,c,d,torang,phaseangle,maxdisk,maxmem,numproc,makecube=None):
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    comfilename= '%s-sp-%d-%d-%d-%d-%03d.com' % (molecprefix,a,b,c,d,round((torang+phaseangle)%360))
    inputname=comfilename.replace('.com','_psi4.dat')
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (mol.GetTotalCharge(), 1))
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')
    if torsppcm==True:
        temp.write('set {'+'\n')
        temp.write(' basis '+poltype.espbasisset.lower()+'\n')
        temp.write(' e_convergence 10 '+'\n')
        temp.write(' d_convergence 10 '+'\n')
        temp.write(' scf_type pk'+'\n')
        temp.write(' pcm true'+'\n')
        temp.write('  pcm_scf_type total '+'\n')
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
    temp.write('memory '+maxmem+'\n')
    temp.write('set_num_threads(%s)'%(numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scratchdir)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write("E, wfn = energy('%s/%s',return_wfn=True)" % (poltype.espmethod.lower(),poltype.espbasisset.lower())+'\n')
    if makecube==True:
       temp.write('oeprop(wfn,"GRID_ESP")'+'\n')
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
    temp.write('%d %d\n' % (mol.GetTotalCharge(), 1))
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')
    if torsppcm==True:
        temp.write('set {'+'\n')
        temp.write(' basis '+poltype.dmabasisset.lower()+'\n')
        temp.write(' e_convergence 10 '+'\n')
        temp.write(' d_convergence 10 '+'\n')
        temp.write(' scf_type pk'+'\n')
        temp.write(' pcm true'+'\n')
        temp.write('  pcm_scf_type total '+'\n')
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
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scratchdir)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write("E, wfn = energy('%s/%s',return_wfn=True)" % (poltype.dmamethod.lower(),poltype.dmabasisset.lower())+'\n')
    temp.write('gdma(wfn)'+'\n')
    temp.write('fchk(wfn, "%s.fchk")'%(comfilename.replace('.com',''))+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    return inputname

def GrabFinalPsi4Energy(poltype,logname):
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
    temp=open('RMSPD.txt','r')
    for line in temp.readlines():
        if 'Root Mean Square Potential Difference :' in line:
            RMSPD=line.split(':')[1].strip()
    temp.close()
    if float(RMSPD)>poltype.maxRMSPD:
        print('Warning: RMSPD of QM and MM optimized structures is high, RMSPD = ',RMSPD)
        poltype.WriteToLog('Warning: RMSPD of QM and MM optimized structures is high, RMSPD = '+ RMSPD+' Tolerance is '+str(poltype.maxRMSPD)+' kcal/mol ') # now report all torsions as bad

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
    optlogfname = os.path.splitext(comfname)[0] + ".log"
    title = "\"" + poltype.molecprefix + " Gaussian SP Calculation on " + gethostname() + "\""
    cmdstr = poltype.babelexe + " --title " + title + " -i g03 " + poltype.gausoptfname + " " + tailfname
    poltype.call_subsystem(cmdstr)

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
        
        opstr="#P %s/%s Sp Density=%s SCF=Save Guess=Huckel MaxDisk=%s\n" % (poltype.espmethod,poltype.espbasisset, densitystring,maxdisk)


    bset=re.search('(?i)(6-31|aug-cc)\S+',opstr)
    if ('I ' in mol.GetSpacedFormula()):
        opstr=re.sub(r'(?i)(6-31|aug-cc)\S+',r'Gen',opstr)
    tmpfh.write(opstr)
    tmpfh.close()
    cmdstr = "tail -n +2 " + tailfname + " | head -n 4 >> " + comfname
    os.system(cmdstr)
    cmdstr = "tail -n +6 " + tailfname + " | sed -e's/ .*//' > tmp1.txt"
    os.system(cmdstr)
    cmdstr = "grep -A " + str(mol.NumAtoms()+4) + " 'Standard orientation' " + poltype.gausoptfname + " | tail -n " + str(mol.NumAtoms()) + " | sed -e's/^.* 0 //' > tmp2.txt"
    os.system(cmdstr)
    cmdstr = "paste tmp1.txt tmp2.txt >> " + comfname
    os.system(cmdstr)
    cmdstr='rm '+'tmp1.txt'
    os.system(cmdstr)
    cmdstr='rm '+'tmp2.txt'
    os.system(cmdstr)
  

def ElectrostaticPotentialFitting(poltype):
    optmpolecmd = poltype.potentialexe + " 6 " + poltype.xyzoutfile + " -k " + poltype.key2fname + " " + poltype.qmesp2fname + " N 0.1"
    poltype.call_subsystem(optmpolecmd,True)

def ElectrostaticPotentialComparison(poltype):
    poltype.WriteToLog("")
    poltype.WriteToLog("=========================================================")
    poltype.WriteToLog("Electrostatic Potential Comparision\n")
    if not os.path.isfile('RMSPD.txt'):
        cmd=poltype.potentialexe + ' 5 ' + poltype.xyzoutfile + ' ' + '-k'+' '+ poltype.key3fname+' '+ poltype.qmesp2fname + ' N > RMSPD.txt'
        poltype.call_subsystem(cmd,True)
    CheckRMSPD(poltype)


def SPForDMA(poltype,optmol,mol):
    if poltype.use_psi4==True or poltype.use_psi4SPonly:
        if os.path.isfile(poltype.chkdmafname):
            os.remove(poltype.chkdmafname)

        gen_comfile(poltype,poltype.comdmafname.replace('.com','_temp.com'),poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkdmafname,poltype.comtmp,mol)
        poltype.WriteToLog("Calling: " + "Psi4 Gradient for DMA")
        term,error=poltype.is_qm_normal_termination(poltype.logdmafname)
        if poltype.use_psi4SPonly:
            save_structfileXYZ(optmol, 'optimized.xyz')
        inputname=CreatePsi4DMAInputFile(poltype,'optimized.xyz',poltype.comdmafname,mol)
        if term==False:
            now = time.strftime("%c",time.localtime())
            poltype.WriteToLog(" Calling: " + "Psi4 GDMA for DMA")
            cmdstr='psi4 '+inputname+' '+poltype.logdmafname
            poltype.call_subsystem(cmdstr,True)
        
    else:
        term,error=is_qm_normal_termination(poltype,poltype.logdmafname)
        if not term:
            if os.path.isfile(poltype.chkdmafname):
                os.remove(poltype.chkdmafname)
            gen_comfile(poltype,poltype.comdmafname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkdmafname,poltype.comtmp,mol)
            cmdstr = 'GAUSS_SCRDIR=' + poltype.scrtmpdir + ' ' + poltype.gausexe + " " + poltype.comdmafname
            poltype.call_subsystem(cmdstr,True)
            cmdstr = poltype.formchkexe + " " + poltype.chkdmafname
            poltype.call_subsystem(cmdstr,True)

def SPForESP(poltype,optmol,mol):
    if poltype.use_psi4 or poltype.use_psi4SPonly:
        shutil.copy(poltype.espgrdfname, 'grid.dat') # same as .grid file (x,y,z) coords
        inputname=CreatePsi4ESPInputFile(poltype,'optimized.xyz',poltype.comespfname,mol,poltype.maxdisk,poltype.maxmem,poltype.numproc,True)
        if CheckNormalTermPsi4ESP(poltype,poltype.logespfname)==False:
            poltype.WriteToLog(" Calling: " + "Psi4 Gradient for ESP")
            cmdstr='psi4 '+poltype.inputname+' '+poltype.logespfname
            poltype.call_subsystem(cmdstr,True)

    else:
        term,error=is_qm_normal_termination(poltype,poltype.logespfname)
        if poltype.espfit and not term:
            if os.path.isfile(poltype.chkespfname):
                os.remove(poltype.chkespfname)
            gen_comfile(poltype,poltype.comespfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkespfname,poltype.comtmp,mol)
            cmdstr = 'GAUSS_SCRDIR=' + poltype.scrtmpdir + ' ' + poltype.gausexe + " " + poltype.comespfname
            poltype.call_subsystem(cmdstr,True)
            cmdstr = poltype.formchkexe + " " + poltype.chkespfname
            poltype.call_subsystem(cmdstr,True)

def is_qm_normal_termination(poltype,logfname): # needs to handle error checking now
    """
    Intent: Checks the *.log file for normal termination
    """
    error=False
    term=False
    if os.path.isfile(logfname):
        if 'psi4' in logfname:
            for line in open(logfname):
                if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line:
                    opt.GrabFinalPsi4XYZStructure(poltype,logfname,logfname.replace('.log','_opt.xyz'))
                    term=True
        else:
            for line in open(logfname):
                if "Normal termination" in line:
                    term=True
                if ('error' in line or 'Error' in line or 'ERROR' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation' in line) and 'DIIS' not in line:
                    error=True

    return term,error


def NormalTerm(poltype,logfname):
    poltype.WriteToLog("Normal termination: %s" % logfname)


def ErrorTerm(poltype,logfname):
    poltype.WriteToLog("ERROR termination: %s" % logfname)


