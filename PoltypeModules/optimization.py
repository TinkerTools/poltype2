import os
import sys
from socket import gethostname
import openbabel
import re
import time
import apicall as call

def CreatePsi4OPTInputFile(poltype,comfilecoords,comfilename,mol):
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
    if poltype.optpcm==True:
        temp.write('set {'+'\n')
        temp.write('  basis '+poltype.optbasisset+'\n')
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
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scratchdir)+'\n')
    temp.write('for _ in range(1):'+'\n')
    temp.write('  try:'+'\n')
    if poltype.optpcm==True:
        temp.write('    set opt_coordinates cartesian'+'\n')
    temp.write("    optimize('%s/%s')" % (poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    if poltype.freq:
        temp.write('    scf_e,scf_wfn=freq(%s/%s,return_wfn=True)'%(poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    temp.write('    break'+'\n')
    temp.write('  except OptimizationConvergenceError:'+'\n')
    temp.write('    try:'+'\n')
    temp.write('      set opt_coordinates cartesian'+'\n')
    temp.write("      optimize('%s/%s')" % (poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    if poltype.freq:
        temp.write('      scf_e,scf_wfn=freq(%s/%s,return_wfn=True)'%(poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    temp.write('      break'+'\n')
    temp.write('    except OptimizationConvergenceError:'+'\n')
    temp.write('      '+'pass'+'\n')

    temp.write('clean()'+'\n')
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname

def NumberInLine(poltype,line):
    numinline=False
    linesplit=line.split()
    
    for e in linesplit:
        try:
            float(e)
            numinline=True
        except:
            continue
            
    return numinline


def GrabFinalXYZStructure(poltype,logname,filename):
    if poltype.use_gaus==False and poltype.use_gausoptonly==False:
        temp=open(logname,'r')
        results=temp.readlines()
        temp.close()
        temp=open(filename,'w')
        temp.write(str(poltype.mol.NumAtoms())+'\n')
        temp.write('\n')
        finalmarker=False
        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Cartesian Geometry' in line:
                lastidx=lineidx
        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Cartesian Geometry' in line and lineidx==lastidx:
                finalmarker=True
                lengthchange=False
            if finalmarker==True and lineidx>lastidx:    
                linesplit=line.split()
                if len(linesplit)!=4:
                    lengthchange=True
                if len(linesplit)==4 and bool(re.search(r'\d', line))==True and 'point' not in line and lengthchange==False:
                    temp.write(line.lstrip())
        temp.close()
    elif poltype.use_gaus==True or poltype.use_gausoptonly==True:
        obConversion = openbabel.OBConversion()
        tempmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(logname)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(tempmol, logname)
        obConversion.SetOutFormat('xyz')
        obConversion.WriteFile(tempmol, filename)


def gen_optcomfile(poltype,comfname,numproc,maxmem,maxdisk,chkname,molecule):
    """
    Intent: Create *.com file for qm opt
    Input:
        comfname: com file name
        numproc: number of processors
        maxmem: max memory size
        chkname: chk file name
        mol: OBMol object
    Output:
        *opt-*.com is written
    Referenced By: run_gaussian
    Description: -
    """
    restraintlist = []
    write_com_header(poltype,comfname,chkname,maxdisk,maxmem,numproc)
    tmpfh = open(comfname, "a")
    optimizeoptlist = ["maxcycle=%s"%(poltype.optmaxcycle),'Loose']
    if restraintlist:
        optimizeoptlist.insert(0,poltype.gausoptcoords)
    optstr=gen_opt_str(poltype,optimizeoptlist)
    if ('I ' in molecule.GetSpacedFormula()):
        optstring="%s HF/Gen freq Guess=INDO MaxDisk=%s\n" % (optstr,maxdisk)
    else:
        if poltype.freq==True:
            if poltype.optpcm==True:
                optstring= "%s %s/%s freq Guess=INDO MaxDisk=%s SCRF=(PCM)\n" % (optstr,poltype.optmethod,poltype.optbasisset,maxdisk)
            else:
                optstring= "%s %s/%s freq Guess=INDO MaxDisk=%s\n" % (optstr,poltype.optmethod,poltype.optbasisset,maxdisk)
        else:
            if poltype.optpcm==True:
                optstring= "%s %s/%s Guess=INDO MaxDisk=%s SCRF=(PCM)\n" % (optstr,poltype.optmethod,poltype.optbasisset,maxdisk)
            else:
                optstring= "%s %s/%s Guess=INDO MaxDisk=%s\n" % (optstr,poltype.optmethod,poltype.optbasisset,maxdisk)
    tmpfh.write(optstring)
    commentstr = poltype.molecprefix + " Gaussian OPT Calculation on " + gethostname()
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d\n' % (molecule.GetTotalCharge(), molecule.GetTotalSpinMultiplicity()))
    tmpfh.close()

    iteratombab = openbabel.OBMolAtomIter(molecule)
    tmpfh = open(comfname, "a")
    etab = openbabel.OBElementTable()
    for atm in iteratombab:
        tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), atm.x(), atm.y(), atm.z()))
    tmpfh.write('\n')

    tmpfh.close()
    
def gen_opt_str(poltype,optimizeoptlist):
    optstr = "#P opt"
    if optimizeoptlist:
        optstr += "=(" + ','.join(optimizeoptlist) + ")"
    return optstr

def write_com_header(poltype,comfname,chkfname,maxdisk,maxmem,numproc):
    """
    Intent: Add header to *.com file
    Referenced By: gen_optcomfile
    """
    tmpfh = open(comfname, "w")
    assert tmpfh, "Cannot create file: " + comfname+' '+os.getcwd()

    tmpfh.write('%RWF=' + poltype.scrtmpdir + '/,' + maxdisk + '\n')
    tmpfh.write("%Nosave\n")
    tmpfh.write("%Chk=" + os.path.splitext(comfname)[0] + ".chk\n")
    tmpfh.write("%Mem=" + maxmem + "\n")
    tmpfh.write("%Nproc=" + str(numproc) + "\n")
    tmpfh.close()

def CheckBondConnectivity(poltype,mol,optmol):
    atomitermol=openbabel.OBMolAtomIter(mol)
    atomiteroptmol=openbabel.OBMolAtomIter(optmol)
    for atm in atomitermol:
       
        atmidxmol=atm.GetIdx()
        atmoptmol=optmol.GetAtom(atmidxmol)
        atmidxoptmol=atmoptmol.GetIdx()
        atmneighbidxlist=[]
        iteratomatommol = openbabel.OBAtomAtomIter(atm)
        for newatm in iteratomatommol:
            atmneighbidxlist.append(newatm.GetIdx())
        atmneighbidxlistoptmol=[]
        iteratomatomoptmol = openbabel.OBAtomAtomIter(atmoptmol)
        for newatm in iteratomatomoptmol:
            atmneighbidxlistoptmol.append(newatm.GetIdx())
        if set(atmneighbidxlist)!=set(atmneighbidxlistoptmol):
            if len(set(atmneighbidxlist))>len(set(atmneighbidxlistoptmol)):
                diff=set(atmneighbidxlist)-set(atmneighbidxlistoptmol)
                idxset='mol'
            else:
                diff=set(atmneighbidxlistoptmol)-set(atmneighbidxlist)
                idxset='optmol'
            RaiseConnectivityError(poltype,diff,idxset)

def RaiseConnectivityError(poltype,diff,idxset):
    print('Error! The bond connectivity before and after structure optimization is different')
    poltype.WriteToLog('Error! The bond connectivity before and after structure optimization is different')
    for atmidx in diff:
        print('The atom index '+str(atmidx)+' from structure '+idxset+' does not have the same connectivity before and after structure optimization')
        poltype.WriteToLog('The atom index '+str(atmidx)+' from structure '+idxset+' does not have the same connectivity before and after structure optimization')
    sys.exit() 

def gen_superposeinfile(poltype):
    """
    Intent: Initialize superpose input file (for tinker's superpose) 
    """
    f = open(poltype.superposeinfile, 'w')
    f.write('\n\n\n\n\n')

# Create file to specify groups of atoms based on molecular symmetry


def CheckRMSD(poltype):
    for line in open(poltype.superposeinfile,'r'):
        if os.stat(poltype.superposeinfile).st_size > 5:
            if 'Root Mean' in line:
                RMSD=''
                for e in line:
                    if e.isdigit() or e=='.':
                        RMSD+=e
            if float(RMSD)>poltype.maxRMSD:
                poltype.WriteToLog('Warning: RMSD of QM and MM optimized structures is high, RMSD = '+ RMSD+' Tolerance is '+str(poltype.maxRMSD)+' kcal/mol ')

                raise ValueError(os.getcwd()+' '+'RMSD of QM and MM optimized structures is high, RMSD = '+str(RMSD))

def StructureMinimization(poltype):
     poltype.WriteToLog("")
     poltype.WriteToLog("=========================================================")
     poltype.WriteToLog("Minimizing structure\n")
    
     cmd='cp ' + poltype.xyzoutfile + ' ' + poltype.tmpxyzfile
     poltype.call_subsystem(cmd)
     cmd='cp ' + poltype.key5fname + ' ' + poltype.tmpkeyfile
     poltype.call_subsystem(cmd)
     cmd = poltype.minimizeexe+' -k '+poltype.tmpkeyfile+' '+poltype.tmpxyzfile+' 0.1 > Minimizedttt.out'
     poltype.call_subsystem(cmd,True)


def GeometryOptimization(poltype,mol):
    poltype.WriteToLog("NEED QM Density Matrix: Executing Gaussian Opt and SP")
    OBOPTname=poltype.molstructfname.replace('.sdf','_OBOPT.pdb')

    if not os.path.isfile(OBOPTname):
        cmdstr=poltype.obminimizeexe+' -ff MMFF94 '+poltype.molstructfname +' '+ '>'+ OBOPTname
        poltype.call_subsystem(cmdstr,True)
    obConversion = openbabel.OBConversion()
    OBOPTmol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(OBOPTname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(OBOPTmol,OBOPTname)
    charge=poltype.totalcharge
    OBOPTmol.SetTotalCharge(charge) # for some reason obminimize does not print charge in output PDB
    

    if poltype.use_gaus==False and poltype.use_gausoptonly==False:
        if not os.path.exists(poltype.scratchdir):
            mkdirstr='mkdir '+poltype.scratchdir
            poltype.call_subsystem(mkdirstr,True)
    else:
        if not os.path.exists(poltype.scrtmpdir):
            mkdirstr='mkdir '+poltype.scrtmpdir
            poltype.call_subsystem(mkdirstr,True)

    if (poltype.use_gaus==True or poltype.use_gausoptonly==True): # try to use gaussian for opt
        term,error=poltype.CheckNormalTermination(poltype.logoptfname)
        if not term:
            mystruct = load_structfile(poltype,poltype.molstructfname)
            gen_optcomfile(poltype,poltype.comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkoptfname,OBOPTmol)
            cmdstr = 'cd '+os.getcwd()+' && '+'GAUSS_SCRDIR=' + poltype.scrtmpdir + ' ' + poltype.gausexe + " " + poltype.comoptfname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logoptfname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdir
            jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix 
            if os.path.isfile(poltype.chkoptfname):
                os.remove(poltype.logoptfname) # if chk point exists just remove logfile, there could be error in it and we dont want WaitForTermination to catch error before job is resubmitted by daemon 
            if poltype.externalapi==None:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,True) # have to skip errors because setting optmaxcycle to low number in gaussian causes it to crash
            else:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilepathprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)

            cmdstr = poltype.formchkexe + " " + poltype.chkoptfname
            poltype.call_subsystem(cmdstr)
        term,error=poltype.CheckNormalTermination(poltype.logoptfname)
        if error and term==False:
            poltype.RaiseOutputFileError(poltype.logoptfname) 
        optmol =  load_structfile(poltype,poltype.logoptfname)
        optmol=rebuild_bonds(poltype,optmol,mol)
        
           
    else:
        gen_optcomfile(poltype,poltype.comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkoptfname,OBOPTmol)
        poltype.WriteToLog("Calling: " + "Psi4 Optimization")
        term,error=poltype.CheckNormalTermination(poltype.logoptfname)
        inputname,outputname=CreatePsi4OPTInputFile(poltype,poltype.comoptfname,poltype.comoptfname,OBOPTmol)
        if term==False:
            cmdstr='cd '+os.getcwd()+' && '+'psi4 '+inputname+' '+poltype.logoptfname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logoptfname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scratchdir
            jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix
            if os.path.isfile(poltype.logoptfname):
                os.remove(poltype.logoptfname)
            if poltype.externalapi==None:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False)
            else:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilepathprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)
            
        term,error=poltype.CheckNormalTermination(poltype.logoptfname) # now grabs final structure when finished with QM if using Psi4
        if error and term==False:
            poltype.RaiseOutputFileError(poltype.logoptfname) 
        GrabFinalXYZStructure(poltype,poltype.logoptfname,poltype.logoptfname.replace('.log','.xyz'))
        optmol =  load_structfile(poltype,poltype.logoptfname.replace('.log','.xyz'))
        optmol=rebuild_bonds(poltype,optmol,mol)


    GrabFinalXYZStructure(poltype,poltype.logoptfname,poltype.logoptfname.replace('.log','.xyz'))
    CheckBondConnectivity(poltype,mol,optmol)
    return optmol


def load_structfile(poltype,structfname):
    """
    Intent: load 'structfname' as an OBMol structure
    Input:
        structfname: structure file name
    Output:
        tmpmol: OBMol object with information loaded from structfname
    Referenced By: run_gaussian, tor_opt_sp, compute_qm_tor_energy, compute_mm_tor_energy 
    Description: -
    """
    strctext = os.path.splitext(structfname)[1]
    tmpconv = openbabel.OBConversion()
    if strctext in '.fchk':
        tmpconv.SetInFormat('fchk')
    elif strctext in '.log':
        tmpconv.SetInFormat('g03')
    else:
        inFormat = openbabel.OBConversion.FormatFromExt(structfname)
        tmpconv.SetInFormat(inFormat)
    tmpmol = openbabel.OBMol()
    tmpconv.ReadFile(tmpmol, structfname)

    return tmpmol

def rebuild_bonds(poltype,newmol, refmol):
    for b in openbabel.OBMolBondIter(refmol):
        beg = b.GetBeginAtomIdx()
        end = b.GetEndAtomIdx()

        if not newmol.GetBond(beg,end):
            newmol.AddBond(beg,end, b.GetBO(), b.GetFlags())
        else:
            newb=newmol.GetBond(beg,end)
            bondorder=newb.GetBondOrder()
            newb.SetBondOrder(b.GetBO())

    return newmol


def PruneBonds(poltype,mol,bondtopology):
    molindexlist=[]
    atomitermol=openbabel.OBMolAtomIter(mol)
    for atom in atomitermol:
        molindexlist.append(atom.GetIdx())
    molidxtonewmolidx={}
    newmol=openbabel.OBMol() # define new OBMol object for the fragment
    atomlist=[] # list of atom objects from mol object
    newatomlist=[] # list of blank atom objects for fragment mol object
    count=1
    for index in molindexlist: # iterate over indexes in torsion
        atom=mol.GetAtom(index) # grab the atom object via index number
        molidx=atom.GetIdx()
        atomlist.append(atom) # append the atom object to list
        newatom=newmol.NewAtom()
        newatom=newatom.Duplicate(atom)
        newatomlist.append(newatom) # just put into blank atom objects
        molidxtonewmolidx[molidx]=count
        count+=1

    bondorderidxdic={} # key=(atom index1 of bond,atom index2 of bond), value=bondorder
    iterbond = openbabel.OBMolBondIter(mol) # iterator for all bond objects in the molecule
    for bond in iterbond:
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        aidx=a.GetIdx()
        bidx=b.GetIdx()
        if aidx in molindexlist and bidx in molindexlist: # check to make sure we want these atoms
            newaidx=molindexlist.index(aidx)+1 # new atom indexes are always in the order of atoms added to molecule via newatomlist above, +1 is because python starts at 0, atom indexes start at 1
            newbidx=molindexlist.index(bidx)+1
            bondorder=bond.GetBondOrder()
            bondorderidxdic[(newaidx,newbidx)]=bondorder

        else:
            continue

    for key in bondorderidxdic: # add back the bond between atoms in original mol object to the fragment mol object
        key=list(key)
        newaidx=key[0]
        newbidx=key[1]
        bondorder=bondorderidxdic[(newaidx,newbidx)]
        bondset=set([newaidx,newbidx])
        if bondset in bondtopology:
            newmol.AddBond(newaidx,newbidx,bondorder)

    return newmol
