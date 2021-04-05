import os
import sys
from socket import gethostname
import openbabel
import re
import time
import apicall as call
import shlex
import numpy as np

def CreatePsi4OPTInputFile(poltype,comfilecoords,comfilename,mol,modred,restraints,skipscferror,chg):
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    inputname=comfilename.replace('.com','.psi4')
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    if chg==None:
        temp.write('%d %d\n' % (mol.GetTotalCharge(), 1))
    else:
        temp.write('%d %d\n' % (chg, 1))

    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')
    if poltype.optpcm==True:
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

    if restraints!=None:
        space='  '
        bondres=[restraints[0]]
        string='frozen_distance'
        temp.write('set optking{'+'\n')
        temp.write('  '+string+' '+'='+' '+'('+'"'+'\n')
        for res in bondres:
            res=[str(i) for i in res]
            resstring=' '.join(res)+'\n'
            temp.write('   '+resstring)
        temp.write('  "'+')'+'\n')
        temp.write('}'+'\n')
       
        anglerestraints=restraints[1:]
        string='frozen_bend'
        temp.write('set optking{'+'\n')
        temp.write('  '+string+' '+'='+' '+'('+'"'+'\n')
        for res in anglerestraints:
            res=[str(i) for i in res]
            resstring=' '.join(res)+'\n'
            temp.write('   '+resstring)
        temp.write('  "'+')'+'\n')
        temp.write('}'+'\n')
   



    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('for _ in range(1):'+'\n')
    temp.write('  try:'+'\n')
    if poltype.optpcm==True:
        temp.write('    set opt_coordinates cartesian'+'\n')
    spacedformulastr=mol.GetSpacedFormula()
    if ('I ' in spacedformulastr):
        temp.write('    basis {'+'\n')
        temp.write('    ['+' '+poltype.optbasissetfile+' '+poltype.iodineoptbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.optbasissetfile,poltype.iodineoptbasissetfile)
        temp.write('    }'+'\n')
        temp.write("    optimize('%s')" % (poltype.optmethod.lower())+'\n')

    else:
        if modred==False:
            temp.write('    set opt_coordinates both'+'\n')
        temp.write("    optimize('%s/%s')" % (poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    if poltype.freq:
        temp.write('    scf_e,scf_wfn=freq("%s/%s",return_wfn=True)'%(poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    temp.write('    break'+'\n')
    temp.write('  except OptimizationConvergenceError:'+'\n')
    temp.write('    break'+'\n')
    if skipscferror==True:
        temp.write('  except SCFConvergenceError:'+'\n')
        temp.write('    pass'+'\n')
   
    temp.write('  else:'+'\n')
    temp.write('    try:'+'\n')
    temp.write('      set opt_coordinates cartesian'+'\n')
    if ('I ' in spacedformulastr):
        temp.write("      optimize('%s')" % (poltype.optmethod.lower())+'\n')
    else:
        temp.write("      optimize('%s/%s')" % (poltype.optmethod.lower(),poltype.optbasisset)+'\n')

    if poltype.freq:
        temp.write('      scf_e,scf_wfn=freq("%s/%s",return_wfn=True)'%(poltype.optmethod.lower(),poltype.optbasisset)+'\n')
    temp.write('      break'+'\n')
    temp.write('    except OptimizationConvergenceError:'+'\n')
    temp.write('      '+'pass'+'\n')

    temp.write('clean()'+'\n')
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname

def ReadInBasisSet(poltype,tmpfh,normalelementbasissetfile,otherelementbasissetfile):
    newtemp=open(poltype.basissetpath+normalelementbasissetfile,'r')
    results=newtemp.readlines()
    newtemp.close()
    for line in results:
        if '!' not in line:
            tmpfh.write('    '+line)


    newtemp=open(poltype.basissetpath+otherelementbasissetfile,'r')
    results=newtemp.readlines()
    newtemp.close()
    for line in results:
        if '!' not in line:
            tmpfh.write('    '+line)
    return tmpfh



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


def GrabFinalXYZStructure(poltype,logname,filename,mol):
    checkifpsi4=CheckIfPsi4Log(poltype,logname)
    if checkifpsi4==True:
        temp=open(logname,'r')
        results=temp.readlines()
        temp.close()
        temp=open(filename,'w')
        temp.write(str(mol.NumAtoms())+'\n')
        temp.write('\n')
        finalmarker=False
        lengthchange=None
        lastsuccessidx=None
        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Successfully symmetrized geometry' in line:
                lastsuccessidx=lineidx
        if lastsuccessidx==None: # sometimes it doesnt print this but converges? 
            lastsuccessidx=len(results)-1
        for lineidx in range(len(results)):
            line=results[lineidx]
            try:
                if lineidx<lastsuccessidx:
                    if 'Geometry (in Angstrom)' in line:
                        lastidx=lineidx
            except:
                if lineidx<lastsuccessidx:
                    if 'Geometry (in Angstrom)' in line:
                        lastidx=lineidx

        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Geometry (in Angstrom)' in line and lineidx==lastidx:
                finalmarker=True
            if finalmarker==True and lineidx>lastidx:    
                linesplit=line.split()
                if (len(linesplit)!=4 and len(linesplit)!=5) and lengthchange==False:
                    lengthchange=True
                    break
                foundfloat=bool(re.search(r'\d', line))
                if (len(linesplit)==4 or len(linesplit)==5) and foundfloat==True and 'point' not in line:
                    temp.write(' '.join(line.lstrip().split()[:3+1])+'\n')
                    lengthchange=False
        temp.close()
    elif checkifpsi4==False:
        obConversion = openbabel.OBConversion()
        tempmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(logname)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(tempmol, logname)
        obConversion.SetOutFormat('xyz')
        obConversion.WriteFile(tempmol, filename)

def gen_optcomfile(poltype,comfname,numproc,maxmem,maxdisk,chkname,molecule,modred=True):
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

    spacedformulastr=molecule.GetSpacedFormula()
    if ('I ' in spacedformulastr):
        modred=False
    if modred==True:
        optimizeoptlist = ["ModRedundant","maxcycles=%s"%(poltype.optmaxcycle),'Loose']
    else:
        optimizeoptlist = ["Cartesian","maxcycles=%s"%(poltype.optmaxcycle)]

    if restraintlist:
        optimizeoptlist.insert(0,poltype.gausoptcoords)
    optstr=gen_opt_str(poltype,optimizeoptlist)
    if ('I ' in spacedformulastr):
        prevoptbasisset=poltype.optbasisset
        poltype.optbasisset='gen'
    if poltype.freq==True:
        if poltype.optpcm==True:
            optstring= "%s %s/%s freq SCRF=(PCM)" % (optstr,poltype.optmethod,poltype.optbasisset)
        else:
            optstring= "%s %s/%s freq" % (optstr,poltype.optmethod,poltype.optbasisset)
    else:
        if poltype.optpcm==True:
            optstring= "%s %s/%s SCRF=(PCM)" % (optstr,poltype.optmethod,poltype.optbasisset)
        else:
            optstring= "%s %s/%s" % (optstr,poltype.optmethod,poltype.optbasisset)
    if ('I ' in spacedformulastr):
        optstring+=' pseudo=read'
    string=' MaxDisk=%s \n'%(maxdisk)
    optstring+=string
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
    
    if ('I ' in spacedformulastr):
        formulalist=spacedformulastr.lstrip().rstrip().split()
        elementtobasissetlines=GenerateElementToBasisSetLines(poltype,poltype.basissetpath+poltype.optbasissetfile)
        for element,basissetlines in elementtobasissetlines.items():
            if element in spacedformulastr:
                for line in basissetlines: 
                    tmpfh.write(line)

        temp=open(poltype.basissetpath+poltype.iodineoptbasissetfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '!' not in line:
                tmpfh.write(line)


        
        tmpfh.write('\n')
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

    tmpfh.write('%RWF=' + poltype.scrtmpdirgau + '/,' + maxdisk + '\n')
    tmpfh.write("%Nosave\n")
    tmpfh.write("%Chk=" + os.path.splitext(comfname)[0] + ".chk\n")
    tmpfh.write("%Mem=" + maxmem + "\n")
    tmpfh.write("%Nproc=" + str(numproc) + "\n")
    tmpfh.close()

def CompareBondLengths(poltype,inioptmol,optmol):
    isnear=True
    tol=.1

    for inib in openbabel.OBMolBondIter(inioptmol):
        beg = inib.GetBeginAtomIdx()
        end = inib.GetEndAtomIdx()
        b=optmol.GetBond(beg,end)
        if b==None:
            isnear=False
            break
        iniblength=inib.GetLength()
        blength=b.GetLength()
        
        diff=np.abs(iniblength-blength)/iniblength
        if diff>=tol:
            isnear=False
    return isnear

def CheckBondConnectivity(poltype,mol,optmol,raiseerror=False):
    atomitermol=openbabel.OBMolAtomIter(mol)
    atomiteroptmol=openbabel.OBMolAtomIter(optmol)
    issame=True

    firstatomnum=mol.NumAtoms()
    secondatomnum=optmol.NumAtoms()
    if firstatomnum!=secondatomnum:
        issame=False
        return issame

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
            issame=False
            if raiseerror==True:
                RaiseConnectivityError(poltype,diff,idxset)
    return issame

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

    cmd = poltype.minimizeexe+' -k '+poltype.tmpkeyfile+' '+poltype.tmpxyzfile+' 0.1 > Minimized_final.out'
    poltype.call_subsystem(cmd, True)


def GeometryOptimization(poltype,mol,checkbonds=True,modred=True,restraints=None,skipscferror=False,charge=None,skiperrors=False): # specify charge instead of reading from mol if charge!=None
    poltype.WriteToLog("NEED QM Density Matrix: Executing Gaussian Opt and SP")
    
    if (poltype.use_gaus==True or poltype.use_gausoptonly==True): # try to use gaussian for opt
        term,error=poltype.CheckNormalTermination(poltype.logoptfname)
        if not term:
            mystruct = load_structfile(poltype,poltype.molstructfname)
            gen_optcomfile(poltype,poltype.comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkoptfname,mol,modred)
            cmdstr = 'cd '+shlex.quote(os.getcwd())+' && '+'GAUSS_SCRDIR=' + poltype.scrtmpdirgau + ' ' + poltype.gausexe + " " + poltype.comoptfname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logoptfname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdirgau
            jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix 
            if os.path.isfile(poltype.chkoptfname) and os.path.isfile(poltype.logoptfname):
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
        gen_optcomfile(poltype,poltype.comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkoptfname,mol,modred)
        poltype.WriteToLog("Calling: " + "Psi4 Optimization")
        term,error=poltype.CheckNormalTermination(poltype.logoptfname,None,skiperrors)
        modred=False

        inputname,outputname=CreatePsi4OPTInputFile(poltype,poltype.comoptfname,poltype.comoptfname,mol,modred,restraints,skipscferror,charge)
        if term==False:
            cmdstr='cd '+shlex.quote(os.getcwd())+' && '+'psi4 '+inputname+' '+poltype.logoptfname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logoptfname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdirpsi4
            jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix
            if os.path.isfile(poltype.logoptfname):
                os.remove(poltype.logoptfname)
            if poltype.externalapi==None:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,False)
            else:
                if len(jobtooutputlog.keys())!=0:
                    call.CallExternalAPI(poltype,jobtolog,jobtologlistfilepathprefix,scratchdir)
                finishedjobs,errorjobs=poltype.WaitForTermination(jobtooutputlog)
        term,error=poltype.CheckNormalTermination(poltype.logoptfname,None,skiperrors) # now grabs final structure when finished with QM if using Psi4
        if error and term==False and skiperrors==False:
            poltype.RaiseOutputFileError(poltype.logoptfname) 
        GrabFinalXYZStructure(poltype,poltype.logoptfname,poltype.logoptfname.replace('.log','.xyz'),mol)
        optmol =  load_structfile(poltype,poltype.logoptfname.replace('.log','.xyz'))
        optmol=rebuild_bonds(poltype,optmol,mol)


    GrabFinalXYZStructure(poltype,poltype.logoptfname,poltype.logoptfname.replace('.log','.xyz'),mol)
    if checkbonds==True:
        issame=CheckBondConnectivity(poltype,mol,optmol,raiseerror=True)
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
