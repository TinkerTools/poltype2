import os
import sys
from socket import gethostname
import openbabel
import re
import time
import apicall as call
import shlex
import numpy as np
import shutil
import torsiongenerator as torgen

def CreatePsi4OPTInputFile(poltype,comfilecoords,comfilename,mol,modred,bondanglerestraints,skipscferror,chg,torsionrestraints=[]):
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
        if len(torsionrestraints)==0:
            temp.write('  g_convergence GAU_LOOSE'+'\n')
        temp.write('  dynamic_level 1'+'\n')
        temp.write('}'+'\n')

    if bondanglerestraints!=None:
        space='  '
        bondres=[bondanglerestraints[0]]
        string='frozen_distance'
        temp.write('set optking{'+'\n')
        temp.write('  '+string+' '+'='+' '+'('+'"'+'\n')
        for res in bondres:
            res=[str(i) for i in res]
            resstring=' '.join(res)+'\n'
            temp.write('   '+resstring)
        temp.write('  "'+')'+'\n')
        temp.write('}'+'\n')
       
        anglerestraints=bondanglerestraints[1:]
        string='frozen_bend'
        temp.write('set optking{'+'\n')
        temp.write('  '+string+' '+'='+' '+'('+'"'+'\n')
        for res in anglerestraints:
            res=[str(i) for i in res]
            resstring=' '.join(res)+'\n'
            temp.write('   '+resstring)
        temp.write('  "'+')'+'\n')
        temp.write('}'+'\n')
    if len(torsionrestraints)!=0:
        temp.write('set optking { '+'\n')
        temp.write('  frozen_dihedral = ("'+'\n')
        for residx in range(len(torsionrestraints)):
            res=torsionrestraints[residx]
            rta,rtb,rtc,rtd=res[:]
            if residx>0:
                temp.write(', %d %d %d %d\n' % (rta,rtb,rtc,rtd))
            else:
                temp.write('    %d %d %d %d\n' % (rta,rtb,rtc,rtd))

        temp.write('  ")'+'\n')
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

def gen_optcomfile(poltype,comfname,numproc,maxmem,maxdisk,chkname,molecule,modred=True,torsionrestraints=[]):
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
    if len(torsionrestraints)!=0:
        tempname=comfname.replace('.com','_temp.com')
        temp=open(comfname,'r')
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
                for res in torsionrestraints:
                    rta,rtb,rtc,rtd=res[:]
                    tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))
                tmpfh.write("\n")
            else:
                tmpfh.write(line)

        tmpfh.close()
        os.remove(comfname)
        shutil.copy(tempname,comfname)





def GenerateElementToBasisSetLines(poltype,basissetfile):
    elementtobasissetlines={}
    temp=open(basissetfile,'r')
    results=temp.readlines()
    temp.close()
    lines=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx==0:  
            linesplit=line.split()
            element=linesplit[0]
            lines=[line]
            elementtobasissetlines[element]=lines
        elif lineidx>0 and  '****' in results[lineidx-1]:
            linesplit=line.split()
            element=linesplit[0]
            lines=[line]
            elementtobasissetlines[element]=lines
        else:
            lines.append(line)
            elementtobasissetlines[element]=lines


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

def AverageBondTableLength(poltype,elementsbondorder):
    elementstobondordertolength={tuple([1,1,1]):.74,tuple([9,9,1]):1.42,tuple([17,17,1]):1.99,tuple([35,35,1]):2.28,tuple([53,53,1]):2.67,tuple([1,6,1]):1.10,tuple([1,7,1]):1.00,tuple([1,8,1]):.97,tuple([1,9,1]):.92,tuple([6,6,1]):1.54,tuple([6,7,1]):1.47,tuple([6,8,1]):1.43,tuple([7,7,1]):1.45,tuple([8,8,1]):1.45,tuple([1,6,1]):1.10,tuple([6,6,2]):1.34,tuple([6,6,3]):1.20,tuple([6,7,2]):1.28,tuple([6,8,2]):1.20,tuple([6,8,3]):1.13,tuple([7,7,2]):1.23,tuple([7,7,3]):1.10,tuple([8,8,2]):1.21,tuple([1,9,1]):.92,tuple([1,17,1]):1.27,tuple([1,35,1]):1.41,tuple([1,53,1]):1.61,tuple([6,16,1]):1.82,tuple([1,6,1]):1.10,tuple([6,9,1]):1.35,tuple([6,17,1]):1.77,tuple([6,35,1]):1.94,tuple([6,53,1]):2.14}
    found=False
    length=None
    if elementsbondorder in elementstobondordertolength.keys():
        length=elementstobondordertolength[elementsbondorder]
        found=True
    elif elementsbondorder[::-1] in elementstobondordertolength.keys():
        length=elementstobondordertolength[elementsbondorder[::-1]]
        found=True
    if found==False:
        tol=.2 # extreme case if any missing above
    else:
        tol=.05 # test this may increase tolerance later
    return tol,length

    

def CompareBondLengths(poltype,inioptmol,optmol):
    isnear=True
    for inib in openbabel.OBMolBondIter(inioptmol):
        beg = inib.GetBeginAtomIdx()
        end = inib.GetEndAtomIdx()
        b=optmol.GetBond(beg,end)
        if b==None:
            isnear=False
            break
        begatom=optmol.GetAtomWithIdx(beg)
        endatom=optmol.GetAtomWithIdx(end)
        begatomicnum=begatom.GetAtomicNum()
        endatomicnum=endatom.GetAtomicNum()
        bondorder=b.GetBondOrder()
        elementsbondorder=[begatomicnum,endatomicnum,bondorder]
        tol,length=AverageBondTableLength(poltype,elementsbondorder)
        blength=b.GetLength()
        diff=np.abs(length-blength)
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

    shutil.copy(poltype.xyzoutfile,poltype.tmpxyzfile)
    shutil.copy(poltype.key5fname,poltype.tmpkeyfile)
    cmd = poltype.minimizeexe+' -k '+poltype.tmpkeyfile+' '+poltype.tmpxyzfile+' 0.1 > Minimized_final.out'
    poltype.call_subsystem(cmd, True)


def FindTorsionRestraints(poltype,mol):
    torsionrestraints=[]
    atomiter=openbabel.OBMolAtomIter(mol)
    atomnum=0
    for atom in atomiter:
        atomnum+=1
    bondnum=0
    for b in openbabel.OBMolBondIter(mol):
        isrot=b.IsRotor()
        if isrot==True:
            bondnum+=1
    if atomnum>=25 or bondnum>=2:
        for b in openbabel.OBMolBondIter(mol):
            isrot=b.IsRotor()
            if isrot==True:
                t2 = b.GetBeginAtom()
                t3 = b.GetEndAtom()
                t1,t4 = torgen.find_tor_restraint_idx(poltype,mol,t2,t3)
                t2idx=t2.GetIdx()
                t3idx=t3.GetIdx()
                t1idx=t1.GetIdx()
                t4idx=t4.GetIdx()
                torsionrestraints.append([t1idx,t2idx,t3idx,t4idx])


    return torsionrestraints

def GeometryOptimization(poltype,mol,checkbonds=True,modred=True,bondanglerestraints=None,skipscferror=False,charge=None,skiperrors=False,overridecheckterm=False): # specify charge instead of reading from mol if charge!=None
    poltype.WriteToLog("NEED QM Density Matrix: Executing Gaussian Opt and SP")
    if bondanglerestraints!=None: # then vdw opt
        pass
    else: # see if need to restrain torsion in extended conformation
        torsionrestraints=FindTorsionRestraints(poltype,mol)
 
    if (poltype.use_gaus==True or poltype.use_gausoptonly==True): # try to use gaussian for opt
        term,error=poltype.CheckNormalTermination(poltype.logoptfname)
        if not term or overridecheckterm==True:
            mystruct = load_structfile(poltype,poltype.molstructfname)
            gen_optcomfile(poltype,poltype.comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkoptfname,mol,modred,torsionrestraints)
            cmdstr = 'GAUSS_SCRDIR=' + poltype.scrtmpdirgau + ' ' + poltype.gausexe + " " + poltype.comoptfname
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

        gen_optcomfile(poltype,poltype.comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,poltype.chkoptfname,mol,modred,torsionrestraints)
        poltype.WriteToLog("Calling: " + "Psi4 Optimization")
        term,error=poltype.CheckNormalTermination(poltype.logoptfname,None,skiperrors)
        modred=False

        inputname,outputname=CreatePsi4OPTInputFile(poltype,poltype.comoptfname,poltype.comoptfname,mol,modred,bondanglerestraints,skipscferror,charge,torsionrestraints)
        if term==False or overridecheckterm==True:
            cmdstr='psi4 '+inputname+' '+poltype.logoptfname
            jobtooutputlog={cmdstr:os.getcwd()+r'/'+poltype.logoptfname}
            jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
            scratchdir=poltype.scrtmpdirpsi4
            jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix
            if os.path.isfile(poltype.logoptfname):
                os.remove(poltype.logoptfname)
            if poltype.externalapi==None:
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,skiperrors)
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
    return optmol,error


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
