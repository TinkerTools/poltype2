import os
import sys
from socket import gethostname
from openbabel import openbabel
import re
import time
import shlex
import numpy as np
import shutil
import torsiongenerator as torgen
import traceback
import warnings
from rdkit import Chem
from rdkit.Chem import rdmolfiles,RingInfo
from PyAstronomy import pyasl

def GeometryOPTWrapper(poltype,molist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    optmolist=[]
    errorlist=[]
    torsionrestraintslist=[]
    totcharge=molist[0].GetTotalCharge()
    redo=False
    for molidx in range(len(molist)):
        mol=molist[molidx]
        suf=str(molidx+1)
        try:
            optmol,error,torsionrestraints = GeometryOptimization(poltype,mol,totcharge,suffix=suf)
        except:
            redo=False
            if poltype.pcm==True:
                 # Gaussian would have been used first if was detected. 
                 poltype.pcm=False
                 poltype.optpcm=False
                 redo=True
            else:
                if poltype.optmethod=='MP2' and (poltype.use_gaus==False and poltype.use_gausoptonly==False) and poltype.foundgauss==True:
                    redo=True
                    poltype.use_gausoptonly=True

            if redo==True:
                shutil.copy(poltype.logoptfname,poltype.logoptfname.replace('.log','_failed.log'))
                optmolist,errorlist,torsionrestraintslist = GeometryOPTWrapper(poltype,molist) # recursive call allows for muliple attempts
            else:
                traceback.print_exc(file=sys.stdout)
                sys.exit()
        if redo==False:
            optmolist.append(optmol)
            errorlist.append(error)
            torsionrestraintslist.append(torsionrestraints)
    return optmolist,errorlist,torsionrestraintslist
 


def CreatePsi4OPTInputFile(poltype,comfilecoords,comfilename,mol,modred,bondanglerestraints,skipscferror,chg,loose,torsionrestraints=[]):
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
    if chg==None:
        temp.write('%d %d\n' % (mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity()))
    else:
        temp.write('%d %d\n' % (chg, mol.GetTotalSpinMultiplicity()))

    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')

    spacedformulastr=mol.GetSpacedFormula()
    if poltype.optpcm==True or (poltype.optpcm==-1 and poltype.pcm):
        temp.write('set {'+'\n')
        temp.write('  MAX_FORCE_G_CONVERGENCE 2.5e-2'+'\n')
        temp.write('  RMS_FORCE_G_CONVERGENCE 1.7e-2'+'\n')
        temp.write('  MAX_DISP_G_CONVERGENCE 1e-1'+'\n')
        temp.write('  RMS_DISP_G_CONVERGENCE 6.7e-2'+'\n')


        temp.write('  scf_type df'+'\n')
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
        if poltype.optconvergence == 'NULL':
            temp.write('  g_convergence CFOUR'+'\n')
            temp.write('  rms_force_g_convergence 1e0'+'\n')
        elif poltype.optconvergence == 'VERYLOOSE':
            temp.write('  MAX_ENERGY_G_CONVERGENCE 2e-4'+'\n')
            temp.write('  RMS_FORCE_G_CONVERGENCE 1.7e-3'+'\n')
            temp.write('  MAX_FORCE_G_CONVERGENCE 2.5e-3'+'\n')
            temp.write('  RMS_DISP_G_CONVERGENCE 1e-1'+'\n')
            temp.write('  MAX_DISP_G_CONVERGENCE 2e-1'+'\n')
        elif poltype.optconvergence == 'LOOSE':
            temp.write('  g_convergence GAU_LOOSE'+'\n')
        else:
            temp.write('  g_convergence GAU'+'\n')

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
        if len(anglerestraints)!=0:
            temp.write('set optking{'+'\n')
            temp.write('  '+string+' '+'='+' '+'('+'"'+'\n')
            for res in anglerestraints:
                res=[str(i) for i in res]
                resstring=' '.join(res)+'\n'
                temp.write('   '+resstring)
            temp.write('  "'+')'+'\n')
            temp.write('}'+'\n')
    temp.write('geometric_keywords = {'+'\n')
    temp.write(" 'coordsys' : 'tric',"+'\n')
    temp.write(" 'convergence_cmax' : 1e0,"+'\n')

    if poltype.optconvergence == 'NULL':
        poltype.WriteToLog('Warning: optconvergence=NULL is not supported by GeomeTRIC. Please set use_psi4_geometric_opt=False to use OPTKING.')
        temp.write(" 'convergence_set' : 'GAU_LOOSE',"+'\n')
        temp.write(" 'convergence_energy' : 2e-4,"+'\n')
        temp.write(" 'convergence_grms' : 1.0e0,"+'\n')
        temp.write(" 'convergence_gmax' : 1.0e0,"+'\n')
        temp.write(" 'convergence_drms' : 1.0e-1,"+'\n')
        temp.write(" 'convergence_dmax' : 2.0e-1,"+'\n')
    elif poltype.optconvergence == 'VERYLOOSE':
        temp.write(" 'convergence_set' : 'GAU_LOOSE',"+'\n')
        temp.write(" 'convergence_energy' : 2e-4,"+'\n')
        temp.write(" 'convergence_grms' : 1.7e-3,"+'\n')
        temp.write(" 'convergence_gmax' : 2.5e-3,"+'\n')
        temp.write(" 'convergence_drms' : 1.0e-1,"+'\n')
        temp.write(" 'convergence_dmax' : 2.0e-1,"+'\n')
    elif poltype.optconvergence == 'LOOSE':
        temp.write(" 'convergence_set' : 'GAU_LOOSE',"+'\n')
        temp.write(" 'convergence_energy' : 1e-4,"+'\n')

    if len(torsionrestraints)!=0:
        temp.write(" 'constraints' : {\n 'set' : [\n")
        geometric_list = []
        for res in torsionrestraints:
            angle = mol.GetTorsion(res[0], res[1], res[2], res[3]) % 360
            _str = "{'type'    : 'dihedral', 'indices' : [ %d , %d , %d , %d ], "%tuple([_-1 for _ in res[0:4]]) \
               + "'value' : %.4f } \n"%(angle)
            geometric_list.append(_str)
        temp.write("       " + "     , ".join(geometric_list) + "    ]\n  }\n")
    temp.write("}"+'\n')

    if len(torsionrestraints)!=0:
        temp.write("""set optking { \n  frozen_dihedral = ("\n""")
        temp.write(''.join(['   %d %d %d %d\n'%tuple(res[0:4]) for res in torsionrestraints]))
        temp.write("""  ")\n}\n""")

    if poltype.allowradicals==True:
        temp.write('set reference uhf '+'\n')
    
    temp.write('memory '+poltype.maxmem+'\n')
    temp.write('set_num_threads(%s)'%(poltype.numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('basis {'+'\n')
    temp.write('assign '+poltype.optbasisset+'\n')
    if ('I ' in spacedformulastr):
        temp.write('assign I '+poltype.iodineoptbasisset+'\n')
    temp.write('}'+'\n')

    temp.write("opt_finished = False\n")
    if poltype.use_psi4_geometric_opt:
        temp.write("if not opt_finished:\n")
        temp.write("    try:\n")
        temp.write("        ener, opt_hist = optimize('%s',engine='%s',optimizer_keywords=geometric_keywords, return_history=True)\n" % (poltype.optmethod.lower(),'geometric'))
        temp.write("        opt_finished = len(opt_hist['energy']) < %d\n"%(poltype.optmaxcycle))
        temp.write("    except Exception as e:\n")
        temp.write("        core.print_out('Exception: %s'%(e))\n")
    temp.write("if not opt_finished:\n")
    temp.write("    try:\n")
    temp.write("        optimize('%s')\n" % (poltype.optmethod.lower()))
    temp.write("    except:\n")
    temp.write("        set opt_coordinates both\n")
    temp.write("        optimize('%s')\n" % (poltype.optmethod.lower()))
    temp.write("    opt_finished = True\n")

    if poltype.freq:
        temp.write('    scf_e,scf_wfn=freq("%s/%s",return_wfn=True)'%(poltype.optmethod.lower(),poltype.optbasisset)+'\n')

    temp.write('clean()'+'\n')
    temp.write("assert opt_finished\n\n")
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname


def CheckIfPsi4Log(poltype,outputlog):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
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
        bohrtoang=.529177
        scale=1
        detectbohr=False # PCM module uses Bohr for some reason even though specify angstroms??
        for lineidx in range(len(results)):
            line=results[lineidx]
            if 'Successfully symmetrized geometry' in line:
                lastsuccessidx=lineidx
        if lastsuccessidx==None: # sometimes it doesnt print this but converges? 
            lastsuccessidx=len(results)-1
        for lineidx in range(len(results)):
            line=results[lineidx]
            if lineidx<lastsuccessidx:
                if 'Geometry (in Angstrom)' in line or 'Geometry (in Bohr)' in line:
                    lastidx=lineidx
                if 'Geometry (in Bohr)' in line:
                    detectbohr=True
                    scale=bohrtoang

        for lineidx in range(len(results)):
            line=results[lineidx]
            if ('Geometry (in Angstrom)' in line or 'Geometry (in Bohr)' in line) and lineidx==lastidx:
                finalmarker=True
                if 'Geometry (in Bohr)' in line:
                    detectbohr=True
                    scale=bohrtoang
            if finalmarker==True and lineidx>lastidx:    
                linesplit=line.split()
                if (len(linesplit)!=4 and len(linesplit)!=5) and lengthchange==False:
                    lengthchange=True
                    break
                foundfloat=bool(re.search(r'\d', line))
                if (len(linesplit)==4 or len(linesplit)==5) and foundfloat==True and 'point' not in line:
                    coords=line.lstrip().split()[:3+1]
                    element=coords[0].replace('BR','Br').replace('CL','Cl')
                    coords=coords[1:]
                    coords=[float(i) for i in coords]
                    coords=[i*scale for i in coords]
                    coords=[str(i) for i in coords]
                    temp.write(' '.join([element]+coords)+'\n')
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
    optconvergence = poltype.optconvergence
    if modred==True:
        optimizeoptlist = ["ModRedundant","maxcycles=%s"%(str(poltype.optmaxcycle)), optconvergence]
    else:
        optimizeoptlist = ["Cartesian","maxcycles=%s"%(str(poltype.optmaxcycle))]

    if restraintlist:
        optimizeoptlist.insert(0,poltype.gausoptcoords)
    optstr=gen_opt_str(poltype,optimizeoptlist)
    if ('I ' in spacedformulastr):
        prevoptbasisset=poltype.optbasisset
        if (poltype.use_gaus==True or poltype.use_gausoptonly==True):
            poltype.optbasisset='gen'
    if poltype.freq==True:
        if poltype.optpcm==True or (poltype.optpcm==-1 and poltype.pcm):
            optstring= "%s %s/%s freq SCRF=(PCM)" % (optstr,poltype.optmethod,poltype.optbasisset)
        else:
            optstring= "%s %s/%s freq" % (optstr,poltype.optmethod,poltype.optbasisset)
    else:
        if poltype.optpcm==True or (poltype.optpcm==-1 and poltype.pcm):
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
    an = pyasl.AtomicNo()
    for atm in iteratombab:
        tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (an.getElSymbol(atm.GetAtomicNum()), atm.x(), atm.y(), atm.z()))
    tmpfh.write('\n')

    # Add a torsion restraint for Aniline(-like) molecules
    # Chengwen Liu
    # Apr. 2025
    # this is to add a restraint for improper torsion
    # it sometimes does not lead to the geometry we wanted
    # use with caution!!

    if poltype.sp2aniline: 
      rdkitmol = Chem.MolFromMolFile(poltype.molstructfname,removeHs=False)
      pattern = Chem.MolFromSmarts('[NH2]([H])([H])[a]')
      matches = rdkitmol.GetSubstructMatches(pattern)
      for match in matches:
        n,h1,h2,x = match
        tmpfh.write(f"{n+1} {x+1} {h1+1} {h2+1}   =0.0       B\n")
        tmpfh.write(f"{n+1} {x+1} {h1+1} {h2+1}   F\n")
    
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
            else:
                tmpfh.write(line)

        tmpfh.close()
        os.remove(comfname)
        shutil.copy(tempname,comfname)





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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    optstr = "#P opt"
    if optimizeoptlist:
        optstr += "=(" + ','.join(optimizeoptlist) + ")"
    return optstr

def write_com_header(poltype,comfname,chkfname,maxdisk,maxmem,numproc):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpfh = open(comfname, "w")
    assert tmpfh, "Cannot create file: " + comfname+' '+os.getcwd()

    tmpfh.write('%RWF=' + poltype.scrtmpdirgau + '/,' + maxdisk + '\n')
    tmpfh.write("%Nosave\n")
    tmpfh.write("%Chk=" + os.path.splitext(comfname)[0] + ".chk\n")
    tmpfh.write("%Mem=" + maxmem + "\n")
    tmpfh.write("%Nproc=" + str(numproc) + "\n")
    tmpfh.close()







def gen_superposeinfile(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltype.WriteToLog("\n")
    poltype.WriteToLog("=========================================================\n")
    poltype.WriteToLog("Structure RMSD Comparison\n\n")
    cmd = poltype.superposeexe + ' ' + poltype.xyzoutfile + ' ' + poltype.tmpxyzfile + '_2'+' 1 N M N 0  > '+ poltype.superposeinfile
    poltype.call_subsystem([cmd],wait=True)


def CheckRMSD(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    RMSD=None
    for line in open(poltype.superposeinfile,'r'):
        if 'Root Mean' in line:
            RMSD=''
            for e in line:
                if e.isdigit() or e=='.':
                    RMSD+=e
    if RMSD!=None:   
        if float(RMSD)>poltype.maxRMSD:
            poltype.WriteToLog('Warning: RMSD of QM and MM optimized structures is high, RMSD = '+ RMSD+' Tolerance is '+str(poltype.maxRMSD))

            warnings.warn(os.getcwd()+' '+'RMSD of QM and MM optimized structures is high, RMSD = '+str(RMSD))
        else:
            poltype.WriteToLog('RMSD = '+ RMSD+' Tolerance is '+str(poltype.maxRMSD))

def StructureMinimization(poltype,torsionrestraints):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    poltype.WriteToLog("")
    poltype.WriteToLog("=========================================================")
    poltype.WriteToLog("Minimizing structure\n")
    AddTorsionRestraints(poltype,poltype.key7fname,torsionrestraints)
    cmd = poltype.minimizeexe+' -k '+poltype.tmpkeyfile+' '+poltype.tmpxyzfile+' 0.1 > Minimized_final.out'
    poltype.call_subsystem([cmd], True)

    torgen.RemoveStringFromKeyfile(poltype,poltype.key7fname,'restrain-torsion')
    torgen.RemoveStringFromKeyfile(poltype,poltype.tmpkeyfile,'restrain-torsion')


def AddTorsionRestraints(poltype,keyname,torsionrestraints):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    tmpfh=open(keyname,'a')
    for res in torsionrestraints:
        a,b,c,d=res[:]
        tmpfh.write('restrain-torsion %d %d %d %d %f\n' % (a,b,c,d,poltype.torsionrestraint))
    tmpfh.close()

def FindTorsionRestraints(poltype):
    """
    Intent: find non-aromatic torsion
    Input: poltype object
    Output: torsion restraint list
    Referenced By: GeometryOptimization 
    Description: Find the torsion atoms we want to restraint during geometry optimization.
    1. Only non-cyclic rotatable bond torsion
    2. Four atoms must be heavy atoms
    3. Two torsions MAX on one rotatable bond
       a. Heaviest two atoms on one side (a,b)
       b. Heaviest two atoms on the other side (c,d)
       c. Two torsions will be (a,x,y,d) and (b,x,y,c)
    """

    torsionrestraints=[]
    mol = Chem.MolFromMolFile(poltype.molstructfname,removeHs=False)

    # find all diherals in the molecule
    rotbond_torsions = {}
    pattern = Chem.MolFromSmarts('[*]~!#[*]~!#[*]~!#[*]')
    matches = mol.GetSubstructMatches(pattern, uniquify=True)
    for match in matches:
      match = list(match)
      a, x, y, b = match
      sameRing = RingInfo.AreAtomsInSameRing(mol.GetRingInfo(), x, y)
      if not sameRing:
        if x < y:
          comb = f"{x}-{y}"
          tor = match
        else:
          comb = f"{y}-{x}"
          tor = [b,y,x,a]
        
        if comb not in rotbond_torsions.keys():
          rotbond_torsions[comb] = [tor]
        else:
          rotbond_torsions[comb].append(tor)
    
    # record atomic numbers of each atom
    atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    torsionrestraints_1 = [] 
    for key, tors in rotbond_torsions.items():
      first_atoms = []
      last_atoms = []
      for tor in tors:
        a, x, y, b = tor
        if a not in first_atoms:
          first_atoms.append(a)
        if b not in last_atoms:
          last_atoms.append(b)
      
      first_atoms_atomic = [atomic_numbers[a] for a in first_atoms]
      last_atoms_atomic = [atomic_numbers[a] for a in last_atoms]
      
      # sort atoms based on their atomic numbers
      first = [val for _, val in sorted(zip(first_atoms_atomic, first_atoms), reverse=True)]
      last = [val for _, val in sorted(zip(last_atoms_atomic, last_atoms), reverse=True)]
     
      # select top2 atoms
      tmp = []
      for f in first:
        if len(tmp) < 2 and (f not in tmp):
          tmp.append(f)
      first = tmp
      
      # select top2 atoms
      tmp = []
      for f in last:
        if len(tmp) < 2 and (f not in tmp):
          tmp.append(f)
      last = tmp
    
      # add torsion to torsionrestraints_1
      if (len(first) == 2) and (len(last) == 2):
        torsionrestraints_1.append([first[0], x, y, last[1]])
        torsionrestraints_1.append([first[1], x, y, last[0]])
      elif (len(first) == 2) and (len(last) == 1):
        torsionrestraints_1.append([first[0], x, y, last[0]])
        torsionrestraints_1.append([first[1], x, y, last[0]])
      elif (len(first) == 1) and (len(last) == 2):
        torsionrestraints_1.append([first[0], x, y, last[0]])
        torsionrestraints_1.append([first[0], x, y, last[1]])
      elif (len(first) == 1) and (len(last) == 1):
        torsionrestraints_1.append([first[0], x, y, last[0]])
      else:
        continue

    # change to 1-Indexed, required by Poltype
    for tor in torsionrestraints_1:
      tor_plus_1 = [i+1 for i in tor]
      torsionrestraints.append(tor_plus_1)
    return torsionrestraints

def GeometryOptimization(poltype,mol,totcharge,suffix='1',loose=False,checkbonds=True,modred=True,bondanglerestraints=None,skipscferror=False,charge=None,skiperrors=False,overridecheckterm=False): # specify charge instead of reading from mol if charge!=None
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    NATOM_SMALL = 25
    if bondanglerestraints!=None or \
        (poltype.generateextendedconf==False and poltype.userconformation==False) or \
        (poltype.isfragjob and mol.NumAtoms() < NATOM_SMALL and not poltype.generate_symm_frag_conf): # then vdw opt
        torsionrestraints=[]
    else: # see if need to restrain torsion in extended conformation
        torsionrestraints=FindTorsionRestraints(poltype)

    logoptfname=poltype.logoptfname.replace('_1','_'+suffix)
    comoptfname=poltype.comoptfname.replace('_1','_'+suffix)
    chkoptfname=poltype.chkoptfname.replace('_1','_'+suffix)

    # Skip optimization if userxyzgeometry is provided
    if poltype.userxyzgeometry != '':
        error = False
        optmol =  load_structfile(poltype, f'../{poltype.userxyzgeometry}')
        optmol=rebuild_bonds(poltype,optmol,mol)
        optxyz = poltype.logoptfname.replace('.log','.xyz')
        os.system(f'cp ../{poltype.userxyzgeometry} {optxyz}')
    else:
        # Define which software will be used for the QM optimization
        
        # Only use Gaussian when users require
        if poltype.use_gaus or poltype.use_gausoptonly:
          Soft = 'Gaussian'
        # Otherwise use Psi4 
        else:
          Soft = 'Psi4' 
        
        # Replace Psi4 with PySCF when pcm is needed
        if (poltype.optpcm==True or (poltype.optpcm==-1 and poltype.pcm)) and (Soft == 'Psi4'):
          Soft = 'PySCF'
        
        if Soft == 'Gaussian': # try to use gaussian for opt
            term,error=poltype.CheckNormalTermination(logoptfname,errormessages=None,skiperrors=True)
            if not term or overridecheckterm==True:
                
                poltype.WriteToLog("QM Geometry Optimization")
                poltype.WriteToLog(f'The software {Soft} will be used to do the QM optimization')

                gen_optcomfile(poltype,comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,chkoptfname,mol,modred,torsionrestraints)
                cmdstr = 'GAUSS_SCRDIR=' + poltype.scrtmpdirgau + ' ' + poltype.gausexe + " " + comoptfname
                jobtooutputlog={cmdstr:os.getcwd()+r'/'+logoptfname}
                jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
                scratchdir=poltype.scrtmpdirgau
                jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix 
                inputfilepath=os.path.join(os.getcwd(),comoptfname)
                jobtoinputfilepaths={cmdstr:[inputfilepath]}
                jobtooutputfiles={cmdstr:[logoptfname]}
                jobtoabsolutebinpath={cmdstr:poltype.which(poltype.gausexe)}
                if poltype.checkinputonly==True:
                    sys.exit()
                if os.path.isfile(chkoptfname) and os.path.isfile(logoptfname):
                    os.remove(logoptfname) # if chk point exists just remove logfile, there could be error in it and we dont want WaitForTermination to catch error before job is resubmitted by daemon 
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,True) # have to skip errors because setting optmaxcycle to low number in gaussian causes it to crash
                cmdstr = poltype.formchkexe + " " + chkoptfname
                poltype.call_subsystem([cmdstr],True)
            term,error=poltype.CheckNormalTermination(logoptfname,errormessages=None,skiperrors=True)
            if error and term==False and skiperrors==False:
                poltype.RaiseOutputFileError(logoptfname) 
            optmol =  load_structfile(poltype,logoptfname)
            optmol=rebuild_bonds(poltype,optmol,mol)
            
        if Soft == 'PySCF':

            import importlib
                
            try:
                importlib.import_module('pyscf')
                poltype.WriteToLog("PySCF is installed.\n")
            except ImportError:
                poltype.WriteToLog("PySCF is not installed.\n Make sure that pyscf version 2.7.0 is installed along with dft3, dft4 and pyscf-dispersion.\n You can install using: \n        `pip install pyscf==2.7.0 dft3 dft4 pyscf-dispersion`")
                raise ImportError("PySCF is not installed.\n Make sure that pyscf version 2.7.0 is installed along with dft3, dft4 and pyscf-dispersion.\n You can install using: \n        `pip install pyscf==2.7.0 dft3 dft4 pyscf-dispersion`")
        

            gen_optcomfile(poltype,comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,chkoptfname,mol,modred,torsionrestraints)
            term,error=poltype.CheckNormalTermination(logoptfname,errormessages=None,skiperrors=True)
            modred=False

            # Import both PySCF setup class
            from pyscf_setup import PySCF_init_setup
            from pyscf_setup import PySCF_post_run

            # Instantiate the initial pyscf setup object
            Opt_prep = PySCF_init_setup(poltype,os.getcwd(),comoptfname,mol,\
                                   skipscferror,charge)
        
            # Read the gaussian input coordinate
            Opt_prep.read_Gauss_inp_coord()

            # Write the initial xyz file to be used by PySCF
            Opt_prep.write_xyz('init')

            # Write the torsion constraints to be used 
            # during the QM optimization
            Opt_prep.write_geometric_tor_const('init',torsionrestraints)

            # Write the PySCF input file
            Opt_prep.write_PySCF_input(True,False)

            # Define the cmd to run PySCF
            cmd = f'python {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_inp_file} &> {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'

            # Create the dictionnary with the commnd to be used by CallJobsSeriallyLocalHost
            cmd_to_run = {cmd: f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'}

            # If the user only want to create input file
            if poltype.checkinputonly==True:
                sys.exit()

            # Make sure to delete the PySCF output file 
            # before running it
            if os.path.isfile(f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'):
                os.remove(f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}')

            finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(cmd_to_run,skiperrors)

            # Check normal termination of PySCF job
            term,error=poltype.CheckNormalTermination(f'{Opt_prep.PySCF_out_file}',None,skiperrors)

            if error and term==False and skiperrors==False:
                # This poltype method does not exists...
                poltype.RaiseOutputFileError(f'{Opt_prep.PySCF_out_file}')

            # Instantiate the post process PySCF class
            Opt_post = PySCF_post_run(poltype,os.getcwd(),f'{Opt_prep.PySCF_out_file}',mol)

            # Read the output file from PySCF
            Opt_post.read_out()

            # Use OpenBabel to build a molecule out of the optimized structure
            optmol =  load_structfile(poltype,Opt_post.final_opt_xyz)
            optmol=rebuild_bonds(poltype,optmol,mol)
        
            # Here we redefine logoptfname and poltype.logoptfname
            # to make sure that hereafter, the pyscf output is used!
            logoptfname = Opt_prep.PySCF_out_file
            poltype.logoptfname = Opt_prep.PySCF_out_file 

        if Soft == 'Psi4':
            gen_optcomfile(poltype,comoptfname,poltype.numproc,poltype.maxmem,poltype.maxdisk,chkoptfname,mol,modred,torsionrestraints)
            term,error=poltype.CheckNormalTermination(logoptfname,errormessages=None,skiperrors=True)
            modred=False

            inputname,outputname=CreatePsi4OPTInputFile(poltype,comoptfname,comoptfname,mol,modred,bondanglerestraints,skipscferror,charge,loose,torsionrestraints)
            if term==False or overridecheckterm==True:
                poltype.WriteToLog("QM Geometry Optimization")
                poltype.WriteToLog(f'The software {Soft} will be used to do the QM optimization')
                cmdstr=' '.join(['psi4 ',poltype.psi4_args,inputname,logoptfname])
                jobtooutputlog={cmdstr:os.getcwd()+r'/'+logoptfname}
                jobtolog={cmdstr:os.getcwd()+r'/'+poltype.logfname}
                scratchdir=poltype.scrtmpdirpsi4
                jobtologlistfilepathprefix=os.getcwd()+r'/'+'optimization_jobtolog_'+poltype.molecprefix
                inputfilepath=os.path.join(os.getcwd(),inputname)
                jobtoinputfilepaths={cmdstr:[inputfilepath]}
                jobtooutputfiles={cmdstr:[logoptfname]}
                jobtoabsolutebinpath={cmdstr:poltype.which('psi4')}
        
                if poltype.checkinputonly==True:
                    sys.exit()

                if os.path.isfile(logoptfname):
                    os.remove(logoptfname)
                finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(jobtooutputlog,True) # have to skip errors because setting optmaxcycle to low number in gaussian causes it to crash

            term,error=poltype.CheckNormalTermination(logoptfname,None,skiperrors) # now grabs final structure when finished with QM if using Psi4
            if error and term==False and skiperrors==False:
                poltype.RaiseOutputFileError(logoptfname) 
            GrabFinalXYZStructure(poltype,logoptfname,logoptfname.replace('.log','.xyz'),mol)
            optmol =  load_structfile(poltype,logoptfname.replace('.log','.xyz'))
            optmol=rebuild_bonds(poltype,optmol,mol)
        optmol.SetTotalCharge(totcharge)

        if Soft != 'PySCF':
            GrabFinalXYZStructure(poltype,logoptfname,logoptfname.replace('.log','.xyz'),mol)

    return optmol,error,torsionrestraints


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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    for b in openbabel.OBMolBondIter(refmol):
        beg = b.GetBeginAtomIdx()
        end = b.GetEndAtomIdx()
        if not newmol.GetBond(beg,end):
            newmol.AddBond(beg,end, b.GetBondOrder(), b.GetFlags())
        else:
            newb=newmol.GetBond(beg,end)
            bondorder=newb.GetBondOrder()
            newb.SetBondOrder(b.GetBondOrder())

    return newmol


def PruneBonds(poltype,mol,bondtopology):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
