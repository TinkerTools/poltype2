

##################################################################
#
# Title: Poltype
# Description: Atomic typer for the polarizable AMOEBA force field
#
# Copyright:            Copyright (c) Johnny C. Wu,
#                   Gaurav Chattree, & Pengyu Ren 2010-2011
#
# Poltype is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3
# as published by the Free Software Foundation.
#
# Poltype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# if not, write to:
# Free Software Foundation, Inc.
# 59 Temple Place, Suite 330
# Boston, MA 02111-1307  USA
#
##################################################################

import os
import sys
#sys.settrace
import heapq
import time
from collections import deque
from types import *
import re
import string
import getopt
import tempfile
import shutil
import pprint
import subprocess
import inspect
from socket import gethostname
import traceback
from math import *
import itertools

import numpy
from scipy import optimize
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import openbabel
import string as st
import valence
import pybel as pb
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)
from scipy.optimize import fmin
import collections



# Implementation Notes
# 1) Minimize Structure
# 2) Run SP Calculations
# 3) Run GDMA                             DONE (electrostatic-param.pl)
# 4) Run poledit (to extract atom types)  DONE (electrostatic-param.pl)
# 4) Run avgmpoles.pl

# Default starting index of atom type
prmstartidx = 401

# These values can be edited. Affect speed of QM calculations.
numproc = 1
maxmem = "700MB"
maxdisk = "100GB"

# Initilize directories
gausdir = None
gdmadir = None
tinkerdir = None
scratchdir = "/scratch/bdw2292"
paramfname = sys.path[0] + "/amoeba_v2_new.prm"
paramhead = sys.path[0] + "/amoeba_v2_new_head.prm"
obdatadir = sys.path[0] + "/datadir"
amoeba_conv_spec_fname = "amoeba_t5_t4.txt"

# Initilize executables
babelexe = "babel"
gausexe =  "g09"
formchkexe =  "formchk"
cubegenexe =  "cubegen"
gdmaexe = "gdma"
eleparmexe = sys.path[0] + "/electrostatic-param.pl"
groupsymexe = sys.path[0] + "/groupsym.exe"
avgmpolesexe = sys.path[0] + "/avgmpoles.pl"
peditexe = "poledit.x"
potentialfit = sys.path[0] + "/potential-fit.py"
potentialexe = "potential.x"
valenceexe = "valence.x"
minimizeexe = "minimize.x"
analyzeexe = "analyze.x"
superposeexe = "superpose.x"
#peditexe = "poledit"
#potentialexe = "potential"
#valenceexe = "valence"
#minimizeexe = "minimize"
#analyzeexe = "analyze"
#superposeexe = "superpose"

# Initialize constants, basis sets
defopbendval = 0.20016677990819662
Hartree2kcal_mol = 627.5095
optbasisset = "6-31G*"               # This is good; bigger is better but slower
toroptbasisset = "6-31G*"            # Same as above
dmabasisset = "cc-pVTZ"              # This can not be changed
popbasisset = "6-31G*"               # We typically dont do pop analysis unless we need bond orders from Gaussain (not sure PSI4 has this).
espbasisset = "aug-cc-pVTZ"          # This can be from 6-311++G**, 6-311++G(2d,2p), or aug-cc-pvtz. For small molecules like yours, use the last two.
torspbasisset = "6-311++G**"         # This can be 6-31+G* or bigger, depending on the cost (we have to do this for many conformations)


#Initialize Methods
optmethod='MP2'                      # MP2 for small molecules, DFT for large (wB97XD)
toroptmethod='HF'                    # HF or DFT in gas. HF may be needed for large molecules unless we break down to fragments for torsion fitting.
torspmethod='MP2'                    # DFT or MP2 in PCM. DFT may be good enough (comparison would be nice). 
dmamethod='MP2'                      # MP2
espmethod='MP2'                      # MP2


# Initialize some global variables such as arrays and booleans
qmonly = False
espfit = True
parmtors = True
uniqidx = False
torkeyfname = None
foldnum=3
nfoldlist = range(1,foldnum+1)
#nfoldlist = range(1,7)
foldoffsetlist = [ 0.0, 180.0, 0.0, 0.0, 0.0, 0.0 ]
torlist = []
rotbndlist = []
fitrotbndslist=[]
mm_tor_count = 1
output_format = 5
toromit_list = []
omittorsion2 = False
do_tor_qm_opt = True
gencomonly=False
fitrotbnds=False # list of rotatable bonds to be fitted with multiple torsions around bond (without hydrogen), boolean that sets option to read from this list as False first
solventname="water"
maxRMSD=.1
maxRMSPD=1
tordatapointsnum=None
gentorsion=False # boolean to specify when program is running gaussian for torsion
gaustorerror=False # boolean to specify when gaussian crashes during gaussian opt or SP for torsion, this way can exclude some points for fitting 
# Poltype begins with the 'main' method which is found towards the bottom of the program
torsionrestraint=.1
onlyrotbndlist=[] # initial empty list default, do not spin only one bond

rotalltors=False
dontdotor=False
dontdotorfit=False
onlyrotbnd=False
toroptpcm=False
optmaxcycle=400
optpcm=False
sppcm=False
torsppcm=False
use_psi4=False
use_psi4SPonly=False
freq=False
postfit=False # just a boolean specifying where in torsion fitting process is


def call_subsystem(cmdstr, iscritical=False):
    """
    Intent: Run 'cmdstr' on the command line
    Input:
        cmdstr: command string to be run on command line
    Output:
    Description: -
    """
    now = time.strftime("%c",time.localtime())
    logfh.write(now + " Calling: " + cmdstr + "\n")
    logfh.flush()
    os.fsync(logfh.fileno())
    p = subprocess.Popen(cmdstr, shell=True,
       stdout=logfh, stderr=logfh)
    if p.wait() != 0:
        now = time.strftime("%c",time.localtime())
        logfh.write(now + " ERROR: " + cmdstr + "\n")
        logfh.flush()
        os.fsync(logfh.fileno())
        if iscritical and gentorsion==True:
            gaustorerror=True
            sys.exit(1) # just exit until torsion spinning is parrallel not sequential
        elif iscritical and gentorsion==False:
            sys.exit(1)
    return p.wait()

def which(program,pathlist=os.environ["PATH"]):
    """
    Intent: Check if the 'program' is in the user's path
    Input:
        program: Program name
        pathlist: members of the user's path
    Output:
       program: path to program
       exe_file: path to program
       None: program was not found
    Referenced By: initialize 
    Description: -
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in pathlist.split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def rotate_list(l1):
    deq = deque(l1)
    deq.rotate(-1)
    return list(deq)

def parse_options(argv):
    """
    Intent: Set up variables based on arguments supplied by user
    Input:
        argv: command line arguments supplied by user
    Output:
    Referenced By: main
    Description:
    """
    global molecprefix
    global qmonly
    global molstructfname
    global chkname
    global fname
    global gausfname
    global gausoptfname
    global gdmafname
    global comfname
    global prmstartidx
    global numproc
    global maxmem
    global maxdisk
    global espfit
    global parmtors
    global uniqidx
    global torkeyfname
    global optbasisset
    global dmabasisset
    global popbasisset
    global espbasisset
    global comoptfname
    global chkoptfname
    global fckoptfname
    global logoptfname
    global compopfname
    global chkpopfname
    global fckpopfname
    global logpopfname
    global comdmafname
    global chkdmafname
    global fckdmafname
    global logdmafname
    global comespfname
    global chkespfname
    global fckespfname
    global logespfname
    global output_format
    global paramhead
    global omittorsion2
    global do_tor_qm_opt
    global optmethod
    global gencomonly
    global solventname
    global toroptmethod
    global torspmethod
    global dmamethod
    global espmethod
    global optbasisset
    global toroptbasisset
    global dmabasisset
    global popbasisset
    global espbasisset
    global torspbasisset
    global fitrotbnds
    global maxRMSD
    global maxRMSPD
    global fitrotbnds
    global tordatapointsnum
    global gentorsion
    global gaustorerror
    global foldnum
    global postfit
    global torsionrestraint
    global onlyrotbnd
    global rotalltors
    global optmaxcycle
    global optpcm
    global use_psi4
    global use_psi4SPonly
    global sppcm
    global torsppcm
    global freq
    global scratchdir


    try:
        opts, xargs = getopt.getopt(argv[1:],'hqn:m:M:a:s:p:d:u:',['torspmethod=','espmethod=','toroptmethod=','freq','torsppcm','sppcm','toroptpcm','use_psi4','use_psi4SPonly','dontdotor','dontdotorfit','optpcm','optmaxcycle=','rotalltors',"onlyrotbnd=","torsionrestraint=","foldnum=","tordatapointsnum=","help","qmonly","optbasisset=","dmabasisset=","popbasisset=","espbasisset=","torspbasisset=","optlog=","dmalog=","esplog=","dmafck=","espfck=","numproc=","maxmem=","maxdisk=","atmidx=","structure=","prefix=","gdmaout=","gbindir=","qm-scratch-dir=","omit-espfit","omit-torsion","test-tor-key=","uniqidx","tinker4format","omit-torsion2","dont-tor-qm-opt","optmethod=","gencomonly","solventname=","fitrotbnds"])
    except Exception as e: 
        print(e)
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("--rotalltors"):
            rotalltors = True
            continue
        elif o in ("-s", "--structure"):
            molstructfname = a
            continue
        elif o in ("--torsppcm"):
            torsppcm = True
            continue
        elif o in ("--freq"):
            freq = True
            continue
        elif o in ("--sppcm"):
            sppcm = True
            continue
        elif o in ("--optpcm"):
            optpcm = True
            continue
        elif o in ("--toroptpcm"):
            optpcm = True
            continue
        elif o in ("--use_psi4"):
            use_psi4 = True
            continue
        elif o in ("--use_psi4SPonly"):
            use_psi4SPonly = True
            continue
        elif o in ("--dontdotor"):
            dontdotor = True
            continue
        elif o in ("--dontdotorfit"):
            dontdotorfit = True
            continue
        elif o in ("--optmaxcycle"):
            optmaxcycle = a
            continue

        elif o in ("--torsionrestraint"):
            torsionrestraint=float(a)
            continue
        elif o in ("--foldnum"):
            foldnum=int(a)
            continue
        elif o in ("--tordatapointsnum"):
            tordatapointsnum=int(a)
            continue
        elif o in ("--fitrotbnds"):
            fitrotbnds=True
            continue
        elif o in ("--solventname"):
            solventname = a
            continue
        elif o in ("--optmethod"):
            optmethod = a
            continue
        elif o in ("--espmethod"):
            espmethod = a
            continue
        elif o in ("--torspmethod"):
            torspmethod = a
            continue
        elif o in ("--toroptmethod"):
            toroptmethod = a
            continue
        elif o in ("--gencomonly"):
            gencomonly = True
            continue
        elif o in ("-n", "--numproc"):
            numproc = a
            continue
        elif o in ("-m", "--maxmem"):
            maxmem = a
            continue
        elif o in ("-M", "--maxdisk"):
            maxdisk = a
            continue
        elif o in ("-a", "--atmidx"):
            prmstartidx = int(a)
            continue
        elif o in ("--optbasisset"):
            optbasisset = a
            continue
        elif o in ("--dmabasisset"):
            dmabasisset = a
            continue
        elif o in ("--popbasisset"):
            popbasisset = a
            continue
        elif o in ("--espbasisset"):
            espbasisset = a
            continue
        elif o in ("--torspbasisset"):
            torspbasisset = a
            continue
        elif o in ("--optlog"):
            logoptfname = a
            continue
        elif o in ("--dmalog"):
            logdmafname = a
            continue
        elif o in ("--esplog"):
            logespfname = a
            continue
        elif o in ("--dmafck"):
            fckdmafname = a
            continue
        elif o in ("--espfck"):
            fckespfname = a
            continue
        elif o in ("--dmachk"):
            chkdmafname = a
            continue
        elif o in ("--espchk"):
            chkespfname = a
            continue
        elif o in ("-f", "--formchk"):
            fname = a
            continue
        elif o in ("-d", "--gdmaout"):
            gdmafname = a
            continue
        elif o in ("-u", "--gbindir"):
            gausdir = a
            continue
        elif o in ("-q", "--qmonly"):
            qmonly = True
            continue
        elif o in ("--omit-espfit"):
            espfit = False
            continue
        elif o in ("--omit-torsion"):
            parmtors = False
            continue
        elif o in ("--omit-torsion2"):
            omittorsion2 = True
            continue
        elif o in ("--dont-tor-qm-opt"):
            do_tor_qm_opt = False
            continue
        elif o in ("--test-tor-key"):
            torkeyfname = a
            continue
        elif o in ("--uniqidx"):
            uniqidx = True
            continue
        elif o in ("-h", "--help"):
            usage()
            sys.exit(2)
            continue
        elif o in ("--version"):
            printversion()
            sys.exit(2)
            continue
        elif o in ("--tinker4format"):
            output_format = 4 
            paramhead = sys.path[0] + "/amoeba_v2_new_head_tinker_4.prm"
            continue
        elif o in ("--onlyrotbnd"): # comma seperated list 
            onlyrotbnd=[i.lstrip().rstrip() for i in a.split(',')]
            continue
        else:
            assert False, "unhandled option"
            

    if "molstructfname" not in globals():
        print("Molecule structure file (-s) needs to be defined.")
        sys.exit()  


class PrettyFloat(float):
    def __repr__(self):
        return '%.3f' % self

def pretty_floats(obj):
    if isinstance(obj, float):
        return PrettyFloat(obj)
    if isinstance(obj, openbabel.OBAtom):
        return obj.GetIdx()
    elif isinstance(obj, dict):
        return dict((k, pretty_floats(v)) for k, v in obj.items())
    elif isinstance(obj, (list, tuple, numpy.ndarray)):
        return map(pretty_floats, obj)
    return obj

def dispvar(label, *vars):
    """
    Intent: Output variable to stdout
    """
    print(label + ': ')
    for var in vars:
        print(pretty_floats(var))
    print("")

def print_error(errstrarr,kill=None):
    now = time.strftime("%c",time.localtime())
    sys.stderr.write('ERROR (%s): ' % now)
    logfh.write('ERROR (%s): ' % now)
    if isinstance(errstrarr, (list,tuple)):
        for errstr in errstrarr:
            sys.stderr.write(str(pretty_floats(errstr)) + '\n')
            logfh.write(str(pretty_floats(errstr)) + '\n')
    else:
        sys.stderr.write(str(pretty_floats(errstrarr)) + '\n')
        logfh.write(str(pretty_floats(errstrarr)) + '\n')


    if kill is not None:
        if isinstance(kill,int):
            sys.exit(kill)
        else:
            sys.exit(255)

def write_arr_to_file(fname, array_list):
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
            outfh.write("%10.4f" % ele)
        outfh.write("\n")

def initialize ():
    """
    Intent: Initialize all paths to needed executables
    Input:
    Output:
    Referenced By: main
    Description: -
    """
    global gausdir
    global gausexe
    global formchkexe
    global cubegenexe
    global tinkerdir
    global peditexe
    global potentialexe
    global valenceexe
    global minimizeexe
    global analyzeexe
    global superposeexe
    global gdmaexe

    if not use_psi4:
        if (gausdir is not None):
            if which(os.path.join(gausdir,"g09")) is not None:
                gausexe    = os.path.join(gausdir,"g09")
                formchkexe = os.path.join(gausdir,formchkexe)
                cubegenexe = os.path.join(gausdir,cubegenexe)
            elif which(os.path.join(gausdir,"g03")) is not None:
                gausexe    = os.path.join(gausdir,"g03")
                formchkexe = os.path.join(gausdir,formchkexe)
                cubegenexe = os.path.join(gausdir,cubegenexe)
            else:
                print("ERROR: Invalid Gaussian directory: ", gausdir)
                sys.exit(1)
        else:
            if which("g09") is not None:
                gausexe    = "g09"
            elif which("g03") is not None:
                gausexe    = "g03"
            else:
                print("ERROR: Cannot find Gaussian executable in $PATH. Please install Gaussian or specify Gaussian directory with --gbindir flag.")
                sys.exit(1)


            


    cmdstr=analyzeexe+' '+sys.path[0]+r'/'+' '+'water.xyz'+' '+'-k'+' '+sys.path[0]+r'/'+'water.key'+' '+'e'+'>'+' '+'version.out'
    try:
        print('Calling: '+cmdstr) 
        returned_value = subprocess.call(cmdstr, shell=True)
    except:
        pass
    temp=open('version.out','r')
    results=temp.readlines()
    temp.close()
    latestversion = False
    for line in results:
        if "Version" in line:
            linesplit=line.split()
            versionnum=float(linesplit[2])
            if versionnum>=8.7:
                latestversion = True
                break
    print('versionnum ',versionnum)   
    if(not latestversion):
        print("Notice: Not latest version of tinker (>=8.7)")
  
    if ("TINKERDIR" in os.environ):
        tinkerdir = os.environ["TINKERDIR"]
        peditexe = os.path.join(tinkerdir,peditexe)
        potentialexe = os.path.join(tinkerdir,potentialexe)
        valenceexe = os.path.join(tinkerdir,valenceexe)
        minimizeexe = os.path.join(tinkerdir,minimizeexe)
        analyzeexe = os.path.join(tinkerdir,analyzeexe)
        superposeexe = os.path.join(tinkerdir,superposeexe)

    if (not which(analyzeexe)):
        print("ERROR: Cannot find TINKER analyze executable")
        sys.exit(2)

    if not use_psi4 and not use_psi4SPonly:
        if ("GDMADIR" in os.environ):
            gdmadir = os.environ["GDMADIR"]
            gdmaexe = os.path.join(gdmadir,gdmaexe)

        if (not which(gdmaexe)):
            print("ERROR: Cannot find GDMA executable")
            sys.exit(2)

        if ("GAUSS_SCRDIR" in os.environ):
            scratchdir = os.environ["GAUSS_SCRDIR"]
            if not os.path.isdir(scratchdir):
                os.mkdir(scratchdir)
            vfs = os.statvfs(scratchdir)
            gbfree = (vfs.f_bavail * vfs.f_frsize) / (1024*1024*1024)
            if(float(maxdisk[:-2]) > gbfree):
                print("ERROR: maxdisk greater than free space in scratch directory")
                sys.exit(2)

        if (not which(scratchdir)):
            print("ERROR: Cannot find Gaussian scratch directory")
            sys.exit(2)


    #os.putenv('BABEL_DATADIR',obdatadir)

def init_filenames ():
    """
    Intent: Initialize file names
    Input:
    Output:
    Referenced By: main
    Description: -
    """
    global logfname
    global molecprefix
    global molstructfname
    global chkname
    global fname
    global gausfname
    global gausoptfname
    global gdmafname
    global keyfname
    global xyzfname
    global peditinfile
    global valinfile
    global superposeinfile
    global espgrdfname
    global qmespfname
    global qmesp2fname
    global grpfname
    global key2fname
    global key3fname
    global key4fname
    global key5fname
    global xyzoutfile
    global valoutfname
    global scrtmpdir
    global tmpxyzfile
    global tmpkeyfile

    head, molstructfname = os.path.split(molstructfname)
    molecprefix =  os.path.splitext(molstructfname)[0]
    logfname = assign_filenames ( "logfname" , "-poltype.log")
    chkname = assign_filenames ( "chkname" , ".chk")
    fname = assign_filenames ( "fname" , ".fchk")
    gausfname = assign_filenames ( "gausfname" , ".log")
    gausoptfname = assign_filenames ( "gausoptfname" , "-opt.log")
    gdmafname = assign_filenames ( "gdmafname" , ".gdmaout")
    keyfname = assign_filenames ( "keyfname" , ".key")
    xyzfname = assign_filenames ( "xyzfname" , ".xyz")
    peditinfile = assign_filenames ( "peditinfile" , "-peditin.txt")
    valinfile = assign_filenames ( "valinfile" , "-valin.txt")
    superposeinfile = assign_filenames ( "superposeinfile" , "-superin.txt")
    espgrdfname = assign_filenames ( "espgrdfname" , ".grid")
    qmespfname = assign_filenames ( "qmespfname" , ".cube")
    #qmesp2fname = assign_filenames ( "qmesp2fname" , ".cube_2")
    qmesp2fname = assign_filenames ( "qmesp2fname" , ".pot")
    grpfname = assign_filenames ( "grpfname" , "-groups.txt")
    key2fname = assign_filenames ( "key2fname" , ".key_2")
    key3fname = assign_filenames ( "key3fname" , ".key_3")
    key4fname = assign_filenames ( "key4fname" , ".key_4")
    key5fname = assign_filenames ( "key5fname" , ".key_5")
    xyzoutfile = assign_filenames ( "xyzoutfile" , ".xyz_2")
    valoutfname = assign_filenames ( "valoutfname" , "sp.valout")
    scrtmpdir = scratchdir + '/Gau-' + molecprefix
    tmpxyzfile = 'ttt.xyz'
    tmpkeyfile = 'ttt.key'

def assign_filenames (filename,suffix):
    global molecprefix
    if filename in globals():
        return eval(filename)
    else:
        return molecprefix + suffix

def printversion ():
    print(os.path.basename(sys.argv[0]) + ''': Version 2.0.0''')

def copyright ():
    print('''
              Poltype -- Polarizable atom typer of small molecules for the
                                 AMOEBA polarizable force field
   Please cite:
   Wu, J.C.; Chattree, G.; Ren, P.Y.; Automation of AMOEBA polarizable force field 
   parameterization for small molecules. Theor Chem Acc. (Accepted).
                                  Version 2.0.0 June 2015
 Copyright (c)  Johnny Wu, Gaurav Chattree, Brandon Walker, Matthew Harger and Pengyu Ren 2019
                                      All Rights Reserved
 ######################################################################################################
 '''
  )

def usage ():
    print('''poltype:
    -h, --help      -- displays this help message
    -s, --structure -- Gaussian checkpoint file
    -a, --atmidx     -- Starting index of atom type
    -m MEMORY , --maxmem    -- 
    -M, --maxdisk
    -q, --qmonly
    --optbasisset
    --dmabasisset
    --espbasisset
    --torspbasisset
    --omit-espfit
    --omit-torsion
    --version       -- displays version of script''')

def load_structfile(structfname):
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

#def verify_topology (newmol, refmol):
    #for atm in openbabel.OBMolAtomIter(refmol):
    #    newatm = newmol.

    #return True

def rebuild_bonds(newmol, refmol):
    for b in openbabel.OBMolBondIter(refmol):
        beg = b.GetBeginAtomIdx()
        end = b.GetEndAtomIdx()
        if not newmol.GetBond(beg,end):
            newmol.AddBond(beg,end, b.GetBO(), b.GetFlags())

    return True 

def get_class_number(idx):
    """
    Intent: Given an atom idx, return the atom's class number
    """
    maxidx =  max(symmetryclass)
    return prmstartidx + (maxidx - symmetryclass[idx - 1])

def get_class_key(a, b, c, d):
    """
    Intent: Given a set of atom idx's, return the class key for the set (the class numbers of the atoms appended together)
    """
    cla = get_class_number(a)
    clb = get_class_number(b)
    clc = get_class_number(c)
    cld = get_class_number(d)

    if ((clb > clc) or (clb == clc and cla > cld)):
        return '%d %d %d %d' % (cld, clc, clb, cla)
    return '%d %d %d %d' % (cla, clb, clc, cld)

def get_uniq_rotbnd(a, b, c, d):
    """
    Intent: Return the atom idx's defining a rotatable bond in the order of the class key
    found by 'get_class_key'
    """
    cla = get_class_number(a)
    clb = get_class_number(b)
    clc = get_class_number(c)
    cld = get_class_number(d)

    tmpkey = '%d %d %d %d' % (cla,clb,clc,cld)
    if (get_class_key(a,b,c,d) == tmpkey):
        return (a, b, c, d)
    return (d, c, b, a)

def obatom2idx(obatoms):
    """
    Intent: Given a list of OBAtom objects, return a list of their corresponding idx's
    Referenced By: get_torlist_opt_angle
    """
    atmidxlist = []
    for obatm in obatoms:
        atmidxlist.append(obatm.GetIdx())
    return atmidxlist

def save_structfile(molstruct, structfname):
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
        maxidx =  max(symmetryclass)
        iteratom = openbabel.OBMolAtomIter(molstruct)
        etab = openbabel.OBElementTable()
        tmpfh.write('%6d   %s\n' % (molstruct.NumAtoms(), molstruct.GetTitle()))
        for ia in iteratom:
            tmpfh.write( '%6d %2s %13.6f %11.6f %11.6f %5d' % (ia.GetIdx(), etab.GetSymbol(ia.GetAtomicNum()), ia.x(), ia.y(), ia.z(), prmstartidx + (maxidx - symmetryclass[ia.GetIdx() - 1])))
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

def save_structfileXYZ(molstruct, structfname):
        
    tmpconv = openbabel.OBConversion()
    tmpconv.SetOutFormat('xyz')
    return tmpconv.WriteFile(molstruct, structfname)

def rads(degrees):
    """
    Intent: Convert degrees to radians
    """
    if type(degrees) == type(list):
        return [ deg * pi / 180 for deg in degrees ]
    else:
        return degrees * pi / 180

def getCan(x):
    """
    Intent: For a given atom, output it's canonical label if it exists
    Input: 
        x: atom in question
    Output:
        canonical label if it exists, -1 if not
    """
    if (str(type(mol.GetAtom(x))).find("OBAtom") >= 0):
        return canonicallabel[x-1]
    else:
        return -1

def get_symm_class(x):
    """
    Intent: For a given atom, output it's symmetry class if it exists
    Input: 
        x: atom in question
    Output:
        symmetry class if it exists, -1 if not
    """
    if (str(type(mol.GetAtom(x))).find("OBAtom") >= 0):
        return symmetryclass[x-1]
    else:
        return -1

def prepend_keyfile(keyfilename):
    """
    Intent: Adds a header to the key file given by 'keyfilename'
    """
    tmpfname = keyfilename + "_tmp"
    tmpfh = open(tmpfname, "w")
    keyfh = open(keyfilename, "r")

    tmpfh.write("parameters " + paramhead + "\n")
    tmpfh.write("bondterm none\n")
    tmpfh.write("angleterm none\n")
    tmpfh.write("torsionterm none\n")
    tmpfh.write("vdwterm none\n")
    tmpfh.write("fix-monopole\n")
    tmpfh.write("potential-offset 1.0\n\n")
    #tmpfh.write('openmp-threads 4'+'\n')
    for line in keyfh:
        tmpfh.write(line)
    shutil.move(tmpfname, keyfilename)

def scale_multipoles (symmclass, mpolelines,scalelist):
    """
    Intent: Scale multipoles based on value in scalelist
    """
    symmclass = int(symmclass)
    if scalelist[symmclass][2]:
        qp1 = float(mpolelines[2])
        qp1 = map(float, mpolelines[2].split())
        qp1 = [ scalelist[symmclass][2] * x for x in qp1 ]
        mpolelines[2] = '%46.5f\n' % qp1[0]
        qp2 = map(float, mpolelines[3].split())
        qp2 = [ scalelist[symmclass][2] * x for x in qp2 ]
        mpolelines[3] = '%46.5f %10.5f\n' % tuple(qp2)
        qp3 = map(float, mpolelines[4].split())
        qp3 = [ scalelist[symmclass][2] * x for x in qp3 ]
        mpolelines[4] = '%46.5f %10.5f %10.5f\n' % tuple(qp3)
    return mpolelines

def rm_esp_terms_keyfile(keyfilename):
    """
    Intent: Remove unnecessary terms from the key file
    """
    tmpfname = keyfilename + "_tmp"
    cmdstr = "sed -e \'/\(.* none$\|^#\|^fix\|^potential-offset\)/d\' " +  keyfilename + " > " + tmpfname
    call_subsystem(cmdstr)
    shutil.move(tmpfname, keyfilename)
    # Removes redundant multipole definitions created after potential fitting
    keyfh = open(keyfilename)
    lines = keyfh.readlines()
    newlines = []
    tmpfh = open(tmpfname, "w")
    mpolelines = 0
    for ln1 in range(len(lines)):
        if mpolelines > 0:
            mpolelines -= 1
            continue
        if 'multipole' in lines[ln1]:
            nllen = len(newlines)
            for ln2 in range(len(newlines)):
                if len(lines[ln1].split()) <= 1 or len(lines[ln2].split()) <= 1:
                    continue
                if lines[ln1].split()[0] in newlines[ln2].split()[0] and \
                lines[ln1].split()[0] in 'multipole' and \
                lines[ln1].split()[1] in newlines[ln2].split()[1]:
                    newlines[ln2] = lines[ln1]
                    newlines[ln2+1] = lines[ln1+1]
                    newlines[ln2+2] = lines[ln1+2]
                    newlines[ln2+3] = lines[ln1+3]
                    newlines[ln2+4] = lines[ln1+4]
                    mpolelines = 4
                    break
            if mpolelines == 0:
                newlines.append(lines[ln1])
        else:
            newlines.append(lines[ln1])

    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()
    shutil.move(tmpfname, keyfilename)

def post_proc_localframes(keyfilename, lfzerox,atomindextoremovedipquad,atomindextoremovedipquadcross):
    """
    Intent: This method runs after the tinker tool Poledit has run and created an
    initial *.key file. The local frames for each multipole are "post processed". 
    Zeroed out x-components of the local frame are set back to their original values
    If certain multipole values were not zeroed out by poltype this method zeroes them out
    Input:
       keyfilename: string containing file name of *.key file
       lfzerox: array containing a boolean about whether the x-component of the local frame
                for a given atom should be zeroed or not. Filled in method 'gen_peditin'.
                'lfzerox' is true for atoms that are only bound to one other atom (valence = 1)
                that have more than one possible choice for the x-component of their local frame
    Output: *.key file is edited
    Referenced By: main
    Description:
    1. Move the original *.key file to *.keyb
    2. Create new empty file *.key
    3. Iterate through the lines of *.keyb. 
        a. If poledit wrote out the local frame with an x-component missing or as 0
        Then rewrite it with the original x-component (lf2) found in gen_peditin
        b. If poledit did not zero out the local frame x-component for an atom but 
        lfzerox is true, zero out the necessary multipole components manually
    4.  Write the edited multipole information to *.key 
    """
    # mv *.key to *.keyb
    tmpfname = keyfilename + "b"
    shutil.move(keyfilename, tmpfname)
    keyfh = open(tmpfname)
    lines = keyfh.readlines()
    newlines = []
    # open *.key which is currently empty
    tmpfh = open(keyfilename, "w")
    mpolelines = 0
    # iterate over lines in the *.key file

    for ln1 in range(len(lines)):
        # skip over the lines containing multipole values
        if mpolelines > 0:
            mpolelines -= 1
            continue
        elif 'multipole' in lines[ln1]:
            # Check what poledit wrote as localframe2
            tmplst = lines[ln1].split()
            if len(tmplst) == 5:
                (keywd,atmidx,lf1,lf2,chg) = tmplst
            elif len(tmplst) == 4:
                (keywd,atmidx,lf1,chg) = tmplst
                lf2 = '0'
            elif len(tmplst) == 6:
                (keywd,atmidx,lf1,lf2,lf3,chg) = tmplst
            """
            # If poledit set lf2 to 0, then replace it with the lf2 found in gen_peditin
            if int(lf2) == 0:
                #lf2 = localframe2[int(atmidx) - 1]
                lines[ln1] = '%s %5s %4s %4d %21s\n' % (keywd,
                    atmidx, lf1, lf2, chg)
            """
            # manually zero out components of the multipole if they were not done by poledit
            if lfzerox[int(atmidx) - 1]:
                tmpmp = list(map(float, lines[ln1+1].split()))
                tmpmp[0] = 0
                lines[ln1+1] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+3].split()))
                tmpmp[0] = 0
                lines[ln1+3] = '%46.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)

            if int(atmidx) in atomindextoremovedipquad.keys():
                tmpmp = list(map(float, lines[ln1+1].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                tmpmp[2] = 0
                lines[ln1+1] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+2].split()))
                tmpmp[0] = 0
                lines[ln1+2] = '%46.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+3].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                lines[ln1+3] = '%46.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                tmpmp[2] = 0
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)

            if int(atmidx) in atomindextoremovedipquadcross.keys():

                tmpmp = list(map(float, lines[ln1+3].split()))
                tmpmp[0] = 0
                lines[ln1+3] = '%46.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)

            if len(tmplst) == 4:
                linesplit=re.split(r'(\s+)', lines[ln1])
                newtmplist=linesplit[:len(linesplit)-4]
                newtmplist.append('    0')
                newtmplist.append(linesplit[-4])
                newtmplist.append(linesplit[-3])
                newtmplist.append(linesplit[-2])
                templine=''.join(newtmplist)
                newlines.extend(templine)
                newlines.extend(lines[ln1+1:ln1+5])
            else:   
                newlines.extend(lines[ln1:ln1+5])
            mpolelines = 4
        # append any lines unrelated to multipoles as is
        else:
            newlines.append(lines[ln1])
    # write out the new lines to *.key
    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()


def post_process_mpoles(keyfilename, scalelist):
    """
    Intent: Iterate through multipoles in 'keyfilename' and scale them if necessary
    Calls 'scale_multipoles'
    Input:
        keyfilename: file name for key file
        scalelist: structure containing scaling information. Found in process_types
    Output: new *.key file
    Referenced By: main
    Description: 
    1. Move old keyfile to tmpfname
    2. open an empty key file
    3. iterate through tmpfname scaling the multipoles if necessary
    4. output new multipoles to the new key file
    """
    tmpfname = keyfilename + "b"
    shutil.move(keyfilename, tmpfname)
    keyfh = open(tmpfname)
    lines = keyfh.readlines()
    newlines = []
    tmpfh = open(keyfilename, "w")
    mpolelines = 0
    for ln1 in range(len(lines)):
        if mpolelines > 0:
            mpolelines -= 1
            continue
        elif 'multipole' in lines[ln1] and len(lines[ln1])==5:
            (keywd,symcls,lf1,lf2,chg) = lines[ln1].split()
            newlines.extend(scale_multipoles(symcls,lines[ln1:ln1+5],scalelist))
            mpolelines = 4
        elif 'multipole' in lines[ln1] and len(lines[ln1])==6:
            (keywd,symcls,lf1,lf2,lf3,chg) = lines[ln1].split()
            newlines.extend(scale_multipoles(symcls,lines[ln1:ln1+5],scalelist))
            mpolelines = 4
        else:
            newlines.append(lines[ln1])

    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()

def append_basisset (comfname, spacedformulastr,basissetstr):
    """
    Intent: Append terms to *.com file if necessary 
    Input:
        comfname: com file name
        spacedformulastr: stoichoimetric molecular formula in spaced form (e.g. C 4 H 6 O 1)
        basissetstr: basis set string
    Output:
        comfname is possibly appended to
    Referenced By: gen_optcomfile, gen_comfile, tor_opt_sp 
    Description:
    """
    if ('I ' in spacedformulastr):
        fh = open(comfname, "a")
        fh.write(' ' + re.sub(r'(\d+|I)\s+',r'', spacedformulastr) + '0\n')
        fh.write(' %s\n' % basissetstr)
        fh.write(' ****\n')

#      if re.search('6-31',basissetstr,re.I):
        fh.write(' I     0 \n')
        fh.write(' S   5   1.00\n')
        fh.write('       444750.000              0.0008900        \n')
        fh.write('       66127.0000              0.0069400        \n')
        fh.write('       14815.0000              0.0360900        \n')
        fh.write('       4144.90000              0.1356800        \n')
        fh.write('       1361.20000              0.3387800        \n')
        fh.write(' S   2   1.00\n')
        fh.write('       508.440000              0.4365900        \n')
        fh.write('       209.590000              0.1837500        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       81.9590000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       36.8050000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       13.4950000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       6.88590000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       2.55200000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       1.20880000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       0.27340000              1.0000000        \n')
        fh.write(' S   1   1.00\n')
        fh.write('       0.10090000              1.0000000        \n')
        fh.write(' P   4   1.00\n')
        fh.write('       2953.60000              0.0122100        \n')
        fh.write('       712.610000              0.0858700        \n')
        fh.write('       236.710000              0.2949300        \n')
        fh.write('       92.6310000              0.4784900        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       39.7320000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       17.2730000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       7.95700000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       3.15290000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       1.33280000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       0.49470000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       0.21600000              1.0000000        \n')
        fh.write(' P   1   1.00\n')
        fh.write('       0.08293000              1.0000000        \n')
        fh.write(' D   3   1.00\n')
        fh.write('       261.950000              0.0314400        \n')
        fh.write('       76.7340000              0.1902800        \n')
        fh.write('       27.5510000              0.4724700        \n')
        fh.write(' D   1   1.00\n')
        fh.write('       10.6060000              1.0000000        \n')
        fh.write(' D   1   1.00\n')
        fh.write('       3.42170000              1.0000000        \n')
        fh.write(' D   1   1.00\n')
        fh.write('       1.13700000              1.0000000        \n')
        fh.write(' D   1   1.00\n')
        fh.write('       0.30200000              1.0000000    \n')
        fh.write(' ****\n\n')
        fh.close()

def is_qm_normal_termination(logfname):
    """
    Intent: Checks the *.log file for normal termination
    """
    if os.path.isfile(logfname):
        for line in open(logfname):
            if "Normal termination" in line:
                logfh.write("Normal termination: %s\n" % logfname)
                return True
    return False

def gen_opt_str(optimizeoptlist):
    optstr = "#P opt"
    if optimizeoptlist:
        optstr += "=(" + ','.join(optimizeoptlist) + ")"
    return optstr

def gen_function_specific_restraints(mol):
    """
    Intent: Specify restraints for specific functional groups.
    Input:
        mol: An openbabel molecule structure
    Output:
        restraintlist: a list atom indices and its restraint value.
    Referenced By: tor_opt_sp
    Description: -
    """

    restraintlist = []
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,"CC(=O)[OH][#1]")
    sp.Match(mol)
    for ia in sp.GetUMapList():
        restraintlist.append((ia[1],ia[2],1.250))
        restraintlist.append((ia[1],ia[3],1.250))
        restraintlist.append((ia[3],ia[4],1.041))
        restraintlist.append((ia[0],ia[1],ia[2],116.166))
        restraintlist.append((ia[0],ia[1],ia[3],116.166))
        restraintlist.append((ia[1],ia[3],ia[4],113.709))
    return restraintlist

def append_restraint(restraintlist,comfname):
    """
    Intent: Add restraint terms (if they exist) to *.com file
    Referenced By: tor_opt_sp
    """
    fh = open(comfname, "a")
    if restraintlist:
        for restraint in restraintlist:
            for ele in restraint:
                if isinstance(ele,int):
                    fh.write("%4d" % ele)
                elif isinstance(ele,float):
                    fh.write("%9.3f" % ele)
            fh.write(" F\n")

def write_com_header(comfname,chkfname):
    """
    Intent: Add header to *.com file
    Referenced By: gen_optcomfile
    """
    tmpfh = open(comfname, "w")
    assert tmpfh, "Cannot create file: " + comfname

    tmpfh.write('%RWF=' + scrtmpdir + '/,' + maxdisk + '\n')
#   tmpfh.write('%Int=' + scrtmpdir + '/,' + maxdisk + '\n')
#   tmpfh.write('%D2E=' + scrtmpdir + '/,' + maxdisk + '\n')
    tmpfh.write("%Nosave\n")
    tmpfh.write("%Chk=" + os.path.splitext(comfname)[0] + ".chk\n")
    tmpfh.write("%Mem=" + maxmem + "\n")
    tmpfh.write("%Nproc=" + str(numproc) + "\n")
    tmpfh.close()

def gen_optcomfile (comfname,numproc,maxmem,chkname,mol):
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
    write_com_header(comfname,chkname)
    tmpfh = open(comfname, "a")
    optimizeoptlist = ["maxcycle=%s"%(optmaxcycle)]
    if restraintlist:
        optimizeoptlist.insert(0,"modred")
    optstr=gen_opt_str(optimizeoptlist)
    if ('I ' in mol.GetSpacedFormula()):
        tmpfh.write("%s HF/Gen freq Guess=INDO MaxDisk=%s\n" % (optstr,maxdisk))
    else:
        if freq==True:
            if optpcm==True:
                tmpfh.write("%s %s/%s freq Guess=INDO MaxDisk=%s SCRF=(PCM)\n" % (optstr,optmethod,optbasisset,maxdisk))
            else:
                tmpfh.write("%s %s/%s freq Guess=INDO MaxDisk=%s\n" % (optstr,optmethod,optbasisset,maxdisk))
        else:
            if optpcm==True:
                tmpfh.write("%s %s/%s Guess=INDO MaxDisk=%s SCRF=(PCM)\n" % (optstr,optmethod,optbasisset,maxdisk))
            else:
                tmpfh.write("%s %s/%s Guess=INDO MaxDisk=%s\n" % (optstr,optmethod,optbasisset,maxdisk))
     
    commentstr = molecprefix + " Gaussian SP Calculation on " + gethostname()
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d\n' % (mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity()))
    iteratom = openbabel.OBMolAtomIter(mol)
    etab = openbabel.OBElementTable()
    for atm in iteratom:
        tmpfh.write('%2s %11.6f %11.6f %11.6f\n' % (etab.GetSymbol(atm.GetAtomicNum()), atm.x(), atm.y(), atm.z()))
    tmpfh.write('\n')

    tmpfh.close()
    #NOTE: Restraints need to be specified before
    #       extended basis sets for Gaussian 03/09.
    append_restraint(restraintlist,comfname)
    tmpfh = open(comfname, "a")
    tmpfh.write("\n")
    tmpfh.close()
    append_basisset(comfname, mol.GetSpacedFormula(), optbasisset)

    tempname=comfname.replace('.com','temp.com')
    temp=open(comfname,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempname,'w')
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx!=9:
            temp.write(line)
        else:
            temp.write('%d %d\n' % (mol.GetTotalCharge(), 1))

    temp.close()
    os.remove(comfname)
    os.rename(tempname,comfname)
    if gencomonly==True:
        sys.exit()

def gen_comfile (comfname,numproc,maxmem,chkname,tailfname,mol):
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
    title = "\"" + molecprefix + " Gaussian SP Calculation on " + gethostname() + "\""
    cmdstr = babelexe + " --title " + title + " -i g03 " + gausoptfname + " " + tailfname
    call_subsystem(cmdstr)

    write_com_header(comfname,chkname)
    tmpfh = open(comfname, "a")
    #NOTE: Need to pass parameter to specify basis set
    if ('dma' in comfname):
        opstr="#P MP2/%s Sp Density=MP2 MaxDisk=%s\n" % (dmabasisset, maxdisk)
    elif ('pop' in comfname):
        opstr="#P HF/%s MaxDisk=%s Pop=SaveMixed\n" % (popbasisset, maxdisk)
    else:
        if sppcm==True:
            opstr="#P %s/%s Sp Density=MP2 SCF=Save Guess=Huckel MaxDisk=%s SCRF=(PCM)\n" % (espmethod,espbasisset, maxdisk)
        else:
            opstr="#P %s/%s Sp Density=MP2 SCF=Save Guess=Huckel MaxDisk=%s\n" % (espmethod,espbasisset, maxdisk)


    bset=re.search('(?i)(6-31|aug-cc)\S+',opstr)
    if ('I ' in mol.GetSpacedFormula()):
        opstr=re.sub(r'(?i)(6-31|aug-cc)\S+',r'Gen',opstr)
    tmpfh.write(opstr)
    tmpfh.close()
    cmdstr = "tail -n +2 " + tailfname + " | head -n 4 >> " + comfname
    os.system(cmdstr)
    cmdstr = "tail -n +6 " + tailfname + " | sed -e's/ .*//' > tmp1.txt"
    os.system(cmdstr)
    cmdstr = "grep -A " + str(mol.NumAtoms()+4) + " 'Standard orientation' " + gausoptfname + " | tail -n " + str(mol.NumAtoms()) + " | sed -e's/^.* 0 //' > tmp2.txt"
    os.system(cmdstr)
    cmdstr = "paste tmp1.txt tmp2.txt >> " + comfname
    os.system(cmdstr)
    if ('dma' in comfname):
        append_basisset(comfname, mol.GetSpacedFormula(), dmabasisset)
    else:
        append_basisset(comfname, mol.GetSpacedFormula(), espbasisset)

def gen_torcomfile (comfname,numproc,maxmem,prevstruct,xyzf):
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
    write_com_header(comfname,os.path.splitext(comfname)[0] + ".chk")
    tmpfh = open(comfname, "a")

    optimizeoptlist = ["modred"]
    optimizeoptlist.append("maxcycle=400")
    optstr=gen_opt_str(optimizeoptlist)

    if ('-opt-' in comfname):
        if toroptpcm==True:
            print('actaully here ')
            operationstr = "%s %s/%s MaxDisk=%s SCRF=(PCM)\n" % (optstr,toroptmethod,toroptbasisset, maxdisk)
        else:
            print('right here')
            operationstr = "%s %s/%s MaxDisk=%s\n" % (optstr,toroptmethod,toroptbasisset, maxdisk)
        commentstr = molecprefix + " Rotatable Bond Optimization on " + gethostname()
    else:
#        operationstr = "#m06L/%s SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=%s\n" % (torspbasisset, maxdisk)
        if torsppcm==True:
            operationstr = "#P %s/%s SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=%s SCRF=(PCM)\n" % (torspmethod,torspbasisset, maxdisk)
        else:       
            operationstr = "#P %s/%s SP SCF=(qc,maxcycle=800) Guess=Indo MaxDisk=%s\n" % (torspmethod,torspbasisset, maxdisk)

        commentstr = molecprefix + " Rotatable Bond SP Calculation on " + gethostname()   

    bset=re.search('6-31\S+',operationstr)
    if ('I ' in mol.GetSpacedFormula()):
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
        tmpfh.close()

def CheckNormalTermPsi4Opt(logname):
    finished=False
    if not os.path.isfile(logname):
        return finished
    temp=open(logname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if "Final optimized geometry" in line or 'Psi4 exiting successfully' in line:
            finished=True
    return finished 

def CheckNormalTermPsi4ESP(logname):
    finished=False
    if not os.path.isfile(logname):
        return finished
    temp=open(logname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line:
            finished=True
    return finished 

def run_psi4(mol):
    """
    Intent: QM calculations are done within this method. 
        The calculations are skipped if log files exist, 
        suggesting that QM calculations have already been done.
    Input: 
        mol: OBMol object with pre optimized geometry
    Output: 
        optmol: OBMol object with post gaussian optimized geometry
    Referenced By: main
    Description: 
    1. Babel is used to convert the input file ('molstructfname') to a tmp 
       .com file ('comtmp')
    2. A temporary scratch directory is created for Gaussian if one doesn't already exist
    3. Molecule is optimized using Gaussian
        b. 'molstructfname' is loaded in through the 'load_structfile' method to create
           an OBMol object; thrun_psi4is object is stored as 'mystruct'
        c. an opt .com file is created ('comoptfname') using the 'gen_optcomfile' method;
           basically this method just adds the proper header to the tmp .com file ('comtmp')
        d. Gaussian is run
        e. formchk utility is used to convert the *-opt.chk file to a formatted *-opt.fchk file
        f. the information from the optimization logfile is loaded in to create the OBMol object
           'optmol'
        g. bond information that was in the 'mol' object is given to the 'optmol' object 
           through the method 'rebuild_bonds'
    4. Gaussian is run using the Density=MP2 keyword to find the electron density matrix 
       The density matrix info is in *-dma.fchk
    5. Gaussian is run using the 'Density=MP2 SCF=Save Guess=Huckel' keywords to
       find information that will be used to find the electrostatic potential grid
    """
    global molecprefix
    global numproc
    global maxmem
    global optbasisset
    global dmabasisset
    global gausoptfname
    global comoptfname
    global chkoptfname
    global fckoptfname
    global logoptfname
    global compopfname
    global chkpopfname
    global fckpopfname
    global logpopfname
    global comdmafname
    global chkdmafname
    global fckdmafname
    global logdmafname
    global comespfname
    global chkespfname
    global fckespfname
    global logespfname

    logfh.write("NEED QM Density Matrix: Executing Gaussian Opt and SP\n")
    comtmp = assign_filenames ( "comtmp" , "-tmp.com")
    comoptfname = assign_filenames ( "comoptfname" , "-opt.com")
    chkoptfname = assign_filenames ( "chkoptfname" , "-opt.chk")
    fckoptfname = assign_filenames ( "fckoptfname" , "-opt.fchk")
    logoptfname = assign_filenames ( "logoptfname" , "-opt.log")
    compopfname = assign_filenames ( "compopfname" , "-pop.com")
    chkpopfname = assign_filenames ( "chkpopfname" , "-pop.chk")
    fckpopfname = assign_filenames ( "fckpopfname" , "-pop.fchk")
    logpopfname = assign_filenames ( "logpopfname" , "-pop.log")
    comdmafname = assign_filenames ( "comdmafname" , "-dma.com")
    chkdmafname = assign_filenames ( "chkdmafname" , "-dma.chk")
    fckdmafname = assign_filenames ( "fckdmafname" , "-dma.fchk")
    logdmafname = assign_filenames ( "logdmafname" , "-dma.log")
    comespfname = assign_filenames ( "comespfname" , "-esp.com")
    chkespfname = assign_filenames ( "chkespfname" , "-esp.chk")
    fckespfname = assign_filenames ( "fckespfname" , "-esp.fchk")
    logespfname = assign_filenames ( "logespfname" , "-esp.log")

    title = "\"" + molecprefix + " Gaussian Optimization on " + gethostname() + "\""
    cmdstr = babelexe + " --title " + title + " "+ molstructfname+ " " + comtmp
    call_subsystem(cmdstr)

    assert os.path.getsize(comtmp) > 0, "Error: " + \
       os.path.basename(babelexe) + " cannot create .com file."

    tempname=comtmp.replace('.com','temp.com')
    temp=open(comtmp,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempname,'w')
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx!=4:
            temp.write(line)
        else:
            temp.write('%d %d\n' % (mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity()))

    temp.close()
    os.remove(comtmp)
    os.rename(tempname,comtmp)
    if not os.path.isdir(scrtmpdir):
        os.mkdir(scrtmpdir)

    if use_psi4SPonly==True: # try to use gaussian for opt


        if not is_qm_normal_termination(logoptfname):
            mystruct = load_structfile(molstructfname)
            if os.path.isfile(chkoptfname):
                os.remove(chkoptfname)
            gen_optcomfile(comoptfname,numproc,maxmem,chkoptfname,mol)

            cmdstr = 'GAUSS_SCRDIR=' + scrtmpdir + ' ' + gausexe + " " + comoptfname
            call_subsystem(cmdstr,iscritical=True)
            cmdstr = formchkexe + " " + chkoptfname
            call_subsystem(cmdstr)

        optmol =  load_structfile(logoptfname)
        rebuild_bonds(optmol,mol)
        gen_comfile(comdmafname.replace('.com','_temp.com'),numproc,maxmem,chkdmafname,comtmp,mol)
        
    else:
        gen_optcomfile(comoptfname,numproc,maxmem,chkoptfname,mol)
        if os.path.isfile(chkdmafname):
            os.remove(chkdmafname)

        now = time.strftime("%c",time.localtime())
        logfh.write(now + " Calling: " + "Psi4 Optimization" + "\n")
        finished=CheckNormalTermPsi4Opt(logoptfname)
        inputname=CreatePsi4OPTInputFile(comoptfname,comoptfname,mol)
        if finished==False:
            cmdstr='psi4 '+inputname+' '+logoptfname
            call_subsystem(cmdstr,iscritical=True)
        GrabFinalPsi4XYZStructure(logoptfname,'optimized.xyz')

        optmol =  load_structfile('optimized.xyz')
        rebuild_bonds(optmol,mol)

    import psi4
    #
    # DMA run
    #
    now = time.strftime("%c",time.localtime())
    logfh.write(now + " Calling: " + "Psi4 Gradient for DMA" + "\n")
    finished=CheckNormalTermPsi4ESP(logdmafname)
    if use_psi4SPonly:
        save_structfileXYZ(optmol, 'optimized.xyz')
    inputname=CreatePsi4DMAInputFile('optimized.xyz',comdmafname,mol)
    if finished==False:
        now = time.strftime("%c",time.localtime())
        logfh.write(now + " Calling: " + "Psi4 GDMA for DMA" + "\n")
        cmdstr='psi4 '+inputname+' '+logdmafname
        call_subsystem(cmdstr,iscritical=True)
      

    return optmol

def run_gaussian(mol):
    """
    Intent: QM calculations are done within this method. 
        The calculations are skipped if log files exist, 
        suggesting that QM calculations have already been done.
    Input: 
        mol: OBMol object with pre optimized geometry
    Output: 
        optmol: OBMol object with post gaussian optimized geometry
    Referenced By: main
    Description: 
    1. Babel is used to convert the input file ('molstructfname') to a tmp 
       .com file ('comtmp')
    2. A temporary scratch directory is created for Gaussian if one doesn't already exist
    3. Molecule is optimized using Gaussian
        b. 'molstructfname' is loaded in through the 'load_structfile' method to create
           an OBMol object; this object is stored as 'mystruct'
        c. an opt .com file is created ('comoptfname') using the 'gen_optcomfile' method;
           basically this method just adds the proper header to the tmp .com file ('comtmp')
        d. Gaussian is run
        e. formchk utility is used to convert the *-opt.chk file to a formatted *-opt.fchk file
        f. the information from the optimization logfile is loaded in to create the OBMol object
           'optmol'
        g. bond information that was in the 'mol' object is given to the 'optmol' object 
           through the method 'rebuild_bonds'
    4. Gaussian is run using the Density=MP2 keyword to find the electron density matrix 
       The density matrix info is in *-dma.fchk
    5. Gaussian is run using the 'Density=MP2 SCF=Save Guess=Huckel' keywords to
       find information that will be used to find the electrostatic potential grid
    """
    global molecprefix
    global numproc
    global maxmem
    global gausoptfname
    global comoptfname
    global chkoptfname
    global fckoptfname
    global logoptfname
    global compopfname
    global chkpopfname
    global fckpopfname
    global logpopfname
    global comdmafname
    global chkdmafname
    global fckdmafname
    global logdmafname
    global comespfname
    global chkespfname
    global fckespfname
    global logespfname

    logfh.write("NEED QM Density Matrix: Executing Gaussian Opt and SP\n")
    comtmp = assign_filenames ( "comtmp" , "-tmp.com")
    comoptfname = assign_filenames ( "comoptfname" , "-opt.com")
    chkoptfname = assign_filenames ( "chkoptfname" , "-opt.chk")
    fckoptfname = assign_filenames ( "fckoptfname" , "-opt.fchk")
    logoptfname = assign_filenames ( "logoptfname" , "-opt.log")
    compopfname = assign_filenames ( "compopfname" , "-pop.com")
    chkpopfname = assign_filenames ( "chkpopfname" , "-pop.chk")
    fckpopfname = assign_filenames ( "fckpopfname" , "-pop.fchk")
    logpopfname = assign_filenames ( "logpopfname" , "-pop.log")
    comdmafname = assign_filenames ( "comdmafname" , "-dma.com")
    chkdmafname = assign_filenames ( "chkdmafname" , "-dma.chk")
    fckdmafname = assign_filenames ( "fckdmafname" , "-dma.fchk")
    logdmafname = assign_filenames ( "logdmafname" , "-dma.log")
    comespfname = assign_filenames ( "comespfname" , "-esp.com")
    chkespfname = assign_filenames ( "chkespfname" , "-esp.chk")
    fckespfname = assign_filenames ( "fckespfname" , "-esp.fchk")
    logespfname = assign_filenames ( "logespfname" , "-esp.log")

    title = "\"" + molecprefix + " Gaussian Optimization on " + gethostname() + "\""
    cmdstr = babelexe + " --title " + title + " "+ molstructfname+ " " + comtmp
    call_subsystem(cmdstr)
    time.sleep(2)

    assert os.path.getsize(comtmp) > 0, "Error: " + \
       os.path.basename(babelexe) + " cannot create .com file."

    tempname=comtmp.replace('.com','temp.com')
    temp=open(comtmp,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempname,'w')
    for lineidx in range(len(results)):
        line=results[lineidx]
        if lineidx!=4:
            temp.write(line)
        else:
            temp.write('%d %d\n' % (mol.GetTotalCharge(), 1))

    temp.close()
    os.remove(comtmp)
    os.rename(tempname,comtmp)

    if not os.path.isdir(scrtmpdir):
        os.mkdir(scrtmpdir)

    if not is_qm_normal_termination(logoptfname):
        mystruct = load_structfile(molstructfname)
        if os.path.isfile(chkoptfname):
            os.remove(chkoptfname)
        gen_optcomfile(comoptfname,numproc,maxmem,chkoptfname,mol)
        cmdstr = 'GAUSS_SCRDIR=' + scrtmpdir + ' ' + gausexe + " " + comoptfname
        call_subsystem(cmdstr,iscritical=True)
        cmdstr = formchkexe + " " + chkoptfname
        call_subsystem(cmdstr)
    optmol =  load_structfile(logoptfname)
    rebuild_bonds(optmol,mol)

    if not is_qm_normal_termination(logdmafname):
        if os.path.isfile(chkdmafname):
            os.remove(chkdmafname)
        gen_comfile(comdmafname,numproc,maxmem,chkdmafname,comtmp,mol)
        cmdstr = 'GAUSS_SCRDIR=' + scrtmpdir + ' ' + gausexe + " " + comdmafname
        call_subsystem(cmdstr,iscritical=True)
        cmdstr = formchkexe + " " + chkdmafname
        call_subsystem(cmdstr,iscritical=True)

    if espfit and not is_qm_normal_termination(logespfname):
        if os.path.isfile(chkespfname):
            os.remove(chkespfname)
        gen_comfile(comespfname,numproc,maxmem,chkespfname,comtmp,mol)
        cmdstr = 'GAUSS_SCRDIR=' + scrtmpdir + ' ' + gausexe + " " + comespfname
        call_subsystem(cmdstr,iscritical=True)
        cmdstr = formchkexe + " " + chkespfname
        call_subsystem(cmdstr)

    return optmol

def gen_canonicallabels(mol):
    global symmetryclass
    symmetryclass = [ 0 ] * mol.NumAtoms()
    """
    Intent: Find the symmetry class that each atom belongs to
    Input: 
        mol: OBMol object 
    Output: 
        The global variable 'symmetryclass' is altered
    Referenced By: main
    Description:
    1. An empty bit vector is created, 'frag_atoms'
    2. OBMol.FindLargestFragment is called to fill in the 'frag_atoms' bit vector (the
    vector is filled with a 1 or 0 depending on whether the atom is part of the largest
    fragment or not)
    3. 'CalculateSymmetry' method is called to find initial symmetry classes
    4. Terminal atoms of the same element are collapsed to one symmetry class
    5. Possibly renumber the symmetry classes
    """
    # Returns symmetry classes for each atom ID
    frag_atoms = openbabel.OBBitVec()
    symmclasslist = []
    mol.FindLargestFragment(frag_atoms)
    CalculateSymmetry(mol, frag_atoms, symmclasslist)
    for ii in range(len(symmetryclass)):
        #print(len(symmetryclass),ii)
        symmetryclass[ii] = symmclasslist[ii][1]

    # Collapse terminal atoms of same element to one type
    for a in openbabel.OBMolAtomIter(mol):
        for b in openbabel.OBAtomAtomIter(a):
            if b.GetValence() == 1:
                for c in openbabel.OBAtomAtomIter(a):
                    if ((b is not c) and
                        (c.GetValence() == 1) and
                        (b.GetAtomicNum() == c.GetAtomicNum()) and
                        (symmetryclass[b.GetIdx()-1] !=
                            symmetryclass[c.GetIdx()-1])):
                        symmetryclass[c.GetIdx()-1] = \
                            symmetryclass[b.GetIdx()-1]

    # Renumber symmetry classes
    allcls=list(set(symmetryclass))
    allcls.sort()
    for ii in range(len(symmetryclass)):
        symmetryclass[ii] = allcls.index(symmetryclass[ii]) + 1


#scaling of multipole values for certain atom types
def process_types (mol):
    """
    Intent: Set up scalelist array for scaling certain multipole values 
    """
    scalelist = {}
    for atm in openbabel.OBMolAtomIter(mol):
        if get_class_number(atm.GetIdx()) not in scalelist:
            scalelist[get_class_number(atm.GetIdx())] = []
            scalelist[get_class_number(atm.GetIdx())].append(None)
            scalelist[get_class_number(atm.GetIdx())].append(None)
            scalelist[get_class_number(atm.GetIdx())].append(None)
            multipole_scale_dict = {}

    for (sckey, scval) in multipole_scale_dict.items():
        sp = openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,sckey)
        sp.Match(mol)
        for ia in sp.GetUMapList():
            scalelist[get_class_number(ia[0])][scval[0]] = scval[1]

    return scalelist


def gen_gdmain(gdmainfname,molecprefix,fname):
    """
    Intent: Generate GDMA input file for the molecule
    Input: 
        gdmainfname: file name for the gdma input file
        molecprefix: name of the molecule
        fname: file name for the *-dma.fchk file
    Output: *.gdmain file is created
    Referenced By: run_gdma
    Description:
        1. create pointer file dma.fchk to *-dma.fchk
        2. create *.gdmain file and write in all the necessary information for the gdma run
    """
    fnamesym = "dma.fchk"
    try:
        os.symlink(fname,fnamesym)
    except Exception as e: 
        print(e)
        
        
    dmamethod='MP2'
    #punfname = os.path.splitext(fname)[0] + ".punch"
    punfname = "dma.punch"

    try:
        tmpfh = open(gdmainfname, "w")
    except(IOError, (errno, strerror)):
        print("I/O error({0}): {1}".format(errno, strerror))
        sys.exit(errno)

    tmpfh.write("Title " + molecprefix + " gdmain\n")
    tmpfh.write("\n")
    if use_psi4 or use_psi4SPonly:
        if dmamethod=='MP2':
            dmamethod='CC' # for some reason fchk outputs CC for MP2 density
        tmpfh.write("File " + fnamesym  + " density %s\n"%(dmamethod))
    else:
        tmpfh.write("File " + fnamesym  + " density %s\n"%(dmamethod))
    tmpfh.write("Angstrom\n")
    tmpfh.write("AU\n")
    tmpfh.write("Multipoles\n")
    tmpfh.write("Switch 0\n") # v1.3, comment out for v2.2
    tmpfh.write("Limit 2\n")
    tmpfh.write("Punch " + punfname + "\n")
    tmpfh.write("Radius H 0.65\n") # v1.3, use 0.31 for v2.2
    tmpfh.write("Radius S 0.80\n") # v1.3, use 0.XX for v2.2
    tmpfh.write("Radius P 0.75\n") # v1.3, use 0.XX for v2.2
    tmpfh.write("Radius Cl 1.10\n") # v1.3, use 0.XX for v2.2
    tmpfh.write("\n")
    tmpfh.write("Start\n")
    tmpfh.write("\n")
    tmpfh.write("Finish\n")
    tmpfh.close()

def run_gdma():
    """
    Intent: Runs GDMA to find multipole information
    The GDMA program carries out distributed multipole analysis of the wavefunctions
    calculated by Gaussian (using the *-dma.fchk file)
    Input:
    Output: *.gdmaout is created; this is used as input by poledit, a tinker program 
    Referenced By: main
    Description: 
    1. Generates the gdma input file by calling 'gen_gdmain'
    2. Runs the following command: gdma < *.gdmain > *.gdmaout
    """
    global molecprefix
    global fckdmafname

    logfh.write("NEED DMA: Executing GDMA\n")

    if not os.path.isfile(fckdmafname):
        fckdmafname = os.path.splitext(fckdmafname)[0]

    assert os.path.isfile(fckdmafname), "Error: " + fckdmafname + " does not exist."
    gdmainfname = assign_filenames ( "gdmainfname" , ".gdmain")
    gen_gdmain(gdmainfname,molecprefix,fckdmafname)

    cmdstr = gdmaexe + " < " + gdmainfname + " > " + gdmafname
    call_subsystem(cmdstr)

    
    
    while not os.path.isfile(gdmafname):
        time.sleep(1)

    assert os.path.getsize(gdmafname) > 0, "Error: " + \
       os.path.basename(gdmaexe) + " cannot create .gdmaout file."

def is_in_polargroup(mol, smarts, bond, f):
    """
    Intent: Check if a given bond is in the group defined by the smarts string
    Input:
        mol: OBMol object
        smarts: Smarts String defining the group (e.g. O-[*])
        bond: bond in question
        f: output file (not used)
    Output:
        True if the bond is in the group, False if not
    Referenced By: gen_peditinfile
    Description: -
    """
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    sp.Match(mol)
    for i in sp.GetUMapList():
        if ((bond.GetBeginAtomIdx() in i) and \
             (bond.GetEndAtomIdx() in i)):
            return True
    return False

# Create file to define local frames for multipole
def gen_peditinfile (mol):
    """
    Intent: Create a file with local frame definitions for each multipole
    These frame definitions are given as input into tinker's poledit
    Input:
        mol: OBMol object
    Output: *-peditin.txt is created
    Referenced By: main
    Description:
    1. Initialize lfzerox, or local-frame-zero-x-component, array
    2. For each atom a
            a. find the two heaviest atoms that a is bound to
               (heaviest defined by their symmetry class; 
               the more atoms an atom is bound to, the higher its symmetry class)
               These two atoms define the local frame for atom a
            b. This information is stored in the localframe1 and localframe2 arrays
    3. If atom a is bound to two atoms of the same symm class (that aren't Hydrogens),
       use these two atoms to define the frame
    4. Find the atoms that are only bound to one atom and define their local frame 
       based on the local frame of the one atom they are bound to. 
    5. Zero out the x-component of the local frame if there is more than one choice for the 
       x-component
       Note: If the atom is only bound to one atom, then the array lfzerox is altered
       If the atom is bound to more than one atom, then the array lf2write is altered
       Related to how poledit handles the local frames of atoms with valence 1
    6. Define bisectors; i.e. local frames where both components belong to the same sym class
    7. Write out frames to *-peditin.txt
    8. Write out polarizabilities for certain atom types
    9. Define polarizability groups by 'cutting' certain bonds
       The openbabel method, IsRotor() is used to decide whether a bond is cut or not
    """
    lfzerox = [ False ] * mol.NumAtoms()
    atomindextoremovedipquad={} # for methane need to make dipole and quadupole on the carbon zeroed out, will return this for post proccesing the keyfile after poledit is run
    atomindextoremovedipquadcross={}
    atomtypetospecialtrace={} # for H on CH4 need to make sure Qxx=Qyy=-1/2*Qzz
    # Assign local frame for atom a based on the symmetry classes of the atoms it is bound to
    for a in openbabel.OBMolAtomIter(mol):
        # iterate over the atoms that a is bound to
        for b in openbabel.OBAtomAtomIter(a):
            lf1 = localframe1[a.GetIdx() - 1] # defined to be zero at first?
            lf2 = localframe2[a.GetIdx() - 1] # defined to be zero at first?
            # Sort list based on symmetry class, largest first
            a1 = sorted((lf1, lf2, b.GetIdx()), key=get_symm_class, reverse=True) # the lf1 and lf2 return -1 via get_symm_class
            while a1[0] == 0: # redundant if you got the order (reverse or not ) in previous line
                a1 = rotate_list(a1)
            # Set localframe1 and localframe2 for atom a to be the first two atoms
            # of the above sorted list 'a1'
            localframe1[a.GetIdx() - 1] = a1[0] # highest neighboring symmetry class gets stored as local frame 1
            localframe2[a.GetIdx() - 1] = a1[1] # next highest is 2

    # if a is bound to two atoms of the same symmetry class that aren't hydrogens
    # use these two atoms to define the local frame
    for a in openbabel.OBMolAtomIter(mol): # it seems that this case is not handled by above case
        classlist = {}
        for b in openbabel.OBAtomAtomIter(a):
            if b.GetAtomicNum() != 1:
                clsidx = symmetryclass[b.GetIdx() - 1]
                if clsidx not in classlist:
                    classlist[clsidx] = []
                classlist[clsidx].append(b.GetIdx())
        for clstype in classlist.values():
            if len(clstype) > 1:
                localframe1[a.GetIdx() - 1] = clstype[0]
                localframe2[a.GetIdx() - 1] = clstype[1]

    # Find atoms bonded to only one atom
    iteratom = openbabel.OBMolAtomIter(mol)
    for a in iteratom:
        lfa1 = localframe1[a.GetIdx() - 1]
        lfa2 = localframe2[a.GetIdx() - 1]
        lfb1 = localframe1[lfa1 - 1]
        lfb2 = localframe2[lfa1 - 1]
        if a.GetValence() == 1:
            # Set lfa2 to the other atom (the atom that isn't 'a') in the local frame of atom 'b'
            if lfb1 != a.GetIdx():
                localframe2[a.GetIdx() - 1] = lfb1
            else:
                localframe2[a.GetIdx() - 1] = lfb2

    # Zero out x-component if more than one possible choice
    # for x-axis local frame
    # Note: Zeroed x-components need to be added to LF2 after writing to
    # poledit in file.

    # lf2write is to write out localframe2 for poledit
    # Since poledit will zero out x-comp if x-axis(i) = 0
    lf2write = list(localframe2)
    for a in openbabel.OBMolAtomIter(mol):
        lf1 = localframe1[a.GetIdx() - 1]
        lf2 = localframe2[a.GetIdx() - 1]

        # Check LF1 instead if LF2 is not bonded to atom "a"
        # If a is only bound to one other atom
        if a.GetValence() == 1:
            center = mol.GetAtom(localframe1[a.GetIdx() - 1])
            for b in openbabel.OBAtomAtomIter(center):
                # lf2 and b are not the same atom, but the same sym class, set lfzerox to true
                if (lf2 != b.GetIdx() and
                    (get_symm_class(lf2) == get_symm_class(b.GetIdx()))):
                    pass
                    #lfzerox[a.GetIdx() - 1] = True
        else:
            for b in openbabel.OBAtomAtomIter(a):
                # if lf2 and b are not the same atom, but the same sym class, set lf2write to 0. Handles analine case. 
                if (lf1 != b.GetIdx() and
                    lf2 != b.GetIdx() and
                    (get_symm_class(lf2) == get_symm_class(b.GetIdx()))):
                    lf2write[a.GetIdx() - 1] = 0

    # Define bisectors
    for a in openbabel.OBMolAtomIter(mol):
        lfa1 = localframe1[a.GetIdx() - 1]
        lfa2 = localframe2[a.GetIdx() - 1]
        if (a.IsConnected(mol.GetAtom(lfa1)) and
            a.IsConnected(mol.GetAtom(lfa2)) and
            get_symm_class(lfa1) == get_symm_class(lfa2)):
            localframe2[a.GetIdx() - 1] *= -1
            lf2write[a.GetIdx() - 1] *= -1

    print('symmetryclass',symmetryclass)
    for a in openbabel.OBMolAtomIter(mol):
        if a.GetValence() == 4:
            aidx=a.GetIdx()
            atype=get_symm_class(aidx)
            neighbtypelist=[]
            neighbidxlist=[]
            for b in openbabel.OBAtomAtomIter(a):
                bidx=b.GetIdx()
                btype=get_symm_class(bidx)
                if btype not in neighbtypelist:
                    neighbtypelist.append(btype)
                    neighbidxlist.append(bidx)
            if len(neighbtypelist)==1 and atype not in neighbtypelist: # then this is like Methane, one type surrounded by four of the same type (different then middle)
                atomindextoremovedipquad[aidx]=True
                atomtypetospecialtrace[bidx]=True
                
                lf2write[b.GetIdx() - 1] = 0
                for b in openbabel.OBAtomAtomIter(a):
                    bidx=b.GetIdx()
                    atomindextoremovedipquadcross[bidx]=True

    #iterate through molecule and search for case such as OP(O)(O)(O) where three of the atoms have the same symmetry class and the fourth atom has a different class. The middle atom is SP3 and the frame needs to be z-only
    for a in openbabel.OBMolAtomIter(mol):
        if a.GetValence() == 4:
            idxtosymmnum={}
            for b in openbabel.OBAtomAtomIter(a):
                bidx=b.GetIdx()
                symmclass=get_symm_class(bidx)
                idxtosymmnum[bidx]=symmclass
            symset=set(idxtosymmnum.values())
            if len(symset)==2: # only two classes, now see if three are the same
                symcls1=list(symset)[0]
                symcls2=list(symset)[1]
                counter=collections.Counter(idxtosymmnum.values())
                if 3 in counter.values(): # then there must also be a 1 and then this is the case we want
                    for key in counter.keys():
                        count=counter[key]
                        if count==3:
                            pass
                        elif count==1:
                            neededclass=key
                            for idx in idxtosymmnum.keys():
                                symmnum=idxtosymmnum[idx]
                                if symmnum==neededclass:
                                    localframe1[a.GetIdx() - 1]=idx
                                    lf2write[a.GetIdx() - 1] = 0
    idxtobisecthenzbool={}
    idxtobisectidxs={}
    # now iterate through molecule and search for case like  OP(O)(O)(O) where the three Oxygens on the end need bisector then z
    for a in openbabel.OBMolAtomIter(mol):
        idxtobisecthenzbool[a.GetIdx()]=False
        if a.GetValence() == 1:
            for b in openbabel.OBAtomAtomIter(a): # possibly like the middle P atom
                bidx=b.GetIdx()
                if b.GetValence() == 4:
                    idxtosymmnum={}
                    for c in openbabel.OBAtomAtomIter(b):
                        cidx=c.GetIdx()
                        symmclass=get_symm_class(cidx)
                        idxtosymmnum[cidx]=symmclass
                    symset=set(idxtosymmnum.values())
                    if len(symset)==2: # only two classes, now see if three are the same
                        symcls1=list(symset)[0]
                        symcls2=list(symset)[1]
                        counter=collections.Counter(idxtosymmnum.values())
                        if 3 in counter.values(): # then there must also be a 1 and then this is the case we want 
                            for key in counter.keys():
                                count=counter[key]
                                if count==3: # class I want but now grab the other two indexes (not a.GetIdx()) for the bisector, make sure localframe1 is now the P atom(for example)
                                    localframe1[a.GetIdx() - 1]=b.GetIdx()
                                    idxtobisecthenzbool[a.GetIdx()]=True
                                    neededclass=key
                                    bisectidxs=[]
                                    for idx in idxtosymmnum.keys():
                                        symmnum=idxtosymmnum[idx]
                                        if symmnum==neededclass and idx!=a.GetIdx():
                                            bisectidxs.append(idx)
                                    idxtobisectidxs[a.GetIdx()]=bisectidxs

    # now iterate through molecule and search for case N(C)(H)(H), we want this to be a z then bisector N C -H -H
    for a in openbabel.OBMolAtomIter(mol):
        if a.GetValence() == 3 and a.GetAtomicNum()==7:
            hydcount=0
            aliphatcarbcount=0
            hydindexes=[]
            for b in openbabel.OBAtomAtomIter(a):
                bidx=b.GetIdx()
                if b.GetAtomicNum()==1:
                    hydcount+=1
                    hydindexes.append(b.GetIdx())
                if b.GetAtomicNum()==6 and b.IsAromatic()==False:
                    aliphatcarbcount+=1
                    carbidx=b.GetIdx()
            if hydcount==2 and aliphatcarbcount==1:
                localframe1[a.GetIdx() - 1]= carbidx # what?
                idxtobisecthenzbool[a.GetIdx()]=True
                idxtobisectidxs[a.GetIdx()]=hydindexes
                   
    # write out the local frames
    iteratom = openbabel.OBMolAtomIter(mol)
    f = open (peditinfile, 'w')

    for a in iteratom:
        if idxtobisecthenzbool[a.GetIdx()]==False:
            f.write(str(a.GetIdx()) + " " + str(localframe1[a.GetIdx() - 1]) + " " + str(lf2write[a.GetIdx() - 1]) + "\n")
        else:
            bisectidxs=idxtobisectidxs[a.GetIdx()]
            f.write(str(a.GetIdx()) + " " + str(localframe1[a.GetIdx() - 1]) + " -" + str(bisectidxs[0])+ " -" + str(bisectidxs[1]) + "\n")
    
    f.write("\n")
    f.write('A'+'\n')

    #Find aromatic carbon, halogens, and bonded hydrogens to correct polarizability
    iteratom = openbabel.OBMolAtomIter(mol)
    for a in iteratom:
        if (a.GetAtomicNum() == 6 and a.IsAromatic()):
            f.write(str(a.GetIdx()) + " " + str(1.750) + "\n")
        elif (a.GetAtomicNum() == 9):
            f.write(str(a.GetIdx()) + " " + str(0.507) + "\n")
        elif (a.GetAtomicNum() == 17):
            f.write(str(a.GetIdx()) + " " + str(2.500) + "\n")
        elif (a.GetAtomicNum() == 35):
            f.write(str(a.GetIdx()) + " " + str(3.595) + "\n")
        elif (a.GetAtomicNum() == 1):
            iteratomatom = openbabel.OBAtomAtomIter(a)
            for b in iteratomatom:
                if (b.GetAtomicNum() == 6 and b.IsAromatic()):
                    f.write(str(a.GetIdx()) + " " + str(0.696) + "\n")

    # Carboxylate ion O-
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,'[OD1]~C~[OD1]')
    sp.Match(mol)
    for ia in sp.GetMapList():
        f.write(str(ia[0]) + " " + str(0.921) + "\n")

    #f.write("\n")

    #Define polarizable groups by cutting bonds
    iterbond = openbabel.OBMolBondIter(mol)
    for b in iterbond:
        if (b.IsRotor()):
            cut_bond = True
            # If in this group, then don't cut bond
            cut_bond = cut_bond and (not is_in_polargroup(mol,'a[CH2][*]', b,f))
            cut_bond = cut_bond and (not is_in_polargroup(mol,'[#6]O-[#1]', b,f))
            #cut_bond = cut_bond and (not is_in_polargroup(mol,'[#6][#6]~O', b,f))
            # Formamide RC=O
            cut_bond = cut_bond and (not is_in_polargroup(mol,'[*][CH]=O', b,f))
            # Amide
            cut_bond = cut_bond and (not is_in_polargroup(mol,'N(C=O)', b,f))
            cut_bond = cut_bond and (not is_in_polargroup(mol,'C(C=O)', b,f))
            cut_bond = cut_bond and (not is_in_polargroup(mol,'aN', b,f))
            cut_bond = cut_bond and (not is_in_polargroup(mol,'O-[*]', b,f))
            #cut_bond = cut_bond and (not is_in_polargroup(mol,'C[NH2]', b,f))
            if (cut_bond):
                f.write( str(b.GetBeginAtomIdx()) + " " + str(b.GetEndAtomIdx()) + "\n")

    f.write('\n')
    f.write('\n')
    f.write("N\n")
    f.close()
    return lfzerox,atomindextoremovedipquad,atomtypetospecialtrace,atomindextoremovedipquadcross

def CalculateSymmetry(pmol, frag_atoms, symmetry_classes):
    """
    Intent: Uses and builds on openbabel's 'GetGIVector' method which,
    "Calculates a set of graph invariant indexes using the graph theoretical distance,
    number of connected heavy atoms, aromatic boolean,
    ring boolean, atomic number, and summation of bond orders connected to the atom", to
    find initial symmetry classes.
    Input: 
        pmol: OBMol object 
        frag_atoms: OBBitVec object containing information about the largest fragment in 'pmol' 
        symmetry_classes: the symmetry_classes array which will be filled in
    Output: 
		nclasses: # of classes
        symmetry_classes: array is filled in
    Referenced By: gen_canonicallabels
    Description:
    1. vectorUnsignedInt object is created, 'vgi'
    2. It is filled in with the call to GetGIVector
    3. The 'symmetry_classes' array is initially filled in to match with 'vgi'
    4. These initial invariant classes do not suit our needs perfectly,
       so the ExtendInvariants method is called to find the more 
       refined classes that we need
    """

    vgi = openbabel.vectorUnsignedInt()
    natoms = pmol.NumAtoms()
    nfragatoms = frag_atoms.CountBits()
    pmol.GetGIVector(vgi)
    iteratom = openbabel.OBMolAtomIter(pmol)
    for atom in iteratom:
        idx = atom.GetIdx()
        if(frag_atoms.BitIsOn(idx)):
            symmetry_classes.append([atom, vgi[idx-1]])
    nclasses = ExtendInvariants(pmol, symmetry_classes,frag_atoms,nfragatoms,natoms)
    return nclasses

def ExtendInvariants(pmol, symmetry_classes,frag_atoms,nfragatoms,natoms):
    """
    Intent: Refine the invariants found by openbabel's GetGIVector
    Input: 
        pmol: OBMol object
        symmetry_classes: symmetry classes array
        frag_atoms: vector containing information about which atoms belong to the largest fragment
        nfragatoms: # of atoms belonging to the largest fragment
        natoms: # of atoms belonging to the molecule
    Output: 
        nclasses1: # of symmetry classes found 
        symmetry_classes: this array is updated
    Referenced By: CalculateSymmetry
    Description:
    1. Find the # of current classes found by openbabel, nclasses1, 
       and renumber (relabel) the classes to 1, 2, 3, ...
    2. Begin loop
       a. CreateNewClassVector is called which fills in the 'tmp_classes' array with
          a new set of classes by considering bonding information
       b. The number of classes in tmp_classes is found, 'nclasses2'
       c. If there was no change, nclasses1 == nclasses2, break
       d. If the number of classes changed, set nclasses1 to nclasses2, then continue loop
       e. The loop is continued because now that the symmetry classes have changed,
          CreateNewClassVector may find more classes
    3. Return # of classes found
    """
    nclasses1 = CountAndRenumberClasses(symmetry_classes)
    tmp_classes = []
    if(nclasses1 < nfragatoms):
        #stops when number of classes don't change
        for i in range(100):
            CreateNewClassVector(pmol,symmetry_classes, tmp_classes, frag_atoms, natoms)
            nclasses2 = CountAndRenumberClasses(tmp_classes)
            del symmetry_classes[:]
            symmetry_classes.extend(tmp_classes)
            if(nclasses1 == nclasses2):
                break
            nclasses1 = nclasses2
    return nclasses1

def CountAndRenumberClasses(symmetry_classes):
    """
    Intent: Counts the number of symmetry classes and renumbers them to 1, 2, 3, ...
    Input: 
        symmetry_classes: Array of symmetry classes
    Output: 
        count: # of symmetry classes
        symmetry_classes array is updated
    Referenced By: ExtendInvariants
    Description: -
    """
    count = 1
    symmetry_classes = sorted(symmetry_classes, key=lambda sym: sym[1])
    if(len(symmetry_classes) > 0):
        idatom = symmetry_classes[0][1]
        symmetry_classes[0][1] = 1
        for i in range(1,len(symmetry_classes)):
            if(symmetry_classes[i][1] != idatom):
                idatom = symmetry_classes[i][1]
                count = count + 1
                symmetry_classes[i][1] = count
            else:
                symmetry_classes[i][1] = count
    return count

def CreateNewClassVector(pmol,symmetry_classes, tmp_classes, frag_atoms, natoms):
    """
    Intent: Find new symmetry classes if possible
    If two atoms were originally of the same sym class but are bound to atoms of differing
    sym classes, then these two atoms will now belong to two different sym classes
    Input:
        pmol: OBMol object
        symmetry_classes: previous set of symmetry classes
        tmp_classes: tmp array of new symmetry classes
        frag_atoms: atoms in largest fragment
        natoms: number of atoms
    Ouptut:
        tmp_classes is edited
    Referenced By: ExtendInvariants
    Description:
    1. dict idx2index is created which maps atom idx's to an index ranging
       from 0 to # of symmetry classes - 1
    2. For each atom a:
           a. For each bond b that a belongs to:
               i. Find the idx of the other atom, nbratom, in the bond b
               ii. Find and append the symmetry class that nbratom belongs to to vtmp
           b. Using vtmp, create a label for atom a
              i. This label contains information about the symmetry classes of the atoms
              that a is bound to
              ii. This label will be different for two atoms that were originally the same 
              symmetry class but are bound to atoms of differing symmetry classes
    """
    idx2index = dict()
    index = 0
    del tmp_classes[:]
    for s in symmetry_classes:
        idx2index.update({s[0].GetIdx() : index})
        index = index + 1
    for s in symmetry_classes:
        iterbond = openbabel.OBMolBondIter(pmol)
        atom = s[0]
        idatom = s[1]
        nbridx =  0
        vtmp = []
        for b in iterbond:
            #if atom belongs to bond b
            if atom.GetIdx() == b.GetEndAtomIdx() or atom.GetIdx() == b.GetBeginAtomIdx():
                if(atom.GetIdx() == b.GetEndAtomIdx()):
                    nbridx = b.GetBeginAtomIdx()
                elif(atom.GetIdx() == b.GetBeginAtomIdx()):
                    nbridx = b.GetEndAtomIdx()
                if(frag_atoms.BitIsOn(nbridx)):
                    vtmp.append(symmetry_classes[idx2index[nbridx]][1])
        vtmp.sort()
        m = 100
        for v in vtmp:
            idatom = idatom + v * m
            m = 100 * m
        tmp_classes.append([atom,idatom])

def find_tor_restraint_idx(mol,b1,b2):
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
    maxb1class = max(b1nbridx, key=get_symm_class)


    iteratomatom = openbabel.OBAtomAtomIter(b2)
    b2nbridx = list(map(lambda x: x.GetIdx(), iteratomatom))
    del b2nbridx[b2nbridx.index(b1idx)]    # Remove b1 from list
    assert(b2nbridx is not [])
    maxb2class = max(b2nbridx, key=get_symm_class)


    t1 = mol.GetAtom(maxb1class)
    t4 = mol.GetAtom(maxb2class)


    return t1,t4

def get_torlist(mol):
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
    for bond in iterbond:
        # is the bond rotatable
        t2 = bond.GetBeginAtom()
        t3 = bond.GetEndAtom()
        t2idx=t2.GetIdx()
        t3idx=t3.GetIdx()
        if (bond.IsRotor()) or (str(t2idx) in onlyrotbndlist and str(t3idx) in onlyrotbndlist):
            skiptorsion = False
            t1,t4 = find_tor_restraint_idx(mol,t2,t3)
            # is the torsion in toromitlist
            if(omittorsion2 and sorttorsion([t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx()]) in toromit_list):
                skiptorsion = True

            #Check to see if the torsion was found in the look up table or not
            #This decides whether it needs to be scanned for or not
            v1 = valence.Valence(output_format,logfname)
            idxtoclass = []
            for i in range(mol.NumAtoms()):
                idxtoclass.append(i+1)
            v1.setidxtoclass(idxtoclass)
            # dorot is set as false in valence.py
            v1.torguess(mol,False,[])
            missed_torsions = v1.get_mt()
            #pause=input('pausing...')
            if(not sorttorsion([t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx()]) in missed_torsions):
                skiptorsion = True

            if [t2.GetIdx(),t3.GetIdx()] in fitrotbndslist or [t3.GetIdx(),t2.GetIdx()] in fitrotbndslist:
                skiptorsion = False # override previous conditions if in list
            if str(t2idx) in onlyrotbndlist and str(t3idx) in onlyrotbndlist and onlyrotbnd==True:
                skiptorsion = False
            if rotalltors==True:
                skiptorsion=False

            rotbndkey = '%d %d' % (t2.GetIdx(), t3.GetIdx())
            rotbndlist[rotbndkey] = []
            if (not skiptorsion):
                # store the torsion and temporary torsion value found by openbabel in torlist
                tor = mol.GetTorsion(t1,t2,t3,t4)
                torlist.append([t1,t2,t3,t4,tor % 360])
                # store torsion in rotbndlist
                rotbndlist[rotbndkey].append(get_uniq_rotbnd(
                        t1.GetIdx(),t2.GetIdx(),
                        t3.GetIdx(),t4.GetIdx()))
                # write out rotatable bond to log
                logfh.write('Rotatable bond found about %s\n' %
                str(rotbndlist[rotbndkey][0]))
            else:
                continue
            
            #Find other possible torsions about this rotatable bond
            iteratomatom = openbabel.OBAtomAtomIter(bond.GetBeginAtom())
            for iaa in iteratomatom:
                iteratomatom2 = openbabel.OBAtomAtomIter(bond.GetEndAtom())
                for iaa2 in iteratomatom2:
                    a = iaa.GetIdx()
                    b = t2.GetIdx()
                    c = t3.GetIdx()
                    d = iaa.GetIdx()
                    skiptorsion = False
                    if ((iaa.GetIdx() != t3.GetIdx() and \
                             iaa2.GetIdx() != t2.GetIdx()) \
                        and not (iaa.GetIdx() == t1.GetIdx() and \
                             iaa2.GetIdx() == t4.GetIdx())):
                        rotbndlist[rotbndkey].append(get_uniq_rotbnd(
                            iaa.GetIdx(),t2.GetIdx(),
                            t3.GetIdx(),iaa2.GetIdx()))

    return (torlist ,rotbndlist)

def get_torlist2(mol,bondlist,scanlist):
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
    v1 = valence.Valence(output_format,logfname)
    idxtoclass = []
    for i in range(mol.NumAtoms()):
        idxtoclass.append(i+1)
    v1.setidxtoclass(idxtoclass)
    # dorot is set as false in valence.py
    v1.torguess(mol,False,[])
    for bond in iterbond:
        t2 = bond.GetBeginAtom()
        t3 = bond.GetEndAtom()
        # is the bond rotatable
        if (bond.IsRotor())and( [t2.GetIdx(),t3.GetIdx()] in bondlist or [t3.GetIdx(),t2.GetIdx()] in bondlist):
            skiptorsion = False

            t1,t4 = find_tor_restraint_idx(mol,t2,t3)
            # is the torsion in toromitlist
            if(omittorsion2 and sorttorsion([t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx()]) in toromit_list):
                skiptorsion = True
            #Check to see if the torsion was found in the look up table or not
            #This decides whether it needs to be scanned for or not

            test1=([t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx()]) in scanlist
            test2= ([t4.GetIdx(),t3.GetIdx(),t2.GetIdx(),t1.GetIdx()]) in scanlist
            if test1==False and test2==False:
                skiptorsion = True
             
            sortedtor=sorttorsion(tor)
            rotbndkey = '%d %d' % (t2.GetIdx(), t3.GetIdx())
            rotbndlist[rotbndkey] = []
            if (not skiptorsion):
                # store the torsion and temporary torsion value found by openbabel in torlist
                tor = mol.GetTorsion(t1,t2,t3,t4)
                torlist.append([t1,t2,t3,t4,tor % 360])
                # store torsion in rotbndlist
                rotbndlist[rotbndkey].append(get_uniq_rotbnd(
                        t1.GetIdx(),t2.GetIdx(),
                        t3.GetIdx(),t4.GetIdx()))
                # write out rotatable bond to log
                logfh.write('Rotatable bond found about %s\n' %
                str(rotbndlist[rotbndkey][0]))
            else:
                continue
            
            #Find other possible torsions about this rotatable bond
            iteratomatom = openbabel.OBAtomAtomIter(bond.GetBeginAtom())
            for iaa in iteratomatom:
                iteratomatom2 = openbabel.OBAtomAtomIter(bond.GetEndAtom())
                for iaa2 in iteratomatom2:
                    a = iaa.GetIdx()
                    b = t2.GetIdx()
                    c = t3.GetIdx()
                    d = iaa2.GetIdx()
                    skiptorsion = False
                    tor=[a,b,c,d]
                    if ((iaa.GetIdx() != t3.GetIdx() and \
                             iaa2.GetIdx() != t2.GetIdx()) \
                        and not (iaa.GetIdx() == t1.GetIdx() and \
                             iaa2.GetIdx() == t4.GetIdx())):
                        rotbndlist[rotbndkey].append(get_uniq_rotbnd(
                            iaa.GetIdx(),t2.GetIdx(),
                            t3.GetIdx(),iaa2.GetIdx()))

    return (torlist ,rotbndlist)

def get_torlist_opt_angle(optmol, torlist):
    tmplist = []
    for tor in torlist:
        a,b,c,d=obatom2idx(tor[0:4])
        e = optmol.GetTorsion(a,b,c,d)
        tmplist.append([a,b,c,d,e % 360])
    return tmplist

def gen_tinker5_to_4_convert_input(mol, amoeba_conv_spec_fname):
    outfh = open(amoeba_conv_spec_fname,'w')
    class_numH_dict = {}
    for atm in openbabel.OBMolAtomIter(mol):
        num_hydrogen = atm.GetValence() - atm.GetHvyValence()
        if not atm.IsHydrogen() and atm.GetValence() > 1:
            class_numH_dict[get_class_number(atm.GetIdx())] = num_hydrogen

    for (cls, nh) in class_numH_dict.items():
        outfh.write( str(cls) + " " + str(nh) + "\n")

def ConvertTinktoXYZ(filename):
    print(os.getcwd())
    print(time.strftime("%c",time.localtime()))
    temp=open(os.getcwd()+r'/'+filename,'r')
    tempwrite=open(os.getcwd()+r'/'+filename.replace('.xyz_2','-t.xyz'),'w')
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
    return filename.replace('.xyz_2','-t.xyz')


def CreatePsi4TorOPTInputFile(comfilecoords,comfilename,mol,b,c,consttorlist,optmol):
    tempread=open(comfilecoords,'r')
    results=tempread.readlines()
    tempread.close()
    inputname=comfilename.replace('.com','_psi4.dat')
    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (mol.GetTotalCharge(),1))
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==4 and '#' not in line:
            temp.write(line)
    temp.write('}'+'\n')
    # Fix all torsions around the rotatable bond b-c
    temp.write('set optking { '+'\n')
    temp.write('  frozen_dihedral = ("'+'\n')
    firsttor=True
    for resttors in rotbndlist[' '.join([str(b),str(c)])]:
        rta,rtb,rtc,rtd = resttors
        rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
        if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
           (optmol.GetAtom(rtd).GetAtomicNum() != 1):
            if not firsttor:
                temp.write(', %d %d %d %d\n' % (rta,rtb,rtc,rtd))
            else:
                temp.write('    %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                firsttor=True
            
    # Leave all torsions around other rotatable bonds fixed
    for constangle in consttorlist:
        csa,csb,csc,csd,csangle = constangle
        for resttors in rotbndlist[' '.join([str(csb),str(csc)])]:
            rta,rtb,rtc,rtd = resttors
            rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
            if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
               (optmol.GetAtom(rtd).GetAtomicNum() != 1):
                if not firsttor:
                    temp.write(', %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                else:
                    temp.write('    %d %d %d %d\n' % (rta,rtb,rtc,rtd))
                    firsttor=True
    temp.write('  ")'+'\n')
    temp.write('}'+'\n')

    if toroptpcm==True:
        temp.write('set {'+'\n')
        temp.write(' basis '+toroptbasisset.lower()+'\n')
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
    temp.write('psi4_io.set_default_path("%s")'%(scratchdir)+'\n')   
    #temp.write("psi4_io.set_specific_path('%s', './')"%(comfilename.replace('.com','.chk'))+'\n')
    #temp.write("psi4_io.set_specific_retention('%s', True)"%(comfilename.replace('.com','.chk'))+'\n') 
    temp.write("optimize('%s/%s')" % (toroptmethod.lower(),toroptbasisset.lower())+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    return inputname


def CreatePsi4OPTInputFile(comfilecoords,comfilename,mol):
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
    if optpcm==True:
        temp.write('set {'+'\n')
        temp.write(' basis '+optbasisset.lower()+'\n')
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
    temp.write('psi4_io.set_default_path("%s")'%(scratchdir)+'\n')
    #temp.write("psi4_io.set_specific_path('%s', './')"%(comfilename.replace('.com','.chk'))+'\n')
    #temp.write("psi4_io.set_specific_retention('%s', True)"%(comfilename.replace('.com','.chk'))+'\n')
    temp.write("optimize('%s/%s')" % (optmethod.lower(),optbasisset.lower())+'\n')
    if freq==True:
        temp.write('scf_e,scf_wfn=freq(%s/%s,return_wfn=True)'%(optmethod.lower(),optbasisset.lower())+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    return inputname


def CreatePsi4ESPInputFile(comfilecoords,comfilename,mol,makecube=None):
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
    if sppcm==True or torsppcm==True:
        temp.write('set {'+'\n')
        temp.write(' basis '+espbasisset.lower()+'\n')
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
    temp.write('psi4_io.set_default_path("%s")'%(scratchdir)+'\n')
    #temp.write("psi4_io.set_specific_path('%s', './')"%(comfilename.replace('.com','.chk'))+'\n')
    #temp.write("psi4_io.set_specific_retention('%s', True)"%(comfilename.replace('.com','.chk'))+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write("E, wfn = energy('%s/%s',return_wfn=True)" % (espmethod.lower(),espbasisset.lower())+'\n')
    if makecube==True:
       temp.write('oeprop(wfn,"GRID_ESP")'+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    return inputname

def CreatePsi4DMAInputFile(comfilecoords,comfilename,mol):
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
    if sppcm==True or torsppcm==True:
        temp.write('set {'+'\n')
        temp.write(' basis '+dmabasisset.lower()+'\n')
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
    temp.write('psi4_io.set_default_path("%s")'%(scratchdir)+'\n')
    #temp.write("psi4_io.set_specific_path('%s', './')"%(comfilename.replace('.com','.chk'))+'\n')
    #temp.write("psi4_io.set_specific_retention('%s', True)"%(comfilename.replace('.com','.chk'))+'\n')    
    temp.write('set freeze_core True'+'\n')
    temp.write("E, wfn = energy('%s/%s',return_wfn=True)" % (dmamethod.lower(),dmabasisset.lower())+'\n')
    temp.write('gdma(wfn)'+'\n')
    temp.write('fchk(wfn, "%s.fchk")'%(comfilename.replace('.com',''))+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    return inputname

def GrabFinalPsi4Energy(logname):
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
            
        

def GrabFinalPsi4XYZStructure(logname,filename):
    finalmarker=False
    temp=open(logname,'r')
    results=temp.readlines()
    temp.close()
    temp=open(filename,'w')
    temp.write(str(mol.NumAtoms())+'\n')
    temp.write('\n')
    for line in results:
        if 'Final optimized geometry and' in line:
            finalmarker=True
        if finalmarker==True:
            linesplit=line.split()
            if len(linesplit)==4 and NumberInLine(line)==True:
                temp.write(line.lstrip())
    temp.close()

def NumberInLine(line):
    numinline=False
    linesplit=line.split()
    
    for e in linesplit:
        try:
            float(e)
            numinline=True
        except:
            continue
            
    return numinline
        

def tor_opt_sp(molecprefix,a,b,c,d,optmol,consttorlist,phaseangle,prevstrctfname,torsionrestraint):
    """
    Intent: Restrain the torsion to the dihedral angle given (using tinker Minimize tool). 
    Use Gaussian SP calculation to find the new energy. If wanted, Gaussian optimization is done
    on the molecule before the energy is found, but after the torsion is restrained.
    Input:
        molecprefix: molecule name
        a: atom 1 in the torsion of interest
        b: atom 2 in the torsion of interest
        c: atom 3 in the torsion of interest
        d: atom 4 in the torsion of interest
        optmol: optimized OBMol object
        consttorlist: all other rotatable bonds
        phaseangle: what torsion value to restrain it by is given by torang+phaseangle
        prevstrctfname: file containing the current coordinates of the molecule
                        i.e. the coordinates of the molecule fixed at the previous torsion value
                        This is done so that the restraining is done only 30 degrees at a time
    Output:
        prevstrctfname: file name containing the latest coordinates (post torsion restraint)
        many *.log, *.com, and *.chk files are generated for and by Gaussian
    Referenced By: gen_torsion
    Description:
    1. find 'torang', current torsion value
    2. If the do_tor_qm_opt flag was set as true:
        a. Use tinker minimize to optimize the molecule with torsion restrained as 
           'torang+phaseangle'
        b. Coordinates are in *.xyz_2
        c. Generate the *-opt*.com file using the coordinates in *.xyz_2
        d. Append various restraints to *-opt.com, i.e. fix all the rotatable bonds at their
           current values; also fix certain parameter values for certain functional groups
        e. Run Gaussian and optimize the molecule
        f. Run Gaussian SP on the newly optimized molecule to find the energy
        g. Set 'prevstrctfname' to the file containing the latest coordinates (*-opt*.log)
    3. If do_tor_qm_opt flag is false (default):
        a. Use tinker minimize to optimize the molecule with torsion restrained as 'torang'
        b. Coordinates are in *.xyz_2
        c. Generate the *sp*.com file using coordinates *.xyz_2
        d. Run Gaussian SP to find energy
        e. Set 'prevstrctfname' to the file containing the latest coordinates (*sp*.log)
    """
    torang = optmol.GetTorsion(a,b,c,d)
    toroptcomfname = ""
    strctfname = ""
    gaustorerror=False # if there was an error before, reset this to False for next datapoint
    # come here if you would like to run Gaussian Optimization post restraint but before SP
    torspcomfname = '%s-m06lsp-%d-%d-%d-%d-%03d.com' % (molecprefix,a,b,c,d,round((torang+phaseangle)%360))
    torsplogfname = os.path.splitext(torspcomfname)[0] + '.log'
    if do_tor_qm_opt:
        # make the opt com file
        toroptcomfname = '%s-opt-%d-%d-%d-%d-%03d.com' % (molecprefix,a,b,c,d,round((torang+phaseangle)%360))
        strctfname = os.path.splitext(toroptcomfname)[0] + '.log'
    if do_tor_qm_opt and not is_qm_normal_termination(strctfname):
        # make the opt chk
        toroptchkfname = os.path.splitext(toroptcomfname)[0] + '.chk'
        if os.path.isfile(toroptchkfname):
            os.remove(toroptchkfname)
        # load prevstruct
        if use_psi4:            
            prevstrctfname=prevstrctfname.replace('.log','_opt.xyz')
        prevstruct = load_structfile(prevstrctfname)
        rebuild_bonds(prevstruct,optmol)       
        # create xyz and key and write restraint then minimize, getting .xyz_2
        torxyzfname = '%s-opt-%d-%d-%d-%d' % (molecprefix,a,b,c,d)
        torxyzfname+='-%03d-t.xyz' % (round((torang+phaseangle)%360))
        tmpkeyfname = 'tmp-%d-%d-%d-%d' % (a,b,c,d)
        tmpkeyfname+='-%03d-t.key' % (round((torang+phaseangle)%360))
        save_structfile(prevstruct,torxyzfname)
        shutil.copy('../'+key4fname, tmpkeyfname)
        tmpkeyfh = open(tmpkeyfname,'a')
        tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (a,b,c,d,torsionrestraint,round((torang+phaseangle)%360),round((torang+phaseangle)%360)))
        for key in rotbndlist:
            torlist=rotbndlist[key]
            for res in torlist:
                resa,resb,resc,resd = res[0:4]
                if (a,b,c,d) != (resa,resb,resc,resd) and (a,b,c,d) != (resd,resc,resb,resa):
                    if (b==resb and c==resc) or (b==resc and c==resb):
                        secondang = optmol.GetTorsion(resa,resb,resc,resd)
                        tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (resa,resb,resc,resd,torsionrestraint,round((secondang+phaseangle)%360),round((secondang+phaseangle)%360)))
                    else:
                        tmpkeyfh.write('restrain-torsion %d %d %d %d %f\n' % (resa,resb,resc,resd,torsionrestraint))
        tmpkeyfh.close()
        mincmdstr = minimizeexe+' -k '+tmpkeyfname+' '+torxyzfname+' 0.01'
        call_subsystem(mincmdstr)

        # generate the com file using *.xyz_2 which has the restraint
        gen_torcomfile(toroptcomfname,numproc,maxmem,prevstruct,torxyzfname+'_2')
        cmdstr='cp '+torxyzfname+'_2 '+torxyzfname.replace('.xyz','prefitpreqmopt.xyz')
        call_subsystem(cmdstr)
        cmdstr='cp '+tmpkeyfname+' '+tmpkeyfname.replace('.key','prefitpreqmopt.key')
        call_subsystem(cmdstr)
        # remove unnecessary files
        rmcmdstr = 'rm '+tmpkeyfname+'; rm '+torxyzfname+'*'
        call_subsystem(rmcmdstr)

        # Append restraints to *opt*.com file
        tmpfh = open(toroptcomfname, "a")
        
        # Fix all torsions around the rotatable bond b-c 
        for resttors in rotbndlist[' '.join([str(b),str(c)])]:
            rta,rtb,rtc,rtd = resttors
            rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
            if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
               (optmol.GetAtom(rtd).GetAtomicNum() != 1):
                tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))

        # Leave all torsions around other rotatable bonds fixed
        for constangle in consttorlist:
            csa,csb,csc,csd,csangle = constangle
            for resttors in rotbndlist[' '.join([str(csb),str(csc)])]:
                rta,rtb,rtc,rtd = resttors
                rtang = optmol.GetTorsion(rta,rtb,rtc,rtd)
                if (optmol.GetAtom(rta).GetAtomicNum() != 1) and \
                   (optmol.GetAtom(rtd).GetAtomicNum() != 1):
                    tmpfh.write('%d %d %d %d F\n' % (rta,rtb,rtc,rtd))
        tmpfh.close()

        # Append restraint values for specific functional groups
        restraintlist = gen_function_specific_restraints(prevstruct)
        append_restraint(restraintlist,toroptcomfname)

        tmpfh = open(toroptcomfname, "a")
        tmpfh.write("\n")
        tmpfh.close()
        append_basisset(toroptcomfname,prevstruct.GetSpacedFormula(),optbasisset)
        if use_psi4==True and use_psi4SPonly==False:
            logtoroptfname=toroptcomfname.replace('.com','.log')
            now = time.strftime("%c",time.localtime())
            logfh.write(now + " Calling: " + "Psi4 Optimization" + "\n")
            finished=CheckNormalTermPsi4Opt(logtoroptfname)
            
            inputname=CreatePsi4TorOPTInputFile(toroptcomfname,toroptcomfname,mol,b,c,consttorlist,optmol)
            if finished==False:
                cmdstr='psi4 '+inputname+' '+logtoroptfname
                call_subsystem(cmdstr,iscritical=True)
            GrabFinalPsi4XYZStructure(logtoroptfname,strctfname.replace('.log','_opt.xyz'))
            


        else:
            cmdstr = 'GAUSS_SCRDIR='+scrtmpdir+' '+gausexe+' '+toroptcomfname
            call_subsystem(cmdstr,iscritical=True)

        

    if do_tor_qm_opt:
        # prevstrct becomes the opt log found above
        
        prevstrctfname = strctfname

    # Start here if do_tor_qm_opt is False
    # generate initial sp files
    torspcomfname = '%s-m06lsp-%d-%d-%d-%d-%03d.com' % (molecprefix,a,b,c,d,round((torang+phaseangle)%360))
    torsplogfname = os.path.splitext(torspcomfname)[0] + '.log'

    if not is_qm_normal_termination(torsplogfname):
        torspchkfname = os.path.splitext(torspcomfname)[0] + '.chk'
        if os.path.isfile(torspchkfname):
            os.remove(torspchkfname)
        # load previous *.log file
        if use_psi4==False:
            prevstruct = load_structfile(prevstrctfname)
        else:
            prevstruct =  load_structfile(strctfname.replace('.log','_opt.xyz'))
            rebuild_bonds(prevstruct,optmol)

        if not do_tor_qm_opt:
             # create xyz and key and write restraint then minimize, getting .xyz_2
            torxyzfname = '%s-opt-%d-%d-%d-%d' % (molecprefix,a,b,c,d)
            torxyzfname+='-%03d-t.xyz' % (round((torang+phaseangle)%360))
            tmpkeyfname = 'tmp-%d-%d-%d-%d' % (a,b,c,d)
            tmpkeyfname+='-%03d-t.key' % (round((torang+phaseangle)%360))
            save_structfile(prevstruct,torxyzfname)
            shutil.copy('../'+key4fname, tmpkeyfname)
            tmpkeyfh = open(tmpkeyfname,'a')
            tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (a,b,c,d,torsionrestraint,round((torang+phaseangle)%360),round((torang+phaseangle)%360)))
            for key in rotbndlist:
                torlist=rotbndlist[key]
                for res in torlist:
                    resa,resb,resc,resd = res[0:4]
                    if (a,b,c,d) != (resa,resb,resc,resd) and (a,b,c,d) != (resd,resc,resb,resa):
                        if (b==resb and c==resc) or (b==resc and c==resb):
                            secondang = optmol.GetTorsion(resa,resb,resc,resd)
                            tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (resa,resb,resc,resd,torsionrestraint,round((secondang+phaseangle)%360),round((secondang+phaseangle)%360)))
                        else:
                            tmpkeyfh.write('restrain-torsion %d %d %d %d %f\n' % (resa,resb,resc,resd,torsionrestraint))
            tmpkeyfh.close()

            # minimize the structure to the restraint
            mincmdstr = minimizeexe+' -k '+tmpkeyfname+' '+torxyzfname+' 0.01'
            call_subsystem(mincmdstr)

            # generate the *.com file using the minimized *.xyz, *.xyz_2
            gen_torcomfile(torspcomfname,numproc,maxmem,prevstruct,torxyzfname+'_2')
        else:
            gen_torcomfile(torspcomfname,numproc,maxmem,prevstruct,"non")

        # append the proper basis set to the *.com file
        append_basisset(torspcomfname,prevstruct.GetSpacedFormula(),torspbasisset)

        # run Gaussian SP on *.com file
        if use_psi4==False:
            cmdstr = 'GAUSS_SCRDIR='+scrtmpdir+' '+gausexe+' '+torspcomfname
            call_subsystem(cmdstr,iscritical=True)
        else:
            inputname=CreatePsi4ESPInputFile(torspcomfname,torspcomfname,mol)  
            if CheckNormalTermPsi4ESP(torsplogfname)==False:
                now = time.strftime("%c",time.localtime())
                logfh.write(now + " Calling: " + "Psi4 Gradient for ESP" + "\n")
                cmdstr='psi4 '+inputname+' '+torsplogfname
                call_subsystem(cmdstr,iscritical=True)
        if use_psi4 or use_psi4SPonly:
            cmdstr='cp '+strctfname.replace('.log','_opt.xyz')+' '+torsplogfname.replace('.log','_opt.xyz')
            call_subsystem(cmdstr)
            #cmdstr='cp '+'optimized.xyz'+' '+strctfname.replace('.log','_opt.xyz')
            #call_subsystem(cmdstr)
        # prevstrctfname is set to the new log file created (if do_tor_qm_opt is false)
        if not do_tor_qm_opt:
            prevstrctfname = torsplogfname
    return prevstrctfname

def gen_torsion(mol,torsionrestraint):
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
    if use_psi4:
        cmdstr='cp '+'optimized.xyz'+' '+'qm-torsion/.'
        call_subsystem(cmdstr)
    os.chdir('qm-torsion')
    files=os.listdir(os.getcwd())

    global gentorsion
    gentorsion=True
    global foundmin
    foundmin=False

    for tor in torlist:
        
        a,b,c,d = tor[0:4]
        torang = mol.GetTorsion(a,b,c,d)
        key=str(b)+' '+str(c)
        anginc=rotbndtoanginc[key]
        #print('key',key,'anginc ',anginc)
        phaselist=range(0,360,anginc)
        #print('normal phaselist',range(0,360,30),'anglist',phaselist,'len(phaselist)/2 ',len(phaselist)/2)
        clock=phaselist[1:int(len(phaselist)/2)]
        counterclock=[i-360 for i in phaselist[int(len(phaselist)/2) :][::-1]]
        #print('normal clock ',range(30,180,30),'generalized clock',clock)
        #print('normal counter clock ',range(-30,-210,-30),'generalized counterclock',counterclock)
        # create a list of all the other rotatable bonds besides this one
        consttorlist = list(torlist)
        consttorlist.remove(tor)
        #pause=input('pausing...')
        prevstrctfname = '%s-opt-%d-%d-%d-%d-%03d.log' % (molecprefix,a,b,c,d,round(torang % 360))

        minstrctfname = prevstrctfname
        
        # copy *-opt.log found early by Gaussian to 'prevstrctfname'
        cmd = 'cp ../%s %s' % (logoptfname,prevstrctfname)
        call_subsystem(cmd)
        time.sleep(5)

        prevstrctfname = minstrctfname

        # run Gaussian SP on logoptme
        prevstrctfname = tor_opt_sp(molecprefix,a,b,c,d,mol,consttorlist,0,prevstrctfname,torsionrestraint)
        for phaseangle in clock:
            prevstrctfname = tor_opt_sp(molecprefix,a,b,c,d,mol,consttorlist,phaseangle,prevstrctfname,torsionrestraint)

        # Rotate torsion counterclockwise, running Gaussian SP at each rotation
        prevstrctfname = minstrctfname
        for phaseangle in counterclock:
            prevstrctfname = tor_opt_sp(molecprefix,a,b,c,d,mol,consttorlist,phaseangle,prevstrctfname,torsionrestraint)


    os.chdir('..')



def opbset (smarts, opbval, opbhash, mol):
    """
    Intent: Set out-of-plane bend (opbend) parameters using Smarts Patterns
    Referenced By: gen_valinfile
    """
    opbsetstr = ""
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    sp.Match(mol)
    for ia in sp.GetMapList():
        iteratomatom = openbabel.OBAtomAtomIter(mol.GetAtom(ia[1]))
        for ib in iteratomatom:
            opbkey = '%d %d 0 0' % (ib.GetIdx(), ia[1])
            if ib.GetIdx() == ia[0]:
                if ((opbkey not in opbhash) or (opbhash[opbkey][1] == False)):
                    opbhash[opbkey] = [opbval, True]
            else:
                if opbkey not in opbhash:
                    opbhash[opbkey] = [defopbendval, False]

def gen_valinfile (mol):
    """
    Intent: Find aromatic carbon and bonded hydrogens to correct polarizability
    Find out-of-plane bend values using a look up table
    Output a list of rotatable bonds (found in get_torlist) for valence.py
    Input:
        mol: OBMol molecule object
    Output:
        opbhash.items(): list of opbend values
        rotbndprmlist: rotatable bonds for valence.py. 
                       This lets valence.py know to set certain torsions as 0.
    Referenced By: main
    Description: -
    """
    #Find aromatic carbon and bonded hydrogens to correct polarizability
    f = open (valinfile, 'w')
    iteratom = openbabel.OBMolAtomIter(mol)
    for a in iteratom:
        if (a.GetAtomicNum() == 6 and a.IsAromatic()):
            f.write(str(a.GetIdx()) + " " + str(1) + "\n")
        elif (a.GetAtomicNum() == 1):
            iteratomatom = openbabel.OBAtomAtomIter(a)
            for b in iteratomatom:
                if (b.GetAtomicNum() == 6 and b.IsAromatic()):
                    f.write(str(a.GetIdx()) + " " + str(1) + "\n")
    f.write("\n")

    #Search for structures that require opbend parameters
    # OP-Bend parameters is between first and second atom in search string
    opbhash = {}      #Create empty dictionary
    # 51 54 0 0 amide C(-N)(=O)
    opbset('O=CN', 1.7, opbhash, mol)
    opbset('C(-N)(=O)', 1.0, opbhash, mol)
    opbset('[OX1]=N[Oh1]', 0.2, opbhash, mol)
    opbset('[Oh1]N=[OX1]', 0.2, opbhash, mol)

    #opbset('C-C(O)', 1.7, 1, opbhash, mol)

    # 37 51 0 0 acetamide [CH3](-C=O)
    # 37 60 0 0 acetaldehyde [CH3](-C=O)
    opbset('[CH3](-C=O)', 0.590, opbhash, mol)
    # 37 54 0 0 methylformamide [CH3](-NC=O)
    opbset('[CH3](-NC=O)', 0.180, opbhash, mol)

    # 52 51 0 0 formamide HC=O
    opbset('[#1]C=O', 1.950, opbhash, mol)

    # 53 51 0 0 amide O=C
    # 53 60 0 0 carboxylic acid/aldehyde O=C
    opbset('O=C', 0.650, opbhash, mol)

    # 54 51 0 0 amide NC=O
    opbset('NC=O', 1.500, opbhash, mol)
    # 55 54 0 0 amide HN
    opbset('[#1]NC=O', 0.080, opbhash, mol)
    # 58 60 0 0 carboxylic acid [OH](C=O)
    opbset('[OH](C=O)', 1.500, opbhash, mol)
    # 61 51 0 0 aldehyde HC=O for 3-formylindole
    # 61 60 0 0 formic acid/aldehyde HC=O
    opbset('[#1]C=O', 1.950, opbhash, mol)
    # 76 89 0 0 pyridinium cnc
    opbset('cnc', 0.150, opbhash, mol)
    # 77 76 0 0 benzene Hc
    # 79 78 0 0 ethylbenzene/phenol/toluene/p-cresol Hc1aaaaa1
    opbset('[#1]c', 0.210, opbhash, mol)

    # 80 83 0 0 3-formylindole n1caaa1
    opbset('n1caaa1', 0.200, opbhash, mol)
    # 83 83 0 0 3-formylindole c1caaa1
    opbset('c1caaa1', 0.200, opbhash, mol)
    # 84 83 0 0 3-ethylindole [CD4]c1aaaa1
    opbset('[CD4]c1aaaa1', 0.200, opbhash, mol)
    # 78 88 0 0 benzamidine C(c)(N)(=N)
    opbset('cC(N)(=N)', 0.020, opbhash, mol)
    # 83 51 0 0 3-formylindole cC=O
    opbset('c[CH]=O', 0.590, opbhash, mol)
    # 85 51 0 0 3-formylindole O=[CH]c
    opbset('O=[CH]c', 0.650, opbhash, mol)
    # 86 88 0 0 benzamidine N~[C](~N)(c)
    opbset('[NH1,NH2]-,=[C]([NH1,NH2])(c)', 0.020, opbhash, mol)
    # 87 86 0 0 benzamidine HN(~C~N)
    opbset('[#1]N(~C~N)', 0.180, opbhash, mol)
    # 88 78 0 0 benzamidine C(c)(-N)(=N)
    opbset('C(c)(-N)(=N)', 0.100, opbhash, mol)
    # 88 86 0 0 benzamidine C(~N)(~N)c
    opbset('C(~N)(~N)c', 0.050, opbhash, mol)
    # 89 76 0 0 pyridinium nc
    opbset('nc', 0.250, opbhash, mol)
    # 90 89 0 0 Hn
    opbset('[#1]n', 0.150, opbhash, mol)
    # 30 78 0 0 alkane C - aromatic C [CH2;X4]c
    # 37 78 0 0 ethylbenzene [CH2;X4]c
    opbset('[CH2;X4]c', 0.200, opbhash, mol)
    # 35 78 0 0 phenol [OH]c
    opbset('[OH]c', 0.200, opbhash, mol)
    # 51 83 0 0 [CH](c)=O
    opbset('[CH](c)=O', 0.200, opbhash, mol)
    # 74 73 0 0 tricyanomethide C(C(C#N)(C#N))#N
    opbset('C(C(C#N)(C#N))#N', 0.200, opbhash, mol)
    # 76 76 0 0 benzene cc
    # 78 78 0 0 ethylbenzene cc
    opbset('cc', 0.200, opbhash, mol)
    #opbset('[#6D3][*]', 0.200, opbhash, mol)

    #print('hash items ',opbhash.items())

    for (opbkey, opbval) in opbhash.items():
        f.write('%s %.5f %d\n' % (opbkey, opbval[0], opbval[1]))
    f.write("\n")

    rotbndprmlist = []
    for rotbnd in rotbndlist.values():
        for rotbndprm in rotbnd:
            rotlist=list(rotbndprm)
            rotbndprmlist.append(rotlist)
            f.write('%2d %2d %2d %2d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n' % (rotbndprm[0], rotbndprm[1], rotbndprm[2],rotbndprm[3], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    f.write("\n")
    f.close()
    return opbhash.items(),rotbndprmlist

def gen_superposeinfile():
    """
    Intent: Initialize superpose input file (for tinker's superpose) 
    """
    f = open(superposeinfile, 'w')
    f.write('\n\n\n\n\n')

# Create file to specify groups of atoms based on molecular symmetry
def gen_avgmpole_groups_file():
    """
    Intent: Print out *-groups.txt which is a map from symm class to idx
    Also, symm class labels are altered from 1, 2, 3, ... to 401, 402, ...
    (or some other labeling system dependent on 'prmstartidx')
    Input:
    Output: 
        *-groups.txt: map from symm group to atom idx
    Referenced By: main
    Description: -
    """
    symgroups = [None] * max(symmetryclass)
    for i in range(0,len(symgroups)):
        symgroups[i] = []
        symgroups[i].append(prmstartidx + (max(symmetryclass) - i - 1))
    for symclsidx in range(0,len(symmetryclass)):
        symgroups[symmetryclass[symclsidx]-1].append(symclsidx+1)
    f = open(grpfname,"w")
    symgroups.sort()
    for ii in symgroups:
        #print(ii)
        for kk in ii:
            f.write(str(kk) + " " )
        f.write("\n")
    f.close()

def gen_esp_grid(mol):
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
    if not os.path.isfile(espgrdfname):
        gengridcmd = potentialexe + " 1 " + xyzfname+' -k '+keyfname
        call_subsystem(gengridcmd)
    #    shutil.move(xyzoutfile,espgrdfname)
    # Run CUBEGEN
    if use_psi4 or use_psi4SPonly:
        import psi4
        shutil.copy(espgrdfname, 'grid.dat') # same as .grid file (x,y,z) coords
        inputname=CreatePsi4ESPInputFile('optimized.xyz',comespfname,mol,True)
        if CheckNormalTermPsi4ESP(logespfname)==False:
            now = time.strftime("%c",time.localtime())
            logfh.write(now + " Calling: " + "Psi4 Gradient for ESP" + "\n")
            cmdstr='psi4 '+inputname+' '+logespfname
            call_subsystem(cmdstr,iscritical=True)
        temp=open('grid_esp.dat','r')
        results=temp.readlines()
        temp.close()
        Vvals=[]
        for line in results:
            Vvals.append(line.split()[0])
            
        now = time.strftime("%c",time.localtime())
        logfh.write(now + " Calling: " + "Generating CUBE File from PSI4" + "\n")
        with open('grid.dat', 'r') as fp:
            gridpts = fp.readlines()
        if len(gridpts) != len(Vvals):
            raise Exception('Dimension error in potential calculation!')
        # Generate a "cube" file.  I have no idea what the format should be (it's not a
        # regular cube file) so I reverse engineered one by looking at TINKER source.
        with open(qmespfname, 'w') as fp:
            fp.write(" %s/%s potential calculation\n\n" % (espmethod,espbasisset))
            fp.write("%5d\n%5d\n\n\n" % (0,len(Vvals)))
            for xyz,v in zip(gridpts, Vvals):
                fp.write("%s %s\n" % (xyz.rstrip(), v))
        
    elif not os.path.isfile(qmespfname):
        fckfname = fckespfname
        if not espfit:
            fckfname = fckdmafname

        if not os.path.isfile(fckfname):
            fckfname = os.path.splitext(fckfname)[0]

        assert os.path.isfile(fckfname), "Error: " + fckfname + " does not exist."
        gencubecmd = cubegenexe + " 0 potential=MP2 " + fckfname + " " + \
                     qmespfname + " -5 h < " + espgrdfname
        call_subsystem(gencubecmd,iscritical=True)

    # Run potential
    if not os.path.isfile(qmesp2fname):
        genqmpotcmd = potentialexe + " 2 " + qmespfname
        call_subsystem(genqmpotcmd,iscritical=True)

def insert_torprmdict_angle(angle, angledict):
    """
    Intent: Increase the count of this angle by one
    """
    anglekey = int(angle)
    if anglekey in angledict:
        angledict[anglekey] += 1
    else:
        angledict[anglekey] = 1

def tor_func_term (parms, x, nfold, C, angle, offset):
    """
    Intent: Returns energy vs dihedral angle profile of the torsion
    Input: 
        parms: torsion parameter estimate
        x: angle list in radians 
        nfold: fold #
        C: number of times this torsion exists
        angle: current dihedral angle
        offset: an offset for the current fold if it exists
    Output:
        energy profile
    Referenced By: fitfunc
    Description: -
    """
    return C*parms/2.0*(1+numpy.cos(nfold*(x+angle)+offset))


def fitfunc (parms, x, torprmdict, debug = False):
    """
    Intent: Gives energy due to torsion around one rotatable bond 
    (torsion energy vs. dihedral angle) given a set of parameters
    This function is used to make the callable function 'errfunc' used by leastsq
    Input:
        parms: current parameter estimate
        x: angle list in radians 
        torprmdict: contains information about the torsions (like nfolds)
    Output:
        tor_energy: energy due to torsion at various dihedral angles (found using 'parms')
                    (for one rotatable bond)
    Referenced By: fit_rot_bond_tors 
    Description: Loops over torsions about this rotatable bond, then further loops about nfolds,
    making multiple calls to 'tor_func_term', summing its results in 'tor_energy'
    """
    tor_energy = [ 0.0 ] * len(x)
    offset =[ 0.0 ]* len(x)
    #print('Fitting Torsions Around Rotatable Bond')
    if debug:
        pass
        #dispvar('PARMS',parms)
        #dispvar('FFTPD',torprmdict)
    for torprmkey in torprmdict.keys():
        torprm=torprmdict[torprmkey]
        # for each torsion about the same rotatable bond (clskey)
        for nfold in torprm['prmdict']:
            # for each nfold for this torsion (# of parameters)
            for clsangle, clscnt in torprm['phasedict'].items():
                # current dihedral angles and how many torsions are this angle 

                # current parameter for this 'fold'
                prm = torprm['prmdict'][nfold]
               
                # not called by 'eval'
                if parms is not 'eval':
                    # get prm from parms array
                    prm = parms[torprm['prmdict'][nfold]]
                #print('Torsion Parameter Key ',torprmkey,'Torsion Parameter',prm,'Fold Number',nfold,'Number of Times Torsion Occurs Around Bond',clscnt,'Phase Angle To Main Torsion',clsangle,'Fold Offset',foldoffsetlist[nfold-1],'Dihedral Angle List',x)
                # TBC. Why is this being summed for each torsion
                tor_energy += tor_func_term (
                    prm, x, nfold, clscnt, rads(clsangle),
                    rads(foldoffsetlist[nfold-1]))

                if debug:
                    pass
                    #dispvar('TERM', prm, nfold,clsangle,clscnt,tor_energy)

        if parms is 'eval' and 'offset' in torprm:
            offset = torprm['offset']


    if parms is not 'eval':
        offset = parms[-1]
    tor_energy += offset
    return tor_energy

def compute_qm_tor_energy(a,b,c,d,startangle,phase_list = None):
    """
    Intent: Store the QM Energies (vs. Dihedral Angle) found in 'gen_torsion' in a list
    Input:
        a: atom 1 in the torsion of interest
        b: atom 2 in the torsion of interest
        c: atom 3 in the torsion of interest
        d: atom 4 in the torsion of interest
        startangle: current or initial dihedral angle
        phase_list: list of phase offsets. default is 0-360 in intervals of 30
    Output:
        list(rows[1]): QM energies
        list(rows[0]): Dihedral angles
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: Read in the *-m06lsp-*.log files created in 'gen_torsion', find the energy
    values, and store them in a list.
    """
    if phase_list is None:
        key=str(b)+' '+str(c)
        anginc=rotbndtoanginc[key]
        #print('key',key,'anginc ',anginc)
        phase_list=range(0,360,anginc) 
        #phase_list = range(0,360,30)

    energy_list = []
    angle_list = []
    energy_dict = {}
    for phaseangle in phase_list:
        angle = (startangle + phaseangle) % 360
        minstrctfname = '%s-m06lsp-%d-%d-%d-%d-%03d.log' % (molecprefix,a,b,c,d,round(angle))
        if not os.path.exists(minstrctfname): # if optimization failed then SP file will not exist
            tor_energy=None
        else:
            tmpstrct = load_structfile(minstrctfname)
            tmpfh = open(minstrctfname, 'r')
            tor_energy = None
            if use_psi4 or use_psi4SPonly:
                mengi=GrabFinalPsi4Energy(minstrctfname)
                tor_energy = float(mengi) * Hartree2kcal_mol
            else:
                for line in tmpfh:
#                m = re.search(r'E\(RM06L\) =\s+(\-*\d+\.\d+)',line)
                    m = re.search(r'EUMP2 =\s+(\-*\d+\.\d+D\+\d+)',line)
    # this change below will work for both M06L (has no D+ in results) and MP2
                    if not m is None:
                        mengi = m.group(1).replace('D+', 'E+')
                        tor_energy = float(mengi) * Hartree2kcal_mol
        tmpfh.close()
        energy_list.append(tor_energy)
        angle_list.append(angle)
        energy_dict[angle] = tor_energy
    rows = zip(*[angle_list, energy_list])
    rows=sorted(rows)
    rows0=list([i[0] for i in rows])
    rows1=list([i[1] for i in rows])
    return rows1,rows0

def compute_mm_tor_energy(mol,a,b,c,d,startangle,phase_list = None,keyfile = None):
    """
    Intent: Use tinker analyze to find the Pre-fit MM Energy vs. Dihedral Angle profile
    Input:
        a: atom 1 in the torsion of interest
        b: atom 2 in the torsion of interest
        c: atom 3 in the torsion of interest
        d: atom 4 in the torsion of interest
        startangle: current or initial dihedral angle
        phase_list: list of phase offsets. default is 0-360 in intervals of 30
        keyfile: keyfile for tinker analyze
    Output:
        list(rows[1]): MM energies
        list(rows[0]): Dihedral angles
        list(rows[2]): Energy just due to torsion 
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description:
    1. For each phase offset
        a. Restrain the dihedral angle at (startangle + phaseangle)
        b. Run tinker analyze (for this new restraint)
        c. Read in and store the energy
    """
    if phase_list is None:
        key=str(b)+' '+str(c)
        anginc=rotbndtoanginc[key]
        #print('key',key,'anginc ',anginc)
        phase_list=range(0,360,anginc)
        #phase_list = range(0,360,30)
    energy_list = []
    torse_list = []
    angle_list = []

    for phaseangle in phase_list:
        angle = (startangle + phaseangle) % 360
        minstrctfname = '%s-m06lsp-%d-%d-%d-%d-%03d.log' % (molecprefix,a,b,c,d,round(angle))
        if use_psi4 or use_psi4SPonly:
            minstrctfname=minstrctfname.replace('.log','_opt.xyz')
        tmpstrct = load_structfile(minstrctfname)
        torxyzfname = '%s-%d-%d-%d-%d-%03d.xyz' % (molecprefix,a,b,c,d,round(angle))
        tmpkeyfname = 'tmp-%d-%d-%d-%d-%03d_%d.key' % (a,b,c,d,round(angle),mm_tor_count)
        result = save_structfile(tmpstrct, torxyzfname)
        toralzfname = os.path.splitext(torxyzfname)[0] + '.alz'
        if keyfile:
            shutil.copy(keyfile, tmpkeyfname)
            tmpkeyfh = open(tmpkeyfname,'a')
            # fix the current dihedral at the angle found above (init + phase_offset)
            # fix the other rotatable bonds where they are currently at
             # create xyz and key and write restraint then minimize, getting .xyz_2
        
            tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (a,b,c,d,torsionrestraint,round((startangle+phaseangle)%360),round((startangle+phaseangle)%360)))
            for key in rotbndlist:
                torlist=rotbndlist[key]     
                for res in torlist:
                    resa,resb,resc,resd = res[0:4]
                    if (a,b,c,d) != (resa,resb,resc,resd) and (a,b,c,d) != (resd,resc,resb,resa):
                        if (b==resb and c==resc) or (b==resc and c==resb):
                            secondang = mol.GetTorsion(resa,resb,resc,resd)
                            tmpkeyfh.write('restrain-torsion %d %d %d %d %f %6.2f %6.2f\n' % (resa,resb,resc,resd,torsionrestraint,round((secondang+phaseangle)%360),round((secondang+phaseangle)%360)))
                        else:
                            tmpkeyfh.write('restrain-torsion %d %d %d %d %f\n' % (resa,resb,resc,resd,torsionrestraint))
            tmpkeyfh.close()

        if not keyfile:
            mincmdstr=minimizeexe+' '+torxyzfname+' 0.01' # make sure other .xyz files were removed
            alzcmdstr=analyzeexe+' '+torxyzfname+'_'+str(2)+' ed > %s' % toralzfname
        else:
            mincmdstr=minimizeexe+' -k '+tmpkeyfname+' '+torxyzfname+' 0.01' # make sure other .xyz files were removed
            alzcmdstr=analyzeexe+' -k '+tmpkeyfname+' '+torxyzfname+'_'+str(2)+' ed > %s' % toralzfname
        print(mm_tor_count,torxyzfname)
        if not do_tor_qm_opt:
            if postfit==False:
                newxyzname=torxyzfname.replace('.xyz','')+'prefit.xyz'
                newxyzminname=torxyzfname.replace('.xyz','')+'prefitmin.xyz'
                newalzname=toralzfname.replace('.alz','')+'prefit.alz'
                newkeyname=tmpkeyfname.replace('_'+str(mm_tor_count)+'.key','prefit.key')
            else:
                newxyzname=torxyzfname.replace('.xyz','')+'postfit.xyz'
                newxyzminname=torxyzfname.replace('.xyz','')+'postfitmin.xyz'
                newalzname=toralzfname.replace('.alz','')+'postfit.alz'
                newkeyname=tmpkeyfname.replace('_'+str(mm_tor_count)+'.key','postfit.key')
        else:
            if postfit==False:
                newxyzname=torxyzfname.replace('.xyz','')+'prefitpostqmopt.xyz'
                newxyzminname=torxyzfname.replace('.xyz','')+'prefitminpostqmopt.xyz'
                newalzname=toralzfname.replace('.alz','')+'prefitpostqmopt.alz'
                newkeyname=tmpkeyfname.replace('_'+str(mm_tor_count)+'.key','prefit.key')
            else:
                newxyzname=torxyzfname.replace('.xyz','')+'postfit.xyz'
                newxyzminname=torxyzfname.replace('.xyz','')+'postfitmin.xyz'
                newalzname=toralzfname.replace('.alz','')+'postfit.alz'
                newkeyname=tmpkeyfname.replace('_'+str(mm_tor_count)+'.key','postfit.key')

        call_subsystem(mincmdstr)
        call_subsystem(alzcmdstr)
        #print('mm_tor_count',mm_tor_count)
        #print('torxyzfname',torxyzfname,'newxyzname',newxyzname)
        renamecmdstr1 = 'mv ' + torxyzfname+ ' '+newxyzname
        renamecmdstr2 = 'mv ' + torxyzfname+'_2'+ ' '+newxyzminname
        copycmdstr3 = 'cp ' + toralzfname+ ' '+newalzname
        copycmdstr4 = 'cp ' + tmpkeyfname+ ' '+newkeyname
        call_subsystem(renamecmdstr1)
        call_subsystem(renamecmdstr2)
        call_subsystem(copycmdstr3)
        call_subsystem(copycmdstr4)
        tmpfh = open(toralzfname, 'r')
        tot_energy = None
        tor_energy = None
        for line in tmpfh:
            m = re.search(r'Potential Energy :\s+(\-*\d+\.\d+)',line)
            if not m is None:
                tot_energy = float(m.group(1))
            m = re.search(r'Torsional Angle\s+(\-*\d+\.\d+)',line)
            if not m is None:
                tor_energy = float(m.group(1))
        tmpfh.close()
        energy_list.append(tot_energy)
        torse_list.append(tor_energy)
        angle_list.append(angle)
    if None in energy_list:
        string='Cannot analyze XYZ file for torsion %d-%d-%d-%d'%(a,b,c,d)
        errstr = [string, energy_list,angle_list]

    rows = zip(*[angle_list, energy_list, torse_list])
    rows=sorted(rows)
    rows0=list([i[0] for i in rows])
    rows1=list([i[1] for i in rows])
    rows2=list([i[2] for i in rows])

    return rows1,rows0,rows2

#def find_mme_error(mme_list,current_ang_list,cumul_ang_list):
#    dup_list = []
#    del_list = []
#    for listidx in range(0,len(mme_list)-1):
#        mm_eng = mme_list[listidx]
#        if mm_eng is None:
#            if current_ang_list[listidx] in cumul_ang_list:
#                dup_list.append(listidx)
#            else:
#                del_list.append(listidx)
#    return dup_list,del_list

def find_del_list(mme_list,current_ang_list):
    """
    Intent: Run through 'mme_list' and remove None objects;
            remove the corresponding angle as well
    Input:
        mme_list: list of MM energies (vs. angle)
        current_ang_list: list of angles
    Output: 
        del_ang_list: List of angles to remove since the MM energy was not found for that angle
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms 
    Description: -
    """
    del_ang_list = []
    for listidx in range(0,len(mme_list)-1):
        mm_eng = mme_list[listidx]
        if mm_eng is None:
            del_ang_list.append(current_ang_list[listidx])
    return del_ang_list

def sum_xy_list(x1,y1,x2,y2):
    for xx in x1:
        if xx in x2:
            idx1 = x1.index(xx)
            idx2 = x2.index(xx)
            y2[idx2] = y1[idx1]

def del_tor_from_fit(dellist, torprmdict):
    #dispvar("TPD",torprmdict)
    for delitem in dellist:
        #dispvar("DEL",delitem[0],torprmdict[delitem[0]]['prmdict'],
        #    torprmdict[delitem[0]]['prmdict'][delitem[1]])
        del torprmdict[delitem[0]]['prmdict'][delitem[1]]
        if not torprmdict[delitem[0]]['prmdict']:
            torprmdict[delitem[0]]['count'] = 0
            torprmdict[delitem[0]]['phasedict'] = {}
    #dispvar('TPDint', torprmdict)

    # Renumber parameter indices
    idx = 0
    allidx = []
#   dispvar("TPDa", toraboutbnd)
    for toraboutbnd in torprmdict:
        for fold in torprmdict[toraboutbnd]['prmdict']:
            allidx.append(torprmdict[toraboutbnd]['prmdict'][fold])

    # Only keep unique indices
    allidx=list(set(allidx))
    allidx.sort()
    for toraboutbnd in torprmdict:
        for fold in torprmdict[toraboutbnd]['prmdict']:
            torprmdict[toraboutbnd]['prmdict'][fold] = \
                allidx.index(torprmdict[toraboutbnd]['prmdict'][fold])

    #dispvar("TPDa", torprmdict, idx)
    return len(allidx) + 1

def find_least_connected_torsion(torprmdict):
    """
    Find least connected torsion (i.e. the two outer atoms have the highest summed class number)
    """

    least_connected_tor = None
    highest_clssum = 0
    keylist = torprmdict.keys()

    for chkclskey in keylist:
        a,b,c,d = chkclskey.split()
        cur_clssum = int(a) + int(d)
        if (least_connected_tor is None or cur_clssum > highest_clssum):
            least_connected_tor = chkclskey
            highest_clssum = cur_clssum
    return least_connected_tor

def prune_mme_error(del_ang_list,*arr_list):
    """
    Intent: Delete the given ids (del_ang_list) in every list in *arr_list
    Input:
        del_ang_list: list of angles to delete
        *arr_list: list of lists; the corresponding elements of each list will be removed
    Output: -
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: -
    """
    arr_list_idx = 0
    x_list = arr_list[0]
    for del_ang in del_ang_list:
        if del_ang in x_list:
            del_idx = x_list.index(del_ang)
            for a_list in arr_list:
                del a_list[del_idx]
    return arr_list

def prune_qme_error(del_ang_list,*arr_list):
    """
    Intent: Delete the given ids (del_ang_list) in every list in *arr_list
    Input:
        del_ang_list: list of angles to delete
        *arr_list: list of lists; the corresponding elements of each list will be removed
    Output: -
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: -
    """
    arr_list_idx = 0
    x_list = arr_list[0]
    for del_ang in del_ang_list:
        if del_ang in x_list:
            del_idx = x_list.index(del_ang)
            for a_list in arr_list:
                del a_list[del_idx]
    return arr_list

def get_qmmm_rot_bond_energy(mol,tmpkey1basename):
    """
    Intent: Form dicts for each torsion in torlist, mapping the torsion class key ('clskey') to 
    an energy profile (dihedral angle vs. energy). 'cls_mm_engy_dict' maps 'clskey' to pre-fit MM 
    calculated energy profiles, 'cls_qm_engy_dict' maps 'clskey' to QM calculated energy profiles
    Input:
        mol: OBMol Structure
        anglist: phase list. default: 0 - 360 in increments of 30
        tmpkey1basename: key file name for tinker
    Output:
        cls_mm_engy_dict: given a class key, this will provide a list of mm energies (vs. angles)
        cls_qm_engy_dict: given a class key, this will provide a list of qm energies (vs. angles)
        cls_angle_dict: given a class key, this provides the angles that the energies
        above are based on
    Referenced By: process_rot_bond_tors
    Description:
    1. For each rotatable bond
        a. Find QM Energy profile
        b. Find MM Energy profile
        c. delete parts of the list where MM energy was not able to be found
        d. Add the profiles found above to 'cls_qm_engy_dict' and 'cls_mm_engy_dict' for this
        clskey
    2. If a clskey has already been looked at more than once, find the average energies of all
       the torsions having the same clskey
    """
    cls_mm_engy_dict = {}
    cls_qm_engy_dict = {}
    cls_angle_dict = {}
    clscount_dict = {}
    for tor in torlist:
        a,b,c,d = tor[0:4]
        torang = mol.GetTorsion(a,b,c,d)
        key=str(b)+' '+str(c)
        anginc=rotbndtoanginc[key]
        #print('key',key,'anginc ',anginc)
        anglist=range(0,360,anginc)
        # create clskey
        clskey = get_class_key(a,b,c,d)
        # initialize dict-values (in this case lists)
        if clskey not in clscount_dict:
            clscount_dict[clskey] = 0
            cls_mm_engy_dict[clskey] = [0]*len(anglist)
            cls_qm_engy_dict[clskey] = [0]*len(anglist)
            cls_angle_dict[clskey] = [0]*len(anglist)

        clscount_dict[clskey] += 1
        mme_list = []  # MM Energy before fitting to QM torsion energy
        qme_list = []  # QM torsion energy
        initangle = mol.GetTorsion(a,b,c,d)

        # find qm, then mm energies of the various torsion values found for 'tor'
        qme_list,qang_list = compute_qm_tor_energy(a,b,c,d,initangle,anglist)
        mme_list,mang_list,tor_e_list = compute_mm_tor_energy(mol,
            a,b,c,d,initangle,anglist,tmpkey1basename)

        # delete members of the list where the energy was not able to be found 
        del_ang_list = find_del_list(mme_list,mang_list)
        (cls_angle_dict[clskey],cls_mm_engy_dict[clskey])=prune_mme_error(del_ang_list,cls_angle_dict[clskey],cls_mm_engy_dict[clskey])
        (mang_list,mme_list,qme_list,qang_list,tor_e_list)=prune_mme_error(del_ang_list,mang_list,mme_list,qme_list,qang_list,tor_e_list)
        del_ang_list = find_del_list(qme_list,qang_list)
        (cls_angle_dict[clskey],cls_mm_engy_dict[clskey])=prune_qme_error(del_ang_list,cls_angle_dict[clskey],cls_mm_engy_dict[clskey])
        (mang_list,mme_list,qme_list,qang_list,tor_e_list)=prune_qme_error(del_ang_list,mang_list,mme_list,qme_list,qang_list,tor_e_list)

        # add the corresponding energy lists found above to cls_*_dict[clskey]
        #print('cls_qm_engy_dict[clskey]',cls_qm_engy_dict[clskey])
        #print('qme_list',qme_list)
        cls_qm_engy_dict[clskey] = [ runsum+eng for runsum,
            eng in zip (cls_qm_engy_dict[clskey], qme_list)]
        cls_mm_engy_dict[clskey] = [ runsum+eng for runsum,
            eng in zip (cls_mm_engy_dict[clskey], mme_list)]
        cls_angle_dict[clskey] = mang_list

    # if multiple class keys, take the average
    for clskey in clscount_dict:
        cnt = clscount_dict[clskey]
        cls_mm_engy_dict[clskey] = [eng/cnt for eng in cls_mm_engy_dict[clskey]]
        cls_qm_engy_dict[clskey] = [eng/cnt for eng in cls_qm_engy_dict[clskey]]

    return cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict

def check_cooperative_tor_terms(clskey, nfold, angle, xvals, torprmdict):
    """ 
    Intent: Check if a cosine term for a torsion is cooperative with the corresponding cosine term of another torsion
    i.e. do they have equal or opposite phase shifts
    Input:
        clskey: class key for the torsion
        nfold: the 'fold' of the cosine term for the torsion
        angle: angle that the torsion is at
        xvals: angle list
        tormprmdict: dictionary containing information about the torsions about a rotatable bond
    Output:
        ck: class key of the torsion that has a cooperative term
    Referenced By: insert_torprmdict
    Description:
    """
    # for each class key
    for (ck, tp) in torprmdict.items():
        # for each angle for the class key
        for phase in tp['phasedict']:
            # if this torsion has the same number of folds
            if nfold in tp['prmdict'] and \
               not isinstance(tp['prmdict'][nfold], tuple):
                #tor_func_term is a cosine function
                outer_tor_e = tor_func_term(1.0,xvals,nfold,1.0,rads(phase),
                                  rads(foldoffsetlist[nfold-1]))
                inner_tor_e = tor_func_term(1.0,xvals,nfold,1.0,rads(angle),
                                  rads(foldoffsetlist[nfold-1]))
                # sum the profiles, then subtract the min from the list
                sum_e = outer_tor_e + inner_tor_e
                sum_e -= min(sum_e)
                # diff the profiles, subtract the min from the list
                diff_e = outer_tor_e - inner_tor_e
                diff_e -= min(diff_e)
                # if when summed, or diffed, the profile becomes close to constant
                # i.e. the difference from the max and the min is < 1e-10
                # then return the class key for this similar torsion
                if max(sum_e) < 1e-10 or max(diff_e) < 1e-10:
                    return ck
    return None

def insert_torphasedict (mol, toraboutbnd, torprmdict, initangle,
    write_prm_dict, keyfilter = None):
    """
    Intent: Adds torsion to be fitted to torprmdict.
    Input:
        mol: An openbabel molecule structure
        toraboutbnd: A list containing the quadruplet of atoms to add
        write_prm_dict: A dict of torsion parameters that don't need fitting
        tormprmdict: dictionary containing information about the torsions about a rotatable bond
        initangle: The torsion angle of the torsion in optimized geometry
        keyfilter: A list of class numbers to allow adding to torprmdict
    Output:
        Modifies appends toraboutbnd to torprmdict or appends to
        write_prm_dict (with 0s for parameters).
    Referenced By: fit_rot_bond_tors
    Description: Adds the torsion 'toraboutbnd' to 'torprmdict'.
    If a torsion with this class key already exists in torprmdict, increase the count
    Adds the current dihedral angle of the torsion to the phasedict of torprmdict
    """
    # quadruplet
    a2,b2,c2,d2 = toraboutbnd
    # create atom structures
    obaa = mol.GetAtom(a2)
    obab = mol.GetAtom(b2)
    obac = mol.GetAtom(c2)
    obad = mol.GetAtom(d2)
    # create a key
    # because it is using symmetry classes instead of atom id's, tpdkey can repeat
    tpdkey = get_class_key(a2, b2, c2, d2)

    # if the key passes the keyfilter or if the keyfilter does exist
    # and, the end atoms are not hydrogens
    if (keyfilter is None or keyfilter == tpdkey): # MODIFY THIS LINE TO ALLOW HYDROGEN TORSION
        # current torsion value (normalized by initangle)
        torabangle = round(
            mol.GetTorsion(obaa,obab,obac,obad) -
            initangle) % 360
        #print('Dihedral Angle of Main Torsion ',initangle,'Other Torsion About Bond',a2,b2,c2,d2,'Phase Angle',torabangle)
        if tpdkey in torprmdict:
            # increase the count for this tpdkey
            torprmdict[tpdkey]['count'] += 1
            insert_torprmdict_angle(
                torabangle, torprmdict[tpdkey]['phasedict'])
        else:
            # set up dict
            torprmdict[tpdkey] = {}
            torprmdict[tpdkey]['count'] = 1
            torprmdict[tpdkey]['phasedict'] = {}
            torprmdict[tpdkey]['prmdict'] = {}
            if len(torprmdict) == 1:
                torprmdict[tpdkey]['offset'] = 1.
            # alter count for this current angle (torabangle)
            insert_torprmdict_angle(
                torabangle,torprmdict[tpdkey]['phasedict'])

    # else, force constants are set to 0
    else:
        write_prm_dict[tpdkey] = {1:0., 2:0., 3:0.}


def insert_torprmdict (mol, torprmdict):
    """
    Intent: Initialize the prmdicts in torprmdict 
    Give each torsion intially 3 folds
    Input:
        mol: OBMol object
        tormprmdict: dictionary containing information about the torsions about a rotatable bond
    Output:
        torprmdict is modified
        prmidx: the number of parameters to be fitted for
    Referenced By: fit_rot_bond_tors
    Description: -
    """
    tmpx = numpy.arange ( 0.0, 360.0, 10)
    prmidx = 0
    # for each cls key
    for (chkclskey, torprm) in torprmdict.items():
        # for each parameter in the energy equation for this torsion
        # nfoldlist = [1,2,3] 
        for nfold in nfoldlist:
            # init array
            test_tor_energy = numpy.zeros(len(tmpx))
            # for each dihedral angle about this rotatable bond and the number of time it occurs 
            # sum up test_tor_energy
            for (angle, scale) in torprm['phasedict'].items():
                # create a test energy list (0 - 360, 10) starting from this angle
                test_tor_energy += tor_func_term(
                    1.0, tmpx, nfold, scale, rads(angle),
                    rads(foldoffsetlist[nfold-1]))

                # check if another torsion has a similar profile 
                basetorkeys = check_cooperative_tor_terms(
                    chkclskey, nfold, angle, tmpx, torprmdict)

            # normalize
            test_tor_energy -= min(test_tor_energy)
            #print('nfold in dict setup ',nfold)
            #print(nfold,test_tor_energy)
            if torprm['phasedict']:
                if basetorkeys is not None:
                    # if the energy profile for this torsion/fold is not so dissimilar from
                    # another torsion/fold, then take its parameter 
                    # in this case, prmidx is not increased by one
                   # print('in here ')
                    torprmdict[chkclskey]['prmdict'][nfold] = torprmdict[basetorkeys]['prmdict'][nfold]
                # will come here first
                
                elif max(test_tor_energy) > 1e-10:
                    #print('actually in here',nfold)
                    #rint(test_tor_energy,max(test_tor_energy))
                    torprm['prmdict'][nfold] = prmidx
                    prmidx += 1
        #print('prmdict',torprmdict[chkclskey]['prmdict'])
        if not torprmdict[chkclskey]['prmdict']:
            # make empty
            torprmdict[chkclskey]['count'] = 0
            torprmdict[chkclskey]['phasedict'] = {}

    prmidx += 1
    return prmidx

def is_torprmdict_all_empty (torprmdict):
    """
    Intent: Determines if all torsion parameters for fitting have been eliminated.
    Input:
        torprmdict: a dict containing list of torsions to be fit
    Ouptput:
        results: True if empty, False if not
    Referenced By: fit_rot_bond_tors
    Description: -
    """
    result = True
    for clskey in torprmdict:
        #dispvar("CHK",torprmdict[clskey])
        if torprmdict[clskey]['prmdict']:
            return False
    return result

def fit_rot_bond_tors(mol,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict):
    """
    Intent: Uses scipy's optimize.leastsq function to find estimates for the torsion 
    parameters based on energy values found at various angles using qm and mm
    Each rotatable bond is fit for one at a time
    Input:
        mol: OBMol structure
        cls_mm_engy_dict: given a class key, this will provide a list of mm energies (vs. angles)
        cls_qm_engy_dict: given a class key, this will provide a list of qm energies (vs. angles)
        cls_angle_dict: given a class key, this provides the angles that the energies
        above are based on
    Output:
        write_prm_dict: map from class key to parameter information. 
                        Used to write out new key file
        fitfunc_dict: energy profile (using the new parameters) to be plotted
    Referenced By: process_rot_bond_tors
    Description:
    1. Initialize 'fitfunc_dict' and 'write_prm_dict'
    2. For each tor in torlist (essentially, for each rotatable bond) 
    (the fit is done for each rotatable bond one at a time):
        a. Initialize 'torprmdict'
            i. For each torsion about the current rotatable bond, 'torprmdict' maps the torsion
               class key to a set of dictionaries containing information about that torsion.
               The three dictionaries containing information are: 
               *count: # times a torsion with this class key exists 
               (there can be multiple torsions about the rotatable bond that have the same class
               key, for example two H-CC-O 's can exist)
               *prmdict: parameters for this torsion
                         The energy equation for the torsion is a sum of cosines
                         'nfolds' is the number of cosines in the sum 
                         each 'fold' or cosine has a coefficient
                         prmdict contains information on the number of folds and the parameters
                         for each fold
               *phasedict: the various dihedral angles that this torsion currently exists at,
                and how many torsions are at this angle
        b. Get the atoms involved in the main torsion around this rotatable bond
        c. Get the current dihedral angle and the class key
        d. Fill in the phasedict portion of 'torprmdict' by calling method 'insert_torphasedict'
           on each torsion about the current rotatable bond 
        e. Edit 'torprmdict' by calling 'insert_torprmdict'
           'prmidx' now equals: number of parameters to be fit
        f. 'angle_list', 'mm_energy_list', 'qm_energy_list' are all initialized for the current
            classkey
        g. normalize the qm_energy_list and mm_energy_list by subtracting all values by the min 
        h. 'tor_energy_list" is created by subtracting mme from qme and is written to a file
        i. 'nfolds' list (of lists) is found. Should at least initially be [[1,2,3],[1,2,3],...]
            nfolds are the number of force constants per torsion in question
        j. 'max_amp', max - min of tor_energy_list, is found 
        k. 'pzero' is initialized
        l. Remove parameters while # of parameters > # of data points
           This can happen if two torsions are very similar so their parameters can be combined
           into one set
        m. now run optimize.leastsq. Keep rerunning it until the parameters no longer have to be 
           'sanitized', meaning that none of the parameter estimates found by leastsq
           are greater than max_amp
        n. If all of the parameter estimates ended up being deleted, leastsq is rerun, 
           this time fitting for only the main torsion
        o. fill in 'torprmdict' with parameter estimates found by leastsq
        p. write out a plot of the fit
        q. write out the parameter estimates
    """
    global postfit
    postfit=True
    fitfunc_dict = {}
    write_prm_dict = {}
    # For each rotatable bond 
    for tor in torlist:
        torprmdict = {}
        # get the atoms in the main torsion about this rotatable bond
        a,b,c,d = tor[0:4]
        # current torsion value
        torang = mol.GetTorsion(a,b,c,d)
        #print('Fitting Torsions around '+str(b)+'-'+str(c))
        # class key; ie symmetry classes key
        clskey = get_class_key(a,b,c,d)
        
        # new list, post fitting
        mm_energy_list2 = [] # MM Energy after fitting

        rotbndkey = '%d %d' % (b, c)
        initangle = mol.GetTorsion(a,b,c,d)
        #print(' torprmdict before insert torphase',torprmdict)
        # Identify all torsion parameters involved with current rotatable bond.
        #print('rotbndlist[rotbndkey] ',rotbndlist[rotbndkey])
        for toraboutbnd in rotbndlist[rotbndkey]:
            # toraboutbnd: some torsion about the current rotatable bond (rotbndkey)
            # However, initangle is the current angle for 'tor' not for 'toraboutbnd'
            insert_torphasedict(mol, toraboutbnd, torprmdict, initangle, write_prm_dict)
        #print(' torprmdict after insert torphase',torprmdict)
        #dispvar('TPDa', torprmdict)
        # create default torsion force constants
        prmidx = insert_torprmdict(mol, torprmdict)
        #dispvar('TPDb', torprmdict)
        #print('Number of Parameters Being Fit '+str(prmidx))
        #print(' torprmdict after insert torprmdict',torprmdict)
        # get all the lists for the current clskey
        angle_list = cls_angle_dict[clskey]  # Torsion angle for each corresponding energy
        mm_energy_list = cls_mm_engy_dict[clskey]  # MM Energy before fitting to QM torsion energy
        qm_energy_list = cls_qm_engy_dict[clskey]  # QM torsion energy
        # 'normalize'
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]

        # Parameterize each group of rotatable bond (identified by
        #  atoms restrained during restrained rotation.
        # tor_energy_list is set as qm - mm
        tor_energy_list = [qme - mme for qme,mme in zip(qm_energy_list,mm_energy_list)]
        Tx = numpy.arange ( 0.0, 360.0, 30)
        txtfname = "%s-fit-%d-%d-%d-%d.txt" % (molecprefix, a, b, c, d)
        # create initial fit file, initially it seems to be 2d instead of 3d
        write_arr_to_file(txtfname,[Tx,tor_energy_list])

        #pzero = []
        #pzero = [len(torprm['prmdict']) for torprm in torprmdict.values()]
        # array of arrays. Usually: [1,2,3], [1,2,3]
        nfolds = [torprm['prmdict'].keys() for torprm in torprmdict.values()]
        

        # max amplitude of function
        max_amp = max(tor_energy_list) - min(tor_energy_list)
        pzero = [ max_amp ] * prmidx
        #print('torprmdict before removing torsions to prevent overfitting',torprmdict)
        # Remove parameters while # of parameters > # data points
        while prmidx > len(mm_energy_list):
            dellist= []
            least_conn_tor = find_least_connected_torsion(torprmdict)
            for nfold in torprmdict[least_conn_tor]['prmdict']:
                dellist.append((least_conn_tor,nfold))
            prmidx = del_tor_from_fit(dellist,torprmdict)
        pzero = [ max_amp ] * prmidx
        #print('torprmdict after removing torsions to prevent overfitting',torprmdict)

        # run leastsq until all the parameter estimates are reasonable
        parm_sanitized = False
        while not parm_sanitized:
            parm_sanitized = True
            dellist = []
            keylist = list(torprmdict.keys())
            keylist.reverse()
            
            # creating a new function, errfunc
            # p: parameters
            # x: angle list
            # torprmdict: torsion information 
            # y: tor_energy_list 
            errfunc = lambda p, x, torprmdict, y: fitfunc(p, x, torprmdict) - y

            # optimize.leastsq is run, found in the scipy library
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html
            # Inputs:
            # errfunc : the callable function that we supply, based on fitfunc
            # This is the function for which the squares are attempted to be minimized
            # pzero : starting parameter estimates
            # args=(rads(numpy.array(angle_list)), torprmdict, tor_energy_list) :
            # These are the arguments that errfunc needs aside from 'p' (which will be supplied by
            # leastsq)
            # full_output : set to True, leastsq returns all optional outputs
            # Outputs:
            # p1 : parameter estimates returned by optimize.leastsq
            # covx : matrix that can be used to find covariance of parameter estimates
            # idict : dictionary containing info about the run
            # msg : string giving cause of failure if one exists
            # ier : an int. if ier is 1,2,3 or 4, a solution was found, else no solution was found
            #p1,covx,idict,msg,ier = optimize.leastsq(errfunc, pzero, args=(rads(numpy.array(angle_list)),torprmdict, tor_energy_list), full_output = True)
            array=optimize.least_squares(errfunc, pzero, jac='2-point', bounds=(-max_amp, max_amp), args=(rads(numpy.array(angle_list)),torprmdict, tor_energy_list))
            p1=array['x']
            #print('Fitted Parameters ',array)

            # Remove parameters found by least.sq that aren't reasonable; 
            # remove parameters found that are greater than max_amp
            for chkclskey in keylist:
                for nfold in torprmdict[chkclskey]['prmdict']:
                    #dispvar("P1",nfold,p1[torprmdict[chkclskey]['prmdict'][nfold]])
                    #dispvar("TPD",max_amp,torprmdict[chkclskey]['prmdict'][nfold])
                    # if the param value is greater than max_amp, remove it
                    if abs(p1[torprmdict[chkclskey]['prmdict'][nfold]]) > max_amp:
                        dellist.append((chkclskey,nfold))
                        parm_sanitized = False
                        break

            dellist = list(set(dellist))
            dispvar("List of torsion parameters being removed that are greater than max amplitude",dellist)
            prmidx = del_tor_from_fit(dellist,torprmdict)
            # new parameter array since prm size may have changed due to deletions
            pzero = [ max_amp ] * prmidx

        # Attempts to insert main torsion type if all are removed
        # Rerun leastsq, this time fitting for the force constants of the main torsion
        if is_torprmdict_all_empty(torprmdict):
            toraboutbnd = rotbndlist[rotbndkey][0]
            insert_torphasedict(
                mol, toraboutbnd, torprmdict,
                initangle, write_prm_dict,keyfilter = clskey)

            prmidx = insert_torprmdict(mol, torprmdict)
            pzero = [ max_amp ] * prmidx

            errfunc = lambda p, x, torprmdict, y: fitfunc(p, x, torprmdict) - y
            array=optimize.least_squares(errfunc, pzero, jac='2-point', bounds=(-max_amp, max_amp), args=(rads(numpy.array(angle_list)),torprmdict, tor_energy_list))
            p1=array['x']
            #print('All parameters were removed, running least squares on main torsion type ')
            #print('Fitted Parameters ',array)

        # fill in torprmdict with the parameter estimates
        for chkclskey in torprmdict:
            for nfold in torprmdict[chkclskey]['prmdict']:
                parm  = p1[torprmdict[chkclskey]['prmdict'][nfold]]
                torprmdict[chkclskey]['prmdict'][nfold] = parm
            write_prm_dict[chkclskey] = torprmdict[chkclskey]['prmdict']
            # if not found, set as 0
            if write_prm_dict[chkclskey] == {}:
                write_prm_dict[chkclskey] = {1:0., 2:0., 3:0.}
            # Check if no arguments were fitted.
            if 'offset' in torprmdict[chkclskey]:
                if isinstance(p1,numpy.ndarray):
                    torprmdict[chkclskey]['offset'] = p1[-1]
                else:
                    torprmdict[chkclskey]['offset'] = p1
        #print('Final TorsionPrmDict For '+clskey,torprmdict)
        Sx = numpy.array(cls_angle_dict[clskey])
        fitfunc_dict[clskey] = fitfunc('eval',rads(Sx),torprmdict,debug=False)
        deriv_qm=numpy.gradient(tor_energy_list)
        weight=numpy.add(.6,numpy.absolute(deriv_qm))
        if len(fitfunc_dict[clskey])==len(tor_energy_list):
            def RMSDW(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.divide(numpy.subtract(fitfunc_dict[clskey],tor_energy_list),weight),c))))
            resultW=fmin(RMSDW,.5)
            minRMSDW=RMSDW(resultW[0])
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.subtract(fitfunc_dict[clskey],tor_energy_list),c))))
            result=fmin(RMSD,.5)
            minRMSD=RMSD(result[0])
        #else:
        #    print('lengths are not equal!',len(fitfunc_dict[clskey]),len(tor_energy_list))
        figfname = '%s-fit-%d-%d' % (molecprefix,b,c)
        figfname+='.png'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #Sx = numpy.array(cls_angle_dict[clskey])
        #fitfunc_dict[clskey] = fitfunc('eval',rads(Sx),torprmdict,debug=False)
        #print('fitfunc_dict[clskey]',fitfunc_dict[clskey],clskey)
        #print('tor_energy_list',tor_energy_list)
        l1, = ax.plot(Sx,fitfunc_dict[clskey],'r',label='Fit')
        l2, = ax.plot(Sx,tor_energy_list,'b',label='QM-MM1')
        plt.legend(handles=[l1,l2],loc='best')
        numprms=1 # offset parameter incldued with torsion force constant parameters
        for classkey in torprmdict:
           numprms+= len(torprmdict[classkey]['prmdict'].keys())
        # plot figure
        string=' , '.join(list(torprmdict.keys()))
        ax.text(0.05, 1.1, 'Torsions Being Fit =%s'%(string), transform=ax.transAxes, fontsize=10,verticalalignment='top')
        ax.text(0, -0.1, 'FoldNum=%s NumPrms=%s DataPts=%s RMSD(fit,QM-MM1) Wei=%s,Abs=%s'%(str(len(nfoldlist)),str(numprms),str(len(mm_energy_list)),round(minRMSDW,2),round(minRMSD,2)), transform=ax.transAxes, fontsize=10,verticalalignment='bottom')
        fig.savefig(figfname)
        # write parameter estimates to file
        write_arr_to_file(txtfname,[Sx,fitfunc_dict[clskey],tor_energy_list])
        #print "\n\n\n"
    return write_prm_dict,fitfunc_dict

def write_key_file(write_prm_dict,tmpkey1basename,tmpkey2basename):
    """
    Intent: Output the new key file based on parameters in write_prm_dict
    """
    tmpfh1 = open(tmpkey1basename, "r")
    tmpfh2 = open(tmpkey2basename, "w")
    for line in tmpfh1:
        m = re.search(r'torsion',line)
        if m is None:
            tmpfh2.write(line)
        else:
            linarr = line.split()
            cl = linarr[1:5]
            clskey = ' '.join(cl) # Order is fine (read from *.prm file)
            torline = line
            torvals = [float(ele) for ele in linarr[5:24:3]]
#            if clskey in write_prm_dict and \
#               torvals == [0.]*len(torvals):
            if clskey in write_prm_dict:
                torline = ' torsion %7s %4s %4s %4s   ' % (cl[0],cl[1],cl[2],cl[3])
                print(write_prm_dict[clskey])
                for (nfold, prm) in write_prm_dict[clskey].items():
                    torline += ' %7.3f %.1f %d' % (prm,foldoffsetlist[nfold - 1], nfold)
                torline += '\n'
            tmpfh2.write(torline)
    tmpfh1.close()
    tmpfh2.close()

def eval_rot_bond_parms(mol,fitfunc_dict,tmpkey1basename,tmpkey2basename):
    """
    Intent: 
    For each torsion whose parameters were fit for:
        Using the new parameters, find the new MM Energy vs. Dihedral Angle Profile
        Output the MM Energy (Pre-fit), MM Energy (Post-fit), and QM Energy profiles as 
        plots in the file: *energy*.png
    Ideally the profiles of MM Energy (Post-fit) will be much closer to the QM Energy profiles
    than the MM Energy (Pre-fit) profiles were. Look at the *png post running poltype to confirm.
    Input:
        mol: OBMol structure
        anglelist: dihedral angles, default 0-360, increments of 30
        fitfunc_dict: energy profile
        tmpkey1basename: Old key file, with old torsion parameters
        tmpkey2basename: New key file, with new torsion parameters
    Output:
        *energy*.png:
    Referenced By: process_rot_bond_tors
    Description:
    1. For each torsion whose parameters have been fit for (for each tor in torlist):
        a. Get each energy profile (MM pre, MM post, QM)
        b. Plot the profiles
    """
    global mm_tor_count
    # for each main torsion
    for tor in torlist:
        a,b,c,d = tor[0:4]
        key=str(b)+' '+str(c)
        anginc=rotbndtoanginc[key]
        #print('key',key,'anginc ',anginc)
        anglelist=range(0,360,anginc)
        torang = mol.GetTorsion(a,b,c,d)
        atmnuma = mol.GetAtom(a).GetAtomicNum()
        atmnumd = mol.GetAtom(d).GetAtomicNum()


        # clskey
        clskey = get_class_key(a, b, c, d)

        mm_energy_list = []
        mm_energy_list2 = []
        qm_energy_list = []

        # get the qm energy profile
        qm_energy_list,qang_list = compute_qm_tor_energy(a,b,c,d,torang,anglelist)
        tmpkeyfname = 'tmp.key'
        shutil.copy(tmpkey1basename, tmpkeyfname)
        # get the original mm energy profile
        mm_energy_list,mang_list,tor_e_list = compute_mm_tor_energy(mol,a,b,c,d,torang,anglelist,tmpkeyfname)
        mm_tor_count += 1
        # get the new mm energy profile (uses new parameters to find energies)
        mm2_energy_list,m2ang_list,tor_e_list2 = compute_mm_tor_energy(mol,a,b,c,d,torang,anglelist,tmpkey2basename)

        # remove angles for which energy was unable to be found
        del_ang_list = find_del_list(mm_energy_list,mang_list)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2)=prune_mme_error(del_ang_list,mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2)
        #dispvar("MMEb", mm_energy_list)
        #dispvar("MME2b", mm2_energy_list)
        del_ang_list = find_del_list(qm_energy_list,qang_list)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2)=prune_qme_error(del_ang_list,mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2)

        # normalize profiles
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]
        mm2_energy_list = [en - min(mm2_energy_list) for en in mm2_energy_list]

        # find the difference between the two energy due to torsion profiles 
        tordif_list = [e2-e1 for (e1,e2) in zip(tor_e_list,tor_e_list2)]
        # normalize
        tordif_list = [en - min(tordif_list) for en in tordif_list]
        # find the difference between the two mm energy profiles
        tordifmm_list = [e1+e2 for (e1,e2) in zip (tordif_list,mm_energy_list)]
        tordifmm_list = [en - min(tordifmm_list) for en in tordifmm_list]
        # TBC
        ff_list = [aa+bb for (aa,bb) in zip(mm_energy_list,fitfunc_dict[clskey])]


        deriv_qm=numpy.gradient(qm_energy_list)
        weight=numpy.add(.6,numpy.absolute(deriv_qm))
        if len(ff_list)==len(mm2_energy_list):
            def RMSDW(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.divide(numpy.subtract(mm2_energy_list,qm_energy_list),weight),c))))
            resultW=fmin(RMSDW,.5)
            minRMSDW=RMSDW(resultW[0])
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.subtract(mm2_energy_list,qm_energy_list),c))))
            result=fmin(RMSD,.5)
            minRMSD=RMSD(result[0])
        # output the profiles as plots
        figfname = "%s-energy-%d-%d" % (molecprefix,b,c)
        figfname+='.png'
        fig = plt.figure()
        #fig.subplots_adjust(right=0.75,left=0.05,top=0.95,bottom=0.05)
        ax = fig.add_subplot(111)
        # energy profiles: mm (pre-fit), mm (post-fit), qm
        line1, =ax.plot(mang_list,mm_energy_list,'g',label='MM1 (prefit)')
        line2, =ax.plot(m2ang_list,mm2_energy_list,'r',label='MM2 (postfit)')
        line3, =ax.plot(qang_list,qm_energy_list,'b',label='QM')

        ax.text(0.05, 1.1, 'RMSD(MM2,QM) Wei=%s,Abs=%s'%(round(minRMSDW,2),round(minRMSD,2)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
        # mm + fit
        line4, =ax.plot(mang_list,ff_list,'md-',label='MM1+Fit')

        plt.xlabel('Dihedral Angle')
        plt.ylabel('SP Energy (kcal/mol)')
        plt.legend(handles=[line1,line2,line3,line4],loc='best')
        fig = plt.gcf()
        plt.show()
        fig.savefig(figfname)
        txtfname = "%s-energy-%d-%d-%d-%d" % (molecprefix, a, b, c, d)


        txtfname = "%s-energy-%d-%d-%d-%d.txt" % (molecprefix, a, b, c, d)
        write_arr_to_file(txtfname,[mang_list,mm_energy_list,mm2_energy_list,qm_energy_list,tordif_list])

def gen_toromit_list():
    """
    Intent: if 'omittorsion2' is True, read in the *.toromit file to see which torsions 
    should not be scanned for
    Input: *.toromit is read in
    Output: 
        toromit_list: list of torsions that scanning should be omitted for
    Referenced By: main
    Description: Read in file, append information to toromit_list
    """
    global toromit_list
    toromitf = open(molecprefix+".toromit")
    for l in toromitf:
        toromit_list.append(sorttorsion([int(l.split()[0]), int(l.split()[1]), int(l.split()[2]), int(l.split()[3])]))
    toromitf.close()

def gen_fitrotbnds_list():
    """
    Intent: if 'fitrotbnds' is True, read in the *.fitrotbnds file to see which rotatable
    should be fitted with multiple torsions per bond
    Input: *.fitrotbnd is read in
    Output: 
        fitrotbndslist: list of rotatable bonds to fit torsions to, does not include hydrogens by default
    Referenced By: main
    Description: Read in file, append information to toromit_list
    """
    global fitrotbndslist
    tempf= open(molecprefix+".fitrotbnds")
    for l in tempf:
        fitrotbndslist.append([int(l.split()[0]), int(l.split()[1])])
    tempf.close()

def sorttorsion(keylist):
    """
    Intent: Sort the torsion key by sorting the two outer terms and then two inner terms
    (e.g. 4-2-1-3 -> 3-1-2-4)
    Input: 
        keylist: torsion key to be sorted
    Output:
        keylist is updated
    Referenced By: get_torlist, gen_toromit_list
    Description: -
    """
    if(keylist[1] > keylist[2] or (keylist[1] == keylist[2] and keylist[0] > keylist[3])):
        temp1 = keylist[1]
        keylist[1] = keylist[2]
        keylist[2] = temp1 
        temp2 = keylist[3]
        keylist[3] = keylist[0]
        keylist[0] = temp2
    return keylist

# Fit torsion parameters for rotatable bonds
def process_rot_bond_tors(mol):
    """
    Intent: Fit torsion parameters for torsions about rotatable bonds 
    Input:
        mol: OBMol structure
    Output:
        *.key_5 is written out, with updated torsion parameters
    Referenced By: main
    Description:
    1. Get the QM Energy vs. Dihedral Angle and (Initial/Pre-Fit) MM Energy vs. Dihedral Angle 
       profiles for each rotatable bond. 
       Store these profiles in 'cls_qm_engy_dict' and 'cls_mm_engy_dict'.
    2. Use these profiles to fit for the torsion parameters by calling 'fit_rot_bond_tors'
    3. Evaluate the new parameters output by the fitting fuction and output informational plots 
       by calling 'eval_rot_bond_parms'
    4. Write out the new keyfile (*.key_5) with these new torsion parameters
    """
    global mm_tor_count

    #create list from 0 - 360 in increments of 30
    #anglist = range(0,360,30)
    tordir = 'qm-torsion'
    tmpkey1basename = 'tinker.key'
    tmpkey2basename = 'tinker.key_2'
    tmpkey1fname = tordir + '/' + tmpkey1basename
    assert os.path.isdir(tordir), \
       "ERROR: Directory '%s' does not exist" % tordir
    # copy *.key_4 to the directory qm-torsion
    shutil.copy(key4fname, tmpkey1fname)
    # change directory to qm-torsion
    os.chdir(tordir)

    # Group all rotatable bonds with the same classes and identify
    # the torsion parameters that need to be fitted.

    # For each rotatable bond, get torsion energy profile from QM
    # and MM (with no rotatable bond torsion parameters)
    # Get QM and MM (pre-fit) energy profiles for torsion parameters
    mm_tor_count += 1
    cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict = get_qmmm_rot_bond_energy(mol,tmpkey1basename)

    # if the fit has not been done already
    if torkeyfname is None:
        # do the fit
        write_prm_dict,fitfunc_dict = fit_rot_bond_tors(
            mol,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict)
        mm_tor_count += 1
        # write out new keyfile
        write_key_file(write_prm_dict,tmpkey1basename,tmpkey2basename)
    else:
        shutil.copy('../' + torkeyfname,tmpkey2basename)
    # evaluate the new parameters
    eval_rot_bond_parms(
        mol,fitfunc_dict,tmpkey1basename,tmpkey2basename)
    shutil.copy(tmpkey2basename,'../' + key5fname)
    os.chdir('..')



def CheckRMSD():
    for line in open(superposeinfile,'r'):
        if os.stat(superposeinfile).st_size > 5:
            if 'Root Mean' in line:
                RMSD=''
                for e in line:
                    if e.isdigit() or e=='.':
                        RMSD+=e
            if float(RMSD)>maxRMSD:
                print('Warning: RMSD of QM and MM optimized structures is high, RMSD = ',RMSD)
                logfh.write('Warning: RMSD of QM and MM optimized structures is high, RMSD = '+ RMSD+' Tolerance is '+str(maxRMSD)+' kcal/mol '+'\n') # now report all torsions as bad
                logfh.flush()
                os.fsync(logfh.fileno())
                

def CheckRMSPD():
    temp=open('RMSPD.txt','r')
    for line in temp.readlines():
        if 'Root Mean Square Potential Difference :' in line:
            RMSPD=line.split(':')[1].strip()
    temp.close()
    if float(RMSPD)>maxRMSPD:
        print('Warning: RMSPD of QM and MM optimized structures is high, RMSPD = ',RMSPD)
        logfh.write('Warning: RMSPD of QM and MM optimized structures is high, RMSPD = '+ RMSPD+' Tolerance is '+str(maxRMSPD)+' kcal/mol '+'\n') # now report all torsions as bad
        logfh.flush()
        os.fsync(logfh.fileno())

def CheckBondConnectivity(mol,optmol):
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
            RaiseConnectivityError(diff,idxset)

def RaiseConnectivityError(diff,idxset):
    print('Error! The bond connectivity before and after structure optimization is different')
    logfh.write('Error! The bond connectivity before and after structure optimization is different'+'\n')
    for atmidx in diff:
        print('The atom index '+str(atmidx)+' from structure '+idxset+' does not have the same connectivity before and after structure optimization')
        logfh.write('The atom index '+str(atmidx)+' from structure '+idxset+' does not have the same connectivity before and after structure optimization'+'\n')
    sys.exit() 

              
#POLTYPE BEGINS HERE
def main():

    # log file object
    global logfh
    # An array that maps atom ids to the symmetry class they belong to
    global symmetryclass
    # the OBMol object
    global mol
    global localframe1
    global localframe2
    global canonicallabel
    global rotbndlist
    global torlist
    global use_psi4
    global use_psi4SPonly
    global nfoldlist
    global dmamethod
    dmamethod="MP2"
    nfoldlist=range(1,foldnum+1)

    # Initialization. 
    # Setting flags, setting up directories, setting up files
    parse_options(sys.argv)
    copyright()
    initialize()
    init_filenames()
    
    # Use openbabel to create a 'mol' object from the input molecular structure file. 
    # Openbabel does not play well with certain molecular structure input files,
    # such as tinker xyz files. (normal xyz are fine)
    assert os.path.isfile(molstructfname), "Error: Cannot open " + molstructfname

    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(molstructfname)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(mol, molstructfname) 


    chg=mol.GetTotalCharge()
    print('total charge ',chg)
    if chg!=0:
        toroptpcm=True
        optpcm=True
        torsppcm=True
        sppcm=True

    # Begin log. *-poltype.log
    logfh = open(logfname,"a")
    if use_psi4==True:
        try:
            import psi4
            print("Psi4 detected; it will be used as the QM engine")
            use_psi4 = True
        except:
            print("Psi4 not detected; attempting to use Gaussian instead")
            use_psi4 = False


    if use_psi4:
        psi4.set_num_threads(int(numproc))
        psi4.set_memory(maxmem)
    logfh.write("Running on host: " + gethostname() + "\n")

    # QM calculations are done here
    # First the molecule is optimized. (-opt) 
    # This optimized molecule is stored in the structure optmol
    # Then the electron density matrix is found (-dma)
    # This is used by GDMA to find multipoles
    # Then information for generating the electrostatic potential grid is found (-esp)
    # This information is used by cubegen
    if use_psi4 or use_psi4SPonly:
        optmol = run_psi4(mol)
    else:
        optmol = run_gaussian(mol)

    CheckBondConnectivity(mol,optmol)

    # End here if qm calculations were all that needed to be done 
    if qmonly:
        now = time.strftime("%c",time.localtime())
        logfh.write(now + " poltype QM-only complete.\n")
        logfh.close()
        sys.exit(0)

    # Initializing arrays
    
    canonicallabel = [ 0 ] * mol.NumAtoms()
    localframe1 = [ 0 ] * mol.NumAtoms()
    localframe2 = [ 0 ] * mol.NumAtoms()

    # Finds the symmetry class for each atom
    # For example in the molecule ethanol: CH3-CH2-OH
    # The 3 Hydrogens bound to the first carbon all belong to the same symmetry class
    # The 2 Hydrogens bound to the second carbon all belong to a second symmetry class
    # The rest of the atoms all belong to their own individual symmetry classes
    # Many babel tools are used in finding the symmetry classes
    gen_canonicallabels(mol)
   
    # scaling of multipole values for certain atom types
    # checks if the molecule contains any atoms that should have their multipole values scaled
    scalelist = process_types(mol)
    
    # if the omittorsion2 flag has been selected, poltype will know not to scan for the torsions
    # of certain rotatble bonds
    if(omittorsion2):
        gen_toromit_list()


    if(fitrotbnds):
        gen_fitrotbnds_list()
    # Find rotatable bonds for future torsion scans
    (torlist, rotbndlist) = get_torlist(mol)
    torlist = get_torlist_opt_angle(optmol, torlist)
    #print('torsion list is ',torlist)
    #print('rotbndlist ',rotbndlist)
    oblist, rotbndlist_forvalence = gen_valinfile(mol)
    #print('rotbndlist_forvalence ',rotbndlist_forvalence)
    global rotbndtoanginc
    rotbndtoanginc={}
    # if there are multiple torsions to fit per rotatable bond, make sure there are enough angles for QM profile to do fitting
    for key in rotbndlist:
        keylist=rotbndlist[key]
        nonhydlist=[]
        for tor in keylist:
            a2,b2,c2,d2=tor[0:4]
            obaa = mol.GetAtom(a2)
            obab = mol.GetAtom(b2)
            obac = mol.GetAtom(c2)
            obad = mol.GetAtom(d2)
            tpdkey = get_class_key(a2, b2, c2, d2)
            #print('Key',tpdkey,obaa.GetAtomicNum(),obab.GetAtomicNum(),obac.GetAtomicNum(),obad.GetAtomicNum())
            append=False
            if obaa.GetAtomicNum() != 1 and obad.GetAtomicNum() != 1:
                append=True
                #print('appending because of no H')
            if append==True and tpdkey not in nonhydlist:
                nonhydlist.append(tpdkey)
        #print('Number of folds',len(nfoldlist),'Number of Parameters to Fit',len(nonhydlist)+1)
        if tordatapointsnum==None:
            ang=360/(len(nfoldlist)*(len(nonhydlist)+1)) # offset parameter is the +1
            if ang> 30:
                ang=30
                #print('Angle Increment is greater than 30, so default is 30 ')
        else:
            ang=360/tordatapointsnum
        rotbndtoanginc[key]=ang
  
    
    # Obtain multipoles from Gaussian fchk file using GDMA

    if not os.path.isfile(gdmafname):
        run_gdma()

    # Set up input file for poledit
    # find multipole local frame definitions 
    lfzerox,atomindextoremovedipquad,atomtypetospecialtrace,atomindextoremovedipquadcross = gen_peditinfile(mol)
    
    
    if (not os.path.isfile(xyzfname) or not os.path.isfile(keyfname)):
        # Run poledit
        cmdstr = peditexe + " 1 " + gdmafname +' '+paramhead+ " < " + peditinfile
        call_subsystem(cmdstr)
        # Add header to the key file output by poledit
        prepend_keyfile(keyfname)
    # post process local frames written out by poledit
    post_proc_localframes(keyfname, lfzerox,atomindextoremovedipquad,atomindextoremovedipquadcross)
    # generate the electrostatic potential grid used for multipole fitting
    gen_esp_grid(optmol)

    # Average multipoles based on molecular symmetry
    # Does this using the script avgmpoles.pl which is found in the poltype directory
    # Atoms that belong to the same symm class will now have only one common multipole definition
    if uniqidx:
        shutil.copy(keyfname, key2fname)
        prepend_keyfile(key2fname)
    elif ((not os.path.isfile(xyzoutfile) or
            not os.path.isfile(key2fname)) and
            not uniqidx):
        # gen input file
        gen_avgmpole_groups_file()
        # call avgmpoles.pl
        avgmpolecmdstr = avgmpolesexe + " " + keyfname + " " + xyzfname + " " + grpfname + " " + key2fname + " " + xyzoutfile + " " + str(prmstartidx)
        call_subsystem(avgmpolecmdstr)
        prepend_keyfile(key2fname)

    if espfit:
        # Optimize multipole parameters to QM ESP Grid (*.cube_2)
        # tinker's potential utility is called, with option 6.
        # option 6 reads: 'Fit Electrostatic Parameters to a Target Grid'
        
        if not os.path.isfile(key3fname):

            #optmpolecmd = potentialfit + " 6 " + xyzoutfile + " -k " + key2fname + " " + qmesp2fname + " N 1"
            optmpolecmd = potentialexe + " 6 " + xyzoutfile + " -k " + key2fname + " " + qmesp2fname + " N 0.1"
            call_subsystem(optmpolecmd)
    else:
        shutil.copy(key2fname, key3fname)
    # Remove header terms from the keyfile
    rm_esp_terms_keyfile(key3fname)

    logfh.write("\n")
    logfh.write("=========================================================\n")
    logfh.write("Electrostatic Potential Comparision\n\n")
    logfh.flush()
    os.fsync(logfh.fileno())
    cmd=potentialexe + ' 5 ' + xyzoutfile + ' ' + qmesp2fname + ' N > RMSPD.txt'
    call_subsystem(cmd)
    CheckRMSPD()
    if not os.path.isfile(key4fname):
        shutil.copy(key3fname, key4fname)
        
        # Multipoles are scaled if needed using the scale found in process_types
        post_process_mpoles(key4fname, scalelist)
        
        # Now that multipoles have been found
        # Other parameters such as opbend, vdw, etc. are found here using a look up table
        # Part of the look up table is here in poltype.py 
        # Most of it is in the file valence.py found in the poltype directory

        # Finds aromatic carbons and associated hydrogens and corrects polarizability
        # Find opbend values using a look up table
        # Outputs a list of rotatable bonds (found in get_torlist) in a form usable by valence.py
        oblist, rotbndlist_forvalence = gen_valinfile(mol)
        print('rotbndlist_forvalence type number',rotbndlist_forvalence)

        # Map from idx to symm class is made for valence.py
        idxtoclass=[]
        for i in range(mol.NumAtoms()):
            idxtoclass.append(get_class_number(i+1))
        v = valence.Valence(output_format,logfname)
        v.setidxtoclass(idxtoclass)
        dorot = True
        # valence.py method is called to find parameters and append them to the keyfile
        v.appendtofile(key4fname, optmol, oblist, dorot, rotbndlist_forvalence) 

        

    # Torsion scanning then fitting. *.key_5 will contain updated torsions
    # default, parmtors = True
    if dontdotor==True:
        sys.exit()
    if (parmtors):
        # torsion scanning
        gen_torsion(optmol,torsionrestraint)
    if (parmtors):
        # torsion fitting
        if dontdotorfit==True:
            sys.exit()
        process_rot_bond_tors(optmol)
    else:
        shutil.copy(key4fname,key5fname)

    #If the output format is set to tinker 4, a key_6
    #is created with parameters in the tinker 4 format
    if output_format == 4:
        key5 = open(key5fname)
        key6 = open(keyfname+"_6","w")
        for line in key5:
            if 'polarize' in line:
                ln = ""
                for i in range(len(line.split())):
                    if i != 3:
                        if i == 1:
                            ln += line.split()[i] + "                          "
                        else:
                            ln += line.split()[i] + "  "
                ln += "\n"
                key6.write(ln)
            else:
                key6.write(line)

    gen_tinker5_to_4_convert_input(mol, amoeba_conv_spec_fname)



    # A series of tests are done so you one can see whether or not the parameterization values
    # found are acceptable and to what degree
    logfh.write("\n")
    logfh.write("=========================================================\n")
    logfh.write("Minimizing structure\n\n")
    logfh.flush()
    os.fsync(logfh.fileno())

    cmd='cp ' + xyzoutfile + ' ' + tmpxyzfile
    os.system(cmd)
    cmd='cp ' + key5fname + ' ' + tmpkeyfile
    os.system(cmd)
    cmd = minimizeexe+' -k '+tmpkeyfile+' '+tmpxyzfile+' 0.1'
    call_subsystem(cmd)
    logfh.flush()
    os.fsync(logfh.fileno())

    logfh.write("\n")
    logfh.write("=========================================================\n")
    logfh.write("QM Dipole moment\n\n")
    logfh.flush()
    os.fsync(logfh.fileno())
    grepcmd = 'grep -A7 "Dipole moment" ' + logespfname
    call_subsystem(grepcmd)
    logfh.flush()
    os.fsync(logfh.fileno())

    logfh.write("\n")
    logfh.write("=========================================================\n")
    logfh.write("MM Dipole moment\n\n")
    logfh.flush()
    os.fsync(logfh.fileno())
    cmd=analyzeexe + ' ' +  xyzoutfile + ' em | grep -A11 Charge'
    call_subsystem(cmd,iscritical=True)

    gen_superposeinfile()
    logfh.write("\n")
    logfh.write("=========================================================\n")
    logfh.write("Structure RMSD Comparison\n\n")
    logfh.flush()
    os.fsync(logfh.fileno())
    cmd = superposeexe + ' ' + xyzoutfile + ' ' + tmpxyzfile + '_2' + ' < ' + superposeinfile
    call_subsystem(cmd)
    CheckRMSD()

    logfh.write("\n")
    logfh.write("=========================================================\n")
    logfh.write("Electrostatic Potential Comparision\n\n")
    logfh.flush()
    os.fsync(logfh.fileno())
    cmd=potentialexe + ' 5 ' + xyzoutfile + ' ' + qmesp2fname + ' N > RMSPD.txt'
    call_subsystem(cmd)
    CheckRMSPD()
    logfh.write('Poltype Job Finished'+'\n')
    logfh.flush()
    os.fsync(logfh.fileno())

if __name__ == '__main__':
    main()
