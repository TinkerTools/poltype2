
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
##################################################################

import os
import sys
from socket import gethostname
import subprocess
import openbabel
import shutil
import time
import getopt
from PoltypeModules import valence
from PoltypeModules import torsiongenerator as torgen
from PoltypeModules import modifiedresidues as modres
from PoltypeModules import symmetry as symm
from PoltypeModules import torsionfit as torfit
from PoltypeModules import optimization as opt
from PoltypeModules import electrostaticpotential as esp
from PoltypeModules import multipole as mpole



class PolarizableTyper():

    def __init__(self):

        self.prmstartidx = 401
        self.numproc = 1
        self.maxmem = "700MB"
        self.maxdisk = "100GB"
        self.gausdir = None
        self.gdmadir = None
        self.tinkerdir = None
        self.scratchdir = "/scratch"
        self.paramhead = sys.path[0] + "/amoebabio18_header.prm"
        self.babelexe = "babel"
        self.gausexe =  "g09"
        self.formchkexe =  "formchk"
        self.cubegenexe =  "cubegen"
        self.gdmaexe = "gdma"
        self.avgmpolesexe = sys.path[0] + "/avgmpoles.pl"
        self.peditexe = "poledit.x"
        self.potentialexe = "potential.x"
        self.valenceexe = "valence.x"
        self.minimizeexe = "minimize.x"
        self.analyzeexe = "analyze.x"
        self.superposeexe = "superpose.x"
        self.defopbendval = 0.20016677990819662
        self.Hartree2kcal_mol = 627.5095
        self.optbasisset = "6-31G*"               # This is good; bigger is better but slower
        self.toroptbasisset = "6-31G*"            # Same as above
        self.dmabasisset = "6-311G**"              # This can not be changed
        self.popbasisset = "6-31G*"               # We typically dont do pop analysis unless we need bond orders from Gaussain (not sure PSI4 has this).
        self.espbasisset = "6-311++G(2d,2p)"          # This can be from 6-311++G**, 6-311++G(2d,2p), or aug-cc-pvtz. For small molecules like yours, use the last two.
        self.torspbasisset = "6-311++G**"         # This can be 6-31+G* or bigger, depending on the cost (we have to do this for many conformations)
        self.optmethod='MP2'                      # MP2 for small molecules, DFT for large (wB97XD)
        self.toroptmethod='HF'                    # HF or DFT in gas. HF may be needed for large molecules unless we break down to fragments for torsion fitting.
        self.torspmethod='MP2'                    # DFT or MP2 in PCM. DFT may be good enough (comparison would be nice). 
        self.dmamethod='MP2'                      # MP2
        self.espmethod='MP2'                      # MP2
        self.qmonly = False
        self.espfit = True
        self.parmtors = True
        self.foldnum=3
        self.nfoldlist = range(1,self.foldnum+1)
        self.foldoffsetlist = [ 0.0, 180.0, 0.0, 0.0, 0.0, 0.0 ]
        self.torlist = []
        self.rotbndlist = []
        self.fitrotbndslist=[]
        self.mm_tor_count = 1
        self.output_format = 5
        self.do_tor_qm_opt = True
        self.fitrotbnds=False # list of rotatable bonds to be fitted with multiple torsions around bond (without hydrogen), boolean that sets option to read from this list as False first
        self.maxRMSD=.1
        self.maxRMSPD=1
        self.tordatapointsnum=None
        self.gentorsion=False # boolean to specify when program is running gaussian for torsion
        self.gaustorerror=False # boolean to specify when gaussian crashes during gaussian opt or SP for torsion, this way can exclude some points for fitting 
        self.torsionrestraint=.1
        self.onlyrotbndlist=[] # initial empty list default, do not spin only one bond
        self.rotalltors=False
        self.dontdotor=False
        self.dontdotorfit=False
        self.onlyrotbnd=False
        self.toroptpcm=False
        self.optpcm=False
        self.torsppcm=False
        self.use_psi4=False
        self.use_psi4SPonly=False
        self.freq=False
        self.postfit=False # just a boolean specifying where in torsion fitting process is
        self.bashrcpath=None
        self.amoebabioprmpath=None
        self.libpath=sys.path[0] + "/lib.bio18_conv1.txt"
        self.SMARTSToTypelibpath=sys.path[0]+'/SMARTSToTypeLib.txt'
        self.ModifiedResiduePrmPath=sys.path[0]+'/ModifiedResidue.prm'
        self.modifiedproteinpdbname=None
        self.unmodifiedproteinpdbname=None
        self.mutatedsidechain=None
        self.mutatedresiduenumber=None
        self.modifiedresiduepdbcode=None
        self.optmaxcycle=400
        self.torkeyfname=None
        self.gausoptcoords=''
        self.uniqidx=False 
        self.helpfile=sys.path[0]+'/README.HELP'
        self.versionfile=sys.path[0]+'/README.VERSION' 
        opts, xargs = getopt.getopt(sys.argv[1:],'h',["help"])


        for o, a in opts:
            if o in ("-h", "--help"):
                self.copyright()
                self.usage()
                sys.exit(2)
        
        
        temp=open(os.getcwd()+r'/'+'poltype.ini','r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '#' not in line:
                if '=' in line:
                    linesplit=line.split('=',1)
                    a=linesplit[1].replace('\n','').rstrip().lstrip()
                    newline=linesplit[0]
                else:
                    newline=line

                if "rotalltors" in newline:
                    self.rotalltors = True
                elif "externalapi" in newline:
                    self.externalapi=a
                elif "gausoptcoords" in newline:
                    self.gausoptcoords = a
                elif "toroptbasisset" in newline:
                    self.toroptbasisset = a
                elif "modifiedresiduepdbcode" in newline:
                    self.modifiedresiduepdbcode = a
                elif "mutatedsidechain" in newline:
                    self.mutatedsidechain = a
                elif "mutatedresiduenumber" in newline:
                    self.mutatedresiduenumber = a
                elif "unmodifiedproteinpdbname" in newline:
                    self.unmodifiedproteinpdbname = a
                    self.molstructfname='ModifiedRes.sdf'
                elif "dmamethod" in newline:
                    self.dmamethod = a
                elif "bashrcpath" in newline:
                    self.bashrcpath = a
                elif "modifiedproteinpdbname" in newline:
                    self.modifiedproteinpdbname = a
                    self.molstructfname='ModifiedRes.sdf'
                elif "amoebabioprmpath" in newline:
                    self.amoebabioprmpath = a
                elif "structure" in newline:
                    self.molstructfname = a
                elif "torsppcm" in newline:
                    self.torsppcm = True
                elif "freq" in newline:
                    self.freq = True
                elif "optpcm" in newline:
                    self.optpcm = True
                elif "toroptpcm" in newline:
                    self.optpcm = True
                elif "use_psi4" in newline:
                    self.use_psi4 = True
                elif "use_psi4SPonly" in newline:
                    self.use_psi4SPonly = True
                elif "dontdotor" in newline:
                    self.dontdotor = True
                elif "dontdotorfit" in newline:
                    self.dontdotorfit = True
                elif "optmaxcycle" in newline:
                    self.optmaxcycle = a
                elif "torsionrestraint" in newline:
                    self.torsionrestraint=float(a)
                elif "foldnum" in newline:
                    self.foldnum=int(a)
                elif "tordatapointsnum" in newline:
                    self.tordatapointsnum=int(a)
                elif "fitrotbndslist" in newline:
                    self.fitrotbndslist=a.split(',')
                    templist=[]
                    for ele in self.fitrotbndslist:
                        nums=ele.lstrip().rstrip().split()
                        temp=[]
                        for e in nums:
                            temp.append(int(e))
                        templist.append(temp)
                    self.fitrotbndslist=templist
                    self.fitrotbnds=True
                elif "optmethod" in newline:
                    self.optmethod = a
                elif "espmethod" in newline:
                    self.espmethod = a
                elif "torspmethod" in newline:
                    self.torspmethod = a
                elif "toroptmethod" in newline:
                    self.toroptmethod = a
                elif "numproc" in newline:
                    self.numproc = a
                elif "maxmem" in newline:
                    self.maxmem = a
                elif "maxdisk" in newline:
                    self.maxdisk = a
                elif "atmidx" in newline:
                    self.prmstartidx = int(a)
                elif "optbasisset" in newline:
                    self.optbasisset = a
                elif "dmabasisset" in newline:
                    self.dmabasisset = a
                elif "popbasisset" in newline:
                    self.popbasisset = a
                elif "espbasisset" in newline:
                    self.espbasisset = a
                elif "torspbasisset" in newline:
                    self.torspbasisset = a
                elif "optlog" in newline:
                    self.logoptfname = a
                elif "dmalog" in newline:
                    self.logdmafname = a
                elif "esplog" in newline:
                    self.logespfname = a
                elif "dmafck" in newline:
                    self.fckdmafname = a
                elif "espfck" in newline:
                    self.fckespfname = a
                elif "dmachk" in newline:
                    self.chkdmafname = a
                elif "espchk" in newline:
                    self.chkespfname = a
                elif "formchk" in newline:
                    self.fname = a
                elif "gdmaout" in newline:
                    self.gdmafname = a
                elif "gbindir" in newline:
                    self.gausdir = a
                elif "qmonly" in newline:
                    self.qmonly = True
                elif "omit-espfit" in newline:
                    self.espfit = False
                elif "dont-tor-qm-opt" in newline:
                    self.do_tor_qm_opt = False
                elif "uniqidx" in newline:
                    self.uniqidx = True
                elif "help" in newline:
                    self.copyright()
                    self.usage()
                    sys.exit(2)
                elif "onlyrotbnd" in newline: # comma seperated list 
                    self.onlyrotbndlist=[i.lstrip().rstrip() for i in a.split(',')]
                else:
                    print('Unrecognized '+line)
                    self.usage()
                    sys.exit()
        self.copyright()
        self.initialize()
        self.init_filenames() 
        # Use openbabel to create a 'mol' object from the input molecular structure file. 
        # Openbabel does not play well with certain molecular structure input files,
        # such as tinker xyz files. (normal xyz are fine)
    
        obConversion = openbabel.OBConversion()
        self.mol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(self.molstructfname)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(self.mol, self.molstructfname) 
    
        # Begin log. *-poltype.log
        self.logfh = open(self.logfname,"a")
    
        self.mol=self.CheckIsInput2D(self.mol,obConversion)
    
        chg=self.mol.GetTotalCharge()
        if chg!=0:
            self.toroptpcm=True
            self.optpcm=True
            self.torsppcm=True
     
        self.main()

    def WriteToLog(self,string):
        now = time.strftime("%c",time.localtime())
        self.logfh.write(now+' '+string+'\n')
        self.logfh.flush()
        os.fsync(self.logfh.fileno())

        
    def which(self,program,pathlist=os.environ["PATH"]):
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


    def initialize(self):
        """
        Intent: Initialize all paths to needed executables
        Input:
        Output:
        Referenced By: main
        Description: -
        """
    
        if not self.use_psi4:
            if (self.gausdir is not None):
                if self.which(os.path.join(self.gausdir,"g09")) is not None:
                    self.gausexe    = os.path.join(self.gausdir,"g09")
                    self.formchkexe = os.path.join(self.gausdir,self.formchkexe)
                    self.cubegenexe = os.path.join(self.gausdir,self.cubegenexe)
                elif self.which(os.path.join(self.gausdir,"g03")) is not None:
                    self.gausexe    = os.path.join(self.gausdir,"g03")
                    self.formchkexe = os.path.join(self.gausdir,self.formchkexe)
                    self.cubegenexe = os.path.join(self.gausdir,self.cubegenexe)
                else:
                    print("ERROR: Invalid Gaussian directory: ", self.gausdir)
                    sys.exit(1)
            else:
                if self.which("g09") is not None:
                    self.gausexe    = "g09"
                elif self.which("g03") is not None:
                    self.gausexe    = "g03"
                else:
                    print("ERROR: Cannot find Gaussian executable in $PATH. Please install Gaussian or specify Gaussian directory with --gbindir flag.")
                    sys.exit(1)
    
    
                
    
    
        cmdstr=self.analyzeexe+' '+sys.path[0]+r'/'+'water.xyz'+' '+'-k'+' '+sys.path[0]+r'/'+'water.key'+' '+'e'+'>'+' '+'version.out'
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
            self.tinkerdir = os.environ["TINKERDIR"]
            self.peditexe = os.path.join(tinkerdir,self.peditexe)
            self.potentialexe = os.path.join(tinkerdir,self.potentialexe)
            self.valenceexe = os.path.join(tinkerdir,self.valenceexe)
            self.minimizeexe = os.path.join(tinkerdir,self.minimizeexe)
            self.analyzeexe = os.path.join(tinkerdir,self.analyzeexe)
            self.superposeexe = os.path.join(tinkerdir,self.superposeexe)
    
        if (not self.which(self.analyzeexe)):
            print("ERROR: Cannot find TINKER analyze executable")
            sys.exit(2)
    
        if not self.use_psi4 and not self.use_psi4SPonly:
            if ("GDMADIR" in os.environ):
                self.gdmadir = os.environ["GDMADIR"]
                self.gdmaexe = os.path.join(self.gdmadir,self.gdmaexe)
    
            if (not self.which(self.gdmaexe)):
                print("ERROR: Cannot find GDMA executable")
                sys.exit(2)
    
            if ("GAUSS_SCRDIR" in os.environ):
                self.scratchdir = os.environ["GAUSS_SCRDIR"]
                if not os.path.isdir(self.scratchdir):
                    os.mkdir(self.scratchdir)
    
    
            if (not self.which(self.scratchdir)):
                print("ERROR: Cannot find Gaussian scratch directory")
                sys.exit(2)
    
    
    
    def init_filenames (self):
        """
        Intent: Initialize file names
        Input:
        Output:
        Referenced By: main
        Description: -
        """
    
    
        if ("GAUSS_SCRDIR" in os.environ):
            self.scratchdir = os.environ["GAUSS_SCRDIR"].rstrip('//')
            if not os.path.isdir(self.scratchdir):
                os.mkdir(self.scratchdir)
    
        head, self.molstructfname = os.path.split(self.molstructfname)
        self.molecprefix =  os.path.splitext(self.molstructfname)[0]
        self.logfname = self.assign_filenames ( "logfname" , "-poltype.log")
        self.chkname = self.assign_filenames ( "chkname" , ".chk")
        self.fname = self.assign_filenames ( "fname" , ".fchk")
        self.gausfname = self.assign_filenames ( "gausfname" , ".log")
        self.gausoptfname = self.assign_filenames ( ".gausoptfname" , "-opt.log")
        self.gdmafname = self.assign_filenames ( "gdmafname" , ".gdmaout")
        self.keyfname = self.assign_filenames ( "keyfname" , ".key")
        self.xyzfname = self.assign_filenames ( "xyzfname" , ".xyz")
        self.peditinfile = self.assign_filenames ( "peditinfile" , "-peditin.txt")
        self.superposeinfile = self.assign_filenames ( "superposeinfile" , "-superin.txt")
        self.espgrdfname = self.assign_filenames ( "espgrdfname" , ".grid")
        self.qmespfname = self.assign_filenames ( "qmespfname" , ".cube")
        self.qmesp2fname = self.assign_filenames ( "qmesp2fname" , ".pot")
        self.grpfname = self.assign_filenames ( "grpfname" , "-groups.txt")
        self.key2fname = self.assign_filenames ( "key2fname" , ".key_2")
        self.key3fname = self.assign_filenames ( "key3fname" , ".key_3")
        self.key4fname = self.assign_filenames ( "key4fname" , ".key_4")
        self.key5fname = self.assign_filenames ( "key5fname" , ".key_5")
        self.xyzoutfile = self.assign_filenames ( "xyzoutfile" , ".xyz_2")
        self.scrtmpdir = self.scratchdir.rstrip('//') + '/Gau-' + self.molecprefix
        self.tmpxyzfile = 'ttt.xyz'
        self.tmpkeyfile = 'ttt.key'
        self.comtmp = self.assign_filenames ( "comtmp" , "-tmp.com")
        self.comoptfname = self.assign_filenames ( "comoptfname" , "-opt.com")
        self.chkoptfname = self.assign_filenames ( "chkoptfname" , "-opt.chk")
        self.fckoptfname = self.assign_filenames ( "fckoptfname" , "-opt.fchk")
        self.logoptfname = self.assign_filenames ( "logoptfname" , "-opt.log")
        self.compopfname = self.assign_filenames ( "compopfname" , "-pop.com")
        self.chkpopfname = self.assign_filenames ( "chkpopfname" , "-pop.chk")
        self.fckpopfname = self.assign_filenames ( "fckpopfname" , "-pop.fchk")
        self.logpopfname = self.assign_filenames ( "logpopfname" , "-pop.log")
        self.comdmafname = self.assign_filenames ( "comdmafname" , "-dma.com")
        self.chkdmafname = self.assign_filenames ( "chkdmafname" , "-dma.chk")
        self.fckdmafname = self.assign_filenames ( "fckdmafname" , "-dma.fchk")
        self.logdmafname = self.assign_filenames ( "logdmafname" , "-dma.log")
        self.comespfname = self.assign_filenames ( "comespfname" , "-esp.com")
        self.chkespfname = self.assign_filenames ( "chkespfname" , "-esp.chk")
        self.fckespfname = self.assign_filenames ( "fckespfname" , "-esp.fchk")
        self.logespfname = self.assign_filenames ( "logespfname" , "-esp.log")


    def assign_filenames (self,filename,suffix):
        if filename in globals():
            return eval(filename)
        else:
            return self.molecprefix + suffix
    
    
    def copyright (self):
        temp=open(self.versionfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            print(line)
    
    def usage (self):
        temp=open(self.helpfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            print(line)
    
    
    
    def CheckIsInput2D(self,mol,obConversion):
        is2d=True
        for atom in openbabel.OBMolAtomIter(mol):
            zcoord=atom.GetZ()
            if zcoord!=0:
                is2d=False
            
        if is2d==True: 
            molprefix=self.molstructfname.split('.')[0]
            newname=molprefix+'_3D'+'.sdf'
            ext=self.molstructfname.split('.')[1]
            cmdstr='babel'+' '+'-i'+ext+' '+self.molstructfname+' '+'--gen3d'+' '+'-osdf'+' '+'-O'+' '+newname
            os.system(cmdstr)
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(newname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol,newname)
    
        return mol
    
    

    def call_subsystem(self,cmdstr,wait=False):
        curdir=os.getcwd()
        self.WriteToLog(" Calling: " + cmdstr)
        p = subprocess.Popen(cmdstr, shell=True,stdout=self.logfh, stderr=self.logfh)
        if wait==True:
            (output,err)=p.communicate()
            rcode = p.wait()



                  
    #POLTYPE BEGINS HERE
    def main(self):
    
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
            knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check=modres.GenerateModifiedProteinPoltypeInput(self)
            self.molstructfname=molname
            head, self.molstructfname = os.path.split(self.molstructfname)
            self.molecprefix =  os.path.splitext(self.molstructfname)[0]

            assert os.path.isfile(self.molstructfname), "Error: Cannot open " + self.molstructfname
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None): # if already have core parameters in modified prm database then dont regenerate parameters
            if check==False:
                self.GenerateParameters(self.mol)
        else:
            self.GenerateParameters(self.mol)
    
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
            modres.GenerateModifiedProteinXYZAndKey(self,knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check)
    
    
        self.WriteToLog('Poltype Job Finished'+'\n')
    
    
    def GenerateParameters(self,mol):
        if self.use_psi4==True:
            try:
                import psi4
                self.WriteToLog("Psi4 detected; it will be used as the QM engine")
                self.use_psi4 = True
            except:
                self.WriteToLog("Psi4 not detected; attempting to use Gaussian instead")
                self.use_psi4 = False
    
    
        if self.use_psi4:
            psi4.set_num_threads(int(self.numproc))
            psi4.set_memory(self.maxmem)
        self.WriteToLog("Running on host: " + gethostname())
    
        # QM calculations are done here
        # First the molecule is optimized. (-opt) 
        # This optimized molecule is stored in the structure optmol
        # Then the electron density matrix is found (-dma)
        # This is used by GDMA to find multipoles
        # Then information for generating the electrostatic potential grid is found (-esp)
        # This information is used by cubegen
        optmol = opt.GeometryOptimization(self,mol)
        opt.CheckBondConnectivity(self,mol,optmol)
        esp.SPForDMA(self,optmol,mol)
        esp.SPForESP(self,optmol,mol) 
            
        # End here if qm calculations were all that needed to be done 
        if self.qmonly:
            self.WriteToLog("poltype QM-only complete.")
            sys.exit(0)
    
        # Initializing arrays
        
        self.canonicallabel = [ 0 ] * mol.NumAtoms()
        self.localframe1 = [ 0 ] * mol.NumAtoms()
        self.localframe2 = [ 0 ] * mol.NumAtoms()
    
        symm.gen_canonicallabels(self,mol)
       
        # scaling of multipole values for certain atom types
        # checks if the molecule contains any atoms that should have their multipole values scaled
        scalelist = mpole.process_types(self,mol)
        
            
        if(self.fitrotbnds):
            torgen.gen_fitrotbnds_list(self)
        # Find rotatable bonds for future torsion scans
        (self.torlist, self.rotbndlist) = torgen.get_torlist(self,mol)
        self.torlist = torgen.get_torlist_opt_angle(self,optmol, self.torlist)
           
    
        self.rotbndtoanginc=torgen.DetermineAngleIncrementForEachTorsion(self,mol,self.rotbndlist)
    
        # Obtain multipoles from Gaussian fchk file using GDMA
    
        if not os.path.isfile(self.gdmafname):
            mpole.run_gdma(self)
    
        # Set up input file for poledit
        # find multipole local frame definitions 
        lfzerox,atomindextoremovedipquad,atomtypetospecialtrace,atomindextoremovedipquadcross = mpole.gen_peditinfile(self,mol)
        
        
        if (not os.path.isfile(self.xyzfname) or not os.path.isfile(self.keyfname)):
            # Run poledit
            cmdstr = self.peditexe + " 1 " + self.gdmafname +' '+self.paramhead+ " < " + self.peditinfile
            self.call_subsystem(cmdstr,True)
            # Add header to the key file output by poledit
            mpole.prepend_keyfile(self,self.keyfname)
        # post process local frames written out by poledit
        mpole.post_proc_localframes(self,self.keyfname, lfzerox,atomindextoremovedipquad,atomindextoremovedipquadcross)
        # generate the electrostatic potential grid used for multipole fitting
        esp.gen_esp_grid(self,optmol)
    
        # Average multipoles based on molecular symmetry
        # Does this using the script avgmpoles.pl which is found in the poltype directory
        # Atoms that belong to the same symm class will now have only one common multipole definition
        if self.uniqidx:
            shutil.copy(self.keyfname, self.key2fname)
            mpole.prepend_keyfile(self,self.key2fname)
        elif ((not os.path.isfile(self.xyzoutfile) or not os.path.isfile(self.key2fname)) and not self.uniqidx):
            mpole.AverageMultipoles(self)
            mpole.prepend_keyfile(self,self.key2fname)    
        if self.espfit and not os.path.isfile(self.key3fname):
            # Optimize multipole parameters to QM ESP Grid (*.cube_2)
            # tinker's potential utility is called, with option 6.
            # option 6 reads: 'Fit Electrostatic Parameters to a Target Grid'
            
            esp.ElectrostaticPotentialFitting(self) 
        else:
            shutil.copy(self.key2fname, self.key3fname)
        # Remove header terms from the keyfile
        mpole.rm_esp_terms_keyfile(self,self.key3fname)
        esp.ElectrostaticPotentialComparison(self)  
        if not os.path.isfile(self.key4fname):
            shutil.copy(self.key3fname, self.key4fname)
            
            # Multipoles are scaled if needed using the scale found in process_types
            mpole.post_process_mpoles(self,self.key4fname, scalelist)
            
            # Now that multipoles have been found
            # Other parameters such as opbend, vdw, etc. are found here using a look up table
            # Part of the look up table is here in poltype.py 
            # Most of it is in the file valence.py found in the poltype directory
    
            # Finds aromatic carbons and associated hydrogens and corrects polarizability
            # Find opbend values using a look up table
            # Outputs a list of rotatable bonds (found in get_torlist) in a form usable by valence.py
    
    
            # Map from idx to symm class is made for valence.py
            idxtoclass=[]
            for i in range(mol.NumAtoms()):
                idxtoclass.append(symm.get_class_number(self,i+1))
            v = valence.Valence(self.output_format,self.logfname)
            v.setidxtoclass(idxtoclass)
            dorot = True
            # valence.py method is called to find parameters and append them to the keyfile
            v.appendtofile(self.key4fname, optmol, dorot,self.rotbndlist) 
    
               
        # Torsion scanning then fitting. *.key_5 will contain updated torsions
        # default, parmtors = True
        if self.dontdotor==True:
            sys.exit()
        if (self.parmtors):
            # torsion scanning
            torgen.gen_torsion(self,optmol,self.torsionrestraint)
        if (self.parmtors):
            # torsion fitting
            if self.dontdotorfit==True:
                sys.exit()
            torfit.process_rot_bond_tors(self,optmol)
        else:
            shutil.copy(self.key4fname,self.key5fname)
        
        # A series of tests are done so you one can see whether or not the parameterization values
        # found are acceptable and to what degree
        opt.StructureMinimization(self)
        opt.gen_superposeinfile(self)
        opt.CheckRMSD(self)
        self.WriteToLog("")
        self.WriteToLog("=========================================================")
        self.WriteToLog("QM Dipole moment\n")
        grepcmd = 'grep -A7 "Dipole moment" ' + self.logespfname
        self.call_subsystem(grepcmd,True)
    
        self.WriteToLog("")
        self.WriteToLog("=========================================================")
        self.WriteToLog("MM Dipole moment\n")
        cmd=self.analyzeexe + ' ' +  self.xyzoutfile + ' em | grep -A11 Charge'
        self.call_subsystem(cmd,True)
obj=PolarizableTyper() 

