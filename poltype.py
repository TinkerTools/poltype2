
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
from PoltypeModules import fragmenter as frag
from parmed.tinker import parameterfile
from rdkit.Chem import rdmolfiles


class PolarizableTyper():

    def __init__(self,poltypepath=os.path.split(sys.argv[0])[0],WBOtol=.01,wholexyz=None,wholemol=None,dontfrag=False,isfragjob=False,dipoletol=.1,externalapi=None,printoutput=False,poltypeini=True,structure=None,prmstartidx=401,numproc=1,maxmem="700MB",maxdisk="100GB",gausdir=None,gdmadir=None,tinkerdir=None,scratchdir="/scratch",paramhead=sys.path[0] + "/amoebabio18_header.prm",babelexe="babel",gausexe='g09',formchkexe='formchk',cubegenexe='cubegen',gdmaexe='gdma',avgmpolesexe=sys.path[0] + "/avgmpoles.pl",peditexe='poledit.x',potentialexe='potential.x',minimizeexe='minimize.x',analyzeexe='analyze.x',superposeexe='superpose.x',defopbendval=0.20016677990819662,Hartree2kcal_mol=627.5095,optbasisset='6-31G*',toroptbasisset='6-31G*',dmabasisset='6-311G**',espbasisset="6-311++G(2d,2p)",torspbasisset="6-311++G**",optmethod='wB97X-D',toroptmethod='wB97X-D',torspmethod='MP2',dmamethod='MP2',espmethod='MP2',qmonly = False,espfit = True,parmtors = True,foldnum=3,foldoffsetlist = [ 0.0, 180.0, 0.0, 0.0, 0.0, 0.0 ],torlist = None,rotbndlist = None,fitrotbndslist=None,maxRMSD=.1,maxRMSPD=1,maxtorRMSPD=2,tordatapointsnum=None,gentorsion=False,gaustorerror=False,torsionrestraint=.1,onlyrotbndlist=None,rotalltors=False,dontdotor=False,dontdotorfit=False,toroptpcm=False,optpcm=False,torsppcm=False,use_gaus=False,use_gausoptonly=False,freq=False,postfit=False,bashrcpath=None,amoebabioprmpath=None,libpath=sys.path[0] + "/lib.bio18_conv1.txt",SMARTSToTypelibpath=sys.path[0]+'/SMARTSToTypeLib.txt',ModifiedResiduePrmPath=sys.path[0]+'/ModifiedResidue.prm',modifiedproteinpdbname=None,unmodifiedproteinpdbname=None,mutatedsidechain=None,mutatedresiduenumber=None,modifiedresiduepdbcode=None,optmaxcycle=400,torkeyfname=None,gausoptcoords='',uniqidx=False ,helpfile=sys.path[0]+'/README.HELP',versionfile=sys.path[0]+'/README.VERSION'): 
        self.WBOtol=WBOtol
        self.isfragjob=isfragjob
        self.wholexyz=wholexyz
        self.wholemol=wholemol
        self.dontfrag=dontfrag
        self.dipoletol=dipoletol
        self.externalapi=externalapi
        self.printoutput=printoutput
        self.poltypepath=poltypepath
        self.poltypeparentdir=os.path.split(self.poltypepath)[0] 
        self.molstructfname=structure
        self.poltypeini=poltypeini
        self.prmstartidx = prmstartidx
        self.numproc = numproc
        self.maxmem = maxmem
        self.maxdisk = maxdisk
        self.gausdir = gausdir
        self.gdmadir = gdmadir
        self.tinkerdir = tinkerdir
        self.scratchdir = scratchdir
        self.paramhead = paramhead
        self.babelexe = babelexe
        self.gausexe =  gausexe
        self.formchkexe =  formchkexe
        self.cubegenexe =  cubegenexe
        self.gdmaexe = gdmaexe
        self.avgmpolesexe = avgmpolesexe
        self.peditexe = peditexe
        self.potentialexe = potentialexe
        self.minimizeexe = minimizeexe
        self.analyzeexe = analyzeexe
        self.superposeexe = superposeexe
        self.defopbendval = defopbendval
        self.Hartree2kcal_mol = Hartree2kcal_mol
        self.optbasisset = optbasisset             
        self.toroptbasisset = toroptbasisset         
        self.dmabasisset = dmabasisset             
        self.espbasisset = espbasisset         
        self.torspbasisset = torspbasisset
        self.optmethod=optmethod                     
        self.toroptmethod=toroptmethod                  
        self.torspmethod=torspmethod                    
        self.dmamethod=dmamethod                      
        self.espmethod=espmethod                  
        self.qmonly = qmonly
        self.espfit = espfit
        self.parmtors = parmtors
        self.foldnum=foldnum
        self.nfoldlist =  list(range(1,self.foldnum+1))
        self.foldoffsetlist = foldoffsetlist
        if torlist==None:
            self.torlist = []
        else:
            self.torlist=torlist
        if rotbndlist==None:
            self.rotbndlist = []
        else:
            self.rotbndlist=rotbndlist
        if fitrotbndslist==None:
            self.fitrotbndslist=[]
        else:
           
            self.fitrotbndslist=fitrotbndslist.split(',')
            templist=[]
            for ele in self.fitrotbndslist:
                nums=ele.lstrip().rstrip().split()
                temp=[]
                for e in nums:
                    temp.append(int(e))
                templist.append(temp)
            self.fitrotbndslist=templist

        self.maxRMSD=maxRMSD
        self.maxRMSPD=maxRMSPD
        self.maxtorRMSPD=maxtorRMSPD
        self.tordatapointsnum=tordatapointsnum
        self.gentorsion=gentorsion
        self.gaustorerror=gaustorerror
        self.torsionrestraint=torsionrestraint
        if onlyrotbndlist==None:
            self.onlyrotbndlist=[]
        else:
            self.onlyrotbndlist=[i.lstrip().rstrip() for i in onlyrotbndlist.split(',')]
        self.rotalltors=rotalltors
        self.dontdotor=dontdotor
        self.dontdotorfit=dontdotorfit
        self.toroptpcm=toroptpcm
        self.optpcm=optpcm
        self.torsppcm=torsppcm
        self.use_gaus=use_gaus
        self.use_gausoptonly=use_gausoptonly
        self.freq=freq
        self.postfit=postfit
        self.bashrcpath=bashrcpath
        self.amoebabioprmpath=amoebabioprmpath
        self.libpath=libpath
        self.SMARTSToTypelibpath=SMARTSToTypelibpath
        self.ModifiedResiduePrmPath=ModifiedResiduePrmPath
        self.modifiedproteinpdbname=modifiedproteinpdbname
        self.unmodifiedproteinpdbname=unmodifiedproteinpdbname
        self.mutatedsidechain=mutatedsidechain
        self.mutatedresiduenumber=mutatedresiduenumber
        self.modifiedresiduepdbcode=modifiedresiduepdbcode
        self.optmaxcycle=optmaxcycle
        self.torkeyfname=torkeyfname
        self.gausoptcoords=gausoptcoords
        self.uniqidx=uniqidx
        self.helpfile=helpfile
        self.versionfile=versionfile 
        opts, xargs = getopt.getopt(sys.argv[1:],'h:u',["help","unittest"])


        for o, a in opts:
            if o in ("-h", "--help"):
                self.copyright()
                self.usage()
                sys.exit(2)
            elif o in("-u","--unittest"):
                cmdstr='python'+' '+self.poltypeparentdir+'/'+'PoltypeModules/'+'test_poltype.py'
                os.system(cmdstr)
                sys.exit()
        
        if self.poltypeini==True:
            temp=open(os.getcwd()+r'/'+'poltype.ini','r')
            results=temp.readlines()
            temp.close()
            for line in results:
                if '#' not in line and line!='\n':
                    if '=' in line:
                        linesplit=line.split('=',1)
                        a=linesplit[1].replace('\n','').rstrip().lstrip()
                        newline=linesplit[0]
                    else:
                        newline=line

                    if "rotalltors" in newline:
                        self.rotalltors = True
                    elif "dontfrag" in newline:
                        self.dontfrag=True
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
                    elif "optpcm" in newline and 'tor' not in line:
                        self.optpcm = True
                    elif "toroptpcm" in newline:
                        self.optpcm = True
                    elif "use_gaus" in newline and 'opt' not in newline:
                        self.use_gaus = True
                    elif "use_gausoptonly" in newline:
                        self.use_gausoptonly = True
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
                    elif "optmethod" in newline and 'tor' not in newline:
                        self.optmethod = a
                    elif "espmethod" in newline and 'tor' not in newline:
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
                    elif "optbasisset" in newline and 'tor' not in newline:
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
    
        
        if not __name__ == '__main__':
            params=self.main()
 

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
    
        if self.use_gaus:
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
    
    
                
    
    
        cmdstr=self.analyzeexe+' '+self.poltypepath+r'/'+'water.xyz'+' '+'-k'+' '+self.poltypepath+r'/'+'water.key'+' '+'e'+'>'+' '+'version.out'
        try:
            if self.printoutput==True:
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
                self.versionnum=float(linesplit[2])
                if self.versionnum>=8.7:
                    latestversion = True
                    break
           
        if(not latestversion):
            
            raise ValueError("Notice: Not latest version of tinker (>=8.7)")
      
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
    
        if self.use_gaus or self.use_gausoptonly:
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
        self.key6fname= self.assign_filenames ( "key6fname" , ".key_6")
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
        if self.printoutput==True:
            print("Calling: " + cmdstr)
        self.WriteToLog(" Calling: " + cmdstr)
        p = subprocess.Popen(cmdstr, shell=True,stdout=self.logfh, stderr=self.logfh)
        
        if p.wait() != 0:
            self.WriteToLog("ERROR: " + cmdstr)
            if wait==True:
                sys.exit(1)

    def WriteOutLiteratureReferences(self,keyfilename): # to use ParmEd on key file need Literature References delimited for parsing
        temp=open(keyfilename,'r')
        results=temp.readlines()
        temp.close()
        tempname=keyfilename.replace('.key','_temp.key')
        temp=open(tempname,'w')
        foundatomblock=False
        for i in range(len(results)):
            line=results[i]
            if 'atom' in line and foundatomblock==False:
                foundatomblock=True
                temp.write('#############################'+'\n')
                temp.write('##                         ##'+'\n')
                temp.write('##  Literature References  ##'+'\n')
                temp.write('##                         ##'+'\n')
                temp.write('#############################'+'\n')
                temp.write('\n')
                temp.write('Wu, J.C.; Chattree, G.; Ren, P.Y.; Automation of AMOEBA polarizable force field'+'\n')
                temp.write('parameterization for small molecules. Theor Chem Acc.'+'\n')
                temp.write('\n')
                temp.write(line)
            else:
                temp.write(line)
        os.remove(keyfilename)
        os.replace(tempname,keyfilename)    

    def RaiseOutputFileError(self,logname):
        raise ValueError('An error occured for '+logname) 
     
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
                self.GenerateParameters()
        else:
           params= self.GenerateParameters()
           self.WriteToLog('Poltype Job Finished'+'\n')
           return params

    
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
            modres.GenerateModifiedProteinXYZAndKey(self,knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check)
    
    
            
    def GenerateParameters(self):
        obConversion = openbabel.OBConversion()
        mol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(self.molstructfname)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(mol, self.molstructfname)
        obConversion.SetOutFormat('mol')
        self.molstructfnamemol=self.molstructfname.replace('.sdf','.mol')
        obConversion.WriteFile(mol,self.molstructfnamemol)
        self.mol=mol
    
        # Begin log. *-poltype.log
        self.logfh = open(self.logfname,"a")
    
        self.mol=self.CheckIsInput2D(self.mol,obConversion)
    
        self.totalcharge=self.mol.GetTotalCharge()
        if self.totalcharge!=0:
            self.toroptpcm=True
            self.optpcm=True
            self.torsppcm=True

     
        self.WriteToLog("Running on host: " + gethostname())
        # Initializing arrays
        
        self.canonicallabel = [ 0 ] * mol.NumAtoms()
        self.localframe1 = [ 0 ] * mol.NumAtoms()
        self.localframe2 = [ 0 ] * mol.NumAtoms()
    
        symm.gen_canonicallabels(self,mol)
 
        # QM calculations are done here
        # First the molecule is optimized. (-opt) 
        # This optimized molecule is stored in the structure optmol
        # Then the electron density matrix is found (-dma)
        # This is used by GDMA to find multipoles
        # Then information for generating the electrostatic potential grid is found (-esp)
        # This information is used by cubegen
        optmol = opt.GeometryOptimization(self,mol)
        esp.SPForDMA(self,optmol,mol)
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
            mpole.prepend_keyfile(self,self.keyfname,optmol)
        # post process local frames written out by poledit
        mpole.post_proc_localframes(self,self.keyfname, lfzerox,atomindextoremovedipquad,atomindextoremovedipquadcross)

        esp.SPForESP(self,optmol,mol) 
            
        # End here if qm calculations were all that needed to be done 
        if self.qmonly:
            self.WriteToLog("poltype QM-only complete.")
            sys.exit(0)
    
               
        # scaling of multipole values for certain atom types
        # checks if the molecule contains any atoms that should have their multipole values scaled
        scalelist = mpole.process_types(self,mol)
        
        
        # generate the electrostatic potential grid used for multipole fitting
        esp.gen_esp_grid(self,optmol)
    
        # Average multipoles based on molecular symmetry
        # Does this using the script avgmpoles.pl which is found in the poltype directory
        # Atoms that belong to the same symm class will now have only one common multipole definition
        if self.uniqidx:
            shutil.copy(self.keyfname, self.key2fname)
            mpole.prepend_keyfile(self,self.key2fname,optmol,True)
        elif ((not os.path.isfile(self.xyzoutfile) or not os.path.isfile(self.key2fname)) and not self.uniqidx):
            mpole.AverageMultipoles(self,optmol)
            mpole.prepend_keyfile(self,self.key2fname,optmol,True)    
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
        # Find rotatable bonds for future torsion scans
        (self.torlist, self.rotbndlist) = torgen.get_torlist(self,mol)
        self.torlist = torgen.get_torlist_opt_angle(self,optmol, self.torlist)
        self.torlist=torgen.RemoveDuplicateRotatableBondTypes(self) # this only happens in very symmetrical molecules
           
        torlist=[i[:4] for i in self.torlist]
    
        self.rotbndtoanginc=torgen.DetermineAngleIncrementForEachTorsion(self,mol,self.rotbndlist)
 
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
            v = valence.Valence(self.versionnum,self.logfname)
            v.setidxtoclass(idxtoclass)
            dorot = True
            # valence.py method is called to find parameters and append them to the keyfile
            v.appendtofile(self.key4fname, optmol, dorot,self.rotbndlist)
            if self.totalcharge!=0:
                torgen.PrependStringToKeyfile(self,self.key4fname,'solvate GK')
        if self.isfragjob==False and not os.path.isfile(self.key5fname) and self.dontfrag==False:

            self.rdkitmol=rdmolfiles.MolFromMolFile(poltype.molstructfnamemol,sanitize=True,removeHs=False)
            WBOmatrix,outputname=frag.GenerateWBOMatrix(poltype,self.rdkitmol,self.logoptfname.replace('.log','.xyz'))
            highlightbonds=[]
            for tor in self.torlist:
                rotbnd=[tor[1]-1,tor[2]-1]
                highlightbonds.append(rotbnd)
            frag.Draw2DMoleculeWithWBO(self,WBOmatrix,self.molstructfname.replace('.sdf',''),self.molstructfnamemol,bondindexlist=highlightbonds)        
            rotbndindextofragindexmap,rotbndindextofragment=frag.GenerateFragments(self,self.mol,torlist,WBOmatrix) # returns list of bond indexes that need parent molecule to do torsion scan for (fragment generated was same as the parent0
            frag.SpawnPoltypeJobsForFragments(self,rotbndindextofragindexmap,rotbndindextofragment)

        if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key5fname):
            frag.GrabTorsionParametersFromFragments(poltype,torlist) # just dump to key_5 since does not exist for parent molecule
            sys.exit()

       
               
        # Torsion scanning then fitting. *.key_5 will contain updated torsions
        if os.path.isfile(self.key5fname):
            self.parmtors=False
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
        self.WriteOutLiteratureReferences(self.key5fname) 
        # A series of tests are done so you one can see whether or not the parameterization values
        # found are acceptable and to what degree
        opt.StructureMinimization(self)
        opt.gen_superposeinfile(self)
        opt.CheckRMSD(self)
        if self.totalcharge!=0:
            torgen.RemoveStringFromKeyfile(self,self.key5fname,'solvate GK')
        esp.CheckDipoleMoments(self,optmol)
        keyfilecopyname=self.key5fname.replace('.key','_copy.key')
        shutil.copy(self.key5fname,keyfilecopyname)
        torgen.RemoveStringFromKeyfile(self,keyfilecopyname,'SOLUTE')
        torgen.RemoveStringFromKeyfile(self,keyfilecopyname,'TARGET-DIPOLE')
        os.system('rm *.chk')
        os.remove(self.tmpxyzfile+'_2') 
        param = parameterfile.AmoebaParameterSet(keyfilecopyname)
        if self.dontfrag==False and self.isfragjob==True and not os.path.isfile(self.key6fname):
            wholemolidxtofragidx=frag.ConvertFragIdxToWholeIdx(self,self.molstructfname,torlist,rotbndindextofragindexmap,WBOmatrix)
        return param

if __name__ == '__main__':
    poltype=PolarizableTyper() 
    params=poltype.main()
