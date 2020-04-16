
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
import valence
import torsiongenerator as torgen
import modifiedresidues as modres
import symmetry as symm
import torsionfit as torfit
import optimization as opt
import electrostaticpotential as esp
import multipole as mpole
import fragmenter as frag
from parmed.tinker import parameterfile
from rdkit import Chem
from rdkit.Chem import rdmolfiles


class PolarizableTyper():

    def __init__(self,readinionly=False,suppressdipoleerr=False,topologylib='residue_connect.txt',poltypepath=os.path.split(__file__)[0],WBOtol=.01,dontfrag=True,isfragjob=False,dipoletol=.1,externalapi=None,printoutput=False,poltypeini=True,structure=None,prmstartidx=401,numproc=1,maxmem="700MB",maxdisk="100GB",gausdir=None,gdmadir=None,tinkerdir=None,scratchdir="/scratch",paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/amoebabio18_header.prm",babelexe="babel",gausexe='g09',formchkexe='formchk',cubegenexe='cubegen',gdmaexe='gdma',avgmpolesexe=os.path.abspath(os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), os.pardir)) + "/avgmpoles.pl",peditexe='poledit',potentialexe='potential',minimizeexe='minimize',analyzeexe='analyze',superposeexe='superpose',defopbendval=0.20016677990819662,Hartree2kcal_mol=627.5095,optbasisset='6-31G*',toroptbasisset='6-31G*',dmabasisset='6-311G**',espbasisset="6-311++G(2d,2p)",torspbasisset="6-311++G**",optmethod='wB97X-D',toroptmethod='wB97X-D',torspmethod='MP2',dmamethod='MP2',espmethod='MP2',qmonly = False,espfit = True,parmtors = True,foldnum=3,foldoffsetlist = [ 0.0, 180.0, 0.0, 0.0, 0.0, 0.0 ],torlist = None,rotbndlist = None,fitrotbndslist=None,maxRMSD=.1,maxRMSPD=1,maxtorRMSPD=2,tordatapointsnum=None,gentorsion=False,gaustorerror=False,torsionrestraint=.1,onlyrotbndlist=None,rotalltors=False,dontdotor=False,dontdotorfit=False,toroptpcm=False,optpcm=False,torsppcm=False,use_gaus=False,use_gausoptonly=False,freq=False,postfit=False,bashrcpath=None,amoebabioprmpath=None,libpath=sys.path[0] + "/lib.bio18_conv1.txt",SMARTSToTypelibpath=sys.path[0]+'/SMARTSToTypeLib.txt',ModifiedResiduePrmPath=sys.path[0]+'/ModifiedResidue.prm',modifiedproteinpdbname=None,unmodifiedproteinpdbname=None,mutatedsidechain=None,mutatedresiduenumber=None,modifiedresiduepdbcode=None,optmaxcycle=400,torkeyfname=None,gausoptcoords='',helpfile='README.HELP',versionfile='README.VERSION'): 
        self.readinionly=readinionly
        self.suppressdipoleerr=suppressdipoleerr
        self.use_gaus=use_gaus
        self.use_gausoptonly=use_gausoptonly
        self.topologylibpath=os.path.abspath(os.path.join(poltypepath, '..'))+r'/'+topologylib
        self.WBOtol=WBOtol
        self.isfragjob=isfragjob
        self.dontfrag=dontfrag
        self.dipoletol=dipoletol
        self.externalapi=externalapi
        self.printoutput=printoutput
        self.poltypepath=poltypepath
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
        self.optmethod=self.SanitizeQMMethod(optmethod,True)                 
        self.toroptmethod=self.SanitizeQMMethod(toroptmethod,True)                  
        self.torspmethod=self.SanitizeQMMethod(torspmethod,False)                    
        self.dmamethod=self.SanitizeQMMethod(dmamethod,False)                      
        self.espmethod=self.SanitizeQMMethod(espmethod,False)                  
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
        self.helpfile=helpfile
        self.versionfile=versionfile 
        opts, xargs = getopt.getopt(sys.argv[1:],'h',["help"])

        for o, a in opts:
            if o in ("-h", "--help"):
                self.copyright()
                self.usage()
                sys.exit(2)
                            
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
                        if '=' not in line:
                            self.rotalltors = True
                        else:
                            self.rotalltors=self.GrabBoolValue(a)
                    elif 'poltypepath' in newline:
                        self.poltypepath=a
                    elif 'printoutput' in newline:
                        if '=' not in line:
                            self.printoutput=True
                        else:
                            self.printoutput=self.GrabBoolValue(a)
                    elif 'suppressdipoleerr' in newline:
                        if '=' not in line:
                            self.suppressdipoleerr=True
                        else:
                            self.suppressdipoleerr=self.GrabBoolValue(a)
                    elif 'isfragjob' in newline:
                        if '=' not in line:
                            self.isfragjob=True
                        else:
                            self.isfragjob=self.GrabBoolValue(a)
                    elif "dontfrag" in newline:
                        if '=' not in line:
                            self.dontfrag=True
                        else:
                            self.dontfrag=self.GrabBoolValue(a)
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
                        self.dmamethod =a
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
                        if '=' not in line:
                            self.torsppcm = True
                        else: 
                            self.torsppcm=self.GrabBoolValue(a)
                    elif "freq" in newline:
                        if '=' not in line:
                            self.freq = True
                        else:
                            self.freq=self.GrabBoolValue(a)
                    elif "optpcm" in newline and 'tor' not in line:
                        if '=' not in line:
                            self.optpcm = True
                        else:
                            self.optpcm=self.GrabBoolValue(a)
                    elif "toroptpcm" in newline:
                        if '=' not in line:
                            self.optpcm = True
                        else:
                            self.optpcm=self.GrabBoolValue(a)
                    elif "use_gaus" in newline and 'opt' not in newline:
                        if '=' not in line:
                            self.use_gaus = True
                        else:
                            self.use_gaus=self.GrabBoolValue(a)
                    elif "use_gausoptonly" in newline:
                        if '=' not in line:
                            self.use_gausoptonly = True
                        else:
                            self.use_gausoptonly=self.GrabBoolValue(a)
                    elif "dontdotor" in newline:
                        if '=' not in line:
                            self.dontdotor = True
                        else:
                            self.dontdotor=self.GrabBoolValue(a)

                    elif "dontdotorfit" in newline:
                        if '=' not in line:
                            self.dontdotorfit = True
                        else:
                            self.dontdotorfit=self.GrabBoolValue(a)
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
                        if '=' not in line:
                            self.qmonly = True
                        else:
                            self.qmonly = self.GrabBoolValue(a)
                    elif "help" in newline:
                        self.copyright()
                        self.usage()
                        sys.exit(2)
                    elif "onlyrotbnd" in newline: # comma seperated list 
                        self.onlyrotbndlist=[i.lstrip().rstrip() for i in a.split(',')]
                    else:
                        print('Unrecognized '+line)
                        self.usage()
                        print('Unrecognized '+line)
                        sys.exit()
        self.optmethod=self.SanitizeQMMethod(optmethod,True)                 
        self.toroptmethod=self.SanitizeQMMethod(toroptmethod,True)                  
        self.torspmethod=self.SanitizeQMMethod(torspmethod,False)                    
        self.dmamethod=self.SanitizeQMMethod(dmamethod,False)                      
        self.espmethod=self.SanitizeQMMethod(espmethod,False)                  
        if self.readinionly==True:
            return
        self.SanitizeMMExecutables()
        self.copyright()
        self.initialize()
        self.init_filenames()
 
        # Use openbabel to create a 'mol' object from the input molecular structure file. 
        # Openbabel does not play well with certain molecular structure input files,
        # such as tinker xyz files. (normal xyz are fine)
    
        
        if not __name__ == '__main__':
            params=self.main()

    def GrabBoolValue(self,value): 
        if value=='true' or value=='True' or value=='TRUE':
            boolvalue=True
        elif value=='false' or value=='False' or value=='FALSE':
            boolvalue=False
        return boolvalue

    def SanitizeQMMethod(self,method,optmethodbool):
        if method[-1]=='D': # assume DFT, gaussian likes D, PSI4 likes -D
            if method[-2]=='-':
                if self.use_gaus or self.use_gausoptonly and optmethodbool==True:
                    method=method.replace('-D','D')
                if self.use_gaus and optmethodbool==False:
                    method=method.replace('-D','D')

            else:
                if not self.use_gaus and not self.use_gausoptonly and optmethodbool==True:
                    method=method.replace('D','-D')
                if not (self.use_gaus) and optmethodbool==False:
                    method=method.replace('D','-D')

        return method
 

    def WriteToLog(self,string):
        now = time.strftime("%c",time.localtime())
        self.logfh.write(now+' '+string+'\n')
        self.logfh.flush()
        os.fsync(self.logfh.fileno())

        
    def SanitizeMMExecutables(self):
        path=self.which(self.peditexe)
        if path==None:
            self.peditexe='poledit.x'
            self.potentialexe='potential.x'
            self.minimizeexe='minimize.x'
            self.analyzeexe='analyze.x'
            self.superposeexe='superpose.x'

    def which(self,program):
        def is_exe(fpath):
            try:
                 return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
            except:
                 return None
    
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
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
    
    
                
    
        cmdstr=self.analyzeexe+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/'+'water.xyz'+' '+'-k'+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/'+'water.key'+' '+'e'+'>'+' '+'version.out'
        try:
            print('Calling: '+cmdstr) 
            returned_value = subprocess.call(cmdstr, shell=True)
        except:
            raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())      
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
            
            raise ValueError("Notice: Not latest working version of tinker (8.7)"+' '+os.getcwd())
      
        if ("TINKERDIR" in os.environ):
            self.tinkerdir = os.environ["TINKERDIR"]
            self.peditexe = os.path.join(self.tinkerdir,self.peditexe)
            self.potentialexe = os.path.join(self.tinkerdir,self.potentialexe)
            self.minimizeexe = os.path.join(self.tinkerdir,self.minimizeexe)
            self.analyzeexe = os.path.join(self.tinkerdir,self.analyzeexe)
            self.superposeexe = os.path.join(self.tinkerdir,self.superposeexe)
    
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
        temp=open(os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/'+self.versionfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            print(line)
    
    def usage (self):
        temp=open(os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/'+self.helpfile,'r')
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
    
    def CallJobsSeriallyLocalHost(self,fulljobtooutputlog,skiperrors=False):
       for job in fulljobtooutputlog.keys():
           self.call_subsystem(job,True,skiperrors)
       finishedjobs,errorjobs=self.WaitForTermination(fulljobtooutputlog)
       return finishedjobs,errorjobs

    def CallJobsLocalHost(self,jobtooutputlog):
        for job in jobtooutputlog.keys():
           self.call_subsystem(job)
        finishedjobs,errorjobs=self.WaitForTermination(jobtooutputlog)
        return finishedjobs,errorjobs


    def WaitForTermination(self,jobtooutputlog):
        finishedjobs=[]
        errorjobs=[]
        sleeptime=1
        errormessages=[]
        while len(finishedjobs)!=len(jobtooutputlog.keys()):
            for job in jobtooutputlog.keys():
                outputlog=jobtooutputlog[job]
                if os.path.isfile(outputlog):
                    statinfo=os.stat(outputlog)
                    size=statinfo.st_size
                    if size==0:
                        continue
                finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages)
                if finished==True and error==False: # then check if SP has been submitted or not
                    if outputlog not in finishedjobs:
                        self.NormalTerm(outputlog)
                        finishedjobs.append(outputlog)
                elif finished==False and error==True:
                    if outputlog not in finishedjobs:
                        self.ErrorTerm(outputlog)
                        finishedjobs.append(outputlog)
                        errorjobs.append(outputlog)
                elif finished==False and error==False:
                    if not os.path.isfile(outputlog):
                        self.WriteToLog('Waiting on '+outputlog+' '+'to begin')
                    else:
                        self.WriteToLog('Waiting on '+outputlog+' '+'for termination ')
                else: # this case is finshed=True and error=True because there stupid quotes sometimes have word error in it                  
                    if outputlog not in finishedjobs:
                        error=False
                        self.NormalTerm(outputlog)
                        finishedjobs.append(outputlog)
    
            #string='Sleeping for %d '%(sleeptime)+' minute '
            #self.WriteToLog(string)
            time.sleep(sleeptime*60) # check logs every minute
            self.WriteToLog('**********************************************')
        self.WriteToLog('All jobs have terminated ')
        return finishedjobs,errorjobs

    def CheckNormalTermination(self,logfname,errormessages=None): # needs to handle error checking now
        """
        Intent: Checks the *.log file for normal termination
        """
        error=False
        term=False
        if os.path.isfile(logfname):
            head,tail=os.path.split(logfname)
            for line in open(logfname):
                if 'poltype' in tail:
                    if 'Poltype Job Finished' in line:
                        term=True
                else:
                    if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line or "LBFGS  --  Normal Termination due to SmallGrad" in line or "Normal termination" in line:
                        term=True
                if ('error' in line or 'Error' in line or 'ERROR' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation, address not mapped to object' in line or 'galloc:  could not allocate memory' in line or 'Erroneous write.' in line) and 'DIIS' not in line and 'mpi' not in line:
                    if 'segmentation violation' in line and 'address not mapped to object' not in line:
                        error=False
                        continue
                    error=True
                if error==True:
                    message='Error '+line+ 'logpath='+logfname
                    
            if error==True and term==False:
                if errormessages!=None:
                    if message not in errormessages:
                        self.WriteToLog(message) 
                        errormessages.append(message)
                else:
                    self.WriteToLog(message) 


        if errormessages!=None:
            return term,error,errormessages
        else:
            return term,error
       
    
        
    def NormalTerm(self,logfname):
        self.WriteToLog("Normal termination: logfile=%s path=%s" % (logfname,os.getcwd()))
    
    
    def ErrorTerm(self,logfname):
        self.WriteToLog("ERROR termination: logfile=%s path=%s" % (logfname,os.getcwd()))


    def call_subsystem(self,cmdstr,wait=False,skiperrors=False):
        curdir=os.getcwd()
        error=False
        if self.printoutput==True:
            print("Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
        self.WriteToLog(" Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
        p = subprocess.Popen(cmdstr, shell=True,stdout=self.logfh, stderr=self.logfh)

        if wait==True:
            p.wait()
            if p.returncode != 0 and skiperrors==False:
                self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
        return error

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

    def WritePoltypeInitializationFile(self,poltypeinput):
        inifilepath=os.getcwd()+r'/'+'poltype.ini'
        temp=open(inifilepath,'w')
        for key,value in poltypeinput.items():
            line=key+'='+str(value)+'\n'
            temp.write(line)
        temp.close()
        return inifilepath

    #POLTYPE BEGINS HERE
    def main(self):
    
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
            knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference=modres.GenerateModifiedProteinPoltypeInput(self)
            self.molstructfname=molname
            head, self.molstructfname = os.path.split(self.molstructfname)
            self.molecprefix =  os.path.splitext(self.molstructfname)[0]

        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None): # if already have core parameters in modified prm database then dont regenerate parameters
            if check==False:
                self.GenerateParameters()
        else:
           params= self.GenerateParameters()
           return params

    
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
            modres.GenerateModifiedProteinXYZAndKey(self,knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference)
    
    
            
    def GenerateParameters(self):
        if os.path.isfile(self.tmpxyzfile+'_2'):
            os.remove(self.tmpxyzfile+'_2') 

        obConversion = openbabel.OBConversion()
        mol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(self.molstructfname)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(mol, self.molstructfname)
        obConversion.SetOutFormat('mol')
        self.molstructfnamemol=self.molstructfname.replace('.sdf','.mol')
        obConversion.WriteFile(mol,self.molstructfnamemol)
        self.mol=mol
        m=Chem.MolFromMolFile(self.molstructfnamemol,removeHs=False)
    
        # Begin log. *-poltype.log
        if os.path.isfile(self.logfname):
            os.remove(self.logfname)
        self.logfh = open(self.logfname,"a",buffering=1)
    
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
    
        symm.CalculateSymmetry(self,m)

 
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
        if os.path.isfile(self.key2fname):
            os.remove(self.key2fname)
        if os.path.isfile(self.xyzoutfile):
            os.remove(self.xyzoutfile)
        mpole.AverageMultipoles(self,optmol)
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
            v = valence.Valence(self.versionnum,self.logfname,self.dontfrag)
            v.setidxtoclass(self.idxtosymclass)
            dorot = True
            # valence.py method is called to find parameters and append them to the keyfile
            v.appendtofile(self.key4fname, optmol, dorot,self.rotbndlist)
            if self.totalcharge!=0:
                torgen.PrependStringToKeyfile(self,self.key4fname,'solvate GK')
        if self.isfragjob==False and not os.path.isfile(self.key5fname) and self.dontfrag==False:

            self.rdkitmol=rdmolfiles.MolFromMolFile(self.molstructfnamemol,sanitize=True,removeHs=False)
            WBOmatrix,outputname=frag.GenerateWBOMatrix(self,self.rdkitmol,self.logoptfname.replace('.log','.xyz'))
            highlightbonds=[]
            for tor in self.torlist:
                rotbnd=[tor[1]-1,tor[2]-1]
                highlightbonds.append(rotbnd)
            frag.Draw2DMoleculeWithWBO(self,WBOmatrix,self.molstructfname.replace('.sdf',''),self.molstructfnamemol,bondindexlist=highlightbonds)        
            rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentfragmentsarray,equivalentrotbndindexarrays=frag.GenerateFragments(self,self.mol,torlist,WBOmatrix) # returns list of bond indexes that need parent molecule to do torsion scan for (fragment generated was same as the parent0
            frag.SpawnPoltypeJobsForFragments(self,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,torlist,equivalentfragmentsarray,equivalentrotbndindexarrays)

        if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key5fname):
            frag.GrabTorsionParametersFromFragments(self,torlist,rotbndindextofragmentfilepath) # just dump to key_5 since does not exist for parent molecule
            return

       
               
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
        self.WriteToLog('Poltype Job Finished'+'\n')
        keyfilecopyname=self.key5fname.replace('.key','_copy.key')
        shutil.copy(self.key5fname,keyfilecopyname)
        torgen.RemoveStringFromKeyfile(self,keyfilecopyname,'SOLUTE')
        torgen.RemoveStringFromKeyfile(self,keyfilecopyname,'TARGET-DIPOLE')
        for fname in os.listdir():
            if fname.endswith('.chk'):
                os.remove(fname)
        if os.path.isdir('qm-torsion'):
            os.chdir('qm-torsion')
            for fname in os.listdir():
                if fname.endswith('.chk'):
                    os.remove(fname)
            os.chdir('..')
        param = parameterfile.AmoebaParameterSet(keyfilecopyname)
        return param



if __name__ == '__main__':
    poltype=PolarizableTyper() 
    params=poltype.main()
