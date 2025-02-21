################################################################
#
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
import warnings
import math
import traceback
import os
import sys
from socket import gethostname
import subprocess
from openbabel import openbabel
import shutil
import time
import copy
import getopt
from collections import OrderedDict
import torsiondatabaseparser
import dimorphite_dl
import torsiongenerator as torgen
import symmetry as symm
import torsionfit as torfit
import optimization as opt
import electrostaticpotential as esp
import multipole as mpole
from distributed_multipole import get_dma_default
from parmmod import modify_key
import fragmenter as frag
import rings
from packaging import version
from rdkit import Chem
from rdkit.Chem import rdmolfiles,AllChem,rdmolops,Descriptors
from rdkit.Geometry import Point3D
import vdwfit
import numpy as np
import itertools
from rdkit.Chem import rdMolTransforms
import psutil
import multiprocessing
import time
import re
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdFMCS
from PyAstronomy import pyasl
from dataclasses import dataclass,field
from pathlib import Path
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdDistGeom
from itertools import product,combinations
import random
import ldatabaseparser
import lmodifytinkerkey

@dataclass
class PolarizableTyper():
        maxtorresnitrogen:int=2
        xtbtorresconstant:float=5
        torfit:bool=True
        toroptmethodlist:list=field(default_factory=lambda : [])
        torspmethodlist:list=field(default_factory=lambda : [])
        anienvname:str='ani'
        xtbenvname:str='xtbenv'
        anifmax:float=.05
        anipath:str=os.path.join(os.path.abspath(os.path.split(__file__)[0]),'ani.py')
        xtbmethod:int=2
        optloose:bool=True
        optconvergence:str="LOOSE"
        inputkeyfile:None=None
        writeoutmultipole:bool=True
        writeoutbond:bool=True
        writeoutangle:bool=True
        writeoutpolarize:bool=True
        writeouttorsion:bool=True
        dontrotbndslist:list=field(default_factory=lambda : [])
        indextompoleframefile:None=None
        potentialoffset:float=1.0
        indextotypefile:None=None
        usesymtypes: bool = True
        ldatabaseparserpath:str=os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'lDatabaseParser', 'lAssignAMOEBAplusPRM.py')
        ldatabaseparserprmdir:str=os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'lDatabaseParser', 'prm')
        numproc:None=None
        analyzepath:str='analyze'
        prmfilepath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18.prm"
        fennixmodeldir:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir, 'FennixModels'))
        fennixmodelname:str='ani2x_model0'
        fixvdwtyperadii:list=field(default_factory=lambda : [])
        maxjobsatsametime:float=10
        onlyrottortorlist:list=field(default_factory=lambda : [])
        numespconfs:int=1
        fitred:bool=False
        vdwprmtypestofit:list=field(default_factory=lambda : ['S','T'])
        lastlogfileupdatetime:float=1
        addlonepairvdwsites:bool=False
        quickdatabasesearch:bool=False
        genprotstatesonly:bool=False
        generateextendedconf:bool=True
        generate_symm_frag_conf:bool=False
        onlyvdwatomlist:None=None
        parentjobsatsametime:int=1
        coresperjob:int=2
        addhydrogens:bool=False
        maximizejobsatsametime:bool=True
        consumptionratio:float=.8
        nonaroringtor1Dscan:bool=False
        skipespfiterror:bool=False
        vdwmaxqmstartingpointspertype:int=1
        vdwmaxtinkergridpoints:int=50
        smallmoleculefragmenter:bool=False
        psi4_args:str=" --loglevel 30 "
        toroptdebugmode:bool=False
        debugmode:bool=False
        fragmenterdebugmode:bool=False
        jobsatsametime:int=0
        usepoleditframes:bool=True
        databasematchonly:bool=False
        setupfragjobsonly:bool=False
        allowradicals:bool=False
        checkinputonly:bool=False
        esprestweight:float=1
        espgrad:float=.1
        issane:bool=True
        deletedfiles:bool=False
        onlyfittorstogether:list=field(default_factory=lambda : [])
        parentname:None=None
        addhydrogentocharged:bool=True
        accuratevdwsp:bool=False
        firstoptfinished:bool=False
        optonly:bool=False
        onlyvdwatomindex:None=None
        use_qmopt_vdw:bool=False
        use_gau_vdw:bool=False
        pcm_auto:bool=True
        deleteallnonqmfiles:bool=False
        totalcharge:None=None
        torspbasissethalogen:str="6-311G*"
        homodimers:bool=False
        tortormissingfilename:str='tortormissing.txt'
        tordebugmode:bool=False
        parameterfilespath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir, 'ParameterFiles'))
        prmmodlist:list=field(default_factory=lambda : [])
        smartstosoluteradiimap:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'SMARTsToSoluteRadiiMap.txt'
        latestsmallmoleculepolarizeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarize.prm'
        updatedsmallmoleculepolarizeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] ))+'/lDatabaseParser/prm/'+'polarize.prm'
        amoebaplussmallmoleculenonbonded_prm:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] ))+'/lDatabaseParser/prm/'+'amoebaplusNonbonded.prm'
        amoebaplussmallmoleculenonbonded_dat:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] ))+'/lDatabaseParser/dat/'+'amoebaplusNonbondedType.dat'
        latestsmallmoleculesmartstotypespolarize:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarcommenttoparameters.txt'
        latestsmallmoleculesmartstotinkerclass:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21smartstoclass.txt'
        latestsmallmoleculeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21.prm'
        boltzmantemp:float=8
        dovdwscan:bool=False
        vdwprobepathname:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/Examples/Assets/'
        vdwprobenames:list=field(default_factory=lambda : ['water'])
        maxtorRMSPDRel:float=.1
        vdwmissingfilename:str='missingvdw.txt'
        databaseprmfilename:str='database.prm'
        tortor:bool=False
        torfit2Drotonly:bool=False
        torfit1Drotonly:bool=False
        externalparameterdatabase:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'externalparameterdatabase.txt'
        fitfirsttorsionfoldphase:bool=False
        skipgridsearch:bool=True
        torsionprmguessfilename:str='torsionprmguess.txt'
        defaultmaxtorsiongridpoints:int=40
        torsionsmissingfilename:str='torsionsmissing.txt'
        smallmoleculemm3prmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'mm3.prm'
        smallmoleculesmartstomm3descrip:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstomm3typedescrip.txt'
        absdipoletol:float=.5
        mmbondtol:float=.5
        mmangletol:float=.5
        transferanyhydrogentor:bool=True
        smallmoleculesmartstotinkerdescrip:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstoamoebatypedescrip.txt'
        smallmoleculeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba09.prm'
        torspbasissetfile:str='6-311+g_st_.0.gbs'
        toroptbasissetfile:str='6-31g_st_.0.gbs'
        optbasissetfile:str='6-31g_st_.0.gbs'
        dmabasissetfile:str='6-311g_st__st_.0.gbs'
        espbasissetfile:str='aug-cc-pvtz.1.gbs'
        iodinetorspbasissetfile:str='def2-svp.1.gbs'
        iodinetorspbasisset:str='def2-svp'
        iodinetoroptbasissetfile:str='def2-svp.1.gbs'
        iodinetoroptbasisset:str='def2-svp'
        iodineoptbasissetfile:str='def2-svp.1.gbs'
        iodineoptbasisset:str='def2-svp'
        iodinedmabasissetfile:str='def2-svp.1.gbs'
        iodineespbasissetfile:str='def2-tzvpp.1.gbs'
        iodineespbasisset:str='def2-tzvpp'
        basissetpath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/'+'BasisSets/'
        refinenonaroringtors:bool=False
        maxgrowthcycles:int=4
        use_gauPCM:bool=False
        fitqmdipole:bool=False
        scfmaxiter:int=500
        suppresstorfiterr:bool=False
        obminimizeexe:str='obminimize'
        suppressdipoleerr:bool=False
        poltypepath:str=os.path.abspath(os.path.split(__file__)[0])
        WBOtol:float=.05
        dontfrag:bool=False
        isfragjob:bool=False
        dipoletol:float=.5
        printoutput:bool=False
        poltypeini:bool=True
        structure:None=None
        prmstartidx:int=401
        maxmem:None=None
        maxdisk:None=None
        gausdir:None=None
        gdmadir:None=None
        tinkerdir:None=None
        scratchdir:str="/scratch"
        paramhead:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18_header.prm"
        gausexe:None=None
        formchkexe:str='formchk'
        cubegenexe:str='cubegen'
        gdmaexe:str='gdma'
        new_gdma:bool=False
        scaleandfixdipole:bool=False
        scalebigmultipole:bool=False
        fragbigmultipole:bool=False
        chargethreshold:float=1.5
        dipolethreshold:float=1.5
        quadrupolethreshold:float=2.4
        atomidsfordmafrag:list=field(default_factory=lambda : [])
        sameleveldmaesp:bool=False
        adaptiveespbasisset:bool=False
        avgmpolesexe:str=os.path.abspath(os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), os.pardir)) + "/PoltypeModules/avgmpoles.pl"
        peditexe:str='poledit.x'
        potentialexe:str='potential.x'
        minimizeexe:str='minimize.x'
        analyzeexe:str='analyze.x'
        superposeexe:str='superpose.x'
        defopbendval:float=0.20016677990819662
        Hartree2kcal_mol:float=627.5095
        optbasisset:str='6-31G*'
        toroptbasisset:str='6-311G'
        dmabasisset:str='6-311G**'
        iodinedmabasisset:str='def2-svp'
        espbasisset:str="aug-cc-pVTZ"
        torspbasisset:str="6-311+G*"
        optmethod:str='MP2'
        pyscf_opt_meth:str = 'wb97x_d3'
        pyscf_sol_imp:str = 'IEF-PCM' # C-PCM, SS(V)PE, COSMO
        pyscf_sol_eps:float = 78.3553 # Water
        toroptmethod:str='xtb'
        torspmethod:str='wB97X-D'
        dmamethod:str='MP2'
        espmethod:str='MP2'
        qmonly:bool = False
        espfit:bool = True
        foldnum:int=3
        foldoffsetlist:list = field(default_factory=lambda : [ 0.0, 180.0, 0.0, 180.0, 0.0, 180.0 ])
        torlist:None = None
        rotbndlist:None = None
        maxRMSD:float=1
        maxRMSPD:float=1
        maxtorRMSPD:float=1
        tordatapointsnum:None=None
        torsionrestraint:float=.5*3282.80354574
        torsionprmrestraintfactorL1:float=0.1
        torsionprmrestraintfactorL2:float=0.1
        onlyrotbndslist:list=field(default_factory=lambda : [])
        rotalltors:bool=False
        dontdotor:bool=False
        dontdotorfit:bool=False
        toroptpcm:int=-1
        optpcm:int=-1
        torsppcm:int=-1
        use_gaus:bool=False
        use_gausoptonly:bool=False
        use_psi4_geometric_opt:bool=True
        freq:bool=False
        postfit:bool=False
        bashrcpath:None=None
        optmaxcycle:int=400
        gausoptcoords:str=''
        forcefield:str="AMOEBA"
        sleeptime:float=.1
        structure:None=None
        espextraconflist:list=field(default_factory=lambda : [])
        userconformation:bool=False
        def __post_init__(self): 
            """
            Intent: Post initialization variables (things you want internal variables but not necesarrily user input). Also for reading input poltype.ini file and changing variable defaults.
            Input: Default variables from self object
            Output: Variables modified by input file and other variables internally used
            Referenced By: N/A 
            Description: 
            1. Initialize some internal variables not used by user input
            2. Read input poltype.ini file and change variable defaults
            """
            self.nfoldlist =  list(range(1,self.foldnum+1))
            self.parentdir=os.getcwd()
            if self.torlist==None:
                self.torlist = []
            else:
                self.torlist=self.torlist
            if self.rotbndlist==None:
                self.rotbndlist = []
            else:
                self.rotbndlist=self.rotbndlist
            self.torsettototalqmtime={}
            self.torsettooptqmtime={}
            self.torsettospqmtime={}
            self.torsettonumpoints={}
            gdma_kws = []
                    
            opts, xargs = getopt.getopt(sys.argv[1:],'h',["help"])

            for o, a in opts:
                if o in ("-h", "--help"):
                    self.copyright()
                    sys.exit(2)
                                
            if self.poltypeini==True:
                temp=open(os.getcwd()+r'/'+'poltype.ini','r')
                results=temp.readlines()
                temp.close()
                for line in results:
                    linesplit=line.split()
                    if line!='\n' and len(linesplit)>0:
                        if len(linesplit)>0:
                            if linesplit[0][0]=='#':
                                continue
                            else:
                                if '#' in line:
                                    for eidx in range(len(linesplit)):
                                        e=linesplit[eidx]
                                        if '#' in e:
                                            theidx=eidx
                                            break
                                    linesplit=linesplit[:theidx]
                                    line=' '.join(linesplit)
                        a=None
                        if '=' in line:
                            linesplit=line.split('=',1)
                            a=linesplit[1].replace('\n','').rstrip().lstrip()
                            newline=linesplit[0]
                            if a=='None':
                                continue
                            commalist=a.split(',')
                            commalist=[k.strip() for k in commalist]
                        else:
                            newline=line.strip()

                        if 'uncomplexedproteinpdbname' in newline:
                            self.uncomplexedproteinpdbname=a
                        elif 'inputkeyfile' in newline:
                            safekeyname='inputkey.key'
                            shutil.copy(a,safekeyname)
                            self.inputkeyfile=safekeyname
                        elif 'listofligands' in newline:
                            self.listofligands=a.split(',')
                            self.listofligands=[i.strip() for i in self.listofligands]
                        
                        elif 'toroptmethodlist' in newline:
                            self.toroptmethodlist=a.split(',')
                            self.toroptmethodlist=[i.strip() for i in self.toroptmethodlist]
                        elif 'torspmethodlist' in newline:
                            self.torspmethodlist=a.split(',')
                            self.torspmethodlist=[i.strip() for i in self.torspmethodlist]
                        elif 'indextotypefile' in newline:
                            self.indextotypefile=a
                        elif 'indextompoleframefile' in newline:
                            self.indextompoleframefile=a
                        elif 'xtbtorresconstant' in newline:
                            self.xtbtorresconstant=float(a)
                        elif 'targetenthalpyerror' in newline:
                            self.targetenthalpyerror=float(a)
                        elif 'targetdensityerror' in newline:
                            self.targetdensityerror=float(a)
                        elif 'potentialoffset' in newline:
                            self.potentialoffset=float(a)
                        elif 'enthalpyrelativeweight' in newline:
                            self.enthalpyrelativeweight=float(a)
                        elif 'densityrelativeweight' in newline:
                            self.densityrelativeweight=float(a)
                        elif 'templateligandxyzfilename' in newline:
                            self.templateligandxyzfilename=a
                        elif 'barinterval' in newline:
                            self.barinterval=int(a)
                        elif 'maxtorresnitrogen' in newline:
                            self.maxtorresnitrogen=int(a)
                        elif 'xtbmethod' in newline:
                            self.xtbmethod=int(a)
                        elif 'fennixmodelname' in newline:
                            self.fennixmodelname=str(a)
                        elif 'templateligandfilename' in newline:
                            self.templateligandfilename=a
                        elif 'prmfilepath' in newline:
                            self.prmfilepath=a
                        elif "useuniquefilenames" in newline:
                            self.useuniquefilenames=self.SetDefaultBool(line,a,True)
                        elif "writeoutmultipole" in newline:
                            self.writeoutmultipole=self.SetDefaultBool(line,a,True)
                        elif "writeoutbond" in newline:
                            self.writeoutbond=self.SetDefaultBool(line,a,True)
                        elif "writeoutangle" in newline:
                            self.writeoutangle=self.SetDefaultBool(line,a,True)
                        elif "optloose" in newline:
                            self.optloose=self.SetDefaultBool(line,a,True)
                            if self.optloose:
                                self.optconvergence = "LOOSE"
                        elif newline.startswith("optconvergence"):
                            self.optconvergence = a.upper()
                        elif "writeoutpolarize" in newline:
                            self.writeoutpolarize=self.SetDefaultBool(line,a,True)
                        elif "writeouttorsion" in newline:
                            self.writeouttorsion=self.SetDefaultBool(line,a,True)
                        elif "usesymtypes" in newline:
                            self.usesymtypes=self.SetDefaultBool(line,a,True)
                        elif "checktraj" in newline:
                            self.checktraj=self.SetDefaultBool(line,a,True)
                        elif "generateinputfilesonly" in newline:
                            self.generateinputfilesonly=self.SetDefaultBool(line,a,True)
                        elif 'bashrcpath' in newline:
                            self.bashrcpath=a
                        elif 'espextraconflist' in newline:
                            self.espextraconflist=a.split(',')
                            self.espextraconflist=[i.strip() for i in self.espextraconflist]
                        elif ("polareps") in newline:
                            self.polareps = a

                        elif "rotalltors" in newline:
                            self.rotalltors=self.SetDefaultBool(line,a,True)

                        elif "vdwprmtypestofit" in newline:
                            self.vdwprmtypestofit=a.split(',')
                            self.vdwprmtypestofit=[i.strip() for i in self.vdwprmtypestofit]

                        elif "fixvdwtyperadii" in newline:
                            self.fixvdwtyperadii=a.split(',')
                            self.fixvdwtyperadii=[i.strip() for i in self.fixvdwtyperadii]

                        elif "addlonepairvdwsites" in newline:
                            self.addlonepairvdwsites=self.SetDefaultBool(line,a,True)
                        elif "genprotstatesonly" in newline:
                            self.genprotstatesonly=self.SetDefaultBool(line,a,True)
                        elif "quickdatabasesearch" in newline:
                            self.quickdatabasesearch=self.SetDefaultBool(line,a,True)
                        elif "fitred" in newline:
                            self.fitred=self.SetDefaultBool(line,a,True)
                        elif "onlyvdwatomlist" in newline:
                            self.onlyvdwatomlist=a.split(',')
                            self.onlyvdwatomlist=[i.strip() for i in self.onlyvdwatomlist]
                            self.onlyvdwatomlist=[int(i) for i in self.onlyvdwatomlist]

                        elif "usepoleditframes" in newline:
                            self.usepoleditframes=self.SetDefaultBool(line,a,True)
                        elif "generateextendedconf" in newline:
                            self.generateextendedconf=self.SetDefaultBool(line,a,True)
                        elif "userconformation" in newline:
                            # for those who would like to use the conf in the input file direcly with torsion constrianed during opt 
                            self.userconformation=self.SetDefaultBool(line,a,True)
                        elif "generate_symm_frag_conf" in newline:
                            self.generate_symm_frag_conf=self.SetDefaultBool(line,a,True)
                        elif "addhydrogens" in newline:
                            self.addhydrogens=self.SetDefaultBool(line,a,True)
                        elif "nonaroringtor1Dscan" in newline:
                            self.nonaroringtor1Dscan=self.SetDefaultBool(line,a,True)
                        elif "fragmenterdebugmode" in newline:
                            self.fragmenterdebugmode=self.SetDefaultBool(line,a,True)
                        elif "skipespfiterror" in newline:
                            self.skipespfiterror=self.SetDefaultBool(line,a,True)
                        elif "smallmoleculefragmenter" in newline:
                            self.smallmoleculefragmenter=self.SetDefaultBool(line,a,True)
                        elif "debugmode" in newline and 'fragmenterdebugmode' not in line and 'tordebugmode' not in line and 'toroptdebugmode' not in line:
                            self.debugmode=self.SetDefaultBool(line,a,True)
                        elif "toroptdebugmode" in newline:
                            self.toroptdebugmode=self.SetDefaultBool(line,a,True)
                        elif "totalcharge" in newline:
                            self.totalcharge=int(a)
                        elif "lastlogfileupdatetime" in newline:
                            self.lastlogfileupdatetime=int(a)
                        elif "numespconfs" in newline:
                            self.numespconfs=int(a)
                        elif "consumptionratio" in newline:
                            self.consumptionratio=float(a)
                        elif "jobsatsametime" in newline and 'max' not in newline and 'parentjobsatsametime' not in newline:
                            self.jobsatsametime=int(a)
                        elif "esprestweight" in newline:
                            self.esprestweight=float(a)
                        elif "espgrad" in newline:
                            self.espgrad=a
                        elif "checkinputonly" in newline:
                            self.checkinputonly=True
                        elif "setupfragjobsonly" in newline:
                            self.setupfragjobsonly=self.SetDefaultBool(line,a,True)
                        elif "databasematchonly" in newline:
                            self.databasematchonly=self.SetDefaultBool(line,a,True)
                        elif "onlyvdwatomindex" in newline:
                            self.onlyvdwatomindex=int(a)
                        elif "parentjobsatsametime" in newline:
                            self.parentjobsatsametime=int(a)
                        elif "coresperjob" in newline:
                            self.coresperjob=int(a)
                        elif "deleteallnonqmfiles" in newline:
                            self.deleteallnonqmfiles=self.SetDefaultBool(line,a,True)
                        elif "addhydrogentocharged" in newline:
                            self.addhydrogentocharged=self.SetDefaultBool(line,a,True)
                        elif "firstoptfinished" in newline:
                            self.firstoptfinished=self.SetDefaultBool(line,a,True)
                        elif "accuratevdwsp" in newline:
                            self.accuratevdwsp=self.SetDefaultBool(line,a,True)
                        elif "homodimers" in newline:
                            self.homodimers=self.SetDefaultBool(line,a,True)
                        elif "optonly" in newline and 'gaus' not in newline:
                            self.optonly=self.SetDefaultBool(line,a,True)
                        elif "use_qmopt_vdw" in newline:
                            self.use_qmopt_vdw=self.SetDefaultBool(line,a,True)
                        elif "use_gau_vdw" in newline:
                            self.use_gau_vdw=self.SetDefaultBool(line,a,True)
                        elif "tortor" in newline and 'only' not in newline:
                            self.tortor=self.SetDefaultBool(line,a,True)
                        elif "tordebugmode" in newline:
                            self.tordebugmode=self.SetDefaultBool(line,a,True)
                        elif "refinenonaroringtors" in newline:
                            self.refinenonaroringtors=self.SetDefaultBool(line,a,True)
                        elif "fitfirsttorsionfoldphase" in newline:
                            self.fitfirsttorsionfoldphase=self.SetDefaultBool(line,a,True)
                        elif "fitqmdipole" in newline:
                            self.fitqmdipole=self.SetDefaultBool(line,a,True)
                        elif "maxgrowthcycles" in newline:
                            self.maxgrowthcycles=int(a)
                        elif 'maxtorRMSPDRel' in newline:
                            self.maxtorRMSPDRel=float(a)
                        elif "boltzmantemp" in newline:
                            self.boltzmantemp=float(a)
                        elif 'absdipoletol' in newline:
                            self.absdipoletol=float(a)
                        elif 'mmbondtol' in newline:
                            self.mmbondtol=float(a)
                        elif 'mmangletol' in newline:
                            self.mmangletol=float(a)
                        elif 'dipoletol' in newline:
                            self.dipoletol=float(a)
                        elif 'maxRMSD' in newline:
                            self.maxRMSD=float(a)
                        elif 'maxRMSPD' in newline:
                            self.maxRMSPD=float(a)
                        elif 'maxtorRMSPD' in newline:
                            self.maxtorRMSPD=float(a)
                        elif "use_gauPCM" in newline:
                            self.use_gauPCM=self.SetDefaultBool(line,a,True)
                        elif "espfit" in newline:
                            self.espfit=self.SetDefaultBool(line,a,True)
                        elif 'poltypepath' in newline and 'poltypepathlist' not in newline:
                            self.poltypepath=a
                        elif 'paramhead' in newline:
                            self.paramhead = a
                        elif 'WBOtol' in newline:
                            self.WBOtol=float(a)
                        elif 'scfmaxiter' in newline:
                            self.scfmaxiter=a
                        elif 'printoutput' in newline:
                            self.printoutput=self.SetDefaultBool(line,a,True)
                        elif 'suppressdipoleerr' in newline:
                            self.suppressdipoleerr=self.SetDefaultBool(line,a,True)
                        elif 'isfragjob' in newline:
                            self.isfragjob=self.SetDefaultBool(line,a,True)
                        elif "dontfrag" in newline:
                            self.dontfrag=self.SetDefaultBool(line,a,True)
                        elif "gausoptcoords" in newline:
                            self.gausoptcoords = a
                        elif "suppresstorfiterr" in newline:
                            self.suppresstorfiterr=self.SetDefaultBool(line,a,True)
                        elif "toroptbasisset" in newline:
                            self.toroptbasisset = a
                        elif "dmamethod" in newline:
                            self.dmamethod =a
                        elif 'new_gdma' in newline:
                            self.new_gdma=self.SetDefaultBool(line,a,True)
                        elif 'scaleandfixdipole' in newline:
                            self.scaleandfixdipole=self.SetDefaultBool(line,a,True)
                        elif 'scalebigmultipole' in newline:
                            self.scalebigmultipole=self.SetDefaultBool(line,a,True)
                        elif 'fragbigmultipole' in newline:
                            self.fragbigmultipole=self.SetDefaultBool(line,a,True)
                        elif newline.startswith("gdmacommand_"):
                            self.__dict__[newline] = a
                            gdma_kws.append((newline[len('gdmacommand_'):], a))
                        elif 'sameleveldmaesp' in newline:
                            self.sameleveldmaesp=self.SetDefaultBool(line,a,True)
                        elif 'adaptiveespbasisset' in newline:
                            self.adaptiveespbasisset=self.SetDefaultBool(line,a,True)
                        elif "bashrcpath" in newline and a!='None':
                            self.bashrcpath = a
                        elif "structure" in newline:
                            self.molstructfname = a
                        elif "dontusepcm" in newline:
                            self.pcm_auto= (not self.SetDefaultBool(line,a,True))
                        elif "pcm_auto" in newline:
                            self.pcm_auto=self.SetDefaultBool(line,a,True)
                        elif "freq" in newline:
                            self.freq=self.SetDefaultBool(line,a,True)
                        elif newline.startswith("optpcm"):
                            self.optpcm=self.GrabSwitchValue(a, 1)
                        elif newline.startswith("toroptpcm"):
                            self.toroptpcm=self.GrabSwitchValue(a, 1)
                        elif newline.startswith("torsppcm"):
                            self.torsppcm=self.GrabSwitchValue(a, 1)
                        elif "use_gaus" in newline and 'opt' not in newline:
                            self.use_gaus=self.SetDefaultBool(line,a,True)
                        elif "use_gausoptonly" in newline:
                            self.use_gausoptonly=self.SetDefaultBool(line,a,True)
                        elif "use_psi4_geometric_opt" in newline:
                            self.use_psi4_geometric_opt=self.SetDefaultBool(line,a,True)
                        elif "dontdotor" in newline:
                            self.dontdotor=self.SetDefaultBool(line,a,True)
                        elif "dovdwscan" in newline:
                            self.dovdwscan=self.SetDefaultBool(line,a,True)
                        elif "dontdotorfit" in newline:
                            self.dontdotorfit=self.SetDefaultBool(line,a,True)
                        elif "optmaxcycle" in newline:
                            self.optmaxcycle = int(a)
                        elif "torsionrestraint" in newline:
                            self.torsionrestraint=float(a)
                        elif "torsionprmrestraintfactorL1" in newline:
                            self.torsionprmrestraintfactorL1=float(a)
                        elif "chargethreshold" in newline:
                            self.chargethreshold=float(a)
                        elif "dipolethreshold" in newline:
                            self.dipolethreshold=float(a)
                        elif "quadrupolethreshold" in newline:
                            self.quadrupolethreshold=float(a)
                        elif "torsionprmrestraintfactorL2" in newline:
                            self.torsionprmrestraintfactorL2=float(a)
                        elif 'maxtorRMSPDRel' in newline:
                            self.maxtorRMSPDRel=float(a)
                        elif "foldnum" in newline:
                            self.foldnum=int(a)
                            self.nfoldlist =  list(range(1,self.foldnum+1))

                        elif "tordatapointsnum" in newline:
                            self.tordatapointsnum=int(a)
                        elif "onlyrotbndslist" in newline:
                            self.onlyrotbndslist=a.split(',')
                            templist=[]
                            for ele in self.onlyrotbndslist:
                                nums=ele.lstrip().rstrip().split()
                                temp=[]
                                for e in nums:
                                    temp.append(int(e))
                                templist.append(temp)
                            self.onlyrotbndslist=templist
                        elif "atomidsfordmafrag" in newline:
                            self.atomidsfordmafrag=a.split(',')
                        elif "dontrotbndslist" in newline:
                            self.dontrotbndslist=a.split(',')
                            templist=[]
                            for ele in self.dontrotbndslist:
                                nums=ele.lstrip().rstrip().split()
                                temp=[]
                                for e in nums:
                                    temp.append(int(e))
                                templist.append(temp)
                            self.dontrotbndslist=templist

                        elif "onlyrottortorlist" in newline:
                            self.onlyrottortorlist=a.split(',')
                            templist=[]
                            for ele in self.onlyrottortorlist:
                                nums=ele.lstrip().rstrip().split()
                                temp=[]
                                for e in nums:
                                    temp.append(int(e))
                                templist.append(temp)
                            self.onlyrottortorlist=templist

                        elif "onlyfittorstogether" in newline:
                            self.onlyfittorstogether=a.split(',')
                            templist=[]
                            for ele in self.onlyfittorstogether:
                                nums=ele.lstrip().rstrip().split()
                                temp=[]
                                for e in nums:
                                    temp.append(int(e))
                                templist.append(temp)
                            self.onlyfittorstogether=templist
                            self.onlyfittorstogether=[tuple(i) for i in self.onlyfittorstogether]

                        elif newline.startswith('prmmodfile'):
                            self.prmmodlist = []
                            for fpath1 in a.split(','):
                                fpath2list = (fpath1, os.path.join(self.parameterfilespath, fpath1),
                                    fpath1+'.mod', os.path.join(self.parameterfilespath, fpath1+'.mod'))
                                fpath3list = list(OrderedDict.fromkeys([os.path.abspath(fpath2) for fpath2 in fpath2list if os.path.isfile(fpath2)]))
                                if len(fpath3list) > 0:
                                    if len(fpath3list) > 1:
                                        warnings.warn(f"Multiple paths found ({len(fpath3list)}) for prmmod file '{fpath1}'; Using '{fpath3list[0]}'")
                                    self.prmmodlist.append(fpath3list[0])
                                else:
                                    warnings.warn(f"Could not locate prmmod file '{fpath1}'")
                        elif newline.startswith("optmethod"):
                            self.optmethod = a
                        elif newline.startswith("pyscf_opt_met"):
                            self.pyscf_opt_met = a
                        elif newline.startswith("pyscf_sol_imp"):
                            self.pyscf_sol_imp = a
                        elif newline.startswith("pyscf_sol_eps"):
                            self.pyscf_sol_eps = a
                        elif newline.startswith("espmethod"):
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
                        elif "prmstartidx" in newline:
                            self.prmstartidx = int(a)
                        elif "atmidx" in newline:
                            self.prmstartidx = int(a)
                        elif 'defaultmaxtorsiongridpoints' in newline:
                            self.defaultmaxtorsiongridpoints=int(a)
                        elif newline.startswith("optbasisset"):
                            self.optbasisset = a
                        elif "dmabasisset" in newline:
                            self.dmabasisset = a
                        elif "parentname" in newline:
                            self.parentname = a
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
                        elif "forcefield" in newline:
                            self.forcefield = a
                        elif "qmonly" in newline:
                            self.qmonly=self.SetDefaultBool(line,a,True)
                        elif "sleeptime" in newline:
                            self.sleeptime = float(a)
                        else:
                            print('Unrecognized '+line)
                            sys.exit()
            
            if self.new_gdma:
                self.gdmainp = get_dma_default('dma4')
                # automatically turn on special treatment if dma4 is used
                # this is the treatment applied to COO atoms
                self.scaleandfixdipole=True
            else:
                self.gdmainp = get_dma_default('dma0')
            self.gdmainp.update(gdma_kws)

            # Downgrade ESP basis set for molecules of 20 or more heavy atoms 
            if self.adaptiveespbasisset:
                m = Chem.MolFromMolFile(self.molstructfname,removeHs=False)
                num_of_heavy_atoms = Descriptors.HeavyAtomCount(m)
                if (num_of_heavy_atoms >= 20) and (self.espbasisset.upper() == "AUG-CC-PVTZ"):
                  self.espbasisset = "AUG-CC-PVDZ"
        
            # For iodine-containing molecule, keep using different level of DMA and ESP
            # for use with psi4
            if self.sameleveldmaesp and (not self.use_gaus):
                m = Chem.MolFromMolFile(self.molstructfname,removeHs=False)
                for i in range(m.GetNumAtoms()):
                  atomi = m.GetAtomWithIdx(i)
                  if atomi.GetAtomicNum() == 53:
                    self.sameleveldmaesp=False
                    
            files=os.listdir()
            foundfinal=False
            foundfinaltor=False
            for f in files:
                if 'final.key' in f:
                    foundfinal=True
                elif 'postfittorsion.key' in f:
                    foundfinaltor=True
            if foundfinal==True and foundfinaltor==True:
                self.deleteallnonqmfiles=False # keep this on during development phase        
            if len(self.toroptmethodlist)==0:
                self.toroptmethodlist.append(self.toroptmethod)
            if len(self.torspmethodlist)==0:
                self.torspmethodlist.append(self.torspmethod)
            self.temptoroptmethodlist=self.toroptmethodlist[:]
            self.temptorspmethodlist=self.torspmethodlist[:]

            if self.jobsatsametime!=0:
                self.maximizejobsatsametime=False
            else:
                if self.isfragjob==True and self.jobsatsametime==0:
                    self.jobsatsametime=1
            if self.molstructfname!=None:
                head, self.molstructfname = os.path.split(self.molstructfname)
                self.molecprefix =  os.path.splitext(self.molstructfname)[0]
            self.SanitizeAllQMMethods()
            self.SanitizeMMExecutables()
            self.copyright()
            self.initialize()
            self.init_filenames()
            p = Path(self.scratchdir)
            self.scratchpath=p.parent.absolute()
            if self.maxdisk==None and self.molstructfname!=None:
                stat= os.statvfs(self.scratchpath) 
                gb=stat.f_bfree*stat.f_bsize*10**-9
                gb=str(int(int(gb*self.consumptionratio)/self.parentjobsatsametime))
                self.maxdisk=gb+'GB'
            if self.maxmem==None:
                b=psutil.virtual_memory().available
                gb=b*10**-9
                gb=str(int(int(gb*self.consumptionratio)/self.parentjobsatsametime))
                self.maxmem=gb+'GB'
               
            if self.numproc==None:
                cpu=multiprocessing.cpu_count()
                cpu=str(int(int(cpu*self.consumptionratio)/self.parentjobsatsametime))
                self.numproc=cpu
            if self.maximizejobsatsametime==True and self.isfragjob==False:
                self.jobsatsametime=math.floor(int(self.numproc)/self.coresperjob)
                if self.jobsatsametime>self.maxjobsatsametime:
                    self.jobsatsametime=self.maxjobsatsametime
            self.partition=False
            self.tempjobsatsametime=self.jobsatsametime
            self.tempmaxmem=self.maxmem
            self.tempmaxdisk=self.maxdisk
            self.tempnumproc=self.numproc
            self.firsterror=False
            if self.sameleveldmaesp==True:
                self.dmamethod=self.espmethod
                self.dmabasisset=self.espbasisset
                self.iodinedmabasissetfile=self.iodineespbasissetfile
                self.iodinedmabasisset=self.iodineespbasisset
            if (self.dmamethod.upper() == self.espmethod.upper()) and (self.dmabasisset.upper() == self.espbasisset.upper()):
                self.sameleveldmaesp=True
            if self.debugmode==True:
                self.optmethod="HF"      
                self.toroptmethod="HF"         
                self.torspmethod="HF"                    
                self.dmamethod="HF"                   
                self.espmethod="HF"             
                self.optbasisset = 'MINIX'         
                self.toroptbasisset = 'MINIX'    
                self.dmabasisset = 'MINIX'             
                self.espbasisset = 'MINIX'        
                self.torspbasisset = 'MINIX'

            self.cmdtopid={} # for killing pids that have jobs stalled too long
            
            if self.isfragjob==False:
                self.parentname=self.molecprefix
            else:
                self.parentname=str(self.parentname)
            
 
            # Use openbabel to create a 'mol' object from the input molecular structure file. 
            # Openbabel does not play well with certain molecular structure input files,
            # such as tinker xyz files. (normal xyz are fine)
    
            
            if not __name__ == '__main__':
                params=self.main()


        def GrabIonizationStates(self,m):
            """
            Intent: Determine ionization-protonation states from input molecule using Dimorphite-DL.
            Input: Rdkit molecule object.
            Output: MOL and SDF files for dominant ionization states.
            Referenced By: GenerateParameters 
            Description: 
            1. Convert input molecules to SMILES.
            2. Convert back to MOL object.
            3. Call dimorphite_dl with input SMILES and pH range.
            4. For each MOL object geneated by Dimorphite, now convert to MOL and SDF files. 
            """
            smi=Chem.MolToSmiles(m)
            # STEP 1
            smiles=[smi]
            # STEP 2
            mols = [Chem.MolFromSmiles(s) for s in smiles]
            for i, mol in enumerate(mols):
                mol.SetProp("msg","Orig SMILES: " + smiles[i])
                
            # STEP 3 
            protonated_mols = dimorphite_dl.run_with_mol_list(
                mols,
                min_ph=6.9,
                max_ph=7.2,
            )
            
            finalsmi=[Chem.MolToSmiles(m) for m in protonated_mols]
            obConversion = openbabel.OBConversion()
            # STEP 4 
            for i,mol in enumerate(protonated_mols):
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                smi=finalsmi[i]
                name='IonizationState_'+str(i)+'.mol'
                rdmolfiles.MolToMolFile(mol,name)
                sdfmol=name.replace('.mol','.sdf')
                mol = openbabel.OBMol()
                obConversion.SetInFormat('mol')
                obConversion.ReadFile(mol,name)
                obConversion.SetOutFormat('sdf')
                obConversion.WriteFile(mol,sdfmol)

        def GrabTautomers(self,m):
            """
            Intent: Tautomers are another possible protonation state not handled by Dimorphite-DL. This is experimental (there may be better tools) using rdkit method to determine "canonical tautomer" or what it considers to be most probably tautomer. The first one generated (TautomerState_0..) is the canonical tautomer and the rest are just enumerated possilbilities.
            Input: Rdkit MOL object. 
            Output: MOL and SDF Files for tautomer states. 
            Referenced By: GenerateParameters 
            Description: 
            1. Call rdkit tautomer enumerator. 
            2. Convert molecule to SMILES and back to rdkit mol object. 
            3. Determine canonical tautormer.
            4. Enumerate tautomers from rdkit mol obkect. 
            5. Convert tautomers to SMILES.
            6. Sort SMILES so that the canonical tautomer is first.
            7. Iterate over SMILES and then output MOL and SDF files.
            """
            # STEP 1
            enumerator = rdMolStandardize.TautomerEnumerator()
            # STEP 2
            smi=Chem.MolToSmiles(m)
            m=Chem.MolFromSmiles(smi)
            # STEP 3
            canon = enumerator.Canonicalize(m)
            csmi = Chem.MolToSmiles(canon)
            res = [canon]
            # STEP 4
            tauts = enumerator.Enumerate(m)
            # STEP 5
            smis = [Chem.MolToSmiles(x) for x in tauts]
            # STEP 6
            stpl = sorted((x,y) for x,y in zip(smis,tauts) if x!=csmi)
            res += [y for x,y in stpl]
            obConversion = openbabel.OBConversion()
            # STEP 7
            for i in range(len(res)):
                taut=res[i]
                taut = Chem.AddHs(taut)
                AllChem.EmbedMolecule(taut)
                name='TautomerState_'+str(i)+'.mol'
                rdmolfiles.MolToMolFile(taut,name)
                sdfmol=name.replace('.mol','.sdf')
                mol = openbabel.OBMol()
                obConversion.SetInFormat('mol')
                obConversion.ReadFile(mol,name)
                obConversion.SetOutFormat('sdf')
                obConversion.WriteFile(mol,sdfmol)


        def SanitizeAllQMMethods(self):
            """
            Intent: Some qm methods between Gaussian and Psi4 have different syntax (like wB97XD vs wB97X-D), so just fix it depending on what user inputs and whatever user wants for QM package (Gaussian or Psi4 etc).
            Input: Internal variables to the self object.
            Output: Changed internal variables (relating to torspmethod,toroptmethod, optmethod etc). 
            Referenced By: _post_init, GenerateParameters()
            Description: 
            1. Call SanitizeQMMethod for optmethod, toroptmethod, torspmethod, dmamethod, espmethod.
            """
            self.optmethod=self.SanitizeQMMethod(self.optmethod,True)                 
            self.toroptmethod=self.SanitizeQMMethod(self.toroptmethod,True)                  
            self.torspmethod=self.SanitizeQMMethod(self.torspmethod,False)                    
            self.dmamethod=self.SanitizeQMMethod(self.dmamethod,False)                      
            self.espmethod=self.SanitizeQMMethod(self.espmethod,False)     


        def GrabBoolValue(self, value):
            """
            Intent: Poltype uses can change variables in input file, this function reads the boolean string.
            Input: String from poltype input file.
            Output: Boolean value for variable.
            Referenced By: _post_init
            Description: 
            1. Convert string to all lower case, if true, set to True, if false set to False.
            """
            # STEP 1
            if value.lower() == 'true':
                return True
            if value.lower() == 'false':
                return False
            raise ValueError('Could not convert "{}" into a boolean!'.format(value))

        @staticmethod
        def GrabSwitchValue(value, default=1):
            """ Convert switch option into integer

            0 <- False, no, off
            1 <- True, yes, on
            -1 <- auto
            default <- None
            """
            if isinstance(value, str):
                if value.lower() in ('false', 'no', 'off', '0'):
                    return 0
                elif value.lower() in ('true', 'yes', 'on', '1'):
                    return 1
                elif value.lower() in ('auto', '-1'):
                    return -1
                else:
                    return default
            elif isinstance(value, int) or isinstance(value, bool):
                return int(value)
            else:
                return default

        def SetDefaultBool(self,line,a,truthvalue):
            """
            Intent: Set default boolean value for variables if =True or = False is not given in the input file.
            Input: Line from input file, default truth value. 
            Output: Boolean value.
            Referenced By: _post_init
            Description:
            1. If truth value not in input line, then set to default truth value.
            2. If truth value in input line, then set to that truth value.
            """
            # STEP 1
            if '=' not in line:
                value = truthvalue
            # STEP 2
            else:
                value=self.GrabBoolValue(a)
            return value


        def SanitizeQMMethod(self,method,optmethodbool):
            """
            Intent: Some qm methods between Gaussian and Psi4 have different syntax (like wB97XD vs wB97X-D), so just fix it depending on what user inputs and whatever user wants for QM package (Gaussian or Psi4 etc).
            Input: QM method string, boolean to determine how to output new QM method.  
            Output: New QM method string.
            Referenced By: SanitizeAllQMMethods
            Description:
            1. If DFT method.
            2. If - in method.
            3. If using Gaussian replace -D with D.
            4. If using Psi4, eSanitizeAllQMMethodsce D with -D
            """
            # STEP 1
            if method[-1]=='D': # assume DFT, gaussian likes D, PSI4 likes -D
                # STEP 2
                if method[-2]=='-':
                    # STEP 3
                    if self.use_gaus or self.use_gausoptonly and optmethodbool==True:
                        method=method.replace('-D','D')
                    if self.use_gaus and optmethodbool==False:
                        method=method.replace('-D','D')

                else:
                    # STEP 4
                    if not self.use_gaus and not self.use_gausoptonly and optmethodbool==True:
                        method=method.replace('D','-D')
                    if not (self.use_gaus) and optmethodbool==False:
                        method=method.replace('D','-D')
            return method
 

        def WriteToLog(self,string,prin=False):
            """
            Intent: Write string to log file.   
            Input: String.
            Output: N/A
            Referenced By: Many functions throughout poltype. 
            Description: 
            1. Grab currrent time.
            2. Write out current time and input string to log file.
            3. Flush to disk.
            """
            # STEP 1
            now = time.strftime("%c",time.localtime())
            if not isinstance(string, str):
                string=string.decode("utf-8")
            # STEP 2
            self.logfh.write(now+' '+string+'\n')
            # STEP 3
            self.logfh.flush()
            os.fsync(self.logfh.fileno())
            if prin==True:
                print(now+' '+ string + "\n")
           

        def SanitizeMMExecutables(self):
            """
            Intent: Sometimes tinker installations have .x in tinker binary executable (analyze vs analyze.x), so detect which one is needed.
            Input: Self object variables (containing executable names).
            Output: Updated self object variables with new executable names. 
            Referenced By: _post_init
            Description:
            1. Call SanitizeMMExecutable for all tinker executables used. 
            """
            self.peditexe=self.SanitizeMMExecutable(self.peditexe)
            self.potentialexe=self.SanitizeMMExecutable(self.potentialexe)
            self.minimizeexe=self.SanitizeMMExecutable(self.minimizeexe)
            self.analyzeexe=self.SanitizeMMExecutable(self.analyzeexe)
            self.superposeexe=self.SanitizeMMExecutable(self.superposeexe)
            self.analyzepath=self.SanitizeMMExecutable(self.analyzepath)

        

        def initialize(self):
            """
            Intent: Initialize all paths to needed executables
            Input: Self object with variable paths
            Output: Modified self object with variable paths
            Referenced By: main
            Description: 
            1. Determine if Gaussian is in PATH  
            2. Find cubegen,formchk executables in PATH
            3. Determine tinker version and check if it satisfies the current minimum tinker version requirement
            4. If using AMOEBA+ forcefield, then switch parameter header file to AMOEBA+ parameter header file.
            5. For tinker executables, join exectuable name with bin directory from PATH.
            6. If cant find tinker exeutable, then crash program.
            7. Search for GDMA in PATH, if cant find then crash program.
            8. Search for Psi4 in PATH, if cant find then crash program.
            9. If Gaussian scratch directory is defined and folder doesnt exist, make a folder.
            """
            # STEP 1
            self.foundgauss=False
            if self.gausdir is None:
                self.gausdir = ""
            if self.gausexe is None:
                self.gausexe = ""

            if self.which(os.path.join(self.gausdir, self.gausexe)) is None:
                for g in ["g16", "g09", "g03"]:
                    gausexe = self.which(os.path.join(self.gausdir, g))
                    if gausexe is not None:
                        self.gausexe = gausexe
                        self.gausdir = os.path.dirname(gausexe)
                        self.foundgauss=True
                        break
            # STEP 2
            if self.use_gaus or self.use_gausoptonly:
                if self.which(os.path.join(self.gausdir, self.gausexe)) is None:
                    print("ERROR: Invalid Gaussian directory: ", self.gausdir)
                    sys.exit(1)
                self.cubegenexe = os.path.join(self.gausdir,self.cubegenexe)
                self.formchkexe = os.path.join(self.gausdir,self.formchkexe)

            # STEP 3
            cmdstr=self.analyzeexe+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/Examples/Assets/'+'water.xyz'+' '+'-k'+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/Examples/Assets/'+'water.key'+' '+'e'+'>'+' '+'version.out'
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
                    linesplit = line.split()
                    self.versionnum = linesplit[2]
                    if version.parse(self.versionnum) >= version.parse("8.10.2"):
                        latestversion = True
                        break

            if not latestversion:
                #raise ValueError("Notice: Not latest working version of tinker (8.10.2)"+' '+os.getcwd())
                print("Notice: Not latest working version of tinker (8.10.2)"+' '+os.getcwd())
           
            # STEP 4
            if self.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
                self.paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebaplus24_header.prm"
            # STEP 5
            if ("TINKERDIR" in os.environ):
                self.tinkerdir = os.environ["TINKERDIR"]
                self.peditexe = os.path.join(self.tinkerdir,self.peditexe)
                self.potentialexe = os.path.join(self.tinkerdir,self.potentialexe)
                self.minimizeexe = os.path.join(self.tinkerdir,self.minimizeexe)
                self.analyzeexe = os.path.join(self.tinkerdir,self.analyzeexe)
                self.superposeexe = os.path.join(self.tinkerdir,self.superposeexe)
            # STEP 6 
            if (not self.which(self.analyzeexe)):
                print("ERROR: Cannot find TINKER analyze executable")
                sys.exit(2)
                
                
                
            # STEP 7
            if self.gdmadir is None:
                self.gdmadir = os.getenv("GDMADIR", default="")
            self.gdmaexe = os.path.join(self.gdmadir, self.gdmaexe) 
            if (not self.which(self.gdmaexe)):
                print("ERROR: Cannot find GDMA executable")
                sys.exit(2)

            # STEP 8
            if (not self.which('psi4')):
                print("ERROR: Cannot find PSI4 executable")
                sys.exit(2)
             
            # STEP 9
            if self.use_gaus or self.use_gausoptonly:
                if ("GAUSS_SCRDIR" in os.environ):
                    self.scratchdir = os.environ["GAUSS_SCRDIR"]
                    if not os.path.isdir(self.scratchdir):
                        os.mkdir(self.scratchdir)
        
        
        
        
        
        def init_filenames (self):
            """
            Intent: Initialize file names and scrath folders.
            Input: Self object.
            Output: Updated self object. 
            Referenced By: main
            Description: 
            1. If gaussian or psi4, then create scratch folder if folder does not already exist.
            2. Initialize all filenames. 
            """
        
            # STEP 1 
            if ("GAUSS_SCRDIR" in os.environ):
                self.scratchdir = os.environ["GAUSS_SCRDIR"].rstrip('//')
                if not os.path.isdir(self.scratchdir):
                    os.mkdir(self.scratchdir)

            else:
                if ("PSI_SCRATCH" in os.environ):
                    self.scratchdir = os.environ["PSI_SCRATCH"]
                    if not os.path.isdir(self.scratchdir):
                        os.mkdir(self.scratchdir)
            # STEP 2
            if self.molstructfname!=None: 
                self.FileNames()

        def FileNames(self):
            """
            Intent: Initialize filenames.
            Input: Self object.
            Output: Updated self object. 
            Referenced By: init_filenames
            Description:
            1. For each file name (logs, keys, xyz files etc...), call assign_filenames taking variable name and suffix as inputs.  
            """
            # STEP 1
            self.logfname = self.assign_filenames ( "logfname" , "-poltype.log")
            self.chkname = self.assign_filenames ( "chkname" , ".chk")
            self.fname = self.assign_filenames ( "fname" , ".fchk")
            self.gausfname = self.assign_filenames ( "gausfname" , ".log")
            self.firstgausoptfname = self.assign_filenames ( ".gausoptfname" , "-opt_1.log")
            self.gdmafname = self.assign_filenames ( "gdmafname" , ".gdmaout")
            self.keyfname = self.assign_filenames ( "keyfname" , "_prefitmultipole.key")
            self.keyfnamefrompoledit = self.assign_filenames ( "keyfname" , ".key")
            self.xyzfname = self.assign_filenames ( "xyzfname" , ".xyz")
            self.peditinfile = self.assign_filenames ( "peditinfile" , "-peditin.txt")
            self.superposeinfile = self.assign_filenames ( "superposeinfile" , "-superin.txt")
            self.espgrdfname = self.assign_filenames ( "espgrdfname" , ".grid")
            self.qmespfname = self.assign_filenames ( "qmespfname" , "_fortinker.cube")
            self.qmesp2fname = self.assign_filenames ( "qmesp2fname" , "_fortinker.pot")
            self.grpfname = self.assign_filenames ( "grpfname" , "-groups.txt")
            self.key2fname = self.assign_filenames ( "key2fname" , "_avgprefitmultipole.key")
            self.key2fnamefromavg = self.assign_filenames ( "key2fname" , ".key_2")
            self.key3fnamefrompot = self.assign_filenames ( "key3fname" , ".key_3")
            self.key3fname = self.assign_filenames ( "key3fname" , "_postfitmultipole.key")
            self.key4fname = self.assign_filenames ( "key4fname" , "_prefitvdw.key")
            self.key4bfname = self.assign_filenames ( "key4fname" , "_prefitvdw.keyb")
            self.key5fname = self.assign_filenames ( "key5fname" , "_postfitvdw.key")
            self.key6fname= self.assign_filenames ( "key6fname" , "_prefittorsion.key")
            self.key7fname= self.assign_filenames ( "key7fname" , "_postfittorsion.key")
            self.xyzoutfile = self.assign_filenames ( "xyzoutfile" , ".xyz_2")
            if self.isfragjob==False:
                self.scrtmpdirgau = self.scratchdir.rstrip('//') + '/Gau-' + self.molecprefix
                self.scrtmpdirpsi4 = self.scratchdir.rstrip('//') + '/Psi4-' + self.molecprefix
            else:
                self.scrtmpdirgau = self.scratchdir.rstrip('//') + '/Gau-' + self.molecprefix+'-'+self.parentname
                self.scrtmpdirpsi4 = self.scratchdir.rstrip('//') + '/Psi4-' + self.molecprefix+'-'+self.parentname


            self.tmpxyzfile = 'final.xyz'
            self.tmpkeyfile = 'final.key'
            self.comtmp = self.assign_filenames ( "comtmp" , "-tmp.com")
            self.firstcomoptfname = self.assign_filenames ( "comoptfname" , "-opt_1.com")
            self.firstchkoptfname = self.assign_filenames ( "chkoptfname" , "-opt_1.chk")
            self.firstfckoptfname = self.assign_filenames ( "fckoptfname" , "-opt_1.fchk")
            self.firstlogoptfname = self.assign_filenames ( "logoptfname" , "-opt_1.log")

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
            self.comoptfname=self.firstcomoptfname 
            self.chkoptfname=self.firstchkoptfname 
            self.fckoptfname=self.firstfckoptfname
            self.logoptfname=self.firstlogoptfname 
            self.gausoptfname=self.firstgausoptfname

        def assign_filenames (self,filename,suffix):
            """
            Intent: Assign filename
            Input: filename string and suffix
            Output: filename variable
            Referenced By: FileNames 
            Description:
            1. Return variablename of molecule prefix and input suffix
            """
            # STEP 1
            if filename in globals():
                return eval(filename)
            else:
                return self.molecprefix + suffix
        
        def copyright(self):
            """
            Intent: Print software version and copyright info
            """
            print("Poltype -- Polarizable atom typer of small molecules for the AMOEBA polarizable force field")
            print("Please cite:")
            print("B. Walker, C. Liu, E. Wait, P. Ren, J. Comput. Chem. 2022, 43(23), 1530")
            print("Version 2.3.1 Feb 2025")
            print("Copyright (c)  Johnny Wu, Gaurav Chattree, Brandon Walker, Matthew Harger and Pengyu Ren 2019-2025")
            print("All Rights Reserved")
            print("##############################################################################################################")
             
            
        def CheckIsInput2D(self,mol,obConversion,rdkitmol):
            """
            Intent: Check if input structure has only 2D coordinates, if so then generate 3D coordinates
            Input: Openbabel mol object, openbabel conversion object, rdkit mol object
            Output: Updated mol objects (if they needed 3D coordinates generated) 
            Referenced By: GenerateParameters() 
            Description:
            1. Determine if 2D only by checking if all Z coordinates are 0.
            2. If detect 2D coordinates, generate new filename for 3D coordinates and update filenames.
            3. Call EmbedMolecule to generate 3D coordinates
            4. Write molecule to new filename.
            """
            # STEP 1
            is2d=True
            for atom in openbabel.OBMolAtomIter(mol):
                zcoord=atom.GetZ()
                if zcoord!=0:
                    is2d=False
            if is2d==True: 
                # STEP 2
                molprefix=self.molstructfname.split('.')[0]
                newname=molprefix+'_3D'+'.mol'
                self.molstructfname=newname
                self.molecprefix =  os.path.splitext(self.molstructfname)[0]
                self.FileNames()
                rdmolfiles.MolToMolFile(rdkitmol,'test.mol',kekulize=True)
                # STEP 3
                AllChem.EmbedMolecule(rdkitmol)
                # STEP 4
                rdmolfiles.MolToMolFile(rdkitmol,newname,kekulize=True)
                newmol = openbabel.OBMol()
                inFormat = obConversion.FormatFromExt(newname)
                obConversion.SetInFormat(inFormat)
                obConversion.ReadFile(newmol,newname)
            else:
                molprefix=self.molstructfname.split('.')[0]
                newname=molprefix+'.mol'
                rdmolfiles.MolToMolFile(rdkitmol,newname,kekulize=True)
                newmol = openbabel.OBMol()
                inFormat = obConversion.FormatFromExt(newname)
                obConversion.SetInFormat(inFormat)
                obConversion.ReadFile(newmol,newname)

        
            return newmol,rdkitmol


        def KillJob(self,job,outputlog):
            """
            Intent: Sometimes even after log file has termination, the process can still be zombie so kill the job.
            Input: self object, command string (job)
            Output: N/A
            Referenced By: CallJobsSeriallyLocalHost, WaitForTermination
            Description:
            1. Optain process from dictionary of command string to process
            2. Try to terminate with p.wait
            3. If that fails, try sending kill signal to all children of process
            """
            # STEP 1
            if job in self.cmdtopid.keys():
                p=self.cmdtopid[job]
                try:
                    # STEP 2
                    p.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    # STEP 3
                    process = psutil.Process(p.pid)
                    for proc in process.children(recursive=True):
                        proc.kill()
                    process.kill()

        
        def CallJobsSeriallyLocalHost(self,fulljobtooutputlog,skiperrors,wait=False):
           """
           Intent: Call jobs serially if jobsatsametime=1, else submit in parralel. 
           Input: Dictionary of command string to output log, boolean if should skip errors, boolean if should wai Dictionary of command string to output log, boolean if should skip errors, boolean if should wait for process to finish before moving on.
           Output: List of finished jobs, list of errorjobs. 
           Referenced By: Many functions in poltype
           Description: 
           1. Initialize empty arrays for fininshed and error jobs and submitted jobs.
           2. Remove any log files from previous runs that failed.
           3. While the number of finished jobs are not the same as number of input jobs
           4. Iterative over jobs, 
           5. If job is not finished (in finished array), if job is finished call KillJob, append to finishedjobs array, remove from submittedjobs array. Write to a log file that job finished successfully.
           6. If job has an error, call KillJob, append to errorjobs and finishedjobs array, remove from submittedjobs array. Write to log file that job had an error.
           7. If job isnt finished and submitted and number of current submitted jobs is not greater than number of allowed jobs at same time, then submit job.
           """
           thepath=os.path.join(os.getcwd(),'Fragments')
           # STEP 1
           finishedjobs=[]
           errorjobs=[]
           submittedjobs=[]
           errormessages=[]
           # STEP 2
           for job,outputlog in fulljobtooutputlog.items():
               finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
               if error==True: # remove log before resubmitting
                   os.remove(outputlog)
           # STEP 3
           while len(finishedjobs)!=len(list(fulljobtooutputlog.keys())):
               # STEP 4
               for job,outputlog in fulljobtooutputlog.items():
                   if job not in finishedjobs:
                      
                      finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
                      # STEP 5
                      if finished==True:
                          self.KillJob(job,outputlog)
                          if outputlog not in finishedjobs:
                              finishedjobs.append(outputlog)
                              self.NormalTerm(outputlog)
                              if job in submittedjobs:
                                  submittedjobs.remove(job)
                      # STEP 6
                      if error==True and job in submittedjobs:
                          self.KillJob(job,outputlog) # sometimes tinker jobs or qm just hang and dont want zombie process 
                          if outputlog not in finishedjobs:
                              errorjobs.append(outputlog)
                              finishedjobs.append(outputlog) 
                              self.ErrorTerm(outputlog,skiperrors)
                              submittedjobs.remove(job)
                      # STEP 7
                      if job not in submittedjobs and len(submittedjobs)<self.jobsatsametime and finished==False and outputlog not in finishedjobs:
                          count=len(finishedjobs)
                          jobsleft=len(fulljobtooutputlog.keys())-count
                          ratio=round((100*count)/len(fulljobtooutputlog.keys()),2)
                          if len(finishedjobs)!=0:
                              if 'poltype.py' in job:
                                  self.ETAQMFinish(thepath,len(fulljobtooutputlog.keys()))
                          if len(list(fulljobtooutputlog.keys()))>1: 
                              self.WriteToLog('Percent of jobs finished '+str(ratio)+', jobs left = '+str(jobsleft))
                          self.call_subsystem([job],wait,skiperrors)
                          submittedjobs.append(job)
                       

               time.sleep(self.sleeptime)
           return finishedjobs,errorjobs

        def TabulateLogs(self,thepath):
            """
            Intent: Collect log files into arrays for later evaluating how much time QM took
            Input: Path to where QM log files are
            Output: listoffilepaths
            Referenced By: ETAQMFinish
            Description:
            1. List all files in directory
            2. Iterate over files, 
            3. If its a directory, then call GrabLogFilePaths and append file paths to array listoffilepaths.
            """
            listoffilepaths=[]
            os.chdir(thepath)
            # STEP 1
            files=os.listdir()
            # STEP 2
            for f in files:
                path=os.path.join(thepath,f)
                # STEP 3
                if os.path.isdir(path):
                    logfilepaths=self.GrabLogFilePaths(path)
                    if len(logfilepaths)!=0:
                        listoffilepaths.append(logfilepaths)
            return listoffilepaths
        
        def GrabLogFilePaths(self,path):
            """
            Intent: Given input directory, find all sub directories and QM log files within them.
            Input: Directory
            Output: QM log filepaths
            Referenced By: TabulateLogs
            Description:
            1. Iterate over all sub directories of input directory
            2. Determine if its a finished poltype job (has final.key), if not skip sub directory
            3. For any QM log files in top directory, append to array
            4. For any QM log files in subdirectories if current directory append to array.
            """
            filepaths=[]
            # STEP 1
            for root, subdirs, files in os.walk(path):
                # STEP 2
                finishedpoltype=False
                for f in files:
                    if 'final.key' in f:
                        finishedpoltype=True
                if finishedpoltype==False:
                    continue
                # STEP 3
                for f in files:
                    if '.log' in f and '_frag' not in f and 'poltype' not in f:
                        filepath=os.path.join(path, f)  
                        filepaths.append(filepath)
                # STEP 4
                for d in subdirs:
                    curdir=os.getcwd()
                    path=os.path.join(root, d)
                    os.chdir(path)
                    newfiles=os.listdir()
                    for f in newfiles:
                        if '.log' in f:
                            filepath=os.path.join(path, f)  
                            filepaths.append(filepath) 
            return filepaths
        
        def ConvertTimeString(self,timestring):
            """
            Intent: Convert time string from log file into one number in units of hours.
            Input: timestring
            Output: totaltime in hours
            Referenced By: TabulateCPUTimes
            Description:
            1. Grab hours string
            2. Grab minutes string, convert to hours
            3. Grab seconds string, convert to hours
            4. Add all times together
            """
            timestringsplit=timestring.split(':')
            # STEP 1
            hours=float(timestringsplit[0])
            # STEP 2
            minutes=float(timestringsplit[1])*(1/60)
            # STEP 3
            seconds=float(timestringsplit[2])*(1/60)*(1/60)
            # STEP 4
            totaltime=hours+minutes+seconds
            return totaltime


        def ConvertTimeStringXTB(self,timestring):
            """
            Intent: Convert time string from log file into one number in units of hours for XTB.
            Input: timestring
            Output: totaltime in hours
            Referenced By: TabulateCPUTimes
            Description:
            1. Grab days string, convert to hours
            2. Grab hours string
            3. Grab minutes string, convert to hours
            4. Grab seconds string, convert to hours
            5. Add all times together
            """
            timestringsplit=timestring.split(',')
            # STEP 1
            days=float(timestringsplit[0].strip().split()[0])*(24)
            # STEP 2
            hours=float(timestringsplit[1].strip().split()[0])
            # STEP 3
            minutes=float(timestringsplit[2].strip().split()[0])*(1/60)
            # STEP 4
            seconds=float(timestringsplit[3].strip().split()[0])*(1/60)*(1/60)
            # STEP 5
            totaltime=days+hours+minutes+seconds
            return totaltime


        
        def ConvertTimeStringGaussian(self,line):
            """
            Intent: Convert time string from log file into one number in units of hours for Gaussian.
            Input: timestring
            Output: walltime in hours
            Referenced By: TabulateCPUTimes
            Description: 
            1. Grab days string, convert to hours
            2. Grab hours string
            3. Grab minutes string, convert to hours
            4. Grab seconds string, convert to hours
            5. Add all times together
            """
            linesplit=line.split()
            # STEP 1
            days=float(linesplit[3])*24
            # STEP 2
            hours=float(linesplit[5])
            # STEP 3
            minutes=float(linesplit[7])*(1/60)
            # STEP 4
            seconds=float(linesplit[9])*(1/60)*(1/60)
            # STEP 5
            walltime=days+hours+minutes+seconds
            return walltime
        
        def TabulateCPUTimes(self,listoffilepaths):
            """
            Intent: Iterate over list of filepaths, determine wall time and append to array of wall times.
            Input: listoffilepaths
            Output: listoflistofcputimes
            Referenced By: ETAQMFinish
            Description:
            1. Initalize array of CPU wall times
            2. Iterate over array of array of filepaths
            3. Iterate over filepaths
            4. Read the timestring appropriately for each log file depending on which program was used (Gaussian,Psi4, XTB) 
            """
            # STEP 1
            listoflistofcputimes=[]
            # STEP 2
            for filepathsidx in range(len(listoffilepaths)):
                filepaths=listoffilepaths[filepathsidx]
                cputimearray=[]
                # STEP 3
                for logfile in filepaths:
                    if os.path.isfile(logfile):
                        with open(logfile) as frb:
                            # STEP 4
                            for line in frb:
                                if 'wall time' in line:
                                    timestring=line.split()[-1]
                                    walltime=self.ConvertTimeString(timestring) 
                                    cputimearray.append(walltime)
                                    break
                                elif "Job cpu time:" in line:
                                    cputime=self.ConvertTimeStringGaussian(line) 
                                    walltime=4*cputime # 4 processors
                                    cputimearray.append(walltime)
        
                                    break
                                elif 'wall-time:' in line:
                                    timestring=' '.join(line.split()[2:])
                                    walltime=self.ConvertTimeStringXTB(timestring)
                                    cputimearray.append(walltime)
        
                listoflistofcputimes.append(cputimearray)
            return listoflistofcputimes
        
        
        def ETAQMFinish(self,thepath,numberofmolstotal):
            """
            Intent: Attempt to give an ETA for QM log files to finish based on averaging time from previously finished log file. 
            Input: Fragments directory with subdirectory of QM log files, total QM jobs
            Output: Updated poltype log file with ETA for QM to finish
            Referenced By: CallJobsSeriallyLocalHost
            Description:
            1. Call TabulateLogs and grab listoffilepaths
            2. Call TabulateCPUTimes with listoffilepaths
            3. Add up the total time for each array in listoflistofcputimes
            4. Compute the min,max and average for QM times.
            5. Based on the jobs left for current jobs being submitted, and average QM time, compute avevage time to complete remaining jobs.
            6. Write the result to poltype log file
            """
            # STEP 1
            listoffilepaths=self.TabulateLogs(thepath)
            # STEP 2
            listoflistofcputimes=self.TabulateCPUTimes(listoffilepaths)
            # STEP 3
            listofcputimes=[sum(ls) for ls in listoflistofcputimes]
            if len(listofcputimes)==0:
                return 
            # STEP 4
            average=np.mean(listofcputimes)
            mintime=min(listofcputimes)
            maxtime=max(listofcputimes)
            ls=[mintime,average,maxtime]
            string='Min, Average, Max time for jobs to finish '+str(ls)+' hours'
            self.WriteToLog(string)
            # STEP 5
            jobsleft=int(numberofmolstotal)-len(listoffilepaths)
            timeleft=average*jobsleft
            # STEP 6
            string='There are '+str(jobsleft)+' jobs left, with ETA='+str(timeleft)+',hours for all jobs to finish'
            self.WriteToLog(string)

        def WaitForTermination(self,jobtooutputlog,skiperrors):
            """
            Intent: Wait for jobs to finish before executing other parts of program
            Input: dictionary of command strings to outputlofs, boolean if should skip errors (dont crash if there is error)
            Output: finishedjobs, errorjobs
            Referenced By: Many functions in poltype 
            Description:
            1. Initialize empty arrays
            2. Remove any previously submitted jobs that failed (from previous poltype runs)
            3. While the number of finished jobs is not the same as the number of input jobs
            4. Iterate over input jobs
            5. If job is finished without errors, then call KillJob and append to finishedjobs array, call NormalTerm
            6. If job is finished with errors, then call KillJob and append to finishedjobsarray and errorjobsarray, if skiperrors is False, then append to errorjobs and call ErrorTerm
            7. If job is not finished and has no errors, write to poltype log that waiting on job to finish
            8. After all jobs are finished write all jobs are finished to poltype log
            """
            # STEP 1
            finishedjobs=[]
            errorjobs=[]
            errormessages=[]
            outputStatusDict = copy.deepcopy(jobtooutputlog)
            # STEP 2
            for job,outputlog in jobtooutputlog.items():
               finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
               if error==True: # remove log before resubmitting
                   self.WriteToLog('Removing old outputlog that failed '+outputlog)
                   os.remove(outputlog)
            time.sleep(1)
            # STEP 3
            while len(finishedjobs)!=len(jobtooutputlog.keys()):
                # STEP 4
                for job in jobtooutputlog.keys():
                    outputlog=jobtooutputlog[job]
                    finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
                    # STEP 5
                    if finished==True and error==False: # then check if SP has been submitted or not
                        self.KillJob(job,outputlog)
                        if outputlog not in finishedjobs:
                            self.NormalTerm(outputlog)
                            finishedjobs.append(outputlog)
                    # STEP 6
                    elif finished==False and error==True:
                        self.KillJob(job,outputlog)
                        if skiperrors==False:
                            if outputlog not in finishedjobs:
                                self.ErrorTerm(outputlog,skiperrors)
                                finishedjobs.append(outputlog)
                                errorjobs.append(outputlog)
                        else:
                            finishedjobs.append(outputlog)
                    # STEP 7
                    elif finished==False and error==False:
                        if not os.path.isfile(outputlog):
                            printStr='Waiting on '+outputlog+' '+'to begin'
                        else:
                            printStr = 'Waiting on '+outputlog+' '+'for termination '
                        if (printStr != outputStatusDict[job]):
                            self.WriteToLog(printStr)
                            outputStatusDict[job] = printStr
                    else: # this case is finshed=True and error=True because there stupid quotes sometimes have word error in it                  
                        if outputlog not in finishedjobs:
                            if skiperrors==True:
                                error=False
                                self.NormalTerm(outputlog)
                                finishedjobs.append(outputlog)
        
                time.sleep(self.sleeptime) # how often to check logs (default every 30 s)
            # STEP 8
            self.WriteToLog('All jobs have terminated '+str(finishedjobs))
            return finishedjobs,errorjobs


        @staticmethod
        def DetectFileType(logfname):
            """
            Intent: check file type
            Input: log file path
            Referenced By: CheckLogfileStatus
            Description: 
            1. Define the patterns for each file type, and the detection threshold for number of patterns 
            2. Count the number of matched patterns, and return the file type of the threshold is reached 
            """
            file_type = 'other'
            patt_thr = {'psi4':3, 'gau':3}
            patt_count = {_:set() for _ in patt_thr}
            patt_list = []
            patt_list.append(('psi4', 'header', re.compile('^\s+Psi4\: An Open-Source Ab Initio Electronic Structure Package')))
            patt_list.append(('psi4', 'date', re.compile('^\s+Psi4 (started|stopped) on: ')))
            patt_list.append(('psi4', 'tstart', re.compile('^\*\*\* (tstart|tstop)\(\) called on')))
            patt_list.append(('psi4', 'time', re.compile('^\s+total time  =\s+\S+ seconds =\s+\S+ minutes')))
            patt_list.append(('psi4', 'basis', re.compile('^\s+=> Loading Basis Set <=')))
            patt_list.append(('psi4', 'normal_term', re.compile('^\*\*\* Psi4 exiting successfully.')))
            patt_list.append(('psi4', 'error_term', re.compile('^\*\*\* Psi4 encountered an error.')))
            patt_list = [_ for _ in patt_list if _[0] in patt_thr]
            with open(logfname) as fh:
                for line in fh:
                    for prog, entry, patt in patt_list:
                        if patt.match(line) is not None:
                            patt_count[prog].add(entry)
                            if len(patt_count[prog]) >= patt_thr[prog]:
                                file_type = prog
                                return file_type
            return file_type

        def CheckLogfileStatus(self, logfname, file_type=None):
            """
            Intent: check for normal terminatino in output files from psi4
            Input: log file path
            Referenced By: CheckNormalTermination
            Description: 
            1. Define the status of several items by using regex
            2. Go through each line of input file, and update the status of each item
            3. Determine whether the job is terminated and wether there is error
            """
            if file_type is None:
                file_type = self.DetectFileType(logfname)
            patt_map = {'psi4':[]}
            patt_map['psi4'].append(('normal_term', True, re.compile('^\*\*\* Psi4 exiting successfully.')))
            patt_map['psi4'].append(('error_term', True, re.compile('^\*\*\* Psi4 encountered an error.')))
            patt_map['psi4'].append(('error_term', True, re.compile('^Traceback \(most recent call last\):')))
            patt_map['psi4'].append(('opt', 'error', re.compile('(^\s+Optimization failed to converge!\s+~|Could not converge geometry optimization in \d+ iterations.)')))
            patt_map['psi4'].append(('opt', 'success', re.compile('(^\s+Optimization converged!\s+~|^\s+Final optimized geometry)')))
            log_status = {}
            if file_type in patt_map:
                with open(logfname) as fh:
                    for line in fh:
                        for entry, status, patt in patt_map[file_type]:
                            if patt.match(line) is not None:
                                log_status[entry] = status
            flag_error = log_status.get('error_term', False)
            flag_term = log_status.get('normal_term', False)
            if flag_term and log_status.get('opt', '') == 'error':
                flag_error = True
            return flag_term, flag_error

        def CheckNormalTermination(self,logfname,errormessages=None,skiperrors=False): 
            """
            Intent: Check for normal termination or errors in output files from tinker, psi4, gaussian, xtb and ANI. 
            Input: Output logname, boolean if should skip any errors. For torsion jobs, sometimes want to skip over optimization errors and remove point from dihedral energy surface so the whole program doesnt crash if a single point crashes.
            Output: Boolean if job terminated and boolean if job has error
            Referenced By: CallJobsSeriallyLocalHost
            Description: 
            1. Check the last update time of files in directory and store in dictionary.
            2. Determine latest updated filetime from all files
            3. Determine updatetime based on if torsion QM log file or other file. Sometimes QM jobs will crash but not send error signal or any output in outputlog, so then need to have update time if not terminated and not updated in updatetime then the job is determined to have failed.
            4. Read lines from output file and check for normal termination text strings or if there was an error.
            5. If error occured, write out line containing error to poltype log file.
            6. xtb only outputs the same filename for optimization jobs, so be sure to copy to desired filename after job finishes. 
            """
            
            error=False
            term=False
            lastupdatetofilename={}
            curdir=os.getcwd()
            logpath=os.path.split(logfname)[0]
            if logpath!='':
                os.chdir(logpath)
            # STEP 1
            files=os.listdir()
            for f in files:
                try:
                    Ftime=os.path.getmtime(f)
                    reltime=time.time()-Ftime
                    htime=reltime*0.000277778
                    lastupdatetofilename[htime]=f
                except:
                    pass
            os.chdir(curdir)
            # STEP 2
            htime=min(lastupdatetofilename.keys())
            if os.path.isfile(logfname):
                head,tail=os.path.split(logfname)
                # STEP 3
                if ('opt-' in logfname or 'sp-' in logfname):
                    updatetime=self.lastlogfileupdatetime
                else:
                    updatetime=2.5 # parent QM can take longer dependeing on resources, poltype log can take a while as well to finish QM
                foundendgau=False # sometimes gaussian comments have keyword Error, ERROR in them
                foundhf=False
                # STEP 4
                for line in open(logfname):
                    if 'poltype' in tail:
                        if 'Poltype Job Finished' in line:
                            term=True
                        elif 'Poltype has crashed!' in line:
                            error=True
                    else:
                        if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line or "LBFGS  --  Normal Termination due to SmallGrad" in line or "Normal termination" in line or 'Normal Termination' in line or 'Total Potential Energy' in line or 'Psi4 stopped on' in line or 'finished run' in line or ('Converged! =D' in line):
                            term=True
                        if ('Tinker is Unable to Continue' in line or 'error' in line or ' Error ' in line or ' ERROR ' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation, address not mapped to object' in line or 'galloc:  could not allocate memory' in line or 'Erroneous write.' in line or 'Optimization failed to converge!' in line or 'Geometry optimization is not converged' in line or 'Error' in line) and 'DIIS' not in line and 'mpi' not in line and 'RMS Error' not in line:
                            error=True
                            errorline=line
                        if 'segmentation violation' in line and 'address not mapped to object' not in line or 'Waiting' in line or ('OptimizationConvergenceError' in line and 'except' in line) or "Error on total polarization charges" in line or 'Erroneous write' in line:
                            error=False
                            continue
                        if ('Error termination request processed by link 9999' in line or 'Error termination via Lnk1e in' in line) or ('OptimizationConvergenceError' in line and 'except' not in line) or 'Could not converge geometry optimization' in line or 'SCFConvergenceError' in line or 'Incomplete Convergence due to' in line or 'Induced Dipoles are not Converged' in line:
                            error=True
                            errorline=line
                        if 'l9999.exe' in line and foundendgau==False:
                            foundendgau=True
                            preverror=error # did find error before comment?
                        if foundendgau==True:
                            if error==True and preverror==False: # then caused by Gaussian comment
                                error=False
                        if ('hf/MINIX' in line or 'HF/gen' in line) and '_frag' not in logfname:
                            foundhf=True
                file_type = self.DetectFileType(logfname)
                if file_type == 'psi4':
                    term, error = self.CheckLogfileStatus(logfname, file_type=file_type)
                    
                if self.debugmode==False and foundhf==True:
                    term=False
                        
                # STEP 5
                if error==True:
                    term=False # sometimes psi4 geometry opt not fully converge but says successfully exiting etc..
                    message='Error '+errorline+ 'logpath='+logfname
                if error==False and term==False and htime>=updatetime:
                    error=True
                    message='Error '+'Job has not been updated in '+str(updatetime)+' hours'+' last update time = '+str(htime)+' hours'+' logname='+logfname
                if error==True and term==False:
                    if errormessages!=None:
                        if message not in errormessages:
                            self.WriteToLog(message) 
                            errormessages.append(message)
                    else:
                        self.WriteToLog(message) 
            # STEP 6
            if self.toroptmethod=='xtb' and 'xtb-opt' in logfname and term==True:
                time.sleep(1) 
                cartxyz=logfname.replace('.log','.xyz')
                if not os.path.isfile(cartxyz) and os.path.isfile('xtbopt.xyz'):
                    shutil.copy('xtbopt.xyz',cartxyz)
            if errormessages!=None:
                return term,error,errormessages
            else:
                return term,error        


        
            
        def NormalTerm(self,logfname):
            """
            Intent: If job completes sucessfully write that it finished to poltype log file.
            Input: output job filename
            Output: String written to poltype log file.
            Referenced By: CallJobsSeriallyLocalHost, CheckNormalTermination
            Description: 
            1. Write termination signal to poltype log file 
            """
            # STEP 1
            self.WriteToLog("Normal termination: logfile=%s path=%s" % (logfname,os.getcwd()))
        
        
        def ErrorTerm(self,logfname,skiperrors):
            """
            Intent: If job has an error, write that the job has an error to poltype log file.
            Input: output job filename and boolean to crash program or not
            Output: String written to poltype log file.
            Referenced By: CallJobsSeriallyLocalHost, CheckNormalTermination
            Description: 
            1. Write error to poltype log file
            2. If boolean to skip errors is False, then crash the program
            """
            string="ERROR termination: logfile=%s path=%s" % (logfname,os.getcwd())
            self.WriteToLog(string)
            if skiperrors==False:
                raise ValueError(string)

        def call_subsystem(self,cmdstrs,wait=False,skiperrors=False):
            """
            Intent: Submit job on local node.
            Input: List of job command strings, boolean if should wait for job to complete before continuing program and boolean to crash or not crash program if receive error signal from subprocess.
            Output: - 
            Referenced By: CallJobsSeriallyLocalHost, many other functions throughout poltype
            Description: 
            1. Write to poltype log that will submit job
            2. Submit job on local host and save PID. Sometimes after job finishes zombie process remains so need PID to kill after seeing termination signal in output log file.
            3. If wait boolean is True then wait for job to finish executing
            4. If received an error code and boolean for skipping crashed jobs is False, then crash program. 
            """
            if self.printoutput==True:
                for cmdstr in cmdstrs:
                    print("Submitting: " + cmdstr+' '+'path'+' = '+os.getcwd())
            procs=[]
            for cmdstr in cmdstrs:
                if cmdstr!='':

                    # STEP 1
                    self.WriteToLog("Submitting: " + cmdstr+' '+'path'+' = '+os.getcwd())
                    # STEP 2
                    if 'pyscf' in cmdstr:
                        out = cmdstr.split()[-1]
                        out_pyscf = open(out, 'w')
                        p = subprocess.Popen(cmdstr,shell=True,stdout=out_pyscf, stderr=out_pyscf)
                    else:
                        p = subprocess.Popen(cmdstr,shell=True,stdout=self.logfh, stderr=self.logfh)

                    procs.append(p)
                    self.cmdtopid[cmdstr]=p

            if wait==True:
                # STEP 3
                exit_codes=[p.wait() for p in procs]
                for i in range(len(procs)):
                    exitcode=exit_codes[i]
                    cmdstr=cmdstrs[i] 
                    # STEP 4
                    if skiperrors==False:
                        if exitcode != 0:
                            if 'pyscf' in cmdstr:
                                out_pyscf.close()
                            self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                            raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                    else:
                        if exitcode != 0:
                            self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
            if 'pyscf' in cmdstr:
                out_pyscf.close()

        def WriteOutLiteratureReferences(self,keyfilename): 
            """
            Intent: Write references to key file
            Input: Keyfilename
            Output: Updated keyfile with references
            Referenced By: GenerateParameters 
            Description: 
            1. Read in key file contents to array
            2. Open temporary keyfile 
            3. Iterate over original keyfile lines and after passing atom block then write out references
            4. Replace temporary keyfile with original keyfile 
            """
            # STEP 1
            temp=open(keyfilename,'r')
            results=temp.readlines()
            temp.close()
            # STEP 2
            tempname=keyfilename.replace('.key','_temp.key')
            temp=open(tempname,'w')
            foundatomblock=False
            for i in range(len(results)):
                line=results[i]
                # STEP 3
                if 'atom' in line and foundatomblock==False:
                    foundatomblock=True
                    temp.write('#############################'+'\n')
                    temp.write('##                         ##'+'\n')
                    temp.write('##  Literature References  ##'+'\n')
                    temp.write('##                         ##'+'\n')
                    temp.write('#############################'+'\n')
                    temp.write('\n')
                    temp.write('Walker, B., Liu, C., Wait, E., Ren, P., J. Comput. Chem. 2022, 1. https://doi.org/10.1002/jcc.26954'+'\n')

                    temp.write('\n')
                    temp.write('Wu, J.C.; Chattree, G.; Ren, P.Y.; Automation of AMOEBA polarizable force field'+'\n')
                    temp.write('parameterization for small molecules. Theor Chem Acc.'+'\n')
                    temp.write('\n')
                    temp.write(line)
                else:
                    temp.write(line)
            # STEP 4
            os.remove(keyfilename)
            os.replace(tempname,keyfilename)    

        

        def WritePoltypeInitializationFile(self,poltypeinput):
            """
            Intent: Write poltype.ini file for fragmenter jobs
            Input: Dictionary of poltype variable strings and variable inputs
            Output: poltype.ini file
            Referenced By: FragmentJobSetup in fragmenter.py
            Description:
            1. Open filename poltype.ini
            2. Iterate over dictionary of poltype inputs and write out values to poltype.ini file 
            """
            # STEP 1
            inifilepath=os.getcwd()+r'/'+'poltype.ini'
            temp=open(inifilepath,'w')
            # STEP 2
            for key,value in poltypeinput.items():
                if value!=None:
                    line=key+'='+str(value)+'\n'
                    temp.write(line)
            temp.close()
            return inifilepath

        def CheckTorsionParameters(self,keyfilename,torsionsmissing): 
            """
            Intent: Sanity check to see if torsions parameters that wanted to be parameterized are all zero (something went wrong with fragmenter or fitting).
            Input: Keyfilename and array of torsions that were parameterized
            Output: -
            Referenced By: GenerateParameters 
            Description:
            1. Read contents of keyfile to array
            2. Iterate over array and search for lines with torsion parameters
            3. If torsion in line is one of the torsions that parameters were derived for check the force constant coefficents to see if all of them are zero
            4. If all of them are 0, then crash the program  
            """
            # STEP 1
            temp=open(keyfilename,'r')
            results=temp.readlines()
            temp.close()
            # STEP 2
            for lineidx in range(len(results)):
                line=results[lineidx]
                # STEP 3
                if line.strip().startswith('torsion') and '#' not in line and 'Missing' not in line and 'none' not in line and 'unit' not in line:
                    allzero=True
                    linesplit=line.split()
                    ls=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
                    revls=ls[::-1]
                    if (ls in torsionsmissing or revls in torsionsmissing):
                        pass
                    else:
                        continue
                    prms=linesplit[5:]
                    newprms=prms[0::3]
                    newprms=[float(i) for i in newprms]
                    for prm in newprms:
                        if prm!=0:
                            allzero=False
                    # STEP 4
                    if allzero==True:
                        self.WriteToLog("Warning: torsion parameters are all zero "+line+' path ='+os.getcwd())
                        #raise ValueError("torsion parameters are all zero "+line+' path ='+os.getcwd())

                        

        def main(self):
            """
            Intent: Call GenerateParameters.
            Input: Self object
            Output: - 
            Referenced By: RunPoltype
            Description: 
            1. Call GenerateParameters
            """   
            # STEP 1
            params= self.GenerateParameters()

        
        
        def GrabIndexToCoordinates(self,mol):
            """
            Intent: Rdkit by default doesnt always use input conformation, so read coordinates from openbabel MOL object into dictionary to later be used to force rdkit first conformor to have those coordinates.
            Input: Openbabel mol object
            Output: dictionary of atom index to coordinates
            Referenced By: GenerateParameters 
            Description: 
            1. Iterate over atoms in mol objcet
            2. Determine atom index, coordinates and save into dictionary
            """
            indextocoordinates={}
            # STEP 1
            iteratom = openbabel.OBMolAtomIter(mol)
            for atom in iteratom:
                # STEP 2
                atomidx=atom.GetIdx()
                rdkitindex=atomidx-1
                coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
                indextocoordinates[rdkitindex]=coords
            return indextocoordinates

        def AddInputCoordinatesAsDefaultConformer(self,m,indextocoordinates):
            """
            Intent: Set default conformor coordinates for rdkit
            Input: Rdkit mol object and dictionary of rdkit atom index to coordinates from input structure
            Output: Updated rdkit mol object
            Referenced By: GenerateParameters 
            Description:
            1. Initialize conformor object
            2. Iterate over rdkit atoms and set coordinates to coordinates from input dictionary 
            """
            # STEP 1
            conf = m.GetConformer()
            # STEP 2
            for i in range(m.GetNumAtoms()):
                x,y,z = indextocoordinates[i]
                conf.SetAtomPosition(i,Point3D(x,y,z))

            return m 

        def CheckIfCartesianXYZ(self,f):
            """
            Intent: Sometimes user will have cartesian XYZ for input structure in same folder but tinker generates same filename extension .xyz so need to delete cartesian XYZ file.
            Input: filename
            Output: boolean if cartesian XYZ or not
            Referenced By: RemoveCartesianXYZFiles
            Description:
            1. Assume file is cartesian XYZ
            2. If number of items in line are greater than 4, then it is not cartesian XYZ file 
            """
            # STEP 1
            check=True
            temp=open(f,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                linesplit=line.split()
                # STEP 2
                if len(linesplit)>4:
                    check=False
            return check

        def RemoveCartesianXYZFiles(self):
            """
            Intent: Remove cartesian XYZ files so that they dont interfere with tinker XYZ files that have same filename extension.
            Input: Self object
            Output: - 
            Referenced By: GenerateParameters
            Description: 
            1. Generate list of filesnames in current directory
            2. Iterate over files and check filename extension
            3. If filename extension is xyz and its not output from QM optimzation
            4. If filename is cartesian xyz (as opposed to tinker) then remove file 
            """
            # STEP 1
            files=os.listdir()
            for f in files:
                # STEP 2
                filename, file_extension = os.path.splitext(f)
                # STEP 3
                if file_extension=='.xyz' and 'opt' not in filename:
                    # STEP 4
                    check=self.CheckIfCartesianXYZ(f)
                    if check==True:
                        os.remove(f)

               
        
        def CheckInputCharge(self,molecule,verbose=False):
            """
            Intent: Assign formal charges to mol object so that users dont have to
            Input: rdkit mol object
            Output: Updated rdkit mol object and dictionary of rdkit atom index to formal charge 
            Referenced By: GenerateParameters 
            Description: 
            1. Initialize dictionary of atomic number to dictionary of valence to formal charge
            2. Iterate over rdkit atoms and grab atomic number, explicit valence and determine what the formal charge should be from dictionary of valence to formal charge for that atomic number.
            3. Special case if atom is carbon or nitrogen and if neighbors contain nitrogen, oyxgen or sulfur (polarizable atoms) then if carbon and explicit valence only 3, give formal charge of +1 (more stable then -1 case)
            4. If carbon has valence less than 4 and doesnt have polarizable neighbors then add hydrogens 
            5. If nitorgen has valence of 2 and no polariable neighbors, then add hydrogens
            6. If oxygen has valence of 1 and user specifies to protonate unfilled valence states with "addhydrogens" then add hydrogen to oxygen
            7. If nitrogen with valence of 2 and a radical is detected set the formal charge to 0 
            8. If oxygen with valence of 2 and a radical is detected set formal charge to 1 
            9. If oxygen with valence of 1 and a radical is detected set formal charge to 0
            10. Otherwise set atom formal charge to formal charge detected earlier 
            """
            array=[]
            totchg=0
            atomindextoformalcharge={}
            # STEP 1
            atomicnumtoformalchg={1:{2:1},5:{4:1},6:{3:-1},7:{2:-1,4:1},8:{1:-1,3:1},15:{4:1},16:{1:-1,3:1,5:-1},17:{0:-1,4:3},9:{0:-1},35:{0:-1},53:{0:-1}}
            for atom in molecule.GetAtoms():
                # STEP 2
                atomidx=atom.GetIdx()
                atomnum=atom.GetAtomicNum()
                val=atom.GetExplicitValence()
                valtochg=atomicnumtoformalchg[atomnum]
                radicals=atom.GetNumRadicalElectrons()
                if val not in valtochg.keys(): # then assume chg=0
                    chg=0
                else:
                    chg=valtochg[val]
                # STEP 3
                polneighb=False
                if atomnum==6 or atomnum==7:
                    for natom in atom.GetNeighbors():
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==7 or natomicnum==8 or natomicnum==16:
                            polneighb=True
                    if polneighb and val==3 and atomnum==6:
                        chg=1
                string='Atom index = '+str(atomidx+1)+' Atomic Number = ' +str(atomnum)+ ' Valence = '+str(val)+ ' Formal charge = '+str(chg)
                array.append(string)
                # STEP 4
                if atomnum==6 and val<4 and (self.addhydrogentocharged==True or self.addhydrogens==True)  and radicals==0 and polneighb==False:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Carbon, will assume missing hydrogens and add'+string) 
                        self.WriteToLog('WARNING! Strange valence for Carbon, will assume missing hydrogens and add '+string)
                    atom.SetNumRadicalElectrons(0)
                    chg=0
                    atom.SetFormalCharge(chg)

                # STEP 5
                elif atomnum==7 and val==2 and (self.addhydrogentocharged==True or self.addhydrogens==True) and radicals==0 and polneighb==False:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Nitrogen, will assume missing hydrogens and add'+string) 
                        self.WriteToLog('WARNING! Strange valence for Nitrogen, will assume missing hydrogens and add '+string)
                    atom.SetNumRadicalElectrons(0)
                    chg=0
                    atom.SetFormalCharge(chg)
                # STEP 6
                elif atomnum==8 and val==1 and (self.addhydrogens==True) and radicals==0:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Oxygen, will assume missing hydrogens and add'+string) 
                        self.WriteToLog('WARNING! Strange valence for Oxygen, will assume missing hydrogens and add '+string)
                    atom.SetNumRadicalElectrons(0)
                    chg=0
                    atom.SetFormalCharge(chg)

                # STEP 7
                elif atomnum==7 and val==2 and radicals==1:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Nitrogen, will assume radical and set charge to zero') 
                        self.WriteToLog('WARNING! Strange valence for Nitrogen, will assume radical and set charge to zero')
                    self.allowradicals=True

                    atom.SetFormalCharge(0)
                    self.addhydrogentocharged=False
                # STEP 8
                elif atomnum==8 and val==2 and radicals==1:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1') 
                        self.WriteToLog('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1')
                    self.allowradicals=True

                    atom.SetFormalCharge(1)
                # STEP 9
                elif atomnum==8 and val==1 and radicals==1:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1') 
                        self.WriteToLog('WARNING! Strange valence for Oxygen, will assume radical and set charge to +0')
                    self.allowradicals=True
                    atom.SetFormalCharge(0)
                    self.addhydrogentocharged=False

                else:
                    # STEP 10
                    atom.SetFormalCharge(chg)
                    if self.allowradicals==False:
                        atom.SetNumRadicalElectrons(0)
 
                    totchg+=chg
                atomindextoformalcharge[atomidx]=chg
            if self.totalcharge!=None:
                if self.totalcharge!=totchg:
                    for row in array:
                        print(row,flush=True)
                    raise ValueError('Valence is not consistent with input total charge')
            else:
                self.totalcharge=totchg 
            return molecule,atomindextoformalcharge


        def DetermineSymmetryMetric(self,rdkitmol,prin):
            rotbnd=self.onlyrotbndslist[0]
            a=rotbnd[0]-1
            b=rotbnd[1]-1
            # need to choose best neighbor to define a plane to project 3D coordinates onto
            allneighbs=[]
            for atomindex in rotbnd:
                atom=rdkitmol.GetAtomWithIdx(atomindex-1)
                neighbs=atom.GetNeighbors()
                for natom in neighbs:
                    natomindex=natom.GetIdx()
                    if natomindex!=a and natomindex!=b:
                        allneighbs.append(natom)
            foundc=False
            for natom in allneighbs:
                hyb=natom.GetHybridization()
                if hyb==Chem.HybridizationType.SP2:
                    c=natom.GetIdx()
                    foundc=True
                    break
            if foundc==False:
                c=allneighbs[0].GetIdx() # if not a ring, just pick any neighbor for plane right now. If see more examples can refine later
            
            apos = rdkitmol.GetConformer().GetAtomPosition(a) 
            A=np.array([float(apos.x),float(apos.y),float(apos.z)])
            bpos = rdkitmol.GetConformer().GetAtomPosition(b) 
            B=np.array([float(bpos.x),float(bpos.y),float(bpos.z)])
            cpos = rdkitmol.GetConformer().GetAtomPosition(c) 
            C=np.array([float(cpos.x),float(cpos.y),float(cpos.z)])
            BA=B-A
            CA=C-A
            u1=BA/np.linalg.norm(BA)
            uCA=CA/np.linalg.norm(CA)
            u2=uCA-(np.dot(uCA,u1)/np.dot(u1,u1))*u1
            dotprod=np.dot(u1,u2)
            atomindextoprojcoords={}
            for atom in rdkitmol.GetAtoms():
                atomidx=atom.GetIdx()
                pos = rdkitmol.GetConformer().GetAtomPosition(atomidx) 
                r=np.array([float(pos.x),float(pos.y),float(pos.z)])-A # shift coordinate system to be on top of A
                a1=np.dot(u1,r)
                a2=np.dot(u2,r)
                atomindextoprojcoords[atomidx]=[a1,a2]
            perpindicularsum=0 # make sure all projected distances in perpinducular direction are as close to 0 as possible
            for atomindex,coords in atomindextoprojcoords.items():
                    perpindicularsum+=coords[1]
            perpindicularsum=np.abs(perpindicularsum)
            return perpindicularsum

        def GenerateConformers(self,rdkitmol,mol, maxconfs=1000):
            AllChem.EmbedMolecule(rdkitmol,randomSeed=10) 
            conformer=rdkitmol.GetConformer(0)
            m2=copy.deepcopy(rdkitmol)
            mp = AllChem.MMFFGetMoleculeProperties(m2)
            mp.SetMMFFOopTerm(False)
            ffm = AllChem.MMFFGetMoleculeForceField(m2, mp)
            confid=0
            bondcutoff=2
            ndim=1
            rotbnds=[]
            frobnds=[]
            distmat=Chem.rdmolops.GetDistanceMatrix(rdkitmol) # bond distance matrix
            rotbnd=self.onlyrotbndslist[0]
            t2idx=rotbnd[0]
            t3idx=rotbnd[1]
            bond=rdkitmol.GetBondBetweenAtoms(t2idx-1,t3idx-1)
            bondidx=bond.GetIdx()
            bondrow=distmat[bondidx]
            for bondidx in range(len(bondrow)):
                bonddist=bondrow[bondidx]
                bond=bondrow[bondidx]
                try:
                    bond=rdkitmol.GetBondWithIdx(bondidx)
                except:
                    continue
                bgnatom=bond.GetBeginAtom()
                endatom=bond.GetEndAtom()
                bgnatomidx=bgnatom.GetIdx()
                endatomidx=endatom.GetIdx()
                key='%d %d' % (bgnatomidx+1,endatomidx+1)
                revkey='%d %d' % (endatomidx+1,bgnatomidx+1)
                found=False
                if key in self.rotbndlist.keys():
                    thekey=key
                    found=True
                elif revkey in self.rotbndlist.keys():
                    thekey=revkey
                    found=True
                if found==True:
                    if bonddist>=bondcutoff:
                        frobnds.append(thekey)
                    else:
                        rotbnds.append(thekey)
                
            key='%d %d' % (t2idx,t3idx)
            revkey='%d %d' % (t3idx,t2idx)
            angles=range(0,360,1)
            courseangles=range(0,370,10)
            phaselists=[]
            rotkeytoindex={}
            count=0
            for rotbndkeyidx in range(len(self.rotbndlist.keys())): # first one is b-c that will derive torsion
                rotkey=list(self.rotbndlist.keys())[rotbndkeyidx]
                if count>=ndim:
                    break
                if rotbndkeyidx==0:
                    phaselists.append(angles)
                    rotkeytoindex[rotkey]=rotbndkeyidx
                    count+=1
                else:
                    if rotkey in rotbnds:
                        phaselists.append(courseangles)
                        rotkeytoindex[rotkey]=count
                        count+=1
            phaselists=phaselists[:ndim]
            phaselist=(list(product(*phaselists)))
            if len(phaselist) > maxconfs:
                phaselist = random.sample(phaselist, k=maxconfs)
            phaselist = np.asarray(phaselist)
            t2=mol.GetAtom(t2idx)
            t3=mol.GetAtom(t3idx)
            t1,t4 = torgen.find_tor_restraint_idx(self,mol,t2,t3)
            torset=tuple([[t1.GetIdx(),t2.GetIdx(),t3.GetIdx(),t4.GetIdx()]])
            phaseangles=[0]
            rotbndtorescount={}
            maxrotbnds=1
            restlist=[]
            variabletorlist=[]
            newtorset=[]
            torsiontophaseangle={}
            torsiontomaintor={}
            rottors,rotbndtorescount,restlist,rotphases,torsiontophaseangle,torsiontomaintor=torgen.RotatableBondRestraints(self,torset,variabletorlist,rotbndtorescount,maxrotbnds,mol,restlist,phaseangles,torsiontophaseangle,torsiontomaintor,crash=False)
            frotors,rotbndtorescount,restlist,torsiontophaseangle,torsiontomaintor=torgen.FrozenBondRestraints(self,torset,variabletorlist,rotbndtorescount,maxrotbnds,mol,restlist,phaseangles,torsiontophaseangle,torsiontomaintor)
            tortophase={}
            tortoangle={}
            alltors=[]
            alltors.extend(rottors)
            alltors.extend(frotors)
            for tor in alltors:
                tor=[i-1 for i in tor]
                a=rdMolTransforms.GetDihedralDeg(conformer, tor[0],tor[1],tor[2],tor[3])
                if a<0:
                    a+=360
                tortoangle[tuple(tor)]=a
            firsttor=rottors[0]
            firsttor=[i-1 for i in firsttor] 
            firstangle=tortoangle[tuple(firsttor)]
            for i in range(len(rottors)):
                rottor=rottors[i]
                rottor=[j-1 for j in rottor]
                if i!=0:
                    a=tortoangle[tuple(rottor)]
                    phase=firstangle-a
                else:
                    phase=0
                tortophase[tuple(rottor)]=phase
            for angletup in phaselist:
                    
                confid+=1
                ff2 = AllChem.MMFFGetMoleculeForceField(m2, mp)
                for j in range(len(rottors)):
                    rottor=rottors[j]
                    rottor=[i-1 for i in rottor]
                    phase=tortophase[tuple(rottor)]
                    rotkey='%d %d' % (rottor[1]+1,rottor[2]+1)
                    if rotkey in rotkeytoindex.keys():
                        angleidx=rotkeytoindex[rotkey]
                        angle=angletup[angleidx]
                        ff2.MMFFAddTorsionConstraint(rottor[0],rottor[1],rottor[2],rottor[3], False, angle+phase - .1, angle+phase + .1, 100.0)

                for j in range(len(frotors)):
                    frotor=frotors[j]
                    frotor=[i-1 for i in frotor]
                    currentangle=tortoangle[tuple(frotor)]
                    rotkey='%d %d' % (frotor[1]+1,frotor[2]+1)
                    if rotkey in rotkeytoindex.keys():
                        angleidx=rotkeytoindex[rotkey]
                        angle=angletup[angleidx]
                    else:
                        angle=currentangle

                    ff2.MMFFAddTorsionConstraint(frotor[0],frotor[1],frotor[2],frotor[3], False, angle - .1, angle + .1, 100.0)



                ff2.Minimize()
                xyz=ff2.Positions()
                new_conf = Chem.Conformer(rdkitmol.GetNumAtoms())
                for i in range(rdkitmol.GetNumAtoms()):
                    new_conf.SetAtomPosition(i, (m2.GetConformer(-1).GetAtomPosition(i)))
                new_conf.SetId(confid)
                rdkitmol.AddConformer(new_conf) 

            return rdkitmol

        def GenerateMaxSymmetryConformer(self,rdkitmol,mol):
            rdkitmol=self.GenerateConformers(rdkitmol,mol) 
            confs=rdkitmol.GetConformers()
            symmetrictoconf={}
            for i in range(len(confs)):
                conf=confs[i]
                name="conftest"+".mol"
                rdmolfiles.MolToMolFile(rdkitmol,name,confId=i)
                mol=rdmolfiles.MolFromMolFile(name,removeHs=False)
                prin=False
                symmetric=self.DetermineSymmetryMetric(mol,prin)
                symmetrictoconf[symmetric]=i
            minsymmetric=min(symmetrictoconf.keys())
            confindex=symmetrictoconf[minsymmetric]
            name="bestconf.mol"
            rdmolfiles.MolToMolFile(rdkitmol,name,confId=confindex)
            confslist=[confindex]
            indextocoordslist=[]
            for confindex in confslist:
                indextocoordinates={}
                rdmolfiles.MolToMolFile(rdkitmol,name,confId=confindex)
                mol=rdmolfiles.MolFromMolFile(name,removeHs=False)
                for i in range(len(mol.GetAtoms())):
                    pos = mol.GetConformer().GetAtomPosition(i) 
                    vec=np.array([float(pos.x),float(pos.y),float(pos.z)])
                    indextocoordinates[i]=vec
                indextocoordslist.append(indextocoordinates)
            return indextocoordslist


        def GenerateExtendedConformer(self):
            """
            Intent: When deriving multipoles, dont want to have folded conformation where densities from other parts of the molecule overlap. So instead try to find most extended conformor and then will freeze dihedrals during QM opt.
            Input: Rdkit mol and openbabel mol
            Output: Dictionary of atomic index to coordinates.
            Referenced By: GenerateParameters
            Description: 
            1. Call EmbedMultipleConfs to add many conformors to rdkit object.
            2. Grab energy of all confomrors with MMFF94 force field
            3. Iterate over all conformors and determine max distance in molecule and energy of molecule and save in dictionaries.
            4. Make the first confomor the maxdistance conformor, 
            5. If user input wants multiple confomors for multipole parameterization then use the next lowest energy confomors (not necearrily maximaly extended conformation) as well.
            6. Save coordinates of conformations in dictionary for later use in QM optimzation  
            """
            pythonpath=self.which('python')
            head,tail=os.path.split(pythonpath)
            pythonpath=Path(head) 
            envdir=pythonpath.parent.absolute()
            envpath=Path(envdir)
            allenvs=envpath.parent.absolute()
            xtbenvpath=os.path.join(allenvs,self.xtbenvname)
            pythonpath=os.path.join(xtbenvpath,'bin')
            xtbpath=os.path.join(pythonpath,'xtb')
            cmdstr = f"python \"{os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'lConformerGenerator.py')}\" -i {self.molstructfname} -p {xtbpath}"
            self.WriteToLog('Calling: '+cmdstr) 
            os.system(cmdstr)
            name = "conftest.mol" 
            indextocoordslist=[]
            indextocoordinates={}
            mol=rdmolfiles.MolFromMolFile(name,removeHs=False)
            for i in range(len(mol.GetAtoms())):
                pos = mol.GetConformer().GetAtomPosition(i) 
                vec=np.array([float(pos.x),float(pos.y),float(pos.z)])
                indextocoordinates[i]=vec
            indextocoordslist.append(indextocoordinates)
            # STEP 6
            if len(self.espextraconflist)!=0:
                curdir=os.getcwd()
                os.chdir('..')
                for filename in self.espextraconflist:
                    shutil.copy(filename,os.path.join('Temp',filename))
                os.chdir(curdir)
                for filename in self.espextraconflist:
                    indextocoordinates=self.ReadCoordinates(filename)
                    indextocoordslist.append(indextocoordinates)

            return indextocoordslist


        def ReadCoordinates(self,filename):
            """
            Intent: Read coordinates from .mol file generated by GenerateExtendedConformer, then save in dictionary to later be used by QM geometry optimization
            Input: .mol filename
            Output: Dictionary of atomic index to coordinates
            Referenced By: GenerateExtendedConformer
            Description:
            1. Read in .mol file to openbabel mol object
            2. Iterate over mol object atoms
            3. Save atom index and coordinates in dictionary
            """
            # STEP 1
            indextocoordinates={}
            obConversion = openbabel.OBConversion()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(filename) # self.molstructfname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, filename) # self.molstructfname)
            # STEP 2
            iteratombab = openbabel.OBMolAtomIter(mol)
            for atm in iteratombab:
                atmindex=atm.GetIdx()
                coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
                # STEP 3
                indextocoordinates[atmindex-1]=coords

            return indextocoordinates


        def FindLongestDistanceInMolecule(self,mol):
            """
            Intent: For the given conformation, compute the longest pairwise atom-atom distance
            Input: Rdkit mol object
            Output: Maximum pairwise distance
            Referenced By: GenerateExtendedConformer
            Description: 
            1. Store all atom positions into an array
            2. Generate all pairwise combinations of coordinate vectors
            3. Iterate over all pairwise combinations of coordinate vectors and compute the distance, then store in array
            4. Determine the maximum pairwise distance from array of pairwise distances
            """
            # STEP 1
            veclist=[]
            for i in range(len(mol.GetAtoms())):
                pos = mol.GetConformer().GetAtomPosition(i) 
                vec=np.array([float(pos.x),float(pos.y),float(pos.z)])
                veclist.append(vec)
            # STEP 2
            pairs=list(itertools.combinations(veclist, 2))
            distlist=[]
            # STEP 3
            for pairidx in range(len(pairs)):
                pair=pairs[pairidx]
                dist=np.linalg.norm(np.array(pair[0])-np.array(pair[1]))
                distlist.append(dist)
            # STEP 4
            maxdist=np.amax(np.array(distlist))
            return maxdist
        
        def SetDefaultCoordinatesBabel(self,mol,indextocoordinates):
            """
            Intent: Need to update openbabel mol object coordinates with coordinates of the extended conformation generated.
            Input: openbabel mol object, dictionary of atomic index to coordinates
            Output: Updated mol object
            Referenced By: GenerateParameters, GenerateListOfMols
            Description: 
            1. Iterate over atom index to coordinates dictionary
            2. Save coordinates from dictionary to mol atom oject
            """
            # STEP 1
            for index,coords in indextocoordinates.items():
                atom=mol.GetAtom(index+1)
                x=coords[0]
                y=coords[1]
                z=coords[2]
                # STEP 2
                atom.SetVector(x,y,z)
            return mol

        def CheckForConcentratedFormalCharges(self,m,atomindextoformalcharge):
            """
            Intent: Often if a zwitterion is present or highly concentrated charge like phosphate group is present with nearby hydrogens, the hydrogens can migrate to the highly charged regions during QM optimization. So need to determine when there is concentrated formal charge and then turn on PCM for qm geometry optimization. 
            Input: Rdkit mol object, dictionary of atomic index to formal charge 
            Output: 0, no charge; 1, at least one charged atom; 2, postively and negatively charged atoms are close to each other
            Referenced By: GenerateParameters 
            Description: 
            1. If the molecule contains hydrogens, then continue, else dont use PCM
            2. Iterate over dictionary of atomic index to formal charges and save atomic indices that have formal charges to an array.
            3. check the bond distance between charged atoms

            """
            distmat=Chem.rdmolops.GetDistanceMatrix(m) # bond distance matrix
            positive_indices = []
            negative_indices = []
            charge_state = 0
            BOND_DIST_UPPER = 4 # upper bound for bond distance for determining zwitterion
            BOND_DIST_LOWER = 2
            # STEP 1
            if self.hashyd==True:
                # STEP 2
                for atomindex,chg in atomindextoformalcharge.items():
                    if chg > 0:
                        positive_indices.append(atomindex)
                        charge_state = 1
                    elif chg < 0:
                        negative_indices.append(atomindex)
                        charge_state = 1
                # STEP 3
                for atomi in positive_indices:
                    for atomj in negative_indices:
                        if distmat[atomi][atomj] <= BOND_DIST_UPPER and \
                            distmat[atomi][atomj] >= BOND_DIST_LOWER:
                            charge_state = 2
                            return charge_state
            return charge_state 
                    
                        
        def DeleteAllNonQMFiles(self,folderpath=None):
            """
            Intent: Sometimes if there is issue with non-qm (tinker files etc) then easier to autodelete instead of user manually needing to know which files to delete before starting parameterization, so this will auto delete non-qm log files.
            Input: -  
            Output: -
            Referenced By: GenerateParameters 
            Description: 
            1. In current directory call  DeleteNonQMFiles
            2. If qm-torsion folder exists, navigate to that folder and call DeleteNonQMFiles
            3. If vdw folder exists, navigate to that folder and call DeleteNonQMFiles
            """
            # STEP 1
            tempdir=os.getcwd()
            if folderpath!=None:
                os.chdir(folderpath)
            self.DeleteNonQMFiles(os.getcwd()) 
            # STEP 2
            if os.path.isdir('qm-torsion'):
                os.chdir('qm-torsion')
                self.DeleteNonQMFiles(os.getcwd()) 
                os.chdir('..')
            # STEP 3
            if os.path.isdir('vdw'):
                os.chdir('vdw')
                self.DeleteNonQMFiles(os.getcwd()) 
                os.chdir('..')

            os.chdir(tempdir)
           
        def DeleteNonQMFiles(self,directory):
            """
            Intent: Delete all non-QM log files in directory
            Input: directory
            Output: - 
            Referenced By: DeleteAllNonQMFiles
            Description:
            1. Change to new directory
            2. If user inputs certain files in poltype.ini, dont delete them by appending to array of files to not delete.
            3. Iterate over files and target files that are not related to qm log files or their associated outputs (.chk,.dat). Dont delete output AMOEBA, xtb or ANI output files also. Save files to be deleted to an array.
            4. Iterate over array of files to delete and delete them. 
            """
            # STEP 1
            tempdir=os.getcwd()
            os.chdir(directory)
            deletearray=[]
            files=os.listdir()
            filestonotdelete=[]
            # STEP 2
            if self.indextotypefile!=None:
                filestonotdelete.append(self.indextotypefile)
            if self.inputkeyfile!=None:
                filestonotdelete.append(self.inputkeyfile)
            # STEP 3
            for f in files:
                if not os.path.isdir(f) and 'nohup' not in f and f[0]!='.' and f!='parentvdw.key':
                    fsplit=f.split('.')
                    if len(fsplit)>1:
                        end=fsplit[1]
                        if 'AMOEBA-opt' not in f and 'ANI-opt' not in f and 'xtb-opt' not in f and 'FENNIX-opt' not in f and 'log' not in end and 'sdf' not in end and 'ini' not in end and 'chk' not in end and 'dat' not in end and 'mol' not in end and 'txt' not in end and f not in filestonotdelete:
                            deletearray.append(f)
            # STEP 4
            for f in deletearray:
                os.remove(f)

            os.chdir(tempdir) 

        
        def CheckIfAtomsAreAllowed(self,m):
            """
            Intent: AMOEBA only supports certain atom types (not all elements possible have atom types), so check if the input molecule contains any elements outside current types allowed and if it doesnt exist, then crash the program.
            Input: Rdkit mol object
            Output: -
            Referenced By: GenerateParameters 
            Description: 
            1. Hard code list of allowed elements.
            2. Iterate over mol object atoms and get the atomic number.
            3. If the atomic number is not in list of allowed atomic numbers, then raise error and crash the program.
            """
            # STEP 1
            listofallowedatoms=[1,5,6,7,8,15,16,17,35,53,9]
            # STEP 2
            for atom in m.GetAtoms():
                atomicnum=atom.GetAtomicNum()
                # STEP 3
                if atomicnum not in listofallowedatoms:
                    raise ValueError('Element not allowed! '+str(atomicnum)) 


        def CheckIfAtomsAreAllowedANI(self,m):
            """
            Intent: ANI-2 doesnt allow certain elements (Phosphorous, Iodine, Bromine etc...), so if the user sets the torsion OPT methods to use ANI-2, ensure that the input molecule has allowed elements. If not, then crash the program.
            Input: Rdkit mol object
            Output: -
            Referenced By: GenerateParameters 
            Description: 
            1. If user specifies to use ANI-2 then continue, else return
            2. Hard code list of allowed elements for ANI-2
            3. Iterate over mol object atoms
            4. Get the atomic number and if doesnt exist in list of allowed atomic numbers, then raise error and crash the program. 
            """
            # STEP 1
            if self.toroptmethod=='ANI' or self.torspmethod=='ANI':
                # STEP 2
                listofallowedatoms=[1,6,7,8,9,17,16]
                # STEP 3
                for atom in m.GetAtoms():
                    # STEP 4
                    atomicnum=atom.GetAtomicNum()
                    if atomicnum not in listofallowedatoms:
                        raise ValueError('Element not allowed for ANI! '+str(atomicnum)) 


        def CheckIfFENNIXCanBeUsed(self,m):
            """
            Currently the NN model in FENNIX can only be used on elements: H,C,N,O,F,S,Cl and neutral molecules
            Modify the method to default ones if this is the case
            """
            atomsAreCovered = True
            if self.toroptmethod=='FENNIX' or self.torspmethod=='FENNIX':
              # FENNIX requires total charge must be ZERO 
              isNeutral = True
              chg = Chem.GetFormalCharge(m) 
              if chg != 0:
                isNeutral = False
                self.WriteToLog(f'Total Charge of the Molecule is not Zero: {chg}') 
                raise ValueError('Total Charge of the Molecule is not Zero') 
              
              listofallowedatoms=[1,6,7,8,9,17,16]
              for atom in m.GetAtoms():
                atomicnum=atom.GetAtomicNum()
                if atomicnum not in listofallowedatoms:
                  self.WriteToLog('Element not allowed for FENNIX! '+str(atomicnum)) 
                  atomsAreCovered = False
            
            if not atomsAreCovered:
              raise ValueError('Not All Elements are allowed for FENNIX! ') 

            

        def CheckForHydrogens(self,m):
            """
            Intent: Sometimes user may download a structure that does not contain hydrogens, so check for hydrogens and print warning if any hydrogens are missing from the input structure. 
            Input: Rdkit mol object
            Output: -
            Referenced By: GenerateParameters
            Description:
            1. Assume there are no hydrogens
            2. Iterate over mol object atoms
            3. If atomic number is a hydrogen then assumption is wrong
            4. If no hydrogens were found, then print warning in log file and to standard output 
            """
            # STEP 1
            self.hashyd=False
            # STEP 2
            for atom in m.GetAtoms():
                atomicnum=atom.GetAtomicNum()
                # STEP 3
                if atomicnum==1:
                    self.hashyd=True
            # STEP 4
            if self.hashyd==False:
                string='No hydrogens detected in input file!'
                warnings.warn(string)
                self.WriteToLog(string)


        def GrabAtomicSymbols(self,rdkitmol):
            """
            Intent: Save dictionary of atomic index to atomic symbol for later use when generating Gaussian and Psi4 input files
            Input: Rdkit mol object
            Output: dictionary of atomic index to atomic symbol
            Referenced By: GenerateParameters 
            Description:
            1. Iterate over mol object
            2. Get atomic number
            3. Use pyastronomy lib to convert atomic number to atomic symbol and then save in dictionary 
            """
            indextoatomicsymbol={}
            an = pyasl.AtomicNo()
            # STEP 1
            for atm in rdkitmol.GetAtoms():
                # STEP 2
                atmnum=atm.GetAtomicNum()
                atmidx=atm.GetIdx()
                # STEP 3
                sym=an.getElSymbol(atmnum)
                indextoatomicsymbol[atmidx+1]=sym

            return indextoatomicsymbol


        def CheckIfInputIsTinkerXYZ(self,molstructfname):
            """
            Intent: Determine if the input structure is tinker XYZ or not
            Input: Input structure file
            Output: boolean specifying if the structure is tinker XYZ format or not
            Referenced By: ConvertInputStructureToSDFFormat
            Description: 
            1. Assume that the input structure is not tinker XYZ
            2. Iterate over lines of file
            3. If there are more than 4 entries in line (cartesian XYZ has 4), then it must be a tinker XYZ file. 
            """
            # STEP 1
            istinkxyz=False
            temp=open(molstructfname,'r')
            results=temp.readlines()
            temp.close()
            # STEP 2
            for line in results:
                linesplit=line.split()
                if len(linesplit)>1:
                    # STEP 3
                    if len(linesplit)>4:
                        istinkxyz=True

            return istinkxyz


        def GenerateIndexToTypeFile(self,indextotype):
            """
            Intent: If user gives tinker XYZ as input, need to save type information to a file and give it to poltype to read in and use instead of default symmetry type detection. 
            Input: Dictionary of atom index to type numbers
            Output: Filename to be processed by poltype later
            Referenced By: ConvertInputStructureToSDFFormat
            Description: 
            1. Create file handle for filename
            2. Iterate over dictionary of atom index to type number
            3. Write out atomic index and type number to the filename 
            """
            filename='indextotype.txt'
            # STEP 1
            temp=open(filename,'w')
            # STEP 2
            for index,typenum in indextotype.items():
                # STEP 3
                temp.write(str(index)+' '+str(typenum)+'\n')
            temp.close()
            return filename


        def GrabIndexToType(self,molstrucfname):
            """
            Intent: Read tinker XYZ file and save atomic index to type information into a dictionary
            Input: Tinker XYZ file
            Output: Dictionary of atomic index to type number
            Referenced By: ConvertInputStructureToSDFFormat
            Description: 
            1. Iterate over lines of tinker XYZ file
            2. Grab atomic index and type number
            3. Save in dictionary
            """
            indextotype={}
            temp=open(molstrucfname,'r')
            results=temp.readlines()
            temp.close()
            # STEP 1
            for line in results:
                linesplit=line.split()
                if len(linesplit)>=5:
                    # STEP 2
                    index=int(linesplit[0])
                    typenum=int(linesplit[5])
                    # STEP 3
                    indextotype[index]=typenum

            return indextotype



        def ConvertInputStructureToSDFFormat(self,molstructfname):
            """
            Intent: Users may give many different possible file formats for the input structure.
            Input: Input structure file
            Output: Converted structure file
            Referenced By: GenerateParameters
            Description: 
            1. If not .mol or .sdf file format (.mol is defauled for fragment job inputs)
            2. Check if input is tinker XYZ file
            3. If file is tinker XYZ, then need to convert to cartesian XYZ and save the type information for later use to ensure that poltype doesnt use default symmetry typing. Otherwise read in structure file to openbabel mol object.
            4. Convert mol object to output SDF file format. 
            """
            obConversion = openbabel.OBConversion()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(molstructfname)
            split=molstructfname.split('.')
            ext=split[-1]
            # STEP 1
            if ext!='sdf' and ext!='mol':
                # STEP 2
                istinkerxyz=False
                if ext=='xyz':
                    istinkerxyz=self.CheckIfInputIsTinkerXYZ(molstructfname)

                obConversion.SetInFormat(ext)
                # STEP 3
                if istinkerxyz==False:
                    obConversion.ReadFile(mol, molstructfname)
                else:
                    newname=self.ConvertTinkerXYZToCartesianXYZ(molstructfname)
                    obConversion.ReadFile(mol, newname)
                    indextotype=self.GrabIndexToType(molstructfname)
                    filename=self.GenerateIndexToTypeFile(indextotype) 
                    self.indextotypefile=filename
                # STEP 4
                obConversion.SetOutFormat('sdf')
                molstructfname=molstructfname.replace('.'+ext,'.sdf')
                obConversion.WriteFile(mol,molstructfname)


            return molstructfname


        def GenerateEntireTinkerXYZ(self,atmindextocoordinates,atmindextotypenum,atmindextoconnectivity,atmindextoelement,filename):
            """
            Intent: If user chooses to skip multipole derivation (which requires QM optimization and running GDMA, poledit, which produces the first tinker XYZ file), then need to generate a tinker XYZ file using input mol object coordinates
            Input: Dictionary of atomic index to coordinates, dictionary of atomic index to type number, dictionary of atomic index to connectivity, dictionary of atomic index to element, filename of tinker XYZ to be made
            Output: - 
            Referenced By: GenerateParameters
            Description: 
            1. Write the total atom number at the top of the tinker XYZ file
            2. Iterate over dictionary of atomic index to coordinates, grab type number, connectivity and element symbol from other dictionaries.
            3. Construct line in format of tinker XYZ and write line to file. 
            """
            temp=open(filename,'w')
            # STEP 1
            temp.write(str(len(atmindextocoordinates.keys()))+'\n')
            for index,coords in atmindextocoordinates.items():
                # STEP 2
                typenum=atmindextotypenum[index]
                conns=atmindextoconnectivity[index]
                element=atmindextoelement[index]
                x=coords[0]
                y=coords[1]
                z=coords[2]
                # STEP 3
                newline='    '+str(index)+'  '+element+'     '+str(x)+'   '+str(y)+'   '+str(z)+'    '+str(typenum)+'     '
                for con in conns:
                    newline+=str(con)+'     '
                newline+='\n'
                temp.write(newline)

            temp.close()


        def GenerateParameters(self):
            """
            Intent: Assign parameters from database and/or derive them from ab initio QM computations
            Input: Self object
            Output: - 
            Referenced By: main  
            Description: 
            1. Convert input structure to SDF format
            2. If not a fragment job, then make a folder called Temp, copy input files into Temp and perform poltype calculations within this folder. This way when return final.xyz and final.key to the user it is a lot less messy.
            3. If user specifies to delete all non QM files, then delete them.
            4. Read in SDF file to openbabel mol object
            5. Initialize poltype.log file handle
            6. Convert SDF to MOL file, so rdkit mol object can be generated (has different functions available than openbabel)
            7. Check if input atoms are allowed (by current types in AMOEBA), if ANI-2 is used, check if input atoms are allowed. If there are no hydrogens in the molecule just print a warning to log file and standard output.
            8. Determine formal charges and total charge from input SDF file, in some cases such as nitrogen or carbon unprotonated, then will automatically add hydrogen. 
            9. Check for multiple fragments in input file, if detected, then crash program.
            10. Determine if need to use PCM based on if any hydrogens and concentrated formal charges. 
            11. Initialize array containing atomic indices to coordinates for each conformation that will be used for multipole fitting (can fit many conformations)
            12. Determine if the whole molecule is one giant ring (like host etc), then dont want to generate extended conformations (cause issues in rdkit)
            13. If molecule input is 2D coordinates, generate 3D coordinates
            14. Compute symmetry type numbers for input molecule. 
            
            15. Find partial double bonds in molecule, for later use in deciding which torsions to transfer and fit from database, also generate list of torsions for rotatable bonds for generating Max symmetry conformer. 
            16. Generate extended conformation via rdkit and also if user specifies, add more conformations for multipole fitting. Generate max symmetry conformer for fragment job. 
            17. Update rdkit object with extended conformation coordinates
            18. If scratch directories dont exist, then make them.
            19. Generate ionization states and enumerate tautomers, if user only wants to generate protonation states, then exit. 
            20. Iodine requires def2 basis sets with ECP and ECP-MP2 analytical gradients are not implemented in psi4, so if Iodine is detected in molecule, then change default basis set from MP2 to wB97X-D
            21. The default basis set for torsion SP, doesnt have Br, so remove the + and use slightly lower basis set.
            22. Remove any cartesian XYZ files, as they interfere with tinker XYZ files (same extension) 
            23. For each conformation needed to be fit for multipoles, compute the QM geometry optimization. Then afterwords check to ensure no bonds were created/broken if so, then crash the program.
            24. If user is providing multipole and polarize parameter to poltype, then dont need to run QM geometry optimization so generate a fake geometry opt XYZ output file from input coordinates.
            25. If user wants to only do geometry optimzation then quit program
            26. Search database for parameters to be later appended to keyfile.
            27. Compute QM needed for GDMA to determine initial multipole parameters
            28. Run GDMA from output of QM in previous step
            29. Grab polarize parameters from database, to be used as input for poledit to generate first tinker XYZ and key files
            30. Generate poledit input file, using polarize parameters as inputs. Multipole symmetry frame detection is used also. 
            31. Run poledit to generate first tinker XYZ and key file
            32. Append header to keyfile, replace empty spaces in multipole frame definition with zeroes, alzo ensure x components of multipole that need to be zeroed out are zeroed out prior to averaging then electrostatic potential fitting.
            33. If the user requested multiple conformations for multipole fitting, then use the first generated tinker XYZ as a template and generate the other tinker XYZs using optimized coordinates from each conformation.
            34. Compute high level QM for generating electrostatic potential grid for multipole fitting later
            35. Generate the electrostatic potential grid using output from QM computed in last step
            36. Average multipole and polarize parameters via atom types 
            37. Add comments to key file for polarize parameters section
            38. Electrostatic potential fitting
            39. Remove header from output key file
            40. Append database parameters to keyfile
            41. Determine list of rotatable bonds that need to have parameters derived based on what was found from database matching, list will be used by fragmenter later. Determine how many points on dihedral surface need to sample based on number of torsions around rotatable bond that need to be derived and how many cosine terms used for fitting. 
            42. Turn off fragmenter if number of atoms < 25
            43. If there are missing vdw parameters and the user specifies to fit vdw parameters, then add relevant atom indices to be probed into array for fragmenter to handle later. 
            44. If user specifies to fit vdw parameters, then call the fragmenter and generate fragment jobs, then run the fragment poltype jobs. 
            45. If user specifies to fit vdw parameters, now parse the output key files from vdw fragmentjobs and update parent key with vdw parameters derived in fragment.
            46. If user specifies to do tor-tor, then add list of rotatble bonds to array for fragmenter to later fragment.
            47. If user specifies to refine non-aromatic ring torsions, then determine the torsions in non-aromatic rings that need to be puckered and add to array that fragmenter will use to fragment later. 
            48. Create and run torsion fragment poltype jobs
            49. Grab torsion parameters from fragment jobs and add to parent key file
            50. Sanity check, ensure torsion parameters are not zeroed out (indicates some fragmenter or parsing issue)
            51. Minimize tinker XYZ with final parameters and restrain same dihedrals restrained during QM geometry optimzation, compute RMSD to qm geometry optimized structure and ensure its not too high, otherwise crash program. 
            52. Compute the QM dipole MM dipole moment and compare them, ensure the difference in magnitude is not too different, otherwise crash program.   
            53. Check for any missing vdw or multipole parameters in final key file, tinker will not crash if missing vdw parameters.
            54. Log how much time torsion opt and SP jobs took in poltype log file.
            55. If the fragmenter was not used but torsion parameters were derived, write out database.prm file of parameters in format needed to be added to database.
            56. Delete scratch files used during QM calculations
            57. Copy final XYZ and key files to top directory, where user submits jobs. Also copy torsion plots to OPENME folder.
            58. Apply any modifications or FLAG to the final.key file 
            """
            # STEP 1 
            self.molstructfname=self.ConvertInputStructureToSDFFormat(self.molstructfname)
            # STEP 2
            if self.isfragjob==False:
                foldername='Temp'
                if not os.path.exists(foldername):
                    os.mkdir(foldername)
                shutil.copy(self.molstructfname,os.path.join(foldername,self.molstructfname))
                shutil.copy('poltype.ini',os.path.join(foldername,'poltype.ini'))

                if self.indextotypefile!=None:
                    shutil.copy(self.indextotypefile,os.path.join(foldername,self.indextotypefile))
                if self.indextompoleframefile!=None:
                    shutil.copy(self.indextompoleframefile,os.path.join(foldername,self.indextompoleframefile))
                if self.inputkeyfile!=None:
                    shutil.copy(self.inputkeyfile,os.path.join(foldername,self.inputkeyfile))

                os.chdir(foldername)
            self.startdir=os.getcwd()
            self.totalcharge=None
            # STEP 3
            if self.deleteallnonqmfiles==True:
                self.DeleteAllNonQMFiles()
            # STEP 4
            obConversion = openbabel.OBConversion()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(self.molstructfname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, self.molstructfname)
            self.atomnum=mol.NumAtoms()
            # STEP 5 
            self.logfh = open(self.logfname,"w",buffering=1)
            # STEP 6
            self.molstructfnamemol=self.molstructfname.replace('.sdf','.mol')
            if '.mol' not in self.molstructfname: 
                obConversion.SetOutFormat('mol')
                obConversion.WriteFile(mol,self.molstructfnamemol)
            indextocoordinates=self.GrabIndexToCoordinates(mol)
            m=Chem.MolFromMolFile(self.molstructfnamemol,removeHs=False,sanitize=False)
            # STEP 7
            self.CheckIfAtomsAreAllowed(m)
            self.CheckIfAtomsAreAllowedANI(m)
            self.CheckIfFENNIXCanBeUsed(m)
            self.CheckForHydrogens(m)
            # STEP 8
            m,atomindextoformalcharge=self.CheckInputCharge(m,verbose=True)
            if self.allowradicals==True:
                self.dontfrag=True # Psi4 doesnt allow UHF and properties (like compute WBO) for fragmenter, so need to turn of fragmenter if radical detected
            m.UpdatePropertyCache()
            if self.addhydrogentocharged==True and self.isfragjob==False:
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m)
            self.indextoatomicsymbol=self.GrabAtomicSymbols(m)
            Chem.SanitizeMol(m)
            # STEP 9
            smarts=rdmolfiles.MolToSmarts(m)
            if '.' in smarts:
                raise ValueError('Multiple fragments detectected in input molecule')
            # STEP 10
            charge_state=self.CheckForConcentratedFormalCharges(m,atomindextoformalcharge)
            self.pcm = (charge_state==2) and (self.pcm_auto)

            self.temptorsppcm=self.torsppcm
            self.temptoroptpcm=self.toroptpcm
            
            if self.use_gauPCM==True:
                self.use_gausoptonly=False
                self.use_gaus=True
                self.SanitizeAllQMMethods()

            # STEP 11 
            indextocoordslist=[indextocoordinates]
            # STEP 12
            iswholemoleculering=self.CheckIfWholeMoleculeIsARing(mol)
            if iswholemoleculering==True:
                self.generateextendedconf=False
            # STEP 13
            mol,m=self.CheckIsInput2D(mol,obConversion,m)

            # STEP 14
            self.canonicallabel = [ 0 ] * mol.NumAtoms()
            self.localframe1 = [ 0 ] * mol.NumAtoms()
            self.localframe2 = [ 0 ] * mol.NumAtoms()
            self.WriteToLog("Atom Type Classification")
            self.idxtosymclass,self.symmetryclass=symm.gen_canonicallabels(self,mol,None,self.usesymtypes,True)
            # STEP 15
            torgen.FindPartialDoubleBonds(self,m,mol) # need to find something to hardcode transfer for partial double amide/acid, currently will derive torsion parameters if doesnt find "good" match in torsion database
            torgen.get_all_torsions(self,mol)
            (torlist, self.rotbndlist,nonaroringtorlist,self.nonrotbndlist) = torgen.get_torlist(self,mol,[],[],allmissing=True) # need to call this to get self.rotbndlist to generate restraints for N-dimensional scan in GenerateMaxSymmetryConformer
            # STEP 16
            if self.firstoptfinished==False and self.isfragjob==False and self.generateextendedconf==True:
                indextocoordslist=self.GenerateExtendedConformer()
                indextocoordinates=indextocoordslist[0]
            if self.isfragjob==True and self.generate_symm_frag_conf and len(self.onlyrotbndslist)!=0:
                self.WriteToLog("Generating Max Symmetry Conformer")
                indextocoordslist=self.GenerateMaxSymmetryConformer(m,mol)
                indextocoordinates=indextocoordslist[0]
            # STEP 17
            Chem.GetSymmSSSR(m)
            m.GetRingInfo().NumRings() 
            try:
                m=self.AddInputCoordinatesAsDefaultConformer(m,indextocoordinates)
            except:
                pass
            if self.generateextendedconf==True:
                rdmolfiles.MolToMolFile(m,'extendedconf.mol')
            mol=self.SetDefaultCoordinatesBabel(mol,indextocoordinates)
            # STEP 18 
            if not os.path.exists(self.scrtmpdirpsi4):
                os.mkdir(self.scrtmpdirpsi4)
            if not os.path.exists(self.scrtmpdirgau):
                os.mkdir(self.scrtmpdirgau)

            molist=self.GenerateListOfMols(mol,indextocoordslist)
            self.mol=mol
            self.rdkitmol=m
            self.mol.SetTotalCharge(self.totalcharge)
            # STEP 19
            self.GrabIonizationStates(m)
            self.GrabTautomers(m)
            if self.genprotstatesonly==True:
                sys.exit()
            # STEP 20
            if ('I ' in self.mol.GetSpacedFormula()):
                self.optmethod=self.SanitizeQMMethod(self.optmethod,True)
            # STEP 21
            if ('Br ' in self.mol.GetSpacedFormula()):
                self.torspbasisset=self.torspbasissethalogen
            
            # STEP 22
            self.RemoveCartesianXYZFiles()
            self.WriteToLog("Running on host: " + gethostname())
             
            
             

                        
            # STEP 23    
            if (self.writeoutpolarize==True and self.writeoutmultipole==True):
                optmolist,errorlist,torsionrestraintslist = opt.GeometryOPTWrapper(self,molist)
                optmol=optmolist[0]
                error=errorlist[0]
                torsionrestraints=torsionrestraintslist[0]
                for omolidx in range(len(optmolist)):
                    bgnmol=molist[omolidx]
                    omol=optmolist[omolidx]
                    bondtopoopt=torgen.GenerateBondTopology(self,omol)
                    bondtopoopt=[list(i) for i in bondtopoopt]
                    bondtopo=torgen.GenerateBondTopology(self,bgnmol)
                    bondtopo=[list(i) for i in bondtopo]
                    for bond in bondtopo:
                        if bond in bondtopoopt or bond[::-1] in bondtopoopt:
                            pass
                        else:
                            raise ValueError('Bond does not exist after optimization !'+str(bond))

                    for bond in bondtopoopt:
                        if bond in bondtopo or bond[::-1] in bondtopo:
                            pass
                        else:
                            raise ValueError('Bond created after optimization !'+str(bond))


                    


            else:
                # STEP 24
                optmol=mol
                tmpconv = openbabel.OBConversion()
                tmpconv.SetOutFormat('xyz')
                tmpconv.WriteFile(mol, self.logoptfname.replace('.log','.xyz'))
                atmindextocoordinates,atmindextoconnectivity,atmindextoelement=self.ExtractMOLInfo(mol)
                self.GenerateEntireTinkerXYZ(atmindextocoordinates,self.idxtosymclass,atmindextoconnectivity,atmindextoelement,self.xyzfname)
            # STEP 25
            if self.optonly==True:
                sys.exit()
            # STEP 26
            if not os.path.isfile(self.key4fname) or not os.path.isfile(self.torsionsmissingfilename) or not os.path.isfile(self.torsionprmguessfilename):
                self.WriteToLog('Searching Database')
                torsionprmstotransferinfo,torsionsmissing,classkeytotorsionparametersguess,tortorprmstotransferinfo,tortorsmissing=torsiondatabaseparser.GrabSmallMoleculeAMOEBAParameters(self,optmol,mol,m)
            if os.path.isfile(self.torsionsmissingfilename):
                # Read missing torsions
                torsionsmissing=torsiondatabaseparser.ReadTorsionList(self,self.torsionsmissingfilename)
            if os.path.isfile(self.torsionprmguessfilename):
                classkeytotorsionparametersguess=torsiondatabaseparser.ReadDictionaryFromFile(self,self.torsionprmguessfilename)
            if os.path.isfile(self.vdwmissingfilename):
                missingvdwatomindices=torsiondatabaseparser.ReadVdwList(self,self.vdwmissingfilename)
            if self.onlyvdwatomlist!=None:
                missingvdwatomindices=self.onlyvdwatomlist[:]
            if os.path.isfile(self.tortormissingfilename):
                tortorsmissing=torsiondatabaseparser.ReadTorTorList(self,self.tortormissingfilename)

            # STEP 27
            if (self.writeoutpolarize==True and self.writeoutmultipole==True):
                esp.SPForDMA(self,optmol,mol)
            

                # STEP 28 
                if not os.path.isfile(self.gdmafname):
                    mpole.run_gdma(self)
                # STEP 30
                lfzerox = [ False ] * mol.NumAtoms()
                if not os.path.isfile(self.peditinfile):
                    lfzerox=mpole.gen_peditinfile(self,mol)
            
                if (not os.path.isfile(self.xyzfname) or not os.path.isfile(self.keyfname)):
                    # STEP 31
                    cmdstr = self.peditexe + " 1 " + self.gdmafname + " < " + self.peditinfile
                    self.call_subsystem([cmdstr],True)
                    # Add header to the key file output by poledit
                    while not os.path.isfile(self.keyfnamefrompoledit):
                        time.sleep(1)
                        self.WriteToLog('Waiting for '+self.keyfnamefrompoledit)
                    # STEP 32 
                    mpole.prepend_keyfile(self,self.keyfnamefrompoledit,optmol)
                    mpole.SanitizeMultipoleFrames(self,self.keyfnamefrompoledit)
                    mpole.post_proc_localframes(self,self.keyfnamefrompoledit, lfzerox)
                    shutil.copy(self.keyfnamefrompoledit,self.keyfname)
                # STEP 33
                xyzfnamelist,keyfnamelist=self.GenerateDuplicateXYZsFromOPTs(self.xyzfname,self.keyfname,optmolist)
            # STEP 34
            if self.atomnum!=1 and (self.writeoutpolarize==True and self.writeoutmultipole==True): 
                 gridnamelist,espnamelist,fchknamelist,cubenamelist=esp.SPForESP(self,optmolist,molist,xyzfnamelist,keyfnamelist) 

            if self.qmonly:
                self.WriteToLog("poltype QM-only complete.")
                sys.exit(0)
        
                   
            
            if (self.writeoutpolarize==True and self.writeoutmultipole==True):
                if self.atomnum!=1:
 
                # STEP 35
                    if not os.path.isfile(self.key3fname):
                        potnamelist=esp.gen_esp_grid(self,optmol,gridnamelist,espnamelist,fchknamelist,cubenamelist)
                # STEP 36  
                if not os.path.isfile(self.key2fnamefromavg):
                    self.WriteToLog("Average Multipoles Via Symmetry")
                    mpole.AverageMultipoles(self,optmol)
                    # STEP 37
                    # the comments and polarize params do not match anymore
                    # so we comment out this line
                    #mpole.AddPolarizeCommentsToKey(self,self.key2fnamefromavg,polartypetotransferinfo)
                    if self.forcefield.upper() in ['APLUS', 'AMOEBA+', 'AMOEBAPLUS']: 
                      self.WriteToLog("Assign Charge Penetration Parameters using DatabaseParser")
                      ldatabaseparser.assign_chgpen_params(self) 
                # STEP 38
                fit=False
                if self.espfit and not os.path.isfile(self.key3fname) and self.atomnum!=1:
                    xyzfnamelist,keyfnamelist=self.GenerateDuplicateXYZsFromOPTs(self.xyzoutfile,self.key2fnamefromavg,optmolist)   
                    if self.scaleandfixdipole:
                      lmodifytinkerkey.modkey2_COO(self)
                      self.WriteToLog("Special Treatment for COO Functional Group ")
                    if self.scalebigmultipole:
                      lmodifytinkerkey.modkey2_scalempole(self)
                    if self.fragbigmultipole:
                      lmodifytinkerkey.modkey2_fragmpole(self)
                    self.WriteToLog("Electrostatic Potential Optimization")
                    combinedxyz,combinedpot=esp.ElectrostaticPotentialFitting(self,xyzfnamelist,keyfnamelist,potnamelist) 
                    shutil.copy(self.key3fnamefrompot,self.key3fname)
                    fit=True
                elif self.atomnum==1 or self.espfit==False:
                    shutil.copy(self.key2fnamefromavg, self.key3fname)
                # Remove header terms from the keyfile
                # STEP 39
                mpole.rm_esp_terms_keyfile(self,self.key3fname)
                if fit==True:
                    if self.atomnum!=1: 
                        esp.ElectrostaticPotentialComparison(self,combinedxyz,combinedpot)
            # STEP 40
            if not os.path.exists(self.key4fname):
                torsiondatabaseparser.appendtofile(self,self.key3fname,self.key4fname, torsionprmstotransferinfo,tortorprmstotransferinfo)
                # Write Zero parameters for special torsion
                ldatabaseparser.zero_special_torsions(self)
                if self.writeoutangle==True:
                    # assign valence term here using new script
                    # Chengwen Liu
                    # Feb 2024
                    self.WriteToLog('Assign bonded parameters using DatabaseParser')
                    ldatabaseparser.assign_bonded_params(self)
                    if self.forcefield.upper() in ['APLUS', 'AMOEBA+', 'AMOEBAPLUS']:
                      self.WriteToLog('Assign Charge flux, Charge transfer, and Van der Waals parameters using DatabaseParser')
                    else:
                      self.WriteToLog('Assign Van der Waals and GK parameters using DatabaseParser')
                    ldatabaseparser.assign_nonbonded_params(self) 
                    torsiondatabaseparser.StiffenZThenBisectorAngleConstants(self,self.key4fname)
                    torsiondatabaseparser.TestBondAngleEquilValues(self)
                self.AddIndicesToKey(self.key4fname)
                if self.databasematchonly==True:
                    sys.exit()
                # make copy of key file before patches
                shutil.copy(self.key4fname,self.key4bfname)
                if len(self.prmmodlist) > 0:
                    # Apply prmmod patches
                    self.WriteToLog('Apply prm patches\n%s'%('\n'.join(self.prmmodlist)))
                    modify_key(self.xyzoutfile, self.key4fname, self.key4fname, sdffile=self.molstructfname, inpfile=self.prmmodlist)
                # make copy of key file before patches
                shutil.copy(self.key4fname,self.key4bfname+"_2")
            
            # STEP 40-41
            # Try to match torsions with lAssignAMOEBAplusPRM.py
            # here it can be thought as a patch to the current torsion database
            # only the missing torsions are replaced in key4 and beyond
            # Chengwen Liu
            # June 2024

            tmpkey = self.key4fname
            tmpsdf = self.molstructfname
            tmpxyz = self.xyzoutfile 
            # Match for AMOEBA FF
            if self.forcefield.upper() == 'AMOEBA':
              cmd = f'python {self.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {tmpsdf} -potent TORSION'
              self.call_subsystem([cmd], True)  
            # Match for AMOEBAplus FF
            if self.forcefield.upper() in ['AMOEBA+', 'AMOEBAPLUS']:
              cmd = f'python {self.ldatabaseparserpath} -xyz {tmpxyz} -key {tmpkey} -sdf {tmpsdf} -potent TORSION+'
              self.call_subsystem([cmd], True)  
            # find torsions to be added
            torsions_toadd = {}
            matched_torlist = []
            if os.path.isfile(tmpkey + '_torsion'):
              tmplines = open(tmpkey + '_torsion').readlines()
              for line in tmplines:
                if (line[0] != '#') and (line.split()[0].upper() == 'TORSION'):
                  s = line.split()
                  tor = [int(ss) for ss in s[1:5]]
                  tor_r = [int(ss) for ss in s[4:0:-1]] 
                  torsions_toadd['-'.join([str(i) for i in tor])] = line
                  torsions_toadd['-'.join([str(i) for i in tor_r])] = line
                  matched_torlist.append(tor)
                  matched_torlist.append(tor_r)

            # remove matched torsion types from torsionsmissing.txt
            if matched_torlist != []:
              shutil.move(self.torsionsmissingfilename, self.torsionsmissingfilename + '.back')
              tmp = []
              for tor in torsionsmissing:
                if tor not in matched_torlist:
                  tmp.append(tor)
              torsionsmissing = tmp
              with open(self.torsionsmissingfilename, 'w') as f:
               for torsionmiss in torsionsmissing:
                f.write(str(torsionmiss) + '\n')

              # write the matched torsions to key4 
              keylines = open(tmpkey).readlines()
              shutil.move(tmpkey, tmpkey+'b')
              with open(tmpkey, 'w') as f:
                for keyline in keylines:
                  if ('torsion ' in keyline) and len(keyline.split()) > 5:
                    s = keyline.split()
                    tor_str = '-'.join(s[1:5])
                    tor_str_r = '-'.join(s[4:0:-1])
                    if tor_str in torsions_toadd.keys():
                      f.write(torsions_toadd[tor_str] + '\n')
                    elif tor_str_r in torsions_toadd.keys():
                      f.write(torsions_toadd[tor_str_r] + '\n')
                    else:
                      f.write(keyline)
                  else:
                    f.write(keyline)
            # STEP 41
            (torlist, self.rotbndlist,nonaroringtorlist,self.nonrotbndlist) = torgen.get_torlist(self,optmol,torsionsmissing,self.onlyrotbndslist)
            torlist,self.rotbndlist=torgen.RemoveDuplicateRotatableBondTypes(self,torlist) # this only happens in very symmetrical molecules
            
            # STEP 41-42
            # remove matched torsion indices from torlist
            atom2type = {}
            xyzlines = open(tmpxyz).readlines()
            for line in xyzlines[1:]:
              s = line.split()
              atom2type[int(s[0])] = s[5]
            
            tmptorlist = []
            for tor in torlist:
              (a1, a2, a3, a4) = tor
              t1 = atom2type[a1]
              t2 = atom2type[a2]
              t3 = atom2type[a3]
              t4 = atom2type[a4]
              comb = '-'.join([t1,t2,t3,t4])
              comb_r = '-'.join([t4,t3,t2,t1])
              if not ((comb in torsions_toadd.keys()) or (comb_r in torsions_toadd.keys())):
                tmptorlist.append(tor)
            torlist = tmptorlist
            
            torlist=[tuple(i) for i in torlist]
            torlist=[tuple([i]) for i in torlist]
            self.torsettovariabletorlist={}
            for torset in torlist:
                self.torsettovariabletorlist[tuple(torset)]=[]
            nonaroringtorlist=[tuple(i) for i in nonaroringtorlist]
            nonaroringtorlist=[tuple([i]) for i in nonaroringtorlist]
            self.rotbndtoanginc=torgen.DetermineAngleIncrementAndPointsNeededForEachTorsionSet(self,mol,self.rotbndlist)
            if self.dontdotor==True:
                torlist=[]
            
              
            self.torsettofilenametorset={}
            self.torsettotortorindex={}
            self.torsettotortorphaseindicestokeep={}
            self.nonaroringtors=[]
            self.nonaroringtorsets=[]
            self.classkeytoinitialprmguess={}
            self.nonarotortotorsbeingfit={}
            # STEP 42
            if self.atomnum<25 and len(nonaroringtorlist)==0 and self.smallmoleculefragmenter==False: 
                self.dontfrag=True
            if self.dontfrag==True and self.toroptmethod!='xtb' and 'xtb' not in self.toroptmethodlist: # if fragmenter is turned off, parition resources by jobs at sametime for parent,cant parralelize xtb since coords always written to same filename
                self.maxmem,self.maxdisk,self.numproc=self.PartitionResources()
                self.partition=True


            
            # STEP 43
            missingvdwatomsets=[]
            if self.isfragjob==False and self.dovdwscan==True:
                for vdwatomindex in missingvdwatomindices:
                    ls=tuple([tuple([vdwatomindex])])
                    missingvdwatomsets.append(ls)
                    self.torlist.append(ls)
  
            
            

            # STEP 44
            if self.isfragjob==False and not os.path.isfile(self.key5fname) and self.dontfrag==False and (self.dovdwscan==True):
                self.WriteToLog('Create and Run vdW Fragment Poltype Jobs')

                WBOmatrix,outputname,error=frag.GenerateWBOMatrix(self,self.rdkitmol,self.mol,self.logoptfname.replace('.log','.xyz'))
                rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor=frag.GenerateFragments(self,self.mol,self.torlist,WBOmatrix,missingvdwatomsets,nonaroringtorlist) # returns list of bond indexes that need parent molecule to do torsion scan for (fragment generated was same as the parent0
                equivalentrotbndindexarrays,rotbndindextoringtor,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray=frag.SpawnPoltypeJobsForFragments(self,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor)
            # STEP 45
            if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key5fname) and (self.dovdwscan==True):
                frag.GrabVdwAndTorsionParametersFromFragments(self,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor,self.key4fname,self.key5fname,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray) # just dump to key_5 since does not exist for parent molecule
            else:         
                if self.dovdwscan==True:
                    if self.dontfrag==False: 
                        if self.isfragjob==True:
                            vdwfit.VanDerWaalsOptimization(self,missingvdwatomindices)    
                    else:
                        vdwfit.VanDerWaalsOptimization(self,missingvdwatomindices)
                else:
                    shutil.copy(self.key4fname,self.key5fname)

            # STEP 46
            shutil.copy(self.key5fname,self.key6fname)
            self.torlist=torlist[:]
            if self.tortor==True and self.dontdotor==False:
                torgen.PrepareTorsionTorsion(self,optmol,mol,tortorsmissing)
            torgen.DefaultMaxRange(self,self.torlist)
            # STEP 47
            if self.refinenonaroringtors==True and self.dontfrag==False:
                rings.RefineNonAromaticRingTorsions(self,mol,optmol,classkeytotorsionparametersguess)
            # STEP 48
            if self.isfragjob==False and not os.path.isfile(self.key7fname) and self.dontfrag==False and (self.dontdotor==False) and len(self.torlist)!=0:
                self.WriteToLog('Create and Run Torsion Fragment Poltype Jobs')

                WBOmatrix,outputname,error=frag.GenerateWBOMatrix(self,self.rdkitmol,self.mol,self.logoptfname.replace('.log','.xyz'))
                highlightbonds=[]
                for torset in self.torlist:
                    for tor in torset:
                        if len(tor)>1:
                            rotbnd=[tor[1]-1,tor[2]-1]
                            highlightbonds.append(rotbnd)
                frag.Draw2DMoleculeWithWBO(self,WBOmatrix,self.molstructfname.replace('.sdf',''),self.rdkitmol,bondindexlist=highlightbonds,imgsize=1500)       
                rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor=frag.GenerateFragments(self,self.mol,self.torlist,WBOmatrix,missingvdwatomsets,nonaroringtorlist) # returns list of bond indexes that need parent molecule to do torsion scan for (fragment generated was same as the parent0
                equivalentrotbndindexarrays,rotbndindextoringtor,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray=frag.SpawnPoltypeJobsForFragments(self,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor)
            # STEP 49
            if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key7fname) and (self.dontdotor==False) and len(self.torlist)!=0:
                frag.GrabVdwAndTorsionParametersFromFragments(self,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor,self.key6fname,self.key7fname,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray) # just dump to key_6 since does not exist for parent molecule
            else:
                # Torsion scanning then fitting. *.key_7 will contain updated torsions
                if not os.path.isfile(self.key7fname):
                    if len(self.torlist)!=0:
                        # torsion scanning
                        for r in range(len(self.toroptmethodlist)):
                            self.toroptmethod=self.toroptmethodlist[r]
                            self.torspmethod=self.torspmethodlist[r]
                            torgen.gen_torsion(self,optmol,self.torsionrestraint,mol)
                        # torsion fitting
                        if self.dontdotorfit==True:
                            shutil.copy(self.key6fname,self.key7fname)
                            sys.exit()
                        torfit.process_rot_bond_tors(self,optmol)
                    else:
                        shutil.copy(self.key6fname,self.key7fname)           
           

            if self.torfit==False:
                sys.exit()
            # STEP 50 
            if self.isfragjob==False and self.dontdotor==False:
                self.CheckTorsionParameters(self.key7fname,torsionsmissing)
            self.WriteOutLiteratureReferences(self.key7fname) 
            if self.writeoutpolarize==False or self.writeoutmultipole==False:
                shutil.copy(self.xyzfname,self.xyzoutfile)
            shutil.copy(self.xyzoutfile,self.tmpxyzfile)
            shutil.copy(self.key7fname,self.tmpkeyfile)

            shutil.copy(self.key7fname,'TEST_tim.key')
            
            # STEP 51
            if self.writeoutpolarize and self.writeoutmultipole==True:
                opt.StructureMinimization(self,torsionrestraints)
                if self.atomnum != 1:
                    opt.gen_superposeinfile(self)
                    opt.CheckRMSD(self)
            # STEP 52
            if self.atomnum!=1: 
                 esp.CheckDipoleMoments(self,optmol)
            # STEP 53
            if os.path.exists(self.tmpkeyfile):
                self.FinalVDWMultipoleCheck(self.tmpkeyfile)
            # STEP 54
            for torset,totaltime in self.torsettototalqmtime.items():
                sptime=self.torsettospqmtime[torset]
                opttime=self.torsettooptqmtime[torset]
                numpoints=self.torsettonumpoints[torset]
                self.WriteToLog('Total torsion QM time for '+str(torset)+' is '+str(round(totaltime,3))+' hours'+' and '+str(numpoints)+' conformations')
                self.WriteToLog('OPT torsion QM time for '+str(torset)+' is '+str(round(opttime,3))+' hours'+' and '+str(numpoints)+' conformations')
                self.WriteToLog('SP torsion QM time for '+str(torset)+' is '+str(round(sptime,3))+' hours'+' and '+str(numpoints)+' conformations')
            # STEP 55
            if self.isfragjob==False and self.dontfrag==True and len(self.torlist)!=0:
                self.WriteOutDatabaseParameterLines()
            # STEP 56
            if os.path.exists(self.scrtmpdirgau):
                shutil.rmtree(self.scrtmpdirgau)
            if os.path.exists(self.scrtmpdirpsi4):
                shutil.rmtree(self.scrtmpdirpsi4)
            # STEP 58
            if self.isfragjob==False:
                previousdir=os.path.abspath(os.path.join(os.getcwd(), os.pardir))
                # copy the optimized xyz file
                if os.path.exists(self.tmpxyzfile + '_2'):
                     shutil.copy(self.tmpxyzfile + '_2',os.path.join(previousdir,self.tmpxyzfile))
                if os.path.exists(self.tmpkeyfile):
                     shutil.copy(self.tmpkeyfile,os.path.join(previousdir,self.tmpkeyfile))
            self.CopyFitPlots()
            os.chdir('..')
            
            if os.path.isfile(self.tmpxyzfile):
              cmdstr = f"python \"{os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'lFormatTXYZ.py')}\" {self.tmpxyzfile}"
              self.WriteToLog(cmdstr)
              os.system(cmdstr)

            # STEP 59
            # apply any modifications to the final.key
            lmodifytinkerkey.mod_final_key(self)
            
            self.WriteToLog('Poltype Job Finished'+'\n')


            


        def GenerateDuplicateXYZsFromOPTs(self,xyzfname,keyfname,optmolist):
            """
            Intent: If user gives multipole parameters as input, then need to skip calling poledit to generate tinker XYZ file, so this method makes tinker XYZ file for you. 
            Input: Tinker XYZ name, tinker key name, list of optimized mol objects (need many when fitting several conformations for multipole).
            Output: list of tinker XYZs and list of tinker keys used as input for generating QM input for electrostatic surface potential fitting for multipole parameters. 
            Referenced By: GenerateParameters
            Description:
            1.  Iterate over all mol objects
            2.  Extract the atom indices and coordinates
            3.  Tinker XYZ may already exist for first opt conformation, so append to list, else need to generate tinker XYZ fysing coordinates from optmol (other conformorations) and use initial tinker XYZ as template 
            """
            xyzfnamelist=[]
            keyfnamelist=[]
            # STEP 1
            for optmolidx in range(len(optmolist)):
                optmol=optmolist[optmolidx]
                # STEP 2
                indextocoordinates=self.GrabIndexToCoordinates(optmol)
                # STEP 3
                if optmolidx==0:
                    xyzfnamelist.append(xyzfname)
                    keyfnamelist.append(keyfname)
                    continue
                suffix='_'+str(optmolidx+1)
                newxyzfname=xyzfname.replace('.xyz',suffix+'.xyz') 
                self.GenerateTinkerXYZ(xyzfname,newxyzfname,indextocoordinates)
                xyzfnamelist.append(newxyzfname)
                keyfnamelist.append(keyfname)

            return xyzfnamelist,keyfnamelist


        def ExtractMOLInfo(self,mol):
            """
            Intent: Used for when need to extract information from mol object and then can later be used for generating tinker XYZ
            Input: mol object
            Output: Dictionaries mapping atom indices to coordinates,connectivity and elements
            Referenced By: GenerateParameters 
            Description: 
            1. Iterate over mol objec atoms
            2. Extract index, element symbol, coordinates
            3. Save information to dictionaries
            """
            atmindextocoordinates={}
            atmindextoconnectivity={}
            atmindextoelement={}
            iteratom = openbabel.OBMolAtomIter(mol)
            an = pyasl.AtomicNo()
            # STEP 1
            for atom in iteratom:
                # STEP 2
                index=atom.GetIdx()
                atomicnum=atom.GetAtomicNum()
                coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
                element=an.getElSymbol(atomicnum)
                # STEP 3
                atmindextocoordinates[index]=coords
                connectivity=self.GrabConnectivity(mol,index)
                atmindextoconnectivity[index]=connectivity
                atmindextoelement[index]=element


            return atmindextocoordinates,atmindextoconnectivity,atmindextoelement


        def GrabConnectivity(self,mol,index):
            """
            Intent: Grab connecting atom indices to input atom. Useful for reconstructing tinker XYZs.
            Input: mol object and atom of interest
            Output: List of atom indices connected to atom index of interest
            Referenced By: ExtractMOLInfo
            Description: 
            1. Iterate over all bonds
            2. Extract bond atom indices (a,b)
            3. If a is atom index of interest, append b, else if b is atom index of interest. append a
            """
            conn=[]
            # STEP 1
            bonditer=openbabel.OBMolBondIter(mol)
            for bond in bonditer:
                # STEP 2
                oendidx = bond.GetEndAtomIdx()
                obgnidx = bond.GetBeginAtomIdx()
                # STEP 3
                if oendidx==index:
                    if obgnidx not in conn:
                        conn.append(obgnidx)
                elif obgnidx==index:
                    if oendidx not in conn:
                        conn.append(oendidx)
            return conn

        def GenerateTinkerXYZ(self,xyzfname,newxyzfname,indextocoordinates):
            """
            Intent: Generate Tinker XYZ for other conformations needed for multipole fitting
            Input: Template tinker XYZ name, new tinker XYZ name, dictionary of atom indices to coordinates
            Output: -
            Referenced By: GenerateDuplicateXYZsFromOPTs
            Description: 
            1. Read in lines of template tinker XYZ into array (results)
            2. Iterate over array (results)
            3. If line contains coordinate information, extract current atom index from line then extract coordinates from dictionary
            4. Replace old coordinates with new coordinates
            5. Write new line to new file
            """
            # STEP 1
            temp=open(xyzfname,'r')
            results=temp.readlines()
            temp.close()
            temp=open(newxyzfname,'w')
            # STEP 2
            for line in results:
                linesplit=line.split()
                if len(linesplit)>1:
                    # STEP 3
                    index=int(linesplit[0])
                    coords=indextocoordinates[index-1]
                    # STEP 4
                    linesplit[2]=str(coords[0])
                    linesplit[3]=str(coords[1]) 
                    linesplit[4]=str(coords[2])
                    line=' '.join(linesplit)+'\n'
                # STEP 5
                temp.write(line)
            temp.close()


        def GenerateListOfMols(self,mol,indextocoordslist):
            """
            Intent: Need many mol objects for each QM geometry optimization, when you want to fit multipoles with many conformations
            Input: mol object, list of dictionaries mapping atom indices to coordinates for different conformations
            Output: List of mol objects
            Referenced By: GenerateParameters 
            Description: 
            1. First item in list needs to be the original mol object (extended conformation)
            2. Iterate over list of dictionaries
            3. For all other conformations in list, read in initial SDF file, then change coordinates based on dictionary information
            4. Append updated mol object to list of mol objects 
            """
            # STEP 1
            molist=[mol]
            obConversion = openbabel.OBConversion()
            # STEP 2
            for i in range(len(indextocoordslist)):
                if i!=0:
                    # STEP 3
                    indextocoordinates=indextocoordslist[i] 
                    othermol = openbabel.OBMol()
                    inFormat = obConversion.FormatFromExt(self.molstructfname)
                    obConversion.SetInFormat(inFormat)
                    obConversion.ReadFile(othermol, self.molstructfname)
                    othermol=self.SetDefaultCoordinatesBabel(othermol,indextocoordinates)
                    othermol.SetTotalCharge(self.totalcharge)
                    # STEP 4
                    molist.append(othermol)

            return molist


        def FinalVDWMultipoleCheck(self,keyfile):
            """
            Intent: Sanity check to ensure no vdw parameters or multipoles are missing from final tinker key file. Tinker will not give undefined parameter error if vdw parameters are missing!
            Input: Keyfile
            Output: - 
            Referenced By: GenerateParameters 
            Description:
            1. Read in keyfile contents to an array (results)
            2. Extract type numbers from dictionary idxtosymclass 
            3. Initialize arrays to not found any vdw/mpole parameters
            4. Iterate over results array and if detect vdw or multipole parameters, set array values to True
            5. If any values in array are False (missing) for vdw or multipole, then raise Error 
            """
            # STEP 1
            temp=open(keyfile,'r')
            results=temp.readlines()
            temp.close()
            vdwtypetofound={}
            mpoletypetofound={}
            # STEP 2
            typenums=list(self.idxtosymclass.values())
            # STEP 3
            for typenum in typenums:
                typenum=str(typenum)
                vdwtypetofound[typenum]=False
                mpoletypetofound[typenum]=False
            # STEP 4
            for line in results:
                if '#' not in line:
                    if 'vdw' in line or 'multipole' in line:
                        linesplit=line.split()
                        typenum=linesplit[1]
                        if 'vdw' in line:
                            vdwtypetofound[typenum]=True
                        elif 'multipole' in line:
                            mpoletypetofound[typenum]=True
            # STEP 5
            missingvdw=[]
            missingmpole=[]
            for typenum,found in vdwtypetofound.items():
                if found==False:
                    missingvdw.append(typenum)

            for typenum,found in mpoletypetofound.items():
                if found==False:
                    missingmpole.append(typenum)

            if len(missingvdw)!=0:
                raise ValueError('Missing vdw parameters '+str(missingvdw))


            if len(missingmpole)!=0:
                raise ValueError('Missing multipole parameters '+str(missingmpole))


        def AddIndicesToKey(self,keyfilename):
            """
            Intent: Its confusing to look at type numbers and compare to structure visualization, so add comments to keyfile that inform user of which atom indices correspond to which atom types
            Input: Keyfile
            Output: -
            Referenced By: GenerateParameters 
            Description: 
            1. Read in input keyfile to array
            2. For each parameter line, determine the indices (1,2,3 etc..) and extract type numbers from the line
            3. For each type number, determine all possible atom indices for that type number. 
            4. Construct a string with all possible atom indices for each type number and add as a comment to new keyfile.
            """
            # STEP 1
            temp=open(keyfilename,'r')
            results=temp.readlines()
            temp.close()
            tempname=keyfilename.replace('.key','_TEMP.key')
            temp=open(tempname,'w')
           
            for line in results:
                # STEP 2
                if '#' not in line and 'none' not in line:
                    indices=[]
                    if 'multipole ' in line or 'polarize ' in line or 'vdw ' in line:
                        indices=[1]
                    elif 'bond ' in line or 'opbend ' in line:
                        indices=[1,2]
                    elif 'angle ' in line or 'strbnd ' in line:
                        indices=[1,2,3]
                    elif 'torsion ' in line:
                        indices=[1,2,3,4]
                    if len(indices)!=0:
                        # STEP 3
                        allindices=[]
                        linesplit=line.split()
                        typenums=[linesplit[i] for i in indices]
                        typenums=[int(i) for i in typenums]
                        for typenum in typenums: 
                            indexes=self.GrabKeysFromValue(self.idxtosymclass,typenum)
                            allindices.append(indexes)
                        # STEP 4
                        string='# '+str(typenums) +' = '+str(allindices)+'\n'
                        temp.write(string)
                temp.write(line)


            temp.close()
            os.remove(keyfilename)
            os.replace(tempname,keyfilename)


        def GrabKeysFromValue(self,dic,thevalue):
            """
            Intent: For example, want to determine atomic indices that map to type number starting with information only of type number. Then need to find all possible keys (atomic indices) that map to value (type number).
            Input: Dictionary of keys to values, the desired value that want to search all keys that map to 
            Output: List of keys that map to desired value
            Referenced By: AddIndicesToKey
            Description: 
            1. Iterate over key, value pairs in dictionary
            2. If the current value is equivalent to the desired value, then append key to array
            3. Return array of keys
            """
            keylist=[]
            # STEP 1
            for key,value in dic.items():
                # STEP 2
                if value==thevalue:
                    keylist.append(key)
            # STEP 3
            return keylist

        
        def CollectElectrostaticDipoleFitLines(self):
            """
            Intent: Users are too lazy to read poltype log file for information about electrostatic fitting etc.., so copy relevent information to a textfile and put in the OPENME folder, in the hope that they will read it there.
            Input: - 
            Output: Array of lines containing relevant information
            Referenced By: CopyFitPlots
            Description: 
            1. Iterate over lines of poltype log
            2. If results for electrostatic potential fitting or comparison to dipole moments between QM and MM is printed then grab those lines and append to an array.
            """
            # STEP 1
            files=os.listdir()
            for f in files:
                if 'poltype.log' in f:
                    name=f
            temp=open(name,'r')
            results=temp.readlines()
            temp.close()  
            # STEP 2
            fitlines=[]
            for line in results:
                if 'RMSPD =' in line or 'QMDipole' in line or 'RMSD of QM and MM' in line:
                    fitlines.append(line)


            return fitlines


        def Instructions(self):
            """
            Intent: Attempt to make it easier for users to understand plots and fitting information in OPENME, without them having to read the documentation README.MD . 
            Input: -
            Output: Array of lines containing instructions of how to interpret results in OPENME folder
            Referenced By: CopyFitPlots
            Description: 
            1. Hard code instructions and then append to array
            """
            # STEP 1
            instructions=[]
            line='Please ensure that the Root Mean Square Deviation (RMSD) between QM (red curve) and MM2 post-fit torsion (blue curve) has a decent fit.'+'\n'
            instructions.append(line) 
            line='Plots with *energy* contain the QM energy vs dihedral angle vs MM2 energy vs dihedral angle'+'\n'
            instructions.append(line)
            line='Plots with *fit* contain the QM-MM1 (prefit torsion energy) vs dihedral angle and the fit spline'+'\n'
            instructions.append(line)
            line='The first two numbers in plot are the rotatable bond in the parent, last two are the rotatable bond in the fragment'+'\n' 
            instructions.append(line)
            line='Plots with _use_weights in the filename are Boltzman fitted plots'+'\n'
            instructions.append(line)
            line='Please ensure that the QM dipole and MM dipole match reasonably well'+'\n'
            instructions.append(line)
          
            line='Please ensure that the RMSD between average QM potential and average MM potential is small'+'\n'
            instructions.append(line)

            return instructions

        def CopyFitPlots(self):
            """
            Intent: Users will likely be unfamiliar find commands or where to find plots or too lazy to find, so grab all plots in one folder called OPENME in top directory. 
            Input: -
            Output: -
            Referenced By: GenerateParameters 
            Description: 
            1. Walk along all sub directories and look for any png files, csv files then append to array
            2. Collect electrostatic fitting results from poltype log file
            3. Collect instructions for how to interpret OPENME results
            4. Create folder OPENME (if doesnt exist)
            5. Copy plots found to OPENME folder
            6. Create text file in OPENME folder and write out results from fitting and instructions on how to interpret fitting results and also how to interpret plot results. 
            """
            plots=[]
            fold='OPENME'
            thecurdir=os.getcwd()
            # STEP 1
            for root, subdirs, files in os.walk(os.getcwd()):
                for d in subdirs:
                    curdir=os.getcwd()
                    path=os.path.join(root, d)
                    os.chdir(path)
                    if fold in path:
                        continue
                    files=os.listdir()
                    for f in files:
                        if '.png' in f and ('energy' in f or 'fit' in f or 'water' in f) or '.csv' in f:
                            plots.append(os.path.join(path,f))
            os.chdir(thecurdir)
            # STEP 2 
            fitlines=self.CollectElectrostaticDipoleFitLines()
            # STEP 3
            instructions=self.Instructions()
            os.chdir('..')
            # STEP 4
            if not os.path.isdir(fold):
                os.mkdir(fold)
            # STEP 5
            for path in plots:
                shutil.copy(path,fold)
            os.chdir(fold)
            # STEP 6
            temp=open("README.txt",'w')
            for line in instructions:
                temp.write(line)
            for line in fitlines:
                temp.write(line)
            temp.close()
            os.chdir(thecurdir)


        def CheckFloat(self,string):
            """
            Intent: Just quick check to see if character in file is float or not, for example when want to determine if line contains box size information or not
            Input: String
            Output: Boolean determining if string is float or not
            Referenced By: Various functions for modifying/reading arc/xyz files 
            Description: 
            1. try to convert string to a float, if it works return True
            2. If try fails, then return False
            """
            try:
                float(string)
                return True
            except:
                return False

            
        def SanitizeMMExecutable(self, executable):
            """
            Intent: Sometimes tinker binaries are installed with program.x rather than just program, so need a function that chooses correctly based on what is in PATH.
            Input: Executable name
            Output: New executable name
            Referenced By: GenerateParameters, _post_init_
            Description: 
            1. If input executable is found in PATH, then return executable 
            2. If not, check for executable with .x in TINKERDIR and return that executable 
            """
            # STEP 1
            if self.which(executable)!=None:
                return executable
            # STEP 2
            if self.tinkerdir is None:
                self.tinkerdir = os.getenv("TINKERDIR", default="")
            exe = os.path.join(self.tinkerdir, executable)
            if self.which(exe) is None:
                exe = exe[:-2] if exe.endswith('.x') else exe + '.x'
                if self.which(exe) is None:
                    print("ERROR: Cannot find Tinker {} executable".format(executable))
                    sys.exit(2)
            return exe

        def which(self,program):
            """
            Intent: Look for executable in PATH
            Input: executable name
            Output: Executable name (if found in PATH), else return None
            Referenced By: SanitizeMMExecutable
            Description: 
            1. Define function that tries to see if function is in PATH
            2. First call function is_exe to check if program is found in PATH
            3. If not found, then itetate over all paths in PATH, and check if the program can be found in any of those paths 
            """
            # STEP 1
            def is_exe(fpath):
                try:
                     return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
                except:
                     return None
            # STEP 2
            fpath, fname = os.path.split(program)
            if fpath:
                result=is_exe(program)
                if result:
                    return program
            else:
                # STEP 3
                for path in os.environ["PATH"].split(os.pathsep):
                    exe_file = os.path.join(path, program)
                    result=is_exe(exe_file)
                    if result:
                        return exe_file
        
            return None

            
        def ConvertTinkerXYZToCartesianXYZ(self,tinkerxyzfilename):
            """
            Intent: Convert tinker XYZ to cartesian XYZ, which then enables creating mol object for processing topology information easier.
            Input: Tinker XYZ file
            Output: Cartesian XYZ file
            Referenced By: Various functions
            Description: 
            1. Create cartesian XYZ filename
            2. Loop over lines of tinker XYZ
            3. Write out total atom number
            4. Write out only element, and coordinates in new line to cartesian XYZ file 
            """
            # STEP 1
            xyzfilename=tinkerxyzfilename.replace('.xyz','_cart.xyz')
            temp=open(tinkerxyzfilename,'r')
            tempwrite=open(xyzfilename,'w')
            results=temp.readlines()
            # STEP 2
            for lineidx in range(len(results)):
                line=results[lineidx]
                if lineidx==0:
                    # STEP 3
                    linesplit=line.split()
                    tempwrite.write(linesplit[0]+'\n')
                    tempwrite.write('\n')
                    tempwrite.flush()
                    os.fsync(tempwrite.fileno())
                else:
                    linesplit=line.split()
                    if len(linesplit)>1:
                        # STEP 4
                        newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
                        tempwrite.write(newline)
                        tempwrite.flush()
                        os.fsync(tempwrite.fileno())
            temp.close()
            tempwrite.close()
            return xyzfilename


        def ConvertTinktoXYZ(self,filename,newfilename):
            """
            Intent: For when user wants to align a new ligand tinker XYZ to a template tinker XYZ. Need to convert tinker XYZ to cartesian XYZ as an intermediate for processing with toolkits such as rdkit.
            Input: Tinker XYZ, new filename for cartesian XYZ
            Output: Cartesian xyz filename generated
            Referenced By: AlignLigandXYZToTemplateXYZ
            Description: 
            1. Iterate over tinker XYZ
            2. If line with atom number, then write that line to new file
            3. If line with coordiantes and type information and connectivity, extract only the element and coordiantes and write to the new file
            """
            temp=open(os.getcwd()+r'/'+filename,'r')
            tempwrite=open(os.getcwd()+r'/'+newfilename,'w')
            results=temp.readlines()
            # STEP 1
            for lineidx in range(len(results)):
                line=results[lineidx]
                # STEP 2
                if lineidx==0:
                    linesplit=line.split()
                    tempwrite.write(linesplit[0]+'\n')
                    tempwrite.write('\n')
                    tempwrite.flush()
                    os.fsync(tempwrite.fileno())
                else:
                    # STEP 3
                    linesplit=line.split()
                    if len(linesplit)>1:
                        newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
                        tempwrite.write(newline)
                        tempwrite.flush()
                        os.fsync(tempwrite.fileno())
            temp.close()
            tempwrite.close()
            return newfilename
       
        
        def ExtractResource(self,string):
            """
            Intent: Exract resource value from string (if need to change or convert units) 
            Input: String with value + memory unit appended
            Output: value, memory unit
            Referenced By: PartitionResources
            Description: 
            1. Parse for cases of MB and GB
            2. Extract memory string and convert memory to float 
            """
            # STEP 1
            if 'MB' in string:
                split=string.split('MB')
                memstring='MB'
            elif 'GB' in string:
                split=string.split('GB')
                memstring='GB'
            # STEP 2
            mem=float(split[0])
        
            return mem,memstring

        def PartitionResources(self):
            """
            Intent: For running QM jobs in parralel, need to take total resources allocated to parent job and split amongst many fragment poltype jobs (or if fragmenter turned off, then parralelize to many QM jobs of parent at same time)
            Input: -
            Output: Modified resources, paritioned by number of jobs that can be run at same time
            Referenced By: GenerateParameters 
            Description:
            1. Extract the memory values of RAM
            2. Divide that by number of jobs at sametime
            3. Repeat for Disk and number of processors 
            """
            # STEP 1
            maxmem,memstring=self.ExtractResource(self.maxmem)
            # STEP 2
            maxmem=int(maxmem/self.jobsatsametime)
            tempmaxmem=str(maxmem)+memstring
            # STEP 3
            maxdisk,diskstring=self.ExtractResource(self.maxdisk)
            maxdisk=int(maxdisk/self.jobsatsametime)
            tempmaxdisk=str(maxdisk)+diskstring
            numproc=math.floor(int(self.numproc)/self.jobsatsametime)
            tempnumproc=str(numproc)
        
            return tempmaxmem,tempmaxdisk,tempnumproc

        def CheckIfWholeMoleculeIsARing(self,mol):
            """
            Intent: Dont want to generate extended conformers if whole molecule is one giant ring
            Input: babel mol object 
            Output: Boolean specifying is whole molecule is one ring
            Referenced By: GenerateParameters
            Description: 
            1. Call function RingAtomicIndices to find which atoms belong to rings
            2. For each ring found, determinte its length and add to dictionary
            3. If find a ring length greater than 7 (lazy detection, dont find large than 7 atom rings), then whole molecule is a ring. 
            """
            iswholemoleculering=False
            # STEP 1
            atomindices=torsiondatabaseparser.RingAtomicIndices(self,mol)
            # STEP 2
            lengthtoring={}
            for ring in atomindices:
                lengthtoring[len(ring)]=ring
            # STEP 3
            if len(lengthtoring.keys())>0:
                maxlength=max(lengthtoring.keys())
                maxring=lengthtoring[maxlength]
                if len(maxring)>7:
                    iswholemoleculering=True
            return iswholemoleculering


        def RepresentsInt(self,s):
            """
            Intent: Detect integers in tinker XYZ/keys for modification purposes
            Input: String
            Output: Boolean specifying if an integer or not
            Referenced By: ShiftParameterTypesComplexXYZ,ShiftParameterTypes
            Description:
            1. Try to convert type to integer if it works return True
            2. If it fails, return False 
            """
            # STEP 1
            try: 
                int(s)
                return True
            # STEP 2
            except ValueError:
                return False


        def GrabXYZInfo(self,xyzfile):
            """
            Intent: Grab information from tinker XYZ file.  
            Input: Tinker xyz file
            Output: array of coordinates, types and connectivity, dictionary of index to type number , tinker xyz atom number, dictionary of index to atomic coordinates, dictionary of index to connected indices , dictionary of index to symbol
            Referenced By: ExtractLigandIndicesFromComplexXYZ, MolecularDynamics
            Description:
            1. Iterate over lines of tinker XYZ
            2. If first line, grab total atom number
            3. Otherwise, grab index, type number, coordinates, connected indices, atomic symbol
            4. Add all the information into respecive arrays/dictionaries
            """
            temp=open(xyzfile,'r')
            xyzfileresults=temp.readlines()
            temp.close()
            xyzatominfo=[]
            indextotypeindex={}
            indextocoords={}
            indextoneighbs={}
            indextosym={}
            # STEP 1
            for lineidx in range(len(xyzfileresults)):
                line=xyzfileresults[lineidx]
                linesplit=line.split()
                if lineidx==0:
                    # STEP 2
                    xyzatomnum=int(linesplit[0])
                else:
                    if len(linesplit)>1 and '90.00' not in line:
                        # STEP 3
                        index=int(linesplit[0])
                        typenum=int(linesplit[5])
                        coords=[linesplit[2],linesplit[3],linesplit[4]]
                        # STEP 4
                        indextotypeindex[index]=typenum
                        xyzatominfo.append(linesplit[2:])
                        indextocoords[index]=coords
                        indextosym[index]=linesplit[1]
                        if len(linesplit)>=6:
                            neighbs=linesplit[6:]
                            neighbs=[int(i) for i in neighbs]
                        else:
                            neighbs=[]
                        indextoneighbs[index]=neighbs
            return xyzatominfo,indextotypeindex,xyzatomnum,indextocoords,indextoneighbs,indextosym



        def WriteOutDatabaseParameterLines(self):
            """
            Intent: Need a way to collect torsion parameters and put into format poltype database can read (this is only called when fragmenter is not used).
            Input: - 
            Output: -
            Referenced By: GenerateParameters 
            Description:
            1. Grab array of atom indices that correpond to each atom in SMARTS (fragidxarray) 
            2. Grab fragment SMARTS
            3. Make dictionary of atom index to atom order in SMARTS (starting from 1,2,..)
            4. Iterate over all torsions
            5. Generate class key (type numbers) for the torsion
            6. Generate smilesposstring (torsion indices in SMARTS atom order starting from 1,..) and fragtorstring (atom indices of torsion)  
            7. Iterate over final key file
            8. Extract torsion class key (all torsion types) from each line that has torsion key word and check if in dictionary of torsions that were fit
            9. Save parameters for classkeys found that were fit
            10.If vdw keyword in line extract type and find corresponding atom index to that type. 
            11.If atom index was a vdw type being fit, then
            12.Construct comment line with vdw SMARTS, atom index in SMARTS, and parameters
            13.If comment containing quality of fit information is detected for classkey, then save it to add to database file later
            14.Iterate over final key file
            15.Extract torsion class key (all torsion types) from each line that has torsion key word and check if in dictionary of torsions that were fit
            16.If torsions were fit, then using information contstructed earlier make comment line with torsion SMARTS, positions of torsion in SMARTS and parameters to add to database file.
            17. Write out all generated lines to database file
            """
            newkey=self.tmpkeyfile.replace('.key','_TEMP.key')
            temp=open(self.tmpkeyfile,'r')
            results=temp.readlines()
            temp.close()
            temp=open(newkey,'w')
            smartstovdwlinelist={}
            valenceprmlist={}
            valkeytosmarts={}
            classkeytoparameters={}
            classkeytosmartsposarray={}
            classkeytosmarts={}
            classkeytotorsionindexes={}
            trueotherparentindextofragindex={}
            classkeytofragmentfilename={}
            classkeytofitresults={}
            for idx,symclass in self.idxtosymclass.items():
                trueotherparentindextofragindex[idx-1]=idx-1
            m=frag.mol_with_atom_index(self,self.rdkitmol)
            fragsmirks=rdmolfiles.MolToSmarts(m)
            # STEP 1
            fragidxarray=frag.GrabAtomOrder(self,fragsmirks)
            tempmol=frag.mol_with_atom_index_removed(self,self.rdkitmol) 
            # STEP 2
            fragsmarts=rdmolfiles.MolToSmarts(tempmol)
            # STEP 3
            fragindextosmartspos=frag.GenerateAtomIndexToSMARTSPosition(self,fragidxarray)
            # STEP 4
            for rotbnd,tors in self.rotbndlist.items():
                for torsion in tors:
                    # STEP 5
                    classkey=torgen.get_class_key(self,torsion[0],torsion[1],torsion[2],torsion[3])
                    # STEP 6
                    smilesposstring,fragtorstring=frag.GenerateSMARTSPositionStringAndAtomIndices(self,torsion,trueotherparentindextofragindex,fragidxarray)
                    classkeytosmartsposarray[classkey]=smilesposstring
                    classkeytosmarts[classkey]=fragsmarts
                    classkeytotorsionindexes[classkey]=fragtorstring
                    classkeytofragmentfilename[classkey]=self.molstructfname
            # STEP 7
            for lineidx in range(len(results)):
                line=results[lineidx]
                newline=line.strip()
                linesplit=newline.split()
                if line.strip().startswith('torsion') and '#' not in line and 'Missing' not in line:
                    # STEP 8
                    typea=int(linesplit[1])
                    typeb=int(linesplit[2])
                    typec=int(linesplit[3])
                    typed=int(linesplit[4])
                    prms=linesplit[5:]
                    tor=[typea,typeb,typec,typed]
                    torkey='%d %d %d %d' % (typea, typeb, typec, typed)
                    revtorkey='%d %d %d %d' % (typed, typec, typeb, typea)
                    if torkey in classkeytofragmentfilename.keys():
                        classkey=torkey
                    elif revtorkey in classkeytofragmentfilename.keys():
                        classkey=revtorkey
                    else:
                        continue
                    # STEP 9
                    classkeytoparameters[classkey]=prms
                elif 'vdw' in line and '#' not in line:
                    for clskey,smrts in classkeytosmarts.items():
                        pass
                    # STEP 10
                    linesplit=line.split() 
                    fragclasskey=linesplit[1]
                    for fragidx,symclass in self.idxtosymclass.items():
                        if symclass==int(fragclasskey):
                            break
                    fragsymclass=int(fragclasskey)
                    prms=linesplit[2:]
                    fragidx=int(fragidx)-1
                    # STEP 11
                    if fragidx in fragindextosmartspos.keys():
                        # STEP 12
                        smartspos=fragindextosmartspos[fragidx]
                        smilesposarray=[smartspos]
                        smilesposarray=[str(i) for i in smilesposarray]
                        smilespos=','.join(smilesposarray)
                        valencestring='vdw'+' % '+smrts+' % '+smilespos+' % '
                        for prm in prms:
                            valencestring+=prm+','
                        valencestring=valencestring[:-1]
                        valencestring+='\n'
                        if smrts not in smartstovdwlinelist.keys():
                            smartstovdwlinelist[smrts]=[]
                        smartstovdwlinelist[smrts].append(valencestring)

                elif 'RMSD(MM2,QM)' in line:
                    # STEP 13
                    typea=int(linesplit[2])
                    typeb=int(linesplit[3])
                    typec=int(linesplit[4])
                    typed=int(linesplit[5])
                    tor=[typea,typeb,typec,typed]
                    torkey='%d %d %d %d' % (typea, typeb, typec, typed)
                    revtorkey='%d %d %d %d' % (typed, typec, typeb, typea)
                    if torkey in classkeytofragmentfilename.keys():
                        classkey=torkey
                    elif revtorkey in classkeytofragmentfilename.keys():
                        classkey=revtorkey
                    else:
                        continue

                    classkeytofitresults[classkey]=' '.join(linesplit)+'\n'

            # STEP 14
            for line in results:
                fitline="# Fitted from Fragment "
                linesplit=line.split()
                if line.strip().startswith('torsion') and '#' not in line and 'Missing' not in line:
                    # STEP 15
                    typea=int(linesplit[1])
                    typeb=int(linesplit[2])
                    typec=int(linesplit[3])
                    typed=int(linesplit[4])
                    torkey='%d %d %d %d' % (typea, typeb, typec, typed)
                    rev='%d %d %d %d' % (typed,typec,typeb,typea)
                    if typeb<typec:
                        valkey=tuple([typeb,typec])
                    else:
                        valkey=tuple([typec,typeb])
                    # STEP 16
                    if torkey in classkeytoparameters.keys():
                        valenceprmlist,valkeytosmarts=frag.ConstructTorsionLineFromFragment(self,torkey,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarray,classkeytosmarts,classkeytotorsionindexes,temp,valenceprmlist,fitline,classkeytofitresults,valkey,valkeytosmarts)

                    elif rev in classkeytoparameters.keys():
                        valenceprmlist,valkeytosmarts=frag.ConstructTorsionLineFromFragment(self,rev,classkeytofragmentfilename,classkeytoparameters,classkeytosmartsposarray,classkeytosmarts,classkeytotorsionindexes,temp,valenceprmlist,fitline,classkeytofitresults,valkey,valkeytosmarts)

                    else:
                        temp.write(line)


            temp.close()
            # STEP 17
            if len(self.torlist)!=0:
                frag.WriteOutDatabaseLines(self,valenceprmlist,valkeytosmarts,smartstovdwlinelist)


if __name__ == '__main__':
    def RunPoltype():
        """
        Intent: Call main poltype function, if error print error stack trace.
        Input: -
        Output: -
        Referenced By:  -
        Description: 
        1. Try to run poltype
        2. If program crashes
        3. If scratch directories exist, remove them
        """
        poltype=PolarizableTyper() 
        # STEP 1
        try:
            poltype.main()
        except:
            # STEP 2
            traceback.print_exc(file=sys.stdout)
            text = str(traceback.format_exc())
            poltype.WriteToLog(text)
            # STEP 3
            if os.path.exists(poltype.scrtmpdirgau):
                shutil.rmtree(poltype.scrtmpdirgau)
            if os.path.exists(poltype.scrtmpdirpsi4):
                shutil.rmtree(poltype.scrtmpdirpsi4)
            raise ValueError('Houston, we have a problem. Buy a developer some coffee!')
    RunPoltype()
