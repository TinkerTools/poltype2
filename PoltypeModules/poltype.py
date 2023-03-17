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
import csv
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
import databaseparser
import dimorphite_dl
import torsiongenerator as torgen
import symmetry as symm
import docking
import torsionfit as torfit
import optimization as opt
import electrostaticpotential as esp
import multipole as mpole
import fragmenter as frag
import rings
from packaging import version
from rdkit import Chem
from rdkit.Chem import rdmolfiles,AllChem,rdmolops
from rdkit.Geometry import Point3D
import vdwfit
import numpy as np
import itertools
from rdkit.Chem import rdMolTransforms
import smtplib
import textwrap
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
import psutil
import multiprocessing
import annihilation as ann
import tables 
import boxsetup as box
import time
import pdbxyz
import restraints
import plots
import submitjobs as submit
import keyfilemodifications as keymods
import re
from scipy.optimize import fmin
import pylab as plt
from scipy.interpolate import interp1d
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdFMCS
from PyAstronomy import pyasl
from dataclasses import dataclass,field
from pathlib import Path
from rdkit.Chem.MolStandardize import rdMolStandardize
from itertools import groupby
from operator import itemgetter
from rdkit.Chem import rdDistGeom
from itertools import product,combinations
import random
import productiondynamics as prod

@dataclass
class PolarizableTyper():
        numbergpus:int=1 # for estimating ETA for dynamics
        estimatedynamictimeonly:bool=False
        etafilename:str='ETA.csv'
        estimatedynamictime:bool=True
        xyzedittranslatevalue:str=''
        xyzeditstrayvalue:str=''
        xyzeditappendvalue:str=''
        xyzeditperiodicvalue:str=''
        xyzeditsoakvalue:str=''
        xyzeditionvalue:str=''
        xyzedittranslatestring:str='Translate All Atoms by an X,Y,Z-Vector'
        xyzeditstraystring:str='Move Stray Molecules into Periodic Box'
        xyzeditappendstring:str='Append a Second XYZ File to Current One'
        xyzeditperiodicstring:str='Create and Fill a Periodic Boundary Box'
        xyzeditsoakstring:str='Soak Current Molecule in Box of Solvent'
        xyzeditionstring:str='Place Monoatomic Ions around a Solute'
        needrot:bool=False
        heavyhyd:bool=False
        maxtorresnitrogen:int=2
        skipchargecheck:bool=False
        useuniquefilenames:bool=False # if users want to have unique filenames for molecular dynamics/BAR otherwise keep same filename to make copying easier from folder to folder
        xtbtorresconstant:float=5
        alignPDB:bool=False
        torfit:bool=True
        makexyzonly:bool=False
        toroptmethodlist:list=field(default_factory=lambda : [])
        torspmethodlist:list=field(default_factory=lambda : [])
        anienvname:str='ani'
        xtbenvname:str='xtbenv'
        anifmax=.05
        anipath:str=os.path.join(os.path.abspath(os.path.split(__file__)[0]),'ani.py')
        complexationonly:bool=False
        removesolventpdbtraj:bool=True
        generatepdbtrajs:bool=True
        alchemical:bool=True
        covalentdock:bool=False
        xtbmethod:int=2
        optloose:bool=True
        inputkeyfile:None=None
        writeoutmultipole:bool=True
        writeoutbond:bool=True
        writeoutangle:bool=True
        writeoutstrbnd:bool=True
        writeoutopbend:bool=True
        writeoutvdw:bool=True
        writeoutpolarize:bool=True
        writeouttorsion:bool=True
        dontrotbndslist:list=field(default_factory=lambda : [])
        relaxFBbox:bool=False
        ecrexpect:float=10
        listofligands:list=field(default_factory=lambda : [])
        targetenthalpyerror:float=.24 # kcal/mol
        targetdensityerror:float=10 # kg/m^3
        gridspacing:float=.4
        nposes:int=10
        vinaexhaustiveness:int=32
        dockgridsize:list=field(default_factory=lambda : [20, 20, 20])
        dockgridcenter:None=None
        usead4:bool=False
        usegold:bool=False
        usevina:bool=False
        usevinardo:bool=False
        goldbin:str='gold_auto'
        dockingenvname:str='dockingprep'
        prepdockscript:str=os.path.join(os.path.split(__file__)[0],'preparedockingfiles.py')
        indextompoleframefile:None=None
        qmrelativeweight:float=1
        liqrelativeweight:float=1.5
        enthalpyrelativeweight:float=1
        densityrelativeweight:float=.01
        indextotypefile:None=None
        pdbcode:None=None
        usesymtypes: bool = True
        barinterval: int =5
        visfolder: str ='VisualizationFiles'
        xyzpdbpath: str ='xyzpdb'
        extractinterforbinding: bool =False
        equiltimeionNPT: float =1
        templateligandxyzfilename: None =None
        templateligandfilename: None=None
        salthfe:bool=False
        runjobslocally:bool=True
        gpucardnumber:int=0
        density:None=None
        neatliquidsim:bool=False
        amoeba09prmfilepath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoeba09.prm"
        endstatestructurefile:None=None
        redobar:bool=False
        rotateframes:bool=False
        perturbedkeyfilename:None=None
        writeinputfiles:bool=False
        printjobsleft:bool=False
        cpujobsonly:bool=False
        minfinished:bool=False
        fep:bool=False
        submitlocally:None=None
        didinputsubmitlocally:bool=False
        numinterpolsforperturbkeyfiles:None=None
        checktraj:bool=False
        expfreeenergy:None=None
        pathtosims:None=None
        simulationstostopfolderpath:None=None
        equilfinished:bool=False
        barfilesfinished:bool=False
        perturbkeyfilelist:None=None
        boxonly:bool=False
        preequilboxpath:str=os.path.join(os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir)),os.path.join('PreEquilibriatedBoxes','waterhuge.xyz'))
        usepreequilibriatedbox:bool=False
        lastNVTequiltime:float=.5
        norotrestrainsphereradius:float=2
        aaxis:None=None
        baxis:None=None
        caxis:None=None
        prodmdfinished:bool=False
        numproc:None=None
        usetinkerforthermoprops:bool=False
        productiondynamicsNVT:bool=False
        generateinputfilesonly:bool=False
        annihilatorpath:str=os.path.abspath(os.path.split(__file__)[0])
        changegasphaseintegrator:bool=False
        annihilatevdw:bool=True
        endstatekey:None=None
        bgnstatekey:None=None
        endstatexyz:None=None
        bgnstatexyz:None=None
        restrainreceptorligand:bool=True
        minonly:bool=False
        usegpu:bool=False
        truedynamicpath:None=None
        truebarpath:None=None
        equilonly:bool=False
        binding:bool=False
        proddyngrprests:bool=True
        equilrestrainsphereradius:float=2
        restrainpositionconstant:float=5.0
        ligandfilename:None=None
        tightmincriteria:float=1
        loosemincriteria:float=10
        rescorrection:float=0
        anglerestraintconstant:float=0.003046
        pdbxyzpath:str='pdbxyz'
        distancerestraintconstant:float=10.0
        poleditpath:str='poledit'
        minimizepath:str='minimize'
        analyzepath:str='analyze'
        potentialpath:str='potential'
        averageenergies:bool=False
        complexedproteinpdbname:None=None
        uncomplexedproteinpdbname:None=None
        addphysioions:bool=True
        equilibriatescheme:list=field(default_factory=lambda : [50,100,150,200,300,300])
        equilibriaterestscheme:list=field(default_factory=lambda : [5,2,1,.5,.1,0])
        prmfilepath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18.prm"
        keyfilenamelist:None=None
        xyzfilename:None=None
        restrainatomsduringminimization:bool=True
        restrainatomgroup1:None=None
        restrainatomgroup2:None=None
        ligandxyzfilenamelist:None=None
        annihilateligandxyzfilenamelist:None=None
        receptorligandxyzfilename:None=None
        xyzeditpath:str='xyzedit'
        lowerperf:float=7
        upperperf:float=12
        torsionrestlist:None=None
        estatlambdascheme:list=field(default_factory=lambda : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
        vdwlambdascheme:list=field(default_factory=lambda : [0,.45,.52,.56,.58,.6,.62,.64,.67,.7,.75,.8,.85,.9,.95,1,1,1,1,1,1,1,1,1,1,1])
        restlambdascheme:list=field(default_factory=lambda : [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0])
        torsionrestscheme:list=field(default_factory=lambda : [.1,.95,.09,.08,.075,.07,.06,.055,.05,.04,.035,.03,.02,.015,.01,0,0,0,0,0,0,0,0,0,0,0])
        waitingtime:float=5
        receptorcharge:float=0
        ligandcharge:None=None
        barpath:str='bar'
        dynamicpath:str='dynamic'
        analyzeommpath:str='analyze9'
        barommpath:str='bar9'
        dynamicommpath:str='dynamic9'
        minimizeommpath:str='minimize9'
        complexation:bool=False
        solvation:bool=False
        flatbotrest:bool=True
        logname:str='TINKER.log'
        equilwritefreq:float=100
        proddynwritefreq:float=2
        equiltimeNVT:int=4
        equiltimeNPT:float=1
        equiltimestep:float=2
        proddyntimestep:float=2
        proddyntime:float=5
        compproddyntime:float=10
        solvproddyntime:float=5
        pressure:float=1
        NVTensem:int=2
        NPTensem:int=4
        vdwcutoff:float=12
        ewaldcutoff:float=7
        polareps:float=0.00001
        barostatmethod:str='montecarlo'
        integrator:str='RESPA'
        thermostat:str='BUSSI'
        listofsaltcons:str='[KCl]=100'
        fixvdwtyperadii:list=field(default_factory=lambda : [])
        maxjobsatsametime:float=10
        onlyrottortorlist:list=field(default_factory=lambda : [])
        numespconfs:int=1
        fitred:bool=False
        vdwtypestoeval:list=field(default_factory=lambda : [])
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
        fragmentjobslocal:bool=False
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
        email:None=None
        firstoptfinished:bool=False
        optonly:bool=False
        onlyvdwatomindex:None=None
        use_qmopt_vdw:bool=False
        use_gau_vdw:bool=False
        dontusepcm:bool=False
        deleteallnonqmfiles:bool=False
        totalcharge:None=None
        torspbasissethalogen:str="6-311G*"
        homodimers:bool=False
        tortormissingfilename:str='tortormissing.txt'
        tordebugmode:bool=False
        amoebapluscfsmartstocommentmap:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebapluscfsmartstocomment.txt'
        amoebapluscfprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'cfprmlib.txt'
        amoebaplusnonbondedprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebaplusnonbonded.prm'
        amoebaplusnonbondedsmartstocommentmap:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebaplusnonbonded.txt'
        smartstosoluteradiimap:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'SMARTsToSoluteRadiiMap.txt'
        latestsmallmoleculepolarizeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarize.prm'
        updatedsmallmoleculepolarizeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba_and_amoebaplus_polarize.prm'
        latestsmallmoleculesmartstotypespolarize:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarcommenttoparameters.txt'
        latestsmallmoleculesmartstotinkerclass:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21smartstoclass.txt'
        latestsmallmoleculeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21.prm'
        boltzmantemp:float=40
        dovdwscan:bool=False
        vdwprobepathname:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/VdwProbes/'
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
        transferanyhydrogentor:bool=True
        smallmoleculesmartstotinkerdescrip:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstoamoebatypedescrip.txt'
        smallmoleculeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba09.prm'
        torspbasissetfile:str='6-311+g_st_.0.gbs'
        toroptbasissetfile:str='6-311g_st_.0.gbs'
        optbasissetfile:str='6-311g_st_.0.gbs'
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
        externalapi:None=None
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
        toroptmethod:str='xtb'
        torspmethod:str='wB97X-D'
        dmamethod:str='MP2'
        espmethod:str='MP2'
        qmonly:bool = False
        espfit:bool = True
        parmtors:bool = True
        foldnum:int=3
        foldoffsetlist:list = field(default_factory=lambda : [ 0.0, 180.0, 0.0, 180.0, 0.0, 180.0 ])
        torlist:None = None
        rotbndlist:None = None
        maxRMSD:float=1
        maxRMSPD:float=1
        maxtorRMSPD:float=1
        tordatapointsnum:None=None
        gentorsion:bool=False
        gaustorerror:bool=False
        torsionrestraint:float=.5*3282.80354574
        torsionprmrestraintfactor:float=1.0
        onlyrotbndslist:list=field(default_factory=lambda : [])
        rotalltors:bool=False
        dontdotor:bool=False
        dontdotorfit:bool=False
        toroptpcm:bool=False
        optpcm:bool=False
        torsppcm:bool=False
        use_gaus:bool=False
        use_gausoptonly:bool=False
        use_psi4_geometric_opt:bool=True
        freq:bool=False
        postfit:bool=False
        bashrcpath:None=None
        optmaxcycle:int=400
        torkeyfname:None=None
        gausoptcoords:str=''
        forcefield:str="AMOEBA"
        helpfile:str='README.md'
        versionfile:str=os.path.join('VersionFiles','version.md')
        sleeptime:float=.1
        structure:None=None
        espextraconflist:list=field(default_factory=lambda : [])
        usepdb2pqr:bool=False
        inputproddyntime:bool=False
        inputequiltimeNVT:bool=False
        inputequiltimeNPT:bool=False
        inputlastNVTequiltime:bool=False
        indicestorestrain:list=field(default_factory=lambda : [])
        hetatmindices:list=field(default_factory=lambda : [])
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
            self.knownresiduesymbs=['A', 'G', 'U', 'A3', 'A5','DA', 'DC', 'DG', 'DT', 'G3', 'G5', 'U3', 'U5','DA3', 'DA5', 'DC3', 'DC5', 'DG3', 'DG5', 'DT3', 'DT5','ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIS','ILE', 'LEU', 'LYD', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYD', 'TYR', 'VAL', 'CYD', 'CYS']
            self.molstructfname=self.structure 
            self.simpath=os.getcwd()        
            self.estatlambdascheme=self.estatlambdascheme[::-1]
            self.vdwlambdascheme=self.vdwlambdascheme[::-1]
            self.restlambdascheme=self.restlambdascheme[::-1]
            self.elementsymtocharge={'K':1,'Cl':-1,'Mg':2,'Li':1,'Na':1,'Rb':1,'Cs':1,'Be':2,'Ca':2,'Zn':2}
            self.elementsymtomass={'H':1.00794,'HN':1.00794,'HO':1.00794,'O':15.9994,'OH':15.9994,'N':14.0067,'C':12.0107,'CA':12.0107,'F':18.9984032,'Cl':35.453,'Cl-':35.453,'S':32.065,'P':30.9737,'Na':22.98976,'K':39.0983,'Ca':40.078,'Mg':24.305,'Mg+':24.305,'K+':39.0983}
            if self.usepreequilibriatedbox==True:
                self.equiltimeNVT=2
                self.equiltimeNPT=1
                self.lastNVTequiltime=.5

            
            
            
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
                            newline=line

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
                        elif 'pdbcode' in newline:
                            self.pdbcode=a
                        elif 'indextotypefile' in newline:
                            self.indextotypefile=a
                        elif 'indextompoleframefile' in newline:
                            self.indextompoleframefile=a
                        elif 'xtbtorresconstant' in newline:
                            self.xtbtorresconstant=float(a)
                        elif 'qmrelativeweight' in newline:
                            self.qmrelativeweight=float(a)
                        elif 'ecrexpect' in newline:
                            self.ecrexpect=float(a)
                        elif 'targetenthalpyerror' in newline:
                            self.targetenthalpyerror=float(a)
                        elif 'targetdensityerror' in newline:
                            self.targetdensityerror=float(a)
                        elif 'liqrelativeweight' in newline:
                            self.liqrelativeweight=float(a)
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
                        elif 'templateligandfilename' in newline:
                            self.templateligandfilename=a
                        elif 'prmfilepath' in newline:
                            self.prmfilepath=a
                        elif 'gpucardnumber' in newline:
                            self.gpucardnumber=int(a)
                        elif 'dockgridsize' in newline:
                            self.dockgridsize=a.split(',')
                            self.dockgridsize=[i.strip() for i in self.dockgridsize]
                        elif 'dockgridcenter' in newline:
                            self.dockgridcenter=a.split(',')
                            self.dockgridcenter=[i.strip() for i in self.dockgridcenter]
                        elif 'gridspacing' in newline:
                            self.gridspacing=int(a)
                        elif 'nposes' in newline:
                            self.nposes=int(a)
                        elif 'vinaexhaustiveness' in newline:
                            self.vinaexhaustiveness=int(a)
                        elif 'dockingenvname' in newline:
                            self.dockingenvname=a
                        elif "estimatedynamictime" in newline and 'estimatedynamictimeonly' not in newline:
                            self.estimatedynamictime=self.SetDefaultBool(line,a,True)
                        elif "estimatedynamictimeonly" in newline:
                            self.estimatedynamictimeonly=self.SetDefaultBool(line,a,True)
                        elif "heavyhyd" in newline:
                            self.heavyhyd=self.SetDefaultBool(line,a,True)
                        elif "skipchargecheck" in newline:
                            self.skipchargecheck=self.SetDefaultBool(line,a,True)
                        elif "useuniquefilenames" in newline:
                            self.useuniquefilenames=self.SetDefaultBool(line,a,True)
                        elif "makexyzonly" in newline:
                            self.makexyzonly=self.SetDefaultBool(line,a,True)
                        elif "complexationonly" in newline:
                            self.complexationonly=self.SetDefaultBool(line,a,True)
                        elif "generatepdbtrajs" in newline:
                            self.generatepdbtrajs=self.SetDefaultBool(line,a,True)
                        elif "removesolventpdbtraj" in newline:
                            self.removesolventpdbtraj=self.SetDefaultBool(line,a,True)
                        elif "writeoutmultipole" in newline:
                            self.writeoutmultipole=self.SetDefaultBool(line,a,True)
                        elif "alchemical" in newline:
                            self.alchemical=self.SetDefaultBool(line,a,True)
                        elif "covalentdock" in newline:
                            self.covalentdock=self.SetDefaultBool(line,a,True)
                        elif "writeoutbond" in newline:
                            self.writeoutbond=self.SetDefaultBool(line,a,True)
                        elif "writeoutangle" in newline:
                            self.writeoutangle=self.SetDefaultBool(line,a,True)
                        elif "optloose" in newline:
                            self.optloose=self.SetDefaultBool(line,a,True)
                        elif "writeoutstrbnd" in newline:
                            self.writeoutstrbnd=self.SetDefaultBool(line,a,True)
                        elif "writeoutopbend" in newline:
                            self.writeoutopbend=self.SetDefaultBool(line,a,True)
                        elif "writeoutvdw" in newline:
                            self.writeoutvdw=self.SetDefaultBool(line,a,True)
                        elif "writeoutpolarize" in newline:
                            self.writeoutpolarize=self.SetDefaultBool(line,a,True)
                        elif "writeouttorsion" in newline:
                            self.writeouttorsion=self.SetDefaultBool(line,a,True)
                        elif "usevinardo" in newline:
                            self.usevinardo=self.SetDefaultBool(line,a,True)
                        elif "usevina" in newline:
                            self.usevina=self.SetDefaultBool(line,a,True)
                        elif "usegold" in newline:
                            self.usegold=self.SetDefaultBool(line,a,True)
                        elif "usead4" in newline:
                            self.usead4=self.SetDefaultBool(line,a,True)
                        elif "addphysioions" in newline:
                            self.addphysioions=self.SetDefaultBool(line,a,True)
                        elif "submitlocally" in newline:
                            self.didinputsubmitlocally=True
                            self.submitlocally=self.SetDefaultBool(line,a,True)
                        elif "usesymtypes" in newline:
                            self.usesymtypes=self.SetDefaultBool(line,a,True)
                        elif "printjobsleft" in newline:
                            self.printjobsleft=self.SetDefaultBool(line,a,True)
                        elif "salthfe" in newline:
                            self.salthfe=self.SetDefaultBool(line,a,True)
                        elif "relaxFBbox" in newline:
                            self.relaxFBbox=self.SetDefaultBool(line,a,True)
                        elif "extractinterforbinding" in newline:
                            self.extractinterforbinding=self.SetDefaultBool(line,a,True)
                        elif "redobar" in newline:
                            self.redobar=self.SetDefaultBool(line,a,True)
                        elif "rotateframes" in newline:
                            self.rotateframes=self.SetDefaultBool(line,a,True)
                        elif "checktraj" in newline:
                            self.checktraj=self.SetDefaultBool(line,a,True)
                        elif "cpujobsonly" in newline:
                            self.cpujobsonly=self.SetDefaultBool(line,a,True)
                        elif "writeinputfiles" in newline:
                            self.writeinputfiles=self.SetDefaultBool(line,a,True)
                        elif 'expfreeenergy' in newline:
                            self.expfreeenergy=a
                        elif 'pathtosims' in newline:
                            self.pathtosims=a
                        elif 'perturbedkeyfilename' in newline:
                            self.perturbedkeyfilename=a
                        elif 'simulationstostopfolderpath' in newline:
                            self.simulationstostopfolderpath=a
                        elif 'aaxis' in newline:
                            self.aaxis=float(a)
                        elif 'energycutoff' in newline:
                            self.energycutoff=float(a)
                        elif 'baxis' in newline:
                            self.baxis=float(a)
                        elif 'caxis' in newline:
                            self.caxis=float(a)
                        elif 'density' in newline:
                            self.density=float(a)
                        elif "usepreequilibriatedbox" in newline:
                            self.usepreequilibriatedbox=self.SetDefaultBool(line,a,True)
                        elif "equilfinished" in newline:
                            self.equilfinished=self.SetDefaultBool(line,a,True)
                        elif "barfilesfinished" in newline:
                            self.barfilesfinished=self.SetDefaultBool(line,a,True)
                        elif "minfinished" in newline:
                            self.minfinished=self.SetDefaultBool(line,a,True)
                        elif "boxonly" in newline:
                            self.boxonly=self.SetDefaultBool(line,a,True)
                        elif "fep" in newline:
                            self.fep=self.SetDefaultBool(line,a,True)
                        elif "generateinputfilesonly" in newline:
                            self.generateinputfilesonly=self.SetDefaultBool(line,a,True)
                        elif "productiondynamicsNVT" in newline:
                            self.productiondynamicsNVT=self.SetDefaultBool(line,a,True)
                        elif "usetinkerforthermoprops" in newline:
                            self.usetinkerforthermoprops=self.SetDefaultBool(line,a,True)
                        elif "changegasphaseintegrator" in newline:
                            self.changegasphaseintegrator=self.SetDefaultBool(line,a,True)
                        elif "prodmdfinished" in newline:
                            self.prodmdfinished=self.SetDefaultBool(line,a,True)
                        elif 'complexedproteinpdbname' in newline:
                            self.complexedproteinpdbname=a
                        elif 'dynamicpath' in newline:
                            self.dynamicpath=a
                        elif 'dynamicommpath' in newline:
                            self.dynamicommpath=a
                        elif 'externalapi' in newline:
                            self.externalapi=a
                        elif 'bashrcpath' in newline:
                            self.bashrcpath=a
                        elif 'restrainatomgroup1' in newline:
                            self.restrainatomgroup1=[int(i) for i in commalist]
                        elif 'espextraconflist' in newline:
                            self.espextraconflist=a.split(',')
                            self.espextraconflist=[i.strip() for i in self.espextraconflist]
                        elif 'perturbkeyfilelist' in newline:
                            self.perturbkeyfilelist=commalist
                        elif 'numinterpolsforperturbkeyfiles' in newline:
                            self.numinterpolsforperturbkeyfiles=commalist
                        elif 'restrainatomgroup2' in newline:
                            self.restrainatomgroup2=[int(i) for i in commalist]
                        elif 'torsionrestlist' in newline:
                            self.torsionrestlist=commalist
                            templist=[]
                            for ele in self.torsionrestlist:
                                paths=ele.lstrip().rstrip().split()
                                temp=[]
                                for e in paths:
                                    temp.append(e)
                                templist.append(temp)
                            self.torsionrestlist=templist
                        elif "flatbotrest" in newline:
                            self.flatbotrest=self.SetDefaultBool(line,a,True)
                        elif "useproddyngrprests" in newline:
                            self.useproddyngrprests=self.SetDefaultBool(line,a,True)
                        elif "restrainreceptorligand" in newline:
                            self.restrainreceptorligand=self.SetDefaultBool(line,a,True)
                        elif "annihilatevdw" in newline:
                            self.annihilatevdw=self.SetDefaultBool(line,a,True)
                        elif "proddynensem" in newline:
                            self.proddynensem = a
                        elif ("keyfilenamelist") in newline:
                            self.keyfilenamelist=commalist
                        elif ("ligandcharge") in newline:
                            self.ligandcharge = int(a)
                        elif ("tightmincriteria") in newline:
                            self.tightmincriteria = float(a)
                        elif ("restrainpositionconstant") in newline:
                            self.restrainpositionconstant = float(a)
                        elif ("loosemincriteria") in newline:
                            self.loosemincriteria = float(a)
                        elif ("equilrestrainsphereradius") in newline:
                            self.equilrestrainsphereradius = float(a)
                        elif ("barostatmethod") in newline:
                            self.barostatmethod = a
                        elif ("receptorcharge") in newline:
                            self.receptorcharge = float(a)
                        elif ("waitingtime") in newline:
                            self.waitingtime = float(a)
                        elif ("listofsaltcons") in newline:
                            self.listofsaltcons = a
                        elif ("ligandxyzfilenamelist") in newline and "receptor" not in newline and 'annihilate' not in newline:
                            self.ligandxyzfilenamelist=commalist
                        elif ("annihilateligandxyzfilenamelist") in newline and "receptor" not in newline:
                            self.annihilateligandxyzfilenamelist=commalist
                        elif ("ligandfilename") in newline:
                            self.ligandfilename = a
                        elif ("receptorligandxyzfilename") in newline:
                            self.receptorligandxyzfilename = a
                        elif ("integrator") in newline:
                            self.integrator = a
                        elif ("thermostat") in newline:
                            self.thermostat = a
                        elif ("vdwcutoff") in newline:
                            self.vdwcutoff = a
                        elif ("ewaldcutoff") in newline:
                            self.ewaldcutoff = a
                        elif ("polareps") in newline:
                            self.polareps = a
                        elif ("bgnstatexyz") in newline:
                            self.bgnstatexyz = a
                        elif ("endstatexyz") in newline:
                            self.endstatexyz = a
                        elif ("bgnstatekey") in newline:
                            self.bgnstatekey = a
                        elif ("endstatekey") in newline:
                            self.endstatekey = a
                        elif ("distancerestraintconstant") in newline:
                            self.distancerestraintconstant= float(a)
                        elif ("anglerestraintconstant") in newline:
                            self.anglerestraintconstant= float(a)
                        elif ("equilibriatescheme") in newline:
                            self.equilibriatescheme=commalist
                        elif ("restlambdascheme") in newline:
                            self.restlambdascheme=commalist
                        elif ("vdwlambdascheme") in newline:
                            self.vdwlambdascheme=commalist
                        elif ("estatlambdascheme") in newline:
                            self.estatlambdascheme=commalist
                        elif ("torsionrestscheme") in newline:
                            self.torsionrestscheme=commalist
                        elif ("equilwritefreq") in newline:
                            self.equilwritefreq= a
                        elif ("equiltimestep") in newline:
                            self.equiltimestep= int(a)
                        elif ("proddyntimestep") in newline:
                            self.proddyntimestep= int(a)
                        elif ("numbergpus") in newline:
                            self.numbergpus= int(a)
                        elif ("equiltimeNPT") in newline:
                            self.equiltimeNPT= float(a)
                            self.inputequiltimeNPT=True
                        elif ("equiltimeionNPT") in newline:
                            self.equiltimeionNPT= float(a)
                        elif ("proddynwritefreq") in newline:
                            self.proddynwritefreq= a
                        elif ("proddyntime") in newline and 'comp' not in newline and 'solv' not in newline:
                            self.proddyntime= float(a)
                            self.inputproddyntime=True
                        elif ("compproddyntime") in newline:
                            self.compproddyntime= float(a)
                        elif ("solvproddyntime") in newline:
                            self.solvproddyntime= float(a)
                        elif ("complexation") in newline:
                            self.complexation=True
                        elif ("equiltimeNVT") in newline:
                            self.equiltimeNVT= float(a)
                            self.inputequiltimeNVT=True
                        elif ("lastNVTequiltime") in newline:
                            self.lastNVTequiltime= float(a)
                            self.lastNVTequiltime=True
                        elif ("equilonly") in newline:
                            self.equilonly=True
                        elif ("minonly") in newline:
                            self.minonly=True
                        elif ("equilrestlambdascheme") in newline:
                            self.equilrestlambdascheme=commalist

                        elif "rotalltors" in newline:
                            self.rotalltors=self.SetDefaultBool(line,a,True)

                        elif "vdwprmtypestofit" in newline:
                            self.vdwprmtypestofit=a.split(',')
                            self.vdwprmtypestofit=[i.strip() for i in self.vdwprmtypestofit]

                        elif "fixvdwtyperadii" in newline:
                            self.fixvdwtyperadii=a.split(',')
                            self.fixvdwtyperadii=[i.strip() for i in self.fixvdwtyperadii]

                        elif "vdwtypestoeval" in newline:
                            self.vdwtypestoeval=a.split(',')
                            self.vdwtypestoeval=[i.strip() for i in self.vdwtypestoeval]
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
                        elif "fragmentjobslocal" in newline:
                            self.fragmentjobslocal=self.SetDefaultBool(line,a,True)
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
                        elif "email" in newline:
                            self.email=a
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
                        elif "externalapi" in newline and a!='None':
                            self.externalapi=a
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
                        elif "bashrcpath" in newline and a!='None':
                            self.bashrcpath = a
                        elif "structure" in newline:
                            self.molstructfname = a
                        elif "dontusepcm" in newline:
                            self.dontusepcm=self.SetDefaultBool(line,a,True)
                        elif "torsppcm" in newline:
                            self.torsppcm=self.SetDefaultBool(line,a,True)
                        elif "freq" in newline:
                            self.freq=self.SetDefaultBool(line,a,True)
                        elif "optpcm" in newline and 'tor' not in line:
                            self.optpcm=self.SetDefaultBool(line,a,True)
                        elif "toroptpcm" in newline:
                            self.toroptpcm=self.SetDefaultBool(line,a,True)
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
                        elif "torsionprmrestraintfactor" in newline:
                            self.torsionprmrestraintfactor=float(a)
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
                        elif 'defaultmaxtorsiongridpoints' in newline:
                            self.defaultmaxtorsiongridpoints=int(a)
                        elif "optbasisset" in newline and 'tor' not in newline:
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
                        elif "gbindir" in newline:
                            self.gausdir = a
                        elif "forcefield" in newline:
                            self.forcefield = a
                        elif "qmonly" in newline:
                            self.qmonly=self.SetDefaultBool(line,a,True)
                        elif "sleeptime" in newline:
                            self.sleeptime = float(a)
                        else:
                            print('Unrecognized '+line)
                            sys.exit()

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
               self.tempjobsatsametime=self.jobsatsametime
               self.tempmaxmem=self.maxmem
               self.tempmaxdisk=self.maxdisk
               self.tempnumproc=self.numproc
               self.partition=False
            self.firsterror=False
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
            if self.covalentdock==True:
                self.usevina=False # doesnt have covalent
                self.usevinardo=False # doesnt have covalent
                self.usead4=False
            if (self.usead4==True or self.usegold==True or self.usevina==True or self.usevinardo==True):
                self.docking=True
            else:
                self.docking=False
            if (self.complexedproteinpdbname!=None or self.receptorligandxyzfilename!=None) and self.docking==False:
                self.binding=True
                self.solvation=True
                self.complexation=True

            if self.binding==False:
                if self.density!=None:
                    self.neatliquidsim=True
            if self.binding==False and self.neatliquidsim==False and self.docking==False and self.molstructfname==None:
                self.solvation=True
            if self.uncomplexedproteinpdbname!=None and self.complexedproteinpdbname==None:
                self.usepdb2pqr=True
            if self.uncomplexedproteinpdbname!=None and self.complexedproteinpdbname!=None:
                self.alignPDB=True

            
            if self.ligandxyzfilenamelist!=None and (self.binding==True or self.solvation==True or self.neatliquidsim==True) or self.pdbcode!=None or self.usepdb2pqr!=False or self.alignPDB==True:
                self.MolecularDynamics()
                sys.exit()
            if (self.docking==True) and self.complexedproteinpdbname!=None:
                docking.DockingWrapper(self,self.complexedproteinpdbname,self.dockgridcenter,self.dockgridsize,self.vinaexhaustiveness,self.nposes,self.goldbin,self.usevina,self.usead4,self.usevinardo,self.usegold,self.dockingenvname,self.prepdockscript,self.gridspacing,self.listofligands,self.ecrexpect,self.covalentdock,self.knownresiduesymbs)
                sys.exit()
            
            if self.isfragjob==False:
                self.parentname=self.molecprefix
            else:
                self.parentname=str(self.parentname)
            
 
            # Use openbabel to create a 'mol' object from the input molecular structure file. 
            # Openbabel does not play well with certain molecular structure input files,
            # such as tinker xyz files. (normal xyz are fine)
    
            
            if not __name__ == '__main__':
                params=self.main()


        def MolecularDynamics(self):
            """
            Intent: Handles initializing more variables for MD simulations/FEP
            Input: Variables from poltype.ini file and other default variables
            Output: Variables needed for molecular dynamics/FEP then calls functions for molecular dynamics/FEP
            Referenced By: __post_init__ 
            Description:
            0. Determine xyzedit options, numbers keep changing but the string doesnt so parse the text to get options.
            1. Change default production dynamics time if no previous input was given and running binding free energy simulations.
            2. Copy prmfilepath (default amoebabio18) to current directory (so dont have to have full path in tinker key file for parameters).
            3. If user inputs for stopping simulations, then call function to stop simulations (generates .end files in appropriate folders) then quit poltype.
            4. Copy ligandxyzfilenamelist to annihilateligandxyzfilenamelist if annihilateligandxyzfilenamelist is not specified by default.
            5. Read in tinker XYZ files in ligandxyzfilenamelist and convert to SMILES strings to be later used for detection in input PDB files.
            6. Create log file for simulations/FEP
            7. Change tinker executable names to reflect what is in the PATH (sometimes people use analyze.x vs analyze etc)
            8. Determine if should use GPU executable names (rather than CPU executables) based on if GPU executable is in PATH.
            9. Determine whether to submit jobs locally or not based on if keywords are enabled such as submitlocally and externalapi.
            10. Remove old intermediate files from previous run such as *.xyz_* etc..)
            11. Check if using a tinker version that is up to date enough to be consistent with poltype
            12. If the appropriate inputs are given for aligning a ligand to template ligand, then perform alignment in pocket and quit poltype. 
            13. Initialize variables to be filled in later in the program stack. Read water/Ions parameter type numbers from prmfilepath.
            14. If receptorligandxyzfilename is given as input, determine ligand indices for complexation box from input file ( and what the indices should be for solvation box as well).
            15. If there are overlapping types in input ligandxyzfilenamelist/keyfilenamelist then shift the types so that they are no longer overlapping. Do the same for receptorligandxyzfilename if this is provided as input.
            16. If pdbcode for crystal is given, download and fill in missing residues via Modeller and quit program.
            17. If user has inputs to use pdb2pqr to protonate protein (via proPka) then call pdb2pqr and quit program.
            18. Append all keys from keyfilenamelist into one key file. Check to ensure net charge is an integer, if not then pick some multipoles on host (protein) and add residual charge onto the protein multipole to ensure net charge is an integer.  
            19. If inputs for host-guest tinker XYZ are not given, then generate from complexed PDB.
            20. If alchemical=False (user doesnt wasnt to run FEP, just dynamics) then set internal lambda schedules to only include one step.
            21. If user wants neat liquid simulation, adjust defaults for equilibriation and production dynamics length. Only incllude dynamics for 100% electrostatics and 100% vdW.
            22. Check if need solvation gas phase regrowth for hydration free energy computations. Some molecules are small enough so that the vdW interactions dont exist. If small enough dont need gas phase vdW regrowth step.
            23. If user gives complexed host-guest tinker XYZ file and box info is on the top of the file, remove the box info.
            24. Ensure no undefined parameters for all key files, and make sure net charge of receptor is an integer, if not add residual charge to one of receptor multipoles and will append modified multipoles to complexation tinker key file. 
            25. If need to add gas phase vdW regrowth, adjust variables for electrostatic and vdW lambda schedules. 
            26. Set up variables such as system total charge, number of ions needed, filenames, variables used for free energy, enthalpy, entropy computations.
            27. If pathtosims is specified, then combine simulation table data from other simulation output tables (different simulation folders) and add into combined table. Generate combined plots as well.
            28. If perturbkeyfilelist is specified, then set up variables for determining filenames, folders and number of interpolations needed later.
            29. Append all keys from keyfilenamelist into keys for solvation/complexation. Append any additional modified multipoles needed to ensure net charge of receptor is an integer. 
            30. Setup keyfile names to be used later on.
            31. Compute step number for dynamics based on dynamic time and time step. Setup box filenames and output files for each dynamic step.
            32. Setup folder names for production dynamics and production dynamics outputfile names. 
            33. Setup production dynamic ARC file paths (for BAR input) and BAR output filepath names.
            34. Rotate into inertial frame. This will make it easier to compute longest length for determining box size later.
            35. Write relevant variables to internal dictionary to be later printed out in CSV files.
            36. Grab properties of ligand XYZ files to later be used when writing out PDB trajectories for equilbriated ARC and production dynamics ARC.
            37. Call annihilator.py (which calls minimization, equilbritation, production dynamics and BAR etc..)
            """
            # STEP 0
            self.DetermineXYZEditOptions()
            # STEP 1
            if self.neatliquidsim==True: 
                self.usepreequilibriatedbox=False # then make new neat liquid box from scratch
            if self.inputproddyntime==True:
                self.compproddyntime=self.proddyntime
                self.solvproddyntime=self.proddyntime

           
            # STEP 2

            head,tail=os.path.split(self.prmfilepath)
            newpath=os.path.join(os.getcwd(),tail)
            if self.prmfilepath!=newpath or head!=None:
                checkpath=os.path.join(os.getcwd(),self.prmfilepath)
                if checkpath!=newpath:
                    shutil.copy(self.prmfilepath,newpath)
            self.prmfilepath=tail
            # STEP 3

            if self.simulationstostopfolderpath!=None:
                self.StopSimulations(self.simulationstostopfolderpath)
                sys.exit()
            # STEP 4

            if self.annihilateligandxyzfilenamelist==None and self.ligandxyzfilenamelist!=None:
                self.annihilateligandxyzfilenamelist=self.ligandxyzfilenamelist.copy()
            if int(self.estatlambdascheme[0])!=1 or int(self.vdwlambdascheme[0])!=1 and self.solvation==True:
                raise ValueError('Please start with 1 on left side of lambdascheme and 0 on right. This way for solvation can add extra steps correctly')
            # STEP 5

            obConversion = openbabel.OBConversion()
            self.ligandsmileslist=[]
            self.annihilateligandsmileslist=[]
            self.ligandxyztosmiles={}
            if self.ligandxyzfilenamelist!=None:
                for xyz in self.ligandxyzfilenamelist:
                    ligmol = openbabel.OBMol()
                    obConversion.SetInFormat('xyz')
                    newname=self.ConvertTinkerXYZToCartesianXYZ(xyz)
                    obConversion.ReadFile(ligmol, newname)
                    obConversion.SetOutFormat('mol')
                    molname=xyz.replace('.xyz','.mol')
                    obConversion.WriteFile(ligmol,molname)
                    ligm=Chem.MolFromMolFile(molname,removeHs=False,sanitize=False)
                    smi=Chem.MolToSmarts(ligm)
                    smi=smi.replace('@','')
                    self.ligandsmileslist.append(smi)
                    self.ligandxyztosmiles[xyz]=smi
                    if xyz in self.annihilateligandxyzfilenamelist:
                        self.annihilateligandsmileslist.append(smi)

            # STEP 6

            self.logname=os.path.basename(os.getcwd())+'_'+self.logname 
            self.outputpath=os.path.join(os.getcwd(),'')
            self.logfh=open(self.outputpath+self.logname,'a+')
            # STEP 7

            self.SanitizeMMExecutables()
            # STEP 8

            foundgpukey=False
            if (self.which(self.dynamicommpath)):
                self.usegpu=True
                foundgpukey=True
            if self.externalapi!=None:
                self.runjobslocally=False 
            # STEP 9

            if self.externalapi==None and self.submitlocally==None:
                self.submitlocally=True
            if self.submitlocally!=True:
                self.submitlocally=False
            if not (self.which(self.dynamicommpath)) and foundgpukey==False:
                self.localdynamicpath=self.dynamicpath
                self.localbarpath=self.barpath
            else:
                self.localdynamicpath=self.dynamicommpath
                self.localbarpath=self.barommpath
            if self.usegpu==True:
                self.truedynamicpath=self.dynamicommpath
                self.truebarpath=self.barommpath
                self.trueanalyzepath=self.analyzeommpath
                self.trueminimizepath=self.minimizeommpath
            else:
                self.truedynamicpath=self.dynamicpath
                self.truebarpath=self.barpath
                self.trueanalyzepath=self.analyzepath
                self.trueminimizepath=self.minimizepath
            # STEP 10

            self.CleanUpFiles()
            # STEP 11

            self.CheckTinkerVersion()


            # STEP 13

            if self.ligandxyzfilenamelist!=None and self.templateligandxyzfilename!=None:
                self.AlignLigandXYZToTemplateXYZ()
                sys.exit()

            if self.alignPDB==True:
                self.AlignComplexedPDBToUnComplexedPDB()
                sys.exit()


            # STEP 14

            self.simfoldname=''
            self.boxfilename='box.xyz'
            self.masterdict={}
            self.masterdict['energy']={}
            self.masterdict['freeenergy']={}
            self.masterdict['summary']={}
            self.masterdict['lambda']={}

            if self.complexation==True and self.solvation==False:
                self.ligandindices=[[]]
                self.allligandindices=[[]]
            elif self.solvation==True and self.complexation==False:
                self.ligandindices=[[]]
                self.allligandindices=[[]]
            elif self.solvation==True and self.complexation==True:
                self.ligandindices=[[],[]]
                self.allligandindices=[[],[]]
            elif self.solvation==False and self.complexation==False and self.neatliquidsim==True:
                self.ligandindices=[[]]
                self.allligandindices=[[]]

            self.elementsymtotinktype={'K':box.GrabTypeNumber(self,'Potassium Ion K+'),'Cl':box.GrabTypeNumber(self,'Chloride Ion Cl-'),'Mg':box.GrabTypeNumber(self,'Magnesium Ion Mg+2'),'Li':box.GrabTypeNumber(self,'Lithium Ion Li+'),'Na':box.GrabTypeNumber(self,'Sodium Ion Na+'),'Rb':box.GrabTypeNumber(self,'Rubidium Ion Rb+'),'Cs':box.GrabTypeNumber(self,'Cesium Ion Cs+'),'Be':box.GrabTypeNumber(self,'Beryllium Ion Be+2'),'Ca':box.GrabTypeNumber(self,'Calcium Ion Ca+2'),'Zn':box.GrabTypeNumber(self,'Zinc Ion Zn+2'),'Mg+':box.GrabTypeNumber(self,'Magnesium Ion Mg+2')}
            self.ReadWaterFromPRMFile()

            # STEP 15


            if self.receptorligandxyzfilename!=None:
                oldtypelist,ligatomnums=self.GrabLigandTypesInfoForIndicesExtraction(self.ligandxyzfilenamelist)
                indextooldtype,complexligands,solvligands=self.ExtractLigandIndicesFromComplexXYZ(self.receptorligandxyzfilename,oldtypelist,ligatomnums)
                self.allligands=[]
                self.allligands.extend(complexligands)
                self.allligands.extend(solvligands)
                oldtypelist,ligatomnums=self.GrabLigandTypesInfoForIndicesExtraction(self.annihilateligandxyzfilenamelist)
                indextooldtype,complexligands,solvligands=self.ExtractLigandIndicesFromComplexXYZ(self.receptorligandxyzfilename,oldtypelist,ligatomnums)
                self.ligands=[]
                self.ligands.extend(complexligands)
                self.ligands.extend(solvligands)


            # STEP 16


            needtoshifttypes=self.CheckIfNeedToShiftTypes(self.keyfilenamelist)
            if needtoshifttypes==True:
                oldtypetonewtypelist=self.GenerateTypeMaps(self.keyfilenamelist)
                for i in range(len(oldtypetonewtypelist)):
                    oldindextonewindex=oldtypetonewtypelist[i]
                    key=self.keyfilenamelist[i]
                    xyz=self.ligandxyzfilenamelist[i]
                    self.ShiftParameterTypes(key,oldindextonewindex)
                    self.ShiftParameterTypes(xyz,oldindextonewindex)
                    
                if self.receptorligandxyzfilename!=None:
                    indextonewtype,complexligands,solvligands=self.ExtractLigandIndicesFromComplexXYZ(self.receptorligandxyzfilename,oldtypetonewtypelist,ligatomnums)
                    self.ShiftParameterTypesComplexXYZ(self.receptorligandxyzfilename,indextonewtype)

            # STEP 17
              

            if self.pdbcode!=None:
                pdbxyz.FillInMissingResidues(self,self.pdbcode)
                sys.exit()

            # STEP 18


            if self.usepdb2pqr==True and self.uncomplexedproteinpdbname!=None:
                pdbxyz.CallPDB2PQR(self,self.uncomplexedproteinpdbname)
                sys.exit()


            # STEP 19

            self.originalkeyfilename=self.AppendKeys(self.keyfilenamelist,'ligands.key')
            mpolearrays=self.CheckInputXYZKeyFiles(ligandonly=True) # check xyz/keys of ligands before making complex XYZ


            # STEP 20

            if self.complexation==True and self.complexedproteinpdbname!=None: 

                pdbxyz.GenerateProteinTinkerXYZFile(self)

            elif self.complexation==True  and self.complexedproteinpdbname==None and self.receptorligandxyzfilename==None: 
                raise ValueError('Missing complexedproteinpdbname, need ligand in complex')


            # STEP 21


            if self.alchemical==False:
                self.estatlambdascheme=[1]
                self.vdwlambdascheme=[1]
                self.restlambdascheme=[0]


            # STEP 22

            self.addgas=False
            if self.neatliquidsim==True:
                self.estatlambdascheme=[1]
                self.vdwlambdascheme=[1]
                self.restlambdascheme=[0]
                self.restrainatomsduringminimization=False
                if self.inputproddyntime==False:
                    self.proddyntime=3
                if self.inputequiltimeNVT==False:
                    self.equiltimeNVT=2
                if self.inputequiltimeNPT==False:
                    self.equiltimeNPT=.5
                if self.inputlastNVTequiltime==False:
                    self.lastNVTequiltime=.5
                self.addgas=True
            self.originalestatlambdascheme=self.estatlambdascheme[:]
            self.originalvdwlambdascheme=self.vdwlambdascheme[:]
            self.originalrestlambdascheme=self.restlambdascheme[:]
            self.originaltorsionrestscheme=self.torsionrestscheme[:]

            

            # STEP 23


            self.checkneedregrow=self.CheckIfNeedRegrowForSolvation()

            # STEP 24

            self.CheckReceptorLigandXYZFile() # sometimes user have box info on top
            
            # STEP 25
            mpolearrays=self.CheckInputXYZKeyFiles()

            # STEP 26

            if (self.extractinterforbinding==True and self.binding==True) or self.checkneedregrow==True:
                self.addgas=True
            if self.addgas:
                regrowelelambdascheme=self.estatlambdascheme[::-1]
                regrowvdwlambdascheme=self.vdwlambdascheme[::-1]
                regrowrestlambdascheme=self.restlambdascheme[::-1]
                self.estatlambdascheme=[self.estatlambdascheme,regrowelelambdascheme]
                self.vdwlambdascheme=[self.vdwlambdascheme,regrowvdwlambdascheme]
                self.restlambdascheme=[self.restlambdascheme,regrowrestlambdascheme]
                self.torsionrestscheme=[self.torsionrestscheme,self.torsionrestscheme]
            else:
                self.estatlambdascheme=[self.estatlambdascheme]
                self.vdwlambdascheme=[self.vdwlambdascheme]
                self.restlambdascheme=[self.restlambdascheme]
                self.torsionrestscheme=[self.torsionrestscheme]

            if self.complexationonly==True:
                self.solvation=False

            # STEP 27 
            if self.useuniquefilenames==True:
                head,self.foldername=os.path.split(os.getcwd())
            else:
                self.foldername='MD'
            if self.complexation==True and self.solvation==False:
                self.bufferlen=[20]

                self.proddynsteps=[str(int((self.compproddyntime*1000000)/self.proddyntimestep))]
                self.proddynframenum=[str(int(float(self.compproddyntime)/(float(self.proddynwritefreq)*0.001)))]
                self.proddyntimearray=[self.compproddyntime]
                self.systemcharge=[self.receptorcharge+self.ligandcharge+self.otherligcharge]
                self.ligandchargelist=[self.ligandcharge]
                self.xyzfilename=[[self.receptorligandxyzfilename]]
                self.iontypetoionnumberneut=[{}]
                self.iontypetoionnumberphysio=[{}]
                self.lambdafolderlist=[[]] 
                self.proddynoutfilepath=[[]]
                self.baroutputfilepath=[[]]
                self.barfilepath=[[]]
                self.thermooutputfilepath=[[]]
                self.secondarcpaths=[[]]
                self.firstarcpaths=[[]]
                self.tabledictkeysused=[[]]
                self.tabledict=[{}]
                self.simfoldname=['Comp'+self.simfoldname]
                self.boxfilename=[self.foldername+'_comp'+self.boxfilename]
                self.newkeyfilename=[self.foldername+'_comp'+'.key']
                self.WriteToLog('Complexation job')
                self.tightminoutput=[self.foldername+'_comp'+'tightmin.out']
                self.looseminoutput=[self.foldername+'_comp'+'loosemin.out']
                self.freeenergy=[[]]
                self.freeenergylist=[[]]
                self.overlaplist=[[]]
                self.freeenergyerrorlist=[[]]
                self.freeenergyfwd=[[]]
                self.freeenergylistfwd=[[]]
                self.freeenergyerrorlistfwd=[[]]
                self.freeenergybwd=[[]]
                self.freeenergylistbwd=[[]]
                self.freeenergyerrorlistbwd=[[]]
                self.freeenergyviabariter=[[]]
                self.freeenergylistviabariter=[[]]
                self.freeenergyerrorlistviabariter=[[]]
                self.freeenergyviabootstrap=[[]]
                self.freeenergylistviabootstrap=[[]]
                self.freeenergyerrorlistviabootstrap=[[]]
                self.enthalpy=[[]]
                self.enthalpylist=[[]]
                self.enthalpyerrorlist=[[]]
                self.enthalpyerrorlisttotal=[[]]
                self.entropy=[[]]
                self.entropylist=[[]]
                self.entropyerrorlist=[[]]
                self.entropyerrorlisttotal=[[]]
                self.freeenergyerror=[[]]
                self.enthalpyerror=[[]]
                self.entropyerror=[[]]
                self.masterdict['boxinfo']=[{}]

            elif self.solvation==True and self.complexation==False:
                self.proddynsteps=[str(int((self.solvproddyntime*1000000)/self.proddyntimestep))]
                self.proddynframenum=[str(int(float(self.solvproddyntime)/(float(self.proddynwritefreq)*0.001)))]
                self.proddyntimearray=[self.solvproddyntime]
                self.bufferlen=[2*float(self.vdwcutoff)+6]
                self.systemcharge=[self.ligandcharge]
                self.ligandchargelist=[self.ligandcharge]
                self.xyzfilename=[self.annihilateligandxyzfilenamelist]
                self.iontypetoionnumberneut=[{}]
                self.iontypetoionnumberphysio=[{}]
                self.lambdafolderlist=[[]]
                self.proddynoutfilepath=[[]]
                self.baroutputfilepath=[[]]
                self.barfilepath=[[]]
                self.thermooutputfilepath=[[]]
                self.secondarcpaths=[[]]
                self.firstarcpaths=[[]]
                self.tabledictkeysused=[[]]
                self.tabledict=[{}]
                self.simfoldname=['Solv'+self.simfoldname]
                self.boxfilename=[self.foldername+'_solv'+self.boxfilename]
                self.newkeyfilename=[self.foldername+'_solv'+'.key']
                self.WriteToLog('Solvation job')
                self.tightminoutput=[self.foldername+'_solv'+'tightmin.out']
                self.looseminoutput=[self.foldername+'_solv'+'loosemin.out']
                self.freeenergy=[[]]
                self.freeenergylist=[[]]
                self.overlaplist=[[]]
                self.freeenergyerrorlist=[[]]
                self.freeenergyfwd=[[]]
                self.freeenergylistfwd=[[]]
                self.freeenergyerrorlistfwd=[[]]
                self.freeenergybwd=[[]]
                self.freeenergylistbwd=[[]]
                self.freeenergyerrorlistbwd=[[]]
                self.freeenergyviabariter=[[]]
                self.freeenergylistviabariter=[[]]
                self.freeenergyerrorlistviabariter=[[]]
                self.freeenergyviabootstrap=[[]]
                self.freeenergylistviabootstrap=[[]]
                self.freeenergyerrorlistviabootstrap=[[]]
                self.enthalpy=[[]]
                self.enthalpylist=[[]]
                self.enthalpyerrorlist=[[]]
                self.enthalpyerrorlisttotal=[[]]
                self.entropy=[[]]
                self.entropylist=[[]]
                self.entropyerrorlist=[[]]
                self.entropyerrorlisttotal=[[]]
                self.freeenergyerror=[[]]
                self.enthalpyerror=[[]]
                self.entropyerror=[[]]
                self.masterdict['boxinfo']=[{}]


            elif self.solvation==False and self.complexation==False and self.neatliquidsim==True:
                self.proddynsteps=[str(int((self.proddyntime*1000000)/self.proddyntimestep))]
                self.proddynframenum=[str(int(float(self.proddyntime)/(float(self.proddynwritefreq)*0.001)))]
                self.proddyntimearray=[self.proddyntime]
                self.bufferlen=[2*float(self.vdwcutoff)+6]
                self.systemcharge=[self.ligandcharge]
                self.xyzfilename=[self.annihilateligandxyzfilenamelist]
                self.iontypetoionnumberneut=[{}]
                self.iontypetoionnumberphysio=[{}]
                self.lambdafolderlist=[[]]
                self.proddynoutfilepath=[[]]
                self.baroutputfilepath=[[]]
                self.barfilepath=[[]]
                self.thermooutputfilepath=[[]]
                self.secondarcpaths=[[]]
                self.firstarcpaths=[[]]
                self.tabledictkeysused=[[]]
                self.tabledict=[{}]
                self.simfoldname=['NeatLiq'+self.simfoldname]
                self.boxfilename=[self.foldername+'_neatliq'+self.boxfilename]
                self.newkeyfilename=[self.foldername+'_neatliq'+'.key']
                self.WriteToLog('Neat Liquid job')
                self.tightminoutput=[self.foldername+'_neatliq'+'tightmin.out']
                self.looseminoutput=[self.foldername+'_neatliq'+'loosemin.out']
                self.masterdict['boxinfo']=[{}]

            elif self.solvation==True and self.complexation==True:
                self.proddynsteps=[str(int((self.compproddyntime*1000000)/self.proddyntimestep)),str(int((self.solvproddyntime*1000000)/self.proddyntimestep))]
                self.proddynframenum=[str(int(float(self.compproddyntime)/(float(self.proddynwritefreq)*0.001))),str(int(float(self.solvproddyntime)/(float(self.proddynwritefreq)*0.001)))]
                self.proddyntimearray=[self.compproddyntime,self.solvproddyntime]
                self.bufferlen=[20,2*float(self.vdwcutoff)+6]
                self.systemcharge=[self.receptorcharge+self.ligandcharge+self.otherligcharge,self.ligandcharge]
                self.ligandchargelist=[self.ligandcharge,self.ligandcharge]
                self.xyzfilename=[[self.receptorligandxyzfilename],self.annihilateligandxyzfilenamelist]
                self.iontypetoionnumberneut=[{},{}]
                self.iontypetoionnumberphysio=[{},{}]
                self.lambdafolderlist=[[],[]]
                self.proddynoutfilepath=[[],[]]
                self.baroutputfilepath=[[],[]]
                self.barfilepath=[[],[]]
                self.thermooutputfilepath=[[],[]]
                self.secondarcpaths=[[],[]]
                self.firstarcpaths=[[],[]]
                self.tabledictkeysused=[[],[]]
                self.tabledict=[{},{}]
                self.simfoldname=['Comp'+self.simfoldname,'Solv'+self.simfoldname]
                self.boxfilename=[self.foldername+'_comp'+self.boxfilename,self.foldername+'_solv'+self.boxfilename]
                self.newkeyfilename=[self.foldername+'_comp'+'.key',self.foldername+'_solv'+'.key']
                self.WriteToLog('Binding job')
                self.tightminoutput=[self.foldername+'_comp'+'tightmin.out',self.foldername+'_solv'+'tightmin.out']
                self.looseminoutput=[self.foldername+'_comp'+'loosemin.out',self.foldername+'_solv'+'loosemin.out']
                self.freeenergy=[[],[]]
                self.freeenergylist=[[],[]]
                self.overlaplist=[[],[]]
                self.freeenergyerrorlist=[[],[]]
                self.freeenergyfwd=[[],[]]
                self.freeenergylistfwd=[[],[]]
                self.freeenergyerrorlistfwd=[[],[]]
                self.freeenergybwd=[[],[]]
                self.freeenergylistbwd=[[],[]]
                self.freeenergyerrorlistbwd=[[],[]]
                self.freeenergyviabariter=[[],[]]
                self.freeenergylistviabariter=[[],[]]
                self.freeenergyerrorlistviabariter=[[],[]]
                self.freeenergyviabootstrap=[[],[]]
                self.freeenergylistviabootstrap=[[],[]]
                self.freeenergyerrorlistviabootstrap=[[],[]]
                self.enthalpy=[[],[]]
                self.enthalpylist=[[],[]]
                self.enthalpyerrorlist=[[],[]]
                self.enthalpyerrorlisttotal=[[],[]]
                self.entropy=[[],[]]
                self.entropylist=[[],[]]
                self.entropyerrorlist=[[],[]]
                self.entropyerrorlisttotal=[[],[]]
                self.freeenergyerror=[[],[]]
                self.enthalpyerror=[[],[]]
                self.entropyerror=[[],[]]
                self.masterdict['boxinfo']=[{},{}]

            # STEP 28

            if self.pathtosims!=None:
                tables.GrabSimDataFromPathList(self)
                plots.PlotEnergyData(self)
                sys.exit()

            # STEP 29

            if self.perturbkeyfilelist!=None:
                ls=[i+1 for i in range(len(self.perturbkeyfilelist))]
                self.perturbkeyfilelisttokeyindex=dict(zip(self.perturbkeyfilelist,ls)) 
                if self.numinterpolsforperturbkeyfiles==None:
                    self.numinterpolsforperturbkeyfiles=[1 for i in range(len(self.perturbkeyfilelist))]
                else:
                    self.numinterpolsforperturbkeyfiles=[int(i) for i in self.numinterpolsforperturbkeyfiles]
                self.perturbkeyfilelisttointerpolations=dict(zip(self.perturbkeyfilelist,self.numinterpolsforperturbkeyfiles))
                if self.fep==True: # if fep keyword is used then dont need to run dynamics only BAR, so submit jobs on current node
                    self.submitlocally=True

            # STEP 30

            for i in range(len(self.newkeyfilename)):
                newkeyfilename=self.newkeyfilename[i]
                self.AppendKeys(self.keyfilenamelist,newkeyfilename)
                if i==0:
                    self.ModifyCharge(newkeyfilename,mpolearrays)

            # STEP 31

            if self.complexation==True and self.solvation==False:
                self.keyfilename=[[self.foldername+'_comp'+'.key']]
                if self.perturbkeyfilelist!=None:
                    for keyname,index in self.perturbkeyfilelisttokeyindex.items():
                        numinterpols=self.perturbkeyfilelisttointerpolations[keyname]
                        for k in range(numinterpols):
                            newkeyname=self.foldername+'_'+str(k)+'_comp'+'.key'
                            shutil.copy(keyname,newkeyname)
                            self.keyfilename[0].append(newkeyname)
                confignamelist=[i.replace('.key','_config.key') for i in self.keyfilename[0]]
                lambdanamelist=[i.replace('.key','_lambda.key') for i in self.keyfilename[0]]
                self.configkeyfilename=[confignamelist]
                self.lambdakeyfilename=[lambdanamelist]

            elif self.solvation==True and self.complexation==False:
                self.keyfilename=[[self.foldername+'_solv'+'.key']]
                if self.perturbkeyfilelist!=None:
                    for keyname,index in self.perturbkeyfilelisttokeyindex.items():
                        numinterpols=self.perturbkeyfilelisttointerpolations[keyname]
                        for k in range(numinterpols):
                            newkeyname=self.foldername+'_'+str(k)+'_solv'+'.key'
                            shutil.copy(keyname,newkeyname)
                            self.keyfilename[0].append(newkeyname)
                confignamelist=[i.replace('.key','_config.key') for i in self.keyfilename[0]]
                lambdanamelist=[i.replace('.key','_lambda.key') for i in self.keyfilename[0]]
                self.configkeyfilename=[confignamelist]
                self.lambdakeyfilename=[lambdanamelist]


            elif self.solvation==False and self.complexation==False and self.neatliquidsim==True:
                self.keyfilename=[[self.foldername+'_neatliq'+'.key']]
                if self.perturbkeyfilelist!=None:
                    for keyname,index in self.perturbkeyfilelisttokeyindex.items():
                        numinterpols=self.perturbkeyfilelisttointerpolations[keyname]
                        for k in range(numinterpols):
                            newkeyname=self.foldername+'_'+str(k)+'_neatliq'+'.key'
                            shutil.copy(keyname,newkeyname)
                            self.keyfilename[0].append(newkeyname)
                confignamelist=[i.replace('.key','_config.key') for i in self.keyfilename[0]]
                lambdanamelist=[i.replace('.key','_lambda.key') for i in self.keyfilename[0]]
                self.configkeyfilename=[confignamelist]
                self.lambdakeyfilename=[lambdanamelist]


            elif self.solvation==True and self.complexation==True:
                self.keyfilename=[[self.foldername+'_comp'+'.key'],[self.foldername+'_solv'+'.key']]
                if self.perturbkeyfilelist!=None:
                    for keyname,index in self.perturbkeyfilelisttokeyindex.items():
                        numinterpols=self.perturbkeyfilelisttointerpolations[keyname]
                        for k in range(numinterpols):
                            newkeyname=self.foldername+'_'+str(k)+'_comp'+'.key'
                            shutil.copy(keyname,newkeyname)
                            self.keyfilename[0].append(newkeyname)
                            newkeyname=self.foldername+'_solv'+'.key'
                            shutil.copy(keyname,newkeyname)
                            self.keyfilename[1].append(newkeyname)
                confignamelistcomp=[i.replace('.key','_config.key') for i in self.keyfilename[0]]
                lambdanamelistcomp=[i.replace('.key','_lambda.key') for i in self.keyfilename[0]]
                confignamelistsolv=[i.replace('.key','_config.key') for i in self.keyfilename[1]]
                lambdanamelistsolv=[i.replace('.key','_lambda.key') for i in self.keyfilename[1]]
                self.configkeyfilename=[confignamelistcomp,confignamelistsolv]
                self.lambdakeyfilename=[lambdanamelistcomp,lambdanamelistsolv]

            # STEP 32 
            
            self.equilstepsNVT=str(int((self.equiltimeNVT*1000000)/self.equiltimestep/len(self.equilibriatescheme)-1)) # convert ns to fs divide by the length of temperature scheme
            self.lastequilstepsNVT=str(int((self.lastNVTequiltime*1000000)/self.equiltimestep))
            self.equilstepsNPT=str(int((self.equiltimeNPT*1000000)/self.equiltimestep))
            self.equilstepsionNPT=str(int((self.equiltimeionNPT*1000000)/self.equiltimestep))
            self.equilframenumNPT=int((float(self.equiltimeNPT))/(float(self.equilwritefreq)*0.001))
            self.equilframenumNVT=int((float(self.equiltimeNVT+self.lastNVTequiltime))/(float(self.equilwritefreq)*0.001))
            self.equilframenum=self.equilframenumNPT+self.equilframenumNVT
            self.minboxfilename=[i.replace('.xyz','min.xyz') for i in self.boxfilename]
            self.minboxfilenamepymol=[i.replace('.xyz','_pymol.xyz') for i in self.minboxfilename]
            self.equilboxfilename=[i.replace('.xyz','equil.xyz') for i in self.boxfilename]
            self.proddynboxfilename=[i.replace('.xyz','proddyn.xyz') for i in self.boxfilename]
            self.proddynboxfilenamepymol=[i.replace('.xyz','_pymol.xyz') for i in self.proddynboxfilename]
            self.equilarcboxfilename=[i.replace('.xyz','.arc') for i in self.equilboxfilename]
            self.proddynarcboxfilename=[i.replace('.xyz','.arc') for i in self.proddynboxfilename]
            self.proddynboxkeyfilename=[i.replace('.xyz','.key') for i in self.proddynboxfilename]
            self.equiljobsfilename=self.foldername+'_equiljobs.txt'
            self.proddynjobsfilename=self.foldername+'_proddynamicsjobs.txt'
            self.barjobsfilename=self.foldername+'_barjobs.txt'
            self.freeenergyjobsfilename=self.foldername+'_freeenergyjobs.txt'
            self.tightminjobsfilename=self.foldername+'_tightboxminjobs.txt'
            self.looseminjobsfilename=self.foldername+'_looseboxminjobs.txt'
            self.analyzejobsfilename=self.foldername+'_analyzejobs.txt'
            self.solviontocount=self.DetermineIonsForChargedSolvation()
            if len(self.solviontocount.keys())!=0 and self.salthfe==True: # then need to add counter ions free energy sims
                self.estatlambdascheme.append(self.originalestatlambdascheme)
                self.vdwlambdascheme.append(self.originalvdwlambdascheme)
                self.restlambdascheme.append(self.originalrestlambdascheme)
                self.ionboxfilename='ionbox.xyz'
                self.ionequilboxfilename='ionequilbox.xyz'
                self.ionproddynboxfilename='ionproddynbox.xyz'
                self.ionequilarcfilename=self.ionequilboxfilename.replace('.xyz','.arc')
                self.ionkeyfilename='ion'+self.configkeyfilename[0][0]
                self.systemcharge.append(self.ligandcharge)

            self.equiloutputarray=[]
            self.equiloutputstepsarray=[]
            for i in range(len(self.simfoldname)):
                simfold=self.simfoldname[i]
                templist=[]
                tempsteps=[]
                for temp in self.equilibriatescheme:
                    templist.append(self.outputpath+self.foldername+'_'+simfold+'_'+str(temp)+'_'+self.equilstepsNVT+'_NVT.out')
                    tempsteps.append(self.equilstepsNVT)

                templist.append(self.outputpath+self.foldername+'_'+simfold+'_'+str(temp)+'_'+self.equilstepsNPT+'_NPT.out')
                tempsteps.append(self.equilstepsNPT)

                if i==0 and self.solvation==True and self.addsolvionwindows==True:
                    templist.append(self.outputpath+'SolvIon'+'_'+simfold+'_'+str(temp)+'_'+self.equilstepsionNPT+'_NPT.out')
                    tempsteps.append(self.equilstepsionNPT)
                self.equiloutputarray.append(templist)
                self.equiloutputstepsarray.append(tempsteps)


            # STEP 33


            self.nextfiletofinish=self.equiloutputarray[0]
            if self.productiondynamicsNVT==True:
                self.proddynensem=self.NPVTensem
            else:
                self.proddynensem=self.NPTensem
            self.foldernametolambdakeyfilename={}
            self.foldernametonuminterpols={}
            self.foldernametointerpolindex={}
            for simfoldidx in range(len(self.simfoldname)):
                simfold=self.simfoldname[simfoldidx]
                templambdafolderlistoflist=[]
                tempproddynoutfilepathlistoflist=[]
                lambdakeyfilelist=self.lambdakeyfilename[simfoldidx]
                for k in range(len(self.vdwlambdascheme)):
                    vdwlambdascheme=self.vdwlambdascheme[k]
                    estatlambdascheme=self.estatlambdascheme[k]
                    restlambdascheme=self.restlambdascheme[k]
                    tempproddynoutfilepath=[]
                    templambdafolderlist=[]
                    for i in range(len(vdwlambdascheme)):
                        elelamb=estatlambdascheme[i]
                        vdwlamb=vdwlambdascheme[i]
                        rest=restlambdascheme[i]
                        lambdakeyfilename=lambdakeyfilelist[0]
                        if 'Comp' in simfold:
                            if k==0:
                                fold=simfold+"E%s_V%s_R%s"%(elelamb,vdwlamb,rest)
                            elif k==1 and self.addgas==True:
                                fold=simfold+"E%s_V%s_R%s_Gas"%(elelamb,vdwlamb,rest)
                        else:
                            if k==0:
                                fold=simfold+"E%s_V%s"%(elelamb,vdwlamb)
                            elif k==1 and self.addgas==True:
                                fold=simfold+"E%s_V%s_Gas"%(elelamb,vdwlamb)
                            elif k==2 and self.addgas==True:
                                fold=simfold+"E%s_V%s_Ion"%(elelamb,vdwlamb)
                            elif k==1 and self.addgas==False:
                                fold=simfold+"E%s_V%s_Ion"%(elelamb,vdwlamb)

                        self.foldernametolambdakeyfilename[fold]=lambdakeyfilename
                        templambdafolderlist.append(fold)
                        outputfilepath=os.getcwd()+'/'+simfold+'/'+fold+'/'
                        outputfilepath+=self.foldername+'_'+fold+'.out'
                        if k==1 and (self.extractinterforbinding==True and self.binding==True):
                            pass
                        else:
                            tempproddynoutfilepath.append(outputfilepath)
                    
                    if self.perturbkeyfilelist!=None:
                        if k==0:
                            elelamb=estatlambdascheme[0]
                            vdwlamb=vdwlambdascheme[0]
                            rest=restlambdascheme[0]
                        elif k==1 and self.addgas==True:
                            elelamb=estatlambdascheme[-1]
                            vdwlamb=vdwlambdascheme[-1]
                            rest=restlambdascheme[-1]
                        elif k==1 and self.addgas==False:
                            elelamb=estatlambdascheme[0]
                            vdwlamb=vdwlambdascheme[0]
                            rest=restlambdascheme[0]
                        elif k==2 and self.addgas==True:
                            elelamb=estatlambdascheme[0]
                            vdwlamb=vdwlambdascheme[0]
                            rest=restlambdascheme[0]

                        count=1
                        for keyfilename,index in self.perturbkeyfilelisttokeyindex.items():
                            numinterpols=self.perturbkeyfilelisttointerpolations[keyfilename]
                            for r in range(numinterpols):
                                lambdakeyfilename=lambdakeyfilelist[count]
                                split=keyfilename.split('.')
                                keyprefix=split[0]
                                count+=1
                                interpolindex=r+1
                                interpolrate=int(((interpolindex/numinterpols)*100))
                                if 'Comp' in simfold:
                                   fold=simfold+"E%s_V%s_R%s_Key_%s_Ipl%s"%(elelamb,vdwlamb,rest,str(keyprefix),str(interpolrate))
                                else:
                                   if k==0:
                                       fold=simfold+"E%s_V%s_Key_%s_Ipl%s"%(elelamb,vdwlamb,str(keyprefix),str(interpolrate))
                                   elif k==1:
                                       fold=simfold+"E%s_V%s_Gas_Key_%s_Ipl%s"%(elelamb,vdwlamb,str(keyprefix),str(interpolrate))
                                   else:
                                       continue
                                self.foldernametolambdakeyfilename[fold]=lambdakeyfilename
                                self.foldernametonuminterpols[fold]=numinterpols
                                self.foldernametointerpolindex[fold]=interpolindex
                                if k==0:
                                    self.estatlambdascheme[k].insert(0,elelamb)
                                    self.vdwlambdascheme[k].insert(0,vdwlamb)
                                    self.restlambdascheme[k].insert(0,rest)
                                    templambdafolderlist.insert(0,fold) 
                                elif k==1 and self.addgas==True:
                                    self.estatlambdascheme[k].append(elelamb)
                                    self.vdwlambdascheme[k].append(vdwlamb)
                                    self.restlambdascheme[k].append(rest)
                                    templambdafolderlist.append(fold) 
                                elif k==1 and self.addgas==False:
                                    self.estatlambdascheme[k].insert(0,elelamb)
                                    self.vdwlambdascheme[k].insert(0,vdwlamb)
                                    self.restlambdascheme[k].insert(0,rest)
                                    templambdafolderlist.insert(0,fold) 
                                elif k==2 and self.addgas==True:
                                    self.estatlambdascheme[k].insert(0,elelamb)
                                    self.vdwlambdascheme[k].insert(0,vdwlamb)
                                    self.restlambdascheme[k].insert(0,rest)
                                    templambdafolderlist.insert(0,fold) 


                                if self.fep==False:
                                    outputfilepath=os.getcwd()+'/'+simfold+'/'+fold+'/'
                                    outputfilepath+=self.foldername+'_'+fold+'.out'
                                    tempproddynoutfilepath.append(outputfilepath)


                    templambdafolderlistoflist.append(templambdafolderlist)
                    tempproddynoutfilepathlistoflist.append(tempproddynoutfilepath)
                self.proddynoutfilepath[simfoldidx]=tempproddynoutfilepathlistoflist
                self.lambdafolderlist[simfoldidx]=templambdafolderlistoflist
                
            
           
            # STEP 34

            
 
            for k in range(len(self.lambdafolderlist)): # stop before last one because using i+1 for grabbing next index
                listoffolderlist=self.lambdafolderlist[k]
                foldname=self.simfoldname[k]
                proddynarcboxfilename=self.proddynarcboxfilename[k]
                secondarcpathslistoflist=[] 
                firstarcpathslistoflist=[] 
                barfilepathlistoflist=[]
                baroutputfilepathlistoflist=[]
                thermooutputfilepathlistoflist=[]
                for j in range(len(listoffolderlist)):
                    folderlist=listoffolderlist[j]
                    secondarcpathslist=[] 
                    firstarcpathslist=[] 
                    barfilepathlist=[]
                    baroutputfilepathlist=[]
                    thermooutputfilepathlist=[]

                    for i in range(len(folderlist)-1):
                        firstfoldname=folderlist[i]
                        secondfoldname=folderlist[i+1]

                        firstarcpath=self.outputpath+foldname+r'/'+firstfoldname+'/'+proddynarcboxfilename.replace(self.foldername+'_',self.foldername+'_'+firstfoldname+'_').replace('0.','0-')


                        secondarcpath=self.outputpath+foldname+r'/'+secondfoldname+'/'+proddynarcboxfilename.replace(self.foldername+'_',self.foldername+'_'+secondfoldname+'_').replace('0.','0-')


                        secondarcpathslist.append(secondarcpath)
                        firstarcpathslist.append(firstarcpath)


                        baroutputfilepath=self.outputpath+foldname+r'/'+secondfoldname+'/'+self.foldername+'_'+firstfoldname+secondfoldname+'BAR1.out'
                        thermooutputfilepath=baroutputfilepath.replace('BAR1.out','BAR2.out')
                        barfilepath=baroutputfilepath.replace('BAR1.out','.bar')
                        barfilepathlist.append(barfilepath)
                        baroutputfilepathlist.append(baroutputfilepath)
                        thermooutputfilepathlist.append(thermooutputfilepath)
                    secondarcpathslistoflist.append(secondarcpathslist)
                    firstarcpathslistoflist.append(firstarcpathslist)
                    barfilepathlistoflist.append(barfilepathlist)
                    baroutputfilepathlistoflist.append(baroutputfilepathlist)
                    thermooutputfilepathlistoflist.append(thermooutputfilepathlist)
                self.secondarcpaths[k]=secondarcpathslistoflist
                self.firstarcpaths[k]=firstarcpathslistoflist
                self.barfilepath[k]=barfilepathlistoflist
                self.baroutputfilepath[k]=baroutputfilepathlistoflist
                self.thermooutputfilepath[k]=thermooutputfilepathlistoflist
            
            # STEP 35

            self.RotateTranslateInertialFrame()

            # STEP 36

            for i in range(len(self.tabledict)):
                self.tabledict[i]['Prod MD Ensemb']=self.proddynensem
                self.tabledict[i]['Solv Prod MD Time']=self.solvproddyntime
                self.tabledict[i]['Dynamic Writeout Frequency (ps)']=self.proddynwritefreq
                self.tabledict[i]['Dynamic Time Step (fs)']=self.equiltimestep
                self.tabledict[i]['Equil Time NPT']=self.equiltimeNPT
                self.tabledict[i]['Equil Time NVT']=self.equiltimeNVT
                self.tabledict[i]['Ligand Charge']=self.ligandcharge
                for idx in range(len(self.annihilateligandxyzfilenamelist)):
                    xyz=self.annihilateligandxyzfilenamelist[idx]
                    self.tabledict[i]['Ligand Name '+str(idx+1)]=xyz.replace('.xyz','')
                if self.expfreeenergy!=None:
                    self.tabledict[i][u'ΔGᵉˣᵖ']=self.expfreeenergy
                if self.complexation==True:
                    self.tabledict[i]['Comp Prod MD Time']=self.compproddyntime
                    self.tabledict[i]['Comp Prod MD Steps']=self.proddynsteps[0]
                    if self.complexationonly==False: 
                        self.tabledict[i]['Solv Prod MD Steps']=self.proddynsteps[1]
                    self.tabledict[i]['Receptor Charge']=self.receptorcharge
                    if self.uncomplexedproteinpdbname!=None:
                        self.tabledict[i]['Receptor Name']=self.uncomplexedproteinpdbname.replace('.pdb','')
                    elif self.receptorligandxyzfilename!=None and self.uncomplexedproteinpdbname==None: 
                        self.tabledict[i]['Receptor Name']=self.receptorligandxyzfilename.replace('.xyz','')
                else:
                    self.tabledict[i]['Solv Prod MD Steps']=self.proddynsteps[0]

            
            tables.WriteTableUpdateToLog(self)

            # STEP 37

            self.ligandindextoneighbslist=[]
            self.ligandindextosymlist=[]
            self.ligandindextotypenumlist=[]

            for xyz in self.annihilateligandxyzfilenamelist:
                statexyzatominfo,stateindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=self.GrabXYZInfo(xyz)
                self.ligandindextoneighbslist.append(indextoneighbs)
                self.ligandindextosymlist.append(indextosym)
                self.ligandindextotypenumlist.append(stateindextotypeindex)

            # STEP 38
            ann.main(self)





        def CheckIfNeedToShiftTypes(self,keyfilenamelist):
            """
            Intent: Sometimes keyfiles may have overlapping types but user wants to treat two ligands differently (dissapear only one etc). So then need to check if type numbers are overlapping and change them if so.
            Input: List of keyfile names
            Output: Boolean determining if need to shift type numbers or not.
            Referenced By: MolecularDynamics 
            Description:
            1. If keyfilenamelist is given as input.
            2. Iterate over list of keyfilenames and determine the max and min type number.
            3. Determine the range of types for the max and min type numbers. 
            4. Iterate over the types and save new types to an array. If encounter existing type already in array then need to shift type numbers.  
            """

            needtoshifttypes=False
            # STEP 1
            if keyfilenamelist!=None: 
                newlist=keyfilenamelist.copy()
                #newlist.append(self.prmfilepath) can only go to 1000 for atom classes in tinker so for now dont shift by max of prm file
                alltypes=[]
                # STEP 2
                for keyfilename in newlist:
                    maxnumberfromkey=self.GrabMaxTypeNumber(keyfilename)
                    minnumberfromkey=self.GrabMinTypeNumber(keyfilename)
                    # STEP 3
                    types=list(range(minnumberfromkey,maxnumberfromkey+1))
                    # STEP 4
                    for typenum in types:
                        if typenum in alltypes:
                            needtoshifttypes=True
                            break
                        alltypes.append(typenum)


            return needtoshifttypes



        def AppendKeys(self,keyfilenamelist,thekey):
            """
            Intent: Append all keys from keyfilenamelist into a single key file.
            Input: List of key files and output keyfilename.
            Output: Final appended keyfile.
            Referenced By: MolecularDynamics 
            Description: 
            1. For each key in keyfilenamelist, read the lines into an array. 
            2. For each line in array, append to final keyfile.
            """
            temp=open(thekey,'w')
            # STEP 1
            for key in keyfilenamelist:
                othertemp=open(key,'r')
                results=othertemp.readlines()
                othertemp.close()
                # STEP 2
                for line in results:
                    temp.write(line)

            temp.close()
            return thekey

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




        def ReturnListOfList(self,string):
            """
            Intent: Some input keywords require a list of list datatype. Example onlyrottortorlist=b1 c1 d1 , b2 c2 d2 . Outer list is comma seperated and inner list is space seperated.
            Input: String from reading poltype.ini
            Output: List of list
            Referenced By: _post_init
            Description:
            1. Split string by comma into a list.
            2. For each item in list, strip front and end extra spaces. Then split by space.
            3. For each item in spaced list, append to inner list.
            4. Append inner list to outerlist.
            """
            # STEP 1
            newlist=string.split(',')
            templist=[]
            for ele in newlist:
                # STEP 2
                nums=ele.lstrip().rstrip().split()
                temp=[]
                # STEP 3
                for e in nums:
                    temp.append(int(e))
                # STEP 4
                templist.append(temp)
            newlist=templist
            return newlist



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
            Referenced By: _post_init, MolecularDynamics 
            Description:
            1. Call SanitizeMMExecutable for all tinker executables used. 
            """
            self.peditexe=self.SanitizeMMExecutable(self.peditexe)
            self.potentialexe=self.SanitizeMMExecutable(self.potentialexe)
            self.minimizeexe=self.SanitizeMMExecutable(self.minimizeexe)
            self.analyzeexe=self.SanitizeMMExecutable(self.analyzeexe)
            self.superposeexe=self.SanitizeMMExecutable(self.superposeexe)
            self.xyzeditpath=self.SanitizeMMExecutable(self.xyzeditpath)
            self.barpath=self.SanitizeMMExecutable(self.barpath)
            self.dynamicpath=self.SanitizeMMExecutable(self.dynamicpath)
            self.minimizepath=self.SanitizeMMExecutable(self.minimizepath)
            self.pdbxyzpath=self.SanitizeMMExecutable(self.pdbxyzpath)
            self.analyzepath=self.SanitizeMMExecutable(self.analyzepath)
            self.potentialpath=self.SanitizeMMExecutable(self.potentialpath)
            self.poleditpath=self.SanitizeMMExecutable(self.poleditpath)

        

        

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
            cmdstr=self.analyzeexe+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/VersionFiles/'+'water.xyz'+' '+'-k'+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/VersionFiles/'+'water.key'+' '+'e'+'>'+' '+'version.out'
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
                raise ValueError("Notice: Not latest working version of tinker (8.10.2)"+' '+os.getcwd())
           
            # STEP 4
            if self.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
                self.paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebaplus21_header.prm"
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
        
        def printfile(self,filename):
            """
            Intent: Print filename contents
            Input: filename
            Output: filename contents written to standard output.
            Referenced By: copyright 
            Description:
            1. Open filename and print file contents
            """
            # STEP 1
            with open(os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/'+filename,'r') as f:
                print(f.read(), end='')
       

        def copyright(self):
            """
            Intent: Print version file and copyright info
            Input: Self object
            Output: Contents of version file printed to standard output.
            Referenced By: GenerateParameters() 
            Description: 
            1. Call printfile on poltype version file
            """
            # STEP 1
            self.printfile(self.versionfile)
        
            
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
                        if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line or "LBFGS  --  Normal Termination due to SmallGrad" in line or "Normal termination" in line or 'Normal Termination' in line or 'Total Potential Energy' in line or 'Psi4 stopped on' in line or 'finished run' in line:
                            term=True
                        if ('Tinker is Unable to Continue' in line or 'error' in line or ' Error ' in line or ' ERROR ' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation, address not mapped to object' in line or 'galloc:  could not allocate memory' in line or 'Erroneous write.' in line or 'Optimization failed to converge!' in line) and 'DIIS' not in line and 'mpi' not in line and 'RMS Error' not in line:
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
                            self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                            raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                    else:
                        if exitcode != 0:
                            self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())

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
                        self.WriteToLog("torsion parameters are all zero "+line+' path ='+os.getcwd())
                        raise ValueError("torsion parameters are all zero "+line+' path ='+os.getcwd())

                        

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
            cmdstr = f"python {os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'lConformerGenerator.py')} {self.molstructfname}"
            print('Calling: '+cmdstr) 
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
            inFormat = obConversion.FormatFromExt(self.molstructfname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, self.molstructfname)
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
            Output: Boolean specifying if should use PCM or not.
            Referenced By: GenerateParameters 
            Description: 
            1. If the molecule contains hydrogens, then continue, else dont use PCM
            2. Iterate over dictionary of atomic index to formal charges and save atomic indices that have formal charges to an array.
            3. Iterate over array of charged atomic indices
            4. Iterate over neighbors of each charged atomic index, if any neighbor also contains a formal charge, then need to use PCM
            5. Iterate over neighbors of the first neighbors, if atomindex is not original charged index but has formal charge then need to use PCM.
            6. Iterate over neighbors of the second neighbors, if atomindex is not the first neighbor but has formal charge then need to use PCM.

            """
            chargedindices=[]
            pcm=False
            # STEP 1
            if self.hashyd==True:
                # STEP 2
                for atomindex,chg in atomindextoformalcharge.items():
                    if chg!=0:
                        chargedindices.append(atomindex)
                # STEP 3
                for atomindex in chargedindices:
                    atom=m.GetAtomWithIdx(atomindex)
                    # STEP 4
                    for atm in atom.GetNeighbors():
                        atmidx=atm.GetIdx()
                        if atmidx in chargedindices:
                            pcm=True
                        # STEP 5
                        for natm in atm.GetNeighbors():
                            natmidx=natm.GetIdx()
                            if natmidx!=atomindex:
                                if natmidx in chargedindices:
                                    pcm=True
                                # STEP 6
                                for nnatm in natm.GetNeighbors():
                                    nnatmidx=nnatm.GetIdx()
                                    if nnatmidx!=atmidx:
                                        if nnatmidx in chargedindices:
                                            pcm=True
            return pcm 
                    
                        
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
                        if 'AMOEBA-opt' not in f and 'ANI-opt' not in f and 'xtb-opt' not in f and 'log' not in end and 'sdf' not in end and 'ini' not in end and 'chk' not in end and 'dat' not in end and 'mol' not in end and 'txt' not in end and f not in filestonotdelete:
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
            listofallowedatoms=[1,6,7,8,15,16,17,35,53,9]
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
            57. If user specifies to send email about job being finished, then send email
            58. Copy final XYZ and key files to top directory, where user submits jobs. Also copy torsion plots to OPENME folder.
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
            pcm=self.CheckForConcentratedFormalCharges(m,atomindextoformalcharge)
            self.pcm=False
            if pcm==True and self.dontusepcm==False:
                self.pcm=True
                self.toroptpcm=True
                self.optpcm=True
                self.torsppcm=True

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
                self.optmethod='wB97X-D'
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
                bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,torsionsmissing,classkeytotorsionparametersguess,missingvdwatomindextoneighbors,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo,tortorsmissing=databaseparser.GrabSmallMoleculeAMOEBAParameters(self,optmol,mol,m)
            if os.path.isfile(self.torsionsmissingfilename):
                torsionsmissing=databaseparser.ReadTorsionList(self,self.torsionsmissingfilename)
            if os.path.isfile(self.torsionprmguessfilename):
                classkeytotorsionparametersguess=databaseparser.ReadDictionaryFromFile(self,self.torsionprmguessfilename)
            if os.path.isfile(self.vdwmissingfilename):
                missingvdwatomindices=databaseparser.ReadVdwList(self,self.vdwmissingfilename)
            if self.onlyvdwatomlist!=None:
                missingvdwatomindices=self.onlyvdwatomlist[:]
            if os.path.isfile(self.tortormissingfilename):
                tortorsmissing=databaseparser.ReadTorTorList(self,self.tortormissingfilename)

            # STEP 27
            if (self.writeoutpolarize==True and self.writeoutmultipole==True):
                esp.SPForDMA(self,optmol,mol)
            

                # STEP 28 
                if not os.path.isfile(self.gdmafname):
                    mpole.run_gdma(self)
        
                # STEP 29
                polarindextopolarizeprm,polartypetotransferinfo=databaseparser.GrabSmallMoleculeAMOEBAParameters(self,optmol,mol,m,polarize=True)
                # STEP 30
                lfzerox=mpole.gen_peditinfile(self,mol,polarindextopolarizeprm)
            
                if (not os.path.isfile(self.xyzfname) or not os.path.isfile(self.keyfname)):
                    # STEP 31
                    cmdstr = self.peditexe + " 1 " + self.gdmafname +' '+self.paramhead+ " < " + self.peditinfile
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
                    mpole.AddPolarizeCommentsToKey(self,self.key2fnamefromavg,polartypetotransferinfo)
                # STEP 38
                fit=False
                if self.espfit and not os.path.isfile(self.key3fname) and self.atomnum!=1:
                    xyzfnamelist,keyfnamelist=self.GenerateDuplicateXYZsFromOPTs(self.xyzoutfile,self.key2fnamefromavg,optmolist)   
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
                databaseparser.appendtofile(self,self.key3fname,self.key4fname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo)
                if self.writeoutangle==True:
                    databaseparser.StiffenZThenBisectorAngleConstants(self,self.key4fname)
                    databaseparser.TestBondAngleEquilValues(self)
                self.AddIndicesToKey(self.key4fname)
                if self.databasematchonly==True:
                    sys.exit()
            # STEP 41
            (torlist, self.rotbndlist,nonaroringtorlist,self.nonrotbndlist) = torgen.get_torlist(self,optmol,torsionsmissing,self.onlyrotbndslist)
            torlist,self.rotbndlist=torgen.RemoveDuplicateRotatableBondTypes(self,torlist) # this only happens in very symmetrical molecules
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
            self.WriteToLog('Poltype Job Finished'+'\n')
            # STEP 55
            if self.isfragjob==False and self.dontfrag==True and len(self.torlist)!=0:
                self.WriteOutDatabaseParameterLines()
            # STEP 56
            if os.path.exists(self.scrtmpdirgau):
                shutil.rmtree(self.scrtmpdirgau)
            if os.path.exists(self.scrtmpdirpsi4):
                shutil.rmtree(self.scrtmpdirpsi4)
            # STEP 57
            if self.email!=None:
                moleculename=self.molstructfname.replace('.sdf','')
                password='amoebaisbest'
                fromaddr = 'poltypecrashreportnoreply@gmail.com'
                toaddr = self.email
                TEXT='Molecule has finished parameterization'
                try: 
                    self.SendReportEmail(TEXT,fromaddr,toaddr,password,moleculename,'Poltype Finished Report ')
                except:
                    pass
            # STEP 58
            if self.isfragjob==False:
                previousdir=os.path.abspath(os.path.join(os.getcwd(), os.pardir))
                if os.path.exists(self.tmpxyzfile):
                     shutil.copy(self.tmpxyzfile,os.path.join(previousdir,self.tmpxyzfile))
                if os.path.exists(self.tmpkeyfile):
                     shutil.copy(self.tmpkeyfile,os.path.join(previousdir,self.tmpkeyfile))
            self.CopyFitPlots()
            os.chdir('..')
            


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

        
        def SendReportEmail(self,TEXT,fromaddr,toaddr,password,filename,subject):
            """
            Intent: If user gives email as input and poltype crashes, then send email report that crashed. 
            Input: Email text, from address, to address, password and filename that has error.
            Output: -
            Referenced By: RunPoltype
            Description:
            1. Initialize email object
            2. Enter from address, to address, subject and message into email object 
            3. Initialize server object that contacts gmail, enter from address and password, then login
            4. Send message and close server object
            """
            # STEP 1
            msg = MIMEMultipart()
            # STEP 2
            msg['From'] = fromaddr
            msg['To'] = toaddr
            msg['Subject'] = subject+filename
            message = TEXT
            msg.attach(MIMEText(message, 'plain'))
            # STEP 3
            s = smtplib.SMTP_SSL('smtp.gmail.com')
            s.ehlo()
            s.login(fromaddr,password)
            text = msg.as_string()
            # STEP 4
            s.sendmail(fromaddr, [toaddr],text)
            s.quit()


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

        def StopSimulations(self,simulationstostopfolderpath):
            """
            Intent: If user wants to stop all simulations, then go into each folder and add a *.end file, tinker will then terminate dynamics safely without corrupting arc or dyn files accidentally. 
            Input: Folder path with subdirectories containing simulation folders (with arcs etc...)
            Output: - 
            Referenced By: MolecularDynamics
            Description: 
            1. Walk over all subdirectories from input filepath
            2. If .arc file is detected within current folder or subfolder, add a file with .end into the folder 
            """
            # STEP 1
            for root, subdirs, files in os.walk(simulationstostopfolderpath):
                for d in subdirs:
                    curdir=os.getcwd()
                    path=os.path.join(root, d)
                    os.chdir(path)
                    files=os.listdir()
                    for f in files:
                        # STEP 2
                        if '.arc' in f:
                            split=f.split('.')
                            prefix=split[0]
                            endname=prefix+'.end'
                            endname=os.path.join(path,endname)
                            with open(endname, 'w') as fp:
                                pass
                        # STEP 2
                        if os.path.isdir(f):
                            os.chdir(f)
                            newfiles=os.listdir()
                            for newf in newfiles:
                                if '.arc' in newf:
                                    split=newf.split('.')
                                    prefix=split[0]
                                    endname=prefix+'.end'
                                    with open(endname, 'w') as fp:
                                        pass
                        os.chdir('..')
                    os.chdir(curdir)


        def ReadWaterFromPRMFile(self):
            """
            Intent: Grab water parameters from prmfile, if user specifices older prmfile, the type numbers may be different, so need to parse via tinker type description rather than hard coded type numbers.
            Input: -
            Output: -
            Referenced By: MolecularDynamics 
            Description: 
            1. Read lines of prmfilepath into an array
            2. Iterate over lines in array and if detect line with water hydrogen or water oxygen description, then save the type numbers to variable names for later use
            """
            # STEP 1
            temp=open(self.prmfilepath,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                # STEP 2
                if 'atom' in line:
                    linesplit=line.split()
                    oxygen='"AMOEBA Water O"'
                    hydrogen='"AMOEBA Water H"'
                    typenum=int(linesplit[1])
                    if oxygen in line:
                        self.waterOtypenum=typenum
                    elif hydrogen in line:
                        self.waterHtypenum=typenum

        def RotateTranslateInertialFrame(self):
            """
            Intent: Rotate into inertial frame to make it easier to compute appropriate dimensions for box based on max length of the molecule. 
            Input: - 
            Output: - 
            Referenced By: MolecularDynamics 
            Description: 
            1. Create filename with inputs for xyzedit program
            2. Call xyzedit with input file created
            3. Rename output filename (which may have _2 etc) after to the original intended filename  
            """
            if self.receptorligandxyzfilename!=None:
                # STEP 1
                filename='xyzedit.in'
                temp=open(filename,'w')
                temp.write('14'+'\n')
                temp.write('\n')
                temp.close()
                cmdstr=self.xyzeditpath+' '+self.receptorligandxyzfilename+' '+'-k'+' '+self.originalkeyfilename+' '+'<'+' '+filename
                # STEP 2 
                submit.call_subsystem(self,cmdstr,wait=True)   
                # STEP 3 
                endsplit=self.receptorligandxyzfilename.split('.')
                end=endsplit[1]
                if '_' in end:
                    newsplit=end.split('_')
                    count=int(newsplit[1])
                else:
                    count=1
                newcount=str(count+1)
                filename=self.receptorligandxyzfilename
                if filename[-2]=='_':
                    filename=filename[:-2]
                newfilename=filename+'_'+newcount
                os.rename(newfilename,self.receptorligandxyzfilename)


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

        def CheckReceptorLigandXYZFile(self):
            """
            Intent: If user has box information on the top in input file (copied from some arc etc), then remove it.
            Input: -
            Output: -
            Referenced By: MolecularDynamics 
            Description:
            1. Read lines of xyz file into array
            2. Iterate over lines of array, then check if there is any box information (6 items, all floats)
            3. Skip over box information and write all other lines to new file, then rename to original filename 
            """
            if self.receptorligandxyzfilename!=None:
                # STEP 1
                temp=open(self.receptorligandxyzfilename,'r')
                results=temp.readlines()
                temp.close()
                tempname=self.receptorligandxyzfilename.replace('.xyz','-t.xyz')
                temp=open(tempname,'w')
                for line in results:
                    linesplit=line.split() 
                    isboxline=True
                    # STEP 2
                    if len(linesplit)==6:
                        for e in linesplit:
                            if self.CheckFloat(e)==False:
                                isboxline=False
                    else:
                        isboxline=False
                    # STEP 3
                    if isboxline==True:
                        continue
                    else:
                        temp.write(line)
                temp.close()
                os.remove(self.receptorligandxyzfilename)
                os.rename(tempname,self.receptorligandxyzfilename)
                        
                         

 
        def CleanUpFiles(self):
            """
            Intent: Remove intermediate files from tinker computations (usually the last one is copied and/or renamed to desired filename)
            Input: - 
            Output: -
            Referenced By: MolecularDynamics
            Description: 
            1. Iterate over all files in directory, remove any .xyz_ files or key_ files 
            """
            files=os.listdir()
            # STEP 1
            for f in files:
                if self.ligandxyzfilenamelist!=None:
                    for xyz in self.ligandxyzfilenamelist:
                        if xyz+'_' in f:
                            os.remove(f)
                if 'water.xyz_' in f or 'key_' in f or 'uncomplexed.' in f:
                    os.remove(f)
                if self.receptorligandxyzfilename!=None:
                    if self.receptorligandxyzfilename+'_' in f:
                        os.remove(f)

        
 
            
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

        

        def ReadLigandOBMol(self,structfname):
            """
            Intent: Convert xyz->mol object to easily probe molecular properties. For example when converting tinker XYZ-> cartesian xyz, then want to convert cartesian xyz to openbabel mol object. 
            Input: Filename
            Output: openbabel mol object
            Referenced By: ReadLigandRdkitMol
            Description:
            1. Initialize openbabel conversion object
            2. Set input format based on file extension
            3. Create new openbabel mol object and read in filename to mol object 
            """
            # STEP 1
            tmpconv = openbabel.OBConversion()
            # STEP 2
            inFormat = openbabel.OBConversion.FormatFromExt(structfname)
            tmpconv.SetInFormat(inFormat)
            # STEP 3
            tmpmol = openbabel.OBMol()
            tmpconv.ReadFile(tmpmol, structfname)
            return tmpmol

        def ReadLigandRdkitMol(self,ligandfilename):
            """
            Intent: Convert xyz->mol object to easily probe molecular properties. For example when converting tinker XYZ-> cartesian xyz, then want to convert cartesian xyz to rdkit mol object. 
            Input: Filename
            Output: rdkit mol object
            Referenced By: Parameter comparison module 
            Description:
            1. First generate openbabel object
            2. Convert openbabel object to .mol file
            3. Read .mol file into rdkit object (they have more limited ways of reading in filetypes than openbabel) 
            """
            # STEP 1
            tmpmol=self.ReadLigandOBMol(ligandfilename)
            # STEP 2
            tmpconv = openbabel.OBConversion()
            tmpconv.SetOutFormat('mol')
            temp=ligandfilename.replace('.sdf','.mol')
            tmpconv.WriteFile(tmpmol, temp)
            # STEP 3
            mol=rdmolfiles.MolFromMolFile(temp,removeHs=False,sanitize=False)

            return mol


        def PDBCoordinate(self,coords,index):
            """
            Intent: Need to format PDB files very strictly based on correct column number. This formats coordinate string for PDB files exactly. 
            Input: Input coordinates, Index for coordinates (for x,y, or z coord)
            Output: String formatted for coordinate to be added to PDB file
            Referenced By: GeneratePDBFromARC
            Description: 
            1. Extract coordinate via index
            2. Convert coordinate to string
            3. Split string via decimal place
            4. Save the string after the decimal
            5. Determine the number of 0's that need to be appended after words based on length of the string after decimal
            6. Append the appropriate number of 0's after decimal place
            7. Determine the number of empty spaces that need to be appended to front of string based on the current total length of string (only 8 characters allowed max). 
            8. Append empty spaces to front of string based on length determined in last step. 
            9. Return string 
            """
            # STEP 1
            x=round(float(coords[index]),3)
            # STEP 2
            xstring=str(x)
            # STEP 3
            xstringsplit=xstring.split('.')
            # STEP 4
            xtail=xstringsplit[1]
            # STEP 5
            diff = 3-len(xtail)
            # STEP 6
            for i in range(1,diff+1):
                xstring+='0'
            # STEP 7
            odiff=8-len(xstring)
            # STEP 8
            space=''
            for i in range(1,odiff+1):
                space+=' '
            xstring=space+xstring
            # STEP 9
            return xstring
       
 
        def WritePDBString(self,string,line,firstindex,lastindex):
            """
            Intent: Need to input string exactly in right columns for line that will go into PDB file.
            Input: String to add to line, original line, first and last indices of where string will go into line
            Output: Modified line
            Referenced By: GeneratePDBFromARC
            Description: 
            1. Iterate over indices in between firstindex and last index
            2. For each character in input string, assign charcter to index in line 
            """
            # STEP 1
            for i in range(firstindex,lastindex+1):
                # STEP 2
                xitem=string[i-firstindex]
                line[i]=xitem

            return line


        def GrabPDBLines(self):
            """
            Intent: Generate dictionary of pdbindex to pdb line, for use in generating PDB trajectory from tinker ARC
            Input: -
            Output:
            Referenced By: GeneratePDBFromARC
            Description: 
            1. Iterate over lines of PDB file
            2. If ATOM or HETATM in line, grab index
            3. Store index and line into dictionary 
            """
            # STEP 1
            temp=open(self.complexedproteinpdbname,'r')
            results=temp.readlines()
            temp.close()
            pdbindextopdbline={}
            # STEP 2
            for line in results:
                if 'ATOM' in line or 'HETATM' in line:
                    index=int(line[6:10+1].strip())
                    # STEP 3
                    pdbindextopdbline[index]=line
            return pdbindextopdbline

        def WritePDBCONECTRecord(self,temp,indextoconn):
            """
            Intent: Add CONECT records to PDB file 
            Input: filehandle, dictionary of atomindex to connected atom indices
            Output: - 
            Referenced By: GeneratePDBFromARC
            Description: 
            1. Iterate over dictionary of index to connected indices
            2. Generate string by calling CONECTRecord
            3. Write line to filehandle 
            """
            # STEP 1
            for index,conns in indextoconn.items():
                if len(conns)!=0:
                    # STEP 2
                    string=self.CONECTRecord(conns,index)
                    # STEP 3
                    temp.write(string)
            temp.write('ENDMDL'+'\n')

        def CONECTRecord(self,conns,index):
            """
            Intent: Add CONECT records to PDB file 
            Input: filehandle, dictionary of atomindex to connected atom indices
            Output: - 
            Referenced By: GeneratePDBFromARC
            Description: 
            1. For each index, create a string and append empty spaces in front of string such that the total length is 5
            2. Write the index string to the CONECT line string
            3. For all connected indices, do the same, create string of length 5 and add in write place to CONECT string
            """

            shift=0
            string='CONECT                                                                       '
            string=list(string)
            # STEP 1
            indexstring=str(index)
            odiff=5-len(indexstring)
            space=''
            for i in range(1,odiff+1):
                space+=' '
            indexstring=space+indexstring
            # STEP 2
            self.WritePDBString(indexstring,string,6+shift,10+shift)
            shift+=5
            # STEP 3
            for idx in conns:
                indexstring=str(idx)
                odiff=5-len(indexstring)
                space=''
                for i in range(1,odiff+1):
                    space+=' '
                indexstring=space+indexstring
                self.WritePDBString(indexstring,string,6+shift,10+shift)
                shift+=5
            string=''.join(string)+'\n'
            return string


        def GeneratePDBFromARC(self,tinkerxyzfilename,pdbfilename):
            """
            Intent: Convert tinker XYZ or tinker ARC into PDB trajectory
            Input: tinker file, desired pdbfilename
            Output: -
            Referenced By: PymolReadableFile
            Description: 
            1. Iterate over lines of tinker file and save in array
            2. Grab dictionary of PDB index to pdb line from pdbfilename
            3. Save generic HETATM string for use in converting lines in tinker XYZ that were not originally in the input PDB (such as waters or ions etc).
            4. Create dictionary of tinker XYZ index to PDB index
            5. Append current model number to top of PDB file
            6. Iterate over tinker XYZ lines 
            7. Extact current element, XYZ index, and connected XYZ indices and pdbline from dictionary of tinker xyz index to PDB line
            8. Save connected indices in dictionry for later
            9. Add index formatted corectly to PDB line
            10. Add element, formatted correctly to PDB line
            11. Add residue, formatted correctly to PDB line 
            12. Add tinker coordinates, formatted correctly to PDB line
            13. Write connect record and then update model number and write model line to PDB file everytime box information is detected in tinker file
            14. Write final connect record and ENDMDL after looping over all tinker XYZs in ARC  
            """
            # STEP 1
            temp=open(tinkerxyzfilename,'r')
            results=temp.readlines()
            temp.close()
            # STEP 2
            pdbindextopdbline=self.GrabPDBLines()
            temp=open(pdbfilename,'w')
            # STEP 3 
            hetatmstring='HETATM22011  H   UNK  1334      -6.447 -12.617  22.452'+'\n'
            # STEP 4
            xyzindextopdbindex={v: k for k, v in self.pdbindextoxyzindex.items()}
            indextoconn={}
            # STEP 5
            first=True
            modelnum=1
            modelstring='MODEL                                            '+'\n'
            modelstring=list(modelstring)
            modelnumstring=str(modelnum)
            odiff=4-len(modelnumstring)
            space=''
            for i in range(1,odiff+1):
                space+=' '
            modelnumstring=space+modelnumstring
            self.WritePDBString(modelnumstring,modelstring,10,13)
            modelstring=''.join(modelstring)
            temp.write(modelstring)
            # STEP 6
            for line in results:
                linesplit=line.split()
                cont=False
                if '90.000' not in line and len(linesplit)>1:
                    first=False
                    # STEP 7
                    xyzindex=int(linesplit[0])
                    element=linesplit[1]
                    conns=linesplit[6:]
                    conns=[int(i) for i in conns]
                    if xyzindex in xyzindextopdbindex.keys():
                        pdbindex=xyzindextopdbindex[xyzindex]
                        pdbline=pdbindextopdbline[pdbindex]
                        residuenum=int(pdbline[22:25+1].strip())
                        lastres=residuenum
                    else:
                        if self.removesolventpdbtraj==True:
                            cont=True
                            continue
                        pdbline=hetatmstring
                        residuenum=lastres+1
                    if cont==False:
                        # STEP 8
                        indextoconn[xyzindex]=conns
                        # STEP 9
                        fullsplit=re.split(r'(\s+)', pdbline)
                        pdbline=list(pdbline)
                        indexstring=str(xyzindex)
                        odiff=5-len(indexstring)
                        space=''
                        for i in range(1,odiff+1):
                            space+=' '
                        indexstring=space+indexstring
                        self.WritePDBString(indexstring,pdbline,6,10)
                        # STEP 10
                        elementstring=str(element)
                        odiff=4-len(elementstring)
                        space=''
                        for i in range(1,odiff+1):
                            space+=' '
                        elementstring=elementstring+space
                        self.WritePDBString(elementstring,pdbline,12,15)
                        # STEP 11
                        residuestring=str(residuenum)
                        odiff=4-len(residuestring)
                        space=''
                        for i in range(1,odiff+1):
                            space+=' '
                        residuestring=space+residuestring
                        self.WritePDBString(residuestring,pdbline,22,25)
                        # STEP 12
                        coords=[float(linesplit[2]),float(linesplit[3]),float(linesplit[4])]        
                        xstring=self.PDBCoordinate(coords,0)
                        ystring=self.PDBCoordinate(coords,1)
                        zstring=self.PDBCoordinate(coords,2)
                        self.WritePDBString(xstring,pdbline,30,37)
                        self.WritePDBString(ystring,pdbline,38,45)
                        self.WritePDBString(zstring,pdbline,46,53)
                        pdbline=''.join(pdbline)
                        temp.write(pdbline)
                elif '90.000' in line and first==False:
                    # STEP 13
                    modelnum+=1
                    self.WritePDBCONECTRecord(temp,indextoconn)
                    modelstring='MODEL                                            '+'\n'
                    modelstring=list(modelstring)
                    modelnumstring=str(modelnum)
                    odiff=4-len(modelnumstring)
                    space=''
                    for i in range(1,odiff+1):
                        space+=' '
                    modelnumstring=space+modelnumstring
                    self.WritePDBString(modelnumstring,modelstring,10,13)
                    modelstring=''.join(modelstring)
                    temp.write(modelstring)
            # STEP 14
            self.WritePDBCONECTRecord(temp,indextoconn)
            temp.write('ENDMDL'+'\n')
            temp.close()


        def PymolReadableFile(self,tinkerxyzfilename,outputname):
            """
            Intent: Convert tinker XYZ/ARC -> PDB trajectory, so can visualize in popular visualization software
            Input: Tinker xyz filename, output name for pymol readable file
            Output:
            Referenced By: Various python modules for when converting tinker XYZ-> PDB, such as minimized XYZ, final frame of equilibriated XYZ and even ARC for production dynamics with full electrostatics and van der waals interactions.  
            Description: 
            1. Convert tinker XYZ to cartesian XYZ
            2. If complexed XYZ/ARC file then generate PDB trajectory
            """
            # STEP 1
            xyzfilename=self.ConvertTinkerXYZToCartesianXYZ(tinkerxyzfilename)
            os.rename(xyzfilename,outputname)
            pdbfilename=outputname.replace('.xyz','.pdb')
            # STEP 2
            if self.complexedproteinpdbname!=None and 'comp' in tinkerxyzfilename:
                self.GeneratePDBFromARC(tinkerxyzfilename,pdbfilename)   


        def GrabIndexToCoordinatesPymol(self,mol):
            """
            Intent: Grab coordinates from mol object
            Input: rdkit mol object
            Output: Dictionary of atom index to coordinates
            Referenced By: pdbxyz module
            Description:
            1. Iterate over atoms in mol object
            2. Grab atom position from default conformer
            3. Save the coordinates into dictionary
            """
            indextocoords={}
            # STEP 1
            for atom in mol.GetAtoms():
                atomidx=atom.GetIdx()
                # STEP 2
                pos = mol.GetConformer(0).GetAtomPosition(atomidx)
                X=pos.x
                Y=pos.y
                Z=pos.z
                atomindex=atomidx+1 
                # STEP 3
                indextocoords[atomindex]=[X,Y,Z]
        
            return indextocoords

            
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


        def CheckIfNeedRegrowForSolvation(self):
            """
            Intent: If molecule is small enough, vdw interactions will not be needed (need at least 4 atoms across to feel vdw interactions). If small enough, then for hydration free energy, dont need to use vdw-annihilate or regrow vdw with gas phase. 
            Input: - 
            Output: Boolean to see if need to regrow vdw interactions in gas phase for HFE
            Referenced By: MolecularDynamics
            Description: 
            1. Assume do not need to regrow
            2. Iterate over neighors
            3. Iterate over neighbors of neighbors
            4. Iterate over neighbors of neighbors of neighbors
            5. Iterate over neighbors of previous neighbors again...
            6. If atom index is not repeated, then molecule is big enough to require gas phase vdw regrowth for HFE  
            """
            # STEP 1
            needregrow=False
            if self.solvation==True:
                if self.binding==False and self.annihilatevdw==True:
                    for xyz in self.ligandxyzfilenamelist:
                        xyzfilename=self.ConvertTinkerXYZToCartesianXYZ(xyz)
                        tmpmol=self.ReadLigandOBMol(xyzfilename)
                        # STEP 2
                        atomiter=openbabel.OBMolAtomIter(tmpmol)
                        for atom in atomiter:
                            atomidx=atom.GetIdx()
                            iteratomatom = openbabel.OBAtomAtomIter(atom)
                            # STEP 3
                            for natom in iteratomatom:
                                natomidx=natom.GetIdx()
                                iteratomatomatom = openbabel.OBAtomAtomIter(natom)
                                # STEP 4
                                for nnatom in iteratomatomatom:
                                    nnatomidx=nnatom.GetIdx()
                                    if nnatomidx!=atomidx: 
                                        # STEP 5
                                        iteratomatomatomatom = openbabel.OBAtomAtomIter(nnatom)
                                        for nnnatom in iteratomatomatomatom:
                                            nnnatomidx=nnnatom.GetIdx()
                                            # STEP 6
                                            if nnnatomidx!=natomidx:
                                                needregrow=True


            return needregrow 


        def CheckTinkerVersion(self):
            """
            Intent: Need to ensure user doesnt have a tinker version that is too old and incompatible with poltype
            Input: -
            Output: -
            Referenced By: MolecularDynamics
            Description: 
            1. Use water.xyz and water.key in VersionFiles folder to construct a command for the analyze program
            2. Call analyze with command
            3. Read output file and parse for version number
            4. If the version number is less than required value, then just crash poltype 
            """
            # STEP 1
            cmdstr=self.analyzepath+' '+os.path.abspath(os.path.join(self.annihilatorpath, os.pardir))+r'/VersionFiles/'+'water.xyz'+' '+'-k'+' '+os.path.abspath(os.path.join(self.annihilatorpath, os.pardir))+r'/VersionFiles/'+'water.key'+' '+'e'+'>'+' '+'version.out'
            
            self.WriteToLog(cmdstr)
            # STEP 2
            try:
                returned_value = subprocess.call(cmdstr, shell=True)
            except:
                raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())      
            # STEP 3
            temp=open('version.out','r')
            results=temp.readlines()
            temp.close()
            latestversion = False
            for line in results:
                if "Version" in line:
                    linesplit = line.split()
                    self.versionnum = linesplit[2]
                    if version.parse(self.versionnum) >= version.parse("8.9.4"):
                        latestversion = True
                        break
            # STEP 4
            if not latestversion:
                raise ValueError("Notice: Not latest working version of tinker (8.9.4)"+' '+os.getcwd())
          

        def GrabIndicesWithTypeNumber(self,xyzfilename,ligandtypes):
            """
            Intent: Find all atom indices corresponding to input type numbers, used when comparing types between input ligand tinker XYZ files as well as receptor-ligand tinker XYZ files (type numbers may be inconsistent or out of order etc). 
            Input: Tinker xyz filename, array of type numbers
            Output: array of atomic indices
            Referenced By: CheckInputXYZKeyFiles
            Description:
            1. Iterate over lines of tinker XYZ file
            2. Determine if line contains box information
            3. If not a box line, extract type number and index, if type number is in input array of type numbers then append index to array of indices. 
            """
            indices=[]
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            # STEP 1
            for line in results:
                linesplit=line.split()
                # STEP 2
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if len(linesplit)>1 and isboxline==False:
                    # STEP 3
                    typenum=int(linesplit[5])
                    index=int(linesplit[0])
                    if typenum in ligandtypes:
                        indices.append(index)
            return indices
            
        def GrabTypeNumbers(self,xyzfilename,indices=None):  
            """
            Intent: Grab all type numbers from tinker XYZ file
            Input: tinker xyz filename
            Output: array of type numbers
            Referenced By: CheckInputXYZKeyFiles
            Description: 
            1. Iterate over lines of tinker XYZ file
            2. If line does not contain box information
            3. Extract type number and append to array
            """
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            typenums=[]
            # STEP 1
            for line in results:
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                # STEP 2
                if len(linesplit)>1 and isboxline==False:
                    # STEP 3
                    typenum=int(linesplit[5])
                    index=int(linesplit[0])
                    if indices!=None:
                        if index in indices:
                            typenums.append(typenum)
                    else:
                        typenums.append(typenum)
            return typenums 

        def GrabCoordinates(self,xyzfilename,indices=None):
            """
            Intent: Want to grab coordinates from complexed XYZ, want to make sure solvation XYZ has same initial coordinates.
            Input: Tinker xyz filename, optionally input array of indices
            Output: Array of coordinates
            Referenced By: CheckInputXYZKeyFiles, minimization module - FindTwoAtomsForRestrainingRotation
            Description:
            1. Iterate over tinker XYZ
            2. If not box line, then grab coordinates
            3. Append coordinates to array of coordinates
            """
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            coords=[]
            # STEP 1
            for line in results:
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                # STEP 2
                if len(linesplit)>1 and isboxline==False:
                    typenum=int(linesplit[5])
                    index=int(linesplit[0])
                    # STEP 3
                    coord=[float(linesplit[2]),float(linesplit[3]),float(linesplit[4])]
                    try:
                        if index in indices:
                            coords.append(coord)
                    except:
                        coords.append(coord)
            return coords


        def CompareTypes(self,receptorligandtypes,ligandtypes):
            """
            Intent: Ensure type numbers in receptor-ligand XYZ are same as from input ligand tinker XYZ file.
            Input: Array of receptorligand types, array of ligand types
            Output: -
            Referenced By: CheckInputXYZKeyFiles
            Description:
            1. Iterate over ligandtypes
            2. Extract corresponding receptorligand type
            3. If the types are not the same, raise error and crash program
            """
            for i in range(len(ligandtypes)):
                receptorligandtype=receptorligandtypes[i]
                ligandtype=ligandtypes[i]
                if receptorligandtype!=ligandtype:
                    raise ValueError('Ligand XYZ and Receptor-Ligand XYZ dont have the same type number order!')

        def RewriteCoordinates(self,xyzfilename,coords):
            """
            Intent: Update coordinates in tinker XYZ with new coordinates, for example if want solvation XYZ to have same coordinates as ligand in complexed-ligand tinker XYZ
            Input: Tinker xyz filename, array of coordinates
            Output:
            Referenced By: CheckInputXYZKeyFiles
            Description: 
            1. Iterate over tinker XYZ
            2. If not box line, then grab coordinates from input array of coordinates
            3. Add coordinates to line
            4. Write out new line with updated coordinates
            """
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            tempname=xyzfilename.replace('.xyz','-t.xyz')
            temp=open(tempname,'w')
            count=0
            # STEP 1
            for line in results:
                linesplit=line.split()
                # STEP 2
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                # STEP 3
                if len(linesplit)>1 and isboxline==False:
                    coord=coords[count]
                    count+=1
                    linesplit[2]=str(coord[0])
                    linesplit[3]=str(coord[1])
                    linesplit[4]=str(coord[2])
                    line=' '.join(linesplit)+'\n'
                if isboxline==False:
                    # STEP 4
                    temp.write(line)
            temp.close()
            os.remove(xyzfilename)
            os.rename(tempname,xyzfilename)


        def ReadCharge(self,output,checkresid=True):
            """
            Intent: Read charge from output using tinker analyze. For example, want to ensure total charge is integer or for some systems need to ensure net charge is 0.
            Input: Output from tinker analyze, boolean specifying to check for residual charge and crash if not integer
            Output: Total charge and residual charge (left over from nearest integer).
            Referenced By: CheckNetChargeIsZero,CheckInputXYZKeyFiles
            Description: 
            1. If boolean specifies to skip checking charge, then dont crash program (for debugging purposes) 
            2. Read in results from analyze output
            3. Iterate over results and extract total charge 
            4. Check for residual charge
            5. If boolean specifies to crash for residual charge and there is residual charge, crash program
            """
            # STEP 1
            if self.skipchargecheck==True:
                checkresid=False
            # STEP 2
            temp=open(output,'r')
            results=temp.readlines()
            temp.close()
            chg=0
            # STEP 3
            for line in results:
                if 'Total Electric Charge :' in line:
                    linesplit=line.split()
                    chg=float(linesplit[-2])
            # STEP 4
            resid=int(round(chg))-chg
            # STEP 5
            if checkresid==True:
                if resid!=0:
                    raise ValueError('Charge is not integer! '+output+' '+str(chg)+' '+'residual charge is '+str(resid))

            return chg,resid


        def CheckEnergies(self,output):
            """
            Intent: Used after minimizing box to determine if box setup or minimization did not perform well (bond energy per bond is too high)
            Input: Outputfile from Tinker analyze
            Output: -
            Referenced By: minimization.py - CheapMinimizationProtocol
            Description: 
            1. Define tolerance of 5 kcal/mol/bond
            2. Iterate over results of tinker analyze output
            3. Extract Bond Stretching Energies
            4. Divide Bond energies by total number of bonds
            5. If bond energy per bond is greater than tolerance, then something is definitely wrong with setup.  
            """
            temp=open(output,'r')
            results=temp.readlines()
            temp.close()
            # STEP 1
            tol=5
            # STEP 2
            for line in results:
                if 'Bond Stretching' in line and 'Parameters' not in line:
                    linesplit=line.split()
                    intnum=int(linesplit[-1])
                    # STEP 3
                    energy=float(linesplit[-2])
                    # STEP 4
                    energyperbond=energy/intnum
                    # STEP 5
                    if energyperbond>tol:
                        raise ValueError('Bad starting box structure, bond energy per bond is way to high '+str(energyperbond)+' kcal/mol')





        def CheckInputXYZKeyFiles(self,ligandonly=False):
            """
            Intent: Append all keys from keyfilenamelist into one key file. Check to ensure net charge is an integer, if not then pick some multipoles on host (protein) and add residual charge onto the protein multipole to ensure net charge is an integer.
            Input: Boolean to only check ligand XYZs/keys (want to check ligand tinker files before making complex XYZ and key) 
            Output: Array of multipole parameters that may need to be added if there is residual charge (not integer net charge) 
            Referenced By: MolecularDynamics 
            Description: 
            1. Remove typical header keywords from keyfile input. Sometimes users copy from other source with header information.
            2. Add parameter file key word to header 
            3. Iterate over all input ligand XYZ/keys 
            """
            keymods.RemoveKeyWords(self,self.originalkeyfilename,['parameters','axis','ewald','pme-grid','pme-order','cutoff','thermostat','integrator','ligand','verbose','archive','neighbor-list','polar-eps','polar-predict','heavy-hydrogen','omp-threads','OPENMP-THREADS'])
            head,tail=os.path.split(self.prmfilepath)
            if head=='':
                self.prmfilepath=os.path.join(os.getcwd(),self.prmfilepath)
            string='parameters '+self.prmfilepath+'\n'
            keymods.AddKeyWord(self,self.originalkeyfilename,string)
            if self.receptorligandxyzfilename!=None and self.ligandxyzfilenamelist!=None:
                for xyz in self.ligandxyzfilenamelist:
                    ligandtypes=self.GrabTypeNumbers(xyz) 
                    indices=self.GrabIndicesWithTypeNumber(self.receptorligandxyzfilename,ligandtypes)
                    receptorligandtypes=self.GrabTypeNumbers(self.receptorligandxyzfilename,indices=indices)
                    self.CompareTypes(receptorligandtypes,ligandtypes)
                receptortypes=list(range(1,348+1))
                receptorindices=self.GrabIndicesWithTypeNumber(self.receptorligandxyzfilename,receptortypes)
                self.totalreceptornumber=len(receptorindices)
                self.receptorindices=list(range(1,self.totalreceptornumber+1))
            if self.ligandxyzfilenamelist!=None and self.originalkeyfilename!=None: 
                self.otherligcharge=0
                self.ligandcharge=0
                for xyz in self.ligandxyzfilenamelist:
                    if xyz not in self.annihilateligandxyzfilenamelist: # more background charge
                        cmdstr=self.trueanalyzepath+' '+xyz+' '+'-k'+' '+self.originalkeyfilename+' '+'e'
                        submit.call_subsystem(self,cmdstr,wait=True)    
                        cmdstr=self.trueanalyzepath+' '+xyz+' '+'-k'+' '+self.originalkeyfilename+' '+'m'+'> alz.out'
                        submit.call_subsystem(self,cmdstr,wait=True)    
                        chg,resid=self.ReadCharge('alz.out')
                        self.otherligcharge+=chg
                for xyz in self.annihilateligandxyzfilenamelist:
                    cmdstr=self.trueanalyzepath+' '+xyz+' '+'-k'+' '+self.originalkeyfilename+' '+'e'
                    submit.call_subsystem(self,cmdstr,wait=True)    
                    cmdstr=self.trueanalyzepath+' '+xyz+' '+'-k'+' '+self.originalkeyfilename+' '+'m'+'> alz.out'
                    submit.call_subsystem(self,cmdstr,wait=True)    
                    chg,resid=self.ReadCharge('alz.out')
                    self.ligandcharge+=chg

            if ligandonly==True:
                return []
            mpolearrays=[]

            if self.receptorligandxyzfilename!=None and self.originalkeyfilename!=None: 
                
                cmdstr=self.trueanalyzepath+' '+self.receptorligandxyzfilename+' '+'-k'+' '+self.originalkeyfilename+' '+'e'
                submit.call_subsystem(self,cmdstr,wait=True)    
                cmdstr=self.trueanalyzepath+' '+self.receptorligandxyzfilename+' '+'-k'+' '+self.originalkeyfilename+' '+'m'+'> alz.out'
                submit.call_subsystem(self,cmdstr,wait=True)    
                self.complexcharge,resid=self.ReadCharge('alz.out',checkresid=False)
                if resid!=0:
                    index,typenum=self.GrabFirstReceptorIndexAndType(self.receptorligandxyzfilename)
                    mpolearrays=self.GrabMultipoleParameters(self.prmfilepath,typenum)
                    mpolearrays=self.ModifyIndexAndChargeMultipole(mpolearrays,resid,index)
                    self.ModifyCharge(self.originalkeyfilename,mpolearrays)
                    cmdstr=self.trueanalyzepath+' '+self.receptorligandxyzfilename+' '+'-k'+' '+self.originalkeyfilename+' '+'m'+'> alz.out'
                    submit.call_subsystem(self,cmdstr,wait=True)    
                    self.complexcharge,resid=self.ReadCharge('alz.out')
                else:
                    mpolearrays=[]

            

            if self.receptorligandxyzfilename!=None:
                self.receptorcharge=self.complexcharge-self.ligandcharge-self.otherligcharge # this includes the charges from ions in pocket initially
            if self.receptorligandxyzfilename!=None and self.ligandxyzfilenamelist!=None:
                for xyz in self.ligandxyzfilenamelist:
                    ligandtypes=self.GrabTypeNumbers(xyz)
                    indices=self.GrabIndicesWithTypeNumber(self.receptorligandxyzfilename,ligandtypes)
                    coords=self.GrabCoordinates(self.receptorligandxyzfilename,indices)
                    self.RewriteCoordinates(xyz,coords)
            return mpolearrays


        def GrabFirstReceptorIndexAndType(self,xyzfilename):
            """
            Intent: When want to deposit residual charge on protein to ensure net charge is an integer. Just pick an atom (hydrogen) on protein and then grab that multipole later on for adding residual charge. This funciton only grabs index and type number.
            Input: Tinker XYZ 
            Output: index and type number of atom to deposit residual charge on
            Referenced By: CheckInputXYZKeyFiles
            Description: 
            1. Iterate over results of tinker XYZ file.
            2. If encounter a hydrogen atom return index and type number 
            """
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            # STEP 1
            for line in results:
                linesplit=line.split()
                if len(linesplit)>1 and '90.0' not in line:
                    index=int(linesplit[0])
                    element=linesplit[1]
                    typenum=int(linesplit[5])
                    # STEP 2
                    if element=='H':
                        return index,typenum



        def GrabMultipoleParameters(self,prmfilepath,typenum):
            """
            Intent: Grab the multipole parameters that will be modified to deposit residual charge upon.
            Input: Tinker parameter file, type number specifying which multipole parameters to grab.
            Output: Array of multipole parameters
            Referenced By: CheckInputXYZKeyFiles
            Description: 
            1. Iterate over parameter file
            2. If multipole keyword is detected and the desired type number corresponds to that multipole parameter
            3. Grab all the associated lines then append to an array
            """
            mpolearrays=[]
            temp=open(prmfilepath,'r')
            results=temp.readlines()
            temp.close()
            # STEP 1
            for lineidx in range(len(results)):
                line=results[lineidx]
                # STEP 2
                if 'multipole' in line:
                    linesplit=line.split()
                    if linesplit[1]==str(typenum):
                        # STEP 3
                        mpolearrays.append(line)
                        dip=results[lineidx+1]
                        firstquad=results[lineidx+2]
                        secondquad=results[lineidx+3]
                        thirdquad=results[lineidx+4]
                        mpolearrays.append(dip)
                        mpolearrays.append(firstquad)
                        mpolearrays.append(secondquad)
                        mpolearrays.append(thirdquad)

            return mpolearrays


        def ModifyIndexAndChargeMultipole(self,mpolearrays,resid,index):
            """
            Intent: Modify multipole parameters to ensure there is no net residual charge (net charge is integer).
            Input: Array of multipole parameters, residual charge to add, atom index of atom that is being modified (tinker allows putting two copies of mlutipole parameters but if you specify atom index instead of type only affects that atom).
            Output: Modifed array of multipole parameters.
            Referenced By: CheckInputXYZKeyFiles
            Description:
            1. Iterate over array of multipole parameters
            2. If line contains multipole keyword (also has charge in that line)
            3. Then modify charge by adding residual and add atom index instead of type number
            4. Append back to array
            """
            # STEP 1
            for idx in range(len(mpolearrays)):
                line=mpolearrays[idx]
                # STEP 2
                if 'multipole' in line:
                    chargelinesplit=line.split()
                    # STEP 3
                    chargelinesplit[1]='-'+str(index)
                    charge=float(chargelinesplit[-1])
                    charge=str(charge+resid)
                    chargelinesplit[-1]=charge
                    chargeline=' '.join(chargelinesplit)+'\n'
                    # STEP 4
                    mpolearrays[idx]=chargeline 

            return mpolearrays


        def ModifyCharge(self,keyfilename,mpolearrays):
            """
            Intent: Append array of modified multipole parameters to keyfile.
            Input: Tinker key file, array of modified multipole parameters
            Output: -
            Referenced By: CheckInputXYZKeyFiles , MolecularDynamics
            Description:
            1. Iterate over lines in array of multipole parameters 
            2. Append line to key file
            """
            temp=open(keyfilename,'a')
            # STEP 1
            for line in mpolearrays:
                # STEP 2
                temp.write(line)
            temp.close()


        def DetermineIonsForChargedSolvation(self):
            """
            Intent: If user specifies to compute the "salt" hydration free energy, charged HFE + HFE of ions with that net charge, then this function computes the ions that need to be added for the HFE of ions box. 
            Input: - 
            Output: Dictionary of ion type number to number of ions needed to be added to box  
            Referenced By: MolecularDynamics 
            Description:
            1. If user wants to compute the salt hydration free energy
            2. If there is a positive net charge, grab Cl type number, negative net charge, grab K type number.
            3. Store type number in dictionary along with magnitude of charge (represents number of ions needed to be added).  
            """
            solviontocount={}
            self.addsolvionwindows=False
            # STEP 1
            if self.binding==False and self.solvation==True and self.salthfe==True:
                for i in range(len(self.systemcharge)):
                    chg=int(self.systemcharge[i])
                    # STEP 2
                    if chg>0:
                        iontypenumber=self.elementsymtotinktype['Cl']
                    elif chg<0:
                        iontypenumber=self.elementsymtotinktype['K']
                    else:
                        continue
                    count=np.abs(chg)
                    # STEP 3
                    solviontocount[iontypenumber]=count
                    self.addsolvionwindows=True
            return solviontocount    

        def ChangeTypeNumbers(self,xyzfile,elementtotype):
            """
            Intent: Preequilbriated box may be using water type numbers from another parameter file, so need to change type numbers to ensure its consistent with input parameter file.
            Input: Tinker XYZ file, dictionary of element to tinker type number
            Output: Updated XYZ filename
            Referenced By: TrimPreEquilibriatedBox
            Description: 
            1. Iterate over XYZ file
            2. If line is not a line containing boxinformation and its not the first line (atom number)
            3. Extract the element from the line
            4. Determine type number from the input dictionary and extracted element
            5. Insert new type number into new line and write to file
            """
            temp=open(xyzfile,'r')
            results=temp.readlines()
            temp.close()
            tempname=xyzfile.replace('.xyz','_new.xyz')
            temp=open(tempname,'w')
            # STEP 1
            for lineidx in range(len(results)):
                line=results[lineidx]
                linesplit=line.split()
                isboxline=True
                # STEP 2
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if isboxline==False and lineidx!=0:
                    # STEP 3
                    element=linesplit[1] 
                    # STEP 4
                    typenum=elementtotype[element]
                    linesplit[5]=str(typenum)
                    # STEP 5
                    line=' '.join(linesplit)+'\n'
                temp.write(line)
            temp.close()
            return tempname
        
        def ModifyBoxSizeInXYZFile(self,xyzfile,boxsize):
            """
            Intent: Change the box size dimensions, needed for when trimming the pre-equilibriated box.
            Input: Tinker XYZ box file, new box size
            Output: -
            Referenced By: TrimPreEquilibriatedBox
            Description: 
            1. Iterate over box file
            2. If line contains box information
            3. Modify the line to contain new box size information
            4. Write new line to file
            """
            temp=open(xyzfile,'r')
            results=temp.readlines()
            temp.close()
            tempname=xyzfile.replace('.xyz','_new.xyz')
            temp=open(tempname,'w')
            # STEP 1
            for lineidx in range(len(results)):
                line=results[lineidx]
                linesplit=line.split()
                # STEP 2
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if isboxline==True:
                    # STEP 3
                    linesplit[0]=str(boxsize[0])
                    linesplit[1]=str(boxsize[1])
                    linesplit[2]=str(boxsize[2])
                    line=' '.join(linesplit)+'\n'
                # STEP 4
                temp.write(line)
            temp.close()
            os.remove(xyzfile)
            os.rename(tempname,xyzfile)

        def TrimPreEquilibriatedBox(self,boxsize):
            """
            Intent: Instead of minimizing new box and equilbriating for long period, it is sometimes easier to trim a larger pre-equibriated box to new box sixe and equilbriate for a shorter time.
            Input: New desired box size
            Output: - 
            Referenced By: BoxSetupProtocol in boxsetup.py
            Description:
            1. Copy pre-equilbriated box from poltype source to current directory 
            2. Change type numbers of water to refelect current parameter file
            3. Change the size of the boxfile to desired box size
            4. Increase box size slightly to accomdate water that will be removed outside of input box size (some may be across original boundary) 
            5. Remove waters outside of original box size
            6. Modify box file to new box size again
            """
            xyzfile=self.preequilboxpath
            split=os.path.split(xyzfile)
            xyzfilename=split[1]
            # STEP 1
            shutil.copy(xyzfile,os.path.join(self.simpath,xyzfilename))    
            elementtotype={'O':self.waterOtypenum,'H':self.waterHtypenum}
            # STEP 2
            filename=self.ChangeTypeNumbers(xyzfilename,elementtotype)
            # STEP 3
            self.ModifyBoxSizeInXYZFile(filename,boxsize)
            # STEP 4
            newboxsize=[2+i for i in boxsize]
            # STEP 5
            self.RemoveWaterOutsideBox(filename,boxsize)
            # STEP 6
            self.ModifyBoxSizeInXYZFile(filename,newboxsize)
            os.rename(filename,'water.xyz_2') # for xyzedit lib

        def RemoveWaterOutsideBox(self,xyzfilename,boxsize):
            """
            Intent: Needed when trimming pre-equilbriated box. Remove waters outside of desired box size. 
            Input: XYZ box file, desired box size.
            Output: -
            Referenced By: TrimPreEquilibriatedBox
            Description:
            1. Create input file to call xyzedit
            2. Write the appropriate options for removing water outside box size
            3. Generate command string to call xyzedit with filename generated as input
            4. Submit command to system
            5. Rename the new filename to original input filename
            """
            # STEP 1
            tempfile=open('xyzedit.in','w')
            # STEP 2
            tempfile.write('17'+'\n')
            tempfile.write(str(boxsize[0])+','+str(boxsize[1])+','+str(boxsize[2])+'\n')
            tempfile.write('\n')
            tempfile.close()
            # STEP 3
            cmdstr=self.xyzeditpath +' '+xyzfilename+' '+self.prmfilepath+' '+' < xyzedit.in '
            # STEP 4
            submit.call_subsystem(self,cmdstr,wait=True)   
            # STEP 5
            newname=xyzfilename+'_2'
            os.rename(newname,xyzfilename)


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
       

        def AlignComplexedPDBToUnComplexedPDB(self):
            """
            Intent: After adding missing residues and assigning protonation state, may need to add ligands back in original PDB. So this requires an allignment procedure.
            Input: -
            Output: -
            Referenced By: MolecularDynamics 
            Description:
            1. Initialize a reference and current mol object
            2. Read in complexed protein PDB to current mol and uncomplexed PDB to reference mol object
            3. Call openbabel aligment object and set the reference and target mol objects
            4. Align the mol objects
            5. Update coordiantes in target mol file
            6. Grab HETATM indices from mol object
            7. Add those indices along with coordinates from mol object (that was aligned to ref) to the uncomplexed PDB
            8. Write out new PDB with uncomplex structure that has aligned ligand from complexed structure.
            """
            # STEP 1
            obConversion = openbabel.OBConversion()
            ref = openbabel.OBMol()
            mol = openbabel.OBMol()
            # STEP 2
            inFormat = obConversion.FormatFromExt(self.uncomplexedproteinpdbname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, self.complexedproteinpdbname)
            obConversion.ReadFile(ref, self.uncomplexedproteinpdbname)
            # STEP 3
            aligner = openbabel.OBAlign(False, False)
            aligner.SetRefMol(ref)
            aligner.SetTargetMol(mol)
            # STEP 4
            aligner.Align()
            # STEP 5
            aligner.UpdateCoords(mol)
            # STEP 6
            hetatmids=self.GrabHETATMS(mol)
            # STEP 7
            ref=self.AddHETATMSToPDB(ref,mol,hetatmids)
            # STEP 8
            obConversion.SetOutFormat('pdb')
            newmolfile=self.uncomplexedproteinpdbname.replace('.pdb','_aligned.pdb')
            obConversion.WriteFile(ref,newmolfile)


        def AddHETATMSToPDB(self,ref,mol,hetatmids):
            """
            Intent: Grab atoms from one mol and add to another (when adding aligned ligand in complexed PDB to uncomplexed PDB)
            Input: reference and current mol objects, array of ligand atom indices
            Output: updated reference mol object
            Referenced By: AlignComplexedPDBToUnComplexedPDB
            Description: 
            1. Iterate over array of ligand atom indices 
            2. Grab atom from current mol object with given atom index in array
            3. Add that atom to reference mol object
            """
            for atomindex in hetatmids:
                atom=mol.GetAtom(atomindex)
                ref.AddAtom(atom)


            return ref


        def GrabHETATMS(self,mol):
            """
            Intent: Determine ligand indices from PDB mol object, needed for when adding aligned ligand to complexed PDB to uncomplexed PDB
            Input: mol object
            Output: Array of PDB ligand indices
            Referenced By: AlignComplexedPDBToUnComplexedPDB
            Description:
            1. Iterate over atoms of mol object
            2. Grab residue information
            3. If residue is not in canoncical residue symbols, then it must be ligand
            4. Append atom index to array
            """
            hetatmids=[]
            iteratom = openbabel.OBMolAtomIter(mol)
            # STEP 1
            for atm in iteratom:
                atmindex=atm.GetIdx()
                # STEP 2
                res=atm.GetResidue()
                reskey=res.GetName()
                # STEP 3
                if reskey not in self.knownresiduesymbs:
                    # STEP 4
                    hetatmids.append(atmindex)


            return hetatmids

            
        def AlignLigandXYZToTemplateXYZ(self):
            """
            Intent: User may wish to align a ligand XYZ to a template ligand XYZ in protein pocket for example.
            Input: -
            Output: -
            Referenced By: MolecularDynamics 
            Description: 
            1. Convert input ligand XYZ and template ligand XYZ into cartesian XYZ files ( to be processed by chemoinformatic tookkits)
            2. Generate reference and current mol objects
            3. Read in cartesian XYZs into mol objects
            4. Convert cartesian XYZs into MOL format with babel (to be read in by rdkit, rdkit doesnt read cartesian) 
            5. Read in MOL files into rdkit mol objects
            6. Remove bond topology information to make matching easier (sometimes cartesian->MOL creates wrong bond information that doesnt match so remove)
            7. Sanitize the rdkit mol objects (assign formal charges otherwise rdkit complain)
            8. Embed many conformations into target mol object (this will allow finding best conformation that matches to other mol)
            9. Compute RMSD between template XYZ and each conformor of input XYZ, then choose conformer with minimum RMSD 
            10.Update XYZ coordinates in tinker files with conformer coordinates
            11.Convert tinker XYZ coordaintes to cartesian XYZ  
            """
            # STEP 1
            ligandcartxyz=self.ConvertTinktoXYZ(self.ligandxyzfilenamelist[0],self.ligandxyzfilenamelist[0].replace('.xyz','_cart.xyz'))
            templateligandcartxyz=self.ConvertTinktoXYZ(self.templateligandxyzfilename,self.templateligandxyzfilename.replace('.xyz','_cart.xyz'))
            # STEP 2
            obConversion = openbabel.OBConversion()
            ref = openbabel.OBMol()
            mol = openbabel.OBMol()
            # STEP 3
            inFormat = obConversion.FormatFromExt(ligandcartxyz)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, ligandcartxyz)
            obConversion.ReadFile(ref, templateligandcartxyz)
            # STEP 4
            obConversion.SetOutFormat('mol')
            molfile=ligandcartxyz.replace('_cart.xyz','.mol')
            reffile=templateligandcartxyz.replace('_cart.xyz','.mol')
            obConversion.WriteFile(mol,molfile)
            obConversion.WriteFile(ref,reffile)
            # STEP 5
            mol=Chem.MolFromMolFile(molfile,removeHs=False,sanitize=False)
            ref=Chem.MolFromMolFile(reffile,removeHs=False,sanitize=False)
            # STEP 6
            for bond in mol.GetBonds():
                if bond.IsInRing()==False:
                    bond.SetBondType(Chem.BondType.SINGLE)
            for bond in ref.GetBonds():
                if bond.IsInRing()==False:
                    bond.SetBondType(Chem.BondType.SINGLE)
            # STEP 7
            self.ligandcharge=None
            mol,atomindextoformalcharge=self.CheckInputCharge(mol)
            self.ligandcharge=None
            ref,atomindextoformalcharge=self.CheckInputCharge(ref)
            self.ligandcharge=None
            Chem.SanitizeMol(mol)
            Chem.SanitizeMol(ref)
            # STEP 8
            confnum=1000
            AllChem.EmbedMultipleConfs(mol, confnum)
            mcs = rdFMCS.FindMCS([ref,mol])
            patt = Chem.MolFromSmarts(mcs.smartsString)
            atomnum=mcs.numAtoms
            refMatch = ref.GetSubstructMatch(patt)
            mv = mol.GetSubstructMatch(patt)
            # STEP 9
            rmstocid={}
            for cid in range(confnum):
                rms = AllChem.AlignMol(mol,ref,prbCid=cid,atomMap=list(zip(mv,refMatch)))
                rmstocid[rms]=cid
            minrms=min(rmstocid.keys())
            mincid=rmstocid[minrms]
            # STEP 10
            indextocoords={}
            for atom in mol.GetAtoms():
                atomidx=atom.GetIdx()
                pos = mol.GetConformer(mincid).GetAtomPosition(atomidx)
                X=pos.x
                Y=pos.y
                Z=pos.z
                atomindex=atomidx+1 
                indextocoords[atomindex]=[X,Y,Z]
            name=self.ligandxyzfilenamelist[0].replace('.xyz','_aligned.xyz')
            self.ReplaceXYZCoords(self.ligandxyzfilenamelist[0],indextocoords,name,replace=False)
            # STEP 11
            alignedligandcartxyz=self.ConvertTinktoXYZ(name,name.replace('.xyz','_cart.xyz'))
 
        
        def CallAnalyze(self,statexyz,statekey,alzout,analyzepath,option):
            """
            Intent: Call analyze for many purposes, checking parameters are defined, checking energy components or total charge.
            Input: XYZ file, key file, analyze output file, analyze bin path, analyze option
            Output: - 
            Referenced By: Many functions 
            Description: 
            1. Construct command string for analyze
            2. Submit command to call_subsystem  
            """
            # STEP 1
            cmdstr=analyzepath+' '+statexyz+' '+'-k'+' '+statekey+' '+option+' > '+alzout
            # STEP 2
            submit.call_subsystem(self,cmdstr,True,False,alzout) 

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
            atomindices=databaseparser.RingAtomicIndices(self,mol)
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


        def ExtractLigand(self,ligandreceptorfilename,coordinates=None,indicestokeep=[]):
            """
            Intent: Used for docking module and pdbxyz module to extract the ligand from PDB into seperate files.
            Input: ligand-receptor PDB file name, optional coordinates
            Output: ligand only PDB, receptor only PDB
            Referenced By:  docking.py - DockingWrapper , pdbxyz.py - DeterminePocketGrid
            Description:
            1. Call ExtractMOLObject for receptor, using known residue symbols only
            2. Call ExtractMOLObject for ligand, grabbing any atoms with unknown residues
            3. Convert residue name UNL to LIG, for other parts of program to read ligand easier in PDB
            """
            ligandpdbfilename='ligand.pdb'
            receptorpdbfilename=ligandreceptorfilename.replace('.pdb','_receptoronly.pdb')
            # STEP 1
            receptormol=self.ExtractMOLObject(ligandreceptorfilename,receptorpdbfilename,None,True,indicestokeep)
            # STEP 2
            ligandmol=self.ExtractMOLObject(ligandreceptorfilename,ligandpdbfilename,coordinates,False,indicestokeep)
            # STEP 3
            self.ConvertUNLToLIG(ligandpdbfilename)
            return ligandpdbfilename,receptorpdbfilename


        def ExtractMOLObject(self,ligandreceptorfilename,newpdbfilename,coordinates,receptor,indicestokeep):
            """
            Intent: Extract either ligand or receptor from PDB file, then create new PDB file.
            Input: ligand-receptor PDB filename, new PDB filename, array of coordinates to only include certain atoms, boolean if grabbing receptor or ligand, array of indices to keep 
            Output: Extracted mol object
            Referenced By: ExtractLigand
            Description:
            1. Define common HETATM symbols such as water
            2. Generate mol object from ligand-receptor PDB
            3. Iterate over atoms
            4. Grab residue information from atom
            5. If not in known residue symbols and trying to extract receptor, then add atom to array of atoms to remove
            6. If extracting ligand, if coordinates array is given, delete any atoms not in coordinates array. Otherwise if residue is a known residue symbol or common HETATM symbol, then append atoms to array to delete. 
            7. Remove atoms that should be deleted
            8. Write out new PDB file
            """
            # STEP 1
            commonhetatmsymbs=['HOH']
            # STEP 2
            pdbmol=self.GenerateMOLObject(ligandreceptorfilename)
            iteratom = openbabel.OBMolAtomIter(pdbmol)
            obConversion = openbabel.OBConversion()
            obConversion.SetOutFormat('pdb')
            atmindicestodelete=[]
            # STEP 3
            for atm in iteratom:
                atmindex=atm.GetIdx()
                # STEP 4
                res=atm.GetResidue()
                reskey=res.GetName()
                coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
                # STEP 5
                if reskey not in self.knownresiduesymbs and receptor==True:
                    atmindicestodelete.append(atmindex)
                elif receptor==False:
                    # STEP 6
                    if coordinates!=None:
                        if coords not in coordinates:
                            atmindicestodelete.append(atmindex)
                    else:
                        if (reskey in self.knownresiduesymbs or reskey in commonhetatmsymbs) and atmindex not in indicestokeep:
                            atmindicestodelete.append(atmindex)
            # STEP 7
            atmindicestodelete.sort(reverse=True)
            for atmindex in atmindicestodelete:
                atm=pdbmol.GetAtom(atmindex)
                pdbmol.DeleteAtom(atm)
            # STEP 8
            obConversion.WriteFile(pdbmol,newpdbfilename)
            return pdbmol

        def GenerateMOLObject(self,pdbfilename):
            """
            Intent: Convert PDB file into mol object
            Input: PDB filename
            Output: PDB mol object
            Referenced By: GrabLigandCentroid,ExtractMOLObject
            Description:
            1. Call conversion object from openbabel
            2. Create empty mol object
            3. Set input format for converter
            4. Read in PDB file into mol file via converter
            """
            # STEP 1
            obConversion = openbabel.OBConversion()
            # STEP 2
            pdbmol = openbabel.OBMol()
            # STEP 3
            obConversion.SetInFormat('pdb')
            # STEP 4
            obConversion.ReadFile(pdbmol,pdbfilename)
            return pdbmol


        def ConvertUNLToLIG(self,filename):
            """
            Intent: Make it more clear in PDB file what is ligand by changing residue name to LIG
            Input: PDB filename
            Output: Modified PDB filename
            Referenced By: ExtractLigand
            Description: 
            1. Iterate over results of PDB filename
            2. If UNL detected in spot for residue name or UNK detected, then replace with LIG
            """
            temp=open(filename,'r')
            results=temp.readlines()
            temp.close()
            tempname=filename.replace('.pdb','_TEMP.pdb')
            temp=open(tempname,'w')
            # STEP 1
            for line in results:
                linesplit=re.split(r'(\s+)', line)
                # STEP 2
                if 'UNL' in line:
                    lineindex=linesplit.index('UNL')
                    linesplit[lineindex]='LIG'
                    line=''.join(linesplit)
                if 'UNK' in line:
                    lineindex=linesplit.index('UNK')
                    linesplit[lineindex]='LIG'
                    line=''.join(linesplit)
        
                temp.write(line)
            temp.close()
            os.rename(tempname,filename)


        def GrabLigandCentroid(self,ligandpdbfilename):
            """
            Intent: Grab centroid of ligand, used for defining a pocket grid to detect any ions/waters near ligand that should be kept when running simulations.  
            Input: Ligand PDB filename
            Output: Coordinates of ligand centroid
            Referenced By: pdbxyz.py - DeterminePocketGrid , docking.py - DockingWrapper 
            Description:
            1. Call GenerateMOLObject to grab mol from ligand PDB
            2. Call GrabAtomPositions to determine atomic coordinates
            3. Compute centroid of atomic positions
            """
            # STEP 1
            pdbmol=self.GenerateMOLObject(ligandpdbfilename)
            # STEP 2
            atomvecls=self.GrabAtomPositions(pdbmol)
            # STEP 3
            centroid=np.array([0.0,0.0,0.0])
            for vec in atomvecls:
                centroid+=np.array(vec)
            if len(atomvecls)==0:
                raise ValueError('Ligand in PDB file is not labeled as LIG')
            centroid=centroid/len(atomvecls)
        
            return centroid


        def GrabAtomPositions(self,pdbmol):
            """
            Intent: Grab atomic coordinates from mol object
            Input: mol object
            Output: Array of atomic coordinates
            Referenced By: GrabLigandCentroid
            Description:
            1. Iterate over atoms in mol object
            2. Grab atom coordinates
            3. Append to array of atomic coordinates
            """
            atomvecls=[]
            iteratombab = openbabel.OBMolAtomIter(pdbmol)
            # STEP 1
            for atm in iteratombab:
                atmindex=atm.GetIdx()
                # STEP 2
                coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
                # STEP 3
                atomvecls.append(coords)
        
         
            return atomvecls


        

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


        def ExtractLigandIndicesFromComplexXYZ(self,receptorligandxyzfilename,oldtypetonewtypelist,ligatomnums):
            """
            Intent: Need indices of ligands that will be annihilated and indices of ligands that will not be annihilated. Need to keep track of indices for restraints of knowing which atoms to define with "ligand" keyword. Also used for shifting type numbers if any types are overlapping in input ligand keys (generating index->new type dictionary)  
            Input: receptor-ligand tinker XYZ, array of old type to new type dictionaries.
            Output: dictionary of atom index to new type number, ligand indices in receptor-ligand complex, ligand indices not in complex 
            Referenced By: MolecularDynamics 
            Description: 
            1. Extract old index to type index from input tinker XYZ
            2. Define list of index to old types to keep track of previous indices used to assist with only adding correct ligand indices for each ligand into seperate arrays (reading all ligand indices from one complexed XYZ file) 
            3. Iterate over array of dictionaries old type -> new type (one for each ligand)
            4. Grab current old type -> new type dictionary
            5. Grab atom number of current ligand
            6. Define empty dictionary indextooldtypenum to keep track of number of atoms adding (stop when reach atom number of current ligand)
            7. Get array of previous atom indices used to assist in making sure only correct atom indices are added for each seperate ligand in their own arrays.
            8. Iterate over dictionary of old index -> type index (from tinker XYZ)
            9. If type index in current ligand dictionary of old type -> new type and current length of indextooldtypenum hasnt reached current atom number for current ligand and indices havent been previously used then add information to indextooldtypenum (contains information of ligand indices seperated into individual dictionaries for each ligand) and indextonewtype for shifting paramteer types later.
            10. Extract ligand indices from indextooldtypenum into seperate arrays and save in complexligands array
            11. Generate the equivalent ligand indices from in complex but for solvation phase (assume indices are on top of box file and solvent after) 
            """
            # STEP 1
            statexyzatominfo,oldindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=self.GrabXYZInfo(receptorligandxyzfilename)
            # STEP 2
            listofindextooldtypes=[]
            indextonewtype={}
            complexligands=[]
            # STEP 3
            for oldtypetonewtypeidx in range(len(oldtypetonewtypelist)):
                # STEP 4
                oldtypetonewtype=oldtypetonewtypelist[oldtypetonewtypeidx]
                # STEP 5
                atomnum=ligatomnums[oldtypetonewtypeidx]
                # STEP 6
                indextooldtypenum={}
                # STEP 7
                if oldtypetonewtypeidx!=0:
                    previousindextooldtypes=listofindextooldtypes[:oldtypetonewtypeidx-1+1]
                else:
                    previousindextooldtypes=[]
                previousindices=[]
                for oindextooldtypenum in previousindextooldtypes:
                    previousindices.extend(oindextooldtypenum.keys()) 
                # STEP 8
                for oldindex,typeindex in oldindextotypeindex.items():
                    # STEP 9
                    if typeindex in oldtypetonewtype.keys() and len(indextooldtypenum.keys())<atomnum and oldindex not in previousindices:
                        indextooldtypenum[oldindex]=typeindex
                        newtype=oldtypetonewtype[typeindex]
                        indextonewtype[oldindex]=newtype
                # STEP 10
                listofindextooldtypes.append(indextooldtypenum)
                complexligands.append(list(indextooldtypenum.keys()))
            # STEP 11
            solvligands=[]
            count=1
            for ligand in complexligands:
                lig=[]
                for i in range(len(ligand)):
                    lig.append(count)
                    count+=1
                solvligands.append(lig)
            

            return indextonewtype,complexligands,solvligands


        def GrabLigandTypesInfoForIndicesExtraction(self,ligandxyzfilenamelist):
            """
            Intent: Grab type numbers and atom numbers for each ligand in input tinker XYZ list. This will be used for generating list of complexed ligand indices and solvation ligand indices. Need to keep track of indices for restraints of knowing which atoms to define with "ligand" keyword. 
            Input: List of tinker XYZs for each ligand
            Output: List of type numbers, list of atom numbers
            Referenced By: MolecularDynamics 
            Description: 
            1. Iterate over list of tinker XYZ
            2. Call GrabTypeNumbers to grab type numbers
            3. Append information about type numbers to array
            4. Call GrabXYZInfo to extract total atom number
            5. Save atom number information to array 
            """
            oldtypelist=[]
            ligatomnums=[]
            # STEP 1
            for xyz in ligandxyzfilenamelist:
                # STEP 2
                ligandtypes=self.GrabTypeNumbers(xyz) 
                oldtypetooldtype=dict(zip(ligandtypes,ligandtypes))
                # STEP 3
                oldtypelist.append(oldtypetooldtype)
                # STEP 4
                statexyzatominfo,oldindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=self.GrabXYZInfo(xyz)
                # STEP 5
                ligatomnums.append(stateatomnum)


            return oldtypelist,ligatomnums


        def ShiftParameterTypesComplexXYZ(self,receptorligandxyzfilename,indextonewtype):
            """
            Intent: If duplicate ligands have overlapping types, then shift them to be different (in case users wants to only annihilate one and not a duplicate). So need to shift types in tinker XYZ also. 
            Input: Receptor-ligand XYZ file, dictionary of index to new type (replaces old type) 
            Output: - 
            Referenced By: MolecularDynamics
            Description: 
            1. Iterate over lines of tinker XYZ
            2. Iterate over strings in each line
            3. If string is integer, grab atom index
            4. Grab new type number from dicionary indextonewtype
            5. Replace new type in line
            """
            tempname=receptorligandxyzfilename+"_TEMP"
            temp=open(receptorligandxyzfilename,'r')
            results=temp.readlines()
            temp.close()
            temp=open(tempname,'w')
            # STEP 1
            for line in results:
                linesplitall=re.split(r'(\s+)', line)
                linesplit=line.split()
                # STEP 2
                for i in range(len(linesplitall)):
                    element=linesplitall[i]
                    if '.xyz' in receptorligandxyzfilename and (i!=12):
                        continue
                    if self.RepresentsInt(element):
                        oldtypenum=np.abs(int(element))
                        # STEP 3
                        index=int(linesplit[0])
                        if index in indextonewtype.keys():
                            # STEP 4
                            typenum=indextonewtype[index]
                            # STEP 5
                            linesplitall[i]=str(typenum)
                line=''.join(linesplitall)
                temp.write(line)
            temp.close()
            os.rename(tempname,receptorligandxyzfilename)



        def ShiftParameterTypes(self,filename,oldtypetonewtype):
            """
            Intent: Shift key/XYZ types if overlapping (only for ligand XYZ/keys not complex)
            Input: Tinker key or XYZ, dictionary of old type to new type
            Output: - 
            Referenced By: MolecularDynamics 
            Description: 
            1. Iterate over lines of tinker XYZ
            2. Iterate over strings in each line
            3. If string is integer, grab type from line
            4. Grab new type number from dicionary oldtypetonewtype
            5. If - in front of type (multipoles), then add back in front of new type
            6. Update line with new type number 
            """
            tempname=filename+"_TEMP"
            temp=open(filename,'r')
            results=temp.readlines()
            temp.close()
            temp=open(tempname,'w')
            # STEP 1
            for line in results:
                linesplitall=re.split(r'(\s+)', line)
                # STEP 2
                for i in range(len(linesplitall)):
                    element=linesplitall[i]
                    if '.xyz' in filename and (i!=12):
                        continue
                    if self.RepresentsInt(element):
                        # STEP 3
                        oldtypenum=np.abs(int(element))
                        if oldtypenum in oldtypetonewtype.keys():
                            # STEP 4
                            newtypenum=oldtypetonewtype[oldtypenum]
                            typenum=newtypenum
                        else:
                            typenum=oldtypenum
                        # STEP 5
                        if '-' in element:
                            typenum=-typenum
                        # STEP 6
                        linesplitall[i]=str(typenum)
                line=''.join(linesplitall)
                temp.write(line)
            temp.close()
            os.rename(tempname,filename)


        def GrabMaxTypeNumber(self,parameterfile):
            """
            Intent: Need max type number for shifting type numbers (dont want overlapping types for same duplicate ligands if want to disappear one and not the duplicate etc). 
            Input: parameter or key file
            Output: max type number
            Referenced By: CheckIfNeedToShiftTypes, GenerateTypeMaps
            Description: 
            1. Assume max type number is 1
            2. Iterate over lines of file
            3. If atom keyword in line and its not a comment, then grab atom type number
            4. If atom type is greater than current max, update the max type 
            """
            # STEP 1
            maxnumberfromprm=1
            temp=open(parameterfile,'r')
            results=temp.readlines()
            temp.close()
            # STEP 2
            for line in results:
                # STEP 3
                if 'atom' in line and '#' not in line:
                    linesplit=line.split()
                    atomtype=int(linesplit[1])
                    # STEP 4
                    if atomtype>maxnumberfromprm:
                        maxnumberfromprm=atomtype
            return maxnumberfromprm
        
        def GrabMinTypeNumber(self,parameterfile):
            """
            Intent: Need min type number for shifting type numbers (dont want overlapping types for same duplicate ligands if want to disappear one and not the duplicate etc). 
            Input: parameter or key file
            Output: min type number
            Referenced By: CheckIfNeedToShiftTypes, GenerateTypeMaps
            Description: 
            1. Assume min type number is 10000
            2. Iterate over lines of file
            3. If atom keyword in line and its not a comment, then grab atom type number
            4. If atom type is smaller than current min, update the min type 
            """
            minnumberfromprm=10000
            temp=open(parameterfile,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                if 'atom' in line and '#' not in line:
                    linesplit=line.split()
                    atomtype=int(linesplit[1])
                    if atomtype<minnumberfromprm:
                        minnumberfromprm=atomtype
            return minnumberfromprm

        def GenerateTypeMaps(self,keyfilelist):
            """
            Intent: This creates old type -> new type dictionaries for each ligand, for shifting parameter types. 
            Input: List of ligand key files
            Output: List of old type -> new type dictionaries for each ligand. 
            Referenced By: MolecularDynamics
            Description: 
            1. Iterate over array of key files
            2. Grab min and max type number from key file
            3. Define shift as previous max type - previous min type + 1
            4. If min from key is smaller than current min, then update current min. 
            5. If max from key is greater than current max, then update current max. 
            6. Generate all types from min type from key and max type from key
            7. Shift all the types by the shift value defined earlier
            8. Create dictionary of old types -> shifted types
            9. Append dictionary to list of dictionaries
            10.Update previousmaxnumberfromkey and previousminnumberfromkey 
            """
            newkeyfilelist=keyfilelist.copy()
            oldtypetonewtypelist=[]
            currentmin=100000
            currentmax=0
            shift=0
            prevmaxnumberfromkey=0
            prevminnumberfromkey=0
            #originalshift=self.GrabMaxTypeNumber(self.prmfilepath)
            originalshift=0 # cant shift because max 1000 for atom classes in tinker
            # STEP 1
            for keyfilename in newkeyfilelist:
                # STEP 2
                maxnumberfromkey=self.GrabMaxTypeNumber(keyfilename)
                minnumberfromkey=self.GrabMinTypeNumber(keyfilename)
                # STEP 3
                shift+=prevmaxnumberfromkey-prevminnumberfromkey+1+originalshift
                # STEP 4
                if minnumberfromkey<currentmin:
                    currentmin=minnumberfromkey
                # STEP 5
                if maxnumberfromkey>currentmax:
                    currentmax=maxnumberfromkey
                # STEP 6
                types=np.arange(minnumberfromkey,maxnumberfromkey+1,1)
                # STEP 7
                shiftedtypes=types+shift
                maxtype=max(shiftedtypes)
                currentmax=maxtype
                # STEP 8
                oldtypetonewtype=dict(zip(types,shiftedtypes))
                temp={}
                for oldtype,newtype in oldtypetonewtype.items():
                    negold=-oldtype # for when multipoles have - in front of type number
                    negnew=-newtype
                    temp[negold]=negnew
                oldtypetonewtype.update(temp)
                # STEP 9
                oldtypetonewtypelist.append(oldtypetonewtype)
                # STEP 10
                prevmaxnumberfromkey=maxnumberfromkey
                prevminnumberfromkey=minnumberfromkey
                
            return oldtypetonewtypelist


        def CheckNetChargeIsZero(self,xyzpath,keypath,alzout):
            """
            Intent: Binding simulations require net charge to be 0, otherwise ewald artifcacts. Just use this as sanity check.  
            Input: Tinker XYZ, tinker key, name for tinker analyze output
            Output: - 
            Referenced By: boxsetup.py - BoxSetUpProtocol , productiondynamics.py - SetupProductionDynamics 
            Description:
            1. Call CallAnalyze with xyz,key and name of analyze output. Use option m for checking net charge. 
            2. Read the charge and any residual charge by calling ReadCharge
            3. If the charge is not 0, then raise error
            """
            # STEP 1
            self.CallAnalyze(xyzpath,keypath,alzout,self.trueanalyzepath,'m')
            # STEP 2
            charge,resid=self.ReadCharge(alzout)
            if charge!=0:
                # STEP 3
                raise ValueError('Net charge is not zero! Net Charge = '+str(charge)+' '+keypath)


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

        def DetermineXYZEditOptions(self):
            """
            Intent: Parse xyzedit options via hardcoded text to determine which values to use for xyzedit
            Input: -
            Output: -
            Referenced By:  MolecularDynamics
            Description: 
            1. Call xyzedit
            2. Parse output for keywords and save options
            """
            xyzout='xyzout.txt'
            blanktxt='blankxyzeditinput.txt'
            temp=open(blanktxt,'w')
            temp.write('\n')
            temp.close()
            cmdstr=self.xyzeditpath+' '+os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/VersionFiles/'+'water.xyz'+' '+'<'+' '+blanktxt+' '+'>'+' '+xyzout
            returned_value = subprocess.call(cmdstr, shell=True)
            temp=open(xyzout,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                if self.xyzedittranslatestring in line:
                    self.xyzedittranslatevalue=line.split()[0].replace('(','').replace(')','')
                if self.xyzeditstraystring in line:
                    self.xyzeditstrayvalue=line.split()[0].replace('(','').replace(')','')
                if self.xyzeditappendstring in line:
                    self.xyzeditappendvalue=line.split()[0].replace('(','').replace(')','')
                if self.xyzeditperiodicstring in line:
                    self.xyzeditperiodicvalue=line.split()[0].replace('(','').replace(')','')
                if self.xyzeditsoakstring in line:
                    self.xyzeditsoakvalue=line.split()[0].replace('(','').replace(')','')
                if self.xyzeditionstring in line:
                    self.xyzeditionvalue=line.split()[0].replace('(','').replace(')','')


        def ReadGPUPerformance(self,outputpath):
            temp=open(outputpath,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                if 'Performance' in line:
                    linesplit=line.split()
                    return float(linesplit[-1])

        def ETAString(self,ETA):
            day = ETA // (24 * 3600)
            ETA = ETA % (24 * 3600)
            hour = ETA // 3600
            ETA %= 3600
            minutes = ETA // 60
            ETA %= 60
            seconds = ETA
            ETA_string="d:h:m:s-> %d:%d:%d:%d" % (day, hour, minutes, seconds)
            return ETA_string


        def GrabOutputFileTypeAndTime(self,filename,jobtype,index,outputlist):
            if '_NVT.out' in filename:
                if index==len(outputlist)-2: # last NVT
                    totaltime=self.lastNVTequiltime
                    outputfiletype=jobtype+'_'+'LastNVTEquilibriation'
                    freq=1
                    serial='True'
                else:
                    totaltime=self.equiltimeNVT
                    outputfiletype=jobtype+'_'+'InitialNVTEquilibriation'
                    freq=1
                    serial='True'

            elif '_NPT.out' in filename:
                totaltime=self.equiltimeNPT
                outputfiletype=jobtype+'_'+'NPTEquilibriation'
                freq=1
                serial='True'

            elif 'E' in filename and 'V' in filename and 'Gas' not in filename:
                totaltime=self.proddyntime
                outputfiletype=jobtype+'_'+'ProductionDynamics'
                freq=len(self.estatlambdascheme[0])+len(self.vdwlambdascheme[0])
                serial='False'

            else:
                totaltime=None
                outputfiletype=None
                freq=None
                serial=None


            return totaltime,outputfiletype,freq,serial


        def ReportETA(self):
            if os.path.isfile(self.etafilename):
                with open(self.etafilename) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    line_count = 0
                    totaleta=0
                    for row in csv_reader:
                        if line_count == 0:
                            pass
                        else:
                            jobtype,eta,freq,serial=row
                            tempeta=float(eta)*int(freq)
                            if serial=='False':
                                tempeta=tempeta/self.numbergpus
                            ETA_string=self.ETAString(tempeta)
                            self.WriteToLog('JobType: '+jobtype+' ETA: '+ETA_string)
                            totaleta+=tempeta
                        line_count += 1
                ETA_string=self.ETAString(totaleta)
                self.WriteToLog('Number of GPUs '+str(self.numbergpus))
                self.WriteToLog('Total ETA: '+ETA_string)


        def EstimateDynamicTime(self):
            if self.externalapi!=None:
                dynamicpath=self.dynamicommpath
            else:
                dynamicpath=self.truedynamicpath
            nvt=self.productiondynamicsNVT
            ensemble=self.proddynensem
            proddyntimestep=self.proddyntimestep # fs
            totaltime=.01 * 1000000 # fs
            self.proddynwritefreq=(totaltime*.001) / 10
            steps=int(totaltime/proddyntimestep)
            outputpaths=[]
            for i in range(len(self.minboxfilename)):
                minboxfilename=self.minboxfilename[i]
                key=self.configkeyfilename[i][0]
                outputpath='QuickMD_'+str(i)+'.out'
                outputpaths.append(outputpath)
                if not os.path.isfile(outputpath):
                    cmdstr=prod.ProductionDynamicsCommand(self,minboxfilename,key,steps,ensemble,outputpath,dynamicpath,proddyntimestep,nvt)
                    submit.call_subsystem(self,cmdstr,wait=True)
            jobtypetoperf={}
            for i,outputpath in enumerate(outputpaths):
                if i==0 and self.binding is True:
                    jobtype='Comp'
                elif i==1 and self.binding is True:
                    jobtype='Solv'
                else: # solvation
                    jobtype='Solv'
                perf=self.ReadGPUPerformance(outputpath)
                jobtypetoperf[jobtype]=perf

            fields=['JobType','ETA','Occurance','Serial']
            with open(self.etafilename, 'w') as csvfile:  
                csvwriter = csv.writer(csvfile)  
                csvwriter.writerow(fields)
                for i,outputlist in enumerate(self.equiloutputarray):
                    if i==0 and self.binding is True:
                        jobtype='Comp'
                    elif i==1 and self.binding is True:
                        jobtype='Solv'
                    else: # solvation
                        jobtype='Solv'
                    perf=jobtypetoperf[jobtype]
                    outputfiletypetoeta={}
                    outputfiletypetofreq={}
                    outputfiletypetoserial={}
                    for j,outputfile in enumerate(outputlist):
                        head,tail=os.path.split(outputfile)
                        totaltime,outputfiletype,freq,serial=self.GrabOutputFileTypeAndTime(tail,jobtype,j,outputlist)
                        ETA=(totaltime/perf)*86400 # seconds
                        ETA_string=self.ETAString(ETA)
                        outputfiletypetoeta[outputfiletype]=ETA
                        outputfiletypetofreq[outputfiletype]=freq
                        outputfiletypetoserial[outputfiletype]=serial
                for j in range(len(self.simfoldname)):
                    simfoldname=self.simfoldname[j]
                    perf=jobtypetoperf[simfoldname]
                    proddynoutfilepathlistoflist=self.proddynoutfilepath[j]
                    for k in range(len(proddynoutfilepathlistoflist)):
                        proddynoutfilepath=proddynoutfilepathlistoflist[k]
                        for i in range(len(proddynoutfilepath)):
                            outputfilepath=proddynoutfilepath[i]
                            head,tail=os.path.split(outputfilepath)
                            totaltime,outputfiletype,freq,serial=self.GrabOutputFileTypeAndTime(tail,jobtype,i,proddynoutfilepath)
                            if totaltime!=None:
                                ETA=(totaltime/perf)*86400 # seconds
                                ETA_string=self.ETAString(ETA)
                                outputfiletypetoeta[outputfiletype]=ETA
                                outputfiletypetofreq[outputfiletype]=freq
                                outputfiletypetoserial[outputfiletype]=serial
                
                for outputfiletype,ETA in outputfiletypetoeta.items():
                    freq=outputfiletypetofreq[outputfiletype]
                    serial=outputfiletypetoserial[outputfiletype]
                    row=[outputfiletype,ETA,freq,serial]
                    csvwriter.writerow(row)
            
            



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
        4. If email is given as input send email report of crash 
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
            # STEP 4
            if poltype.email!=None:
                password='amoebaisbest'
                fromaddr = 'poltypecrashreportnoreply@gmail.com'
                toaddr = poltype.email
                filename=poltype.logfname
                poltype.WriteToLog(text)
                poltype.WriteToLog('Poltype has crashed!')
                try:
                    poltype.SendReportEmail(text,fromaddr,toaddr,password,filename,'Poltype Crash Report ')
                except:
                    pass
            raise ValueError('Houston, we have a problem. Buy a developer some coffee!')
    RunPoltype()
