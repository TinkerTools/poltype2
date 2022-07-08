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
import modifiedresidues as modres
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
import forcebalancepoltypeinterface as fb
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
import plotFBresults
import annihilation as ann
import tables 
import boxsetup as box
import time
import pdbxyz
import restraints
import plots
import mutation as mutate
import submitjobs as submit
import keyfilemodifications as keymods
import re
from scipy.optimize import fmin
import pylab as plt
from scipy.interpolate import interp1d
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdFMCS
import parametercomparison
from PyAstronomy import pyasl
from dataclasses import dataclass,field
from pathlib import Path
from rdkit.Chem.MolStandardize import rdMolStandardize
from itertools import groupby
from operator import itemgetter


@dataclass
class PolarizableTyper():
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
        modifiedproteinpdbname:None=None
        unmodifiedproteinpdbname:None=None
        mutatedsidechain:None=None
        mutatedresiduenumber:None=None
        modifiedresiduepdbcode:str='MOD'
        ModifiedResiduePrmPath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/ModifiedResidue.prm'
        SMARTSToTypelibpath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ModifiedResidueLibraries/SMARTSToTypeLib.txt'
        topologylibpath:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ModifiedResidueLibraries/residue_connect.txt"
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
        pymolenvname:str='oldpy3'
        prepdockscript:str=os.path.join(os.path.split(__file__)[0],'preparedockingfiles.py')
        pocketscript:str=os.path.join(os.path.split(__file__)[0],'pocketdetection.py')
        indextompoleframefile:None=None
        qmrelativeweight:float=.5
        liqrelativeweight:float=2.5
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
        bgnstatestructurefile:None=None
        endstatestructurefile:None=None
        redobar:bool=False
        rotateframes:bool=False
        perturbedkeyfilename:None=None
        writeinputfiles:bool=False
        compareparameters:bool=False
        energycutoff:float=.5
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
        mutlambdascheme:list=field(default_factory=lambda : [])
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
        forcebalancejobsdir:None=None
        fixvdwtyperadii:list=field(default_factory=lambda : [])
        maxjobsatsametime:float=10
        liquid_equ_time:float=.5
        gas_equ_time:float=.5
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
        onlyvdwatomlist:None=None
        poltypepathlist:None=None
        vdwtypeslist:None=None
        fittypestogether:None=None
        csvexpdatafile:None=None
        liquid_equ_steps:float=500000
        liquid_prod_steps:float=5000000
        liquid_timestep:float=1.0
        liquid_interval:float=1
        gas_equ_steps:float=500000
        gas_prod_steps:float=5000000
        gas_timestep:float=1.0
        gas_interval:float=1
        md_threads:int=4
        liquid_prod_time:float=5
        gas_prod_time:float=5
        WQ_PORT:None=None
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
        toroptdebugmode:bool=False
        debugmode:bool=False
        fragmenterdebugmode:bool=False
        jobsatsametime:int=0
        usepoleditframes:bool=False
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
        addhydrogentononcharged:bool=True
        accuratevdwsp:bool=False
        inputmoleculefolderpaths:None=None
        email:None=None
        firstoptfinished:bool=False
        optonly:bool=False
        onlyvdwatomindex:None=None
        use_qmopt_vdw:bool=False
        use_gau_vdw:bool=False
        dontusepcm:bool=False
        deleteallnonqmfiles:bool=True
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
        latestsmallmoleculesmartstotypespolarize:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarcommenttoparameters.txt'
        latestsmallmoleculesmartstotinkerclass:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21smartstoclass.txt'
        latestsmallmoleculeprmlib:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21.prm'
        boltzmantemp:float=8
        dovdwscan:bool=False
        vdwprobepathname:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/VdwProbes/'
        vdwprobenames:list=field(default_factory=lambda : ['water'])
        use_gausgeomoptonly:bool=False
        maxtorRMSPDRel:float=.2
        vdwmissingfilename:str='missingvdw.txt'
        databaseprmfilename:str='database.prm'
        tortor:bool=False
        torfit2Drotonly:bool=False
        torfit1Drotonly:bool=False
        externalparameterdatabase:str=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'externalparameterdatabase.txt'
        fitfirsttorsionfoldphase:bool=False
        keyfiletoaddtodatabase:None=None
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
        iodinetoroptbasissetfile:str='def2-svp.1.gbs'
        iodineoptbasissetfile:str='def2-svp.1.gbs'
        iodinedmabasissetfile:str='def2-svp.1.gbs'
        iodineespbasissetfile:str='def2-tzvpp.1.gbs'
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
        espbasisset:str="aug-cc-pVTZ"
        torspbasisset:str="6-311+G*"
        optmethod:str='MP2'
        toroptmethod:str='wB97X-D'
        torspmethod:str='wB97X-D'
        dmamethod:str='MP2'
        espmethod:str='MP2'
        qmonly:bool = False
        espfit:bool = True
        parmtors:bool = True
        foldnum:int=6
        foldoffsetlist:list = field(default_factory=lambda : [ 0.0, 180.0, 0.0, 180.0, 0.0, 180.0 ])
        torlist:None = None
        rotbndlist:None = None
        maxRMSD:float=1
        maxRMSPD:float=1
        maxtorRMSPD:float=1.8
        tordatapointsnum:None=None
        gentorsion:bool=False
        gaustorerror:bool=False
        torsionrestraint:float=.5*3282.80354574
        onlyrotbndslist:list=field(default_factory=lambda : [])
        rotalltors:bool=False
        dontdotor:bool=False
        dontdotorfit:bool=False
        toroptpcm:bool=False
        optpcm:bool=False
        torsppcm:bool=False
        use_gaus:bool=False
        use_gausoptonly:bool=False
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
        def __post_init__(self): 
        
            self.knownresiduesymbs=['A', 'G', 'U', 'A3', 'A5','DA', 'DC', 'DG', 'DT', 'G3', 'G5', 'U3', 'U5','DA3', 'DA5', 'DC3', 'DC5', 'DG3', 'DG5', 'DT3', 'DT5','ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIS','ILE', 'LEU', 'LYD', 'LYS', 'MET','PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYD', 'TYR', 'VAL', 'CYD', 'CYS']
            self.molstructfname=self.structure 
            self.simpath=os.getcwd()        
            self.simname=os.path.basename(self.simpath)
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
                        elif 'pdbcode' in newline:
                            self.pdbcode=a
                        elif 'indextotypefile' in newline:
                            self.indextotypefile=a
                        elif 'indextompoleframefile' in newline:
                            self.indextompoleframefile=a
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
                        elif "writeoutmultipole" in newline:
                            if '=' not in line:
                                self.writeoutmultipole = True
                            else:
                                self.writeoutmultipole=self.GrabBoolValue(a)
                        elif "writeoutbond" in newline:
                            if '=' not in line:
                                self.writeoutbond = True
                            else:
                                self.writeoutbond=self.GrabBoolValue(a)
                        elif "writeoutangle" in newline:
                            if '=' not in line:
                                self.writeoutangle = True
                            else:
                                self.writeoutangle=self.GrabBoolValue(a)
                        
                        elif "optloose" in newline:
                            if '=' not in line:
                                self.optloose = True
                            else:
                                self.optloose=self.GrabBoolValue(a)

                        elif "writeoutstrbnd" in newline:
                            if '=' not in line:
                                self.writeoutstrbnd = True
                            else:
                                self.writeoutstrbnd=self.GrabBoolValue(a)
                        elif "writeoutopbend" in newline:
                            if '=' not in line:
                                self.writeoutopbend = True
                            else:
                                self.writeoutopbend=self.GrabBoolValue(a)
                        elif "writeoutvdw" in newline:
                            if '=' not in line:
                                self.writeoutvdw = True
                            else:
                                self.writeoutvdw=self.GrabBoolValue(a)
                        elif "writeoutpolarize" in newline:
                            if '=' not in line:
                                self.writeoutpolarize = True
                            else:
                                self.writeoutpolarize=self.GrabBoolValue(a)
                        elif "writeouttorsion" in newline:
                            if '=' not in line:
                                self.writeouttorsion = True
                            else:
                                self.writeouttorsion=self.GrabBoolValue(a)


                        elif "usevinardo" in newline:
                            if '=' not in line:
                                self.usevinardo = True
                            else:
                                self.usevinardo=self.GrabBoolValue(a)
                        elif "usevina" in newline:
                            if '=' not in line:
                                self.usevina = True
                            else:
                                self.usevina=self.GrabBoolValue(a)

                        elif "usegold" in newline:
                            if '=' not in line:
                                self.usegold = True
                            else:
                                self.usegold=self.GrabBoolValue(a)
                        elif "usead4" in newline:
                            if '=' not in line:
                                self.usead4 = True
                            else:
                                self.usead4=self.GrabBoolValue(a)


                        elif "addphysioions" in newline:
                            if '=' not in line:
                                self.addphysioions = True
                            else:
                                self.addphysioions=self.GrabBoolValue(a)
                        elif "usepdb2pqr" in newline:
                            if '=' not in line:
                                self.usepdb2pqr = True
                            else:
                                self.usepdb2pqr=self.GrabBoolValue(a)
                        elif "submitlocally" in newline:
                            if '=' not in line:
                                self.submitlocally = True
                            else:
                                self.submitlocally=self.GrabBoolValue(a)
                            self.didinputsubmitlocally=True

                        elif "usesymtypes" in newline:
                            if '=' not in line:
                                self.usesymtypes = True
                            else:
                                self.usesymtypes=self.GrabBoolValue(a)

                        elif "printjobsleft" in newline:
                            if '=' not in line:
                                self.printjobsleft = True
                            else:
                                self.printjobsleft=self.GrabBoolValue(a)
                        elif "salthfe" in newline:
                            if '=' not in line:
                                self.salthfe = True
                            else:
                                self.salthfe=self.GrabBoolValue(a)
                        elif "relaxFBbox" in newline:
                            if '=' not in line:
                                self.relaxFBbox = True
                            else:
                                self.relaxFBbox=self.GrabBoolValue(a)

                        elif "extractinterforbinding" in newline:
                            if '=' not in line:
                                self.extractinterforbinding = True
                            else:
                                self.extractinterforbinding=self.GrabBoolValue(a)


                        elif "neatliquidsim" in newline:
                            if '=' not in line:
                                self.neatliquidsim = True
                            else:
                                self.neatliquidsim=self.GrabBoolValue(a)

                        elif "redobar" in newline:
                            if '=' not in line:
                                self.redobar = True
                            else:
                                self.redobar=self.GrabBoolValue(a)
                        elif "rotateframes" in newline:
                            if '=' not in line:
                                self.rotateframes = True
                            else:
                                self.rotateframes=self.GrabBoolValue(a)
                        elif "checktraj" in newline:
                            if '=' not in line:
                                self.checktraj = True
                            else:
                                self.checktraj=self.GrabBoolValue(a)
                        elif "cpujobsonly" in newline:
                            if '=' not in line:
                                self.cpujobsonly = True
                            else:
                                self.cpujobsonly=self.GrabBoolValue(a)
                        elif "writeinputfiles" in newline:
                            if '=' not in line:
                                self.writeinputfiles = True
                            else:
                                self.writeinputfiles=self.GrabBoolValue(a)
                        elif "compareparameters" in newline:
                            if '=' not in line:
                                self.compareparameters = True
                            else:
                                self.compareparameters=self.GrabBoolValue(a)

                        elif 'expfreeenergy' in newline:
                            self.expfreeenergy=a
                        elif 'pathtosims' in newline:
                            self.pathtosims=a
                        elif 'perturbedkeyfilename' in newline:
                            self.perturbedkeyfilename=a
                        elif 'bgnstatestructurefile' in newline:
                            self.bgnstatestructurefile=a
                        elif 'endstatestructurefile' in newline:
                            self.endstatestructurefile=a
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
                            if '=' not in line:
                                self.usepreequilibriatedbox = True
                            else:
                                self.usepreequilibriatedbox=self.GrabBoolValue(a)
                        elif "equilfinished" in newline:
                            if '=' not in line:
                                self.equilfinished = True
                            else:
                                self.equilfinished=self.GrabBoolValue(a)
                        elif "barfilesfinished" in newline:
                            if '=' not in line:
                                self.barfilesfinished = True
                            else:
                                self.barfilesfinished=self.GrabBoolValue(a)

                        elif "minfinished" in newline:
                            if '=' not in line:
                                self.minfinished = True
                            else:
                                self.minfinished=self.GrabBoolValue(a)
                        elif "boxonly" in newline:
                            if '=' not in line:
                                self.boxonly = True
                            else:
                                self.boxonly=self.GrabBoolValue(a)
                        elif "fep" in newline:
                            if '=' not in line:
                                self.fep = True
                            else:
                                self.fep=self.GrabBoolValue(a)
                        elif "generateinputfilesonly" in newline:
                            if '=' not in line:
                                self.generateinputfilesonly = True
                            else:
                                self.generateinputfilesonly=self.GrabBoolValue(a)
                        elif "productiondynamicsNVT" in newline:
                            if '=' not in line:
                                self.productiondynamicsNVT = True
                            else:
                                self.productiondynamicsNVT=self.GrabBoolValue(a)
                        elif "usetinkerforthermoprops" in newline:
                            if '=' not in line:
                                self.usetinkerforthermoprops = True
                            else:
                                self.usetinkerforthermoprops=self.GrabBoolValue(a)
                        elif 'binding' in newline:
                            self.binding=True
                            self.solvation=True
                            self.complexation=True
                        elif "changegasphaseintegrator" in newline:
                            if '=' not in line:
                                self.changegasphaseintegrator = True
                            else:
                                self.changegasphaseintegrator=self.GrabBoolValue(a)
                        elif "prodmdfinished" in newline:
                            if '=' not in line:
                                self.prodmdfinished = True
                            else:
                                self.prodmdfinished=self.GrabBoolValue(a)
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
                            if '=' not in line:
                                self.flatbotrest = True
                            else:
                                self.flatbotrest=self.GrabBoolValue(a)

                        elif "useproddyngrprests" in newline:
                            if '=' not in line:
                                self.useproddyngrprests = True
                            else:
                                self.useproddyngrprests=self.GrabBoolValue(a)
                        elif "restrainreceptorligand" in newline:
                            if '=' not in line:
                                self.restrainreceptorligand = True
                            else:
                                self.restrainreceptorligand=self.GrabBoolValue(a)
                        elif "restrainreceptorligand" in newline:
                            if '=' not in newline:
                                self.restrainreceptorligand = True
                            else:
                                self.restrainreceptorligand=self.GrabBoolValue(a)
                        elif "annihilatevdw" in newline:
                            if '=' not in line:
                                self.annihilatevdw = True
                            else:
                                self.annihilatevdw=self.GrabBoolValue(a)
                        elif "proddynensem" in newline:
                            self.proddynensem = a
                        elif ("keyfilenamelist") in newline:
                            self.keyfilenamelist=a.split(',')
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
                            self.ligandxyzfilenamelist=a.split(',')
                        elif ("annihilateligandxyzfilenamelist") in newline and "receptor" not in newline:
                            self.annihilateligandxyzfilenamelist=a.split(',')
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
                        elif ("mutlambdascheme") in newline:
                            self.mutlambdascheme=commalist
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
                        elif ("equiltimeNPT") in newline:
                            self.equiltimeNPT= float(a)
                        elif ("equiltimeionNPT") in newline:
                            self.equiltimeionNPT= float(a)
                        elif ("proddynwritefreq") in newline:
                            self.proddynwritefreq= a
                        elif ("proddyntime") in newline:
                            self.proddyntime= float(a)
                            self.inputproddyntime=True
                        elif ("complexation") in newline:
                            self.complexation=True
                        elif ("solvation") in newline:
                            self.solvation=True
                        elif ("equiltimeNVT") in newline:
                            self.equiltimeNVT= float(a)
                        elif ("lastNVTequiltime") in newline:
                            self.lastNVTequiltime= float(a)
                        elif ("equilonly") in newline:
                            self.equilonly=True
                        elif ("minonly") in newline:
                            self.minonly=True
                        elif ("equilrestlambdascheme") in newline:
                            self.equilrestlambdascheme=commalist

                        elif "rotalltors" in newline:
                            if '=' not in line:
                                self.rotalltors = True
                            else:
                                self.rotalltors=self.GrabBoolValue(a)


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
                            if '=' not in line:
                                self.addlonepairvdwsites = True
                            else:
                                self.addlonepairvdwsites=self.GrabBoolValue(a)


                        elif "genprotstatesonly" in newline:
                            if '=' not in line:
                                self.genprotstatesonly = True
                            else:
                                self.genprotstatesonly=self.GrabBoolValue(a)

                        elif "quickdatabasesearch" in newline:
                            if '=' not in line:
                                self.quickdatabasesearch = True
                            else:
                                self.quickdatabasesearch=self.GrabBoolValue(a)

                        elif "fitred" in newline:
                            if '=' not in line:
                                self.fitred = True
                            else:
                                self.fitred=self.GrabBoolValue(a)

                        elif "onlyvdwatomlist" in newline:
                            self.onlyvdwatomlist=a.split(',')
                            self.onlyvdwatomlist=[i.strip() for i in self.onlyvdwatomlist]
                            self.onlyvdwatomlist=[int(i) for i in self.onlyvdwatomlist]

                        elif "usepoleditframes" in newline:
                            if '=' not in line:
                                self.usepoleditframes = True
                            else:
                                self.usepoleditframes=self.GrabBoolValue(a)

                        elif "generateextendedconf" in newline:
                            if '=' not in line:
                                self.generateextendedconf = True
                            else:
                                self.generateextendedconf=self.GrabBoolValue(a)


                        elif "addhydrogens" in newline:
                            if '=' not in line:
                                self.addhydrogens = True
                            else:
                                self.addhydrogens=self.GrabBoolValue(a)

                        elif "nonaroringtor1Dscan" in newline:
                            if '=' not in line:
                                self.nonaroringtor1Dscan = True
                            else:
                                self.nonaroringtor1Dscan=self.GrabBoolValue(a)

                        elif "fragmenterdebugmode" in newline:
                            if '=' not in line:
                                self.fragmenterdebugmode = True
                            else:
                                self.fragmenterdebugmode=self.GrabBoolValue(a)
                        elif "skipespfiterror" in newline:
                            if '=' not in line:
                                self.skipespfiterror = True
                            else:
                                self.skipespfiterror=self.GrabBoolValue(a)

                        elif "fragmentjobslocal" in newline:
                            if '=' not in line:
                                self.fragmentjobslocal = True
                            else:
                                self.fragmentjobslocal=self.GrabBoolValue(a)
                        elif "smallmoleculefragmenter" in newline:
                            if '=' not in line:
                                self.smallmoleculefragmenter = True
                            else:
                                self.smallmoleculefragmenter=self.GrabBoolValue(a)


                        elif "debugmode" in newline and 'fragmenterdebugmode' not in line and 'tordebugmode' not in line and 'toroptdebugmode' not in line:
                            if '=' not in line:
                                self.debugmode = True
                            else:
                                self.debugmode=self.GrabBoolValue(a)

                        elif "toroptdebugmode" in newline:
                            if '=' not in line:
                                self.toroptdebugmode = True
                            else:
                                self.toroptdebugmode=self.GrabBoolValue(a)

                        elif "totalcharge" in newline:
                            self.totalcharge=int(a)
                        elif "lastlogfileupdatetime" in newline:
                            self.lastlogfileupdatetime=int(a)

                        elif "numespconfs" in newline:
                            self.numespconfs=int(a)

                        elif "consumptionratio" in newline:
                            self.consumptionratio=float(a)
                        elif "liquid_equ_time" in newline:
                            self.liquid_equ_time=float(a)
                        elif "gas_equ_time" in newline:
                            self.gas_equ_time=float(a)

                        elif "jobsatsametime" in newline and 'max' not in newline and 'parentjobsatsametime' not in newline:
                            self.jobsatsametime=int(a)
                        elif "esprestweight" in newline:
                            self.esprestweight=float(a)
                        elif "espgrad" in newline:
                            self.espgrad=a
                        elif "forcebalancejobsdir" in newline:
                            self.forcebalancejobsdir=a
                        elif "checkinputonly" in newline:
                            self.checkinputonly=True
                        elif "setupfragjobsonly" in newline:
                            if '=' not in line:
                                self.setupfragjobsonly = True
                            else:
                                self.setupfragjobsonly=self.GrabBoolValue(a)

                        elif "databasematchonly" in newline:
                            if '=' not in line:
                                self.databasematchonly = True
                            else:
                                self.databasematchonly=self.GrabBoolValue(a)

                        elif "onlyvdwatomindex" in newline:
                            self.onlyvdwatomindex=int(a)
                        elif "parentjobsatsametime" in newline:
                            self.parentjobsatsametime=int(a)
                        elif "coresperjob" in newline:
                            self.coresperjob=int(a)
                        elif "deleteallnonqmfiles" in newline:
                            if '=' not in line:
                                self.deleteallnonqmfiles = True
                            else:
                                self.deleteallnonqmfiles=self.GrabBoolValue(a)

                        elif "addhydrogentononcharged" in newline:
                            if '=' not in line:
                                self.addhydrogentononcharged = True
                            else:
                                self.addhydrogentononcharged=self.GrabBoolValue(a)


                        elif "firstoptfinished" in newline:
                            if '=' not in line:
                                self.firstoptfinished = True
                            else:
                                self.firstoptfinished=self.GrabBoolValue(a)
                        elif "accuratevdwsp" in newline:
                            if '=' not in line:
                                self.accuratevdwsp = True
                            else:
                                self.accuratevdwsp=self.GrabBoolValue(a)

                        elif "homodimers" in newline:
                            if '=' not in line:
                                self.homodimers = True
                            else:
                                self.homodimers=self.GrabBoolValue(a)

                        elif "optonly" in newline and 'gaus' not in newline:
                            if '=' not in line:
                                self.optonly = True
                            else:
                                self.optonly=self.GrabBoolValue(a)

                        elif "use_qmopt_vdw" in newline:
                            if '=' not in line:
                                self.use_qmopt_vdw = True
                            else:
                               self.use_qmopt_vdw=self.GrabBoolValue(a)


                        elif "use_gausgeomoptonly" in newline:
                            if '=' not in line:
                                self.use_gausgeomoptonly = True
                            else:
                                self.use_gausgeomoptonly=self.GrabBoolValue(a)
                        elif "use_gau_vdw" in newline:
                            if '=' not in line:
                                self.use_gau_vdw = True
                            else:
                                self.use_gau_vdw=self.GrabBoolValue(a)


                        elif "tortor" in newline and 'only' not in newline:
                            if '=' not in line:
                                self.tortor = True
                            else:
                                self.tortor=self.GrabBoolValue(a)

                        elif "tordebugmode" in newline:
                            if '=' not in line:
                                self.tordebugmode = True
                            else:
                                self.tordebugmode=self.GrabBoolValue(a)
                        elif "keyfiletoaddtodatabase" in newline:
                            self.keyfiletoaddtodatabase=a
                        elif "email" in newline:
                            self.email=a
                        elif "refinenonaroringtors" in newline:
                            if '=' not in line:
                                self.refinenonaroringtors = True
                            else:
                                self.refinenonaroringtors=self.GrabBoolValue(a)
                        elif "fitfirsttorsionfoldphase" in newline:
                            if '=' not in line:
                                self.fitfirsttorsionfoldphase = True
                            else:
                                self.fitfirsttorsionfoldphase=self.GrabBoolValue(a)

                        elif "fitqmdipole" in newline:
                            if '=' not in line:
                                self.fitqmdipole = True
                            else:
                                self.fitqmdipole=self.GrabBoolValue(a)
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
                            if '=' not in line:
                                self.use_gauPCM = True
                            else:
                                self.use_gauPCM=self.GrabBoolValue(a)

                        elif "espfit" in newline:
                            if '=' not in line:
                                self.espfit = True
                            else:
                                self.espfit=self.GrabBoolValue(a)


                        elif 'poltypepath' in newline and 'poltypepathlist' not in newline:
                            self.poltypepath=a
                        elif 'paramhead' in newline:
                            self.paramhead = a
                        elif 'inputmoleculefolderpaths' in newline:
                            self.inputmoleculefolderpaths=a
                        elif 'WBOtol' in newline:
                            self.WBOtol=float(a)

                        elif 'scfmaxiter' in newline:
                            self.scfmaxiter=a
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
                        elif "externalapi" in newline and a!='None':
                            self.externalapi=a
                        elif "gausoptcoords" in newline:
                            self.gausoptcoords = a
                        elif "suppresstorfiterr" in newline:
                            if '=' not in line:
                                self.suppresstorfiterr=True
                            else:
                                self.suppresstorfiterr=self.GrabBoolValue(a)
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
                        elif "bashrcpath" in newline and a!='None':
                            self.bashrcpath = a
                        elif "modifiedproteinpdbname" in newline:
                            self.modifiedproteinpdbname = a
                            self.molstructfname='ModifiedRes.sdf'
                        elif "structure" in newline:
                            self.molstructfname = a
                        elif "dontusepcm" in newline:
                            if '=' not in line:
                                self.dontusepcm = True
                            else: 
                                self.dontusepcm=self.GrabBoolValue(a)

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

                        elif "dovdwscan" in newline:
                            if '=' not in line:
                                self.dovdwscan = True
                            else:
                                self.dovdwscan=self.GrabBoolValue(a)


                        elif "dontdotorfit" in newline:
                            if '=' not in line:
                                self.dontdotorfit = True
                            else:
                                self.dontdotorfit=self.GrabBoolValue(a)
                        elif "optmaxcycle" in newline:
                            self.optmaxcycle = int(a)
                        elif "torsionrestraint" in newline:
                            self.torsionrestraint=float(a)
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
                            if '=' not in line:
                                self.qmonly = True
                            else:
                                self.qmonly = self.GrabBoolValue(a)
                        elif "sleeptime" in newline:
                            self.sleeptime = float(a)
                        elif  "poltypepathlist" in newline:
                            self.poltypepathlist=a.split(',')
                        elif "csvexpdatafile" in newline:
                            self.csvexpdatafile=a
                        elif "fittypestogether" in newline:
                            self.fittypestogether=self.ReturnListOfList(a)
                        elif "vdwtypeslist" in newline:
                            self.vdwtypeslist=self.ReturnListOfList(a)
                        elif "liquid_equ_steps" in newline:
                            self.liquid_equ_steps=int(a)
                        elif 'liquid_prod_steps' in newline:
                            self.liquid_prod_steps=int(a) 
                        elif 'liquid_timestep' in newline:
                            self.liquid_timestep=int(a)
                        elif 'liquid_interval' in newline:
                            self.liquid_interval=float(a)
                        elif 'gas_equ_steps' in newline:
                            self.gas_equ_steps=int(a)
                        elif 'gas_prod_steps' in newline:
                            self.gas_prod_steps=int(a)
                        elif 'gas_timestep' in newline:
                            self.gas_timestep=float(a)
                        elif 'gas_interval' in newline:
                            self.gas_interval=float(a)
                        elif 'md_threads' in newline:
                            self.md_threads=int(a)
                        elif 'liquid_prod_time' in newline:
                            self.liquid_prod_time=float(a)
                        elif 'gas_prod_time' in newline:
                            self.gas_prod_time=float(a)
                        elif 'WQ_PORT' in newline:
                            self.WQ_PORT=a
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


            if self.jobsatsametime!=0:
                self.maximizejobsatsametime=False
            else:
                if self.isfragjob==True:
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
           


            if self.poltypepathlist!=None:
                fb.GenerateForceBalanceInputs(self,self.poltypepathlist,self.vdwtypeslist,self.liquid_equ_steps,self.liquid_prod_steps,self.liquid_timestep,self.liquid_interval,self.gas_equ_steps,self.gas_prod_steps,self.gas_timestep,self.gas_interval,self.md_threads,self.liquid_prod_time,self.gas_prod_time,self.WQ_PORT,self.csvexpdatafile,self.fittypestogether,self.vdwprmtypestofit,self.vdwtypestoeval,self.liquid_equ_time,self.gas_equ_time,self.qmrelativeweight,self.liqrelativeweight,self.enthalpyrelativeweight,self.densityrelativeweight,self.relaxFBbox)
                sys.exit()
            if self.forcebalancejobsdir!=None:
                plotFBresults.PlotForceBalanceResults(self.forcebalancejobsdir,self.targetdensityerror,self.targetenthalpyerror)
                sys.exit()
            if self.compareparameters==True:
                self.MolecularDynamics()
                sys.exit()
            if self.ligandxyzfilenamelist!=None and (self.binding==True or self.solvation==True or self.neatliquidsim==True) or self.pdbcode!=None or self.usepdb2pqr!=False:
                self.MolecularDynamics()
                sys.exit()
            if (self.usead4==True or self.usegold==True or self.usevina==True or self.usevinardo==True) and self.complexedproteinpdbname!=None:
                docking.DockingWrapper(self,self.complexedproteinpdbname,self.dockgridcenter,self.dockgridsize,self.vinaexhaustiveness,self.nposes,self.goldbin,self.usevina,self.usead4,self.usevinardo,self.usegold,self.dockingenvname,self.prepdockscript,self.gridspacing,self.listofligands,self.ecrexpect)
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
            if self.neatliquidsim==True and self.preequilboxpath==preequilboxpath: # if default water box
                self.usepreequilibriatedbox=False # then make new neat liquid box from scratch
            if self.inputproddyntime==False and self.binding==True:
                self.proddyntime=10
            head,tail=os.path.split(self.prmfilepath)
            newpath=os.path.join(os.getcwd(),tail)
            if self.prmfilepath!=newpath or head!=None:
                shutil.copy(self.prmfilepath,newpath)
            self.prmfilepath=tail
            if self.simulationstostopfolderpath!=None:
                self.StopSimulations(self.simulationstostopfolderpath)
                sys.exit()
            if self.annihilateligandxyzfilenamelist==None:
                self.annihilateligandxyzfilenamelist=self.ligandxyzfilenamelist.copy()
            if int(self.estatlambdascheme[0])!=1 or int(self.vdwlambdascheme[0])!=1 and self.solvation==True:
                raise ValueError('Please start with 1 on left side of lambdascheme and 0 on right. This way for solvation can add extra steps correctly')
            obConversion = openbabel.OBConversion()
            self.ligandsmileslist=[]
            self.annihilateligandsmileslist=[]
            self.ligandxyztosmiles={}
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


            self.logname=os.path.basename(os.getcwd())+'_'+self.logname 
            self.outputpath=os.path.join(os.getcwd(),'')
            self.logfh=open(self.outputpath+self.logname,'a+')
            self.SanitizeMMExecutables()

            foundgpukey=False
            if (self.which(self.dynamicommpath)):
                self.usegpu=True
                foundgpukey=True

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
            self.CleanUpFiles()
            self.CheckTinkerVersion()
            if self.perturbedkeyfilename!=None:
                self.compareparameters=True

            if self.bgnstatexyz!=None and self.bgnstatekey!=None and self.endstatexyz!=None and self.endstatekey!=None and self.compareparameters==True: 
                parametercomparison.CompareBgnEndParameters(self)
                if self.perturbedkeyfilename==None:
                    sys.exit()


            if self.ligandxyzfilenamelist!=None and self.templateligandxyzfilename!=None:
                self.AlignLigandXYZToTemplateXYZ()
                sys.exit()

            if self.binding==False and self.solvation==False and self.neatliquidsim==False and self.usepdb2pqr==False and self.pdbcode==None:
                raise ValueError('Please choose either solvation or binding, or neat liquid simulation mode')

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
            needtoshifttypes=self.CheckIfNeedToShiftTypes(self.keyfilenamelist)
            if needtoshifttypes==True:
                oldtypetonewtypelist=self.GenerateTypeMaps(self.keyfilenamelist)
                for i in range(len(oldtypetonewtypelist)):
                    oldindextonewindex=oldtypetonewtypelist[i]
                    key=self.keyfilenamelist[i]
                    xyz=self.ligandxyzfilenamelist[i]
                    self.ShiftParameterTypes(key,oldindextonewindex)
                    self.ShiftParameterTypes(xyz,oldindextonewindex)
            if self.pdbcode!=None:
                pdbxyz.FillInMissingResidues(self,self.pdbcode)
                sys.exit()

            if self.usepdb2pqr==True and self.uncomplexedproteinpdbname!=None:
                pdbxyz.CallPDB2PQR(self,self.uncomplexedproteinpdbname)
                sys.exit()

            self.originalkeyfilename=self.AppendKeys(self.keyfilenamelist,'ligands.key')
            mpolearrays=self.CheckInputXYZKeyFiles(ligandonly=True) # check xyz/keys of ligands before making complex XYZ
            if self.complexation==True and self.complexedproteinpdbname!=None: 

                pdbxyz.GenerateProteinTinkerXYZFile(self)
            elif self.complexation==True  and self.complexedproteinpdbname==None and self.receptorligandxyzfilename==None: 
                raise ValueError('Missing complexedproteinpdbname, need ligand in complex')
            if self.ligandfilename!=None:
                self.ReadLigandCharge()

            if self.neatliquidsim==True:
                self.estatlambdascheme=[1]
                self.vdwlambdascheme=[1]
                self.restlambdascheme=[1]
                self.restrainatomsduringminimization=False

            self.originalestatlambdascheme=self.estatlambdascheme[:]
            self.originalvdwlambdascheme=self.vdwlambdascheme[:]
            self.originalrestlambdascheme=self.restlambdascheme[:]
            self.originaltorsionrestscheme=self.torsionrestscheme[:]

            


            self.checkneedregrow=self.CheckIfNeedRegrowForSolvation()
            self.CheckReceptorLigandXYZFile() # sometimes user have box info on top
            mpolearrays=self.CheckInputXYZKeyFiles()
            self.addgas=False
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

            if self.ligandcharge==None: # assume since no input given
                self.ligandcharge=0

            head,self.foldername=os.path.split(os.getcwd())
            if self.complexation==True and self.solvation==False:
                self.bufferlen=[20]
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

            if self.pathtosims!=None:
                tables.GrabSimDataFromPathList(self)
                plots.PlotEnergyData(self)
                sys.exit()
            if self.perturbkeyfilelist!=None:
                ls=[i+1 for i in range(len(self.perturbkeyfilelist))]
                self.perturbkeyfilelisttokeyindex=dict(zip(self.perturbkeyfilelist,ls)) 
                if self.numinterpolsforperturbkeyfiles==None:
                    self.numinterpolsforperturbkeyfiles=[1 for i in range(len(self.perturbkeyfilelist))]
                else:
                    self.numinterpolsforperturbkeyfiles=[int(i) for i in self.numinterpolsforperturbkeyfiles]
                self.perturbkeyfilelisttointerpolations=dict(zip(self.perturbkeyfilelist,self.numinterpolsforperturbkeyfiles))
                if self.fep==True:
                    self.submitlocally=True

            for i in range(len(self.newkeyfilename)):
                newkeyfilename=self.newkeyfilename[i]
                self.AppendKeys(self.keyfilenamelist,newkeyfilename)
                if i==0:
                    self.ModifyCharge(newkeyfilename,mpolearrays)
            if self.complexation==True and self.solvation==False:
                self.keyfilename=[[self.foldername+'_comp'+'.key']]
                if self.perturbkeyfilelist!=None:
                    for keyname,index in self.perturbkeyfilelisttokeyindex.items():
                        numinterpols=self.perturbkeyfilelisttointerpolations[keyname]
                        for k in range(numinterpols):
                            newkeyname=self.foldername+'_comp'+'.key'
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
                            newkeyname=self.foldername+'_solv'+'.key'
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
                            newkeyname=self.foldername+'_neatliq'+'.key'
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
                            newkeyname=self.foldername+'_comp'+'.key'
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

            
            
            self.equilstepsNVT=str(int((self.equiltimeNVT*1000000)/self.equiltimestep/len(self.equilibriatescheme)-1)) # convert ns to fs divide by the length of temperature scheme
            self.lastequilstepsNVT=str(int((self.lastNVTequiltime*1000000)/self.equiltimestep))
            self.equilstepsNPT=str(int((self.equiltimeNPT*1000000)/self.equiltimestep))
            self.equilstepsionNPT=str(int((self.equiltimeionNPT*1000000)/self.equiltimestep))
            self.proddynframenum=str(int(float(self.proddyntime)/(float(self.proddynwritefreq)*0.001)))
            self.proddynsteps=str(int((self.proddyntime*1000000)/self.proddyntimestep))
            self.equilframenumNPT=int((float(self.equiltimeNPT))/(float(self.equilwritefreq)*0.001))
            self.equilframenumNVT=int((float(self.equiltimeNVT+self.lastNVTequiltime))/(float(self.equilwritefreq)*0.001))
            self.equilframenum=self.equilframenumNPT+self.equilframenumNVT
            self.proddynframenum=int((self.proddyntime)/(float(self.proddynwritefreq)*0.001))
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
                    templist.append(self.outputpath+self.simname+'_'+simfold+'_'+str(temp)+'_'+self.equilstepsNVT+'_NVT.out')
                    tempsteps.append(self.equilstepsNVT)

                templist.append(self.outputpath+self.simname+'_'+simfold+'_'+str(temp)+'_'+self.equilstepsNPT+'_NPT.out')
                tempsteps.append(self.equilstepsNPT)

                if i==0 and self.solvation==True and self.addsolvionwindows==True:
                    templist.append(self.outputpath+'SolvIon'+'_'+simfold+'_'+str(temp)+'_'+self.equilstepsionNPT+'_NPT.out')
                    tempsteps.append(self.equilstepsionNPT)
                self.equiloutputarray.append(templist)
                self.equiloutputstepsarray.append(tempsteps)


            self.nextfiletofinish=self.equiloutputarray[0]
            if self.productiondynamicsNVT==True:
                self.proddynensem=self.NPVTensem
            else:
                self.proddynensem=self.NPTensem
            mut=False
            self.foldernametolambdakeyfilename={}
            self.foldernametonuminterpols={}
            self.foldernametointerpolindex={}
            if len(self.mutlambdascheme)==0:
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
                    
            else:
                mut=True
                index=0
                for i in range(len(self.mutlambdascheme)):
                    mutlambda=self.mutlambdascheme[i]
                    fold=self.simfoldname+"Mut%s"%(mutlambda)
                    self.lambdafolderlist[index].append(fold)
                    outputfilepath=os.getcwd()+'/'+self.simfoldname+'/'+fold+'/'
                    outputfilepath+=self.foldername+'_'+fold+'.out'
                    self.proddynoutfilepath[index].append(outputfilepath)
            if mut==True:
                mutate.SingleTopologyMutationProtocol(self)
           

            
 
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
            

            self.RotateTranslateInertialFrame()

            for i in range(len(self.tabledict)):
                self.tabledict[i]['Prod MD Ensemb']=self.proddynensem
                self.tabledict[i]['Prod MD Time']=self.proddyntime
                self.tabledict[i]['Prod MD Steps']=self.proddynsteps
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
                    self.tabledict[i]['Receptor Charge']=self.receptorcharge
                    if self.uncomplexedproteinpdbname!=None:
                        self.tabledict[i]['Receptor Name']=self.uncomplexedproteinpdbname.replace('.pdb','')
                    elif self.receptorligandxyzfilename!=None and self.uncomplexedproteinpdbname==None: 
                        self.tabledict[i]['Receptor Name']=self.receptorligandxyzfilename.replace('.xyz','')

            
            tables.WriteTableUpdateToLog(self)
            self.ligandindextoneighbslist=[]
            self.ligandindextosymlist=[]
            self.ligandindextotypenumlist=[]

            for xyz in self.annihilateligandxyzfilenamelist:
                statexyzatominfo,stateindextotypeindex,stateatomnum,indextocoords,indextoneighbs,indextosym=parametercomparison.GrabXYZInfo(self,xyz)
                self.ligandindextoneighbslist.append(indextoneighbs)
                self.ligandindextosymlist.append(indextosym)
                self.ligandindextotypenumlist.append(stateindextotypeindex)

            ann.main(self)





        def CheckIfNeedToShiftTypes(self,keyfilenamelist):
            newlist=keyfilenamelist.copy()
            newlist.append(self.prmfilepath)
            needtoshifttypes=False
            alltypes=[]
            for keyfilename in newlist:
                maxnumberfromkey=self.GrabMaxTypeNumber(keyfilename)
                minnumberfromkey=self.GrabMinTypeNumber(keyfilename)
                types=list(range(minnumberfromkey,maxnumberfromkey+1))
                for typenum in types:
                    if typenum in alltypes:
                        needtoshifttypes=True
                    alltypes.append(typenum)


            return needtoshifttypes



        def AppendKeys(self,keyfilenamelist,thekey):
            temp=open(thekey,'w')
            for key in keyfilenamelist:
                othertemp=open(key,'r')
                results=othertemp.readlines()
                othertemp.close()
                for line in results:
                    temp.write(line)

            temp.close()

            return thekey




        def GrabIonizationStates(self,m):
            smi=Chem.MolToSmiles(m)
            smiles=[smi]
            mols = [Chem.MolFromSmiles(s) for s in smiles]
            for i, mol in enumerate(mols):
                mol.SetProp("msg","Orig SMILES: " + smiles[i])
                
            
            protonated_mols = dimorphite_dl.run_with_mol_list(
                mols,
                min_ph=6.9,
                max_ph=7.2,
            )
            
            finalsmi=[Chem.MolToSmiles(m) for m in protonated_mols]
            obConversion = openbabel.OBConversion()
            
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

            enumerator = rdMolStandardize.TautomerEnumerator()
            smi=Chem.MolToSmiles(m)
            m=Chem.MolFromSmiles(smi) 
            canon = enumerator.Canonicalize(m)
            csmi = Chem.MolToSmiles(canon)
            res = [canon]
            tauts = enumerator.Enumerate(m)
            smis = [Chem.MolToSmiles(x) for x in tauts]
            stpl = sorted((x,y) for x,y in zip(smis,tauts) if x!=csmi)
            res += [y for x,y in stpl]
            obConversion = openbabel.OBConversion()
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
            newlist=string.split(',')
            templist=[]
            for ele in newlist:
                nums=ele.lstrip().rstrip().split()
                temp=[]
                for e in nums:
                    temp.append(int(e))
                templist.append(temp)
            newlist=templist
            return newlist



        def SanitizeAllQMMethods(self):
            self.optmethod=self.SanitizeQMMethod(self.optmethod,True)                 
            self.toroptmethod=self.SanitizeQMMethod(self.toroptmethod,True)                  
            self.torspmethod=self.SanitizeQMMethod(self.torspmethod,False)                    
            self.dmamethod=self.SanitizeQMMethod(self.dmamethod,False)                      
            self.espmethod=self.SanitizeQMMethod(self.espmethod,False)     
           
        def GrabBoolValue(self, value):
            if value.lower() == 'true':
                return True
            if value.lower() == 'false':
                return False
            raise ValueError('Could not convert "{}" into a boolean!'.format(value))

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
 

        def WriteToLog(self,string,prin=False):
            now = time.strftime("%c",time.localtime())
            if not isinstance(string, str):
                string=string.decode("utf-8")
            self.logfh.write(now+' '+string+'\n')
            self.logfh.flush()
            os.fsync(self.logfh.fileno())
            if prin==True:
                print(now+' '+ string + "\n")
            
        def SanitizeMMExecutables(self):
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
            Input:
            Output:
            Referenced By: main
            Description: -
            """
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

            if self.use_gaus or self.use_gausoptonly:
                if self.which(os.path.join(self.gausdir, self.gausexe)) is None:
                    print("ERROR: Invalid Gaussian directory: ", self.gausdir)
                    sys.exit(1)
                self.cubegenexe = os.path.join(self.gausdir,self.cubegenexe)
                self.formchkexe = os.path.join(self.gausdir,self.formchkexe)


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
            
            if self.forcefield.upper() in ["AMOEBAPLUS", "APLUS", "AMOEBA+"]:
                self.paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebaplus21_header.prm"
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
                
                
                

            if self.gdmadir is None:
                self.gdmadir = os.getenv("GDMADIR", default="")
            self.gdmaexe = os.path.join(self.gdmadir, self.gdmaexe) 
            if (not self.which(self.gdmaexe)):
                print("ERROR: Cannot find GDMA executable")
                sys.exit(2)


            if (not self.which('psi4')):
                print("ERROR: Cannot find PSI4 executable")
                sys.exit(2)
             

            if self.use_gaus or self.use_gausoptonly:
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

            else:
                if ("PSI_SCRATCH" in os.environ):
                    self.scratchdir = os.environ["PSI_SCRATCH"]
                    if not os.path.isdir(self.scratchdir):
                        os.mkdir(self.scratchdir)
        
            if self.molstructfname!=None: 
                self.FileNames()

        def FileNames(self):
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


        def assign_filenames (self,filename,suffix):
            if filename in globals():
                return eval(filename)
            else:
                return self.molecprefix + suffix
        
        def printfile(self,filename):
            with open(os.path.abspath(os.path.join(self.poltypepath, os.pardir))+r'/'+filename,'r') as f:
                print(f.read(), end='')
        
        def copyright(self):
            self.printfile(self.versionfile)
        
            
        def CheckIsInput2D(self,mol,obConversion,rdkitmol):
            is2d=True
            for atom in openbabel.OBMolAtomIter(mol):
                zcoord=atom.GetZ()
                if zcoord!=0:
                    is2d=False
            if is2d==True: 
                molprefix=self.molstructfname.split('.')[0]
                newname=molprefix+'_3D'+'.mol'
                self.molstructfname=newname
                self.molecprefix =  os.path.splitext(self.molstructfname)[0]
                self.FileNames()

                rdmolfiles.MolToMolFile(rdkitmol,'test.mol',kekulize=True)
                AllChem.EmbedMolecule(rdkitmol)
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
            p=self.cmdtopid[job]
            try:
                p.wait(timeout=1)
            except subprocess.TimeoutExpired:
                process = psutil.Process(p.pid)
                for proc in process.children(recursive=True):
                    proc.kill()
                process.kill()

        
        def CallJobsSeriallyLocalHost(self,fulljobtooutputlog,skiperrors,wait=False):
           thepath=os.path.join(os.getcwd(),'Fragments')
           finishedjobs=[]
           errorjobs=[]
           submittedjobs=[]
           errormessages=[]
           for job,outputlog in fulljobtooutputlog.items():
               finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
               if error==True: # remove log before resubmitting
                   os.remove(outputlog)

           while len(finishedjobs)!=len(list(fulljobtooutputlog.keys())):
               for job,outputlog in fulljobtooutputlog.items():
                   if job not in finishedjobs:
                      
                      finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
                      
                      if finished==True:
                          self.KillJob(job,outputlog)
                          if outputlog not in finishedjobs:
                              finishedjobs.append(outputlog)
                              self.NormalTerm(outputlog)
                              if job in submittedjobs:
                                  submittedjobs.remove(job) 
                      if error==True and job in submittedjobs:
                          self.KillJob(job,outputlog) # sometimes tinker jobs or qm just hang and dont want zombie process 
                          if outputlog not in finishedjobs:
                              errorjobs.append(outputlog)
                              finishedjobs.append(outputlog) 
                              self.ErrorTerm(outputlog,skiperrors)
                              submittedjobs.remove(job)
                      if job not in submittedjobs and len(submittedjobs)<self.jobsatsametime and finished==False and outputlog not in finishedjobs:
                          count=len(finishedjobs)
                          ratio=round((100*count)/len(fulljobtooutputlog.keys()),2)
                          if len(finishedjobs)!=0:
                              if 'poltype.py' in job:
                                  self.ETAQMFinish(thepath,len(fulljobtooutputlog.keys()))
                          if len(list(fulljobtooutputlog.keys()))>1: 
                              self.WriteToLog('Percent of jobs finished '+str(ratio))
                          self.call_subsystem([job],wait,skiperrors)
                          submittedjobs.append(job)
                       

               time.sleep(self.sleeptime)
           return finishedjobs,errorjobs

        def TabulateLogs(self,thepath):
            listoffilepaths=[]
            os.chdir(thepath)
            files=os.listdir()
            for f in files:
                path=os.path.join(thepath,f)
                if os.path.isdir(path):
                    logfilepaths=self.GrabLogFilePaths(path)
                    if len(logfilepaths)!=0:
                        listoffilepaths.append(logfilepaths)
            return listoffilepaths
        
        def GrabLogFilePaths(self,path):
            filepaths=[]
            for root, subdirs, files in os.walk(path):
                finishedpoltype=False
                for f in files:
                    if 'final.key' in f:
                        finishedpoltype=True
                if finishedpoltype==False:
                    continue
                for f in files:
                    if '.log' in f and '_frag' not in f and 'poltype' not in f:
                        filepath=os.path.join(path, f)  
                        filepaths.append(filepath)
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
            timestringsplit=timestring.split(':')
            hours=float(timestringsplit[0])
            minutes=float(timestringsplit[1])*(1/60)
            seconds=float(timestringsplit[2])*(1/60)*(1/60)
            totaltime=hours+minutes+seconds
            return totaltime
        
        def ConvertTimeStringGaussian(self,line):
            linesplit=line.split()
            days=float(linesplit[3])*24
            hours=float(linesplit[5])
            minutes=float(linesplit[7])*(1/60)
            seconds=float(linesplit[9])*(1/60)*(1/60)
            walltime=days+hours+minutes+seconds
            return walltime
        
        def TabulateCPUTimes(self,listoffilepaths):
            n=5
            listoflistofcputimes=[]
            for filepathsidx in range(len(listoffilepaths)):
                filepaths=listoffilepaths[filepathsidx]
                cputimearray=[]
                for logfile in filepaths:
                    if os.path.isfile(logfile):
                        with open(logfile) as frb:
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
        
                listoflistofcputimes.append(cputimearray)
            return listoflistofcputimes
        
        
        def ETAQMFinish(self,thepath,numberofmolstotal):
            listoffilepaths=self.TabulateLogs(thepath)
            listoflistofcputimes=self.TabulateCPUTimes(listoffilepaths)
            listofcputimes=[sum(ls) for ls in listoflistofcputimes]
            if len(listofcputimes)==0:
                return 
            average=np.mean(listofcputimes)
            mintime=min(listofcputimes)
            maxtime=max(listofcputimes)
            ls=[mintime,average,maxtime]
            string='Min, Average, Max time for jobs to finish '+str(ls)+' hours'
            self.WriteToLog(string)

            jobsleft=int(numberofmolstotal)-len(listoffilepaths)
            timeleft=average*jobsleft
            string='There are '+str(jobsleft)+' jobs left, with ETA='+str(timeleft)+',hours for all jobs to finish'
            self.WriteToLog(string)

        def CallJobsLocalHost(self,fulljobtooutputlog,skiperrors):
           for job in fulljobtooutputlog.keys():
               self.call_subsystem([job],False,skiperrors)
           finishedjobs,errorjobs=self.WaitForTermination(fulljobtooutputlog,skiperrors)
           return finishedjobs,errorjobs


        def WaitForTermination(self,jobtooutputlog,skiperrors):
            finishedjobs=[]
            errorjobs=[]
            errormessages=[]
            outputStatusDict = copy.deepcopy(jobtooutputlog)
            for job,outputlog in jobtooutputlog.items():
               finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
               if error==True: # remove log before resubmitting
                   os.remove(outputlog)
            while len(finishedjobs)!=len(jobtooutputlog.keys()):
                for job in jobtooutputlog.keys():
                    outputlog=jobtooutputlog[job]
                    finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
                    if finished==True and error==False: # then check if SP has been submitted or not
                        self.KillJob(job,outputlog)
                        if outputlog not in finishedjobs:
                            self.NormalTerm(outputlog)
                            finishedjobs.append(outputlog)
                    elif finished==False and error==True:
                        self.KillJob(job,outputlog)
                        if skiperrors==False:
                            if outputlog not in finishedjobs:
                                self.ErrorTerm(outputlog,skiperrors)
                                finishedjobs.append(outputlog)
                                errorjobs.append(outputlog)
                        else:
                            finishedjobs.append(outputlog)

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
            self.WriteToLog('All jobs have terminated '+str(finishedjobs))
            return finishedjobs,errorjobs

        def CycleCount(self,logname):
            temp=open(logname,'r')
            results=temp.readlines()
            temp.close()
            count=0
            for line in results:
                if 'Converged' in line or 'Convergence Check' in line:
                    count+=1
            return count


        def CheckNormalTermination(self,logfname,errormessages=None,skiperrors=False): # needs to handle error checking now
            """
            Intent: Checks the *.log file for normal termination
            """
            error=False
            term=False
            lastupdatetofilename={}
            curdir=os.getcwd()
            logpath=os.path.split(logfname)[0]
            if logpath!='':
                os.chdir(logpath)
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
            htime=min(lastupdatetofilename.keys())
            if os.path.isfile(logfname):
                head,tail=os.path.split(logfname)
                if ('opt-' in logfname or 'sp-' in logfname):
                    updatetime=self.lastlogfileupdatetime
                else:
                    updatetime=2.5 # parent QM can take longer dependeing on resources, poltype log can take a while as well to finish QM
                foundendgau=False # sometimes gaussian comments have keyword Error, ERROR in them
                foundhf=False
                for line in open(logfname):
                    if 'poltype' in tail:
                        if 'Poltype Job Finished' in line:
                            term=True
                        elif 'Poltype has crashed!' in line:
                            error=True
                    else:
                        if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line or "LBFGS  --  Normal Termination due to SmallGrad" in line or "Normal termination" in line or 'Normal Termination' in line or 'Total Potential Energy' in line or 'Psi4 stopped on' in line:
                            term=True
                        if ('Tinker is Unable to Continue' in line or 'error' in line or ' Error ' in line or ' ERROR ' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation, address not mapped to object' in line or 'galloc:  could not allocate memory' in line or 'Erroneous write.' in line) and 'DIIS' not in line and 'mpi' not in line and 'RMS Error' not in line:
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
                if self.debugmode==False and foundhf==True:
                    term=False
                        

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


            if errormessages!=None:
                return term,error,errormessages
            else:
                return term,error        


        
            
        def NormalTerm(self,logfname):
            self.WriteToLog("Normal termination: logfile=%s path=%s" % (logfname,os.getcwd()))
        
        
        def ErrorTerm(self,logfname,skiperrors):
            string="ERROR termination: logfile=%s path=%s" % (logfname,os.getcwd())
            self.WriteToLog(string)
            if skiperrors==False:
                raise ValueError(string)

        def call_subsystem(self,cmdstrs,wait=False,skiperrors=False):
            if self.printoutput==True:
                for cmdstr in cmdstrs:
                    print("Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
            procs=[]
            for cmdstr in cmdstrs:
                self.WriteToLog("Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
                p = subprocess.Popen(cmdstr,shell=True,stdout=self.logfh, stderr=self.logfh)
                procs.append(p)
                self.cmdtopid[cmdstr]=p

            if wait==True:
                exit_codes=[p.wait() for p in procs]
                for i in range(len(procs)):
                    exitcode=exit_codes[i]
                    cmdstr=cmdstrs[i] 
                    if skiperrors==False:
                        if exitcode != 0:
                            self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                            raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
                    else:
                        if exitcode != 0:
                            self.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())

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
                if value!=None:
                    line=key+'='+str(value)+'\n'
                    temp.write(line)
            temp.close()
            return inifilepath

        def CheckTorsionParameters(self,keyfilename,torsionsmissing,hydtorsions): # dont check torsions skipped due to rule that if a-b-c-d, a or d is hydroten and not all possible a and d around b-c is hydrogen then torsion is skipped, if database transfers all zeros do not check
            temp=open(keyfilename,'r')
            results=temp.readlines()
            temp.close()
            for lineidx in range(len(results)):
                line=results[lineidx]
                if 'torsion' in line and '#' not in line and 'Missing' not in line and 'none' not in line and 'unit' not in line:
                    allzero=True
                    linesplit=line.split()
                    ls=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
                    revls=ls[::-1]
                    if (ls in torsionsmissing or revls in torsionsmissing) and (ls not in hydtorsions or revls not in hydtorsions):
                        pass
                    else:
                        continue
                    prms=linesplit[5:]
                    newprms=prms[0::3]
                    newprms=[float(i) for i in newprms]
                    for prm in newprms:
                        if prm!=0:
                            allzero=False
                    if allzero==True:
                        self.WriteToLog("torsion parameters are all zero "+line+' path ='+os.getcwd())
                        raise ValueError("torsion parameters are all zero "+line+' path ='+os.getcwd())

                        

        def main(self):
            
             
            if (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
                knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference,modligidxs=modres.GenerateModifiedProteinPoltypeInput(self)
                self.molstructfname=molname
                head, self.molstructfname = os.path.split(self.molstructfname)
                self.molecprefix =  os.path.splitext(self.molstructfname)[0]

            if (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None): # if already have core parameters in modified prm database then dont regenerate parameters
                if check==False:
                    params=self.GenerateParameters()
            else:
               params= self.GenerateParameters()
               return params

        
            if (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
                modres.GenerateModifiedProteinXYZAndKey(self,knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference,modligidxs)
        
        
        def GrabIndexToCoordinates(self,mol):
            indextocoordinates={}
            iteratom = openbabel.OBMolAtomIter(mol)
            for atom in iteratom:
                atomidx=atom.GetIdx()
                rdkitindex=atomidx-1
                coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
                indextocoordinates[rdkitindex]=coords
            return indextocoordinates

        def AddInputCoordinatesAsDefaultConformer(self,m,indextocoordinates):
            conf = m.GetConformer()
            for i in range(m.GetNumAtoms()):
                x,y,z = indextocoordinates[i]
                conf.SetAtomPosition(i,Point3D(x,y,z))

            return m 

                   

        def CheckIfCartesianXYZ(self,f):
            check=True
            temp=open(f,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                linesplit=line.split()
                if len(linesplit)>4:
                    check=False
            return check

        def RemoveCartesianXYZFiles(self):
            files=os.listdir()
            for f in files:
                filename, file_extension = os.path.splitext(f)
                if file_extension=='.xyz' and 'opt' not in filename:
                    check=self.CheckIfCartesianXYZ(f)
                    if check==True:
                        os.remove(f)

               
        
        def CheckInputCharge(self,molecule,verbose=False):
            array=[]
            totchg=0
            atomindextoformalcharge={}
            atomicnumtoformalchg={1:{2:1},5:{4:1},6:{3:-1},7:{2:-1,4:1},8:{1:-1,3:1},15:{4:1},16:{1:-1,3:1,5:-1},17:{0:-1,4:3},9:{0:-1},35:{0:-1},53:{0:-1}}
            for atom in molecule.GetAtoms():
                atomidx=atom.GetIdx()
                atomnum=atom.GetAtomicNum()
                val=atom.GetExplicitValence()
                valtochg=atomicnumtoformalchg[atomnum]
                radicals=atom.GetNumRadicalElectrons()
                if val not in valtochg.keys(): # then assume chg=0
                    chg=0
                else:
                    chg=valtochg[val]
                polneighb=False
                if atomnum==6:
                    for natom in atom.GetNeighbors():
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==7 or natomicnum==8 or natomicnum==16:
                            polneighb=True
                    if polneighb and val==3:
                        chg=1
                string='Atom index = '+str(atomidx+1)+' Atomic Number = ' +str(atomnum)+ ' Valence = '+str(val)+ ' Formal charge = '+str(chg)
                array.append(string)
                if atomnum==6 and val==3 and (self.addhydrogentononcharged==True or self.addhydrogens==True)  and radicals==0:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Carbon, will assume missing hydrogens and add'+string) 
                        self.WriteToLog('WARNING! Strange valence for Carbon, will assume missing hydrogens and add '+string)
                    atom.SetNumRadicalElectrons(0)
                    chg=0
                    atom.SetFormalCharge(chg)


                elif atomnum==7 and val==2 and (self.addhydrogentononcharged==True or self.addhydrogens==True) and radicals==0:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Nitrogen, will assume missing hydrogens and add'+string) 
                        self.WriteToLog('WARNING! Strange valence for Nitrogen, will assume missing hydrogens and add '+string)
                    atom.SetNumRadicalElectrons(0)
                    chg=0
                    atom.SetFormalCharge(chg)

                elif atomnum==8 and val==1 and (self.addhydrogens==True) and radicals==0:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Oxygen, will assume missing hydrogens and add'+string) 
                        self.WriteToLog('WARNING! Strange valence for Oxygen, will assume missing hydrogens and add '+string)
                    atom.SetNumRadicalElectrons(0)
                    chg=0
                    atom.SetFormalCharge(chg)


                elif atomnum==7 and val==2 and radicals==1:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Nitrogen, will assume radical and set charge to zero') 
                        self.WriteToLog('WARNING! Strange valence for Nitrogen, will assume radical and set charge to zero')
                    self.allowradicals=True

                    atom.SetFormalCharge(0)
                    self.addhydrogentononcharged=False

                elif atomnum==8 and val==2 and radicals==1:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1') 
                        self.WriteToLog('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1')
                    self.allowradicals=True

                    atom.SetFormalCharge(1)
                elif atomnum==8 and val==1 and radicals==1:
                    if verbose==True:
                        warnings.warn('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1') 
                        self.WriteToLog('WARNING! Strange valence for Oxygen, will assume radical and set charge to +0')
                    self.allowradicals=True
                    atom.SetFormalCharge(0)
                    self.addhydrogentononcharged=False

                else:
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

        

        

        def GenerateExtendedConformer(self,rdkitmol,mol):
            numconf=100 # just try this
            AllChem.EmbedMultipleConfs(rdkitmol, numConfs=numconf,useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
            energies = AllChem.MMFFOptimizeMoleculeConfs(rdkitmol,maxIters=2000, nonBondedThresh=100.0)
            confs=rdkitmol.GetConformers()
            listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=databaseparser.GrabAtomsForParameters(self,mol)
            disttoconf={}
            energytoconf={}
            for i in range(len(confs)):
                conf=confs[i]
                energy=energies[i][1]
                name="conftest.mol"
                rdmolfiles.MolToMolFile(rdkitmol,name,confId=i)
                mol=rdmolfiles.MolFromMolFile(name,removeHs=False)
                maxdist=self.FindLongestDistanceInMolecule(mol)
                disttoconf[maxdist]=i
                energytoconf[energy]=i
            distances=list(disttoconf.keys())
            if len(distances)!=0:
                maxdist=max(distances)
                confindex=disttoconf[maxdist]
            else:
                confindex=0
            minenergy=min(energytoconf.keys())
            newenergytoconf={}
            for e,conf in energytoconf.items():
                newe=e-minenergy
                newenergytoconf[newe]=conf 
            sortedenergy=sorted(newenergytoconf)
            sortedconf=[newenergytoconf[i] for i in sortedenergy]
            sortedenergytoconf=dict(zip(sortedenergy,sortedconf))
            
            confslist=[confindex]
            numextraconfs=self.numespconfs-1
            count=1
            for energy,idx in sortedenergytoconf.items():
                if idx not in confslist:
                    if count>numextraconfs:
                        break
                    confslist.append(idx)
                    count+=1
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
            indextocoordinates={}
            obConversion = openbabel.OBConversion()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(self.molstructfname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, self.molstructfname)
            iteratombab = openbabel.OBMolAtomIter(mol)
            for atm in iteratombab:
                atmindex=atm.GetIdx()
                coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
                indextocoordinates[atmindex-1]=coords

            return indextocoordinates


        def FindLongestDistanceInMolecule(self,mol):
            veclist=[]
            for i in range(len(mol.GetAtoms())):
                pos = mol.GetConformer().GetAtomPosition(i) 
                vec=np.array([float(pos.x),float(pos.y),float(pos.z)])
                veclist.append(vec)
            pairs=list(itertools.combinations(veclist, 2))
            distlist=[]
            for pairidx in range(len(pairs)):
                pair=pairs[pairidx]
                dist=np.linalg.norm(np.array(pair[0])-np.array(pair[1]))
                distlist.append(dist)
            maxdist=np.amax(np.array(distlist))
            return maxdist
        
        def SetDefaultCoordinatesBabel(self,mol,indextocoordinates):
            for index,coords in indextocoordinates.items():
                atom=mol.GetAtom(index+1)
                x=coords[0]
                y=coords[1]
                z=coords[2]
                atom.SetVector(x,y,z)
            return mol

        def CheckForConcentratedFormalCharges(self,m,atomindextoformalcharge):
            chargedindices=[]
            pcm=False
            for atomindex,chg in atomindextoformalcharge.items():
                if chg!=0:
                    chargedindices.append(atomindex)
            for atomindex in chargedindices:
                atom=m.GetAtomWithIdx(atomindex)
                for atm in atom.GetNeighbors():
                    atmidx=atm.GetIdx()
                    if atmidx in chargedindices:
                        pcm=True
                    for natm in atm.GetNeighbors():
                        natmidx=natm.GetIdx()
                        if natmidx!=atomindex:
                            if natmidx in chargedindices:
                                pcm=True
                            for nnatm in natm.GetNeighbors():
                                nnatmidx=nnatm.GetIdx()
                                if nnatmidx!=atmidx:
                                    if nnatmidx in chargedindices:
                                        pcm=True
            return pcm 
                    
                        
        def DeleteAllNonQMFiles(self,folderpath=None):
            tempdir=os.getcwd()
            if folderpath!=None:
                os.chdir(folderpath)
            self.DeleteNonQMFiles(os.getcwd()) 
            if os.path.isdir('qm-torsion'):
                os.chdir('qm-torsion')
                self.DeleteNonQMFiles(os.getcwd()) 
                os.chdir('..')
            if os.path.isdir('vdw'):
                os.chdir('vdw')
                self.DeleteNonQMFiles(os.getcwd()) 
                os.chdir('..')

            os.chdir(tempdir)
           
        def DeleteNonQMFiles(self,directory):
            tempdir=os.getcwd()
            os.chdir(directory)
            deletearray=[]
            files=os.listdir()
            filestonotdelete=[]
            if self.indextotypefile!=None:
                filestonotdelete.append(self.indextotypefile)
            if self.inputkeyfile!=None:
                filestonotdelete.append(self.inputkeyfile)
            for f in files:
                if not os.path.isdir(f) and 'nohup' not in f and f[0]!='.' and f!='parentvdw.key':
                    fsplit=f.split('.')
                    if len(fsplit)>1:
                        end=fsplit[1]
                        if 'log' not in end and 'sdf' not in end and 'ini' not in end and 'chk' not in end and 'dat' not in end and 'mol' not in end and 'txt' not in end and f not in filestonotdelete:
                            deletearray.append(f)
            for f in deletearray:
                os.remove(f)

            os.chdir(tempdir) 

        def ResourceInputs(self,jobpaths,ramalljobs,diskalljobs,numprocalljobs):
            scratchspacelist=[]
            ramlist=[]
            numproclist=[]
            for i in range(len(jobpaths)):
                scratchspacelist.append(diskalljobs)
                ramlist.append(ramalljobs)
                numproclist.append(numprocalljobs)
        
            return scratchspacelist,ramlist,numproclist

        def GenerateDaemonInput(self,joblist,outputlogpath,scratchspacelist,ramlist,numproclist,jobpaths,inputdaemonfilepath):
            head,tail=os.path.split(outputlogpath)
            os.chdir(head)
            temp=open(inputdaemonfilepath,'w')
            for i in range(len(joblist)):
                job=joblist[i]
                scratchspace=scratchspacelist[i]
                ram=ramlist[i]
                numproc=numproclist[i]
                jobpath=jobpaths[i]
                string='--job='+job+' '+'--numproc='+str(1)+' '+'--ram=10GB'+' '+'--disk=0GB'+' '+'--inputfilepaths='+os.path.join(jobpath,'poltype.ini')+' '+'\n'
                temp.write(string)
            temp.close()


        def CreatePoltypeInputFilesMultipleMolecules(self):
            dic={'externalapi':self.externalapi,'accuratevdwsp':self.accuratevdwsp,'email':self.email,'firstoptfinished':self.firstoptfinished,'optonly':self.optonly,'onlyvdwatomindex':self.onlyvdwatomindex,'use_qmopt_vdw':self.use_qmopt_vdw,'use_gau_vdw':self.use_gau_vdw,'dontusepcm':self.dontusepcm,'deleteallnonqmfiles':self.deleteallnonqmfiles,'totalcharge':self.totalcharge,'torspbasissethalogen':self.torspbasissethalogen,'homodimers':self.homodimers,'boltzmantemp':self.boltzmantemp,'dovdwscan':self.dovdwscan,'use_gausgeomoptonly':self.use_gausgeomoptonly,'maxtorRMSPDRel':self.maxtorRMSPDRel,'tortor':self.tortor,'fitfirsttorsionfoldphase':self.fitfirsttorsionfoldphase,'defaultmaxtorsiongridpoints':self.defaultmaxtorsiongridpoints,'absdipoletol':self.absdipoletol,'refinenonaroringtors':self.refinenonaroringtors,'maxgrowthcycles':self.maxgrowthcycles,'use_gauPCM':self.use_gauPCM,'fitqmdipole':self.fitqmdipole,'WBOtol':self.WBOtol,'dontfrag':self.dontfrag,'dipoletol':self.dipoletol,'numproc':self.numproc,'maxmem':self.maxmem,'maxdisk':self.maxdisk,'optbasisset':self.optbasisset,'toroptbasisset':self.toroptbasisset,'dmabasisset':self.dmabasisset,'espbasisset':self.espbasisset,'torspbasisset':self.torspbasisset,'optmethod':self.optmethod,'toroptmethod':self.toroptmethod,'torspmethod':self.torspmethod,'dmamethod':self.dmamethod,'espmethod':self.espmethod,'qmonly' : self.qmonly,'espfit' : self.espfit,'foldnum':self.foldnum,'maxRMSD':self.maxRMSD,'maxRMSPD':self.maxRMSPD,'maxtorRMSPD':self.maxtorRMSPD,'tordatapointsnum':self.tordatapointsnum,'torsionrestraint':self.torsionrestraint,'rotalltors':self.rotalltors,'dontdotor':self.dontdotor,'dontdotorfit':self.dontdotorfit,'toroptpcm':self.toroptpcm,'optpcm':self.optpcm,'torsppcm':self.torsppcm,'use_gaus':self.use_gaus,'use_gausoptonly':self.use_gausoptonly,'freq':self.freq,'optmaxcycle':self.optmaxcycle,'forcefield':self.forcefield}
            os.chdir(self.inputmoleculefolderpaths)
            files=os.listdir()
            jobpaths=[]
            joblist=[]
            for f in files:
                if os.path.isdir(f):
                    os.chdir(f)
                    subfiles=os.listdir()
                    for subf in subfiles:
                        if '.sdf' in subf:
                            dic['structure']=subf
                            inifilepath=self.WritePoltypeInitializationFile(dic)
                            curdir=os.getcwd()
                            jobpaths.append(curdir)
                            poltypefilepath=os.path.join(self.poltypepath,'poltype.py')
                            joblist.append('cd '+curdir+' '+'&&'+' '+'python '+poltypefilepath)

                    os.chdir('..')
            if self.externalapi!=None:
                outputlogfilepath=os.path.join(self.inputmoleculefolderpaths,'outputlog.txt')
                inputdaemonfilepath=os.path.join(self.inputmoleculefolderpaths,'jobinfo.txt')
                scratchspacelist,ramlist,numproclist=self.ResourceInputs(jobpaths,self.maxmem,self.maxdisk,self.numproc)
                self.GenerateDaemonInput(joblist,outputlogfilepath,scratchspacelist,ramlist,numproclist,jobpaths,inputdaemonfilepath)

            sys.exit()

        def CheckIfAtomsAreAllowed(self,m):
            listofallowedatoms=[1,6,7,8,15,16,17,35,53,9]
            for atom in m.GetAtoms():
                atomicnum=atom.GetAtomicNum()
                if atomicnum not in listofallowedatoms:
                    raise ValueError('Element not allowed! '+str(atomicnum)) 


        def CheckMP2OptFailed(self):
            files=os.listdir()
            mp2failed=False
            found=False
            for f in files:
                if 'opt' in f and ('.com' in f or '.psi4' in f) and '_temp' not in f:
                    found=True
                    break
            if found==True:
                temp=open(f,'r')
                results=temp.readlines()
                temp.close()
                foundHF=False
                foundminix=False
                for line in results:
                    if 'optimize' in line or 'opt' in line:
                        if 'hf' in line or 'HF' in line:
                            foundHF=True
                        if 'minix' in line or 'MINIX' in line:
                            foundminix=True 
                if foundHF==True and foundminix==False:
                    mp2failed=True
            return mp2failed


        def GrabAtomicSymbols(self,rdkitmol):
            indextoatomicsymbol={}
            an = pyasl.AtomicNo()
            for atm in rdkitmol.GetAtoms():
                atmnum=atm.GetAtomicNum()
                atmidx=atm.GetIdx()
                sym=an.getElSymbol(atmnum)
                indextoatomicsymbol[atmidx+1]=sym

            return indextoatomicsymbol


        def CheckIfInputIsTinkerXYZ(self,molstructfname):
            istinkxyz=False
            temp=open(molstructfname,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                linesplit=line.split()
                if len(linesplit)>1:
                    if len(linesplit)>4:
                        istinkxyz=True

            return istinkxyz


        def GenerateIndexToTypeFile(self,indextotype):
            filename='indextotype.txt'
            temp=open(filename,'w')
            for index,typenum in indextotype.items():
                temp.write(str(index)+' '+str(typenum)+'\n')
            temp.close()
            return filename


        def GrabIndexToType(self,molstrucfname):
            indextotype={}
            temp=open(molstrucfname,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                linesplit=line.split()
                if len(linesplit)>=5:
                    index=int(linesplit[0])
                    typenum=int(linesplit[5])
                    indextotype[index]=typenum



            return indextotype



        def ConvertInputStructureToSDFFormat(self,molstructfname):
            obConversion = openbabel.OBConversion()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(molstructfname)
            split=molstructfname.split('.')
            ext=split[-1]
            if ext!='sdf':
                istinkerxyz=False
                if ext=='xyz':
                    istinkerxyz=self.CheckIfInputIsTinkerXYZ(molstructfname)

                obConversion.SetInFormat(ext)
                if istinkerxyz==False:
                    obConversion.ReadFile(mol, molstructfname)
                else:
                    newname=self.ConvertTinkerXYZToCartesianXYZ(molstructfname)
                    obConversion.ReadFile(mol, newname)
                    indextotype=self.GrabIndexToType(molstructfname)
                    filename=self.GenerateIndexToTypeFile(indextotype) 
                    self.indextotypefile=filename

                obConversion.SetOutFormat('sdf')
                molstructfname=molstructfname.replace('.'+ext,'.sdf')
                obConversion.WriteFile(mol,molstructfname)


            return molstructfname


        def GenerateEntireTinkerXYZ(self,atmindextocoordinates,atmindextotypenum,atmindextoconnectivity,atmindextoelement,filename):
            temp=open(filename,'w')
            temp.write(str(len(atmindextocoordinates.keys()))+'\n')
            for index,coords in atmindextocoordinates.items():
                typenum=atmindextotypenum[index]
                conns=atmindextoconnectivity[index]
                element=atmindextoelement[index]
                x=coords[0]
                y=coords[1]
                z=coords[2]
                newline='    '+str(index)+'  '+element+'     '+str(x)+'   '+str(y)+'   '+str(z)+'    '+str(typenum)+'     '
                for con in conns:
                    newline+=str(con)+'     '
                newline+='\n'
                temp.write(newline)

            temp.close()


        def GenerateParameters(self):
            
            temp=open(os.getcwd()+r'/'+'poltype.ini','r')
            results=temp.readlines()
            temp.close()
            for line in results:
                if '#' not in line and line!='\n':
                    if '=' in line:
                        linesplit=line.split('=',1)
                        a=linesplit[1].replace('\n','').rstrip().lstrip()
                        newline=linesplit[0]
                        if a=='None':
                            continue
                    else:
                        newline=line

                    if "structure" in newline:
                        self.molstructfname = a
            self.molstructfname=self.ConvertInputStructureToSDFFormat(self.molstructfname)
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
            if self.deleteallnonqmfiles==True:
                self.DeleteAllNonQMFiles()
            if self.inputmoleculefolderpaths!=None:
                self.CreatePoltypeInputFilesMultipleMolecules() 
            if self.optmaxcycle>=100:
                self.fullopt=True
            else:
                self.fullopt=False
            obConversion = openbabel.OBConversion()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(self.molstructfname)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, self.molstructfname)
            self.atomnum=mol.NumAtoms() 
            self.logfh = open(self.logfname,"w",buffering=1)

            self.molstructfnamemol=self.molstructfname.replace('.sdf','.mol')
            if '.mol' not in self.molstructfname: 
                obConversion.SetOutFormat('mol')
                obConversion.WriteFile(mol,self.molstructfnamemol)
            indextocoordinates=self.GrabIndexToCoordinates(mol)
            m=Chem.MolFromMolFile(self.molstructfnamemol,removeHs=False,sanitize=False)
            self.CheckIfAtomsAreAllowed(m)
            m,atomindextoformalcharge=self.CheckInputCharge(m,verbose=True)
            if self.allowradicals==True:
                self.dontfrag=True # Psi4 doesnt allow UHF and properties (like compute WBO) for fragmenter, so need to turn of fragmenter if radical detected
            m.UpdatePropertyCache()
            if self.addhydrogentononcharged==True and self.isfragjob==False:
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m)

            self.indextoatomicsymbol=self.GrabAtomicSymbols(m)
            Chem.SanitizeMol(m)
            smarts=rdmolfiles.MolToSmarts(m)
            if '.' in smarts:
                raise ValueError('Multiple fragments detectected in input molecule')
            pcm=self.CheckForConcentratedFormalCharges(m,atomindextoformalcharge)
            cpm = copy.deepcopy(m)
            indextocoordslist=[indextocoordinates]
            iswholemoleculering=self.CheckIfWholeMoleculeIsARing(mol)
            if iswholemoleculering==True:
                self.generateextendedconf=False
            if self.firstoptfinished==False and self.isfragjob==False and self.generateextendedconf==True:
                indextocoordslist=self.GenerateExtendedConformer(m,mol)
                indextocoordinates=indextocoordslist[0]
            Chem.GetSymmSSSR(m)
            m.GetRingInfo().NumRings() 
            try:
                m=self.AddInputCoordinatesAsDefaultConformer(m,indextocoordinates)
            except:
                pass
            if self.generateextendedconf==True:
                rdmolfiles.MolToMolFile(m,'extendedconf.mol')
            mol,m=self.CheckIsInput2D(mol,obConversion,m)
            if not os.path.exists(self.scrtmpdirpsi4):
                os.mkdir(self.scrtmpdirpsi4)
            if not os.path.exists(self.scrtmpdirgau):
                os.mkdir(self.scrtmpdirgau)

            mol=self.SetDefaultCoordinatesBabel(mol,indextocoordinates)
            molist=self.GenerateListOfMols(mol,indextocoordslist)
            self.mol=mol

            self.rdkitmol=m
            self.mol.SetTotalCharge(self.totalcharge)
            if self.keyfiletoaddtodatabase!=None:
                databaseparser.AddKeyFileParametersToParameterFile(self,m)   
                sys.exit()
            self.GrabIonizationStates(m)
            self.GrabTautomers(m)
            if self.genprotstatesonly==True:
                sys.exit()
            if ('Br ' in self.mol.GetSpacedFormula()):
                self.torspbasisset=self.torspbasissethalogen
            self.pcm=False
            if pcm==True and self.dontusepcm==False:
                self.pcm=True
                self.toroptpcm=True
                self.optpcm=True
                self.torsppcm=True

            
            if self.use_gauPCM==True:
                self.use_gausoptonly=False
                self.use_gaus=True
                self.SanitizeAllQMMethods()

            atomiter=openbabel.OBMolAtomIter(self.mol)
            atomnum=0
            for atom in atomiter:
                atomnum+=1
                atomidx=atom.GetIdx()
                atomicnum=atom.GetAtomicNum()


            self.RemoveCartesianXYZFiles()
 
            self.WriteToLog("Running on host: " + gethostname())
            # Initializing arrays
            
            self.canonicallabel = [ 0 ] * mol.NumAtoms()
            self.localframe1 = [ 0 ] * mol.NumAtoms()
            self.localframe2 = [ 0 ] * mol.NumAtoms()
            self.WriteToLog("Atom Type Classification")

            self.idxtosymclass,self.symmetryclass=symm.gen_canonicallabels(self,mol,None,self.usesymtypes,True) 
 
            # QM calculations are done here
            # First the molecule is optimized. (-opt) 
            # This optimized molecule is stored in the structure optmol
            # Then the electron density matrix is found (-dma)
            # This is used by GDMA to find multipoles
            # Then information for generating the electrostatic potential grid is found (-esp)
            # This information is used by cubegen
            self.comoptfname=self.firstcomoptfname 
            self.chkoptfname=self.firstchkoptfname 
            self.fckoptfname=self.firstfckoptfname
            self.logoptfname=self.firstlogoptfname 
            self.gausoptfname=self.firstgausoptfname 

            if self.use_gausgeomoptonly==True:
                self.use_gausoptonly=True

            torgen.FindPartialDoubleBonds(self,m)
                
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
                optmol=mol
                tmpconv = openbabel.OBConversion()
                tmpconv.SetOutFormat('xyz')
                tmpconv.WriteFile(mol, self.logoptfname.replace('.log','.xyz'))
                atmindextocoordinates,atmindextoconnectivity,atmindextoelement=self.ExtractMOLInfo(mol)
                self.GenerateEntireTinkerXYZ(atmindextocoordinates,self.idxtosymclass,atmindextoconnectivity,atmindextoelement,self.xyzfname)

            if self.optonly==True:
                sys.exit()
            if self.use_gausgeomoptonly==True:
                self.use_gausoptonly=False
                self.use_gaus=False
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
            #try:
            if (self.writeoutpolarize==True and self.writeoutmultipole==True):
                esp.SPForDMA(self,optmol,mol)
            #except:
            #    if self.use_gaus==False: 
            #        self.use_gaus=True
            #        esp.SPForDMA(self,optmol,mol) 
            #        self.use_gaus=False
            #    else:
            #        traceback.print_exc(file=sys.stdout)
            #        sys.exit()

            # Obtain multipoles from Gaussian fchk file using GDMA
        
                if not os.path.isfile(self.gdmafname):
                    mpole.run_gdma(self)
        

                polarindextopolarizeprm,polartypetotransferinfo=databaseparser.GrabSmallMoleculeAMOEBAParameters(self,optmol,mol,m,polarize=True)
                mpole.gen_peditinfile(self,mol,polarindextopolarizeprm)
            
                if (not os.path.isfile(self.xyzfname) or not os.path.isfile(self.keyfname)):
                    # Run poledit
                    cmdstr = self.peditexe + " 1 " + self.gdmafname +' '+self.paramhead+ " < " + self.peditinfile
                    self.call_subsystem([cmdstr],True)
                    # Add header to the key file output by poledit
                    while not os.path.isfile(self.keyfnamefrompoledit):
                        time.sleep(1)
                        self.WriteToLog('Waiting for '+self.keyfnamefrompoledit)
                        
                    mpole.prepend_keyfile(self,self.keyfnamefrompoledit,optmol)
                    mpole.SanitizeMultipoleFrames(self,self.keyfnamefrompoledit)
                    shutil.copy(self.keyfnamefrompoledit,self.keyfname)

                xyzfnamelist,keyfnamelist=self.GenerateDuplicateXYZsFromOPTs(self.xyzfname,self.keyfname,optmolist)
            # post process local frames written out by poledit
            if self.atomnum!=1 and (self.writeoutpolarize==True and self.writeoutmultipole==True): 
                 #try:
                 gridnamelist,espnamelist,fchknamelist,cubenamelist=esp.SPForESP(self,optmolist,molist,xyzfnamelist,keyfnamelist) 
                 #except:
                 #    if self.use_gaus==False: 
                 #        self.use_gaus=True
                 #        gridnamelist,espnamelist,fchknamelist,cubenamelist=esp.SPForESP(self,optmolist,molist,xyzfnamelist,keyfnamelist) 
                 #        self.use_gaus=False
                 #    else:
                 #        traceback.print_exc(file=sys.stdout)
                 #        sys.exit()

            # End here if qm calculations were all that needed to be done 
            if self.qmonly:
                self.WriteToLog("poltype QM-only complete.")
                sys.exit(0)
        
                   
            
            # generate the electrostatic potential grid used for multipole fitting
            if (self.writeoutpolarize==True and self.writeoutmultipole==True):
                if self.atomnum!=1: 
                    if not os.path.isfile(self.key3fname):
                        potnamelist=esp.gen_esp_grid(self,optmol,gridnamelist,espnamelist,fchknamelist,cubenamelist)
        
                # Average multipoles based on molecular symmetry
                # Does this using the script avgmpoles.pl which is found in the poltype directory
                # Atoms that belong to the same symm class will now have only one common multipole definition
                if not os.path.isfile(self.key2fnamefromavg):
                    self.WriteToLog("Average Multipoles Via Symmetry")
                    mpole.AverageMultipoles(self,optmol)
                    mpole.AddPolarizeCommentsToKey(self,self.key2fnamefromavg,polartypetotransferinfo)
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
                mpole.rm_esp_terms_keyfile(self,self.key3fname)
                if fit==True:
                    if self.atomnum!=1: 
                        esp.ElectrostaticPotentialComparison(self,combinedxyz,combinedpot) 
            
            if not os.path.exists(self.key4fname):
                databaseparser.appendtofile(self,self.key3fname,self.key4fname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo)
                if self.writeoutangle==True:
                    databaseparser.StiffenZThenBisectorAngleConstants(self,self.key4fname)
                    databaseparser.TestBondAngleEquilValues(self)
                self.AddIndicesToKey(self.key4fname)
                if self.databasematchonly==True:
                    sys.exit()

            if self.torsppcm:
                torgen.PrependStringToKeyfile(self,self.key4fname,'solvate GK')
            torgen.get_all_torsions(self,mol)
            # Find rotatable bonds for future torsion scans
            (torlist, self.rotbndlist,hydtorsions,nonaroringtorlist) = torgen.get_torlist(self,mol,torsionsmissing)
            if atomnum<25 and len(nonaroringtorlist)==0 and self.smallmoleculefragmenter==False: 
                self.dontfrag=True
            if self.dontfrag==True: # if fragmenter is turned off, parition resources by jobs at sametime for parent
                self.maxmem,self.maxdisk,self.numproc=self.PartitionResources()
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
            
            # add missingvdwindices to torlist (for fragmenter input)
            missingvdwatomsets=[]
            if self.isfragjob==False and self.dovdwscan==True:
                for vdwatomindex in missingvdwatomindices:
                    ls=tuple([tuple([vdwatomindex])])
                    missingvdwatomsets.append(ls)
                    self.torlist.append(ls)
  
            self.torsettofilenametorset={}
            self.torsettotortorindex={}
            self.torsettotortorphaseindicestokeep={}
            self.nonaroringtors=[]
            self.nonaroringtorsets=[]
            self.classkeytoinitialprmguess={}
            self.nonarotortotorsbeingfit={}
            


            if self.isfragjob==False and not os.path.isfile(self.key5fname) and self.dontfrag==False and (self.dovdwscan==True):
                self.WriteToLog('Create and Run vdW Fragment Poltype Jobs')

                WBOmatrix,outputname,error=frag.GenerateWBOMatrix(self,self.rdkitmol,self.mol,self.logoptfname.replace('.log','.xyz'))
                rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor=frag.GenerateFragments(self,self.mol,self.torlist,WBOmatrix,missingvdwatomsets,nonaroringtorlist) # returns list of bond indexes that need parent molecule to do torsion scan for (fragment generated was same as the parent0
                equivalentrotbndindexarrays,rotbndindextoringtor,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray=frag.SpawnPoltypeJobsForFragments(self,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor)
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
            shutil.copy(self.key5fname,self.key6fname)
            self.torlist=torlist[:]
            if self.tortor==True and self.dontdotor==False:
                torgen.PrepareTorsionTorsion(self,optmol,mol,tortorsmissing)
            torgen.DefaultMaxRange(self,self.torlist)
            if self.refinenonaroringtors==True and self.dontfrag==False:
                rings.RefineNonAromaticRingTorsions(self,mol,optmol,classkeytotorsionparametersguess)
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
            if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key7fname) and (self.dontdotor==False) and len(self.torlist)!=0:
                frag.GrabVdwAndTorsionParametersFromFragments(self,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor,self.key6fname,self.key7fname,rotbndindextoparentrotbndindexes,rotbndindextosmartsindexarray) # just dump to key_6 since does not exist for parent molecule
            else:
                # Torsion scanning then fitting. *.key_7 will contain updated torsions
                if not os.path.isfile(self.key7fname):
                    if len(self.torlist)!=0:
                        # torsion scanning
                        torgen.gen_torsion(self,optmol,self.torsionrestraint,mol)
                        # torsion fitting
                        if self.dontdotorfit==True:
                            shutil.copy(self.key6fname,self.key7fname)
                            sys.exit()
                        torfit.process_rot_bond_tors(self,optmol)
                    else:
                        shutil.copy(self.key6fname,self.key7fname)           
           


 
            if self.isfragjob==False and self.dontdotor==False:
                self.CheckTorsionParameters(self.key7fname,torsionsmissing,hydtorsions)
            self.WriteOutLiteratureReferences(self.key7fname) 
            # A series of tests are done so you one can see whether or not the parameterization values
            # found are acceptable and to what degree
            if self.writeoutpolarize==False or self.writeoutmultipole==False:
                shutil.copy(self.xyzfname,self.xyzoutfile)
            shutil.copy(self.xyzoutfile,self.tmpxyzfile)
            shutil.copy(self.key7fname,self.tmpkeyfile)

            if self.writeoutpolarize and self.writeoutmultipole==True:
                opt.StructureMinimization(self,torsionrestraints)
                if self.atomnum != 1:
                    opt.gen_superposeinfile(self)
                    opt.CheckRMSD(self)

            if self.torsppcm:
                torgen.RemoveStringFromKeyfile(self,self.key7fname,'solvate GK')
            if self.atomnum!=1: 
                 esp.CheckDipoleMoments(self,optmol)
            if os.path.exists(self.tmpkeyfile):
                self.FinalVDWMultipoleCheck(self.tmpkeyfile)
            for torset,totaltime in self.torsettototalqmtime.items():
                sptime=self.torsettospqmtime[torset]
                opttime=self.torsettooptqmtime[torset]
                numpoints=self.torsettonumpoints[torset]
                self.WriteToLog('Total torsion QM time for '+str(torset)+' is '+str(round(totaltime,3))+' hours'+' and '+str(numpoints)+' conformations')
                self.WriteToLog('OPT torsion QM time for '+str(torset)+' is '+str(round(opttime,3))+' hours'+' and '+str(numpoints)+' conformations')
                self.WriteToLog('SP torsion QM time for '+str(torset)+' is '+str(round(sptime,3))+' hours'+' and '+str(numpoints)+' conformations')
            self.WriteToLog('Poltype Job Finished'+'\n')
            
            if os.path.exists(self.scrtmpdirgau):
                shutil.rmtree(self.scrtmpdirgau)
            if os.path.exists(self.scrtmpdirpsi4):
                shutil.rmtree(self.scrtmpdirpsi4)
            if self.email!=None:
                moleculename=self.molstructfname.replace('.sdf','')
                password='amoebaisbest'
                fromaddr = 'poltypecrashreportnoreply@gmail.com'
                toaddr = self.email
                TEXT='Molecule has finished parameterization'
                try: 
                    self.SendFinalReportEmail(TEXT,fromaddr,toaddr,password,moleculename)
                except:
                    pass
            if self.isfragjob==False:
                previousdir=os.path.abspath(os.path.join(os.getcwd(), os.pardir))
                if os.path.exists(self.tmpxyzfile):
                     shutil.copy(self.tmpxyzfile,os.path.join(previousdir,self.tmpxyzfile))
                if os.path.exists(self.tmpkeyfile):
                     shutil.copy(self.tmpkeyfile,os.path.join(previousdir,self.tmpkeyfile))

            self.CopyFitPlots()
            if (self.binding==True or self.solvation==True or self.neatliquidsim==True) or self.usepdb2pqr==True or self.pdbcode!=None:
                self.ligandxyzfilenamelist=[self.tmpxyzfile]
                self.keyfilename=self.tmpkeyfile
                newfolder='Sim'
                if not os.path.exists(newfolder):
                    os.mkdir(newfolder)
                os.chdir(newfolder)
                for xyz in self.ligandxyzfilenamelist:
                    shutil.copy(xyz,os.path.join(newfolder,xyz))
                for key in self.keyfilenamelist:
                    shutil.copy(key,os.path.join(newfolder,key))
                self.MolecularDynamics()



        def GenerateDuplicateXYZsFromOPTs(self,xyzfname,keyfname,optmolist):
            xyzfnamelist=[]
            keyfnamelist=[]
            for optmolidx in range(len(optmolist)):
                optmol=optmolist[optmolidx]
                indextocoordinates=self.GrabIndexToCoordinates(optmol)
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
            atmindextocoordinates={}
            atmindextoconnectivity={}
            atmindextoelement={}
            iteratom = openbabel.OBMolAtomIter(mol)
            an = pyasl.AtomicNo()
            for atom in iteratom:
                index=atom.GetIdx()
                atomicnum=atom.GetAtomicNum()
                coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
                element=an.getElSymbol(atomicnum)
                atmindextocoordinates[index]=coords
                connectivity=self.GrabConnectivity(mol,index)
                atmindextoconnectivity[index]=connectivity
                atmindextoelement[index]=element


            return atmindextocoordinates,atmindextoconnectivity,atmindextoelement


        def GrabConnectivity(self,mol,index):
            conn=[]
            bonditer=openbabel.OBMolBondIter(mol)
            for bond in bonditer:
                oendidx = bond.GetEndAtomIdx()
                obgnidx = bond.GetBeginAtomIdx()
                if oendidx==index:
                    if obgnidx not in conn:
                        conn.append(obgnidx)
                elif obgnidx==index:
                    if oendidx not in conn:
                        conn.append(oendidx)
            return conn


        def GenerateTinkerXYZ(self,xyzfname,newxyzfname,indextocoordinates):
            temp=open(xyzfname,'r')
            results=temp.readlines()
            temp.close()
            temp=open(newxyzfname,'w')
            for line in results:
                linesplit=line.split()
                if len(linesplit)>1:
                    index=int(linesplit[0])
                    coords=indextocoordinates[index-1]
                    linesplit[2]=str(coords[0])
                    linesplit[3]=str(coords[1]) 
                    linesplit[4]=str(coords[2])
                    line=' '.join(linesplit)+'\n'
                temp.write(line)
            temp.close()


        def GenerateListOfMols(self,mol,indextocoordslist):
            molist=[mol]
            obConversion = openbabel.OBConversion()
            for i in range(len(indextocoordslist)):
                if i!=0:
                    indextocoordinates=indextocoordslist[i] 
                    othermol = openbabel.OBMol()
                    inFormat = obConversion.FormatFromExt(self.molstructfname)
                    obConversion.SetInFormat(inFormat)
                    obConversion.ReadFile(othermol, self.molstructfname)
                    othermol=self.SetDefaultCoordinatesBabel(othermol,indextocoordinates)
                    othermol.SetTotalCharge(self.totalcharge)
                    molist.append(othermol)

            return molist


        def FinalVDWMultipoleCheck(self,keyfile):
            temp=open(keyfile,'r')
            results=temp.readlines()
            temp.close()
            vdwtypetofound={}
            mpoletypetofound={}
            typenums=list(self.idxtosymclass.values())
            for typenum in typenums:
                typenum=str(typenum)
                vdwtypetofound[typenum]=False
                mpoletypetofound[typenum]=False
            for line in results:
                if '#' not in line:
                    if 'vdw' in line or 'multipole' in line:
                        linesplit=line.split()
                        typenum=linesplit[1]
                        if 'vdw' in line:
                            vdwtypetofound[typenum]=True
                        elif 'multipole' in line:
                            mpoletypetofound[typenum]=True

            missingvdw=[]
            missingmpole=[]
            for typenum,found in vdwtypetofound.items():
                if found==False:
                    missingvdw.append(typenum)

            for typenum,found in mpoletypetofound.items():
                if found==False:
                    missingmpole.append(typenum)

            if len(missingvdw)!=0:
                if self.firsterror==True:
                    raise ValueError('Missing vdw parameters '+str(missingvdw))
                else:
                    self.DeleteFilesWithExtension(['key','xyz','key_2','key_3','key_4','key_5','xyz_2'])
                    self.firsterror=True
                    self.GenerateParameters()


            if len(missingmpole)!=0:
                if self.firsterror==True:
                    raise ValueError('Missing multipole parameters '+str(missingmpole))
                else:
                    self.DeleteFilesWithExtension(['key','xyz','key_2','key_3','key_4','key_5','xyz_2'])
                    self.firsterror=True
                    self.GenerateParameters()


        def AddIndicesToKey(self,keyfilename):
            temp=open(keyfilename,'r')
            results=temp.readlines()
            temp.close()
            tempname=keyfilename.replace('.key','_TEMP.key')
            temp=open(tempname,'w')
           
            for line in results:
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
                        allindices=[]
                        linesplit=line.split()
                        typenums=[linesplit[i] for i in indices]
                        typenums=[int(i) for i in typenums]
                        for typenum in typenums: 
                            indexes=self.GrabKeysFromValue(self.idxtosymclass,typenum)
                            allindices.append(indexes)
                        string='# '+str(typenums) +' = '+str(allindices)+'\n'
                        temp.write(string)
                temp.write(line)


            temp.close()
            os.remove(keyfilename)
            os.replace(tempname,keyfilename)


        def GrabKeysFromValue(self,dic,thevalue):
            keylist=[]
            for key,value in dic.items():
                if value==thevalue:
                    keylist.append(key)
            return keylist


        def DeleteFilesWithExtension(self,ls):
            files=os.listdir()
            for f in files:
                if '.' in f and not os.path.isdir(f):
                    ext=f.split('.')[1]
                    if ext in ls:
                        os.remove(f)

        def DeleteFilesWithString(self,ls):
            files=os.listdir()
            for f in files:
                for string in ls:
                    if string in f:
                        os.remove(f)


        def DeleteAllFiles(self,deletefragments=False):
            files=os.listdir()
            for f in files:
                if os.path.isdir(f):
                    if deletefragments==True:
                        shutil.rmtree(f)
                    else:
                        if f!='Fragments':
                            shutil.rmtree(f)
                elif '.' in f and not os.path.isdir(f) and f!='parentvdw.key':
                    ext=f.split('.')[1]
                    if ext!='txt' and ext!='mol' and ext!='sdf' and ext!='ini' and 'nohup' not in f and f[0]!='.':
                        os.remove(f)

        def CheckFileSanity(self):
            optxyzname=self.logoptfname.replace('.log','.xyz')
            temp=open(optxyzname,'r')
            results=temp.readlines()
            temp.close()
            atomindextocoords={}
            atomcount=1
            for line in results:
                linesplit=line.split()
                if len(linesplit)==1:
                    optatomnum=int(linesplit[0])
                else:
                    if len(linesplit)==4:
                        coords=[float(linesplit[1]),float(linesplit[2]),float(linesplit[3])]
                        atomindextocoords[atomcount]=coords
                        atomcount+=1
            temp=open(self.xyzfname,'r')
            results=temp.readlines()
            temp.close()
            tinkatomindextocoords={}
            atomcount=1
            for line in results:
                linesplit=line.split()
                if len(linesplit)==1:
                    atomnum=int(linesplit[0])
                else:
                    coords=[float(linesplit[2]),float(linesplit[3]),float(linesplit[4])]
                    tinkatomindextocoords[atomcount]=coords
                    atomcount+=1
            issane=True
            if optatomnum!=atomnum:
                warnings.warn('Tinker XYZ file and OPT XYZ file are not same atom number, program changes have been made between runs, need to delete and restart parameterization')
                issane=False
            else:
                for atomnum,optcoord in atomindextocoords.items():
                    tinkcoord=tinkatomindextocoords[atomnum]
                    for i in range(len(optcoord)):
                        optnum=optcoord[i]
                        tinknum=tinkcoord[i]
                        if np.abs(optnum-tinknum)>.01:
                            issane=False
                            break
                if issane==False:
                    warnings.warn('Tinker XYZ file and OPT XYZ file are not same coordinates, program changes have been made between runs, need to delete and restart parameterization')
                    sys.exit()
            return issane

        
        def SendCrashReportEmail(self,TEXT,fromaddr,toaddr,password,filename):
            msg = MIMEMultipart()
            msg['From'] = fromaddr
            msg['To'] = toaddr
            msg['Subject'] = 'Poltype Crash Report '+filename
            message = TEXT
            msg.attach(MIMEText(message, 'plain'))
            s = smtplib.SMTP_SSL('smtp.gmail.com')
            s.ehlo()
            s.login(fromaddr,password)
            text = msg.as_string()
            s.sendmail(fromaddr, [toaddr],text)
            s.quit()

 
        def SendFinalReportEmail(self,TEXT,fromaddr,toaddr,password,moleculename):
            msg = MIMEMultipart()
            msg['From'] = fromaddr
            msg['To'] = toaddr
            msg['Subject'] = 'Poltype Finished Report '+moleculename
            message = TEXT
            msg.attach(MIMEText(message, 'plain'))
            s = smtplib.SMTP_SSL('smtp.gmail.com')
            s.ehlo()
            s.login(fromaddr,password)
            text = msg.as_string()
            s.sendmail(fromaddr, [toaddr],text)
            s.quit()

        def CollectElectrostaticDipoleFitLines(self):
            files=os.listdir()
            for f in files:
                if 'poltype.log' in f:
                    name=f
            temp=open(name,'r')
            results=temp.readlines()
            temp.close()  
            fitlines=[]
            for line in results:
                if 'RMSPD =' in line or 'QMDipole' in line or 'RMSD of QM and MM' in line:
                    fitlines.append(line)


            return fitlines


        def Instructions(self):
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
            plots=[]
            fold='OPENME'
            thecurdir=os.getcwd()
            for root, subdirs, files in os.walk(os.getcwd()):
                for d in subdirs:
                    curdir=os.getcwd()
                    path=os.path.join(root, d)
                    os.chdir(path)
                    if fold in path:
                        continue
                    files=os.listdir()
                    for f in files:
                        if '.png' in f and ('energy' in f or 'fit' in f or 'water' in f):
                            plots.append(os.path.join(path,f))
            os.chdir(thecurdir)
            
            fitlines=self.CollectElectrostaticDipoleFitLines()
            instructions=self.Instructions()
            os.chdir('..')
            if not os.path.isdir(fold):
                os.mkdir(fold)
            for path in plots:
                shutil.copy(path,fold)
            os.chdir(fold)
            temp=open("README.txt",'w')
            for line in instructions:
                temp.write(line)
            for line in fitlines:
                temp.write(line)
            temp.close()
            os.chdir(thecurdir)

        def StopSimulations(self,simulationstostopfolderpath):
            for root, subdirs, files in os.walk(simulationstostopfolderpath):
                for d in subdirs:
                    curdir=os.getcwd()
                    path=os.path.join(root, d)
                    os.chdir(path)
                    files=os.listdir()
                    for f in files:
                        if '.arc' in f:
                            split=f.split('.')
                            prefix=split[0]
                            endname=prefix+'.end'
                            endname=os.path.join(path,endname)
                            with open(endname, 'w') as fp:
                                pass

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
            temp=open(self.prmfilepath,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
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
            if self.receptorligandxyzfilename!=None:
                filename='xyzedit.in'
                temp=open(filename,'w')
                temp.write('14'+'\n')
                temp.write('\n')
                temp.close()
                cmdstr=self.xyzeditpath+' '+self.receptorligandxyzfilename+' '+'-k'+' '+self.originalkeyfilename+' '+'<'+' '+filename
                
                submit.call_subsystem(self,cmdstr,wait=True)    
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
                #os.remove(self.receptorligandxyzfilename)
                os.rename(newfilename,self.receptorligandxyzfilename)


        def CheckFloat(self,string):
            try:
                float(string)
                return True
            except:
                return False

        def CheckReceptorLigandXYZFile(self):
            if self.receptorligandxyzfilename!=None:
                temp=open(self.receptorligandxyzfilename,'r')
                results=temp.readlines()
                temp.close()
                tempname=self.receptorligandxyzfilename.replace('.xyz','-t.xyz')
                temp=open(tempname,'w')
                for line in results:
                    linesplit=line.split() 
                    isboxline=True
                    if len(linesplit)==6:
                        for e in linesplit:
                            if self.CheckFloat(e)==False:
                                isboxline=False
                    else:
                        isboxline=False
                    if isboxline==True:
                        continue
                    else:
                        temp.write(line)
                temp.close()
                os.remove(self.receptorligandxyzfilename)
                os.rename(tempname,self.receptorligandxyzfilename)
                        
                         

 
        def CleanUpFiles(self):
            files=os.listdir()
            for f in files:
                for xyz in self.ligandxyzfilenamelist:
                    if xyz+'_' in f:
                        os.remove(f)
                if 'water.xyz_' in f or 'key_' in f:
                    os.remove(f)
                if self.receptorligandxyzfilename!=None:
                    if self.receptorligandxyzfilename+'_' in f:
                        os.remove(f)

        
 
            
        def SanitizeMMExecutable(self, executable):
            # Try to find Tinker executable with/without suffix
            if self.which(executable)!=None:
                return executable
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
            def is_exe(fpath):
                try:
                     return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
                except:
                     return None
        
            fpath, fname = os.path.split(program)
            if fpath:
                result=is_exe(program)
                if result:
                    return program
            else:
                for path in os.environ["PATH"].split(os.pathsep):
                    exe_file = os.path.join(path, program)
                    result=is_exe(exe_file)
                    if result:
                        return exe_file
        
            return None

        def ReadReceptorCharge(self):
            pdbmol=openbabel.OBMol()
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat('pdb')
            obConversion.ReadFile(pdbmol,self.uncomplexedproteinpdbname)
            chg=pdbmol.GetTotalCharge()
            self.receptorcharge=chg 


        def ReadLigandOBMol(self,structfname):
            tmpconv = openbabel.OBConversion()
            inFormat = openbabel.OBConversion.FormatFromExt(structfname)
            tmpconv.SetInFormat(inFormat)
            tmpmol = openbabel.OBMol()
            tmpconv.ReadFile(tmpmol, structfname)
            return tmpmol

        def ReadLigandRdkitMol(self,ligandfilename):
            tmpmol=self.ReadLigandOBMol(ligandfilename)
            tmpconv = openbabel.OBConversion()
            tmpconv.SetOutFormat('mol')
            temp=ligandfilename.replace('.sdf','.mol')
            tmpconv.WriteFile(tmpmol, temp)
            mol=rdmolfiles.MolFromMolFile(temp,removeHs=False,sanitize=False)


            return mol


        def ReadLigandCharge(self):
            tmpmol=self.ReadLigandRdkitMol(self.ligandfilename)
            tmpmol,atomindextoformalcharge=self.CheckLigandInputCharge(tmpmol)

        def PymolReadableFile(self,tinkerxyzfilename,outputname):
            if os.path.exists(tinkerxyzfilename):
                temp=open(tinkerxyzfilename,'r')
                results=temp.readlines()
                temp.close()
                temp=open(outputname,'w')
                for line in results:
                    linesplit=line.split()
                    writeline=True
                    if len(linesplit)==6:
                        second=linesplit[1]
                        try:
                            conv=float(second)
                            writeline=False
                        except:
                            continue
                    if writeline==True:
                        temp.write(line)
                temp.close()

        def ConvertTinkerXYZToCartesianXYZ(self,tinkerxyzfilename):
            xyzfilename=tinkerxyzfilename.replace('.xyz','_cart.xyz')
            temp=open(tinkerxyzfilename,'r')
            tempwrite=open(xyzfilename,'w')
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
                    if len(linesplit)>1:
                        newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
                        tempwrite.write(newline)
                        tempwrite.flush()
                        os.fsync(tempwrite.fileno())
            temp.close()
            tempwrite.close()
            return xyzfilename


        def CheckIfNeedRegrowForSolvation(self):
            needregrow=False
            if self.solvation==True:
                if self.binding==False and self.annihilatevdw==True:
                    for xyz in self.ligandxyzfilenamelist:
                        xyzfilename=self.ConvertTinkerXYZToCartesianXYZ(xyz)
                        tmpmol=self.ReadLigandOBMol(xyzfilename)
                        atomiter=openbabel.OBMolAtomIter(tmpmol)
                        for atom in atomiter:
                            atomidx=atom.GetIdx()
                            iteratomatom = openbabel.OBAtomAtomIter(atom)
                            for natom in iteratomatom:
                                natomidx=natom.GetIdx()
                                iteratomatomatom = openbabel.OBAtomAtomIter(natom)
                                for nnatom in iteratomatomatom:
                                    nnatomidx=nnatom.GetIdx()
                                    if nnatomidx!=atomidx: 
                                        iteratomatomatomatom = openbabel.OBAtomAtomIter(nnatom)
                                        for nnnatom in iteratomatomatomatom:
                                            nnnatomidx=nnnatom.GetIdx()
                                            if nnnatomidx!=natomidx:
                                                needregrow=True


            return needregrow 


        def CheckLigandInputCharge(self,molecule):
            array=[]
            totchg=0
            atomindextoformalcharge={}
            atomicnumtoformalchg={1:{2:1},5:{4:1},6:{3:-1},7:{2:-1,4:1},8:{1:-1,3:1},15:{4:1},16:{1:-1,3:1,5:-1},17:{0:-1,4:3},9:{0:-1},35:{0:-1},53:{0:-1}}
            for atom in molecule.GetAtoms():
                atomidx=atom.GetIdx()
                atomnum=atom.GetAtomicNum()
                val=atom.GetExplicitValence()
                valtochg=atomicnumtoformalchg[atomnum]
                if val not in valtochg.keys(): # then assume chg=0
                    chg=0
                else:
                    chg=valtochg[val]
                polneighb=False
                if atomnum==6:
                    for natom in atom.GetNeighbors():
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==7 or natomicnum==8 or natomicnum==16:
                            polneighb=True
                    if polneighb and val==3:
                        chg=1
                atom.SetFormalCharge(chg) 
                atom.SetNumRadicalElectrons(0)
                atomindextoformalcharge[atomidx]=chg
                totchg+=chg
                string='Atom index = '+str(atomidx+1)+' Atomic Number = ' +str(atomnum)+ ' Valence = '+str(val)+ ' Formal charge = '+str(chg)
                array.append(string)
            if self.ligandcharge!=None: 
                if self.ligandcharge!=totchg:
                    for row in array:
                        print(row,flush=True)
                    raise ValueError('Valence is not consistent with input total charge')
            else:
                self.ligandcharge=totchg
            return molecule,atomindextoformalcharge


        def CheckTinkerVersion(self):
            cmdstr=self.analyzepath+' '+os.path.abspath(os.path.join(self.annihilatorpath, os.pardir))+r'/VersionFiles/'+'water.xyz'+' '+'-k'+' '+os.path.abspath(os.path.join(self.annihilatorpath, os.pardir))+r'/VersionFiles/'+'water.key'+' '+'e'+'>'+' '+'version.out'
            
            self.WriteToLog(cmdstr)
            try:
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
                    if version.parse(self.versionnum) >= version.parse("8.9.4"):
                        latestversion = True
                        break

            if not latestversion:
                raise ValueError("Notice: Not latest working version of tinker (8.9.4)"+' '+os.getcwd())
          

        def GrabIndicesWithTypeNumber(self,xyzfilename,ligandtypes):
            indices=[]
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if len(linesplit)>1 and isboxline==False:
                    typenum=int(linesplit[5])
                    index=int(linesplit[0])
                    if typenum in ligandtypes:
                        indices.append(index)
            return indices
            
        def GrabTypeNumbers(self,xyzfilename,indices=None):   
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            typenums=[]
            for line in results:
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if len(linesplit)>1 and isboxline==False:
                    typenum=int(linesplit[5])
                    index=int(linesplit[0])
                    if indices!=None:
                        if index in indices:
                            typenums.append(typenum)
                    else:
                        typenums.append(typenum)
            return typenums 

        def GrabCoordinates(self,xyzfilename,indices=None):
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            coords=[]
            for line in results:
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if len(linesplit)>1 and isboxline==False:
                    typenum=int(linesplit[5])
                    index=int(linesplit[0])
                    coord=[float(linesplit[2]),float(linesplit[3]),float(linesplit[4])]
                    try:
                        if index in indices:
                            coords.append(coord)
                    except:
                        coords.append(coord)
            return coords


        def CompareTypes(self,receptorligandtypes,ligandtypes):
            for i in range(len(ligandtypes)):
                receptorligandtype=receptorligandtypes[i]
                ligandtype=ligandtypes[i]
                if receptorligandtype!=ligandtype:
                    raise ValueError('Ligand XYZ and Receptor-Ligand XYZ dont have the same type number order!')

        def RewriteCoordinates(self,xyzfilename,coords):
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            tempname=xyzfilename.replace('.xyz','-t.xyz')
            temp=open(tempname,'w')
            count=0
            for line in results:
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if len(linesplit)>1 and isboxline==False:
                    coord=coords[count]
                    count+=1
                    linesplit[2]=str(coord[0])
                    linesplit[3]=str(coord[1])
                    linesplit[4]=str(coord[2])
                    line=' '.join(linesplit)+'\n'
                if isboxline==False:
                    temp.write(line)
            temp.close()
            os.remove(xyzfilename)
            os.rename(tempname,xyzfilename)


        def ReadCharge(self,output,checkresid=True):
            temp=open(output,'r')
            results=temp.readlines()
            temp.close()
            chg=0
            for line in results:
                if 'Total Electric Charge :' in line:
                    linesplit=line.split()
                    chg=float(linesplit[-2])
            resid=int(round(chg))-chg
            if checkresid==True:
                if resid!=0:
                    raise ValueError('Charge is not integer! '+output+' '+str(chg)+' '+'residual charge is '+str(resid))

            return chg,resid


        def CheckEnergies(self,output):
            temp=open(output,'r')
            results=temp.readlines()
            temp.close()
            tol=5
            for line in results:
                if 'Bond Stretching' in line and 'Parameters' not in line:
                    linesplit=line.split()
                    intnum=int(linesplit[-1])
                    energy=float(linesplit[-2])
                    energyperbond=energy/intnum
                    if energyperbond>tol:
                        raise ValueError('Bad starting box structure, bond energy per bond is way to high '+str(energyperbond)+' kcal/mol')


            for line in results:
                if 'Total Electric Charge :' in line:
                    linesplit=line.split()
                    chg=float(linesplit[-2])



        def CheckInputXYZKeyFiles(self,ligandonly=False):
            keymods.RemoveKeyWords(self,self.originalkeyfilename,['parameters','axis','ewald','pme-grid','pme-order','cutoff','thermostat','integrator','ligand','verbose','archive','neighbor-list','polar-eps','polar-predict','heavy-hydrogen','omp-threads','OPENMP-THREADS'])
            string='parameters '+self.prmfilepath+'\n'
            keymods.AddKeyWord(self,self.originalkeyfilename,string)
            if self.receptorligandxyzfilename!=None and self.ligandxyzfilenamelist!=None:
                for xyz in self.ligandxyzfilenamelist:
                    ligandtypes=self.GrabTypeNumbers(xyz)  
                    indices=self.GrabIndicesWithTypeNumber(self.receptorligandxyzfilename,ligandtypes)
                    receptorligandtypes=self.GrabTypeNumbers(self.receptorligandxyzfilename,indices=indices)
                    self.CompareTypes(receptorligandtypes,ligandtypes)

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
            temp=open(xyzfilename,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                linesplit=line.split()
                if len(linesplit)>1 and '90.0' not in line:
                    index=int(linesplit[0])
                    element=linesplit[1]
                    typenum=int(linesplit[5])
                    if element=='H':
                        return index,typenum



        def GrabMultipoleParameters(self,prmfilepath,typenum):
            mpolearrays=[]
            temp=open(prmfilepath,'r')
            results=temp.readlines()
            temp.close()
            for lineidx in range(len(results)):
                line=results[lineidx]
                if 'multipole' in line:
                    linesplit=line.split()
                    if linesplit[1]==str(typenum):
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
            for idx in range(len(mpolearrays)):
                line=mpolearrays[idx]
                if 'multipole' in line:
                    chargelinesplit=line.split()
                    chargelinesplit[1]='-'+str(index)
                    charge=float(chargelinesplit[-1])
                    charge=str(charge+resid)
                    chargelinesplit[-1]=charge
                    chargeline=' '.join(chargelinesplit)+'\n'
                    mpolearrays[idx]=chargeline 

            return mpolearrays




        def ModifyCharge(self,keyfilename,mpolearrays):
            temp=open(keyfilename,'a')
            for line in mpolearrays:
                temp.write(line)
            temp.close()


        def DetermineIonsForChargedSolvation(self):
            solviontocount={}
            self.addsolvionwindows=False
            if self.binding==False and self.solvation==True and self.salthfe==True:
                for i in range(len(self.systemcharge)):
                    chg=int(self.systemcharge[i])
                    if chg>0:
                        iontypenumber=self.elementsymtotinktype['Cl']
                    elif chg<0:
                        iontypenumber=self.elementsymtotinktype['K']
                    else:
                        continue
                    count=np.abs(chg)
                    solviontocount[iontypenumber]=count
                    self.addsolvionwindows=True
            return solviontocount    

        def ChangeTypeNumbers(self,xyzfile,elementtotype):
            temp=open(xyzfile,'r')
            results=temp.readlines()
            temp.close()
            tempname=xyzfile.replace('.xyz','_new.xyz')
            temp=open(tempname,'w')
            for lineidx in range(len(results)):
                line=results[lineidx]
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if isboxline==False and lineidx!=0:
                    element=linesplit[1] 
                    typenum=elementtotype[element]
                    linesplit[5]=str(typenum)
                    line=' '.join(linesplit)+'\n'
                temp.write(line)
            temp.close()
            return tempname
        
        def ModifyBoxSizeInXYZFile(self,xyzfile,boxsize):
            temp=open(xyzfile,'r')
            results=temp.readlines()
            temp.close()
            tempname=xyzfile.replace('.xyz','_new.xyz')
            temp=open(tempname,'w')
            for lineidx in range(len(results)):
                line=results[lineidx]
                linesplit=line.split()
                isboxline=True
                if len(linesplit)==6:
                    for e in linesplit:
                        if self.CheckFloat(e)==False:
                            isboxline=False
                else:
                    isboxline=False
                if isboxline==True:
                    linesplit[0]=str(boxsize[0])
                    linesplit[1]=str(boxsize[1])
                    linesplit[2]=str(boxsize[2])
                    line=' '.join(linesplit)+'\n'
                temp.write(line)
            temp.close()
            os.remove(xyzfile)
            os.rename(tempname,xyzfile)

        def TrimPreEquilibriatedBox(self,boxsize):
            xyzfile=self.preequilboxpath
            split=os.path.split(xyzfile)
            xyzfilename=split[1]
            shutil.copy(xyzfile,os.path.join(self.simpath,xyzfilename))    
            elementtotype={'O':self.waterOtypenum,'H':self.waterHtypenum}
            filename=self.ChangeTypeNumbers(xyzfilename,elementtotype)
            self.ModifyBoxSizeInXYZFile(filename,boxsize)
            newboxsize=[2+i for i in boxsize]
            self.RemoveWaterOutsideBox(filename,boxsize)
            self.ModifyBoxSizeInXYZFile(filename,newboxsize)
            os.rename(filename,'water.xyz_2') # for xyzedit lib

        def RemoveWaterOutsideBox(self,xyzfilename,boxsize):
            tempfile=open('xyzedit.in','w')
            tempfile.write('17'+'\n')
            tempfile.write(str(boxsize[0])+','+str(boxsize[1])+','+str(boxsize[2])+'\n')
            tempfile.write('\n')
            tempfile.close()
            cmdstr=self.xyzeditpath +' '+xyzfilename+' '+self.prmfilepath+' '+' < xyzedit.in '
            submit.call_subsystem(self,cmdstr,wait=True)    
            newname=xyzfilename+'_2'
            os.rename(newname,xyzfilename)


        def ConvertTinktoXYZ(self,filename,newfilename):
            temp=open(os.getcwd()+r'/'+filename,'r')
            tempwrite=open(os.getcwd()+r'/'+newfilename,'w')
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
                    if len(linesplit)>1:
                        newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
                        tempwrite.write(newline)
                        tempwrite.flush()
                        os.fsync(tempwrite.fileno())
            temp.close()
            tempwrite.close()
            return newfilename
        

        def AlignLigandXYZToTemplateXYZ(self):
            ligandcartxyz=self.ConvertTinktoXYZ(self.ligandxyzfilenamelist[0],self.ligandxyzfilenamelist[0].replace('.xyz','_cart.xyz'))
            templateligandcartxyz=self.ConvertTinktoXYZ(self.templateligandxyzfilename,self.templateligandxyzfilename.replace('.xyz','_cart.xyz'))
            obConversion = openbabel.OBConversion()
            ref = openbabel.OBMol()
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(ligandcartxyz)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, ligandcartxyz)
            obConversion.ReadFile(ref, templateligandcartxyz)
            obConversion.SetOutFormat('mol')
            molfile=ligandcartxyz.replace('_cart.xyz','.mol')
            reffile=templateligandcartxyz.replace('_cart.xyz','.mol')
            obConversion.WriteFile(mol,molfile)
            obConversion.WriteFile(ref,reffile)
            mol=Chem.MolFromMolFile(molfile,removeHs=False,sanitize=False)
            ref=Chem.MolFromMolFile(reffile,removeHs=False,sanitize=False)
            for bond in mol.GetBonds():
                if bond.IsInRing()==False:
                    bond.SetBondType(Chem.BondType.SINGLE)
            for bond in ref.GetBonds():
                if bond.IsInRing()==False:
                    bond.SetBondType(Chem.BondType.SINGLE)

            self.ligandcharge=None
            mol,atomindextoformalcharge=self.CheckLigandInputCharge(mol)
            self.ligandcharge=None
            ref,atomindextoformalcharge=self.CheckLigandInputCharge(ref)
            self.ligandcharge=None
            Chem.SanitizeMol(mol)
            Chem.SanitizeMol(ref)
            confnum=1000
            AllChem.EmbedMultipleConfs(mol, confnum)
            mcs = rdFMCS.FindMCS([ref,mol])
            patt = Chem.MolFromSmarts(mcs.smartsString)
            atomnum=mcs.numAtoms
            refMatch = ref.GetSubstructMatch(patt)
            mv = mol.GetSubstructMatch(patt)
            rmstocid={}
            for cid in range(confnum):
                rms = AllChem.AlignMol(mol,ref,prbCid=cid,atomMap=list(zip(mv,refMatch)))
                rmstocid[rms]=cid
            minrms=min(rmstocid.keys())
            mincid=rmstocid[minrms]
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
            alignedligandcartxyz=self.ConvertTinktoXYZ(name,name.replace('.xyz','_cart.xyz'))
 
        
        def CallAnalyze(self,statexyz,statekey,alzout,analyzepath,option):
            cmdstr=analyzepath+' '+statexyz+' '+'-k'+' '+statekey+' '+option+' > '+alzout
            submit.call_subsystem(self,cmdstr,True,False,alzout) 

        def ExtractResource(self,string):
            if 'MB' in string:
                split=string.split('MB')
                memstring='MB'
            elif 'GB' in string:
                split=string.split('GB')
                memstring='GB'
            mem=float(split[0])
        
            return mem,memstring

        def PartitionResources(self):
            maxmem,memstring=self.ExtractResource(self.maxmem)
            maxmem=int(maxmem/self.jobsatsametime)
            tempmaxmem=str(maxmem)+memstring
            maxdisk,diskstring=self.ExtractResource(self.maxdisk)
            maxdisk=int(maxdisk/self.jobsatsametime)
            tempmaxdisk=str(maxdisk)+diskstring
            numproc=math.floor(int(self.numproc)/self.jobsatsametime)
            tempnumproc=str(numproc)
        
            return tempmaxmem,tempmaxdisk,tempnumproc


        def CheckIfWholeMoleculeIsARing(self,mol):
            iswholemoleculering=False
            atomindices=databaseparser.RingAtomicIndices(self,mol)
            lengthtoring={}
            for ring in atomindices:
                lengthtoring[len(ring)]=ring
            if len(lengthtoring.keys())>0:
                maxlength=max(lengthtoring.keys())
                maxring=lengthtoring[maxlength]
                if len(maxring)>7:
                    iswholemoleculering=True
            return iswholemoleculering


        def ExtractLigand(self,ligandreceptorfilename,coordinates):
            ligandpdbfilename='ligand.pdb'
            receptorpdbfilename=ligandreceptorfilename.replace('.pdb','_receptoronly.pdb')
            receptormol=self.ExtractMOLObject(ligandreceptorfilename,receptorpdbfilename,coordinates=None,receptor=True)
            ligandmol=self.ExtractMOLObject(ligandreceptorfilename,ligandpdbfilename,coordinates=coordinates,receptor=False)
        
            self.ConvertUNLToLIG(ligandpdbfilename)
            return ligandpdbfilename,receptorpdbfilename

        def ExtractMOLObject(self,ligandreceptorfilename,newpdbfilename,coordinates,receptor):
            pdbmol=self.GenerateMOLObject(ligandreceptorfilename)
            iteratom = openbabel.OBMolAtomIter(pdbmol)
            obConversion = openbabel.OBConversion()
            obConversion.SetOutFormat('pdb')
            atmindicestodelete=[]
            for atm in iteratom:
                atmindex=atm.GetIdx()
                res=atm.GetResidue()
                reskey=res.GetName()
                coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
                if reskey not in self.knownresiduesymbs and receptor==True:
                    atmindicestodelete.append(atmindex)
                elif receptor==False:
                    if coords not in coordinates:
                        atmindicestodelete.append(atmindex)
            atmindicestodelete.sort(reverse=True)
            for atmindex in atmindicestodelete:
                atm=pdbmol.GetAtom(atmindex)
                pdbmol.DeleteAtom(atm)
            obConversion.WriteFile(pdbmol,newpdbfilename)
            return pdbmol

        def GenerateMOLObject(self,pdbfilename):
            obConversion = openbabel.OBConversion()
            pdbmol = openbabel.OBMol()
            obConversion.SetInFormat('pdb')
            obConversion.ReadFile(pdbmol,pdbfilename)
            return pdbmol


        def ConvertUNLToLIG(self,filename):
            temp=open(filename,'r')
            results=temp.readlines()
            temp.close()
            tempname=filename.replace('.pdb','_TEMP.pdb')
            temp=open(tempname,'w')
            for line in results:
                linesplit=re.split(r'(\s+)', line)
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
            pdbmol=self.GenerateMOLObject(ligandpdbfilename)
            atomvecls=self.GrabAtomPositions(pdbmol)
            centroid=np.array([0.0,0.0,0.0])
            for vec in atomvecls:
                centroid+=np.array(vec)
            if len(atomvecls)==0:
                raise ValueError('Ligand in PDB file is not labeled as LIG')
            centroid=centroid/len(atomvecls)
        
            
            return centroid


        def GrabAtomPositions(self,pdbmol):
            atomvecls=[]
            iteratombab = openbabel.OBMolAtomIter(pdbmol)
            for atm in iteratombab:
                atmindex=atm.GetIdx()
                coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
                atomvecls.append(coords)
        
         
            return atomvecls


        def FindConsecutiveIndices(self,nums):
            consec=[nums[0]]
            allconsec=[]
            for i in range(1,len(nums)):
                num=nums[i]
                prevnum=nums[i-1]
                if num==prevnum+1:
                    consec.append(num)
                else:
                    allconsec.append(consec)
                    consec=[num]
            allconsec.append(consec)
            return allconsec


        def RepresentsInt(self,s):
            try: 
                int(s)
                return True
            except ValueError:
                return False


        def ShiftParameterTypes(self,filename,oldtypetonewtype):
            tempname=filename+"_TEMP"
            temp=open(filename,'r')
            results=temp.readlines()
            temp.close()
            temp=open(tempname,'w')
            for line in results:
                linesplitall=re.split(r'(\s+)', line)
                for i in range(len(linesplitall)):
                    element=linesplitall[i]
                    if '.xyz' in filename and (i!=12):
                        continue
                    if self.RepresentsInt(element):
                        oldtypenum=np.abs(int(element))
                        if oldtypenum in oldtypetonewtype.keys():
                            newtypenum=oldtypetonewtype[oldtypenum]
                            typenum=newtypenum
                        else:
                            typenum=oldtypenum
                        if '-' in element:
                            typenum=-typenum
                        linesplitall[i]=str(typenum)
                line=''.join(linesplitall)
                temp.write(line)
            temp.close()
            os.rename(tempname,filename)


        def GrabMaxTypeNumber(self,parameterfile):
            maxnumberfromprm=1
            temp=open(parameterfile,'r')
            results=temp.readlines()
            temp.close()
            for line in results:
                if 'atom' in line and '#' not in line:
                    linesplit=line.split()
                    atomtype=int(linesplit[1])
                    if atomtype>maxnumberfromprm:
                        maxnumberfromprm=atomtype
            return maxnumberfromprm
        
        def GrabMinTypeNumber(self,parameterfile):
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
            newkeyfilelist=keyfilelist.copy()
            oldtypetonewtypelist=[]
            currentmin=100000
            currentmax=0
            shift=0
            prevmaxnumberfromkey=0
            prevminnumberfromkey=0
            originalshift=self.GrabMaxTypeNumber(self.prmfilepath)
            for keyfilename in newkeyfilelist:
                maxnumberfromkey=self.GrabMaxTypeNumber(keyfilename)
                minnumberfromkey=self.GrabMinTypeNumber(keyfilename)
                shift+=prevmaxnumberfromkey-prevminnumberfromkey+1+originalshift
                if minnumberfromkey<currentmin:
                    currentmin=minnumberfromkey
                if maxnumberfromkey>currentmax:
                    currentmax=maxnumberfromkey
                types=np.arange(minnumberfromkey,maxnumberfromkey+1,1)
                shiftedtypes=types+shift
                maxtype=max(shiftedtypes)
                currentmax=maxtype
                oldtypetonewtype=dict(zip(types,shiftedtypes))
                temp={}
                for oldtype,newtype in oldtypetonewtype.items():
                    negold=-oldtype
                    negnew=-newtype
                    temp[negold]=negnew
                oldtypetonewtype.update(temp)
                oldtypetonewtypelist.append(oldtypetonewtype)
                prevmaxnumberfromkey=maxnumberfromkey
                prevminnumberfromkey=minnumberfromkey
                
            return oldtypetonewtypelist


        def CheckNetChargeIsZero(self,xyzpath,keypath,alzout):
            self.CallAnalyze(xyzpath,keypath,alzout,self.trueanalyzepath,'m')
            charge,resid=self.ReadCharge(alzout)
            if charge!=0:
                raise ValueError('Net charge is not zero! Net Charge = '+str(charge)+' '+keypath)



if __name__ == '__main__':
    def RunPoltype():
        poltype=PolarizableTyper() 
        try:
            poltype.main()
        except:
            traceback.print_exc(file=sys.stdout)
            text = str(traceback.format_exc())
            if os.path.exists(poltype.scrtmpdirgau):
                shutil.rmtree(poltype.scrtmpdirgau)
            if os.path.exists(poltype.scrtmpdirpsi4):
                shutil.rmtree(poltype.scrtmpdirpsi4)
            if poltype.email!=None:
                password='amoebaisbest'
                fromaddr = 'poltypecrashreportnoreply@gmail.com'
                toaddr = poltype.email
                filename=poltype.logfname
                poltype.WriteToLog(text)
                poltype.WriteToLog('Poltype has crashed!')
                try:
                    poltype.SendCrashReportEmail(text,fromaddr,toaddr,password,filename)
                except:
                    pass
            raise ValueError('Houston, we have a problem. Buy a developer some coffee!')
    RunPoltype()
