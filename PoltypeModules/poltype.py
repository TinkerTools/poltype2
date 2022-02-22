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
import openbabel
import shutil
import time
import copy
import getopt
import databaseparser
import dimorphite_dl
import forcebalancepoltypeinterface as fb
import torsiongenerator as torgen
import modifiedresidues as modres
import symmetry as symm
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



class PolarizableTyper():
    def __init__(self,fitred=False,vdwtypestoeval=['S','T','D'],vdwprmtypestofit=['S','T'],lastlogfileupdatetime=1,addlonepairvdwsites=False,quickdatabasesearch=False,genprotstatesonly=False,generateextendedconf=True,onlyvdwatomlist=None,poltypepathlist=None,vdwtypeslist=None,fittypestogether=None,csvexpdatafile=None,liquid_equ_steps=10000,liquid_prod_steps=5000000,liquid_timestep=1.0,liquid_interval=0.1,gas_equ_steps=500000,gas_prod_steps=1000000,gas_timestep=1.0,gas_interval=0.1,md_threads=4,liquid_prod_time=5,gas_prod_time=5,WQ_PORT=None,parentjobsatsametime=1,coresperjob=2,addhydrogens=False,maximizejobsatsametime=True,consumptionratio=.8,scratchpath='/scratch',nonaroringtor1Dscan=False,skipespfiterror=False,vdwmaxqmstartingpointspertype=1,vdwmaxtinkergridpoints=50,smallmoleculefragmenter=False,fragmentjobslocal=False,toroptdebugmode=False,debugmode=False,fragmenterdebugmode=False,jobsatsametime=0,usepoleditframes=False,databasematchonly=False,setupfragjobsonly=False,allowradicals=False,checkinputonly=False,username=None,esprestweight=1,espgrad=.1,issane=True,deletedfiles=False,onlyfittorstogether=[],parentname=None,addhydrogentononcharged=True,accuratevdwsp=False,inputmoleculefolderpaths=None,email=None,firstoptfinished=False,optonly=False,onlyvdwatomindex=None,use_qmopt_vdw=False,use_gau_vdw=False,dontusepcm=False,deleteallnonqmfiles=True,totalcharge=None,torspbasissethalogen="6-311G*",homodimers=False,tortormissingfilename='tortormissing.txt',tordebugmode=False,amoebapluscfsmartstocommentmap=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebapluscfsmartstocomment.txt',amoebapluscfprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'cfprmlib.txt',amoebaplusnonbondedprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebaplusnonbonded.prm',amoebaplusnonbondedsmartstocommentmap=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebaplusnonbonded.txt',smartstosoluteradiimap=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'SMARTsToSoluteRadiiMap.txt',latestsmallmoleculepolarizeprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarize.prm',latestsmallmoleculesmartstotypespolarize=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarcommenttoparameters.txt',latestsmallmoleculesmartstotinkerclass=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21smartstoclass.txt',latestsmallmoleculeprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21.prm',boltzmantemp=8,dovdwscan=False,vdwprobepathname=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/VdwProbes/',vdwprobenames=['water'],use_gausgeomoptonly=False,maxtorRMSPDRel=.2,vdwmissingfilename='missingvdw.txt',databaseprmfilename='database.prm',tortor=False,torfit2Drotonly=False,torfit1Drotonly=False,externalparameterdatabase=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'externalparameterdatabase.txt',fitfirsttorsionfoldphase=False,keyfiletoaddtodatabase=None,skipgridsearch=True,torsionprmguessfilename='torsionprmguess.txt',defaultmaxtorsiongridpoints=40,torsionsmissingfilename='torsionsmissing.txt',smallmoleculemm3prmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'mm3.prm',smallmoleculesmartstomm3descrip=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstomm3typedescrip.txt',absdipoletol=.5,transferanyhydrogentor=True,smallmoleculesmartstotinkerdescrip=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstoamoebatypedescrip.txt',smallmoleculeprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba09.prm',torspbasissetfile='6-311+g_st_.0.gbs',toroptbasissetfile='6-311g_st_.0.gbs',optbasissetfile='6-311g_st_.0.gbs',dmabasissetfile='6-311g_st__st_.0.gbs',espbasissetfile='aug-cc-pvtz.1.gbs',iodinetorspbasissetfile='def2-svp.1.gbs',iodinetoroptbasissetfile='def2-svp.1.gbs',iodineoptbasissetfile='def2-svp.1.gbs',iodinedmabasissetfile='def2-svp.1.gbs',iodineespbasissetfile='def2-tzvpp.1.gbs',basissetpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/'+'BasisSets/',refinenonaroringtors=False,maxgrowthcycles=4,use_gauPCM=False,fitqmdipole=False,scfmaxiter=500,suppresstorfiterr=False,obminimizeexe='obminimize',readinionly=False,suppressdipoleerr=False,topologylib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ModifiedResidueLibraries/residue_connect.txt",poltypepath=os.path.abspath(os.path.split(__file__)[0]),WBOtol=.05,dontfrag=False,isfragjob=False,dipoletol=.5,externalapi=None,printoutput=False,poltypeini=True,structure=None,prmstartidx=401,numproc=None,maxmem=None,maxdisk=None,gausdir=None,gdmadir=None,tinkerdir=None,scratchdir="/scratch",paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18_header.prm",gausexe=None,formchkexe='formchk',cubegenexe='cubegen',gdmaexe='gdma',avgmpolesexe=os.path.abspath(os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), os.pardir)) + "/PoltypeModules/avgmpoles.pl",peditexe='poledit.x',potentialexe='potential.x',minimizeexe='minimize.x',analyzeexe='analyze.x',superposeexe='superpose.x',defopbendval=0.20016677990819662,Hartree2kcal_mol=627.5095,optbasisset='6-31G*',toroptbasisset='6-311G',dmabasisset='6-311G**',espbasisset="aug-cc-pVTZ",torspbasisset="6-311+G*",optmethod='MP2',toroptmethod='wB97X-D',torspmethod='wB97X-D',dmamethod='MP2',espmethod='MP2',qmonly = False,espfit = True,parmtors = True,foldnum=3,foldoffsetlist = [ 0.0, 180.0, 0.0, 180.0, 0.0, 180.0 ],torlist = None,rotbndlist = None,maxRMSD=1,maxRMSPD=1,maxtorRMSPD=1.8,tordatapointsnum=None,gentorsion=False,gaustorerror=False,torsionrestraint=.1*3282.80354574,onlyrotbndslist=None,rotalltors=False,dontdotor=False,dontdotorfit=False,toroptpcm=False,optpcm=False,torsppcm=False,use_gaus=False,use_gausoptonly=False,freq=False,postfit=False,bashrcpath=None,amoebabioprmpath=None,libpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ModifiedResidueLibraries/lib.bio18_conv1.txt",SMARTSToTypelibpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ModifiedResidueLibraries/SMARTSToTypeLib.txt',ModifiedResiduePrmPath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/ModifiedResidue.prm',modifiedproteinpdbname=None,unmodifiedproteinpdbname=None,mutatedsidechain=None,mutatedresiduenumber=None,modifiedresiduepdbcode=None,optmaxcycle=400,torkeyfname=None,gausoptcoords='',forcefield="AMOEBA",helpfile='README.md',versionfile='version.md',sleeptime=5):
        self.fitred=fitred
        self.vdwtypestoeval=vdwtypestoeval
        self.vdwprmtypestofit=vdwprmtypestofit 
        self.lastlogfileupdatetime=lastlogfileupdatetime
        self.addlonepairvdwsites=addlonepairvdwsites
        self.quickdatabasesearch=quickdatabasesearch
        self.genprotstatesonly=genprotstatesonly
        self.generateextendedconf=generateextendedconf
        self.onlyvdwatomlist=onlyvdwatomlist
        self.poltypepathlist=poltypepathlist
        self.vdwtypeslist=vdwtypeslist
        self.fittypestogether=fittypestogether
        self.csvexpdatafile=csvexpdatafile
        self.liquid_equ_steps=liquid_equ_steps
        self.liquid_prod_steps=liquid_prod_steps
        self.liquid_timestep=liquid_timestep
        self.liquid_interval=liquid_interval
        self.gas_equ_steps=gas_equ_steps
        self.gas_prod_steps=gas_prod_steps
        self.gas_timestep=gas_timestep
        self.gas_interval=gas_interval
        self.md_threads=md_threads
        self.liquid_prod_time=liquid_prod_time
        self.gas_prod_time=gas_prod_time
        self.WQ_PORT=WQ_PORT
        self.parentjobsatsametime=parentjobsatsametime
        self.coresperjob=coresperjob
        self.addhydrogens=addhydrogens
        self.maximizejobsatsametime=maximizejobsatsametime 
        self.consumptionratio=consumptionratio
        self.scratchpath=scratchpath
        self.nonaroringtor1Dscan=nonaroringtor1Dscan
        self.skipespfiterror=skipespfiterror 
        self.vdwmaxqmstartingpointspertype=vdwmaxqmstartingpointspertype
        self.vdwmaxtinkergridpoints=vdwmaxtinkergridpoints

        self.smallmoleculefragmenter=smallmoleculefragmenter 
        self.fragmentjobslocal=fragmentjobslocal
        self.toroptdebugmode=toroptdebugmode
        self.debugmode=debugmode
        self.fragmenterdebugmode=fragmenterdebugmode
        self.jobsatsametime=jobsatsametime
        self.usepoleditframes=usepoleditframes
        self.databasematchonly=databasematchonly
        self.setupfragjobsonly=setupfragjobsonly
        self.allowradicals=allowradicals
        self.checkinputonly=checkinputonly
        self.username=username
        self.esprestweight=esprestweight
        self.espgrad=espgrad
        self.issane=issane
        self.deletedfiles=deletedfiles
        self.onlyfittorstogether=onlyfittorstogether
        self.parentname=parentname
        self.addhydrogentononcharged=addhydrogentononcharged
        self.accuratevdwsp=accuratevdwsp
        self.inputmoleculefolderpaths=inputmoleculefolderpaths
        self.email=email
        self.firstoptfinished=firstoptfinished
        self.onlyvdwatomindex=onlyvdwatomindex
        self.optonly=optonly
        self.use_qmopt_vdw=use_qmopt_vdw
        self.use_gau_vdw=use_gau_vdw
        self.dontusepcm=dontusepcm
        self.deleteallnonqmfiles=deleteallnonqmfiles
        self.totalcharge=totalcharge
        self.torspbasissethalogen=torspbasissethalogen
        self.homodimers=homodimers
        self.tortormissingfilename=tortormissingfilename
        self.tordebugmode=tordebugmode
        self.amoebapluscfprmlib=amoebapluscfprmlib
        self.amoebapluscfsmartstocommentmap=amoebapluscfsmartstocommentmap
        self.amoebaplusnonbondedsmartstocommentmap=amoebaplusnonbondedsmartstocommentmap
        self.amoebaplusnonbondedprmlib=amoebaplusnonbondedprmlib
        self.smartstosoluteradiimap=smartstosoluteradiimap
        self.latestsmallmoleculepolarizeprmlib=latestsmallmoleculepolarizeprmlib
        self.latestsmallmoleculesmartstotypespolarize=latestsmallmoleculesmartstotypespolarize
        self.latestsmallmoleculesmartstotinkerclass=latestsmallmoleculesmartstotinkerclass
        self.latestsmallmoleculeprmlib=latestsmallmoleculeprmlib
        self.boltzmantemp=boltzmantemp
        self.dovdwscan=dovdwscan
        self.vdwprobepathname=vdwprobepathname
        self.vdwprobenames=vdwprobenames
        self.use_gausgeomoptonly=use_gausgeomoptonly
        self.maxtorRMSPDRel=maxtorRMSPDRel
        self.vdwmissingfilename=vdwmissingfilename
        self.databaseprmfilename=databaseprmfilename
        self.tortor=tortor
        self.torfit2Drotonly=torfit2Drotonly
        self.torfit1Drotonly=torfit1Drotonly
        self.externalparameterdatabase=externalparameterdatabase
        self.fitfirsttorsionfoldphase=fitfirsttorsionfoldphase
        self.keyfiletoaddtodatabase=keyfiletoaddtodatabase
        self.skipgridsearch=skipgridsearch
        self.torsionprmguessfilename=torsionprmguessfilename
        self.defaultmaxtorsiongridpoints=defaultmaxtorsiongridpoints
        self.torsionsmissingfilename=torsionsmissingfilename
        self.smallmoleculemm3prmlib=smallmoleculemm3prmlib
        self.smallmoleculesmartstomm3descrip=smallmoleculesmartstomm3descrip
        self.transferanyhydrogentor=transferanyhydrogentor

        self.absdipoletol=absdipoletol
        self.torspbasissetfile=torspbasissetfile
        self.toroptbasissetfile=toroptbasissetfile
        self.optbasissetfile=optbasissetfile
        self.dmabasissetfile=dmabasissetfile
        self.espbasissetfile=espbasissetfile
        self.iodinetorspbasissetfile=iodinetorspbasissetfile
        self.iodinetoroptbasissetfile=iodinetoroptbasissetfile
        self.iodineoptbasissetfile=iodineoptbasissetfile
        self.iodinedmabasissetfile=iodinedmabasissetfile
        self.iodineespbasissetfile=iodineespbasissetfile
        self.basissetpath=basissetpath
        self.smallmoleculeprmlib=smallmoleculeprmlib
        self.smallmoleculesmartstotinkerdescrip=smallmoleculesmartstotinkerdescrip
        self.refinenonaroringtors=refinenonaroringtors
        self.fitqmdipole=fitqmdipole
        self.maxgrowthcycles=maxgrowthcycles
        self.use_gauPCM=use_gauPCM
        self.scfmaxiter=scfmaxiter
        self.suppresstorfiterr=suppresstorfiterr
        self.obminimizeexe=obminimizeexe
        self.readinionly=readinionly
        self.suppressdipoleerr=suppressdipoleerr
        self.use_gaus=use_gaus
        self.use_gausoptonly=use_gausoptonly
        self.topologylibpath=topologylib
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
        self.qmonly = qmonly
        self.sleeptime = sleeptime
        self.espfit = espfit
        self.parmtors = parmtors
        self.foldnum=foldnum
        self.nfoldlist =  list(range(1,self.foldnum+1))
        self.foldoffsetlist = foldoffsetlist
        self.parentdir=os.getcwd()
        if torlist==None:
            self.torlist = []
        else:
            self.torlist=torlist
        if rotbndlist==None:
            self.rotbndlist = []
        else:
            self.rotbndlist=rotbndlist
           
            
        self.maxRMSD=maxRMSD
        self.maxRMSPD=maxRMSPD
        self.maxtorRMSPD=maxtorRMSPD
        self.tordatapointsnum=tordatapointsnum
        self.gentorsion=gentorsion
        self.gaustorerror=gaustorerror
        self.torsionrestraint=torsionrestraint
        if onlyrotbndslist==None:
            self.onlyrotbndslist=[]
        else:
           
            self.onlyrotbndslist=onlyrotbndslist.split(',')
            templist=[]
            for ele in self.onlyrotbndslist:
                nums=ele.lstrip().rstrip().split()
                temp=[]
                for e in nums:
                    temp.append(int(e))
                templist.append(temp)
            self.onlyrotbndslist=templist
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
        self.forcefield=forcefield
        self.helpfile=helpfile
        self.versionfile=versionfile 
        self.optmethod=optmethod               
        self.toroptmethod=toroptmethod              
        self.torspmethod=torspmethod                   
        self.dmamethod=dmamethod                    
        self.espmethod=espmethod                 
        
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
                        if a=='None':
                            continue
                    else:
                        newline=line

                    if "rotalltors" in newline:
                        if '=' not in line:
                            self.rotalltors = True
                        else:
                            self.rotalltors=self.GrabBoolValue(a)


                    elif "vdwprmtypestofit" in newline:
                        self.vdwprmtypestofit=a.split(',')
                        self.vdwprmtypestofit=[i.strip() for i in self.vdwprmtypestofit]


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


                    elif "consumptionratio" in newline:
                        self.consumptionratio=float(a)
                    elif "scratchpath" in newline:
                        self.scratchpath=a
                    elif "jobsatsametime" in newline and 'max' not in newline and 'parentjobsatsametime' not in newline:
                        self.jobsatsametime=int(a)
                    elif "esprestweight" in newline:
                        self.esprestweight=float(a)
                    elif "espgrad" in newline:
                        self.espgrad=a
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

                    elif "username" in newline:
                        self.username=a
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


                    elif "tortor" in newline:
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
                    elif "amoebabioprmpath" in newline:
                        self.amoebabioprmpath = a
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
                    elif "help" in newline:
                        self.copyright()
                        self.usage()
                        sys.exit(2)
                    else:
                        print('Unrecognized '+line)
                        self.usage()
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


        self.startdir=os.getcwd()
        self.SanitizeAllQMMethods()
        if self.readinionly==True:
            return
        self.SanitizeMMExecutables()
        self.copyright()
        if self.poltypepathlist!=None:
            fb.GenerateForceBalanceInputs(self.poltypepathlist,self.vdwtypeslist,self.liquid_equ_steps,self.liquid_prod_steps,self.liquid_timestep,self.liquid_interval,self.gas_equ_steps,self.gas_prod_steps,self.gas_timestep,self.gas_interval,self.md_threads,self.liquid_prod_time,self.gas_prod_time,self.WQ_PORT,self.csvexpdatafile,self.fittypestogether,self.vdwprmtypestofit,self.vdwtypestoeval)
            sys.exit()
        head, self.molstructfname = os.path.split(self.molstructfname)
        self.molecprefix =  os.path.splitext(self.molstructfname)[0]
        self.initialize()
        if self.isfragjob==False:
            self.parentname=self.molecprefix
        else:
            self.parentname=str(self.parentname)
        
        self.init_filenames()
 
        # Use openbabel to create a 'mol' object from the input molecular structure file. 
        # Openbabel does not play well with certain molecular structure input files,
        # such as tinker xyz files. (normal xyz are fine)
    
        
        if not __name__ == '__main__':
            params=self.main()


    def GrabProtStates(self,m):
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
        for i,mol in enumerate(protonated_mols):
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            smi=finalsmi[i]
            name='ProtonationState_'+str(i)+'.mol'
            rdmolfiles.MolToMolFile(mol,name)

        return finalsmi



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
 

    def WriteToLog(self,string):
        now = time.strftime("%c",time.localtime())
        self.logfh.write(now+' '+string+'\n')
        self.logfh.flush()
        os.fsync(self.logfh.fileno())

        
    def SanitizeMMExecutables(self):
        self.peditexe=self.SanitizeMMExecutable(self.peditexe)
        self.potentialexe=self.SanitizeMMExecutable(self.potentialexe)
        self.minimizeexe=self.SanitizeMMExecutable(self.minimizeexe)
        self.analyzeexe=self.SanitizeMMExecutable(self.analyzeexe)
        self.superposeexe=self.SanitizeMMExecutable(self.superposeexe)

    def SanitizeMMExecutable(self, executable):
        # Try to find Tinker executable with/without suffix
        if self.tinkerdir is None:
            self.tinkerdir = os.getenv("TINKERDIR", default="")
        exe = os.path.join(self.tinkerdir, executable)
        if self.which(exe) is None:
            exe = exe[:-2] if exe.endswith('.x') else exe + '.x'
            if self.which(exe) is None:
                print("ERROR: Cannot find Tinker {} executable".format(exe))
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
                if version.parse(self.versionnum) >= version.parse("8.9.4"):
                    latestversion = True
                    break

        if not latestversion:
            raise ValueError("Notice: Not latest working version of tinker (8.9.4)"+' '+os.getcwd())
        
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
    
    def usage(self):
        self.printfile(self.helpfile)
    
        
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

    

    
    def CallJobsSeriallyLocalHost(self,fulljobtooutputlog,skiperrors,wait=False):
       thepath=os.path.join(os.getcwd(),'Fragments')
       finishedjobs=[]
       errorjobs=[]
       submittedjobs=[]
       errormessages=[]
           
       while len(finishedjobs)!=len(list(fulljobtooutputlog.keys())):
           for job,outputlog in fulljobtooutputlog.items():
               if job not in finishedjobs:
                  finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
                  if finished==True:
                      if outputlog not in finishedjobs:
                          finishedjobs.append(outputlog)
                          self.NormalTerm(outputlog)
                          if job in submittedjobs:
                              submittedjobs.remove(job) 
                  if error==True and job in submittedjobs:
                      if outputlog not in finishedjobs:
                          errorjobs.append(outputlog)
                          finishedjobs.append(outputlog) 
                          self.ErrorTerm(outputlog,skiperrors)
                          submittedjobs.remove(job)
                  if job not in submittedjobs and len(submittedjobs)<self.jobsatsametime and finished==False and outputlog not in finishedjobs:
                      count=len(finishedjobs)
                      ratio=(100*count)/len(fulljobtooutputlog.keys())
                      if len(finishedjobs)!=0:
                          if 'poltype.py' in job:
                              self.ETAQMFinish(thepath,len(fulljobtooutputlog.keys()))
                      
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
        while len(finishedjobs)!=len(jobtooutputlog.keys()):
            for job in jobtooutputlog.keys():
                outputlog=jobtooutputlog[job]
                finished,error,errormessages=self.CheckNormalTermination(outputlog,errormessages,skiperrors)
                if finished==True and error==False: # then check if SP has been submitted or not
                    if outputlog not in finishedjobs:
                        self.NormalTerm(outputlog)
                        finishedjobs.append(outputlog)
                elif finished==False and error==True:
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
                    if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line or "LBFGS  --  Normal Termination due to SmallGrad" in line or "Normal termination" in line or 'Normal Termination' in line or 'Total Potential Energy' in line:
                        term=True
                    if ('Tinker is Unable to Continue' in line or 'error' in line or 'Error' in line or 'ERROR' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation, address not mapped to object' in line or 'galloc:  could not allocate memory' in line or 'Erroneous write.' in line) and 'DIIS' not in line and 'mpi' not in line and 'RMS Error' not in line:
                        error=True
                        errorline=line
                    if 'segmentation violation' in line and 'address not mapped to object' not in line or 'Waiting' in line or ('OptimizationConvergenceError' in line and 'except' in line) or "Error on total polarization charges" in line or 'Erroneous write' in line:
                        error=False
                        continue
                    if ('Error termination request processed by link 9999' in line or 'Error termination via Lnk1e in' in line) or ('OptimizationConvergenceError' in line and 'except' not in line) or 'Could not converge geometry optimization' in line or 'SCFConvergenceError' in line or 'Incomplete Convergence due to' in line:
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
        if skiperrors==False:
            self.WriteToLog("ERROR termination: logfile=%s path=%s" % (logfname,os.getcwd()))


    def call_subsystem(self,cmdstrs,wait=False,skiperrors=False):
        if self.printoutput==True:
            for cmdstr in cmdstrs:
                print("Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
        procs=[]
        for cmdstr in cmdstrs:
            self.WriteToLog("Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
            p = subprocess.Popen(cmdstr,shell=True,stdout=self.logfh, stderr=self.logfh)
            procs.append(p)

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
            if 'torsion' in line and '#' not in line and 'Missing' not in line and 'none' not in line:
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
        
         
        if self.amoebabioprmpath!=None and (self.modifiedproteinpdbname!=None or self.unmodifiedproteinpdbname!=None):
            knownresiduesymbs,modproidxs,proboundidxs,boundaryatomidxs,proOBmol,molname,modresiduelabel,proidxtoligidx,ligidxtoproidx,modmol,smarts,check,connectedatomidx,backboneindexesreference,modligidxs=modres.GenerateModifiedProteinPoltypeInput(self)
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

    def CheckBondTopology(self,outputlog,rdkitmol):
        bondtoposame=True
        if self.use_gaus==False and self.use_gausoptonly==False:
            fname=outputlog.replace('.log','.xyz')
        else:
            fname=outputlog
        outputname='preQMopt.mol'
        rdmolfiles.MolToMolFile(rdkitmol,outputname)
        optmol = opt.load_structfile(self,fname)
        inioptmol = opt.load_structfile(self,outputname)
        issame=opt.CheckBondConnectivity(self,inioptmol,optmol,outputlog.replace('.log','.xyz'))
        if issame==False:
            bondtoposame=False
        if self.fullopt==False: 
            isnear=opt.CompareBondLengths(self,inioptmol,optmol,outputlog)
        else:
            isnear=True
        if isnear==False:
            bondtoposame=False
    
        return bondtoposame

    def CheckIfTorsionUndefined(self,listoftorsionsforprm,conf): # sometimes rdkit conformation has 3 linear atoms
        isitsafe=True
        for tor in listoftorsionsforprm:
            middle=[tor[1],tor[2]]
            firstangle=rdMolTransforms.GetAngleDeg(conf,tor[0],tor[1],tor[2])
            secondangle=rdMolTransforms.GetAngleDeg(conf,tor[1],tor[2],tor[3])

            if firstangle<0:
                firstangle=firstangle+360
            if secondangle<0:
                secondangle=secondangle+360
            angletol=2
            if np.abs(180-firstangle)<=2 or np.abs(180-secondangle)<=2:
                isitsafe=False
                break

        return isitsafe



    def GenerateExtendedConformer(self,rdkitmol,mol):
        numconf=100 # just try this
        AllChem.EmbedMultipleConfs(rdkitmol, numConfs=numconf)
        confs=rdkitmol.GetConformers()
        listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=databaseparser.GrabAtomsForParameters(self,mol)
        disttoconf={}
        for i in range(len(confs)):
            conf=confs[i]
            name="conftest.mol"
            rdmolfiles.MolToMolFile(rdkitmol,name,confId=i)
            mol=rdmolfiles.MolFromMolFile(name,removeHs=False)
            isitsafe=self.CheckIfTorsionUndefined(listoftorsionsforprm,conf)
            if isitsafe==False:
                continue
            maxdist=self.FindLongestDistanceInMolecule(mol)

            disttoconf[maxdist]=i
        distances=list(disttoconf.keys())
        if len(distances)!=0:
            maxdist=max(distances)
            confindex=disttoconf[maxdist]
        else:
            confindex=0
        indextocoordinates={}
        rdmolfiles.MolToMolFile(rdkitmol,name,confId=confindex)
        mol=rdmolfiles.MolFromMolFile(name,removeHs=False)
        for i in range(len(mol.GetAtoms())):
            pos = mol.GetConformer().GetAtomPosition(i) 
            vec=np.array([float(pos.x),float(pos.y),float(pos.z)])
            indextocoordinates[i]=vec
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
        atleastonehashydrogen=False
        for atomindex in chargedindices:
            atom=m.GetAtomWithIdx(atomindex)
            for natom in atom.GetNeighbors():
                natomicnum=natom.GetAtomicNum()
                if natomicnum==1:
                    atleastonehashydrogen=True
        if atleastonehashydrogen==True:         
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
        for f in files:
            if not os.path.isdir(f) and 'nohup' not in f and f[0]!='.' and f!='parentvdw.key':
                fsplit=f.split('.')
                if len(fsplit)>1:
                    end=fsplit[1]
                    if 'log' not in end and 'sdf' not in end and 'ini' not in end and 'chk' not in end and 'dat' not in end and 'mol' not in end and 'txt' not in end: 
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
            string='--job='+job+' '+'--numproc='+str(1)+' '+'--ram=10GB'+' '+'--disk=0GB'+' '+'--inputfilepaths='+os.path.join(jobpath,'poltype.ini')+' '+'--username='+self.username+'\n'
            temp.write(string)
        temp.close()


    def CreatePoltypeInputFilesMultipleMolecules(self):
        dic={'username':self.username,'externalapi':self.externalapi,'accuratevdwsp':self.accuratevdwsp,'email':self.email,'firstoptfinished':self.firstoptfinished,'optonly':self.optonly,'onlyvdwatomindex':self.onlyvdwatomindex,'use_qmopt_vdw':self.use_qmopt_vdw,'use_gau_vdw':self.use_gau_vdw,'dontusepcm':self.dontusepcm,'deleteallnonqmfiles':self.deleteallnonqmfiles,'totalcharge':self.totalcharge,'torspbasissethalogen':self.torspbasissethalogen,'homodimers':self.homodimers,'boltzmantemp':self.boltzmantemp,'dovdwscan':self.dovdwscan,'use_gausgeomoptonly':self.use_gausgeomoptonly,'maxtorRMSPDRel':self.maxtorRMSPDRel,'tortor':self.tortor,'fitfirsttorsionfoldphase':self.fitfirsttorsionfoldphase,'defaultmaxtorsiongridpoints':self.defaultmaxtorsiongridpoints,'absdipoletol':self.absdipoletol,'refinenonaroringtors':self.refinenonaroringtors,'maxgrowthcycles':self.maxgrowthcycles,'use_gauPCM':self.use_gauPCM,'fitqmdipole':self.fitqmdipole,'WBOtol':self.WBOtol,'dontfrag':self.dontfrag,'dipoletol':self.dipoletol,'numproc':self.numproc,'maxmem':self.maxmem,'maxdisk':self.maxdisk,'optbasisset':self.optbasisset,'toroptbasisset':self.toroptbasisset,'dmabasisset':self.dmabasisset,'espbasisset':self.espbasisset,'torspbasisset':self.torspbasisset,'optmethod':self.optmethod,'toroptmethod':self.toroptmethod,'torspmethod':self.torspmethod,'dmamethod':self.dmamethod,'espmethod':self.espmethod,'qmonly' : self.qmonly,'espfit' : self.espfit,'foldnum':self.foldnum,'maxRMSD':self.maxRMSD,'maxRMSPD':self.maxRMSPD,'maxtorRMSPD':self.maxtorRMSPD,'tordatapointsnum':self.tordatapointsnum,'torsionrestraint':self.torsionrestraint,'rotalltors':self.rotalltors,'dontdotor':self.dontdotor,'dontdotorfit':self.dontdotorfit,'toroptpcm':self.toroptpcm,'optpcm':self.optpcm,'torsppcm':self.torsppcm,'use_gaus':self.use_gaus,'use_gausoptonly':self.use_gausoptonly,'freq':self.freq,'optmaxcycle':self.optmaxcycle,'forcefield':self.forcefield}
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
        Chem.SanitizeMol(m)
        smarts=rdmolfiles.MolToSmarts(m)
        if '.' in smarts:
            raise ValueError('Multiple fragments detectected in input molecule')
        pcm=self.CheckForConcentratedFormalCharges(m,atomindextoformalcharge)
        cpm = copy.deepcopy(m)
        if self.firstoptfinished==False and self.isfragjob==False and self.generateextendedconf==True:
            indextocoordinates=self.GenerateExtendedConformer(m,mol)
        Chem.GetSymmSSSR(m)
        m.GetRingInfo().NumRings() 
        m=self.AddInputCoordinatesAsDefaultConformer(m,indextocoordinates)
        rdmolfiles.MolToMolFile(m,'test.mol')
        mol,m=self.CheckIsInput2D(mol,obConversion,m)
        if not os.path.exists(self.scrtmpdirpsi4):
            os.mkdir(self.scrtmpdirpsi4)
        if not os.path.exists(self.scrtmpdirgau):
            os.mkdir(self.scrtmpdirgau)

        mol=self.SetDefaultCoordinatesBabel(mol,indextocoordinates)
        self.mol=mol

        self.rdkitmol=m
        self.mol.SetTotalCharge(self.totalcharge)
        if self.keyfiletoaddtodatabase!=None:
            databaseparser.AddKeyFileParametersToParameterFile(self,m)   
            sys.exit()
        self.GrabProtStates(m)
        if self.genprotstatesonly==True:
            sys.exit()
        if ('I ' in self.mol.GetSpacedFormula()):
            if self.foundgauss==True:
                self.use_gaus=True
        if ('Br ' in self.mol.GetSpacedFormula()):
            self.torspbasisset=self.torspbasissethalogen
        self.pcm=False
        if pcm==True and self.dontusepcm==False:
            if self.foundgauss==True:
                self.use_gauPCM=True
                self.SanitizeAllQMMethods()
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
        self.idxtosymclass,self.symmetryclass=symm.gen_canonicallabels(self,mol) 
 
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
            

        if self.firstoptfinished==False:
            optmol,error,torsionrestraints = opt.GeometryOPTWrapper(self,mol)
            finished,error=self.CheckNormalTermination(self.firstlogoptfname)

            bondtopoopt=torgen.GenerateBondTopology(self,optmol)
            bondtopoopt=[list(i) for i in bondtopoopt]
            bondtopo=torgen.GenerateBondTopology(self,mol)
            bondtopo=[list(i) for i in bondtopo]
            for bond in bondtopo:
                if bond in bondtopoopt or bond[::-1] in bondtopoopt:
                    pass
                else:
                    if self.deletedfiles==True:
                        raise ValueError('Bond does not exist after optimization !'+str(bond))
                    else:
                        self.deletedfiles=True
                        self.DeleteAllFiles()
                        self.GenerateParameters()

            for bond in bondtopoopt:
                if bond in bondtopo or bond[::-1] in bondtopo:
                    pass
                else:
                    if self.deletedfiles==True:
                        raise ValueError('Bond created after optimization !'+str(bond))
                    else:
                        self.deletedfiles=True
                        self.DeleteAllFiles()
                        self.GenerateParameters()


            if finished==False:
                bondtoposame=self.CheckBondTopology(self.firstlogoptfname,self.rdkitmol)
            else:
                bondtoposame=True
            attempts=0
            maxiter=4
            cartxyz=self.firstlogoptfname.replace('.log','.xyz')
            opt.GrabFinalXYZStructure(self,self.firstlogoptfname,cartxyz,mol)

            inioptmol = opt.load_structfile(self,cartxyz)
            inioptmol.SetTotalCharge(mol.GetTotalCharge())

            while bondtoposame==False:
                if attempts>=maxiter or finished==True:
                    break
                try:           
                    optmol,error = opt.GeometryOptimization(self,inioptmol,loose=False,checkbonds=True,modred=True,bondanglerestraints=None,skipscferror=False,charge=None,skiperrors=True,overridecheckterm=True)
                    finished,error=self.CheckNormalTermination(self.firstlogoptfname)
                    cartxyz=self.firstlogoptfname.replace('.log','.xyz')
                    opt.GrabFinalXYZStructure(self,self.firstlogoptfname,cartxyz,mol)
                    inioptmol = opt.load_structfile(self,cartxyz)
                    inioptmol.SetTotalCharge(mol.GetTotalCharge())
                except:
                    pass
                bondtoposame=self.CheckBondTopology(self.firstlogoptfname,self.rdkitmol)
                attempts+=1
        else:
            cartxyz=self.firstlogoptfname.replace('.log','.xyz')
            optmol = opt.load_structfile(self,cartxyz)
            optmol.SetTotalCharge(mol.GetTotalCharge())

        optatomnums=optmol.NumAtoms()
        molatomnums=mol.NumAtoms()
        if optatomnums!=molatomnums: # then program behaviour changed need to restart
            self.deletedfiles=True
            self.DeleteAllFiles()
            self.GenerateParameters()

        checkmp2optfailed=self.CheckMP2OptFailed()
        if checkmp2optfailed==True:
            self.deletedfiles=True
            self.DeleteAllFiles()
            self.GenerateParameters()


        if self.optonly==True:
            sys.exit()
        if self.use_gausgeomoptonly==True:
            self.use_gausoptonly=False
            self.use_gaus=False
        if not os.path.isfile(self.key4fname) or not os.path.isfile(self.torsionsmissingfilename) or not os.path.isfile(self.torsionprmguessfilename):
            self.WriteToLog('Searching database for parameters')
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
        try:
            esp.SPForDMA(self,optmol,mol)
        except:
            if self.use_gaus==False: 
                self.use_gaus=True
                esp.SPForDMA(self,optmol,mol) 
                self.use_gaus=False
            else:
                traceback.print_exc(file=sys.stdout)
                sys.exit()

        # Obtain multipoles from Gaussian fchk file using GDMA
    
        if not os.path.isfile(self.gdmafname):
            mpole.run_gdma(self)
    
        # Set up input file for poledit
        # find multipole local frame definitions

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
        self.issane=self.CheckFileSanity()
        if self.issane==False:
            self.deletedfiles=True
            self.DeleteAllFiles()
            self.GenerateParameters()

            
        # post process local frames written out by poledit
        if self.atomnum!=1: 
             try:
                 esp.SPForESP(self,optmol,mol) 
             except:
                 if self.use_gaus==False: 
                     self.use_gaus=True
                     esp.SPForESP(self,optmol,mol) 
                     self.use_gaus=False
                 else:
                     traceback.print_exc(file=sys.stdout)
                     sys.exit()

        # End here if qm calculations were all that needed to be done 
        if self.qmonly:
            self.WriteToLog("poltype QM-only complete.")
            sys.exit(0)
    
               
        
        # generate the electrostatic potential grid used for multipole fitting
        if self.atomnum!=1: 
            esp.gen_esp_grid(self,optmol)
    
        # Average multipoles based on molecular symmetry
        # Does this using the script avgmpoles.pl which is found in the poltype directory
        # Atoms that belong to the same symm class will now have only one common multipole definition
        if not os.path.isfile(self.key2fnamefromavg):
            mpole.AverageMultipoles(self,optmol)
            mpole.AddPolarizeCommentsToKey(self,self.key2fnamefromavg,polartypetotransferinfo)
        if self.espfit and not os.path.isfile(self.key3fname) and self.atomnum!=1:
            # Optimize multipole parameters to QM ESP Grid (*.cube_2)
            # tinker's potential utility is called, with option 6.
            # option 6 reads: 'Fit Electrostatic Parameters to a Target Grid'
            
            esp.ElectrostaticPotentialFitting(self) 
            shutil.copy(self.key3fnamefrompot,self.key3fname)
        elif self.atomnum==1 or self.espfit==False:
            shutil.copy(self.key2fnamefromavg, self.key3fname)
        # Remove header terms from the keyfile
        mpole.rm_esp_terms_keyfile(self,self.key3fname)
        if self.atomnum!=1: 
            esp.ElectrostaticPotentialComparison(self) 
        
        if not os.path.exists(self.key4fname):
            databaseparser.appendtofile(self,self.key3fname,self.key4fname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo)
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
        if self.isfragjob==False and not os.path.isfile(self.key7fname) and self.dontfrag==False and (self.dontdotor==False):
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
        if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key7fname) and (self.dontdotor==False):
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
        try:
            opt.StructureMinimization(self,torsionrestraints)
            if self.atomnum != 1:
                opt.gen_superposeinfile(self)
                opt.CheckRMSD(self)
        except: # in case old key_4,key_5 files not working delete and restart
                #self.DeleteFilesWithExtension(['key_4','key_5'])
                #self.GenerateParameters()
                pass

        if self.torsppcm:
            torgen.RemoveStringFromKeyfile(self,self.key7fname,'solvate GK')
        if self.atomnum!=1: 
             esp.CheckDipoleMoments(self,optmol)
        self.FinalVDWMultipoleCheck(self.tmpkeyfile)
        self.WriteToLog('Poltype Job Finished'+'\n')
        
        if os.path.exists(self.scrtmpdirgau):
            shutil.rmtree(self.scrtmpdirgau)
        if os.path.exists(self.scrtmpdirpsi4):
            shutil.rmtree(self.scrtmpdirpsi4)
        string='Poltype has completed successfully, but this software is still under active development. It is your responsibility to check your own final parameters, while we are still in development phase.'
        warnings.warn(string)
        self.WriteToLog(string)
        self.CopyFitPlots()
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
        if not os.path.isdir(fold):
            os.mkdir(fold)
        for path in plots:
            shutil.copy(path,fold)
        fitlines=self.CollectElectrostaticDipoleFitLines()
        instructions=self.Instructions()
        os.chdir(fold)
        temp=open("README.txt",'w')
        for line in instructions:
            temp.write(line)
        for line in fitlines:
            temp.write(line)
        temp.close()

if __name__ == '__main__':
    def RunPoltype():
        poltype=PolarizableTyper() 
        try:
            poltype.main()
        except:
            traceback.print_exc(file=sys.stdout)
            text = str(traceback.format_exc())
            #if os.path.exists(poltype.scrtmpdirgau):
            #    shutil.rmtree(poltype.scrtmpdirgau)
            #if os.path.exists(poltype.scrtmpdirpsi4):
            #    shutil.rmtree(poltype.scrtmpdirpsi4)
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
