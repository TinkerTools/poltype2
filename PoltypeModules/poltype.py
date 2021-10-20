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


class PolarizableTyper():
    def __init__(self,fragmentjobslocal=False,toroptdebugmode=False,debugmode=False,fragmenterdebugmode=False,jobsatsametime=1,usepoleditframes=False,databasematchonly=False,setupfragjobsonly=False,allowradicals=False,checkinputonly=False,username=None,esprestweight=1,espgrad=.1,issane=True,deletedfiles=False,onlyfittorstogether=[],parentname=None,addhydrogentononcharged=True,accuratevdwsp=False,inputmoleculefolderpaths=None,email=None,firstoptfinished=False,optonly=False,onlyvdwatomindex=None,use_qmopt_vdw=False,use_gau_vdw=False,dontusepcm=False,deleteallnonqmfiles=False,totalcharge=None,torspbasissethalogen="6-311G*",homodimers=False,tortormissingfilename='tortormissing.txt',tordebugmode=False,amoebapluscfsmartstocommentmap=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebapluscfsmartstocomment.txt',amoebapluscfprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'cfprmlib.txt',amoebaplusnonbondedprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebaplusnonbonded.prm',amoebaplusnonbondedsmartstocommentmap=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoebaplusnonbonded.txt',smartstosoluteradiimap=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'SMARTsToSoluteRadiiMap.txt',latestsmallmoleculepolarizeprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarize.prm',latestsmallmoleculesmartstotypespolarize=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21polarcommenttoparameters.txt',latestsmallmoleculesmartstotinkerclass=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21smartstoclass.txt',latestsmallmoleculeprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba21.prm',boltzmantemp=8,dontdovdwscan=False,vdwprobepathname=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/VdwProbes/',vdwprobenames=['water'],use_gausgeomoptonly=False,maxtorRMSPDRel=.2,vdwmissingfilename='missingvdw.txt',databaseprmfilename='database.prm',tortor=False,torfit2Drotonly=False,torfit1Drotonly=False,externalparameterdatabase=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'externalparameterdatabase.txt',fitfirsttorsionfoldphase=False,keyfiletoaddtodatabase=None,skipgridsearch=True,torsionprmguessfilename='torsionprmguess.txt',defaultmaxtorsiongridpoints=40,torsionsmissingfilename='torsionsmissing.txt',smallmoleculemm3prmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'mm3.prm',smallmoleculesmartstomm3descrip=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstomm3typedescrip.txt',absdipoletol=.5,transferanyhydrogentor=True,smallmoleculesmartstotinkerdescrip=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'smartstoamoebatypedescrip.txt',smallmoleculeprmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba09.prm',torspbasissetfile='6-311+g_st_.0.gbs',toroptbasissetfile='6-311g_st_.0.gbs',optbasissetfile='6-311g_st_.0.gbs',dmabasissetfile='6-311g_st__st_.0.gbs',espbasissetfile='aug-cc-pvtz.1.gbs',iodinetorspbasissetfile='def2-svp.1.gbs',iodinetoroptbasissetfile='def2-svp.1.gbs',iodineoptbasissetfile='def2-svp.1.gbs',iodinedmabasissetfile='def2-svp.1.gbs',iodineespbasissetfile='def2-tzvpp.1.gbs',basissetpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/'+'BasisSets/',refinenonaroringtors=False,maxgrowthcycles=4,use_gauPCM=False,fitqmdipole=False,scfmaxiter=500,suppresstorfiterr=False,obminimizeexe='obminimize',readinionly=False,suppressdipoleerr=False,topologylib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ModifiedResidueLibraries/residue_connect.txt",poltypepath=os.path.abspath(os.path.split(__file__)[0]),WBOtol=.05,dontfrag=False,isfragjob=False,dipoletol=.5,externalapi=None,printoutput=False,poltypeini=True,structure=None,prmstartidx=401,numproc="1",maxmem="700MB",maxdisk="100GB",gausdir=None,gdmadir=None,tinkerdir=None,scratchdir="/scratch",paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18_header.prm",gausexe=None,formchkexe='formchk',cubegenexe='cubegen',gdmaexe='gdma',avgmpolesexe=os.path.abspath(os.path.join(os.path.abspath(os.path.join(__file__, os.pardir)), os.pardir)) + "/PoltypeModules/avgmpoles.pl",peditexe='poledit.x',potentialexe='potential.x',minimizeexe='minimize.x',analyzeexe='analyze.x',superposeexe='superpose.x',defopbendval=0.20016677990819662,Hartree2kcal_mol=627.5095,optbasisset='6-31G*',toroptbasisset='6-311G',dmabasisset='6-311G**',espbasisset="aug-cc-pVTZ",torspbasisset="6-311+G*",optmethod='MP2',toroptmethod='wB97X-D',torspmethod='wB97X-D',dmamethod='MP2',espmethod='MP2',qmonly = False,espfit = True,parmtors = True,foldnum=3,foldoffsetlist = [ 0.0, 180.0, 0.0, 180.0, 0.0, 180.0 ],torlist = None,rotbndlist = None,maxRMSD=1,maxRMSPD=1,maxtorRMSPD=1.8,tordatapointsnum=None,gentorsion=False,gaustorerror=False,torsionrestraint=.1*3282.80354574,onlyrotbndslist=None,rotalltors=False,dontdotor=False,dontdotorfit=False,toroptpcm=False,optpcm=False,torsppcm=False,use_gaus=False,use_gausoptonly=False,freq=False,postfit=False,bashrcpath=None,amoebabioprmpath=None,libpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ModifiedResidueLibraries/lib.bio18_conv1.txt",SMARTSToTypelibpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ModifiedResidueLibraries/SMARTSToTypeLib.txt',ModifiedResiduePrmPath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/ModifiedResidue.prm',modifiedproteinpdbname=None,unmodifiedproteinpdbname=None,mutatedsidechain=None,mutatedresiduenumber=None,modifiedresiduepdbcode=None,optmaxcycle=400,torkeyfname=None,gausoptcoords='',forcefield="AMOEBA",helpfile='README.md',versionfile='version.md',sleeptime=.1):
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
        self.dontdovdwscan=dontdovdwscan
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

                    elif "usepoleditframes" in newline:
                        if '=' not in line:
                            self.usepoleditframes = True
                        else:
                            self.usepoleditframes=self.GrabBoolValue(a)
                    elif "fragmenterdebugmode" in newline:
                        if '=' not in line:
                            self.fragmenterdebugmode = True
                        else:
                            self.fragmenterdebugmode=self.GrabBoolValue(a)
                    elif "fragmentjobslocal" in newline:
                        if '=' not in line:
                            self.fragmentjobslocal = True
                        else:
                            self.fragmentjobslocal=self.GrabBoolValue(a)




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
                    elif "jobsatsametime" in newline:
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

                    elif "optonly" in newline:
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


                    elif 'poltypepath' in newline:
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

                    elif "dontdovdwscan" in newline:
                        if '=' not in line:
                            self.dontdovdwscan = True
                        else:
                            self.dontdovdwscan=self.GrabBoolValue(a)


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
                        self.sleeptime = int(a)
                    elif "help" in newline:
                        self.copyright()
                        self.usage()
                        sys.exit(2)
                    else:
                        print('Unrecognized '+line)
                        self.usage()
                        print('Unrecognized '+line)
                        sys.exit()


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



        self.SanitizeAllQMMethods()
        if self.readinionly==True:
            return
        self.SanitizeMMExecutables()
        self.copyright()
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
            if self.forcefield.upper() != "AMOEBA+":  #allow old version for AMOEBA+
                raise ValueError("Notice: Not latest working version of tinker (8.9.4)"+' '+os.getcwd())
      
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
        self.keyfname = self.assign_filenames ( "keyfname" , ".key")
        self.xyzfname = self.assign_filenames ( "xyzfname" , ".xyz")
        self.peditinfile = self.assign_filenames ( "peditinfile" , "-peditin.txt")
        self.superposeinfile = self.assign_filenames ( "superposeinfile" , "-superin.txt")
        self.espgrdfname = self.assign_filenames ( "espgrdfname" , ".grid")
        self.qmespfname = self.assign_filenames ( "qmespfname" , "_fortinker.cube")
        self.qmesp2fname = self.assign_filenames ( "qmesp2fname" , "_fortinker.pot")
        self.grpfname = self.assign_filenames ( "grpfname" , "-groups.txt")
        self.key2fname = self.assign_filenames ( "key2fname" , ".key_2")
        self.key3fname = self.assign_filenames ( "key3fname" , ".key_3")
        self.key4fname = self.assign_filenames ( "key4fname" , ".key_4")
        self.key5fname = self.assign_filenames ( "key5fname" , ".key_5")
        self.key6fname= self.assign_filenames ( "key6fname" , ".key_6")
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

    def Chunks(self,lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    
    def CallJobsSeriallyLocalHost(self,fulljobtooutputlog,skiperrors):
       jobchunks=list(self.Chunks(list(fulljobtooutputlog.keys()),self.jobsatsametime))   
       thepath=os.path.join(os.getcwd(),'Fragments')
       for jobidx in range(len(jobchunks)):
           jobs=jobchunks[jobidx]
           count=jobidx+len(jobs)
           ratio=(100*count)/len(fulljobtooutputlog.keys())
           self.WriteToLog('Percent of jobs submitted '+str(ratio))   
           if jobidx!=0:
               if 'poltype.py' in job:
                   self.ETAQMFinish(thepath,len(fulljobtooutputlog.keys()))
           self.call_subsystem(jobs,False,skiperrors)
           temp={}
           for job in jobs:
               temp[job]=fulljobtooutputlog[job]
           finishedjob,errorjob=self.WaitForTermination(temp,skiperrors)
       finishedjobs,errorjobs=self.WaitForTermination(fulljobtooutputlog,skiperrors)
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
                if os.path.isfile(outputlog):
                    statinfo=os.stat(outputlog)
                    size=statinfo.st_size
                    if size==0:
                        continue
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
        if os.path.isfile(logfname):
            head,tail=os.path.split(logfname)
            Ftime=os.path.getmtime(logfname)
            reltime=time.time()-Ftime
            htime=reltime*0.000277778
            updatetime=4 # hours. sometimes psi4 gives return code 0 even though program crashes
            foundendgau=False # sometimes gaussian comments have keyword Error, ERROR in them
            for line in open(logfname):
                if 'poltype' in tail:
                    if 'Poltype Job Finished' in line:
                        term=True
                else:
                    if "Final optimized geometry" in line or "Electrostatic potential computed" in line or 'Psi4 exiting successfully' in line or "LBFGS  --  Normal Termination due to SmallGrad" in line or "Normal termination" in line or 'Normal Termination' in line or 'Total Potential Energy' in line:
                        term=True
                    if ('error' in line or 'Error' in line or 'ERROR' in line or 'impossible' in line or 'software termination' in line or 'segmentation violation, address not mapped to object' in line or 'galloc:  could not allocate memory' in line or 'Erroneous write.' in line) and 'DIIS' not in line and 'mpi' not in line:
                        error=True
                        errorline=line
                    if 'segmentation violation' in line and 'address not mapped to object' not in line or 'Waiting' in line or ('OptimizationConvergenceError' in line and 'except' in line) or "Error on total polarization charges" in line or 'Erroneous write' in line:
                        error=False
                        continue
                    if ('Error termination request processed by link 9999' in line or 'Error termination via Lnk1e in' in line) or ('OptimizationConvergenceError' in line and 'except' not in line) or 'Could not converge geometry optimization' in line:
                        error=True
                        errorline=line
                    if 'l9999.exe' in line and foundendgau==False:
                        foundendgau=True
                        preverror=error # did find error before comment?
                    if foundendgau==True:
                        if error==True and preverror==False: # then caused by Gaussian comment
                            error=False
                    

            if error==True:
                term=False # sometimes psi4 geometry opt not fully converge but says successfully exiting etc..
                message='Error '+errorline+ 'logpath='+logfname
            #if error==False and term==False and htime>=updatetime:
            #    error=True
            #    message='Error '+'Job has not been updated in '+str(updatetime)+' hours'+' last update time = '+str(htime)+' hours'+' logname='+logfname
            if error==True and term==False and skiperrors==False:
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
            if 'torsion' in line and '#' not in line and 'Missing' not in line:
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

    def CheckMemorySettings(self):
        proc=int(self.numproc)
        if proc>8:
            raise ValueError('Too many input processors, lower the numproc value to 8 or less')
           
    
    def CheckInputCharge(self,molecule):
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
            if atomnum==6 and val==3 and self.addhydrogentononcharged==True and radicals==0:
                warnings.warn('WARNING! Strange valence for Carbon, will assume missing hydrogens and add'+string) 
                self.WriteToLog('WARNING! Strange valence for Carbon, will assume missing hydrogens and add '+string)
                atom.SetNumRadicalElectrons(0)
                chg=0

            elif atomnum==7 and val==2 and self.addhydrogentononcharged==True and radicals==0:
                warnings.warn('WARNING! Strange valence for Nitrogen, will assume missing hydrogens and add'+string) 
                self.WriteToLog('WARNING! Strange valence for Nitrogen, will assume missing hydrogens and add '+string)
                atom.SetNumRadicalElectrons(0)
                chg=0

            elif atomnum==7 and val==2 and radicals==1:
                warnings.warn('WARNING! Strange valence for Nitrogen, will assume radical and set charge to zero') 
                self.WriteToLog('WARNING! Strange valence for Nitrogen, will assume radical and set charge to zero')
                self.allowradicals=True

                atom.SetFormalCharge(0)
                self.addhydrogentononcharged=False

            elif atomnum==8 and val==2 and radicals==1:
                warnings.warn('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1') 
                self.WriteToLog('WARNING! Strange valence for Oxygen, will assume radical and set charge to +1')
                self.allowradicals=True

                atom.SetFormalCharge(1)
            elif atomnum==8 and val==1 and radicals==1:
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
            if not os.path.isdir(f) and 'nohup' not in f:
                fsplit=f.split('.')
                if len(fsplit)>1:
                    end=fsplit[1]
                    if 'log' not in end and 'sdf' not in end and 'ini' not in end and 'chk' not in end and 'dat' not in end: 
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
        dic={'username':self.username,'externalapi':self.externalapi,'accuratevdwsp':self.accuratevdwsp,'email':self.email,'firstoptfinished':self.firstoptfinished,'optonly':self.optonly,'onlyvdwatomindex':self.onlyvdwatomindex,'use_qmopt_vdw':self.use_qmopt_vdw,'use_gau_vdw':self.use_gau_vdw,'dontusepcm':self.dontusepcm,'deleteallnonqmfiles':self.deleteallnonqmfiles,'totalcharge':self.totalcharge,'torspbasissethalogen':self.torspbasissethalogen,'homodimers':self.homodimers,'boltzmantemp':self.boltzmantemp,'dontdovdwscan':self.dontdovdwscan,'use_gausgeomoptonly':self.use_gausgeomoptonly,'maxtorRMSPDRel':self.maxtorRMSPDRel,'tortor':self.tortor,'fitfirsttorsionfoldphase':self.fitfirsttorsionfoldphase,'defaultmaxtorsiongridpoints':self.defaultmaxtorsiongridpoints,'absdipoletol':self.absdipoletol,'refinenonaroringtors':self.refinenonaroringtors,'maxgrowthcycles':self.maxgrowthcycles,'use_gauPCM':self.use_gauPCM,'fitqmdipole':self.fitqmdipole,'WBOtol':self.WBOtol,'dontfrag':self.dontfrag,'dipoletol':self.dipoletol,'numproc':self.numproc,'maxmem':self.maxmem,'maxdisk':self.maxdisk,'optbasisset':self.optbasisset,'toroptbasisset':self.toroptbasisset,'dmabasisset':self.dmabasisset,'espbasisset':self.espbasisset,'torspbasisset':self.torspbasisset,'optmethod':self.optmethod,'toroptmethod':self.toroptmethod,'torspmethod':self.torspmethod,'dmamethod':self.dmamethod,'espmethod':self.espmethod,'qmonly' : self.qmonly,'espfit' : self.espfit,'foldnum':self.foldnum,'maxRMSD':self.maxRMSD,'maxRMSPD':self.maxRMSPD,'maxtorRMSPD':self.maxtorRMSPD,'tordatapointsnum':self.tordatapointsnum,'torsionrestraint':self.torsionrestraint,'rotalltors':self.rotalltors,'dontdotor':self.dontdotor,'dontdotorfit':self.dontdotorfit,'toroptpcm':self.toroptpcm,'optpcm':self.optpcm,'torsppcm':self.torsppcm,'use_gaus':self.use_gaus,'use_gausoptonly':self.use_gausoptonly,'freq':self.freq,'optmaxcycle':self.optmaxcycle,'forcefield':self.forcefield}
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




    def GenerateParameters(self):
        if self.deleteallnonqmfiles==True:
            self.DeleteAllNonQMFiles()
        self.CheckMemorySettings()       
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
        # Begin log. *-poltype.log
        self.logfh = open(self.logfname,"w",buffering=1)

        obConversion.SetOutFormat('mol')
        self.molstructfnamemol=self.molstructfname.replace('.sdf','.mol')
        obConversion.WriteFile(mol,self.molstructfnamemol)
        indextocoordinates=self.GrabIndexToCoordinates(mol)
        m=Chem.MolFromMolFile(self.molstructfnamemol,removeHs=False,sanitize=False)
        m,atomindextoformalcharge=self.CheckInputCharge(m)
        if self.allowradicals==True:
            self.dontfrag=True # Psi4 doesnt allow UHF and properties (like compute WBO) for fragmenter, so need to turn of fragmenter if radical detected
        m.UpdatePropertyCache()
        if self.addhydrogentononcharged==True:
            m = Chem.AddHs(m)
            AllChem.EmbedMolecule(m)
        Chem.SanitizeMol(m)

        pcm=self.CheckForConcentratedFormalCharges(m,atomindextoformalcharge)
        cpm = copy.deepcopy(m)
        if self.firstoptfinished==False:
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
            RunPoltype()

        if self.optonly==True:
            sys.exit()
        if self.use_gausgeomoptonly==True:
            self.use_gausoptonly=False
            self.use_gaus=False

        if not os.path.isfile(self.key4fname) or not os.path.isfile(self.torsionsmissingfilename) or not os.path.isfile(self.torsionprmguessfilename):
            bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,torsionsmissing,classkeytotorsionparametersguess,missingvdwatomindextoneighbors,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo,tortorsmissing=databaseparser.GrabSmallMoleculeAMOEBAParameters(self,optmol,mol,m)
        if os.path.isfile(self.torsionsmissingfilename):
            torsionsmissing=databaseparser.ReadTorsionList(self,self.torsionsmissingfilename)
        if os.path.isfile(self.torsionprmguessfilename):
            classkeytotorsionparametersguess=databaseparser.ReadDictionaryFromFile(self,self.torsionprmguessfilename)
        if os.path.isfile(self.vdwmissingfilename):
            missingvdwatomindices=databaseparser.ReadVdwList(self,self.vdwmissingfilename)

        if os.path.isfile(self.tortormissingfilename):
            tortorsmissing=databaseparser.ReadTorTorList(self,self.tortormissingfilename)
        try:
            esp.SPForDMA(self,optmol,mol)
        except:
            if self.use_gaus==True: # if gaussian failed try psi4
                self.use_gaus=False
                esp.SPForDMA(self,optmol,mol)
                self.use_gaus=True
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
            while not os.path.isfile(self.keyfname):
                time.sleep(1)
                self.WriteToLog('Waiting for '+self.keyfname)
                
            mpole.prepend_keyfile(self,self.keyfname,optmol)
            mpole.SanitizeMultipoleFrames(self,self.keyfname)
            

        self.issane=self.CheckFileSanity()
        if self.issane==False:
            self.deletedfiles=True
            self.DeleteAllFiles()
            RunPoltype()
            
        # post process local frames written out by poledit
        if self.atomnum!=1: 
             try:
                 esp.SPForESP(self,optmol,mol) 
             except:
                 if self.use_gaus==True: # if gaussian failed try psi4
                     self.use_gaus=False
                     esp.SPForESP(self,optmol,mol) 
                     self.use_gaus=True
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
        if not os.path.isfile(self.key2fname):
            mpole.AverageMultipoles(self,optmol)
            mpole.AddPolarizeCommentsToKey(self,self.key2fname,polartypetotransferinfo)
        if self.espfit and not os.path.isfile(self.key3fname) and self.atomnum!=1:
            # Optimize multipole parameters to QM ESP Grid (*.cube_2)
            # tinker's potential utility is called, with option 6.
            # option 6 reads: 'Fit Electrostatic Parameters to a Target Grid'
            
            esp.ElectrostaticPotentialFitting(self) 
        elif self.atomnum==1 or self.espfit==False:
            shutil.copy(self.key2fname, self.key3fname)
        # Remove header terms from the keyfile
        mpole.rm_esp_terms_keyfile(self,self.key3fname)
        if self.atomnum!=1: 
            esp.ElectrostaticPotentialComparison(self) 
            #if self.failedrmspd==True and self.deletedfiles==False:
            #    self.DeleteFilesWithExtension(['pot','grid','key','xyz','key_2','key_3','key_4','key_5','xyz_2','cube'])
            #    self.DeleteFilesWithString(['esp','dma'])
            #    self.deletedfiles=True
            #    self.GenerateParameters()

        
        
        # Now that multipoles have been found
        # Other parameters such as opbend, vdw, etc. are found here using a look up table
        # Part of the look up table is here in poltype.py 
        # Most of it is in the file databaseparser.py found in the poltype directory
    
        # Finds aromatic carbons and associated hydrogens and corrects polarizability
        # Find opbend values using a look up table
        # Outputs a list of rotatable bonds (found in get_torlist) in a form usable by databaseparser.py
    
        # Map from idx to symm class is made for databaseparser.py
        # databaseparser.py method is called to find parameters and append them to the keyfile
        if not os.path.exists(self.key4fname):
            databaseparser.appendtofile(self,self.key3fname,self.key4fname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo)
            databaseparser.TestBondAngleEquilValues(self)
            if self.databasematchonly==True:
                sys.exit()

        if self.torsppcm:
            torgen.PrependStringToKeyfile(self,self.key4fname,'solvate GK')
        torgen.get_all_torsions(self,mol)
        # Find rotatable bonds for future torsion scans
        (self.torlist, self.rotbndlist,hydtorsions,nonaroringtorlist) = torgen.get_torlist(self,mol,torsionsmissing)
        if atomnum<25 and len(nonaroringtorlist)==0: 
            self.dontfrag=True
        self.torlist,self.rotbndlist=torgen.RemoveDuplicateRotatableBondTypes(self) # this only happens in very symmetrical molecules
        self.torlist=[tuple(i) for i in self.torlist]
        self.torlist=[tuple([i]) for i in self.torlist]
        self.torsettovariabletorlist={}
        for torset in self.torlist:
            self.torsettovariabletorlist[tuple(torset)]=[]
        nonaroringtorlist=[tuple(i) for i in nonaroringtorlist]
        nonaroringtorlist=[tuple([i]) for i in nonaroringtorlist]
        self.rotbndtoanginc=torgen.DetermineAngleIncrementAndPointsNeededForEachTorsionSet(self,mol,self.rotbndlist)
        torgen.DefaultMaxRange(self,self.torlist)
        if self.dontdotor==True:
            self.torlist=[]

        # add missingvdwindices to torlist (for fragmenter input)
        missingvdwatomsets=[]
        if self.isfragjob==False and self.dontfrag==False and self.dontdovdwscan==False:
            for vdwatomindex in missingvdwatomindices:
                ls=tuple([tuple([vdwatomindex])])
                missingvdwatomsets.append(ls)
                self.torlist.append(ls)
        if self.dontdotor==True and self.dontdovdwscan==True:
            shutil.copy(self.key4fname,self.key5fname)
        self.torsettofilenametorset={}
        self.torsettotortorindex={}
        self.torsettotortorphaseindicestokeep={}
        self.nonaroringtors=[]
        self.nonaroringtorsets=[]
        self.classkeytoinitialprmguess={}
        if self.tortor==True and self.dontdotor==False:
            torgen.PrepareTorsionTorsion(self,optmol,mol,tortorsmissing)
        self.nonarotortotorsbeingfit={}
        if self.refinenonaroringtors==True and self.dontfrag==False:
            rings.RefineNonAromaticRingTorsions(self,mol,optmol,classkeytotorsionparametersguess)


        if self.isfragjob==False and not os.path.isfile(self.key5fname) and self.dontfrag==False and (self.dontdotor==False or self.dontdovdwscan==False):

            WBOmatrix,outputname,error=frag.GenerateWBOMatrix(self,self.rdkitmol,self.mol,self.logoptfname.replace('.log','.xyz'))
            highlightbonds=[]
            for torset in self.torlist:
                for tor in torset:
                    if len(tor)>1:
                        rotbnd=[tor[1]-1,tor[2]-1]
                        highlightbonds.append(rotbnd)
            frag.Draw2DMoleculeWithWBO(self,WBOmatrix,self.molstructfname.replace('.sdf',''),self.rdkitmol,bondindexlist=highlightbonds,imgsize=1500)        
            rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor=frag.GenerateFragments(self,self.mol,self.torlist,WBOmatrix,missingvdwatomsets,nonaroringtorlist) # returns list of bond indexes that need parent molecule to do torsion scan for (fragment generated was same as the parent0
            equivalentrotbndindexarrays,rotbndindextoringtor=frag.SpawnPoltypeJobsForFragments(self,rotbndindextoparentindextofragindex,rotbndindextofragment,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor)
        if self.dontfrag==False and self.isfragjob==False and not os.path.isfile(self.key5fname) and (self.dontdotor==False or self.dontdovdwscan==False):
            frag.GrabVdwAndTorsionParametersFromFragments(self,rotbndindextofragmentfilepath,equivalentrotbndindexarrays,rotbndindextoringtor) # just dump to key_5 since does not exist for parent molecule
        else:         
            # Torsion scanning then fitting. *.key_5 will contain updated torsions
            if not os.path.isfile(self.key5fname):
                if len(self.torlist)!=0:
                    # torsion scanning
                    torgen.gen_torsion(self,optmol,self.torsionrestraint,mol)
                    # torsion fitting
                    if self.dontdotorfit==True:
                        shutil.copy(self.key4fname,self.key5fname)
                        sys.exit()
                    torfit.process_rot_bond_tors(self,optmol)
                else:
                    shutil.copy(self.key4fname,self.key5fname)
        if self.isfragjob and len(self.onlyrotbndslist)!=0:
            self.dontdovdwscan=True
        if self.dontdovdwscan==False:
            if self.dontfrag==False: 
                if self.isfragjob==True:
                    vdwfit.VanDerWaalsOptimization(self,missingvdwatomindices)    
            else:
                vdwfit.VanDerWaalsOptimization(self,missingvdwatomindices)       
        
        if self.isfragjob==False and self.dontdotor==False:
            self.CheckTorsionParameters(self.key5fname,torsionsmissing,hydtorsions)
        self.WriteOutLiteratureReferences(self.key5fname) 
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
            torgen.RemoveStringFromKeyfile(self,self.key5fname,'solvate GK')
        if self.atomnum!=1: 
             esp.CheckDipoleMoments(self,optmol)
        self.WriteToLog('Poltype Job Finished'+'\n')
        
        if os.path.exists(self.scrtmpdirgau):
            shutil.rmtree(self.scrtmpdirgau)
        if os.path.exists(self.scrtmpdirpsi4):
            shutil.rmtree(self.scrtmpdirpsi4)
        string='Poltype has completed successfully, but this software is still under active development. It is your responsibility to check your own final parameters, while we are still in development phase.'
        warnings.warn(string)
        self.WriteToLog(string)
        if self.email!=None:
            moleculename=self.molstructfname.replace('.sdf','')
            password='amoebaisbest'
            fromaddr = 'poltypecrashreportnoreply@gmail.com'
            toaddr = self.email
            TEXT='Molecule has finished parameterization'
            self.SendFinalReportEmail(TEXT,fromaddr,toaddr,password,moleculename)


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


    def DeleteAllFiles(self):
        files=os.listdir()
        for f in files:
            if os.path.isdir(f):
                shutil.rmtree(f)
            elif '.' in f and not os.path.isdir(f):
                ext=f.split('.')[1]
                if ext!='sdf' and ext!='ini' and ext!='out':
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
                poltype.SendCrashReportEmail(text,fromaddr,toaddr,password,filename)
            raise ValueError('Houston, we have a problem. Buy a developer some coffee!')
    RunPoltype()
