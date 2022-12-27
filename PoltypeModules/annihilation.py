import os
import sys
import time
import shutil
import numpy as np
import boxsetup as box
import keyfilemodifications as keymods
import minimization as mini     
import bar
import restraints as res
import plots
import tables
import productiondynamics as prod
import equilbriation as equil


def main(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for i in range(len(poltype.keyfilename)):
        keylist=poltype.keyfilename[i]
        configkeylist=poltype.configkeyfilename[i]
        for k in range(len(keylist)):
            firstkey=keylist[k]
            secondkey=configkeylist[k]
            shutil.copy(firstkey,secondkey)
    for keylist in poltype.configkeyfilename:
        for key in keylist:
            keymods.RemoveKeyWords(poltype,key,['parameters','TARGET-DIPOLE','OPENMP-THREADS'])
            keymods.InsertKeyfileHeader(poltype,key)    

    box.BoxSetupProtocol(poltype)
    if poltype.generateinputfilesonly==True:
        sys.exit()
    mini.CheapMinimizationProtocol(poltype)
    if poltype.estimatedynamictime==True:
        poltype.EstimateDynamicTime()
    poltype.ReportETA()
    if poltype.estimatedynamictimeonly==True:
        sys.exit()
    if poltype.equiltimeNVT!=0 and poltype.equiltimeNPT!=0:
        equil.EquilibriationProtocol(poltype)
    if poltype.proddyntime!=0:
        prod.ProductionDynamicsProtocol(poltype) 
        bar.BARProtocol(poltype)  
        tables.GenerateSimInfoTable(poltype)
        plots.PlotEnergyData(poltype)

